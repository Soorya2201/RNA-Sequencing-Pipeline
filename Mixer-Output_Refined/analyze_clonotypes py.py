
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
analyze_clonotypes.py — Merge MiXCR clonotype TSVs into a single CSV and generate visualizations.

Usage examples:
  python analyze_clonotypes.py --input-glob "/mnt/data/patient_1_mixcr_clonotypes_*.tsv" --outdir /mnt/data --sample patient_1
  python analyze_clonotypes.py --input-glob "./data/*clonotypes_*.tsv" --outdir ./out --sample S01 --topn 25

What it does:
- Reads MiXCR clonotype TSV(s), auto-detects columns (CDR3 AA, V/J, count/fraction).
- Normalizes V/J names (drop *alleles and (score) suffixes).
- Computes clone fractions (prefers cloneFraction; else derived from cloneCount).
- Writes combined CSV and diversity-by-locus CSV.
- Saves PNG plots: top clonotypes, V/J usage per locus, CDR3 length per locus, cumulative frequency curve (overall).

Notes:
- This script uses matplotlib only (no seaborn) and one chart per PNG, no custom colors.
- Supported loci auto-detected from filenames: TRA/TRB/TRD/TRG/IGH/IGK/IGL (fallback "UNK").
"""

import argparse
import glob
import os
import re
from pathlib import Path
from collections import Counter

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")  # headless-safe
import matplotlib.pyplot as plt


SUPPORTED_LOCI = {"TRA","TRB","TRD","TRG","IGH","IGK","IGL"}


def parse_args():
    ap = argparse.ArgumentParser(description="Merge MiXCR clonotype TSVs and visualize.")
    ap.add_argument("--input-glob", required=True, help="Glob for input TSV files (e.g., './data/*clonotypes_*.tsv')") 
    ap.add_argument("--outdir", required=True, help="Output directory for CSV and PNG files") 
    ap.add_argument("--sample", default="sample", help="Sample name to annotate rows (default: sample)")
    ap.add_argument("--topn", type=int, default=20, help="Top N clonotypes for the bar chart (default: 20)")
    return ap.parse_args()


# ------------------ Column helpers ------------------
def first_present(df, candidates, default=None):
    for c in candidates:
        if c in df.columns:
            return c
    return default

def parse_gene(gene_val):
    if pd.isna(gene_val):
        return None
    s = str(gene_val).strip()
    s = s.split(",")[0]          # first hit if multi
    s = s.split("(")[0]          # drop score
    s = s.split("*")[0]          # drop allele
    return s.strip() or None

def compute_fraction(df):
    """Pick cloneFraction if numeric, else derive from cloneCount, else uniform."""
    frac_col = first_present(df, ["cloneFraction", "clone.freq", "Fraction", "freq"]) 
    if frac_col is not None:
        try:
            return pd.to_numeric(df[frac_col], errors="coerce").fillna(0.0).clip(lower=0.0)
        except Exception:
            pass
    count_col = first_present(df, ["cloneCount", "clone.count", "Count", "count"]) 
    if count_col is not None:
        counts = pd.to_numeric(df[count_col], errors="coerce").fillna(0.0)
        total = counts.sum() if counts.sum() > 0 else 1.0
        return counts / total
    n = len(df)
    return pd.Series(np.full(n, 1.0 / max(n,1)), index=df.index)

def get_cdr3_aa(df):
    aa_col = first_present(df, ["aaCDR3", "CDR3.aa", "cdr3aa", "AA.JUNCTION", "aaSeqCDR3", "cdr3"]) 
    if aa_col is None:
        return pd.Series([None]*len(df))
    return df[aa_col].astype(str)

def get_v_col(df):
    return first_present(df, ["allVHitsWithScore", "bestVHit", "V.name", "v.segm", "V", "allVHits"]) 

def get_j_col(df):
    return first_present(df, ["allJHitsWithScore", "bestJHit", "J.name", "j.segm", "J", "allJHits"]) 

def compute_cdr3_len(aa_series):
    lens = []
    for x in aa_series:
        if pd.isna(x) or str(x) in ("None", "nan"):
            lens.append(np.nan)
        else:
            aa = str(x).replace('"',"").replace("'","").strip()
            lens.append(len(aa))
    return pd.Series(lens, index=aa_series.index)

def locus_from_path(p):
    s = Path(p).stem
    for tag in sorted(SUPPORTED_LOCI, key=len, reverse=True):
        if tag in s:
            return tag
    # Try regex fallback
    m = re.search(r"_(TR[ABDG]|IG[HKL])", s)
    return m.group(1) if m else "UNK"

def short_clone_label(row):
    cdr3 = row.get("cdr3_aa", None)
    v = row.get("V_gene", None)
    j = row.get("J_gene", None)
    parts = [p for p in [v, (None if (cdr3 in (None, "None", "nan")) else cdr3), j] if p]
    lab = "/".join(parts) if parts else "clonotype"
    if len(lab) > 40:
        lab = lab[:37] + "…"
    return lab


# ------------------ Diversity metrics ------------------
def shannon(p):
    p = np.asarray(p, dtype=float)
    p = p[p > 0]
    return float(-np.sum(p * np.log(p))) if p.size else float("nan") 

def simpson(p):
    p = np.asarray(p, dtype=float)
    return float(1.0 - np.sum(p**2)) if p.size else float("nan") 


# ------------------ Main ------------------
def main():
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    files = sorted(glob.glob(args.input_glob))
    if not files:
        raise SystemExit(f"No input files matched: {args.input_glob}")

    all_rows = []
    for fp in files:
        try:
            df = pd.read_csv(fp, sep="\t")
        except Exception:
            df = pd.read_csv(fp, sep="\t", engine="python")
        locus = locus_from_path(fp)
        aa = get_cdr3_aa(df)
        v_col = get_v_col(df)
        j_col = get_j_col(df)
        V_gene = df[v_col].map(parse_gene) if v_col else pd.Series([None]*len(df))
        J_gene = df[j_col].map(parse_gene) if j_col else pd.Series([None]*len(df))
        frac = compute_fraction(df)
        cdr3_len = compute_cdr3_len(aa)
        clone_id_col = first_present(df, ["cloneId", "clone.id", "ID", "index"]) 
        clone_id = df[clone_id_col] if clone_id_col else pd.Series(np.arange(len(df)))

        sub = pd.DataFrame({
            "sample": args.sample,
            "locus": locus,
            "clone_id": clone_id,
            "cdr3_aa": aa,
            "V_gene": V_gene,
            "J_gene": J_gene,
            "clone_fraction": frac.astype(float),
            "cdr3_aa_len": cdr3_len.astype(float),
            "source_file": os.path.basename(fp),
        })
        cnt_col = first_present(df, ["cloneCount", "clone.count", "Count", "count"]) 
        if cnt_col:
            sub["clone_count"] = pd.to_numeric(df[cnt_col], errors="coerce")
        all_rows.append(sub)

    combined = pd.concat(all_rows, ignore_index=True)

    # Exports
    combined_csv = outdir / "combined_clonotypes.csv"
    combined.to_csv(combined_csv, index=False)

    # Diversity by locus
    div = []
    for locus, g in combined.groupby("locus"):
        pf = g["clone_fraction"].fillna(0).values
        div.append({
            "locus": locus,
            "shannon": shannon(pf),
            "simpson": simpson(pf),
            "n_clonotypes": len(g),
        })
    div_df = pd.DataFrame(div)
    div_csv = outdir / "diversity_by_locus.csv"
    div_df.to_csv(div_csv, index=False)

    # ---------- Plots (one chart per PNG, default styles) ----------
    # 1) Top-N overall clonotypes
    topn = max(1, int(args.topn))
    top = combined.sort_values("clone_fraction", ascending=False).head(topn).copy()
    if len(top) > 0:
        top["label"] = top.apply(short_clone_label, axis=1)
        plt.figure(figsize=(10, 6))
        plt.barh(top["label"][::-1], top["clone_fraction"][::-1])
        plt.xlabel("Clone fraction")
        plt.ylabel(f"Top {topn} clonotypes")
        plt.title("Top clonotypes (overall)")
        plt.tight_layout()
        plt.savefig(outdir / "top_overall.png")
        plt.close()

    # 2) V usage per locus
    for locus, g in combined.groupby("locus"):
        vc = g["V_gene"].dropna()
        if len(vc) == 0: 
            continue
        counts = vc.value_counts().head(15)
        plt.figure(figsize=(10,6))
        plt.barh(counts.index[::-1], counts.values[::-1])
        plt.xlabel("Count"); plt.ylabel("V gene"); plt.title(f"V gene usage ({locus})")
        plt.tight_layout()
        plt.savefig(outdir / f"V_usage_{locus}.png"); plt.close()

    # 3) J usage per locus
    for locus, g in combined.groupby("locus"):
        jc = g["J_gene"].dropna()
        if len(jc) == 0: 
            continue
        counts = jc.value_counts().head(15)
        plt.figure(figsize=(10,6))
        plt.barh(counts.index[::-1], counts.values[::-1])
        plt.xlabel("Count"); plt.ylabel("J gene"); plt.title(f"J gene usage ({locus})")
        plt.tight_layout()
        plt.savefig(outdir / f"J_usage_{locus}.png"); plt.close()

    # 4) CDR3 length per locus
    for locus, g in combined.groupby("locus"):
        lengths = g["cdr3_aa_len"].dropna().values
        if lengths.size == 0: 
            continue
        plt.figure(figsize=(8,5))
        bins = range(int(np.nanmin(lengths)), int(np.nanmax(lengths)) + 2)
        plt.hist(lengths, bins=bins)
        plt.xlabel("CDR3 length (aa)"); plt.ylabel("Count"); plt.title(f"CDR3 length distribution ({locus})")
        plt.tight_layout()
        plt.savefig(outdir / f"cdr3_len_{locus}.png"); plt.close()

    # 5) Cumulative fraction overall
    freq = combined["clone_fraction"].fillna(0).sort_values(ascending=False).to_numpy()
    if freq.size > 0:
        cum = np.cumsum(freq)
        x = np.arange(1, len(cum)+1)
        plt.figure(figsize=(8,5))
        plt.plot(x, cum)
        plt.xlabel("Top N clonotypes"); plt.ylabel("Cumulative fraction"); plt.title("Cumulative clonotype frequency (overall)")
        plt.tight_layout()
        plt.savefig(outdir / "cumulative_fraction_overall.png"); plt.close()

    print(f"\nWrote: {combined_csv}") 
    print(f"Wrote: {div_csv}")
    print(f"PNG files in: {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
