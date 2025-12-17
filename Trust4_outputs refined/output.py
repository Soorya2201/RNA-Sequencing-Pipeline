#!/usr/bin/env python3
import argparse
import re
from pathlib import Path
from typing import Optional, Tuple
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(
        description="Produce two files: final.out parsed and assembled_reads parsed (header, sequence, length, count)."
    )
    p.add_argument("--final", type=Path, required=True, help="Path to *_final.out")
    p.add_argument("--assembled", type=Path, required=True, help="Path to *_assembled_reads.fa")
    p.add_argument("--out_prefix", type=Path, default=Path("Patient_1"),
                   help="Output prefix (default: Patient_1)")
    # If you want to keep the raw header (incl. trailing count token) in the assembled output, flip this flag.
    p.add_argument("--keep-raw-header", action="store_true",
                   help="Do NOT strip trailing numeric count from assembled header column")
    return p.parse_args()


# -------------------- helpers --------------------

def _is_numeric_matrix_line(s: str) -> bool:
    """True if line is just digits/space/dots (seen in TRUST4 .out coverage matrices)."""
    return bool(s) and re.fullmatch(r"[0-9\s\.]+", s) is not None


def _looks_like_dna(s: str) -> bool:
    """Heuristic: treat as sequence line if it contains DNA letters (A/C/G/T/N)."""
    return re.search(r"[ACGTNacgtn]", s) is not None


def _parse_count_from_header(header: str) -> Tuple[Optional[int], str]:
    """
    Extract the LAST integer token from the header as 'count'.
    Returns (count, header_without_that_token).
    If none found, returns (None, original_header).
    """
    tokens = header.strip().split()
    # Find index of last all-digit token
    idx = None
    for i in range(len(tokens) - 1, -1, -1):
        if tokens[i].isdigit():
            idx = i
            break
    if idx is None:
        return None, header
    count = int(tokens[idx])
    cleaned = " ".join(tokens[:idx] + tokens[idx+1:])
    # Trim extra whitespace
    cleaned = cleaned.strip()
    return count, cleaned if cleaned else header


# -------------------- parsers --------------------

def parse_fasta_multi(path: Path) -> pd.DataFrame:
    """
    Generic FASTA parser: accumulates multi-line sequences until next '>'.
    Returns columns: header (raw), sequence (upper).
    """
    if not path.exists():
        return pd.DataFrame(columns=["header", "sequence"])

    records = []
    with path.open("r", errors="ignore") as fh:
        header = None
        seq_chunks = []
        for raw in fh:
            line = raw.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                # flush previous
                if header is not None:
                    seq = "".join(seq_chunks).replace(" ", "").upper()
                    records.append({"header": header, "sequence": seq})
                header = line[1:].strip()
                seq_chunks = []
            else:
                seq_chunks.append(line.strip())
        # last
        if header is not None:
            seq = "".join(seq_chunks).replace(" ", "").upper()
            records.append({"header": header, "sequence": seq})

    return pd.DataFrame.from_records(records)


def parse_final_out_fasta_like(path: Path) -> pd.DataFrame:
    """
    TRUST4 final.out may interleave header/sequence with numeric matrices.
    We capture headers (>) and then accumulate subsequent DNA-looking lines,
    skipping numeric-only blocks, until the next header.
    Returns: header, sequence (upper).
    """
    if not path.exists():
        return pd.DataFrame(columns=["header", "sequence"])

    with path.open("r", errors="ignore") as fh:
        lines = [ln.rstrip("\n") for ln in fh]

    recs = []
    i, n = 0, len(lines)
    while i < n:
        line = lines[i]
        if line.startswith(">"):
            header = line[1:].strip()
            i += 1
            seq_chunks = []
            while i < n and not lines[i].startswith(">"):
                s = lines[i].strip()
                # skip empty & numeric matrix lines
                if s and not _is_numeric_matrix_line(s) and _looks_like_dna(s):
                    seq_chunks.append(s)
                i += 1
            seq = "".join(seq_chunks).replace(" ", "").upper()
            recs.append({"header": header, "sequence": seq})
        else:
            i += 1

    return pd.DataFrame.from_records(recs)


# -------------------- main --------------------

def main():
    args = parse_args()

    # 1) final.out → final_sequences.csv
    df_final = parse_final_out_fasta_like(args.final)
    if df_final.empty:
        df_final["header"] = []
        df_final["sequence"] = []
    # useful extras for final view
    df_final["length"] = df_final["sequence"].fillna("").map(len)
    # Try to add contig_id, gene_hint, chain for context
    df_final["contig_id"] = df_final["header"].str.extract(r"(assemble\d+)", expand=False)
    df_final["gene_hint"] = df_final["header"].str.extract(r"^\s*(?:assemble\d+)?\s*([A-Za-z0-9\-]+)", expand=False)

    def _infer_chain(header: str, gene_hint: Optional[str]) -> Optional[str]:
        m = re.search(r"chain=([A-Za-z]+)", header or "")
        if m:
            return m.group(1).upper()
        gh = (gene_hint or "").upper()
        for ch in ("TRA", "TRB", "TRG", "TRD", "IGH", "IGK", "IGL"):
            if gh.startswith(ch) or ch in (header or "").upper():
                return ch
        return None

    df_final["chain"] = [
        _infer_chain(h, g) for h, g in zip(df_final["header"], df_final["gene_hint"])
    ]

    final_cols = [c for c in ("header", "sequence", "length", "contig_id", "gene_hint", "chain") if c in df_final.columns]
    final_csv = Path(str(args.out_prefix) + "_final_sequences.csv")
    df_final[final_cols].to_csv(final_csv, index=False)

    # 2) assembled_reads.fa → assembled_reads.csv with EXACT columns header, sequence, length, count
    df_assm = parse_fasta_multi(args.assembled)
    if df_assm.empty:
        df_assm["header"] = []
        df_assm["sequence"] = []

    # derive length
    df_assm["length"] = df_assm["sequence"].fillna("").map(len)

    # pull count out of header and optionally strip it from the header field
    counts = []
    cleaned_headers = []
    for h in df_assm["header"].astype(str):
        cnt, cleaned = _parse_count_from_header(h)
        counts.append(cnt)
        cleaned_headers.append(cleaned if not args.keep_raw_header else h)
    count_series = pd.to_numeric(pd.Series(counts, dtype="object"), errors="coerce")

    try:
        # Prefer nullable Int64 if your pandas supports it
        df_assm["count"] = count_series.astype("Int64")
    except (TypeError, ValueError, AttributeError):
        # Fallback for older pandas: keep float if there are NaNs, else cast to int
        if count_series.isna().any():
            df_assm["count"] = count_series.astype(float)
        else:
            df_assm["count"] = count_series.astype(int)
        # enforce exactly the requested columns/order
    df_assm_out = df_assm[["header", "sequence", "length", "count"]]

    assm_csv = Path(str(args.out_prefix) + "_assembled_reads.csv")
    df_assm_out.to_csv(assm_csv, index=False)

    print("Wrote:")
    print(f"  - {final_csv}")
    print(f"  - {assm_csv}")


if __name__ == "__main__":
    main()
