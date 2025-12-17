#!/usr/bin/env python3
import argparse, gzip, logging
from collections import defaultdict, Counter
from pathlib import Path
import re

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
log = logging.getLogger("mixcer-manual")

# ===== 1) PASTE / EDIT YOUR MANUAL REFERENCES HERE (IGH + IGK) =====
# Use UPPERCASE DNA (ACGT). Add as many entries as you want.
# Keys are gene names; values are nucleotide sequences (no spaces/newlines).
MANUAL_REFS = {
    "IGH": {
        "V": {
            # EXAMPLE from your snippet:
            "IGHV1-18*01": (
                "CAGGTGCAGCTGGTGCAGTCTGGGGCTGAGGTGAAGAAGCCTGGGGCCTCAGTGAAGGTCTCCTGCAAGGCTTCTGGATACACCTTCACCAACTATGCTATGCACTGGGTGCGACAGGCCCCTGGACAAGGGCTTGAGTGGATGGGATGGATCAACCCTAACAGTGGTGGCACAAACTATGCACAGAAGTTTCAGGGCAGGGTCACCATGACCAGGGACACGTCCATCAGCACAGCCTACATGGAGCTGAGCAGGCTGAGATCTGACGACACGGCCGTGTATTACTGTGCGAGAGA"
            ),
            # Add more IGH V genes here...
        },
        "D": {
            "IGHD2-2*01": "AGGATATTGTAGTAGTACCAGCTGCTATGGAC",
            # Add more IGH D genes here...
        },
        "J": {
            "IGHJ1*01": "CTGAATACTTCCAGCACTGGGGCCAGGGCACCCTGGTCACCGTCTCCTCAG",
            # Add more IGH J genes here...
        },
        "C": {
            # Optional: IGH constant region snips for quick chain routing (not required)
            # "IGHA1": "..." 
        }
    },
    "IGK": {
        "V": {
            "IGKV1-5*03": (
                "GACATCCAGATGACCCAGTCTCCATCCTCACTGTCTGCATCTGTAGGAGACAGAGTCACCATCACTTGCCGGGCAAGTCAGGGCATTAGAAATTATTTAGCCTGGTATCAGCAGAAACCAGGGAAAGCCCCTAAGCTCCTGATCTATGCTGCATCCAGCTTGCAAAGTGGGGTCCCATCAAGGTTCAGTGGCAGTGGATCTGGGACAGATTTCACTCTCACCATCAGCAGTCTGCAACCTGAAGATTTTGCAACTTACTATTGTCAACAG"
            ),
            # Add more IGK V genes here...
        },
        "D": {},  # light chain has no D
        "J": {
            "IGKJ2*01": "TTTATACCCCGTACACTTTTGGCCAGGGGACCAAGCTGGAGATCAAA",
            # Add more IGK J genes here...
        },
        "C": {
            # Optional IGKC snippet(s)
        }
    }
}
# ===== END MANUAL REFS =====

# ----- CDR3 motif finders (DNA) -----
CDR3_V_CYS_DNA = re.compile(r"TGT", re.I)  # cysteine codon
# very tolerant J motif (FGXG-like, WGXG-like)
CDR3_J_MOTIF_DNA = re.compile(r"(TT[TC][ACGT]GG|TGGGG[ACGT]GG)", re.I)

def translate(nt: str) -> str:
    table = {
        'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
        'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
        'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
        'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
        'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
        'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
        'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
        'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G'
    }
    nt = nt.upper()
    aa = []
    for i in range(0, len(nt) - 2, 3):
        aa.append(table.get(nt[i:i+3], 'X'))
    return "".join(aa)

def min_phred_in_slice(qual: str, s: int, e: int) -> int:
    if not qual or e <= s: return 0
    return min(max(ord(c)-33, 0) for c in qual[s:e])

def find_cdr3_nt_bounds(read_nt: str):
    v_cys = None
    for m in CDR3_V_CYS_DNA.finditer(read_nt):
        v_cys = m.start()
    if v_cys is None:
        return None, None
    jm = CDR3_J_MOTIF_DNA.search(read_nt, pos=v_cys + 3)
    if not jm:
        return None, None
    start = v_cys
    end = jm.start() + 3
    rem = (end - start) % 3
    if rem != 0:
        end -= rem
    if end <= start or (end-start) < 9:
        return None, None
    return start, end

# ----- very simple k-mer Jaccard scoring -----
def build_kmer_index(refs_for_seg: dict, k: int):
    idx = defaultdict(set)
    for name, seq in refs_for_seg.items():
        s = seq.upper()
        if len(s) < k: 
            continue
        for i in range(len(s)-k+1):
            idx[s[i:i+k]].add(name)
    return idx

def jaccard_hits(seq: str, idx: dict, k: int, topk: int = 3):
    s = seq.upper()
    if len(s) < k: 
        return []
    q = {s[i:i+k] for i in range(len(s)-k+1)}
    scores = Counter()
    for mer in q:
        if mer in idx:
            for name in idx[mer]:
                scores[name] += 1
    if not scores:
        return []
    denom = len(q)
    out = [(name, cnt/denom) for name, cnt in scores.items()]
    out.sort(key=lambda x: x[1], reverse=True)
    return out[:topk]

def write_tsv(path: Path, rows: list, header: list):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as w:
        w.write("\t".join(header) + "\n")
        for r in rows:
            w.write("\t".join(str(r.get(h,"")) for h in header) + "\n")

def main():
    p = argparse.ArgumentParser(description="Manual IGH/IGK clonotype caller (no downloads).")
    p.add_argument("--fastq", required=True, help="FASTQ(.gz)")
    p.add_argument("--outdir", required=True, help="Output directory")
    p.add_argument("--k", type=int, default=15, help="k-mer size for V/J/D scoring")
    p.add_argument("--min-cdr3-aa", type=int, default=5)
    args = p.parse_args()

    outdir = Path(args.outdir)
    fastq = Path(args.fastq)

    # Build indices only for IGH and IGK
    chains = ["IGH", "IGK"]
    refs = {ch: MANUAL_REFS[ch] for ch in chains}

    log.info("Building manual k-mer indices...")
    kidx = {ch: {seg: build_kmer_index(refs[ch][seg], args.k) for seg in refs[ch]} for ch in chains}

    header_cols = [
        "cloneId","readCount","readFraction","targetSequences","targetQualities",
        "allVHitsWithScore","allDHitsWithScore","allJHitsWithScore","allCHitsWithScore",
        "allVAlignments","allDAlignments","allJAlignments","allCAlignments",
        "nSeqCDR3","minQualCDR3","aaSeqCDR3","refPoints"
    ]

    per_chain_reads = {ch: [] for ch in chains}

    openf = gzip.open if str(fastq).endswith(".gz") else open
    with openf(fastq, "rt") as fh:
        total_reads = 0
        while True:
            h = fh.readline()
            if not h: break
            s = fh.readline().strip().upper()
            fh.readline()  # +
            q = fh.readline().strip()
            total_reads += 1

            cdr3_s, cdr3_e = find_cdr3_nt_bounds(s)
            if cdr3_s is None:
                continue
            nseq = s[cdr3_s:cdr3_e]
            aaseq = translate(nseq)
            if "*" in aaseq or len(aaseq) < args.min_cdr3_aa:
                continue
            minq = min_phred_in_slice(q, cdr3_s, cdr3_e)

            for ch in chains:
                Vh = jaccard_hits(s, kidx[ch]["V"], args.k) if kidx[ch]["V"] else []
                Jh = jaccard_hits(s, kidx[ch]["J"], args.k) if kidx[ch]["J"] else []
                # Require at least a V or a J match to consider this chain
                if not Vh and not Jh:
                    continue
                Dh = jaccard_hits(s, kidx[ch]["D"], args.k) if "D" in kidx[ch] and kidx[ch]["D"] else []
                Ch = jaccard_hits(s, kidx[ch]["C"], args.k) if "C" in kidx[ch] and kidx[ch]["C"] else []

                per_chain_reads[ch].append({
                    "readId": h[1:].strip(),
                    "sequence": s,
                    "quality": q,
                    "nSeqCDR3": nseq,
                    "aaSeqCDR3": aaseq,
                    "minQualCDR3": minq,
                    "V_hits": Vh, "D_hits": Dh, "J_hits": Jh, "C_hits": Ch
                })

    log.info("Collapsing to clonotypes...")
    for ch in chains:
        reads = per_chain_reads[ch]
        if not reads:
            log.info(f"{ch}: no reads matched manual refs.")
            continue

        # Clone key: aaCDR3 + top V + top J
        clones = defaultdict(list)
        for r in reads:
            v_top = r["V_hits"][0][0] if r["V_hits"] else "Unknown"
            j_top = r["J_hits"][0][0] if r["J_hits"] else "Unknown"
            clones[(r["aaSeqCDR3"], v_top, j_top)].append(r)

        total = sum(len(v) for v in clones.values())
        rows = []
        for i, ((aa, v, j), bucket) in enumerate(sorted(clones.items(), key=lambda kv: len(kv[1]), reverse=True), 1):
            rep = bucket[0]
            rows.append({
                "cloneId": i,
                "readCount": float(len(bucket)),
                "readFraction": (len(bucket)/total) if total else 0.0,
                "targetSequences": rep["sequence"],
                "targetQualities": rep["quality"],
                "allVHitsWithScore": ",".join(f"{n}({sc:.3f})" for n, sc in rep["V_hits"]),
                "allDHitsWithScore": ",".join(f"{n}({sc:.3f})" for n, sc in rep["D_hits"]),
                "allJHitsWithScore": ",".join(f"{n}({sc:.3f})" for n, sc in rep["J_hits"]),
                "allCHitsWithScore": ",".join(f"{n}({sc:.3f})" for n, sc in rep["C_hits"]),
                "allVAlignments": "", "allDAlignments": "", "allJAlignments": "", "allCAlignments": "",
                "nSeqCDR3": rep["nSeqCDR3"],
                "minQualCDR3": rep["minQualCDR3"],
                "aaSeqCDR3": aa,
                "refPoints": ""
            })

        out = Path(args.outdir) / f"clonotypes.{ch.lower()}.tsv"
        header = [
            "cloneId","readCount","readFraction","targetSequences","targetQualities",
            "allVHitsWithScore","allDHitsWithScore","allJHitsWithScore","allCHitsWithScore",
            "allVAlignments","allDAlignments","allJAlignments","allCAlignments",
            "nSeqCDR3","minQualCDR3","aaSeqCDR3","refPoints"
        ]
        write_tsv(out, rows, header)
        log.info(f"{ch}: wrote {len(rows)} clonotypes → {out}")

if __name__ == "__main__":
    main()
