# Immune Repertoire Profiling from Bulk RNA-seq

**Reconstructs T-cell and B-cell receptor repertoires (TRA/TRB/TRG/TRD, IGH/IGK/IGL) from bulk RNA-seq using MiXCR and TRUST4, then quantifies clonal diversity, V/J gene usage and CDR3 length per locus.**

![License](https://img.shields.io/badge/license-MIT-yellow.svg)

## Highlights

- **Two independent reconstruction tools on the same reads.** MiXCR 4.6 (`align → assemble → exportClones`) and TRUST4 (against `hg38_bcrtcr.fa`) were run on two public paired-end RNA-seq samples from the European Nucleotide Archive (SRR1918236, SRR1918241), so their clonotype calls can be compared.
- **Repertoire analytics** (`Mixer-Output_Refined/analyze_clonotypes py.py`): per-locus Shannon and Simpson diversity, top clonotypes, cumulative clonal fraction, V and J gene usage, and CDR3 amino-acid length distributions across all 7 loci.
- **A lightweight V(D)J caller** (`Fast_Mixcer/faster_Mixcer.py`): CDR3 boundary detection, a Phred-quality gate, and k-mer Jaccard scoring (k = 15) of reads against V/D/J reference segments. It is a fast, transparent way to sanity-check the full aligners.
- **An interactive Colab front end** (`TCR_seq_Pipeline.ipynb`) with ipywidgets for selecting samples and segments, running the analysis and generating a summary report.

## Results (patient 1, MiXCR)

| Locus | Clonotypes | Shannon | Simpson |
|---|---:|---:|---:|
| IGH | 596 | 6.39 | 0.998 |
| IGK | 808 | 6.69 | 0.999 |
| IGL | 732 | 6.60 | 0.999 |
| TRB | 46 | 3.83 | 0.978 |
| TRA | 11 | 2.40 | 0.909 |
| TRD | 1 | 0.00 | 0.000 |
| TRG | 1 | 0.00 | 0.000 |

B-cell loci are highly diverse with no dominant clone. T-cell receptor recovery is sparse, as expected from bulk (not TCR-enriched) RNA-seq: few TCR transcripts are captured, so TRA/TRB/TRD/TRG counts are low. Plots are in `Mixer-Output_Refined/Outputs/` (`V_usage_*.png`, `J_usage_*.png`, `cdr3_len_*.png`, `cumulative_fraction_overall.png`).

## Pipeline

```
ENA FASTQ (paired-end, SRR1918236 / SRR1918241)
  ├── MiXCR  align -s hs -p rna-seq → assemble → exportClones  → per-locus clonotype tables
  ├── TRUST4 -f hg38_bcrtcr.fa -1 R1 -2 R2                       → CDR3 calls, assembled reads, AIRR output
  └── faster_Mixcer.py (k-mer Jaccard V/D/J scoring)             → quick clonotype calls for IGH/IGK
                │
                ▼
analyze_clonotypes: diversity (Shannon, Simpson), V/J usage, CDR3 length, top clones
```

## Repository layout

| Path | Contents |
|---|---|
| `Biotech_data_Reserach (1).ipynb` | End-to-end research notebook: tool installs, ENA downloads, MiXCR and TRUST4 runs |
| `TCR_Sequence_Extracror.ipynb` | Standalone MiXCR extraction run |
| `TCR_seq_Pipeline.ipynb` | Interactive Colab UI (see note below) |
| `Fast_Mixcer/` | Lightweight k-mer V(D)J caller and its IGH/IGK results |
| `Mixer-Output_Refined/` | Clonotype analysis script, diversity table and plots |
| `Trust4_outputs refined/` | TRUST4 post-processing script and cleaned outputs |
| `Outputs/`, `Double Run/` | Raw MiXCR, TRUST4 and Bulk2Space outputs from the runs |

## Run it

```bash
pip install pandas numpy matplotlib biopython
python "Mixer-Output_Refined/analyze_clonotypes py.py" \
    --input-glob "Double Run/patient_1_mixcr_clonotypes_*.tsv" \
    --outdir results/ --sample patient_1                         # diversity, V/J usage, CDR3 plots
python Fast_Mixcer/faster_Mixcer.py --help                        # k-mer V(D)J caller
```

To regenerate clonotypes, open `Biotech_data_Reserach (1).ipynb` in Google Colab. MiXCR needs a free academic license: set it as `MIXCR_LICENSE_KEY` and never commit it.

## Limitations

- Two samples only: this is a methods walk-through, not a cohort study.
- `TCR_seq_Pipeline.ipynb` currently runs in **demo mode**: its `run_trust4` step writes example CDR3 output so the UI can be exercised without a GPU or reference download. The real TRUST4 and MiXCR runs live in the research notebook. Wiring the UI to the real tools is the next step.
- `faster_Mixcer.py` uses a small, manually curated reference set; it is a sanity check, not a replacement for MiXCR or IMGT-based alignment.

## Tools and credits

[MiXCR](https://github.com/milaboratory/mixcr) (MiLaboratories) · [TRUST4](https://github.com/liulab-dfci/TRUST4) (Liu Lab, Dana-Farber) · Python · pandas · matplotlib · Google Colab

## License

MIT for the code in this repository — see [LICENSE](LICENSE). MiXCR and TRUST4 have their own licenses.
