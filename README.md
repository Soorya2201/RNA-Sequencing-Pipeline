# MIXCER-Trust$ RNA-Seq Matching Pipeline

> End-to-end, single-notebook pipeline for RNA-seq–based matching and instant reporting.

---

## Overview

**MIXCER-Trust4** is a fully automated pipeline that takes raw **RNA sequencing runs** (single or multiple “runs of runs”) and turns them into **clean, matched outputs** in a single command.

The pipeline:

- Ingests raw RNA-seq files (e.g., FASTQ from multiple runs / batches)  
- Performs quality control, alignment, and quantification  
- Feeds processed signals into the **MIXCER** and **Trust$** modules  
- Automatically **matches** samples / signatures according to your rules  
- Produces **instant, consolidated reports** ready for downstream analysis  

This entire workflow was **designed and implemented by me** to take a fragmented, manual RNA-seq workflow and turn it into a **robust, reproducible, single-shot pipeline**.

---

## Core Concepts

### MIXCER

`MIXCER` is the module responsible for:

- Combining information from **multiple RNA-seq runs** per sample  
- Handling **mixtures / batches** and normalizing across runs  
- Aggregating expression profiles into a **stable, comparable representation**  

In short, MIXCER makes sure your data is **consistent, de-noised, and comparable**, even when it comes from different runs or sequencing batches.

### Trust$

`Trust$` is the **matching and scoring engine**:

- Takes MIXCER-processed profiles as input  
- Computes **similarity scores, confidence metrics, and rankings**  
- Applies customizable rules to **match samples, conditions, or reference signatures**  
- Outputs interpretable **match tables and summary reports**  

Trust$ turns raw expression profiles into **decisions**: “which sample matches what?” and “how confident are we?”

---

## Pipeline Architecture

The pipeline is structured as:

1. **Input Handling**
   - Accepts single or paired FASTQ (optionally multiple runs per sample)
   - Reads sample metadata / mapping files

2. **Preprocessing & QC**
   - Adapter trimming & quality filtering  
   - Basic QC summaries and logs

3. **Alignment & Quantification**
   - Aligns reads to the reference genome / transcriptome  
   - Generates gene / transcript-level expression matrices

4. **MIXCER Layer**
   - Normalizes expression across runs and batches  
   - Aggregates multiple runs into unified per-sample profiles  
   - Optionally performs deconvolution / mixture modelling (depending on config)

5. **Trust$ Matching Engine**
   - Computes similarity/matching scores  
   - Ranks candidate matches  
   - Flags high-confidence / ambiguous cases

6. **Reporting**
   - Final **match table** (CSV/TSV)  
   - Per-sample summary reports  
   - Run-level logs for full reproducibility  

---

## Key Features

- **Single-Notebook Workflow**  
  From raw reads to final matches by running one notebook top to bottom — no juggling separate scripts per stage.

- **Run-of-Runs Aware**  
  Designed for situations where you have **multiple runs per sample** and need them merged intelligently.

- **Automated Matching**  
  No more manual spreadsheets or ad-hoc scripts. Trust$ handles matching and scoring end-to-end.

- **Reproducible & Config-Driven**  
  Every parameter is controlled via a config file, making the pipeline portable and repeatable.

- **Scalable**  
  Works on small test datasets all the way up to large cohorts and multi-run projects.

- **Instant, Consolidated Output**  
  Final results are delivered as clean tables and reports, ready for plotting, statistics, or integration into downstream workflows.

---

## Repository Structure

This project is implemented as a set of Jupyter notebooks rather than a single CLI script:

- `TCR_seq_Pipeline.ipynb` — main pipeline: QC, alignment, quantification, and MIXCER/Trust4 matching
- `TCR_Sequence_Extracror.ipynb` — TCR sequence extraction from raw reads
- `Biotech_data_Reserach (1).ipynb` — exploratory research and analysis
- `Fast_Mixcer/`, `Mixer-Output_Refined/`, `Double Run/`, `Outputs/`, `Trust4_outputs refined/` — output directories from prior pipeline runs

## Running the Pipeline

```bash
git clone https://github.com/Soorya2201/RNA-Sequencing-Pipeline.git
cd RNA-Sequencing-Pipeline

python -m venv .venv
source .venv/bin/activate

pip install jupyter pandas numpy biopython
jupyter notebook
```

Open `TCR_seq_Pipeline.ipynb` and run the cells top to bottom against your own FASTQ input. Adjust file paths and config values at the top of the notebook for your dataset.

---

## License

MIT — see [LICENSE](LICENSE) for details.
