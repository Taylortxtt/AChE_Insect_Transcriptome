# 📘 AChE Insect Transcriptome  
### Pipeline for Identifying Longest AChE Transcripts Across Insect Species

[![R](https://img.shields.io/badge/R-≥4.0-blue)](https://www.r-project.org/)
![Status](https://img.shields.io/badge/status-active-success)
![License](https://img.shields.io/badge/license-MIT-lightgrey)

This repository contains a complete, reproducible pipeline for retrieving, processing, aligning, and analyzing **acetylcholinesterase (AChE)** transcripts across multiple insect species.

The workflow supports *locusts*, *crickets*, and *mosquitoes* and includes:

- Automated NCBI sequence download  
- Extraction of the **longest transcript per gene**  
- Multi-species FASTA construction  
- MAFFT alignment  
- Phylogenetic reconstruction with  
  - **FastTree** (fast, approximate ML)  
  - **IQ-TREE** (slower, model-tested ML)  
- Tree visualization (FastTree & IQ-TREE)  
- **Mutation screening**  
- AChE copy-number summaries  
- A clean config system (no PATH hacking required)

The goal is to evaluate AChE gene expansion patterns and detect functional or insecticide-resistance–associated variants, while keeping the pipeline simple, readable, and easy for collaborators to run.

---

## ⚡ Quickstart

1. **Clone or download** this repository  
   ```bash
   git clone https://github.com/<your-username>/AChE_Insect_Transcriptome.git
   cd AChE_Insect_Transcriptome
   ```

2. **Install external tools** (MAFFT, FastTree, IQ-TREE)

   **macOS (Homebrew):**
   ```bash
   brew install mafft fasttree iqtree
   ```

   **Ubuntu / Debian:**
   ```bash
   sudo apt-get install mafft fasttree iqtree
   ```

   **Windows:**
   - Install MAFFT  
   - Install FastTree  
   - Install IQ-TREE  
   - Note installation paths

3. **Set tool paths in `config.R`**

   Open `config.R` and set:

   ```r
   MAFFT_PATH    <- "/path/to/mafft"
   FASTTREE_PATH <- "/path/to/FastTree"
   IQTREE_PATH   <- "/path/to/iqtree2"
   ```

4. **Run the full pipeline in R**

   ```r
   source("run_pipeline.R")
   ```

5. Check the outputs in:
   - `data/combined/`
   - `data/alignment/`
   - `data/tree/`
   - any summary tables / plots in `exports/` (if used)

---

## 🧬 Species Included

- **Locusts:** *Schistocerca gregaria*, *S. cancellata*, *S. piceifrons*  
- **Cricket:** *Anabrus simplex*  
- **Mosquitoes:** *Aedes aegypti*, *Anopheles gambiae*, *Culex quinquefasciatus*

These species span high-copy and low-copy AChE lineages, allowing comparisons involving gene family expansion and potential insecticide resistance.

---

## 🗂 Project Structure

```text
AChE_Insect_Transcriptome/
│
├── AChE_Project.Rproj
├── config.R                       # MAFFT, FastTree, IQ-TREE executable paths
├── run_pipeline.R                 # One-click full pipeline
│
├── scripts/
│   ├── 01_download_AChE_sequences.R
│   ├── 02_extract_longest_isoforms.R      # longest transcript per *gene*
│   ├── 03_merge_longest_fastas.R
│   ├── 04_run_alignment.R
│   ├── 05_build_tree_fasttree.R
│   ├── 05b_build_tree_iqtree.R           # optional IQ-TREE workflow
│   ├── 06_plot_tree.R
│   ├── 06b_plot_tree_iqtree.R            # optional IQ-TREE plot
│   ├── 07_summarize_copy_number.R        # AChE gene counts per species
│   ├── 08_screen_mutations.R             # mutation presence/absence
│
├── data/
│   ├── raw/            # Raw FASTAs from NCBI
│   ├── longest/        # Per-species longest-per-gene FASTAs
│   ├── combined/       # Combined FASTAs + metadata CSV
│   ├── alignment/      # MAFFT alignment output
│   └── tree/           # FastTree & IQ-TREE files + PNG plots
│
└── exports/            # Optional: figures, tables, reports
```

---

## 🔧 Configuration (`config.R`)

This file stores the **absolute paths** to the external tools.  
You edit this file **once**, and all scripts use these paths automatically.

```r
# config.R
MAFFT_PATH    <- "/path/to/mafft"
FASTTREE_PATH <- "/path/to/FastTree"
IQTREE_PATH   <- "/path/to/iqtree2"
```

### Example – macOS (Apple Silicon, Homebrew)

```r
MAFFT_PATH    <- "/opt/homebrew/bin/mafft"
FASTTREE_PATH <- "/opt/homebrew/bin/FastTree"
IQTREE_PATH   <- "/opt/homebrew/bin/iqtree2"
```

### Example – Windows

```r
MAFFT_PATH    <- "C:/Program Files/mafft/mafft.bat"
FASTTREE_PATH <- "C:/Program Files/FastTree/FastTree.exe"
IQTREE_PATH   <- "C:/Program Files/iqtree2/iqtree2.exe"
```

If a path is wrong, the corresponding script will stop with a clear error message telling you what to fix.

---

## 🔁 Pipeline Overview

### One-click run

Once `config.R` is set up:

```r
source("run_pipeline.R")
```

By default, `run_pipeline.R` calls, in order:

1. `01_download_AChE_sequences.R` – download raw AChE sequences from NCBI  
2. `02_extract_longest_isoforms.R` – keep the **longest transcript per gene**  
3. `03_merge_longest_fastas.R` – merge all species into one FASTA + metadata CSV  
4. `04_run_alignment.R` – align with MAFFT  
5. `05_build_tree_fasttree.R` – FastTree ML tree  
6. `05b_build_tree_iqtree.R` – IQ-TREE ML + support values  
7. `06_plot_tree.R` – plot the FastTree tree  
8. `06b_plot_tree_iqtree.R` – plot the IQ-TREE tree  
9. `07_summarize_copy_number.R` – AChE copy numbers per species  
10. `08_screen_mutations.R` – mutation presence/absence table  

You can comment out any `run_script()` line inside `run_pipeline.R` if you want to skip a step.

---

## 🔀 Pipeline Diagram

```text
         ┌────────────────────────────┐
         │ 01_download_AChE_sequences │
         └─────────────┬──────────────┘
                       ↓
         ┌────────────────────────────┐
         │ 02_extract_longest_isoforms│  (longest per gene)
         └─────────────┬──────────────┘
                       ↓
         ┌────────────────────────────┐
         │ 03_merge_longest_fastas    │  (multi-species FASTA + metadata)
         └─────────────┬──────────────┘
                       ↓
         ┌────────────────────────────┐
         │ 04_run_alignment (MAFFT)   │
         └─────────────┬──────────────┘
                       ↓
        ┌──────────────┼───────────────┐
        ↓                              ↓
┌─────────────────────┐       ┌─────────────────────┐
│05_build_tree_fasttree│      │05b_build_tree_iqtree│
└─────────────┬────────┘       └─────────────┬──────┘
              ↓                              ↓
  ┌────────────────────┐           ┌────────────────────────┐
  │06_plot_tree        │           │06b_plot_tree_iqtree    │
  └─────────────┬──────┘           └─────────────┬──────────┘
                ↓                              ↓
       ┌───────────────────┐        ┌────────────────────────┐
       │07_summarize_copy_ │        │08_screen_mutations     │
       │   number          │        │ (resistance mutations) │
       └───────────────────┘        └────────────────────────┘
```

---

## 📂 Key Outputs

### Combined FASTA & Metadata

```text
data/combined/ache_longest_by_gene_all_species.fasta
data/combined/ache_longest_by_gene_metadata.csv
```

### Alignment

```text
data/alignment/ache_longest_by_gene_aligned.fasta
```

### FastTree Outputs

```text
data/tree/ache_longest_by_gene_fasttree.nwk
data/tree/ache_longest_by_gene_tree.png
```

### IQ-TREE Outputs

```text
data/tree/ache_longest_by_gene_iqtree.treefile
data/tree/ache_longest_by_gene_iqtree.log
data/tree/ache_longest_by_gene_iqtree.iqtree
data/tree/ache_longest_by_gene_iqtree_plot.png
# plus IQ-TREE support / model files
```

### Summaries

- Copy number table + plots (from `07_summarize_copy_number.R`)  
- Mutation presence/absence table (from `08_screen_mutations.R`)

---

## 🧪 Methods Summary

### Data Download

NCBI nuccore is queried via `{rentrez}` using:

- Organism filters (genus + species names)
- AChE-related keywords (acetylcholinesterase, ache, ace)
- mRNA/CDS-focused search terms

Raw FASTAs are written to `data/raw/`.

### Longest-Per-Gene Extraction

Raw sequences are parsed, gene IDs are inferred from FASTA headers using:

- `gene=...` annotations when present  
- `LOC...` identifiers  
- Accession-root fallback (e.g., `XM_123456789`)

For each gene, the **longest isoform** is retained. Outputs go to `data/longest/`.

### Multi-Species Merging

Species-specific FASTAs are merged into:

- `data/combined/ache_longest_by_gene_all_species.fasta`  
- `data/combined/ache_longest_by_gene_metadata.csv`

Headers are standardized as:

```text
shortname|original_header
```

making it easy to link back to species and original records.

### Alignment (MAFFT)

The combined FASTA is aligned with MAFFT using:

```text
--auto
```

The resulting alignment is saved in `data/alignment/`.

### Tree Inference

Two tree builders are supported:

- **FastTree**: approximate ML tree with GTR + Gamma  
- **IQ-TREE**: model-finding (`-m MFP`), ultrafast bootstrap (`-bb 1000`), and SH-aLRT (`-alrt 1000`)

Trees are output to `data/tree/` in Newick format, plus PNG visualizations generated by `06_*` scripts.

### Mutation Screening

`08_screen_mutations.R` scans aligned AChE sequences for a curated set of resistance-associated mutations reported in the literature.  
Output is a presence/absence matrix that can be used for:

- Heatmaps  
- Annotation layers on phylogenies  
- Simple summary tables

---

## 📜 License & Citation

- **License:** See `LICENSE` file in this repository.  
- **Citation:** See `CITATION.cff` for citation metadata.

If you use or adapt this pipeline, please cite this repository and acknowledge **Taylor Johnson**.

---

## 📬 Contact

Maintained by **Taylor Johnson**.  
For questions, suggestions, or collaboration, feel free to open an issue or pull request.
