# AChE Gene Evolution – Insect Transcriptome Pipeline
### A Multi-Species Workflow for Acetylcholinesterase (AChE) Gene Analysis

This repository contains a complete R-based bioinformatics pipeline for collecting, filtering, aligning, and analyzing acetylcholinesterase (AChE) transcripts across **locusts**, **crickets**, and **mosquitoes**. The workflow identifies the longest isoforms, builds multi-species FASTA files, prepares alignments, constructs phylogenetic trees, and generates figures for research poster or manuscript preparation.

---

## 🧬 Project Overview

Locusts possess unusually high AChE gene copy numbers yet **fail to evolve organophosphate pesticide resistance**, unlike mosquitoes. This project explores:

- AChE copy-number differences across key insects  
- Whether locust AChE expansion is ancestral or derived  
- How gene family size affects resistance evolution  
- Whether AChE copies contain known resistance-associated mutations  

This repository provides the complete computational workflow for that analysis.

---

## ⚙️ Pipeline Summary

### **1. Download AChE sequences**  
**File:** `scripts/01_download_AChE_sequences.R`  
Downloads all AChE CDS/mRNA isoforms from NCBI into `data/raw/`.

### **2. Extract longest isoforms**  
**File:** `scripts/02_extract_longest_isoforms.R`  
Identifies the longest CDS per AChE gene and writes results to `data/longest/`.

### **3. Merge longest isoforms**  
**File:** `scripts/03_merge_longest_fastas.R`  
Builds a unified multi-species FASTA in `data/combined/`.

### **4. Prepare alignment**  
**File:** `scripts/04_run_alignment.R`  
Prepares or runs MAFFT/Clustal alignment and saves files to `data/alignment/`.

### **5. Build phylogenetic tree**  
**File:** `scripts/05_build_tree_fasttree.R`  
Constructs a phylogenetic tree using FastTree or IQ-TREE, saved to `data/tree/`.

### **6. Plot phylogeny**  
**File:** `scripts/06_plot_tree.R`  
Generates polished tree visualizations.

---

## 🗂 Repository Structure (Non-breaking format)

**Top-level files**
- `AChE_Project.Rproj` — RStudio project file  
- `README.md` — documentation  
- `LICENSE` — MIT license  
- `CITATION.cff` — citation metadata  
- `.gitignore` — ignored files  

**Folder: `scripts/`**
- 01_download_AChE_sequences.R  
- 02_extract_longest_isoforms.R  
- 03_merge_longest_fastas.R  
- 04_run_alignment.R  
- 05_build_tree_fasttree.R  
- 06_plot_tree.R  

**Folder: `data/`**
- `raw/` — raw NCBI downloads *(ignored by Git)*  
- `longest/` — longest isoform FASTAs *(ignored by Git)*  
- `combined/` — merged multi-species dataset *(ignored by Git)*  
- `alignment/` — alignment files *(ignored by Git)*  
- `tree/` — final tree outputs *(tracked)*  

---

## 🦗 Species Included

**Locusts (Schistocerca spp.)**  
- *Schistocerca gregaria*  
- *Schistocerca americana*  
- *Schistocerca piceifrons*

**Cricket**  
- *Anabrus simplex* (Mormon cricket)

**Mosquitoes**  
- *Anopheles gambiae*  
- *Aedes aegypti*  
- *Culex quinquefasciatus*

---

## 📦 Dependencies

R packages:
- rentrez  
- seqinr  
- dplyr  
- stringr  
- ape  
- ggtree  
- ggplot2  
- readr  

External tools:
- **MAFFT** — alignment  
- **FastTree / IQ-TREE** — phylogenetics  

---

## 🚀 Quickstart

**1. Clone the repository**
git clone https://github.com/Taylortxtt/AChE_Insect_Transcriptome.git

**2. Open the project**
- Open `AChE_Project.Rproj` in RStudio

**3. Run the pipeline**
Run scripts **in numerical order**:

1 → 2 → 3 → 4 → 5 → 6

Each script outputs files into the corresponding `data/` subfolder.

**4. Add new species**
- Drop FASTA files into `data/raw/`  
- Re-run the pipeline starting at script 02 or 03

---

## 👩‍🔬 Author

**Taylor M. Johnson**  
Department of Biochemistry  
Mississippi State University  

---

## 📄 License

MIT License — see the `LICENSE` file.

## 📣 Citation

Please cite using the included `CITATION.cff` file.

---
