# AChE-Insect-Transcriptome
### **A Multi-Species Pipeline for Analyzing Acetylcholinesterase Gene Evolution**

This repository contains a complete workflow for downloading, processing, and analyzing acetylcholinesterase (AChE) transcripts across **locusts, crickets, and mosquitoes**. The pipeline identifies the **longest isoform per AChE gene**, builds a combined multi-species dataset, prepares sequence alignments, constructs phylogenetic trees, and produces final visualizations for research and poster presentation.

---

## 🧬 **Project Overview**

Locusts possess unusually high AChE gene copy numbers yet **fail to evolve organophosphate pesticide resistance**, unlike fast-resistance species such as mosquitoes. This project investigates:

- **AChE copy-number differences** across insects  
- Whether locust AChE expansion is **ancestral** (shared with crickets) or **derived**  
- How gene family size may constrain resistance evolution  
- Whether any AChE copies carry **known resistance-associated mutations**  

This repository contains the code and structure enabling all computational steps.

---

## 🗂 **Repository Structure**

AChE_Project/
│ AChE_Project.Rproj
│ README.md
│ .gitignore
│
├── scripts/
│ 01_download_AChE_sequences.R
│ 02_extract_longest_isoforms.R
│ 03_merge_longest_fastas.R
│ 04_run_alignment.R
│ 05_build_tree_fasttree.R
│ 06_plot_tree.R
│
└── data/
├── raw/ # raw NCBI downloads (ignored by Git)
├── longest/ # longest isoforms per gene (ignored by Git)
├── combined/ # merged multi-species FASTA (ignored by Git)
├── alignment/ # alignment outputs (ignored by Git)
└── tree/ # final tree files (committed)

The repo is structured for **reproducibility**, **clarity**, and **collaborator usability**.

---

## ⚙️ **Pipeline Summary**

### **1. Download AChE Sequences**  
`scripts/01_download_AChE_sequences.R`  
- Pulls all AChE CDS/mRNA sequences for each species from NCBI  
- Outputs to `data/raw/`  

### **2. Extract Longest Isoforms**  
`scripts/02_extract_longest_isoforms.R`  
- Groups isoforms by gene name  
- Keeps the longest CDS per gene  
- Outputs to `data/longest/`  

### **3. Merge Multi-Species FASTA**  
`scripts/03_merge_longest_fastas.R`  
- Combines longest isoforms across species  
- Produces the dataset used for all downstream analyses  

### **4. Prepare Alignment**  
`scripts/04_run_alignment.R`  
- Runs or prepares input for MAFFT or Clustal Omega  
- Saves alignments to `data/alignment/`  

### **5. Build Phylogenetic Tree**  
`scripts/05_build_tree_fasttree.R`  
- Runs FastTree/IQ-TREE on the alignment  
- Saves output to `data/tree/`  

### **6. Plot Tree**  
`scripts/06_plot_tree.R`  
- Reads the final tree and generates publishable figures  

---

## 🧪 **Species Panel**

### **Locusts (Schistocerca spp.)**
- *S. gregaria*  
- *S. americana*  
- *S. piceifrons*

### **Cricket**
- *Anabrus simplex* (Mormon cricket)

### **Mosquitoes**
- *Anopheles gambiae*  
- *Aedes aegypti*  
- *Culex quinquefasciatus*

---

## 📦 **Dependencies**

- **R ≥ 4.0**
- Packages:  
  `rentrez`, `seqinr`, `dplyr`, `stringr`, `ape`, `ggtree`, `ggplot2`, `readr`

External tools (optional depending on workflow):
- **MAFFT** (multiple sequence alignment)
- **FastTree** or **IQ-TREE** (phylogenetics)

---

## 🧑‍💻 **Usage**

1. Clone the repository  
2. Open `AChE_Project.Rproj` in RStudio  
3. Run scripts **in numerical order (01 → 06)**  
4. Add new species or update raw data simply by adding FASTA files and re-running  

---

## 👩‍🔬 **Author**

**Taylor M. Johnson**  
Department of Biochemistry  
Mississippi State University  

This repository is actively updated as part of an undergraduate research project supervised by Hunter Walt and Dr. Federico Hoffmann.

---
