# AChE Gene Evolution – Insect Transcriptome Pipeline

Pipeline for identifying the longest acetylcholinesterase (AChE) transcripts across locusts, crickets, and mosquitoes, and preparing data for phylogenetic analysis.

## Pipeline Steps

1. **Download sequences**
   - `scripts/01_download_AChE_sequences.R`
   - Downloads AChE CDS/mRNA sequences from NCBI for each species into `data/raw/`.

2. **Extract longest isoforms**
   - `scripts/02_extract_longest_isoforms.R`
   - For each gene, keeps only the longest CDS isoform and writes results to `data/longest/`.

3. **Merge longest isoforms**
   - `scripts/03_merge_longest_fastas.R`
   - Combines species-specific FASTA files into a multi-species FASTA in `data/combined/`.

4. **Run alignment**
   - `scripts/04_run_alignment.R`
   - Prepares or runs multiple sequence alignment (e.g., MAFFT) and saves output to `data/alignment/`.

5. **Build phylogenetic tree**
   - `scripts/05_build_tree_fasttree.R`
   - Uses the alignment to build a phylogenetic tree, saved in `data/tree/`.

6. **Plot tree**
   - `scripts/06_plot_tree.R`
   - Reads tree files and generates figures for the poster/manuscript.

## Directory Structure

- `scripts/` – All R scripts for the pipeline (01–06).
- `data/raw/` – Raw NCBI downloads (ignored by Git).
- `data/longest/` – Longest-isoform FASTA files (ignored by Git).
- `data/combined/` – Combined multi-species FASTA files (ignored by Git).
- `data/alignment/` – Alignment outputs (ignored by Git).
- `data/tree/` – Final tree files and small outputs that are tracked.

## Author

Taylor M. Johnson  
Department of Biochemistry, Mississippi State University
