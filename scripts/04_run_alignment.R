############################################################
# 04_run_alignment.R
#
# Goal
# ----
# Align the multi-species AChE FASTA using MAFFT.
#
# Input
# -----
#  Preferred (if present):
#   data/combined/ache_longest_by_gene_ace12_blast_only.fasta
#
#  Otherwise:
#   data/combined/ache_longest_by_gene_all_species.fasta
#
# Output
# ------
#   data/alignment/ache_longest_by_gene_aligned.fasta
#
# Notes
# -----
# * Uses MAFFT_PATH from config.R
# * Requires MAFFT to be installed and the path correctly
#   specified in config.R.
############################################################

### 0. Load config -----------------------------------------------------------

if (!file.exists("config.R")) {
  stop("config.R not found in project root.\n",
       "Create config.R with MAFFT_PATH defined.")
}
source("config.R")

if (!exists("MAFFT_PATH") || !nzchar(MAFFT_PATH)) {
  stop("MAFFT_PATH not defined in config.R")
}

### 1. Input / output paths --------------------------------------------------

# Prefer BLAST-filtered ace1/ace2 FASTA if it was created
blast_fasta   <- file.path("data", "combined",
                           "ache_longest_by_gene_ace12_blast_only.fasta")
default_fasta <- file.path("data", "combined",
                           "ache_longest_by_gene_all_species.fasta")

if (file.exists(blast_fasta)) {
  in_fasta <- blast_fasta
  cat("Using BLAST-filtered FASTA as input:\n  ",
      normalizePath(in_fasta), "\n")
} else {
  in_fasta <- default_fasta
  cat("BLAST-filtered FASTA not found.\n",
      "Using full combined FASTA as input:\n  ",
      normalizePath(in_fasta), "\n")
}

align_dir <- file.path("data", "alignment")
if (!dir.exists(align_dir)) {
  dir.create(align_dir, recursive = TRUE)
}

out_fasta <- file.path(align_dir,
                       "ache_longest_by_gene_aligned.fasta")

cat("Output alignment will be written to:\n  ",
    normalizePath(out_fasta), "\n\n")

### 2. Check input exists ----------------------------------------------------

if (!file.exists(in_fasta)) {
  stop("Input FASTA not found at: ", in_fasta,
       "\nRun 03_merge_longest_fastas.R (and optionally BLAST) first.")
}

### 3. Run MAFFT (streaming directly to file) --------------------------------

cmd_args <- c(
  "--auto",
  in_fasta
)

cat("Running MAFFT...\n")
cat("Command: ", MAFFT_PATH, " ",
    paste(cmd_args, collapse = " "), "\n\n", sep = "")

# Stream MAFFT output directly to the alignment file
exit_code <- system2(
  command = MAFFT_PATH,
  args    = cmd_args,
  stdout  = out_fasta,
  stderr  = ""
)

if (!identical(exit_code, 0L)) {
  stop("MAFFT finished with a non-zero exit code: ", exit_code)
}

### 4. Done ------------------------------------------------------------------

cat("Alignment complete.\n")
cat("Aligned FASTA saved to:\n  ",
    normalizePath(out_fasta), "\n")
