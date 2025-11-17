############################################################
# 04_run_alignment.R
#
# Goal
# ----
# Run MAFFT on the combined longest-isoform FASTA and
# save the multiple sequence alignment to:
#
#   data/alignment/AChE_alignment_mafft.fasta
#
# This script assumes:
#   * You are running from the project root (AChE_Project).
#   * MAFFT is installed and either:
#       - available on your system PATH as "mafft", or
#       - you edit mafft_path below to point to mafft/mafft.bat.
############################################################

### 1. Paths ---------------------------------------------------------------

combined_fasta <- file.path("data", "combined", "AChE_longest_all_species.fasta")
align_dir      <- file.path("data", "alignment")

if (!file.exists(combined_fasta)) {
  stop(
    "Combined FASTA not found at: ", combined_fasta, "\n",
    "Run 03_merge_longest_fastas.R first."
  )
}

if (!dir.exists(align_dir)) {
  dir.create(align_dir, recursive = TRUE)
}

aligned_fasta <- file.path(align_dir, "AChE_alignment_mafft.fasta")

cat("Input FASTA:  ", normalizePath(combined_fasta), "\n")
cat("Output FASTA: ", normalizePath(aligned_fasta),  "\n")

### 2. MAFFT executable path -----------------------------------------------

# Option 1 (recommended): assume MAFFT is on PATH as "mafft".
mafft_path <- "mafft"

# Option 2: hard-code a local MAFFT path (uncomment and edit if needed).
# mafft_path <- "path/to/mafft.bat"        # Windows example
# mafft_path <- "/usr/local/bin/mafft"     # macOS / Linux example

mafft_found <- Sys.which(mafft_path)

if (mafft_found == "") {
  warning(
    "MAFFT executable not found for '", mafft_path, "'.\n",
    "Make sure MAFFT is installed and on your system PATH,\n",
    "or edit mafft_path in 04_run_alignment.R to point to it."
  )
} else {
  cat("Using MAFFT at:", mafft_found, "\n")
}

### 3. Build MAFFT command --------------------------------------------------

# --auto     : MAFFT chooses a sensible strategy
# --thread 2 : use 2 CPU threads (edit if needed)
mafft_args <- c(
  "--auto",
  "--thread", "2",
  combined_fasta
)

cat("Running MAFFT...\n")

### 4. Run MAFFT via system2 -----------------------------------------------

# Capture MAFFT's stdout and write it directly to the alignment file.
system2(
  command = mafft_path,
  args    = mafft_args,
  stdout  = aligned_fasta,
  stderr  = ""   # "" = send stderr to the R console
)

cat("Done.\nAlignment saved to:\n", normalizePath(aligned_fasta), "\n")
