############################################################
# 03_merge_longest_fastas.R
#
# Goal
# ----
# Merge all *_longest.fasta files in data/longest/ into
# a single FASTA file for downstream alignment and tree
# building.
#
# Output:
#   data/combined/AChE_longest_all_species.fasta
############################################################

### 1. Directories -----------------------------------------------------------

longest_dir  <- file.path("data", "longest")
combined_dir <- file.path("data", "combined")

if (!dir.exists(combined_dir)) {
  dir.create(combined_dir, recursive = TRUE)
}

cat("Longest dir:  ", normalizePath(longest_dir),  "\n")
cat("Combined dir: ", normalizePath(combined_dir), "\n")

### 2. Find all longest FASTA files -----------------------------------------

longest_files <- list.files(
  longest_dir,
  pattern    = "_longest\\.fasta$",
  full.names = TRUE
)

cat("Found", length(longest_files), "longest FASTA files:\n")
print(longest_files)

if (length(longest_files) == 0) {
  stop("No *_longest.fasta files found in ", longest_dir)
}

### 3. Read and concatenate --------------------------------------------------

all_lines <- character(0)

for (f in longest_files) {
  cat("Adding file:", f, "\n")
  lines <- readLines(f)
  all_lines <- c(all_lines, lines)
}

cat("Total lines in merged FASTA:", length(all_lines), "\n")

### 4. Write combined FASTA --------------------------------------------------

out_path <- file.path(combined_dir, "AChE_longest_all_species.fasta")

writeLines(all_lines, con = out_path)

cat("Saved combined FASTA to:", out_path, "\n")
