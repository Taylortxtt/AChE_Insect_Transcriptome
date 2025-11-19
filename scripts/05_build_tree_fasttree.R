############################################################
# 05_build_tree_fasttree.R
#
# Goal
# ----
# Build a phylogenetic tree from the MAFFT alignment using
# FastTree and save it as a Newick file.
#
# Input
# -----
#   data/alignment/ache_longest_by_gene_aligned.fasta
#
# Output
# ------
#   data/tree/ache_longest_by_gene_fasttree.nwk
#
# Notes
# -----
# * Uses FASTTREE_PATH from config.R
# * Requires FastTree (or FastTreeMP) to be installed and
#   the path correctly specified in config.R.
############################################################

### 0. Load config -----------------------------------------------------------

if (!file.exists("config.R")) {
  stop("config.R not found in project root.\n",
       "Create config.R with FASTTREE_PATH defined.")
}

source("config.R")   # loads FASTTREE_PATH

if (!exists("FASTTREE_PATH")) {
  stop("FASTTREE_PATH not defined in config.R.")
}

### 1. Paths -----------------------------------------------------------------

align_dir <- file.path("data", "alignment")
tree_dir  <- file.path("data", "tree")

if (!dir.exists(align_dir)) {
  stop("Directory 'data/alignment' does not exist. ",
       "Run 04_run_alignment.R first.")
}
if (!dir.exists(tree_dir)) {
  dir.create(tree_dir, recursive = TRUE)
}

in_align <- file.path(
  align_dir,
  "ache_longest_by_gene_aligned.fasta"
)

out_tree <- file.path(
  tree_dir,
  "ache_longest_by_gene_fasttree.nwk"
)

cat("Input alignment: ", normalizePath(in_align), "\n")
cat("Output tree:     ", normalizePath(out_tree), "\n\n")

if (!file.exists(in_align)) {
  stop("Aligned FASTA not found at: ", in_align,
       "\nRun 04_run_alignment.R first.")
}

### 2. Verify FastTree path --------------------------------------------------

fasttree_path <- FASTTREE_PATH

if (!file.exists(fasttree_path)) {
  stop(
    "FastTree not found at path: ", fasttree_path, "\n",
    "Edit FASTTREE_PATH in config.R to point to your FastTree executable.\n",
    "Example paths:\n",
    "  Mac (compiled): /usr/local/bin/FastTree\n"
  )
}

cat("Using FastTree at: ", fasttree_path, "\n\n")

### 3. Run FastTree ----------------------------------------------------------

# Recommended model: GTR + Gamma
cmd_args <- c(
  "-gtr",
  "-gamma",
  in_align
)

cat("Running FastTree...\n")
cat("Command: ", fasttree_path, " ", paste(cmd_args, collapse = " "), "\n\n")

# Stream FastTree output directly to the .nwk file
exit_code <- system2(
  command = fasttree_path,
  args    = cmd_args,
  stdout  = out_tree,  # write tree directly to file
  stderr  = ""         # progress messages go to console
)

if (!identical(exit_code, 0L)) {
  stop("FastTree finished with a non-zero exit code: ", exit_code)
}

### 4. Done ------------------------------------------------------------------

cat("Tree building complete.\n")
cat("Tree saved to:\n  ", normalizePath(out_tree), "\n")
