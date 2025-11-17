############################################################
# 05_build_tree_fasttree.R
#
# Goal
# ----
# Use FastTree to build a phylogenetic tree from
# the MAFFT alignment and save it as Newick.
#
# Input:
#   data/alignment/AChE_alignment_mafft.fasta
#
# Output:
#   data/tree/AChE_tree_fasttree.nwk
#
# Assumes:
#   * You are running from the project root (AChE_Project).
#   * FastTree is installed and either:
#       - available on PATH as "FastTree", or
#       - you edit fasttree_path below.
############################################################

### 1. Paths -----------------------------------------------------------------

align_fasta <- file.path("data", "alignment", "AChE_alignment_mafft.fasta")
tree_dir    <- file.path("data", "tree")

if (!file.exists(align_fasta)) {
  stop(
    "Alignment file not found at: ", align_fasta, "\n",
    "Run 04_run_alignment.R first."
  )
}

if (!dir.exists(tree_dir)) {
  dir.create(tree_dir, recursive = TRUE)
}

tree_file <- file.path(tree_dir, "AChE_tree_fasttree.nwk")

cat("Alignment file:", normalizePath(align_fasta), "\n")
cat("Tree output:   ", normalizePath(tree_file),   "\n")

### 2. FastTree executable path ---------------------------------------------

# Option 1 (recommended): assume FastTree is on PATH as "FastTree" or "fasttree".
fasttree_path <- "FastTree"

# Option 2: hard-code a local FastTree path (uncomment and edit if needed).
# fasttree_path <- "path/to/FastTree.exe"      # Windows example
# fasttree_path <- "/usr/local/bin/FastTree"   # macOS / Linux example

fasttree_found <- Sys.which(fasttree_path)

if (fasttree_found == "") {
  warning(
    "FastTree executable not found for '", fasttree_path, "'.\n",
    "Make sure FastTree is installed and on your system PATH,\n",
    "or edit fasttree_path in 05_build_tree_fasttree.R to point to it."
  )
} else {
  cat("Using FastTree at:", fasttree_found, "\n")
}

### 3. Build FastTree arguments ---------------------------------------------
# -nt : nucleotide sequences
# -gtr: GTR model (reasonable default)

fasttree_args <- c(
  "-nt",
  "-gtr",
  align_fasta
)

cat("Running FastTree...\n")

### 4. Run FastTree via system2 and capture output ---------------------------

tree_output <- system2(
  command = fasttree_path,
  args    = fasttree_args,
  stdout  = TRUE,   # capture all output lines
  stderr  = ""      # show log in console
)

cat("FastTree finished. Total output lines:",
    length(tree_output), "\n")

### 5. Extract only the Newick tree line ------------------------------------

# FastTree's actual tree line starts with "("
tree_line <- tree_output[grepl("^\\(", tree_output)][1]

if (is.na(tree_line)) {
  stop("Could not find a Newick line starting with '(' in FastTree output.")
}

### 6. Write the clean tree to file -----------------------------------------

writeLines(tree_line, con = tree_file)

cat("Done.\nTree saved to:\n", normalizePath(tree_file), "\n")
