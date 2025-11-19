############################################################
# 05b_build_tree_iqtree.R
#
# Goal
# ----
# Build a phylogenetic tree from the MAFFT alignment using
# IQ-TREE and save it (and its support values) in standard
# output files.
#
# Input
# -----
#   data/alignment/ache_longest_by_gene_aligned.fasta
#
# Outputs (all in data/tree/)
# ---------------------------
#   ache_longest_by_gene_iqtree.treefile   # main Newick tree
#   ache_longest_by_gene_iqtree.log        # run log
#   ache_longest_by_gene_iqtree.iqtree     # detailed report
#   ache_longest_by_gene_iqtree.*          # other IQ-TREE files
#
# Notes
# -----
# * Uses IQTREE_PATH from config.R
# * Requires IQ-TREE (iqtree or iqtree2) to be installed.
# * This is an extra, slower but more thorough alternative
#   to the FastTree tree in 05_build_tree_fasttree.R.
############################################################

### 0. Load config -----------------------------------------------------------

if (!file.exists("config.R")) {
  stop("config.R not found in project root.\n",
       "Create config.R with IQTREE_PATH defined.")
}

source("config.R")   # loads IQTREE_PATH

if (!exists("IQTREE_PATH")) {
  stop("IQTREE_PATH not defined in config.R.")
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

# IQ-TREE uses a prefix instead of a single output file path
# All output files will start with this prefix:
out_prefix <- file.path(
  tree_dir,
  "ache_longest_by_gene_iqtree"
)

cat("Input alignment:  ", normalizePath(in_align),   "\n")
cat("Output prefix:    ", normalizePath(out_prefix), "\n\n")

if (!file.exists(in_align)) {
  stop("Aligned FASTA not found at: ", in_align,
       "\nRun 04_run_alignment.R first.")
}

### 2. Verify IQ-TREE path ---------------------------------------------------

iqtree_path <- IQTREE_PATH

if (!file.exists(iqtree_path)) {
  stop(
    "IQ-TREE not found at path: ", iqtree_path, "\n",
    "Edit IQTREE_PATH in config.R to point to your IQ-TREE executable.\n",
    "Example paths:\n",
    "  Mac (brew):   /opt/homebrew/bin/iqtree2\n",
    "  Linux:        /usr/bin/iqtree2\n",
    "  Windows:      C:/Program Files/iqtree2/iqtree2.exe\n"
  )
}

cat("Using IQ-TREE at: ", iqtree_path, "\n\n")

### 3. Run IQ-TREE -----------------------------------------------------------

# Common, reasonable defaults:
#   -s      : alignment file
#   -m MFP  : model finder (automatically chooses best model)
#   -bb 1000: ultrafast bootstrap (1000 replicates)
#   -alrt 1000: SH-aLRT support (1000 replicates)
#   -nt AUTO: automatically decide number of threads
#   -pre    : prefix for output files
cmd_args <- c(
  "-s", in_align,
  "-m", "MFP",
  "-bb", "1000",
  "-alrt", "1000",
  "-nt", "AUTO",
  "-pre", out_prefix
)

cat("Running IQ-TREE...\n")
cat("Command: ", iqtree_path, " ", paste(cmd_args, collapse = " "), "\n\n")

# We don't capture stdout here because IQ-TREE writes its own files.
# We just run it and let it write into data/tree/.
exit_code <- system2(
  command = iqtree_path,
  args    = cmd_args
)

if (!identical(exit_code, 0L)) {
  stop("IQ-TREE finished with a non-zero exit code: ", exit_code)
}

### 4. Summary ---------------------------------------------------------------

treefile <- paste0(out_prefix, ".treefile")

if (file.exists(treefile)) {
  cat("IQ-TREE run complete.\n")
  cat("Main tree file:\n  ", normalizePath(treefile), "\n")
} else {
  warning("IQ-TREE run finished, but .treefile was not found.\n",
          "Check the log at:\n  ",
          normalizePath(paste0(out_prefix, ".log")))
}
