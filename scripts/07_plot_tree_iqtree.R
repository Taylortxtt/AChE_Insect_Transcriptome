############################################################
# 06b_plot_tree_iqtree.R
#
# Goal
# ----
# Plot the IQ-TREE phylogeny and save it as a PNG.
#
# Input
# -----
#   data/tree/ache_longest_by_gene_iqtree.treefile
#
# Output
# ------
#   data/tree/ache_longest_by_gene_iqtree_plot.png
#
# Notes
# -----
# * Assumes 05b_build_tree_iqtree.R has been run.
# * Uses the same color scheme as 06_plot_tree.R.
############################################################

### 0. Packages -------------------------------------------------------------

if (!requireNamespace("ape", quietly = TRUE)) {
  install.packages("ape")
}
library(ape)

### 1. Paths ---------------------------------------------------------------

tree_dir <- file.path("data", "tree")

if (!dir.exists(tree_dir)) {
  stop("Directory 'data/tree' does not exist. ",
       "Run 05b_build_tree_iqtree.R first.")
}

tree_file <- file.path(
  tree_dir,
  "ache_longest_by_gene_iqtree.treefile"
)

png_file <- file.path(
  tree_dir,
  "ache_longest_by_gene_iqtree_plot.png"
)

cat("IQ-TREE file: ", normalizePath(tree_file), "\n")
cat("PNG file:     ", normalizePath(png_file),  "\n\n")

if (!file.exists(tree_file)) {
  stop("IQ-TREE tree file not found at: ", tree_file,
       "\nRun 05b_build_tree_iqtree.R first.")
}

### 2. Read tree -----------------------------------------------------------

tree <- read.tree(tree_file)

cat("Number of tips in IQ-TREE tree:", length(tree$tip.label), "\n\n")

### 3. Define groups & colors ----------------------------------------------

# Tip labels look like:
#   shortname|XM_...
# We use the shortname (before "|") to assign groups.

get_short_name <- function(tip_label) {
  strsplit(tip_label, "\\|")[[1]][1]
}

short_names <- vapply(tree$tip.label, get_short_name, character(1))

# Map short names to broader groups for coloring
group <- ifelse(short_names %in% c("gregaria", "cancellata", "piceifrons"),
                "locust",
                ifelse(short_names %in% c("anabrus"),
                       "cricket",
                       ifelse(short_names %in% c("anopheles", "aedes", "culex"),
                              "mosquito",
                              "other")))

group_factor <- factor(group,
                       levels = c("locust", "cricket", "mosquito", "other"))

group_colors <- c(
  locust   = "forestgreen",
  cricket  = "goldenrod",
  mosquito = "firebrick",
  other    = "gray50"
)

tip_cols <- group_colors[as.character(group_factor)]

### 4. Plot and save --------------------------------------------------------

png(
  filename = png_file,
  width  = 2000,
  height = 2000,
  res    = 300
)

par(mar = c(2, 2, 2, 2))

plot(
  tree,
  show.tip.label = TRUE,
  cex            = 0.4,
  tip.color      = tip_cols
)
title("AChE Longest-Per-Gene Tree (IQ-TREE)")

legend(
  "topleft",
  legend = levels(group_factor),
  col    = group_colors[levels(group_factor)],
  pch    = 19,
  bty    = "n",
  cex    = 0.7,
  title  = "Species group"
)

dev.off()

cat("IQ-TREE plot saved to:\n  ", normalizePath(png_file), "\n")
