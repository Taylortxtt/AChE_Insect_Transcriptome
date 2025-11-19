############################################################
# 06_plot_tree.R
#
# Goal
# ----
# Plot the FastTree phylogeny and save:
#   1) A clean overview tree (no tip labels).
#   2) A detailed tree with tiny tip labels (for zooming).
#
# Input
# -----
#   data/tree/ache_longest_by_gene_fasttree.nwk
#
# Output
# ------
#   data/tree/ache_longest_by_gene_tree_overview.png
#   data/tree/ache_longest_by_gene_tree_labeled.png
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
       "Run 05_build_tree_fasttree.R first.")
}

tree_file <- file.path(
  tree_dir,
  "ache_longest_by_gene_fasttree.nwk"
)

png_overview <- file.path(
  tree_dir,
  "ache_longest_by_gene_tree_overview.png"
)

png_labeled <- file.path(
  tree_dir,
  "ache_longest_by_gene_tree_labeled.png"
)
cat("FastTree file: ", normalizePath(tree_file), "\n")
cat("Overview PNG:  ", png_overview, "\n")
cat("Labeled PNG:   ", png_labeled,  "\n\n")

if (!file.exists(tree_file)) {
  stop("FastTree tree file not found at: ", tree_file,
       "\nRun 05_build_tree_fasttree.R first.")
}

### 2. Read tree -----------------------------------------------------------

tree <- read.tree(tree_file)

cat("Number of tips in FastTree tree:", length(tree$tip.label), "\n\n")

### 3. Define groups & colors ----------------------------------------------

get_short_name <- function(tip_label) {
  strsplit(tip_label, "\\|")[[1]][1]
}

short_names <- vapply(tree$tip.label, get_short_name, character(1))

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

### 4. Overview tree (no tip labels) ---------------------------------------

png(
  filename = png_overview,
  width  = 2000,
  height = 2000,
  res    = 300
)

par(mar = c(2, 2, 2, 2))
plot(
  tree,
  show.tip.label = FALSE,  # key change: no labels
  cex            = 0.4,
  tip.color      = tip_cols
)
title("AChE Longest-Per-Gene Tree (FastTree) – Overview")

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

cat("Overview FastTree plot saved to:\n  ",
    normalizePath(png_overview), "\n\n")

### 5. Detailed tree (tiny labels, for zooming) -----------------------------

png(
  filename = png_labeled,
  width  = 3000,
  height = 3000,
  res    = 400
)

par(mar = c(2, 2, 2, 2))
plot(
  tree,
  show.tip.label = TRUE,
  cex            = 0.2,    # much smaller labels
  tip.color      = tip_cols
)
title("AChE Longest-Per-Gene Tree (FastTree) – Labeled")

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

cat("Labeled FastTree plot saved to:\n  ",
    normalizePath(png_labeled), "\n")
