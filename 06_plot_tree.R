############################################################
# 06_plot_tree.R
#
# Goal
# ----
# Load the FastTree Newick tree and plot it with
# common-name tip labels only, plus simple color coding.
#
# Notes
# -----
# * The underlying .nwk file is not modified.
# * Only the in-memory tree object is relabeled for plotting.
############################################################

library(ape)

### 1. Load the tree ---------------------------------------------------------

tree_dir  <- file.path("data", "tree")
tree_path <- file.path(tree_dir, "AChE_tree_fasttree.nwk")

if (!file.exists(tree_path)) {
  stop("Tree file not found at: ", tree_path)
}

tree <- read.tree(tree_path)
cat("Loaded tree with", length(tree$tip.label), "tips.\n")

### 2. Use common names as tip labels ----------------------------------------
# Tree tip labels look like:
#   gregaria_AChE_XM_049998832
#   aedes_AChE_XM_021851332.1
#
# For the figure we only want the common name:
#   gregaria, aedes, anabrus, etc.
# → take everything before the first underscore.

original_labels <- tree$tip.label
common_names    <- sub("_.*$", "", original_labels)

cat("Original tip labels:\n")
print(original_labels)

cat("Common-name labels used for plotting:\n")
print(common_names)

tree$tip.label <- common_names
tips <- tree$tip.label

### 3. Define groups and colors ----------------------------------------------

locusts  <- c("gregaria", "cancellata", "piceifrons")
cricket  <- c("anabrus")                       # Mormon cricket
mosquito <- c("anopheles", "aedes", "culex")

tip_group <- ifelse(tips %in% locusts,  "locust",
                    ifelse(tips %in% cricket,  "cricket",
                           ifelse(tips %in% mosquito, "mosquito", "other")))

tip_col <- ifelse(tip_group == "locust",   "forestgreen",
                  ifelse(tip_group == "cricket",  "steelblue",
                         ifelse(tip_group == "mosquito", "darkorange", "black")))

### 4. Save PNG to data/tree/ -----------------------------------------------

fig_path <- file.path(tree_dir, "AChE_tree_fasttree_clean.png")

png(fig_path, width = 960, height = 768, res = 120)

# Slightly larger outer top margin so the title is not clipped
par(
  mar = c(5, 4, 4, 2) + 0.1,   # inner margins
  oma = c(0, 0, 3, 0)          # outer margins (top used for title)
)

plot(
  tree,
  type            = "phylogram",
  use.edge.length = TRUE,
  tip.color       = tip_col,
  cex             = 0.9,
  font            = 2,
  no.margin       = FALSE,
  main            = ""          # title added with mtext() below
)

axisPhylo()

# Title in the outer margin so it stays fully visible
mtext(
  "AChE FastTree Phylogeny (Longest Isoforms)",
  side  = 3,
  outer = TRUE,
  line  = 1.2,
  cex   = 1.6,
  font  = 2
)

dev.off()

cat("Tree figure saved to:", fig_path, "\n")
