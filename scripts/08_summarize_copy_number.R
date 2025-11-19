############################################################
# 07_summarize_copy_number.R
#
# Goal
# ----
# Summarize AChE copy number (longest-per-gene) per species
# using the metadata table generated in 03_merge_longest_fastas.R.
#
# Input
# -----
#   data/combined/ache_longest_by_gene_metadata.csv
#
# Outputs
# -------
#   data/combined/ache_copy_number_summary.csv
#   exports/ache_copy_number_barplot.png   (if ggplot2 is available)
#
# Notes
# -----
# * Must be run after 03_merge_longest_fastas.R.
# * Each row in the metadata table represents the longest
#   transcript for one gene in one species.
############################################################

### 0. Packages (optional for plotting) --------------------------------------

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  message("Package 'ggplot2' not found. ",
          "Plots will be skipped unless you install it with install.packages('ggplot2').")
}

### 1. Paths -----------------------------------------------------------------

combined_dir <- file.path("data", "combined")
exports_dir  <- file.path("exports")

if (!dir.exists(combined_dir)) {
  stop("Directory 'data/combined' does not exist. ",
       "Run 03_merge_longest_fastas.R first.")
}

if (!dir.exists(exports_dir)) {
  dir.create(exports_dir, recursive = TRUE)
}

meta_file <- file.path(
  combined_dir,
  "ache_longest_by_gene_metadata.csv"
)

out_summary_csv <- file.path(
  combined_dir,
  "ache_copy_number_summary.csv"
)

out_barplot_png <- file.path(
  exports_dir,
  "ache_copy_number_barplot.png"
)

cat("Metadata file:  ", normalizePath(meta_file),        "\n")
cat("Summary CSV:    ", normalizePath(out_summary_csv),  "\n")
cat("Barplot (PNG):  ", normalizePath(out_barplot_png),  "\n\n")

if (!file.exists(meta_file)) {
  stop("Metadata file not found at: ", meta_file,
       "\nRun 03_merge_longest_fastas.R first.")
}

### 2. Read metadata ---------------------------------------------------------

meta <- read.csv(meta_file, stringsAsFactors = FALSE)

# Expecting at least: species_short, gene_id, header, seq_len
required_cols <- c("species_short", "gene_id")
missing_cols  <- setdiff(required_cols, names(meta))

if (length(missing_cols) > 0) {
  stop("Metadata file is missing required columns: ",
       paste(missing_cols, collapse = ", "))
}

cat("Rows in metadata: ", nrow(meta), "\n")

### 3. Calculate copy number per species ------------------------------------

# Just in case, deduplicate by species + gene_id
unique_genes <- unique(meta[, c("species_short", "gene_id")])

cat("Unique species/gene combinations: ", nrow(unique_genes), "\n\n")

# Count genes per species
copy_summary <- aggregate(
  gene_id ~ species_short,
  data = unique_genes,
  FUN  = length
)

names(copy_summary)[names(copy_summary) == "gene_id"] <- "n_genes"

cat("Copy-number summary:\n")
print(copy_summary)
cat("\n")

# Save summary CSV
write.csv(copy_summary, out_summary_csv, row.names = FALSE)

cat("Copy-number summary written to:\n  ",
    normalizePath(out_summary_csv), "\n\n")

### 4. Optional barplot ------------------------------------------------------

if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  
  p <- ggplot(copy_summary,
              aes(x = species_short, y = n_genes)) +
    geom_col() +
    theme_minimal() +
    labs(
      title = "AChE Copy Number per Species (Longest-Per-Gene)",
      x     = "Species (short name)",
      y     = "Number of AChE genes"
    )
  
  ggsave(
    filename = out_barplot_png,
    plot     = p,
    width    = 6,
    height   = 4,
    dpi      = 300
  )
  
  cat("Barplot saved to:\n  ",
      normalizePath(out_barplot_png), "\n")
} else {
  cat("Skipping barplot because 'ggplot2' is not installed.\n")
}
