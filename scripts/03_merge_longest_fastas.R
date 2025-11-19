############################################################
# 03_merge_longest_fastas.R
#
# Merge all species' longest-by-gene FASTAs into:
#   data/combined/ache_longest_by_gene_all_species.fasta
# and metadata:
#   data/combined/ache_longest_by_gene_metadata.csv
############################################################

library(stringr)
library(dplyr)

gene_dir <- file.path("data", "longest_by_gene")
out_dir  <- file.path("data", "combined")

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

out_fasta <- file.path(out_dir, "ache_longest_by_gene_all_species.fasta")
out_meta  <- file.path(out_dir, "ache_longest_by_gene_metadata.csv")

# Start fresh
if (file.exists(out_fasta)) file.remove(out_fasta)
if (file.exists(out_meta))  file.remove(out_meta)

cat("Longest-by-gene directory:", normalizePath(gene_dir), "\n")
cat("Combined FASTA:           ", normalizePath(out_fasta), "\n")
cat("Metadata CSV:             ", normalizePath(out_meta),  "\n\n")

# ---- Helper: read FASTA -----------------------------------------------
read_fasta <- function(path) {
  lines   <- readLines(path)
  headers <- grep("^>", lines)
  if (length(headers) == 0) {
    stop("No FASTA headers (>) found in file: ", path)
  }
  starts  <- headers + 1
  ends    <- c(headers[-1] - 1, length(lines))
  
  hdr <- sub("^>", "", lines[headers])
  seq <- character(length(headers))
  for (i in seq_along(headers)) {
    seq[i] <- paste(lines[starts[i]:ends[i]], collapse = "")
  }
  
  data.frame(
    header = hdr,
    seq    = seq,
    stringsAsFactors = FALSE
  )
}

# ---- Main --------------------------------------------------------------
gene_files <- list.files(
  gene_dir,
  pattern = "_longest_by_gene\\.fasta$",
  full.names = TRUE
)

if (length(gene_files) == 0) {
  stop("No longest_by_gene FASTAs found in ", gene_dir)
}

all_meta <- list()

for (f in gene_files) {
  sp <- sub("_longest_by_gene\\.fasta$", "", basename(f))
  cat("Species:", sp, "\n")
  cat("  Input FASTA:", f, "\n")
  
  df <- read_fasta(f)
  
  # qseqid is what BLAST will see (first token of header)
  df$qseqid <- sub(" .*", "", df$header)
  df$species <- sp
  df$len <- nchar(df$seq)
  
  cat("  Sequences:", nrow(df), "\n")
  
  # Append to combined FASTA (using full header for readability)
  con <- file(out_fasta, open = if (file.exists(out_fasta)) "a" else "w")
  for (i in seq_len(nrow(df))) {
    writeLines(paste0(">", df$header[i]), con)
    writeLines(df$seq[i], con)
  }
  close(con)
  
  all_meta[[length(all_meta) + 1]] <- df[, c("species", "qseqid", "header", "len")]
}

meta_df <- bind_rows(all_meta)

write.csv(meta_df, out_meta, row.names = FALSE)

cat("\nCombined FASTA written to:\n  ", normalizePath(out_fasta), "\n")
cat("Metadata CSV written to:\n  ", normalizePath(out_meta),  "\n")
cat("Done merging.\n")
