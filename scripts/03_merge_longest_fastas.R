############################################################
# 03_merge_longest_fastas.R
#
<<<<<<< HEAD
# Merges all species' "longest-by-gene" FASTAs into a single
# multi-species FASTA, and builds a metadata table.
#
# Inputs:
#   data/longest_by_gene/<species>_longest_by_gene.fasta
#
# Outputs:
#   data/combined/ache_longest_by_gene_all_species.fasta
=======
# Merge all species' longest-by-gene FASTAs into:
#   data/combined/ache_longest_by_gene_all_species.fasta
# and metadata:
>>>>>>> 600b4705d5255729d9071965d3f21768e31ffdc9
#   data/combined/ache_longest_by_gene_metadata.csv
############################################################

library(stringr)
library(dplyr)
library(readr)

<<<<<<< HEAD
# 1. Directories -----------------------------------------------------------

by_gene_dir   <- file.path("data", "longest_by_gene")
combined_dir  <- file.path("data", "combined")

if (!dir.exists(by_gene_dir)) {
  stop("No longest_by_gene directory found: ", by_gene_dir)
}
if (!dir.exists(combined_dir)) {
  dir.create(combined_dir, recursive = TRUE)
}

out_fasta <- file.path(combined_dir, "ache_longest_by_gene_all_species.fasta")
out_meta  <- file.path(combined_dir, "ache_longest_by_gene_metadata.csv")

# 2. Helper to read a FASTA into a data.frame -----------------------------

read_fasta_df <- function(path) {
  lines   <- readLines(path)
  header_i <- grep("^>", lines)
  if (length(header_i) == 0) {
    return(data.frame(header = character(), seq = character(), stringsAsFactors = FALSE))
  }
  
  starts <- header_i + 1
  ends   <- c(header_i[-1] - 1, length(lines))
  
  df <- data.frame(
    header = character(length(header_i)),
    seq    = character(length(header_i)),
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(header_i)) {
    df$header[i] <- lines[header_i[i]]
    df$seq[i]    <- paste(lines[starts[i]:ends[i]], collapse = "")
  }
  
  df
}

# 3. Loop over species files ----------------------------------------------

by_gene_files <- list.files(by_gene_dir, pattern = "_longest_by_gene\\.fasta$", full.names = TRUE)

if (length(by_gene_files) == 0) {
  stop("No longest_by_gene FASTAs found in ", by_gene_dir)
}

cat("Found", length(by_gene_files), "longest-by-gene FASTA files in", by_gene_dir, "\n\n")

# Start clean output FASTA
if (file.exists(out_fasta)) file.remove(out_fasta)

all_meta <- list()
first_write <- TRUE

for (f in by_gene_files) {
  sp <- basename(f)
  sp <- sub("_longest_by_gene.fasta$", "", sp)
  cat("Species:", sp, "\n")
  cat("  Input FASTA:", f, "\n")
  
  df <- read_fasta_df(f)
  n  <- nrow(df)
  cat("  Sequences:", n, "\n\n")
  
  if (n == 0) next
  
  # qseqid = what BLAST will see (first token after '>')
  df$qseqid  <- sub("^>", "", df$header)
  df$qseqid  <- sub(" .*", "", df$qseqid)
  
  # species column
  df$species <- sp
  
  # Append to combined FASTA
  con <- file(out_fasta, open = if (first_write) "w" else "a")
  for (i in seq_len(nrow(df))) {
    writeLines(df$header[i], con)
    writeLines(df$seq[i],    con)
  }
  close(con)
  first_write <- FALSE
  
  # Build metadata chunk
  meta_chunk <- df %>%
    transmute(
      species = species,
      qseqid  = qseqid,
      header  = header,
      seq_len = nchar(seq)
    )
  
  all_meta[[length(all_meta) + 1]] <- meta_chunk
}

# 4. Write metadata -------------------------------------------------------

meta_df <- bind_rows(all_meta)

write_csv(meta_df, out_meta)

cat("Combined FASTA written to:\n  ", out_fasta, "\n")
cat("Metadata CSV written to:\n  ", out_meta,  "\n")
=======
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
>>>>>>> 600b4705d5255729d9071965d3f21768e31ffdc9
cat("Done merging.\n")
cat("Next step: run scripts/04_classify_genes_blast.R\n")
