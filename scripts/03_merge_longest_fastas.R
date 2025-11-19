############################################################
# 03_merge_longest_fastas.R
#
# Goal
# ----
# Take the per-species FASTA files that already contain the
# **longest transcript per gene** (from 02_extract_longest_isoforms.R)
# and merge them into a single multi-species FASTA.
#
# Outputs
# -------
# 1) Combined FASTA:
#      data/combined/ache_longest_by_gene_all_species.fasta
#
#    Headers are prefixed with the species short name, e.g.:
#      >gregaria|XM_012345678.1 Schistocerca gregaria ...
#
# 2) Metadata table (CSV) for downstream summaries:
#      data/combined/ache_longest_by_gene_metadata.csv
#
#    Columns:
#      species_short, header, gene_id, seq_len
#
# Notes
# -----
# * Must be run after 02_extract_longest_isoforms.R.
# * Assumes input files:
#      data/longest/<short_name>_longest_by_gene.fasta
############################################################

### 0. Project folders --------------------------------------------------------

longest_dir  <- file.path("data", "longest")
combined_dir <- file.path("data", "combined")

if (!dir.exists(longest_dir)) {
  stop("Directory 'data/longest' does not exist. ",
       "Run 02_extract_longest_isoforms.R first.")
}

if (!dir.exists(combined_dir)) {
  dir.create(combined_dir, recursive = TRUE)
}

cat("Longest-isoform directory:", normalizePath(longest_dir),  "\n")
cat("Combined output directory:", normalizePath(combined_dir), "\n\n")

### 1. Species list (must match scripts 01 & 02) ------------------------------

species_list <- list(
  list(short_name = "gregaria",   organism_term = "Schistocerca gregaria"),
  list(short_name = "cancellata", organism_term = "Schistocerca cancellata"),
  list(short_name = "piceifrons", organism_term = "Schistocerca piceifrons"),
  list(short_name = "anabrus",    organism_term = "Anabrus simplex"),
  list(short_name = "anopheles",  organism_term = "Anopheles gambiae"),
  list(short_name = "aedes",      organism_term = "Aedes aegypti"),
  list(short_name = "culex",      organism_term = "Culex quinquefasciatus")
)

### 2. Small helpers (FASTA + gene ID) ---------------------------------------

read_fasta <- function(path) {
  # Simple FASTA reader that returns a data frame with:
  #   header  = full header line (no ">")
  #   seq     = sequence (one string per record)
  lines <- readLines(path)
  is_header <- startsWith(lines, ">")
  
  if (!any(is_header)) {
    stop("No FASTA headers (>) detected in file: ", path)
  }
  
  header_idx <- which(is_header)
  header_idx_end <- c(header_idx[-1] - 1, length(lines))
  
  headers <- character(length(header_idx))
  seqs    <- character(length(header_idx))
  
  for (i in seq_along(header_idx)) {
    h <- lines[header_idx[i]]
    s <- lines[(header_idx[i] + 1):header_idx_end[i]]
    headers[i] <- sub("^>", "", h)
    seqs[i]    <- paste0(s, collapse = "")
  }
  
  data.frame(
    header = headers,
    seq    = seqs,
    stringsAsFactors = FALSE
  )
}

get_gene_id_from_header <- function(header) {
  # Same logic as in 02_extract_longest_isoforms.R
  # to keep gene IDs consistent.
  
  # 1) gene= pattern
  gene_match <- regexpr("gene=([^ ]+)", header)
  if (gene_match != -1) {
    gene <- regmatches(header, gene_match)
    gene <- sub("^gene=", "", gene)
    return(gene)
  }
  
  # 2) LOC-style IDs
  loc_match <- regexpr("LOC[0-9]+", header)
  if (loc_match != -1) {
    loc_id <- regmatches(header, loc_match)
    return(loc_id)
  }
  
  # 3) Fallback: accession without version
  first_token <- strsplit(header, " ")[[1]][1]
  accession_root <- sub("\\.[0-9]+$", "", first_token)
  
  return(accession_root)
}

write_fasta <- function(fasta_df, path) {
  lines <- character(2 * nrow(fasta_df))
  j <- 1
  for (i in seq_len(nrow(fasta_df))) {
    lines[j]   <- paste0(">", fasta_df$header[i])
    lines[j+1] <- fasta_df$seq[i]
    j <- j + 2
  }
  writeLines(lines, con = path)
}

### 3. Main: combine all species ---------------------------------------------

combined_fasta   <- list()
combined_meta_df <- list()

cat("Merging longest-by-gene FASTAs into a combined file...\n\n")

for (sp in species_list) {
  short_name <- sp$short_name
  
  in_file <- file.path(
    longest_dir,
    paste0(short_name, "_longest_by_gene.fasta")
  )
  
  cat("--------------------------------------------------\n")
  cat("Species:", short_name, "\n")
  cat("Input:", in_file, "\n")
  
  if (!file.exists(in_file)) {
    warning("Longest-by-gene FASTA not found for ", short_name,
            ". Run 02_extract_longest_isoforms.R first.")
    next
  }
  
  # Read FASTA for this species
  df <- read_fasta(in_file)
  cat("  Records:", nrow(df), "\n")
  
  # Build metadata: species_short, gene_id, seq_len
  df$species_short <- short_name
  df$gene_id       <- vapply(
    df$header,
    get_gene_id_from_header,
    FUN.VALUE = character(1)
  )
  df$seq_len       <- nchar(df$seq)
  
  combined_meta_df[[short_name]] <- df[, c(
    "species_short", "header", "gene_id", "seq_len"
  )]
  
  # For the combined FASTA, prefix headers with species short name
  # so it's always clear where each sequence came from.
  df$header <- paste0(short_name, "|", df$header)
  
  combined_fasta[[short_name]] <- df[, c("header", "seq")]
}

# Bind all species together
combined_fasta_df <- do.call(rbind, combined_fasta)
combined_meta     <- do.call(rbind, combined_meta_df)

# Output paths
out_fasta <- file.path(
  combined_dir,
  "ache_longest_by_gene_all_species.fasta"
)
out_meta  <- file.path(
  combined_dir,
  "ache_longest_by_gene_metadata.csv"
)

# Write combined FASTA and CSV
write_fasta(combined_fasta_df, out_fasta)
write.csv(combined_meta, out_meta, row.names = FALSE)

cat("--------------------------------------------------\n")
cat("Wrote combined FASTA to:\n  ", normalizePath(out_fasta), "\n")
cat("Wrote metadata CSV to:\n  ", normalizePath(out_meta),  "\n")
cat("Done.\n")
