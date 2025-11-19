############################################################
# 02_extract_longest_isoforms.R
#
# From each species' raw FASTA, keep the *longest transcript
# per gene* (no ace1/ace2 filtering here).
#
# Outputs, per species:
#   data/longest_isoforms/<sp>_longest_isoforms.fasta
#   data/longest_by_gene/<sp>_longest_by_gene.fasta
############################################################

library(stringr)
library(dplyr)

raw_dir  <- file.path("data", "raw")
iso_dir  <- file.path("data", "longest_isoforms")
gene_dir <- file.path("data", "longest_by_gene")

if (!dir.exists(iso_dir))  dir.create(iso_dir,  recursive = TRUE)
if (!dir.exists(gene_dir)) dir.create(gene_dir, recursive = TRUE)

cat("Raw FASTAs directory:     ", normalizePath(raw_dir),  "\n")
cat("Longest-isoform directory:", normalizePath(iso_dir),  "\n")
cat("Longest-by-gene directory:", normalizePath(gene_dir), "\n\n")

# Species short names (must match the *_raw.fasta files)
species_list <- c(
  "gregaria",
  "piceifrons",
  "nitens",
  "cancellata",
  "serialis",
  "anabrus",
  "anopheles",
  "aedes",
  "culex"
)

# ---- Helper: read FASTA into a data.frame -------------------------------
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

# ---- Helper: get gene ID from header -----------------------------------
# Strategy:
#   1) If there's a LOCxxxxx in the header, use that as gene_id.
#   2) Otherwise, use the first token (accession) before the first space.
get_gene_id <- function(header_vec) {
  loc <- str_extract(header_vec, "LOC\\d+")
  first_token <- sub(" .*", "", header_vec)
  ifelse(!is.na(loc), loc, first_token)
}

# ---- Main --------------------------------------------------------------
for (sp in species_list) {
  cat("Species:", sp, "\n")
  
  in_file   <- file.path(raw_dir,  paste0(sp, "_raw.fasta"))
  iso_file  <- file.path(iso_dir,  paste0(sp, "_longest_isoforms.fasta"))
  gene_file <- file.path(gene_dir, paste0(sp, "_longest_by_gene.fasta"))
  
  if (!file.exists(in_file)) {
    cat("  Raw FASTA missing – skipping.\n\n")
    next
  }
  
  df <- read_fasta(in_file)
  
  df$gene_id <- get_gene_id(df$header)
  df$len     <- nchar(df$seq)
  
  cat("  Raw sequences:   ", nrow(df), "\n")
  cat("  Unique gene IDs: ", length(unique(df$gene_id)), "\n")
  
  # One longest transcript per gene
  longest <- df %>%
    group_by(gene_id) %>%
    slice_max(order_by = len, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  cat("  Longest-by-gene: ", nrow(longest), " sequences\n")
  
  # ---- Write FASTA: longest isoforms (for this project, same as by-gene) ----
  write_fasta <- function(df_sub, path) {
    con <- file(path, "w")
    on.exit(close(con), add = TRUE)
    for (i in seq_len(nrow(df_sub))) {
      writeLines(paste0(">", df_sub$header[i]), con)
      writeLines(df_sub$seq[i], con)
    }
  }
  
  write_fasta(longest, iso_file)
  write_fasta(longest, gene_file)
  
  cat("  Saved longest isoforms FASTA to:  ", iso_file,  "\n")
  cat("  Saved longest-by-gene FASTA to:   ", gene_file, "\n\n")
}

cat("Done extracting longest isoforms / longest-by-gene.\n")
