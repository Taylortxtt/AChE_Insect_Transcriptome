############################################################
# 09_classify_genes_blast.R
#
# Use BLAST to classify each longest-per-gene sequence as
# ace-1, ace-2, or other based on a curated reference DB,
# and create a BLAST-confirmed ace1/ace2 FASTA.
############################################################

### 0. Load config -----------------------------------------------------------

if (!file.exists("config.R")) {
  stop("config.R not found in project root.")
}
source("config.R")

if (!exists("BLASTN_PATH") || !nzchar(BLASTN_PATH)) {
  stop("BLASTN_PATH not defined in config.R")
}
if (!exists("BLAST_DB") || !nzchar(BLAST_DB)) {
  stop("BLAST_DB not defined in config.R")
}

### 1. Input files -----------------------------------------------------------

combined_fasta <- file.path("data", "combined",
                            "ache_longest_by_gene_all_species.fasta")
meta_file      <- file.path("data", "combined",
                            "ache_longest_by_gene_metadata.csv")

if (!file.exists(combined_fasta)) {
  stop("Combined FASTA not found at: ", combined_fasta,
       "\nRun 03_merge_longest_fastas.R first.")
}
if (!file.exists(meta_file)) {
  stop("Metadata file not found at: ", meta_file,
       "\nRun 03_merge_longest_fastas.R first.")
}

query_path <- normalizePath(combined_fasta)

cat("Combined FASTA: ", query_path, "\n")
cat("Metadata CSV:   ", normalizePath(meta_file), "\n")
cat("BLASTN:         ", BLASTN_PATH, "\n")
cat("BLAST DB:       ", BLAST_DB, "\n\n")

### 2. Load metadata + FASTA -------------------------------------------------

meta <- read.csv(meta_file, stringsAsFactors = FALSE)
meta$qseqid <- sub(" .*", "", meta$header)

read_fasta <- function(path) {
  lines <- readLines(path)
  is_header <- startsWith(lines, ">")
  if (!any(is_header)) {
    stop("No FASTA headers found in ", path)
  }
  header_idx     <- which(is_header)
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

fasta_df <- read_fasta(query_path)
fasta_df$qseqid <- sub(" .*", "", fasta_df$header)

### 3. Run BLAST with outfmt 6 ----------------------------------------------

blast_args <- c(
  "-query", query_path,
  "-db", BLAST_DB,
  "-outfmt", "6",          # standard tabular, 12 columns
  "-max_target_seqs", "10"
)

cat("Running BLASTN with:\n")
cat("  ", BLASTN_PATH, paste(blast_args, collapse = " "), "\n\n")

blast_res <- system2(
  command = BLASTN_PATH,
  args    = blast_args,
  stdout  = TRUE,
  stderr  = TRUE
)

# If BLAST produced absolutely nothing, that usually just
# means: no hits, no warnings. Treat that as "no_hit".
if (length(blast_res) == 0) {
  cat("BLASTN produced no output (no hits, no warnings).\n")
  cat("Interpreting this as: no sequences had BLAST hits to the\n",
      "current reference_ache database.\n\n")
  
  meta_blast <- meta
  meta_blast$ace_type_blast <- "no_hit"
  
  out_meta_blast <- file.path(
    "data", "combined",
    "ache_longest_by_gene_metadata_blast.csv"
  )
  write.csv(meta_blast, out_meta_blast, row.names = FALSE)
  cat("Metadata with BLAST classification (all 'no_hit') saved to:\n  ",
      normalizePath(out_meta_blast), "\n\n")
  
  summary_tab <- as.data.frame(
    table(meta_blast$species_short, meta_blast$ace_type_blast),
    stringsAsFactors = FALSE
  )
  colnames(summary_tab) <- c("species_short", "ace_type_blast", "count")
  
  out_summary <- file.path(
    "data", "combined",
    "ache_blast_summary_by_species.csv"
  )
  write.csv(summary_tab, out_summary, row.names = FALSE)
  cat("BLAST summary by species saved to:\n  ",
      normalizePath(out_summary), "\n\n")
  
  cat("NOTE: No BLAST hits were found with the current reference_ache.fasta.\n",
      "You may want to supply real ace1/ace2 reference sequences and rerun.\n\n")
  
  cat("Done (no BLAST hits).\n")
  quit(save = "no")  # stop this script cleanly if sourced via Rscript
}

### 4. Parse BLAST table (normal case with hits) -----------------------------

is_table_line <- grepl("^[^#].*\\t", blast_res)
table_lines   <- blast_res[is_table_line]
msg_lines     <- blast_res[!is_table_line]

if (length(msg_lines) > 0) {
  cat("BLAST messages:\n")
  cat(paste0("  ", msg_lines, collapse = "\n"), "\n\n")
}

if (length(table_lines) == 0) {
  cat("BLASTN produced only messages but no tabular hits.\n",
      "Treating all sequences as 'no_hit'.\n\n")
  
  meta_blast <- meta
  meta_blast$ace_type_blast <- "no_hit"
  
  out_meta_blast <- file.path(
    "data", "combined",
    "ache_longest_by_gene_metadata_blast.csv"
  )
  write.csv(meta_blast, out_meta_blast, row.names = FALSE)
  cat("Metadata with BLAST classification (all 'no_hit') saved to:\n  ",
      normalizePath(out_meta_blast), "\n\n")
  
  summary_tab <- as.data.frame(
    table(meta_blast$species_short, meta_blast$ace_type_blast),
    stringsAsFactors = FALSE
  )
  colnames(summary_tab) <- c("species_short", "ace_type_blast", "count")
  
  out_summary <- file.path(
    "data", "combined",
    "ache_blast_summary_by_species.csv"
  )
  write.csv(summary_tab, out_summary, row.names = FALSE)
  cat("BLAST summary by species saved to:\n  ",
      normalizePath(out_summary), "\n\n")
  
  cat("Done (no BLAST hits).\n")
  quit(save = "no")
}

blast_df <- read.table(
  text            = paste(table_lines, collapse = "\n"),
  sep             = "\t",
  header          = FALSE,
  stringsAsFactors = FALSE
)

colnames(blast_df) <- c(
  "qseqid", "sseqid", "pident", "aln_length",
  "mismatch", "gapopen",
  "qstart", "qend", "sstart", "send",
  "evalue", "bitscore"
)

cat("Total BLAST hits:", nrow(blast_df), "\n")

blast_df <- blast_df[order(blast_df$qseqid,
                           -blast_df$bitscore,
                           blast_df$evalue), ]
best_idx  <- !duplicated(blast_df$qseqid)
best_hits <- blast_df[best_idx, ]

cat("Unique queries with ≥1 hit:", nrow(best_hits), "\n\n")

### 5. Classify ace-1 / ace-2 / other ---------------------------------------

classify_from_subject <- function(sid) {
  s <- tolower(sid)
  if (grepl("ace1", s) || grepl("ache1", s)) return("ace1")
  if (grepl("ace2", s) || grepl("ache2", s)) return("ace2")
  return("other")
}

best_hits$ace_type_blast <- vapply(
  best_hits$sseqid,
  classify_from_subject,
  FUN.VALUE = character(1)
)

### 6. Merge back into metadata ---------------------------------------------

meta_blast <- merge(
  meta,
  best_hits,
  by = "qseqid",
  all.x = TRUE
)

meta_blast$ace_type_blast[is.na(meta_blast$ace_type_blast)] <- "no_hit"

out_meta_blast <- file.path(
  "data", "combined",
  "ache_longest_by_gene_metadata_blast.csv"
)
write.csv(meta_blast, out_meta_blast, row.names = FALSE)
cat("Metadata with BLAST classification saved to:\n  ",
    normalizePath(out_meta_blast), "\n\n")

### 7. Summary by species ----------------------------------------------------

summary_tab <- as.data.frame(
  table(meta_blast$species_short, meta_blast$ace_type_blast),
  stringsAsFactors = FALSE
)
colnames(summary_tab) <- c("species_short", "ace_type_blast", "count")

out_summary <- file.path(
  "data", "combined",
  "ache_blast_summary_by_species.csv"
)
write.csv(summary_tab, out_summary, row.names = FALSE)
cat("BLAST summary by species saved to:\n  ",
    normalizePath(out_summary), "\n\n")

### 8. BLAST-confirmed ace1/ace2 FASTA --------------------------------------

keep_qseqid <- meta_blast$qseqid[meta_blast$ace_type_blast %in% c("ace1", "ace2")]
keep_fasta  <- fasta_df[fasta_df$qseqid %in% keep_qseqid, ]

cat("Sequences classified as ace1/ace2:",
    nrow(keep_fasta), "of", nrow(fasta_df), "\n")

if (nrow(keep_fasta) > 0) {
  out_fasta <- file.path(
    "data", "combined",
    "ache_longest_by_gene_ace12_blast_only.fasta"
  )
  lines <- character(2 * nrow(keep_fasta))
  j <- 1
  for (i in seq_len(nrow(keep_fasta))) {
    lines[j]   <- paste0(">", keep_fasta$header[i])
    lines[j+1] <- keep_fasta$seq[i]
    j <- j + 2
  }
  writeLines(lines, con = out_fasta)
  cat("BLAST-confirmed ace1/ace2 FASTA written to:\n  ",
      normalizePath(out_fasta), "\n")
} else {
  cat("NOTE: No sequences classified as ace1/ace2 by BLAST.\n")
}

cat("\nDone.\n")
