############################################################
# 04_classify_genes_blast.R
#
# Use BLASTN against a small AChE reference panel to classify
# each "longest-by-gene" transcript as ace1 / ace2 / other.
#
# Inputs
#   data/combined/ache_longest_by_gene_all_species.fasta
#   data/combined/ache_longest_by_gene_metadata.csv
#   BLAST DB: data/blast/ache_ref  (built from reference_ache.fasta)
#
# Outputs
#   data/combined/ache_longest_by_gene_metadata_blast.csv
#   data/combined/ache_longest_by_gene_ace12_blast_only.fasta
############################################################

source("config.R")   # provides BLASTN_PATH and BLAST_DB

library(dplyr)
library(stringr)

# ---- Paths ---------------------------------------------------------

combo_fasta <- file.path("data", "combined", "ache_longest_by_gene_all_species.fasta")
meta_csv    <- file.path("data", "combined", "ache_longest_by_gene_metadata.csv")
out_meta    <- file.path("data", "combined", "ache_longest_by_gene_metadata_blast.csv")
out_fasta   <- file.path("data", "combined", "ache_longest_by_gene_ace12_blast_only.fasta")

if (!file.exists(combo_fasta)) {
  stop("Combined FASTA not found: ", combo_fasta)
}
if (!file.exists(meta_csv)) {
  stop("Metadata CSV not found: ", meta_csv)
}
if (!file.exists(BLASTN_PATH)) {
  stop("BLASTN_PATH does not exist: ", BLASTN_PATH)
}

if (!dir.exists("data/blast")) dir.create("data/blast", recursive = TRUE)

# Check that the BLAST DB is there (nhr file is enough)
if (!file.exists(paste0(BLAST_DB, ".nhr"))) {
  stop(
    "BLAST DB files not found for '", BLAST_DB,
    "'.\nBuild them with makeblastdb using data/blast/reference_ache.fasta."
  )
}

cat("Combined FASTA: ", normalizePath(combo_fasta), "\n")
cat("Metadata CSV:   ", normalizePath(meta_csv), "\n")
cat("BLASTN:         ", BLASTN_PATH, "\n")
cat("BLAST DB:       ", BLAST_DB, "\n\n")

# ---- Run BLASTN ----------------------------------------------------

blast_args <- c(
  "-query", combo_fasta,
  "-db", BLAST_DB,
  "-outfmt", "6 qseqid sseqid pident length evalue bitscore",
  "-max_target_seqs", "5"
)

cat("Running BLASTN with:\n  ",
    BLASTN_PATH, paste(blast_args, collapse = " "), "\n\n")

blast_out <- system2(
  BLASTN_PATH,
  args   = blast_args,
  stdout = TRUE,
  stderr = TRUE
)

if (length(blast_out) == 0) {
  stop("BLASTN returned no output at all. Check BLAST installation and database.")
}

# Keep only tab-delimited result lines (stdout may contain warnings)
tab_lines <- blast_out[grepl("\t", blast_out)]
if (length(tab_lines) == 0) {
  stop(
    "No tab-delimited BLAST results found.\n",
    "Full BLAST output:\n", paste(blast_out, collapse = "\n")
  )
}

blast_df <- read.table(
  text = paste(tab_lines, collapse = "\n"),
  sep  = "\t",
  header = FALSE,
  stringsAsFactors = FALSE
)
colnames(blast_df) <- c("qseqid", "sseqid", "pident", "align_len", "evalue", "bitscore")

# ---- Pick best hit per query & classify ----------------------------

best_hits <- blast_df %>%
  group_by(qseqid) %>%
  slice_max(order_by = bitscore, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    ace_type_blast = case_when(
      str_detect(sseqid, "ace1", ignore_case = TRUE) ~ "ace1",
      str_detect(sseqid, "ace2", ignore_case = TRUE) ~ "ace2",
      TRUE ~ "other"
    )
  )

meta <- read.csv(meta_csv, stringsAsFactors = FALSE)

meta2 <- meta %>%
  left_join(
    best_hits %>%
      select(qseqid, sseqid, ace_type_blast, pident, align_len, evalue, bitscore),
    by = "qseqid"
  )

write.csv(meta2, out_meta, row.names = FALSE)
cat("Metadata with BLAST classification written to:\n  ",
    normalizePath(out_meta), "\n\n")

# ---- Write FASTA with BLAST-confirmed ace1/ace2 --------------------

keep_ids <- meta2 %>%
  filter(ace_type_blast %in% c("ace1", "ace2")) %>%
  pull(qseqid)

cat("Sequences classified as ace1/ace2 by BLAST: ",
    length(keep_ids), " of ", nrow(meta2), "\n", sep = "")

if (length(keep_ids) > 0) {
  lines <- readLines(combo_fasta)
  header_idx <- grep("^>", lines)
  header_idx <- c(header_idx, length(lines) + 1)
  
  keep_lines <- character(0)
  
  for (i in seq_len(length(header_idx) - 1)) {
    h <- lines[header_idx[i]]
    id <- sub("^>", "", h)
    id <- strsplit(id, "\\s+")[[1]][1]
    
    if (id %in% keep_ids) {
      seq_chunk <- lines[(header_idx[i] + 1):(header_idx[i + 1] - 1)]
      keep_lines <- c(keep_lines, h, seq_chunk)
    }
  }
  
  writeLines(keep_lines, out_fasta)
  cat("BLAST-confirmed ace1/ace2 FASTA written to:\n  ",
      normalizePath(out_fasta), "\n")
} else {
  cat("NOTE: No sequences were classified as ace1/ace2 by BLAST.\n")
}

cat("\nDone.\n")
