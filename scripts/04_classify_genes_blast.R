############################################################
# 04_classify_genes_blast.R
#
# Uses BLASTN against a reference AChE database to classify
# each "longest-by-gene" transcript as ace1 / ace2 / other.
#
# Inputs  (from script 03):
#   data/combined/ache_longest_by_gene_all_species.fasta
#   data/combined/ache_longest_by_gene_metadata.csv
#
# Outputs:
#   data/combined/ache_longest_by_gene_metadata_blast.csv
#   data/blast/ache_longest_by_gene_blast6.tsv
############################################################

# Load BLAST paths (BLASTN_PATH, BLAST_DB)
source("config.R")

library(readr)
library(dplyr)
library(stringr)

## 1. Paths -------------------------------------------------------------

query_fasta <- file.path("data", "combined",
                         "ache_longest_by_gene_all_species.fasta")
meta_csv    <- file.path("data", "combined",
                         "ache_longest_by_gene_metadata.csv")
blast_dir   <- file.path("data", "blast")
if (!dir.exists(blast_dir)) dir.create(blast_dir, recursive = TRUE)

blast_out_tsv <- file.path(blast_dir, "ache_longest_by_gene_blast6.tsv")

if (!file.exists(query_fasta)) {
  stop("Combined FASTA not found: ", normalizePath(query_fasta))
}
if (!file.exists(meta_csv)) {
  stop("Metadata CSV not found: ", normalizePath(meta_csv))
}

cat("Combined FASTA: ", query_fasta, "\n")
cat("Metadata CSV:   ", meta_csv, "\n")
cat("BLASTN:         ", BLASTN_PATH, "\n")
cat("BLAST DB:       ", BLAST_DB, "\n\n")

## 2. Run BLASTN --------------------------------------------------------

# NOTE: outfmt fields are passed as ONE string so BLAST parses them correctly.
blast_args <- c(
  "-query", query_fasta,
  "-db",    BLAST_DB,
  "-outfmt", "6 qseqid sseqid pident length evalue bitscore",
  "-max_target_seqs", "5",
  "-out", blast_out_tsv
)

cat("Running BLASTN with command:\n ",
    BLASTN_PATH, paste(blast_args, collapse = " "), "\n\n")

status <- system2(
  BLASTN_PATH,
  args = blast_args
)

if (status != 0) {
  stop("BLASTN exited with non-zero status code: ", status,
       "\nCheck that BLASTN_PATH and BLAST_DB in config.R are correct.")
}

if (!file.exists(blast_out_tsv) || file.size(blast_out_tsv) == 0) {
  stop("BLAST output file is missing or empty: ",
       normalizePath(blast_out_tsv))
}

cat("Raw BLAST table written to:\n  ",
    normalizePath(blast_out_tsv), "\n\n")

## 3. Read BLAST results ------------------------------------------------

blast_df <- read_tsv(
  blast_out_tsv,
  col_names = c("qseqid", "sseqid", "pident", "length", "evalue", "bitscore"),
  show_col_types = FALSE
)

if (nrow(blast_df) == 0) {
  stop("BLAST table exists but has 0 rows: ",
       normalizePath(blast_out_tsv))
}

# Best hit per query (highest bitscore, then lowest evalue)
blast_best <- blast_df %>%
  arrange(qseqid, desc(bitscore), evalue) %>%
  group_by(qseqid) %>%
  slice(1) %>%
  ungroup()

# Classify ace1 / ace2 / other based on reference ID
blast_best <- blast_best %>%
  mutate(
    ace_type_blast = case_when(
      str_detect(sseqid, regex("ace1", ignore_case = TRUE)) ~ "ace1",
      str_detect(sseqid, regex("ace2", ignore_case = TRUE)) ~ "ace2",
      TRUE ~ "other"
    )
  )

## 4. Join with metadata ------------------------------------------------

meta_df <- read_csv(meta_csv, show_col_types = FALSE)

# meta_df$qseqid comes from script 03; should match BLAST qseqid
meta_joined <- meta_df %>%
  left_join(
    blast_best %>%
      select(qseqid, sseqid, pident, length, evalue, bitscore, ace_type_blast),
    by = "qseqid"
  )

out_meta_blast <- file.path("data", "combined",
                            "ache_longest_by_gene_metadata_blast.csv")
write_csv(meta_joined, out_meta_blast)

cat("BLAST-annotated metadata written to:\n  ",
    normalizePath(out_meta_blast), "\n\n")

## 5. Simple summary ----------------------------------------------------

summary_tab <- meta_joined %>%
  group_by(species, ace_type_blast) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(species, ace_type_blast)

cat("Sequences classified by BLAST (per species):\n")
print(summary_tab)
cat("\nDone.\n")
