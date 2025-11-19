############################################################
# 04_classify_genes_blast.R
#
# Uses BLASTN to classify longest isoforms as ace1 / ace2
# based on a reference database.
#
# Output → data/longest_by_gene/
############################################################

library(dplyr)
library(stringr)

iso_dir  <- file.path("data", "longest_isoforms")
out_dir  <- file.path("data", "longest_by_gene")
combined_dir <- file.path("data", "combined")

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

BLASTN_PATH <- "/opt/homebrew/bin/blastn"  # adjust if needed
BLAST_DB    <- "data/blast/ache_ref"       # no file extension

iso_files <- list.files(iso_dir, full.names = TRUE)

# Read combined metadata from script 03
meta_path <- file.path(combined_dir, "ache_longest_isoforms_metadata.csv")
meta <- read.csv(meta_path, stringsAsFactors = FALSE)

# ---- RUN BLAST AND CLASSIFY -------------------------------------------------

results <- list()

for (f in iso_files) {
  sp <- basename(f) |> str_replace("_longest_isoforms.fasta", "")
  out_blast <- tempfile(fileext = ".txt")
  
  cmd <- c(
    "-query", f,
    "-db", BLAST_DB,
    "-outfmt", "6 qseqid sseqid pident length bitscore",
    "-max_target_seqs", "5"
  )
  
  out <- system2(BLASTN_PATH, cmd, stdout = TRUE, stderr = TRUE)
  
  if (!length(out)) {
    cat("BLAST returned no output for", sp, "\n")
    next
  }
  
  hits <- read.table(text = out, sep = "\t", header = FALSE,
                     col.names = c("qseqid", "sseqid", "pident", "length", "bitscore"))
  
  if (nrow(hits) == 0) {
    cat("No BLAST hits for", sp, "\n")
    next
  }
  
  # Choose best hit
  best <- hits |> arrange(desc(bitscore)) |> slice(1)
  
  # Determine gene assignment from reference header
  assigned <- ifelse(str_detect(best$sseqid, "ace1"), "ace1", "ace2")
  
  results[[sp]] <- assigned
}

# ---- Write final BLAST-classified FASTA -------------------------------------

combined_fasta <- file.path(combined_dir, "ache_longest_isoforms_all_species.fasta")
df <- readLines(combined_fasta)

# Simple FASTA reader
idx <- grep("^>", df)
starts <- c(idx, length(df)+1)
parts <- data.frame(header=character(), seq=character(), stringsAsFactors=FALSE)

for (i in seq_along(idx)) {
  h <- df[idx[i]]
  s <- paste(df[(idx[i]+1):(starts[i+1]-1)], collapse="")
  parts[nrow(parts)+1,] <- c(h, s)
}

final <- data.frame()
for (sp in names(results)) {
  gene <- results[[sp]]
  keep <- parts[str_detect(parts$header, sp) & str_detect(parts$header, gene), ]
  final <- rbind(final, keep)
}

out_fasta <- file.path(out_dir, "ache_longest_by_gene_blast_classified.fasta")
con <- file(out_fasta, "w")
for (i in seq_len(nrow(final))) {
  writeLines(final$header[i], con)
  writeLines(final$seq[i], con)
}
close(con)

cat("BLAST-classified longest-by-gene FASTA saved to:\n  ", out_fasta, "\n")
cat("Done.\n")
