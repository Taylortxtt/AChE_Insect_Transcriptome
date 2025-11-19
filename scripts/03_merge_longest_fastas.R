############################################################
# 03_merge_longest_fastas.R
#
# Merges all species’ longest isoform FASTAs into one dataset.
############################################################

library(stringr)
library(dplyr)

iso_dir <- file.path("data", "longest_isoforms")
combined_dir <- file.path("data", "combined")
if (!dir.exists(combined_dir)) dir.create(combined_dir, recursive = TRUE)

out_fasta <- file.path(combined_dir, "ache_longest_isoforms_all_species.fasta")
out_meta  <- file.path(combined_dir, "ache_longest_isoforms_metadata.csv")

iso_files <- list.files(iso_dir, full.names = TRUE)

# ---- Helper ---------------------------------------------------
read_fasta <- function(path) {
  lines <- readLines(path)
  headers <- grep("^>", lines)
  starts <- c(headers, length(lines) + 1)
  df <- data.frame(header = character(), seq = character(), stringsAsFactors = FALSE)
  for (i in seq_along(headers)) {
    df[nrow(df)+1, ] <- c(
      lines[headers[i]],
      paste(lines[(headers[i]+1):(starts[i+1]-1)], collapse = "")
    )
  }
  df
}

# ---- Merge -----------------------------------------------------
all_df <- data.frame()
writeLines("", out_fasta)

for (f in iso_files) {
  sp <- basename(f) |> str_replace("_longest_isoforms.fasta", "")
  df <- read_fasta(f)
  
  df$species <- sp
  df$gene <- ifelse(
    str_detect(df$header, "ace1", ignore_case = TRUE),
    "ace1", "ace2"
  )
  
  all_df <- rbind(all_df, df)
  
  # Append to combined FASTA
  for (i in seq_len(nrow(df))) {
    write(df$header[i], file = out_fasta, append = TRUE)
    write(df$seq[i], file = out_fasta, append = TRUE)
  }
}

write.csv(all_df, out_meta, row.names = FALSE)

cat("Combined FASTA:", out_fasta, "\n")
cat("Metadata CSV:", out_meta, "\n")
cat("Done merging.\n")
