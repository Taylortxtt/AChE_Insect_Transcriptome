############################################################
# 02_extract_longest_isoforms.R
#
# Takes each species’ raw FASTA and extracts the *longest* 
# transcript per gene (ace1 / ace2).
#
# Output → data/longest_isoforms/
############################################################

library(stringr)

raw_dir <- file.path("data", "raw")
iso_dir <- file.path("data", "longest_isoforms")
if (!dir.exists(iso_dir)) dir.create(iso_dir, recursive = TRUE)

cat("Raw FASTAs directory:", normalizePath(raw_dir), "\n")
cat("Longest-isoform directory:", normalizePath(iso_dir), "\n\n")

species_list <- c(
  "gregaria", "piceifrons", "nitens", "cancellata",
  "serialis", "anabrus", "anopheles", "aedes", "culex"
)

# ---- Helper: read FASTA -----------------------------------------------------
read_fasta <- function(path) {
  lines <- readLines(path)
  headers <- grep("^>", lines)
  starts <- c(headers, length(lines) + 1)
  fasta_df <- data.frame(header = character(), seq = character(), stringsAsFactors = FALSE)
  
  for (i in seq_along(headers)) {
    h <- lines[headers[i]]
    s <- paste(lines[(headers[i] + 1):(starts[i + 1] - 1)], collapse = "")
    fasta_df[nrow(fasta_df) + 1, ] <- c(h, s)
  }
  return(fasta_df)
}

# ---- Main -----------------------------------------------------
for (sp in species_list) {
  cat("Species:", sp, "\n")
  
  in_file  <- file.path(raw_dir, paste0(sp, "_raw.fasta"))
  out_file <- file.path(iso_dir, paste0(sp, "_longest_isoforms.fasta"))
  
  if (!file.exists(in_file)) {
    cat("  Raw FASTA missing — skipping.\n\n")
    next
  }
  
  df <- read_fasta(in_file)
  
  # Mark ace1/ace2 by header text
  df$gene <- case_when(
    str_detect(df$header, regex("ace1", ignore_case = TRUE)) ~ "ace1",
    str_detect(df$header, regex("ace2", ignore_case = TRUE)) ~ "ace2",
    TRUE ~ NA
  )
  
  df <- df[!is.na(df$gene), ]
  if (nrow(df) == 0) {
    cat("  No ace1/ace2 transcripts found.\n\n")
    next
  }
  
  # Keep LONGEST transcript per gene
  df$len <- nchar(df$seq)
  longest <- df |>
    dplyr::group_by(gene) |>
    dplyr::slice_max(order_by = len, n = 1) |>
    dplyr::ungroup()
  
  # Write output
  con <- file(out_file, "w")
  for (i in seq_len(nrow(longest))) {
    writeLines(longest$header[i], con)
    writeLines(longest$seq[i], con)
  }
  close(con)
  
  cat("  Saved:", out_file, "\n\n")
}

cat("Done extracting longest isoforms.\n")
