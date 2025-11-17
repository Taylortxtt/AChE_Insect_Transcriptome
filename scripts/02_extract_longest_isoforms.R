############################################################
# 02_extract_longest_isoforms.R
#
# Goal
# ----
# For each species' raw FASTA in data/raw/:
#   1) Write all (filtered) raw sequences into a single
#      combined FASTA with standardized headers.
#   2) For each gene type ("AChE1", "AChE2", or "other"),
#      keep only the longest sequence per species and save
#      it to data/longest/.
#
# Header format:
#   >CommonName_GeneName_GeneID
#
# Example:
#   >gregaria_AChE_XM_049998832
############################################################

### 0. Helper: read a FASTA file into a data.frame ---------------------------

# Base-R parser: returns one row per sequence with
# columns: header, seq
read_fasta_to_df <- function(path) {
  lines      <- readLines(path)
  header_idx <- which(startsWith(lines, ">"))
  
  if (length(header_idx) == 0) {
    stop("No FASTA headers ('>') found in file: ", path)
  }
  
  # Add an artificial "end" index so each header knows where it stops
  header_idx_end <- c(header_idx[-1] - 1, length(lines))
  
  headers <- character(length(header_idx))
  seqs    <- character(length(header_idx))
  
  for (i in seq_along(header_idx)) {
    h_start <- header_idx[i]
    h_end   <- header_idx_end[i]
    
    headers[i] <- lines[h_start]
    seq_lines  <- lines[(h_start + 1):h_end]
    seqs[i]    <- paste(seq_lines, collapse = "")
  }
  
  data.frame(
    header = headers,
    seq    = seqs,
    stringsAsFactors = FALSE
  )
}

### 1. Classify gene type from header ---------------------------------------

# Very simple pattern-based classifier:
#   ace-1 / AChE1 / "acetylcholinesterase 1" → "AChE1"
#   ace-2 / AChE2 / "acetylcholinesterase 2" → "AChE2"
#   everything else                          → "other"
classify_gene_type <- function(header) {
  h <- tolower(header)
  
  if (grepl("ace-1|ache1|acetylcholinesterase 1", h)) {
    return("AChE1")
  }
  
  if (grepl("ace-2|ache2|acetylcholinesterase 2", h)) {
    return("AChE2")
  }
  
  # If we can't clearly tell 1 vs 2, call it "other"
  "other"
}

### 1B. Rename header to CommonName_GeneName_GeneID -------------------------

# CommonName = species_short (e.g. "gregaria")
# GeneName   = "AChE1", "AChE2", or generic "AChE"
# GeneID     = accession like XM_049998832 (no version suffix)
rename_header <- function(header, species_name, gene_type) {
  # Strip leading ">"
  id_raw <- sub("^>", "", header)
  
  # gene_id = first token before a space or dot
  #  e.g. XM_049998832 from "XM_049998832.1 something something"
  gene_id <- sub("^([^ .]+).*", "\\1", id_raw)
  
  # Collapse gene_type for labeling
  gene_label <- ifelse(gene_type %in% c("AChE1", "AChE2"),
                       gene_type,
                       "AChE")
  
  paste0(">", species_name, "_", gene_label, "_", gene_id)
}

### 2. Directories -----------------------------------------------------------

raw_dir      <- file.path("data", "raw")
longest_dir  <- file.path("data", "longest")
combined_dir <- file.path("data", "combined")

if (!dir.exists(longest_dir)) {
  dir.create(longest_dir, recursive = TRUE)
}
if (!dir.exists(combined_dir)) {
  dir.create(combined_dir, recursive = TRUE)
}

cat("Raw directory:      ", normalizePath(raw_dir),      "\n")
cat("Longest directory:  ", normalizePath(longest_dir),  "\n")
cat("Combined directory: ", normalizePath(combined_dir), "\n")

# Combined RAW FASTA (all species together, renamed headers)
combined_raw_path <- file.path(combined_dir, "AChE_raw_all_species.fasta")
combined_con      <- file(combined_raw_path, open = "w")

### 3. Species short names (match script 01) ---------------------------------

species_short <- c(
  "gregaria",
  "cancellata",
  "piceifrons",
  "anabrus",
  "anopheles",
  "aedes",
  "culex"
)

### 4. Process each species --------------------------------------------------

for (sp in species_short) {
  cat("-------------------------------------------------------------\n")
  cat("Processing species:", sp, "\n")
  
  raw_path <- file.path(raw_dir, paste0(sp, "_raw.fasta"))
  
  if (!file.exists(raw_path)) {
    warning("RAW FASTA not found for ", sp, " at ", raw_path)
    next
  }
  
  fasta_df <- read_fasta_to_df(raw_path)
  cat("  Sequences read:", nrow(fasta_df), "\n")
  
  # Special case: drop the partial Culex record JX292118.1
  if (sp == "culex") {
    before   <- nrow(fasta_df)
    fasta_df <- fasta_df[!grepl("JX292118.1", fasta_df$header), ]
    after    <- nrow(fasta_df)
    
    cat("  Removed partial Culex sequence JX292118.1 (",
        before - after, " record(s) dropped)\n", sep = "")
  }
  
  # If nothing left after filtering, skip this species
  if (nrow(fasta_df) == 0) {
    warning("No sequences left for species ", sp, " after filtering.")
    next
  }
  
  # Annotate sequences
  fasta_df$length    <- nchar(fasta_df$seq)
  fasta_df$gene_type <- vapply(fasta_df$header, classify_gene_type, character(1))
  
  ### 4A. Add all (renamed) raw sequences to the combined FASTA -------------
  
  raw_headers_renamed <- vapply(
    seq_len(nrow(fasta_df)),
    function(i) {
      rename_header(
        header       = fasta_df$header[i],
        species_name = sp,
        gene_type    = fasta_df$gene_type[i]
      )
    },
    character(1)
  )
  
  for (i in seq_len(nrow(fasta_df))) {
    writeLines(raw_headers_renamed[i], combined_con)
    writeLines(fasta_df$seq[i],        combined_con)
  }
  
  ### 4B. For each gene_type, keep the longest sequence ----------------------
  
  # Split by gene_type, then take the row with the max length in each group.
  longest_list <- lapply(
    split(fasta_df, fasta_df$gene_type),
    function(df) df[which.max(df$length), , drop = FALSE]
  )
  
  longest_df <- do.call(rbind, longest_list)
  
  cat("  Gene types found:",
      paste(unique(longest_df$gene_type), collapse = ", "), "\n")
  cat("  Sequences kept:", nrow(longest_df), "\n")
  
  ### 4C. Write per-species 'longest' FASTA with renamed headers -------------
  
  out_path <- file.path(longest_dir, paste0(sp, "_longest.fasta"))
  con      <- file(out_path, open = "w")
  
  for (i in seq_len(nrow(longest_df))) {
    new_header <- rename_header(
      header       = longest_df$header[i],
      species_name = sp,
      gene_type    = longest_df$gene_type[i]
    )
    
    writeLines(new_header,        con)
    writeLines(longest_df$seq[i], con)
  }
  
  close(con)
  cat("  Saved longest isoforms to:", out_path, "\n")
}

### 5. Close combined RAW FASTA and wrap up ---------------------------------

close(combined_con)
cat("Combined RAW FASTA saved to:", combined_raw_path, "\n")
cat("-------------------------------------------------------------\n")
cat("Done! Check data/longest/ for per-species longest-isoform FASTA files.\n")
