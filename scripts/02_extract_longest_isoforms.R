############################################################
# 02_extract_longest_isoforms.R
#
# Goal
# ----
# For each species, take all raw AChE(-like) sequences in
# data/raw/ and keep ONLY the longest transcript per *gene*.
#
# Then, additionally, filter to sequences whose headers
# or gene_ids clearly look like ace-1 / ace-2 and collapse
# to at most one longest ace-1 and one longest ace-2 per
# species (when such labels exist).
#
# Outputs
# -------
# One FASTA per species in data/longest/:
#   data/longest/<short_name>_longest_by_gene.fasta
############################################################

### 0. Packages --------------------------------------------------------------

# base R only

### 1. Project folders -------------------------------------------------------

raw_dir     <- file.path("data", "raw")
longest_dir <- file.path("data", "longest")

if (!dir.exists(raw_dir)) {
  stop("Raw directory 'data/raw' does not exist. ",
       "Run 01_download_AChE_sequences.R first.")
}

if (!dir.exists(longest_dir)) {
  dir.create(longest_dir, recursive = TRUE)
}

cat("Raw FASTAs directory:     ", normalizePath(raw_dir), "\n")
cat("Longest-isoform directory:", normalizePath(longest_dir), "\n\n")

### 2. Species list (must match script 01) -----------------------------------

species_list <- list(
  list(short_name = "gregaria",   organism_term = "Schistocerca gregaria"),
  list(short_name = "cancellata", organism_term = "Schistocerca cancellata"),
  list(short_name = "piceifrons", organism_term = "Schistocerca piceifrons"),
  list(short_name = "anabrus",    organism_term = "Anabrus simplex"),
  list(short_name = "anopheles",  organism_term = "Anopheles gambiae"),
  list(short_name = "aedes",      organism_term = "Aedes aegypti"),
  list(short_name = "culex",      organism_term = "Culex quinquefasciatus")
)

### 3. Helper: parse FASTA file ----------------------------------------------

read_fasta <- function(path) {
  lines <- readLines(path)
  is_header <- startsWith(lines, ">")
  if (!any(is_header)) {
    stop("No FASTA headers (>) detected in file: ", path)
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

### 4. Helper: infer gene ID from a FASTA header -----------------------------

get_gene_id_from_header <- function(header) {
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
  first_token    <- strsplit(header, " ")[[1]][1]
  accession_root <- sub("\\.[0-9]+$", "", first_token)
  return(accession_root)
}

### 5. Helper: choose longest transcript per gene ----------------------------

select_longest_per_gene <- function(fasta_df) {
  fasta_df$gene_id <- vapply(
    fasta_df$header,
    get_gene_id_from_header,
    FUN.VALUE = character(1)
  )
  fasta_df$seq_len <- nchar(fasta_df$seq)
  
  o <- order(fasta_df$gene_id, -fasta_df$seq_len)
  fasta_df <- fasta_df[o, ]
  keep_idx <- !duplicated(fasta_df$gene_id)
  
  fasta_df[keep_idx, c("header", "seq", "gene_id", "seq_len")]
}

### 6. Helper: write FASTA back to disk --------------------------------------

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

### 7. Patterns: ace-1 / ace-2 detection ------------------------------------

# Broad match for any ace-1 / ace-2 mention
ace12_pattern <- "(ace[-_ ]?1|ace1|achE1|AChE1|acetylcholinesterase[- ]?1|
                   ace[-_ ]?2|ace2|achE2|AChE2|acetylcholinesterase[- ]?2)"

is_ace12_like <- function(header, gene_id) {
  h_match <- grepl(ace12_pattern, header,  ignore.case = TRUE)
  g_match <- grepl(ace12_pattern, gene_id, ignore.case = TRUE)
  h_match | g_match
}

# Finer classification into ace1 vs ace2
ace1_pattern_type <- "(ace[-_ ]?1|ace1|achE1|AChE1|acetylcholinesterase[- ]?1)"
ace2_pattern_type <- "(ace[-_ ]?2|ace2|achE2|AChE2|acetylcholinesterase[- ]?2)"

classify_ace_type <- function(header, gene_id) {
  is_ace1 <- grepl(ace1_pattern_type, header,  ignore.case = TRUE) |
    grepl(ace1_pattern_type, gene_id, ignore.case = TRUE)
  is_ace2 <- grepl(ace2_pattern_type, header,  ignore.case = TRUE) |
    grepl(ace2_pattern_type, gene_id, ignore.case = TRUE)
  
  if (is_ace1) return("ace1")
  if (is_ace2) return("ace2")
  return("unknown")
}

### 8. Main loop over species -------------------------------------------------

cat("Selecting longest transcript per gene (ace-1 / ace-2 only)...\n\n")

for (sp in species_list) {
  short_name <- sp$short_name
  
  in_file  <- file.path(raw_dir,     paste0(short_name, "_raw.fasta"))
  out_file <- file.path(longest_dir, paste0(short_name, "_longest_by_gene.fasta"))
  
  cat("--------------------------------------------------\n")
  cat("Species:", short_name, "\n")
  cat("Input FASTA: ", in_file,  "\n")
  cat("Output FASTA:", out_file, "\n")
  
  if (!file.exists(in_file)) {
    warning("Raw FASTA not found for ", short_name,
            ". Run 01_download_AChE_sequences.R first.")
    next
  }
  
  fasta_df <- read_fasta(in_file)
  cat("  Raw sequences:", nrow(fasta_df), "\n")
  
  # Longest per gene
  longest_df <- select_longest_per_gene(fasta_df)
  cat("  Unique genes before filters:", nrow(longest_df), "\n")
  
  # Length sanity filter (remove huge scaffolds)
  max_len    <- 6000
  before_len <- nrow(longest_df)
  longest_df <- longest_df[longest_df$seq_len <= max_len, ]
  after_len  <- nrow(longest_df)
  cat("  Filtered out", before_len - after_len,
      "overly large sequences (> ", max_len, " bp).\n", sep = "")
  
  # ace-1 / ace-2-like filter
  ace_flag <- is_ace12_like(longest_df$header, longest_df$gene_id)
  n_ace    <- sum(ace_flag)
  
  if (n_ace > 0) {
    cat("  Initially keeping", n_ace,
        "genes whose headers look like ace-1 / ace-2.\n")
    longest_df <- longest_df[ace_flag, ]
    
    # Classify into ace1 / ace2 / unknown
    ace_type <- vapply(
      seq_len(nrow(longest_df)),
      function(i) classify_ace_type(longest_df$header[i], longest_df$gene_id[i]),
      FUN.VALUE = character(1)
    )
    longest_df$ace_type <- ace_type
    
    # Collapse to at most one ace1 and one ace2 per species
    has_typed <- ace_type %in% c("ace1", "ace2")
    if (any(has_typed)) {
      keep_idx <- rep(FALSE, nrow(longest_df))
      for (t in c("ace1", "ace2")) {
        idx <- which(longest_df$ace_type == t)
        if (length(idx) > 0) {
          # keep the longest sequence within this type
          best_idx <- idx[which.max(longest_df$seq_len[idx])]
          keep_idx[best_idx] <- TRUE
        }
      }
      longest_df <- longest_df[keep_idx, ]
      cat("  After collapsing to one per ace-1 / ace-2 type:",
          nrow(longest_df), "genes.\n")
    } else {
      cat("  NOTE: ace-1/ace-2-like pattern matched, but unable to",
          "classify into ace1/ace2 types; keeping all ace-like genes.\n")
    }
    
  } else {
    cat("  WARNING: No clear ace-1 / ace-2-like headers found for this species.\n")
    cat("           Keeping ALL genes after length filter instead.\n")
  }
  
  cat("  Genes kept for downstream analysis:", nrow(longest_df), "\n")
  
  write_fasta(longest_df, out_file)
  cat("  Saved longest-per-gene FASTA.\n\n")
}

cat("--------------------------------------------------\n")
cat("Done! Check data/longest/ for *_longest_by_gene.fasta files.\n")
