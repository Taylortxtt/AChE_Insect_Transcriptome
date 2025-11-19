############################################################
# 08_screen_mutations.R
#
# Goal
# ----
# Screen aligned AChE sequences for a set of known or
# hypothesized resistance-associated mutations.
#
# Input
# -----
#   data/alignment/ache_longest_by_gene_aligned.fasta
#
# Outputs
# -------
#   data/combined/ache_mutation_presence_matrix.csv
#
# Notes
# -----
# * Positions are interpreted in terms of the MAFFT alignment
#   (1-based alignment positions, including gaps).
# * mutation_definitions can be edited to reflect curated
#   literature mutations or new hypotheses.
############################################################

### 0. User-editable mutation definitions -----------------------------------
# Each row defines one mutation to screen for, using:
#   mutation_id        - short label for the mutation
#   alignment_position - integer (1-based) in the alignment
#   ref_aa             - expected reference amino acid
#   alt_aa             - "mutant" amino acid we are looking for
#
# You can edit, remove, or add rows here as needed.

mutation_definitions <- data.frame(
  mutation_id        = c("Example_1", "Example_2"),
  alignment_position = c(100, 250),
  ref_aa             = c("A", "G"),
  alt_aa             = c("V", "D"),
  stringsAsFactors   = FALSE
)

### 1. Paths -----------------------------------------------------------------

align_dir    <- file.path("data", "alignment")
combined_dir <- file.path("data", "combined")

if (!dir.exists(align_dir)) {
  stop("Directory 'data/alignment' does not exist. ",
       "Run 04_run_alignment.R first.")
}
if (!dir.exists(combined_dir)) {
  dir.create(combined_dir, recursive = TRUE)
}

align_fasta <- file.path(
  align_dir,
  "ache_longest_by_gene_aligned.fasta"
)

out_matrix_csv <- file.path(
  combined_dir,
  "ache_mutation_presence_matrix.csv"
)

cat("Alignment FASTA: ", normalizePath(align_fasta),   "\n")
cat("Output CSV:      ", normalizePath(out_matrix_csv), "\n\n")

if (!file.exists(align_fasta)) {
  stop("Alignment FASTA not found at: ", align_fasta,
       "\nRun 04_run_alignment.R first.")
}

### 2. FASTA reader + gene ID helper (reuse logic from earlier scripts) ------

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

get_gene_id_from_header <- function(header) {
  # Try gene=, then LOC..., then accession root as fallback.
  gene_match <- regexpr("gene=([^ ]+)", header)
  if (gene_match != -1) {
    gene <- regmatches(header, gene_match)
    gene <- sub("^gene=", "", gene)
    return(gene)
  }
  
  loc_match <- regexpr("LOC[0-9]+", header)
  if (loc_match != -1) {
    loc_id <- regmatches(header, loc_match)
    return(loc_id)
  }
  
  first_token     <- strsplit(header, " ")[[1]][1]
  accession_root  <- sub("\\.[0-9]+$", "", first_token)
  return(accession_root)
}

get_short_name <- function(tip_label) {
  # tip_label looks like: shortname|original_header
  strsplit(tip_label, "\\|")[[1]][1]
}

### 3. Read alignment --------------------------------------------------------

aln_df <- read_fasta(align_fasta)

cat("Number of aligned sequences:", nrow(aln_df), "\n")
cat("Alignment length (characters) of first sequence:",
    nchar(aln_df$seq[1]), "\n\n")

# Basic safety: check all sequences have the same length
aln_lengths <- nchar(aln_df$seq)
if (length(unique(aln_lengths)) != 1) {
  warning("Not all sequences have the same alignment length. ",
          "Some positions may be out of bounds for shorter sequences.")
}

### 4. Build base metadata for matrix ----------------------------------------

aln_df$species_short <- vapply(aln_df$header, get_short_name, character(1))
aln_df$gene_id       <- vapply(aln_df$header, get_gene_id_from_header,
                               FUN.VALUE = character(1))

# Initialize the result data frame with base columns
result_df <- aln_df[, c("species_short", "gene_id", "header")]
rownames(result_df) <- NULL

### 5. For each mutation, compute presence/absence ---------------------------

# For each mutation row, we create a new column in result_df
for (i in seq_len(nrow(mutation_definitions))) {
  mut_id   <- mutation_definitions$mutation_id[i]
  pos      <- mutation_definitions$alignment_position[i]
  ref_aa   <- mutation_definitions$ref_aa[i]
  alt_aa   <- mutation_definitions$alt_aa[i]
  
  cat("Screening mutation:", mut_id,
      "at alignment position", pos,
      "(ref:", ref_aa, "alt:", alt_aa, ")\n")
  
  presence <- rep(NA, nrow(aln_df))  # default NA
  
  for (j in seq_len(nrow(aln_df))) {
    seq_j <- aln_df$seq[j]
    if (nchar(seq_j) < pos) {
      presence[j] <- NA
      next
    }
    
    aa <- substr(seq_j, pos, pos)
    
    if (aa == "-") {
      # gap at this position
      presence[j] <- NA
    } else if (aa == alt_aa) {
      presence[j] <- TRUE
    } else if (aa == ref_aa) {
      presence[j] <- FALSE
    } else {
      # neither ref nor alt (different amino acid)
      presence[j] <- FALSE
    }
  }
  
  # Add this mutation as a new logical column
  col_name <- paste0("mut_", mut_id)
  result_df[[col_name]] <- presence
}

cat("\nPreview of mutation presence matrix:\n")
print(utils::head(result_df))
cat("\n")

### 6. Write output matrix ----------------------------------------------------

write.csv(result_df, out_matrix_csv, row.names = FALSE)

cat("Mutation presence matrix written to:\n  ",
    normalizePath(out_matrix_csv), "\n")
cat("Columns include species_short, gene_id, header, and one column per mutation.\n")
