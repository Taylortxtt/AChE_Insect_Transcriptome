############################################################
# 01_download_AChE_sequences.R
#
# Downloads AChE-related nucleotide sequences from NCBI
# (mRNA + CDS only) for each species in species_list.
#
# Output: data/raw/<short>_raw.fasta
############################################################

if (!requireNamespace("rentrez", quietly = TRUE)) {
  install.packages("rentrez")
}
library(rentrez)

# ---- Directories ------------------------------------------------------------
raw_dir <- file.path("data", "raw")
if (!dir.exists(raw_dir)) dir.create(raw_dir, recursive = TRUE)

cat("Raw data directory:", normalizePath(raw_dir), "\n\n")

# ---- Species list -----------------------------------------------------------
species_list <- list(
  list(short = "gregaria",   org = "Schistocerca gregaria"),
  list(short = "piceifrons", org = "Schistocerca piceifrons"),
  list(short = "nitens",     org = "Schistocerca nitens"),
  list(short = "cancellata", org = "Schistocerca cancellata"),
  list(short = "serialis",   org = "Schistocerca serialis cubense"),
  list(short = "anabrus",    org = "Anabrus simplex"),
  list(short = "anopheles",  org = "Anopheles gambiae"),
  list(short = "aedes",      org = "Aedes aegypti"),
  list(short = "culex",      org = "Culex quinquefasciatus")
)

# ---- Helper: fetch sequence safely -----------------------------------------
fetch_one_fasta <- function(ncbi_id, max_tries = 5, sleep_sec = 0.5) {
  for (i in seq_len(max_tries)) {
    res <- try(
      rentrez::entrez_fetch(
        db = "nucleotide",
        id = ncbi_id,
        rettype = "fasta"
      ),
      silent = TRUE
    )
    if (!inherits(res, "try-error")) return(res)
    Sys.sleep(sleep_sec)
  }
  return(NULL)
}

# ---- Main loop --------------------------------------------------------------
cat("Starting AChE sequence downloads from NCBI...\n")
cat("------------------------------------------------------------\n\n")

for (sp in species_list) {
  cat("Species:", sp$short, "\n")
  
  term <- paste0(
    '("', sp$org, '"[Organism]) AND (acetylcholinesterase[All Fields] ',
    'OR AChE[All Fields] OR ace1[All Fields] OR ace2[All Fields]) AND ',
    '(mRNA[Title] OR cds[Title])'
  )
  
  search_res <- rentrez::entrez_search(
    db = "nucleotide",
    term = term,
    retmax = 200
  )
  
  ids <- search_res$ids
  cat("Hits found:", length(ids), "\n")
  
  out_path <- file.path(raw_dir, paste0(sp$short, "_raw.fasta"))
  writeLines("", out_path)
  
  if (length(ids) == 0) {
    cat("  No records found — skipping.\n\n")
    next
  }
  
  cat("  Downloading sequences...\n")
  total_ok <- 0
  
  for (n in ids) {
    fa <- fetch_one_fasta(n)
    if (!is.null(fa)) {
      write(fa, file = out_path, append = TRUE)
      total_ok <- total_ok + 1
    }
  }
  
  cat("  Completed:", total_ok, "saved.", "\n\n")
}

cat("------------------------------------------------------------\n")
cat("Done!\n")
