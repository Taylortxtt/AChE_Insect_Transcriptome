############################################################
# 01_download_AChE_sequences.R
#
# Goal
# ----
# Download acetylcholinesterase (AChE) mRNA/CDS sequences
# from NCBI for all project species into data/raw/.
#
# Design
# ------
# * Uses rentrez (NCBI E-utilities).
# * One helper function per species.
# * Outputs one FASTA file per species:
#       data/raw/<short_name>_raw.fasta
############################################################

### 0. Packages --------------------------------------------------------------

if (!requireNamespace("rentrez", quietly = TRUE)) {
  install.packages("rentrez")
}
library(rentrez)

### 1. Project folders --------------------------------------------------------

# Assumes you are running this from the project root (AChE_Project).
# All raw NCBI downloads go into: data/raw/

raw_dir <- file.path("data", "raw")

if (!dir.exists(raw_dir)) {
  dir.create(raw_dir, recursive = TRUE)
}

cat("Raw data directory:", normalizePath(raw_dir), "\n")

### 2. Species list -----------------------------------------------------------

# Each entry has:
#   short_name    → used in file names
#   organism_term → how we tell NCBI which organism to search
species_list <- list(
  list(short_name = "gregaria",   organism_term = "Schistocerca gregaria"),
  list(short_name = "cancellata", organism_term = "Schistocerca cancellata"),
  list(short_name = "piceifrons", organism_term = "Schistocerca piceifrons"),
  list(short_name = "anabrus",    organism_term = "Anabrus simplex"),
  list(short_name = "anopheles",  organism_term = "Anopheles gambiae"),
  list(short_name = "aedes",      organism_term = "Aedes aegypti"),
  list(short_name = "culex",      organism_term = "Culex quinquefasciatus")
)

### 3. (Optional) NCBI API key -----------------------------------------------
# If you have an NCBI API key, uncomment and paste it here to speed things up:
#
# entrez_key("YOUR_API_KEY_HERE")

### 4. Helper: download AChE sequences for one species -----------------------

download_AChE_for_species <- function(short_name, organism_term) {
  cat("--------------------------------------------------\n")
  cat("Species:", organism_term, "(", short_name, ")\n")
  
  # Search for AChE using a fairly generous set of name variants.
  search_term <- paste0(
    "(", organism_term, "[Organism]) AND ",
    "(",
    "acetylcholinesterase[All Fields] OR ",
    "AChE[All Fields] OR ",
    "ace-1[All Fields] OR ",
    "ace1[All Fields] OR ",
    "ace-2[All Fields] OR ",
    "ace2[All Fields]",
    ") AND ",
    "(mRNA[Filter] OR cds[Filter])"
  )
  
  cat("NCBI search term:\n  ", search_term, "\n")
  
  # Search the nucleotide database
  search_res <- rentrez::entrez_search(
    db     = "nuccore",
    term   = search_term,
    retmax = 500
  )
  
  n_hits <- search_res$count
  cat("Hits found:", n_hits, "\n")
  
  if (n_hits == 0) {
    warning("No sequences found for: ", organism_term)
    return(invisible(NULL))
  }
  
  # Fetch all matching records as FASTA
  fasta_raw <- rentrez::entrez_fetch(
    db      = "nuccore",
    id      = search_res$ids,
    rettype = "fasta",
    retmode = "text"
  )
  
  # Save to data/raw/<short_name>_raw.fasta
  out_file <- file.path(raw_dir, paste0(short_name, "_raw.fasta"))
  writeLines(fasta_raw, con = out_file)
  
  cat("Saved FASTA to:", out_file, "\n")
  invisible(out_file)
}

### 5. Loop over species ------------------------------------------------------

for (sp in species_list) {
  download_AChE_for_species(
    short_name    = sp$short_name,
    organism_term = sp$organism_term
  )
}

cat("--------------------------------------------------\n")
cat("Done! Check data/raw/ for the downloaded FASTA files.\n")
