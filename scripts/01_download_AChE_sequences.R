############################################################
# 01_download_AChE_sequences.R
#
# Goal
# ----
# Download acetylcholinesterase (AChE) mRNA/CDS sequences
# from NCBI for all project species into data/raw/.
#
# Outputs
# -------
# One FASTA file per species:
#   data/raw/<short_name>_raw.fasta
#
# Notes
# -----
# * Designed to be run from the RStudio project root
#   (the folder that contains AChE_Project.Rproj).
# * Uses the {rentrez} package (NCBI E-utilities).
# * Later scripts will:
#     - choose the longest transcript per
#     - align sequences and build the tree.
############################################################

### 0. Packages --------------------------------------------------------------

if (!requireNamespace("rentrez", quietly = TRUE)) {
  install.packages("rentrez")
}
library(rentrez)

### 1. Project folders --------------------------------------------------------

# All raw NCBI downloads go into: data/raw/
raw_dir <- file.path("data", "raw")

if (!dir.exists(raw_dir)) {
  dir.create(raw_dir, recursive = TRUE)
}

cat("Raw data directory:", normalizePath(raw_dir), "\n\n")

### 2. Species list -----------------------------------------------------------

# Each entry has:
#   short_name    → used in file names
#   organism_term → how we tell NCBI which organism to search
#
# Currently includes:
#   * Three Schistocerca species (locusts)
#   * One cricket
#   * Three mosquitoes
#
# If you want to add more species later, just append another
# list(short_name = "...", organism_term = "Genus species").
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

# If you have an NCBI API key, you can paste it below to speed
# up requests.
#
# Get a key (free) by making an NCBI account and requesting one.
# Then uncomment and replace "YOUR_API_KEY_HERE".
#
# entrez_key("YOUR_API_KEY_HERE")

### 4. Helper: download AChE sequences for one species ------------------------

download_AChE_for_species <- function(short_name, organism_term) {
  cat("--------------------------------------------------\n")
  cat("Species:", organism_term, "(", short_name, ")\n")
  
  # Output FASTA file for this species
  out_file <- file.path(raw_dir, paste0(short_name, "_raw.fasta"))
  
  # If file already exists, skip by default so we don't constantly
  # re-download from NCBI. Delete the file if you need a fresh copy.
  if (file.exists(out_file)) {
    cat("File already exists, skipping download:\n  ", out_file, "\n\n")
    return(invisible(out_file))
  }
  
  # Search term:
  # * Limit to the organism
  # * Look for acetylcholinesterase / ache / ace
  # * Prefer mRNA/CDS records
  search_term <- paste0(
    "(", organism_term, "[Organism]) AND (",
    "acetylcholinesterase[All Fields] OR ",
    "ache[All Fields] OR ",
    "ace[All Fields]",
    ") AND (mRNA[Title] OR cds[Title] OR \"complete cds\"[Title])"
  )
  
  cat("NCBI search term:\n  ", search_term, "\n")
  
  # Search NCBI nuccore
  search_res <- rentrez::entrez_search(
    db         = "nuccore",
    term       = search_term,
    retmax     = 500,
    use_history = TRUE
  )
  
  n_hits <- search_res$count
  cat("Hits found:", n_hits, "\n")
  
  if (n_hits == 0) {
    warning("No sequences found for: ", organism_term)
    return(invisible(NULL))
  }
  
  # Fetch all matching records as FASTA using the web history
  fasta_txt <- rentrez::entrez_fetch(
    db         = "nuccore",
    web_history = search_res$web_history,
    rettype    = "fasta",
    retmode    = "text"
  )
  
  # Basic sanity check
  if (!grepl(">", fasta_txt, fixed = TRUE)) {
    warning("No FASTA headers ('>') found for: ", organism_term)
  }
  
  # Write to disk
  writeLines(fasta_txt, con = out_file)
  cat("Saved FASTA to:\n  ", normalizePath(out_file), "\n\n")
  
  invisible(out_file)
}

### 5. Loop over species ------------------------------------------------------

cat("Starting AChE sequence downloads from NCBI...\n")

for (sp in species_list) {
  download_AChE_for_species(
    short_name    = sp$short_name,
    organism_term = sp$organism_term
  )
}

cat("--------------------------------------------------\n")
cat("Done! Check data/raw/ for the downloaded FASTA files.\n")
