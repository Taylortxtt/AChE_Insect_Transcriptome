############################################################
# config.R
#
# Central configuration for the AChE Insect Transcriptome
# project: paths to external tools etc.
############################################################

# ---- MAFFT ------------------------------------------------
# Multiple sequence alignment
# Find with `which mafft` in Terminal.
MAFFT_PATH <- "/opt/homebrew/bin/mafft"   # <-- change if needed

# ---- FastTree ---------------------------------------------
# Approximate maximum-likelihood tree
# Find with `which FastTree` in Terminal
FASTTREE_PATH <- "/usr/local/bin/FastTree"   # <-- change if needed

# ---- IQ-TREE (optional) -----------------------------------
# If you installed IQ-TREE2; otherwise this can stay unused
# Find with `which iqtree2` in Terminal
IQTREE_PATH <- "/opt/homebrew/bin/iqtree2"   # <-- change if needed

# ---- BLASTN -----------------------------------------------
# NCBI BLAST+ nucleotide search
# Install with:  brew install blast
# Then find with: which blastn
BLASTN_PATH <- "/opt/homebrew/bin/blastn"    # <-- change if needed

# Path PREFIX of the BLAST database you build with makeblastdb
# (don't add .nsq/.nin/.nhr)
BLAST_DB <- "data/blast/ache_ref"
