############################################################
# run_pipeline.R
#
# One-click runner for the AChE Insect Transcriptome
# pipeline. Comment out any steps you don't want.
############################################################

run_script <- function(script_name) {
  path <- file.path("scripts", script_name)
  cat("\n==================================================\n")
  cat("Running", path, "\n")
  cat("==================================================\n\n")
  source(path)
}

# ---- Core pipeline ---------------------------------------------------------
# Comment any of these if you do not want to run them as
# part of the one-click pipeline.

run_script("01_download_AChE_sequences.R")
run_script("02_extract_longest_isoforms.R")
run_script("03_merge_longest_fastas.R")
run_script("09_classify_genes_blast.R")   # BLAST classification
run_script("04_run_alignment.R")
run_script("05_build_tree_fasttree.R")
run_script("05b_build_tree_iqtree.R")    # if you have IQ-TREE
run_script("06_plot_tree.R")
run_script("06b_plot_tree_iqtree.R")     # if you have IQ-TREE
run_script("07_summarize_copy_number.R")
run_script("08_screen_mutations.R")

cat("\nAll selected pipeline steps completed.\n")
cat("Check the data/ and exports/ folders for outputs.\n")
