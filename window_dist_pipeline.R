#!/usr/bin/env Rscript

### this script takes whole genome gds file with SNPs (you can create it from vcf) and create pairwise distance matrix
# for non-overlap windows of given length in bp. It does calculation parallelly. Also you get whole-genome pairwise dist matrix.
# Main functions for this script are in non_overlap_windows_functions.R, so keep it in the same directory as this script.
suppressPackageStartupMessages({
  library(SeqArray)
  library(data.table)
})

# ---- Argument parsing ----
args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag) {
  idx <- which(args == flag)
  if (length(idx) == 0 || idx == length(args)) return(NULL)
  return(args[idx + 1])
}

gds_path <- get_arg("--seq-gds")
threads  <- as.integer(get_arg("--threads"))
win_size <- as.integer(get_arg("--win-size"))

# ---- Validate input ----
if (is.null(gds_path) || is.null(threads) || is.null(win_size)) {
  stop("Usage: script.R --seq-gds <file.gds> --threads <int> --win-size <int>")
}

if (!file.exists(gds_path)) {
  stop("GDS file does not exist: ", gds_path)
}

# ---- Resolve script directory (robust sourcing) ----
get_script_path <- function() {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  fileArg <- "--file="
  match <- grep(fileArg, cmdArgs)

  if (length(match) > 0) {
    return(normalizePath(sub(fileArg, "", cmdArgs[match])))
  } else {
    stop("Cannot determine script path (are you running via Rscript?)")
  }
}

script_path <- get_script_path()
script_dir  <- dirname(script_path)

# ---- Source helper functions ----
source(file.path(script_dir, "non_overlap_windows_functions.R"))

# ---- Prepare names ----
base_name <- sub("\\.gds.*$", "", basename(gds_path))
win_size_f <- format_bp(win_size)

win_dt_file <- paste0(base_name, "_win_dt_", win_size_f, ".rds")
dist_file   <- paste0(base_name, "_cust_dist.rds")

# ---- Main pipeline ----
cat("Opening GDS...\n")
gds <- seqOpen(gds_path)

samples <- seqGetData(gds, "sample.id")
n <- length(samples)

cat("Splitting into windows...\n")
win_dt <- split_to_windows(
  gds,
  window_size = win_size,
  min_n_var = 3000
)

cat("Chunking windows for parallel processing...\n")
win_chunks <- split_win_dt(win_dt, threads)

seqClose(gds)

cat("Running parallel distance calculations...\n")
res <- run_parallel_windows(
  win_chunks,
  gds_path,
  threads,
  samples,
  n
)

cat("Combining results...\n")
win_dt <- combine_win_results(res)

# ---- Save window table ----
cat("Saving window data to:", win_dt_file, "\n")
saveRDS(win_dt, file = win_dt_file)

# ---- Global distance ----
cat("Computing global distance matrix...\n")
global_dist <- calc_global_dist(win_dt, samples, n)

cat("Saving global distance to:", dist_file, "\n")
saveRDS(global_dist[[1]], file = dist_file)

cat("Done.\n")
