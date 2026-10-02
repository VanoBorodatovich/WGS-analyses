library(SeqArray)
library(data.table)
library(ape)

# these functions work with SeqArray gds data format.
### FUNCTIONS ###
split_to_windows <- function(gds, window_size = 1e5, min_n_var = 1000) {
  # extract required vectors
  dt <- data.table(
    chr = seqGetData(gds, "chromosome"),
    pos = seqGetData(gds, "position"),
    vid = seqGetData(gds, "variant.id")
  )
  # ensure proper ordering
  setkey(dt, chr, pos)
  # assign windows per chromosome
  dt[, window_id := floor((pos - 1) / window_size), by = chr]
  # aggregate into window table
  win_dt <- dt[, .(
    wind_st   = window_id[1] * window_size + 1,
    wind_end  = (window_id[1] + 1) * window_size,
    var_start = vid[1],
    var_end   = vid[.N],
    n_var     = .N
  ), by = .(chr, window_id)]
  # filter by minimum number of variants
  win_dt <- win_dt[n_var >= min_n_var]
  # free mem
  rm(dt)
  gc()
  return(win_dt[])
}



calc_window_distances <- function(gds, win_dt, samples, n) {
  win_dt[, dist := vector("list", .N)]
  win_dt[, n_snps := vector("list", .N)]
  for (i in seq_len(nrow(win_dt))) {
    v_start <- win_dt$var_start[i]
    v_end   <- win_dt$var_end[i]
    seqSetFilter(gds, variant.sel = v_start:v_end)
    dist_mat <- matrix(0, n, n)
    dosage_mat <- seqGetData(gds, "$dosage")
    results <- calc_dist_from_matrix(dosage_mat, samples, n)
    win_dt$dist[[i]] <- results[[1]]
    win_dt$n_snps[[i]] <- results[[2]]
    seqResetFilter(gds)
    if (i %% 10 == 0) {
      message("Processed window ", i, "/", nrow(win_dt))
    }
  }
  return(win_dt)
}



calc_dist_from_matrix <- function(dosage_mat, samples, n) {
  # Initialize an empty square matrix
  dist_m <- matrix(0, nrow = n, ncol = n)
  n_snps <- matrix(0, nrow = n, ncol = n)
  colnames(dist_m) <- rownames(dist_m) <- samples
  colnames(n_snps) <- rownames(n_snps) <- samples
  # Iterate through unique pairs (vectorized across SNPs)
  pairs <- combn(n, 2)
  for (i in 1:ncol(pairs)) {
    s1 <- pairs[1, i]
    s2 <- pairs[2, i]
    # Manhattan distance: sum of absolute differences
    # If both are 1 (het), |1-1| = 0.
    # If 0 and 2 (homoz opposites), |0-2| = 2.
    valid_snps <- !is.na(dosage_mat[s1,]) & !is.na(dosage_mat[s2,])
    n_valid <- sum(valid_snps)
    if (n_valid == 0) {
      d <- NA
    } else {
      d <- sum(abs(dosage_mat[s1, ] - dosage_mat[s2, ]), na.rm = TRUE) / (2 * n_valid)
    }
    dist_m[s1, s2] <- d
    dist_m[s2, s1] <- d
    n_snps[s1, s2] <- n_valid
    n_snps[s2, s1] <- n_valid
  }
  res <- list(dist_m, n_snps)
  return(res)
}



split_win_dt <- function(win_dt, n_chunks) {
  stopifnot(data.table::is.data.table(win_dt))
  stopifnot(n_chunks >= 1)
  n <- nrow(win_dt)
  if (n_chunks > n) n_chunks <- n
  # balanced split
  chunk_id <- rep(seq_len(n_chunks), length.out = n)
  win_dt[, chunk_id := chunk_id]
  chunks <- split(win_dt, by = "chunk_id", keep.by = FALSE)
  win_dt[, chunk_id := NULL]
  # optional: name them like tmp_win_dt_1, ...
  names(chunks) <- paste0("tmp_win_dt_", seq_along(chunks))
  return(chunks)
}



run_parallel_windows <- function(win_chunks, gds_file, n_cores = 1, samples, n) {
  stopifnot(file.exists(gds_file))
  # Define worker outside or ensure it's exported
  worker_fun <- function(chunk, gds_path, samples, n) {
    # Create a message
    # Get the ID of the current process to see who is talking
    pid <- Sys.getpid()
    msg <- paste0(Sys.time(), " [PID ", pid, "]: Processing ", nrow(chunk), " windows.\n")
    cat(msg, file = "parallel_log.txt", append = TRUE)
    # It is safer to load libraries inside if not using a formal package
    suppressPackageStartupMessages({
      library(SeqArray)
      library(data.table)
    })
    setDTthreads(1)
    gds <- seqOpen(gds_path)
    # Ensure file closes even if calc_window_distances fails
    res <- tryCatch({
      calc_window_distances(gds, chunk, samples, n)
    }, finally = {
      seqClose(gds)
    })
    return(res)
  }
  if (.Platform$OS.type == "unix") {
    res_list <- parallel::mclapply(
      win_chunks,
      worker_fun,
      gds_path = gds_file,
      samples = samples,
      n = n,
      mc.cores = n_cores
    )
  } else {
    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    # Need to export the distance function to the workers!
    parallel::clusterExport(cl, varlist = c("calc_window_distances", "calc_dist_from_matrix"))
    res_list <- parallel::parLapply(
      cl,
      win_chunks,
      worker_fun,
      gds_path = gds_file,
      samples = samples,
      n = n
    )
  }
  return(res_list)
}



combine_win_results <- function(res_list) {
  stopifnot(is.list(res_list))
  combined <- data.table::rbindlist(res_list, use.names = TRUE, fill = TRUE)
  data.table::setorder(combined, window_id)
  return(combined)
}



fill_trees <- function(win_dt, outgroup = NULL) {
  stopifnot(is.data.table(win_dt))
  win_dt[, tree := lapply(dist, function(m) {
    tr <- bionj(as.dist(m))
    if (!is.null(outgroup)) {
      root(tr, outgroup = outgroup, resolve.root = TRUE)
    } else {
      midpoint(tr)
    }
  })]
  return(win_dt)
}



format_bp <- function(x) {
  ifelse(x >= 1e6,
         paste0(x / 1e6, "m"),
         ifelse(x >= 1e3,
                paste0(x / 1e3, "kb"),
                as.character(x)))
}



calc_global_dist <- function(win_dt, samples, n) {
  # initialize accumulators
  sum_weighted_dist <- matrix(0, n, n)
  sum_snps <- matrix(0, n, n)
  colnames(sum_snps) <- rownames(sum_snps) <- samples
  # iterate over windows
  for (i in seq_len(nrow(win_dt))) {
    d_mat <- win_dt$dist[[i]]
    n_mat <- win_dt$n_snps[[i]]
    # skip empty windows
    if (is.null(d_mat) || is.null(n_mat)) next
    # accumulate
    valid <- !is.na(d_mat) & n_mat > 0
    sum_weighted_dist[valid] <- sum_weighted_dist[valid] + d_mat[valid] * n_mat[valid]
    sum_snps[valid] <- sum_snps[valid] + n_mat[valid]
  }
  # compute final weighted mean
  global_dist <- matrix(NA, n, n)
  valid <- sum_snps > 0
  global_dist[valid] <- sum_weighted_dist[valid] / sum_snps[valid]
  colnames(global_dist) <- rownames(global_dist) <- samples
  # enforce symmetry and zero diagonal
  diag(global_dist) <- 0
  return(list(
    dist = global_dist,
    n_snps = sum_snps
  ))
}
