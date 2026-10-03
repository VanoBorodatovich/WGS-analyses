# This script perform PCA for samples you choose in your vcf, then project on this PCA all other samples using only genotyped position 
# in ech of them. 
library(vcfR)

vcf <- read.vcfR("data/lemmus_autosomes_clean_popan.vcf.gz", verbose = FALSE)
sites <- vcf@fix
gt <- extract.gt(vcf, IDtoRowNames = FALSE)
code <- c("0/0" = 0, "0/1" = 1, "1/0" = 1, "1/1" = 2)
if (any(!is.na(gt) & !(gt %in% names(code)))) {
  stop("Expected diploid GT calls at biallelic SNPs.")
}
G <- matrix(unname(code[as.vector(gt)]), nrow(gt), ncol(gt))
colnames(G) <- colnames(gt)
rownames(G) <- paste(sites[, "CHROM"], sites[, "POS"],
                     sites[, "REF"], sites[, "ALT"], sep = ":")
positions <- paste(sites[, "CHROM"], sites[, "POS"], sep = ":")
if (anyDuplicated(positions)) stop("Duplicate SNP positions in VCF.")
if (ncol(G) == 0 || anyDuplicated(colnames(G))) stop("Invalid sample names.")

### STEP2 subset genotype matrix
samples_bad <- c("Verkhoyansk","Sikhote","6328_Yuribey_Gulf","Ob","Yugorsky2","Indigirka3","GBritain","NUrals")
samples_good <- setdiff(colnames(G), samples_bad)
high = G[, samples_good, drop = FALSE]
low = G[, samples_bad, drop = FALSE]

### STEP3
# Recalculate frequencies AFTER subsetting to reference samples.
min_maf = 0.05
max_missing = 0.05
n_pc = 6
p <- rowMeans(high, na.rm = TRUE) / 2
missing <- rowMeans(is.na(high))
keep <- is.finite(p) & p > 0 & p < 1 &
  pmin(p, 1 - p) >= min_maf & missing <= max_missing
keep[is.na(keep)] <- FALSE
high <- high[keep, , drop = FALSE]
p <- p[keep]
if (nrow(high) <= n_pc) stop("Too few SNPs after reference filtering.")

# Fixed SNPs were removed above, so every scaling factor must be positive.
mu <- 2 * p
sigma <- sqrt(2 * p * (1 - p))
stopifnot(all(is.finite(mu)), all(is.finite(sigma)), all(sigma > 0))

Z <- sweep(high, 1, mu, "-")
Z <- sweep(Z, 1, sigma, "/")
Z[is.na(high)] <- 0  # Mean-impute ONLY missing reference genotypes.
stopifnot(all(is.finite(Z)))

# prcomp expects samples in rows and SNPs in columns.
pca <- prcomp(t(Z), center = FALSE, scale. = FALSE, rank. = n_pc)
if (pca$sdev[n_pc] < pca$sdev[1] * 1e-8) stop("Insufficient rank: reduce n_pc.")
message("Reference PCA: ", ncol(G), " samples, ", nrow(G), " SNPs.")

model <- list(scores = pca$x, loadings = pca$rotation,
            snps = rownames(high), mean = mu, scale = sigma,
            variance_pct = 100 * pca$sdev^2 / sum(pca$sdev^2))



### STPEP4 Fit coordinates for each remaining sample using observed SNPs only.
# Apply the reference SNP panel and its exact SNP order.
min_snps <- 100000
index <- match(model$snps, rownames(low))
if (anyNA(index)) stop("Reference SNPs absent from projection matrix.")
low <- low[index, , drop = FALSE]
if (any(!is.na(G) & !(G %in% 0:2))) stop("Invalid genotype dosage.")

n_pc <- ncol(model$loadings)
scores <- matrix(NA_real_, ncol(low), n_pc,
                 dimnames = list(colnames(low), colnames(model$loadings)))
qc <- data.frame(sample = colnames(low), n_observed = colSums(!is.na(low)),
                 status = rep("too_few_SNPs", ncol(low)), row.names = NULL)

for (i in seq_len(ncol(low))) {
  observed <- which(!is.na(low[, i]))
  if (length(observed) < max(min_snps, n_pc)) next
  
  y <- (low[observed, i] - model$mean[observed]) / model$scale[observed]
  L <- model$loadings[observed, , drop = FALSE]
  stopifnot(all(is.finite(y)), all(is.finite(L)))
  
  # Check whether these SNPs allow a stable fit of all requested PCs.
  q <- qr(L)
  condition <- kappa(L, exact = TRUE)
  if (q$rank < n_pc || !is.finite(condition) || condition > 1e6) {
    qc$status[i] <- "unstable_fit"
    next
  }
  
  # Least squares, no intercept; equivalent to qr.solve(L, y).
  scores[i, ] <- qr.coef(q, y)
  qc$status[i] <- "ok"
}
if (any(qc$status != "ok")) warning("Some projections failed; inspect $qc.")
projected <- list(scores = scores, qc = qc)
