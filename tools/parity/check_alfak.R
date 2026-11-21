# Validate files written by the Python ALFA-K exporter using R readRDS.
args <- commandArgs(TRUE)
if (length(args) < 2L) stop("Usage: Rscript tools/parity/check_alfak.R R_LIBRARY OUTPUT_DIR...")
library(karst, lib.loc = normalizePath(args[[1L]]))

for (folder in args[-1L]) {
  read <- function(name) readRDS(file.path(folder, paste0(name, ".Rds")))
  landscape <- read("landscape")
  samples <- read("landscape_posterior_samples")
  bootstrap <- read("bootstrap_res")
  validation <- read("xval")
  models <- read("landscape_data")
  stopifnot(
    identical(names(landscape), c("k", "mean", "median", "sd", "fq", "nn")),
    identical(colnames(bootstrap$final_fitness), landscape$k[landscape$fq]),
    identical(colnames(bootstrap$nn_fitness), landscape$k[landscape$nn]),
    all(abs(rowSums(bootstrap$final_frequencies) - 1) < 1e-12),
    is.data.frame(attr(bootstrap, "bootstrap_retry_log")),
    length(attr(bootstrap, "bootstrap_retry_summary")) == 5L,
    identical(colnames(validation$tmp), c("test_f", "est_f")),
    isTRUE(all.equal(unname(rowMeans(samples)), landscape$mean, tolerance = 1e-12)),
    isTRUE(all.equal(unname(apply(samples, 1L, sd)), landscape$sd, tolerance = 1e-12)),
    isTRUE(all.equal(as.numeric(predict(models$mean, landscape$k)), landscape$mean, tolerance = 1e-12))
  )
  cat(folder, ": read and validated five RDS files\n", sep = "")
}
