test_est_threshold <- function(d, r, m) {
  x <- list('d' = d, 'm' = m)
  axes <- (d - r + 1):d
  arl <- 500
  est_tpca_threshold(x, axes, arl)
}
