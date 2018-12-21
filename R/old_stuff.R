cumvars_rev <- function(x) {
  d <- nrow(x)
  m <- ncol(x)
  x <- x[, m:1]
  cumvars <- matrix(0, nrow = d, ncol = m - 1)
  for (i in 2:m) {
    cumvars[, i - 1] <- (i - 1) / i * rowVars(x[, 1:i])
  }
  cumvars
}

mixture_log_liks2 <- function(sums, x, m, t, w, p0) {
  s_full <- sums$s[, t + 1]
  d <- length(s_full)
  ks <- max(0, t - (w + 1)):(t - 2)
  log_liks <- matrix(0, nrow = d, ncol = length(ks))
  for (k in ks) {
    s_pre <- sums$s[, k + 1]
    s_post <- (t - k - 1) / (t - k) * rowVars(x[, (k + 1):t])
    if (t == 1000 && k == 998) print(s_post)
    i <- k - min(ks) + 1
    log_liks[, i] <- (m + t) * log(s_full) -
                     (m + k) * log(s_pre) -
                     (t - k) * log(s_post)

  }
  c <- bartlett_corr(m, t, ks)
  corr_log_liks <- log_liks %*% diag(1 / c, nrow = length(c))
  global_log_liks <- colSums(log(1 - p0 + p0 * exp(1 / 2 * corr_log_liks)))
  global_log_liks
}
