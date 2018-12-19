test_est_threshold <- function(d, r, m) {
  mu <- rep(0, d)
  Sigma <- tpca::rcov_mat(d, range_sd = c(0.5, 2))
  x <- t(MASS::mvrnorm(m, mu = mu, Sigma = Sigma))
  axes <- (d - r + 1):d
  n <- 500
  alpha <- 0.05
  est_tpca_threshold(x, axes, n, alpha)
}

test_tpca_arl <- function(data_dim) {
  threshold <- 200
  m <- data_dim * 5
  mu <- rnorm(data_dim, mean = 0, sd = 1 / sqrt(m))
  Sigma <- tpca::rcor_mat_est(data_dim, n = m)
  axes <- round(c(0.8, 1) * data_dim)
  n <- 500
  w <- 200
  n_sim <- 100
  plot(tpca_arl(threshold, mu, Sigma, axes, n, w, n_sim))
}

benchmark_tpca_arl <- function(data_dim) {
  threshold <- 200
  m <- data_dim * 5
  mu <- rnorm(data_dim, mean = 0, sd = 1 / sqrt(m))
  Sigma <- tpca::rcor_mat_est(data_dim, n = m)
  axes <- round(c(0.8, 0.9, 0.95, 1) * data_dim)
  n <- 500
  w <- 200
  n_sim <- 1
  microbenchmark::microbenchmark(
    tpca_arl(threshold, mu, Sigma, axes, n, w, n_sim)
  )
  # 500 milliseconds on data_dim = 100.
}
