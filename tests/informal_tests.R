test_est_threshold <- function(d, r, m) {
  mu <- rep(0, d)
  Sigma <- tpca::rcov_mat(d, range_sd = c(0.5, 2))
  x <- t(MASS::mvrnorm(m, mu = mu, Sigma = Sigma))
  axes <- (d - r + 1):d
  n <- 500
  alpha <- 0.05
  tpca_threshold(x, axes, n, alpha)
}

test_tpca_arl <- function(data_dim) {
  threshold <- 200
  m <- data_dim * 5
  mu <- rnorm(data_dim, mean = 0, sd = 1 / sqrt(m))
  Sigma <- tpca::rcor_mat_est(data_dim, n = m)
  axes <- round(c(0.8, 1) * data_dim)
  n <- 500
  w <- 200
  n_sim <- 1
  tpca_arl_obj <- tpca_arl(threshold, mu, Sigma, axes, n, w, n_sim)
  plot(colMeans(tpca_arl_obj))
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
  # 750 milliseconds on data_dim = 100 after isolating work in functions.
  # 236 milliseconds on data_dim = 100 with Rcpp log_lik
  # 223 milliseconds on data_dim = 100 with Rcpp log_lik and update_sums
}

test_llC <- function(data_dim) {
  set.seed(30)
  m <- data_dim * 5
  n <- data_dim * 5
  w <- round(n / 2)
  mu <- rnorm(data_dim, mean = 0, sd = 1 / sqrt(m))
  Sigma <- tpca::rcor_mat_est(data_dim, n = m)
  x_train <- t(MASS::mvrnorm(m, mu = mu, Sigma = Sigma))
  x <- t(MASS::mvrnorm(n, mu = mu, Sigma = Sigma))

  sums <- init_sums(x_train, x[, 1], n)
  sumsC <- init_sums(x_train, x[, 1], n)
  for(t in 2:n) {
    sums <- update_sums(sums, x[, t], m, t)
    sumsC <- update_sumsC(sums, x[, t], m, t)
    print(c(pryr::address(sum2C), pryr::refs(sum2C)))
    log_liks <- mixture_log_liks(sums, m, t, w, 1)
    log_liksC <- mixture_log_liksC(sumsC, m, t, w, 1)
  }
  # print(cumvars_rev(cbind(x_train, x)))
  print(log_liks)
  print(log_liksC)
}

benchmark_llC <- function(data_dim) {
  set.seed(20)
  m <- 200
  n <- 300
  w <- 200
  mu <- rnorm(data_dim, mean = 0, sd = 1 / sqrt(m))
  Sigma <- tpca::rcor_mat_est(data_dim, n = m)
  x_train <- t(MASS::mvrnorm(m, mu = mu, Sigma = Sigma))
  x <- t(MASS::mvrnorm(n, mu = mu, Sigma = Sigma))

  microbenchmark::microbenchmark({
    sums <- init_sums(x_train, x[, 1], n)
    for(t in 2:n) {
      sums <- update_sums(sums, x[, t], m, t)
      log_liks <- mixture_log_liks(sums, m, t, w, 1)
    }
  }, {
    sums <- init_sums(x_train, x[, 1], n)
    for(t in 2:n) {
      sums <- update_sums(sums, x[, t], m, t)
      log_liks <- mixture_log_liksC(sums, m, t, w, 1)
    }
  }, {
    sums <- init_sums(x_train, x[, 1], n)
    for(t in 2:n) {
      sums <- update_sums(sums, x[, t], m, t)
      log_liks <- mixture_log_liks2(sums, x, m, t, w, 1)
    }
  },
  times = 10)
}

test_llC2 <- function(data_dim, seed) {
  set.seed(seed)
  m <- 200
  n <- 1000
  w <- 200
  mu <- rnorm(data_dim, mean = 0, sd = 1 / sqrt(m))
  Sigma <- tpca::rcor_mat_est(data_dim, n = m)
  x_train <- t(MASS::mvrnorm(m, mu = mu, Sigma = Sigma))
  x <- t(MASS::mvrnorm(n, mu = mu, Sigma = Sigma))

  sums <- init_sums(x_train, x[, 1], n)
  for (t in 2:n) {
  sums <- update_sums(sums, x[, t], m, t)
  log_liks <- mixture_log_liks(sums, m, t, w, 1)
  ks <- max(0, t - (w + 1)):(t - 2)
  log_liks2 <- mixture_log_liks2(sums, x,  m, t, w, 1)
  log_liksC <- mixture_log_liksC(sums, m, t, w, 1)
  }
  print(1/2 * rowVars(x[, (n-1):n]))
}
