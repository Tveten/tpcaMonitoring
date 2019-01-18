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

update_sums <- function(sums, next_x, m, t) {
  prev_u <- sums$u[, t]
  prev_v <- sums$v[, t]
  dev2 <- (next_x - prev_u / (m + t - 1))^2
  sums$u[, t + 1] <- prev_u + next_x
  sums$v[, t + 1] <- prev_v + (m + t - 1) / (m + t) * dev2
  sums
}

mixture_log_liks <- function(sums, m, t, w, p0) {

  log_lik_ratio <- function(m, t, k, s_full, s_pre, s_post) {
    if (any(s_post <= 0)) {
      message(paste0('Numerical error in estimating post-change variance, t = ',
                     t, ', t - k = ', t - k, ', s_post = ',
                     paste(s_post[s_post <= 0], collapse = ', '), '\n'))
      log_liks[, i] <- 0
    } else
      log_liks[, i] <- (m + t) * log(s_full) -
                       (m + k) * log(s_pre) -
                       (t - k) * log(s_post)
  }

  u_full <- sums$u[, t + 1]
  v_full <- sums$v[, t + 1]
  s_full <- v_full / (m + t)
  d <- length(u_full)
  ks <- max(0, t - (w + 1)):(t - 2)
  log_liks <- matrix(0, nrow = d, ncol = length(ks))
  for (k in ks) {
    u_pre <- sums$u[, k + 1]
    v_pre <- sums$v[, k + 1]
    s_pre <- v_pre / (m + k)
    avg_pre <- u_pre / (m + k)

    avg_post <- (u_full - u_pre) / (t - k)
    v_post <- v_full - v_pre - ((m + k) * (t - k)) / (m + t) * (avg_pre - avg_post)^2
    s_post <- 1 / (t - k) * v_post
    if (t == 1000 && k == 998) print(s_post)
    i <- k - min(ks) + 1
    log_liks[, i] <- log_lik_ratio(m, t, k, s_full, s_pre, s_post)
  }

  c <- bartlett_corr(m, t, ks)
  corr_log_liks <- log_liks %*% diag(1 / c, nrow = length(c))
  global_log_liks <- colSums(log(1 - p0 + p0 * exp(1 / 2 * corr_log_liks)))
  global_log_liks
}
