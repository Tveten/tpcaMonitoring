init_sums <- function(x_train, x_1, t_max) {
  d <- nrow(x_train)
  m <- ncol(x_train)
  sums <- list('u' = matrix(0, nrow = d, ncol = t_max + 1),
               'v' = matrix(0, nrow = d, ncol = t_max + 1),
               's' = matrix(0, nrow = d, ncol = t_max + 1))
  sums$u[, 1] <- rowSums(x_train)
  sums$s[, 1] <- (m - 1) / m * rowVars(x_train)
  sums$v[, 1] <- sums$s[, 1] * m
  sums <- update_sums(sums, x_1, m, 1)
}

update_sums <- function(sums, next_x, m, t) {
  prev_u <- sums$u[, t]
  prev_v <- sums$v[, t]
  next_u <- prev_u + next_x

  dev <- next_x - 1 / (m + t - 1) * prev_u
  next_v <- prev_v + (m + t - 1) / (m + t) * dev^2
  next_s <- next_v / (m + t)

  sums$u[, t + 1] <- next_u
  sums$v[, t + 1] <- next_v
  sums$s[, t + 1] <- next_s
  sums
}

bartlett_corr <- function(m, t, k) {
  1 / 2 * (-(m + t) * log(m + t) +
            (m + k) * log(m + k) +
            (t - k) * log(t - k) +
            (m + t) * digamma((m + t - 1) / 2) -
            (m + k) * digamma((m + k - 1) / 2) -
            (t - k) * digamma((t - k - 1) / 2))
}

mixture_log_liks <- function(sums, m, t, w, p0) {

  log_lik_ratio <- function(m, t, k, s_full, s_pre, s_post) {
    if (any(s_post <= 0)) {
      message(paste0('Numerical error in estimating post-change variance, t = ',
                     t, ', t - k = ', t - k, '. v_full = ', 'v_pre = ', v_pre,
                     'last term = ', ((m + k) * (t - k)) / (m + t) * (avg_pre - avg_post)^2))
      log_liks[, i] <- 0
    } else
      log_liks[, i] <- (m + t) * log(s_full) -
                       (m + k) * log(s_pre) -
                       (t - k) * log(s_post)
  }

  u_full <- sums$u[, t + 1]
  v_full <- sums$v[, t + 1]
  s_full <- sums$s[, t + 1]
  d <- length(u_full)
  ks <- max(0, t - (w + 1)):(t - 2)
  log_liks <- matrix(0, nrow = d, ncol = length(ks))
  for (k in ks) {
    u_pre <- sums$u[, k + 1]
    v_pre <- sums$v[, k + 1]
    s_pre <- sums$s[, k + 1]
    avg_pre <- u_pre / (m + k)

    avg_post <- (u_full - u_pre) / (t - k)
    v_post <- v_full - v_pre - ((m + k) * (t - k)) / (m + t) * (avg_pre - avg_post)^2
    s_post <- 1 / (t - k) * v_post
    i <- k - min(ks) + 1
    log_liks[, i] <- log_lik_ratio(m, t, k, s_full, s_pre, s_post)
  }

  c <- bartlett_corr(m, t, ks)
  corr_log_liks <- log_liks %*% diag(1 / c, nrow = length(c))
  global_log_liks <- colSums(log(1 - p0 + p0 * exp(1 / 2 * corr_log_liks)))
  global_log_liks
}
