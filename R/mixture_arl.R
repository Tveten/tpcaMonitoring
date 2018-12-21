#' Estimates average run length (ARL) for monitoring by the mixture statistic.
#'
#' Description
#'
#' Details
#'
#' @param threshold A numeric specifying the threshold value for when a change
#'   is declared.
#' @param mu_x A mean vector estimated from training data, representing the null distribution.
#' @param Sigma_x A covariance matrix estimated from training data, representing the null distribution.
#' @param n The number of observations to monitor for an estimate of the ARL.
#'   See details.
#' @param w The window size. Number of recent time-points to consider for a
#'   change.
#' @param n_sim The number of simulations to base the estimate on.
#'
#' @return An estimate of the average run length.
#'
#' @export

mixture_arl <- function(threshold, m, d, p0, n, w, n_sim) {
  # TODO: Parallelize.
  n_cores <- parallel::detectCores()
  c1 <- parallel::makeCluster(n_cores - 1, outfile = '', type = 'PSOCK')
  doParallel::registerDoParallel(c1)
  mu_x <- rep(0, d)
  Sigma_x <- diag(1, nrow = d)
  x_train <- gen_norm_data(m, mu_x, Sigma_x)
  `%dopar%` <- foreach::`%dopar%`
  run_lengths <- foreach::foreach(b = 1:n_sim, .combine = 'c') %dopar% {
    x <- gen_norm_data(n, mu_x, Sigma_x)
    t <- 1
    sums <- init_sums(x_train, x[, t], n)
    detection_stat <- 0
    while ((detection_stat < threshold) & (t < n)) {
      t <- t + 1
      sums <- update_sumsC(sums, x[, t], m, t)
      log_liks <- mixture_log_liksC(sums, m, t, w, 1)
      detection_stat <- max(log_liks)
    }
    t
  }
  parallel::stopCluster(c1)
  arl_est <- n / mean(as.numeric(run_lengths < n))
  arl_est
}
