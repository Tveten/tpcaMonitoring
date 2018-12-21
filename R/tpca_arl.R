#' Estimates average run length (ARL) for monitoring by tpca
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

tpca_arl <- function(threshold, mu_x, Sigma_x, axes, n, w, n_sim) {
  m <- attr(Sigma_x, 'n_obs')
  d <- length(mu_x)
  r <- length(axes)
  # TODO: Parallelize.
  n_cores <- parallel::detectCores()
  c1 <- parallel::makeCluster(n_cores - 1, outfile = '', type = 'PSOCK')
  doParallel::registerDoParallel(c1)
  `%dopar%` <- foreach::`%dopar%`
  run_lengths <- foreach::foreach(b = 1:n_sim, .combine = 'c') %dopar% {
    z_train <- boot_z_train(m, mu_x, Sigma_x, axes)
    z <- gen_norm_data(n, mu = z_train$mu, Sigma = diag(z_train$sigma2))
    t <- 1
    sums <- init_sums(z_train$data, z[, t], n)
    detection_stat <- 0
    while ((detection_stat < threshold) & (t < n)) {
      t <- t + 1
      sums <- update_sumsC(sums, z[, t], m, t)
      log_liks <- mixture_log_liksC(sums, m, t, w, 1)
      detection_stat <- max(log_liks)
    }
    t
  }
  parallel::stopCluster(c1)
  arl_est <- n / mean(as.numeric(run_lengths < n))
  arl_est
}
