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
  mu_x <- rep(0, d)
  Sigma_x <- diag(1, nrow = d)
  x_train <- gen_norm_data(m, mu_x, Sigma_x)

  comp_cluster <- setup_parallel()
  `%dopar%` <- foreach::`%dopar%`
  run_lengths <- foreach::foreach(b = 1:n_sim, .combine = 'c') %dopar% {
    x <- gen_norm_data(n, mu_x, Sigma_x)
    mixture_rl(threshold, x_train, x, p0, kappa = Inf, w)
  }
  stop_parallel(comp_cluster)
  est_arl(run_lengths, n, n_sim)
}

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

  comp_cluster <- setup_parallel()
  `%dopar%` <- foreach::`%dopar%`
  run_lengths <- foreach::foreach(b = 1:n_sim, .combine = 'c') %dopar% {
    z_train <- boot_z_train(m, mu_x, Sigma_x, axes)
    z <- gen_norm_data(n, mu = z_train$mu, Sigma = diag(z_train$sigma2, nrow = r))
    mixture_rl(threshold, z_train$data, z, 1, kappa = Inf, w)
  }
  stop_parallel(comp_cluster)
  est_arl(run_lengths, n, n_sim)
}
