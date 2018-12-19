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
  run_lengths <- rep(0, n_sim)
  for (b in 1:n_sim) {
    x_train <- t(MASS::mvrnorm(m, mu = mu_x, Sigma = Sigma_x))
    mu_hat <- rowMeans(x_train)
    sigma_hat <- rowSds(x_train)
    x_train <- (x_train - mu_hat) / sigma_hat

    cor_mat_hat <- 1 / (m - 1) * x_train%*% t(x_train)
    pca_obj <- tpca::pca(cor_mat_hat, axes = axes)
    V <- pca_obj$vectors
    lambda <- pca_obj$values

    z_train <- V %*% x_train / sqrt(lambda)
    mu_z <- 1 / sqrt(lambda) * V %*% ((mu_x - mu_hat) / sigma_hat)
    D <- diag(1 / sigma_hat)
    sigma2_z <- 1 / lambda * diag(V %*% (D %*% Sigma_x %*% D) %*% t(V))

    z <- t(MASS::mvrnorm(n, mu = mu_z, Sigma = diag(sigma2_z)))
    t <- 1
    sums <- init_sums(z_train, z[, t], n)
    detection_stat <- 0
    while ((detection_stat < threshold) & (t < n)) {
      t <- t + 1
      sums <- update_sums(sums, z[, t], m, t)
      log_liks <- mixture_log_liks(sums, m, t, w, 1)
      detection_stat <- max(log_liks, na.rm = TRUE)
    }
    run_lengths[b] <- t
  }
  arl_est <- n / mean(as.numeric(run_lengths < n))
  arl_est
}

#' Estimating the threshold for tpca changepoint detection
#'
#' Description
#'
#' \code{n} and \code{alpha} governs the false alarm rate by the relation
#'   P(T < n | H_0) <= \code{alpha}.
#' The corresponding average run length (\code{arl}) is approximately given by n / alpha.
#'
#' \code{rel_tol} and \code{thresh_alpha} governs the number of simulations used in
#' each step of the algorithm towards a more and more certain estimate.
#' At each step, the number of simulations is chosen so that
#' \code{[(1 - rel_tol)arl, (1 + rel_tol)arl]} approximately covers the true
#' average run length at confidence level \code{thresh_alpha}.
#' For example, when the last \code{rel_tol} is 0.025, it means that the final
#' estimated threshold corresponds to an average run length of approximately
#' \code{arl} +- 0.025 * \code{arl} at confidence level \code{thresh_alpha}.
#' The algorithm should start with a large relative error tolerance,
#' and then narrow it down for the quickest convergence.
#'
#' @param x The d x m training data matrix, where d is the dimension of the data
#'   and m the number of training samples
#' @param axes Indices of the principal axes to be used in simulations.
#' @param n The length of the segment to monitor for false alarms. See details.
#' @param alpha Probability of type I error (false alarm) within the time window . See details.
#' @param w The window size (an integer).
#' @param rel_tol A vector with the sequence of relative error tolerances
#'   allowed at each step towards convergence in the algorithm. See details.
#' @param thresh_alpha 1 - \code{thresh_alpha} is the confidence level. See details.
#' @param init_thresh Use if a custom initial value for the threshold is wanted.
#' @param learning_coef The learning rate is defined as \code{rel_tol} /
#' \code{learning_coef}. The default \code{learning_coef} is 3, which has
#' been found a good choice after a lot of experimenting.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{\code{threshold}}{The final estimate }
#' }
#' \code{est_tpca_threshold} also creates a .txt file with each line showing
#' summaries of each step in the estimation procedure. The last line corresponds
#' to the final estimate.
#'
#' @export

est_tpca_threshold <- function(x, axes, n, alpha,
                               w             = 200,
                               rel_tol       = c(0.2, 0.1, 0.05, 0.025),
                               thresh_alpha  = 0.05,
                               init_thresh   = NULL,
                               learning_coef = NULL) {

  set_log_name <- function() {
    files_in_wdir <- list.files()
    root_name <- paste0('est_tpca_threshold_log_',
                        'm', as.character(m),
                        'd', as.character(d),
                        'r', as.character(r))
    n_logs <- sum(grepl(root_name, files_in_wdir))
    if (n_logs > 0)
      root_name <- paste0(root_name, '(', as.character(n_logs + 1), ')')
    paste0(root_name, '.txt')
  }

  set_results_name <- function() {
    files_in_wdir <- list.files()
    root_name <- paste0('tpca_threshold_results_',
                        'm', as.character(m),
                        'd', as.character(d),
                        'r', as.character(r))
    n_logs <- sum(grepl(root_name, files_in_wdir))
    if (n_logs > 0)
      root_name <- paste0(root_name, '(', as.character(n_logs + 1), ')')
    paste0(root_name, '.txt')
  }

  log_next_stage <- function(rel_tol) {
    write(sprintf('=========================================================='),
          file = log_file, append = TRUE)
    write(sprintf('Relative error tolerance: %.3f',
                  rel_tol),
          file = log_file, append = TRUE)
    write(sprintf(' '),
          file = log_file, append = TRUE)
  }

  log_current <- function(n_sim, threshold) {
    write(sprintf('Currently: threshold = %.4f with %.0f simulations',
                  threshold, n_sim),
          file = log_file, append = TRUE)
  }

  log_results <- function(arl_est, threshold) {
    write(sprintf('Result:    ARL       = %.0f', arl_est), file = log_file, append = TRUE)
  }

  init_log_file <- function() {
    write(sprintf('Data dimension             = %.0f', d), file = log_file, append = TRUE)
    write(sprintf('Reduced dimension          = %.0f', r), file = log_file, append = TRUE)
    write(paste0('Axes                       = ', paste(axes, collapse = ' ')), file = log_file, append = TRUE)
    write(sprintf('Number of training samples = %.0f', m), file = log_file, append = TRUE)
    write(sprintf('Window length              = %.0f', w), file = log_file, append = TRUE)
  }

  init_results_file <- function() {
    header_items <- c('data_dim', 'reduced_dim', 'n_training',
                      'window_size', 'n_sim', 'threshold',
                      'arl_est', 'arl_lower', 'arl_upper')
    write(paste(header_items, collapse = ' '),
          file = results_file, append = TRUE, ncolumns = length(header_items))
  }

  store_results <- function(n_sim, threshold, arl_est, arl_conf_int) {
    stored_values <- c(d, r, m, w, n_sim, threshold, arl_est, arl_conf_int)
    write(stored_values, file = results_file,
          append = TRUE, ncolumns = length(stored_values))

  }

  geom_conf_int <- function(mean_est, n, thresh_alpha) {
    p_est <- 1 / mean_est
    sd_est <- sqrt((1 - p_est) / p_est^2)
    quartile <- qnorm((1 - thresh_alpha / 2))
    lower <- mean_est - quartile * sd_est / sqrt(n)
    upper <- mean_est + quartile * sd_est / sqrt(n)
    conf_int <- c(lower, upper)
    return(conf_int)
  }

  get_n_sim <- function(thresh_alpha, rel_tol) {
    z <- qnorm(1 - thresh_alpha / 2)
    n_sim <- (z / rel_tol)^2
    return(round(n_sim))
  }

  init_threshold <- function(r) {
    if (is.null(init_thresh)) return(9 + 1.2 * r)
    else return(init_thresh)
  }

  update_threshold <- function(threshold, rel_tol, arl, arl_est) {
    if (is.null(learning_coef)) learning_rate <- rel_tol / 3
    else learning_rate <- rel_tol / learning_coef
    threshold_change <- (arl - arl_est) / arl * learning_rate
    threshold * (1 + threshold_change)
  }

  d <- nrow(x)
  m <- ncol(x)
  mu_x <- rowMeans(x)
  Sigma_x <- 1 / (m - 1) * x %*% t(x)
  attr(Sigma_x, 'n_obs') <- m

  r <- length(axes)
  n_sim <- get_n_sim(thresh_alpha, rel_tol)
  threshold <- init_threshold(r)
  arl <- n / alpha

  log_file <- set_log_name()
  results_file <- set_results_name()
  init_log_file()
  init_results_file()

  # When rel_tol is low, the algorithm might jump back and forth over the true
  # value. Thus, it should be stopped at some point. Too many runs on the final
  # stage is not worth it.
  max_retries <- 5
  n_retries <- 1
  for (k in 1:length(n_sim)) {
    arl_conf_int <- c(0, 0)
    log_next_stage(rel_tol[k])
    while (!is_in_interval(arl, arl_conf_int) && (n_retries <= max_retries)) {
      log_current(n_sim[k], threshold)
      arl_est <- tpca_arl(threshold, mu_x, Sigma_x, axes, n, w, n_sim[k])
      arl_conf_int <- geom_conf_int(arl_est, n_sim[k], thresh_alpha)
      log_results(arl_est, threshold)
      store_results(n_sim[k], threshold, arl_est, arl_conf_int)

      if (!is_in_interval(arl, arl_conf_int))
        threshold <- update_threshold(threshold, rel_tol[k], arl, arl_est)
      if (k == length(n_sim)) n_retries <- n_retries + 1
    }
  }
  threshold
}
