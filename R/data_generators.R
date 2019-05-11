gen_norm_data <- function(n, mu, Sigma) {
  t(MASS::mvrnorm(n, mu = mu, Sigma = Sigma))
}

gen_changed_data <- function(d, n, mu0, Sigma0,
                             kappa = 0, p = 0,
                             mu = NULL, sigma = NULL, rho_scale = NULL) {

  change_mu <- function(mu0, mu, D) {
    if (!is.null(mu)) {
      mu1 <- mu0
      mu1[D] <- mu
      return(mu1)
    } else return(mu0)
  }

  change_Sigma <- function(Sigma0, sigma, rho_scale, D) {
    if (!is.null(sigma) || !is.null(rho_scale)) {
      draw_sigma <- NULL
      draw_rho <- NULL
      if (!is.null(sigma)) draw_sigma <- function(n) rep(sigma, n)
      if (!is.null(rho_scale)) {
        draw_rho <- function(n) rep(rho_scale, n)
      }
      Sigma1 <- tpca::change_cor_mat(Sigma0, D,
                                     draw_cor = draw_rho,
                                     draw_sd  = draw_sigma)
      return(Sigma1)
    } else return(Sigma0)

  }

  assert_cor_mat_attribute(Sigma0)

  mu0 <- rep(0, d)
  k <- round(p * d)

  if (k == 0) {
    x <- gen_norm_data(n, mu0, Sigma0)
    warning(paste0('p = ', p, 'is so small that no dimensions are changed.'))
  } else if (k > 0 && k <= d) {
    affected_dims <- sample(1:d, k)
    mu1 <- change_mu(mu0, mu, affected_dims)
    Sigma1 <- change_Sigma(Sigma0, sigma, rho_scale, affected_dims)

    if (kappa == 0)
      x <- gen_norm_data(n, mu1, Sigma1)
    else {
      x0 <- gen_norm_data(kappa, mu0, Sigma0)
      x1 <- gen_norm_data(n - kappa, mu1, Sigma1)
      x <- cbind(x0, x1)
    }
  } else stop('Invalid proportion of affected streams p. No data generated.')

  return(x)
}

#' Generate multivariate normal training data.
#'
#' @param d The data dimension
#' @param m The number of training samples.
#' @param seed The seed to be used.
#' @param alphad The value of alpha_d in \code{\link{tpca::rcor_mat}}. A high
#' value means low correlations and vice versa.
#' @param return_all A logical value. If FALSE: Only returns the data matrix x.
#' If TRUE: Returns a list with 'x', 'mu', 'Sigma' and 'nr', where
#' 'mu' is the mean vector used, 'Sigma' is the covariance matrix and 'nr' indicates
#' which seed that was used.
#'
#' @return See the argument return_all.
#'
#' @export
gen_train <- function(d, m, seed, alphad = 1, return_all = FALSE) {
  set.seed(as.numeric(seed))
  mu <- rep(0, d)
  Sigma <- tpca::rcor_mat(d, alphad = alphad)
  x <- gen_norm_data(m, mu, Sigma)
  if (return_all) return(list('x' = x, 'mu' = mu, 'Sigma' = Sigma,
                              'nr' = as.numeric(substr_right(seed, 2))))
  else return(x)
}

#' Generate multivariate normal training data.
#'
#' @param n_sets The number of training sets to be generated.
#' @param d The data dimension.
#' @param m The number of training samples.
#' @param return_all A logical value. If FALSE: Only returns the data matrix x.
#' If TRUE: Returns a list with 'x', 'mu', 'Sigma' and 'nr', where
#' 'mu' is the mean vector used, 'Sigma' is the covariance matrix and 'nr' indicates
#' which seed that was used.
#'
#' @return A list of training sets from \code{\link{gen_train}}.
#' The argument return_all determines the type of each list element.
#'
#' @export
get_training_sets <- function(n_sets, d, m, return_all = TRUE) {
  seed_seq <- paste0(300, sprintf('%.2d', 1:n_sets))
  n_large_cor <- round(n_sets / 2)
  n_small_cor <- n_sets - n_large_cor
  alphad_seq <- c(seq(0.05, 0.95, length.out = n_large_cor),
                  seq(1, 50, length.out = n_small_cor))
  lapply(1:n_sets, function(i) gen_train(d, m, seed_seq[i], alphad_seq[i], return_all))
}
