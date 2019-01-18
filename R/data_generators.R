generate_cor_data <- function(d, n, mu0, Sigma0,
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

gen_norm_data <- function(n, mu, Sigma) {
  t(MASS::mvrnorm(n, mu = mu, Sigma = Sigma))
}
