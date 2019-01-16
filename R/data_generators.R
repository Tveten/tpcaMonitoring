  }
  return(X)
}

generate_cor_data <- function(N, stream.length, K0 = N,
                              K = 0, kappa = 0,
                              mu = 0, sigma = 1, rho.scale = 1,
                              return.cor.mat = FALSE) {

  mu0 <- rep(0, N)
  Sigma0 <- generate_cor_mat(N, K0 = K0)

  if (K == 0)
    X <- t(mvrnorm(stream.length, mu = mu0, Sigma = Sigma0))
  else if (K > 0 && K <= N) {
    affected.dims <- sample(1:N, K)

    mu1 <- rep(0, N)
    Sigma1 <- Sigma0
    if (mu != 0) mu1[affected.dims] <- mu

    if (sigma != 1 || rho.scale != 1) {
      draw_sigma <- NULL
      draw_rho <- NULL
      if (sigma != 1) draw_sigma <- function(n) rep(sigma, n)
      if (rho.scale != 1) draw_rho <- function(n) rep(rho.scale, n)
      Sigma1 <- change_cor_mat(Sigma0, affected.dims, draw_rho, draw_sigma)
    }

    if (kappa == 0)
      X <- t(mvrnorm(stream.length, mu = mu1, Sigma = Sigma1))
    else {
      X0 <- t(mvrnorm(kappa, mu = mu0, Sigma = Sigma0))
      X1 <- t(mvrnorm(stream.length - kappa, mu = mu1, Sigma = Sigma1))
      X <- cbind(X0, X1)
    }
  } else stop('Invalid number of affected streams K. No data generated.')

  if (return.cor.mat)
    return(list('X' = X, 'cor.mat' = Sigma0))
  else
    return(X)
}

gen_norm_data <- function(n, mu, Sigma) {
  t(MASS::mvrnorm(n, mu = mu, Sigma = Sigma))
}
