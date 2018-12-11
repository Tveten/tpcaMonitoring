library(MASS)
library(clusterGeneration)
library(Matrix)

gen_ind_data <- function(N, stream.length,
                         K = 0, kappa = 0, mu = 0, sigma2 = 1) {
  # Generates all relevant different data streams.
  #
  # Args:
  #   N: Dimension of the original data stream.
  #   stream.length:  Wanted length of the data stream.
  #   n.affected = 0: Number of affected streams.
  #   mu = 0: Size of the change in all entries in the n.affected-dimensional
  #           mean vector.
  #   sigma2 = 1: Variance of each indiviual stream post-change.
  #
  # Return:
  #   X: The specified (N x stream.length)-dimensional data stream

  #set.seed(1)
  length.affected <- stream.length - kappa
  X <- matrix(rnorm(N * stream.length), ncol = stream.length)
  if (K > 0) {
    X.affected <- t(mvrnorm(length.affected,
                            mu = rep(mu, K),
                            Sigma = diag(sigma2, nrow = K)))
    X[1:K, (kappa + 1):stream.length] <- X.affected
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

generate_cor_mat <- function(N, K0 = N, all.positive = FALSE) {
  # k : Sparsity level, number of correlated dimensions.

  # set.seed(652)
  Sigma <- genPositiveDefMat(K0, covMethod = 'onion', rangeVar = c(1, 1))$Sigma
  if (K0 != N) {
    identity.mat <- diag(rep(1, N - K0))
    zero.mat <- matrix(0, ncol = N - K0, nrow = K0)
    Sigma <- cbind(Sigma, zero.mat)
    Sigma <- rbind(Sigma, cbind(t(zero.mat), identity.mat))
  }
  if (all.positive) {
    Sigma <- sqrt(Sigma^2)
  }
  return(Sigma)
}

change_cor_mat <- function(R, affected.dims,
                               draw_rho = NULL, draw_sigma = NULL) {
  # At least one of the NULL-arguments must be supplied:
  #   functions draw_rho or draw_sigma
  #
  # Returns:
  #   Sigma2: The change covariance matrix.

  change_cor <- function(R.sub, draw_rho, n.affected) {
    affected.cor.dims <- sample(1:K0, min(n.affected, K0))
    ind <- combn(affected.cor.dims, 2)
    change.factor <- draw_rho(sum(1:length(ind)))
    R.changed <- R.sub
    for (i in 1:ncol(ind)) {
      R.changed[ind[1, i], ind[2, i]] <- change.factor[i] * R.changed[ind[1, i], ind[2, i]]
      R.changed[ind[2, i], ind[1, i]] <- change.factor[i] * R.changed[ind[2, i], ind[1, i]]
    }
    R.changed <- as.matrix(nearPD(R.changed,
                                  corr = TRUE,
                                  eig.tol = 1e-10,
                                  maxit = 10^3)$mat)
  }

  N <- ncol(R)
  K0 <- sum(R[, 1] != 0)
  K <- length(affected.dims)
  Sigma2 <- R

  if (all(is.null(c(draw_rho, draw_sigma))))
    stop('ERROR: Either a variance or a correlation change distribution must be specified')

  # Correlation change handling
  if (!is.null(draw_rho)) {
    R.sub <- R[1:K0, 1:K0, drop = FALSE]
    Sigma2[1:K0, 1:K0] <- change_cor(R.sub, draw_rho, K)
  }

  # Variance change handling
  if (!is.null(draw_sigma)) {
    sigma2.vec <- rep(1, N)
    sigma2.vec[affected.dims] <- draw_sigma(K)
    Sigma2 <- diag(sigma2.vec) %*% Sigma2 %*% diag(sigma2.vec)
  }

  return(Sigma2)
}

## TEST
# N <- 100
# Sigma1 <- GenerateCorrMatrix(N)
# Sigma2 <- ChangeCorrMatrix(Sigma1, 10, 1, 2, method = 'svd')
