boot_z_train <- function(m, mu_x, Sigma_x, axes) {
  d <- length(mu_x)
  r <- length(axes)
  x_train <- gen_norm_data(m, mu_x, Sigma_x)
  mu_hat <- rowMeans(x_train)
  sigma_hat <- rowSds(x_train)
  x_train <- (x_train - mu_hat) / sigma_hat

  cor_mat_hat <- 1 / (m - 1) * x_train %*% t(x_train)
  pca_obj <- tpca::pca(cor_mat_hat, axes = axes)
  V <- pca_obj$vectors
  lambda <- pca_obj$values

  z_train <- V %*% x_train / sqrt(lambda)
  mu_z <- 1 / sqrt(lambda) * V %*% ((mu_x - mu_hat) / sigma_hat)
  D <- diag(1 / sigma_hat, nrow = d)
  sigma2_z <- 1 / lambda * diag(V %*% (D %*% Sigma_x %*% D) %*% t(V))
  return(list('data'   = z_train,
              'mu'     = as.vector(mu_z),
              'sigma2' = as.vector(sigma2_z)))
}

boot_z_train_np <- function(x, axes) {
  m <- ncol(x)
  x_train <- x[, sample(1:m, m, replace = TRUE)]
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
  return(list('data'   = z_train,
              'mu'     = mu_z,
              'sigma2' = sigma2_z))
}

est_arl <- function(run_lengths, n, n_sim) {
  arl_est <- n / mean(as.numeric(run_lengths < n))
  if (is.infinite(arl_est)) arl_est <- 5 * n / (1 / n_sim)
  round(arl_est)
}

conf_int <- function(x, alpha = 0.05) {
  ci <- quantile(x, probs = c(alpha / 2, (1 - alpha) / 2))
}

geom_conf_int <- function(mean_est, n, thresh_alpha) {
  p_est <- 1 / mean_est
  sd_est <- sqrt((1 - p_est) / p_est^2)
  percentile <- qnorm((1 - thresh_alpha / 2))
  lower <- mean_est - percentile * sd_est / sqrt(n)
  upper <- mean_est + percentile * sd_est / sqrt(n)
  conf_int <- c(floor(lower), ceiling(upper))
  return(conf_int)
}
