pca_arl <- function(threshold, d, r, m, n_sim, ) {
  x_cor_mat <- rcor_mat(d)
  x_mean <- rep(0, d)
  x_train <- t(mvrnorm(m, x_mean, x_cor_mat))
  
  cor_mat_est <- standardize_cov_mat(1 / (m - 1) * x_train %*% t(x_train))
  V <- pca(cor_mat_est, axes = (d - r + 1):d)$vectors
  y_train
}