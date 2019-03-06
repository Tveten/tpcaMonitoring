gen_train_dense <- function(d, m, seed, return_all = FALSE) {
  set.seed(as.numeric(seed))
  mu <- rep(0, d)
  Sigma <- tpca::rcor_mat(d)
  x <- gen_norm_data(m, mu, Sigma)
  if (return_all) return(list('x' = x, 'mu' = mu, 'Sigma' = Sigma,
                              'nr' = as.numeric(substr_right(seed, 2))))
  else return(x)
}

gen_train_halfsparse <- function(d, m, return_all = FALSE) {
  set.seed(102)
  mu <- rep(0, d)
  Sigma <- tpca::rcor_mat(d, k0 = round(d / 2))

  # Separate seed so Sigma can be kept same, but data varied for verification.
  set.seed(202)
  x <- gen_norm_data(m, mu, Sigma)
  if (return_all) return(list('x' = x, 'mu' = mu, 'Sigma' = Sigma,
                              'type' = 'halfsparse'))
  else return(x)
}

gen_train_sparse <- function(d, m, return_all = FALSE) {
  set.seed(103)
  mu <- rep(0, d)
  Sigma <- tpca::rcor_mat(d, k0 = round(d / 10))

  # Separate seed so Sigma can be kept same, but data varied for verification.
  set.seed(203)
  x <- gen_norm_data(m, mu, Sigma)
  if (return_all) return(list('x' = x, 'mu' = mu, 'Sigma' = Sigma,
                              'type' = 'sparse'))
  else return(x)
}

get_train_dense <- function(return_all = FALSE) {
  gen_train_dense(100, 300, return_all)
}

get_train_halfsparse <- function(return_all = FALSE) {
  gen_train_halfsparse(100, 300, return_all)
}

get_train_sparse <- function(return_all = FALSE) {
  gen_train_sparse(100, 300, return_all)
}

get_training_sets <- function(n_sets, d, m, return_all = TRUE) {
  seed_seq <- paste0(300, sprintf('%.2d', 1:n_sets))
  lapply(seed_seq, function(s) gen_train_dense(d, m, s, return_all))
}
