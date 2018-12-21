gen_train_dense <- function(d, m) {
  set.seed(101)
  mu <- rep(0, d)
  Sigma <- tpca::rcor_mat(d)
  print(Sigma)

  # Separate seed so Sigma can be kept same, but data varied for verification.
  set.seed(201)
  gen_norm_data(m, mu, Sigma)
}

gen_train_halfsparse <- function(d, m) {
  set.seed(102)
  mu <- rep(0, d)
  Sigma <- tpca::rcor_mat(d, k0 = round(d / 2))

  # Separate seed so Sigma can be kept same, but data varied for verification.
  set.seed(202)
  gen_norm_data(m, mu, Sigma)
}

gen_train_sparse <- function(d, m) {
  set.seed(103)
  mu <- rep(0, d)
  Sigma <- tpca::rcor_mat(d, k0 = round(d / 10))

  # Separate seed so Sigma can be kept same, but data varied for verification.
  set.seed(203)
  gen_norm_data(m, mu, Sigma)
}

get_train_dense <- function() {
  gen_train_dense(100, 300)
}

get_train_halfsparse <- function() {
  gen_train_halfsparse(100, 300)
}

get_train_sparse <- function() {
  gen_train_sparse(100, 300)
}
