rowVars <- function(x) {
  d <- ncol(x)
  x_mean <- rowMeans(x)
  d / (d - 1) * rowMeans((x - x_mean)^2)
}

rowSds <- function(x) {
  sqrt(rowVar(x))
}

get_training_data <- function(x) {
  if (is.list(x)) {
    x_true_mean <- rep(0, x$d)
    x_true_cor_mat <- tpca::rcor_mat(x$d)
    x_train <- t(mvrnorm(x$m, x_mean_true, x_cor_mat_true))
  } else if (is.matrix(x)) {
    x_train <- x
  } else {
    stop('x must be either a matrix with training data, or a list(d = ., m = .) where d is the dimension of the data and m the number of training samples.')
  }
}
