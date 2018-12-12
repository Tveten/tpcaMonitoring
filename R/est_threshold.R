#' Estimates average run length (ARL) for monitoring by tpca
#'
#' Description
#'
#' Details
#'
#' @param x The d x m training data matrix, where d is the dimension of the data and m the number of training samples, or a list(d = ., m = .) describing the dimensions of normal training data to be drawn.
#' @param threshold A numeric specifying the threshold value for when a change
#'   is declared.
#' @param d The dimension of the original data stream.
#' @param r The dimension of the reduced data stream, i.e.
#'   the number of principal axes to project onto.
#' @param m The size of the training set.
#' @param n The number of observations to monitor for an estimate of the ARL.
#'   See details.
#' @param w The window size. Number of recent time-points to consider for a
#'   change.
#' @param n_sim The number of simulations to base the estimate on.
#'
#' @return An estimate of the average run length.
#'
#' @export

tpca_arl <- function(x, threshold, axes, n, n_sim) {
  x_train <- get_training_data(x)
  x_mean <- rowMeans(x_train)
  x_sd <- rowSds(x_train)
  x_train <- (x_train - x_mean) / x_sd
  m <- ncol(x_train)
  d <- nrow(x_train)
  r <- length(axes)

  x_cor_mat <- 1 / (m - 1) * x_train %*% t(x_train)
  V <- tpca::pca(x_cor_mat, axes = axes)$vectors
  y_train <- V %*% x_train
}

find_tpca_threshold <- function(x, axes, arl) {

}
