rowVars <- function(x) {
  d <- ncol(x)
  x_mean <- rowMeans(x)
  d / (d - 1) * rowMeans((x - x_mean)^2)
}

rowSds <- function(x) {
  sqrt(rowVars(x))
}

setup_parallel <- function() {
  n_cores <- parallel::detectCores()
  c <- parallel::makeCluster(n_cores - 1, outfile = '', type = 'PSOCK')
  doParallel::registerDoParallel(c)
  c
}

stop_parallel <- function(c) {
  parallel::stopCluster(c)
}

is_in_interval <- function(x, interval) {
  x >= interval[1] & x <= interval[2]
}

is_equal_vectors <- function(x, y) {
  if (length(x) != length(y)) return(FALSE)
  is_equal <- function(u, v) isTRUE(all.equal(u, v))
  all(unlist(Map(is_equal, x, y)))
}

assert_cor_mat_attribute <- function(Sigma) {
  if (is.null(attr(Sigma, 'which_dims_cor')))
    stop('Sigma must contain the attribute "which_dims_cor", as when generated from tpca::r_cor_mat or tpca::r_cov_mat')
}
