find_tpca_thresholds <- function(x_train, cov_mat_type) {
  unique_vectors <- function(list_of_vectors) {
    list_of_vectors[!duplicated(lapply(list_of_vectors, sort))]
  }

  log_tpca <- function(change_distr, cutoff, axes) {
    out <- c(change_distr, cutoff, axes)
    tpca_log_file <- paste0(cov_mat_type, '_tpca_axes.txt')
    write(out, file = tpca_log_file, ncolumns = length(out), append = TRUE)
  }

  which_axes <- function(prop_max, keep_prop, max_axes) {
    order_axes <- order(prop_max, decreasing = TRUE)
    cum_prop <- cumsum(prop_max[order_axes])
    n_keep <- min(sum(cum_prop < keep_prop) + 1, max_axes)
    order_axes[1:n_keep]
  }

  get_axes_list <- function() {
    axes_pca_max <- lapply(r, function(j) 1:j)
    axes_pca_min <- lapply(r, function(j) (d - j + 1):d)
    axes_list <- list(axes_pca_max, axes_pca_min)
    change_distrs <- c('full_uniform', 'mean_only', 'sd_only', 'cor_only')
    for (i in seq_along(change_distrs)) {
      prop_max <- tpca::tpca(cov_mat,
                             change_distr = change_distrs[i],
                             cutoff = cutoff[1],
                             max_axes = max_axes)$prop_axes_max
      which_axes_list <- list()
      for (j in seq_along(cutoff)) {
        which_axes_list[[j]] <- which_axes(prop_max, cutoff[j], max_axes)
        log_tpca(change_distrs[i], cutoff[j], which_axes_list[[j]])
      }
      axes_list[[2 + i]] <- which_axes_list
    }
    axes_list <- unlist(axes_list, recursive = FALSE)
    unique_axes <- unique_vectors(axes_list)
    unique_axes
  }

  m <- ncol(x_train)
  d <- nrow(x_train)
  cov_mat <- 1 / (m - 1) * x_train %*% t(x_train)
  r <- c(1, 2, 3, 5, 10, 20)
  cutoff <- c(0.8, 0.9, 0.95, 0.99, 0.995, 0.999)
  max_axes <- round(d / 5)
  axes_list <- get_axes_list()

  n <- 1000
  alpha <- 0.01

  for (i in seq_along(axes_list))
    threshold_finder(x_train, 'tpca', n, alpha, axes = axes_list[[i]],
                     file_id = paste0('_', cov_mat_type))
}

find_mixture_thresholds <- function(x_train) {
  m <- ncol(x_train)
  d <- nrow(x_train)
  # p0 <- c(0.03, 0.1, 0.3, 1)
  p0 <- c(0.03, 0.1, 0.3)
  init_thresh <- c(30, 50, 80)

  n <- 1000
  alpha <- 0.01

  for (i in seq_along(p0))
    threshold_finder(x_train, 'mixture', n, alpha, p0 = p0[i],
                     init_thresh   = init_thresh[i],
                     learning_coef = 10)
}

run_threshold_finder_dense <- function() {
  find_tpca_thresholds(get_train_dense(), 'dense')
}

run_threshold_finder_halfsparse <- function() {
  find_tpca_thresholds(get_train_halfsparse(), 'halfsparse')
}

run_threshold_finder_sparse <- function() {
  find_tpca_thresholds(get_train_sparse(), 'sparse')
}

run_threshold_finder_mixture <- function() {
  # The dimension of all the training sets are the same,
  # so any of them can be input below with the same result.
  find_mixture_thresholds(get_train_dense())
}

run_all_threshold_finders <- function() {
  run_threshold_finder_dense()
  run_threshold_finder_halfsparse()
  run_threshold_finder_sparse()
}

run_selected_threshold_finders <- function() {
  run_threshold_finder_mixture()
  run_threshold_finder_sparse()
}
