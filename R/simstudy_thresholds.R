find_tpca_thresholds <- function(train_obj, n, alpha, rel_tol, r, cutoff) {
  unique_vectors <- function(list_of_vectors) {
    list_of_vectors[!duplicated(lapply(list_of_vectors, sort))]
  }

  log_tpca <- function(change_distr, cutoff, axes) {
    out <- c(change_distr, cutoff, axes)
    dir <- './thresholds/axes/'
    tpca_log_file <- paste0(dir, 'tpca_axes_', cov_mat_nr, '_m', m, 'd', d, '.txt')
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
    change_distrs <- c('mean_only', 'sd_only', 'cor_only')
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

  write_global_log(paste0('Finding tpca thresholds for cov_mat_nr ', train_obj$nr, '.'))

  x_train <- train_obj$x
  cov_mat_nr <- train_obj$nr
  m <- ncol(x_train)
  d <- nrow(x_train)
  cov_mat <- 1 / (m - 1) * x_train %*% t(x_train)
  max_axes <- round(d / 5)
  axes_list <- get_axes_list()

  for (i in seq_along(axes_list))
    threshold_finder(x_train, 'tpca', n, alpha, axes = axes_list[[i]],
                     rel_tol = rel_tol,
                     file_id = paste0(cov_mat_nr, '_'))
}

find_mixture_thresholds <- function(train_obj, n, alpha, rel_tol, p0) {
  write_global_log('Finding mixture thresholds.')
  x_train <- train_obj$x
  m <- ncol(x_train)
  d <- nrow(x_train)
  # init_thresh <- c(30, 50, 80)

  for (i in seq_along(p0))
    threshold_finder(x_train, 'mixture', n, alpha, p0 = p0[i], rel_tol = rel_tol)
    # threshold_finder(x_train, 'mixture', n, alpha, p0 = p0[i],
    #                  init_thresh   = init_thresh[i],
    #                  learning_coef = 10)
}

find_all_thresholds <- function(training_sets, n, alpha, rel_tol, p0, r, cutoff) {
  find_mixture_thresholds(training_sets[[1]], n, alpha, rel_tol, p0)
  lapply(training_sets, find_tpca_thresholds, n = n, alpha = alpha,
         rel_tol = rel_tol, r = r, cutoff = cutoff)
}

extract_final_thresholds <- function(d, m, n, alpha,
                                     mon_type = c('mixture', 'tpca')) {
  extract_p0 <- function(p0_id) {
      p0 <- strsplit(p0_id, '-')[[1]][2]
      p0 <- paste0('0.', p0)
      10 * as.numeric(p0)
  }

  get_final_threshold <- function(curr_file, id) {
    threshold_df <- read.table(paste0(path2, curr_file), header = TRUE, sep = ' ')
    split_id <- strsplit(id, 'n|alpha|m|d')
    n <- as.numeric(split_id[[1]][2])
    alpha <- as.numeric(paste0('0.', split_id[[1]][3]))
    target_arl <- n / alpha
    threshold_sub_df <- subset(threshold_df, n_sim == max(threshold_df$n_sim))
    final_threshold_ind <- which.min(abs(threshold_sub_df$arl_est - target_arl))
    threshold_sub_df$threshold[final_threshold_ind]
  }

  alpha <- strsplit(as.character(alpha), '[.]')[[1]][2]
  path1 <- './thresholds/'
  path2 <- paste0(path1, 'threshold_results/')

  if (any(mon_type == 'mixture')) {
    type <- 'mixture'
    files <- list_files_matching(path2, type, 'results', paste0('d', d),
                                 paste0('m', m), paste0('n', n), paste0('alpha', alpha))
    split_files <- strsplit(files, '_|p0|ax|[.]')
    id_ind <- 4
    simulation_id <- vapply(split_files, function(x) x[id_ind], character(1))
    unique_id <- unique(simulation_id)

    final_files <- paste0(path1, type, '_thresholds_', unique_id, '_FINAL.txt')
    for (i in seq_along(final_files)) {
      which_files <- simulation_id == unique_id[i]
      curr_files <- files[which_files]
      curr_split_files <- split_files[which_files]

      p0_id <- vapply(curr_split_files, function(split_file) split_file[id_ind + 1], character(1))
      p0 <- vapply(p0_id, extract_p0, numeric(1))
      write(c('p0', 'threshold'), final_files[i], ncolumns = 2)
      for (j in seq_along(curr_files)) {
        final_threshold <- get_final_threshold(curr_files[j], unique_id[i])
        write(c(p0[j], final_threshold), final_files[i], ncolumns = 2, append = TRUE)
      }
    }
  }

  if (any(mon_type == 'tpca')) {
    type <- 'tpca'
    files <- list_files_matching(path2, type, 'results', paste0('d', d),
                                 paste0('m', m), paste0('n', n), paste0('alpha', alpha))
    split_files <- strsplit(files, '_|p0|ax|[.]')
    id_ind <- c(3, 5)
    simulation_id <- vapply(split_files, function(x) {
        paste(x[id_ind], collapse = '_')
      }, character(1))
    unique_id <- unique(simulation_id)

    final_files <- paste0(path1, type, '_thresholds_', unique_id, '_FINAL.txt')
    for (i in seq_along(unique_id)) {
      which_files <- simulation_id == unique_id[i]
      curr_files <- files[which_files]
      curr_split_files <- split_files[which_files]

      axes_id <- vapply(curr_split_files, function(split_file) split_file[6], character(1))
      write(c('axes', 'threshold'), final_files[i], ncolumns = 2)
      for (j in seq_along(curr_files)) {
        final_threshold <- get_final_threshold(curr_files[j], unique_id[i])
        write(c(axes_id[j], final_threshold), final_files[i], ncolumns = 2, append = TRUE)
      }
    }
  }
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
