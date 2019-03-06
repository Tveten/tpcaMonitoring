#' Estimates the edd for a given change for all methods.
#'
#' Description
#'
#' Details.
#'
#' @param kappa The d x m training data matrix, where d is the dimension of the data
#'   and m the number of training samples
#'
#' @return A list with the following components:
#' \describe{
#'   \item{\code{threshold}}{The final estimate }
#' }
#' \code{tpca_threshold} also creates a .txt file with each line showing
#' summaries of each step in the estimation procedure. The last line corresponds
#' to the final estimate.
#'
#' @export
est_edd_all_methods <- function(train_obj, n_max, n_sim, edd_file, kappa, p, w,
                                p0s, r_pca, threshold_settings,
                                mu = NULL, sigma = NULL, rho_scale = NULL) {

  get_change_type <- function() {
    n_not_null <- 0
    if (!is.null(mu)) {
      n_not_null <- n_not_null + 1
      change_type <- 'mean'
    }
    if(!is.null(sigma)) {
      n_not_null <- n_not_null + 1
      change_type <- 'sd'
    }
    if (!is.null(rho_scale)) {
      n_not_null <- n_not_null + 1
      change_type <- 'cor'
    }
    if (n_not_null > 1) stop('Only one of mu, sigma and rho_scale can be non-NULL at a time')

    change_type
  }

  set_change_param <- function(change_type) {
    if (change_type == 'mean') return(mu)
    if (change_type == 'sd') return(sigma)
    if (change_type == 'cor') return(rho_scale)
  }

  get_tpca_train_info <- function(cov_mat_nr, change_type) {
    dir <- './thresholds/axes/'
    axes_location <- paste0(dir, 'tpca_axes_', cov_mat_nr, '_m', m, 'd', d, '.txt')
    all_lines_vec <- readLines(axes_location)
    all_lines_split <- strsplit(all_lines_vec, ' ')
    axes_list <- list()
    cutoffs <- list()
    for (i in seq_along(all_lines_split)) {
      line <- all_lines_split[[i]]
      if (line[1] == paste0(change_type, '_only')) {
        axes_list[[length(axes_list) + 1]] <- as.numeric(line[3:length(line)])
        cutoffs[[length(cutoffs) + 1]] <- as.numeric(line[2])
      }
    }
    cutoffs <- unlist(cutoffs)
    names(axes_list) <- as.character(cutoffs)
    list('axes' = axes_list, 'cutoffs' = cutoffs)
  }

  get_mixture_threshold <- function(m, d, cov_mat_nr, p0) {
    dir <- './thresholds/'
    alpha_str <- strsplit(as.character(threshold_settings$alpha), '[.]')[[1]][2]
    threshold_file <- paste0(dir, 'mixture_thresholds_n', threshold_settings$n,
                             'alpha', alpha_str, 'm', m, 'd', d, '_FINAL.txt')
    threshold_df <- read.table(threshold_file, header = TRUE, sep = ' ')
    threshold <- threshold_df$threshold[threshold_df$p0 == p0]
    if (length(threshold) == 0)
      stop(paste0('Threshold for m = ', m, ', d = ', d, ', p0 = ', p0, ' does not exist.'))
    if (length(threshold) > 1) {
      threshold <- threshold[1]
      warning('More than one threshold returned from file. Picking the first one.')
    }
    threshold
  }

  get_tpca_threshold <- function(m, d, cov_mat_nr, axes) {
    dir <- './thresholds/'
    alpha_str <- strsplit(as.character(threshold_settings$alpha), '[.]')[[1]][2]
    threshold_file <- paste0(dir, 'tpca_thresholds_', cov_mat_nr,
                             '_n', threshold_settings$n, 'alpha', alpha_str,
                             'm', m, 'd', d, '_FINAL.txt')
    threshold_df <- read.table(threshold_file, header = TRUE, sep = ' ')
    threshold_df$axes <- as.character(threshold_df$axes)
    axes_in_df <- lapply(strsplit(threshold_df$axes, '-'), as.numeric)
    axes_in_df <- lapply(axes_in_df, sort)
    axes <- sort(axes)
    ind <- which(vapply(axes_in_df, is_equal_vectors, logical(1), y = axes))
    threshold <- threshold_df$threshold[ind]
    if (length(threshold) == 0)
      stop(paste0('Threshold for m = ', m, ', d = ', d, ', axes = (',
                  paste(axes, collapse = ', '), ') does not exist.'))
    if (length(threshold) > 1) {
      threshold <- threshold[1]
      warning('More than one threshold returned from file. Picking the first one.')
    }
    threshold
  }

  get_threshold <- function(m, d, cov_mat_nr, p0 = NULL, axes = NULL) {
    if (!is.null(p0)) return(get_mixture_threshold(m, d, cov_mat_nr, p0))
    if (!is.null(axes)) return(get_tpca_threshold(m, d, cov_mat_nr, axes))
  }

  get_log_name <- function(cov_mat_nr, change_type) {
    dir <- './results/'
    log_file <- list_files_matching(dir, 'log', cov_mat_nr, change_type)
    paste0(dir, log_file)
  }

  log_new_change <- function(log_file) {
    write('===================================', file = log_file, append = TRUE)
  }

  log_current_stage <- function(log_file, change_type, l) {
    time_used <- round(proc.time()[3]/60, 2)
    log_str <- paste0('Sim nr. ', l, ' with p = ', p, ' and ',
                      change_type, ' = ', change_param, '. ',
                      time_used, ' min used.')
    write(log_str, file = log_file, append = TRUE)
  }

  get_rl_mat <- function(all_run_lengths, method) {
    rl_mat <- lapply(all_run_lengths, function(x) x[[method]])
    do.call('rbind', rl_mat)
  }

  store_edd_results <- function(rl_mat, method, edd_file, method_param) {
    for (k in 1:length(method_param)) {
      run_lengths <- rl_mat[, k]
      edd <- mean(run_lengths)
      sd_rl <- sd(run_lengths)
      edd_conf_int <- conf_int(run_lengths, alpha = 0.05)
      stored_values <- round(c(method_param[k], p, change_param, n_sim, edd, sd_rl,
                               edd_conf_int), digits = 3)
      write(c(method, stored_values), file = edd_file,
            append = TRUE, ncolumns = length(stored_values) + 1)
    }
  }

  store_results <- function(all_run_lengths) {
    methods <- names(all_run_lengths[[1]])
    rl_mats <- lapply(methods, function(method) get_rl_mat(all_run_lengths, method))
    edd_file_rep <- rep(edd_file, length(methods))
    if (change_type == 'cor') method_params <- list(r_pca, r_pca, c_tpca)
    else method_params <- list(r_pca, r_pca, c_tpca, p0s)
    Map(store_edd_results, rl_mats, methods, edd_file_rep, method_params)
  }

  x_train <- train_obj$x
  mu0 <- train_obj$mu
  Sigma0 <- train_obj$Sigma
  cov_mat_nr <- train_obj$nr
  d <- nrow(x_train)
  m <- ncol(x_train)

  Sigma0_hat <- 1 / (m - 1) * x_train %*% t(x_train)
  pca_obj <- tpca::pca(Sigma0_hat)
  V <- pca_obj$vectors
  lambda <- pca_obj$values
  z_train <- V %*% x_train / sqrt(lambda)

  change_type <- get_change_type()
  change_param <- set_change_param(change_type)

  # Method parameters
  axes_max_pca <- lapply(r_pca, function(j) 1:j)
  axes_min_pca <- lapply(r_pca, function(j) (d - j + 1):d)
  tpca_train_info <- get_tpca_train_info(cov_mat_nr, change_type)
  c_tpca <- tpca_train_info$cutoffs
  axes_tpca <- tpca_train_info$axes

  # Initialization
  sim_mixture_rl <- function(p0, x) {
    threshold <- get_threshold(m, d, cov_mat_nr, p0 = p0)
    mixture_rl(threshold, x_train, x, p0, kappa, w)
  }

  sim_tpca_rl <- function(axes, z) {
    threshold <- get_threshold(m, d, cov_mat_nr, axes = axes)
    z_train_sub <- z_train[axes, , drop = FALSE]
    z_sub <- z[axes, , drop = FALSE]
    mixture_rl(threshold, z_train_sub, z_sub, 1, kappa, w)
  }

  log_file <- get_log_name(cov_mat_nr, change_type)
  log_new_change(log_file)
  comp_cluster <- setup_parallel()
  `%dopar%` <- foreach::`%dopar%`
  all_run_lengths <- foreach::foreach(l = 1:n_sim) %dopar% {
  # all_run_lengths <- list()
  # for (l in 1:n_sim) {
    log_current_stage(log_file, change_type, l)
    x <- generate_cor_data(d, n_max, mu0, Sigma0, kappa, p, mu, sigma, rho_scale)
    z <- V %*% x / sqrt(lambda)
    rl_max_pca <- vapply(axes_max_pca, sim_tpca_rl, numeric(1), z = z)
    rl_min_pca <- vapply(axes_min_pca, sim_tpca_rl, numeric(1), z = z)
    rl_tpca <- vapply(axes_tpca, sim_tpca_rl, numeric(1), z = z)
    rl_list <- list('max_pca' = rl_max_pca,
                    'min_pca' = rl_min_pca,
                    'tpca'    = rl_tpca)
    if (any(change_type == c('mean', 'sd'))) {
      rl_mix <- vapply(p0s, sim_mixture_rl, numeric(1), x = x)
      rl_list$mix <- rl_mix
    }
    rl_list
    # all_run_lengths[[l]] <- rl_list
  }
  stop_parallel(comp_cluster)

  store_results(all_run_lengths)
}
