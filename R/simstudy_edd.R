edd_sim <- function(train_obj, n_sim, kappa, p, mu, sigma, rho_scale,
                    p0s, r_pca, threshold_settings) {
  set_edd_file_name <- function(change_type) {
    x <- train_obj$x
    d <- nrow(x)
    m <- ncol(x)
    w <- 200
    alpha_str <- strsplit(as.character(threshold_settings$alpha), '[.]')[[1]][2]
    dir <- './results/'
    paste0(dir, 'edd_', change_type, '_', cov_mat_nr,
           '_n', threshold_settings$n,
           'alpha', alpha_str,
           'm', as.character(m),
           'd', as.character(d),
           'w', as.character(w), '.txt')
  }

  init_edd_file <- function(edd_file) {
    header_items <- c('method', 'method_param',  'p', 'change_param', 'n_sim',
                      'edd', 'sd_edd', 'edd_lower', 'edd_upper')
    write(paste(header_items, collapse = ' '),
          file = edd_file, append = TRUE, ncolumns = length(header_items))
  }

  init_log_files <- function(cov_mat_nr) {
    dir <- './results/'
    log_files  <- c(paste0(dir, 'edd_log_mean_', cov_mat_nr),
                    paste0(dir, 'edd_log_sd_', cov_mat_nr),
                    paste0(dir, 'edd_log_cor_', cov_mat_nr))
    lapply(log_files, function(file) {
      write('+++++++++++++++++++++++++++++++++++', file = file)
      write(paste0('cov_mat nr. ', cov_mat_nr), file = file, append = TRUE)
      write('+++++++++++++++++++++++++++++++++++', file = file, append = TRUE)
    })
  }

  run_mean_change_simulations <- function() {
    for (i in 1:length(p)) {
      for (j in 1:length(mu)) {
        write_global_log(paste0('Running EDD sims for changes in MEAN and cov_mat_nr ',
                                train_obj$nr, '. p = ', p[i],
                                ' and mu = ', mu[j], '.'))
        est_edd_all_methods(train_obj, n_max, n_sim, edd_mean_file,
                            kappa, p[i], w, p0s, r_pca, mu = mu[j],
                            threshold_settings = threshold_settings)
      }
    }
  }

  run_sd_change_simulations <- function() {
    for (i in seq_along(p)) {
      for (j in seq_along(sigma)) {
        write_global_log(paste0('Running EDD sims for changes in VAR and cov_mat_nr ',
                                train_obj$nr, '. p = ', p[i],
                                ' and sigma = ', sigma[j], '.'))
        est_edd_all_methods(train_obj, n_max, n_sim, edd_sd_file,
                            kappa, p[i], w, p0s, r_pca, sigma = sigma[j],
                            threshold_settings = threshold_settings)
      }
    }
  }

  run_cor_change_simulations <- function() {
    # Since we only study decreases in correlation, only the correlated
    # dimensions can be changed. Hence, it's no point in running p for larger
    # values than the proportion of correlated dimensions.
    n_dims_cor <- length(attr(train_obj$Sigma, 'which_dims_cor'))
    d <- ncol(train_obj$Sigma)
    prop_cor_dims <- n_dims_cor / d
    p <- p[p <= prop_cor_dims]

    for (i in 1:length(p)) {
      for (j in 1:length(rho_scale)) {
        write_global_log(paste0('Running EDD sims for changes in COR and cov_mat_nr ',
                                train_obj$nr, '. p = ', p[i],
                                ' and rho_scale = ', rho_scale[j], '.'))
        est_edd_all_methods(train_obj, n_max, n_sim, edd_cor_file,
                            kappa, p[i], w, p0s, r_pca, rho_scale = rho_scale[j],
                            threshold_settings = threshold_settings)
      }
    }
  }

  # Simulation setup
  w <- 200
  n_max <- 2000

  # File setup
  cov_mat_nr <- train_obj$nr
  edd_mean_file <- set_edd_file_name('mean')
  edd_sd_file <- set_edd_file_name('sd')
  edd_cor_file <- set_edd_file_name('cor')
  lapply(c(edd_mean_file, edd_sd_file, edd_cor_file), init_edd_file)
  init_log_files(cov_mat_nr)

  # Simulation runs
  run_mean_change_simulations()
  run_sd_change_simulations()
  run_cor_change_simulations()
}

#' Estimates the edd for a given change for all methods.
#'
#' Description
#'
#' Simulates \code{n_sim} run lengths and calculates the EDD for all
#' methods (Mixture, Min PCA, Max PCA and TPCA) together with other summary
#' statistics of the run length distribution. These values are not returned,
#' but printed to the input \code{edd_file}. The input \code{kappa}, \code{p},
#' \code{mu}, \code{sigma} and \code{rho_scale} constitute the change scenario
#' to estimate the EDD for. \code{p0s} is a vector of parameters for the the
#' mixture procedure, while \code{r_pca} is a vector of the number of projections
#' to retain in Max and Min PCA. The chosen projections for TPCA is
#' found from the file axes-file created during the threshold-finding phase.
#'
#' @param train_obj A list returned from \code{\link{get_training_sets}}.
#' @param n_max The maximum number of monitoring samples.
#' @param n_sim The number of run length simulations to estimate the EDD.
#' @param edd_file A string specifying the file to print results to.
#' @param kappa The change-point.
#' @param p A number between 0 and 1. The proportion of affected streams.
#' @param w The window length.
#' @param p0s A vector of p0-values to be used in the mixture procedure.
#' @param r_pca A vector of the numbers of projections to keep in Max and Min PCA.
#' @param threshold_settings A list with named elements 'n' and 'alpha', giving
#' the settings used to control the probability of false alarm.
#' @param mu A numeric specifying the change in mean.
#' @param sigma A numeric specifying the change in standard deviation.
#' @param rho_scale A numeric specifying the change in correlation.
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
    if (length(axes_list) == 0) return(NULL)
    else return(list('axes' = axes_list, 'cutoffs' = cutoffs))
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
      detection_delay <- run_lengths - kappa
      edd <- mean(detection_delay)
      sd_rl <- sd(detection_delay)
      edd_conf_int <- conf_int(detection_delay, alpha = 0.05)
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
    method_params <- list(r_pca, r_pca)
    if (!is.null(tpca_train_info))
      method_params[[length(method_params) + 1]] <- c_tpca
    if (any(change_type == c('mean', 'sd')))
      method_params[[length(method_params) + 1]] <- p0s
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
  if (!is.null(tpca_train_info)) {
    c_tpca <- tpca_train_info$cutoffs
    axes_tpca <- tpca_train_info$axes
  }

  # Initialization
  sim_mixture_rl <- function(p0, x) {
    threshold <- get_threshold(m, d, cov_mat_nr, p0 = p0)
    mixture_rl(threshold, x_train, x, p0, w)
  }

  sim_tpca_rl <- function(axes, z) {
    threshold <- get_threshold(m, d, cov_mat_nr, axes = axes)
    z_train_sub <- z_train[axes, , drop = FALSE]
    z_sub <- z[axes, , drop = FALSE]
    mixture_rl(threshold, z_train_sub, z_sub, 1, w)
  }

  log_file <- get_log_name(cov_mat_nr, change_type)
  log_new_change(log_file)
  comp_cluster <- setup_parallel()
  `%dopar%` <- foreach::`%dopar%`
  all_run_lengths <- foreach::foreach(l = 1:n_sim) %dopar% {
  # all_run_lengths <- list()
  # for (l in 1:n_sim) {
    log_current_stage(log_file, change_type, l)
    x <- gen_changed_data(d, n_max, mu0, Sigma0, kappa, p, mu, sigma, rho_scale)
    z <- V %*% x / sqrt(lambda)
    rl_max_pca <- vapply(axes_max_pca, sim_tpca_rl, numeric(1), z = z)
    rl_min_pca <- vapply(axes_min_pca, sim_tpca_rl, numeric(1), z = z)
    rl_list <- list('max_pca' = rl_max_pca,
                    'min_pca' = rl_min_pca)
    if (!is.null(tpca_train_info))
      rl_list$tpca <- vapply(axes_tpca, sim_tpca_rl, numeric(1), z = z)
    if (any(change_type == c('mean', 'sd')))
      rl_list$mix <- vapply(p0s, sim_mixture_rl, numeric(1), x = x)
    rl_list
    # all_run_lengths[[l]] <- rl_list
  }
  stop_parallel(comp_cluster)

  invisible(store_results(all_run_lengths))
}
