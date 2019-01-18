edd_sim <- function(train_obj) {
  set_edd_file_name <- function(change_type) {
    x <- train_obj$x
    d <- nrow(x)
    m <- ncol(x)
    w <- 200
    paste0('edd_', change_type, '_', cov_mat_type,
           '_m', as.character(m),
           'd', as.character(d),
           'w', as.character(w), '.txt')
  }

  init_edd_file <- function(edd_file) {
    header_items <- c('method', 'method_param',  'p', 'change_param', 'n_sim',
                      'edd', 'sd_edd', 'edd_lower', 'edd_upper')
    write(paste(header_items, collapse = ' '),
          file = edd_file, append = TRUE, ncolumns = length(header_items))
  }

  init_log_files <- function(cov_mat_type) {
    log_files  <- c(paste0('edd_log_mean_', cov_mat_type),
                    paste0('edd_log_sd_', cov_mat_type),
                    paste0('edd_log_cor_', cov_mat_type))
    lapply(log_files, function(file) {
      write('+++++++++++++++++++++++++++++++++++', file = file, append = TRUE)
      write(cov_mat_type, file = file, append = TRUE)
      write('+++++++++++++++++++++++++++++++++++', file = file, append = TRUE)
    })
  }

  run_mean_change_simulations <- function() {
    for (i in 1:length(p)) {
      for (j in 1:length(mu)) {
        est_edd_all_methods(train_obj, n_max, n_sim, edd_mean_file,
                            kappa, p[i], w, mu = mu[j])
      }
    }
  }

  run_sd_change_simulations <- function() {
    for (i in seq_along(p)) {
      for (j in seq_along(sigma)) {
        est_edd_all_methods(train_obj, n_max, n_sim, edd_sd_file,
                            kappa, p[i], w, sigma = sigma[j])
      }
    }
  }

  run_cor_change_simulations <- function() {
    for (i in 1:length(p)) {
      for (j in 1:length(rho_scale)) {
        est_edd_all_methods(train_obj, n_max, n_sim, edd_cor_file,
                            kappa, p[i], w, rho_scale = rho_scale[j])
      }
    }
  }

  # Simulation setup
  w <- 200
  n_max <- 2000
  n_sim <- 500

  # Change scenarios
  kappa <- 0
  p <- rev(c(0.02, 0.05, 0.1, 0.3, 0.5, 0.7, 0.9, 0.95, 0.98))
  mu <- c(0.5, 0.7, 1, 1.3)
  sigma <- c(0.5, 0.75, 1.5, 2)
  rho_scale <- c(0, 0.25, 0.5, 0.75)

  # File setup
  cov_mat_type <- train_obj$type
  edd_mean_file <- set_edd_file_name('mean')
  edd_sd_file <- set_edd_file_name('sd')
  edd_cor_file <- set_edd_file_name('cor')
  lapply(c(edd_mean_file, edd_sd_file, edd_cor_file), init_edd_file)
  init_log_files(cov_mat_type)

  # Simulation runs
  run_mean_change_simulations()
  run_sd_change_simulations()
  run_cor_change_simulations()
}

run_edd_sim_dense <- function() {
  edd_sim(get_train_dense(return_all = TRUE))
}

run_edd_sim_halfsparse <- function() {
  edd_sim(get_train_halfsparse(return_all = TRUE))
}

run_edd_sim_sparse <- function() {
  edd_sim(get_train_sparse(return_all = TRUE))
}

run_all_edd_sims <- function() {
  run_edd_sim_dense()
  run_edd_sim_halfsparse()
  run_edd_sim_sparse()
}

run_selected_edd_sims <- function() {
  run_edd_sim_sparse()
}
