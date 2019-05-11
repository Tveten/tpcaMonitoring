write_global_log <- function(text) {
  # text indicates where the simulations are at.
  # Assumes that a global variable named log_file_G exists with the name of the
  # global log file.
  curr_time <- Sys.time()
  write(paste0(text, ' Current time: ', curr_time, '. Time used so far: ',
               round(difftime(curr_time, start_time_G, units = 'hour'), 2), ' hours.'),
        log_file_G, append = TRUE)
}

make_dirs <- function() {
  # Creates the following directories in the current directory:
  # ./thresholds
  #     /axes
  #     /threshold_results
  # ./results
  #     /figures

  dir.create('thresholds', showWarnings = FALSE)
  dir.create('thresholds/axes', showWarnings = FALSE)
  dir.create('thresholds/threshold_results', showWarnings = FALSE)
  dir.create('results', showWarnings = FALSE)
  dir.create('results/figures', showWarnings = FALSE)
}

#' Reproduce simulation study with D = 100.
#'
#' This function reproduces the entire simulation study in
#' "Online Detection of Sparse Changes in High-Dimensional Data Streams Using Tailored Projections"
#' for D = 100 (Section 4). On an average computer it takes about a week to run on four cores.
#' Thus, the files with the results of running \code{run_simstudy} used in the paper
#' are included in the package or found at github.com/Tveten/tpcaMonitoring.
#'
#' The following directories are created in the current directory:
#' ./thresholds,
#' ./thresholds/axes,
#' .thresholds/threshold_results,
#' ./results,
#' ./results/figures.
#' Files created to save information about thresholds are saved in the threshold
#' directory, while EDD results are stored in "results".
#' The results we obtained by running \code{run_simstudy} are stored in the
#' same structure of directories on Github.
#'
#' @export
run_simstudy <- function() {
  init_global_log <- function() {
    alpha_str <- strsplit(as.character(alpha), '[.]')[[1]][2]
    log_file_G <<- paste0('overall_log_', 'n', n, 'alpha', alpha_str,
                               'm', m, 'd', d, 'txt')
    start_time_G <<- Sys.time()
    write(paste0('Simulation start at ', start_time_G), log_file_G, append = TRUE)

    # Log every
    #   - find_mixture_thresholds
    #   - find_tpca_thresholds
    #   - In edd_sim:
    #      * run_mean_change_simulations
    #      * run_sd_change_simulations
    #      * run_cor_change_simulations
  }

  n_sets <- 30
  d <- 100
  m <- 2 * d
  training_sets <- get_training_sets(n_sets, d = d, m = m)

  # Method parameters
  p0 <- c(0.03, 0.1, 0.3, 1)
  J <- c(1, 2, 3, 5, 10, 20)
  cutoff <- c(0.8, 0.9, 0.95, 0.99, 0.995, 0.999)

  # Threshold handling
  n <- 100
  alpha <- 0.01
  rel_tol <- c(0.2, 0.1, 0.05, 0.025)
  init_global_log()
  make_dirs()
  find_all_thresholds(training_sets, n, alpha, rel_tol, p0, J, cutoff)
  extract_final_thresholds(d, m, n, alpha)

  # Change scenarios
  kappa <- 0
  p <- rev(c(0.02, 0.05, 0.1, 0.3, 0.5, 0.7, 0.9, 0.95, 0.98))
  mu <- c(0.5, 0.7, 1, 1.3)
  sigma <- c(0.5, 0.75, 1.5, 2)
  rho_scale <- c(0, 0.25, 0.5, 0.75)
  n_sim <- 500
  invisible(lapply(training_sets, edd_sim,
                   n_sim     = n_sim,
                   kappa     = kappa,
                   p         = p,
                   mu        = mu,
                   sigma     = sigma,
                   rho_scale = rho_scale,
                   p0        = p0,
                   r         = J,
                   threshold_settings = list('n' = n, 'alpha' = alpha)))
}
