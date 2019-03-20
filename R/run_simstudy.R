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
  m <- 200
  training_sets <- get_training_sets(n_sets, d = d, m = m)

  # Method parameters
  p0 <- c(0.03, 0.1, 0.3, 1)
  # r <- c(1, 2, 3, 5, 10, 20)
  r <- c(1, 2, 3, 5)
  # cutoff <- c(0.8, 0.9, 0.95, 0.99, 0.995, 0.999)
  cutoff <- c(0.8, 0.9, 0.95, 0.99)

  # Threshold handling
  n <- 100
  alpha <- 0.01
  # rel_tol <- c(0.2, 0.1, 0.05, 0.025)
  rel_tol <- c(0.2, 0.1, 0.05, 0.025)
  find_all_thresholds(training_sets, n, alpha, rel_tol, p0, r, cutoff)
  extract_final_thresholds(d, m, n, alpha)

  # Change scenarios
  kappa <- 0
  # p <- rev(c(0.02, 0.05, 0.1, 0.3, 0.5, 0.7, 0.9, 0.95, 0.98))
  p <- c(0.9, 0.7, 0.5, 0.3, 0.1)
  mu <- c(0.5, 0.7, 1)
  sigma <- c(0.5, 1.5, 2)
  rho_scale <- c(0, 0.5, 0.75)
  n_sim <- 500
  invisible(lapply(training_sets, edd_sim,
                   n_sim     = n_sim,
                   kappa     = kappa,
                   p         = p,
                   mu        = mu,
                   sigma     = sigma,
                   rho_scale = rho_scale,
                   p0        = p0,
                   r         = r,
                   threshold_settings = list('n' = n, 'alpha' = alpha)))
}

write_global_log <- function(text) {
  # text indicates where the simulations are at.
  curr_time <- Sys.time()
  write(paste0(text, ' Current time: ', curr_time, '. Time used so far: ',
               round(difftime(curr_time, start_time_G, units = 'hour'), 2), ' hours.'),
        log_file_G, append = TRUE)
