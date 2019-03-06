run_simstudy <- function() {
  n_sets <- 3
  d <- 20
  m <- 40
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
