get_single_edd <- function(change_type, cov_mat_nr, n, alpha, m, d,
                           col_name = 'edd') {
  # Returns a list with each method being an entry: max_pca, min_pca, tpca, mix.
  # Each entry is a 3-array with EDD for  (p, change_size, method_param)

  extract_EDD <- function(method, method.params) {
    a <- length(p)
    b <- length(change_param)
    c <- length(method_params)
    edd_array <- array(NA, dim = c(a, b, c))
    for (i in 1:a) {
      for (j in 1:b) {
        for (k in 1:c) {
          index <- (edd_results$method == method &
                      edd_results$p == p[i] &
                      edd_results$change_param == change_param[j] &
                      edd_results$method_param == method_params[k])
          if (sum(index) >= 1) edd_array[i, j, k] <- edd_results[[col_name]][index][1]
          else {
            if (method == 'mix') next
            else {
              warning(paste('No data on', method, p[i], change_param[j], method_params[k]))
              edd_array[i, j, k] <- NA
            }
          }
        }
      }
    }
    rownames(edd_array) <- sapply(p, function(x) paste('p =', x))
    colnames(edd_array) <- sapply(change_param, function(x) paste('change_size =', x))
    dimnames(edd_array)[[3]] <- sapply(method_params, function(x) paste('method_param =', x))
    return(edd_array)
  }

  path <- './results/'
  alpha_str <- strsplit(as.character(alpha), '[.]')[[1]][2]
  file_name <- paste0('edd_', change_type, '_', cov_mat_nr, '_n', n,
                      'alpha', alpha_str, 'm', m, 'd', d, 'w200.txt')
  edd_results <- read.table(paste0(path, file_name), header = TRUE)
  p <- sort(unique(edd_results$p))
  # p <- c(0.02, 0.05, 0.1, 0.3, 0.5)
  change_param <- sort(unique(edd_results$change_param))
  methods <- as.character(unique(edd_results$method))

  edd_list <- list()
  for (i in 1:length(methods)) {
    method_ind <- edd_results$method == methods[i]
    method_params <- sort(unique(edd_results$method_param[method_ind]))
    edd_list[[methods[i]]] <- extract_EDD(methods[i], method_params)
  }
  return(edd_list)
}

get_summary_edd <- function(change_type, n, alpha, m, d,
                            cov_mat_type = NULL,
                            summary_func = mean,
                            col_name = 'edd') {
  get_n_files <- function() {
    path <- './results/'
    alpha_str <- strsplit(as.character(alpha), '[.]')[[1]][2]
    file_ids <- paste0('edd_', change_type, '\\w+', 'n', n,
                    'alpha', alpha_str, 'm', m, 'd', d, 'w200.txt')
    files <- list.files(path, file_ids)
    length(files)
  }

  get_cov_mat_nr <- function() {
    if (is.null(cov_mat_type)) {
      return(1:n_files)
    } else if (cov_mat_type == 'low') {
      return(16:n_files)
    } else if (cov_mat_type == 'high') {
      return(1:15)
    }
  }

  n_files <- get_n_files()
  cov_mat_nrs <- get_cov_mat_nr()
  edd_lists <- lapply(cov_mat_nrs, function(i) {
    get_single_edd(change_type, i, n, alpha, m, d, col_name)
  })

  methods <- names(edd_lists[[1]])
  edd_arrays <- lapply(methods, function(method) {
    array(0, dim = c(dim(edd_lists[[1]][[method]]), n_files))
  })
  names(edd_arrays) <- methods
  for (i in seq_along(edd_lists)) {
    for (j in seq_along(methods)) {
      edd_arrays[[methods[j]]][, , , i] <- edd_lists[[i]][[methods[j]]]
    }
  }

  avg_edds <- lapply(methods, function(method) {
    avg_edd <- apply(edd_arrays[[method]], 1:3, summary_func)
    dimnames(avg_edd) <- dimnames(edd_lists[[1]][[method]])
    avg_edd
  })
  names(avg_edds) <- methods
  avg_edds
}

get_edd <- function(change_type, n, alpha, m, d,
                    cov_mat_nr = NULL, cov_mat_type = NULL,
                    summary_func = mean, col_name = 'edd') {
  if (is.null(cov_mat_nr))
    edd_list <- get_summary_edd(change_type, n, alpha, m, d, cov_mat_type,
                                summary_func, col_name)
  else
    edd_list <- get_single_edd(change_type, cov_mat_nr, n, alpha, m, d, col_name)
  edd_list
}

get_single_rl_confint <- function(change_type, cov_mat_nr, n, alpha, m, d) {
  # Returns a list with each method being an entry: max_pca, min_pca, tpca, mix.
  # Each entry is a 3-array with EDD for  (p, change_size, method_param)

  extract_EDD <- function(method, method.params) {
    a <- length(p)
    b <- length(change_param)
    c <- length(method_params)
    edd_array <- array(NA, dim = c(a, b, c))
    for (i in 1:a) {
      for (j in 1:b) {
        for (k in 1:c) {
          index <- (edd_results$method == method &
                      edd_results$p == p[i] &
                      edd_results$change_param == change_param[j] &
                      edd_results$method_param == method_params[k])
          if (sum(index) >= 1) {
            edd <- edd_results$edd[index][1]
            edd_lower <- round(edd_results$edd_lower[index][1])
            edd_upper <- round(edd_results$edd_upper[index][1])
            edd_array[i, j, k] <- paste0(edd, ' (', edd_lower, ', ', edd_upper, ')')
          } else {
            if (method == 'mix') next
            else {
              warning(paste('No data on', method, p[i], change_param[j], method_params[k]))
              edd_array[i, j, k] <- NA
            }
          }
        }
      }
    }
    rownames(edd_array) <- sapply(p, function(x) paste('p =', x))
    colnames(edd_array) <- sapply(change_param, function(x) paste('change_size =', x))
    dimnames(edd_array)[[3]] <- sapply(method_params, function(x) paste('method_param =', x))
    return(edd_array)
  }

  path <- './results/'
  alpha_str <- strsplit(as.character(alpha), '[.]')[[1]][2]
  file_name <- paste0('edd_', change_type, '_', cov_mat_nr, '_n', n,
                      'alpha', alpha_str, 'm', m, 'd', d, 'w200.txt')
  edd_results <- read.table(paste0(path, file_name), header = TRUE)
  p <- sort(unique(edd_results$p))
  change_param <- sort(unique(edd_results$change_param))
  methods <- as.character(unique(edd_results$method))

  edd_list <- list()
  for (i in 1:length(methods)) {
    method_ind <- edd_results$method == methods[i]
    method_params <- sort(unique(edd_results$method_param[method_ind]))
    edd_list[[methods[i]]] <- extract_EDD(methods[i], method_params)
  }
  return(edd_list)
}

get_avg_rl_confint <- function(change_type, cov_mat_nr, n, alpha, m, d) {
  # Returns a list with each method being an entry: max_pca, min_pca, tpca, mix.
  # Each entry is a 3-array with EDD for  (p, change_size, method_param)

  EDD <- function(method, method.params) {
    a <- length(p)
    b <- length(change_param)
    c <- length(method_params)
    edd_array <- array(NA, dim = c(a, b, c))
    for (i in 1:a) {
      for (j in 1:b) {
        for (k in 1:c) {
          edd_array[i, j, k] <- paste0(edd, ' (', edd_lower, ', ', edd_upper, ')')
        }
      }
    }
    rownames(edd_array) <- sapply(p, function(x) paste('p =', x))
    colnames(edd_array) <- sapply(change_param, function(x) paste('change_size =', x))
    dimnames(edd_array)[[3]] <- sapply(method_params, function(x) paste('method_param =', x))
    return(edd_array)
  }

  # path <- './results/'
  # alpha_str <- strsplit(as.character(alpha), '[.]')[[1]][2]
  # file_name <- paste0('edd_', change_type, '_', cov_mat_nr, '_n', n,
  #                     'alpha', alpha_str, 'm', m, 'd', d, 'w200.txt')
  # edd_results <- read.table(paste0(path, file_name), header = TRUE)
  # p <- sort(unique(edd_results$p))
  # change_param <- sort(unique(edd_results$change_param))
  # methods <- as.character(unique(edd_results$method))

  avg_edd <- get_summary_edd(change_type, n, alpha, m, d, col_name = 'edd')
  avg_lower_rl <- get_summary_edd(change_type, n, alpha, m, d, col_name = 'edd_lower')
  avg_upper_rl <- get_summary_edd(change_type, n, alpha, m, d, col_name = 'edd_upper')

  rl_summary_list <- avg_edd
  for (i in 1:length(rl_summary_list)) {
    for (i in 1:dim(rl_summary_list)) {
      for (j in 1:b) {
        for (k in 1:c) {
          edd_array[i, j, k] <- paste0(edd, ' (', edd_lower, ', ', edd_upper, ')')
        }
      }
    }
  }
  return(edd_list)
}


##### FIGURES ----------------------------------------------------------------
ggplot_edd <- function(edd_list, change_type, change_param, ylim = c(0, 2000)) {
  get_change_type <- function(file) {
    strsplit(file, '_')[[1]][2]
  }

  get_all_methods <- function() {
    c('mix', 'max_pca', 'min_pca', 'tpca')
  }

  set_col <- function(methods) {
    all_col <- list('darkgreen', 'darkgoldenrod', 'cornflowerblue', 'red')
    names(all_col) <- get_all_methods()
    unlist(all_col[methods])
  }

  set_labels <- function(methods) {
    all_labels <- list('Mixture', 'Max PCA', 'Min PCA', 'TPCA')
    names(all_labels) <- get_all_methods()
    unlist(all_labels[methods])
  }

  set_title <- function(change_type) {
    title_str <- ''
    # if (cov_mat_title)
    #   title_str <- paste0(title_str, first_up(cov_mat_type), ' $\\mathbf{\\Sigma}_0$, ')
    if (change_type == 'mean') {
      title_str <- paste0(title_str, sprintf('$\\mu_d = %.1f$', change_param))
    }
    if (change_type == 'sd') {
      title_str <- paste0(title_str, sprintf('$\\sigma_d = %.2f$', change_param))
    }
    if (change_type == 'cor') {
      title_str <- paste0(title_str, sprintf('$a_{di} = %.2f$', change_param))
      # rho_char <- '\u03C1'
      # title_str <- paste0(title_str, 'cor ', rho_char, ' -> ', change_param, rho_char)
    }
    latex2exp::TeX(title_str)
  }

  extract_double <- function(x, split = '= ') {
    x_temp <- strsplit(as.character(x), split)
    vapply(x_temp, function(x) as.double(x[2]), numeric(1))
  }

  # edd_list <- get_edd(file)
  # change_type <- get_change_type(file)
  methods <- names(edd_list)
  col <- set_col(methods)
  labels <- set_labels(methods)
  title <- set_title(change_type)

  edd_df <- reshape2::melt(edd_list, varnames = c('p', 'change_size', 'method_param'),
                           value.name = 'edd', level = 'method')
  if (change_type == 'cor') edd_df <- subset(edd_df, Lmethod != 'mix')
  edd_df$p <- log(extract_double(edd_df$p))
  edd_df$change_size <- extract_double(edd_df$change_size)
  edd_df$method_param <- as.factor(extract_double(edd_df$method_param))
  edd_subset <- subset(edd_df, change_size == change_param)

  ggplot2::ggplot(edd_subset, ggplot2::aes(x = p, y = edd, color = Lmethod, linetype = method_param)) +
    ggplot2::geom_line() +
    ggplot2::theme_light() +
    ggplot2::labs(title = title, x = 'log(p)', y = 'EDD') +
    ggplot2::coord_cartesian(ylim = ylim) +
    ggplot2::scale_color_manual('', values = col, labels = labels) +
    ggplot2::scale_linetype_manual(values = rep(1, length(unique(edd_df$method_param))), guide = FALSE)
}

#' Plot EDDS results.
#'
#' \code{multiplot_edd_summary} creates a grid of plots with results from
#' running \code{\link{run_simstudy}}
#'
#' @param n The monitoring length the probability of false alarm is controlled for.
#' @param alpha The probability of false alarms
#' @param m The number of training samples.
#' @param d The data dimension.
#' @param cov_mat_nr The number given the covariance matrix when generated through \code{\link{gen_train}}
#' @param cov_mat_type Either 'low' or 'high', to get results only for the low
#' correlation or high correlation matrices.
#' @param summary_func The function to summarize the EDDs over the cov_mat_types.
#' @param ylim A vector with lower and upper bounds on the EDD values shown in the plots.
#'
#' @examples
#' # Plot of the results for low correlation matrices:
#' multiplot_edd_summary(cov_mat_type = 'low')
#'
#' # Plot of the results for low correlation matrices:
#' multiplot_edd_summary(cov_mat_type = 'high')
#'
#' # Plot of the EDD results for a single covariance matrix/training set:
#' multiplot_edd_summary(cov_mat_nr = 5)
#'
#' @export
multiplot_edd_summary <- function(n = 100, alpha = 0.01, m = 200, d = 100,
                                  cov_mat_nr    = NULL,
                                  cov_mat_type  = NULL,
                                  change_params = NULL,
                                  summary_func  = mean,
                                  ylim          = c(0, 300)) {
  #  summary_func: Function to summarize over cov_mats.
  change_types <- c('mean', 'sd', 'cor')
  if (is.null(change_params))
    change_params <- list('mean' = c(0.5, 0.7, 1),
                          'sd' = c(0.5, 1.5, 2),
                          'cor' = c(0, 0.5, 0.75))
  edd_lists <- lapply(change_types, function(x) {
    get_edd(x, n, alpha, m, d, cov_mat_nr, cov_mat_type, summary_func)
  })

  n <- length(edd_lists)
  plots <- list()
  for (i in seq_along(change_types)) {
    l <- length(change_params)
    plots[[i]] <- lapply(change_params[[i]], function(param) {
      ggplot_edd(edd_lists[[i]], change_types[i], param, ylim = ylim)
    })
  }
  plots <- unlist(plots, recursive = FALSE)

  ggpubr::ggarrange(plotlist = plots, ncol = n, nrow = ceiling(length(plots) / n),
                    common.legend = TRUE, legend = "bottom")
}

save_summary_plot <- function(n, alpha, m, d,
                              cov_mat_nr   = NULL,
                              cov_mat_type = NULL,
                              summary_func = mean,
                              ylim         = c(0, 300),
                              extension    = 'png') {
  base_width <- 2.88
  base_height <- 2.9
  n_col <- 3
  n_row <- 3
  width <- n_col * base_width
  height <- n_row * base_height
  plot_obj <- multiplot_edd_summary(n, alpha, m, d, cov_mat_nr,
                                    cov_mat_type, summary_func, ylim)
  show(plot_obj)
  dir <- './results/figures/'
  alpha_str <- strsplit(as.character(alpha), '[.]')[[1]][2]
  if (is.null(cov_mat_nr)) {
    cov_mat_str <- as.character(substitute(summary_func))
    if (!is.null(cov_mat_type))
      cov_mat_str <- paste0(cov_mat_str, '-', cov_mat_type)
  } else cov_mat_str <- as.character(cov_mat_nr)
  file_name <- paste0('edd_summary_', cov_mat_str, '_n', n,
                      'alpha', alpha_str, 'm', m, 'd', d,
                      '_ylim', ylim[1], '-', ylim[2])
  ggplot2::ggsave(paste0(dir, file_name, '.', extension),
                  plot = plot_obj,
                  width = width, height = height, units = 'in')
}
#' Plot EDD results with confidence intervals.
#'
#' See the documentation for \code{\link{multiplot_edd_summary}} for an explanation
#' of the arguments. The additional ones are described here.
#'
#' @param change_type Either 'mean', 'sd' or 'cor'
#' @param change_param A numeric specifying the change size for the chosen change_type.
#'
#' @export
ggplot_edd_conf_int <- function(n, alpha, m, d, change_type, change_param,
                                ci_alpha     = 0.05,
                                cov_mat_nr   = NULL,
                                cov_mat_type = NULL,
                                summary_func = mean,
                                ylim         = c(0, 300),
                                extension    = 'png') {
  get_change_type <- function(file) {
    strsplit(file, '_')[[1]][2]
  }

  get_all_methods <- function() {
    c('mix', 'max_pca', 'min_pca', 'tpca')
  }

  set_col <- function(methods) {
    all_col <- list('darkgreen', 'darkgoldenrod', 'cornflowerblue', 'red')
    names(all_col) <- get_all_methods()
    unlist(all_col[methods])
  }

  set_labels <- function(methods) {
    all_labels <- list('Mixture', 'Max PCA', 'Min PCA', 'TPCA')
    names(all_labels) <- get_all_methods()
    unlist(all_labels[methods])
  }

  set_title <- function(change_type) {
    title_str <- ''
    # if (cov_mat_title)
    #   title_str <- paste0(title_str, first_up(cov_mat_type), ' $\\mathbf{\\Sigma}_0$, ')
    if (change_type == 'mean') {
      title_str <- paste0(title_str, sprintf('$\\mu_d = %.1f$', change_param))
    }
    if (change_type == 'sd') {
      title_str <- paste0(title_str, sprintf('$\\sigma_d = %.2f$', change_param))
    }
    if (change_type == 'cor') {
      title_str <- paste0(title_str, sprintf('$a_{di} = %.2f$', change_param))
      # rho_char <- '\u03C1'
      # title_str <- paste0(title_str, 'cor ', rho_char, ' -> ', change_param, rho_char)
    }
    latex2exp::TeX(title_str)
  }

  extract_double <- function(x, split = '= ') {
    x_temp <- strsplit(as.character(x), split)
    vapply(x_temp, function(x) as.double(x[2]), numeric(1))
  }

  # n <- 100
  # alpha <- 0.01
  # m <- 200
  # d <- 100
  # ci_alpha <- 0.05
  # change_type <- 'mean'
  # cov_mat_nr <- 1
  # cov_mat_type <- NULL
  # summary_func <- mean
  # ylim <- c(0, 300)

  edd_list <- get_edd(change_type, n, alpha, m, d, cov_mat_nr, cov_mat_type,
                      summary_func)
  sd_edd_list <- get_edd(change_type, n, alpha, m, d, cov_mat_nr, cov_mat_type,
                         summary_func, col_name = 'sd_edd')

  edd_df <- reshape2::melt(edd_list, varnames = c('p', 'change_size', 'method_param'),
                           value.name = 'edd', level = 'method')
  sd_edd_df <- reshape2::melt(sd_edd_list, varnames = c('p', 'change_size', 'method_param'),
                              value.name = 'sd_edd', level = 'method')
  z_alpha <- qnorm(1 - ci_alpha/2)
  edd_df$ci <- z_alpha * sd_edd_df$sd_edd / sqrt(500)

  if (change_type == 'cor') edd_df <- subset(edd_df, Lmethod != 'mix')
  edd_df$p <- log(extract_double(edd_df$p))
  edd_df$change_size <- extract_double(edd_df$change_size)
  edd_df$method_param <- as.factor(extract_double(edd_df$method_param))
  edd_subset <- subset(edd_df, change_size == change_param)
  edd_subset <- subset(edd_subset, method_param == 2 | method_param == 20 | method_param == 0.03 | method_param == 0.3 | method_param == 0.9 | method_param == 0.995)

  methods <- names(edd_list)
  col <- set_col(methods)
  labels <- set_labels(methods)
  title <- set_title(change_type)
  pd <- suppressWarnings(ggplot2::position_dodge(0.05))

  suppressWarnings(ggplot2::ggplot(edd_subset, ggplot2::aes(x = p, y = edd, color = Lmethod, linetype = method_param)) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = edd - ci,
                                                         ymax = edd + ci),
                                            width = 0.5, position = pd) +
    ggplot2::geom_line(position = pd) +
    ggplot2::theme_light() +
    ggplot2::labs(title = title, x = 'log(p)', y = 'EDD') +
    ggplot2::coord_cartesian(ylim = ylim) +
    ggplot2::scale_color_manual('', values = col, labels = labels) +
    ggplot2::scale_linetype_manual(values = rep(1, length(unique(edd_df$method_param))), guide = FALSE))
}

multiplot_ci <- function(n, alpha, m, d, change_param,
                         ci_alpha     = 0.05,
                         cov_mat_nr   = NULL,
                         cov_mat_type = NULL,
                         summary_func = mean,
                         ylim         = c(0, 300)) {
  change_types <- c('mean', 'sd', 'cor')
  plots <- lapply(seq_along(change_types), function(i) {
    ggplot_edd_conf_int(n, alpha, m, d, change_types[i], change_param[i],
                        ci_alpha, cov_mat_nr, cov_mat_type,
                        summary_func, ylim, extension)
  })
  ggpubr::ggarrange(plotlist = plots, ncol = length(change_types), nrow = 1,
                    common.legend = TRUE, legend = "bottom")
}

save_ci_plot <- function(n, alpha, m, d, change_param,
                         ci_alpha     = 0.05,
                         cov_mat_nr   = NULL,
                         cov_mat_type = NULL,
                         summary_func = mean,
                         ylim         = c(0, 300),
                         extension    = 'png') {
  base_width <- 2.88
  base_height <- 3.3
  n_col <- 3
  n_row <- 1
  width <- n_col * base_width
  height <- n_row * base_height
  plot_obj <- multiplot_ci(n, alpha, m, d, change_param, ci_alpha,
                           cov_mat_nr, cov_mat_type, summary_func, ylim)
  show(plot_obj)
  dir <- './results/figures/'
  alpha_str <- strsplit(as.character(alpha), '[.]')[[1]][2]
  if (is.null(cov_mat_nr)) {
    cov_mat_str <- as.character(substitute(summary_func))
    if (!is.null(cov_mat_type))
      cov_mat_str <- paste0(cov_mat_str, '-', cov_mat_type)
  } else cov_mat_str <- as.character(cov_mat_nr)
  file_name <- paste0('edd_ci_', cov_mat_str, '_n', n,
                      'alpha', alpha_str, 'm', m, 'd', d,
                      '_ylim', ylim[1], '-', ylim[2])
  ggplot2::ggsave(paste0(dir, file_name, '.', extension),
                  plot = plot_obj,
                  width = width, height = height, units = 'in')
}

ggplot_compare_edd_dim <- function(n, alpha, change_type, change_param,
                                   ylim = c(0, 300)) {
  get_change_type <- function(file) {
    strsplit(file, '_')[[1]][2]
  }

  get_all_methods <- function() {
    c('mix', 'max_pca', 'min_pca', 'tpca')
  }

  set_col <- function(methods) {
    all_col <- list('darkgreen', 'darkgoldenrod', 'cornflowerblue', 'red')
    names(all_col) <- get_all_methods()
    unlist(all_col[methods])
  }

  set_labels <- function(methods) {
    all_labels <- list('Mixture', 'Max PCA', 'Min PCA', 'TPCA')
    names(all_labels) <- get_all_methods()
    unlist(all_labels[methods])
  }

  set_title <- function(change_type) {
    title_str <- ''
    # if (cov_mat_title)
    #   title_str <- paste0(title_str, first_up(cov_mat_type), ' $\\mathbf{\\Sigma}_0$, ')
    if (change_type == 'mean') {
      title_str <- paste0(title_str, sprintf('$\\mu_d = %.1f$', change_param))
    }
    if (change_type == 'sd') {
      title_str <- paste0(title_str, sprintf('$\\sigma_d = %.2f$', change_param))
    }
    if (change_type == 'cor') {
      title_str <- paste0(title_str, sprintf('$a_{di} = %.2f$', change_param))
      # rho_char <- '\u03C1'
      # title_str <- paste0(title_str, 'cor ', rho_char, ' -> ', change_param, rho_char)
    }
    latex2exp::TeX(title_str)
  }

  extract_double <- function(x, split = '= ') {
    x_temp <- strsplit(as.character(x), split)
    vapply(x_temp, function(x) as.double(x[2]), numeric(1))
  }

  m <- c(200, 1000)
  d <- c(100, 500)

  edd_list100 <- get_edd(change_type, n, alpha, m[1], d[1], 15)
  edd_list500 <- get_edd(change_type, n, alpha, m[2], d[2], 1)

  edd_df100 <- reshape2::melt(edd_list100, varnames = c('p', 'change_size', 'method_param'),
                              value.name = 'edd', level = 'method')
  edd_df100$D <- rep(d[1], nrow(edd_df100))
  edd_df500 <- reshape2::melt(edd_list500, varnames = c('p', 'change_size', 'method_param'),
                              value.name = 'edd', level = 'method')
  edd_df500$D <- rep(d[2], nrow(edd_df500))
  edd_df <- rbind(edd_df100, edd_df500)

  if (change_type == 'cor') edd_df <- subset(edd_df, Lmethod != 'mix')
  edd_df$p <- log(extract_double(edd_df$p))
  edd_df$change_size <- extract_double(edd_df$change_size)
  edd_df$method_param <- as.factor(extract_double(edd_df$method_param))
  edd_df$D <- as.factor(edd_df$D)
  edd_subset <- subset(edd_df, change_size == change_param)
  edd_subset <- subset(edd_subset, method_param == 5 | method_param == 0.03 | method_param == 0.99)

  methods <- names(edd_list)
  col <- set_col(methods)
  labels <- set_labels(methods)
  title <- set_title(change_type)

  ggplot2::ggplot(edd_subset, ggplot2::aes(x = p, y = edd, color = Lmethod, linetype = D)) +
    ggplot2::geom_line() +
    ggplot2::theme_light() +
    ggplot2::labs(title = title, x = 'log(p)', y = 'EDD') +
    ggplot2::coord_cartesian(ylim = ylim) +
    ggplot2::scale_color_manual('', values = col, labels = labels)
    # ggplot2::scale_linetype_manual(values = rep(1, length(unique(edd_df$method_param))), guide = FALSE)
}

multiplot_compare_dim <- function(n, alpha, change_param = c(0.5, 1.5, 0),
                                 ylim = c(0, 300)) {
  change_types <- c('mean', 'sd', 'cor')
  plots <- lapply(seq_along(change_types), function(i) {
    ggplot_compare_edd_dim(n, alpha, change_types[i], change_param[i], ylim)
  })
  ggpubr::ggarrange(plotlist = plots, ncol = length(change_types), nrow = 1,
                    common.legend = TRUE, legend = "bottom")
}

save_compare_plot <- function(n, alpha, change_param = c(0.5, 1.5, 0),
                              ylim = c(0, 300), extension = 'png') {
  base_width <- 2.88
  base_height <- 3.3
  n_col <- 3
  n_row <- 1
  width <- n_col * base_width
  height <- n_row * base_height
  plot_obj <- multiplot_compare_dim(n, alpha, change_param, ylim)
  show(plot_obj)
  dir <- './results/figures/'
  alpha_str <- strsplit(as.character(alpha), '[.]')[[1]][2]
  file_name <- paste0('edd_compare_dim100-500', '_n', n,
                      'alpha', alpha_str,
                      '_ylim', ylim[1], '-', ylim[2])
  ggplot2::ggsave(paste0(dir, file_name, '.', extension),
                  plot = plot_obj,
                  width = width, height = height, units = 'in')
}

##### TABLES -----------------------------------------------------------------
summarize_edd <- function(change_type, n, alpha, m, d,
                          summary_func = mean, cov_mat_nr = NULL,
                          cov_mat_type = NULL, ...) {
  #  summary_func: Function to summarize over change size and change sparsity.

  edd_list <- get_edd(change_type, n, alpha, m, d, cov_mat_nr, cov_mat_type, summary_func)
  method <- names(edd_list)
  n_params <- lapply(method, function(i) dim(edd_list[[i]])[3])
  names(n_params) <- method

  edd_summary <- list()
  for (i in seq_along(method)) {
    if (!is.null(edd_list[[method[i]]])) {
      edd_summary[[method[i]]] <- apply(edd_list[[method[i]]], 3, summary_func, ...)
    } else {
      edd_summary[[method[i]]] <- rep(NA, n_params[[method[i]]])
    }
  }
  edd_summary
}

min_edd_summary <- function(n, alpha, m, d, summary_func = mean,
                            cov_mat_nr = NULL, cov_mat_type = NULL,
                            show = FALSE, ...) {
  # Summarizes EDD for each change_type for a given cov_mat_type.
  # Picks out the minimizer among method_parameters.
  # TODO: Combined table with all cov_mat_types separately, with method_param in parenthesis.

  which_method_approx_min <- function(y) {
    # Finds the minimizing arguments of y within a tolerance of 0.5.
    # Here: Minimum of average EDDs.
    tol <- 1
    if (length(which.min(y)) != 0) {
      ind <- which(abs(y - min(y)) <= tol)
      params_temp <- names(y)[ind]
      params_temp <- strsplit(params_temp, ' = ')
      params <- vapply(params_temp, `[`, character(1), 2)
      min_char <- paste0('method_param = ', paste(params, collapse = ', '))
      return(min_char)
    } else return(NA)
  }

  change_types <- c('mean', 'sd', 'cor')
  mean_edd <- lapply(change_types, function(x) {
    summarize_edd(x, n, alpha, m, d, summary_func, cov_mat_nr, cov_mat_type)
  })
  names(mean_edd) <- change_types
  if (show) print(mean_edd)

  min_mean_edd <- lapply(mean_edd, function(x) as.data.frame(lapply(x, min)))
  min_mean_edd$cor$mix <- NA
  min_mean_edd <- do.call('rbind', min_mean_edd)

  which_min_mean_edd <- lapply(mean_edd, function(x) {
    as.data.frame(lapply(x, which_method_approx_min), stringsAsFactors = FALSE)
  })
  which_min_mean_edd$cor$mix <- NA
  which_min_mean_edd <- do.call('rbind', which_min_mean_edd)
  return(list('min' = min_mean_edd, 'which_min' = which_min_mean_edd))
}

edd_summary_table <- function(n, alpha, m, d, summary_func = mean,
                              cov_mat_nr = NULL, cov_mat_type = NULL, ...) {
  min_edd_list <- min_edd_summary(n, alpha, m, d, summary_func = summary_func,
                                  cov_mat_nr = cov_mat_nr,
                                  cov_mat_type = cov_mat_type,  ...)
  min_edd <- min_edd_list$min
  which_min_edd <- min_edd_list$which_min

  caption <- 'Average EDD per change type for each method\'s best method parameters (in parenthesis). '
  caption <- paste0(caption, 'Each listed method parameter is within 1 of the method\'s minimum average EDD.')
  caption <- paste0(caption, 'The average is taken over change sparsity, change size and the 30 training sets.')
  label <- 'tab:EDD_summary_dense'
  date_line <- paste0('% ', date())
  begin_table <- paste('\\begin{table}[htb]',
                       paste0('\\caption{', caption, '}'),
                       paste0('\\label{', label ,'}'),
                       '\\centering',
                       '\\begin{tabular}{lllll}',
                       '\\toprule', sep = ' \n')
  end_table <- paste('\\bottomrule',
                     '\\end{tabular}',
                     '\\end{table}', sep = ' \n')

  headings1 <- paste('\\multirow{2}{*}{Change type}',
                    '\\multicolumn{4}{c}{EDD} \\\\', sep = ' & ')
  headings1 <- paste0(headings1, '\n\\cmidrule{2-5}')
  shown_headings <- list('mix'     = 'Mixture($p_0$)',
                         'max_pca' = 'Max PCA($J$)',
                         'min_pca' = 'Min PCA($J$)',
                         'tpca'    = 'TPCA$(c)$')
  headings2 <- paste(unlist(shown_headings[colnames(min_edd)]), collapse = ' & ')
  headings2 <- paste0(' & ', headings2, '\\\\')
  mid_sep <- ' \n\\midrule'

  shown_type_names <- list('mean' = 'Mean', 'sd' = 'Variance', 'cor' = 'Correlation')
  latex_table <- paste(date_line, begin_table, headings1, headings2, sep = ' \n')
  latex_table <- paste0(latex_table, mid_sep)
  for (i in 1:nrow(min_edd)) {
    type <- shown_type_names[[rownames(min_edd)[i]]]
    table_line <- type
    for (j in 1:ncol(min_edd)) {
      edd <- round(min_edd[i, j], 1)
      which_param <- unlist(strsplit(which_min_edd[i, j], ' = '))[2]
      table_line <- paste0(table_line, ' & ', edd, ' (', which_param, ')')
      if (j == ncol(min_edd))
        table_line <- paste0(table_line, ' \\\\')
    }
    latex_table <- paste(latex_table, table_line, sep = ' \n')
  }
  latex_table <- paste0(latex_table, ' \n', end_table)
  cat(latex_table)
}

##### TESTS ------------------------------------------------------------------
test_multiplot <- function(cov_mat_type, change_type) {
  multiplot_edd(cov_mat_type, change_type)
}

test_all_plots <- function() {
  cov_mat_types <- c('dense', 'halfsparse', 'sparse')
  change_types <- c('mean', 'sd', 'cor')

  for (i in seq_along(cov_mat_types)) {
    print(cov_mat_types[i])
    for (j in seq_along(change_types)) {
      print(change_types[j])
      multiplot_edd(cov_mat_types[i], change_types[j])
      Sys.sleep(5)
    }
  }
}

test_plot <- function(cov_mat_type, change_type, change_params) {
  multiplot_edd(cov_mat_type, change_type, change_params)
}

test_get_summary_edd <- function(summary_func) {
  change_type <- 'mean'
  n <- 100
  alpha <- 0.01
  m <- 40
  d <- 20
  summary_edd <- get_summary_edd(change_type, n, alpha, m, d, summary_func)
  ggplot_edd(summary_edd, change_type, 0.5)
}
