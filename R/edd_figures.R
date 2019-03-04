get_edd <- function(file) {
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
          if (sum(index) >= 1) edd_array[i, j, k] <- edd_results$edd[index][1]
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

  path <- 'results/'
  edd_results <- read.table(paste0(path, file), header = TRUE)
  p <- sort(unique(edd_results$p))
  change_param <- sort(unique(edd_results$change_param))
  methods <- as.character(unique(edd_results$method))

  edd_list <- list()
  for (i in 1:length(methods)) {
    method_ind <- edd_results$method == methods[i]
    method_params <- sort(unique(edd_results$method_param[method_ind]))
    edd_list[[methods[i]]] <- extract_EDD(methods[i], method.params)
  }
  return(edd_list)
}


##### FIGURES ----------------------------------------------------------------
ggplot_edd <- function(file, change_param, cov_mat_title = FALSE,
                       ylim = c(0, 2000)) {
  # edd.list: An object received from get.EDD
  # change.size.str: For instance: change.size.str = "mu = 0.75".
  # ...: Other subset araguments

  get_change_type <- function(file) {
    strsplit(file, '_')[[1]][2]
  }

  get_cov_mat_type <- function(file) {
    cov_mat_type <- strsplit(file, '_')[[1]][3]
    if (cov_mat_type == 'halfsparse') cov_mat_type <- 'half-sparse'
    cov_mat_type
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

  set_title <- function(change_type, cov_mat_type) {
    title_str <- ''
    if (cov_mat_title)
      title_str <- paste0(title_str, first_up(cov_mat_type), ' $\\mathbf{\\Sigma}_0$, ')
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

  edd_list <- get_edd(file)
  change_type <- get_change_type(file)
  cov_mat_type <- get_cov_mat_type(file)
  methods <- names(edd_list)
  col <- set_col(methods)
  labels <- set_labels(methods)
  title <- set_title(change_type, cov_mat_type)

  edd_df <- reshape2::melt(edd_list, varnames = c('p', 'change_size', 'method_param'),
                           value.name = 'edd', level = 'method')
  if (change_type == 'cor') edd_df <- subset(edd_df, Lmethod != 'mix')
  edd_df$p <- extract_double(edd_df$p)
  edd_df$change_size <- extract_double(edd_df$change_size)
  edd_df$method_param <- as.factor(extract_double(edd_df$method_param))
  edd_subset <- subset(edd_df, change_size == change_param)

  ggplot2::ggplot(edd_subset, ggplot2::aes(x = p, y = edd, color = Lmethod, linetype = method_param)) +
    ggplot2::geom_line() +
    ggplot2::theme_light() +
    ggplot2::labs(title = title, x = 'p', y = 'EDD') +
    ggplot2::coord_cartesian(ylim = ylim) +
    ggplot2::scale_color_manual('', values = col, labels = labels) +
    ggplot2::scale_linetype_manual(values = rep(1, length(unique(edd_df$method_param))), guide = FALSE)
    # scale_linetype_manual('Method parameters', values = c(1, 5, 2, 4, 1, 5, 2, 4, 3),
    #                       labels = c('p0 = 0.03',
    #                                  'p0 = 0.1',
    #                                  'p0 = 0.3',
    #                                  'p0 = 1',
    #                                  'M = 2',
    #                                  'M = 3',
    #                                  'M = 5',
    #                                  'M = 10',
    #                                  'M = 20'))
}

multiplot_edd <- function(cov_mat_type, change_type,
                          change_params = NULL, ylim = c(0, 2000)) {
  path <- 'results/'
  file <- paste0('edd_', change_type, '_', cov_mat_type, '_m300d100w200.txt')
  if (is.null(change_params)) {
    edd_results <- read.table(paste0(path, file), header = TRUE)
    change_params <- sort(unique(edd_results$change_param))
  }
  l <- length(change_params)
  ggplots <- lapply(change_params, function(x) ggplot_edd(file, x, ylim))

  invisible(gridExtra::grid.arrange(grobs = ggplots,
                                    nrow = ceiling(l / 2),
                                    ncol = min(l, 2)))
}

multiplot_edd_cov_mats <- function(change_type, change_params = NULL,
                              ylim = c(0, 2000)) {
  path <- 'results/'
  cov_mat_types <- c('dense', 'halfsparse', 'sparse')
  files <- paste0('edd_', change_type, '_', cov_mat_types, '_m300d100w200.txt')
  n <- length(files)
  plots <- list()
  for (i in seq_along(files)) {
    if (is.null(change_params)) {
      edd_results <- read.table(paste0(path, files[i]), header = TRUE)
      change_params <- sort(unique(edd_results$change_param))
    }
    l <- length(change_params)
    # plots <- c(plots, lapply(change_params, function(x) ggplot_edd(files[i], x, ylim)))
    plots_temp <- lapply(change_params, function(x) ggplot_edd(files[i], x, ylim))

    # ggarrange fills up row by row, so some rearrangement is needed.
    for (j in seq_along(change_params)) {
      plots[[(j - 1) * n + i]] <- plots_temp[[j]]
    }
  }

  ggpubr::ggarrange(plotlist = plots, ncol = n, nrow = ceiling(length(plots) / n),
                    common.legend = TRUE, legend = "bottom")
}

multiplot_edd_summary <- function(cov_mat_type, ylim = c(0, 2000)) {
  path <- 'results/'
  change_types <- c('mean', 'sd', 'cor')
  change_params <- list('mean' = c(0.5, 0.7, 1),
                        'sd' = c(0.5, 1.5, 2),
                        'cor' = c(0, 0.5, 0.75))
  files <- paste0('edd_', change_types, '_', cov_mat_type, '_m300d100w200.txt')
  n <- length(files)
  plots <- list()
  for (i in seq_along(files)) {
    edd_results <- read.table(paste0(path, files[i]), header = TRUE)
    l <- length(change_params)
    plots[[i]] <- lapply(change_params[[i]], function(param) {
        ggplot_edd(files[i], param, ylim = ylim)
      })
  }
  plots <- unlist(plots, recursive = FALSE)

  ggpubr::ggarrange(plotlist = plots, ncol = n, nrow = ceiling(length(plots) / n),
                    common.legend = TRUE, legend = "bottom")
}

save_all_plots <- function(ylim = c(0, 2000)) {
  cov_mat_types <- c('dense', 'halfsparse', 'sparse')
  change_types <- c('mean', 'sd', 'cor')
  for (i in seq_along(cov_mat_types)) {
    for (j in seq_along(change_types)) {
      fig <- multiplot_edd(cov_mat_types[i], change_types[j], ylim = ylim)
      file_name <- paste0('edd_', change_types[j], '_', cov_mat_types[i],
                          '_ylim', ylim[1], '-', ylim[2], '.png')
      ggplot2::ggsave(file_name, plot = fig, path = './results/figures/',
                      width = 20, height = 14, units = 'cm')
    }
  }

}

save_change_type_plots <- function(change_type, extension = 'png',
                                   ylim = c(0, 200)) {
  base_width <- 2.88
  base_height <- 2.3
  n_col <- 3
  n_row <- 4
  width <- n_col * base_width
  height <- n_row * base_height
  plot_obj <- multiplot_edd_all(change_type, ylim = ylim)
  show(plot_obj)
  dir <- './results/figures/'
  file_name <- paste0('edd_', change_type, '_ylim', ylim[1], '-', ylim[2])
  ggplot2::ggsave(paste0(dir, file_name, '.', extension),
                  plot = plot_obj,
                  width = width, height = height, units = 'in')
}

save_all_change_type_plots <- function(ylim = c(0, 2000), extension = 'png') {
  change_types <- c('mean', 'sd', 'cor')
  lapply(change_types, save_change_type_plots, ylim = ylim, extension = extension)
}

save_summary_plot <- function(cov_mat_type = 'dense', ylim = c(0, 200), extension = 'png') {
  base_width <- 2.88
  base_height <- 2.9
  n_col <- 3
  n_row <- 3
  width <- n_col * base_width
  height <- n_row * base_height
  plot_obj <- multiplot_edd_summary(cov_mat_type, ylim = ylim)
  show(plot_obj)
  dir <- './results/figures/'
  file_name <- paste0('edd_', cov_mat_type, '_ylim', ylim[1], '-', ylim[2])
  ggplot2::ggsave(paste0(dir, file_name, '.', extension),
                  plot = plot_obj,
                  width = width, height = height, units = 'in')
}

##### TABLES -----------------------------------------------------------------
best_method_param <- function(cov_mat_type, change_type) {
  file <- paste0('edd_', change_type, '_', cov_mat_type, '_m300d100w200.txt')

  # Returns list with 3-array (p, size, method_param) for each method.
  edd <- get_edd(file)
  best_method_param <- list()
  method <- names(edd)
  for (i in seq_along(edd)) {
    best_method_param[[method[i]]] <- apply(edd[[method[i]]], 1:2, function(x) {
      dimnames(edd[[method[i]]])[[3]][which.min(x)]
    })
  }
  best_method_param
}

summarize_edd <- function(change_type, cov_mat_type, summary_func, ...) {
  file <- paste0('edd_', change_type, '_', cov_mat_type, '_m300d100w200.txt')
  edd <- get_edd(file)
  method <- c('mix', 'max_pca', 'min_pca', 'tpca')
  n_params <- list(4, 6, 6, 6)
  names(n_params) <- method
  edd_summary <- list()
  for (i in seq_along(method)) {
    if (!is.null(edd[[method[i]]])) {
      edd_summary[[method[i]]] <- apply(edd[[method[i]]], 3, summary_func, ...)
    } else {
      edd_summary[[method[i]]] <- rep(NA, n_params[[method[i]]])
    }
  }
  edd_summary
}

mean_edd_type <- function(change_type, show = FALSE) {
  # Averages over cov_mat_types

  cov_mat_types <- c('dense', 'halfsparse', 'sparse')
  mean_edd <- summarize_edd(change_type, cov_mat_types[1], mean)
  method <- names(mean_edd)
  for (i in 2:length(cov_mat_types)) {
    next_edd_mean <- summarize_edd(change_type, cov_mat_types[i], mean)
    for (j in seq_along(method)) {
      mean_edd[[method[j]]] <- rbind(mean_edd[[method[j]]], next_edd_mean[[method[j]]])
    }
  }
  for (j in seq_along(mean_edd)) {
    rownames(mean_edd[[j]]) <- cov_mat_types
  }
  if (show) print(mean_edd)
  mean_edd_type <- lapply(mean_edd, function(edd_mat) apply(edd_mat, 2, mean))
  mean_edd_type
}

min_edd_summary <- function(cov_mat_type, summary_func, show = FALSE, ...) {
  # Summarizes EDD for each change_type for a given cov_mat_type.
  # Picks out the minimizer among method_parameters.
  # TODO: Combined table with all cov_mat_types separately, with method_param in parenthesis.

  which_method_approx_min <- function(y) {
    # Finds the minimizing arguments of y within a tolerance of 0.5.
    # Here: Minimum of average EDDs.
    tol <- 0.5
    if (length(which.min(y)) != 0) {
      ind <- which(abs(y - min(y)) < tol)
      params_temp <- names(y)[ind]
      params_temp <- strsplit(params_temp, ' = ')
      params <- vapply(params_temp, `[`, character(1), 2)
      min_char <- paste0('method_param = ', paste(params, collapse = ', '))
      return(min_char)
    } else return(NA)
  }

  change_types <- c('mean', 'sd', 'cor')
  if (cov_mat_type == 'all')
    mean_edd <- lapply(change_types, mean_edd_type)
  else
    mean_edd <- lapply(change_types, function(type) summarize_edd(type, cov_mat_type, summary_func, ...))
  names(mean_edd) <- change_types
  if (show) print(mean_edd)

  min_mean_edd <- lapply(mean_edd, function(x) as.data.frame(lapply(x, min)))
  min_mean_edd <- do.call('rbind', min_mean_edd)

  which_min_mean_edd <- lapply(mean_edd, function(x) {
    as.data.frame(lapply(x, which_method_approx_min), stringsAsFactors = FALSE)})
  which_min_mean_edd <- do.call('rbind', which_min_mean_edd)
  return(list('min' = min_mean_edd, 'which_min' = which_min_mean_edd))
}

edd_summary_table <- function(summary_func, ...) {
  caption <- 'Average EDD per type of pre-change covariance matrix and change type for each method\'s best method parameter (in parenthesis). The average is taken over change sparsity and change size.'
  label <- 'tab:EDD_summary'
  date_line <- paste0('% ', date())
  begin_table <- paste('\\begin{table}[htb]',
                       paste0('\\caption{', caption, '}'),
                       paste0('\\label{', label ,'}'),
                       '\\centering',
                       '\\begin{tabular}{llllll}',
                       '\\toprule', sep = ' \n')
  end_table <- paste('\\bottomrule',
                     '\\end{tabular}',
                     '\\end{table}', sep = ' \n')
  headings1 <- paste('\\multirow{2}{*}{Sparsity of $\\bSigma_0$}',
                    '\\multirow{2}{*}{Change type}',
                    '\\multicolumn{4}{c}{EDD} \\\\', sep = ' & ')
  headings1 <- paste0(headings1, '\n\\cmidrule{3-6}')
  headings2 <- ' & & Mixture($p_0$) & Max PCA($J$) & Min PCA($J$) & TPCA$(c)$ \\\\'
  mid_sep <- ' \n\\midrule'
  latex_table <- paste(date_line, begin_table, headings1, headings2, sep = ' \n')

  shown_type_names <- list('mean' = 'Mean', 'sd' = 'Var', 'cor' = 'Cor')
  cov_mat_types <- c('dense', 'halfsparse', 'sparse')
  for (l in seq_along(cov_mat_types)) {
    min_edd_list <- min_edd_summary(cov_mat_types[l], summary_func, ...)
    min_edd <- min_edd_list$min
    which_min_edd <- min_edd_list$which_min
    latex_table <- paste0(latex_table, mid_sep)
    for (i in 1:nrow(min_edd)) {
      type <- shown_type_names[[rownames(min_edd)[i]]]
      table_line <- paste0(' & ', type)
      if (i == 1)
        table_line <- paste0('\\multirow{3}{*}{', cov_mat_types[l],'}', table_line)
      for (j in 1:ncol(min_edd)) {
        edd <- round(min_edd[i, j], 1)
        which_param <- unlist(strsplit(which_min_edd[i, j], ' = '))[2]
        table_line <- paste0(table_line, ' & ', edd, ' (', which_param, ')')
        if (j == ncol(min_edd))
          table_line <- paste0(table_line, ' \\\\')
      }
      latex_table <- paste(latex_table, table_line, sep = ' \n')
    }
  }
  latex_table <- paste0(latex_table, ' \n', end_table)
  cat(latex_table)
}

edd_summary_table_dense <- function(summary_func, ...) {
  caption <- 'Average EDD per change type for each method\'s best method parameters (in parenthesis). '
  caption <- paste0(caption, 'Each listed method parameter is within 0.5 of the method\'s minimum average EDD.')
  caption <- paste0(caption, 'The average is taken over change sparsity and change size.')
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
  headings2 <- ' & Mixture($p_0$) & Max PCA($J$) & Min PCA($J$) & TPCA$(c)$ \\\\'
  mid_sep <- ' \n\\midrule'
  latex_table <- paste(date_line, begin_table, headings1, headings2, sep = ' \n')

  cov_mat_type <- 'dense'
  shown_type_names <- list('mean' = 'Mean', 'sd' = 'Variance', 'cor' = 'Correlation')
  min_edd_list <- min_edd_summary(cov_mat_type, summary_func, ...)
  min_edd <- min_edd_list$min
  which_min_edd <- min_edd_list$which_min
  latex_table <- paste0(latex_table, mid_sep)
  for (i in 1:nrow(min_edd)) {
    type <- shown_type_names[[rownames(min_edd)[i]]]
    table_line <- type
    # if (i == 1)
    #   table_line <- paste0('\\multirow{3}{*}{', cov_mat_type,'}', table_line)
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

