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

  edd_results <- read.table(paste0('./results/', file), header = TRUE)
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

ggplot_edd <- function(file, change_param,
                       ylim = c(0, 2000)) {
  # edd.list: An object received from get.EDD
  # change.size.str: For instance: change.size.str = "mu = 0.75".
  # ...: Other subset araguments

  get_change_type <- function(file) {
    strsplit(file, '_')[[1]][2]
  }

  get_all_methods <- function() {
    c('mix', 'max_pca', 'min_pca', 'tpca')
  }

  set_col <- function(methods) {
    all_col <- list('darkgreen', 'darkgoldenrod', 'darkblue', 'red')
    names(all_col) <- get_all_methods()
    unlist(all_col[methods])
  }

  set_labels <- function(methods) {
    all_labels <- list('Raw data', 'Max PCA', 'Min PCA', 'TPCA')
    names(all_labels) <- get_all_methods()
    unlist(all_labels[methods])
  }

  set_title <- function(change_type) {
    if (change_type == 'mean') {
      title_str <- paste0('Change in mean: 0 -> ', change_param)
    }
    if (change_type == 'sd') {
      title_str <- paste0('Change in sd: 1 -> ', change_param)
    }
    if (change_type == 'cor') {
      rho_char <- '\u03C1'
      title_str <- paste0('Change in cor: ', rho_char, ' -> ', change_param)
    }
    return(title_str)
  }

  extract_double <- function(x, split = '= ') {
    x_temp <- strsplit(as.character(x), split)
    vapply(x_temp, function(x) as.double(x[2]), numeric(1))
  }

  edd_list <- get_edd(file)
  change_type <- get_change_type(file)
  methods <- names(edd_list)
  col <- set_col(methods)
  labels <- set_labels(methods)
  title <- set_title(change_type)

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
  file <- paste0('./results/edd_', change_type, '_', cov_mat_type, '_m300d100w200.txt')
  if (is.null(change_params)) {
    edd_results <- read.table(file, header = TRUE)
    change_params <- sort(unique(edd_results$change_param))
  }
  l <- length(change_params)
  ggplots <- lapply(change_params, function(x) ggplot_edd(file, x, ylim))

  invisible(gridExtra::grid.arrange(grobs = ggplots,
                                    nrow = ceiling(l / 2),
                                    ncol = min(l, 2)))
}

test_plot <- function(cov_mat_type, change_type, change_params) {
  multiplot_edd(cov_mat_type, change_type, change_params)
}

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

save_all_plots <- function(ylim = c(0, 2000)) {
  cov_mat_types <- c('dense', 'halfsparse', 'sparse')
  change_types <- c('mean', 'sd', 'cor')
  for (i in seq_along(cov_mat_types)) {
    for (j in seq_along(change_types)) {
      fig <- multiplot_edd(cov_mat_types[i], change_types[j], ylim = ylim)
      file_name <- paste0('edd_', change_types[j], '_', cov_mat_types[i],
                          '_ylim', ylim[1], '-', ylim[2], '.png')
      ggplot2::ggsave(file_name, plot = fig, path = './results/',
                      width = 20, height = 14, units = 'cm')
    }
  }

}
