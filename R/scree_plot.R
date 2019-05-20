scree_plots <- function() {
  n_sets <- 30
  d <- 100
  m <- 200
  eigen_values <- lapply(training_sets, function(train_obj) {
    x <- train_obj$x
    Sigma_hat <- 1 / (ncol(x) - 1) * (x %*% t(x))
    tpca::pca(Sigma_hat)$value
  })
  eigen_values <- do.call('cbind', eigen_values)

  eigen_values_I <- lapply(1:15, function(i) {
    I_hat <- 1 / (m - 1) * rWishart(1, m - 1, diag(rep(1, d)))[, , 1]
    tpca::pca(I_hat)$value
  })
  eigen_values_I <- do.call('cbind', eigen_values_I)
  eigen_values <- cbind(eigen_values, eigen_values_I)

  col <- c(rep('red', 15), rep('darkgreen', 15), rep('black', 15))
  legend_breaks <- c(1, 16, n_sets + 1)
  legend_labels <- c('High correlation', 'Low correlation', 'No correlation')
  line_sizes <- c(rep(0.2, n_sets + 15))
  eigen_value_df <- reshape2::melt(eigen_values,
                                   varnames   = c('component_nr', 'cov_mat_nr'),
                                   value.name = 'eigen_value')
  eigen_value_df$cov_mat_nr <- as.factor(eigen_value_df$cov_mat_nr)
  ggplot2::ggplot(eigen_value_df, ggplot2::aes(x     = component_nr,
                                               y     = eigen_value,
                                               color = cov_mat_nr)) +
    ggplot2::geom_line(ggplot2::aes(size = cov_mat_nr)) +
    ggplot2::theme_light() +
    ggplot2::labs(x = 'Component nr.', y = 'Eigen value') +
    ggplot2::scale_color_manual('',
                                values = col,
                                breaks = legend_breaks,
                                labels = legend_labels) +
    ggplot2::scale_size_manual(values = line_sizes, guide = FALSE)
}

cor_density <- function() {
  n_sets <- 30
  d <- 100
  m <- 200
  cors <- lapply(training_sets, function(train_obj) {
    x <- train_obj$x
    Sigma_hat <- 1 / (ncol(x) - 1) * (x %*% t(x))
    sd_x <- rowSds(x)
    Sigma_hat <- diag(1 / sd_x) %*% Sigma_hat %*% diag(1 / sd_x)
    Sigma_hat[lower.tri(Sigma_hat)]
  })
  cors <- do.call('cbind', cors)

  cors_I <- lapply(1:15, function(i) {
    I_hat <- 1 / (m - 1) * rWishart(1, m - 1, diag(rep(1, d)))[, , 1]
    I_hat[lower.tri(I_hat)]
  })
  cors_I <- do.call('cbind', cors_I)
  cors <- cbind(cors, cors_I)

  col <- c(rep('red', 15), rep('darkgreen', 15), rep('black', 15))
  legend_breaks <- c(1, 16, n_sets + 1)
  legend_labels <- c('High correlation', 'Low correlation', 'No correlation')
  line_sizes <- c(rep(0.2, n_sets + 15))
  cor_df <- reshape2::melt(cors,
                           varnames   = c('cor_nr', 'cov_mat_nr'),
                           value.name = 'cors')
  cor_df$cov_mat_nr <- as.factor(cor_df$cov_mat_nr)
  ggplot2::ggplot(cor_df, ggplot2::aes(cors,
                                       color = cov_mat_nr)) +
    ggplot2::geom_density(ggplot2::aes(size = cov_mat_nr)) +
    ggplot2::theme_light() +
    ggplot2::coord_cartesian(xlim = c(-0.5, 0.5)) +
    ggplot2::labs(x = 'Correlation', y = 'Density') +
    ggplot2::scale_color_manual('',
                                values = col,
                                breaks = legend_breaks,
                                labels = legend_labels) +
    ggplot2::scale_size_manual(values = line_sizes, guide = FALSE)
}

scree_cor_plot <- function() {
  plots <- list()
  plots[[1]] <- scree_plots()
  plots[[2]] <- cor_density()
  ggpubr::ggarrange(plotlist = plots, ncol = 2, nrow = 1,
                    common.legend = TRUE, legend = "bottom")
}

save_scree_cor <- function() {
  ggplot_obj <- scree_cor_plot()
  width <- 7
  height <- 2.8
  dir <- './results/figures/'
  file_name <- 'scree_cor_plot_d100nsets30'
  ggplot2::ggsave(paste0(dir, file_name, '.png'), ggplot_obj,
                  width = width, height = height, units = 'in')
}
