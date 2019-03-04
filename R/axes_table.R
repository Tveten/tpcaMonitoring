get_tpca_axes <- function(cov_mat_type, change_type) {
  axes_location <- paste0('./thresholds/axes/', cov_mat_type, '_tpca_axes.txt')
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

make_axes_table <- function(cov_mat_type) {
  # x: A list of vectors

  begin_table <- paste('\\begin{table}[htb]', '\\caption{Add caption}',
                       '\\label{tab:add_label}', '\\centering',
                       '\\begin{tabular}{l@{\\qquad}l@{\\qquad}l}', '\\toprule',
                       sep = ' \n')
  end_table <- paste('\\bottomrule', '\\end{tabular}', '\\end{table}', sep = ' \n')

  headings <- paste('Change distribution', 'Cutoff',
                    '$\\mathcal{J}$ in order of sensitivity \\\\', sep = ' & ')
  headings <- paste0(headings, '\n\\midrule')

  latex_table <- paste0(begin_table, ' \n', headings)

  change_types <- c('mean', 'sd', 'cor')
  for (change_type in change_types) {
    tpca_axes <- get_tpca_axes(cov_mat_type, change_type)
    axes <- tpca_axes$axes
    cutoffs <- tpca_axes$cutoffs
    for (j in seq_along(cutoffs)) {
      axes_string <- paste(axes[[j]], collapse = ', ')
      table_line <- paste0('& ', cutoffs[j], ' & ', axes_string, ' \\\\')
      if (j == 1)
        table_line <- paste0('\\multirow{6}{*}{', change_type, ' only} ', table_line)
      latex_table <- paste0(latex_table, ' \n', table_line)
      if (j == length(cutoffs) && change_type != 'cor')
        latex_table <- paste0(latex_table, ' \n', '\\midrule')
    }
  }
  latex_table <- paste0(latex_table, ' \n', end_table)
  cat(latex_table)
}
