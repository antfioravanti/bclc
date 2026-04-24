# Needed packages:
library(dplyr)
library(ggplot2)
library(scales)
library(viridisLite)

# Helper to make an "N" tag for filenames
.n_tag <- function(N_vals) {
  N_vals <- sort(unique(N_vals))
  paste(N_vals, collapse = "-")
}

# Create a palette function with your specific breakpoints
color_vals <- rescale(c(0, 0.5, 0.75, 0.85, 0.95, 1, 1.05, 1.15, 1.25, 1.5, 1.75, 2))

palette_colors <- c(
  "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#F7F7F7",
  "white",
  "#FEE0D2", "#FCBBA1", "#FC9272", "#FB6A4A", "#EF3B2C", "#CB181D"
)

# Create palette function
my_palette <- gradient_n_pal(palette_colors, values = color_vals)

#-----------------------------
# 1) Largest m with MSEReduced
#-----------------------------
plot_max_m_reduced <- function(summary_df,
                               lambda_filter,
                               title_bool = TRUE,          
                               out_dir = ".",
                               filename = NULL,            # auto if NULL
                               width = 8, height = 6, dpi = 300,
                               label_size = 3,
                               palette_option = "plasma",
                               palette_direction = -1) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  # helper for N tag (in case it's not defined globally)
  .n_tag <- function(Nvec) if (length(Nvec)) paste(sort(unique(Nvec)), collapse = "-") else "NA"
  
  df_reduced <- summary_df %>%
    dplyr::filter(MSEReduced, lambda == lambda_filter)
  
  # Build default filename using unique N in filtered data
  if (is.null(filename)) {
    n_tag <- .n_tag(df_reduced$N)
    if(title_bool){
      filename <- sprintf("largest_m_MSEReduced_lambda%s_N%s.png", as.character(lambda_filter), n_tag)
    }else{
      filename <- sprintf("non_largest_m_MSEReduced_lambda%s_N%s.png", as.character(lambda_filter), n_tag)
    }
  }
  
  df_max <- df_reduced %>%
    dplyr::group_by(h1, h2) %>%
    dplyr::summarise(max_m = max(m, na.rm = TRUE), .groups = "drop")
  
  pal_fun <- scales::viridis_pal(option = palette_option, direction = palette_direction)
  
  df_max_mutated <- df_max %>%
    dplyr::mutate(
      scale_pos = scales::rescale(max_m, to = c(0, 1), from = range(max_m, na.rm = TRUE)),
      idx       = pmin(101, pmax(1, ceiling(scale_pos * 100) + 1)),
      fill_hex  = pal_fun(101)[idx],
      luminance = {
        rgb_mat <- as.data.frame(t(col2rgb(fill_hex)))
        (0.2126 * rgb_mat$red + 0.7152 * rgb_mat$green + 0.0722 * rgb_mat$blue) / 255
      },
      text_col = ifelse(luminance < 0.5, "white", "black")
    )
  
  # labs layer depends on title_bool
  labs_layer <- if (isTRUE(title_bool)) {
    ggplot2::labs(
      title = paste("Largest m with MSEReduced = TRUE, lambda =", lambda_filter),
      x = "h1", y = "h2"
      # add subtitle here if you ever want one
    )
  } else {
    ggplot2::labs(x = "h1", y = "h2")
  }
  
  p <- ggplot2::ggplot(df_max_mutated, ggplot2::aes(x = factor(h1), y = factor(h2), fill = max_m)) +
    ggplot2::geom_tile(colour = "white") +
    ggplot2::geom_text(ggplot2::aes(label = max_m, colour = text_col), size = label_size) +
    ggplot2::scale_colour_identity() +
    ggplot2::scale_fill_viridis_c(option = palette_option, direction = palette_direction, name = "Max m") +
    labs_layer +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), legend.position = "none")
  
  out_path <- file.path(out_dir, filename)
  ggplot2::ggsave(out_path, plot = p, width = width, height = height, dpi = dpi, bg = "white")
  
  if (interactive()) print(p)
  return(p)
}

#-------------------------------------------------------------
# 2) Best m = argmin(MSEChatCorr) among rows with MSEReduced=1
#-------------------------------------------------------------
plot_best_m_min_MSE <- function(summary_df,
                                lambda_filter,
                                title_bool = TRUE,          # <— added
                                out_dir = ".",
                                filename = NULL,            # auto if NULL
                                width = 8, height = 6, dpi = 300,
                                label_size = 3,
                                palette_option = "plasma",
                                palette_direction = -1) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  # helper for N tag (in case not defined globally)
  .n_tag <- function(Nvec) if (length(Nvec)) paste(sort(unique(Nvec)), collapse = "-") else "NA"
  
  df_reduced <- summary_df %>%
    dplyr::filter(MSEReduced, lambda == lambda_filter)
  
  # Build default filename using unique N in filtered data
  if (is.null(filename)) {
    n_tag <- .n_tag(df_reduced$N)
    if(title_bool){
      filename <- sprintf("best_m_minMSEChatCorr_MSEReduced_lambda%s_N%s.png",
                          as.character(lambda_filter), n_tag)
    }else{
      filename <- sprintf("non_best_m_minMSEChatCorr_MSEReduced_lambda%s_N%s.png",
                          as.character(lambda_filter), n_tag)
    }
  }
  
  df_min <- df_reduced %>%
    dplyr::group_by(h1, h2) %>%
    dplyr::slice_min(order_by = MSEChatCorr, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::select(h1, h2, m, MSEChatCorr) %>%
    dplyr::rename(best_m = m, best_MSE = MSEChatCorr)
  
  pal_fun <- scales::viridis_pal(option = palette_option, direction = palette_direction)
  
  df_min_mutated <- df_min %>%
    dplyr::mutate(
      scale_pos = scales::rescale(best_m, to = c(0, 1), from = range(best_m, na.rm = TRUE)),
      idx       = pmin(101, pmax(1, ceiling(scale_pos * 100) + 1)),
      fill_hex  = pal_fun(101)[idx],
      luminance = {
        rgb_mat <- as.data.frame(t(col2rgb(fill_hex)))
        (0.2126 * rgb_mat$red + 0.7152 * rgb_mat$green + 0.0722 * rgb_mat$blue) / 255
      },
      text_col = ifelse(luminance < 0.5, "white", "black")
    )
  
  # labs layer controls title/subtitle
  labs_layer <- if (isTRUE(title_bool)) {
    ggplot2::labs(
      title = paste("Best m (min MSEChatCorr) with MSEReduced = TRUE, lambda =", lambda_filter),
      x = "h1", y = "h2"
    )
  } else {
    ggplot2::labs(x = "h1", y = "h2")
  }
  
  p <- ggplot2::ggplot(df_min_mutated, ggplot2::aes(x = factor(h1), y = factor(h2), fill = best_m)) +
    ggplot2::geom_tile(colour = "white") +
    ggplot2::geom_text(ggplot2::aes(label = best_m, colour = text_col), size = label_size) +
    ggplot2::scale_colour_identity() +
    ggplot2::scale_fill_viridis_c(option = palette_option, direction = palette_direction, name = "Best m") +
    labs_layer +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), legend.position = "none")
  
  out_path <- file.path(out_dir, filename)
  ggplot2::ggsave(out_path, plot = p, width = width, height = height, dpi = dpi, bg = "white")
  
  if (interactive()) print(p)
  return(p)
}


#-------------------------------------------------------------
# 2) Best m = argmin(MSEChatCorr) among rows with MSEReduced=1 and MSEZred=0
#-------------------------------------------------------------
# plot_best_m_min_MSE_flags <- function(summary_df,
#                                       lambda_filter,
#                                       title_bool = TRUE,
#                                       out_dir = ".",
#                                       filename = NULL,            # auto if NULL
#                                       width = 8, height = 6, dpi = 300,
#                                       label_size = 3,
#                                       palette_option = "plasma",
#                                       palette_direction = -1,
#                                       msered_col = "MSEZRed",      # second-criterion column
#                                       show_legend = TRUE) {
#   
#   n_i = unique(summary_df$N)
#   .n_tag <- function(Nvec) if (length(Nvec)) paste(sort(unique(Nvec)), collapse = "-") else "NA"
#   if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
#   
#   df_lam <- dplyr::filter(summary_df, .data$lambda == !!lambda_filter)
#   if (nrow(df_lam) == 0) stop("No rows for the requested lambda_filter.")
#   
#   if (is.null(filename)) {
#     n_tag <- .n_tag(df_lam$N)
#     if(title_bool){
#       filename <- sprintf("bestm_minMSE_flags_lambda%s_N%s.png",
#                           as.character(lambda_filter), n_tag)
#     }else{
#       filename <- sprintf("non_bestm_minMSE_flags_lambda%s_N%s.png",
#                           as.character(lambda_filter), n_tag)
#     }
#   }
#   out_path <- file.path(out_dir, filename)
#   
#   if (!msered_col %in% names(df_lam)) {
#     warning(sprintf("Column '%s' not found; X overlay will be skipped.", msered_col))
#     df_lam[[msered_col]] <- NA
#   }
#   
#   grid_hh <- df_lam %>% dplyr::distinct(h1, h2)
#   df_reduced <- df_lam %>% dplyr::filter(.data$MSEReduced)
#   
#   df_best <- df_reduced %>%
#     dplyr::group_by(h1, h2) %>%
#     dplyr::arrange(.data$MSEChatCorr, .data$m, .by_group = TRUE) %>%
#     dplyr::slice(1) %>%
#     dplyr::ungroup() %>%
#     dplyr::transmute(
#       h1, h2,
#       best_m = .data$m,
#       best_MSE = .data$MSEChatCorr,
#       zero_reduced = .data[[msered_col]]
#     )
#   
#   df_plot <- dplyr::left_join(grid_hh, df_best, by = c("h1","h2"))
#   
#   pal_fun <- scales::viridis_pal(option = palette_option, direction = palette_direction)
#   
#   if (any(!is.na(df_plot$best_m))) {
#     rng <- range(df_plot$best_m, na.rm = TRUE)
#     df_plot <- df_plot %>%
#       dplyr::mutate(
#         scale_pos = ifelse(is.na(.data$best_m), NA_real_,
#                            scales::rescale(.data$best_m, to = c(0,1), from = rng)),
#         idx       = ifelse(is.na(.data$scale_pos), NA_integer_,
#                            pmin(101L, pmax(1L, as.integer(ceiling(.data$scale_pos * 100) + 1)))),
#         fill_hex  = dplyr::if_else(is.na(.data$idx), NA_character_, pal_fun(101)[.data$idx]),
#         luminance = dplyr::if_else(
#           is.na(.data$fill_hex), NA_real_,
#           {
#             rgb_mat <- as.data.frame(t(col2rgb(.data$fill_hex)))
#             (0.2126 * rgb_mat$red + 0.7152 * rgb_mat$green + 0.0722 * rgb_mat$blue) / 255
#           }
#         ),
#         text_col  = dplyr::if_else(is.na(.data$luminance), "black",
#                                    dplyr::if_else(.data$luminance < 0.5, "white", "black")),
#         label     = dplyr::if_else(is.na(.data$best_m), "", as.character(.data$best_m))
#       )
#   } else {
#     df_plot <- df_plot %>%
#       dplyr::mutate(fill_hex = NA_character_, text_col = "black", label = "")
#   }
#   
#   df_plot <- df_plot %>%
#     dplyr::mutate(
#       x = as.numeric(factor(.data$h1, levels = sort(unique(.data$h1)))),
#       y = as.numeric(factor(.data$h2, levels = sort(unique(.data$h2))))
#     )
#   
#   # --- minimal change: build labs layer depending on title_bool ---
#   labs_layer <- if (isTRUE(title_bool)) {
#     ggplot2::labs(
#       title = paste0("m with smallest MSE of Corrected Estimator with \u03BB = ",
#                      lambda_filter, ", N = ", n_i, " \u00D7 ", n_i),
#       subtitle = "Grey = No MSE Reduction at (h1,h2); X = Zero Estimation has smaller MSE",
#       x = "h1", y = "h2"
#     )
#   } else {
#     ggplot2::labs(x = "h1", y = "h2")
#   }
#   
#   p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = factor(h1), y = factor(h2), fill = best_m)) +
#     ggplot2::geom_tile(color = "white", linewidth = 0.3, na.rm = FALSE) +
#     ggplot2::geom_text(ggplot2::aes(label = label, colour = I(text_col)), size = label_size) +
#     ggplot2::geom_segment(
#       data = df_plot %>% dplyr::filter(!is.na(best_m), !is.na(zero_reduced), zero_reduced == FALSE),
#       ggplot2::aes(x = x - 0.45, xend = x + 0.45, y = y - 0.45, yend = y + 0.45),
#       inherit.aes = FALSE, linewidth = 0.6
#     ) +
#     ggplot2::geom_segment(
#       data = df_plot %>% dplyr::filter(!is.na(best_m), !is.na(zero_reduced), zero_reduced == FALSE),
#       ggplot2::aes(x = x - 0.45, xend = x + 0.45, y = y + 0.45, yend = y - 0.45),
#       inherit.aes = FALSE, linewidth = 0.6
#     ) +
#     ggplot2::scale_fill_viridis_c(
#       option = palette_option,
#       direction = palette_direction,
#       name = "Best m",
#       na.value = "white"
#     ) +
#     labs_layer +
#     ggplot2::theme_minimal() +
#     ggplot2::theme(
#       panel.grid = ggplot2::element_blank(),
#       legend.position = if (isTRUE(show_legend)) "right" else "none"
#     )
#   
#   # Extra safety: remove legend when show_legend = FALSE
#   if (!isTRUE(show_legend)) {
#     p <- p + ggplot2::guides(fill = "none")
#   }
#   
#   ggplot2::ggsave(out_path, plot = p, width = width, height = height, dpi = dpi, bg = "white")
#   if (interactive()) print(p)
#   invisible(list(plot = p, file = out_path))
# }
#-------------------------------------------------------------
# 4) Bias vs Lag Norm for different m values
#-------------------------------------------------------------
plot_bias_vs_norm <- function(summary_df,
                              lambda_filter,
                              N_filter,
                              title_bool = TRUE,
                              out_dir = ".",
                              filename = NULL,
                              width = 10, height = 6, dpi = 300,
                              palette_option = "plasma",
                              show_legend = TRUE,
                              line_size = 0.8,
                              point_size = 2,
                              benchmark_color = "black",
                              benchmark_size = 0.9,
                              benchmark_linetype = "dashed") {
  
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  # Filter data
  df_filtered <- summary_df %>%
    dplyr::filter(lambda == lambda_filter, N == N_filter)
  
  if (nrow(df_filtered) == 0) {
    stop("No data found for the specified lambda and N filters.")
  }
  
  # Build default filename
  if (is.null(filename)) {
    if (title_bool) {
      filename <- sprintf("bias_vs_norm_lambda%s_N%s.png", 
                          as.character(lambda_filter), N_filter)
    } else {
      filename <- sprintf("non_bias_vs_norm_lambda%s_N%s.png", 
                          as.character(lambda_filter), N_filter)
    }
  }
  
  # Calculate lag norm and aggregate by unique norm values
  df_plot <- df_filtered %>%
    dplyr::mutate(
      lag_norm = sqrt(h1^2 + h2^2)
    ) %>%
    dplyr::group_by(m, lag_norm) %>%
    dplyr::summarise(
      bias = mean(bias, na.rm = TRUE),
      biasCorr = mean(biasCorr, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Benchmark: uncorrected bias (bias column, same for all m at each lag_norm)
  df_benchmark <- df_plot %>%
    dplyr::group_by(lag_norm) %>%
    dplyr::summarise(
      benchmark_bias = first(bias),
      .groups = "drop"
    ) %>%
    dplyr::mutate(estimator_type = "Benchmark")
  
  # Data for corrected estimators (all m values with biasCorr)
  df_corrected <- df_plot %>%
    dplyr::mutate(estimator_type = paste0("m = ", m))
  
  # Create color palette for m values
  m_values <- sort(unique(df_corrected$m))
  n_colors <- length(m_values)
  
  # Check if palette_option is a vector of colors or a viridis name
  if (length(palette_option) > 1) {
    # Custom palette provided
    colors <- colorRampPalette(palette_option)(n_colors)
  } else {
    # Viridis palette name provided
    colors <- viridisLite::viridis(n_colors, option = palette_option)
  }
  names(colors) <- paste0("m = ", m_values)
  
  # Add benchmark color
  all_colors <- c("Benchmark" = benchmark_color, colors)
  
  # Build title/subtitle layer
  labs_layer <- if (isTRUE(title_bool)) {
    ggplot2::labs(
      title = sprintf("Bias vs Lag Norm: λ = %s, N = %d × %d", 
                      lambda_filter, N_filter, N_filter),
      x = "Lag Norm ||h||",
      y = "Bias",
      color = "Estimator"
    )
  } else {
    ggplot2::labs(
      x = "Lag Norm ||h||",
      y = "Bias",
      color = "Estimator"
    )
  }
  
  # Create plot
  p <- ggplot2::ggplot() +
    # Benchmark line (uncorrected bias from 'bias' column)
    ggplot2::geom_line(
      data = df_benchmark,
      ggplot2::aes(x = lag_norm, y = benchmark_bias, color = estimator_type),
      linetype = benchmark_linetype,
      linewidth = benchmark_size
    ) +
    ggplot2::geom_point(
      data = df_benchmark,
      ggplot2::aes(x = lag_norm, y = benchmark_bias, color = estimator_type),
      size = point_size
    ) +
    # Corrected estimator lines (biasCorr for all m)
    ggplot2::geom_line(
      data = df_corrected,
      ggplot2::aes(x = lag_norm, y = biasCorr, color = estimator_type, group = m),
      linewidth = line_size
    ) +
    ggplot2::geom_point(
      data = df_corrected,
      ggplot2::aes(x = lag_norm, y = biasCorr, color = estimator_type),
      size = point_size
    ) +
    ggplot2::scale_color_manual(
      values = all_colors,
      name = "Estimator",
      breaks = c("Benchmark", paste0("m = ", m_values))
    ) +
    labs_layer +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = if (isTRUE(show_legend)) "right" else "none",
      panel.grid.minor = ggplot2::element_blank()
    )
  
  # Add horizontal line at y = 0
  p <- p + ggplot2::geom_hline(yintercept = 0, linetype = "dotted", 
                               color = "gray40", linewidth = 0.5)
  
  # Save plot
  out_path <- file.path(out_dir, filename)
  ggplot2::ggsave(out_path, plot = p, width = width, height = height, 
                  dpi = dpi, bg = "white")
  
  if (interactive()) print(p)
  
  invisible(list(plot = p, file = out_path, 
                 data = list(corrected = df_corrected, benchmark = df_benchmark)))
}
#-------------------------------------------------------------------------------
plot_mse_vs_norm <- function(summary_df,
                             lambda_filter,
                             N_filter,
                             title_bool = TRUE,
                             out_dir = ".",
                             filename = NULL,
                             width = 10, height = 6, dpi = 300,
                             palette_option = "plasma",
                             show_legend = TRUE,
                             line_size = 0.8,
                             point_size = 2,
                             benchmark_color = "red",
                             benchmark_size = 1.2,
                             benchmark_linetype = "dashed") {
  
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  # Filter data
  df_filtered <- summary_df %>%
    dplyr::filter(lambda == lambda_filter, N == N_filter)
  
  if (nrow(df_filtered) == 0) {
    stop("No data found for the specified lambda and N filters.")
  }
  
  # Build default filename
  if (is.null(filename)) {
    if (title_bool) {
      filename <- sprintf("mse_vs_norm_lambda%s_N%s.png", 
                          as.character(lambda_filter), N_filter)
    } else {
      filename <- sprintf("non_mse_vs_norm_lambda%s_N%s.png", 
                          as.character(lambda_filter), N_filter)
    }
  }
  
  # Calculate lag norm and aggregate by unique norm values
  df_plot <- df_filtered %>%
    dplyr::mutate(
      lag_norm = sqrt(h1^2 + h2^2)
    ) %>%
    dplyr::group_by(m, lag_norm) %>%
    dplyr::summarise(
      MSEChat = mean(MSEChat, na.rm = TRUE),
      MSEChatCorr = mean(MSEChatCorr, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Benchmark: uncorrected MSE (MSEChat column, same for all m at each lag_norm)
  df_benchmark <- df_plot %>%
    dplyr::group_by(lag_norm) %>%
    dplyr::summarise(
      benchmark_mse = first(MSEChat),
      .groups = "drop"
    ) %>%
    dplyr::mutate(estimator_type = "Naive")
  
  # Data for corrected estimators (all m values with MSEChatCorr)
  df_corrected <- df_plot %>%
    dplyr::mutate(estimator_type = paste0("m = ", m))
  
  # Create color palette for m values
  m_values <- sort(unique(df_corrected$m))
  n_colors <- length(m_values)
  colors <- viridisLite::viridis(n_colors, option = palette_option)
  names(colors) <- paste0("m = ", m_values)
  
  # Add benchmark color
  all_colors <- c("Naive" = benchmark_color, colors)
  
  # Build title/subtitle layer
  labs_layer <- if (isTRUE(title_bool)) {
    ggplot2::labs(
      title = sprintf("MSE vs Lag Norm: λ = %s, N = %d × %d", 
                      lambda_filter, N_filter, N_filter),
      subtitle = "",
      x = "Lag Norm ||h||",
      y = "MSE",
      color = "Estimator"
    )
  } else {
    ggplot2::labs(
      x = "Lag Norm ||h||",
      y = "MSE",
      color = "Estimator"
    )
  }
  
  # Create plot
  p <- ggplot2::ggplot() +
    # Benchmark line (uncorrected MSE from 'MSEChat' column)
    ggplot2::geom_line(
      data = df_benchmark,
      ggplot2::aes(x = lag_norm, y = benchmark_mse, color = estimator_type),
      linetype = benchmark_linetype,
      linewidth = benchmark_size
    ) +
    ggplot2::geom_point(
      data = df_benchmark,
      ggplot2::aes(x = lag_norm, y = benchmark_mse, color = estimator_type),
      size = point_size
    ) +
    # Corrected estimator lines (MSEChatCorr for all m)
    ggplot2::geom_line(
      data = df_corrected,
      ggplot2::aes(x = lag_norm, y = MSEChatCorr, color = estimator_type, group = m),
      linewidth = line_size
    ) +
    ggplot2::geom_point(
      data = df_corrected,
      ggplot2::aes(x = lag_norm, y = MSEChatCorr, color = estimator_type),
      size = point_size
    ) +
    ggplot2::scale_color_manual(
      values = all_colors,
      name = "Estimator",
      breaks = c("Naive", paste0("m = ", m_values))
    ) +
    labs_layer +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = if (isTRUE(show_legend)) "right" else "none",
      panel.grid.minor = ggplot2::element_blank()
    )
  
  # Save plot
  out_path <- file.path(out_dir, filename)
  ggplot2::ggsave(out_path, plot = p, width = width, height = height, 
                  dpi = dpi, bg = "white")
  
  if (interactive()) print(p)
  
  invisible(list(plot = p, file = out_path, 
                 data = list(corrected = df_corrected, benchmark = df_benchmark)))
}


#-------------------------------------------------------------------------------
plot_best_m_min_MSE_flags_general <- function(summary_df,
                                      lambda_filter,
                                      metric = "MSE",           # "MSE", "bias", "Var"
                                      use_N_version = FALSE,    # TRUE for _N columns
                                      title_bool = TRUE,
                                      out_dir = ".",
                                      filename = NULL,            # auto if NULL
                                      width = 8, height = 6, dpi = 300,
                                      label_size = 3,
                                      palette_option = "plasma",
                                      palette_direction = -1,
                                      msered_col = NULL,         # auto-determined if NULL
                                      show_legend = TRUE) {
  
  n_i = unique(summary_df$N)
  .n_tag <- function(Nvec) if (length(Nvec)) paste(sort(unique(Nvec)), collapse = "-") else "NA"
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  df_lam <- dplyr::filter(summary_df, .data$lambda == !!lambda_filter)
  if (nrow(df_lam) == 0) stop("No rows for the requested lambda_filter.")
  
  # Construct column names based on metric and use_N_version
  n_suffix <- if (use_N_version) "_N" else ""
  
  metric_col <- paste0(metric, "ChatCorr", n_suffix)
  reduced_col <- paste0(metric, "Reduced_bool", n_suffix)
  
  # Auto-determine msered_col if not provided
  if (is.null(msered_col)) {
    msered_col <- paste0(metric, "ZeroRed_bool", n_suffix)
  }
  
  # Check if columns exist
  if (!metric_col %in% names(df_lam)) {
    stop(sprintf("Column '%s' not found in data.", metric_col))
  }
  if (!reduced_col %in% names(df_lam)) {
    stop(sprintf("Column '%s' not found in data.", reduced_col))
  }
  
  if (is.null(filename)) {
    n_tag <- .n_tag(df_lam$N)
    metric_tag <- paste0(metric, if (use_N_version) "_N" else "")
    if(title_bool){
      filename <- sprintf("bestm_min%s_flags_lambda%s_N%s.png",
                          metric_tag, as.character(lambda_filter), n_tag)
    }else{
      filename <- sprintf("non_bestm_min%s_flags_lambda%s_N%s.png",
                          metric_tag, as.character(lambda_filter), n_tag)
    }
  }
  out_path <- file.path(out_dir, filename)
  
  if (!msered_col %in% names(df_lam)) {
    warning(sprintf("Column '%s' not found; X overlay will be skipped.", msered_col))
    df_lam[[msered_col]] <- NA
  }
  
  grid_hh <- df_lam %>% dplyr::distinct(h1, h2)
  
  # Use the dynamically determined reduced_col
  df_reduced <- df_lam %>% dplyr::filter(.data[[reduced_col]])
  
  df_best <- df_reduced %>%
    dplyr::group_by(h1, h2) %>%
    dplyr::arrange(.data[[metric_col]], .data$m, .by_group = TRUE) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(
      h1, h2,
      best_m = .data$m,
      best_metric = .data[[metric_col]],
      zero_reduced = .data[[msered_col]]
    )
  
  df_plot <- dplyr::left_join(grid_hh, df_best, by = c("h1","h2"))
  
  pal_fun <- scales::viridis_pal(option = palette_option, direction = palette_direction)
  
  if (any(!is.na(df_plot$best_m))) {
    rng <- range(df_plot$best_m, na.rm = TRUE)
    df_plot <- df_plot %>%
      dplyr::mutate(
        scale_pos = ifelse(is.na(.data$best_m), NA_real_,
                           scales::rescale(.data$best_m, to = c(0,1), from = rng)),
        idx       = ifelse(is.na(.data$scale_pos), NA_integer_,
                           pmin(101L, pmax(1L, as.integer(ceiling(.data$scale_pos * 100) + 1)))),
        fill_hex  = dplyr::if_else(is.na(.data$idx), NA_character_, pal_fun(101)[.data$idx]),
        luminance = dplyr::if_else(
          is.na(.data$fill_hex), NA_real_,
          {
            rgb_mat <- as.data.frame(t(col2rgb(.data$fill_hex)))
            (0.2126 * rgb_mat$red + 0.7152 * rgb_mat$green + 0.0722 * rgb_mat$blue) / 255
          }
        ),
        text_col  = dplyr::if_else(is.na(.data$luminance), "black",
                                   dplyr::if_else(.data$luminance < 0.5, "white", "black")),
        label     = dplyr::if_else(is.na(.data$best_m), "", as.character(.data$best_m))
      )
  } else {
    df_plot <- df_plot %>%
      dplyr::mutate(fill_hex = NA_character_, text_col = "black", label = "")
  }
  
  df_plot <- df_plot %>%
    dplyr::mutate(
      x = as.numeric(factor(.data$h1, levels = sort(unique(.data$h1)))),
      y = as.numeric(factor(.data$h2, levels = sort(unique(.data$h2))))
    )
  
  # Build title with metric information
  metric_name <- paste0(metric, if (use_N_version) " (N-version)" else "")
  
  labs_layer <- if (isTRUE(title_bool)) {
    ggplot2::labs(
      title = paste0("m with smallest ", metric_name, " of Corrected Estimator with \u03BB = ",
                     lambda_filter, ", N = ", n_i, " \u00D7 ", n_i),
      subtitle = paste0("White = No ", metric, " Reduction at (h1,h2); X = Zero Estimation has smaller ", metric),
      x = "h1", y = "h2"
    )
  } else {
    ggplot2::labs(x = "h1", y = "h2")
  }
  
  p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = factor(h1), y = factor(h2), fill = best_m)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.3, na.rm = FALSE) +
    ggplot2::geom_text(ggplot2::aes(label = label, colour = I(text_col)), size = label_size) +
    ggplot2::geom_segment(
      data = df_plot %>% dplyr::filter(!is.na(best_m), !is.na(zero_reduced), zero_reduced == FALSE),
      ggplot2::aes(x = x - 0.45, xend = x + 0.45, y = y - 0.45, yend = y + 0.45),
      inherit.aes = FALSE, linewidth = 0.6
    ) +
    ggplot2::geom_segment(
      data = df_plot %>% dplyr::filter(!is.na(best_m), !is.na(zero_reduced), zero_reduced == FALSE),
      ggplot2::aes(x = x - 0.45, xend = x + 0.45, y = y + 0.45, yend = y - 0.45),
      inherit.aes = FALSE, linewidth = 0.6
    ) +
    ggplot2::scale_fill_viridis_c(
      option = palette_option,
      direction = palette_direction,
      name = "Best m",
      na.value = "white"
    ) +
    labs_layer +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      legend.position = if (isTRUE(show_legend)) "right" else "none"
    )
  
  # Extra safety: remove legend when show_legend = FALSE
  if (!isTRUE(show_legend)) {
    p <- p + ggplot2::guides(fill = "none")
  }
  
  ggplot2::ggsave(out_path, plot = p, width = width, height = height, dpi = dpi, bg = "white")
  if (interactive()) print(p)
  invisible(list(plot = p, file = out_path))
}

#-------------------------------------------------------------------------------
plot_bias_vs_norm_general <- function(summary_df,
                                      lambda_filter,
                                      N_filter,
                                      metric = "bias",          # "bias", "MSE", "Var"
                                      use_N_version = FALSE,    # TRUE for _N columns
                                      title_bool = TRUE,
                                      out_dir = ".",
                                      filename = NULL,
                                      width = 10, height = 6, dpi = 300,
                                      palette_option = "plasma",
                                      show_legend = TRUE,
                                      line_size = 0.8,
                                      point_size = 2,
                                      benchmark_color = "black",
                                      benchmark_size = 0.9,
                                      benchmark_linetype = "solid") {
  
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  # Filter data
  df_filtered <- summary_df %>%
    dplyr::filter(lambda == lambda_filter, N == N_filter)
  
  if (nrow(df_filtered) == 0) {
    stop("No data found for the specified lambda and N filters.")
  }
  
  # Construct column names based on metric and use_N_version
  n_suffix <- if (use_N_version) "_N" else ""
  
  # Determine column names based on metric
  if (metric == "bias") {
    col_original <- paste0("bias", n_suffix)
    col_corrected <- paste0("biasCorr", n_suffix)
  } else if (metric == "MSE") {
    col_original <- paste0("MSEChat", n_suffix)
    col_corrected <- paste0("MSEChatCorr", n_suffix)
  } else if (metric == "Var") {
    col_original <- paste0("VarChat", n_suffix)
    col_corrected <- paste0("VarChatCorr", n_suffix)
  } else {
    stop("metric must be one of: 'bias', 'MSE', 'Var'")
  }
  
  # Check if columns exist
  if (!col_original %in% names(df_filtered)) {
    stop(sprintf("Column '%s' not found in data.", col_original))
  }
  if (!col_corrected %in% names(df_filtered)) {
    stop(sprintf("Column '%s' not found in data.", col_corrected))
  }
  
  # Build default filename
  if (is.null(filename)) {
    metric_tag <- paste0(metric, if (use_N_version) "_N" else "")
    if (title_bool) {
      filename <- sprintf("%s_vs_norm_lambda%s_N%s.png", 
                          tolower(metric_tag), as.character(lambda_filter), N_filter)
    } else {
      filename <- sprintf("non_%s_vs_norm_lambda%s_N%s.png", 
                          tolower(metric_tag), as.character(lambda_filter), N_filter)
    }
  }
  
  # Calculate lag norm and aggregate by unique norm values
  df_plot <- df_filtered %>%
    dplyr::mutate(
      lag_norm = sqrt(h1^2 + h2^2)
    ) %>%
    dplyr::group_by(m, lag_norm) %>%
    dplyr::summarise(
      original_metric = mean(.data[[col_original]], na.rm = TRUE),
      corrected_metric = mean(.data[[col_corrected]], na.rm = TRUE),
      .groups = "drop"
    )
  
  # Benchmark: uncorrected metric (same for all m at each lag_norm)
  df_benchmark <- df_plot %>%
    dplyr::group_by(lag_norm) %>%
    dplyr::summarise(
      benchmark_value = first(original_metric),
      .groups = "drop"
    ) %>%
    dplyr::mutate(estimator_type = "Benchmark")
  
  # Data for corrected estimators (all m values with corrected metric)
  df_corrected <- df_plot %>%
    dplyr::mutate(estimator_type = paste0("m = ", m))
  
  # Create color palette for m values
  m_values <- sort(unique(df_corrected$m))
  n_colors <- length(m_values)
  
  # Check if palette_option is a vector of colors or a viridis name
  if (length(palette_option) > 1) {
    # Custom palette provided
    colors <- colorRampPalette(palette_option)(n_colors)
  } else {
    # Viridis palette name provided
    colors <- viridisLite::viridis(n_colors, option = palette_option)
  }
  names(colors) <- paste0("m = ", m_values)
  
  # Add benchmark color
  all_colors <- c("Benchmark" = benchmark_color, colors)
  
  # Determine metric display name
  metric_display <- switch(metric,
                           "bias" = "Bias",
                           "MSE" = "MSE",
                           "Var" = "Variance")
  metric_name <- paste0(metric_display, if (use_N_version) " (N-version)" else "")
  
  # Build title/subtitle layer
  labs_layer <- if (isTRUE(title_bool)) {
    ggplot2::labs(
      title = sprintf("%s vs Lag Norm: λ = %s, N = %d × %d", 
                      metric_name, lambda_filter, N_filter, N_filter),
      x = "Lag Norm ||h||",
      y = metric_name,
      color = "Estimator"
    )
  } else {
    ggplot2::labs(
      x = "Lag Norm ||h||",
      y = metric_name,
      color = "Estimator"
    )
  }
  
  # Create plot
  p <- ggplot2::ggplot() +
    # Benchmark line (uncorrected metric)
    ggplot2::geom_line(
      data = df_benchmark,
      ggplot2::aes(x = lag_norm, y = benchmark_value, color = estimator_type),
      linetype = benchmark_linetype,
      linewidth = benchmark_size
    ) +
    ggplot2::geom_point(
      data = df_benchmark,
      ggplot2::aes(x = lag_norm, y = benchmark_value, color = estimator_type),
      size = point_size
    ) +
    # Corrected estimator lines (corrected metric for all m)
    ggplot2::geom_line(
      data = df_corrected,
      ggplot2::aes(x = lag_norm, y = corrected_metric, color = estimator_type, group = m),
      linewidth = line_size
    ) +
    ggplot2::geom_point(
      data = df_corrected,
      ggplot2::aes(x = lag_norm, y = corrected_metric, color = estimator_type),
      size = point_size
    ) +
    ggplot2::scale_color_manual(
      values = all_colors,
      name = "Estimator",
      breaks = c("Benchmark", paste0("m = ", m_values))
    ) +
    labs_layer +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = if (isTRUE(show_legend)) "right" else "none",
      panel.grid.minor = ggplot2::element_blank()
    )
  
  # Add horizontal line at y = 0 (only for bias)
  if (metric == "bias") {
    p <- p + ggplot2::geom_hline(yintercept = 0, linetype = "dotted", 
                                 color = "gray40", linewidth = 0.5)
  }
  
  # Save plot
  out_path <- file.path(out_dir, filename)
  ggplot2::ggsave(out_path, plot = p, width = width, height = height, 
                  dpi = dpi, bg = "white")
  
  if (interactive()) print(p)
  
  invisible(list(plot = p, file = out_path, 
                 data = list(corrected = df_corrected, benchmark = df_benchmark)))
}

plot_metric_vs_norm_grid <- function(summary_df,
                                     N_filter,
                                     lambda_values = c(2, 4, 6, 8),
                                     metric = "bias",
                                     use_N_version = FALSE,
                                     out_dir = ".",
                                     filename = NULL,
                                     width = 14, height = 10, dpi = 300,
                                     palette_option = "plasma",
                                     line_size = 0.8,
                                     point_size = 2,
                                     benchmark_color = "black",
                                     benchmark_size = 0.9,
                                     benchmark_linetype = "solid",
                                     ncol = 2,
                                     axis_text_size = 11) {  # New parameter for axis text size
  
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  # Check if we have exactly 4 lambda values for 2x2 grid
  if (length(lambda_values) != 4) {
    warning("Expected 4 lambda values for 2x2 grid. Adjusting layout.")
  }
  
  # Construct column names based on metric and use_N_version
  n_suffix <- if (use_N_version) "_N" else ""
  
  # Determine column names based on metric
  if (metric == "bias") {
    col_original <- paste0("bias", n_suffix)
    col_corrected <- paste0("biasCorr", n_suffix)
  } else if (metric == "MSE") {
    col_original <- paste0("MSEChat", n_suffix)
    col_corrected <- paste0("MSEChatCorr", n_suffix)
  } else if (metric == "Var") {
    col_original <- paste0("VarChat", n_suffix)
    col_corrected <- paste0("VarChatCorr", n_suffix)
  } else {
    stop("metric must be one of: 'bias', 'MSE', 'Var'")
  }
  
  # Filter data for all lambda values at once
  df_filtered <- summary_df %>%
    dplyr::filter(N == N_filter, lambda %in% lambda_values)
  
  if (nrow(df_filtered) == 0) {
    stop("No data found for the specified N and lambda filters.")
  }
  
  # Check if columns exist
  if (!col_original %in% names(df_filtered)) {
    stop(sprintf("Column '%s' not found in data.", col_original))
  }
  if (!col_corrected %in% names(df_filtered)) {
    stop(sprintf("Column '%s' not found in data.", col_corrected))
  }
  
  # Prepare data for all subplots
  df_plot <- df_filtered %>%
    dplyr::mutate(
      lag_norm = sqrt(h1^2 + h2^2),
      lambda_label = paste0("λ = ", lambda)
    ) %>%
    dplyr::group_by(lambda, lambda_label, m, lag_norm) %>%
    dplyr::summarise(
      original_metric = mean(.data[[col_original]], na.rm = TRUE),
      corrected_metric = mean(.data[[col_corrected]], na.rm = TRUE),
      .groups = "drop"
    )
  
  # Benchmark data (uncorrected metric)
  df_benchmark <- df_plot %>%
    dplyr::group_by(lambda, lambda_label, lag_norm) %>%
    dplyr::summarise(
      benchmark_value = first(original_metric),
      .groups = "drop"
    ) %>%
    dplyr::mutate(estimator_type = "Naive")
  
  # Corrected estimator data
  df_corrected <- df_plot %>%
    dplyr::mutate(estimator_type = paste0("m = ", m))
  
  # Create color palette for m values
  m_values <- sort(unique(df_corrected$m))
  n_colors <- length(m_values)
  
  # Check if palette_option is a vector of colors or a viridis name
  if (length(palette_option) > 1) {
    colors <- colorRampPalette(palette_option)(n_colors)
  } else {
    colors <- viridisLite::viridis(n_colors, option = palette_option)
  }
  names(colors) <- paste0("m = ", m_values)
  
  # Add benchmark color
  all_colors <- c("Naive" = benchmark_color, colors)
  
  # Determine metric display name
  metric_display <- switch(metric,
                           "bias" = "Bias",
                           "MSE" = "MSE",
                           "Var" = "Variance")
  metric_name <- paste0(metric_display, if (use_N_version) " (N-version)" else "")
  
  # Create the plot
  p <- ggplot2::ggplot() +
    # Benchmark lines
    ggplot2::geom_line(
      data = df_benchmark,
      ggplot2::aes(x = lag_norm, y = benchmark_value, color = estimator_type),
      linetype = benchmark_linetype,
      linewidth = benchmark_size
    ) +
    ggplot2::geom_point(
      data = df_benchmark,
      ggplot2::aes(x = lag_norm, y = benchmark_value, color = estimator_type),
      size = point_size
    ) +
    # Corrected estimator lines
    ggplot2::geom_line(
      data = df_corrected,
      ggplot2::aes(x = lag_norm, y = corrected_metric, color = estimator_type, group = m),
      linewidth = line_size
    ) +
    ggplot2::geom_point(
      data = df_corrected,
      ggplot2::aes(x = lag_norm, y = corrected_metric, color = estimator_type),
      size = point_size
    ) +
    ggplot2::scale_color_manual(
      values = all_colors,
      name = "Estimator",
      breaks = c("Naive", paste0("m = ", m_values))
    ) +
    ggplot2::facet_wrap(~ lambda_label, ncol = ncol, scales = "free_y") +
    ggplot2::labs(
      title = sprintf("%s vs Lag Norm for N = %d × %d", metric_name, N_filter, N_filter),
      x = "Lag Norm ||h||",
      y = metric_name
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "right",
      panel.grid.minor = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = 11, face = "bold"),
      strip.background = ggplot2::element_blank(),  # Remove grey background
      axis.text = ggplot2::element_text(size = axis_text_size),  # Larger axis numbers
      axis.text.x = ggplot2::element_text(size = axis_text_size),
      axis.text.y = ggplot2::element_text(size = axis_text_size)
    )
  
  # Add horizontal line at y = 0 (only for bias)
  if (metric == "bias") {
    p <- p + ggplot2::geom_hline(yintercept = 0, linetype = "dotted", 
                                 color = "gray40", linewidth = 0.5)
  }
  
  # Build default filename
  if (is.null(filename)) {
    metric_tag <- paste0(metric, if (use_N_version) "_N" else "")
    lambda_tag <- paste(lambda_values, collapse = "_")
    filename <- sprintf("%s_vs_norm_grid_N%s_lambda%s.png", 
                        tolower(metric_tag), N_filter, lambda_tag)
  }
  
  # Save plot
  out_path <- file.path(out_dir, filename)
  ggplot2::ggsave(out_path, plot = p, width = width, height = height, 
                  dpi = dpi, bg = "white")
  
  if (interactive()) print(p)
  
  invisible(list(
    plot = p, 
    file = out_path,
    data = list(corrected = df_corrected, benchmark = df_benchmark)
  ))
}

#-------------------------------------------------------------------------------
plot_metric_vs_norm_publication <- function(
    summary_df,
    N_filter,
    lambda_values = c(2, 4, 6, 8),
    metric = "bias",
    normalization = "Nh",
    max_m = NULL,
    best_m_only = FALSE,
    show_CH = FALSE,
    CH_color = "#E41A1C",
    naive_color = "black",
    out_dir = ".",
    filename = NULL,
    width = 7,
    height = 6,
    dpi = 600,  # High DPI for publication
    base_size = 11,
    font_scale = 1,  # Multiplier for all text sizes (axis labels, titles, legend, facet strips)
    line_size = 0.6,
    point_size = 1.8,
    use_bw = FALSE,  # Black & white friendly with shapes
    show_title = FALSE,  # Often not needed in papers
    legend_position = "right",
    palette_option = "viridis",  # viridis palette: "viridis", "plasma", "magma", "inferno", "cividis"
    CH_clip_factor = 1.5  # For MSE: drop CH points exceeding this factor times the naive value
) {
  
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  # Validate normalization parameter
  if (!normalization %in% c("Nh", "N")) {
    stop("normalization must be either 'Nh' or 'N'")
  }
  
  # Compute scaled font sizes
  sz_axis_title  <- base_size * font_scale        # Axis names (x, y labels)
  sz_axis_text   <- (base_size - 1) * font_scale  # Axis tick labels (numbers/letters)
  sz_legend_title <- base_size * font_scale        # Legend title
  sz_legend_text  <- (base_size - 1) * font_scale  # Legend entries
  sz_strip_text   <- base_size * font_scale        # Facet strip labels
  
  # Construct column names based on metric and normalization
  n_suffix <- if (normalization == "N") "_N" else ""
  
  if (metric == "bias") {
    col_original <- paste0("bias", n_suffix)
    col_corrected <- paste0("biasCorr", n_suffix)
    col_CH <- "biasCH"
    y_label <- "Bias"
  } else if (metric == "MSE") {
    col_original <- paste0("MSEChat", n_suffix)
    col_corrected <- paste0("MSEChatCorr", n_suffix)
    col_CH <- "MSECH"
    y_label <- "MSE"
  } else if (metric == "Var") {
    col_original <- paste0("VarChat", n_suffix)
    col_corrected <- paste0("VarChatCorr", n_suffix)
    col_CH <- "VarCH"
    y_label <- "Variance"
  } else {
    stop("metric must be one of: 'bias', 'MSE', 'Var'")
  }
  
  #if (normalization == "N") y_label <- paste0(y_label, " (N-normalized)")
  
  # Filter data
  df_filtered <- summary_df %>%
    dplyr::filter(N == N_filter, lambda %in% lambda_values)
  
  # Filter by max_m if specified
  if (!is.null(max_m)) {
    df_filtered <- df_filtered %>%
      dplyr::filter(m <= max_m)
  }
  
  if (nrow(df_filtered) == 0) {
    stop("No data found for the specified N and lambda filters.")
  }
  
  # Prepare data
  df_plot <- df_filtered %>%
    dplyr::mutate(
      lag_norm = sqrt(h1^2 + h2^2),
      # Use expression-friendly labels
      lambda_label = factor(
        lambda,
        levels = lambda_values,
        labels = paste0("lambda == ", lambda_values)
      )
    ) %>%
    dplyr::group_by(lambda, lambda_label, m, lag_norm) %>%
    dplyr::summarise(
      original_metric = mean(.data[[col_original]], na.rm = TRUE),
      corrected_metric = mean(.data[[col_corrected]], na.rm = TRUE),
      CH_metric = if (col_CH %in% names(df_filtered)) mean(.data[[col_CH]], na.rm = TRUE) else NA_real_,
      .groups = "drop"
    )
  
  # Benchmark data
  df_benchmark <- df_plot %>%
    dplyr::group_by(lambda, lambda_label, lag_norm) %>%
    dplyr::summarise(
      benchmark_value = first(original_metric),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      estimator = "Naïve",
      m_factor = factor("Naïve")
    )
  
  # Corrected estimator data
  m_values <- sort(unique(df_plot$m))
  
  # If best_m_only, select the m that minimizes the corrected metric for each (lambda, lag_norm)
  if (best_m_only) {
    df_corrected <- df_plot %>%
      dplyr::group_by(lambda, lambda_label, lag_norm) %>%
      dplyr::slice_min(order_by = abs(corrected_metric), n = 1, with_ties = FALSE) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        estimator = "Best m",
        m_factor = factor("Best m")
      )
    m_values <- "best"  # Placeholder for legend
  } else {
    df_corrected <- df_plot %>%
      dplyr::mutate(
        estimator = paste0("m = ", m),
        m_factor = factor(paste0("m = ", m), levels = paste0("m = ", m_values))
      )
  }
  
  # Cressie-Hawkins data (if requested)
  if (show_CH && col_CH %in% names(df_filtered)) {
    df_CH <- df_plot %>%
      dplyr::group_by(lambda, lambda_label, lag_norm) %>%
      dplyr::summarise(
        CH_value = first(CH_metric),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        estimator = "CH",
        m_factor = factor("CH")
      )
    
    # For MSE: clip CH points that exceed CH_clip_factor * naive value
    if (metric == "MSE") {
      df_CH <- df_CH %>%
        dplyr::left_join(
          df_benchmark %>% dplyr::select(lambda, lag_norm, benchmark_value),
          by = c("lambda", "lag_norm")
        ) %>%
        dplyr::filter(CH_value <= CH_clip_factor * benchmark_value) %>%
        dplyr::select(-benchmark_value)
    }
  }
  
  # Combine for unified legend
  df_combined <- dplyr::bind_rows(
    df_benchmark %>% 
      dplyr::rename(y_value = benchmark_value) %>%
      dplyr::select(lambda, lambda_label, lag_norm, y_value, estimator, m_factor),
    df_corrected %>% 
      dplyr::rename(y_value = corrected_metric) %>%
      dplyr::select(lambda, lambda_label, lag_norm, y_value, estimator, m_factor, m)
  )
  
  # Add CH if requested
  if (show_CH && col_CH %in% names(df_filtered)) {
    df_combined <- dplyr::bind_rows(
      df_combined,
      df_CH %>%
        dplyr::rename(y_value = CH_value) %>%
        dplyr::select(lambda, lambda_label, lag_norm, y_value, estimator, m_factor)
    )
  }
  
  # Create estimator factor with proper ordering
  if (best_m_only) {
    all_levels <- c("Naïve", "Best m")
  } else {
    all_levels <- c("Naïve", paste0("m = ", m_values))
  }
  if (show_CH && col_CH %in% names(df_filtered)) {
    all_levels <- c(all_levels, "CH")
  }
  df_combined$estimator <- factor(df_combined$estimator, levels = all_levels)
  
  # Define aesthetics based on use_bw
  n_estimators <- length(all_levels)
  
  if (best_m_only) {
    # Simplified palette for best_m_only mode
    if (use_bw) {
      colors <- c(naive_color, "gray40")
      linetypes <- c("solid", "dashed")
      shapes <- c(16, 17)
    } else {
      colors <- c(naive_color, viridis(1, option = palette_option, begin = 0.5))
      linetypes <- c("solid", "solid")
      shapes <- c(16, 17)
    }
    
    # Add CH styling if present
    if (show_CH && "CH" %in% all_levels) {
      colors <- c(colors, if (use_bw) "gray20" else CH_color)
      linetypes <- c(linetypes, if (use_bw) "dotdash" else "dashed")
      shapes <- c(shapes, 4)
    }
    
    names(colors) <- all_levels
    names(linetypes) <- all_levels
    names(shapes) <- all_levels
    
  } else if (use_bw) {
    # Black & white friendly: use shapes and line types
    colors <- c(naive_color, rep("gray40", length(m_values)))
    linetypes <- c("solid", rep("dashed", length(m_values)))
    shapes <- c(16, 17, 15, 18, 8, 4, 3, 1, 2, 0)[1:(1 + length(m_values))]
    
    # Add CH styling if present
    if (show_CH && "CH" %in% all_levels) {
      colors <- c(colors, "gray20")
      linetypes <- c(linetypes, "dotdash")
      shapes <- c(shapes, 4)  # X shape for CH
    }
    
    names(colors) <- all_levels
    names(linetypes) <- all_levels
    names(shapes) <- all_levels
    
  } else {
    # Color version with viridis palette
    colors <- c(naive_color, viridis(length(m_values), option = palette_option, end = 0.9))
    linetypes <- c("solid", rep("solid", length(m_values)))
    shapes <- c(16, rep(16, length(m_values)))
    
    # Add CH styling if present
    if (show_CH && "CH" %in% all_levels) {
      colors <- c(colors, CH_color)
      linetypes <- c(linetypes, "dashed")
      shapes <- c(shapes, 17)  # Triangle for CH
    }
    
    names(colors) <- all_levels
    names(linetypes) <- all_levels
    names(shapes) <- all_levels
  }
  
  # Build the plot
  p <- ggplot(df_combined, aes(x = lag_norm, y = y_value, 
                               color = estimator, 
                               linetype = estimator,
                               shape = estimator)) +
    # Lines
    geom_line(aes(group = estimator), linewidth = line_size) +
    # Points
    geom_point(size = point_size) +
    # Facets with parsed labels (for Greek lambda)
    facet_wrap(~ lambda_label, ncol = 2, scales = "free_y",
               labeller = label_parsed) +
    # Scales
    scale_color_manual(values = colors, name = "Estimator") +
    scale_linetype_manual(values = linetypes, name = "Estimator") +
    scale_shape_manual(values = shapes, name = "Estimator") +
    # Labels
    labs(
      x = "Lag norm",
      y = y_label
    ) +
    # Publication theme
    theme_bw(base_size = base_size) +
    theme(
      # Panel appearance
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.3, color = "gray85"),
      panel.border = element_rect(linewidth = 0.5, color = "black"),
      
      # Strip (facet labels)
      strip.background = element_blank(),
      strip.text = element_text(size = sz_strip_text, face = "plain", 
                                margin = margin(b = 5, t = 5)),
      
      # Legend
      legend.position = legend_position,
      legend.title = element_text(size = sz_legend_title, face = "plain"),
      legend.text = element_text(size = sz_legend_text),
      legend.key.width = unit(1.5, "cm"),
      legend.background = element_blank(),
      legend.box.background = element_blank(),
      
      # Axis
      axis.title = element_text(size = sz_axis_title),
      axis.text = element_text(size = sz_axis_text, color = "black"),
      axis.ticks = element_line(linewidth = 0.3, color = "black"),
      
      # Margins
      plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
    )
  
  # Add horizontal line at 0 for bias
  if (metric == "bias") {
    p <- p + geom_hline(yintercept = 0, linetype = "dotted", 
                        color = "gray50", linewidth = 0.4)
  }
  
  # Add title if requested
  if (show_title) {
    p <- p + ggtitle(bquote(.(y_label) ~ "vs lag norm for" ~ N == .(N_filter) %*% .(N_filter)))
  }
  
  # Adjust legend layout based on position
  if (legend_position == "bottom") {
    p <- p + guides(
      color = guide_legend(nrow = 1, byrow = TRUE),
      linetype = guide_legend(nrow = 1, byrow = TRUE),
      shape = guide_legend(nrow = 1, byrow = TRUE)
    )
  } else if (legend_position == "right") {
    p <- p + guides(
      color = guide_legend(ncol = 1),
      linetype = guide_legend(ncol = 1),
      shape = guide_legend(ncol = 1)
    )
  }
  
  # Build filename
  if (is.null(filename)) {
    metric_tag <- paste0(metric, "_", normalization)
    if (best_m_only) metric_tag <- paste0(metric_tag, "_bestm")
    if (show_CH) metric_tag <- paste0(metric_tag, "_withCH")
    filename <- sprintf("%s_vs_norm_N%d_publication.png", 
                        tolower(metric_tag), N_filter)
  }
  
  # Ensure filename ends with .png
  if (!grepl("\\.png$", filename, ignore.case = TRUE)) {
    filename <- sub("\\.[^.]+$", "", filename)  # strip any existing extension
    filename <- paste0(filename, ".png")
  }
  
  out_path_png <- file.path(out_dir, filename)
  
  ggsave(out_path_png, plot = p, width = width, height = height, dpi = dpi, bg = "white")
  
  message("Saved: ", out_path_png)
  
  invisible(list(
    plot = p,
    png_file = out_path_png,
    data = df_combined
  ))
}

plot_mse_ratio_heatmap_grid <- function(
    summary_df,
    N_filter,
    lambda_values = c(2, 4, 6, 8),
    normalization = "Nh",
    best_m_only = TRUE,
    m_filter = NULL,
    show_numbers = TRUE,
    out_dir = ".",
    filename = NULL,
    width = 8,
    height = 8,
    dpi = 600,
    label_size = 3,
    show_legend = TRUE,
    palette_option = "diverging",
    palette_direction = 1,
    fixed_aspect = FALSE,
    ncol = 2
) {
  
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  # Construct column names
  n_suffix <- if (normalization == "N") "_N" else ""
  col_naive <- paste0("MSEChat", n_suffix)
  col_corrected <- paste0("MSEChatCorr", n_suffix)
  col_zero <- "MSEZeroEst"
  
  # Filter data by N and lambda values
  df_filtered <- summary_df %>%
    dplyr::filter(N == N_filter, lambda %in% lambda_values)
  
  if (nrow(df_filtered) == 0) stop("No rows for the requested N_filter and lambda values.")
  
  n_i <- N_filter
  
  # Get best m or filter by specific m for each lambda
  if (best_m_only) {
    df_best <- df_filtered %>%
      dplyr::group_by(lambda, h1, h2) %>%
      dplyr::slice_min(order_by = .data[[col_corrected]], n = 1, with_ties = FALSE) %>%
      dplyr::ungroup()
    m_label <- "best m"
  } else {
    if (is.null(m_filter)) stop("Must specify m_filter when best_m_only = FALSE")
    df_best <- df_filtered %>%
      dplyr::filter(m == m_filter)
    m_label <- paste0("m = ", m_filter)
  }
  
  # Calculate ratio and percentage
  df_plot <- df_best %>%
    dplyr::mutate(
      mse_ratio = .data[[col_corrected]] / .data[[col_naive]],
      pct_change = (mse_ratio - 1) * 100,
      # Zero beats the best available (whichever is smaller: naive or corrected)
      best_available = pmin(.data[[col_naive]], .data[[col_corrected]], na.rm = TRUE),
      zero_better = !is.na(.data[[col_zero]]) & !is.na(best_available) &
        (.data[[col_zero]] < best_available),
      # Format label: integer percentage, empty string for NaN
      label = dplyr::if_else(
        is.nan(mse_ratio) | is.na(mse_ratio),
        "",
        sprintf("%+.0f%%", pct_change)
      ),
      lambda_label = factor(lambda, levels = lambda_values,
                            labels = paste0("lambda == ", lambda_values))
    )
  
  # Remove NaN rows
  df_plot <- df_plot %>%
    dplyr::filter(!is.nan(mse_ratio) & !is.na(mse_ratio))
  
  # If show_numbers is FALSE, clear all labels
  if (!show_numbers) {
    df_plot <- df_plot %>%
      dplyr::mutate(label = "")
  }
  
  # Text color based on ratio
  df_plot <- df_plot %>%
    dplyr::mutate(
      text_col = dplyr::case_when(
        mse_ratio < 0.85 ~ "white",
        mse_ratio > 1.15 ~ "white",
        TRUE ~ "black"
      )
    )
  
  # Prepare coordinates for X marks (need per-facet coordinates)
  df_plot <- df_plot %>%
    dplyr::group_by(lambda_label) %>%
    dplyr::mutate(
      x_num = as.numeric(factor(h1, levels = sort(unique(h1)))),
      y_num = as.numeric(factor(h2, levels = sort(unique(h2))))
    ) %>%
    dplyr::ungroup()
  
  df_cross <- df_plot %>%
    dplyr::filter(zero_better == TRUE)
  
  # Determine scale limits
  max_deviation <- max(abs(df_plot$mse_ratio - 1), na.rm = TRUE)
  limit <- ceiling(max_deviation * 10) / 10 + 0.1
  scale_limits <- c(max(0.5, 1 - limit), min(1.5, 1 + limit))
  
  # Build filename
  if (is.null(filename)) {
    m_tag <- if (best_m_only) "bestm" else paste0("m", m_filter)
    norm_tag <- if (normalization == "N") "_N" else "_Nh"
    filename <- sprintf("mse_ratio_heatmap_grid_%s%s_N%d.pdf", m_tag, norm_tag, N_filter)
  }
  
  # Create plot
  p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = factor(h1), y = factor(h2), fill = mse_ratio)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.2)
  
  # Add text labels if show_numbers is TRUE
  if (show_numbers) {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = label, color = I(text_col)),
      size = label_size
    )
  }
  
  # Add X marks for zero_better
  if (nrow(df_cross) > 0) {
    p <- p +
      ggplot2::geom_segment(
        data = df_cross,
        ggplot2::aes(x = x_num - 0.4, xend = x_num + 0.4, 
                     y = y_num - 0.4, yend = y_num + 0.4),
        inherit.aes = FALSE, linewidth = 0.5, color = "black"
      ) +
      ggplot2::geom_segment(
        data = df_cross,
        ggplot2::aes(x = x_num - 0.4, xend = x_num + 0.4, 
                     y = y_num + 0.4, yend = y_num - 0.4),
        inherit.aes = FALSE, linewidth = 0.5, color = "black"
      )
  }
  
  p <- p + ggplot2::facet_wrap(~ lambda_label, ncol = ncol, labeller = ggplot2::label_parsed)
  
  # Add fill scale based on palette option
  if (palette_option == "diverging") {
    p <- p + ggplot2::scale_fill_gradient2(
      low = "#2166AC",
      mid = "#F7F7F7",
      high = "#B2182B",
      midpoint = 1,
      limits = scale_limits,
      name = "MSE Ratio",
      na.value = "grey70",
      oob = scales::squish
    )
  } else {
    p <- p + ggplot2::scale_fill_viridis_c(
      option = palette_option,
      direction = palette_direction,
      limits = scale_limits,
      name = "MSE Ratio",
      na.value = "grey70",
      oob = scales::squish
    )
  }
  
  p <- p +
    ggplot2::labs(
      x = expression(h[1]),
      y = expression(h[2])
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = 12),
      legend.position = if (isTRUE(show_legend)) "right" else "none",
      legend.key.height = ggplot2::unit(2, "cm"),
      legend.key.width = ggplot2::unit(0.4, "cm"),
      axis.text = ggplot2::element_text(size = 9),
      panel.spacing = ggplot2::unit(1, "lines")
    )
  
  # Add coord_fixed only if requested
  if (fixed_aspect) {
    p <- p + ggplot2::coord_fixed()
  }
  
  # Save
  out_path_pdf <- file.path(out_dir, filename)
  out_path_png <- file.path(out_dir, sub("\\.pdf$", ".png", filename))
  
  ggplot2::ggsave(out_path_pdf, plot = p, width = width, height = height, device = cairo_pdf)
  ggplot2::ggsave(out_path_png, plot = p, width = width, height = height, dpi = dpi, bg = "white")
  
  message("Saved: ", out_path_pdf)
  message("Saved: ", out_path_png)
  
  if (interactive()) print(p)
  
  invisible(list(plot = p, pdf_file = out_path_pdf, png_file = out_path_png, data = df_plot))
}


# UPDATED VERSION OF THE BEST HEATMAP m
plot_best_m_heatmap <- function(summary_df,
                                lambda_filter,
                                N_filter,
                                normalization = c("Nh", "N"),
                                title_bool = TRUE,
                                out_dir = ".",
                                filename = NULL,
                                width = 7, height = 7, dpi = 300,
                                label_size = 4,
                                palette_option = "plasma",
                                palette_direction = -1,
                                show_legend = TRUE,
                                tick_every = NULL) {
  # tick_every: show axis labels every N lags (NULL = auto based on grid size).
  #             Always shows the extreme lags and zero.
  
  norm_type <- match.arg(normalization)
  
  # Column mapping based on normalization type
  if (norm_type == "Nh") {
    mse_corr_col    <- "MSEChatCorr"
    mse_naive_col   <- "MSEChat"
    mse_reduced_col <- "MSEReduced_bool"
    norm_label      <- "Nh"
  } else {
    mse_corr_col    <- "MSEChatCorr_N"
    mse_naive_col   <- "MSEChat_N"
    mse_reduced_col <- "MSEReduced_bool_N"
    norm_label      <- "N"
  }
  mse_zero_col <- "MSEZeroEst"
  
  n_i <- N_filter
  .n_tag <- function(Nvec) if (length(Nvec)) paste(sort(unique(Nvec)), collapse = "-") else "NA"
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  df_lam <- dplyr::filter(summary_df,
                          .data$lambda == !!lambda_filter,
                          .data$N == !!N_filter)
  if (nrow(df_lam) == 0) stop("No rows for the requested lambda_filter and N_filter.")
  
  if (is.null(filename)) {
    n_tag <- .n_tag(df_lam$N)
    filename <- sprintf("%sbestm_minMSE_flags_%s_lambda%s_N%s.png",
                        if (title_bool) "" else "non_",
                        norm_label,
                        as.character(lambda_filter), as.character(N_filter))
  }
  out_path <- file.path(out_dir, filename)
  
  # Check MSEZeroEst column exists
  if (!mse_zero_col %in% names(df_lam)) {
    warning(sprintf("Column '%s' not found; X overlay will be skipped.", mse_zero_col))
    df_lam[[mse_zero_col]] <- NA
  }
  
  grid_hh <- df_lam %>% dplyr::distinct(h1, h2)
  
  # Filter to rows where correction reduces MSE
  df_reduced <- df_lam %>% dplyr::filter(.data[[mse_reduced_col]])
  
  # Find best m (smallest MSE of corrected estimator) per (h1, h2)
  df_best <- df_reduced %>%
    dplyr::group_by(h1, h2) %>%
    dplyr::arrange(.data[[mse_corr_col]], .data$m, .by_group = TRUE) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(
      h1, h2,
      best_m = .data$m,
      best_MSE = .data[[mse_corr_col]]
    )
  
  # For each (h1, h2): get the naive MSE and zero MSE (constant across m, take first)
  df_ref <- df_lam %>%
    dplyr::group_by(h1, h2) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(
      h1, h2,
      naive_MSE = .data[[mse_naive_col]],
      zero_MSE  = .data[[mse_zero_col]]
    )
  
  df_plot <- grid_hh %>%
    dplyr::left_join(df_best, by = c("h1", "h2")) %>%
    dplyr::left_join(df_ref, by = c("h1", "h2"))
  
  # Determine the "best available" MSE per cell:
  #   - If correction helped: min(best_MSE, naive_MSE)  (best_MSE should be < naive)
  #   - If no correction helped (grey): naive_MSE
  # Then X-mark if zero_MSE < best_available
  df_plot <- df_plot %>%
    dplyr::mutate(
      best_available_MSE = dplyr::if_else(
        is.na(best_m),
        naive_MSE,                               # grey tile: naive is best
        pmin(best_MSE, naive_MSE, na.rm = TRUE)  # corrected tile
      ),
      zero_beats = !is.na(zero_MSE) & !is.na(best_available_MSE) &
        (zero_MSE < best_available_MSE)
    )
  
  # Color mapping
  pal_fun <- scales::viridis_pal(option = palette_option, direction = palette_direction)
  
  if (any(!is.na(df_plot$best_m))) {
    rng <- range(df_plot$best_m, na.rm = TRUE)
    df_plot <- df_plot %>%
      dplyr::mutate(
        scale_pos = ifelse(is.na(.data$best_m), NA_real_,
                           scales::rescale(.data$best_m, to = c(0, 1), from = rng)),
        idx       = ifelse(is.na(.data$scale_pos), NA_integer_,
                           pmin(101L, pmax(1L, as.integer(ceiling(.data$scale_pos * 100) + 1)))),
        fill_hex  = dplyr::if_else(is.na(.data$idx), NA_character_, pal_fun(101)[.data$idx]),
        luminance = dplyr::if_else(
          is.na(.data$fill_hex), NA_real_,
          {
            rgb_mat <- as.data.frame(t(col2rgb(.data$fill_hex)))
            (0.2126 * rgb_mat$red + 0.7152 * rgb_mat$green + 0.0722 * rgb_mat$blue) / 255
          }
        ),
        text_col  = dplyr::if_else(is.na(.data$luminance), "black",
                                   dplyr::if_else(.data$luminance < 0.5, "white", "black")),
        label     = dplyr::if_else(is.na(.data$best_m), "", as.character(.data$best_m))
      )
  } else {
    df_plot <- df_plot %>%
      dplyr::mutate(fill_hex = NA_character_, text_col = "black", label = "")
  }
  
  df_plot <- df_plot %>%
    dplyr::mutate(
      x = as.numeric(factor(.data$h1, levels = sort(unique(.data$h1)))),
      y = as.numeric(factor(.data$h2, levels = sort(unique(.data$h2))))
    )
  
  # --- Axis label thinning ---
  # Helper: given sorted unique lag values, return which to display
  thin_labels <- function(lag_vals, every) {
    n_lags <- length(lag_vals)
    if (is.null(every)) {
      # Auto: show all if <= 15, otherwise thin
      if (n_lags <= 15) return(lag_vals)
      every <- max(2, round(n_lags / 12))
    }
    if (every <= 1) return(lag_vals)
    
    # Always keep: min, max, and 0 (if present)
    keep <- c(lag_vals[1], lag_vals[n_lags])
    if (0 %in% lag_vals) keep <- c(keep, 0)
    
    # From zero outward, keep every `every`-th lag symmetrically
    pos_vals <- lag_vals[lag_vals > 0]
    neg_vals <- lag_vals[lag_vals < 0]
    
    if (length(pos_vals) > 0) {
      keep <- c(keep, pos_vals[seq(every, length(pos_vals), by = every)])
    }
    if (length(neg_vals) > 0) {
      # Reverse so we go outward from zero
      neg_sorted <- rev(neg_vals)
      keep <- c(keep, neg_sorted[seq(every, length(neg_sorted), by = every)])
    }
    
    sort(unique(keep))
  }
  
  h1_vals <- sort(unique(df_plot$h1))
  h2_vals <- sort(unique(df_plot$h2))
  
  h1_show <- thin_labels(h1_vals, tick_every)
  h2_show <- thin_labels(h2_vals, tick_every)
  
  # Build label vectors: show number for kept lags, empty string for others
  h1_labels <- ifelse(h1_vals %in% h1_show, as.character(h1_vals), "")
  h2_labels <- ifelse(h2_vals %in% h2_show, as.character(h2_vals), "")
  
  # Title / labels
  labs_layer <- if (isTRUE(title_bool)) {
    ggplot2::labs(
      title = paste0("Best m (min MSE, ", norm_label, " norm.) | \u03BB = ",
                     lambda_filter, ", N = ", n_i, " \u00D7 ", n_i),
      subtitle = paste0("Grey = No MSE reduction; X = Zero estimation has smaller MSE (",
                        norm_label, ")"),
      x = "h1", y = "h2"
    )
  } else {
    ggplot2::labs(x = "h1", y = "h2")
  }
  
  # X markers: any cell where the zero estimator beats the best available
  df_x <- df_plot %>% dplyr::filter(zero_beats)
  
  p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = factor(h1), y = factor(h2), fill = best_m)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.3, na.rm = FALSE) +
    ggplot2::geom_text(ggplot2::aes(label = label, colour = I(text_col)), size = label_size) +
    ggplot2::geom_segment(
      data = df_x,
      ggplot2::aes(x = x - 0.45, xend = x + 0.45, y = y - 0.45, yend = y + 0.45),
      inherit.aes = FALSE, linewidth = 0.6
    ) +
    ggplot2::geom_segment(
      data = df_x,
      ggplot2::aes(x = x - 0.45, xend = x + 0.45, y = y + 0.45, yend = y - 0.45),
      inherit.aes = FALSE, linewidth = 0.6
    ) +
    ggplot2::scale_x_discrete(labels = h1_labels) +
    ggplot2::scale_y_discrete(labels = h2_labels) +
    ggplot2::scale_fill_viridis_c(
      option = palette_option,
      direction = palette_direction,
      name = "Best m",
      na.value = "white"
    ) +
    labs_layer +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      legend.position = if (isTRUE(show_legend)) "right" else "none"
    )
  
  if (!isTRUE(show_legend)) {
    p <- p + ggplot2::guides(fill = "none")
  }
  
  ggplot2::ggsave(out_path, plot = p, width = width, height = height, dpi = dpi, bg = "white")
  if (interactive()) print(p)
  invisible(list(plot = p, file = out_path))
}
