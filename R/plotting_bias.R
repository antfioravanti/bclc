#------------------------------------------------------------------------------
# Author: Antonio Fioravanti
# Result heatmaps and line plots over a summary_df keyed by N, lambda, m, h1, h2.
#------------------------------------------------------------------------------

# Needed packages:
library(dplyr)
library(ggplot2)
library(scales)
library(viridisLite)


#------------------------------------------------------------------------------
# SHARED HELPERS AND PALETTES
#------------------------------------------------------------------------------

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

#------------------------------------------------------------------------------
# 1) Largest m with MSEReduced
#------------------------------------------------------------------------------
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

#------------------------------------------------------------------------------
# 2) Best m = argmin(MSEChatCorr) among rows with MSEReduced=1
#------------------------------------------------------------------------------
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


#------------------------------------------------------------------------------
# 2) Best m = argmin(MSEChatCorr) among rows with MSEReduced=1 and MSEZred=0
#------------------------------------------------------------------------------
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
#------------------------------------------------------------------------------
# 4) Bias vs Lag Norm for different m values
#------------------------------------------------------------------------------
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
#------------------------------------------------------------------------------
# MSE vs lag norm across m values
#------------------------------------------------------------------------------
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


#------------------------------------------------------------------------------
# Best m heatmap with reduction / zero-estimator flags (generic metric)
#------------------------------------------------------------------------------
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

#------------------------------------------------------------------------------
# Metric vs lag norm across m values (generic metric)
#------------------------------------------------------------------------------
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

#------------------------------------------------------------------------------
# grade metric vs lag norm 
#------------------------------------------------------------------------------

# Colours for the ordered m-family. `palette_option` may be
#   * a viridisLite option name ("viridis", "magma", ..., or "A".."H"),
#   * a paletteer spec "package::palette" (e.g. "scico::berlin",
#     "wesanderson::Zissou1"), 
.m_palette <- function(palette_option, n, begin = 0, end = 1) {
  if (n < 1) return(character(0))

  if (length(palette_option) > 1) {
    ramp <- palette_option                       # explicit colour vector
  } else if (grepl("::", palette_option, fixed = TRUE)) {
    if (!requireNamespace("paletteer", quietly = TRUE)) {
      stop("palette_option \"", palette_option, "\" looks like a paletteer ",
           "palette, but the paletteer package is not installed.\n",
           "  install.packages(\"paletteer\")  - or pass a vector of colours ",
           "instead, e.g. palette_option = c(\"#2166AC\", \"#B2182B\").",
           call. = FALSE)
    }

    ramp <- tryCatch(
      as.character(paletteer::paletteer_c(palette_option, n = 256)),
      error = function(e) tryCatch(
        as.character(paletteer::paletteer_d(palette_option)),
        error = function(e2) stop("paletteer does not know the palette \"",
                                  palette_option, "\" as either continuous or ",
                                  "discrete. See paletteer::palettes_c_names ",
                                  "and paletteer::palettes_d_names.",
                                  call. = FALSE)
      )
    )
  } else {

    known <- c("magma", "inferno", "plasma", "viridis", "cividis",
               "rocket", "mako", "turbo", LETTERS[1:8])
    if (!palette_option %in% known) {
      stop("Unknown palette_option \"", palette_option, "\".\n",
           "  viridisLite names: ", paste(known[1:8], collapse = ", "),
           " (or \"A\"..\"H\").\n",
           "  For a palette from another package use the paletteer form ",
           "\"package::palette\", e.g. \"colorBlindness::", palette_option,
           "\".\n  Or pass a vector of colours, e.g. c(\"#DEEBF7\", \"#08306B\").",
           call. = FALSE)
    }
    
    return(viridisLite::viridis(n, option = palette_option,
                                begin = begin, end = end))
  }

  # Trim, then interpolate to the number of m values
  dense <- grDevices::colorRampPalette(ramp)(256)
  lo <- max(1, floor(begin * 255) + 1)
  hi <- min(256, ceiling(end * 255) + 1)
  if (hi < lo) hi <- lo
  grDevices::colorRampPalette(dense[lo:hi])(n)
}

# Column names + display label of one metric, for a given normalization.
.norm_metric_cols <- function(metric, normalization) {
  n_suffix <- if (normalization == "N") "_N" else ""
  if (metric == "bias") {
    list(original = paste0("bias", n_suffix),
         corrected = paste0("biasCorr", n_suffix),
         CH = "biasCH", label = "Bias")
  } else if (metric == "MSE") {
    list(original = paste0("MSEChat", n_suffix),
         corrected = paste0("MSEChatCorr", n_suffix),
         CH = "MSECH", label = "MSE")
  } else if (metric == "Var") {
    list(original = paste0("VarChat", n_suffix),
         corrected = paste0("VarChatCorr", n_suffix),
         CH = "VarCH", label = "Variance")
  } else {
    stop("metric must be one of: 'bias', 'MSE', 'Var'")
  }
}

# Long-format plotting data for ONE metric: Naive benchmark + corrected (all m or
# best m) + Cressie-Hawkins rows, tagged with `metric` / `metric_label`. Shared by
# the single-metric and the stacked (bias row over MSE row) layouts of
# plot_metric_vs_norm_publication().
.norm_metric_long <- function(df_filtered, metric, normalization, lambda_values,
                              best_m_only, show_CH, CH_clip_factor,
                              max_lag_norm = NULL) {
  cols <- .norm_metric_cols(metric, normalization)
  for (nm in c(cols$original, cols$corrected)) {
    if (!nm %in% names(df_filtered)) {
      stop(sprintf("Column '%s' not found in data (metric = '%s', normalization = '%s').",
                   nm, metric, normalization))
    }
  }
  has_CH <- show_CH && cols$CH %in% names(df_filtered)

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
      original_metric = mean(.data[[cols$original]], na.rm = TRUE),
      corrected_metric = mean(.data[[cols$corrected]], na.rm = TRUE),
      CH_metric = if (has_CH) mean(.data[[cols$CH]], na.rm = TRUE) else NA_real_,
      .groups = "drop"
    )

  # Cut the long lags out of the data rather than zooming the axis, so the free y
  # scales are computed from what is actually visible
  if (!is.null(max_lag_norm)) {
    df_plot <- df_plot %>% dplyr::filter(lag_norm <= max_lag_norm)
    if (nrow(df_plot) == 0) {
      stop("max_lag_norm = ", max_lag_norm, " leaves no lags to plot.")
    }
  }

  # Benchmark data
  df_benchmark <- df_plot %>%
    dplyr::group_by(lambda, lambda_label, lag_norm) %>%
    dplyr::summarise(
      benchmark_value = first(original_metric),
      .groups = "drop"
    ) %>%
    dplyr::mutate(estimator = "Naïve")

  m_values <- sort(unique(df_plot$m))

  # Corrected estimator data: one curve per m, or the best m per (lambda, lag)
  if (best_m_only) {
    df_corrected <- df_plot %>%
      dplyr::group_by(lambda, lambda_label, lag_norm) %>%
      dplyr::slice_min(order_by = abs(corrected_metric), n = 1, with_ties = FALSE) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(estimator = "Best m")
  } else {
    df_corrected <- df_plot %>%
      dplyr::mutate(estimator = paste0("m = ", m))
  }

  # Combine for unified legend
  df_combined <- dplyr::bind_rows(
    df_benchmark %>%
      dplyr::rename(y_value = benchmark_value) %>%
      dplyr::select(lambda, lambda_label, lag_norm, y_value, estimator),
    df_corrected %>%
      dplyr::rename(y_value = corrected_metric) %>%
      dplyr::select(lambda, lambda_label, lag_norm, y_value, estimator, m)
  )

  # Cressie-Hawkins data (if requested and available)
  if (has_CH) {
    df_CH <- df_plot %>%
      dplyr::group_by(lambda, lambda_label, lag_norm) %>%
      dplyr::summarise(
        CH_value = first(CH_metric),
        .groups = "drop"
      ) %>%
      dplyr::mutate(estimator = "CH")

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

    df_combined <- dplyr::bind_rows(
      df_combined,
      df_CH %>%
        dplyr::rename(y_value = CH_value) %>%
        dplyr::select(lambda, lambda_label, lag_norm, y_value, estimator)
    )
  }

  df_combined$metric       <- metric
  df_combined$metric_label <- cols$label

  list(df = df_combined, m_values = m_values, has_CH = has_CH,
       y_label = cols$label)
}

plot_metric_vs_norm_publication <- function(
    summary_df,
    N_filter,
    lambda_values = c(2, 4, 6, 8),
    metric = "bias",
    normalization = "Nh",
    max_m = NULL,
    max_lag_norm = NULL,  
    show_CH = FALSE,
    CH_color = "#E41A1C",
    naive_color = "black",
    out_dir = ".",
    filename = NULL,
    width = 7,
    height = 6,
    dpi = 600,  
    base_size = 11,
    font_scale = 1,
    line_size = 0.6,
    point_size = 1.8,
    show_points = TRUE,
    m_points = TRUE,
    m_line_size = NULL,
    m_point_size = NULL,
    naive_point_size = NULL,
    CH_point_size = NULL,
    m_show = NULL,
    label_m_extremes = FALSE,
    palette_begin = 0,
    palette_end = 0.9,
    use_bw = FALSE,
    show_title = FALSE,
    legend_position = "right",
    facet_ncol = 2,
    legend_nrow = NULL,
    palette_option = "viridis",
    CH_clip_factor = 1.5,
    stack_metrics = FALSE,
    stack_independent_y = FALSE
) {

  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  # Validate normalization parameter
  if (!normalization %in% c("Nh", "N")) {
    stop("normalization must be either 'Nh' or 'N'")
  }

  # Which metrics become rows: a single panel row (the classic layout) or one row
  # per metric stacked on top of each other, sharing the lambda columns.
  if (stack_metrics) {
    metrics_to_plot <- if (length(metric) > 1) metric else c("bias", "MSE")
  } else {
    metrics_to_plot <- metric[1]
  }

  # Compute scaled font sizes
  sz_axis_title  <- base_size * font_scale        # Axis names (x, y labels)
  sz_axis_text   <- (base_size - 1) * font_scale  # Axis tick labels (numbers/letters)
  sz_legend_title <- base_size * font_scale        # Legend title
  sz_legend_text  <- (base_size - 1) * font_scale  # Legend entries
  sz_strip_text   <- base_size * font_scale        # Facet strip labels
  
  #if (normalization == "N") y_label <- paste0(y_label, " (N-normalized)")

  # Filter data
  df_filtered <- summary_df %>%
    dplyr::filter(N == N_filter, lambda %in% lambda_values)
  
  # Filter by max_m if specified
  if (!is.null(max_m)) {
    df_filtered <- df_filtered %>%
      dplyr::filter(m <= max_m)
  }

  # Keep only selected m values. The Naive and CH columns are constant across m,
  # so this thins the m-family without touching the two reference curves.
  if (!is.null(m_show)) {
    df_filtered <- df_filtered %>%
      dplyr::filter(m %in% m_show)
    if (nrow(df_filtered) == 0) {
      stop("m_show selected no rows; available m: ",
           paste(sort(unique(summary_df$m)), collapse = ", "))
    }
  }

  # Resolve the per-series styling defaults
  if (is.null(m_line_size))      m_line_size      <- line_size
  if (is.null(m_point_size))     m_point_size     <- point_size
  if (is.null(naive_point_size)) naive_point_size <- point_size
  if (is.null(CH_point_size))    CH_point_size    <- point_size

  if (nrow(df_filtered) == 0) {
    stop("No data found for the specified N and lambda filters.")
  }
  
  # Prepare data: one long-format block per metric row, then stack them
  parts <- lapply(metrics_to_plot, function(met) {
    .norm_metric_long(df_filtered, met, normalization, lambda_values,
                      best_m_only = best_m_only, show_CH = show_CH,
                      CH_clip_factor = CH_clip_factor,
                      max_lag_norm = max_lag_norm)
  })
  names(parts) <- metrics_to_plot

  df_combined <- dplyr::bind_rows(lapply(parts, `[[`, "df"))
  has_CH      <- any(vapply(parts, `[[`, logical(1), "has_CH"))
  y_labels    <- vapply(parts, `[[`, character(1), "y_label")
  y_label     <- y_labels[[1]]

  # Metric rows in the requested order (labels are parseable, so label_parsed
  # renders them as plain text next to the lambda expressions)
  df_combined$metric_label <- factor(df_combined$metric_label, levels = y_labels)

  # Drop m values that have no finite entry anywhere in this figure. m = n - 1 makes
  # W*_m singular, so the summary carries it as all-NaN; keeping it would add a
  # legend entry (and a palette slot) for a curve that is never drawn.
  if (!best_m_only && any(!is.na(df_combined$m))) {
    finite_m <- df_combined %>%
      dplyr::filter(!is.na(m)) %>%
      dplyr::group_by(m) %>%
      dplyr::summarise(any_finite = any(is.finite(y_value)), .groups = "drop")
    drop_m <- finite_m$m[!finite_m$any_finite]
    if (length(drop_m)) {
      message("plot_metric_vs_norm_publication(): no finite values for m = ",
              paste(drop_m, collapse = ", "), " - dropping from the figure.")
      df_combined <- df_combined %>%
        dplyr::filter(is.na(m) | !m %in% drop_m)
    }
  }

  # m values present (placeholder level when only the best m is drawn)
  if (best_m_only) {
    m_values <- "best"
  } else {
    m_values <- sort(unique(df_combined$m[!is.na(df_combined$m)]))
  }

  # Create estimator factor with proper ordering
  if (best_m_only) {
    all_levels <- c("Naïve", "Best m")
  } else {
    all_levels <- c("Naïve", paste0("m = ", m_values))
  }
  if (has_CH) {
    all_levels <- c(all_levels, "CH")
  }
  df_combined$estimator <- factor(df_combined$estimator, levels = all_levels)
  df_combined$m_factor  <- df_combined$estimator   # kept for backward compatibility

  # Define aesthetics based on use_bw
  n_estimators <- length(all_levels)
  
  if (best_m_only) {
    # Simplified palette for best_m_only mode
    if (use_bw) {
      colors <- c(naive_color, "gray40")
      linetypes <- c("solid", "dashed")
      shapes <- c(16, 17)
    } else {
      mid <- (palette_begin + palette_end) / 2
      colors <- c(naive_color, .m_palette(palette_option, 1, mid, mid))
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
    # Colour version: viridis name, paletteer spec, or an explicit colour vector
    colors <- c(naive_color, .m_palette(palette_option, length(m_values),
                                        palette_begin, palette_end))
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
  
  # Panel columns: lambda values when stacking metrics, facet_ncol otherwise.
  # When stacking, the metric names label the rows, so the y axis title is dropped.
  n_facet_cols <- if (stack_metrics) length(lambda_values) else facet_ncol
  y_axis_label <- if (stack_metrics) NULL else y_label

  # With markers on the references only, blank the m-family shapes so the legend
  # keys show exactly what is drawn (a shape of NA draws no marker in the key)
  if (show_points && !m_points) {
    shapes[setdiff(all_levels, c("Naïve", "CH"))] <- NA
  }

  # Split the series into the ordered m-family and the two reference estimators, so
  # the references can be drawn last (on top) with their own weight and markers
  # instead of being buried under a dozen overlapping m-curves.
  ref_levels <- intersect(c("CH", "Naïve"), all_levels)   # CH first, Naïve on top
  df_family <- df_combined %>% dplyr::filter(!estimator %in% ref_levels)
  df_CH     <- df_combined %>% dplyr::filter(estimator == "CH")
  df_naive  <- df_combined %>% dplyr::filter(estimator == "Naïve")

  # Direct labels for the lowest and highest m curves, in the left margin. With an
  # ordered family, reading "which line is which" off the curve beats matching a
  # colour against a nine-entry legend. The left edge (lag norm 0) is where the
  # family fans out furthest, so labels there collide least.
  df_label <- NULL
  if (label_m_extremes && !best_m_only && nrow(df_family) &&
      any(!is.na(df_family$m))) {
    keep_m <- range(df_family$m, na.rm = TRUE)
    df_label <- df_family %>%
      # High m can be NaN at the shortest lags, and a label on a missing y is
      # silently dropped - anchor each label to the curve's first finite point
      dplyr::filter(m %in% keep_m, is.finite(y_value)) %>%
      dplyr::group_by(metric_label, lambda_label, estimator) %>%
      dplyr::slice_min(order_by = lag_norm, n = 1, with_ties = FALSE) %>%
      dplyr::ungroup()
  }

  # Build the plot
  p <- ggplot(df_combined, aes(x = lag_norm, y = y_value,
                               color = estimator,
                               linetype = estimator,
                               shape = estimator)) +
    # m-family: thin lines, optional (small) markers. ggplot2 >= 3.5 draws legend
    # keys per layer, so each key ends up showing exactly its own curve's weight.
    geom_line(data = df_family, aes(group = estimator), linewidth = m_line_size) +
    (if (show_points && m_points)
      geom_point(data = df_family, size = m_point_size) else NULL) +
    # Reference estimators, each its own layer so line weight and marker size are
    # tunable per series. CH first, Naive last, so Naive ends up on top.
    (if (nrow(df_CH)) geom_line(data = df_CH, linewidth = line_size) else NULL) +
    (if (show_points && nrow(df_CH))
      geom_point(data = df_CH, size = CH_point_size) else NULL) +
    geom_line(data = df_naive, linewidth = line_size) +
    (if (show_points) geom_point(data = df_naive, size = naive_point_size) else NULL) +
    # Direct labels on the outermost m-curves, just left of the panel content
    (if (!is.null(df_label))
      geom_text(data = df_label, aes(label = paste0("m=", m)),
                hjust = 1.15, vjust = 0.4, show.legend = FALSE,
                size = (base_size - 3) * font_scale / .pt) else NULL) +
    # Scales. `breaks` pins the legend order to Naive -> m ... -> CH; without it the
    # keys would follow the layer order above (family first) instead.
    scale_color_manual(values = colors, breaks = all_levels, name = "Estimator") +
    scale_linetype_manual(values = linetypes, breaks = all_levels, name = "Estimator") +
    scale_shape_manual(values = shapes, breaks = all_levels, name = "Estimator") +
    # Labels
    labs(
      x = "Lag norm",
      y = y_axis_label
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

  # Room on the left for the direct labels so they are not clipped
  if (!is.null(df_label)) {
    p <- p + scale_x_continuous(expand = expansion(mult = c(0.16, 0.03)))
  }

  # Facets with parsed labels (for Greek lambda)
  if (stack_metrics) {
    if (stack_independent_y) {
      # facet_wrap over metric x lambda: every panel gets its own y range, at the
      # price of repeating the metric name in each panel strip
      p <- p + facet_wrap(~ metric_label + lambda_label, ncol = n_facet_cols,
                          scales = "free_y", labeller = label_parsed)
    } else {
      # facet_grid: metric rows, lambda columns. The row strip sits on the left of
      # the axis and so reads as the y-axis title; the y range is shared within a row.
      p <- p + facet_grid(metric_label ~ lambda_label, scales = "free_y",
                          switch = "y",
                          labeller = labeller(lambda_label = label_parsed,
                                              metric_label = label_value)) +
        theme(strip.placement = "outside",
              strip.text.y.left = element_text(size = sz_strip_text, angle = 90,
                                               margin = margin(r = 5)))
    }
  } else {
    p <- p + facet_wrap(~ lambda_label, ncol = n_facet_cols, scales = "free_y",
                        labeller = label_parsed)
  }

  # Add horizontal line at 0 for bias (only in the bias row when metrics are stacked)
  if ("bias" %in% metrics_to_plot) {
    if (stack_metrics) {
      df_zero <- data.frame(
        metric_label = factor("Bias", levels = levels(df_combined$metric_label)),
        y_zero = 0
      )
      p <- p + geom_hline(data = df_zero, aes(yintercept = y_zero),
                          inherit.aes = FALSE, linetype = "dotted",
                          color = "gray50", linewidth = 0.4)
    } else {
      p <- p + geom_hline(yintercept = 0, linetype = "dotted",
                          color = "gray50", linewidth = 0.4)
    }
  }

  # Add title if requested
  if (show_title) {
    title_metric <- paste(y_labels, collapse = " and ")
    p <- p + ggtitle(bquote(.(title_metric) ~ "vs lag norm for" ~ N == .(N_filter) %*% .(N_filter)))
  }

  # Adjust legend layout based on position
  if (legend_position == "bottom") {
    # A single row overflows when many m-levels are shown (e.g. all m plus
    # Naïve/CH on a 4-panel-wide layout). Wrap the legend over several rows so
    # it stays within the figure width. Aim for ~ (#facets + 1) entries per row.
    if (!is.null(legend_nrow)) {
      n_rows <- legend_nrow
    } else {
      per_row <- n_facet_cols + 1
      n_rows <- max(1, ceiling(n_estimators / per_row))
    }
    p <- p + guides(
      color = guide_legend(nrow = n_rows, byrow = TRUE),
      linetype = guide_legend(nrow = n_rows, byrow = TRUE),
      shape = guide_legend(nrow = n_rows, byrow = TRUE)
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
    metric_tag <- paste0(paste(metrics_to_plot, collapse = "_"), "_", normalization)
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

  if (interactive()) print(p)

  invisible(list(
    plot = p,
    png_file = out_path_png,
    data = df_combined
  ))
}

#------------------------------------------------------------------------------
# MSE ratio heatmap grid (corrected / naive, faceted by lambda)
#------------------------------------------------------------------------------
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
    ncol = 2,
    value = c("ratio", "pct"),
    axes_all = FALSE,
    panel_axis_titles = FALSE
) {
  # value    "ratio" fills with MSEChatCorr / MSEChat, centred on 1 (the
  #          historical behaviour). "pct" fills with the signed percent change
  #          (ratio - 1) * 100, centred on 0, so a REDUCTION carries a minus
  #          sign in the legend exactly as it already does in the cell labels.
  # axes_all TRUE draws the h1 / h2 tick labels on every panel instead of only
  #          on the outer margins.
  # panel_axis_titles TRUE gives every panel its own "h1" / "h2" axis TITLES.
  #          A facetted plot has one pair of titles for the whole figure, so
  #          this builds one plot per lambda and lays them out with patchwork;
  #          the result is a patchwork object, still printable and ggsave-able.
  value <- match.arg(value)
  if (panel_axis_titles && !requireNamespace("patchwork", quietly = TRUE))
    stop("panel_axis_titles = TRUE needs the 'patchwork' package: ",
         "install.packages(\"patchwork\").", call. = FALSE)

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
    val_tag <- if (value == "pct") "_pct" else ""
    filename <- sprintf("mse_ratio_heatmap_grid_%s%s%s_N%d.png", m_tag, norm_tag,
                        val_tag, N_filter)
  }
  
  # The fill is the ratio or the signed percent change; the cell labels are the
  # signed percent change either way.
  fill_var <- if (value == "pct") "pct_change" else "mse_ratio"

  # Fill scale: on the ratio scale the neutral point is 1, on the percent scale
  # it is 0 and a reduction is negative. A fresh scale object per panel, so the
  # per-lambda plots never share one.
  if (value == "pct") {
    fill_limits <- (scale_limits - 1) * 100
    fill_mid    <- 0
    fill_name   <- "MSE change (%)"
    fill_labels <- function(x) sprintf("%+.0f%%", x)
  } else {
    fill_limits <- scale_limits
    fill_mid    <- 1
    fill_name   <- "MSE Ratio"
    fill_labels <- ggplot2::waiver()
  }
  fill_scale <- function() {
    if (palette_option == "diverging") {
      ggplot2::scale_fill_gradient2(
        low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
        midpoint = fill_mid, limits = fill_limits, labels = fill_labels,
        name = fill_name, na.value = "grey70", oob = scales::squish)
    } else {
      ggplot2::scale_fill_viridis_c(
        option = palette_option, direction = palette_direction,
        limits = fill_limits, labels = fill_labels,
        name = fill_name, na.value = "grey70", oob = scales::squish)
    }
  }

  base_theme <- ggplot2::theme_minimal(base_size = 11) +
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

  # One heatmap, over whatever subset it is handed. Used for the whole grid
  # (then facetted) and for a single lambda (then laid out by patchwork).
  build_panel <- function(dat, cross_dat) {
    q <- ggplot2::ggplot(dat, ggplot2::aes(x = factor(h1), y = factor(h2),
                                           fill = .data[[fill_var]])) +
      ggplot2::geom_tile(color = "white", linewidth = 0.2)
    if (show_numbers) {
      q <- q + ggplot2::geom_text(
        ggplot2::aes(label = label, color = I(text_col)), size = label_size)
    }
    # X marks where the zero estimator beats the best available
    if (nrow(cross_dat) > 0) {
      q <- q +
        ggplot2::geom_segment(
          data = cross_dat,
          ggplot2::aes(x = x_num - 0.4, xend = x_num + 0.4,
                       y = y_num - 0.4, yend = y_num + 0.4),
          inherit.aes = FALSE, linewidth = 0.5, color = "black") +
        ggplot2::geom_segment(
          data = cross_dat,
          ggplot2::aes(x = x_num - 0.4, xend = x_num + 0.4,
                       y = y_num + 0.4, yend = y_num - 0.4),
          inherit.aes = FALSE, linewidth = 0.5, color = "black")
    }
    q <- q + fill_scale() +
      ggplot2::labs(x = expression(h[1]), y = expression(h[2])) +
      base_theme
    if (fixed_aspect) q <- q + ggplot2::coord_fixed()
    q
  }

  if (panel_axis_titles) {
    # One plot per lambda, so each keeps its own h1 / h2 axis titles; the lambda
    # label becomes the panel title and the shared legend is collected once.
    lvls <- levels(droplevels(df_plot$lambda_label))
    panels <- lapply(lvls, function(lv) {
      build_panel(df_plot[df_plot$lambda_label == lv, , drop = FALSE],
                  df_cross[df_cross$lambda_label == lv, , drop = FALSE]) +
        ggplot2::ggtitle(parse(text = lv)[[1]]) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5,
                                                          size = 12))
    })
    p <- patchwork::wrap_plots(panels, ncol = ncol, guides = "collect")
  } else {
    # axes = "all" repeats the h1 / h2 tick labels (not titles) on every panel.
    p <- build_panel(df_plot, df_cross) +
      ggplot2::facet_wrap(~ lambda_label, ncol = ncol,
                          labeller = ggplot2::label_parsed,
                          axes = if (axes_all) "all" else "margins")
  }

  # Save (PNG only -- this project never produces PDF plots)
  out_path_png <- file.path(out_dir, sub("\\.pdf$", ".png", filename))

  ggplot2::ggsave(out_path_png, plot = p, width = width, height = height, dpi = dpi, bg = "white")

  message("Saved: ", out_path_png)

  if (interactive()) print(p)

  invisible(list(plot = p, png_file = out_path_png, data = df_plot))
}


#------------------------------------------------------------------------------
# Best m heatmap with axis thinning and zero-estimator flags
#------------------------------------------------------------------------------
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
