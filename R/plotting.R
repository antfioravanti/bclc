if (!require(RColorBrewer)) install.packages("RColorBrewer"); library(RColorBrewer)
#-------------------------------------------------------------------------------
# CUSTOM PALETTES
palette_custom <- list(
  # Core colors (ordered for nice gradients)
  trr          = rgb(113, 173, 69, maxColorValue = 255),   # Green
  stairs       = rgb(123, 193, 187, maxColorValue = 255),  # Teal
  boxcolor     = rgb(186, 215, 230, maxColorValue = 255),  # Light blue
  white        = rgb(255, 255, 255, maxColorValue = 255),  # White
  areaa_light  = rgb(251, 220, 212, maxColorValue = 255),  # Light orange
  areaa        = rgb(234, 84, 40, maxColorValue = 255)     # Orange/red
)


# Palette 2: Diverging color scheme (unchanged)
palette2 <- c(
  "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#F7F7F7",
  "white",
  "#FEE0D2", "#FCBBA1", "#FC9272", "#FB6A4A", "#EF3B2C", "#CB181D"
)

#-------------------------------------------------------------------------------
# PLOTTING
#-------------------------------------------------------------------------------
plot_matrix = function(X,
                       main = "Matrix", 
                       labels = FALSE,
                       show_cell_num = FALSE,
                       show_index = FALSE, 
                       show_M_num = FALSE, 
                       show_values = FALSE,
                       col_palette = NA,
                       show_legend = FALSE,
                       asp = 1,
                       mar = c(1, 1, 2, 1),
                       oma = c(0, 0, 0, 0),
                       out_dir = NULL,
                       filename = NULL,
                       width = 7,
                       height = 7,
                       dpi = 300,
                       file_type = "png") {
  
  do_plot = function() {
    old_par = par(no.readonly = TRUE)
    on.exit(par(old_par))
    
    # If legend requested, add extra right margin
    plot_mar = mar
    if (show_legend) {
      plot_mar[4] = max(mar[4], 5)
    }
    
    par(mar = plot_mar, oma = oma, pty = if (asp == 1) "s" else "m")
    
    n1 = nrow(X)
    n2 = ncol(X)
    
    if (length(col_palette) == 1 && is.na(col_palette)) {
      image(t(X)[, nrow(X):1],
            xaxt = "n", yaxt = "n",
            main = main,
            asp = asp)
    } else {
      image(t(X)[, nrow(X):1],
            col = col_palette,
            xaxt = "n", yaxt = "n",
            main = main,
            asp = asp)
    }
    
    if (labels) {
      axis(1, at = seq(0, 1, length.out = ncol(X)), labels = 1:ncol(X))
      axis(2, at = seq(0, 1, length.out = nrow(X)), labels = nrow(X):1)
    }
    
    if (show_values) {
      for (i in seq_len(nrow(X))) {
        for (j in seq_len(ncol(X))) {
          x_coord = (j - 1) / (ncol(X) - 1)
          y_coord = 1 - (i - 1) / (nrow(X) - 1)
          val_str = sprintf("%.4f", X[i, j])
          text(x_coord, y_coord, val_str, cex = 0.7)
        }
      }
    }
    
    if (show_legend) {
      if (!requireNamespace("fields", quietly = TRUE)) {
        stop("Install 'fields' for legend support.")
      }
      fields::image.plot(legend.only = TRUE,
                         zlim = range(X, na.rm = TRUE),
                         col = col_palette,
                         legend.width = 1.2,
                         legend.shrink = 0.8,
                         legend.mar = 4.5)
    }
  }
  
  do_plot()
  
  if (!is.null(out_dir) && !is.null(filename)) {
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    file_path <- file.path(out_dir, paste0(filename, ".", file_type))
    
    if (file_type == "png") {
      png(file_path, width = width, height = height, units = "in", res = dpi)
    } else if (file_type == "pdf") {
      pdf(file_path, width = width, height = height)
    } else {
      stop("file_type must be 'png' or 'pdf'")
    }
    
    do_plot()
    dev.off()
  }
}



plot_matrix_grayscale = function(X){
  nvec = dim(X)
  return(filled.contour(x=1:nvec[1], y=1:nvec[2],
                        X, color.palette=gray.colors))
}
plot_3D_matrix = function(X){
  nvec = dim(X)
  return(persp(x=(1:nvec[1]), y=(1:nvec[2]), X,
               theta=45, phi=35, r=5, expand=0.6, axes=T,
               ticktype="detailed", xlab="t1 ", ylab="t2", zlab="X_t1,t2"))
}
