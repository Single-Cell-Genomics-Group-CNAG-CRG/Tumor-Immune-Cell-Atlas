plot_spaniel <- function(data_df, 
                         grob, 
                         x, 
                         y, 
                         point_colour, 
                         point_size,
                         point_alpha){
  
  # Inverse Y to flip the coordinates
  # data_df$y_inv <- 36 - data_df$y
  # data_df$y_inv <- data_df$y
  
  data_df[, point_size] <- if_else(data_df[, point_size] == 0, NA_real_, data_df[, point_size])
  
  tmp_plt <- ggplot2::ggplot(data_df,
                  ggplot2::aes_string(x, y,
                                      color = point_colour, 
                                      size = point_size,
                                      alpha = point_alpha)) +
    ggplot2::xlim(1, 36) +
    ggplot2::ylim(1, 36) +
    # Layer 1 - Plot image
    ggplot2::annotation_custom(grob,
                               xmin = 1,
                               xmax = 36,
                               ymin = 1,
                               ymax = 36) +
    # Layer 2 - Plot points
    geom_point() +
    # Layer 3 - Join legends all with same name
    labs(color = "Proportion",
         size = "Proportion"
         # alpha = "Proportion"
    ) +
    # Coordinates fixed so x and y keep the proportions
    coord_fixed(1) +
    # Tune colour parameters
    ggplot2::scale_size_continuous(range=c(0, 3), limits = c(0, 1)) +
    ggplot2::scale_color_gradientn(
      colours = heat.colors(10, rev = TRUE),
      limits = c(0, 1)) +
    ggplot2::scale_alpha(range = c(0, 1), limits = c(0, 1)) +
    # Join legends into one
    ggplot2::guides(color = ggplot2::guide_legend(), 
                    size = ggplot2::guide_legend(),
                    alpha = ggplot2::guide_legend()
    )
  
  return(tmp_plt)
}


################################################################################
################################################################################
############################Join Spatial plots##################################
################################################################################
################################################################################

plt_theme <- function(plt, rm_leg) {
  tmp_plt <- plt + 
    scale_fill_gradientn(
      colours = heat.colors(10, rev = TRUE),
      limits = c(0, 1)) +
    scale_alpha(range = c(0, 1)) +
    labs(title = "",
         fill = "Proportion") +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"),
          plot.title = element_text(hjust = 0.5))
  
  if (rm_leg) {
    tmp_plt + theme(legend.position = "none")
  } else {
    tmp_plt + theme(legend.position = "right")
  }
}

draw_title <- function(feat,
                       title_size = 20,
                       title_face = "bold",
                       title_margin_t = 0,
                       title_margin_r = 0,
                       title_margin_b = 0,
                       title_margin_l = 10) {
  title <- cowplot::ggdraw() + 
    cowplot::draw_label(
      feat,
      fontface = 'bold',
      x = 0,
      hjust = 0
    ) +
    ggplot2::theme(
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(0, 0, 0, 7)
    )
  return(title)
}


join_seurat_spatial <- function(se_obj,
                                plt,
                                grp,
                                feat,
                                title_size = 20,
                                title_face = "bold",
                                title_margin_t = 0,
                                title_margin_r = 0,
                                title_margin_b = 0,
                                title_margin_l = 10) {
  
  # Determine which groups are 0
  metadata_ds <- se_obj@meta.data
  slices <- unique(metadata_ds[, grp])
  grp1 <- sum(metadata_ds[metadata_ds$slice == slices[[1]], feat])
  grp2 <- sum(metadata_ds[metadata_ds$slice == slices[[2]], feat])
  
  # Process individual pannels
  plt_1 <- suppressMessages(plt_theme(plt = plt[[1]], rm_leg = TRUE))
  ### Remove dots if 0s everywhere
  if (grp1 == 0) {
    plt_1 <- suppressMessages(plt_1 + scale_alpha(range = c(0,0)))
  }
  
  plt_2 <- suppressMessages(plt_theme(plt = plt[[2]], rm_leg = FALSE))
  ### Remove dots if 0s everywhere
  if (grp2 == 0) {
    plt_2 <- suppressMessages(plt_2 + scale_alpha(range = c(0,0)))
  }
  
  # Recombine plt 1 and 2
  tmp_plt <- cowplot::plot_grid(plt_1,
                                plt_2,
                                align = "hv",
                                axis = "trbl",
                                nrow = 1)
  
  # Add title to plot composition
  title <- draw_title(feat = feat,
                      title_size = title_size,
                      title_face = title_face,
                      title_margin_t = title_margin_t,
                      title_margin_r = title_margin_r,
                      title_margin_b = title_margin_b,
                      title_margin_l = title_margin_l)
  
  plt_arr <- cowplot::plot_grid(
    title, tmp_plt,
    ncol = 1,
    # rel_heights values control vertical title margins
    rel_heights = c(0.1, 1)
  )
  
  return(plt_arr)
  
}

