# For "scenario = 5", run "application.R" until line 210

library("patchwork")

USA <- readRDS(file = "APPLICATION/USA.rds")
USA <- USA[!(USA$shapeName %in% c("Alaska", "American Samoa","Commonwealth of the Northern Mariana Islands", "Guam", "Hawaii", "Puerto Rico", "United States Virgin Islands")), ]
USA <- st_transform(x = USA, crs = "+init=epsg:6345 +units=km +no_defs")
USA_all <- st_union(USA$geom)

bfs <- bfs_pred
coord <- coord_pred

n_base_functions <- length(bfs)

bfs_all <- list()
for (i in 1:n_base_functions) {
  bfs_all[[i]] <- data.frame(x = coord$x[1:n_coord_pr], y = coord$y[1:n_coord_pr], z = bfs[[i]][1:n_coord_pr])
  
  coordinates(bfs_all[[i]]) <- ~ x + y
  gridded(bfs_all[[i]]) <- TRUE
  
  bfs_all[[i]] <- raster(bfs_all[[i]])
}

pps <- list()
for (i in 1:n_base_functions) {
  crs(bfs_all[[i]]) <- "+init=epsg:6345 +units=km +no_defs"
  bfs_plot <- as(bfs_all[[i]], "SpatialPixelsDataFrame")
  bfs_plot <- as.data.frame(bfs_plot)
  colnames(bfs_plot) <- c("value", "x", "y")
  
  pal <- jet.col(n = 100, alpha = 1)
  labs <- seq(0, 1, length.out = 6)
  
  if (i %in% c(13, 14, 15)     ) { x_lab = "" } else { x_lab = "" }
  if (i %in% c(1, 4, 7, 10, 13)) { y_lab = ""  } else { y_lab = "" }
  if (i %in% c(13, 14, 15)) { guide = "colourbar" } else { guide = "none" }
  if (i %in% c(13, 14, 15)) { plot_margin = margin(0, 0, 0, 0, "cm") } else { plot_margin = margin(0, 0, 2.9, 0, "cm") }
  
  pps[[i]] <- ggplot() +
    geom_tile(data = bfs_plot, aes(x = x, y = y, fill = value)) +
    geom_sf(data = USA_all, color = "black", fill = NA, lwd = 0.5) +
    scale_fill_gradientn(guide = guide, 
                         name = bquote(phi*"(x)"), 
                         colors = pal, breaks = labs, 
                         labels = as.character(format(labs, nsmall = 1)), 
                         limits = c(labs[1], tail(labs, 1))) +
    labs(x = x_lab, y = y_lab, title ="") + 
    theme_bw() +
    theme(text = element_text(size = 30, family = "LM Roman 10"), 
          legend.key.height = unit(1.0, "cm"),
          legend.key.width  = unit(2.9, "cm"),
          plot.margin = plot_margin,
          plot.title = element_text(hjust = 0.5),
          plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank()) +
    coord_sf(ndiscr = 0)   
  
  # 
}

pp <- wrap_elements(pps[[ 1]] + pps[[ 2]] + pps[[ 3]] + plot_layout(widths = c(1, 1, 1))) / 
      wrap_elements(pps[[ 4]] + pps[[ 5]] + pps[[ 6]] + plot_layout(widths = c(1, 1, 1))) /
      wrap_elements(pps[[ 7]] + pps[[ 8]] + pps[[ 9]] + plot_layout(widths = c(1, 1, 1))) /
      wrap_elements(pps[[10]] + pps[[11]] + pps[[12]] + plot_layout(widths = c(1, 1, 1))) /
      wrap_elements(pps[[13]] + pps[[14]] + pps[[15]] + 
                      plot_layout(widths = c(1, 1, 1), guides = "collect") & 
                      scale_fill_gradientn(name = "",
                                           limits = c(labs[1], tail(labs, 1)),
                                           colors = pal, breaks = labs, 
                                           labels = as.character(format(labs, nsmall = 1))) & theme(legend.position = "bottom"),
                    )

ggsave(filename = paste("APPLICATION/IMAGES/phi.jpeg", sep = ""), plot = pp, width = 8000, height = 12000, units = c("px"), dpi = 300, bg = "white")
  
  