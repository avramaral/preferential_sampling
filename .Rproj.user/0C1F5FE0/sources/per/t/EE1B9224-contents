# args <- commandArgs(trailingOnly = TRUE)
# scenario <- as.numeric(args[1]) 

library("plot3D")
library("tidyverse")
library("patchwork")

source("header.R")
source("functions_data_generation.R")
source("functions_model_fitting.R")
source("functions_result_processing.R")
source("functions_aux.R")

scenario <- 5

data <- readRDS(file = "APPLICATION/data.rds")
result_obj <- readRDS(file = paste("APPLICATION/FITTED_MODELS/scenario", scenario, ".rds", sep = ""))
USA <- readRDS(file = "APPLICATION/USA.rds")

##################################################
##################################################
USA <- USA[!(USA$shapeName %in% c("Alaska", "American Samoa","Commonwealth of the Northern Mariana Islands", "Guam", "Hawaii", "Puerto Rico", "United States Virgin Islands")), ]
USA <- st_transform(x = USA, crs = "+init=epsg:6345 +units=km +no_defs")
##################################################
##################################################

result <- result_obj$result
bfs_pred <- result_obj$bfs_pred
fs_formula <- result_obj$fs_formula
coord_pred <- result_obj$coord_pred
n_coord_pr <- result_obj$n_coord_pr

idx_lat <- inla.stack.index(fs_formula$full_stack, tag = "pred_y")$data[(n_coord_pr + 1):nrow(coord_pred)]
fitted_latent <- cbind(coord_pred[(n_coord_pr + 1):nrow(coord_pred), ], result$fit$summary.fitted.values[idx_lat, c("mean")])
colnames(fitted_latent) <- c("x", "y", "z")
data$z <- fitted_latent$z

##################################################
# Estimated process in the observation locations #
##################################################
pal <- jet.col(n = 100, alpha = 0.9)
labs <- seq(round(min(data$z) - 0.01, 2), round(max(data$z) + 0.01, 2), length.out = 6)

p <- ggplot() +
  geom_sf(data = USA, color = "black", fill = NA, lwd = 0.5) +
  geom_sf(data = data, aes(color = z), size = 5) + 
  scale_colour_gradientn(name = expression(paste(PM[2.5], " level", sep = "")), colors = pal, breaks = labs, labels = as.character(format(labs, nsmall = 3)), limits = c(labs[1], tail(labs, 1))) +
  labs(x = "Longitude", y = "Latitude", title = "") + 
  theme_bw() +
  theme(text = element_text(size = 24, family = "LM Roman 10"), 
        legend.key.height = unit(2.57, "cm"),
        legend.key.width  = unit(0.75, "cm"),
        plot.margin = margin(0, 0, 0, 0, "cm"))
ggsave(filename = paste("APPLICATION/IMAGES/", sprintf("%02d", scenario), "/estim_obs.jpeg", sep = ""), plot = p, width = 3250, height = 2000, units = c("px"), dpi = 300, bg = "white")
##################################################
##################################################

ss <- (result$fit$summary.fitted.values[idx_lat, c("sd")] ^ 2) + (1 / result$fit$summary.hyperpar["Precision for the Gaussian observations", "mean"])
CRPS_latent  <- mean(CRPS(orig = data$mean, pred = fitted_latent$z, s = sqrt(ss)))
SCRPS_latent <- mean(SCRPS(orig = data$mean, pred = fitted_latent$z, s = sqrt(ss)))
error_latent <- MSE(orig = data$mean, pred = fitted_latent$z)

scores <- list(error = error_latent,
               CRPS  = CRPS_latent,
               SCRPS = SCRPS_latent)

saveRDS(object = scores, file =  paste("APPLICATION/IMAGES/", sprintf("%02d", scenario), "/scores.rds", sep = ""))

if (TRUE) {
  print(paste("Scenario ", scenario, sep = ""))
  print(paste("Error (latent): ", round(x = error_latent, digits = 4), sep = ""))
  print(paste("CRPS  (latent): ", round(x = CRPS_latent,  digits = 4), sep = ""))
  print(paste("SCRPS (latent): ", round(x = SCRPS_latent, digits = 4), sep = ""))
}

##################################################
##################################################

#####################
# Intensity Process #
#####################

if (scenario != 1) {
  idx_int <- inla.stack.index(fs_formula$full_stack, tag = "pred_pp")$data[1:n_coord_pr]
  fitted_intensity <- cbind(coord_pred[1:n_coord_pr, ], exp(result$fit$summary.fitted.values[idx_int, c("mean")]))
  colnames(fitted_intensity) <- c("x", "y", "estimated")
  
  # Create a gridded spatial object from "process"
  coordinates(fitted_intensity) <- ~ x + y
  gridded(fitted_intensity) <- TRUE
  fitted_intensity <- raster(fitted_intensity)
  crs(fitted_intensity) = "+init=epsg:6345 +units=km +no_defs"
  
  fitted_intensity    <- as(fitted_intensity, "SpatialPixelsDataFrame")
  fitted_intensity_df <- as.data.frame(fitted_intensity)
  colnames(fitted_intensity_df) <- c("estimated", "x", "y")
  
  pal <- jet.col(n = 100, alpha = 0.9)
  labs <- seq(round(min(fitted_intensity_df$estimated) - 0.000001, 6), round(max(fitted_intensity_df$estimated) + 0.000001, 6), length.out = 6)
  
  p_1 <- ggplot() +
    geom_tile(data = fitted_intensity_df, mapping = aes(x = x, y = y, fill = estimated)) + 
    geom_sf(data = USA, color = "black", fill = NA, lwd = 0.5) +
    scale_fill_gradientn(name = "Intensity", colors = pal, breaks = labs, labels = as.character(format(labs, nsmall = 3)), limits = c(labs[1], tail(labs, 1))) +
    labs(x = "Longitude", y = "Latitude", title = "") + 
    theme_bw() +
    theme(text = element_text(size = 24, family = "LM Roman 10"), 
          legend.key.height = unit(2.57, "cm"),
          legend.key.width  = unit(0.75, "cm"),
          plot.margin = margin(0, 0, 0, 0, "cm"))
  ggsave(filename = paste("APPLICATION/IMAGES/", sprintf("%02d", scenario), "/intensity.jpeg", sep = ""), plot = p_1, width = 3250, height = 2000, units = c("px"), dpi = 300, bg = "white")
}

##################
# Latent Process #
##################

idx_lat <- inla.stack.index(fs_formula$full_stack, tag = "pred_y")$data[1:n_coord_pr]
fitted_latent <- cbind(coord_pred[1:n_coord_pr, ], result$fit$summary.fitted.values[idx_lat, c("mean")])
colnames(fitted_latent) <- c("x", "y", "estimated")

# Create a gridded spatial object from "process"
coordinates(fitted_latent) <- ~ x + y
gridded(fitted_latent) <- TRUE
fitted_latent <- raster(fitted_latent)
crs(fitted_latent) = "+init=epsg:6345 +units=km +no_defs"

fitted_latent    <- as(fitted_latent, "SpatialPixelsDataFrame")
fitted_latent_df <- as.data.frame(fitted_latent)
colnames(fitted_latent_df) <- c("estimated", "x", "y")

pal <- jet.col(n = 100, alpha = 0.9)
labs <- seq(round(min(fitted_latent_df$estimated) - 0.01, 2), round(max(fitted_latent_df$estimated) + 0.01, 2), length.out = 6)
# labs <- seq(2.28, 14.68, length.out = 6)

p_2 <- ggplot() +
  geom_tile(data = fitted_latent_df, mapping = aes(x = x, y = y, fill = estimated)) + 
  geom_sf(data = USA, color = "black", fill = NA, lwd = 0.5) +
  scale_fill_gradientn(name = expression(paste(PM[2.5], " level", sep = "")), colors = pal, breaks = labs, labels = as.character(format(labs, nsmall = 3)), limits = c(labs[1], tail(labs, 1))) +
  labs(x = "Longitude", y = "Latitude", title = "") + 
  theme_bw() +
  theme(text = element_text(size = 24, family = "LM Roman 10"), 
        legend.key.height = unit(2.57, "cm"),
        legend.key.width  = unit(0.75, "cm"),
        plot.margin = margin(0, 0, 0, 0, "cm"))
ggsave(filename = paste("APPLICATION/IMAGES/", sprintf("%02d", scenario), "/latent.jpeg", sep = ""), plot = p_2, width = 3250, height = 2000, units = c("px"), dpi = 300, bg = "white")

p_2_loc <- p_2 + geom_sf(data = data, color = rgb(red = 0, green = 0, blue = 0, alpha = 1), size = 3, shape = 4)
ggsave(filename = paste("APPLICATION/IMAGES/", sprintf("%02d", scenario), "/usa_pm25.jpeg", sep = ""), plot = p_2_loc, width = 3250, height = 2000, units = c("px"), dpi = 300, bg = "white")

# 95% equal-tailed credible interval
if (TRUE) {
  idx_lat <- inla.stack.index(fs_formula$full_stack, tag = "pred_y")$data[1:n_coord_pr]
  fitted_latent_0025 <- cbind(coord_pred[1:n_coord_pr, ], result$fit$summary.fitted.values[idx_lat, c("0.025quant")])
  fitted_latent_0975 <- cbind(coord_pred[1:n_coord_pr, ], result$fit$summary.fitted.values[idx_lat, c("0.975quant")])
  colnames(fitted_latent_0025) <- c("x", "y", "estimated")
  colnames(fitted_latent_0975) <- c("x", "y", "estimated")
  
  # Create a gridded spatial object from "process"
  coordinates(fitted_latent_0025) <- ~ x + y
  gridded(fitted_latent_0025) <- TRUE
  fitted_latent_0025 <- raster(fitted_latent_0025)
  crs(fitted_latent_0025) = "+init=epsg:6345 +units=km +no_defs"
  
  coordinates(fitted_latent_0975) <- ~ x + y
  gridded(fitted_latent_0975) <- TRUE
  fitted_latent_0975 <- raster(fitted_latent_0975)
  crs(fitted_latent_0975) = "+init=epsg:6345 +units=km +no_defs"
  
  
  fitted_latent_0025    <- as(fitted_latent_0025, "SpatialPixelsDataFrame")
  fitted_latent_0025_df <- as.data.frame(fitted_latent_0025)
  colnames(fitted_latent_0025_df) <- c("estimated", "x", "y")
  
  fitted_latent_0975    <- as(fitted_latent_0975, "SpatialPixelsDataFrame")
  fitted_latent_0975_df <- as.data.frame(fitted_latent_0975)
  colnames(fitted_latent_0975_df) <- c("estimated", "x", "y")
  
  
  pal <- jet.col(n = 100, alpha = 0.9)
  labs <- seq(round(min(c(fitted_latent_0025_df$estimated, fitted_latent_0975_df$estimated)) - 0.01, 2), round(max(c(fitted_latent_0025_df$estimated, fitted_latent_0975_df$estimated)) + 0.01, 2), length.out = 6)
  
  p_0025 <- ggplot() +
    geom_tile(data = fitted_latent_0025_df, mapping = aes(x = x, y = y, fill = estimated)) + 
    geom_sf(data = USA, color = "black", fill = NA, lwd = 0.5) +
    scale_fill_gradientn(name = expression(paste(PM[2.5], " level", sep = "")), colors = pal, breaks = labs, labels = as.character(format(labs, nsmall = 3)), limits = c(labs[1], tail(labs, 1))) +
    labs(x = "Longitude", y = "Latitude", title = expression(paste("2.5" ^ "th", " percentile", sep = ""))) + 
    theme_bw() +
    theme(text = element_text(size = 24, family = "LM Roman 10"), 
          legend.key.height = unit(2.57, "cm"),
          legend.key.width  = unit(0.75, "cm"),
          plot.margin = margin(0, 1, 0, 0, "cm"))
  
  p_0975 <- ggplot() +
    geom_tile(data = fitted_latent_0975_df, mapping = aes(x = x, y = y, fill = estimated)) + 
    geom_sf(data = USA, color = "black", fill = NA, lwd = 0.5) +
    scale_fill_gradientn(name = expression(paste(PM[2.5], " level", sep = "")), colors = pal, breaks = labs, labels = as.character(format(labs, nsmall = 3)), limits = c(labs[1], tail(labs, 1))) +
    labs(x = "Longitude", y = "Latitude", title = expression(paste("97.5" ^ "th", " percentile", sep = ""))) +
    theme_bw() +
    theme(text = element_text(size = 24, family = "LM Roman 10"), 
          legend.key.height = unit(2.57, "cm"),
          legend.key.width  = unit(0.75, "cm"),
          plot.margin = margin(0, 0, 0, 0, "cm"))

  
  p_95 <- p_0025 + p_0975 & theme(legend.position = "right", legend.box.margin = margin(0, 0, 0, 0.5, "cm"))
  p_95 <- p_95 + plot_layout(guides = "collect")
  
  ggsave(filename = paste("APPLICATION/IMAGES/", sprintf("%02d", scenario), "/latent_credible_interval.jpeg", sep = ""), plot = p_95, width = 6000, height = 2000, units = c("px"), dpi = 300, bg = "white")
}


###########################
# Preferentiality Process #
###########################

if (!scenario %in% 1:2) {
  fitted_preferentiality <- compute_preferentiatility(fit = result$fit, coord_pred = coord_pred, bfs_pred = bfs_pred, xlim = xlim, by = by, n_coord_pr = n_coord_pr)
  crs(fitted_preferentiality) = "+init=epsg:6345 +units=km +no_defs"
  
  fitted_preferentiality    <- as(fitted_preferentiality, "SpatialPixelsDataFrame")
  fitted_preferentiality_df <- as.data.frame(fitted_preferentiality)
  colnames(fitted_preferentiality_df) <- c("estimated", "x", "y")
  
  pal <- jet.col(n = 100, alpha = 0.9)
  labs <- seq(round(min(fitted_preferentiality_df$estimated) - 0.001, 3), round(max(fitted_preferentiality_df$estimated) + 0.001, 3), length.out = 6)
  
  p_3 <- ggplot() +
    geom_tile(data = fitted_preferentiality_df, mapping = aes(x = x, y = y, fill = estimated)) + 
    geom_sf(data = USA, color = "black", fill = NA, lwd = 0.5) +
    scale_fill_gradientn(name = "Deg. Pref.", colors = pal, breaks = labs, labels = as.character(format(labs, nsmall = 3)), limits = c(labs[1], tail(labs, 1))) +
    labs(x = "Longitude", y = "Latitude", title = "") + 
    theme_bw() +
    theme(text = element_text(size = 24, family = "LM Roman 10"), 
          legend.key.height = unit(2.57, "cm"),
          legend.key.width  = unit(0.75, "cm"),
          plot.margin = margin(0, 0, 0, 0, "cm"))
  ggsave(filename = paste("APPLICATION/IMAGES/", sprintf("%02d", scenario), "/preferentiality.jpeg", sep = ""), plot = p_3, width = 3250, height = 2000, units = c("px"), dpi = 300, bg = "white")
}

