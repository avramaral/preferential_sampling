args <- commandArgs(trailingOnly = TRUE)
scenario <- as.numeric(args[1]) 
div      <- as.numeric(args[2]) # 1:09

print("Application process CV")

library("plot3D")
library("tidyverse")

source("header.R")
source("functions_data_generation.R")
source("functions_model_fitting.R")
source("functions_result_processing.R")
source("functions_aux.R")

# scenario <- 3
# div      <- 8

data <- readRDS(file = "APPLICATION/FITTED_MODELS/CV/data_CV.rds")
data <- data[[div]]
result_obj <- readRDS(file = paste("APPLICATION/FITTED_MODELS/CV/scenario_", scenario, "_division_", div, ".rds", sep = ""))
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
fitted_latent <- as.data.frame(fitted_latent)
fitted_latent <- st_as_sf(x = fitted_latent, coords = c("x", "y"), crs = "+init=epsg:6345 +units=km +no_defs")

##########################################################
# Estimated process in the missing observation locations #
##########################################################
opt <- "Estimated"

if (opt == "Original") {
  plotting <- data
  plotting <- plotting |> dplyr::select(mean, geometry)
  colnames(plotting) <- c("value", "geometry")
} else if (opt == "Estimated") {
  plotting <- fitted_latent
  plotting <- plotting |> dplyr::select(z, geometry)
  colnames(plotting) <- c("value", "geometry")
}

pal <- jet.col(n = 100, alpha = 0.9)
labs <- seq(round(min(plotting$value) - 0.01, 2), round(max(plotting$value) + 0.01, 2), length.out = 6)

p <- ggplot() +
  geom_sf(data = USA, color = "black", fill = NA, lwd = 0.5) +
  geom_sf(data = plotting, mapping = aes(color = value), size = 5) + 
  scale_colour_gradientn(name = expression(paste(PM[2.5], " level", sep = "")), colors = pal, breaks = labs, labels = as.character(format(labs, nsmall = 3)), limits = c(labs[1], tail(labs, 1))) +
  labs(x = "Longitude", y = "Latitude", title = opt) + 
  theme_bw() +
  theme(text = element_text(size = 24, family = "LM Roman 10"), 
        legend.key.height = unit(2.57, "cm"),
        legend.key.width  = unit(0.75, "cm"),
        plot.margin = margin(0, 0, 0, 0, "cm"))
# p
##################################################
##################################################

ss <- (result$fit$summary.fitted.values[idx_lat, c("sd")] ^ 2) + (1 / result$fit$summary.hyperpar["Precision for the Gaussian observations", "mean"])
CRPS_latent  <- mean(CRPS(orig = data$mean, pred = fitted_latent$z, s = sqrt(ss)))
SCRPS_latent <- mean(SCRPS(orig = data$mean, pred = fitted_latent$z, s = sqrt(ss)))
error_latent <- MSE(orig = data$mean, pred = fitted_latent$z)

scores <- list(error = error_latent,
               CRPS  = CRPS_latent,
               SCRPS = SCRPS_latent)

saveRDS(object = scores, file = paste("APPLICATION/FITTED_MODELS/CV/SCORE/scenario_", scenario, "_division_", div, ".rds", sep = ""))

##################################################
##################################################

if (FALSE) {
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
  p_2_loc <- p_2 + geom_sf(data = data, color = rgb(red = 0, green = 0, blue = 0, alpha = 1), size = 3, shape = 4)
  p_2_loc
  
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
    p_3_loc <- p_3 + geom_sf(data = data, color = rgb(red = 0, green = 0, blue = 0, alpha = 1), size = 3, shape = 4)
    p_3_loc
  }
}









