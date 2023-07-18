library("tidyverse")
library("sf")
library("rgeoboundaries")
library("gstat")
library("plot3D")
library("terra")
library("raster")
library("patchwork")
library("dismo")

#########################
# READ DATA
#########################

year <- 2022
data <- read_csv(file = paste("APPLICATION/usa_filtered_data_", year, ".csv", sep = ""))
data <- data |> filter(`Sample Duration` %in% c("24 HOUR", "24-HR BLK AVG"))
data <- data |> distinct(Latitude, Longitude, .keep_all = TRUE) # Choose the first POC, if applied
data <- data |> dplyr::select("Latitude", "Longitude", "Arithmetic Mean", "Arithmetic Standard Dev")
data <- data |> rename(latt = Latitude, long = Longitude, mean = `Arithmetic Mean`, sd = `Arithmetic Standard Dev`)
data <- data |> filter(latt > 22.5, latt < 50, long > -130) # Remove stations located on Alaska and other small territories
data <- data |> filter(mean > 0) 
# data <- data |> filter(mean < 15) # Values higher than 15 does not make much sense based on other works (10 out of 942: > 15)
data <- st_as_sf(x = data, coords = c("long", "latt"), crs = "WGS84")
data <- st_transform(x = data, crs = "+init=epsg:6345 +units=km +no_defs")
nrow(data)
# saveRDS(data, file = "APPLICATION/data.rds")

USA <- readRDS(file = "APPLICATION/USA.rds")
USA <- USA[!(USA$shapeName %in% c("Alaska", "American Samoa","Commonwealth of the Northern Mariana Islands", "Guam", "Hawaii", "Puerto Rico", "United States Virgin Islands")), ]
USA <- st_transform(x = USA, crs = "+init=epsg:6345 +units=km +no_defs")

pal <- jet.col(n = 100, alpha = 0.75)
labs <- seq(round(min(data$mean) - 0.01, 2), round(max(data$mean) + 0.01, 2), length.out = 6)

p1 <- ggplot() +
  geom_sf(data = USA, color = "black", fill = NA, lwd = 0.5) +
  geom_sf(data = data, aes(color = mean), size = 5) +
  scale_color_gradientn(name = expression(paste(PM[2.5], " level", sep = "")), colors = pal, breaks = labs, labels = as.character(format(labs, nsmall = 3)), limits = c(labs[1], tail(labs, 1))) +
  labs(x = "Longitude", y = "Latitude", title = "") + 
  theme_bw() +
  theme(text = element_text(size = 24, family = "LM Roman 10"), 
        legend.key.height = unit(2.57, "cm"),
        legend.key.width  = unit(0.75, "cm"),
        plot.margin = margin(0, 0, 0, 0, "cm"))

#########################
# FIT MODEL
#########################

coord_pred <- readRDS(file = "APPLICATION/coord_pred.rds")
coord_pred <- st_as_sf(coord_pred, coords = c("x", "y"), crs = st_crs(data))

grid <- terra::rast(USA, nrows = 100, ncols = 100)
xy <- terra::xyFromCell(grid, 1:ncell(grid))
coord_pred <- st_as_sf(as.data.frame(xy), coords = c("x", "y"), crs = st_crs(data))
coord_pred <- st_filter(coord_pred, USA)

res  <- gstat(formula = mean ~ 1, locations = data, nmax = nrow(data), set = list(idp = 2.1))
resp <- predict(res, coord_pred)

resp$x <- st_coordinates(resp)[, 1]
resp$y <- st_coordinates(resp)[, 1]
resp$pred <- resp$var1.pred
pred <- terra::rasterize(resp, grid, field = "pred", fun = "mean")
pred <- raster(pred)
pred <- as(pred, "SpatialPixelsDataFrame")
pred <- as.data.frame(pred)
colnames(pred) <- c("mean", "x", "y")

pal <- jet.col(n = 100, alpha = 0.75)
labs <- seq(round(min(pred$mean) - 0.01, 2), round(max(pred$mean) + 0.01, 2), length.out = 6)

p2 <- ggplot() +
  geom_sf(data = USA, color = "black", fill = NA, lwd = 0.5) +
  geom_tile(data = pred, aes(x = x, y = y, fill = mean)) +
  geom_sf(data = data, color = rgb(red = 0, green = 0, blue = 0, alpha = 1), size = 3, shape = 4, alpha = 1) + 
  scale_fill_gradientn(name = expression(paste(PM[2.5], " level", sep = "")), colors = pal, breaks = labs, labels = as.character(format(labs, nsmall = 3)), limits = c(labs[1], tail(labs, 1))) +
  labs(x = "Longitude", y = "Latitude", title = "") + 
  theme_bw() +
  theme(text = element_text(size = 24, family = "LM Roman 10"), 
        legend.key.height = unit(2.57, "cm"),
        legend.key.width  = unit(0.75, "cm"),
        plot.margin = margin(0, 0, 0, 0, "cm"))

ggsave(filename = paste("APPLICATION/IMAGES/free_assumptions.jpeg", sep = ""), plot = p2, width = 3250, height = 2000, units = c("px"), dpi = 300, bg = "white")
# p_total <- p1 + p2 + plot_annotation(theme = theme(plot.margin = margin()))

#########################
# CROSS-VALIDATION
#########################

set.seed(1)

RMSE <- function (observed, predicted) { sqrt(mean((observed - predicted) ^ 2)) }

idps <- seq(0, 5, 0.1)

k <- 5
KF <- dismo::kfold(nrow(data), k = k)

rmse_list <- list()
for (i in 1:k) {
  rmse_list[[paste("k", i, sep = "")]] <- list()
  
  test  <- data[KF == i, ]
  train <- data[KF != i, ]
  
  for (j in 1:length(idps)) {
    gs <- gstat(formula = mean ~ 1, locations = train, nmax = nrow(train), set = list(idp = idps[j]))
    pred_tmp <- predict(gs, test)$var1.pred
    rmse_list[[paste("k", i, sep = "")]][[paste("idp", j, sep = "")]] <- RMSE(test$mean, pred_tmp)
  }
}

rmse_summ <- rep(0, length(idps))
for (i in 1:k) {
  for (j in 1:length(idps)) {
    rmse_summ[j] <- rmse_summ[j] + rmse_list[[paste("k", i, sep = "")]][[paste("idp", j, sep = "")]]
  }
}
rmse_summ <- rmse_summ / k
idps[which.min(rmse_summ)]


