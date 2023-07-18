args <- commandArgs(trailingOnly = TRUE)
scenario <- as.numeric(args[1]) # 1:11

library("plot3D")
library("tidyverse")
library("sf")
library("rgeoboundaries")

##################################################
# 01: Non-PS
# 02: Constant-PS
# 03: Varying-PS with Wendland (05 functions)
# 04: Varying-PS with Wendland (10 functions)
# 05: Varying-PS with Wendland (15 functions)
# 06: Varying-PS with horiz-triang (05 functions)
# 07: Varying-PS with horiz-triang (10 functions)
# 08: Varying-PS with verti-triang (05 functions)
# 09: Varying-PS with verti-triang (10 functions)
# 10: Varying-PS with partition per region (04 functions)
# 11: Varying-PS with partition per division (09 functions)
##################################################
# scenario <- 04
##################################################

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
saveRDS(data, file = "APPLICATION/data.rds")

if (FALSE) {
  USA <- geoboundaries(country = "USA", adm_lvl = 1)  
} else {
  USA <- readRDS(file = "APPLICATION/USA.rds")
}
USA <- USA[!(USA$shapeName %in% c("Alaska", "American Samoa","Commonwealth of the Northern Mariana Islands", "Guam", "Hawaii", "Puerto Rico", "United States Virgin Islands")), ]
USA <- st_transform(x = USA, crs = "+init=epsg:6345 +units=km +no_defs")

pal <- jet.col(n = 100, alpha = 0.75)
labs <- seq(round(min(data$mean) - 0.01, 2), round(max(data$mean) + 0.01, 2), length.out = 6)

p <- ggplot() +
  geom_sf(data = USA, color = "black", fill = NA, lwd = 0.5) +
  geom_sf(data = data, mapping = aes(color = mean), size = 7, shape = 19) + #geom_sf(data = data, aes(color = mean), size = 5) + 
  # scale_colour_gradientn(name = expression(paste(PM[2.5], " level", sep = "")), colors = pal, breaks = labs, labels = as.character(format(labs, nsmall = 3)), limits = c(labs[1], tail(labs, 1))) +
  scale_colour_gradientn(name = expression(paste(PM[2.5], " level", sep = "")), colors = pal, breaks = labs, labels = as.character(format(labs, nsmall = 3)), limits = c(labs[1], tail(labs, 1))) +
  labs(x = "Longitude", y = "Latitude", title = "") + 
  theme_bw() +
  theme(text = element_text(size = 24, family = "LM Roman 10"), 
        legend.key.height = unit(2.57, "cm"),
        legend.key.width  = unit(0.75, "cm"),
        plot.margin = margin(0, 0, 0, 0, "cm"))
#ggsave(filename = "APPLICATION/usa_pm25.jpeg", plot = p, width = 2850, height = 2000, units = c("px"), dpi = 300, bg = "white")
ggsave(filename = "APPLICATION/usa_pm25.jpeg", plot = p, width = 3250, height = 2000, units = c("px"), dpi = 300, bg = "white")

# https://en.wikipedia.org/wiki/List_of_regions_of_the_United_States
northwest <- st_union(USA[USA$shapeName %in% c("Connecticut", "Maine", "Massachusetts", "New Hampshire", "Rhode Island", "Vermont", "New Jersey", "New York", "Pennsylvania"), ], by_feature = FALSE) %>% st_as_sf() %>% st_cast("POLYGON") %>% mutate(area = st_area(.)) %>% arrange(desc(area)) %>% slice(1)
midwest   <- st_union(USA[USA$shapeName %in% c("Illinois", "Indiana", "Michigan", "Ohio", "Wisconsin", "Iowa", "Kansas", "Minnesota", "Missouri", "Nebraska", "North Dakota", "South Dakota"), ], by_feature = FALSE) %>% st_as_sf() %>% st_cast("POLYGON") %>% mutate(area = st_area(.)) %>% arrange(desc(area)) %>% slice(1)
south     <- st_union(USA[USA$shapeName %in% c("Delaware", "Florida", "Georgia", "Maryland", "North Carolina", "South Carolina", "Virginia", "District of Columbia", "West Virginia", "Alabama", "Kentucky", "Mississippi", "Tennessee", "Arkansas", "Louisiana", "Oklahoma", "Texas"), ], by_feature = FALSE) %>% st_as_sf() %>% st_cast("POLYGON") %>% mutate(area = st_area(.)) %>% arrange(desc(area)) %>% slice(1)
west      <- st_union(USA[USA$shapeName %in% c("Arizona", "Colorado", "Idaho", "Montana", "Nevada", "New Mexico", "Utah", "Wyoming", "Alaska", "California", "Hawaii", "Oregon", "Washington"), ], by_feature = FALSE) %>% st_as_sf() %>% st_cast("POLYGON") %>% mutate(area = st_area(.)) %>% arrange(desc(area)) %>% slice(1)

division_1 <- st_union(USA[USA$shapeName %in% c("Connecticut", "Maine", "Massachusetts", "New Hampshire", "Rhode Island", "Vermont"), ], by_feature = FALSE) %>% st_as_sf() %>% st_cast("POLYGON") %>% mutate(area = st_area(.)) %>% arrange(desc(area)) %>% slice(1)
division_2 <- st_union(USA[USA$shapeName %in% c("New Jersey", "New York", "Pennsylvania"), ], by_feature = FALSE) %>% st_as_sf() %>% st_cast("POLYGON") %>% mutate(area = st_area(.)) %>% arrange(desc(area)) %>% slice(1)
division_3 <- st_union(USA[USA$shapeName %in% c("Illinois", "Indiana", "Michigan", "Ohio", "Wisconsin"), ], by_feature = FALSE) %>% st_as_sf() %>% st_cast("POLYGON") %>% mutate(area = st_area(.)) %>% arrange(desc(area)) %>% slice(1)
division_4 <- st_union(USA[USA$shapeName %in% c("Iowa", "Kansas", "Minnesota", "Missouri", "Nebraska", "North Dakota", "South Dakota"), ], by_feature = FALSE) %>% st_as_sf() %>% st_cast("POLYGON") %>% mutate(area = st_area(.)) %>% arrange(desc(area)) %>% slice(1)
division_5 <- st_union(USA[USA$shapeName %in% c("Delaware", "Florida", "Georgia", "Maryland", "North Carolina", "South Carolina", "Virginia", "District of Columbia", "West Virginia"), ], by_feature = FALSE) %>% st_as_sf() %>% st_cast("POLYGON") %>% mutate(area = st_area(.)) %>% arrange(desc(area)) %>% slice(1)
division_6 <- st_union(USA[USA$shapeName %in% c("Alabama", "Kentucky", "Mississippi", "Tennessee"), ], by_feature = FALSE) %>% st_as_sf() %>% st_cast("POLYGON") %>% mutate(area = st_area(.)) %>% arrange(desc(area)) %>% slice(1)
division_7 <- st_union(USA[USA$shapeName %in% c("Arkansas", "Louisiana", "Oklahoma", "Texas"), ], by_feature = FALSE) %>% st_as_sf() %>% st_cast("POLYGON") %>% mutate(area = st_area(.)) %>% arrange(desc(area)) %>% slice(1)
division_8 <- st_union(USA[USA$shapeName %in% c("Arizona", "Colorado", "Idaho", "Montana", "Nevada", "New Mexico", "Utah", "Wyoming"), ], by_feature = FALSE) %>% st_as_sf() %>% st_cast("POLYGON") %>% mutate(area = st_area(.)) %>% arrange(desc(area)) %>% slice(1)
division_9 <- st_union(USA[USA$shapeName %in% c("Alaska", "California", "Hawaii", "Oregon", "Washington"), ], by_feature = FALSE) %>% st_as_sf() %>% st_cast("POLYGON") %>% mutate(area = st_area(.)) %>% arrange(desc(area)) %>% slice(1)

USA_regions   <- st_cast(dplyr::bind_rows(list(northwest, midwest, south, west))$x, "MULTIPOLYGON")
USA_divisions <- st_cast(dplyr::bind_rows(list(division_1, division_2, division_3, division_4, division_5, division_6, division_7, division_8, division_9))$x, "MULTIPOLYGON")

USA_regions2   <- USA_regions   |> sf::st_buffer(dist = 100) # Buffered by 100 km
USA_divisions2 <- USA_divisions |> sf::st_buffer(dist = 100) # Buffered by 100 km

##################################################
##################################################
source("header.R")
source("functions_model_fitting.R")
source("functions_aux.R")
source("functions_result_processing.R")
##################################################
##################################################

loct <- matrix(data = NA, nrow = length(data$geometry), ncol = 2)
colnames(loct) <- c("x", "y")
for (p in 1:length(data$geometry)) { loct[p, ] <- sf::st_coordinates(data$geometry[[p]]) }

loct  <- list(loct = loct)
if (FALSE) {
  USA0 <- geoboundaries(country = "USA", adm_lvl = 0) %>% st_as_sf() %>% st_cast("POLYGON") %>% mutate(area = st_area(.)) %>% arrange(desc(area)) %>% slice(1) # Main land  
} else {
  USA0 <- readRDS(file = "APPLICATION/USA0.rds")
}
USA0  <- st_transform(x = USA0, crs = "+init=epsg:6345 +units=km +no_defs")
USA02 <- USA0$geometry |> sf::st_buffer(dist = 100) # Buffered by 100 km

if (FALSE) {
  mesh_wgt <- mesh_weights(latt = USA0, max_edge_mf = c(300, 3000), offset_mf = c(300, 1500), PLOT = TRUE)  
} else {
  mesh_wgt <- readRDS(file = "APPLICATION/mesh_wgt.rds")
}

bandwidth <- 2500
map       <-  USA02 
if (scenario ==  1) {
  smoothing_kernel <- "Wendland" # Irrelevant
  center_pts <- rbind(c(-90, 40)) 
} else if (scenario ==  2) { 
  smoothing_kernel <- "Wendland" # Irrelevant
  center_pts <- rbind(c(-90, 40)) 
} else if (scenario ==  3) { 
  smoothing_kernel <- "Wendland"
  center_pts <- rbind(c(-118, 43.5), c(-95, 40), c(-80, 38), c(-105, 36), c(-90, 45)) 
} else if (scenario ==  4) { 
  smoothing_kernel <- "Wendland"
  center_pts <- rbind(c(-118, 35), c(-118, 44), c(-105, 45), c(-98, 34), c(-90, 45), c(-85, 37), c(-75, 42), c(-95, 38.5), c(-108, 38), c(-83, 32.5))
} else if (scenario ==  5) { 
  smoothing_kernel <- "Wendland"
  center_pts <- rbind(c(-118, 35), c(-118, 44), c(-105, 45), c(-98, 34), c(-90, 45), c(-85, 37), c(-75, 42), c(-95, 38.5), c(-112, 38), c(-83, 32.5), c(-105, 32.5), c(-109, 42.5), c(-96, 46.5), c(-82, 41), c(-96, 43)) 
} else if (scenario ==  6) {
  smoothing_kernel <- "unidirectional_triangular"
  uni_triang_hor <- range(USA0$geometry[[1]][[1]][, 2])
  uni_triang_hor <- seq(uni_triang_hor[1], uni_triang_hor[2], length.out = 5)
  center_pts <- matrix(data = c(rep(0, 5), uni_triang_hor), nrow = 5, ncol = 2) 
} else if (scenario ==  7) { 
  smoothing_kernel <- "unidirectional_triangular"
  uni_triang_hor <- range(USA0$geometry[[1]][[1]][, 2])
  uni_triang_hor <- seq(uni_triang_hor[1], uni_triang_hor[2], length.out = 10)
  center_pts <- matrix(data = c(rep(0, 10), uni_triang_hor), nrow = 10, ncol = 2)
} else if (scenario ==  8) { 
  smoothing_kernel <- "unidirectional_triangular"
  uni_triang_ver <- range(USA0$geometry[[1]][[1]][, 1])
  uni_triang_ver <- seq(uni_triang_ver[1], uni_triang_ver[2], length.out = 5)
  center_pts <- matrix(data = c(uni_triang_ver, rep(4500, 5)), nrow = 5, ncol = 2)
} else if (scenario ==  9) { 
  smoothing_kernel <- "unidirectional_triangular"
  uni_triang_ver <- range(USA0$geometry[[1]][[1]][, 1])
  uni_triang_ver <- seq(uni_triang_ver[1], uni_triang_ver[2], length.out = 10)
  center_pts <- matrix(data = c(uni_triang_ver, rep(4500, 10)), nrow = 10, ncol = 2)
} else if (scenario == 10) { 
  smoothing_kernel <- "partition_map"
  center_pts <- rbind(c(-90, 40)) # Irrelevant
  map <- USA_regions 
} else if (scenario == 11) { 
  smoothing_kernel <- "partition_map" 
  center_pts <- rbind(c(-90, 40)) # Irrelevant
  map <- USA_divisions
} else if (scenario == 12) {
  smoothing_kernel <- "partition_map" 
  center_pts <- rbind(c(-90, 40)) # Irrelevant
  map <- USA02
}

if (scenario %in% 3:5) {
  colnames(center_pts) <- c("long", "latt")
  center_pts <- st_as_sf(x = as_tibble(as.data.frame(center_pts)), coords = c("long", "latt"), crs = "WGS84")
  center_pts <- st_transform(x = center_pts, crs = "+init=epsg:6345 +units=km +no_defs")
  center_pts <- st_coordinates(center_pts)
  colnames(center_pts) <- NULL
}

# Locate the center points
if (FALSE) {
  plot(USA0$geometry[[1]][[1]], xlab = "Long", ylab = "Latt", main = "", cex = 0.5, pch = 19)
  points(center_pts, cex = 2, pch = 19, col = 2)
  points(loct$loct, cex = 0.5, pch = 19, col = 4)  
}

NS <- ifelse(scenario %in% 8:9, FALSE, TRUE)
bfs <- base_functions(center_pts = center_pts, mesh = mesh_wgt$mesh, loct = loct$loct, smoothing_kernel = smoothing_kernel, bandwidth = bandwidth, NS = NS, map = map)

# If prediction
if (FALSE) {
  resolution <- 25 
  pts_bdy <- USA0$geometry[[1]][[1]]
  pts_bdy_x <- range(pts_bdy[, 1])
  pts_bdy_y <- range(pts_bdy[, 2])
  coord_pred <- expand.grid(x = seq(pts_bdy_x[1], pts_bdy_x[2], by = resolution), y = seq(pts_bdy_y[1], pts_bdy_y[2], by = resolution))
  coordinates(coord_pred) <- ~ x + y
  
  xx <- as(st_as_sf(coord_pred), "sf");    st_crs(xx) <- st_crs(USA0$geometry)
  yy <- as(st_as_sf(USA0$geometry), "sf"); st_crs(yy) <- st_crs(USA0$geometry)
  
  pppts <- st_intersection(x = xx, y = yy)
  
  coord_pred <- matrix(data = NA, nrow = length(pppts$geometry), ncol = 2)
  colnames(coord_pred) <- c("x", "y")
  for (p in 1:length(pppts$geometry)) { coord_pred[p, ] <- sf::st_coordinates(pppts$geometry[[p]]) }
  coord_pred <- data.frame(x = coord_pred[, 1], y = coord_pred[, 2])
} else {
  coord_pred <- readRDS(file = "APPLICATION/coord_pred.rds")
}

n_coord_pr <- nrow(coord_pred)
data_coord <- sf::st_coordinates(data$geometry)
colnames(data_coord) <- c("x", "y")
coord_pred <- rbind(coord_pred, data_coord)

bfs_pred <- base_functions(center_pts = center_pts, loct = coord_pred, smoothing_kernel = smoothing_kernel, bandwidth = bandwidth, NS = NS, map = map)
if (TRUE) { plot_bfs(bfs = bfs_pred, coord = coord_pred, n_coord_pr = n_coord_pr) }

# Projection matrices (PM)
if (scenario == 1) {
  PS_pm <- FALSE
  bfs_pm <- NULL
  bfs_pred_pm <- NULL 
  coord_pred_pm <- coord_pred 
} else if (scenario == 2) {
  PS_pm <- TRUE
  bfs_pm <- NULL
  bfs_pred_pm <- NULL 
  coord_pred_pm <- coord_pred 
} else {
  PS_pm <- TRUE
  bfs_pm <- bfs 
  bfs_pred_pm <- bfs_pred 
  coord_pred_pm <- coord_pred 
}

loct_pm <- list(x = loct$loct[, 1], y = loct$loct[, 2], n = nrow(loct$loct)) 

pm <- projection_matrices(mesh = mesh_wgt$mesh, wgt = mesh_wgt$wgt, loct = loct_pm, PS = PS_pm, bfs = bfs_pm, bfs_pred = bfs_pred_pm, coord_pred = coord_pred_pm)

SPDE <- inla.spde2.pcmatern(mesh = mesh_wgt$mesh, alpha = 2, 
                            prior.range = c(1e3, 0.90), # P(range < 1e3) = 0.90
                            prior.sigma = c(1.0, 0.01)) # P(sigma > 1.0) = 0.01

gaus_prior <- list(prior = "gaussian", param = c(0, 5)) # Prior for the degree of preferentiality

fs_formula <- full_stack_formula(Yobs = data$mean, pm = pm, SPDE = SPDE, gaus_prior = gaus_prior)

##################################################
##################################################
result <- model_fitting(fs_formula = fs_formula, compute_scores = TRUE, verbose = TRUE)
##################################################
##################################################

result_obj <- list(result = result,
                   bfs_pred = bfs_pred,
                   fs_formula = fs_formula,
                   coord_pred = coord_pred,
                   n_coord_pr = n_coord_pr)
saveRDS(object = result_obj, file = paste("APPLICATION/FITTED_MODELS/scenario", scenario, ".rds", sep = ""))

############
############
# Analysis #
############
############

#####################
# Intensity Process #
#####################

idx_int <- inla.stack.index(fs_formula$full_stack, tag = "pred_pp")$data[1:n_coord_pr]
fitted_intensity <- cbind(coord_pred[1:n_coord_pr, ], exp(result$fit$summary.fitted.values[idx_int, c("mean")])[1:n_coord_pr])
colnames(fitted_intensity) <- c("x", "y", "z")

# Create a gridded spatial object from "process"
coordinates(fitted_intensity) <- ~ x + y
gridded(fitted_intensity) <- TRUE
fitted_intensity <- raster(fitted_intensity)

plot(fitted_intensity, main = "Intensity Process")

##################
# Latent Process #
##################

idx_lat <- inla.stack.index(fs_formula$full_stack, tag = "pred_y")$data[1:n_coord_pr]
fitted_latent <- cbind(coord_pred[1:n_coord_pr, ], result$fit$summary.fitted.values[idx_lat, c("mean")][1:n_coord_pr])
colnames(fitted_latent) <- c("x", "y", "z")

# Create a gridded spatial object from "process"
coordinates(fitted_latent) <- ~ x + y
gridded(fitted_latent) <- TRUE
fitted_latent <- raster(fitted_latent)

plot(fitted_latent, main = "Latent Process")

###########################
# Preferentiality Process #
###########################

fitted_preferentiality <- compute_preferentiatility(fit = result$fit, coord_pred = coord_pred, bfs_pred = bfs_pred, xlim = xlim, by = by, n_coord_pr = n_coord_pr)
  
plot(fitted_preferentiality, main = "Preferentiality Process")


