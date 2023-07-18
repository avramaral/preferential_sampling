args <- commandArgs(trailingOnly = TRUE)
scenario <- as.numeric(args[1]) # 1:11
div      <- as.numeric(args[2]) # 1:09

print(paste("Division: ", div, sep = ""))

library("plot3D")
library("tidyverse")
library("sf")
library("rgeoboundaries")
library("tidycensus")

data(fips_codes)

fips_codes <- fips_codes |> dplyr::select(state_code, state_name)
fips_codes <- fips_codes |> distinct(state_code, state_name, .keep_all = TRUE)
fips_codes <- as_tibble(fips_codes)
fips_codes <- fips_codes |> mutate(across(.cols = matches("state_code"), .fns = ~ as.double(.x)))
fips_codes <- fips_codes |> rename(`State Code` = state_code)

##################################################
# scenario <- 04
##################################################

d1 <- c("Connecticut", "Maine", "Massachusetts", "New Hampshire", "Rhode Island", "Vermont")
d2 <- c("New Jersey", "New York", "Pennsylvania")
d3 <- c("Illinois", "Indiana", "Michigan", "Ohio", "Wisconsin")
d4 <- c("Iowa", "Kansas", "Minnesota", "Missouri", "Nebraska", "North Dakota", "South Dakota")
d5 <- c("Delaware", "Florida", "Georgia", "Maryland", "North Carolina", "South Carolina", "Virginia", "District of Columbia", "West Virginia")
d6 <- c("Alabama", "Kentucky", "Mississippi", "Tennessee")
d7 <- c("Arkansas", "Louisiana", "Oklahoma", "Texas")
d8 <- c("Arizona", "Colorado", "Idaho", "Montana", "Nevada", "New Mexico", "Utah", "Wyoming")
d9 <- c("Alaska", "California", "Hawaii", "Oregon", "Washington")

d10 <- c("Montana", "North Dakota", "South Dakota", "Wyoming", "Nebraska")
d11 <- c("Montana", "North Dakota", "South Dakota", "Wyoming")

dd <- list(d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11)

##################################################

year <- 2022
data <- read_csv(file = paste("APPLICATION/usa_filtered_data_", year, ".csv", sep = ""))
data <- data |> filter(`Sample Duration` %in% c("24 HOUR", "24-HR BLK AVG"))
data <- data |> distinct(Latitude, Longitude, .keep_all = TRUE) # Choose the first POC, if applied
data <- data |> left_join(y = fips_codes, by = "State Code")
data <- data |> dplyr::select("state_name", "Latitude", "Longitude", "Arithmetic Mean", "Arithmetic Standard Dev")
data <- data |> rename(latt = Latitude, long = Longitude, mean = `Arithmetic Mean`, sd = `Arithmetic Standard Dev`)
data <- data |> filter(latt > 22.5, latt < 50, long > -130) # Remove stations located on Alaska and other small territories
data <- data |> filter(mean > 0) 
data <- st_as_sf(x = data, coords = c("long", "latt"), crs = "WGS84")
data <- st_transform(x = data, crs = "+init=epsg:6345 +units=km +no_defs")
nrow(data)

data_CV <- list()
data_test <- list()
for (d in 1:11) {
  data_CV[[d]]   <- data |> filter(!state_name %in% dd[[d]]) 
  data_test[[d]] <- data |> filter( state_name %in% dd[[d]]) 
}

saveRDS(data_test, file = "APPLICATION/FITTED_MODELS/CV/data_CV.rds")
data <- data_CV[[div]]

##################################################

if (FALSE) {
  USA <- geoboundaries(country = "USA", adm_lvl = 1)  
} else {
  USA <- readRDS(file = "APPLICATION/USA.rds")
}
USA <- USA[!(USA$shapeName %in% c("Alaska", "American Samoa","Commonwealth of the Northern Mariana Islands", "Guam", "Hawaii", "Puerto Rico", "United States Virgin Islands")), ]
USA <- st_transform(x = USA, crs = "+init=epsg:6345 +units=km +no_defs")

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
  center_pts <- rbind(c(-118, 35), c(-116, 41), c(-105, 45), c(-98, 34), c(-90, 45), c(-85, 37), c(-75, 42), c(-95, 38.5), c(-108, 38), c(-83, 32.5))
} else if (scenario ==  5) { 
  smoothing_kernel <- "Wendland"
  center_pts <- rbind(c(-118, 35), c(-116, 41), c(-105, 45), c(-98, 34), c(-90, 45), c(-85, 37), c(-75, 42), c(-95, 38.5), c(-112, 38), c(-83, 32.5), c(-105, 32.5), c(-109, 42.5), c(-96, 46.5), c(-82, 41), c(-96, 43)) 
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

# If prediction ##################################
coord_pred <- readRDS(file = "APPLICATION/coord_pred.rds")
n_coord_pr <- nrow(coord_pred)

data_coord <- sf::st_coordinates(data_test[[div]]$geometry)
colnames(data_coord) <- c("x", "y")
coord_pred <- rbind(coord_pred, data_coord)
##################################################

bfs_pred <- base_functions(center_pts = center_pts, loct = coord_pred, smoothing_kernel = smoothing_kernel, bandwidth = bandwidth, NS = NS, map = map)

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
result <- model_fitting(fs_formula = fs_formula, compute_scores = TRUE, verbose = FALSE)
##################################################
##################################################

result_obj <- list(result = result,
                   bfs_pred = bfs_pred,
                   fs_formula = fs_formula,
                   coord_pred = coord_pred,
                   n_coord_pr = n_coord_pr)

saveRDS(object = result_obj, file = paste("APPLICATION/FITTED_MODELS/CV/scenario_", scenario, "_division_", div, ".rds", sep = ""))
