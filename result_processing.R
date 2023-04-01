args <- commandArgs(trailingOnly = TRUE)
scenario         <- as.numeric(args[1]) 
sample_size      <- as.numeric(args[2]) 
smoothing_kernel <- args[3] 
partition_const  <- args[4] 
n_center_points  <- args[5] ################################### UPDATE THIS
fitted_model     <- args[6] 
normalize        <- as.logical(as.numeric(args[7]))

source("header.R")
source("functions_data_generation.R")
source("functions_model_fitting.R")
source("functions_result_processing.R")
source("functions_aux.R")

normalize_function <- function (r1, r2, ...) {
  mm <- max(c(values(r1), values(r2)))
  mn <- min(c(values(r1), values(r2)))
  values(r1) <- (values(r1) - mn) / (mm - mn)
  values(r2) <- (values(r2) - mn) / (mm - mn)
  
  list(r1, r2)
}

plotting_true_fitted <- function (true, fitted, title = "Latent", fitted_model = "non_PS", normalize = normalize, ...) {
  par(mfrow = c(1, 2)); 
  len  <- 199
  mn   <- ifelse(normalize, 0, min(c(values(true), values(fitted))) - 0.01) 
  mm   <- ifelse(normalize, 1, max(c(values(true), values(fitted))) + 0.01) 
  brks <- round(seq(from = mn, to = mm, length.out = len), 9)
  arg  <- list(at = seq(from = mn, to = mm, length.out = 5), labels = round(seq(from = mn, to = mm, length.out = 5), digits = 2)) 
  plot(true,   breaks = brks, col = rainbow(n = len, start = 0.05, end = 0.85), axis.arg = arg, main = paste("Original ", title, " Process", sep = ""))
  plot(fitted, breaks = brks, col = rainbow(n = len, start = 0.05, end = 0.85), axis.arg = arg, main = paste("Estimated ", title, " Process (", fitted_model, ")", sep = ""))
  par(mfrow = c(1, 1))
}

set.seed(1)

##################################################
# Initial parameters
##################################################
xlim <- c(0, 10)
ylim <- c(0, 10)
by <- 0.05
mu <- 2
nu <- 1
scl <- 2
sig2 <- 0.5
sig2_error <- 0.1

n_points <- NULL
pref_sampling <- TRUE

NS <- FALSE
normalize <- normalize
plotting <- TRUE
##################################################

N_sim <- 110

data <- readRDS(file = paste("DATA/", sprintf("%02d", scenario), "/", sample_size, ".rds", sep = ""))
results <- readRDS(file = paste("FITTED_MODELS/", sprintf("%02d", scenario), "/", fitted_model, "-", smoothing_kernel, "-", partition_const, "-", n_center_points, "-", sample_size, ".rds", sep = ""))
others <- readRDS(file = paste("FITTED_MODELS/", sprintf("%02d", scenario), "/OTHERS/", fitted_model, "-", smoothing_kernel, "-", partition_const, "-", n_center_points, "-", sample_size, ".rds", sep = ""))

##################################################

N_models <- 100
failed_models <- rep(x = 0, times = N_sim)
for (n in 1:N_sim) { failed_models[n] <- ifelse(test = !results[[n]]$fit$ok, yes = 1, no = 0) }
idxs <- which(failed_models == 0)
if (length(idxs) >= N_models) { idxs <- idxs[1:N_models] } else { stop(paste("Less than ", N_models , " successfully fitted models.", sep = "")) }

data <- data[idxs]
results <- results[idxs]
others <- others[idxs]

##################################################

errors <- list()
for (n in 1:N_models) {
  print(paste("Simulation: ", sprintf("%03d", n), sep = ""))
  
  latt <- data[[n]]$latt
  pref <- data[[n]]$pref
  loct <- data[[n]]$loct
  Yobs <- data[[n]]$Yobs
  orig <- data[[n]]$orig
  
  alpha <- data[[n]]$alpha
  
  result <- results[[n]]
  fs_formula <- others[[n]]$fs_formula
  
  ##############################
  # Intensity Process
  ##############################
  
  if (fitted_model %in% c("const_DP_PS", "var_DP_PS")) { 
    
    idx_int <- inla.stack.index(fs_formula$full_stack, tag = "pred_pp")$data
    fitted_intensity_partial <- exp(result$fit$summary.fitted.values[idx_int, c("mean")])
    
    fitted_intensity <- loct$lambda
    values(fitted_intensity) <- fitted_intensity_partial
    
    MAPE_intensity <- MAPE(orig = loct$lambda, pred = fitted_intensity)
    print(paste("MAPE (intensity): ", round(x = MAPE_intensity, digits = 4), sep = ""))
    
    if (normalize) { r <- normalize_function(r1 = loct$lambda, r2 = fitted_intensity); true_intensity <- r[[1]]; fitted_intensity <- r[[2]] } else { true_intensity <- loct$lambda }
    
    error_intensity <- MSE(orig = true_intensity, pred = fitted_intensity)
    print(paste("Prediction error (intensity): ", round(x = error_intensity, digits = 4), sep = ""))
    
    
    if (plotting) { plotting_true_fitted(true = true_intensity, fitted = fitted_intensity, title = "Intensity", fitted_model = fitted_model, normalize = normalize) }
  } else { error_intensity <- NA ; MAPE_intensity <- NA }
  
  ##############################
  # Latent Process
  ##############################
  
  idx_lat <- inla.stack.index(fs_formula$full_stack, tag = "pred_y")$data
  fitted_latent_partial <- result$fit$summary.fitted.values[idx_lat, c("mean")]
  
  fitted_latent <- latt 
  values(fitted_latent) <- fitted_latent_partial
  
  ss <- (result$fit$summary.fitted.values[idx_lat, c("sd")] ^ 2) + (1 / result$fit$summary.hyperpar["Precision for the Gaussian observations", "mean"])
  CRPS_latent <- mean(CRPS(orig = orig, pred = fitted_latent, s = sqrt(ss)))
  SCRPS_latent <- mean(SCRPS(orig = orig, pred = fitted_latent, s = sqrt(ss)))
  print(paste("CRPS (latent):  ", round(x = CRPS_latent,  digits = 4), sep = ""))
  print(paste("SCRPS (latent): ", round(x = SCRPS_latent, digits = 4), sep = ""))
  
  MAPE_latent <- MAPE(orig = orig, pred = fitted_latent)
  print(paste("MAPE (latent): ", round(x = MAPE_latent, digits = 4), sep = ""))
  
  if (normalize) { r <- normalize_function(r1 = orig, r2 = fitted_latent); true_latent <- r[[1]]; fitted_latent <- r[[2]] } else { true_latent <- orig }
  
  error_latent <- MSE(orig = true_latent, pred = fitted_latent)
  print(paste("Prediction error (latent): ", round(x = error_latent, digits = 4), sep = ""))

  if (plotting) { plotting_true_fitted(true = true_latent, fitted = fitted_latent, title = "Latent", fitted_model = fitted_model, normalize = normalize) }

  ##############################
  # Preferentiality Process
  ##############################
  
  coord_pred <- others[[n]]$coord_pred
  bfs_pred <- others[[n]]$bfs_pred
  
  if (fitted_model == "var_DP_PS") {
    
    fitted_preferentiality <- compute_preferentiatility(fit = result$fit, coord_pred = coord_pred, bfs_pred = bfs_pred, xlim = xlim, by = by)
    
    MAPE_preferentiality <- MAPE(orig = pref, pred = fitted_preferentiality)
    print(paste("MAPE (latent): ", round(x = MAPE_preferentiality, digits = 4), sep = ""))
    
    if (normalize) { r <- normalize_function(r1 = pref, r2 = fitted_preferentiality); true_preferentiality <- r[[1]]; fitted_preferentiality <- r[[2]] } else { true_preferentiality <- pref }
    
    error_preferentiality <- MSE(orig = true_preferentiality, pred = fitted_preferentiality)
    print(paste("Prediction error (preferentiatility): ", round(x = error_preferentiality, digits = 4), sep = ""))
    
    if (plotting) { plotting_true_fitted(true = true_preferentiality, fitted = fitted_preferentiality, title = "Preferentiality", fitted_model = fitted_model, normalize = normalize) }
    
  } else { error_preferentiality <- NA ; MAPE_preferentiality <- NA }
  
  ##############################
  ##############################
  
  errors[[n]] <- list(error_intensity = error_intensity, 
                      error_latent = error_latent, 
                      error_preferentiality = error_preferentiality, 
                      MAPE_intensity = MAPE_intensity,
                      MAPE_latent = MAPE_latent,
                      MAPE_preferentiality = MAPE_preferentiality,
                      DIC_latent = result$fit$dic$dic,
                      WAIC_latent = result$fit$waic$waic,
                      CRPS_latent = CRPS_latent, 
                      SCRPS_latent = SCRPS_latent)
}

if (normalize) {
  saveRDS(errors, file = paste("FITTED_MODELS/", sprintf("%02d", scenario), "/ERRORS/", fitted_model, "-", smoothing_kernel, "-", partition_const, "-", n_center_points, "-", sample_size, ".rds", sep = ""))  
} else {
  saveRDS(errors, file = paste("FITTED_MODELS/", sprintf("%02d", scenario), "/ERRORS/UNNORMALIZED_", fitted_model, "-", smoothing_kernel, "-", partition_const, "-", n_center_points, "-", sample_size, ".rds", sep = ""))  
}

