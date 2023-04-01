scenario         <- 8
sample_size      <- 1000
smoothing_kernel <- "Wendland"
partition_const  <- "non_constant"
n_center_points  <- "manymany"

normalized <- F

source("header.R")
source("functions_data_generation.R")
source("functions_model_fitting.R")
source("functions_result_processing.R")
source("functions_aux.R")

fitted_models <- c("var_DP_PS", "var_DP_PS", "var_DP_PS") 
# fitted_models <- c("non_PS", "const_DP_PS", "var_DP_PS")
N_models <- 100

int_error <- list()
lat_error <- list()
pre_error <- list()
int_MAPE <- list()
lat_MAPE <- list()
pre_MAPE <- list()

DIC_latent <- list()
WAIC_latent <- list()

CRPS_latent <- list()
SCRPS_latent <- list()
for (fitted_model in fitted_models) {
  
  int_error[[fitted_model]] <- rep(x = 0, times = N_models)
  lat_error[[fitted_model]] <- rep(x = 0, times = N_models)
  pre_error[[fitted_model]] <- rep(x = 0, times = N_models)
  
  int_MAPE[[fitted_model]] <- rep(x = 0, times = N_models)
  lat_MAPE[[fitted_model]] <- rep(x = 0, times = N_models)
  pre_MAPE[[fitted_model]] <- rep(x = 0, times = N_models)
  
  DIC_latent[[fitted_model]] <- rep(x = 0, times = N_models)
  WAIC_latent[[fitted_model]] <- rep(x = 0, times = N_models)
  
  CRPS_latent[[fitted_model]] <- rep(x = 0, times = N_models)
  SCRPS_latent[[fitted_model]] <- rep(x = 0, times = N_models)
  
  if (normalized) {
    errors <- readRDS(file = paste("FITTED_MODELS/", sprintf("%02d", scenario), "/ERRORS/", fitted_model, "-", smoothing_kernel, "-", partition_const, "-", n_center_points, "-", sample_size, ".rds", sep = "")) 
  } else {
    errors <- readRDS(file = paste("FITTED_MODELS/", sprintf("%02d", scenario), "/ERRORS/UNNORMALIZED_", fitted_model, "-", smoothing_kernel, "-", partition_const, "-", n_center_points, "-", sample_size, ".rds", sep = ""))
  }
  
  for (n in 1:N_models) {
    int_error[[fitted_model]][n] <- errors[[n]]$error_intensity
    lat_error[[fitted_model]][n] <- errors[[n]]$error_latent
    pre_error[[fitted_model]][n] <- errors[[n]]$error_preferentiality
    
    int_MAPE[[fitted_model]][n] <- errors[[n]]$MAPE_intensity
    lat_MAPE[[fitted_model]][n] <- errors[[n]]$MAPE_latent
    pre_MAPE[[fitted_model]][n] <- errors[[n]]$MAPE_preferentiality
    
    DIC_latent[[fitted_model]][n] <- errors[[n]]$DIC_latent
    WAIC_latent[[fitted_model]][n] <- errors[[n]]$WAIC_latent
    
    CRPS_latent[[fitted_model]][n] <- errors[[n]]$CRPS_latent
    SCRPS_latent[[fitted_model]][n] <- errors[[n]]$SCRPS_latent
  }
}

errors_summary <- list(int_error = int_error, 
                       lat_error = lat_error, 
                       pre_error = pre_error,
                       int_MAPE = int_MAPE, 
                       lat_MAPE = lat_MAPE, 
                       pre_MAPE = pre_MAPE,
                       DIC_latent = DIC_latent, 
                       WAIC_latent = WAIC_latent,
                       CRPS_latent = CRPS_latent, 
                       SCRPS_latent = SCRPS_latent)

if (TRUE) {
  print("Intensity Process")
  print(paste("non_PS:      ", round(mean(errors_summary$int_error$non_PS),      4), " (", round(sd(errors_summary$int_error$non_PS),      4), ")", sep = ""))
  print(paste("const_DP_PS: ", round(mean(errors_summary$int_error$const_DP_PS), 4), " (", round(sd(errors_summary$int_error$const_DP_PS), 4), ")", sep = ""))
  print(paste("var_DP_PS:   ", round(mean(errors_summary$int_error$var_DP_PS),   4), " (", round(sd(errors_summary$int_error$var_DP_PS),   4), ")", sep = ""))
  cat("\n")
  
  print("Latent Process")
  print(paste("non_PS:      ", round(mean(errors_summary$lat_error$non_PS),      4), " (", round(sd(errors_summary$lat_error$non_PS),      4), ")", sep = ""))
  print(paste("const_DP_PS: ", round(mean(errors_summary$lat_error$const_DP_PS), 4), " (", round(sd(errors_summary$lat_error$const_DP_PS), 4), ")", sep = ""))
  print(paste("var_DP_PS:   ", round(mean(errors_summary$lat_error$var_DP_PS),   4), " (", round(sd(errors_summary$lat_error$var_DP_PS),   4), ")", sep = ""))
  cat("\n")

  print("Preferentiality Process")
  print(paste("non_PS:      ", round(mean(errors_summary$pre_error$non_PS),      4), " (", round(sd(errors_summary$pre_error$non_PS),      4), ")", sep = ""))
  print(paste("const_DP_PS: ", round(mean(errors_summary$pre_error$const_DP_PS), 4), " (", round(sd(errors_summary$pre_error$const_DP_PS), 4), ")", sep = ""))
  print(paste("var_DP_PS:   ", round(mean(errors_summary$pre_error$var_DP_PS),   4), " (", round(sd(errors_summary$pre_error$var_DP_PS),   4), ")", sep = ""))
  cat("\n")
  
  cat("---------------------------------------------")
  cat("\n")
  print("Intensity Process")
  print(paste("non_PS:      ", round(mean(errors_summary$int_MAPE$non_PS),      4), " (", round(sd(errors_summary$int_MAPE$non_PS),      4), ")", sep = ""))
  print(paste("const_DP_PS: ", round(mean(errors_summary$int_MAPE$const_DP_PS), 4), " (", round(sd(errors_summary$int_MAPE$const_DP_PS), 4), ")", sep = ""))
  print(paste("var_DP_PS:   ", round(mean(errors_summary$int_MAPE$var_DP_PS),   4), " (", round(sd(errors_summary$int_MAPE$var_DP_PS),   4), ")", sep = ""))
  cat("\n")
  
  print("Latent Process")
  print(paste("non_PS:      ", round(mean(errors_summary$lat_MAPE$non_PS),      4), " (", round(sd(errors_summary$lat_MAPE$non_PS),      4), ")", sep = ""))
  print(paste("const_DP_PS: ", round(mean(errors_summary$lat_MAPE$const_DP_PS), 4), " (", round(sd(errors_summary$lat_MAPE$const_DP_PS), 4), ")", sep = ""))
  print(paste("var_DP_PS:   ", round(mean(errors_summary$lat_MAPE$var_DP_PS),   4), " (", round(sd(errors_summary$lat_MAPE$var_DP_PS),   4), ")", sep = ""))
  cat("\n")
  
  print("Preferentiality Process")
  print(paste("non_PS:      ", round(mean(errors_summary$pre_MAPE$non_PS),      4), " (", round(sd(errors_summary$pre_MAPE$non_PS),      4), ")", sep = ""))
  print(paste("const_DP_PS: ", round(mean(errors_summary$pre_MAPE$const_DP_PS), 4), " (", round(sd(errors_summary$pre_MAPE$const_DP_PS), 4), ")", sep = ""))
  print(paste("var_DP_PS:   ", round(mean(errors_summary$pre_MAPE$var_DP_PS),   4), " (", round(sd(errors_summary$pre_MAPE$var_DP_PS),   4), ")", sep = ""))
  cat("\n")
  
  cat("---------------------------------------------")
  cat("\n")
  print("Latent Process (others)")
  print("DIC")
  print(paste("non_PS:      ", round(mean(errors_summary$DIC_latent$non_PS),      4), " (", round(sd(errors_summary$DIC_latent$non_PS),      4), ")", sep = ""))
  print(paste("const_DP_PS: ", round(mean(errors_summary$DIC_latent$const_DP_PS), 4), " (", round(sd(errors_summary$DIC_latent$const_DP_PS), 4), ")", sep = ""))
  print(paste("var_DP_PS:   ", round(mean(errors_summary$DIC_latent$var_DP_PS),   4), " (", round(sd(errors_summary$DIC_latent$var_DP_PS),   4), ")", sep = ""))
  cat("\n")
  print("WAIC")
  print(paste("non_PS:      ", round(mean(errors_summary$WAIC_latent$non_PS),      4), " (", round(sd(errors_summary$WAIC_latent$non_PS),      4), ")", sep = ""))
  print(paste("const_DP_PS: ", round(mean(errors_summary$WAIC_latent$const_DP_PS), 4), " (", round(sd(errors_summary$WAIC_latent$const_DP_PS), 4), ")", sep = ""))
  print(paste("var_DP_PS:   ", round(mean(errors_summary$WAIC_latent$var_DP_PS),   4), " (", round(sd(errors_summary$WAIC_latent$var_DP_PS),   4), ")", sep = ""))
  cat("\n")
  
  cat("---------------------------------------------")
  cat("\n")
  print("Latent Process (S)CRPS")
  print("CRPS")
  print(paste("non_PS:      ", round(mean(errors_summary$CRPS_latent$non_PS),      4), " (", round(sd(errors_summary$CRPS_latent$non_PS),      4), ")", sep = ""))
  print(paste("const_DP_PS: ", round(mean(errors_summary$CRPS_latent$const_DP_PS), 4), " (", round(sd(errors_summary$CRPS_latent$const_DP_PS), 4), ")", sep = ""))
  print(paste("var_DP_PS:   ", round(mean(errors_summary$CRPS_latent$var_DP_PS),   4), " (", round(sd(errors_summary$CRPS_latent$var_DP_PS),   4), ")", sep = ""))
  cat("\n")
  print("SCRPS")
  print(paste("non_PS:      ", round(mean(errors_summary$SCRPS_latent$non_PS),      4), " (", round(sd(errors_summary$SCRPS_latent$non_PS),      4), ")", sep = ""))
  print(paste("const_DP_PS: ", round(mean(errors_summary$SCRPS_latent$const_DP_PS), 4), " (", round(sd(errors_summary$SCRPS_latent$const_DP_PS), 4), ")", sep = ""))
  print(paste("var_DP_PS:   ", round(mean(errors_summary$SCRPS_latent$var_DP_PS),   4), " (", round(sd(errors_summary$SCRPS_latent$var_DP_PS),   4), ")", sep = ""))
  cat("\n")
}

