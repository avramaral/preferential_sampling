args <- commandArgs(trailingOnly = TRUE)
scenario         <- as.numeric(args[1]) # 1, 2, 3, 4, and 5
sample_size      <- as.numeric(args[2]) # 500, 1000, 2000, 4000
smoothing_kernel <- args[3] # partition", "unidirectional_triangular", "Wendland", "radial_basis"
partition_const  <- args[4] # "non_constant" and "constant" (for "partition")
n_center_points  <- args[5] # "few" and "many"
fitted_model     <- args[6] # "non_PS", "const_DP_PS", and "var_DP_PS"

print(paste("FITTED_MODELS/", sprintf("%02d", scenario), "/", fitted_model, "-", smoothing_kernel, "-", partition_const, "-", n_center_points, "-", sample_size, ".rds", sep = ""))

source("header.R")
source("functions_data_generation.R")
source("functions_model_fitting.R")
source("functions_result_processing.R")
source("functions_aux.R")

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
##################################################

N_sim <- 200
N_fit <- 110 # 110
data <- readRDS(file = paste("DATA/", sprintf("%02d", scenario), "/", sample_size, ".rds", sep = ""))

##################################################

results <- list()
others  <- list()
valid_n <- 0
n       <- 1
while ((n <= N_sim) & (valid_n < N_fit)) {
  print(paste("Count: ", sprintf("%03d", n), sep = ""))
  
  latt <- data[[n]]$latt
  pref <- data[[n]]$pref
  loct <- data[[n]]$loct
  Yobs <- data[[n]]$Yobs
  orig <- data[[n]]$orig

  alpha <- data[[n]]$alpha

  ##################################################
  # INLA fitting
  ##################################################

  mesh_wgt <- mesh_weights(latt = latt, max_edge_mf = c(0.1, 0.25), offset_mf = c(0.1, 0.15), cutoff_mf = 0.02, PLOT = FALSE)

  center_pts <- define_center_pts(smoothing_kernel = smoothing_kernel, partition_const = partition_const, n_center_points = n_center_points, xlim = xlim, NS = NS)

  bandwidth <- 10
  bfs <- base_functions(center_pts = center_pts, mesh = mesh_wgt$mesh, loct = loct$loct, smoothing_kernel = smoothing_kernel, bandwidth = bandwidth, NS = FALSE)

  # If prediction
  coord_pred <- as.data.frame(rasterToPoints(latt)[, c("x", "y")])
  bfs_pred <- base_functions(center_pts = center_pts, loct = coord_pred, smoothing_kernel = smoothing_kernel, bandwidth = bandwidth, NS = FALSE)
  if (FALSE) { plot_bfs(bfs = bfs_pred, coord = coord_pred) }

  # Projection matrices
  # The following parameters define the model that we will fit (out of the 6 options: non_PS, const_DP_PS, and var_DP_PS, all with or without the prediction step)

  if (fitted_model == "non_PS") {
    PS_pm  <- FALSE
    bfs_pm <- NULL
    bfs_pred_pm <- NULL
  } else if (fitted_model == "const_DP_PS") {
    PS_pm  <- TRUE
    bfs_pm <- NULL
    bfs_pred_pm <- NULL
  } else if (fitted_model == "var_DP_PS") {
    PS_pm  <- TRUE
    bfs_pm <- bfs
    bfs_pred_pm <- bfs_pred
  } else { stop("Select a valid model.") }
  coord_pred_pm <- coord_pred

  pm <- projection_matrices(mesh = mesh_wgt$mesh, wgt = mesh_wgt$wgt, loct = loct$loct, PS = PS_pm, bfs = bfs_pm, bfs_pred = bfs_pred_pm, coord_pred = coord_pred_pm)

  SPDE <- inla.spde2.pcmatern(mesh = mesh_wgt$mesh, alpha = 2,
                              prior.range = c(0.05 * diff(xlim), 0.01), # P(range < 0.05 * diff(xlim)) = 0.01
                              prior.sigma = c(1.00, 0.01))              # P(sigma > 1.00) = 0.01

  gaus_prior <- list(prior = "gaussian", param = c(0, 5)) # Prior for DP

  fs_formula <- full_stack_formula(Yobs = Yobs, pm = pm, SPDE = SPDE, gaus_prior = gaus_prior)
  
  partial_result <- model_fitting(fs_formula = fs_formula, compute_scores = TRUE, verbose = FALSE)
  
  if ((partial_result$fit$summary.hyperpar["Stdev for i", "mean"] > sig2_error) & (partial_result$fit$summary.hyperpar["Precision for the Gaussian observations", "mean"]) > 0.1) {
    others[[(valid_n + 1)]] <- list(fs_formula = fs_formula, bfs_pred = bfs_pred, coord_pred = coord_pred)
    results[[(valid_n + 1)]] <- partial_result
    valid_n <- valid_n + 1
  }
  
  print(paste("Number of valid sim.: ", sprintf("%03d", valid_n), sep = ""))
  cat("\n") 
  
  n <- n + 1
}

saveRDS(results, file = paste("FITTED_MODELS/", sprintf("%02d", scenario), "/", fitted_model, "-", smoothing_kernel, "-", partition_const, "-", n_center_points, "-", sample_size, ".rds", sep = ""))
saveRDS(others,  file = paste("FITTED_MODELS/", sprintf("%02d", scenario), "/OTHERS/", fitted_model, "-", smoothing_kernel, "-", partition_const, "-", n_center_points, "-", sample_size, ".rds", sep = ""))
