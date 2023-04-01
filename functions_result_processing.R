# Result processing

# Compute error
MSE <- function (orig, pred, ...) {
  if ((class(orig) == "RasterLayer") & (class(pred) == "RasterLayer")) {
    orig <- values(orig)
    pred <- values(pred)
  }
  (sum((orig - pred) ^ 2) / length(orig))
}

MAPE <- function (orig, pred, ...) {
  if ((class(orig) == "RasterLayer") & (class(pred) == "RasterLayer")) {
    orig <- values(orig)
    pred <- values(pred)
  }
  (sum(abs((orig - pred) / orig)) / length(orig))
}

CRPS <- function(orig, pred, s) {
  if ((class(orig) == "RasterLayer") & (class(pred) == "RasterLayer")) {
    orig <- values(orig)
    pred <- values(pred)
  }
  md <- orig - pred
  s / sqrt(pi) - 2 * s * dnorm(md / s) + md * (1 - 2 * pnorm(md / s))
}

SCRPS <- function(orig, pred, s) {
  if ((class(orig) == "RasterLayer") & (class(pred) == "RasterLayer")) {
    orig <- values(orig)
    pred <- values(pred)
  }
  md <- orig - pred
  -0.5 * log(2 * s / sqrt(pi)) - sqrt(pi) * (s * dnorm(md / s) - md / 2 + md * pnorm(md / s)) / s
}

compute_preferentiatility <- function (fit, coord_pred, bfs_pred, xlim, by, n_coord_pr = NA, ...) {
  
  if (is.na(n_coord_pr)) { n_coord_pr <- length(bfs_pred[[1]]) }
  
  n_base_functions <- length(bfs_pred)
  n_hyperparameter <- nrow(fit$summary.hyperpar)
  
  coeff <- fit$summary.hyperpar[((n_hyperparameter - n_base_functions + 1):n_hyperparameter), "mean"][1:n_coord_pr]
  value <- rep(x = 0, times = n_coord_pr) 
  for (i in 1:n_base_functions) {
    # bfs_pred_tmp <- c(matrix(data = bfs_pred[[i]], nrow = length(seq(from = xlim[1] + (by / 2), to = xlim[2] - (by / 2), by = by)), byrow = TRUE))
    bfs_pred_tmp <- bfs_pred[[i]][1:n_coord_pr]
    value <- value + (bfs_pred_tmp * coeff[i])
  }
  
  pref <- data.frame(x = coord_pred$x[1:n_coord_pr], y = coord_pred$y[1:n_coord_pr], z = value)
  
  # Create a gridded spatial object from "pref"
  coordinates(pref) <- ~ x + y
  gridded(pref) <- TRUE
  
  raster(pref)
}
