library(RandomFields)
library(raster)
library(spatstat)

data_generation = function (xlim, ylim, by, mu, nu, sig2, scl) {
  x = seq(from = xlim[1], to = xlim[2], by = by)
  y = seq(from = ylim[1], to = ylim[2], by = by)
  coord = expand.grid(x = x, y = y)
  
  # Define the covariance structure
  cov = RMmatern(nu = nu, var = sig2, scale = scl) + RMtrend(mean = mu)
  
  # Generate data
  process = RFsimulate(model = cov, x = coord)
  process = data.frame(x = process@coords[, 1], y = process@coords[, 2], z = process$variable1)
  
  # Created a gridded spatial object from 'process'
  coordinates(process) = ~ x + y
  gridded(process) = TRUE
  
  # Return the 'process' as a 'raster object'
  raster(process)
}

intensity = function (process, x, y, alpha, beta) {
  z = raster::extract(process, cbind(x, y))
  exp(alpha + beta * z)
}

sample_locations = function (process, n_points, alpha, beta, pref_sampling = TRUE, xlim = c(0, 1), ylim = c(0, 1)) {
  locations = data.frame(x = c(), y = c())
  pts = rpoispp(exp(alpha + beta * max(process@data@values)), win = as.owin(list(xrange = c(xlim[1], xlim[2]), yrange = c(ylim[1], ylim[2]))))
  # Assuming that there are sufficiently many points to be chosen
  count = 0
  if (pref_sampling) {
    for (i in (1:length(pts$x))) {
      if (rbinom(n = 1, size = 1, prob = (intensity(process, pts$x[i], pts$y[i], alpha = alpha, beta = beta) / exp(alpha + beta * max(process@data@values))))){
        locations = rbind(locations, data.frame(x = pts$x[i], y = pts$y[i]))   
        count = count + 1
      }
      if (count == n_points) { break }
    }
  } else {
    locations = data.frame(x = pts$x[1:n_points], y = pts$y[1:n_points])
  }
  locations
}

observed_values = function (process, locations, n_points, sig2_error) {
  Y = rep(0, n_points)
  for (i in (1:n_points)) {
    Y[i] = rnorm(n = 1, 
                 mean = raster::extract(process, rbind(c(locations$x[i], locations$y[i]))),
                 sd = sqrt(sig2_error))
  }
  Y
}

compute_mse = function (process1, process2, process3) {
  MSE1 = sum((process1@data@values - process2@data@values)^2) / length(process1@data@values)
  MSE2 = sum((process1@data@values - process3@data@values)^2) / length(process1@data@values)
  c(MSE1, MSE2)
}







