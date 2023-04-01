# Data generation functions

# Latent process
latt_generation <- function (xlim, ylim, by, nu, scl, sig2, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed = seed) 
  }
  
  x <- seq(from = xlim[1] + (by / 2), to = xlim[2] - (by / 2), by = by)
  y <- seq(from = ylim[1] + (by / 2), to = ylim[2] - (by / 2), by = by)
  coord <- expand.grid(x = x, y = y)
  
  # Define the covariance structure
  cov <- RMmatern(nu = nu, var = sig2, scale = scl) # + RMtrend(mean = mu)
  
  # Generate data
  process <- RFsimulate(model = cov, x = coord)
  process <- data.frame(x = process@coords[, 1], y = process@coords[, 2], z = process$variable1)
  
  # Create a gridded spatial object from "process"
  coordinates(process) <- ~ x + y
  gridded(process) <- TRUE
  
  raster(process)
}

# Preferentiality process
pref_generation <- function (xlim, ylim, by, f, maximum = 1, ...) {
  x <- seq(from = xlim[1] + (by / 2), to = xlim[2] - (by / 2), by = by)
  y <- seq(from = ylim[1] + (by / 2), to = ylim[2] - (by / 2), by = by)
  coord <- expand.grid(x = x, y = y)
  
  # Compute f((x, y)), for all x, y
  process <- coord
  process$z <- c(my_outer(X = x, Y = y, FUN = f))

  # Create a gridded spatial object from "process"
  coordinates(process) <- ~ x + y
  gridded(process) <- TRUE
  
  process <- raster(process)
  
  # Change the maximum (assuming data in [x, 1], for any x < 1)
  values(process) <- values(process) * maximum
  
  process
}

# Sample locations
sample_location <- function (latt, pref, alpha = 0, n_points = NULL, pref_sampling = T, seed = NULL, ...) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  lims <- extent(latt)
  xlim <- c(lims[1], lims[2])
  ylim <- c(lims[3], lims[4])
  
  # Compute the intensity process
  lambda <- latt
  if (pref_sampling) {
    values(lambda) <- exp(alpha + values(pref) * values(latt))
  } else {
    values(lambda) <- exp(alpha + values(pref)) # Review it (is it necessary to include the preferentiality process?)
  }
  
  # Generate points according to the intensity process
  if (!is.null(n_points)) {
    loct <- rpoint(n = n_points, f = as.im(lambda), win = owin(xrange = xlim, yrange = ylim))
  } else {
    loct <- rpoispp(lambda = as.im(lambda))
  }
  list(lambda = lambda, loct = loct)
}

# Observed process
observed_values = function (latt, loct, mu, sig2_error, seed = NULL, ...) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  values(latt) <- values(latt) + mu
  
  n_points <- loct$n
  
  Y <- rep(0, n_points)
  for (i in (1:n_points)) {
    Y[i] <- rnorm(n = 1, 
                  mean = raster::extract(latt, rbind(c(loct$x[i], loct$y[i]))),
                  sd = sqrt(sig2_error))
  }
  Y
}

# Implemented preferentiality functions 
select_preferentiality_function <- function (x, ...) {
  if (x == 1) { # Non-preferential sampling
    f <- function (x, y, ...) { 0 }
  } else if (x == 2) { # Preferential sampling with fixed degree of preferentiality
    f <- function (x, y, ...) { 1 }
  } else if (x == 3) { # Partition
    f <- function (x, y, ...) {
      xy <- cbind(x, y)
      z <- rep(x = 0, times = nrow(xy))
      for (i in 1:nrow(xy)) {
        if (xy[i, 2] <= mean(ylim)) { z[i] <- 0.50 } else if (xy[i, 2] > mean(ylim)) { z[i] <- 1.00 }
      }
      z
    }
  } else if (x == 4) { # Horizontal increasing
    f <- function (x, y, ...) { (x / diff(xlim)) }
  } else if (x == 5) { # Surface of second order (or second-degree surface)
    f <- function (x, y, ...) { 0.5 * ((x / diff(xlim)) ^ 2 + (y / diff(ylim)) ^ 2) }
    # f <- function (x, y, ...) { 0.5 * (1 + sin(2 * pi * 0.5 * ((x / diff(xlim)) + (y / diff(ylim))))) }
  } else {
    stop("Choose a valid function.")
  }
  f
}
