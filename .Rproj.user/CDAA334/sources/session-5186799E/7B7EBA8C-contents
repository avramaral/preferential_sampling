# Auxiliary functions

my_outer <- function (X, Y, FUN, ...) {
  coord <- expand.grid(x = X, y = Y)
  z <- do.call(what = FUN, args = list(x = coord$x, y = coord$y))
  if (length(z) == 1) {
    z <- rep(x = z, times = length(nrow(coord)))
  }
  z
}

check_points <- function (xy, xlim, ylim) {
  x_test <- (xy[, 1] >= xlim[1]) & (xy[, 1] <= xlim[2]) 
  y_test <- (xy[, 2] >= ylim[1]) & (xy[, 2] <= ylim[2])
  which(x_test & y_test)
}

Euclidean_norm <- function(x, ...) { sqrt(sum(x^2)) }

# Base functions

radial_basis <- function (x, y, center_x, center_y, bandwidth = 0.5, ...) {
  pars <- list(...)
  if (length(pars) != 0) {
    bandwidth <- ifelse(test = is.null(pars$bandwidth), yes = bandwidth, no = pars$bandwidth)
  }
  pts <- cbind((x - center_x), (y - center_y))
  k <- exp(-(apply(X = pts, MARGIN = 1, FUN = Euclidean_norm) ^ 2) / (2 * bandwidth))
  k <- k / max(k)
  k
}

Wendland <- function (x, y, center_x, center_y, bandwidth = 0.5, ...) {
  pars <- list(...)
  if (length(pars) != 0) {
    bandwidth <- ifelse(test = is.null(pars$bandwidth), yes = bandwidth, no = pars$bandwidth)
  }
  pts <- cbind((x - center_x), (y - center_y))
  d <- apply(X = pts, MARGIN = 1, FUN = Euclidean_norm) / bandwidth
  k <- (1 - d) ^ 6 * (35 * d ^ 2 + 18 * d + 3) / 3 * as.numeric(d <= 1)
  k <- k / max(k)
  k
}

triangular_basis <- function (x, y, center_x, center_y, bandwidth = 0.5, ...) {
  pars <- list(...)
  if (length(pars) != 0) {
    bandwidth <- ifelse(test = is.null(pars$bandwidth), yes = bandwidth, no = pars$bandwidth)
  }
  pts <- cbind((x - center_x), (y - center_y))
  k <- pmax(bandwidth - abs(x - center_x), 0) + pmax(bandwidth - abs(y - center_y), 0) 
  k <- k / max(k)
  k
}

unidirectional_triangular <- function (x, y, center_x, center_y, i = i, NS = TRUE, ...) {
  pars <- list(...)
  if (length(pars) != 0) {
    NS <- ifelse(test = is.null(pars$NS), yes = NS, no = pars$NS)
  }

  if (NS) {
    k <- rep(x = 0, times = length(y))
    for (j in 1:length(y)) {
      if ((i > 1) & (i < length(center_y))) {
        if ((center_y[(i - 1)] <= y[j]) & (y[j] <= center_y[i])) {
          k[j] <- (y[j] - center_y[(i - 1)]) / (center_y[i] - center_y[(i - 1)])
        } else if ((center_y[i] <= y[j]) & (y[j] <= center_y[(i + 1)])) {
          k[j] <- (y[j] - center_y[(i + 1)]) / (center_y[i] - center_y[(i + 1)])
        } else {
          k[j] <- 0 
        } 
      } else if (i == 1) {
        if (y[j] <= center_y[(i + 1)]) {
          k[j] <- (y[j] - center_y[(i + 1)]) / (center_y[i] - center_y[(i + 1)])
        } else {
          k[j] <- 0
        }
      } else {
        if (y[j] >= center_y[(i - 1)]) {
          k[j] <- (y[j] - center_y[(i - 1)]) / (center_y[i] - center_y[(i - 1)])
        } else {
          k[j] <- 0
        }
      }
    }
  } else {
    k <- rep(x = 0, times = length(x))
    for (j in 1:length(x)) {
      if ((i > 1) & (i < length(center_x))) {
        if ((center_x[(i - 1)] <= x[j]) & (x[j] <= center_x[i])) {
          k[j] <- (x[j] - center_x[(i - 1)]) / (center_x[i] - center_x[(i - 1)])
        } else if ((center_x[i] <= x[j]) & (x[j] <= center_x[(i + 1)])) {
          k[j] <- (x[j] - center_x[(i + 1)]) / (center_x[i] - center_x[(i + 1)])
        } else {
          k[j] <- 0 
        } 
      } else if (i == 1) {
        if (x[j] <= center_x[(i + 1)]) {
          k[j] <- (x[j] - center_x[(i + 1)]) / (center_x[i] - center_x[(i + 1)])
        } else {
          k[j] <- 0
        }
      } else {
        if (x[j] >= center_x[(i - 1)]) {
          k[j] <- (x[j] - center_x[(i - 1)]) / (center_x[i] - center_x[(i - 1)])
        } else {
          k[j] <- 0
        }
      }
    }
  }
  k
}

diagonal_triangular <- function (x, y, center_x, center_y, i = i, ...) {
  pars <- list(...)
  if (length(pars) != 0) { }
  
  compute_triangles <- function (x, y, center_x, center_y, i = i, NS, ...) {
    if (NS) {
      k <- rep(x = 0, times = length(y))
      for (j in 1:length(y)) {
        if ((i > 1) & (i < length(center_y))) {
          if ((center_y[(i - 1)] <= y[j]) & (y[j] <= center_y[i])) {
            k[j] <- (y[j] - center_y[(i - 1)]) / (center_y[i] - center_y[(i - 1)])
          } else if ((center_y[i] <= y[j]) & (y[j] <= center_y[(i + 1)])) {
            k[j] <- (y[j] - center_y[(i + 1)]) / (center_y[i] - center_y[(i + 1)])
          } else {
            k[j] <- 0 
          } 
        } else if (i == 1) {
          if (y[j] <= center_y[(i + 1)]) {
            k[j] <- (y[j] - center_y[(i + 1)]) / (center_y[i] - center_y[(i + 1)])
          } else {
            k[j] <- 0
          }
        } else {
          if (y[j] >= center_y[(i - 1)]) {
            k[j] <- (y[j] - center_y[(i - 1)]) / (center_y[i] - center_y[(i - 1)])
          } else {
            k[j] <- 0
          }
        }
      }
    } else {
      k <- rep(x = 0, times = length(x))
      for (j in 1:length(x)) {
        if ((i > 1) & (i < length(center_x))) {
          if ((center_x[(i - 1)] <= x[j]) & (x[j] <= center_x[i])) {
            k[j] <- (x[j] - center_x[(i - 1)]) / (center_x[i] - center_x[(i - 1)])
          } else if ((center_x[i] <= x[j]) & (x[j] <= center_x[(i + 1)])) {
            k[j] <- (x[j] - center_x[(i + 1)]) / (center_x[i] - center_x[(i + 1)])
          } else {
            k[j] <- 0 
          } 
        } else if (i == 1) {
          if (x[j] <= center_x[(i + 1)]) {
            k[j] <- (x[j] - center_x[(i + 1)]) / (center_x[i] - center_x[(i + 1)])
          } else {
            k[j] <- 0
          }
        } else {
          if (x[j] >= center_x[(i - 1)]) {
            k[j] <- (x[j] - center_x[(i - 1)]) / (center_x[i] - center_x[(i - 1)])
          } else {
            k[j] <- 0
          }
        }
      }
    }
    k
  }
  
  n_horizontal <- nrow(center_x[[1]])
  n_vertical   <- nrow(center_y[[2]])
  if (i <= n_horizontal) { j <- i; NS <- FALSE } else { j <- i - n_horizontal; NS <- TRUE }
  k <- compute_triangles(x = x, y = y, center_x = center_x[[ifelse(!NS, 1, 2)]][, 1], center_y = center_y[[ifelse(!NS, 1, 2)]][, 2], i = j, NS = NS)
  k
}

partition <- function (x, y, center_x, center_y, ...) {
  pars <- list(...)
  if (length(pars) != 0) {  }
  
  k <- rep(x = 0, times = length(x))
  for (j in 1:length(x)) {
    if (((center_x[1]<= x[j]) & (x[j] <= center_x[2])) & ((center_y[1] <= y[j]) & (y[j] <= center_y[2]))) {
      k[j] <- 1
    } else {
      k[j] <- 0
    }
  }
  k
}

partition_map <- function (x, y, i, ...) {
  pars <- list(...)
  if (length(pars) != 0) {
    map <- pars$map
  }
  
  if (is.na(map)[1]) { stop("Map not found.") }
  
  
  xy <- data.frame(x = x, y = y)
  coordinates(xy) <- ~ x + y
  
  xy <- st_as_sf(xy)
  st_crs(xy) <- st_crs(map)
  xy <- st_transform(xy)  
  
  intersection <- as.integer(st_intersects(xy, map[[i]]))
  
  k <- rep(x = 0, times = length(x))
  for (j in 1:length(x)) {
    if (!is.na(intersection[j])) {
      k[j] <- 1
    } else {
      k[j] <- 0
    }
  }
  k
}

# Plot base functions
plot_bfs <- function (bfs, coord, n_coord_pr = NA) {
  
  if (is.na(n_coord_pr)) { n_coord_pr <- length(bfs[[1]]) }
  
  n_base_functions <- length(bfs)
  
  bfs_all <- list()
  for (i in 1:n_base_functions) {
    bfs_all[[i]] <- data.frame(x = coord$x[1:n_coord_pr], y = coord$y[1:n_coord_pr], z = bfs[[i]][1:n_coord_pr])
    
    # Create a gridded spatial object from "bfs_all[[i]]"
    coordinates(bfs_all[[i]]) <- ~ x + y
    gridded(bfs_all[[i]]) <- TRUE
    
    bfs_all[[i]] <- raster(bfs_all[[i]])
    plot(bfs_all[[i]], main = paste("Base function ", sprintf('%02d', i), sep = ""))
  }

  bfs_all  
}
