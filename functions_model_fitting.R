# Model fitting functions

# Create mesh and compute weights
mesh_weights <- function (latt, max_edge_mf = c(0.1, 0.25), offset_mf = c(0.1, 0.25), cutoff_mf = 0.05, PLOT = FALSE, ...) {
  
  if(!(class(latt) == "sf")[1]) {
    ext_latt <- extent(latt)
    coo_latt <- cbind(c(ext_latt[1], ext_latt[2], ext_latt[2], ext_latt[1], ext_latt[1]), 
                      c(ext_latt[3], ext_latt[3], ext_latt[4], ext_latt[4], ext_latt[3]))
    
    range_ext <- diff(c(ext_latt[1], ext_latt[2]))
    # Create the mesh
    latt_mesh <- inla.mesh.2d(loc.domain = coo_latt, max.edge = c(max_edge_mf[1] * range_ext, max_edge_mf[2] * range_ext), offset = c(offset_mf[1] * range_ext, offset_mf[2] * range_ext), cutoff = cutoff_mf * range_ext)  
  } else {
    coo_latt <- sf::st_coordinates(latt)
    coo_latt <- matrix(c(coo_latt[, 1], coo_latt[, 2]), ncol = 2)
    colnames(coo_latt) <- c("x", "y")
    
    latt_mesh <- inla.mesh.2d(loc.domain = coo_latt, max.edge = max_edge_mf, offset = offset_mf)
  }
  
  if (PLOT) {
    plot(latt_mesh, main = "Mesh")
    if(!(class(latt) == "sf")[1]) {
      plot(st_as_sf(as.owin(as.im(latt))), border = "green", lwd = 3, add = T)
    } else {
      plot(latt$geometry, col = NA, lwd = 3, border = "green", add = T)
    }
  }
  
  # Dual mesh
  tiles <- voronoi.polygons(SpatialPoints(latt_mesh$loc[, 1:2]))
  lattP <- SpatialPolygons(list(Polygons(list(Polygon(coo_latt)), ID = "1")))
  
  # Weights
  wgt <- sapply(1:length(tiles), function(p) {
    aux <- tiles[p, ]  
    if (gIntersects(aux, lattP)) {
      return(gArea(gIntersection(aux, lattP)))
    } else {
      return(0)
    }
  })
  
  
  if (PLOT) {
    plot(tiles, main = "Dual mesh")
    points(latt_mesh$loc, pch = 19)
    if(!(class(latt) == "sf")[1]) {
      plot(st_as_sf(as.owin(as.im(latt))), border = "green", lwd = 3, add = T)
    } else {
      plot(latt$geometry, col = NA, lwd = 3, border = "green", add = T)
    }
  }
  
  list(mesh = latt_mesh, wgt = wgt)
}

# Define the center points (Make it less rigid)
define_center_pts <- function (smoothing_kernel, partition_const = "non_constant", n_center_points = "few", xlim = c(0, 1), NS = TRUE, ...) {
  if (smoothing_kernel == "partition") {
    if (partition_const == "constant") {
      center_pts <- list(list(xlim = c(-0.25, 1.25), ylim = c(-0.25, 1.25)))  
    } else if (partition_const == "non_constant") {
      center_pts <- list(list(xlim = c(-0.25, 1.25), ylim = c(-0.25, 0.50)), list(xlim = c(-0.25, 1.25), ylim = c(0.50, 1.25))) 
    }
  } else if (smoothing_kernel == "unidirectional_triangular") {
    if (NS) {
      center_pts <- rbind(c(0.50, -0.25), c(0.50, 0.25), c(0.50, 0.50), c(0.50, 0.75), c(0.50, 1.25)) # Vertical
    } else {
      center_pts <- rbind(c(-0.25, 0.50), c(0.25, 0.50), c(0.50, 0.50), c(0.75, 0.50), c(1.25, 0.50)) # Horizontal  
    }
  } else if (smoothing_kernel == "Wendland") {
    if (n_center_points == "few") {
      center_pts <- rbind(c(0.25, 0.25), c(0.75, 0.25), c(0.25, 0.75), c(0.75, 0.75), c(0.50, 0.50))
    } else if (n_center_points == "many") {
      center_pts <- rbind(c(0.25, 0.25), c(0.75, 0.25), c(0.25, 0.75), c(0.75, 0.75), c(0.50, 0.50), c(0.25, 0.50), c(0.50, 0.25), c(0.50, 0.75), c(0.75, 0.50))
    } else if (n_center_points == "manymany") {
      center_pts <- rbind(c(0.25, 0.25), c(0.75, 0.25), c(0.25, 0.75), c(0.75, 0.75), c(0.50, 0.50), c(0.25, 0.50), c(0.50, 0.25), c(0.50, 0.75), c(0.75, 0.50), c(0.00, 0.00), c(0.00, 1.00), c(1.00, 0.00), c(1.00, 1.00))
    }
  } else if (smoothing_kernel == "radial_basis") {
    if (n_center_points == "few") {
      center_pts <- rbind(c(0.25, 0.25), c(0.75, 0.25), c(0.25, 0.75), c(0.75, 0.75), c(0.50, 0.50))
    } else if (n_center_points == "many") {
      center_pts <- rbind(c(0.25, 0.25), c(0.75, 0.25), c(0.25, 0.75), c(0.75, 0.75), c(0.50, 0.50), c(0.25, 0.50), c(0.50, 0.25), c(0.50, 0.75), c(0.75, 0.50))
    } else if (n_center_points == "manymany") {
      center_pts <- rbind(c(0.25, 0.25), c(0.75, 0.25), c(0.25, 0.75), c(0.75, 0.75), c(0.50, 0.50), c(0.25, 0.50), c(0.50, 0.25), c(0.50, 0.75), c(0.75, 0.50), c(0.00, 0.00), c(0.00, 1.00), c(1.00, 0.00), c(1.00, 1.00))
    }
  } else {
    stop("Choose a valid smoothing kernel.")
  }
  
  if (is.list(center_pts)) { center_pts <- lapply(X = center_pts, FUN = function (x) { lapply(X = x, FUN = function (y) { y * diff(xlim) }) }) } else { center_pts <- center_pts * diff(xlim) }
  
  center_pts
}

# Compute base functions 
base_functions <- function (center_pts, loct, mesh = NULL, smoothing_kernel = "gaussian_kernel", ...) {
  if (is.null(mesh)) {
    if (!class(loct)[1] == "matrix") {
      all_pts <- cbind(loct$x, loct$y)  
    } else {
      all_pts <- loct
    }
  } else {
    if (!class(loct)[1] == "matrix") {
      all_pts <- rbind(mesh$loc[, 1:2], cbind(loct$x, loct$y))
    } else {
      all_pts <- rbind(mesh$loc[, 1:2], loct)
    }
  }
  
  f <- get(smoothing_kernel)
  b <- list()
  if (smoothing_kernel == "diagonal_triangular") { n_base_functions <- nrow(center_pts[[1]]) + nrow(center_pts[[2]]) } else if (is.list(center_pts)) { n_base_functions <- length(center_pts) } else if (smoothing_kernel == "partition_map") { n_base_functions <- length(map) } else { n_base_functions <- nrow(center_pts) }
  for (i in 1:n_base_functions) {
    if (smoothing_kernel == "unidirectional_triangular") {
      b[[i]] <- f(x = all_pts[, 1], y = all_pts[, 2], center_x = center_pts[ , 1], center_y = center_pts[ , 2], i = i, ...)
    } else if (smoothing_kernel == "diagonal_triangular") {
      b[[i]] <- f(x = all_pts[, 1], y = all_pts[, 2], center_x = center_pts, center_y = center_pts, i = i, ...)
    } else if (smoothing_kernel == "partition") {
      b[[i]] <- f(x = all_pts[, 1], y = all_pts[, 2], center_x = center_pts[[i]]$xlim, center_y = center_pts[[i]]$ylim, ...)
    } 
    else if (smoothing_kernel == "partition_map") {
      b[[i]] <- f(x = all_pts[, 1], y = all_pts[, 2], i = i, ...)
    } else {
      b[[i]] <- f(x = all_pts[, 1], y = all_pts[, 2], center_x = center_pts[i, 1], center_y = center_pts[i, 2], ...)
    }
  }
  b
}

# Compute projection matrices
projection_matrices <- function (mesh, wgt, loct, PS = FALSE, bfs = NULL, bfs_pred = NULL, coord_pred = NULL, ...) {
  
  n_vtx <- mesh$n
  n_pts <- loct$n
  
  y_pp <- rep(x = 0:1, times = c(n_vtx, n_pts))
  e_pp <- c(wgt, rep(x = 0, times = n_pts))
  
  coor <- cbind(loct$x, loct$y)
  imat <- Diagonal(n = n_vtx, x = rep(x = 1, times = n_vtx)) # Projection matrix for the integration points
  ymat <- inla.spde.make.A(mesh = mesh, loc = coor)  # Projection matrix for the observed points 
  
  # Estimation
  pm <- list(y_pp = y_pp, e_pp = e_pp, ymat = ymat, imat = imat)
  if (PS) { # Preferential sampling
    A_pp_base <- rbind(imat, ymat)
    if (is.null(bfs)) { # Constant degree of preferentiality
      A_pp <- A_pp_base
    } else { # Spatially varying degree of preferentiality
      n_base_functions <- length(bfs)
      A_pp <- list()
      for(i in 1:n_base_functions) {
        A_pp[[i]] <- A_pp_base * bfs[[i]] # Multiply each column by the base function values
      }
    }
    pm$A_pp <- A_pp
  }
  
  # Prediction
  if (!is.null(coord_pred)) {
    pmat_base <- inla.spde.make.A(mesh = mesh, loc = as.matrix(coord_pred)) # Projection matrix for the prediction points
    pm$pmat_base <- pmat_base
    pm$pmat <- pmat_base
    if (PS & (!is.null(bfs))) { # Preferential sampling and spatially varying degree of preferentiality
      if (is.null(bfs_pred)) {
        stop("Please, provide the base functions evaluated at the prediction points.")
      } else {
        pmat_aux <- list()
        for (i in 1:n_base_functions) {
          pmat_aux[[i]] <- pmat_base * bfs_pred[[i]]
        }
        pm$pmat <- pmat_aux
      }
    } 
  }
  
  pm
}

# Stack and model formula
full_stack_formula <- function (Yobs, pm, SPDE = SPDE, gaus_prior, ...) { 
  y_pp <- pm$y_pp
  e_pp <- pm$e_pp
  ymat <- pm$ymat
  imat <- pm$imat
  A_pp <- pm$A_pp
  pmat_base <- pm$pmat_base
  pmat <- pm$pmat
  
  n_pts <- length(Yobs)
  n_vtx <- length(pm$y_pp) - n_pts
  if (is.list(pmat)) { n_pts_pred <- nrow(pmat[[1]]) } else { n_pts_pred <- nrow(pmat) }
  
  SPDE_index <- 1:n_vtx
  
  if (is.null(A_pp)) { # Non-preferential sampling
    
    non_ps_stk_e <- inla.stack(data = list(y = Yobs), A = list(1, ymat), effects = list(mu = rep(x = 1, times = n_pts), i = SPDE_index), tag = "est_y") 
    
    if (is.null(pmat)) { # No prediction
      
      full_stack <- inla.stack(non_ps_stk_e)
    
    } else { # Prediction
        
      non_ps_stk_p <- inla.stack(data = list(y = NA), A = list(1, pmat), effects = list(mu = rep(x = 1, times = n_pts_pred), i = SPDE_index), tag = "pred_y")
      full_stack <- inla.stack(non_ps_stk_e, non_ps_stk_p)
          
    }
    
    formula <- "y ~ 0 + mu + f(i, model = SPDE)"
    
  } else { # Preferential sampling
    
    ps_stk_y_e  <- inla.stack(data = list(y = cbind(Yobs, NA), e = rep(NA, n_pts)), A = list(1, ymat), effects = list(mu = rep(x = 1, times = n_pts), i = SPDE_index), tag = "est_y") 
    
    if (!is.list(A_pp)) { # Constant degree of preferentiality
      
      ps_stk_pp_e <- inla.stack(data = list(y = cbind(NA, y_pp), e = e_pp), A = list(1, A_pp), effects = list(alpha = rep(x = 1, times = n_vtx + n_pts), j = SPDE_index), tag = "est_pp")
      
      if (is.null(pmat)) { # No prediction
        
        full_stack <- inla.stack(ps_stk_y_e, ps_stk_pp_e)
        
      } else { # Prediction
        
        ps_stk_y_p  <- inla.stack(data = list(y = cbind(rep(x = NA, times = n_pts_pred), NA), e = rep(x = NA, times = n_pts_pred)), A = list(1, pmat), effects = list(mu = rep(x = 1, times = n_pts_pred), i = SPDE_index), tag = "pred_y") 
        ps_stk_pp_p <- inla.stack(data = list(y = cbind(NA, rep(x = NA, times = n_pts_pred)), e = rep(x = 1,  times = n_pts_pred)), A = list(1, pmat), effects = list(alpha = rep(x = 1, times = n_pts_pred), j = SPDE_index), tag = "pred_pp")
        
        full_stack <- inla.stack(ps_stk_y_e, ps_stk_pp_e, ps_stk_y_p, ps_stk_pp_p)
        
      }
      
      formula <- "y ~ 0 + mu + alpha + f(i, model = SPDE) + f(j, copy = \"i\", fixed = FALSE, hyper = list(beta = gaus_prior))"
      
    } else { # Spatially varying degree of preferentiality
      
      effects <- list(rep(x = 1, times = n_vtx + n_pts))
      latt_id <- list()
      for (i in 1:length(A_pp)) { latt_id <- append(latt_id, list(SPDE_index)) } 
      effects <- append(effects, latt_id)
      names(effects) <- c("alpha", paste("j_", 1:length(A_pp), sep = ""))
      
      ps_stk_pp_e <- inla.stack(data = list(y = cbind(NA, y_pp), e = e_pp), A = append(1, A_pp), effects = effects, tag = "est_pp")
      
      if (is.null(pmat)) { # No prediction
        
        full_stack <- inla.stack(ps_stk_y_e, ps_stk_pp_e)
        
      } else { # Prediction
        
        effects_pred <- effects
        effects_pred[[1]] <- rep(x = 1, times = n_pts_pred)
        
        ps_stk_y_p  <- inla.stack(data = list(y = cbind(rep(x = NA, times = n_pts_pred), NA), e = rep(x = NA, times = n_pts_pred)), A = list(1, pmat_base), effects = list(mu = rep(x = 1, times = n_pts_pred), i = SPDE_index), tag = "pred_y")
        ps_stk_pp_p <- inla.stack(data = list(y = cbind(NA, rep(x = NA, times = n_pts_pred)), e = rep(x =  1, times = n_pts_pred)), A = append(1, pmat), effects = effects_pred, tag = "pred_pp") 
        full_stack <- inla.stack(ps_stk_y_e, ps_stk_pp_e, ps_stk_y_p, ps_stk_pp_p)
      }
      
      random_effects <- paste(paste("f(j_",  1:length(A_pp), ", copy = \"i\", fixed = FALSE, hyper = list(beta = gaus_prior))", sep = ""), collapse = " + ")
      formula <- paste("y ~ 0 + mu + alpha + f(i, model = SPDE) + ", random_effects, sep = "")
      
    }
  }
  
  list(full_stack = full_stack, formula = formula, PS = ifelse(is.null(A_pp), F, T), varying_prefer = ifelse(is.list(A_pp), T, F), prediction = ifelse(is.null(pmat), F, T))
}

# Fitting model
model_fitting <- function (fs_formula, compute_scores = F, verbose = T, ...) {
  
  full_stack <- fs_formula$full_stack
  formula <- eval(parse(text = fs_formula$formula))
  PS <- fs_formula$PS
  
  start_time <- Sys.time()
  
  if (!PS) {
    if (!compute_scores) {
      fit <- inla(formula,
                  family = "gaussian",
                  data = inla.stack.data(full_stack), 
                  control.predictor = list(link = 1, compute = T, A = inla.stack.A(full_stack)),
                  verbose = verbose)
    } else {
      fit <- inla(formula,
                  family = "gaussian",
                  data = inla.stack.data(full_stack), 
                  control.predictor = list(link = 1, compute = T, A = inla.stack.A(full_stack)), 
                  control.compute = list(po = TRUE, cpo = TRUE, dic = TRUE, waic = TRUE),
                  verbose = verbose)
    }
    
  } else {
    if (!compute_scores) {
      fit <- inla(formula, 
                  family = c("gaussian", "poisson"), 
                  data = inla.stack.data(full_stack),
                  E = inla.stack.data(full_stack)$e,
                  control.predictor = list(link = 1, compute = T, A = inla.stack.A(full_stack)), 
                  verbose = verbose) 
    } else {
      fit <- inla(formula, 
                  family = c("gaussian", "poisson"), 
                  data = inla.stack.data(full_stack),
                  E = inla.stack.data(full_stack)$e,
                  control.predictor = list(link = 1, compute = T, A = inla.stack.A(full_stack)), 
                  control.compute = list(po = TRUE, cpo = TRUE, dic = TRUE, waic = TRUE),
                  verbose = verbose)       
    }
  }
  
  end_time <- Sys.time()
  time_taken <- end_time - start_time
  
  list(fit = fit, time_taken = time_taken, PS = PS, varying_pref = fs_formula$varying_pref, prediction = fs_formula$prediction)
}
