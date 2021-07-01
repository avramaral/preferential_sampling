# setwd(dir = dirname(rstudioapi::getActiveDocumentContext()$path)) # It only works on RStudio
rm(list = ls())

library(RandomFields)
library(raster)
library(spatstat)
library(base)
library(xlsx)
library(INLA)
library(prodlim)

source('main_functions.R')

# General parameters
n_iter = 50
possible_parameters = c(1, 2, 3, 4)
betas = c(0, 0.5)
sample_size = c(100, 500)

# Common parameters
sig2_error = 1
alpha      = 10 
xlim       = c(0, 1)
ylim       = c(0, 1)
by         = 0.01

for (beta in betas) {
  for (n_points in sample_size) {
    for (choice_par in possible_parameters) {
        if (choice_par == 1) {
          mu   = 0
          sig2 = 4
          nu   = 1
          scl  = 0.05
        } else if (choice_par == 2) {
          mu   = 0
          sig2 = 4
          nu   = 0.5
          scl  = 0.5
        } else if (choice_par == 3) {
          mu   = 0
          sig2 = 1
          nu   = 1
          scl  = 0.05
        } else if (choice_par == 4) {
          mu   = 0
          sig2 = 1
          nu   = 0.5
          scl  = 0.5
        }
        
        MSE = data.frame(matrix(data = NA, nrow = n_iter, ncol = 2))
        colnames(MSE) = c('MSE INLA', 'MSE INLA PS') 
        
        for(iter in (1:n_iter)) {
          print(paste('beta:', beta, 'npoints:', n_points, 'choicepar:', choice_par, 'iter:', iter))
         
          ######################################
          #           Data Generation          #
          ######################################
          
          original_process = data_generation(xlim, ylim, by, mu, nu, sig2, scl)
          locations = sample_locations(original_process, n_points, alpha = alpha, beta = beta, xlim = xlim, ylim = ylim)
          Y = observed_values(original_process, locations, n_points, sig2_error)
      
          ######################################
          #                INLA                #
          ######################################
      
          par(mfrow = c(1, 1), family = "LM Roman 10")
          x = seq(from = xlim[1], to = xlim[2], by = by)
          y = seq(from = ylim[1], to = ylim[2], by = by)
          coord = expand.grid(x = x, y = y)
      
          ######################################
          # INLA WITHOUT PREFERENTIAL SAMPLING #
          ######################################
          mesh = inla.mesh.2d(loc = as.matrix(coord),
                              offset = c(0.25, 0.25),
                              cutoff = 0.05,
                              max.edge = c(0.1, 0.5))
          
          spde = inla.spde2.matern(mesh = mesh, alpha = 2) # alpha = nu + (d / 2)
          indx = seq(from = 1, to = spde$mesh$n, by = 1) 
      
          A = inla.spde.make.A(mesh = mesh, loc = as.matrix(locations))
          Ap = inla.spde.make.A(mesh = mesh, loc = as.matrix(coord))
      
          stk_e = inla.stack( # Data used to estimate the parameters
            tag = 'est',
            data = list(y = Y),
            A = list(1, A),
            effects = list(data.frame(mu = rep(1, nrow(locations))), x = indx)
          )
      
          stk_p = inla.stack( # Prediction
            tag = 'pred',
            data = list(y = NA),
            A = list(1, Ap),
            effects = list(data.frame(mu = rep(1, nrow(coord))), x = indx)
          )
      
          stk_full = inla.stack(stk_e, stk_p)
      
          # Fit the model
          M2 = inla(y ~ 0 + mu + f(x, model = spde),
                    data = inla.stack.data(stk_full),
                    control.predictor = list(compute = TRUE, A = inla.stack.A(stk_full)))
      
          index = inla.stack.index(stk_full, tag = 'pred')$data
          pred_mean = M2$summary.fitted.values[index, 'mean'] # Predicted values for the given locations
      
          INLAwithoutPS = data.frame(x = coord[, 1], y = coord[, 2], z = pred_mean)
      
          # Created a gridded spatial object from 'INLAwithoutPS'
          coordinates(INLAwithoutPS) = ~ x + y
          gridded(INLAwithoutPS) = TRUE
      
          INLAwithoutPS = raster(INLAwithoutPS)
      
          ######################################
          # INLA WITH    PREFERENTIAL SAMPLING #
          ######################################
          new_coord = unique(rbind(coord, locations))
      
          mesh = inla.mesh.2d(loc = as.matrix(new_coord),
                              offset = c(0.25, 0.25),
                              cutoff = 0.05,
                              max.edge = c(0.1, 0.5))
      
          spde = inla.spde2.matern(mesh = mesh, alpha = 2) # alpha = nu + (d / 2)
      
          i_loc = mesh$idx$loc 
      
          ig = c(i_loc, rep(NA, nrow(new_coord)))
          ip = c(rep(NA, nrow(new_coord)), i_loc)
          mu_inla = c(rep(1, nrow(new_coord)), rep(0, nrow(new_coord)))
          alpha_inla = c(rep(0, nrow(new_coord)), rep(1, nrow(new_coord)))
      
          y = rep(NA, nrow(new_coord)) # Create an empty vector
          samp_rows = row.match(locations, new_coord) # Find the sampling locations in 'new_coord'
          y[samp_rows] = Y # Store the observed values in the 'y' vector
          y = c(y, rep(NA, nrow(new_coord)))
          p = rep(0, nrow(new_coord)) # Create an empty vector
          p[samp_rows] = 1 # 'p' stands for points
          p = c(rep(NA, nrow(new_coord)), p)
          y = matrix(c(y, p), nrow(new_coord) * 2, 2)
      
          data = list(y = y, mu = mu_inla, alpha = alpha_inla, ig = ig, ip = ip)
      
          M3 = inla(formula = y ~ 0 + mu + alpha + f(ig, model = spde) + f(ip, copy = 'ig', fixed = FALSE),
                     family = c('gaussian', 'poisson'),
                     data = data,
                     control.predictor = list(compute = TRUE, link = 1), verbose = FALSE) #
      
          pred_mean = M3$summary.fitted.values[(1:(nrow(coord))), 'mean']
      
          INLAwithPS = data.frame(x = coord[, 1], y = coord[, 2], z = pred_mean)
      
          # Created a gridded spatial object from 'INLAwithPS
          coordinates(INLAwithPS) = ~ x + y
          gridded(INLAwithPS) = TRUE
      
          INLAwithPS = raster(INLAwithPS) 
          
          ######################################
          #                 MSE                #
          ######################################
          
          MSE[iter, ] = compute_mse(original_process, INLAwithoutPS, INLAwithPS)
          print(MSE)
        }
        
      # To save the MSE files, uncomment the following line and set an appropriate file path  
      # write.csv(x = MSE, file = paste0('path/MSE_beta', gsub('\\.', '@', toString(beta)), '_npoints', toString(n_points), '_choicepar', toString(choice_par), '.csv'), row.names = FALSE)  
    }
  }
}
