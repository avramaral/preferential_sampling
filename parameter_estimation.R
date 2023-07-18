library("raster")
library("sp")
library("INLA")

####################
xlim <- c(0, 10)
ylim <- c(0, 10)
by <- 0.05
####################

sample_sizes <- c(100, 200, 500, 1000)
fitted_model <- c("non_PS", "const_DP_PS", "var_DP_PS")

count <- 1
res <- list()
for (s in sample_sizes) {
  res[[as.character(s)]] <- list()
  for (i in 1:3) {
    print(count)
    res[[as.character(s)]][[i]] <- readRDS(file = paste("FITTED_MODELS/01/", fitted_model[[i]], "-partition-constant-many-", s, ".rds", sep = "")) 
    count <- count + 1
  } 
}

# Fixed Effects
for (s in sample_sizes) {
  fixed_effects <- c(0, 0)
  for (i in 1:100) {
    fixed_effects <- fixed_effects + res[[as.character(s)]][[3]][[i]]$fit$summary.fixed[, 1] # c(mu, alpha)  
  }
  fixed_effects <- fixed_effects / 100
  print(round(fixed_effects, 2))
}

# Hyperpar
for (s in sample_sizes) {
  hyperpar <- c(0, 0, 0)
  for (i in 1:100) {
    tmp <- res[[as.character(s)]][[3]][[i]]$fit$summary.hyperpar[1:3, 1] # (sigma_error, range, sigma_zeta)
    tmp[1] <- 1 / tmp[1]
    tmp[3] <- tmp[3] ^ 2    
    hyperpar <- hyperpar + tmp
  }
  hyperpar <- hyperpar / 100
  print(round(hyperpar, 2))
}

# Betas
for (s in sample_sizes) {
  betas <- c(0, 0, 0, 0, 0)
  for (i in 1:100) {
    betas <- + betas + res[[as.character(s)]][[3]][[i]]$fit$summary.hyperpar[4:8, 1] # (betas)
  }
  betas <- betas / 100
  print(round(betas, 2))
}


