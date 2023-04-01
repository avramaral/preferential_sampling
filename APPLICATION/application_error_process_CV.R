
scenarios <- 1:9 # :11
divs <- 1:9 # 1:9

for (scenario in scenarios) {
  error <- 0
  CRPS  <- 0
  SCRPS <- 0
  for (div in divs) {
    scores <- readRDS(file = paste("APPLICATION/FITTED_MODELS/CV/SCORE/scenario_", scenario, "_division_", div, ".rds", sep = ""))
    
    error <- error + scores$error
    CRPS  <- CRPS  + scores$CRPS
    SCRPS <- SCRPS + scores$SCRPS
  }
  
  if (TRUE) {
    print(paste("Scenario ", sprintf("%02d", scenario), sep = ""))
    print(paste("Error (latent): ", round(x = error / length(divs), digits = 4), sep = ""))
    print(paste("CRPS  (latent): ", round(x = CRPS  / length(divs), digits = 4), sep = ""))
    print(paste("SCRPS (latent): ", round(x = SCRPS / length(divs), digits = 4), sep = ""))
    cat("\n")
  }
}

