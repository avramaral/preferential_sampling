
scenarios <- 1:9 # :11

for (scenario in scenarios) {
  scores <- readRDS(file =  paste("APPLICATION/IMAGES/", sprintf("%02d", scenario), "/scores.rds", sep = ""))
  
  if (TRUE) {
    print(paste("Scenario ", sprintf("%02d", scenario), sep = ""))
    print(paste("Error (latent): ", round(x = scores$error, digits = 4), sep = ""))
    print(paste("CRPS  (latent): ", round(x = scores$CRPS,  digits = 4), sep = ""))
    print(paste("SCRPS (latent): ", round(x = scores$SCRPS, digits = 4), sep = ""))
    cat("\n")
  }
}

