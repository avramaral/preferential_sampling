library(raster)
library(sp)
library(INLA)

####################
xlim <- c(0, 10)
ylim <- c(0, 10)
by <- 0.05
####################

sample_sizes <- c(100, 200, 500, 1000)

for (s in sample_sizes) {
    dat[[as.character(s)]] <- readRDS(file = paste("DATA/01/", s, ".rds", sep = "")) 
}

aa <- c()
for (i in 1:100) { aa <- c(aa, dat$`100`[[i]]$alpha) }
print(paste(round(mean(aa), 4), " (",round(sd(aa), 4), ")", sep = ""))

