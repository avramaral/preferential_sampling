library(plot3D)
library(raster)
library(rasterVis)
library(sp)
library(sf) ###
library(lattice)
library(latticeExtra)
library(gridExtra)
library(INLA)
library(tidyverse) ###
library(plotrix)

####################################################################################################
# Figure 01
####################################################################################################

scenario <- 4
sample_size <- 100

data <- readRDS(file = paste("DATA/", sprintf("%02d", scenario), "/", sample_size, ".rds", sep = ""))

for(n in 1:100) {
  # n <- 1
  latt <- data[[n]]$latt; pref <- data[[n]]$pref; loct <- data[[n]]$loct; Yobs <- data[[n]]$Yobs; orig <- data[[n]]$orig;
  
  jpeg(filename = paste("plots/simulated-example_", n, ".jpeg", sep = ""), width = 5000, height = 1100, res = 300)
  par(mfrow = c(1, 4), family = 'LM Roman 10', mar = c(0, 0, 0, 0) + 0.05)
  if (TRUE) {
    font_list <- list(axis.text = list(fontfamily = "LM Roman 10"),
                      par.xlab.text = list(fontfamily = "LM Roman 10"),
                      par.ylab.text = list(fontfamily = "LM Roman 10"),
                      par.main.text = list(fontfamily = "LM Roman 10"),
                      par.sub.text  = list(fontfamily = "LM Roman 10"))
    
    pts <- data.frame(x = loct$loct$x, y = loct$loct$y, name = as.character(1:loct$loct$n))
    coordinates(pts) <- ~ x + y
    
    # Latent Field
    pal <- jet.col(n = 100, alpha = 1)
    cuts <- seq(round(min(values(latt)) - 0.01, 2), round(max(values(latt)) + 0.01, 2), length.out = 100)
    labs <- seq(round(min(values(latt)) - 0.01, 2), round(max(values(latt)) + 0.01, 2), length.out = 5)
    p1 <- levelplot(latt, main = list(label = '(1)', side = 1, line = 0.5, cex = 1.5),
                    col.regions = pal, 
                    at = cuts, 
                    breaks = cuts, 
                    margin = F,
                    scales = list(x = list(at = seq(0, 10, 2), cex = 1.5),
                                  y = list(at = seq(0, 10, 2), cex = 1.5)),
                    colorkey = list(at = cuts, 
                                   labels = list(at = labs, cex = 1.5,
                                                 labels = as.character(format(labs, nsmall = 3)))),
                    par.settings = font_list)
    p1
    # Preferentiality
    pal <- jet.col(n = 100, alpha = 1)
    cuts <- seq(round(min(values(pref)) - 0.01, 2), round(max(values(pref)) + 0.01, 2), length.out = 100)
    labs <- seq(round(min(values(pref)) - 0.01, 2), round(max(values(pref)) + 0.01, 2), length.out = 6)
    p2 <- levelplot(pref, main = list(label = '(2)', side = 1, line = 0.5, cex = 1.5),
                    col.regions = pal, 
                    at = cuts, 
                    breaks = cuts, 
                    margin = F,
                    scales = list(x = list(at = seq(0, 10, 2), cex = 1.5),
                                  y = list(at = seq(0, 10, 2), cex = 1.5)),
                    colorkey = list(at = cuts, 
                                    labels = list(at = labs, cex = 1.5,
                                                  labels = as.character(format(labs, nsmall = 3)))),
                    par.settings = font_list)
    
    # Intensity
    
    pal <- jet.col(n = 100, alpha = 1)
    cuts <- seq(round(min(values(loct$lambda)) - 0.01, 2), round(max(values(loct$lambda)) + 0.01, 2), length.out = 100)
    labs <- seq(round(min(values(loct$lambda)) - 0.01, 2), round(max(values(loct$lambda)) + 0.01, 2), length.out = 5)
    p3 <- levelplot(loct$lambda, main = list(label = '(3)', side = 1, line = 0.5, cex = 1.5), 
                    col.regions = pal, 
                    at = cuts, 
                    breaks = cuts, 
                    margin = F,
                    scales = list(x = list(at = seq(0, 10, 2), cex = 1.5),
                                  y = list(at = seq(0, 10, 2), cex = 1.5)),
                    colorkey = list(at = cuts, 
                                    labels = list(at = labs, cex = 1.5,
                                                  labels = as.character(format(labs, nsmall = 3)))),
                    par.settings = font_list)
    
    # Latent Field
    pal <- jet.col(n = 100, alpha = 1)
    cuts <- seq(round(min(values(orig)) - 0.01, 2), round(max(values(orig)) + 0.01, 2), length.out = 100)
    labs <- seq(round(min(values(orig)) - 0.01, 2), round(max(values(orig)) + 0.01, 2), length.out = 6)
    p4 <- levelplot(orig, main = list(label = '(4)', side = 1, line = 0.5, cex = 1.5),
                    col.regions = pal, 
                    at = cuts, 
                    breaks = cuts, 
                    margin = F,
                    scales = list(x = list(at = seq(0, 10, 2), cex = 1.5),
                                  y = list(at = seq(0, 10, 2), cex = 1.5)),
                    colorkey = list(at = cuts, 
                                    labels = list(at = labs, cex = 1.5,
                                                  labels = as.character(format(labs, nsmall = 3)))),
                    par.settings = font_list) + layer(sp.points(pts, pch = 19, cex = 1, col = 1))
  }
  grid.arrange(p1, p2, p3, p4, ncol = 4)
  dev.off()
}

####################################################################################################
# Figure 02
####################################################################################################

source("functions_result_processing.R")

####################
xlim <- c(0, 10)
ylim <- c(0, 10)
by <- 0.05
####################

dat4 <- readRDS(file = "DATA/04/1000.rds")
dat5 <- readRDS(file = "DATA/05/1000.rds")

res4 <- readRDS(file = "FITTED_MODELS/04/var_DP_PS-unidirectional_triangular-non_constant-many-1000.rds")
res5 <- readRDS(file = "FITTED_MODELS/05/var_DP_PS-Wendland-non_constant-few-1000.rds")

oth4 <- readRDS(file = "FITTED_MODELS/04/OTHERS/var_DP_PS-unidirectional_triangular-non_constant-many-1000.rds")
oth5 <- readRDS(file = "FITTED_MODELS/05/OTHERS/var_DP_PS-Wendland-non_constant-few-1000.rds")

# for (n in 1:20) {
#   est_pref4 <- compute_preferentiatility(fit = res4[[n]]$fit, coord_pred = oth4[[n]]$coord_pred, bfs_pred = oth4[[n]]$bfs_pred, xlim = xlim, by = by)
#   
#   par(mfrow = c(1, 2))
#   plot(dat4[[n]]$pref, main = n)
#   plot(est_pref4, main = n)
#   par(mfrow = c(1, 1))
# }
# 
# for (n in 41:60) {
#   est_pref5 <- compute_preferentiatility(fit = res5[[n]]$fit, coord_pred = oth5[[n]]$coord_pred, bfs_pred = oth5[[n]]$bfs_pred, xlim = xlim, by = by)
#   
#   par(mfrow = c(1, 2))
#   plot(dat5[[n]]$pref, main = n)
#   plot(est_pref5, main = n)
#   par(mfrow = c(1, 1))
# }


font_list <- list(axis.text = list(fontfamily = "LM Roman 10"),
                  par.xlab.text = list(fontfamily = "LM Roman 10"),
                  par.ylab.text = list(fontfamily = "LM Roman 10"),
                  par.main.text = list(fontfamily = "LM Roman 10"),
                  par.sub.text  = list(fontfamily = "LM Roman 10"))

##################################################
##################################################

n4 <- 15

est_pref4 <- compute_preferentiatility(fit = res4[[n4]]$fit, coord_pred = oth4[[n4]]$coord_pred, bfs_pred = oth4[[n4]]$bfs_pred, xlim = xlim, by = by)

min4 <- min(values(dat4[[n4]]$pref), values(est_pref4))
max4 <- max(values(dat4[[n4]]$pref), values(est_pref4))

n5 <- 9

est_pref5 <- compute_preferentiatility(fit = res5[[n5]]$fit, coord_pred = oth5[[n5]]$coord_pred, bfs_pred = oth5[[n5]]$bfs_pred, xlim = xlim, by = by)

min5 <- min(values(dat5[[n5]]$pref), values(est_pref5))
max5 <- max(values(dat5[[n5]]$pref), values(est_pref5))

mn <- min(min4, min5)
mx <- max(max4, max5)

##################################################
##################################################

if (TRUE) {
  #jpeg(filename = "estimated-pref.jpeg", width = 5000, height = 1100, res = 300)
  par(mfrow = c(1, 4), family = 'LM Roman 10', mar = c(0, 0, 0, 0) + 0.05)
  
  pal <- jet.col(n = 100, alpha = 1)
  cuts <- seq(round(min4 - 0.00, 2), round(max4 + 0.00, 2), length.out = 100)
  labs <- seq(round(min4 - 0.00, 2), round(max4 + 0.00, 2), length.out = 6)
  p1 <- levelplot(dat4[[n4]]$pref, main = list(label = 'True (DS 04)', side = 1, line = 0.5, cex = 1.5),
                  col.regions = pal, 
                  at = cuts, 
                  breaks = cuts, 
                  margin = F,
                  scales = list(x = list(at = seq(0, 10, 2), cex = 1.5),
                                y = list(at = seq(0, 10, 2), cex = 1.5)),
                  colorkey = list(at = cuts, 
                                  labels = list(at = labs, cex = 1.5,
                                                labels = as.character(format(labs, nsmall = 3)))),
                  par.settings = font_list)
  
  p2 <- levelplot(est_pref4, main = list(label = 'Estimated (DS 04)', side = 1, line = 0.5, cex = 1.5),
                  col.regions = pal, 
                  at = cuts, 
                  breaks = cuts, 
                  margin = F,
                  scales = list(x = list(at = seq(0, 10, 2), cex = 1.5),
                                y = list(at = seq(0, 10, 2), cex = 1.5)),
                  colorkey = list(at = cuts, 
                                  labels = list(at = labs, cex = 1.5,
                                                labels = as.character(format(labs, nsmall = 3)))),
                  par.settings = font_list)
  
  ##################################################
  ##################################################
  
  pal <- jet.col(n = 100, alpha = 1)
  cuts <- seq(round(min5 - 0.01, 2), round(max4 + 0.00, 2), length.out = 100)
  labs <- seq(round(min5 - 0.01, 2), round(max4 + 0.00, 2), length.out = 6)
  p3 <- levelplot(dat5[[n5]]$pref, main = list(label = 'True (DS 05)', side = 1, line = 0.5, cex = 1.5),
                  col.regions = pal, 
                  at = cuts, 
                  breaks = cuts, 
                  margin = F,
                  scales = list(x = list(at = seq(0, 10, 2), cex = 1.5),
                                y = list(at = seq(0, 10, 2), cex = 1.5)),
                  colorkey = list(at = cuts, 
                                  labels = list(at = labs, cex = 1.5,
                                                labels = as.character(format(labs, nsmall = 3)))),
                  par.settings = font_list)
  
  p4 <- levelplot(est_pref5, main = list(label = 'Estimated (DS 05)', side = 1, line = 0.5, cex = 1.5),
                  col.regions = pal, 
                  at = cuts, 
                  breaks = cuts, 
                  margin = F,
                  scales = list(x = list(at = seq(0, 10, 2), cex = 1.5),
                                y = list(at = seq(0, 10, 2), cex = 1.5)),
                  colorkey = list(at = cuts, 
                                  labels = list(at = labs, cex = 1.5,
                                                labels = as.character(format(labs, nsmall = 3)))),
                  par.settings = font_list)
  
  grid.arrange(p1, p2, p3, p4, ncol = 4)
  #dev.off()
}

####################################################################################################
# Figure 03
####################################################################################################

scenarios <- c(1, 2, 3, 4, 5)
fitted_model <- "var_DP_PS"
smoothing_kernel <- c("partition", "partition", "partition", "unidirectional_triangular", "Wendland")
partition_const  <- c("constant", "constant", "non_constant", "non_constant", "non_constant")
n_center_points  <- c(rep("many", 4), "few")
sample_sizes <- c(100, 200, 500, 1000)


errors_df <- as.data.frame(matrix(data = 0, nrow = 2000, ncol = 3))
errors <- list()
cc <- 1
for (sample_size in sample_sizes) {
  errors[[as.character(sample_size)]] <- list()
  for (i in 1:5) {
    errors[[as.character(sample_size)]][[i]] <- readRDS(file = paste("FITTED_MODELS/", sprintf("%02d", scenarios[i]), "/ERRORS/UNNORMALIZED_", fitted_model, "-", smoothing_kernel[i], "-", partition_const[i], "-", n_center_points[i], "-", sample_size, ".rds", sep = ""))   
    for (j in 1:100) {
      errors_df[cc, ] <- c(sample_size, i, errors[[as.character(sample_size)]][[i]][[j]]$error_preferentiality)
      cc <- cc + 1
    }
  }
}
colnames(errors_df) <- c("sample_size", "data_scenario", "error")
errors_df <- errors_df[errors_df$error < 40, ] # Exclude non-sense outlier

jpeg(filename = "boxplot-pref.jpeg", width = 5000, height = 1100, res = 300)
at <- c(1:4, 5.5:8.5, 10:13, 14.5:17.5, 19:22)
par(mfrow = c(1, 1), family = 'LM Roman 10', mar = c(3.5, 4, 0, 0) + 0.1)
boxplot(error ~ sample_size + data_scenario, data = errors_df, xlim = c(1.3, 21.8), ylim = c(0, 3), axes = T, xaxt = "n", yaxt = "n", at = at, xlab = "", ylab = "Squared Error (SE)", cex.lab = 1.5)
axis(1, at = at, labels = rep(c("100", "200", "500", "1,000"), 5), cex.axis = 1.5)
axis(2, at = 0:5, labels = 0:5, cex.axis = 1.5)
text(2.5,  -0.75, "Data scenario 01", xpd = NA, cex = 1.5)
text(7.0,  -0.75, "Data scenario 02", xpd = NA, cex = 1.5) 
text(11.5, -0.75, "Data scenario 03", xpd = NA, cex = 1.5) 
text(16.0, -0.75, "Data scenario 04", xpd = NA, cex = 1.5) 
text(20.5, -0.75, "Data scenario 05", xpd = NA, cex = 1.5) 
dev.off()

####################################################################################################
# Figure 04
####################################################################################################

scenarios <- c(4, 6, 5, 7)
order_scenarios <- c(1, 2, 3, 4)
fitted_model <- "var_DP_PS"
smoothing_kernel <- c("unidirectional_triangular", "Wendland", "Wendland", "unidirectional_triangular")
partition_const  <- c("non_constant", "non_constant", "non_constant", "non_constant")
n_center_points  <- c("many", "few", "few", "many")
sample_sizes <- c(100, 200, 500, 1000)


errors_df <- as.data.frame(matrix(data = 0, nrow = 1600, ncol = 3))
errors <- list()
cc <- 1
for (sample_size in sample_sizes) {
  errors[[as.character(sample_size)]] <- list()
  for (i in 1:4) {
    errors[[as.character(sample_size)]][[i]] <- readRDS(file = paste("FITTED_MODELS/", sprintf("%02d", scenarios[i]), "/ERRORS/UNNORMALIZED_", fitted_model, "-", smoothing_kernel[i], "-", partition_const[i], "-", n_center_points[i], "-", sample_size, ".rds", sep = ""))   
    for (j in 1:100) {
      errors_df[cc, ] <- c(sample_size, order_scenarios[i], errors[[as.character(sample_size)]][[i]][[j]]$error_preferentiality)
      cc <- cc + 1
    }
  }
}
colnames(errors_df) <- c("sample_size", "data_scenario", "error")

jpeg(filename = "boxplot-pref-miss.jpeg", width = 4175, height = 1100, res = 300)
at <- c(1:4, 5.5:8.5, 11.5:14.5, 16.0:19.0)
par(mfrow = c(1, 1), family = 'LM Roman 10', mar = c(5.5, 4, 0, 0) + 0.1)
boxplot(error ~ sample_size + data_scenario, data = errors_df, xlim = c(1.2, 18.8), ylim = c(0, 3), axes = T, xaxt = "n", yaxt = "n", at = at, xlab = "", ylab = "Squared Error (SE)", cex.lab = 1.5)
axis(1, at = at, labels = rep(c("100", "200", "500", "1,000"), 4), cex.axis = 1.5)
axis(2, at = 0:5, labels = 0:5, cex.axis = 1.5)
text(2.5,  -0.9, "Hor. Uni. Tri.", xpd = NA, cex = 1.5)
text(7.0,  -0.9, "RB Wendland", xpd = NA, cex = 1.5) 
text(13.0, -0.9, "Hor. Uni. Tri.", xpd = NA, cex = 1.5)
text(17.5, -0.9, "RB Wendland", xpd = NA, cex = 1.5) 
text(4.75, -1.3, "Data Scenario 04", xpd = NA, cex = 1.5) 
text(15.25,-1.3, "Data Scenario 05", xpd = NA, cex = 1.5) 
dev.off()

####################################################################################################
# Figure 05
####################################################################################################

source("functions_result_processing.R")

####################
xlim <- c(0, 10)
ylim <- c(0, 10)
by <- 0.05
####################

dat4 <- readRDS(file = "DATA/04/1000.rds")
dat5 <- readRDS(file = "DATA/05/1000.rds")

res4 <- readRDS(file = "FITTED_MODELS/04/var_DP_PS-unidirectional_triangular-non_constant-many-1000.rds")
res5 <- readRDS(file = "FITTED_MODELS/05/var_DP_PS-Wendland-non_constant-few-1000.rds")
res6 <- readRDS(file = "FITTED_MODELS/06/var_DP_PS-Wendland-non_constant-few-1000.rds")
res7 <- readRDS(file = "FITTED_MODELS/07/var_DP_PS-unidirectional_triangular-non_constant-many-1000.rds")

oth4 <- readRDS(file = "FITTED_MODELS/04/OTHERS/var_DP_PS-unidirectional_triangular-non_constant-many-1000.rds")
oth5 <- readRDS(file = "FITTED_MODELS/05/OTHERS/var_DP_PS-Wendland-non_constant-few-1000.rds")
oth6 <- readRDS(file = "FITTED_MODELS/06/OTHERS/var_DP_PS-Wendland-non_constant-few-1000.rds")
oth7 <- readRDS(file = "FITTED_MODELS/07/OTHERS/var_DP_PS-unidirectional_triangular-non_constant-many-1000.rds")

##################################################
##################################################

font_list <- list(axis.text = list(fontfamily = "LM Roman 10"),
                  par.xlab.text = list(fontfamily = "LM Roman 10"),
                  par.ylab.text = list(fontfamily = "LM Roman 10"),
                  par.main.text = list(fontfamily = "LM Roman 10"),
                  par.sub.text  = list(fontfamily = "LM Roman 10"))

##################################################
##################################################

n4 <- 15

est_pref4 <- compute_preferentiatility(fit = res4[[n4]]$fit, coord_pred = oth4[[n4]]$coord_pred, bfs_pred = oth4[[n4]]$bfs_pred, xlim = xlim, by = by)

min4 <- min(values(dat4[[n4]]$pref), values(est_pref4))
max4 <- max(values(dat4[[n4]]$pref), values(est_pref4))

n5 <- 9

est_pref5 <- compute_preferentiatility(fit = res5[[n5]]$fit, coord_pred = oth5[[n5]]$coord_pred, bfs_pred = oth5[[n5]]$bfs_pred, xlim = xlim, by = by)

min5 <- min(values(dat5[[n5]]$pref), values(est_pref5))
max5 <- max(values(dat5[[n5]]$pref), values(est_pref5))

n6 <- 15

est_pref6 <- compute_preferentiatility(fit = res6[[n6]]$fit, coord_pred = oth6[[n6]]$coord_pred, bfs_pred = oth6[[n6]]$bfs_pred, xlim = xlim, by = by)

min6 <- min(values(dat4[[n4]]$pref), values(est_pref6))
max6 <- max(values(dat4[[n4]]$pref), values(est_pref6))


n7 <- 31

est_pref7 <- compute_preferentiatility(fit = res7[[n7]]$fit, coord_pred = oth7[[n7]]$coord_pred, bfs_pred = oth7[[n7]]$bfs_pred, xlim = xlim, by = by)

min7 <- min(values(dat5[[n5]]$pref), values(est_pref7))
max7 <- max(values(dat5[[n5]]$pref), values(est_pref7))

mn <- min(min4, min5, min6, min7)
mx <- max(max4, max5, max6, max7)

##################################################
##################################################

if (TRUE) {
  jpeg(filename = "estimated-pref-miss.jpeg", width = 4000, height = 2200, res = 300)
  par(mfrow = c(1, 4), family = 'LM Roman 10', mar = c(0, 0, 0, 0) + 0.05)
  
  pal <- jet.col(n = 100, alpha = 1)
  cuts <- seq(round(min(min4, min6) - 0.00, 2), round(max(max4, max6) + 0.01, 2), length.out = 100)
  labs <- seq(round(min(min4, min6) - 0.00, 2), round(max(max4, max6) + 0.01, 2), length.out = 6)
  p1 <- levelplot(dat4[[n4]]$pref, main = list(label = 'True "Data Scenario (DS) 04"', side = 1, line = 0.5, cex = 1.5),
                  col.regions = pal, 
                  at = cuts, 
                  breaks = cuts, 
                  margin = F,
                  scales = list(x = list(at = seq(0, 10, 2), cex = 1.5),
                                y = list(at = seq(0, 10, 2), cex = 1.5)),
                  colorkey = list(at = cuts, 
                                  labels = list(at = labs, cex = 1.5,
                                                labels = as.character(format(labs, nsmall = 3)))),
                  par.settings = font_list)
  
  p2 <- levelplot(est_pref4, main = list(label = 'Estimated (Hor. Uni. Tri) DS 04', side = 1, line = 0.5, cex = 1.5),
                  col.regions = pal, 
                  at = cuts, 
                  breaks = cuts, 
                  margin = F,
                  scales = list(x = list(at = seq(0, 10, 2), cex = 1.5),
                                y = list(at = seq(0, 10, 2), cex = 1.5)),
                  colorkey = list(at = cuts, 
                                  labels = list(at = labs, cex = 1.5,
                                                labels = as.character(format(labs, nsmall = 3)))),
                  par.settings = font_list)
  
  p3 <- levelplot(est_pref6, main = list(label = 'Estimated (RB Wendland) DS 04', side = 1, line = 0.5, cex = 1.5),
                  col.regions = pal, 
                  at = cuts, 
                  breaks = cuts, 
                  margin = F,
                  scales = list(x = list(at = seq(0, 10, 2), cex = 1.5),
                                y = list(at = seq(0, 10, 2), cex = 1.5)),
                  colorkey = list(at = cuts, 
                                  labels = list(at = labs, cex = 1.5,
                                                labels = as.character(format(labs, nsmall = 3)))),
                  par.settings = font_list)
  
  ##################################################
  ##################################################
  
  pal <- jet.col(n = 100, alpha = 1)
  cuts <- seq(round(min(min5, min7) - 0.01, 2), round(max(max5, max7) + 0.00, 2), length.out = 100)
  labs <- seq(round(min(min5, min7) - 0.01, 2), round(max(max5, max7) + 0.00, 2), length.out = 6)
  p4 <- levelplot(dat5[[n5]]$pref, main = list(label = 'True "Data Scenario (DS) 05"', side = 1, line = 0.5, cex = 1.5),
                  col.regions = pal, 
                  at = cuts, 
                  breaks = cuts, 
                  margin = F,
                  scales = list(x = list(at = seq(0, 10, 2), cex = 1.5),
                                y = list(at = seq(0, 10, 2), cex = 1.5)),
                  colorkey = list(at = cuts, 
                                  labels = list(at = labs, cex = 1.5,
                                                labels = as.character(format(labs, nsmall = 3)))),
                  par.settings = font_list)
  
  p5 <- levelplot(est_pref5, main = list(label = 'Estimated (RB Wendland) DS 05', side = 1, line = 0.5, cex = 1.5),
                  col.regions = pal, 
                  at = cuts, 
                  breaks = cuts, 
                  margin = F,
                  scales = list(x = list(at = seq(0, 10, 2), cex = 1.5),
                                y = list(at = seq(0, 10, 2), cex = 1.5)),
                  colorkey = list(at = cuts, 
                                  labels = list(at = labs, cex = 1.5,
                                                labels = as.character(format(labs, nsmall = 3)))),
                  par.settings = font_list)
  
  p6 <- levelplot(est_pref7, main = list(label = 'Estimated (Hor. Uni. Tri) DS 05', side = 1, line = 0.5, cex = 1.5),
                  col.regions = pal, 
                  at = cuts, 
                  breaks = cuts, 
                  margin = F,
                  scales = list(x = list(at = seq(0, 10, 2), cex = 1.5),
                                y = list(at = seq(0, 10, 2), cex = 1.5)),
                  colorkey = list(at = cuts, 
                                  labels = list(at = labs, cex = 1.5,
                                                labels = as.character(format(labs, nsmall = 3)))),
                  par.settings = font_list)
  # plot(p1)
  grid.arrange(p1, p2, p3, p4, p5, p6, ncol = 3)
  dev.off()
}


####################################################################################################
# Figure 06
####################################################################################################

scenarios <- c(5, 5, 8)
order_scenarios <- c(1, 2, 3)
fitted_model <- "var_DP_PS"
smoothing_kernel <- c("Wendland", "Wendland", "Wendland")
partition_const  <- c("non_constant", "non_constant", "non_constant")
n_center_points  <- c("few", "many", "manymany")
sample_sizes <- c(100, 200, 500, 1000)

errors_df <- as.data.frame(matrix(data = 0, nrow = 1200, ncol = 3))
errors <- list()
cc <- 1
for (sample_size in sample_sizes) {
  errors[[as.character(sample_size)]] <- list()
  for (i in 1:3) {
    errors[[as.character(sample_size)]][[i]] <- readRDS(file = paste("FITTED_MODELS/", sprintf("%02d", scenarios[i]), "/ERRORS/UNNORMALIZED_", fitted_model, "-", smoothing_kernel[i], "-", partition_const[i], "-", n_center_points[i], "-", sample_size, ".rds", sep = ""))   
    for (j in 1:100) {
      errors_df[cc, ] <- c(sample_size, order_scenarios[i], errors[[as.character(sample_size)]][[i]][[j]]$error_preferentiality)
      cc <- cc + 1
    }
  }
}
colnames(errors_df) <- c("sample_size", "data_scenario", "error")

jpeg(filename = "boxplot-number-basis.jpeg", width = 3800, height = 1100, res = 300)
at <- c(1:4, 5.5:8.5, 10:13)
par(mfrow = c(1, 1), family = 'LM Roman 10', mar = c(3.5, 4, 0, 0) + 0.1)
boxplot(error ~ sample_size + data_scenario, data = errors_df, xlim = c(1.0, 13), ylim = c(0, 2.05), axes = T, xaxt = "n", yaxt = "n", at = at, xlab = "", ylab = "Squared Error (SE)", cex.lab = 1.5)
axis(1, at = at, labels = rep(c("100", "200", "500", "1,000"), 3), cex.axis = 1.5)
axis(2, at = 0:5, labels = 0:5, cex.axis = 1.5)
text(2.5,  -0.55, "K = 5",  xpd = NA, cex = 1.5)
text(7.0,  -0.55, "K = 9",  xpd = NA, cex = 1.5) 
text(11.5, -0.55, "K = 13", xpd = NA, cex = 1.5)
dev.off()

####################################################################################################
# Figure 07
####################################################################################################

source("functions_result_processing.R")

####################
xlim <- c(0, 10)
ylim <- c(0, 10)
by <- 0.05
####################

dat5 <- readRDS(file = "DATA/05/1000.rds")

res5.1 <- readRDS(file = "FITTED_MODELS/05/var_DP_PS-Wendland-non_constant-few-1000.rds")
res5.2 <- readRDS(file = "FITTED_MODELS/05/var_DP_PS-Wendland-non_constant-many-1000.rds")
res5.3 <- readRDS(file = "FITTED_MODELS/08/var_DP_PS-Wendland-non_constant-manymany-1000.rds")

oth5.1 <- readRDS(file = "FITTED_MODELS/05/OTHERS/var_DP_PS-Wendland-non_constant-few-1000.rds")
oth5.2 <- readRDS(file = "FITTED_MODELS/05/OTHERS/var_DP_PS-Wendland-non_constant-many-1000.rds")
oth5.3 <- readRDS(file = "FITTED_MODELS/08/OTHERS/var_DP_PS-Wendland-non_constant-manymany-1000.rds")

##################################################
##################################################

font_list <- list(axis.text = list(fontfamily = "LM Roman 10"),
                  par.xlab.text = list(fontfamily = "LM Roman 10"),
                  par.ylab.text = list(fontfamily = "LM Roman 10"),
                  par.main.text = list(fontfamily = "LM Roman 10"),
                  par.sub.text  = list(fontfamily = "LM Roman 10"))

##################################################
##################################################

nn <- 9

n5.1 <- nn

est_pref5.1 <- compute_preferentiatility(fit = res5.1[[n5.1]]$fit, coord_pred = oth5.1[[n5.1]]$coord_pred, bfs_pred = oth5.1[[n5.1]]$bfs_pred, xlim = xlim, by = by)

min5.1 <- min(values(dat5[[n5.1]]$pref), values(est_pref5.1))
max5.1 <- max(values(dat5[[n5.1]]$pref), values(est_pref5.1))

n5.2 <- nn

est_pref5.2 <- compute_preferentiatility(fit = res5.2[[n5.2]]$fit, coord_pred = oth5.2[[n5.2]]$coord_pred, bfs_pred = oth5.2[[n5.2]]$bfs_pred, xlim = xlim, by = by)

min5.2 <- min(values(dat5[[n5.1]]$pref), values(est_pref5.2))
max5.2 <- max(values(dat5[[n5.1]]$pref), values(est_pref5.2))

n5.3 <- nn

est_pref5.3 <- compute_preferentiatility(fit = res5.3[[n5.3]]$fit, coord_pred = oth5.3[[n5.3]]$coord_pred, bfs_pred = oth5.3[[n5.3]]$bfs_pred, xlim = xlim, by = by)

min5.3 <- min(values(dat5[[n5.1]]$pref), values(est_pref5.3))
max5.3 <- max(values(dat5[[n5.1]]$pref), values(est_pref5.3))

mn <- min(min5.1, min5.2, min5.3)
mx <- max(max5.1, max5.2, max5.3)

##################################################
##################################################

if (TRUE) {
  jpeg(filename = "estimated-pref-number-basis.jpeg", width = 4000, height = 2200, res = 300)
  par(mfrow = c(1, 4), family = 'LM Roman 10', mar = c(0, 0, 0, 0) + 0.05)
  
  pal <- jet.col(n = 100, alpha = 1)
  cuts <- seq(round(mn - 0.01, 2), round(mx + 0.00, 2), length.out = 100)
  labs <- seq(round(mn - 0.01, 2), round(mx + 0.00, 2), length.out = 6)
  p1 <- levelplot(dat5[[n5.1]]$pref, main = list(label = 'True Data Scenario (DS) 05', side = 1, line = 0.5, cex = 1.5),
                  col.regions = pal, 
                  at = cuts, 
                  breaks = cuts, 
                  margin = F,
                  scales = list(x = list(at = seq(0, 10, 2), cex = 1.5),
                                y = list(at = seq(0, 10, 2), cex = 1.5)),
                  colorkey = list(at = cuts, 
                                  labels = list(at = labs, cex = 1.5,
                                                labels = as.character(format(labs, nsmall = 3)))),
                  par.settings = font_list)
  
  p2 <- levelplot(est_pref5.1, main = list(label = 'Estimated DS 05 (K = 5)', side = 1, line = 0.5, cex = 1.5),
                  col.regions = pal, 
                  at = cuts, 
                  breaks = cuts, 
                  margin = F,
                  scales = list(x = list(at = seq(0, 10, 2), cex = 1.5),
                                y = list(at = seq(0, 10, 2), cex = 1.5)),
                  colorkey = list(at = cuts, 
                                  labels = list(at = labs, cex = 1.5,
                                                labels = as.character(format(labs, nsmall = 3)))),
                  par.settings = font_list)
  
  p3 <- levelplot(est_pref5.2, main = list(label = 'Estimated DS 05 (K = 9)', side = 1, line = 0.5, cex = 1.5),
                  col.regions = pal, 
                  at = cuts, 
                  breaks = cuts, 
                  margin = F,
                  scales = list(x = list(at = seq(0, 10, 2), cex = 1.5),
                                y = list(at = seq(0, 10, 2), cex = 1.5)),
                  colorkey = list(at = cuts, 
                                  labels = list(at = labs, cex = 1.5,
                                                labels = as.character(format(labs, nsmall = 3)))),
                  par.settings = font_list)
  

  p4 <- levelplot(est_pref5.3, main = list(label = 'Estimated DS 05 (K = 13)', side = 1, line = 0.5, cex = 1.5),
                  col.regions = pal, 
                  at = cuts, 
                  breaks = cuts, 
                  margin = F,
                  scales = list(x = list(at = seq(0, 10, 2), cex = 1.5),
                                y = list(at = seq(0, 10, 2), cex = 1.5)),
                  colorkey = list(at = cuts, 
                                  labels = list(at = labs, cex = 1.5,
                                                labels = as.character(format(labs, nsmall = 3)))),
                  par.settings = font_list)
  

  # plot(p1)
  grid.arrange(p1, p2, p3, p4, layout_matrix = matrix(c(1, NA, NA, 4:6), ncol = 3, byrow = TRUE))
  dev.off()
}

####################################################################################################
# Figure 08
####################################################################################################

USA <- readRDS(file = "APPLICATION/USA0.rds")
USA  <- st_transform(x = USA, crs = "+init=epsg:6345 +units=km +no_defs")

project_center_pts <- function (center_pts) {
  colnames(center_pts) <- c("long", "latt")
  center_pts <- st_as_sf(x = as_tibble(as.data.frame(center_pts)), coords = c("long", "latt"), crs = "WGS84")
  center_pts <- st_transform(x = center_pts, crs = "+init=epsg:6345 +units=km +no_defs")
  center_pts <- st_coordinates(center_pts)
  colnames(center_pts) <- NULL
  center_pts
}

center_pts_5  <- project_center_pts(rbind(c(-118, 43.5), c(-95, 40), c(-80, 38), c(-105, 36), c(-90, 45)))
center_pts_10 <- project_center_pts(rbind(c(-118, 35), c(-118, 44), c(-105, 45), c(-98, 34), c(-90, 45), c(-85, 37), c(-75, 42), c(-95, 38.5), c(-108, 38), c(-83, 32.5)))
center_pts_15 <- project_center_pts(rbind(c(-118, 35), c(-118, 44), c(-105, 45), c(-98, 34), c(-90, 45), c(-85, 37), c(-75, 42), c(-95, 38.5), c(-112, 38), c(-83, 32.5), c(-105, 32.5), c(-109, 42.5), c(-96, 46.5), c(-82, 41), c(-96, 43)))

jpeg(filename = "center_pts_wendland.jpeg", width = 4600, height = 1100, res = 300)
par(mfrow = c(1, 3), family = 'LM Roman 10', mar = c(0, 0, 1.5, 0) + 0.05)
plot(USA$geometry, main = "(1)", cex.main = 3)
points(center_pts_5,  pch = 4, col = "black", cex = 2.5)
plot(USA$geometry, main = "(2)", cex.main = 3)
points(center_pts_10, pch = 4, col = "black", cex = 2.5)
plot(USA$geometry, main = "(3)", cex.main = 3)
points(center_pts_15, pch = 4, col = "black", cex = 2.5)
par(mfrow = c(1, 1))
dev.off()

####################################################################################################
# Figure 09
####################################################################################################

source("header.R")
library("tidyverse")
library("plot3D")
library("patchwork")

scenarios <- 1:9

data <- readRDS(file = "APPLICATION/data.rds")
USA  <- readRDS(file = "APPLICATION/USA.rds")

USA <- USA[!(USA$shapeName %in% c("Alaska", "American Samoa","Commonwealth of the Northern Mariana Islands", "Guam", "Hawaii", "Puerto Rico", "United States Virgin Islands")), ]
USA <- st_transform(x = USA, crs = "+init=epsg:6345 +units=km +no_defs")

result_obj <- list()
for (scenario in scenarios) { result_obj[[scenario]] <- readRDS(file = paste("APPLICATION/FITTED_MODELS/scenario", scenario, ".rds", sep = "")) }

##################################################
##################################################

result <- list(); bfs_pred <- list(); fs_formula <- list(); coord_pred <- list()

for (scenario in scenarios) {
  result[[scenario]]     <- result_obj[[scenario]]$result
  bfs_pred[[scenario]]   <- result_obj[[scenario]]$bfs_pred
  fs_formula[[scenario]] <- result_obj[[scenario]]$fs_formula
  coord_pred <- result_obj[[scenario]]$coord_pred
  n_coord_pr <- result_obj[[scenario]]$n_coord_pr
}

sce_1 <- 1
sce_2 <- 2
sce_3 <- 5

idx_lat_1 <- inla.stack.index(fs_formula[[sce_1]]$full_stack, tag = "pred_y")$data[1:n_coord_pr]
fitted_latent_1 <- cbind(coord_pred[1:n_coord_pr, ], result[[sce_1]]$fit$summary.fitted.values[idx_lat_1, c("mean")])
colnames(fitted_latent_1) <- c("x", "y", "estimated")

idx_lat_2 <- inla.stack.index(fs_formula[[sce_2]]$full_stack, tag = "pred_y")$data[1:n_coord_pr]
fitted_latent_2 <- cbind(coord_pred[1:n_coord_pr, ], result[[sce_2]]$fit$summary.fitted.values[idx_lat_2, c("mean")])
colnames(fitted_latent_2) <- c("x", "y", "estimated")

idx_lat_3 <- inla.stack.index(fs_formula[[sce_3]]$full_stack, tag = "pred_y")$data[1:n_coord_pr]
fitted_latent_3 <- cbind(coord_pred[1:n_coord_pr, ], result[[sce_3]]$fit$summary.fitted.values[idx_lat_3, c("mean")])
colnames(fitted_latent_3) <- c("x", "y", "estimated")

fitted_latent_diff_1 <- fitted_latent_1
fitted_latent_diff_2 <- fitted_latent_1
fitted_latent_diff_1[, "estimated"] <- fitted_latent_3[, "estimated"] - fitted_latent_1[, "estimated"]
fitted_latent_diff_2[, "estimated"] <- fitted_latent_3[, "estimated"] - fitted_latent_2[, "estimated"]

# Create a gridded spatial object from "process"

coordinates(fitted_latent_diff_1) <- ~ x + y
gridded(fitted_latent_diff_1) <- TRUE
fitted_latent_diff_1 <- raster(fitted_latent_diff_1)
crs(fitted_latent_diff_1) = "+init=epsg:6345 +units=km +no_defs"

fitted_latent_diff_1    <- as(fitted_latent_diff_1, "SpatialPixelsDataFrame")
fitted_latent_diff_df_1 <- as.data.frame(fitted_latent_diff_1)
colnames(fitted_latent_diff_df_1) <- c("estimated", "x", "y")

#####

coordinates(fitted_latent_diff_2) <- ~ x + y
gridded(fitted_latent_diff_2) <- TRUE
fitted_latent_diff_2 <- raster(fitted_latent_diff_2)
crs(fitted_latent_diff_2) = "+init=epsg:6345 +units=km +no_defs"

fitted_latent_diff_2    <- as(fitted_latent_diff_2, "SpatialPixelsDataFrame")
fitted_latent_diff_df_2 <- as.data.frame(fitted_latent_diff_2)
colnames(fitted_latent_diff_df_2) <- c("estimated", "x", "y")

#####

pal <- colorRampPalette(colors = c("blue", "white", "red"))(n = 101)
# pal <- jet.col(n = 100, alpha = 0.9)
mn_1 <- round(min(fitted_latent_diff_df_1$estimated) - 0.01, 2)
mx_1 <- round(max(fitted_latent_diff_df_1$estimated) + 0.01, 2)
mn_2 <- round(min(fitted_latent_diff_df_2$estimated) - 0.01, 2)
mx_2 <- round(max(fitted_latent_diff_df_2$estimated) + 0.01, 2)
mn <- min(mn_1, mn_2)
mx <- max(mx_1, mx_2)
mn <- max(c(abs(mn), abs(mx))) * -1
mx <- mn * -1
labs <- seq(mn, mx, length.out = 5)
# labs <- seq(round(min(fitted_latent_diff_df$estimated) - 0.01, 2), round(max(fitted_latent_diff_df$estimated) + 0.01, 2), length.out = 6)

p_diff_1 <- ggplot() +
  geom_tile(data = fitted_latent_diff_df_1, mapping = aes(x = x, y = y, fill = estimated)) + 
  geom_sf(data = USA, color = "black", fill = NA, lwd = 0.5) +
  scale_fill_gradientn(name = expression(paste("Diff. ", PM[2.5], sep = "")), colors = pal, breaks = labs, labels = as.character(format(labs, nsmall = 3)), limits = c(labs[1], tail(labs, 1))) +
  labs(x = "Longitude", y = "Latitude", title = "Diff. between \"Sp. Var. PS\" and \"Non-PS\"") + 
  theme_bw() +
  theme(text = element_text(size = 24, family = "LM Roman 10"), 
        legend.key.height = unit(2.57, "cm"),
        legend.key.width  = unit(0.75, "cm"),
        plot.margin = margin(0, 1, 0, 0, "cm"))

p_diff_2 <- ggplot() +
  geom_tile(data = fitted_latent_diff_df_2, mapping = aes(x = x, y = y, fill = estimated)) + 
  geom_sf(data = USA, color = "black", fill = NA, lwd = 0.5) +
  scale_fill_gradientn(name = expression(paste("Diff. ", PM[2.5], sep = "")), colors = pal, breaks = labs, labels = as.character(format(labs, nsmall = 3)), limits = c(labs[1], tail(labs, 1))) +
  labs(x = "Longitude", y = "Latitude", title = "Diff. between \"Sp. Var. PS\" and \"Const. PS\"") + 
  theme_bw() +
  theme(text = element_text(size = 24, family = "LM Roman 10"), 
        legend.key.height = unit(2.57, "cm"),
        legend.key.width  = unit(0.75, "cm"),
        plot.margin = margin(0, 0, 0, 0, "cm"))

p_diff <- p_diff_1 + p_diff_2 & theme(legend.position = "right", legend.box.margin = margin(0, 0, 0, 0.5, "cm"))
p_diff <- p_diff + plot_layout(guides = "collect")

ggsave(filename = paste("APPLICATION/IMAGES/diff_latent.jpeg", sep = ""), plot = p_diff, width = 6010, height = 2000, units = c("px"), dpi = 300, bg = "white")





####################################################################################################
# Figure 10
####################################################################################################

x <- seq(0, 1, length.out = 1000)  
y <- x  

constant <- function (x, y, k = 0.1, ...) {
  rep(k, length(x))
} 

unidirectional_triangular <- function (x, y, center_x, center_y, i = i, NS = FALSE, ...) {
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

partition <- function (x, y, center_x, center_y, ...) {
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

Wendland <- function (x, y, center_x, center_y, bandwidth = 0.5, ...) {
  pts <- cbind((x - center_x), (y - center_y))
  d <- apply(X = pts, MARGIN = 1, FUN = Euclidean_norm) / bandwidth
  k <- (1 - d) ^ 6 * (35 * d ^ 2 + 18 * d + 3) / 3 * as.numeric(d <= 1)
  k
}

if (TRUE) {
  jpeg(filename = "basis-functions.jpeg", width = 3600, height = 900, res = 300)
  par(mfrow = c(1, 4), family = 'LM Roman 10', mar = c(0, 0, 2, 0) + 0.1)
  # 1
  z1 <- outer(x, y, constant, k = 0.5); z1[1, 1] <- 0.5 + 1e-6
  
  persp3D(x, y, z1, xlim = c(0, 1), ylim = c(0, 1), zlim = c(0, 1), theta = 45, phi = 30, axes = F, xlab = "x1", ylab = "x2", zlab = "\U03b3(x)", alpha = 0.5, colkey = F, clim = c(0, 1), main = "(1)")
  
  # 2
  z1 <- outer(x, y, partition, center_x = c(0.0, 1.0), center_y = c(0.0, 0.499)) * 0.25
  z2 <- outer(x, y, partition, center_x = c(0.0, 1.0), center_y = c(0.5, 1.000)) * 0.75
  persp3D(x, y, z1, xlim = c(0, 1), ylim = c(0, 1), zlim = c(0, 1), theta = 45, phi = 30, axes = F, xlab = "x1", ylab = "x2", zlab = "\U03b3(x)", alpha = 0.5, colkey = F, clim = c(0, 1), main = "(2)")
  persp3D(x, y, z2, xlim = c(0, 1), ylim = c(0, 1), zlim = c(0, 1), theta = 45, phi = 30, axes = F, xlab = "x1", ylab = "x2", zlab = "\U03b3(x)", alpha = 0.5, colkey = F, clim = c(0, 1), add = T)
  
  # 3
  z1 <- outer(x, y, unidirectional_triangular, rbind(c(0, 0.50), c(0.50, 0.50), c(1, 0.50)), i = 1)
  z2 <- outer(x, y, unidirectional_triangular, rbind(c(0, 0.50), c(0.50, 0.50), c(1, 0.50)), i = 2)
  z3 <- outer(x, y, unidirectional_triangular, rbind(c(0, 0.50), c(0.50, 0.50), c(1, 0.50)), i = 3)
  persp3D(x, y, z1, xlim = c(0, 1), ylim = c(0, 1), zlim = c(0, 1), theta = 45, phi = 30, axes = F, xlab = "x1", ylab = "x2", zlab = "\U03b3(x)", alpha = 0.5, colkey = F, clim = c(0, 1), main = "(3)")
  persp3D(x, y, z2, xlim = c(0, 1), ylim = c(0, 1), zlim = c(0, 1), theta = 45, phi = 30, axes = F, xlab = "x1", ylab = "x2", zlab = "\U03b3(x)", alpha = 0.5, colkey = F, clim = c(0, 1), add = T)
  persp3D(x, y, z3, xlim = c(0, 1), ylim = c(0, 1), zlim = c(0, 1), theta = 45, phi = 30, axes = F, xlab = "x1", ylab = "x2", zlab = "\U03b3(x)", alpha = 0.5, colkey = F, clim = c(0, 1), add = T)
  
  # 4
  z1 <- outer(x, y, Wendland, center_x = 0.25, center_y = 0.25)
  z2 <- outer(x, y, Wendland, center_x = 0.75, center_y = 0.25)
  z3 <- outer(x, y, Wendland, center_x = 0.25, center_y = 0.75)
  z4 <- outer(x, y, Wendland, center_x = 0.75, center_y = 0.75)
  
  persp3D(x, y, z1, xlim = c(0, 1), ylim = c(0, 1), zlim = c(0, 1), theta = 45, phi = 30, axes = F, xlab = "x1", ylab = "x2", zlab = "\U03b3(x)", alpha = 0.5, colkey = F, clim = c(0, 1), main = "(4)")
  persp3D(x, y, z2, xlim = c(0, 1), ylim = c(0, 1), zlim = c(0, 1), theta = 45, phi = 30, axes = F, xlab = "x1", ylab = "x2", zlab = "\U03b3(x)", alpha = 0.5, colkey = F, clim = c(0, 1), add = T)
  persp3D(x, y, z3, xlim = c(0, 1), ylim = c(0, 1), zlim = c(0, 1), theta = 45, phi = 30, axes = F, xlab = "x1", ylab = "x2", zlab = "\U03b3(x)", alpha = 0.5, colkey = F, clim = c(0, 1), add = T)
  persp3D(x, y, z4, xlim = c(0, 1), ylim = c(0, 1), zlim = c(0, 1), theta = 45, phi = 30, axes = F, xlab = "x1", ylab = "x2", zlab = "\U03b3(x)", alpha = 0.5, colkey = F, clim = c(0, 1), add = T)
  
  dev.off()
}
