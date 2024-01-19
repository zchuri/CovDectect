#!/usr/bin/env Rscript

# Clear workspace
rm(list = ls())
# Save original parameters
op <- par()

# Load libraries
if(!require(devtools)) install.packages("devtools"); library(devtools)
if(!require(dplyr)) install.packages("dplyr"); library(dplyr)
if(!require(RecordTest)) install.packages("RecordTest"); library(RecordTest)

# Data directory
geo_df <- read.csv("geo_info.csv")
tmax_df <- readRDS("Y365.rds")
tmax_dim <- dim(tmax_df)

# Draw one station with p.plot
p.plot(tmax_df[,,9], record = c(1,0,0,0), conf.int = F, smooth = F)

# Binarize all
bin_df <- array(data = 0, dim = tmax_dim)
for(nn in 1:tmax_dim[3]){
  # Zaragoza case
  tmax_zgz <- tmax_df[,,nn]
  # Find records
  rcd_zgz <- lapply(1:ncol(tmax_zgz),
                    function(x){
                      records(X = tmax_zgz[,x], plot = F)[[1]][,1]
                    })
  
  # Binarize records
  bin_zgz <- array(data = 0, dim = dim(tmax_zgz))
  for(dd in 1:length(rcd_zgz)) bin_zgz[rcd_zgz[[dd]],dd] <- 1
  
  bin_df[,,nn] <- bin_zgz
}

#############################################################
# GLM with only temporal components
#############################################################

# Extract for the good observatories
qc_idx <- which(geo_df$QC==1)
bin_df <- bin_df[,,qc_idx]
bin_dim <- dim(bin_df)
geo_df <- geo_df[qc_idx,]

# Extract Zaragoza observatory (9)
bin_zgz <- data.frame(y = c(t(bin_df[,,9])),
                      t = rep(x = 1:bin_dim[1], each = bin_dim[2]),
                      day = rep(x = 1:bin_dim[2], times = bin_dim[1]))

# Create daily sin and cos
bin_zgz$sin.l <- sin(2*pi*bin_zgz$day/bin_dim[2])
bin_zgz$cos.l <- cos(2*pi*bin_zgz$day/bin_dim[2])

# Construct GLM
frm_zgz <- as.formula(paste0("y ~  poly(t,2) + sin.l + cos.l + 
                                   poly(t,2):sin.l + poly(t,2):cos.l +
                                   lag(y,1,0) + lag(y,2,0) +
                                   lag(y,1,0):lag(y,2,0) +
                                   t:lag(y,1,0) + t:lag(y,2,0) +
                                   t:lag(y,1,0):lag(y,2,0)"))
fit_zgz <- glm(formula = frm_zgz, data = bin_zgz, family = binomial)
coef_zgz <- summary(fit_zgz)$coef

#############################################################
# Extract GLM coefficients for all observatories
#############################################################

# Initial variables
coef_n <- nrow(coef_zgz)
coef_lab <- rownames(coef_zgz) 
coef_mat <- matrix(data = 0, nrow = nrow(geo_df), ncol = coef_n)

# Extract coefficients
for(oo in 1:nrow(geo_df)){
  
  # Update observatory records
  bin_zgz$y <- c(t(bin_df[,,oo]))
  # Apply GLM
  fit_zgz <- glm(formula = frm_zgz, data = bin_zgz, family = binomial)
  # Save coefficients
  coef_zgz <- summary(fit_zgz)$coef
  coef_mat[oo,] <- coef_zgz[,1]
  
}

# Add labels
rownames(coef_mat) <- as.character(geo_df$ID)
colnames(coef_mat) <- coef_lab


#############################################################
# Find covariate effects with coefficients
#############################################################

# Explore effects
geo_df$b0 <- coef_mat[,1]
geo_df$logdist <- log(geo_df$CoastDist)

# Correlation
cor(cbind(geo_df$b0, geo_df$logdist, geo_df$HGHT))

# Plot log(dist) effects
library(ggplot2)

# Output directory
out_dir <- file.path(getwd(), "Results")
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = T)

# Plot
for(oo in 1:coef_n){
  
  # Update beta coefficient
  geo_df$b0 <- coef_mat[,oo]
  ttl <- paste0("log(dist(s)) vs glm ",coef_lab[oo])
  cor_text <- paste0("rho=",round(cor(geo_df$logdist,geo_df$b0),2))
  g1 <- ggplot(data = geo_df, aes(x=logdist, y=b0, label = abb)) +
    #geom_point() +
    geom_smooth(method=lm , se=FALSE) +
    geom_label(aes(fill = factor(Zona)), #colour = "white",
               fontface = "bold", size = 2) +
    theme_bw() +
    labs(title = ttl,
         subtitle = cor_text,
         caption = paste0("(",nrow(geo_df)," stations)"))+
    ylab(coef_lab[1]) +
    xlab("log(dist(s))") +
    scale_fill_brewer(palette = "Set1")
  show(g1)
  # Save plot
  out_name <- paste0("log-dist_vs_beta",oo,".png")
  ggsave(filename = out_name, plot = g1, device = "png", path = out_dir,
         scale = 1, width = 7, height = 5)

}# for oo

