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
# Compute the stacked data.frame
#############################################################

# Stacked data.frame object
df_ln <- prod(bin_dim) # Series length
st_ln <- df_ln/bin_dim[3]
stk_df <- data.frame(y = rep(0,df_ln),
                     st = rep("",df_ln),
                     lat = rep(0,df_ln),
                     lon = rep(0,df_ln),
                     elev = rep(0,df_ln),
                     dist = rep(0,df_ln),
                     day = rep(0,df_ln),
                     sin.d = rep(0,df_ln),
                     cos.d = rep(0,df_ln),
                     t = rep(0,df_ln),
                     lag1 = rep(0,df_ln),
                     lag2 = rep(0,df_ln))

# Begin to stack
for(ss in 1:bin_dim[3]){
  
  cat(paste0(ss, "."))
  # Stacked indices
  idx <- (1+(ss-1)*st_ln):(ss*st_ln)
  # Extract records
  stk_df$y[idx] <- c(t(bin_df[,,ss]))
  # Station name
  stk_df$st[idx] <- geo_df$abb[ss]
  # Latitude
  stk_df$lat[idx] <- geo_df$LAT[ss]
  # Longitude
  stk_df$lon[idx] <- geo_df$LON[ss]
  # Elevation
  stk_df$elev[idx] <- geo_df$HGHT[ss]
  # Log dist to coast
  stk_df$dist[idx] <- log(geo_df$CoastDist[ss])
  # Day
  stk_df$day[idx] <- rep(x = 1:bin_dim[2], times = bin_dim[1])
  # Year
  stk_df$t[idx] <- rep(x = 1:bin_dim[1], each = bin_dim[2])
  # Lags
  stk_df$lag1[idx] <- lag(stk_df$y[idx],1,0)
  stk_df$lag2[idx] <- lag(stk_df$lag1[idx],1,0)
}
cat("\n")

# Compute the rest of the predictors
stk_df$sin.d <- sin(2*pi*stk_df$day/bin_dim[2])
stk_df$cos.d <- cos(2*pi*stk_df$day/bin_dim[2])
stk_df$trend <- log(stk_df$t-1)
stk_df$st <- factor(stk_df$st)

# Remove infinite values
stk_df <- stk_df[which(!is.infinite(stk_df$trend)),]

#############################################################
# Model fitting - Main table
#############################################################

# Fit temporal single term model
fit_null <- glm(formula = y ~ offset(I(-trend)) - 1,
                family = binomial,
                data = stk_df)
fit_null.L <- glm(formula = y ~ trend,
                  family = binomial,
                  data = stk_df)
fit_null.Q <- glm(formula = y ~ poly(trend,2),
                  family = binomial,
                  data = stk_df)
fit_Q.lags <- glm(formula = y ~ poly(trend,2)*lag1*lag2,
                  family = binomial,
                  data = stk_df)
fit_Q.lags.harm <- glm(formula = y ~ poly(trend,2)*(lag1*lag2+sin.d+cos.d),
                       family = binomial,
                       data = stk_df)
frm_paper <- as.formula("y ~ poly(trend,2) * (dist + sin.d + cos.d) + (trend + dist)*lag1*lag2")
fit_paper <- glm(formula = frm_paper,
                 family = binomial,
                 data = stk_df)
# AICs
AIC(fit_null, fit_null.L, fit_null.Q, fit_Q.lags, fit_Q.lags.harm, fit_paper)

#############################################################
# Model fitting - Forward
#############################################################

# Remove current models to save working memory
rm(list = ls()[grep("fit_",ls())])
gc()

# Fit temporal single term model
fit_null <- glm(formula = y ~ offset(I(-trend)) - 1,
                family = binomial,
                data = stk_df)
fit_null.L <- glm(formula = y ~ trend,
                  family = binomial,
                  data = stk_df)
fit_null.Q <- glm(formula = y ~ poly(trend,2),
                  family = binomial,
                  data = stk_df)

# Fit temporal single term model
frm_tmp_single <- as.formula("y ~ trend + sin.d + cos.d + lag1*lag2")
fit_tmp_single <- glm(formula = frm_tmp_single,
                      family = binomial,
                      data = stk_df)

# Fit temporal single term model
frm_tmp_trend2 <- as.formula("y ~ poly(trend,2) + sin.d + cos.d + lag1*lag2")
fit_tmp_trend2 <- glm(formula = frm_tmp_trend2,
                      family = binomial,
                      data = stk_df)

# Fit temporal single term model
frm_tmp_trend3 <- as.formula("y ~ poly(trend,3) + sin.d + cos.d + lag1*lag2")
fit_tmp_trend3 <- glm(formula = frm_tmp_trend3,
                      family = binomial,
                      data = stk_df)
AIC(fit_tmp_single, fit_tmp_trend2, fit_tmp_trend3)

# Fit temporal single term model plus log-dist
fit_trend2_dist <- update(fit_tmp_trend2, .~.+dist)

# Fit temporal interaction model
frm_tmp_2nd <- as.formula("y ~ poly(trend,2) * (sin.d + cos.d) + trend*lag1*lag2")
fit_tmp_2nd <- glm(formula = frm_tmp_2nd,
                   family = binomial,
                   data = stk_df)

# Fit temporal interaction model plus log-dist
fit_2nd_t.dist <- update(fit_tmp_2nd, .~.+dist*poly(trend,2))
fit_2nd_lag.dist <- update(fit_tmp_2nd, .~.+dist*lag1*lag2)

# Current paper model
frm_paper <- as.formula("y ~ poly(trend,2) * (dist + sin.d + cos.d) + (trend + dist)*lag1*lag2")
fit_paper <- glm(formula = frm_paper,
                 family = binomial,
                 data = stk_df)

# Compare AIC's
AIC(fit_null,
    fit_null.L,
    fit_null.Q,
    fit_tmp_single,
    fit_tmp_trend2,
    fit_tmp_trend3,
    fit_trend2_dist,
    fit_tmp_2nd,
    fit_2nd_t.dist,
    fit_2nd_lag.dist,
    fit_paper)

###############################################################################
# Model fitting intermediate models
###############################################################################

# Remove current models to save working memory
rm(list = ls()[grep("fit_",ls())])
gc()

# Fit temporal single term model
fit_null.Q <- glm(formula = y ~ poly(trend,2),
                  family = binomial,
                  data = stk_df)

# Fit temporal interaction model plus log-dist
fit_Q_lag1 <- update(fit_null.Q, .~.+lag1*trend)

# Compare AIC's
AIC(fit_null.Q, fit_Q_lag1)

###############################################################################
# Model fitting adding extra interactions
###############################################################################

# Remove current models to save working memory
rm(list = ls()[grep("fit_",ls())])
gc()

# Current paper model
frm_paper <- as.formula("y ~ poly(trend,2) * (dist + sin.d + cos.d) + (trend + dist)*lag1*lag2")
fit_paper <- glm(formula = frm_paper,
                 family = binomial,
                 data = stk_df)

# Adding (sin+cos)*lag1*lag2
fit_add_sincos.lag12 <- update(fit_paper, .~.+(sin.d + cos.d)*lag1*lag2)

# Adding (sin+cos)*log(dist)
fit_add_sincos.dist <- update(fit_paper, .~.+(sin.d + cos.d)*dist)

# Compare AIC's
AIC(fit_paper, fit_add_sincos.dist, fit_add_sincos.lag12)

###############################################################################
# Model fitting adding extra interactions
###############################################################################

# Remove current models to save working memory
rm(list = ls()[grep("fit_",ls())])
gc()

# Current paper model
frm_paper <- as.formula("y ~ poly(trend,2) * (dist + sin.d + cos.d) + (trend + dist)*lag1*lag2")
fit_paper <- glm(formula = frm_paper,
                 family = binomial,
                 data = stk_df)

# Selected model without trend2
frm_paper_no_t2 <- as.formula("y ~ trend * (dist + sin.d + cos.d) + (trend + dist)*lag1*lag2")
fit_paper_no_t2 <- glm(formula = frm_paper_no_t2,
                       family = binomial,
                       data = stk_df)

# Selected model without trend2*dist
frm_paper_no_tdist <- as.formula("y ~ poly(trend,2) * (sin.d + cos.d) + (trend + dist)*lag1*lag2")
fit_paper_no_tdist <- glm(formula = frm_paper_no_tdist,
                          family = binomial,
                          data = stk_df)

# Selected model without dist*lag1*lag1
frm_paper_no_distlag12 <- as.formula("y ~ poly(trend,2) * (dist + sin.d + cos.d) + trend*lag1*lag2")
fit_paper_no_distlag12 <- glm(formula = frm_paper_no_distlag12,
                              family = binomial,
                              data = stk_df)

# Compare AIC's
AIC(fit_paper, fit_paper_no_t2, fit_paper_no_tdist, fit_paper_no_distlag12)