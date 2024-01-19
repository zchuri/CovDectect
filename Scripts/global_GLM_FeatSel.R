#!/usr/bin/env Rscript

# Clear workspace
rm(list = ls())
# Save original parameters
op <- par()

# Load libraries
if(!require(devtools)) install.packages("devtools"); library(devtools)
if(!require(dplyr)) install.packages("dplyr"); library(dplyr)
if(!require(glmnet)) install.packages("glmnet"); library(glmnet)
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
                     dist = rep(0,df_ln),
                     day = rep(0,df_ln),
                     sin.d = rep(0,df_ln),
                     cos.d = rep(0,df_ln),
                     t = rep(0,df_ln),
                     t2 = rep(0,df_ln),
                     lag1 = rep(0,df_ln),
                     lag2 = rep(0,df_ln),
                     lag12 = rep(0,df_ln))

# Begin to stack
for(ss in 1:bin_dim[3]){
  
  cat(paste0(ss, "."))
  # Stacked indices
  idx <- (1+(ss-1)*st_ln):(ss*st_ln)
  # Extract records
  stk_df$y[idx] <- c(t(bin_df[,,ss]))
  # Station name
  stk_df$st[idx] <- geo_df$abb[ss]
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
stk_df$t2 <- stk_df$t^2
stk_df$lag12 <- stk_df$lag1*stk_df$lag2
stk_df$st <- factor(stk_df$st)

# Fit single term model
frm_one <- as.formula("y ~ t + t2 + sin.d + cos.d + dist + lag1 + lag2 + lag12")
tic <- Sys.time()
fit_one <- glm(formula = frm_one,
               family = binomial,
               data = stk_df)
toc <- Sys.time()
print(toc-tic)
# Show summary
summary(fit_one)

#############################################################
# Generate order two interactions (manually)
#############################################################

# (Pseudo)all possible two order interactions
frm_two <- as.formula("y ~ t + t2 + sin.d + cos.d + dist + lag1 + lag2 + lag12 +
                      t:sin.d + t:cos.d + t:dist + t:lag1 + t:lag2 + t:lag12 +
                      t2:sin.d + t2:cos.d + t2:dist + t2:lag1 + t2:lag2 + t2:lag12 +
                      dist:sin.d + dist:cos.d + dist:lag1 + dist:lag2 + dist:lag12 +
                      lag1:sin.d + lag1:cos.d +
                      lag2:sin.d + lag2:cos.d +
                      lag12:sin.d + lag12:cos.d")


# Extract Zaragoza observatory (9)
bin_zgz <- data.frame(y = c(t(bin_df[,,9])),
                      t = rep(x = 1:bin_dim[1], each = bin_dim[2]),
                      day = rep(x = 1:bin_dim[2], times = bin_dim[1]))
tic <- Sys.time()
fit_two <- glm(formula = frm_two,
               family = binomial,
               data = stk_df)
toc <- Sys.time()
print(toc-tic)
summary(fit_two)

#############################################################
# Feature Selection via LASSO
#############################################################

# https://www.r-bloggers.com/2020/05/quick-tutorial-on-lasso-regression-with-example/
# Extract matrices
x_vars <- model.matrix(frm_two, stk_df)[,-1]
y_var <- stk_df$y
lambda_seq <- 10^seq(2, -2, by = -.1)

# Perform k-fold cross-validation to find optimal lambda value
# https://www.statology.org/lasso-regression-in-r/
fit_cv <- cv.glmnet(x_vars, y_var, alpha = 1)