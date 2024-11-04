# This code calculates conditional pseudo-R2 values for linear mixed model
# without kinship matrix. We calculate R2 two ways: (1) for main effects, and
# (2) for variable types (main effects + any interactions containing that main
# effect)

# Load libraries
library(tidyverse); library(readr); library(brms); library(abind);

# Read in Rdata object of no kinship matrix model (model that we are doing most
# of our inference from)
brms_lin_nokin <- read_rds("outputs/phenology_nokin_final.rds")

# Extract posterior samples 
draws <- as_draws(brms_lin_nokin)

# Bind samples of all parameters across all chains in one matrix 
samples_mat <- do.call(rbind, map(draws, ~ abind(.x, along=2)))

# Extract posterior samples of linear predictor 
Y_pred <- t(posterior_linpred(brms_lin_nokin))
n_mcmc <- dim(Y_pred)[2]

# Extract fixed effects matrix and Y
data <- make_standata(formula(brms_lin_nokin), brms_lin_nokin$data)
X <- data$X
col_names <- colnames(X)
var_types <- c("density", "gravel", "pc1", "pc2", "site")
p <- length(var_types)
n <- nrow(X)
Y_mat <- matrix(rep(data$Y, n_mcmc), n, n_mcmc)

# Calculate residual sum of squares for each MCMC iteration for the full model 
RSS <- apply(Y_mat - Y_pred, 2, var)*(n-1)

###
### One variable TYPE at a time 
###

# Table for output 
pseudo_R2_type <- 
  tibble(
    Variable = var_types,
    R2 = NA,
    R2_sd = NA,
  )


for(j in 1:p){
  
  # Get index of every fixed effect (main effects and 2 and/or 3 way interactions) that has to do with variable type j 
  idx <- which(str_detect(col_names, var_types[j]))
  X_tmp <- as.matrix(X[,idx])
  # What are the total sum of squares for a model that includes all fixed and random effects except for those in X_tmp 
  TSS <- apply(Y_mat - Y_pred - X_tmp%*%t(samples_mat[,idx]), 2, var)*(n-1) # could remove (n-1) bc it well cancel but keeping for clarity
  R <- 1 - RSS/TSS
  pseudo_R2_type$R2[j] <- mean(R)
  pseudo_R2_type$R2_sd[j] <- sd(R)
  
}

###
### One variable MAIN effect at a time 
###

# Table for output 
pseudo_R2_MAIN <- 
  tibble(
    Variable = var_types,
    R2 = NA,
    R2_sd = NA,
  )

for(j in 1:p){
  
  # Get index of every fixed effect (main effects and 2 and/or 3 way interactions) that has to do with variable type j 
  idx <- which(str_detect(col_names, var_types[j]) & !str_detect(col_names, ":")) # MAIN effects only for each variable type 
  X_tmp <- as.matrix(X[,idx])
  # What are the total sum of squares for a model that includes all fixed and random effects except for those in X_tmp 
  TSS <- apply(Y_mat - Y_pred - X_tmp%*%t(samples_mat[,idx]), 2, var)*(n-1) # could remove (n-1) bc it well cancel but keeping for clarity
  R <- 1 - RSS/TSS
  pseudo_R2_MAIN$R2[j] <- mean(R)
  pseudo_R2_MAIN$R2_sd[j] <- sd(R)
  
}

pseudo_R2_MAIN
pseudo_R2_type
