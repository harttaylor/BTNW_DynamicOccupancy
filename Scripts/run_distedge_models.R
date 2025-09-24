# ===============================================================================
# RUN DISTANCE TO EDGE MODELS
# ===============================================================================
# Load packages 
library(jagsUI)
library(loo)
library(dplyr)

# Set working directory
#setwd()

# ===============================================================================
# LOAD DATA AND ARRAYS
# ===============================================================================

cat("Loading data...\n")

# Load the standardized distance arrays (from standardization script)
distance_arrays <- readRDS("Data/distance_model_arrays.rds")

# Load the prepared model data (detection data, etc.)
model_data <- readRDS("Data/model_data.rds")

# MCMC settings - using your final settings
settings <- list(
  ni = 40000,    # Number of iterations  
  nt = 1,        # Thinning rate
  nb = 20000,    # Burn-in
  nc = 3         # Number of chains
)

cat("Data loaded successfully\n")
cat("Sites:", length(distance_arrays$sites), "Years:", length(distance_arrays$years), "\n")

# ===============================================================================
# INITIALIZATION FUNCTION
# ===============================================================================

inits <- function() { 
  list(z = apply(model_data$y, c(1, 2), max, na.rm = TRUE))
}

# ===============================================================================
# MODEL 1: NULL MODEL
# ===============================================================================

null_data <- list(
  y = model_data$y, 
  nsite = model_data$nsite, 
  nyear = model_data$nyear, 
  nsurv = model_data$nsurv, 
  J = model_data$J, 
  ind = model_data$ind
)

params_null <- c("psi1", "phi", "gamma", "n.occ", "growthr", "turnover", "p", "N", "z", "muZ", "l.score", "y.prob")

system.time({
  null_model <- jags(
    data = null_data, 
    inits = inits, 
    parameters.to.save = params_null, 
    model.file = "Models/NULL.txt", 
    n.chains = settings$nc, 
    n.thin = settings$nt, 
    n.iter = settings$ni, 
    n.burnin = settings$nb, 
    parallel = TRUE
  )
})

print(null_model)

saveRDS(null_model, "Results/01_null_model.rds")

# ===============================================================================
# MODEL 2: DISTANCE EDGE ONLY (LOG TRANSFORMATION)
# ===============================================================================

dist_edge_data_log <- list(
  y = model_data$y,
  nsite = model_data$nsite,
  nyear = model_data$nyear,
  nsurv = model_data$nsurv,
  J = model_data$J,
  x.psi = distance_arrays$first_year_covariates,
  nbeta.psi = ncol(distance_arrays$first_year_covariates),
  x.phi = distance_arrays$log_with_intercept_phi,
  nbeta.phi = dim(distance_arrays$log_with_intercept_phi)[3],
  x.gamma = distance_arrays$log_with_intercept_gamma,
  nbeta.gamma = dim(distance_arrays$log_with_intercept_gamma)[3],
  x.p = model_data$x_p,
  nbeta.p = model_data$nbeta_p,
  ind = model_data$ind
)

params_dist <- c("beta.psi", "beta.phi", "beta.gamma", "beta.p", "psi", "phi", "gamma", "N", "z", "muZ", "l.score", "y.prob")


system.time({
dist_edge_log_model <- jags(
    data = dist_edge_data_log,
    inits = inits,
    parameters.to.save = params_dist,
    model.file = "Models/DistEdge.txt",
    n.chains = settings$nc,
    n.thin = settings$nt,
    n.iter = settings$ni,
    n.burnin = settings$nb,
    parallel = TRUE
  )
}) 

saveRDS(dist_edge_log_model, "Results/02_dist_edge_log_model.rds")

print(dist_edge_log_model)

# ===============================================================================
# MODEL 3: DISTANCE EDGE ONLY (LINEAR TRANSFORMATION)
# ===============================================================================

dist_edge_data_linear <- list(
  y = model_data$y,
  nsite = model_data$nsite,
  nyear = model_data$nyear,
  nsurv = model_data$nsurv,
  J = model_data$J,
  x.psi = distance_arrays$first_year_covariates,
  nbeta.psi = ncol(distance_arrays$first_year_covariates),
  x.phi = distance_arrays$linear_with_intercept_phi,
  nbeta.phi = dim(distance_arrays$linear_with_intercept_phi)[3],
  x.gamma = distance_arrays$linear_with_intercept_gamma,
  nbeta.gamma = dim(distance_arrays$linear_with_intercept_gamma)[3],
  x.p = model_data$x_p,
  nbeta.p = model_data$nbeta_p,
  ind = model_data$ind
)


system.time({
  dist_edge_linear_model <- jags(
    data = dist_edge_data_linear,
    inits = inits,
    parameters.to.save = params_dist,
    model.file = "Models/DistEdge.txt",
    n.chains = settings$nc,
    n.thin = settings$nt,
    n.iter = settings$ni,
    n.burnin = settings$nb,
    parallel = TRUE
  )
})
print(dist_edge_linear_model)
saveRDS(dist_edge_linear_model, "Results/03_dist_edge_linear_model.rds")

# ===============================================================================
# MODEL 4: LINEAR DISTANCE + RANDOM PATCH
# ===============================================================================

linear_patch_data <- list(
  y = model_data$y,
  nsite = model_data$nsite,
  nyear = model_data$nyear,
  nsurv = model_data$nsurv,
  J = model_data$J,
  x.psi = distance_arrays$first_year_covariates,
  nbeta.psi = ncol(distance_arrays$first_year_covariates),
  x.phi = distance_arrays$linear_no_intercept_phi,
  nbeta.phi = dim(distance_arrays$linear_no_intercept_phi)[3],
  x.gamma = distance_arrays$linear_no_intercept_gamma,
  nbeta.gamma = dim(distance_arrays$linear_no_intercept_gamma)[3],
  x.p = model_data$x_p,
  nbeta.p = model_data$nbeta_p,
  ind = model_data$ind,
  patch = model_data$patch,
  npatch = model_data$npatch
)

# Params to save 
params_linear_patch <- c("beta.psi", "beta.phi", "beta.gamma", "beta.p",
                         "patch.phi", "patch.gamma", "phi", "gamma", 
                         "N", "psi", "z", "muZ", "l.score", "y.prob")
system.time({
  linear_patch_model <- jags(
    data = linear_patch_data,
    inits = inits,
    parameters.to.save = params_linear_patch,
    model.file = "Models/DistEdge_RPatch.txt",
    n.chains = settings$nc,
    n.thin = settings$nt,
    n.iter = settings$ni,
    n.burnin = settings$nb,
    parallel = TRUE
  )
})

print(linear_patch_model)

saveRDS(linear_patch_model, "Results/05_linear_patch_model.rds")

# ===============================================================================
# MODEL 5: LOG DISTANCE + RANDOM PATCH
# ===============================================================================

log_patch_data <- list(
  y = model_data$y,
  nsite = model_data$nsite,
  nyear = model_data$nyear,
  nsurv = model_data$nsurv,
  J = model_data$J,
  x.psi = distance_arrays$first_year_covariates,
  nbeta.psi = ncol(distance_arrays$first_year_covariates),
  x.phi = distance_arrays$log_no_intercept_phi,
  nbeta.phi = dim(distance_arrays$log_no_intercept_phi)[3],
  x.gamma = distance_arrays$log_no_intercept_gamma,
  nbeta.gamma = dim(distance_arrays$log_no_intercept_gamma)[3],
  x.p = model_data$x_p,
  nbeta.p = model_data$nbeta_p,
  ind = model_data$ind,
  patch = model_data$patch,
  npatch = model_data$npatch
)


system.time({
  log_patch_model <- jags(
    data = log_patch_data,
    inits = inits,
    parameters.to.save = params_linear_patch,  # Same params as linear
    model.file = "Models/DistEdge_RPatch.txt",
    n.chains = settings$nc,
    n.thin = settings$nt,
    n.iter = settings$ni,
    n.burnin = settings$nb,
    parallel = TRUE
  )
})

print(log_patch_model)

saveRDS(log_patch_model, "Results/06_log_patch_model.rds")




