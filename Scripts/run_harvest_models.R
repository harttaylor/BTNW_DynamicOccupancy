# ===============================================================================
# COMPLETE HARVEST MODELS ANALYSIS - DISTANCE-BASED APPROACH
# ===============================================================================

library(jagsUI)
library(loo)
library(ggplot2)
library(ROCR)
library(reshape2)

# ===============================================================================
# 1. LOAD DATA AND PREPARE FOR HARVEST MODELS
# ===============================================================================

cat("Loading data and preparing harvest models...\n")

# Load the ecological harvest data (from the preparation script)
model_data <- readRDS("Data/model_data_submodel.rds")

# Load existing null model for comparison
null_model <- readRDS("Results/01_null_model.rds")

# MCMC settings
settings <- list(
  ni = 40000,
  nt = 1,
  nb = 20000,
  nc = 3
)

# Parameters to monitor for harvest models
params_harvest <- c("beta.psi", "beta.phi", "beta.gamma", "beta.p",
                    "patch.phi", "patch.gamma", "mu.phi", "mu.gamma", "psi",
                    "phi", "gamma", "N", "z", "muZ", "l.score", "y.prob")

# Initialization function
inits <- function() { 
  list(z = apply(model_data$y, c(1, 2), max, na.rm = TRUE))
}

cat("Data loaded. Running 4 harvest models...\n")

# ===============================================================================
# 2. MODEL 1: HARVEST DISTANCE ONLY
# ===============================================================================

cat("\n=== HARVEST MODEL 1: Distance Only ===\n")

harvest_dist_data <- list(
  y = model_data$y,
  nsite = model_data$nsite,
  nyear = model_data$nyear,
  nsurv = model_data$nsurv,
  J = model_data$J,
  x.psi = model_data$x_psi,
  nbeta.psi = model_data$nbeta_psi,
  x.phi = model_data$x_phi_harv_dist_submodel,
  nbeta.phi = dim(model_data$x_phi_harv_dist_submodel)[3],
  x.gamma = model_data$x_gamma_harv_dist_submodel,
  nbeta.gamma = dim(model_data$x_gamma_harv_dist_submodel)[3],
  x.p = model_data$x_p,
  nbeta.p = model_data$nbeta_p,
  ind = model_data$ind,
  patch = model_data$patch,
  npatch = model_data$npatch
)

# Verify dimensions
cat("Distance model dimensions:\n")
cat("x.phi:", dim(harvest_dist_data$x.phi), "\n")
cat("nbeta.phi:", harvest_dist_data$nbeta.phi, "\n")

system.time({
  harvest_dist_model <- jags(
    data = harvest_dist_data,
    inits = inits,
    parameters.to.save = params_harvest,
    model.file = "Models/DistEdge_RPatch.txt",  # Using same model as distance-to-edge
    n.chains = settings$nc,
    n.thin = settings$nt,
    n.iter = settings$ni,
    n.burnin = settings$nb,
    parallel = TRUE
  )
})

print(harvest_dist_model)
saveRDS(harvest_dist_model, "Results/harvest_01_distance_model.rds")

# ===============================================================================
# 3. MODEL 2: HARVEST AGE ONLY
# ===============================================================================

cat("\n=== HARVEST MODEL 2: Age Only ===\n")

harvest_age_data <- list(
  y = model_data$y,
  nsite = model_data$nsite,
  nyear = model_data$nyear,
  nsurv = model_data$nsurv,
  J = model_data$J,
  x.psi = model_data$x_psi,
  nbeta.psi = model_data$nbeta_psi,
  x.phi = model_data$x_phi_harv_age_submodel,
  nbeta.phi = dim(model_data$x_phi_harv_age_submodel)[3],
  x.gamma = model_data$x_gamma_harv_age_submodel,
  nbeta.gamma = dim(model_data$x_gamma_harv_age_submodel)[3],
  x.p = model_data$x_p,
  nbeta.p = model_data$nbeta_p,
  ind = model_data$ind,
  patch = model_data$patch,
  npatch = model_data$npatch
)

cat("Age model dimensions:\n")
cat("x.phi:", dim(harvest_age_data$x.phi), "\n")

system.time({
  harvest_age_model <- jags(
    data = harvest_age_data,
    inits = inits,
    parameters.to.save = params_harvest,
    model.file = "Models/DistEdge_RPatch.txt",
    n.chains = settings$nc,
    n.thin = settings$nt,
    n.iter = settings$ni,
    n.burnin = settings$nb,
    parallel = TRUE
  )
})

print(harvest_age_model)
saveRDS(harvest_age_model, "Results/harvest_02_age_model.rds")

# ===============================================================================
# 4. MODEL 3: HARVEST MAIN EFFECTS (DISTANCE + AGE)
# ===============================================================================

cat("\n=== HARVEST MODEL 3: Distance + Age (Main Effects) ===\n")

harvest_main_data <- list(
  y = model_data$y,
  nsite = model_data$nsite,
  nyear = model_data$nyear,
  nsurv = model_data$nsurv,
  J = model_data$J,
  x.psi = model_data$x_psi,
  nbeta.psi = model_data$nbeta_psi,
  x.phi = model_data$x_phi_harv_main_submodel,
  nbeta.phi = dim(model_data$x_phi_harv_main_submodel)[3],
  x.gamma = model_data$x_gamma_harv_main_submodel,
  nbeta.gamma = dim(model_data$x_gamma_harv_main_submodel)[3],
  x.p = model_data$x_p,
  nbeta.p = model_data$nbeta_p,
  ind = model_data$ind,
  patch = model_data$patch,
  npatch = model_data$npatch
)

cat("Main effects model dimensions:\n")
cat("x.phi:", dim(harvest_main_data$x.phi), "\n")

system.time({
  harvest_main_model <- jags(
    data = harvest_main_data,
    inits = inits,
    parameters.to.save = params_harvest,
    model.file = "Models/DistEdge_RPatch.txt",
    n.chains = settings$nc,
    n.thin = settings$nt,
    n.iter = settings$ni,
    n.burnin = settings$nb,
    parallel = TRUE
  )
})

print(harvest_main_model)
saveRDS(harvest_main_model, "Results/harvest_03_main_model.rds")

# ===============================================================================
# 5. MODEL 4: HARVEST INTERACTION (DISTANCE × AGE)
# ===============================================================================

cat("\n=== HARVEST MODEL 4: Distance × Age Interaction ===\n")

harvest_interact_data <- list(
  y = model_data$y,
  nsite = model_data$nsite,
  nyear = model_data$nyear,
  nsurv = model_data$nsurv,
  J = model_data$J,
  x.psi = model_data$x_psi,
  nbeta.psi = model_data$nbeta_psi,
  x.phi = model_data$x_phi_harv_interact_submodel,
  nbeta.phi = dim(model_data$x_phi_harv_interact_submodel)[3],
  x.gamma = model_data$x_gamma_harv_interact_submodel,
  nbeta.gamma = dim(model_data$x_gamma_harv_interact_submodel)[3],
  x.p = model_data$x_p,
  nbeta.p = model_data$nbeta_p,
  ind = model_data$ind,
  patch = model_data$patch,
  npatch = model_data$npatch
)

cat("Interaction model dimensions:\n")
cat("x.phi:", dim(harvest_interact_data$x.phi), "\n")
cat("Should have 3 covariates: distance, age, interaction\n")

system.time({
  harvest_interact_model <- jags(
    data = harvest_interact_data,
    inits = inits,
    parameters.to.save = params_harvest,
    model.file = "Models/DistEdge_RPatch.txt",
    n.chains = settings$nc,
    n.thin = settings$nt,
    n.iter = settings$ni,
    n.burnin = settings$nb,
    parallel = TRUE
  )
})

print(harvest_interact_model)
saveRDS(harvest_interact_model, "Results/harvest_04_interaction_model.rds")


 