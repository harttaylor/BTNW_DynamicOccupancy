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
harvest_dist_model <- readRDS("Results/harvest_01_distance_model.rds")
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
additive_harvest_model <- readRDS("Results/harvest_03_main_model.rds")

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

# ===============================================================================
# 6. MODEL COMPARISON USING WAIC
# ===============================================================================

cat("HARVEST MODEL COMPARISON\n")

# Function to calculate WAIC for occupancy models
calculate_waic_for_occupancy <- function(model) {
  if ("y.prob" %in% names(model$sims.list)) {
    y_prob <- model$sims.list$y.prob
    log_lik <- log(y_prob)
    
    n_iter <- dim(log_lik)[1]
    n_sites <- dim(log_lik)[2]
    n_years <- dim(log_lik)[3]
    
    # Skip first year (initial state)
    log_lik_matrix <- matrix(NA, nrow = n_iter, ncol = n_sites * (n_years-1))
    for (i in 1:n_sites) {
      for (j in 2:n_years) {
        log_lik_matrix[, (i-1)*(n_years-1) + (j-1)] <- log_lik[, i, j]
      }
    }
    
    return(waic(log_lik_matrix))
  } else {
    return(NULL)
  }
}

# Load all harvest models
harvest_models <- list(
  "Null" = null_model,
  "Distance Only" = harvest_dist_model,
  "Age Only" = harvest_age_model,
  "Distance + Age" = harvest_main_model,
  "Distance × Age Interaction" = harvest_interact_model
)

# Calculate WAIC for all models
harvest_waic_results <- data.frame(
  Model = character(),
  WAIC = numeric(),
  SE = numeric(),
  dWAIC = numeric(),
  Weight = numeric(),
  stringsAsFactors = FALSE
)

cat("Calculating WAIC for each model...\n")
for(model_name in names(harvest_models)) {
  cat("Processing:", model_name, "\n")
  model <- harvest_models[[model_name]]
  waic_result <- calculate_waic_for_occupancy(model)
  
  if(!is.null(waic_result)) {
    harvest_waic_results <- rbind(harvest_waic_results, data.frame(
      Model = model_name,
      WAIC = waic_result$estimates["waic", "Estimate"],
      SE = waic_result$estimates["waic", "SE"],
      dWAIC = 0,
      Weight = 0
    ))
    cat("✓ WAIC calculated successfully\n")
  } else {
    cat("✗ WAIC calculation failed\n")
  }
}

# Calculate delta WAIC and weights
if(nrow(harvest_waic_results) > 0) {
  min_waic <- min(harvest_waic_results$WAIC)
  harvest_waic_results$dWAIC <- harvest_waic_results$WAIC - min_waic
  harvest_waic_results$Weight <- exp(-0.5 * harvest_waic_results$dWAIC) / 
    sum(exp(-0.5 * harvest_waic_results$dWAIC))
  
  # Order by WAIC (best to worst)
  harvest_waic_results <- harvest_waic_results[order(harvest_waic_results$WAIC), ]
  
  cat("\n=== HARVEST MODEL COMPARISON RESULTS ===\n")
  print(harvest_waic_results)
  
  # Save results
  write.csv(harvest_waic_results, "Results/harvest_model_comparison.csv", row.names = FALSE)
  
  # Best model
  best_model <- harvest_waic_results$Model[1]
  best_weight <- harvest_waic_results$Weight[1]
  
  cat("\nBest supported model:", best_model, "\n")
  cat("Model weight:", round(best_weight, 3), "\n")
  
} else {
  cat("Error: No WAIC values calculated successfully\n")
}

# ===============================================================================
# 7. PREDICTIVE CHECKS
# ===============================================================================
# posterior predictive check
posterior_predictive_annual_occupancy <- function(model, model_name) {
  # Extract true occupancy and predicted occupancy
  z_samples <- model$sims.list$z
  muZ_samples <- model$sims.list$muZ
  
  # Calculate observed annual occupancy (from posterior mean)
  annual_occ_obs <- apply(model$mean$z, 2, mean, na.rm = TRUE)
  
  # Calculate predicted annual occupancy from simulations
  n_samples <- dim(z_samples)[1]
  n_years <- dim(z_samples)[3]
  
  annual_occ_pred <- matrix(NA, nrow = n_samples, ncol = n_years)
  for(s in 1:min(n_samples, 1000)) {  # Use subset for speed
    for(y in 1:n_years) {
      # Simulate occupancy based on predicted probabilities
      z_sim <- rbinom(dim(z_samples)[2], 1, muZ_samples[s, , y])
      annual_occ_pred[s, y] <- mean(z_sim)
    }
  }
  
  # Calculate prediction intervals
  pred_intervals <- apply(annual_occ_pred, 2, quantile, c(0.025, 0.5, 0.975), na.rm = TRUE)
  
  # Create summary
  years <- 1993:(1993 + n_years - 1)
  ppc_data <- data.frame(
    Year = years,
    Observed = annual_occ_obs,
    Predicted_Median = pred_intervals[2, ],
    CI_Lower = pred_intervals[1, ],
    CI_Upper = pred_intervals[3, ]
  )
  
  # Calculate Bayesian p-value for annual occupancy
  p_vals <- rep(NA, n_years)
  for(y in 1:n_years) {
    p_vals[y] <- mean(annual_occ_pred[, y] > annual_occ_obs[y], na.rm = TRUE)
  }
  
  overall_p_val <- mean(p_vals, na.rm = TRUE)
  
  cat(sprintf("Model: %s\n", model_name))
  cat(sprintf("Overall Bayesian p-value: %.3f\n", overall_p_val))
  cat(sprintf("Good fit: %s\n", ifelse(overall_p_val > 0.1 & overall_p_val < 0.9, "Yes", "Check")))
  
  return(list(data = ppc_data, p_value = overall_p_val))
}

# Run for the best harvest model
harvest_interact_ppc <- posterior_predictive_annual_occupancy(harvest_interact_model, "Harvest Interaction")
 