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

# ===============================================================================
# MODEL COMPARISON
# ===============================================================================

# Function to calculate WAIC 
calculate_waic_occupancy <- function(model) {
  tryCatch({
    if ("y.prob" %in% names(model$sims.list)) {
      y_prob <- model$sims.list$y.prob
      log_lik <- log(y_prob)
      
      # Reshape for WAIC (skip first year)
      n_iter <- dim(log_lik)[1]
      n_sites <- dim(log_lik)[2] 
      n_years <- dim(log_lik)[3]
      
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
  }, error = function(e) {
    cat("Error calculating WAIC:", e$message, "\n")
    return(NULL)
  })
}

# Load all models 
model_files <- list(
  "Null" = "Results/01_null_model.rds",
  "Distance Edge Only (Log)" = "Results/02_dist_edge_log_model.rds",
  "Distance Edge Only (Linear)" = "Results/03_dist_edge_linear_model.rds",
  "Linear Distance + Random Patch" = "Results/05_linear_patch_model.rds",
  "Log Distance + Random Patch" = "Results/06_log_patch_model.rds"
)

# Initialize models list
models <- list()

# Load models from files
for(model_name in names(model_files)) {
  if(file.exists(model_files[[model_name]])) {
    models[[model_name]] <- readRDS(model_files[[model_name]])
    cat("✓ Loaded", model_name, "\n")
  } else {
    cat("✗ Missing", model_name, "\n")
  }
}


# Check if models loaded
if(length(models) == 0) {
  cat("✗ No models loaded - check file paths and ensure model files exist\n")
  stop("Cannot proceed without loaded models")
}

cat("\nLoaded", length(models), "models successfully\n")

# Calculate WAIC for all models
waic_results <- data.frame(
  Model = character(),
  WAIC = numeric(),
  SE = numeric(),
  dWAIC = numeric(),
  Weight = numeric(),
  stringsAsFactors = FALSE
)

for(model_name in names(models)) {
  cat("Calculating WAIC for:", model_name, "\n")
  model <- models[[model_name]]
  waic_result <- calculate_waic_occupancy(model)
  
  if(!is.null(waic_result)) {
    waic_results <- rbind(waic_results, data.frame(
      Model = model_name,
      WAIC = waic_result$estimates["waic", "Estimate"],
      SE = waic_result$estimates["waic", "SE"],
      dWAIC = 0,  # Will calculate below
      Weight = 0  # Will calculate below
    ))
    cat("✓ WAIC calculated successfully\n")
  } else {
    cat("✗ WAIC calculation failed\n")
  }
}

# Calculate delta WAIC and weights
if(nrow(waic_results) > 0) {
  min_waic <- min(waic_results$WAIC)
  waic_results$dWAIC <- waic_results$WAIC - min_waic
  waic_results$Weight <- exp(-0.5 * waic_results$dWAIC) / sum(exp(-0.5 * waic_results$dWAIC))
  waic_results <- waic_results[order(waic_results$WAIC), ]
  
  cat("\n=== DISTANCE-TO-EDGE MODEL COMPARISON RESULTS ===\n")
  print(waic_results)
  
  # Determine best functional form for harvest models
  best_model <- waic_results$Model[1]
  if(grepl("Linear", best_model)) {
    best_functional_form <- "linear"
  } else if(grepl("Log", best_model)) {
    best_functional_form <- "log"
  } else if(grepl("Patch", best_model)) {
    # If patch-only wins, default to log (your original preference)
    best_functional_form <- "log"
  } else {
    # If null wins, still default to log
    best_functional_form <- "log"
  }
  
  cat("\n*** RESULTS SUMMARY ***\n")
  cat("Best performing model:", best_model, "\n")
  cat("Model weight:", round(waic_results$Weight[1], 3), "\n")
  cat("Recommended functional form for harvest models:", toupper(best_functional_form), "\n")
  
  # Create Results directory if it doesn't exist
  if(!dir.exists("Results")) {
    dir.create("Results")
    cat("✓ Created Results directory\n")
  }
  
  # Save results
  write.csv(waic_results, "Results/distance_edge_model_comparison.csv", row.names = FALSE)
  saveRDS(best_functional_form, "Results/best_functional_form.rds")
  
  cat("✓ Model comparison results saved\n")
} else {
  cat("✗ No WAIC results calculated - check model outputs\n")
  
  # Debug information
  cat("\nDEBUG INFO:\n")
  cat("Number of models loaded:", length(models), "\n")
  if(length(models) > 0) {
    for(model_name in names(models)) {
      model <- models[[model_name]]
      cat("Model:", model_name, "\n")
      cat("  - Has sims.list:", "sims.list" %in% names(model), "\n")
      if("sims.list" %in% names(model)) {
        cat("  - Has y.prob:", "y.prob" %in% names(model$sims.list), "\n")
        if("y.prob" %in% names(model$sims.list)) {
          cat("  - y.prob dimensions:", dim(model$sims.list$y.prob), "\n")
        }
      }
    }
  }
}
print(waic_results)



