# ===============================================================================
# MODEL DIAGNOSTICS - DISTANCE-TO-EDGE AND HARVEST MODELS
# ===============================================================================
library(jagsUI)
library(loo)
library(ggplot2)
library(ROCR)
library(ape)

# ===============================================================================
# 1. WAIC CALCULATION FUNCTION
# ===============================================================================
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

# ===============================================================================
# 2. BAYESIAN P-VALUES FUNCTION
# ===============================================================================
calculate_bayesian_pvalues <- function(model, model_name) {
  # Extract model components
  if(!"z" %in% names(model$sims.list) || !"muZ" %in% names(model$sims.list)) {
    cat("Cannot calculate Bayesian p-values for", model_name, "\n")
    return(NULL)
  }
  
  z_samples <- model$sims.list$z
  muZ_samples <- model$sims.list$muZ
  
  n_samples <- dim(z_samples)[1]
  n_sites <- dim(z_samples)[2]
  n_years <- dim(z_samples)[3]
  
  # Calculate observed test statistic: Sum of occupancy states by year
  N_observed <- apply(z_samples, c(1, 3), sum)
  
  # Generate replicated datasets
  N_replicated <- matrix(NA, nrow = n_samples, ncol = n_years)
  
  for (s in 1:n_samples) {
    for (y in 1:n_years) {
      z_rep <- rbinom(n_sites, 1, muZ_samples[s, , y])
      N_replicated[s, y] <- sum(z_rep)
    }
  }
  
  # Calculate Bayesian p-values
  p_values <- matrix(NA, nrow = n_samples, ncol = n_years)
  
  for (s in 1:n_samples) {
    for (y in 1:n_years) {
      p_values[s, y] <- ifelse(N_replicated[s, y] > N_observed[s, y], 1, 0)
    }
  }
  
  mean_p_values <- colMeans(p_values)
  overall_p_value <- mean(mean_p_values, na.rm = TRUE)
  
  return(list(
    model = model_name,
    overall_p_value = overall_p_value,
    by_year_p_values = mean_p_values,
    assessment = ifelse(overall_p_value > 0.1 & overall_p_value < 0.9, "Good fit", 
                        ifelse(overall_p_value <= 0.1, "Overfitted", "Underfitted"))
  ))
}

# ===============================================================================
# 3. SPATIAL AUTOCORRELATION FUNCTION
# ===============================================================================
analyze_spatial_autocorr <- function(model, coords, model_name) {
  # Extract Pearson residuals
  if(!"z" %in% names(model$mean) || !"muZ" %in% names(model$mean)) {
    cat("Cannot extract residuals for", model_name, "\n")
    return(NULL)
  }
  
  z_obs <- model$mean$z
  z_exp <- model$mean$muZ
  
  # Calculate Pearson residuals
  variance <- z_exp * (1 - z_exp)
  variance[variance < 1e-6] <- 1e-6
  pearson_residuals <- (z_obs - z_exp) / sqrt(variance)
  
  n_years <- ncol(pearson_residuals)
  
  # Create distance matrix and weights
  dist_matrix <- as.matrix(dist(coords))
  diag(dist_matrix) <- max(dist_matrix) * 10
  w <- 1 / dist_matrix
  w <- w / rowSums(w)
  
  morans_results <- data.frame(
    Year = 1:n_years,
    Moran_I = NA,
    p_value = NA,
    Significant = NA
  )
  
  # Calculate Moran's I for each year
  for(year in 1:n_years) {
    year_residuals <- pearson_residuals[, year]
    
    if(sum(!is.na(year_residuals)) >= 10) {
      valid_sites <- !is.na(year_residuals)
      year_residuals_clean <- year_residuals[valid_sites]
      w_year <- w[valid_sites, valid_sites]
      w_year <- w_year / rowSums(w_year)
      
      moran_result <- try(ape::Moran.I(year_residuals_clean, w_year), silent = TRUE)
      
      if(!inherits(moran_result, "try-error")) {
        morans_results$Moran_I[year] <- moran_result$observed
        morans_results$p_value[year] <- moran_result$p.value
        morans_results$Significant[year] <- moran_result$p.value < 0.05
      }
    }
  }
  
  # Summary statistics
  summary_stats <- list(
    model = model_name,
    mean_moran = mean(abs(morans_results$Moran_I), na.rm = TRUE),
    max_moran = max(abs(morans_results$Moran_I), na.rm = TRUE),
    prop_significant = mean(morans_results$Significant, na.rm = TRUE),
    assessment = ifelse(max(abs(morans_results$Moran_I), na.rm = TRUE) < 0.1, "Acceptable", 
                        ifelse(max(abs(morans_results$Moran_I), na.rm = TRUE) < 0.2, "Weak", "Moderate"))
  )
  
  return(list(
    results = morans_results,
    summary = summary_stats
  ))
}

# ===============================================================================
# 4. COMPREHENSIVE VALIDATION FUNCTION
# ===============================================================================
validate_model <- function(model, model_name, coords = NULL, verbose = TRUE) {
  if(verbose) cat("\nProcessing", model_name, "...\n")
  
  validation_results <- list(model_name = model_name)
  
  # 1. WAIC calculation
  waic_result <- calculate_waic_occupancy(model)
  if(!is.null(waic_result)) {
    validation_results$waic <- waic_result$estimates["waic", "Estimate"]
    validation_results$waic_se <- waic_result$estimates["waic", "SE"]
    if(verbose) cat("  WAIC:", round(validation_results$waic, 2), "\n")
  }
  
  # 2. Bayesian p-values
  pval_result <- calculate_bayesian_pvalues(model, model_name)
  if(!is.null(pval_result)) {
    validation_results$bayesian_pval = pval_result
    if(verbose) cat("  Bayesian p-value:", round(pval_result$overall_p_value, 3), 
                    "-", pval_result$assessment, "\n")
  }
  
  # 3. Spatial autocorrelation (if coordinates provided)
  if(!is.null(coords)) {
    spatial_result <- analyze_spatial_autocorr(model, coords, model_name)
    if(!is.null(spatial_result)) {
      validation_results$spatial = spatial_result
      if(verbose) cat("  Spatial autocorr: Max |Moran's I| =", 
                      round(spatial_result$summary$max_moran, 3), 
                      "-", spatial_result$summary$assessment, "\n")
    }
  }
  
  return(validation_results)
}

# ===============================================================================
# 5. MODEL COMPARISON FUNCTION
# ===============================================================================
create_model_comparison <- function(validation_results) {
  # Extract WAIC values
  waic_data <- data.frame(
    Model = names(validation_results),
    WAIC = sapply(validation_results, function(x) x$waic),
    SE = sapply(validation_results, function(x) x$waic_se),
    stringsAsFactors = FALSE
  )
  
  # Remove rows with missing WAIC
  waic_data <- waic_data[!is.na(waic_data$WAIC), ]
  
  if(nrow(waic_data) > 0) {
    # Calculate delta WAIC and weights
    min_waic <- min(waic_data$WAIC)
    waic_data$dWAIC <- waic_data$WAIC - min_waic
    waic_data$Weight <- exp(-0.5 * waic_data$dWAIC) / sum(exp(-0.5 * waic_data$dWAIC))
    waic_data <- waic_data[order(waic_data$WAIC), ]
    
    return(waic_data)
  }
  return(NULL)
}

# ===============================================================================
# 6. FIT SUMMARY FUNCTION
# ===============================================================================
create_fit_summary <- function(validation_results) {
  fit_data <- lapply(validation_results, function(x) {
    list(
      model = x$model_name,
      bayesian_p = if(!is.null(x$bayesian_pval)) x$bayesian_pval$overall_p_value else NA,
      bayesian_assessment = if(!is.null(x$bayesian_pval)) x$bayesian_pval$assessment else NA,
      spatial_max = if(!is.null(x$spatial)) x$spatial$summary$max_moran else NA,
      spatial_assessment = if(!is.null(x$spatial)) x$spatial$summary$assessment else NA
    )
  })
  
  # Remove entries with all missing values
  fit_data <- fit_data[!sapply(fit_data, function(x) all(is.na(c(x$bayesian_p, x$spatial_max))))]
  
  if(length(fit_data) > 0) {
    fit_summary <- do.call(rbind, lapply(fit_data, function(x) {
      data.frame(
        Model = x$model,
        Bayesian_P_Value = round(x$bayesian_p, 3),
        Fit_Assessment = x$bayesian_assessment,
        Max_Morans_I = round(x$spatial_max, 3),
        Spatial_Assessment = x$spatial_assessment,
        stringsAsFactors = FALSE
      )
    }))
    return(fit_summary)
  }
  return(NULL)
}

# ===============================================================================
# 7. DISTANCE-TO-EDGE MODEL ANALYSIS
# ===============================================================================
cat("=== DISTANCE-TO-EDGE MODEL DIAGNOSTICS ===\n")

# Load distance-to-edge models
distance_models <- list(
  "Null" = readRDS("Results/01_null_model.rds"),
  "Distance Edge Only (Log)" = readRDS("Results/02_dist_edge_log_model.rds"),
  "Distance Edge Only (Linear)" = readRDS("Results/03_dist_edge_linear_model.rds"),
  "Linear Distance + Random Patch" = readRDS("Results/05_linear_patch_model.rds"),
  "Log Distance + Random Patch" = readRDS("Results/06_log_patch_model.rds")
)

# Load site coordinates
sitecoords <- read.csv("Data/thinnedpoints.csv")

# Validate all distance models
distance_results <- list()
for(model_name in names(distance_models)) {
  distance_results[[model_name]] <- validate_model(
    distance_models[[model_name]], 
    model_name, 
    coords = sitecoords[, c("longitude", "latitude")],
    verbose = TRUE
  )
}

# Create comparison tables
distance_waic <- create_model_comparison(distance_results)
distance_fit <- create_fit_summary(distance_results)

if(!is.null(distance_waic)) {
  cat("\n=== DISTANCE-TO-EDGE MODEL COMPARISON ===\n")
  print(distance_waic)
  write.csv(distance_waic, "Results/distance_edge_model_comparison.csv", row.names = FALSE)
}

if(!is.null(distance_fit)) {
  cat("\n=== DISTANCE-TO-EDGE FIT SUMMARY ===\n")
  print(distance_fit)
  write.csv(distance_fit, "Results/distance_edge_fit_summary.csv", row.names = FALSE)
}

# ===============================================================================
# 8. HARVEST MODEL ANALYSIS
# ===============================================================================
cat("\n=== HARVEST MODEL DIAGNOSTICS ===\n")

# Load harvest models
harvest_models <- list(
  "Null" = readRDS("Results/01_null_model.rds"),
  "Distance Only" = readRDS("Results/harvest_01_distance_model.rds"),
  "Age Only" = readRDS("Results/harvest_02_age_model.rds"),
  "Distance + Age" = readRDS("Results/harvest_03_main_model.rds"),
  "Distance × Age Interaction" = readRDS("Results/harvest_04_interaction_model.rds")
)

# Validate all harvest models
harvest_results <- list()
for(model_name in names(harvest_models)) {
  harvest_results[[model_name]] <- validate_model(
    harvest_models[[model_name]], 
    model_name, 
    coords = sitecoords[, c("longitude", "latitude")],
    verbose = TRUE
  )
}

# Create comparison tables
harvest_waic <- create_model_comparison(harvest_results)
harvest_fit <- create_fit_summary(harvest_results)

if(!is.null(harvest_waic)) {
  cat("\n=== HARVEST MODEL COMPARISON ===\n")
  print(harvest_waic)
  write.csv(harvest_waic, "Results/harvest_model_comparison.csv", row.names = FALSE)
}

if(!is.null(harvest_fit)) {
  cat("\n=== HARVEST FIT SUMMARY ===\n")
  print(harvest_fit)
  write.csv(harvest_fit, "Results/harvest_fit_summary.csv", row.names = FALSE)
}

