# ===============================================================================
# MODEL DIAGNOSTICS - DISTANCE-TO-EDGE AND HARVEST MODELS
# ===============================================================================
library(jagsUI)
library(loo)
library(ggplot2)
library(ROCR)
library(ape)

# ===============================================================================
# 1. BAYESIAN P-VALUES FUNCTION
# ===============================================================================

calculate_bayesian_pvalues_comprehensive <- function(model, model_name) {
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

analyze_spatial_autocorr_comprehensive <- function(model, coords, model_name) {
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
# 4. MODEL DIAGNOSTICS FOR DISTANCE-TO-EDGE MODELS
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

# Run comprehensive diagnostics
distance_diagnostics <- list()

for(model_name in names(distance_models)) {
  cat("\nProcessing", model_name, "...\n")
  
  model <- distance_models[[model_name]]
  
  
  # Bayesian p-values
  pval_result <- calculate_bayesian_pvalues_comprehensive(model, model_name)
  if(!is.null(pval_result)) {
    cat("  Bayesian p-value:", round(pval_result$overall_p_value, 3), "-", pval_result$assessment, "\n")
  }
  
  # Spatial autocorrelation
  spatial_result <- analyze_spatial_autocorr_comprehensive(model, 
                                                           sitecoords[, c("longitude", "latitude")], 
                                                           model_name)
  if(!is.null(spatial_result)) {
    cat("  Spatial autocorr: Max |Moran's I| =", round(spatial_result$summary$max_moran, 3), 
        "-", spatial_result$summary$assessment, "\n")
  }
  
  # Store results
  distance_diagnostics[[model_name]] <- list(
    bayesian_pval = pval_result,
    spatial = spatial_result
  )
}
print(distance_diagnostics)

# ===============================================================================
# 5. DIAGNOSTICS FOR HARVEST MODELS
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

# Run comprehensive diagnostics
harvest_diagnostics <- list()

for(model_name in names(harvest_models)) {
  cat("\nProcessing", model_name, "...\n")
  
  model <- harvest_models[[model_name]]
  
  # Bayesian p-values
  pval_result <- calculate_bayesian_pvalues_comprehensive(model, model_name)
  if(!is.null(pval_result)) {
    cat("  Bayesian p-value:", round(pval_result$overall_p_value, 3), "-", pval_result$assessment, "\n")
  }
  
  # Spatial autocorrelation
  spatial_result <- analyze_spatial_autocorr_comprehensive(model, 
                                                           sitecoords[, c("longitude", "latitude")], 
                                                           model_name)
  if(!is.null(spatial_result)) {
    cat("  Spatial autocorr: Max |Moran's I| =", round(spatial_result$summary$max_moran, 3), 
        "-", spatial_result$summary$assessment, "\n")
  }
  
  # Store results
  harvest_diagnostics[[model_name]] <- list(
    bayesian_pval = pval_result,
    spatial = spatial_result
  )
}

# ===============================================================================
# 6. SUMMARY TABLES
# ===============================================================================

# Model Fit Summary Table
create_fit_summary <- function(diagnostics_list) {
  fit_data <- lapply(diagnostics_list, function(x) {
    list(
      model = x$bayesian_pval$model,
      bayesian_p = x$bayesian_pval$overall_p_value,
      bayesian_assessment = x$bayesian_pval$assessment,
      spatial_max = x$spatial$summary$max_moran,
      spatial_assessment = x$spatial$summary$assessment
    )
  })
  
  fit_data <- fit_data[!sapply(fit_data, function(x) any(sapply(x, is.null)))]
  
  if(length(fit_data) > 0) {
    fit_summary <- do.call(rbind, lapply(fit_data, function(x) {
      data.frame(
        Model = x$model,
        Bayesian_P_Value = round(x$bayesian_p, 3),
        Fit_Assessment = x$bayesian_assessment,
        Max_Morans_I = round(x$spatial_max, 3),
        Spatial_Assessment = x$spatial_assessment
      )
    }))
    return(fit_summary)
  }
  return(NULL)
}

# Generate summary tables
cat("\n=== DISTANCE-TO-EDGE MODEL SUMMARY ===\n")
distance_fit_summary <- create_fit_summary(distance_diagnostics)
if(!is.null(distance_fit_summary)) {
  print(distance_fit_summary)
  write.csv(distance_fit_summary, "Results/distance_edge_fit_summary.csv", row.names = FALSE)
}

cat("\n=== HARVEST MODEL SUMMARY ===\n")
harvest_fit_summary <- create_fit_summary(harvest_diagnostics)
if(!is.null(harvest_fit_summary)) {
  print(harvest_fit_summary)
  write.csv(harvest_fit_summary, "Results/harvest_fit_summary.csv", row.names = FALSE)
}

cat("\n✓ Comprehensive model diagnostics complete!\n")
cat("Summary tables saved to Results/ directory\n")