# ===============================================================================
# HARVEST DATA PREPARATION - SUB-MODEL APPROACH
# This script prepares harvest-specific covariates for dynamic occupancy models
# ===============================================================================
# Load packages 
library(jagsUI)
library(loo)
library(dplyr)

# Load model data
model_data <- readRDS("Data/model_data.rds")
yearly_covariates <- read.csv("Data/harvest_yearly_covs.csv") 
yearly_covariates <- yearly_covariates %>% filter(SS %in% model_data$sites)

# ===============================================================================
# CREATE SUB-MODEL STRUCTURE THROUGH DATA
# ===============================================================================

# Step 1: Create harvest indicator (within 565m = has harvest age data)
yearly_covariates$has_harvest_565m <- as.numeric(yearly_covariates$NEAR.DIST.harvest <= 565)

# Step 2: Distance effects for ALL sites (full ecological gradient)
yearly_covariates$harv_dist_std <- scale(yearly_covariates$NEAR.DIST.harvest)[,1]

# Step 3: Age effects - actual age where harvest exists, 0 elsewhere
# This is the key: 0 = "no age effect" (not "youngest age")
yearly_covariates$harv_age_std <- ifelse(
  yearly_covariates$has_harvest_565m == 1,
  scale(yearly_covariates$MEANAGE.565.harvest[yearly_covariates$has_harvest_565m == 1])[,1],
  0  # No age effect when no nearby harvest
)

# Step 4: Interaction automatically handled (distance × age)
yearly_covariates$harv_interaction_std <- yearly_covariates$harv_dist_std * yearly_covariates$harv_age_std

# Print summary of the sub-model structure
cat("SUB-MODEL STRUCTURE CREATED:\n")
cat("Sites with harvest within 565m:", sum(yearly_covariates$has_harvest_565m), "site-years\n")
cat("Sites without harvest within 565m:", sum(1 - yearly_covariates$has_harvest_565m), "site-years\n")
cat("Distance effects applied to: ALL site-years\n")
cat("Age effects applied to:", sum(yearly_covariates$harv_age_std != 0), "site-years\n")
cat("Interaction effects applied to:", sum(yearly_covariates$harv_interaction_std != 0), "site-years\n\n")

# ===============================================================================
# CREATE MODEL ARRAYS (SAME STRUCTURE AS BEFORE)
# ===============================================================================

NSITE <- model_data$nsite
NYEAR <- model_data$nyear
actual_years_map <- if(min(model_data$years) == 0) model_data$years + 1993 else model_data$years

# Function to create harvest arrays (same as before!)
create_harvest_arrays <- function(covariate_names) {
  num_covs <- length(covariate_names)
  
  x_phi <- array(NA, dim = c(NSITE, NYEAR, num_covs))
  x_gamma <- array(NA, dim = c(NSITE, NYEAR, num_covs))
  
  for(i in 1:NSITE) {
    for(j in 1:NYEAR) {
      actual_year <- actual_years_map[j]
      
      site_data <- yearly_covariates[yearly_covariates$SS == model_data$sites[i] & 
                                       yearly_covariates$YEAR == actual_year, ]
      
      if(nrow(site_data) > 0) {
        for(k in 1:num_covs) {
          covariate_value <- site_data[[covariate_names[k]]]
          x_phi[i, j, k] <- covariate_value
          x_gamma[i, j, k] <- covariate_value
        }
      } else {
        # Use 0 for missing values
        for(k in 1:num_covs) {
          x_phi[i, j, k] <- 0
          x_gamma[i, j, k] <- 0
        }
      }
    }
  }
  
  return(list(x_phi = x_phi, x_gamma = x_gamma))
}

# Create model variants (exactly same process as before)
cat("Creating harvest model arrays...\n")

# Distance only
dist_arrays <- create_harvest_arrays("harv_dist_std")
x_phi_harv_dist <- dist_arrays$x_phi
x_gamma_harv_dist <- dist_arrays$x_gamma

# Age only (will be 0 for sites without harvest)
age_arrays <- create_harvest_arrays("harv_age_std")
x_phi_harv_age <- age_arrays$x_phi
x_gamma_harv_age <- age_arrays$x_gamma

# Main effects (distance + age)
main_arrays <- create_harvest_arrays(c("harv_dist_std", "harv_age_std"))
x_phi_harv_main <- main_arrays$x_phi
x_gamma_harv_main <- main_arrays$x_gamma

# Interaction (distance + age + interaction)
interact_arrays <- create_harvest_arrays(c("harv_dist_std", "harv_age_std", "harv_interaction_std"))
x_phi_harv_interact <- interact_arrays$x_phi
x_gamma_harv_interact <- interact_arrays$x_gamma

# ===============================================================================
# VALIDATION
# ===============================================================================

cat("VALIDATION OF SUB-MODEL STRUCTURE:\n")
cat("Distance model - range:", range(x_phi_harv_dist, na.rm = TRUE), "\n")
cat("Age model - range:", range(x_phi_harv_age, na.rm = TRUE), "\n")
cat("Age model - zero values:", sum(x_phi_harv_age == 0, na.rm = TRUE), "out of", length(x_phi_harv_age[!is.na(x_phi_harv_age)]), "\n")
cat("Interaction model - range:", range(x_phi_harv_interact, na.rm = TRUE), "\n")
cat("Interaction model - zero values:", sum(x_phi_harv_interact[,,3] == 0, na.rm = TRUE), "interaction terms\n\n")

# ===============================================================================
# UPDATE MODEL DATA OBJECT
# ===============================================================================

# Add to existing model_data (same structure as before)
model_data$x_phi_harv_dist_submodel <- x_phi_harv_dist
model_data$x_gamma_harv_dist_submodel <- x_gamma_harv_dist
model_data$x_phi_harv_age_submodel <- x_phi_harv_age
model_data$x_gamma_harv_age_submodel <- x_gamma_harv_age
model_data$x_phi_harv_main_submodel <- x_phi_harv_main
model_data$x_gamma_harv_main_submodel <- x_gamma_harv_main
model_data$x_phi_harv_interact_submodel <- x_phi_harv_interact
model_data$x_gamma_harv_interact_submodel <- x_gamma_harv_interact

# Save updated data
saveRDS(model_data, "Data/model_data_submodel.rds")

