# Dynamic Occupancy Models for Black-throated Green Warbler (*Setophaga virens*) at Calling Lake, Alberta

## Overview

This repository contains data and code for analyzing 25 years (1993-2018) of Black-throated Green Warbler occupancy dynamics in response to anthropogenic edge effects in Alberta's boreal forest. The analysis uses hierarchical Bayesian dynamic occupancy models to estimate colonization and persistence probabilities while accounting for imperfect detection.

**Citation:** [manuscript citation here]

## Software Requirements

- **R version:** 4.4.2 or higher
- **JAGS:** 4.3.0 or higher (must be installed separately from https://mcmc-jags.sourceforge.io/)
- **R packages:**
  - `jagsUI` (1.5.2+): Interface for JAGS
  - `loo` (2.5.1+): Model comparison using WAIC
  - `dplyr` (1.1.0+): Data manipulation
  - `ggplot2` (3.4.0+): Visualization
  - `ROCR` (1.0-11+): AUC calculations
  - `ape` (5.7+): Spatial autocorrelation analysis

Install required R packages:
```r
install.packages(c("jagsUI", "loo", "dplyr", "ggplot2", "ROCR", "ape"))
```

## Data Files

### Core Data (Data/)
- `model_data.rds`: Complete data object for distance-to-edge models containing:
  - Detection array (114 sites × 25 years × 4 visits)
  - Standardized covariates for initial occupancy, detection, persistence, and colonization
  - Site-patch assignments for random effects
  - Survey metadata (Julian day, time of day)
  - All arrays properly aligned and ready for JAGS models
- `model_data_submodel.rds`: Enhanced data object including harvest-specific covariates created by prepare_harvest_data.R

### Supporting Data (Data/)
- thinnedpoints.csv: Site coordinates (longitude, latitude) for spatial diagnostics
- raw_yearly_covs.csv: Annual covariate data for harvest model preparation

## Model Files (Models/)

1. **NULL.txt**: Baseline model with no covariates
2. **DistEdge.txt**: Distance-to-edge effects model  
3. **DistEdge_RPatch.txt**: Distance-to-edge with random patch effects


## Analysis Workflow

### 1. Distance-to-Edge Analysis
Run `run_distedge_models.R` to:
- Compare functional forms (linear vs. log transformation)
- Test distance-to-edge effects on colonization/persistence
- Include random patch effects
- Generate model comparison using WAIC

### 2. Harvest Data Preparation
Run `prepare_harvest_data.R` to:
- Create harvest-specific covariates using sub-model structure
- Generate distance and age effects for sites within 565m of harvest
- Save enhanced data object as `model_data_submodel.rds`

**Key feature:** Hierarchical covariate structure where distance effects use all sites but age effects only use sites within 565m of harvest areas.

### 3. Harvest Regeneration Analysis
Run `run_harvest_models.R` to:
- Test harvest regeneration hypothesis
- Evaluate distance × age interaction effects
- Compare models using WAIC
- Generate parameter estimates and predictive checks

**Requires:** `model_data_submodel.rds` from Step 2

### 4. Model Diagnostics (Optional)
Run `model_diagnostics.R` to:
- Perform Bayesian posterior predictive checks
- Assess spatial autocorrelation (Moran's I)
- Generate diagnostic plots and summary tables

## Model Structure

### Dynamic Occupancy Framework
- **Initial occupancy (ψ):** Function of percent white spruce and stand age
- **Detection (p):** Includes Julian day and time of day effects
- **Persistence (φ):** Probability occupied site remains occupied
- **Colonization (γ):** Probability unoccupied site becomes occupied
- **Random effects:** Patch-level effects account for spatial clustering

### Key Features
- Sites thinned to minimum 400m separation to reduce spatial dependence
- All continuous covariates standardized to mean=0, SD=1
- Distance variables tested with both linear and log transformations
- 40,000 MCMC iterations with 20,000 burn-in across 3 chains
- Model selection using WAIC

## Usage Notes

1. **JAGS installation:** Ensure JAGS is properly installed and accessible to R
2. **Results directory:** Will be created automatically by scripts if it doesn't exist
4. **Runtime:** Each model takes 10-80 minutes depending on system specifications

## GitHub Repository

For the most recent code updates and version history, see:
https://github.com/harttaylor/BTNW_DynamicOccupancy.git

## Citation

Data: [This Dryad repository DOI]
Manuscript: [Journal citation]

## Contact

tah@ualberta.ca

