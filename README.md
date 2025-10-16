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

All data files are located in the `Data/` directory.

### Core Model Data Files (RDS format)

**`model_data.rds`** - Base data structure for dynamic occupancy models. This R data object contains:
- **y**: Detection history array [114 sites × 25 years × 4 visits]. Values are 1 (species detected), 0 (not detected), or NA (site not surveyed). Years span 1993-2018 excluding 2004.
- **nsite**: Number of study sites (114)
- **nyear**: Number of survey years (25)
- **nsurv**: Array [114 × 25 × 4] indicating which visits were conducted at each site-year
- **J**: Matrix [114 × 25] indicating number of visits conducted per site-year (ranges 1-4)
- **ind**: Matrix [114 × 25] detection indicator (1 = species detected at least once that year, 0 = never detected)
- **x_psi**: Matrix [114 × 4] of covariates for initial occupancy in 1993, including intercept, standardized percent white spruce, percent white spruce squared, and standardized stand age
- **nbeta_psi**: Number of initial occupancy parameters (4)
- **x_p**: Array [114 × 25 × 4 × 3] of detection covariates (intercept, standardized Julian day, standardized time of day)
- **nbeta_p**: Number of detection parameters (3)
- **patch**: Vector [114] indicating which of 33 spatial patches each site belongs to, for random effects
- **npatch**: Number of spatial patches (33)
- **sites**: Character vector [114] of site identifiers
- **years**: Integer vector [25] of survey years (0-24, corresponding to 1993-2018 excluding 2004)

**`distance_model_arrays.rds`** - Distance-to-edge covariates formatted for colonization and persistence models. This R data object contains:
- **x_phi** and **x_gamma**: Arrays [114 × 24 × 5] of covariates affecting persistence (phi) and colonization (gamma), respectively. The 5 covariates are: intercept, standardized distance to nearest seismic line, standardized distance to nearest road, standardized distance to nearest pipeline, and standardized distance to nearest harvest area. All distances measured in meters, standardized to mean=0, SD=1.
- **Linear vs. log variants**: Contains both linear and log-transformed versions of distance covariates, with and without intercepts, to test different functional forms
- **first_year_covariates**: Matrix of site characteristics in 1993 for initial occupancy modeling
- **sites**: Site identifiers matching model_data.rds
- **years**: Year indices matching model_data.rds

**`model_data_submodel.rds`** - Extended data object including harvest regeneration covariates. Contains all elements from `model_data.rds` plus:
- **x_phi_harv_dist_submodel**: Array [114 × 25 × 1] of standardized distance to harvest edges for persistence
- **x_gamma_harv_dist_submodel**: Array [114 × 25 × 1] of standardized distance to harvest edges for colonization
- **x_phi_harv_age_submodel**: Array [114 × 25 × 1] of standardized harvest age (years since harvest). Values are actual standardized ages for sites within 565m of harvest, 0 for sites beyond 565m (indicating no age effect)
- **x_gamma_harv_age_submodel**: Array [114 × 25 × 1] of standardized harvest age for colonization (same structure as phi)
- **x_phi_harv_main_submodel**: Array [114 × 25 × 2] with both distance and age effects for persistence
- **x_gamma_harv_main_submodel**: Array [114 × 25 × 2] with both distance and age effects for colonization
- **x_phi_harv_interact_submodel**: Array [114 × 25 × 3] with distance, age, and distance×age interaction for persistence
- **x_gamma_harv_interact_submodel**: Array [114 × 25 × 3] with distance, age, and distance×age interaction for colonization

### Supporting Data Files (CSV format)

**`thinnedpoints.csv`** - Spatial coordinates for the 114 study sites. Contains:
- **SS**: Site identifier (matches site IDs in RDS files)
- **longitude**: Site longitude in decimal degrees
- **latitude**: Site latitude in decimal degrees

**`harvest_yearly_covs.csv`** - Annual site-level covariate data used for harvest model preparation. Contains 2,850 rows (114 sites × 25 years). Key variables:
- **SS**: Site identifier
- **YEAR**: Survey year (1993-2018)
- **NEAR.DIST.harvest**: Distance to nearest harvest edge (meters, 0-1000)
- **MEANAGE.565.harvest**: Mean age of harvested areas within 565m radius (years since harvest, 0-25)
- Distance variables: Variables following pattern `NEAR.DIST.*` represent distances to different disturbance types (seismic lines, roads, pipelines) in meters
- Standardized variables: Variables ending in `_std` are covariates standardized to mean=0, SD=1 for use in models

### Model Files (Models/)

1. **NULL.txt**: Baseline model with no covariates
2. **DistEdge.txt**: Distance-to-edge effects model  
3. **DistEdge_RPatch.txt**: Distance-to-edge with random patch effects

## Analysis Workflow

The analysis follows four sequential steps: 

### 1. Distance-to-Edge Analysis
Run `run_distedge_models.R` to:
- Evaluate how distance to anthropogenic edges (seismic lines, roads, pipelines, harvest areas) affects colonization and persistence probabilities
- Fits five competing models: null, log distance only, linear distance only, log distance with random patch effects, linear distance with random patch effects
- Tests whether edge effects follow linear or logarithmic patterns with distance
- Incorporates patch-level random effects to account for spatial clustering
- Saves fitted model objects to `Results/` directory

**Required inputs:** 
- `Data/model_data.rds`
- `Data/distance_model_arrays.rds`
- Model specification files in `Models/` directory

### 2. Harvest Data Preparation
Run `prepare_harvest_data.R` to:
- Create covariate structures for testing harvest regeneration hypothesis using a hierarchical sub-model approach
- Create distance effects for all sites (full ecological gradient)
- Create age effects only for sites within 565m of harvest (0 elsewhere)
- Format interaction terms (distance × age)
- Save enhanced data object as `model_data_submodel.rds`

**Key feature:** Hierarchical covariate structure where distance effects use all sites but age effects only use sites within 565m of harvest areas.

**Required inputs:**
- `Data/model_data.rds`
- `Data/harvest_yearly_covs.csv`

### 3. Harvest Regeneration Analysis
Run `run_harvest_models.R` to:
- Test whether negative effects of harvest edges diminish as forests regenerate over time
- Fit five competing models testing harvest effects: null, distance only, age only, distance + age (additive), distance × age (interaction)
- Uses the same patch random effects structure as distance-to-edge models

**Required inputs:**
- `Data/model_data_submodel.rds` (created in Step 2)
- Model specification files in `Models/` directory

### 4. Model Validation and Comparison (Optional)
Run `model_diagnostics.R` to:
- Compare all fitted models using information criteria and validate top-performing models using multiple diagnostics
- Calculates WAIC (Watanabe-Akaike Information Criterion) for all distance-to-edge and harvest models to identify top models 
- Performs detailed validation on top models:
   - Bayesian posterior predictive checks 
   - Spatial autocorrelation analysis (Moran's I)

**Required inputs:**
- All model objects from Steps 1 and 3
- `Data/thinnedpoints.csv` (for spatial autocorrelation analysis)


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

## Data Dictionary

### Detection Data Array (y)
- **Dimensions:** 114 sites × 25 years × 4 visits
- **Values:** 1 = detected, 0 = not detected, NA = not surveyed
- **Years:** 1993-2018 (excluding 2004 = 25 total years)

### Site-Level Covariates (x_psi)
| Variable | Description | Units | Range |
|----------|-------------|-------|-------|
| SS | Site identifier | - | Unique codes |
| percentconifer_std | Standardized white spruce cover | SD units | -2 to 2 |
| standage_std | Standardized forest stand age | SD units | -2 to 2 |
| patch | Patch assignment for random effects | - | 1-33 |

### Yearly Covariates (x_phi, x_gamma)
| Variable | Description | Units | Range |
|----------|-------------|-------|-------|
| NEAR.DIST.seismic | Distance to nearest seismic line | meters | 0-1000 |
| NEAR.DIST.road | Distance to nearest road | meters | 14-1000 |
| NEAR.DIST.pipeline | Distance to nearest pipeline | meters | 45-1000 |
| NEAR.DIST.harvest | Distance to nearest harvest edge | meters | 0-1000 |
| MEANAGE.565.harvest | Mean age of harvest within 565m | years | 0-25 |

### Detection Covariates (x_p)
| Variable | Description | Units | Range |
|----------|-------------|-------|-------|
| jcen | Standardized Julian day of survey | SD units | -3 to 3 |
| tcen | Standardized time of day for survey start | SD units | -3 to 3 |

### Derived Variables in Model Data
| Variable | Description | Calculation |
|----------|-------------|-------------|
| J | Number of visits per site-year | Count of non-NA values in detection array |
| ind | Detection indicator per site-year | 1 if species detected at least once that year, 0 if never detected |
| nsurv | Visit index array | Indicates which of the four possible visits were conducted |

## Usage Notes

1. **JAGS installation:** Ensure JAGS is properly installed and accessible to R
2. **Run scripts in order:** 
   - `run_distedge_models.R` (requires `model_data.rds` + `distance_model_arrays.rds`)
   - `prepare_harvest_data.R` (creates `model_data_submodel.rds`)
   - `run_harvest_models.R` (requires `model_data_submodel.rds`)
   - `model_diagnostics.R` (validates all fitted models)
3. **Results directory:** Will be created automatically by scripts if it doesn't exist
4. **Runtime:** Each model takes 10-80 minutes depending on system specifications
5. **Model validation:** The diagnostics script provides all WAIC comparisons, spatial autocorrelation checks, and Bayesian p-values for both distance-to-edge and harvest models
6. **Memory requirements:** Models require approximately 8-16 GB RAM for MCMC sampling with 60,000 posterior samples per model

## GitHub Repository

For the most recent code updates and version history, see:
https://github.com/harttaylor/BTNW_DynamicOccupancy.git

## Citation

Data and code: https://doi.org/10.5061/dryad.95x69p8z6
Manuscript: [DOI:10.1093/ornithapp/duaf068]

## Contact

tah@ualberta.ca

