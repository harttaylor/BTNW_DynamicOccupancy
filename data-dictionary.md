# Data Dictionary

## Detection Data Array (y)
- **Dimensions:** 114 sites × 25 years × 4 visits
- **Values:** 1 = detected, 0 = not detected, NA = not surveyed
- **Years:** 1993-2018 (excluding 2004)

## Site-Level Covariates (x_psi)
| Variable | Description | Units | Range |
|----------|-------------|-------|-------|
| SS | Site identifier | - | Unique codes |
| percentconifer_std | Standardized white spruce cover | SD units | -2 to 2 |
| standage_std | Standardized forest stand age | SD units | -2 to 2 |
| patch | Patch assignment for random effects | - | 1-38 |

## Time-Varying Covariates (x_phi, x_gamma)
| Variable | Description | Units | Range |
|----------|-------------|-------|-------|
| NEAR.DIST.seismic | Distance to nearest seismic line | meters | 0-1000 |
| NEAR.DIST.road | Distance to nearest road | meters | 14-1000 |
| NEAR.DIST.pipeline | Distance to nearest pipeline | meters | 45-1000 |
| NEAR.DIST.harvest | Distance to nearest harvest edge | meters | 0-1000 |
| MEANAGE.565.harvest | Mean age of harvest within 565m | years | 0-25 |

## Detection Covariates (x_p)
| Variable | Description | Units | Range |
|----------|-------------|-------|-------|
| jcen | Standardized Julian day | SD units | -3 to 3 |
| tcen | Standardized time of day | SD units | -3 to 3 |

## Derived Variables
| Variable | Description | Calculation |
|----------|-------------|-------------|
| J | Number of visits per site-year | Count of non-NA detections |
| ind | Detection indicator | 1 if detected at least once, 0 otherwise |
| nsurv | Visit indices | Which visits were conducted |