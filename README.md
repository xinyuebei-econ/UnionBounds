# UnionBounds

**UnionBounds** is an R package for calculating confidence intervals of a target parameter $\theta$ whose identified set is a union bound, i.e. $\theta \in \left[\min_{b \in \mathcal{B}} \lambda_{\ell,b}, \max_{b \in \mathcal{B}} \lambda_{u,b}\right]$. It is based on the modified conditional inference proposed by [Bei (2024)](https://xinyuebei-econ.github.io/files/Bei_Xinyue_JMP.pdf).

## Overview

This package includes three main functions:
- `CI_ModifiedConditional`: Works for general union bounds defined in Bei (2024) Equation 1.
- `CI_RM`: Calculates CI of the treatment effect under relative magnitudes relaxation.
- `CI_SDRM`: Calculates CI of the treatment effect under second differences relative magnitudes relaxation.

## Installation

You can install the package directly from GitHub using the `remotes` package.

```r
# Install the remotes package if you haven't already
install.packages("remotes")

# Install UnionBounds from GitHub
remotes::install_github("xinyuebei-econ/UnionBounds")
```

## Examples

```r
# load the library
library('UnionBounds')

# view the documentation
?CI_ModifiedConditional
?CI_RM
?CI_SDRM
```

Example for CI_ModifiedConditional: $\theta \in \left[\min \left\[ \delta_1, \delta_2\right\], \max \left\[ \delta_3, \delta_2\right\]\right]$
```r
set.seed(0)

deltaSigma <- matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 1), nrow = 3, ncol = 3, byrow = TRUE)
deltahat <- MASS::mvrnorm(n = 1, c(0,0, 0), deltaSigma)
deltahat <- matrix(deltahat, nrow = 3, ncol = 1)

Al <- matrix(c(1, 0, 0, 0, 1, 0), nrow = 2, ncol = 3, byrow = TRUE)
Au <- matrix(c(0, 0, 1, 0, 1, 0), nrow = 2, ncol = 3, byrow = TRUE)

CI <- CI_ModifiedConditional(deltahat, deltaSigma, Al, Au, alpha = 0.05)
print(CI$ConfidenceInterval) # This is the confidence interval [ -2.072250  3.007971]
```

Example for CI_RM
```r
# the sample data is from Dustmann et al (2022), see Bei(2024) Section 5
data(example_data, package = "UnionBounds")

betahat <- example_data$betahat                     # Diff in Diff estimator 
betaSigma <- example_data$betaSigma                 # covariance matrix for Diff in Diff estimator
prePeriodIndices <- example_data$prePeriodIndices   # number of prepolicy periods
ell_post <- c(1, 0)                                 # the parameter of interest is the treatment effect at time 1

CI <- CI_RM(betahat, betaSigma, prePeriodIndices, ell_post, M = 1, alpha = 0.05)
print(CI$ConfidenceInterval)                        # This is the RM confidence interval
```

Example for CI_SDRM
```r
# the sample data is from Dustmann et al (2022), see Bei(2024) Section 5
data(example_data, package = "UnionBounds")

betahat <- example_data$betahat                     # Diff in Diff estimator 
betaSigma <- example_data$betaSigma                 # covariance matrix for Diff in Diff estimator
prePeriodIndices <- example_data$prePeriodIndices   # number of prepolicy periods
ell_post <- c(1, 0)                                 # the parameter of interest is the treatment effect at time 1

CI <- CI_SDRM(betahat, betaSigma, prePeriodIndices, ell_post, M = 1, alpha = 0.05)
print(CI$ConfidenceInterval)                        # This is the SDRM confidence interval
```
