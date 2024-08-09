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
