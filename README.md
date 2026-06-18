# Dependent Dirichlet processes via thinning


[![arXiv:2506.18223](https://img.shields.io/badge/arXiv-2506.18223-b31b1b.svg)](https://arxiv.org/abs/2506.18223)
![Last Commit](https://img.shields.io/github/last-commit/laura-dangelo/thinned_DDP)

This repository contains the code to replicate the results in the paper "_Dependent Dirichlet processes via thinning_" by D'Angelo, Nipoti, and Ongaro (available on [arXiv](https://arxiv.org/abs/2506.18223)). It also contains the R package `thinnedDDP`, which implements the Gibbs sampler for posterior inference (see [Installation](#installation)).


## Structure of the Repository

Main structure:  
- `01_Simulation_study` : folder that contains the scripts and outputs of the simulation study;  
- `02_Application` : folder that contains the script and outputs of the real data analysis;   
- `03_Sensitivity_study` : folder that contains the scripts and outputs of the sensitivity study;  
- `thinnedDDP` : R package implementing the Gibbs sampler.   
- `thinnedDDP_project` : R project to be used as main directory for running the code.  



Large files (typically, outputs of the Gibbs sampler algorithm and extracted quantities) are excluded from the repository and are available in the Google Drive [folder](https://drive.google.com/drive/folders/1VVzDpo1fLS5QG8ByCHO6oXbZSlVjNahu?usp=sharing). The name of each subfolder corresponds to the path of the files in the repository. Note that these files are not necessary to produce the plots in the article, which can be obtained starting from the RDS files available in the folders named `output_RDS`.  


## Installation
You can install the package `thinned_DDP` for running the Gibbs sampler from GitHub with:
``` r
# install.packages("devtools")
devtools::install_github("laura-dangelo/thinned_DDP", subdir='thinnedDDP')
```

## Usage
The function `sampler_thinnedDDP()` implements a blocked Gibbs sampler for a thinned Dependent Dirichlet Process (thinned-DDP) mixture model for grouped data.
```r
sampler_thinnedDDP(
  nrep,
  burnin,
  thinning_factor = 2,
  y,
  group,
  trunc = 50,
  m0 = 0,
  tau0 = 0.1,
  gamma0 = 3,
  lambda0 = 2,
  alpha = 1,
  a_beta = 1,
  b_beta = 1,
  mu_start,
  sigma2_start,
  cl_start,
  progressbar = TRUE
)
```

---

### Required Inputs

| Parameter      | Description                                        | Requirements                                                                       |
| -------------- | -------------------------------------------------- | ---------------------------------------------------------------------------------- |
| `y`            | Vector of observations.                            | Numeric vector of length `N`; no missing values.                                   |
| `group`        | Group membership associated with each observation. | Same length as `y`; groups must be coded as consecutive integers `0, 1, ..., G-1`. |
| `nrep`         | Total number of MCMC iterations.                   | Integer > 0 and `nrep > burnin`.                                                   |
| `burnin`       | Number of initial iterations discarded as burn-in. | Integer ≥ 0 and `burnin < nrep`.                                                   |
| `mu_start`     | Initial values of cluster means.                   | Numeric vector of length `trunc`.                                                  |
| `sigma2_start` | Initial values of cluster variances.               | Positive numeric vector of length `trunc`.                                         |
| `cl_start`     | Initial cluster allocation for each observation.   | Integer vector of length `N`; values must belong to `{0, ..., trunc-1}`.           |

#### Input Dimensions

The following conditions must hold:

```r
length(y) == length(group)
length(cl_start) == length(y)

length(mu_start) == trunc
length(sigma2_start) == trunc
```

#### Group Coding

Group labels must be coded as consecutive integers starting from zero: e.g.,

```r
group <- c(0, 0, 0, 1, 1, 2, 2)
```

If the groups are stored as factors or character labels, convert them before calling the sampler:

```r
group <- as.numeric(factor(group)) - 1
```

#### Cluster Initialization

Initial cluster allocations must be integers between `0` and `trunc - 1`:

```r
cl_start <- sample(
  0:(trunc - 1),
  length(y),
  replace = TRUE
)
```


### Optional Parameters

| Parameter         | Default | Description                                                                                      | Constraints       |
| ----------------- | ------- | ------------------------------------------------------------------------------------------------ | ----------------- |
| `thinning_factor` | `2`     | Store one draw every `thinning_factor` MCMC iterations.                                          | Integer ≥ 1       |
| `trunc`           | `50`    | Truncation level of the stick-breaking representation (maximum number of mixture components).    | Integer > 0       |
| `m0`              | `0`     | Prior mean for cluster means.                                                                    | —                 |
| `tau0`            | `0.1`   | Prior precision parameter in (\mu_j \mid \sigma_j^2 \sim N(m_0,\sigma_j^2/\tau_0)).              | `tau0 > 0`        |
| `gamma0`          | `3`     | Shape parameter of the inverse-Gamma prior on cluster variances.                                 | `gamma0 > 0`      |
| `lambda0`         | `2`     | Rate parameter of the inverse-Gamma prior on cluster variances.                                  | `lambda0 > 0`     |
| `alpha`           | `1`     | Dirichlet Process concentration parameter. Larger values typically favor more occupied clusters. | `alpha > 0`       |
| `a_beta`          | `1`     | First shape parameter of the Beta prior on group-specific thinning probabilities.                | `a_beta > 0`      |
| `b_beta`          | `1`     | Second shape parameter of the Beta prior on group-specific thinning probabilities.               | `b_beta > 0`      |
| `progressbar`     | `TRUE`  | Display a progress bar during MCMC execution.                                                    | `TRUE` or `FALSE` |

### Prior Specification

The sampler uses the following prior distributions:

[
\mu_j \mid \sigma_j^2 \sim N\left(m_0,\frac{\sigma_j^2}{\tau_0}\right)
]

[
1/\sigma_j^2 \sim \mathrm{Gamma}(\gamma_0,\lambda_0)
]

[
p_g \sim \mathrm{Beta}(a_\beta,b_\beta)
]

where (p_g) denotes the group-specific thinning probability.

## Example

```r
set.seed(123)

N <- 200
G <- 3

group <- rep(0:(G - 1), length.out = N)

y <- c(
  rnorm(70, -2, 1),
  rnorm(70, 0, 1),
  rnorm(60, 3, 1)
)

trunc <- 20

fit <- sampler_thinnedDDP(
  nrep = 5000,
  burnin = 1000,
  thinning_factor = 5,
  y = y,
  group = group,
  trunc = trunc,
  mu_start = rnorm(trunc),
  sigma2_start = rep(1, trunc),
  cl_start = sample(0:(trunc - 1), N, replace = TRUE)
)
```

---

## Output

The function returns a list containing:

| Object          | Description                                              |
| --------------- | -------------------------------------------------------- |
| `mu`            | Posterior draws of cluster means                         |
| `sigma2`        | Posterior draws of cluster variances                     |
| `cl`            | Posterior draws of cluster allocations                   |
| `pi`            | Group-specific mixture weights                           |
| `ell`           | Group-specific thinning indicators                       |
| `thinning_prob` | Posterior draws of group-specific thinning probabilities |
| `y`             | Original data vector                                     |
| `group`         | Group membership vector                                  |
| `time`          | Total execution time                                     |
| `hyperpar`      | Prior hyperparameters                                    |
| `alpha`         | DP concentration parameter                               |

---

## Important Notes

1. Group labels must be coded as

```r
0, 1, ..., G - 1
```

with no gaps.

2. Cluster labels in `cl_start` must be coded as

```r
0, 1, ..., trunc - 1
```

3. Input dimensions must satisfy:

```r
length(y) == length(group)
length(cl_start) == length(y)

length(mu_start) == trunc
length(sigma2_start) == trunc
```

4. The number of saved posterior samples is approximately

```r
(nrep - burnin) / thinning_factor
```

depending on the burn-in and thinning configuration.
