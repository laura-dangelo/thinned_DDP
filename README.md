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

## Required Inputs

#### `y`
Numeric vector containing the observed data.

Requirements:
* numeric vector;
* length `N`;
* no missing values (`NA`).

---

#### `group` 
Vector indicating the group membership of each observation.

Requirements:
* same length as `y`;
* one group label for each observation;
* groups must be coded as consecutive integers starting from `0`: `(0, 1, 2, ..., G - 1)`, where `G` is the number of groups.

If your groups are stored as factors:
```r
group <- as.numeric(factor(group)) - 1
```

---

#### `nrep`
Total number of MCMC iterations (positive integer). 

---

#### `burnin`
Number of initial iterations discarded as burn-in (must satisfy `burnin < nrep`).

---

#### `mu_start`, `sigma2_start`
Initial values of the cluster means and variances. Numeric vectors of length `trunc` (truncation level of the DP - default is `trunc = 50`). For `sigma2_start`, all entries must be positive.

---

#### `cl_start`
Initial cluster allocation for each observation. Numeric vector of length equal to `length(y)`. Cluster labels must be integers in `0, 1, ..., trunc - 1`.

---

## Optional Inputs

### `thinning_factor`

Controls how often MCMC draws are stored.

```r
thinning_factor = 2
```

Examples:

| Value | Stored draws           |
| ----- | ---------------------- |
| 1     | Every iteration        |
| 2     | Every second iteration |
| 5     | Every fifth iteration  |

Requirements:

```r
thinning_factor >= 1
```

---

### `trunc`

Truncation level of the stick-breaking representation.

```r
trunc = 50
```

Requirements:

* positive integer;
* should be large enough so that posterior mass does not accumulate at the largest cluster index.

---

### Prior Hyperparameters

#### Mean prior

The model assumes

[
\mu_j \mid \sigma_j^2 \sim N\left(m_0,\frac{\sigma_j^2}{\tau_0}\right)
]

Parameters:

```r
m0
tau0
```

Defaults:

```r
m0 = 0
tau0 = 0.1
```

Requirement:

```r
tau0 > 0
```

---

#### Variance prior

The model assumes

[
1/\sigma_j^2 \sim \text{Gamma}(\gamma_0,\lambda_0)
]

Parameters:

```r
gamma0
lambda0
```

Defaults:

```r
gamma0 = 3
lambda0 = 2
```

Requirements:

```r
gamma0 > 0
lambda0 > 0
```

---

#### DP concentration parameter

```r
alpha = 1
```

Requirement:

```r
alpha > 0
```

Larger values generally encourage a larger number of occupied clusters.

---

#### Thinning probability prior

For each group-specific thinning probability,

[
p_g \sim \text{Beta}(a_\beta,b_\beta)
]

Parameters:

```r
a_beta
b_beta
```

Defaults:

```r
a_beta = 1
b_beta = 1
```

Requirements:

```r
a_beta > 0
b_beta > 0
```

---

### `progressbar`

Logical value controlling whether a progress bar is displayed.

```r
progressbar = TRUE
```

---

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
