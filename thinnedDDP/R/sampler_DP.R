#' Gibbs sampler for the DP Gaussian mixture model
#'
#' @param nrep Integer. Total number of MCMC iterations.
#' @param burnin Integer. Number of initial iterations discarded as burn-in.
#' @param thinning_factor Integer. Thinning interval used to store posterior
#'   samples. For example, `thinning_factor = 2` retains every second draw.
#' @param y Numeric vector of observations.
#' @param trunc Integer. Truncation level of the stick-breaking representation,
#'   corresponding to the maximum number of mixture components.
#' @param m0 Prior mean hyperparameter for component means.
#' @param tau0 Prior precision parameter for component means. The prior is
#'   \eqn{\mu_j \mid \sigma_j^2 \sim N(m_0, \sigma_j^2 / \tau_0)}.
#' @param gamma0 Shape hyperparameter of the inverse-Gamma prior on component
#'   variances.
#' @param lambda0 Rate hyperparameter of the inverse-Gamma prior on component
#'   variances. The prior is
#'   \eqn{1/\sigma_j^2 \sim \mathrm{Gamma}(\gamma_0, \lambda_0)}.
#' @param alpha Positive scalar. Concentration parameter of the Dirichlet
#'   process controlling the stick-breaking weights.
#' @param mu_start Numeric vector of length `trunc` containing initial values
#'   for the component means.
#' @param sigma2_start Numeric vector of length `trunc` containing initial
#'   values for the component variances.
#' @param cl_start Numeric or integer vector of initial cluster allocations.
#'   Values should be in `{0, ..., trunc - 1}`.
#' @param progressbar Logical. If `TRUE`, displays a progress bar during MCMC
#'   sampling.
#'
#' @return A list containing:
#' \describe{
#'   \item{mu}{Matrix of posterior samples for component means. Rows correspond
#'   to retained MCMC iterations and columns to mixture components.}
#'   \item{sigma2}{Matrix of posterior samples for component variances.}
#'   \item{cl}{Matrix of posterior samples for cluster allocations.}
#'   \item{pi}{Matrix of posterior samples for the mixture weights.}
#'   \item{y}{Observed data.}
#'   \item{nrep}{Total number of MCMC iterations.}
#'   \item{burnin}{Number of burn-in iterations.}
#'   \item{trunc}{Truncation level of the stick-breaking representation.}
#'   \item{hyperpar}{Vector containing the prior hyperparameters
#'   \code{(m0, tau0, gamma0, lambda0)}.}
#'   \item{alpha}{Dirichlet process concentration parameter.}
#'   \item{mu_start}{Initial values for component means.}
#'   \item{sigma2_start}{Initial values for component variances.}
#'   \item{cl_start}{Initial cluster allocations.}
#'   \item{time}{Elapsed computation time.}
#' }
#' 
#' @export

sampler_DP <- function(nrep, burnin, 
                       thinning_factor = 2,
                       y, 
                       trunc = 50, 
                       m0 = 0, tau0 = 0.1, 
                       gamma0 = 3, lambda0 = 2,
                       alpha = 1, 
                       mu_start, sigma2_start, cl_start, 
                       progressbar = TRUE)
{ 
  start = Sys.time()
  out = sampler_DP_arma(nrep, burnin, 
                        thinning_factor,
                        y, 
                        trunc, 
                        m0, tau0, 
                        gamma0, lambda0,
                        alpha, 
                        mu_start, sigma2_start, cl_start, 
                        progressbar)
  end = Sys.time()
  out$time = end-start
  return(out)
}
