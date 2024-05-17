\name{msNMix}
\alias{msNMix}
\title{Function for Fitting Multi-species N-mixture Models}

\description{
Function for fitting multi-species N-mixture models.
}

\usage{
msNMix(abund.formula, det.formula, data, inits, priors,
       tuning, n.batch, batch.length, accept.rate = 0.43,
       family = 'Poisson', n.omp.threads = 1, verbose = TRUE, n.report = 100,
       n.burn = round(.10 * n.batch * batch.length), n.thin = 1,
       n.chains = 1, ...)
}

\arguments{

\item{abund.formula}{a symbolic description of the model to be fit
  for the abundance portion of the model using R's model syntax. Only
  right-hand side of formula is specified. See example below. Random intercepts
  and slopes are allowed using lme4 syntax (Bates et al. 2015).}

\item{det.formula}{a symbolic description of the model to be fit
  for the detection portion of the model using R's model syntax. Only
  right-hand side of formula is specified. See example below. Random intercepts
  and slopes are allowed using lme4 syntax (Bates et al. 2015).}

\item{data}{a list containing data necessary for model fitting.
  Valid tags are \code{y}, \code{abund.covs}, \code{det.covs}, and \code{offset}. \code{y}
  is a three-dimensional array of observed count data with first dimension equal to the number
  of species, second dimension equal to the number of sites, and third
  dimension equal to the maximum number of replicates at a given site.
  \code{abund.covs} is a matrix or data frame containing the variables
  used in the abundance portion of the model, with \eqn{J}{J} rows for
  each column (variable). \code{det.covs} is a list of variables included
  in the detection portion of the model. Each list element is a different
  detection covariate, which can be site-level or observational-level.
  Site-level covariates are specified as a vector of length \eqn{J}{J} while
  observation-level covariates are specified as a matrix or data frame
  with the number of rows equal to \eqn{J}{J} and number of columns equal
  to the maximum number of replicates at a given site. \code{offset}
  is an offset to use in the abundance model (e.g., an area offset). This can be either a single
  value or a vector with an offset for each site (e.g., if survey area differed in size).}

\item{inits}{a list with each tag corresponding to a parameter name.
  Valid tags are \code{alpha.comm}, \code{beta.comm}, \code{beta},
  \code{alpha}, \code{tau.sq.beta}, \code{tau.sq.alpha}, \code{sigma.sq.mu},
  \code{sigma.sq.p}, \code{kappa}, and \code{N}. \code{sigma.sq.mu} and
  \code{sigma.sq.p} are only relevant when including random effects in the
  abundance and detection portion of the model, respectively. \code{kappa} is
  only relevant when \code{family = 'NB'}.
  The value portion of each tag is the parameter's initial value. See
  \code{priors} description for definition of each parameter name.
  Additionally, the tag \code{fix} can be set to \code{TRUE}
  to fix the starting values across all chains. If \code{fix} is not specified
  (the default), starting values are varied randomly across chains.}

\item{priors}{a list with each tag corresponding to a parameter name.
  Valid tags are \code{beta.comm.normal}, \code{alpha.comm.normal},
  \code{tau.sq.beta.ig}, \code{tau.sq.alpha.ig}, \code{sigma.sq.mu.ig},
  \code{sigma.sq.p.ig}, and \code{kappa.unif}. Community-level abundance
  (\code{beta.comm}) and detection (\code{alpha.comm}) regression coefficients
  are assumed to follow a normal distribution. The hyperparameters of the
  normal distribution are passed as a list of length two with the first
  and second elements corresponding to the mean and variance of the normal
  distribution, which are each specified as vectors of
  length equal to the number of coefficients to be estimated or of length one
  if priors are the same for all coefficients. If not
  specified, prior means are set to 0 and prior variances for the abundance
  coefficients are set to 100 and for the detection coefficients are set to 2.72.
  Community-level variance parameters for abundance (\code{tau.sq.beta}) and
  detection (\code{tau.sq.alpha}) are assumed to follow an inverse Gamma
  distribution. The hyperparameters of the inverse gamma distribution are
  passed as a list of length two with the first and second elements
  corresponding to the shape and scale parameters, which are each specified
  as vectors of length equal to the number of coefficients to be estimated
  or a single value if all parameters are assigned the same prior.
  If not specified, prior shape and scale parameters are set to 0.1. \code{sigma.sq.mu} and
  \code{sigma.sq.p} are the random effect variances for any abundance or
  detection random effects, respectively, and are assumed to follow an inverse
  Gamma distribution. The hyperparameters of the inverse-Gamma distribution
  are passed as a list of length two with first and second elements corresponding
  to the shape and scale parameters, respectively, which are each specified as
  vectors of length equal to the number of random intercepts/slopes or of length one
  if priors are the same for all random effect variances. \code{kappa} is the
  negative binomial dispersion parameter for each species and is assumed to
  follow a uniform distribution. The hyperparameters of the uniform distribution
  are passed as a list of length two with first and second elements corresponding to the
  lower and upper bounds of the uniform distribution, respectively, which are each
  specified as vectors of length equal to the number of species or of length one
  if priors are the same for all species-specific dispersion parameters.}

\item{tuning}{a list with each tag corresponding to a parameter name,
  whose value defines the initial variance of the adaptive sampler.
  Valid tags are \code{beta}, \code{alpha}, \code{beta.star} (the abundance
  random effect values), \code{alpha.star} (the detection random effect values), and
  \code{kappa}. See Roberts and Rosenthal (2009) for details.}

\item{n.batch}{the number of MCMC batches in each chain to run for the Adaptive MCMC
  sampler. See Roberts and Rosenthal (2009) for details.}

\item{batch.length}{the length of each MCMC batch in each chain to run for the Adaptive
  MCMC sampler. See Roberts and Rosenthal (2009) for details.}

\item{accept.rate}{target acceptance rate for Adaptive MCMC. Default is
  0.43. See Roberts and Rosenthal (2009) for details.}

\item{family}{the distribution to use for the latent abundance process. Currently
  supports \code{'NB'} (negative binomial) and \code{'Poisson'}.}

\item{n.omp.threads}{a positive integer indicating the number of threads
  to use for SMP parallel processing. The package must be compiled for
  OpenMP support. For most Intel-based machines, we recommend setting
  \code{n.omp.threads} up to the number of hypterthreaded cores. Note,
  \code{n.omp.threads} > 1 might not work on some systems. Currently only
  relevant for spatial models.}

\item{verbose}{if \code{TRUE}, messages about data preparation,
  model specification, and progress of the sampler are printed to the screen.
  Otherwise, no messages are printed.}

\item{n.report}{the interval to report MCMC progress.}

\item{n.burn}{the number of samples out of the total \code{n.samples} to
  discard as burn-in for each chain. By default, the first 10\% of samples is discarded.}

\item{n.thin}{the thinning interval for collection of MCMC samples. The
  thinning occurs after the \code{n.burn} samples are discarded. Default
  value is set to 1.}

\item{n.chains}{the number of chains to run in sequence.}

\item{...}{currently no additional arguments}

}

\references{
  Bates, Douglas, Martin Maechler, Ben Bolker, Steve Walker (2015).
  Fitting Linear Mixed-Effects Models Using lme4. Journal of
  Statistical Software, 67(1), 1-48. \doi{10.18637/jss.v067.i01}.

  Royle, J. A. (2004). N‐mixture models for estimating population size
  from spatially replicated counts. Biometrics, 60(1), 108-115.

  Yamaura, Y., Royle, J. A., Shimada, N., Asanuma, S., Sato, T.,
  Taki, H., & Makino, S. I. (2012). Biodiversity of man-made open habitats
  in an underused country: a class of multispecies abundance models for
  count data. Biodiversity and Conservation, 21(6), 1365-1380.
}

\author{
  Jeffrey W. Doser \email{doserjef@msu.edu}, \cr
}

\value{
  An object of class \code{msNMix} that is a list comprised of:

  \item{beta.comm.samples}{a \code{coda} object of posterior samples
    for the community level abundance regression coefficients.}

  \item{alpha.comm.samples}{a \code{coda} object of posterior samples
    for the community level detection regression coefficients.}

  \item{tau.sq.beta.samples}{a \code{coda} object of posterior samples
    for the abundance community variance parameters.}

  \item{tau.sq.alpha.samples}{a \code{coda} object of posterior samples
    for the detection community variance parameters.}

  \item{beta.samples}{a \code{coda} object of posterior samples
    for the species level abundance regression coefficients.}

  \item{alpha.samples}{a \code{coda} object of posterior samples
    for the species level detection regression coefficients.}

  \item{kappa.samples}{a \code{coda} object of posterior samples
    for the species level abundance dispersion parameters. Only included
    when \code{family = 'NB'}.}

  \item{N.samples}{a three-dimensional array of posterior samples for
    the latent abundance values for each species.}

  \item{mu.samples}{a three-dimensional array of posterior samples for
    the latent expected abundance values for each species.}

  \item{sigma.sq.mu.samples}{a \code{coda} object of posterior samples
    for variances of random effects included in the abundance portion
    of the model. Only included if random effects are specified in
    \code{abund.formula}.}

  \item{sigma.sq.p.samples}{a \code{coda} object of posterior samples
    for variances of random effects included in the detection portion
    of the model. Only included if random effects are specified in
    \code{det.formula}.}

  \item{beta.star.samples}{a \code{coda} object of posterior samples
    for the abundance random effects. Only included if random effects
    are specified in \code{abund.formula}.}

  \item{alpha.star.samples}{a \code{coda} object of posterior samples
    for the detection random effects. Only included if random effects
    are specified in \code{det.formula}.}

  \item{rhat}{a list of Gelman-Rubin diagnostic values for some of the model
    parameters.}

  \item{ESS}{a list of effective sample sizes for some of the model parameters.}

  \item{run.time}{MCMC sampler execution time reported using \code{proc.time()}.}

  The return object will include additional objects used for
  subsequent prediction and/or model fit evaluation. Note that detection probability
  estimated values are not included in the model object, but can be extracted
  using \code{fitted()}.
}

\examples{
set.seed(408)
J.x <- 8
J.y <- 8
J <- J.x * J.y
n.rep <- sample(5, size = J, replace = TRUE)
n.sp <- 6
# Community-level covariate effects
# Abundance
beta.mean <- c(0, 0.5)
p.abund <- length(beta.mean)
tau.sq.beta <- c(0.2, 1.2)
# Detection
alpha.mean <- c(0, 0.5, 0.8)
tau.sq.alpha <- c(0.2, 1, 1.5)
p.det <- length(alpha.mean)
# Random effects
mu.RE <- list()
p.RE <- list()
# Draw species-level effects from community means.
beta <- matrix(NA, nrow = n.sp, ncol = p.abund)
alpha <- matrix(NA, nrow = n.sp, ncol = p.det)
for (i in 1:p.abund) {
  beta[, i] <- rnorm(n.sp, beta.mean[i], sqrt(tau.sq.beta[i]))
}
for (i in 1:p.det) {
  alpha[, i] <- rnorm(n.sp, alpha.mean[i], sqrt(tau.sq.alpha[i]))
}

dat <- simMsNMix(J.x = J.x, J.y = J.y, n.rep = n.rep, n.sp = n.sp, beta = beta, alpha = alpha,
                 mu.RE = mu.RE, p.RE = p.RE, sp = FALSE, family = 'Poisson')

y <- dat$y
X <- dat$X
X.p <- dat$X.p
X.re <- dat$X.re
X.p.re <- dat$X.p.re

# Package all data into a list
abund.covs <- X
colnames(abund.covs) <- c('int', 'abund.cov.1')
det.covs <- list(det.cov.1 = as.data.frame(X.p[, , 2]),
                 det.cov.2 = as.data.frame(X.p[, , 3]))
data.list <- list(y = y,
                  abund.covs = abund.covs,
                  det.covs = det.covs)
prior.list <- list(beta.comm.normal = list(mean = rep(0, p.abund),
                                           var = rep(100, p.abund)),
                   alpha.comm.normal = list(mean = rep(0, p.det),
                                            var = rep(2.72, p.det)),
                   tau.sq.beta.ig = list(a = 0.1, b = 0.1),
                   tau.sq.alpha.ig = list(a = 0.1, b = 0.1))
inits.list <- list(beta.comm = 0, alpha.comm = 0,
                   beta = 0, alpha = 0,
                   tau.sq.beta = 0.5, tau.sq.alpha = 0.5,
                   N = apply(y, c(1, 2), max, na.rm = TRUE))
tuning.list <- list(beta = 0.5, alpha = 0.5)

n.batch <- 4
batch.length <- 25
n.burn <- 0
n.thin <- 1
n.chains <- 1

out <- msNMix(abund.formula = ~ abund.cov.1,
              det.formula = ~ det.cov.1 + det.cov.2,
              data = data.list,
              n.batch = n.batch,
              inits = inits.list,
              priors = prior.list,
              tuning = tuning.list,
              batch.length = batch.length,
              n.omp.threads = 1,
              verbose = TRUE,
              n.report = 1,
              n.burn = n.burn,
              n.thin = n.thin,
              n.chains = n.chains)

summary(out, level = 'community')
}
