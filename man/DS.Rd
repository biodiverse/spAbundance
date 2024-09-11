\name{DS}
\alias{DS}
\title{Function for Fitting Single-Species Hierarchical Distance Sampling Models}

\usage{
DS(abund.formula, det.formula, data, inits, priors, tuning,
   n.batch, batch.length, accept.rate = 0.43, family = 'Poisson',
   transect = 'line', det.func = 'halfnormal',
   n.omp.threads = 1, verbose = TRUE,
   n.report = 100, n.burn = round(.10 * n.batch * batch.length), n.thin = 1,
   n.chains = 1, ...)
}

\description{
  Function for fitting single-sepcies hierarchical distance sampling models
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
  Valid tags are \code{y}, \code{covs}, \code{dist.breaks}, and \code{offset}. \code{y}
  is a matrix or data frame of the observed count values,
  with first dimension equal to the number of
  sites (\eqn{J}{J}) and second dimension equal to the number of
  distance bins. \code{covs} is a matrix or data frame
  containing the variables used in the the abundance and/or the detection
  portion of the model, with
  \eqn{J}{J} rows for each column (variable). \code{dist.breaks} is a vector of
  distances that denote the breakpoints of the distance bands. \code{dist.breaks} should
  have length equal to the number of columns in \code{y} plus one. \code{offset} is an
  offset that can be used to scale estimates from abundance per transect to density per
  some desired unit of measure. This can be either a single value or a vector with an offset
  value for each site (e.g., if transects differ in length)}

\item{inits}{a list with each tag corresponding to a parameter name.
  Valid tags are \code{N}, \code{beta}, \code{alpha}, \code{kappa},
  \code{sigma.sq.mu}, and \code{sigma.sq.p}. The value portion of each tag is the
  parameter's initial value. \code{sigma.sq.mu} and \code{sigma.sq.p} are
  only relevant when including random effects in the abundance and
  detection portion of the distance sampling model, respectively. \code{kappa} is
  only relevant when \code{family = 'NB'}. See \code{priors}
  description for definition of each parameter name.
  Additionally, the tag \code{fix} can be set to \code{TRUE}
  to fix the starting values across all chains. If \code{fix} is not specified
  (the default), starting values are varied randomly across chains.}

\item{priors}{a list with each tag corresponding to a parameter name.
  Valid tags are \code{beta.normal}, \code{alpha.normal}, \code{kappa.unif},
  \code{sigma.sq.mu.ig}, and \code{sigma.sq.p.ig}.
  Abundance (\code{beta}) and detection (\code{alpha})
  regression coefficients are assumed to follow a normal distribution.
  The hyperparameters of the normal distribution are passed as a list of
  length two with the first and second elements corresponding to the mean
  and variance of the normal distribution, which are each specified as vectors of
  length equal to the number of coefficients to be estimated or of length one
  if priors are the same for all coefficients. If not specified, prior means
  are set to 0 and prior variances set to 100. \code{kappa} is the negative binomial
  dispersion parameter and is assumed to follow a uniform distribution. The
  hyperparameters of the uniform distribution are passed as a vector of length
  two with the first and second elements corresponding to the lower and upper
  bounds of the uniform distribution. \code{sigma.sq.mu} and
  \code{sigma.sq.p} are the random effect variances for any abundance or
  detection random effects, respectively, and are assumed to follow an inverse
  Gamma distribution. The hyperparameters of the inverse-Gamma distribution
  are passed as a list of length two with first and second elements corresponding
  to the shape and scale parameters, respectively, which are each specified as
  vectors of length equal to the number of random intercepts/slopes or of length one
  if priors are the same for all random effect variances.}

\item{tuning}{a list with each tag corresponding to a parameter name, whose value
  defines the initial tuning variance of the adaptive sampler. Valid tags include
  \code{beta}, \code{alpha}, \code{beta.star} (the abundance
  random effect values), \code{alpha.star} (the detection random effect values), and
  \code{kappa}. See Roberts and Rosenthal (2009) for details.}

\item{n.batch}{the number of MCMC batches in each chain to run for the adaptive MCMC
  sampler. See Roberts and Rosenthal (2009) for details.}

\item{batch.length}{the number of MCMC samples in each batch in each chain to run for the Adaptive
  MCMC sampler. See Roberts and Rosenthal (2009) for details.}

\item{accept.rate}{target acceptance rate for Adaptive MCMC. Default is
  0.43. See Roberts and Rosenthal (2009) for details.}

\item{family}{the distribution to use for the latent abundance process. Currently
  supports \code{'NB'} (negative binomial) and \code{'Poisson'}.}

\item{transect}{the type of transect. Currently supports line transects (\code{'line'}) or
  circular transects (i.e., point counts; \code{'point'}).}

\item{det.func}{the detection model used to describe how detection probability varies
  with distance. In other software, this is often referred to as the key function. Currently
  supports two functions: half normal (\code{'halfnormal'}) and negative exponential
  (\code{'negexp'}).}

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

  Royle, J. A., Dawson, D. K., & Bates, S. (2004). Modeling
  abundance effects in distance sampling. Ecology, 85(6), 1591-1597.
}

\author{
  Jeffrey W. Doser \email{doserjef@msu.edu}, \cr
}

\value{
  An object of class \code{DS} that is a list comprised of:

  \item{beta.samples}{a \code{coda} object of posterior samples
    for the abundance regression coefficients.}

  \item{alpha.samples}{a \code{coda} object of posterior samples
    for the detection regression coefficients.}

  \item{kappa.samples}{a \code{coda} object of posterior samples
    for the abundance dispersion parameter. Only included when
    \code{family = 'NB'}.}

  \item{N.samples}{a \code{coda} object of posterior samples
    for the latent abundance values. Note that these values always
    represent transect-level abundance, even when an offset is
    supplied.}

  \item{mu.samples}{a \code{coda} object of posterior samples
    for the latent expected abundance values. When an offset is
    supplied in the \code{data} object, these correspond to expected
    abundance per unit area (i.e., density).}

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

  \item{y.rep.samples}{a three-dimensional array of fitted values.
    Array dimensions correspond to MCMC samples, sites, and distance band.}

  \item{pi.samples}{a three-dimensional array of cell-specific detection
    probabilities. Array dimensions correspond to MCMC samples, sites, and distance band.}

  \item{rhat}{a list of Gelman-Rubin diagnostic values for some of the model
    parameters.}

  \item{ESS}{a list of effective sample sizes for some of the model parameters.}

  \item{run.time}{execution time reported using \code{proc.time()}.}

  The return object will include additional objects used for
  subsequent prediction and/or model fit evaluation.
}

\examples{
set.seed(123)
J.x <- 10
J.y <- 10
J <- J.x * J.y
# Number of distance bins from which to simulate data.
n.bins <- 5
# Length of each bin. This should be of length n.bins
bin.width <- c(.10, .10, .20, .3, .1)
# Abundance coefficients
beta <- c(1.0, 0.2, 0.3, -0.2)
p.abund <- length(beta)
# Detection coefficients
alpha <- c(-1.0, -0.3)
p.det <- length(alpha)
# Detection decay function
det.func <- 'halfnormal'
mu.RE <- list()
p.RE <- list()
sp <- FALSE
family <- 'NB'
kappa <- 0.1
offset <- 1.8
transect <- 'point'

dat <- simDS(J.x = J.x, J.y = J.y, n.bins = n.bins, bin.width = bin.width,
             beta = beta, alpha = alpha, det.func = det.func, kappa = kappa,
             mu.RE = mu.RE, p.RE = p.RE, sp = sp,
             sigma.sq = sigma.sq, phi = phi, nu = nu, family = family,
             offset = offset, transect = transect)

y <- dat$y
X <- dat$X
X.re <- dat$X.re
X.p <- dat$X.p
X.p.re <- dat$X.p.re
dist.breaks <- dat$dist.breaks

covs <- cbind(X, X.p)
colnames(covs) <- c('int.abund', 'abund.cov.1', 'abund.cov.2', 'abund.cov.3',
                    'int.det', 'det.cov.1')

data.list <- list(y = y,
                  covs = covs,
                  dist.breaks = dist.breaks,
                  offset = offset)

# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 10),
                   alpha.normal = list(mean = 0,
                                       var = 10),
                   kappa.unif = c(0, 100))
# Starting values
inits.list <- list(alpha = 0,
                   beta = 0,
                   kappa = 1)
# Tuning values
tuning <- list(beta = 0.1, alpha = 0.1, beta.star = 0.3, alpha.star = 0.1,
               kappa = 0.2)

out <- DS(abund.formula = ~ abund.cov.1 + abund.cov.2 + abund.cov.3,
          det.formula = ~ det.cov.1,
          data = data.list,
          n.batch = 10,
          batch.length = 25,
          inits = inits.list,
          family = 'NB',
          det.func = 'halfnormal',
          transect = 'point',
          tuning = tuning,
          priors = prior.list,
          accept.rate = 0.43,
          n.omp.threads = 1,
          verbose = TRUE,
          n.report = 100,
          n.burn = 100,
          n.thin = 1,
          n.chains = 1)
summary(out)
}
