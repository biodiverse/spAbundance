\name{msDS}
\alias{msDS}
\title{Function for Fitting Multi-Species Hierarchical Distance Sampling Models}

\description{
Function for fitting multi-species hierarchical distance sampling models. 
}

\usage{
msDS(abund.formula, det.formula, data, inits, priors,  
     tuning, n.batch, batch.length, accept.rate = 0.43, 
     family = 'Poisson', transect = 'line', det.func = 'halfnormal', 
     n.omp.threads = 1, verbose = TRUE, n.report = 100, 
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
  Valid tags are \code{y}, \code{covs}, and \code{dist.breaks}, and \code{offset}. \code{y}
  is a three-dimensional array of observed count data with first dimension equal to the number 
  of species, second dimension equal to the number of sites, and third 
  dimension equal to the maximum number of replicates at a given site. 
  \code{covs} is a matrix or data frame containing the variables 
  used in the abundance and/or the detection portion of the model, with
  \eqn{J}{J} rows for each column (variable). \code{dist.breaks} is a vector of 
  distances that denote the breakpoints of the distance bands. \code{dist.breaks} should
  have length equal to the third dimension of \code{y} plus one. \code{offset} is an 
  offset that can be used to scale estimates from abundance per transect to density per 
  some desired unit of measure. This can be either a single value or a vector with an offset
  value for each site (e.g., if transects differ in length)}

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
  specified, prior means are set to 0 and prior variances are set to 100.
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
  vectors of length equal to the number of random effects or of length one 
  if priors are the same for all random effect variances. \code{kappa} is the 
  negative binomial dispersion parameter for each species and is assumed to 
  follow a uniform distribution. The hyperparameters of the uniform distribution
  are passed as a list of length two with first and second elements corresponding to the
  lower and upper bounds of the uniform distribution, respectively, which are each 
  specified as vectors of length equal to the number of species or of length one
  if priors are the same for all species-specific dispersion parameters.}

\item{tuning}{a list with each tag corresponding to a parameter name, whose
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

  Royle, J. A. (2004). N‐mixture models for estimating population size
  from spatially replicated counts. Biometrics, 60(1), 108-115.

  Sollmann, R., Gardner, B., Williams, K. A., Gilbert, A. T., 
  & Veit, R. R. (2016). A hierarchical distance sampling model to 
  estimate abundance and covariate associations of species 
  and communities. Methods in Ecology and Evolution, 7(5), 529-537.

}

\author{
  Jeffrey W. Doser \email{doserjef@msu.edu}, \cr
}

\value{
  An object of class \code{msDS} that is a list comprised of: 

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
    the latent abundance values for each species. Note that these values
    always represent transect-level abundance, even when an offset is 
    supplied. Array dimensions correspond to MCMC sample, species, and site.}

  \item{mu.samples}{a three-dimensional array of posterior samples for
    the latent expected abundance values for each species. When an offset
    is supplied in the \code{data} object, these correspond to expected
    abundance per unit area (i.e., density). Array dimensions correspond to 
    MCMC samples, species, and site.}

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

  \item{y.rep.samples}{a four-dimensional array of fitted values. 
    Array dimensions correspond to MCMC samples, species, sites, and distance band.}

  \item{pi.samples}{a four-dimensional array of cell-specific detection
    probabilities. Array dimensions correspond to MCMC samples, species, 
    sites, and distance band.}

  \item{rhat}{a list of Gelman-Rubin diagnostic values for some of the model
    parameters.}

  \item{ESS}{a list of effective sample sizes for some of the model parameters.}

  \item{run.time}{MCMC sampler execution time reported using \code{proc.time()}.}

  The return object will include additional objects used for 
  subsequent prediction and/or model fit evaluation.
}

\examples{
set.seed(210)
J.x <- 10
J.y <- 10 
J <- J.x * J.y
# Number of distance bins from which to simulate data. 
n.bins <- 5
# Length of each bin. This should be of length n.bins
bin.width <- c(.10, .10, .20, .3, .1)
# Number of species
n.sp <- 5
# Community-level abundance coefficients
beta.mean <- c(-1, 0.2, 0.3, -0.2)
p.abund <- length(beta.mean)
tau.sq.beta <- c(0.2, 0.3, 0.5, 0.4)
# Detection coefficients
alpha.mean <- c(-1.0, -0.3)
p.det <- length(alpha.mean)
tau.sq.alpha <- c(0.1, 0.2)
# Detection decay function
det.func <- 'halfnormal'
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
sp <- FALSE 
family <- 'Poisson'
kappa <- runif(n.sp, 0.3, 3) 
offset <- pi * .8^2
transect <- 'line'
factor.model <- FALSE

dat <- simMsDS(J.x = J.x, J.y = J.y, n.bins = n.bins, bin.width = bin.width,
               n.sp = n.sp, beta = beta, alpha = alpha, det.func = det.func, 
               mu.RE = mu.RE, p.RE = p.RE, sp = sp, cov.model = cov.model,
               sigma.sq = sigma.sq, phi = phi, nu = nu, family = family, 
               offset = offset, transect = transect, factor.model = factor.model)

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
prior.list <- list(beta.comm.normal = list(mean = 0, var = 10),
                   alpha.comm.normal = list(mean = 0,
                                            var = 10), 
                   kappa.unif = list(0, 100), 
                   tau.sq.beta.ig = list(a = 0.1, b = 0.1),
                   tau.sq.alpha.ig = list(a = 0.1, b = 0.1)) 
# Starting values
inits.list <- list(alpha.comm = 0, beta.comm = 0, beta = 0,
                   alpha = 0, kappa = 1)

tuning <- list(beta = 0.1, alpha = 0.1, beta.star = 0.3, alpha.star = 0.1, 
               kappa = 0.8) 

n.batch <- 4 
batch.length <- 25
n.burn <- 0
n.thin <- 1
n.chains <- 1

out <- msDS(abund.formula = ~ abund.cov.1 + abund.cov.2 + abund.cov.3,
            det.formula = ~ det.cov.1,
            data = data.list, 
            n.batch = n.batch, 
            batch.length = batch.length, 
            inits = inits.list, 
            family = 'Poisson',
            det.func = 'halfnormal', 
            transect = transect, 
            tuning = tuning,
            priors = prior.list, 
            accept.rate = 0.43, 
            n.omp.threads = 1, 
            verbose = TRUE, 
            n.report = 10,
            n.burn = n.burn,
            n.thin = n.thin,
            n.chains = n.chains) 
summary(out, level = 'community')
}
