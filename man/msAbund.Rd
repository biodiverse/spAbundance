\name{msAbund}
\alias{msAbund}
\title{Function for Fitting Multivariate Abundance GLMMs}

\description{
  The function \code{msAbund} fits multivariate abundance GLMMs. 
}

\usage{
msAbund(formula, data, inits, priors, tuning, 
        n.batch, batch.length, accept.rate = 0.43, family = 'Poisson',
        n.omp.threads = 1, verbose = TRUE, n.report = 100, 
        n.burn = round(.10 * n.batch * batch.length), n.thin = 1, n.chains = 1,
        save.fitted = TRUE, ...)
}

\arguments{

  \item{formula}{a symbolic description of the model to be fit
  for the model using R's model syntax. Only
  right-hand side of formula is specified. See example below. Random intercepts
  and slopes are allowed using lme4 syntax (Bates et al. 2015).}

  \item{data}{a list containing data necessary for model fitting.
    Valid tags are \code{y}, \code{covs}, \code{z}, and \code{offset}. 
    \code{y} is a two or three-dimensional array of observed count data. The
    first dimension of the array is equal to the 
    number of species and the second dimension is equal to the number of sites. If
    specified as a three-dimensional array, the third dimension corresponds to 
    replicate observations at each site (e.g., sub-samples, repeated sampling 
    over multiple seasons). \code{covs} is a list or data frame
    containing the variables used in the model. If a data frame, each row 
    of \code{covs} is a site and each column is a variable. If a list, 
    each list element is a different covariate, which 
    can be site-level or observation-level. Site-level covariates
    are specified as a vector of length \eqn{J}{J}, while observation-level covariates
    are specified as a matrix or data frame with the number of rows equal to \eqn{J}{J}
    and number of columns equal to the maximum number of replicate observations at a 
    given site. For zero-inflated Gaussian models, the tag \code{z} is used to specify the 
    binary component of the model and should have the same dimensions as \code{y}. \code{offset} 
    is an offset to use in the abundance model (e.g., an area offset). 
    This can be either a single value, a vector with an offset for each site (e.g., if survey area 
    differed in size), or a site x replicate matrix if more than one count is available at a given site.}

  \item{inits}{a list with each tag corresponding to a parameter name.
    Valid tags are \code{beta.comm}, \code{beta}, 
    \code{tau.sq.beta}, \code{sigma.sq.mu}, \code{kappa}, \code{tau.sq}. 
    \code{kappa} is only specified if \code{family = 'NB'}, \code{tau.sq} is only
    specified for Gaussian or zero-inflated Gaussian models, 
    and \code{sigma.sq.mu} is only specified if random effects are included in \code{formula}. 
    The value portion of each tag is 
    the parameter's initial value. See \code{priors} description for definition
    of each parameter name. Additionally, the tag \code{fix} can be set to \code{TRUE} 
    to fix the starting values across all chains. If \code{fix} is not specified
    (the default), starting values are varied randomly across chains.}

  \item{priors}{a list with each tag corresponding to a parameter name. 
  Valid tags are \code{beta.comm.normal}, \code{tau.sq.beta.ig}, \code{sigma.sq.mu}, 
  \code{kappa.unif}, \code{tau.sq.ig}. Community-level (\code{beta.comm}) 
  regression coefficients are assumed to follow a 
  normal distribution. The hyperparameters of the normal distribution
  are passed as a list of length two with the first and second elements 
  corresponding to the mean and variance of the normal distribution, 
  which are each specified as vectors of length equal to the number of 
  coefficients to be estimated or of length one if priors are the same for 
  all coefficients. If not specified, prior means are set 
  to 0 and prior variances to 100. Community-level variance parameters 
  (\code{tau.sq.beta}) are 
  assumed to follow an inverse Gamma distribution. The hyperparameters of 
  the inverse gamma distribution are passed as a list of length two with 
  the first and second elements corresponding to the shape and scale parameters,
  which are each specified as vectors of length equal to the number of 
  coefficients to be estimated or a single value if priors are the same for all 
  parameters. If not specified, prior shape and scale 
  parameters are set to 0.1. \code{sigma.sq.mu} are the random 
  effect variances random effects, respectively, and are assumed to follow an inverse
  Gamma distribution. The hyperparameters of the inverse-Gamma distribution
  are passed as a list of length two with first and second elements corresponding
  to the shape and scale parameters, respectively, which are each specified as
  vectors of length equal to the number of random intercepts or of length one 
  if priors are the same for all random effect variances. \code{kappa} is the 
  negative binomial dispersion parameter for each species and is assumed to 
  follow a uniform distribution. The hyperparameters of the uniform distribution
  are passed as a list of length two with first and second elements corresponding to the
  lower and upper bounds of the uniform distribution, respectively, which are each 
  specified as vectors of length equal to the number of species or of length one
  if priors are the same for all species-specific dispersion parameters. \code{tau.sq} is the 
  species-specific residual variance for Gaussian (or zero-inflated Gaussian) models, and it is assigned 
  an inverse-Gamma prior. The hyperparameters of the inverse-Gamma are passed as a list 
  of length two, with the first and second element corresponding to the shape and 
  scale parameters, respectively, which are each specified as vectors of length
  equal to the number of species or a single value if priors are the same for all species.}

\item{tuning}{a list with each tag corresponding to a parameter name, 
  whose value defines the initial variance of the adaptive sampler.
  Valid tags are \code{beta}, \code{beta.star} (the abundance
  random effect values), and \code{kappa}. See Roberts and Rosenthal (2009) for details. 
  Note that no tuning is necessary for Gaussian or zero-inflated Gaussian models.}

  \item{n.batch}{the number of MCMC batches  in each chain to run for the adaptive MCMC 
    sampler. See Roberts and Rosenthal (2009) for details.}
  
  \item{batch.length}{the length of each MCMC batch to run for the adaptive 
    MCMC sampler. See Roberts and Rosenthal (2009) for details.}
  
  \item{accept.rate}{target acceptance rate for adaptive MCMC. Defaul is 
    0.43. See Roberts and Rosenthal (2009) for details.}

  \item{family}{the distribution to use for the latent abundance process. Currently 
    supports \code{'NB'} (negative binomial), \code{'Poisson'}, \code{'Gaussian'}, 
    and \code{'zi-Gaussian'}.} 
  
  \item{n.omp.threads}{a positive integer indicating
   the number of threads to use for SMP parallel processing. The package must
   be compiled for OpenMP support. For most Intel-based machines, we
   recommend setting \code{n.omp.threads} up to the number of
   hyperthreaded cores. Note, \code{n.omp.threads} > 1 might not
   work on some systems.}
 
  \item{verbose}{if \code{TRUE}, messages about data preparation, 
    model specification, and progress of the sampler are printed to the screen. 
    Otherwise, no messages are printed.}
 
  \item{n.report}{the interval to report Metropolis sampler acceptance and
    MCMC progress. Note this is specified in terms of batches and not overall
    samples for spatial models.}

  \item{n.burn}{the number of samples out of the total \code{n.samples} to 
    discard as burn-in for each chain. By default, the first 10\% of samples is discarded.}
  
  \item{n.thin}{the thinning interval for collection of MCMC samples. The
    thinning occurs after the \code{n.burn} samples are discarded. Default 
    value is set to 1.}

  \item{n.chains}{the number of chains to run in sequence.}

  \item{save.fitted}{logical value indicating whether or not fitted values and likelihood values
    should be saved in the resulting model object. If \code{save.fitted = FALSE}, the components
    \code{y.rep.samples}, \code{mu.samples}, and \code{like.samples} will not be included
    in the model object, and subsequent functions for calculating WAIC, fitted values, and
    posterior predictive checks will not work, although they all can be calculated manually if
    desired. Setting \code{save.fitted = FALSE} can be useful when working with very large 
    data sets to minimize the amount of RAM needed when fitting and storing the model object in 
    memory.}

  \item{...}{currently no additional arguments}
}

\references{
  Roberts, G.O. and Rosenthal J.S. (2009) Examples  of adaptive MCMC.
  \emph{Journal of Computational and Graphical Statistics}, 18(2):349-367.

  Bates, Douglas, Martin Maechler, Ben Bolker, Steve Walker (2015).
  Fitting Linear Mixed-Effects Models Using lme4. Journal of
  Statistical Software, 67(1), 1-48. \doi{10.18637/jss.v067.i01}.
}

\author{
  Jeffrey W. Doser \email{doserjef@msu.edu}, \cr
  Andrew O. Finley \email{finleya@msu.edu}, \cr
}

\value{
  An object of class \code{msAbund} that is a list comprised of: 

  \item{beta.comm.samples}{a \code{coda} object of posterior samples
    for the community level regression coefficients.}
  
  \item{tau.sq.beta.samples}{a \code{coda} object of posterior samples
    for the abundance community variance parameters.}

  \item{beta.samples}{a \code{coda} object of posterior samples
    for the species level abundance regression coefficients.}

  \item{kappa.samples}{a \code{coda} object of posterior samples
    for the species level abundance dispersion parameters. Only included
    when \code{family = 'NB'}.}

  \item{tau.sq.samples}{a \code{coda} object of posterior samples
    for the Gaussian residual variance parameter. Only included when
    \code{family = 'Gaussian'} or \code{family = 'zi-Gaussian'}.}

  \item{y.rep.samples}{a three or four-dimensional array of posterior samples for
    the fitted (replicate) values for each species with dimensions corresponding
    to MCMC sample, species, site, and replicate.}

  \item{mu.samples}{a three or four-dimensional array of posterior samples for
    the expected abundance values for each species with dimensions
    corresponding to MCMC samples, species, site, and replicate.}

  \item{sigma.sq.mu.samples}{a \code{coda} object of posterior samples
    for variances of random effects included in the abundance portion
    of the model. Only included if random effects are specified in 
    \code{abund.formula}.}

  \item{beta.star.samples}{a \code{coda} object of posterior samples
    for the abundance random effects. Only included if random effects 
    are specified in \code{abund.formula}.}

  \item{like.samples}{a three-dimensional array of posterior samples
    for the likelihood value associated with each site and species. 
    Used for calculating WAIC.}

  \item{rhat}{a list of Gelman-Rubin diagnostic values for some of the model
    parameters.}

  \item{ESS}{a list of effective sample sizes for some of the model parameters.}

  \item{run.time}{MCMC sampler execution time reported using \code{proc.time()}.}

  The return object will include additional objects used for 
  subsequent prediction and/or model fit evaluation. 
}

\examples{
set.seed(408)
J.x <- 8
J.y <- 8
J <- J.x * J.y
n.rep <- sample(3, size = J, replace = TRUE)
n.sp <- 6
# Community-level covariate effects
beta.mean <- c(-2, 0.5)
p.abund <- length(beta.mean)
tau.sq.beta <- c(0.2, 1.2)
# Random effects (two random intercepts)
mu.RE <- list(levels = c(10, 15),
              sigma.sq.mu = c(0.43, 0.5))
# Draw species-level effects from community means.
beta <- matrix(NA, nrow = n.sp, ncol = p.abund)
for (i in 1:p.abund) {
  beta[, i] <- rnorm(n.sp, beta.mean[i], sqrt(tau.sq.beta[i]))
}
sp <- FALSE
kappa <- runif(n.sp, 0.1, 1)

dat <- simMsAbund(J.x = J.x, J.y = J.y, n.rep = n.rep, n.sp = n.sp, beta = beta,
                  mu.RE = mu.RE, sp = sp, kappa = kappa, family = 'NB')

y <- dat$y
X <- dat$X
X.re <- dat$X.re

# Package all data into a list
covs <- list(int = X[, , 1],
             abund.cov.1 = X[, , 2],
             abund.factor.1 = X.re[, , 1],
             abund.factor.2 = X.re[, , 2])
data.list <- list(y = y, covs = covs)
prior.list <- list(beta.comm.normal = list(mean = 0, var = 100),
                   kappa.unif = list(a = 0, b = 10),
                   tau.sq.beta.ig = list(a = .1, b = .1))
inits.list <- list(beta.comm = 0,
                   beta = 0,
                   kappa = 0.5,
                   tau.sq.beta = 1)
tuning.list <- list(kappa = 0.3, beta = 0.1, beta.star = 0.1)

# Small
n.batch <- 2
batch.length <- 25
n.burn <- 20
n.thin <- 1
n.chains <- 1

out <- msAbund(formula = ~ abund.cov.1 + (1 | abund.factor.1) + 
                           (1 | abund.factor.2),
               data = data.list,
               n.batch = n.batch,
               inits = inits.list,
               priors = prior.list,
               tuning = tuning.list,
               batch.length = batch.length,
               n.omp.threads = 3,
               verbose = TRUE,
               n.report = 1,
               n.burn = n.burn,
               n.thin = n.thin,
               n.chains = n.chains)
summary(out)
}
