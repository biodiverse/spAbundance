\name{abund}
\alias{abund}
\title{Function for Fitting Univariate Abundance GLMMs}

\usage{
abund(formula, data, inits, priors, tuning,
      n.batch, batch.length, accept.rate = 0.43, family = 'Poisson',
      n.omp.threads = 1, verbose = TRUE,
      n.report = 100, n.burn = round(.10 * n.batch * batch.length), n.thin = 1, 
      n.chains = 1, save.fitted = TRUE, ...)
}

\description{
  Function for fitting univariate abundance generalized linear (mixed) models 
}

\arguments{
\item{formula}{a symbolic description of the model to be fit
  for the model using R's model syntax. Only
  right-hand side of formula is specified. See example below. Random intercepts
  and slopes are allowed using lme4 syntax (Bates et al. 2015).}

\item{data}{a list containing data necessary for model fitting.
  Valid tags are \code{y}, \code{covs}, \code{z}, and \code{offset}. \code{y}
  is a vector, matrix, or data frame of the observed count values. If a vector, 
  the values represent the observed counts at each site. If multiple replicate
  observations are obtained at the sites (e.g., sub-samples, repeated sampling over
  multiple seasons), \code{y} can be specified as a matrix or data frame 
  with first dimension equal to the number of 
  sites (\eqn{J}{J}) and second dimension equal to the maximum number of 
  replicates at a given site. \code{covs} is a list or data frame 
  containing the variables used in the model. Each list element is a different
  covariate, which can be site-level or observation-level. Site-level covariates
  are specified as a vector of length \eqn{J}{J} (or column
  in a data frame), while observation-level covariates
  are specified as a matrix or data frame with the number of rows equal to \eqn{J}{J}
  and number of columns equal to the maximum number of replicate observations at a 
  given site. For zero-inflated Gaussian models, the tag \code{z} is used to specify the 
  binary component of the zero-inflated model and should have the same length as \code{y}. 
  \code{offset} is an offset to use in the abundance model (e.g., an area offset). 
  This can be either a single value, a vector with an offset for each site (e.g., if survey area differed in size), or a site x replicate matrix if more than one count is available at a given site.} 

\item{inits}{a list with each tag corresponding to a parameter name.
  Valid tags are \code{beta}, \code{kappa}, \code{sigma.sq.mu}, and \code{tau.sq}. 
  The value portion of each tag is the 
  parameter's initial value. \code{sigma.sq.mu} is 
  only relevant when including random effects in the model. \code{kappa} is 
  only relevant when \code{family = 'NB'}. \code{tau.sq} is only 
  relevant when \code{family = 'Gaussian'} or \code{family = 'zi-Gaussian'}. 
  See \code{priors} description for definition of each parameter name.
  Additionally, the tag \code{fix} can be set to \code{TRUE} 
  to fix the starting values across all chains. If \code{fix} is not specified
  (the default), starting values are varied randomly across chains.}

\item{priors}{a list with each tag corresponding to a parameter name. 
  Valid tags are \code{beta.normal}, \code{kappa.unif}, \code{sigma.sq.mu.ig}, 
  and \code{tau.sq.ig}. 
  Abundance (\code{beta}) regression coefficients are assumed to follow a normal distribution. 
  The hyperparameters of the normal distribution are passed as a list of 
  length two with the first and second elements corresponding to the mean 
  and variance of the normal distribution, which are each specified as vectors of 
  length equal to the number of coefficients to be estimated or of length one
  if priors are the same for all coefficients. If not specified, prior means 
  are set to 0 and prior variances set to 100. \code{kappa} is the negative binomial 
  over-dispersion parameter and is assumed to follow a uniform distribution. The 
  hyperparameters of the uniform distribution are passed as a vector of length
  two with the first and second elements corresponding to the lower and upper
  bounds of the uniform distribution. \code{sigma.sq.mu}
  are the random effect variances for any abundance
  random effects, respectively, and are assumed to follow an inverse
  Gamma distribution. The hyperparameters of the inverse-Gamma distribution
  are passed as a list of length two with first and second elements corresponding
  to the shape and scale parameters, respectively, which are each specified as
  vectors of length equal to the number of random effects or of length one 
  if priors are the same for all random effect variances. \code{tau.sq} is the 
  residual variance for Gaussian (or zero-inflated Gaussian) models, and it is assigned 
  an inverse-Gamma prior. The hyperparameters of the inverse-Gamma are passed as a vector
  of length two, with the first and second element corresponding to the shape and 
  scale parameters, respectively.}

\item{tuning}{a list with each tag corresponding to a parameter name, 
  whose value defines the initial variance of the adaptive sampler.
  Valid tags are \code{beta}, \code{beta.star} (the abundance
  random effect values), and \code{kappa}. See Roberts and Rosenthal (2009) for details.
  Note that no tuning is necessary for Gaussian or zero-inflated Gaussian models.}

\item{n.batch}{the number of MCMC batches in each chain to run for the adaptive MCMC 
  sampler. See Roberts and Rosenthal (2009) for details.}

\item{batch.length}{the length of each MCMC batch in each chain to run for the Adaptive 
  MCMC sampler. See Roberts and Rosenthal (2009) for details.}

\item{accept.rate}{target acceptance rate for Adaptive MCMC. Default is 
  0.43. See Roberts and Rosenthal (2009) for details.}

\item{family}{the distribution to use for the latent abundance process. Currently 
  supports \code{'NB'} (negative binomial), \code{'Poisson'}, \code{'Gaussian'}, 
  and \code{'zi-Gaussian'}.} 

\item{n.omp.threads}{a positive integer indicating the number of threads
  to use for SMP parallel processing. The package must be compiled for 
  OpenMP support. For most Intel-based machines, we recommend setting 
  \code{n.omp.threads} up to the number of hypterthreaded cores. Note, 
  \code{n.omp.threads} > 1 might not work on some systems. Currently only 
  relevant for spatially-explicit models.}

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

  Bates, Douglas, Martin Maechler, Ben Bolker, Steve Walker (2015).
  Fitting Linear Mixed-Effects Models Using lme4. Journal of
  Statistical Software, 67(1), 1-48. \doi{10.18637/jss.v067.i01}.

}

\author{
  Jeffrey W. Doser \email{doserjef@msu.edu}, \cr
  Andrew O. Finley \email{finleya@msu.edu}, \cr
}

\value{
  An object of class \code{abund} that is a list comprised of: 

  \item{beta.samples}{a \code{coda} object of posterior samples
    for the regression coefficients.}

  \item{kappa.samples}{a \code{coda} object of posterior samples
    for the abundance overdispersion parameter. Only included when
    \code{family = 'NB'}.}
  
  \item{tau.sq.samples}{a \code{coda} object of posterior samples
    for the Gaussian residual variance parameter. Only included when
    \code{family = 'Gaussian'} or \code{family = 'zi-Gaussian'}.}

  \item{y.rep.samples}{a two or three-dimensional array of posterior samples 
    for the abundance replicate (fitted) values with dimensions 
    corresponding to MCMC samples, site, and an optional third dimension of replicate.}

  \item{mu.samples}{a two or three-dimensional array of posterior samples
    for the expected abundance samples with dimensions corresponding 
    to MCMC samples, site, and an optional third dimension of replicate.}

  \item{sigma.sq.mu.samples}{a \code{coda} object of posterior samples
    for variances of random effects included in the model. 
    Only included if random effects are specified in \code{formula}.}

  \item{beta.star.samples}{a \code{coda} object of posterior samples
    for the random effects. Only included if random effects 
    are specified in \code{formula}.}

  \item{like.samples}{a \code{coda} object of posterior samples
    for the likelihood value associated with each site. Used for calculating
    WAIC.}

  \item{rhat}{a list of Gelman-Rubin diagnostic values for some of the model
    parameters.}

  \item{ESS}{a list of effective sample sizes for some of the model parameters.}

  \item{run.time}{execution time reported using \code{proc.time()}.}

  The return object will include additional objects used for 
  subsequent prediction and/or model fit evaluation.

}

\examples{
set.seed(1010)
J.x <- 15
J.y <- 15
J <- J.x * J.y
n.rep <- sample(3, J, replace = TRUE)
beta <- c(0, -1.5, 0.3, -0.8)
p.abund <- length(beta)
mu.RE <- list(levels = c(30),
	       sigma.sq.mu = c(1.3))
kappa <- 0.5
sp <- FALSE 
family <- 'NB'
dat <- simAbund(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, 
                kappa = kappa, mu.RE = mu.RE, sp = sp, family = 'NB')

y <- dat$y
X <- dat$X
X.re <- dat$X.re

covs <- list(int = X[, , 1], 
             abund.cov.1 = X[, , 2], 
             abund.cov.2 = X[, , 3], 
             abund.cov.3 = X[, , 4],
             abund.factor.1 = X.re[, , 1])

data.list <- list(y = y, covs = covs)

# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 100),
                   kappa.unif = c(0.001, 10)) 
# Starting values
inits.list <- list(beta = 0, kappa = kappa)

tuning <- list(kappa = 0.2, beta = 0.1, beta.star = 0.2)
n.batch <- 5
batch.length <- 25
n.burn <- 0
n.thin <- 1
n.chains <- 1

out <- abund(formula = ~ abund.cov.1 + abund.cov.2 + abund.cov.3 + 
                         (1 | abund.factor.1),
                         data = data.list, 
                         n.batch = n.batch, 
                         batch.length = batch.length, 
                         inits = inits.list, 
                         tuning = tuning,
                         priors = prior.list, 
                         accept.rate = 0.43, 
                         n.omp.threads = 1, 
                         verbose = TRUE, 
                         n.report = 1,
                         n.burn = n.burn,
                         n.thin = n.thin,
                         n.chains = n.chains) 
summary(out)
}
