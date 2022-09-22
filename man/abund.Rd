\name{abund}
\alias{abund}
\title{Function for Fitting Single-Species Abundance GLMs}

\usage{
abund(formula, data, inits, priors, tuning,
     n.batch, batch.length, accept.rate = 0.43, family = 'NB',
     n.omp.threads = 1, verbose = TRUE,
     n.report = 100, n.burn = round(.10 * n.batch * batch.length), n.thin = 1, 
     n.chains = 1, k.fold, k.fold.threads = 1, k.fold.seed = 100, 
     k.fold.only = FALSE, ...)
}

\description{
  Function for fitting single-species abundance generalized linear models 
}

\arguments{
\item{formula}{a symbolic description of the model to be fit
  for the model using R's model syntax. Only
  right-hand side of formula is specified. See example below. Random intercepts
  and slopes are allowed using lme4 syntax (Bates et al. 2015).}

\item{data}{a list containing data necessary for model fitting.
  Valid tags are \code{y} and \code{abund.covs}. \code{y}
  is a vector, matrix, or data frame of the observed count values. If a vector, 
  the values represent the observed counts at each site. If multiple replicate
  observations are obtained at the sites (e.g., sub-samples, repeated sampling over
  multiple seasons), \code{y} can be specified as a matrix or data frame 
  with first dimension equal to the number of 
  sites (\eqn{J}{J}) and second dimension equal to the maximum number of 
  replicates at a given site. \code{abund.covs} is a list 
  containing the variables used in the model. Each list element is a different
  covariate, which can be site-level or observation-level. Site-level covariates
  are specified as a vector of length \eqn{J}{J}, while observation-level covariates
  are specified as a matrix or data frame with the number of rows equal to \eqn{J}{J}
  and number of columns equal to the maximum number of replicate observations at a 
  given site.} 

\item{inits}{a list with each tag corresponding to a parameter name.
  Valid tags are \code{beta}, \code{kappa}, and \code{sigma.sq.mu}. 
  The value portion of each tag is the 
  parameter's initial value. \code{sigma.sq.mu} is 
  only relevant when including random effects in the model. \code{kappa} is 
  only relevant when \code{family = 'NB'}. See \code{priors} 
  description for definition of each parameter name.
  Additionally, the tag \code{fix} can be set to \code{TRUE} 
  to fix the starting values across all chains. If \code{fix} is not specified
  (the default), starting values are varied randomly across chains.}

\item{priors}{a list with each tag corresponding to a parameter name. 
  Valid tags are \code{beta.normal}, \code{kappa.unif}, and \code{sigma.sq.mu.ig}. 
  Abundance (\code{beta}) regression coefficients are assumed to follow a normal distribution. 
  The hyperparameters of the normal distribution are passed as a list of 
  length two with the first and second elements corresponding to the mean 
  and variance of the normal distribution, which are each specified as vectors of 
  length equal to the number of coefficients to be estimated or of length one
  if priors are the same for all coefficients. If not specified, prior means 
  are set to 0 and prior variances set to 100. \code{kappa} is the negative binomial 
  dispersion parameter and is assumed to follow a uniform distribution. The 
  hyperparameters of the uniform distribution are passed as a vector of length
  two with the first and second elements corresponding to the lower and upper
  bounds of the uniform distribution. \code{sigma.sq.mu}
  are the random effect variances for any abundance
  random effects, respectively, and are assumed to follow an inverse
  Gamma distribution. The hyperparameters of the inverse-Gamma distribution
  are passed as a list of length two with first and second elements corresponding
  to the shape and scale parameters, respectively, which are each specified as
  vectors of length equal to the number of random intercepts or of length one 
  if priors are the same for all random effect variances.}

\item{tuning}{a single numeric value representing the initial tuning variance of the
  adaptive sampler for \code{kappa}, the negative binomial dispersion parameter. 
  name. See Roberts and Rosenthal (2009) for details.}

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

\item{k.fold}{specifies the number of \emph{k} folds for cross-validation. 
  If not specified as an argument, then cross-validation is not performed
  and \code{k.fold.threads} and \code{k.fold.seed} are ignored. In \emph{k}-fold
  cross-validation, the data specified in \code{data} is randomly
  partitioned into \emph{k} equal sized subsamples. Of the \emph{k} subsamples, 
  \emph{k} - 1 subsamples are used to fit the model and the remaining \emph{k}
  samples are used for prediction. The cross-validation process is repeated 
  \emph{k} times (the folds). As a scoring rule, we use the model deviance 
  as described in Hooten and Hobbs (2015). Cross-validation is performed
  after the full model is fit using all the data. Cross-validation results
  are reported in the \code{k.fold.deviance} object in the return list.}

\item{k.fold.threads}{number of threads to use for cross-validation. If 
  \code{k.fold.threads > 1} parallel processing is accomplished using the 
  \pkg{foreach} and \pkg{doParallel} packages. Ignored if \code{k.fold}
  is not specified.} 

\item{k.fold.seed}{seed used to split data set into \code{k.fold} parts
  for k-fold cross-validation. Ignored if \code{k.fold} is not specified.}

\item{k.fold.only}{logical value indicating whether to only perform
  k-fold cross-validation. Setting this to \code{TRUE} will only return
  the k-fold deviance.}

\item{...}{currently no additional arguments}
}

\note{
 Some of the underlying code used for generating random numbers from the
 Polya-Gamma distribution is taken from the \pkg{pgdraw} package written
 by Daniel F. Schmidt and Enes Makalic. Their code implements Algorithm
 6 in PhD thesis of Jesse Bennett Windle (2013) \url{https://repositories.lib.utexas.edu/handle/2152/21842}.
}

\references{

  Bates, Douglas, Martin Maechler, Ben Bolker, Steve Walker (2015).
  Fitting Linear Mixed-Effects Models Using lme4. Journal of
  Statistical Software, 67(1), 1-48. \doi{10.18637/jss.v067.i01}.

  Polson, N.G., J.G. Scott, and J. Windle. (2013) Bayesian Inference for
  Logistic Models Using Polya-Gamma Latent Variables. 
  \emph{Journal of the American Statistical Association}, 108:1339-1349.

}

\author{
  Jeffrey W. Doser \email{doserjef@msu.edu}, \cr
  Andrew O. Finley \email{finleya@msu.edu}
}

\value{
  An object of class \code{abund} that is a list comprised of: 

  \item{beta.samples}{a \code{coda} object of posterior samples
    for the regression coefficients.}

  \item{kappa.samples}{a \code{coda} object of posterior samples
    for the abundance dispersion parameter. Only included when
    \code{family = 'NB'}.}

  \item{y.rep.samples}{a three-dimensional array of posterior samples 
    for the abundance replicate (fitted) values with dimensions 
    corresponding to MCMC samples, site, and replicate.}

  \item{mu.samples}{a three-dimensional array of posterior samples
    for the expected abundance samples with dimensions corresponding 
    to MCMC samples, site, and replicate.}

  \item{sigma.sq.mu.samples}{a \code{coda} object of posterior samples
    for variances of random intercepts included in the abundance portion
    of the model. Only included if random intercepts are specified in 
    \code{formula}.}

  \item{beta.star.samples}{a \code{coda} object of posterior samples
    for the random effects. Only included if random intercepts
    are specified in \code{formula}.}

  \item{like.samples}{a \code{coda} object of posterior samples
    for the likelihood value associated with each site. Used for calculating
    WAIC.}

  \item{rhat}{a list of Gelman-Rubin diagnostic values for some of the model
    parameters.}

  \item{ESS}{a list of effective sample sizes for some of the model parameters.}

  \item{run.time}{execution time reported using \code{proc.time()}.}

  \item{k.fold.deviance}{scoring rule (deviance) from k-fold cross-validation. 
    Only included if \code{k.fold} is specified in function call.}

  The return object will include additional objects used for 
  subsequent prediction and/or model fit evaluation. Note that detection
  probability estimated values are not included in the model object, but can be 
  extracted using \code{fitted()}. 

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

abund.covs <- list(int = X[, , 1], 
		   abund.cov.1 = X[, , 2], 
		   abund.cov.2 = X[, , 3], 
		   abund.cov.3 = X[, , 4],
		   abund.factor.1 = X.re[, , 1])

data.list <- list(y = y, 
		  covs = abund.covs)

# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 100),
                   kappa.unif = c(0.001, 10)) 
# Starting values
inits.list <- list(beta = 0,
		   kappa = kappa)

tuning <- 0.3
n.batch <- 5
batch.length <- 25
n.burn <- 0
n.thin <- 1
n.chains <- 1

out <- abund(formula = ~ abund.cov.1 + abund.cov.2 + abund.cov.3 + (1 | abund.factor.1),
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
}
