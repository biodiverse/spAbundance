\name{sfMsAbund}
\alias{sfMsAbund}
\title{Function for Fitting Spatial Factor Multivariate Abundance GLMMs}

\description{
  The function \code{sfMsAbund} fits multivariate spatial abundance GLMMs with species correlations (i.e., a spatially-explicit abundace-based joint species distribution model). We use  a spatial factor modeling approach. Currently, models are implemented using a Nearest Neighbor Gaussian Process. Future development may allow for running the models using full Gaussian Processes. 
}

\usage{
sfMsAbund(formula, data, inits, priors,  
          tuning, cov.model = 'exponential', NNGP = TRUE, 
          n.neighbors = 15, search.type = 'cb', n.factors,
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
    Valid tags are \code{y}, \code{covs}, \code{z}, \code{coords}, and \code{offset}. 
    \code{y} is a two or three-dimensional array of observed count data. The
    first dimension of the array is equal to the 
    number of species and the second dimension is equal to the number of sites. If
    specified as a three-dimensional array, the third dimension corresponds to 
    replicate observations at each site (e.g., sub-samples, repeated sampling 
    over multiple seasons). \code{covs} is a list 
    containing the variables used in the model. If a data frame, each row 
    of \code{covs} is a site and each column is a variable. 
    If a list, each list element is a different
    covariate, which can be site-level or observation-level. Site-level covariates
    are specified as a vector of length \eqn{J}{J}, while observation-level covariates
    are specified as a matrix or data frame with the number of rows equal to \eqn{J}{J}
    and number of columns equal to the maximum number of replicate observations at a 
    given site. \code{coords} is a 
    \eqn{J \times 2}{J x 2} matrix of the observation coordinates. Note that 
    \code{spAbundance} assumes coordinates are specified in a projected coordinate system. 
    For zero-inflated Gaussian models, the tag \code{z} is used to specify the 
    binary component of the model and should have the same dimensions as \code{y}. 
    \code{offset} is an offset to use in the abundance model (e.g., an area offset). 
    This can be either a single value, a vector with an offset for each site (e.g., if survey 
    area differed in size), or a site x replicate matrix if more than one count is available 
    at a given site.}

  \item{inits}{a list with each tag corresponding to a parameter name.
    Valid tags are \code{beta.comm}, \code{beta}, 
    \code{tau.sq.beta}, \code{sigma.sq.mu}, \code{kappa},
    \code{phi}, \code{lambda}, \code{nu}, and \code{tau.sq}. \code{nu} is only specified if 
    \code{cov.model = "matern"}, \code{kappa} is only specified if \code{family = 'NB'}, 
    \code{tau.sq} is only specified for Gaussian and zero-inflated Gaussian models,
    and \code{sigma.sq.mu} is only specified if random effects are included in \code{formula}. 
    The value portion of each tag is 
    the parameter's initial value. See \code{priors} description for definition
    of each parameter name. Additionally, the tag \code{fix} can be set to \code{TRUE} 
    to fix the starting values across all chains. If \code{fix} is not specified
    (the default), starting values are varied randomly across chains.}

  \item{priors}{a list with each tag corresponding to a parameter name. 
  Valid tags are \code{beta.comm.normal}, \code{tau.sq.beta.ig}, \code{sigma.sq.mu}, 
  \code{kappa.unif}, \code{phi.unif}, \code{nu.unif}, and \code{tau.sq.ig}. 
  Community-level (\code{beta.comm}) regression coefficients are assumed to follow a 
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
  parameters are set to 0.1. The spatial factor model fits \code{n.factors} independent
  spatial processes. The spatial decay \code{phi} and smoothness \code{nu} parameters 
  for each latent factor are assumed to follow Uniform distributions. 
  The hyperparameters of the Uniform are passed as a list with two elements, 
  with both elements being vectors of length \code{n.factors} corresponding to the lower and 
  upper support, respectively, or as a single value if the same value is assigned
  for all factors. The priors for the factor loadings matrix \code{lambda} are fixed
  following the standard spatial factor model to ensure parameter 
  identifiability (Christensen and Amemlya 2002). The
  upper triangular elements of the \code{n.sp x n.factors} matrix are fixed at 0 and the 
  diagonal elements are fixed at 1. The lower triangular elements are assigned a 
  standard normal prior (i.e., mean 0 and variance 1).
  \code{sigma.sq.mu} are the random 
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

\item{tuning}{a single numeric value representing the initial variance of the
  adaptive sampler for \code{beta}, \code{alpha}, \code{beta.star} (the abundance 
  random effect values), \code{kappa}, \code{phi}, \code{lambda}. 
  See Roberts and Rosenthal (2009) for details. Note that only \code{phi} and \code{nu}
  are tuned for Gaussian or zero-inflated Gaussian models.}

  \item{cov.model}{a quoted keyword that specifies the covariance
    function used to model the spatial dependence structure among the
    observations.  Supported covariance model key words are:
    \code{"exponential"}, \code{"matern"}, \code{"spherical"}, and
    \code{"gaussian"}.}
  
  \item{NNGP}{if \code{TRUE}, model is fit with an NNGP. If \code{FALSE}, 
    a full Gaussian process is used. See Datta et al. (2016) and 
    Finley et al. (2019) for more information. For spatial factor models, only 
    \code{NNGP = TRUE} is currently supported.}
  
  \item{n.neighbors}{number of neighbors used in the NNGP. Only used if 
  \code{NNGP = TRUE}. Datta et al. (2016) showed that 15 neighbors is usually 
  sufficient, but that as few as 5 neighbors can be adequate for certain data
  sets, which can lead to even greater decreases in run time. We recommend
  starting with 15 neighbors (the default) and if additional gains in computation
  time are desired, subsequently compare the results with a smaller number
  of neighbors using WAIC.}
  
  \item{search.type}{a quoted keyword that specifies the type of nearest
    neighbor search algorithm. Supported method key words are: \code{"cb"} and
    \code{"brute"}. The \code{"cb"} should generally be much
    faster. If locations do not have identical coordinate values on
    the axis used for the nearest neighbor ordering then \code{"cb"} 
    and \code{"brute"} should produce identical neighbor sets. 
    However, if there are identical coordinate values on the axis used 
    for nearest neighbor ordering, then \code{"cb"} and \code{"brute"} 
    might produce different, but equally valid, neighbor sets, 
    e.g., if data are on a grid. }

  \item{n.factors}{the number of factors to use in the spatial factor model approach. 
    Typically, the number of factors is set to be small (e.g., 4-5) relative to the 
    total number of species in the community, which will lead to substantial 
    decreases in computation time. However, the value can be anywhere 
    between 1 and the number of species in the modeled community.}

  \item{n.batch}{the number of MCMC batches  in each chain to run for the adaptive MCMC 
    sampler. See Roberts and Rosenthal (2009) for details.}
  
  \item{batch.length}{the length of each MCMC batch to run for the adaptive 
    MCMC sampler. See Roberts and Rosenthal (2009) for details.}
  
  \item{accept.rate}{target acceptance rate for adaptive MCMC. Defaul is 
    0.43. See Roberts and Rosenthal (2009) for details.}

  \item{family}{the distribution to use for the abundance. Currently 
    supports \code{'NB'} (negative binomial), \code{'Poisson'} (Poisson), \code{'Gaussian'} (Gaussian), 
    and \code{'zi-Gaussian'} (zero-inflated Gaussian).}
  
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
  Datta, A., S. Banerjee, A.O. Finley, and A.E. Gelfand. (2016)
  Hierarchical Nearest-Neighbor Gaussian process models for large
  geostatistical datasets. \emph{Journal of the American Statistical
    Association}, \doi{10.1080/01621459.2015.1044091}.

  Finley, A.O., A. Datta, B.D. Cook, D.C. Morton, H.E. Andersen, and
  S. Banerjee. (2019) Efficient algorithms for Bayesian Nearest Neighbor
  Gaussian Processes. \emph{Journal of Computational and Graphical
    Statistics}, \doi{10.1080/10618600.2018.1537924}.

  Roberts, G.O. and Rosenthal J.S. (2009) Examples  of adaptive MCMC.
  \emph{Journal of Computational and Graphical Statistics}, 18(2):349-367.

  Bates, Douglas, Martin Maechler, Ben Bolker, Steve Walker (2015).
  Fitting Linear Mixed-Effects Models Using lme4. Journal of
  Statistical Software, 67(1), 1-48. \doi{10.18637/jss.v067.i01}.

  Christensen, W. F., and Amemiya, Y. (2002). Latent variable analysis 
  of multivariate spatial data. \emph{Journal of the American Statistical Association}, 
  97(457), 302-317.

}

\author{
  Jeffrey W. Doser \email{doserjef@msu.edu}, \cr
  Andrew O. Finley \email{finleya@msu.edu}
}

\value{
  An object of class \code{sfMsAbund} that is a list comprised of: 

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

  \item{theta.samples}{a \code{coda} object of posterior samples
    for the spatial correlation parameters.}

  \item{lambda.samples}{a \code{coda} object of posterior samples
    for the latent spatial factor loadings.}

  \item{y.rep.samples}{a three or four-dimensional array of posterior samples for
    the fitted (replicate) values for each species with dimensions corresponding
    to MCMC sample, species, site, and replicate.}

  \item{mu.samples}{a three or four-dimensional array of posterior samples for
    the expected abundance values for each species with dimensions
    corresponding to MCMC samples, species, site, and replicate.}

  \item{w.samples}{a three-dimensional array of posterior samples for 
    the latent effects for each latent factor. Array dimensions correspond
    to MCMC sample, latent factor, then site.}

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
sp <- TRUE
n.factors <- 2
factor.model <- TRUE
phi <- runif(n.factors, 3/1, 3 / .1)
kappa <- runif(n.sp, 0.1, 1)

dat <- simMsAbund(J.x = J.x, J.y = J.y, n.rep = n.rep, n.sp = n.sp, beta = beta,
                  mu.RE = mu.RE, sp = sp, kappa = kappa, family = 'NB',
                  factor.model = factor.model, phi = phi,
                  cov.model = 'exponential', n.factors = n.factors)

y <- dat$y
X <- dat$X
X.re <- dat$X.re
coords <- dat$coords

# Package all data into a list
covs <- list(int = X[, , 1],
             abund.cov.1 = X[, , 2],
             abund.factor.1 = X.re[, , 1],
             abund.factor.2 = X.re[, , 2])
data.list <- list(y = y, covs = covs, coords = coords)
prior.list <- list(beta.comm.normal = list(mean = 0, var = 100),
                   kappa.unif = list(a = 0, b = 10),
                   phi.unif = list(a = 3 / 1, b = 3 / .1),
                   tau.sq.beta.ig = list(a = .1, b = .1))
inits.list <- list(beta.comm = 0, beta = 0, kappa = 0.5,
                   tau.sq.beta = 1, phi = 3 / 0.5)

# Small
n.batch <- 2
batch.length <- 25
n.burn <- 20
n.thin <- 1
n.chains <- 1

out <- sfMsAbund(formula = ~ abund.cov.1 + (1 | abund.factor.1) + 
                             (1 | abund.factor.2),
                 data = data.list,
                 n.batch = n.batch,
                 inits = inits.list,
                 priors = prior.list,
                 NNGP = TRUE,
                 cov.model = 'exponential',
                 n.neighbors = 5,
                 n.factors = n.factors,
                 batch.length = batch.length,
                 n.omp.threads = 3,
                 verbose = TRUE,
                 n.report = 1,
                 n.burn = n.burn,
                 n.thin = n.thin,
                 n.chains = n.chains)
summary(out)

}
