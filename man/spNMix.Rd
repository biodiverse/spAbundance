\name{spNMix}
\alias{spNMix}
\title{Function for Fitting Single-Species Spatial N-Mixture Models}

\description{
  The function \code{spNMix} fits single-species spatial N-mixture models. Spatial models are fit using Nearest Neighbor Gaussian Processes. 
}

\usage{
spNMix(abund.formula, det.formula, data, inits, priors, tuning,
       cov.model = 'exponential', NNGP = TRUE, 
       n.neighbors = 15, search.type = 'cb',
       n.batch, batch.length, accept.rate = 0.43, family = 'Poisson',
       n.omp.threads = 1, verbose = TRUE, n.report = 100, 
       n.burn = round(.10 * n.batch * batch.length), n.thin = 1, 
       n.chains = 1, ...)
}

\arguments{

  \item{abund.formula}{a symbolic description of the model to be fit
    for the abundance portion of the model using R's model syntax. Only
    right-hand side of formula is specified. See example below. Random intercepts
    and random slopes are allowed using lme4 syntax (Bates et al. 2015).}
  
  \item{det.formula}{a symbolic description of the model to be fit
  for the detection portion of the model using R's model syntax. Only
  right-hand side of formula is specified. See example below. Random intercepts
  and random slopes are allowed using lme4 syntax (Bates et al. 2015).}

  \item{data}{a list containing data necessary for model fitting.
    Valid tags are \code{y}, \code{abund.covs}, \code{det.covs}, \code{offset}, and \code{coords}. 
    \code{y} is the count data matrix or data frame with 
    first dimension equal to the number of sites (\eqn{J}{J}) and second 
    dimension equal to the maximum number of replicates at a given site. 
    \code{abund.covs} is a matrix or data frame containing the variables used 
    in the abundance portion of the model, with \eqn{J}{J} rows for each column 
    (variable). \code{det.covs} is a list of variables included in the 
    detection portion of the model. Each list element is a different detection 
    covariate, which can be site-level or observational-level. Site-level 
    covariates are specified as a vector of length \eqn{J}{J} while 
    observation-level covariates are specified 
    as a matrix or data frame with the number of rows equal to \eqn{J}{J} and 
    number of columns equal to the maximum number of replicates at a given site. 
    \code{coords} is a \eqn{J \times 2}{J x 2} matrix of the observation coordinates. 
    Note that \code{spAbundance} assumes coordinates are specified 
    in a projected coordinate system. \code{offset} is an offset to use in 
    the abundance model (e.g., an area offset). This can be either a single value or a 
    vector with an offset for each site (e.g., if survey area differed in size).}
  
  \item{inits}{a list with each tag corresponding to a parameter name.
    Valid tags are \code{N}, \code{beta}, \code{alpha}, \code{sigma.sq}, 
    \code{phi}, \code{w}, \code{nu}, \code{kappa}, \code{sigma.sq.mu}, \code{sigma.sq.p}. 
    \code{nu} is only specified if \code{cov.model = "matern"}, \code{sigma.sq.p}
    is only specified if there are random effects in \code{det.formula}, \code{sigma.sq.mu}
    is only specified if there are random effects in \code{abund.formula}, and
    \code{kappa} is only specified when \code{family = 'NB'}. 
    The value portion of each tag is the parameter's initial value. See \code{priors}
    description for definition of each parameter name.
    Additionally, the tag \code{fix} can be set to \code{TRUE} 
    to fix the starting values across all chains. If \code{fix} is not specified
    (the default), starting values are varied randomly across chains.}
  
  \item{priors}{a list with each tag corresponding to a parameter name. 
    Valid tags are \code{beta.normal}, \code{alpha.normal}, \code{phi.unif}, 
    \code{sigma.sq.ig}, \code{nu.unif}, \code{kappa.unif},
    \code{sigma.sq.mu.ig}, and \code{sigma.sq.p.ig}. Abundance 
    (\code{beta}) and detection (\code{alpha}) regression coefficients 
    are assumed to follow a normal distribution. The hyperparameters of the 
    normal distribution are passed as a list of length two with the first
    and second elements corresponding to the mean and variance of the normal
    distribution, which are each specified as vectors of 
    length equal to the number of coefficients to be estimated or of length
    one if priors are the same for all coefficients. If not
    specified, prior means are set to 0 and prior variances for abundance coefficients
    are set to 100 and for detection coefficients set to 2.72. The 
    spatial variance parameter, \code{sigma.sq}, is assumed to follow an 
    inverse-Gamma distribution. The spatial decay \code{phi}, spatial 
    smoothness \code{nu}, and negative binomial dispersion \code{kappa}
    parameters are assumed to follow Uniform 
    distributions. The hyperparameters of the inverse-Gamma for \code{sigma.sq} 
    are passed as a vector of length two, with the first and second 
    elements corresponding to the \emph{shape} and \emph{scale}, respectively. 
    The hyperparameters of the Uniform are also passed as a vector of 
    length two with the first and second elements corresponding to 
    the lower and upper support, respectively. \code{sigma.sq.mu} and 
    \code{sigma.sq.p} are the random effect variances for any abundance or 
    detection random effects, respectively, and are assumed to follow an 
    inverse-Gamma distribution. The hyperparameters of
    the inverse-Gamma distribution are passed as a list of length two with the 
    first and second elements corresponding to the shape and scale parameters, 
    respectively, which are each specified as vectors of length equal to the 
    number of random intercepts/slopes or of length one if priors are the same for all
    random effect variances.}
 
  \item{cov.model}{a quoted keyword that specifies the covariance
    function used to model the spatial dependence structure among the
    observations. Supported covariance model key words are:
    \code{"exponential"}, \code{"matern"}, \code{"spherical"}, and
    \code{"gaussian"}.}

  \item{tuning}{a single numeric value representing the initial variance of the
  adaptive sampler for \code{beta}, \code{alpha}, \code{beta.star} (the abundance 
  random effect values), \code{alpha.star} (the detection random effect values),
  \code{kappa}, \code{phi}, \code{nu}, and \code{w}. See Roberts and Rosenthal (2009) for details.}

  \item{NNGP}{if \code{TRUE}, model is fit with an NNGP.  See Datta et al. (2016) and 
    Finley et al. (2019) for more information. Currently only NNGP is supported,
    functionality for a Gaussian Process may be addded in future package development.}
  
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
    e.g., if data are on a grid.}
 
  \item{n.batch}{the number of MCMC batches in each chain to run for the Adaptive MCMC 
    sampler. See Roberts and Rosenthal (2009) for details.}
  
  \item{batch.length}{the length of each MCMC batch in each chain to run for the Adaptive 
    MCMC sampler. See Roberts and Rosenthal (2009) for details.}
  
  \item{accept.rate}{target acceptance rate for Adaptive MCMC. Default is 
    0.43. See Roberts and Rosenthal (2009) for details.}

  \item{family}{the distribution to use for the latent abundance process. Currently 
    supports \code{'NB'} (negative binomial) and \code{'Poisson'}.}
  
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
    MCMC progress.}

  \item{n.burn}{the number of samples out of the total \code{n.batch * batch.length} 
    samples in each chain to discard as burn-in. By default, the first 
    10\% of samples is discarded.}
  
  \item{n.thin}{the thinning interval for collection of MCMC samples. The
    thinning occurs after the \code{n.burn} samples are discarded. Default 
    value is set to 1.}

  \item{n.chains}{the number of MCMC chains to run in sequence.}

  \item{...}{currently no additional arguments}
}

\references{

  Bates, Douglas, Martin Maechler, Ben Bolker, Steve Walker (2015).
  Fitting Linear Mixed-Effects Models Using lme4. Journal of
  Statistical Software, 67(1), 1-48. \doi{10.18637/jss.v067.i01}.

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

  Royle, J. A. (2004). N‐mixture models for estimating population size 
  from spatially replicated counts. Biometrics, 60(1), 108-115.
}

\author{
  Jeffrey W. Doser \email{doserjef@msu.edu}, \cr
  Andrew O. Finley \email{finleya@msu.edu}
}

\value{
  An object of class \code{spNMix} that is a list comprised of: 

  \item{beta.samples}{a \code{coda} object of posterior samples
    for the abundance regression coefficients.}

  \item{alpha.samples}{a \code{coda} object of posterior samples
    for the detection regression coefficients.}

  \item{kappa.samples}{a \code{coda} object of posterior samples
    for the abundance dispersion parameter. Only included when
    \code{family = 'NB'}.}

  \item{N.samples}{a \code{coda} object of posterior samples 
    for the latent abundance values}

  \item{mu.samples}{a \code{coda} object of posterior samples
    for the latent expected abundance values}

  \item{theta.samples}{a \code{coda} object of posterior samples
    for spatial covariance parameters.}

  \item{w.samples}{a \code{coda} object of posterior samples
    for latent spatial random effects.}

  \item{sigma.sq.mu.samples}{a \code{coda} object of posterior samples
    for variances of random intercepts/slopes included in the abundance portion
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

  \item{run.time}{execution time reported using \code{proc.time()}.}

  The return object will include additional objects used for 
  subsequent prediction and/or model fit evaluation. Note that detection
  probability values are not included in the model object, but can be 
  extracted using \code{fitted()}.
}

\examples{
set.seed(350)
# Simulate Data -----------------------------------------------------------
J.x <- 15
J.y <- 15 
J <- J.x * J.y
n.rep <- sample(3, J, replace = TRUE)
beta <- c(0.5, 1.5)
p.abund <- length(beta)
alpha <- c(0.5, 1.2, -0.5)
p.det <- length(alpha)
mu.RE <- list()
p.RE <- list()
phi <- runif(1, 3 / 1, 3 / .1)
sigma.sq <- runif(1, 0.2, 1.5) 
kappa <- 0.5
sp <- TRUE 
cov.model <- 'exponential'
dat <- simNMix(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, alpha = alpha,
               kappa = kappa, mu.RE = mu.RE, p.RE = p.RE, sp = sp, 
               phi = phi, sigma.sq = sigma.sq, cov.model = cov.model, 
               family = 'NB')

y <- dat$y
X <- dat$X
X.re <- dat$X.re
X.p <- dat$X.p
X.p.re <- dat$X.p.re
coords <- dat$coords

abund.covs <- X
colnames(abund.covs) <- c('int', 'abund.cov.1')

det.covs <- list(det.cov.1 = X.p[, , 2], 
                 det.cov.2 = X.p[, , 3]) 

data.list <- list(y = y, 
                  abund.covs = abund.covs,
                  det.covs = det.covs, 
                  coords = coords)

# Priors
prior.list <- list(beta.normal = list(mean = rep(0, p.abund), 
                                      var = rep(100, p.abund)),
                   alpha.normal = list(mean = rep(0, p.det),
                                       var = rep(2.72, p.det)), 
                   kappa.unif = c(0, 10)) 
# Starting values
inits.list <- list(alpha = alpha,
                   beta = beta,
                   kappa = kappa, 
                   w = rep(0, J),
                   phi = 3 / 0.5,
                   sigma.sq = 1,
                   N = apply(y, 1, max, na.rm = TRUE))
# Tuning values 
tuning.list <- list(phi = 0.5, kappa = 0.5, beta = 0.1, alpha = 0.1, w = 0.1)


n.batch <- 4
batch.length <- 25
n.burn <- 0
n.thin <- 1
n.chains <- 1

out <- spNMix(abund.formula = ~ abund.cov.1,
              det.formula = ~ det.cov.1 + det.cov.2, 
              data = data.list, 
              n.batch = n.batch, 
              batch.length = batch.length, 
              inits = inits.list, 
              priors = prior.list, 
              NNGP = TRUE,
              cov.model = 'spherical',
              n.neighbors = 10,
              accept.rate = 0.43, 
              n.omp.threads = 1, 
              verbose = TRUE, 
              n.report = 1,
              n.burn = n.burn,
              n.thin = n.thin,
              n.chains = n.chains)
summary(out)
}
