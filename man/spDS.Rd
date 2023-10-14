\name{spDS}
\alias{spDS}
\title{Function for Fitting Single-Species Spatially-Explicit Hierarchical Distance Sampling Models}

\usage{
spDS(abund.formula, det.formula, data, inits, priors, tuning,
     cov.model = 'exponential', NNGP = TRUE,
     n.neighbors = 15, search.type = 'cb',
     n.batch, batch.length, accept.rate = 0.43, family = 'Poisson',
     transect = 'line', det.func = 'halfnormal',
     n.omp.threads = 1, verbose = TRUE,
     n.report = 100, n.burn = round(.10 * n.batch * batch.length), n.thin = 1, 
     n.chains = 1, ...)
}

\description{
  Function for fitting single-sepcies spatially-explicit hierarchical distance sampling models. Spatial models are fit using Nearest Neighbor Gaussian Processes.
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
  Valid tags are \code{y}, \code{covs}, \code{coords}, \code{dist.breaks}, and \code{offset}. \code{y}
  is a matrix or data frame of the observed count values, 
  with first dimension equal to the number of 
  sites (\eqn{J}{J}) and second dimension equal to the number of 
  distance bins. \code{covs} is a matrix or data frame
  containing the variables used in the abundance and/or the detection
  portion of the model, with
  \eqn{J}{J} rows for each column (variable). \code{dist.breaks} is a vector of 
  distances that denote the breakpoints of the distance bands. \code{dist.breaks} should
  have length equal to the number of columns in \code{y} plus one. \code{offset} is an 
  offset that can be used to scale estimates from abundance per transect to density per 
  some desired unit of measure. This can be either a single value or a vector with an offset
  value for each site (e.g., if transects differ in length). 
  \code{coords} is a \eqn{J \times 2}{J x 2} matrix of the observation coordinates. 
  Note that \code{spAbundance} assumes coordinates are specified 
  in a projected coordinate system. }

\item{inits}{a list with each tag corresponding to a parameter name.
  Valid tags are \code{N}, \code{beta}, \code{alpha}, \code{kappa}, 
  \code{sigma.sq}, \code{phi}, \code{w}, \code{nu},
  \code{sigma.sq.mu}, and \code{sigma.sq.p}. The value portion of each tag is the 
  parameter's initial value. \code{sigma.sq.mu} and \code{sigma.sq.p} are 
  only relevant when including random effects in the abundance and 
  detection portion of the abundance model, respectively. \code{kappa} is 
  only relevant when \code{family = 'NB'}. \code{nu} is only 
  specified if \code{cov.model = "matern"}. See \code{priors} 
  description for definition of each parameter name.
  Additionally, the tag \code{fix} can be set to \code{TRUE} 
  to fix the starting values across all chains. If \code{fix} is not specified
  (the default), starting values are varied randomly across chains.}

\item{priors}{a list with each tag corresponding to a parameter name. 
  Valid tags are \code{beta.normal}, \code{alpha.normal}, \code{kappa.unif}, 
  \code{phi.unif}, \code{sigma.sq.ig}, \code{nu.unif}, 
  \code{sigma.sq.mu.ig}, and \code{sigma.sq.p.ig}. 
  Abundance (\code{beta}) and detection (\code{alpha}) 
  regression coefficients are assumed to follow a normal distribution. 
  The hyperparameters of the normal distribution are passed as a list of 
  length two with the first and second elements corresponding to the mean 
  and variance of the normal distribution, which are each specified as vectors of 
  length equal to the number of coefficients to be estimated or of length one
  if priors are the same for all coefficients. If not specified, prior means 
  are set to 0 and prior variances set to 100. The spatial variance parameter, 
  \code{sigma.sq}, is assumed to follow an 
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
  detection random effects, respectively, and are assumed to follow an inverse
  Gamma distribution. The hyperparameters of the inverse-Gamma distribution
  are passed as a list of length two with first and second elements corresponding
  to the shape and scale parameters, respectively, which are each specified as
  vectors of length equal to the number of random intercepts/slopes or of length one 
  if priors are the same for all random effect variances.}

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

  Datta, A., S. Banerjee, A.O. Finley, and A.E. Gelfand. (2016)
  Hierarchical Nearest-Neighbor Gaussian process models for large
  geostatistical datasets. \emph{Journal of the American Statistical
  Association}, \doi{10.1080/01621459.2015.1044091}.
  
  Finley, A.O., A. Datta, B.D. Cook, D.C. Morton, H.E. Andersen, and
  S. Banerjee. (2019) Efficient algorithms for Bayesian Nearest Neighbor
  Gaussian Processes. \emph{Journal of Computational and Graphical
  Statistics}, \doi{10.1080/10618600.2018.1537924}.

  Royle, J. A., Dawson, D. K., & Bates, S. (2004). Modeling 
  abundance effects in distance sampling. Ecology, 85(6), 1591-1597.
}

\author{
  Jeffrey W. Doser \email{doserjef@msu.edu}, \cr
  Andrew O. Finley \email{finleya@msu.edu}
}

\value{
  An object of class \code{spDS} that is a list comprised of: 

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

  \item{theta.samples}{a \code{coda} object of posterior samples
    for spatial covariance parameters.}

  \item{w.samples}{a \code{coda} object of posterior samples
    for latent spatial random effects.}

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
sp <- TRUE
phi <- 3 / .5
sigma.sq <- 0.8
cov.model <- 'exponential'
family <- 'NB'
kappa <- 0.1
offset <- 1.8
transect <- 'point'

dat <- simDS(J.x = J.x, J.y = J.y, n.bins = n.bins, bin.width = bin.width,
             beta = beta, alpha = alpha, det.func = det.func, kappa = kappa, 
             mu.RE = mu.RE, p.RE = p.RE, sp = sp,
             offset = offset, transect = transect, phi = phi, sigma.sq = sigma.sq,
             cov.model = cov.model)

y <- dat$y
X <- dat$X
X.re <- dat$X.re
X.p <- dat$X.p
X.p.re <- dat$X.p.re
dist.breaks <- dat$dist.breaks
coords <- dat$coords

covs <- cbind(X, X.p)
colnames(covs) <- c('int.abund', 'abund.cov.1', 'abund.cov.2', 'abund.cov.3', 
                    'int.det', 'det.cov.1')

data.list <- list(y = y, 
                  covs = covs,
                  dist.breaks = dist.breaks, 
                  coords = coords,
                  offset = offset)

# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 10),
                   alpha.normal = list(mean = 0,
                                       var = 10), 
                   kappa.unif = c(0, 100), 
                   phi.unif = c(3 / 1, 3 / .1),
                   sigma.sq.ig = c(2, 1)) 
# Starting values
inits.list <- list(alpha = 0,
                   beta = 0,
                   kappa = 1, 
                   phi = 3 / .5, 
                   sigma.sq = 1)
# Tuning values
tuning <- list(beta = 0.1, alpha = 0.1, beta.star = 0.3, alpha.star = 0.1, 
               kappa = 0.2, phi = 1, w = 1) 

out <- spDS(abund.formula = ~ abund.cov.1 + abund.cov.2 + abund.cov.3,
            det.formula = ~ det.cov.1,
            data = data.list, 
            n.batch = 10, 
            batch.length = 25, 
            inits = inits.list, 
            family = 'NB',
            det.func = 'halfnormal', 
            transect = 'point', 
            cov.model = 'exponential', 
            NNGP = TRUE,
	    n.neighbors = 5,
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
