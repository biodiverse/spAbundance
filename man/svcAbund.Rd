\name{svcAbund}
\alias{svcAbund}
\title{Function for Fitting Univariate Spatialy-Varying Coefficient GLMMs}

\description{
  The function \code{svcAbund} fits univariate spatially-varying coefficient GLMMs.  
}

\usage{
svcAbund(formula, data, inits, priors, tuning,
         svc.cols = 1, cov.model = 'exponential', NNGP = TRUE, 
         n.neighbors = 15, search.type = 'cb', n.batch,
         batch.length, accept.rate = 0.43, family = 'Gaussian',
         n.omp.threads = 1, verbose = TRUE, n.report = 100, 
         n.burn = round(.10 * n.batch * batch.length), n.thin = 1, 
         n.chains = 1, ...)
}

\arguments{

  \item{formula}{a symbolic description of the model to be fit
    for the model using R's model syntax. Only
    right-hand side of formula is specified. See example below. Random intercepts
    and slopes are allowed using lme4 syntax (Bates et al. 2015).}
  
  \item{data}{a list containing data necessary for model fitting.
    Valid tags are \code{y}, \code{covs}, \code{coords}, and \code{z} (for 
    \code{family = 'zi-Gaussian'} only. \code{y}
    is a vector of the observed count values, where the values 
    represent the observed values at each site. \code{covs} is a list, matrix, or data
    frame of covariates used in the model, where each column (or list element) represents a 
    different covariate. \code{coords} is a \eqn{J \times 2}{J x 2} matrix of the observation coordinates.
    Note that \code{svcAbundance} assumes coordinates are specified 
    in a projected coordinate system. \code{z} is used for fitting a zero-inflated 
    Gaussian model. It is a vector where each value indicates the binary 
    component of the model. In the context of abundance models, this can be 
    thought of as the component of the model that indicates whether the species is 
    present at each location, and then the supplied values in \code{y} are the 
    observed abundance values at those locations where \code{z = 1}.}
  
  \item{inits}{a list with each tag corresponding to a parameter name.
    Valid tags are \code{beta}, \code{sigma.sq}, 
    \code{phi}, \code{w}, \code{nu}, \code{tau.sq}, \code{sigma.sq.mu}. 
    \code{nu} is only specified if \code{cov.model = "matern"}, \code{sigma.sq.mu}
    is only specified if there are random effects in \code{formula}, and
    The value portion of each tag is the parameter's initial value. See \code{priors}
    description for definition of each parameter name.
    Additionally, the tag \code{fix} can be set to \code{TRUE} 
    to fix the starting values across all chains. If \code{fix} is not specified
    (the default), starting values are varied randomly across chains.}
  
  \item{priors}{a list with each tag corresponding to a parameter name. 
    Valid tags are \code{beta.normal}, \code{phi.unif}, 
    \code{sigma.sq.ig}, \code{nu.unif}, \code{tau.sq.ig},
    \code{sigma.sq.mu.ig}. Abundance (\code{beta}) regression coefficients 
    are assumed to follow a normal distribution. The hyperparameters of the 
    normal distribution are passed as a list of length two with the first
    and second elements corresponding to the mean and variance of the normal
    distribution, which are each specified as vectors of 
    length equal to the number of coefficients to be estimated or of length
    one if priors are the same for all coefficients. If not
    specified, prior means are set to 0 and prior variances are set to 100. The 
    spatial variance parameter, \code{sigma.sq}, and the Gaussian 
    residual variance parameter, \code{tau.sq}, are assumed to follow an 
    inverse-Gamma distribution. The spatial decay \code{phi} and spatial 
    smoothness \code{nu}, parameters are assumed to follow Uniform 
    distributions. The hyperparameters of the inverse-Gamma for \code{sigma.sq} 
    is passed as a list of length two with the first and second elements corresponding
    to the shape and scale parameters of the inverse-Gamma distribution either for
    each spatially-varying coefficient, or a single value if assumign the same values
    for all spatially-varying coefficients. The hyperparameters of the inverse-Gamma for 
    \code{tau.sq} is passed as a vector of length two, with the first and second 
    elements corresponding to the \emph{shape} and \emph{scale}, respectively. 
    The hyperparameters of the Uniform are also passed as a list of 
    length two with the first and second elements corresponding to 
    the lower and upper support, respectively, for each SVC or a single value
    if giving the same prior for each SVC. \code{sigma.sq.mu}
    are the random effect variances for any random effects, and are assumed to follow an 
    inverse-Gamma distribution. The hyperparameters of
    the inverse-Gamma distribution are passed as a list of length two with the 
    first and second elements corresponding to the shape and scale parameters, 
    respectively, which are each specified as vectors of length equal to the 
    number of random effects or of length one if priors are the same for all
    random effect variances.}

  \item{svc.cols}{a vector indicating the variables whose effects will be
    estimated as spatially-varying coefficients. \code{svc.cols} can be an
    integer vector with values indicating the order of covariates specified
    in the model formula (with 1 being the intercept if specified), or it can
    be specified as a character vector with names corresponding to variable
    names in \code{occ.covs} (for the intercept, use '(Intercept)'). \code{svc.cols}
    default argument of 1 results in a univariate spatial GLMM analogous to
    \code{spAbund} (assuming an intercept is included in the model).}

  \item{cov.model}{a quoted keyword that specifies the covariance
    function used to model the spatial dependence structure among the
    observations.  Supported covariance model key words are:
    \code{"exponential"}, \code{"matern"}, \code{"spherical"}, and
    \code{"gaussian"}.}

  \item{tuning}{a single numeric value representing the initial variance of the
  adaptive sampler for \code{phi} and \code{nu}. See Roberts and Rosenthal (2009) for details.}

  \item{NNGP}{if \code{TRUE}, model is fit with an NNGP. See Datta et al. (2016) and 
    Finley et al. (2019) for more information. Currently only NNGP is supported,
    functionality for a full GP may be addded in future package development.}
  
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
 
  \item{n.batch}{the number of MCMC batches in each chain to run for the adaptive MCMC 
    sampler. See Roberts and Rosenthal (2009) for details.}
  
  \item{batch.length}{the length of each MCMC batch in each chain to run for the adaptive 
    MCMC sampler. See Roberts and Rosenthal (2009) for details.}
  
  \item{accept.rate}{target acceptance rate for adaptive MCMC. Default is 
    0.43. See Roberts and Rosenthal (2009) for details.}

  \item{family}{the distribution to use for abundance. Currently, spatially-varying 
    coefficient models are available for \code{family = 'Gaussian'} and 
    \code{family = 'zi-Gaussian'}.}
  
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
}

\author{
  Jeffrey W. Doser \email{doserjef@msu.edu}, \cr
  Andrew O. Finley \email{finleya@msu.edu}
}

\value{
  An object of class \code{svcAbund} that is a list comprised of: 

  \item{beta.samples}{a \code{coda} object of posterior samples
    for the abundance regression coefficients.}

  \item{tau.sq.samples}{a \code{coda} object of posterior samples
    for the residual variance parameter.}

  \item{y.rep.samples}{a \code{coda} object of posterior samples 
    for the abundance replicate (fitted) values with dimensions 
    corresponding to MCMC samples and site.}

  \item{mu.samples}{a \code{coda} object of posterior samples
    for the expected abundance samples with dimensions corresponding 
    to MCMC samples and site.}

  \item{theta.samples}{a \code{coda} object of posterior samples
    for spatial covariance parameters.}

  \item{w.samples}{a three-dimensional array of posterior samples
    for the spatially-varying coefficients with dimensions corresponding
    to MCMC sample, SVC, and site.}

  \item{sigma.sq.mu.samples}{a \code{coda} object of posterior samples
    for variances of random effects included in the model. 
    Only included if random effects are specified in 
    \code{formula}.}

  \item{beta.star.samples}{a \code{coda} object of posterior samples
    for the abundance random effects. Only included if random effects 
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
set.seed(1000)
# Sites
J.x <- 10
J.y <- 10
J <- J.x * J.y
# Abundance ---------------------------
beta <- c(5, 0.5, -0.2, 0.75)
p <- length(beta)
mu.RE <- list()
mu.RE <- list(levels = c(35, 40),
              sigma.sq.mu = c(0.7, 1.5),
              beta.indx = list(1, 1))
# Spatial parameters ------------------
sp <- TRUE
svc.cols <- c(1, 2)
p.svc <- length(svc.cols)
cov.model <- "exponential"
sigma.sq <- runif(p.svc, 0.4, 4)
phi <- runif(p.svc, 3/1, 3/0.6)
tau.sq <- 2
z <- rbinom(J, 1, 0.5)

# Get all the data
dat <- simAbund(J.x = J.x, J.y = J.y, beta = beta, tau.sq = tau.sq,
                mu.RE = mu.RE, sp = sp, svc.cols = svc.cols, 
                family = 'zi-Gaussian', cov.model = cov.model, 
                sigma.sq = sigma.sq, phi = phi, z = z)
# Get data in format for spAbundance --------------------------------------
y <- dat$y
X <- dat$X
X.re <- dat$X.re
coords <- dat$coords

# Package all data into a list
covs <- cbind(X, X.re)
colnames(covs) <- c('int', 'cov.1', 'cov.2', 'cov.3', 'factor.1', 'factor.2')

# Data list bundle
data.list <- list(y = y, covs = covs, coords = coords, z = z)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 1000),
                   sigma.sq.ig = list(a = 2, b = 1), tau.sq = c(2, 1),
                   sigma.sq.mu.ig = list(a = 2, b = 1),
                   phi.unif = list(a = 3 / 1, b = 3 / 0.1))

# Starting values
inits.list <- list(beta = 0, alpha = 0,
                   sigma.sq = 1, phi = 3 / 0.5, 
                   tau.sq = 2, sigma.sq.mu = 0.5)
# Tuning
tuning.list <- list(phi = 1)

n.batch <- 10
batch.length <- 25
n.burn <- 100
n.thin <- 1

out <- svcAbund(formula = ~ cov.1 + cov.2 + cov.3 + 
                            (1 | factor.1) + (1 | factor.2),
                svc.cols = c(1, 2),
                data = data.list,
                n.batch = n.batch,
                batch.length = batch.length,
                inits = inits.list,
                priors = prior.list,
                accept.rate = 0.43,
                family = 'zi-Gaussian',
                cov.model = "exponential",
                tuning = tuning.list,
                n.omp.threads = 1,
                verbose = TRUE,
                NNGP = TRUE,
                n.neighbors = 5,
                n.report = 25,
                n.burn = n.burn,
                n.thin = n.thin,
                n.chains = 3)
}
