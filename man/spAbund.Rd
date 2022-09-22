\name{spAbund}
\alias{spAbund}
\title{Function for Fitting Single-Species Spatial Abundance GLMs}

\description{
  The function \code{spAbund} fits single-species spatial abundance GLMs.  
}

\usage{
spAbund(formula, data, inits, priors, tuning,
       cov.model = 'exponential', NNGP = TRUE, 
       n.neighbors = 15, search.type = 'cb',
       n.batch, batch.length, accept.rate = 0.43, family = 'NB',
       n.omp.threads = 1, verbose = TRUE, n.report = 100, 
        n.burn = round(.10 * n.batch * batch.length), n.thin = 1, 
       n.chains = 1, k.fold, k.fold.threads = 1, k.fold.seed = 100, 
       k.fold.only = FALSE, ...)
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
    given site. \code{coords} is a \eqn{J \times 2}{J x 2} matrix of the observation coordinates. 
    Note that \code{spAbundance} assumes coordinates are specified 
    in a projected coordinate system.}
  
  \item{inits}{a list with each tag corresponding to a parameter name.
    Valid tags are \code{beta}, \code{sigma.sq}, 
    \code{phi}, \code{w}, \code{nu}, \code{kappa}, \code{sigma.sq.mu}. 
    \code{nu} is only specified if \code{cov.model = "matern"}, \code{sigma.sq.mu}
    is only specified if there are random effects in \code{formula}, and
    \code{kappa} is only specified when \code{family = 'NB'}. 
    The value portion of each tag is the parameter's initial value. See \code{priors}
    description for definition of each parameter name.
    Additionally, the tag \code{fix} can be set to \code{TRUE} 
    to fix the starting values across all chains. If \code{fix} is not specified
    (the default), starting values are varied randomly across chains.}
  
  \item{priors}{a list with each tag corresponding to a parameter name. 
    Valid tags are \code{beta.normal}, \code{phi.unif}, 
    \code{sigma.sq.ig}, \code{nu.unif}, \code{kappa.unif},
    \code{sigma.sq.mu.ig}. Abundance 
    (\code{beta}) regression coefficients 
    are assumed to follow a normal distribution. The hyperparameters of the 
    normal distribution are passed as a list of length two with the first
    and second elements corresponding to the mean and variance of the normal
    distribution, which are each specified as vectors of 
    length equal to the number of coefficients to be estimated or of length
    one if priors are the same for all coefficients. If not
    specified, prior means are set to 0 and prior variances are set to 100. The 
    spatial variance parameter, \code{sigma.sq}, is assumed to follow an 
    inverse-Gamma distribution. The spatial decay \code{phi}, spatial 
    smoothness \code{nu}, and negative binomial dispersion \code{kappa}
    parameters are assumed to follow Uniform 
    distributions. The hyperparameters of the inverse-Gamma for \code{sigma.sq} 
    are passed as a vector of length two, with the first and second 
    elements corresponding to the \emph{shape} and \emph{scale}, respectively. 
    The hyperparameters of the Uniform are also passed as a vector of 
    length two with the first and second elements corresponding to 
    the lower and upper support, respectively. \code{sigma.sq.mu}
    are the random effect variances for any random effects, and are assumed to follow an 
    inverse-Gamma distribution. The hyperparameters of
    the inverse-Gamma distribution are passed as a list of length two with the 
    first and second elements corresponding to the shape and scale parameters, 
    respectively, which are each specified as vectors of length equal to the 
    number of random intercepts or of length one if priors are the same for all
    random effect variances.}
 
  \item{cov.model}{a quoted keyword that specifies the covariance
    function used to model the spatial dependence structure among the
    observations.  Supported covariance model key words are:
    \code{"exponential"}, \code{"matern"}, \code{"spherical"}, and
    \code{"gaussian"}.}

  \item{tuning}{a list with each tag corresponding to a parameter
    name. Valid tags are \code{phi}, \code{nu}, and \code{kappa}. 
    The value portion of each tag defines the initial variance of the 
    Adaptive sampler. See Roberts and Rosenthal (2009) for details.}
  
  \item{NNGP}{if \code{TRUE}, model is fit with an NNGP. If \code{FALSE}, 
    a full Gaussian process is used. See Datta et al. (2016) and 
    Finley et al. (2019) for more information. Currently only NNGP is supported,
    functionality for a full GP will be addded in future package development.}
  
  \item{n.neighbors}{number of neighbors used in the NNGP. Only used if 
  \code{NNGP = TRUE}. Datta et al. (2016) showed that 15 neighbors is usually 
  sufficient, but that as few as 5 neighbors can be adequate for certain data
  sets, which can lead to even greater decreases in run time. We recommend
  starting with 15 neighbors (the default) and if additional gains in computation
  time are desired, subsequently compare the results with a smaller number
  of neighbors using WAIC or k-fold cross-validation.}
  
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

  Finley, A. O., Datta, A., and Banerjee, S. (2020). spNNGP R package 
  for nearest neighbor Gaussian process models. \emph{arXiv} preprint arXiv:2001.09111.

  Hooten, M. B., and Hobbs, N. T. (2015). A guide to Bayesian model 
  selection for ecologists. \emph{Ecological Monographs}, 85(1), 3-28.

  Polson, N.G., J.G. Scott, and J. Windle. (2013) Bayesian Inference for
  Logistic Models Using Polya-Gamma Latent Variables. 
  \emph{Journal of the American Statistical Association}, 108:1339-1349.

  Roberts, G.O. and Rosenthal J.S. (2009) Examples  of adaptive MCMC. 
  \emph{Journal of Computational and Graphical Statistics}, 18(2):349-367.
}

\author{
  Jeffrey W. Doser \email{doserjef@msu.edu}, \cr
  Andrew O. Finley \email{finleya@msu.edu}
}

\value{
  An object of class \code{spAbund} that is a list comprised of: 

  \item{beta.samples}{a \code{coda} object of posterior samples
    for the abundance regression coefficients.}

  \item{kappa.samples}{a \code{coda} object of posterior samples
    for the abundance dispersion parameter. Only included when
    \code{family = 'NB'}.}

  \item{y.rep.samples}{a three-dimensional object of posterior samples 
    for the abundance replicate (fitted) values with dimensions 
    corresponding to MCMC samples, site, and replicate.}

  \item{mu.samples}{a three-dimensional array of posterior samples
    for the expected abundance samples with dimensions corresponding 
    to MCMC samples, site, and replicate.}

  \item{theta.samples}{a \code{coda} object of posterior samples
    for spatial covariance parameters.}

  \item{w.samples}{a \code{coda} object of posterior samples
    for latent spatial random effects.}

  \item{sigma.sq.mu.samples}{a \code{coda} object of posterior samples
    for variances of random intercepts included in the model. 
    Only included if random intercepts are specified in 
    \code{formula}.}

  \item{beta.star.samples}{a \code{coda} object of posterior samples
    for the abundance random effects. Only included if random intercepts
    are specified in \code{formula}.}

  \item{like.samples}{a \code{coda} object of posterior samples
    for the likelihood value associated with each site. Used for calculating
    WAIC.}

  \item{rhat}{a list of Gelman-Rubin diagnostic values for some of the model
    parameters.}

  \item{ESS}{a list of effective sample sizes for some of the model parameters.}

  \item{run.time}{execution time reported using \code{proc.time()}.}

  \item{k.fold.deviance}{soring rule (deviance) from k-fold cross-validation. 
    Only included if \code{k.fold} is specified in function call.}

  The return object will include additional objects used for 
  subsequent prediction and/or model fit evaluation. Note that detection
  probability values are not included in the model object, but can be 
  extracted using \code{fitted()}.
}

\examples{
set.seed(888)
J.x <- 8
J.y <- 8
J <- J.x * J.y
n.rep <- sample(3, J, replace = TRUE)
beta <- c(0, -1.5, 0.3, -0.8)
p.abund <- length(beta)
mu.RE <- list(levels = c(50, 45),
              sigma.sq.mu = c(1.3, 0.5), 
              beta.indx = c(1, 2))
phi <- 3/.6
sigma.sq <- 2
kappa <- 0.2
sp <- TRUE 
cov.model <- 'exponential'
family <- 'NB'
dat <- simAbund(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, 
	        kappa = kappa, mu.RE = mu.RE, sp = sp, phi = phi, 
		sigma.sq = sigma.sq, cov.model = cov.model, family = 'NB')

y <- dat$y
X <- dat$X
X.re <- dat$X.re
coords <- dat$coords

abund.covs <- list(int = X[, , 1], 
		   abund.cov.1 = X[, , 2], 
		   abund.cov.2 = X[, , 3], 
		   abund.cov.3 = X[, , 4],
		   abund.factor.1 = X.re[, , 1], 
		   abund.factor.2 = X.re[, , 2])

data.list <- list(y = y, 
		  covs = abund.covs, 
                  coords = coords)

# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 100),
		   phi.unif = c(3 / 1, 3 / .1),
		   sigma.sq.ig = c(2, 1),
                   kappa.unif = c(0.001, 10)) 
# Starting values
inits.list <- list(beta = beta,
		   kappa = kappa, 
                   sigma.sq = sigma.sq,
                   phi = phi)

tuning <- list(phi = 0.3, kappa = 0.05) 
n.batch <- 4
batch.length <- 25
n.burn <- 20 
n.thin <- 1
n.chains <- 1

out <- spAbund(formula = ~ abund.cov.1 + abund.cov.2 + abund.cov.3 +
	                   (1 | abund.factor.1) + (abund.cov.1 | abund.factor.2),
	       data = data.list, 
	       n.batch = n.batch, 
	       batch.length = batch.length, 
	       inits = inits.list, 
	       tuning = tuning,
	       priors = prior.list, 
	       NNGP = TRUE, 
	       cov.model = 'exponential',
	       search.type = 'cb',
	       n.neighbors = 5,
	       accept.rate = 0.43, 
	       n.omp.threads = 1, 
	       verbose = TRUE, 
	       n.report = 1,
	       n.burn = n.burn,
	       n.thin = n.thin,
	       n.chains = n.chains) 
}
