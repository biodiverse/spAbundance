\name{svcTIntAbund}
\alias{svcTIntAbund}
\title{Function for Fitting Multi-Season Univariate Spatially-Varying Coefficient GLMMs with Multiple Data Sources}

\description{
  The function \code{svcTIntAbund} is experimental and should be used with caution. \code{svcTIntAbund} fits univariate spatially-varying coefficient GLMMs with multiple data sources. The data sources are assumed to have data collected over both space and time. One could also have data that follow the robust design, where data are collected over a larger time span (seasons) and then within a season there are multiple replicate surveys of a given site. Models are fit using Nearest Neighbor Gaussian Processes.
}

\usage{
svcTIntAbund(abund.formula, det.formula, data, inits, priors, tuning,
         svc.cols = 1, cov.model = 'exponential', NNGP = TRUE,
         n.neighbors = 15, search.type = 'cb', n.batch,
         batch.length, accept.rate = 0.43, family = 'Poisson',
         n.omp.threads = 1, verbose = TRUE, n.report = 100,
         n.burn = round(.10 * n.batch * batch.length), n.thin = 1,
         n.chains = 1, save.fitted = TRUE, ...)
}

\arguments{

  \item{abund.formula}{a symbolic description of the model to be fit
    for the latent abundance component of the model using R's model syntax. Only
    right-hand side of formula is specified. See example below. Random intercepts
    (but not slopes) are allowed using lme4 syntax (Bates et al. 2015).}

 \item{det.formula}{a list of symbolic descriptions of the models to be fit
  for the observation (detection) portion of the model using R's model syntax for 
  each data set. Each element in the list is a formula for the detection 
  model of a given data set. Only right-hand side of formula is specified. 
  See example below. Random intercepts are allowed using lme4 syntax (Bates et al. 2015).}

  \item{data}{a list containing data necessary for model fitting.
    Valid tags are \code{y}, \code{abund.covs}, \code{det.covs}, \code{coords}, \code{sites}, 
    \code{seasons}, and \code{offset}. 
    \code{y} is a list of three-dimensional arrays with first dimensional equal to the 
    number of sites surveyed in that data set, second dimension equal to the
    number of primary time periods (i.e., years or seasons), and third dimension 
    equal to the maximum number of replicate surveys at a site within a given season. 
    Note that it is not required that there are multiple replicate surveys within a given season.
    \code{abund.covs} is a list of variables hypothesized to influence the latent ecological 
    process (abundance). Each list element is a different covariate, which can be site level
    or site/primary time period level. Site-level covariates are specified as a vector of 
    length \eqn{J}{J} while site/primary time period level covariates are specified as a matrix
    with rows corresponding to sites and columns corresponding to primary time periods.
    \code{det.covs} is a list of variables included in the detection 
    portion of the model for each data source. \code{det.covs} should have the same 
    number of elements as \code{y}, where each element is itself a list. 
    Each element of the list for a given data source is a different detection covariate, 
    which can be site-level , site-season-level, or observation-level. Site-level covariates
    and site/primary time period level covariates are specified in the same manner as 
    \code{abund.covs}. Observation-level
    covariates are specified as a three-dimensional array with first dimension corresponding 
    to sites, second dimension corresponding to primary time period, and third dimension
    corresponding to replicate. \code{sites} is a list of site indices with number of 
    elements equal to the number of data sources being modeled. Each element contains a vector of 
    length equal to the number of sites that specific data source contains. Each value in the 
    vector indicates the corresponding site in \code{abund.covs} covariates that corresponds 
    with the specific row of the detection-nondetection data for the data source. This is used 
    to properly link sites across data sets. Similarly, \code{seasons} is a list of season 
    indices with number of elements equal to the number of data sources being modeled. 
    Each element contains a vector of length equal to the number of seasons that a specific 
    data source is available for. This is used to properly link seasons across data sets. 
    Each value in the vector indicates the corresponding season in \code{abund.covs} covariates 
    that correspond with the specific column of the detection-nondetection data for the 
    given data source. This is used to properly link seasons across data sets, which can 
    have a differing number of seasons surveyed.  \code{coords} is a matrix of the 
    site coordinates. Note that \code{spOccupancy} assumes coordinates are specified 
    in a projected coordinate system.}

  \item{inits}{a list with each tag corresponding to a parameter name.
    Valid tags are \code{beta}, \code{alpha}, \code{sigma.sq},
    \code{phi}, \code{w}, \code{nu}, \code{sigma.sq.mu}, \code{sigma.sq.p}, \code{kappa}.
    \code{nu} is only specified if \code{cov.model = "matern"}, \code{sigma.sq.mu}
    is only specified if there are abundance random effects in \code{abund.formula}, and 
    \code{sigma.sq.p} should only be specified if there are detection random effects in 
    \code{det.formula}. The value portion of each tag is the parameter's initial value. See \code{priors}
    description for definition of each parameter name. The tag \code{alpha} is a list comprised 
    of the initial values for the detection parameters for each data source. 
    Each element of the list should be a vector of initial values for all detection 
    parameters in the given data source or a single value for each data source to 
    assign all parameters for a given data source the same initial value. See \code{priors} 
    description for definition of each parameter name.
    Additionally, the tag \code{fix} can be set to \code{TRUE}
    to fix the starting values across all chains. If \code{fix} is not specified
    (the default), starting values are varied randomly across chains.}

  \item{priors}{a list with each tag corresponding to a parameter name.
    Valid tags are \code{beta.normal}, \code{alpha.normal}, \code{phi.unif},
    \code{sigma.sq.ig}, \code{nu.unif}, \code{sigma.sq.mu.ig}, \code{kappa.unif}, 
    and \code{sigma.sq.p.ig}. Abundance (\code{beta}) and detection (\code{alpha}) 
    regression coefficients are assumed to follow a normal distribution.  
    For \code{beta} hyperparameters of the normal distribution are passed as a list of 
    length two with the first and second elements corresponding to the mean and variance 
    of the normal distribution, which are each specified as vectors of length equal to 
    the number of coefficients to be estimated or of length one if priors are the same 
    for all coefficients. For the detection coefficients \code{alpha}, the mean
    and variance hyperparameters are themselves passed in as lists, with each element of 
    the list corresponding to the specific hyperparameters for the detection parameters 
    in a given data source. If not specified, prior means are set to 0 and prior 
    variances set to 2.72. The spatial variance parameter, \code{sigma.sq} is assumed to follow an
    inverse-Gamma distribution. The spatial decay \code{phi} and spatial
    smoothness \code{nu}, parameters are assumed to follow Uniform
    distributions. The hyperparameters of the inverse-Gamma for \code{sigma.sq}
    is passed as a list of length two with the first and second elements corresponding
    to the shape and scale parameters of the inverse-Gamma distribution either for
    each spatially-varying coefficient, or a single value if assuming the same values
    for all spatially-varying coefficients. 
    The hyperparameters of the Uniform are passed as a list of
    length two with the first and second elements corresponding to
    the lower and upper support, respectively, for each SVC or a single value
    if giving the same prior for each SVC. \code{sigma.sq.mu} and \code{sigma.sq.p}
    are the random effect variances for any abundance and detection random effects, 
    respectively, and are assumed to follow an inverse-Gamma distribution. The hyperparameters of
    the inverse-Gamma distribution are passed as a list of length two with the
    first and second elements corresponding to the shape and scale parameters,
    respectively, which are each specified as vectors of length equal to the
    number of random effects or of length one if priors are the same for all
    random effect variances. The negative binomial dispersion parameter \code{kappa} is 
    assumed to follow a Uniform distribution. The hyperparameters of the Uniform are passed 
    as a list of two elements that correspond to the lower and upper bound of the uniform 
    prior distribution. Each of the two elements should be a vector, with the number of values in 
    each vector equal to one or the number of data sources in the model.}

  \item{svc.cols}{a vector indicating the variables whose effects will be
    estimated as spatially-varying coefficients. \code{svc.cols} can be an
    integer vector with values indicating the order of covariates specified
    in the model formula (with 1 being the intercept if specified), or it can
    be specified as a character vector with names corresponding to variable
    names in \code{abund.covs} (for the intercept, use '(Intercept)').}

  \item{cov.model}{a quoted keyword that specifies the covariance
    function used to model the spatial dependence structure among the
    observations.  Supported covariance model key words are:
    \code{"exponential"}, \code{"matern"}, \code{"spherical"}, and
    \code{"gaussian"}.}

  \item{tuning}{a list with each tag corresponding to a parameter name,
    whose value defines the initial variance of the adaptive sampler.
    Valid tags are \code{phi}, \code{nu}, \code{beta}, \code{alpha}, \code{beta.star}
    (the abundance random effect values), \code{alpha.star} (the detection 
    random effect values), \code{kappa}, \code{w}. See Roberts and Rosenthal (2009) for details.}

  \item{NNGP}{if \code{TRUE}, model is fit with an NNGP. See Datta et al. (2016) and
    Finley et al. (2019) for more information. Currently only NNGP is supported.}

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

  \item{family}{the distribution to use for the latent abundance process. Currently
    supports \code{'NB'} (negative binomial) and \code{'Poisson'}. This can be a single
    value if using the same distribution for all data sets, or it can be a vector with a 
    specific distribution for each data set being modeled.}

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
  Jeffrey W. Doser \email{doserjef@msu.edu}
}

\value{
  An object of class \code{svcTIntAbund} that is a list comprised of:

  \item{beta.samples}{a \code{coda} object of posterior samples
    for the abundance regression coefficients.}

  \item{alpha.samples}{a \code{coda} object of posterior samples
    for the detection regression coefficients for all data sources.}

  \item{kappa.samples}{a \code{coda} object of posterior samples
    for the abundance dispersion parameters. Only included when
    \code{family = 'NB'}.}

  \item{y.rep.samples}{a two or three-dimensional object of posterior samples
    for the abundance replicate (fitted) values with dimensions
    corresponding to MCMC samples, site, and replicate.}

  \item{mu.samples}{a two or -three-dimensional array of posterior samples
    for the expected abundance samples with dimensions corresponding
    to MCMC samples, site, and season (primary time period). Note that this 
    corresponds to only the parameters related to the abundance portion of the model, 
    not the covariates that are designed as observational (detection) covariates.}
  
  \item{lambda.samples}{a two or three-dimensional object of posterior samples
    for the expected relative abundance (including both the abundance and observational 
    variables) with dimensions corresponding to MCMC samples, site, and replicate.}

  \item{theta.samples}{a \code{coda} object of posterior samples
    for spatial covariance parameters.}

  \item{w.samples}{a three-dimensional array of posterior samples
    for the spatially-varying coefficients with dimensions corresponding
    to MCMC sample, SVC, and site.}

  \item{sigma.sq.mu.samples}{a \code{coda} object of posterior samples
    for variances of abundance random effects included in the model.
    Only included if random effects are specified in
    \code{abund.formula}.}

  \item{sigma.sq.p.samples}{a \code{coda} object of posterior samples
    for variances of random intercpets included in the detection portion 
    of the model. Includes random effect variances for all data sources. 
    Only included if random intercepts are specified in 
    \code{det.formula}.}

  \item{beta.star.samples}{a \code{coda} object of posterior samples
    for the abundance random effects. Only included if random effects
    are specified in \code{formula}.}

  \item{alpha.star.samples}{a \code{coda} object of posterior samples
    for the detection random effects in any of the data sources. 
    Only included if random intercepts are specified in at least one of the 
    individual data set detection formulas in \code{det.formula}.}

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
# TODO: 
}
