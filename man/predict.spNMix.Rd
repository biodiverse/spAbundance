\name{predict.spNMix}
\alias{predict.spNMix}
\title{Function for prediction at new locations for single-species spatial N-mixture models}

\description{
  The function \code{predict} collects posterior predictive samples for a set of new locations given an object of class `spNMix`. Prediction is possible for both the latent abundance state as well as detection. 
}

\usage{
\method{predict}{spNMix}(object, X.0, coords.0, n.omp.threads = 1, 
        verbose = TRUE, n.report = 100, ignore.RE = FALSE, 
        type = 'abundance', include.sp = TRUE, ...)
}

\arguments{

  \item{object}{an object of class spNMix}

  \item{X.0}{the design matrix of covariates at the prediction locations. This should include a column of 1s for the intercept if an intercept is included in the model. If random effects are included in the abundance (or detection if \code{type = 'detection'}) portion of the model, the levels of the random effects at the new locations should be included as a column in the design matrix. The ordering of the levels should match the ordering used to fit the data in \code{spNMix}. Columns should correspond to the order of how covariates were specified in the corresponding formula argument of \code{spNMix}. Column names of all variables must match the names of variables used when fitting the model (for the intercept, use \code{'(Intercept)'}).}

  \item{coords.0}{the spatial coordinates corresponding to \code{X.0}. Note that \code{spAbundance} assumes coordinates are specified 
  in a projected coordinate system.}

   \item{n.omp.threads}{a positive integer indicating
    the number of threads to use for SMP parallel processing. The package must
    be compiled for OpenMP support. For most Intel-based machines, we recommend setting
    \code{n.omp.threads} up to the number of hyperthreaded cores. 
    Note, \code{n.omp.threads} > 1 might not work on some systems.}

  \item{verbose}{if \code{TRUE}, model specification and progress of the 
    sampler is printed to the screen. Otherwise, nothing is printed to
    the screen.}
  
  \item{n.report}{the interval to report sampling progress.}

  \item{ignore.RE}{logical value that specifies whether or not to remove random abundance (or detection if \code{type = 'detection'}) effects from the subsequent predictions. If \code{TRUE}, random effects will be included. If \code{FALSE}, random effects will be set to 0 and predictions will only be generated from the fixed effects.}

  \item{type}{a quoted keyword indicating what type of prediction to produce. Valid keywords are 'abundance' to predict latent abundance and expected abundance values (this is the default), or 'detection' to predict detection probability given new values of detection covariates.}

  \item{include.sp}{a logical value used to indicate whether spatial random effects should be included in the predictions. By default, this is set to \code{TRUE}. If set to \code{FALSE}, predictions are given using the covariates and any unstructured random effects in the model. If \code{FALSE}, the \code{coords} argument is not required.}

  \item{...}{currently no additional arguments}
}

\note{
  When \code{ignore.RE = FALSE}, both sampled levels and non-sampled levels of random effects are supported for prediction. For sampled levels, the posterior distribution for the random intercept corresponding to that level of the random effect will be used in the prediction. For non-sampled levels, random values are drawn from a normal distribution using the posterior samples of the random effect variance, which results in fully propagated uncertainty in predictions with models that incorporate random effects. 
}

\author{
  Jeffrey W. Doser \email{doserjef@msu.edu}, \cr
  Andrew O. Finley \email{finleya@msu.edu}, \cr
}

\value{
  A list object of class \code{predict.spNMix}. When \code{type = 'abundance'}, the list consists of: 

  \item{mu.0.samples}{a \code{coda} object of posterior predictive samples for the
    expected abundance values. Note these will be per unit area if an offset
    was used when fitting the model with \code{NMix()}}

  \item{N.0.samples}{a \code{coda} object of posterior predictive samples for the
    latent abundance values. These will be in the same units as \code{mu.0.samples}.}

  \item{w.0.samples}{a \code{coda} object of posterior predictive samples 
    for the latent spatial random effects.}

  When \code{type = 'detection'}, the list consists of: 

  \item{p.0.samples}{a \code{coda} object of posterior predictive samples for the 
    detection probability values.}

  The return object will include additional objects used for standard 
  extractor functions. 
}

\examples{
set.seed(200)
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

# Split into fitting and prediction data set
pred.indx <- sample(1:J, round(J * .5), replace = FALSE)
y <- dat$y[-pred.indx, ]
# Abundance covariates
X <- dat$X[-pred.indx, ]
# Prediction covariates
X.0 <- dat$X[pred.indx, ]
# Detection covariates
X.p <- dat$X.p[-pred.indx, , ]
coords <- as.matrix(dat$coords[-pred.indx, ])
coords.0 <- as.matrix(dat$coords[pred.indx, ])
mu.0 <- dat$mu[pred.indx]
w.0 <- dat$w[pred.indx]

abund.covs <- X
colnames(abund.covs) <- c('int', 'abund.cov.1')

det.covs <- list(det.cov.1 = X.p[, , 2], det.cov.2 = X.p[, , 3])

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

# Predict at new locations ------------------------------------------------
colnames(X.0) <- c('intercept', 'abund.cov')
out.pred <- predict(out, X.0, coords.0)
str(out.pred)
}
