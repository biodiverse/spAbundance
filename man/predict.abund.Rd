\name{predict.abund}
\alias{predict.abund}
\title{Function for prediction at new locations for univariate GLMMs}

\description{
  The function \code{predict} collects posterior predictive samples for a set of new locations given an object of class `abund`.  
}

\usage{
\method{predict}{abund}(object, X.0, ignore.RE = FALSE, z.0.samples, ...)
}

\arguments{

  \item{object}{an object of class abund}

  \item{X.0}{the design matrix of covariates at the prediction locations. This should be a three-dimensional array, with dimensions corresponding to site, replicate, and covariate, respectively. Note that the first covariate should consist of all 1s for the intercept if an intercept is included in the model. If random effects are included in the abundance portion of the model, the levels of the random effects at the new locations/time periods should be included as an element of the three-dimensional array. The ordering of the levels should match the ordering used to fit the data in \code{abund}. The covariates should be organized in the same order as they were specified in the corresponding formula argument of \code{abund}. Names of the third dimension (covariates) of any random effects in X.0 must match the name of the random effects used to fit the model, if specified in the corresponding formula argument of \code{abund}. See example below. If there is only one replicate per location, the design matrix can be a two-dimensional matrix instead of a three-dimensional array.}

  \item{ignore.RE}{logical value that specifies whether or not to remove unstructured random effects from the subsequent predictions. If \code{TRUE}, unstructured random effects will be included. If \code{FALSE}, unstructured random effects will be set to 0 and predictions will only be generated from the fixed effects.}

  \item{z.0.samples}{a matrix with rows corresponding to MCMC samples and columns corresponding to prediction locations containing the full posterior samples of the predicted binary portion of a zero-inflated Gaussian model. In the context of abundance models, this typically corresponds to estimates of the presence or absence of the species at the location. When using \code{spOccupancy} to generate the first stage samples of the zero-inflated Gaussian model, this is the object contained in the \code{z.0.samples} object of the predition function for the spOccupancy object. Ignored for all model types other than zero-inflated Gaussian.}

  \item{...}{currently no additional arguments}
}

\note{
  When \code{ignore.RE = FALSE}, both sampled levels and non-sampled levels of random effects are supported for prediction. For sampled levels, the posterior distribution for the random effect corresponding to that level of the random effect will be used in the prediction. For non-sampled levels, random values are drawn from a normal distribution using the posterior samples of the random effect variance, which results in fully propagated uncertainty in predictions with models that incorporate random effects. 
}

\author{
  Jeffrey W. Doser \email{doserjef@msu.edu}, \cr
}

\value{
  A list object of class \code{predict.abund}. The list consists of: 

  \item{mu.0.samples}{a three-dimensional object of posterior predictive samples for the
    expected abundance values with dimensions corresponding to posterior predictive
    sample, site, and replicate. When there is no replication, this will be a two-dimensional
    matrix. Note if an offset was used when fitting the model with \code{abund}, the abundance
    values are reported per unit of the offset.}

  \item{y.0.samples}{a three-dimensional object of posterior predictive samples for the
    abundance values with dimensions corresponding to posterior predictive sample, site, 
    and replicate. When there is no replication, this will be a two-dimensional matrix.
    These will be in the same units as \code{mu.0.samples}.}

  The return object will include additional objects used for standard 
  extractor functions. 
}

\examples{
set.seed(1010)
J.x <- 15
J.y <- 15
J <- J.x * J.y
n.rep <- sample(1, J, replace = TRUE)
beta <- c(0, -1.5, 0.3, -0.8)
p.abund <- length(beta)
mu.RE <- list()
kappa <- 0.5
sp <- FALSE 
family <- 'NB'
dat <- simAbund(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, 
                kappa = kappa, mu.RE = mu.RE, sp = sp, family = 'NB')

# Split into fitting and prediction data set
pred.indx <- sample(1:J, round(J * .25), replace = FALSE)
y <- dat$y[-pred.indx, ]
# Abundance covariates
X <- dat$X[-pred.indx, , , drop = FALSE]
# Prediction covariates
X.0 <- dat$X[pred.indx, , ]
coords <- as.matrix(dat$coords[-pred.indx, ])
coords.0 <- as.matrix(dat$coords[pred.indx, ])

abund.covs <- list(int = X[, , 1], 
                   abund.cov.1 = X[, , 2], 
                   abund.cov.2 = X[, , 3], 
                   abund.cov.3 = X[, , 4])

data.list <- list(y = y, covs = abund.covs)

# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 100),
                   kappa.unif = c(0.001, 10)) 
# Starting values
inits.list <- list(beta = 0, kappa = kappa)

n.batch <- 5
batch.length <- 25
n.burn <- 0
n.thin <- 1
n.chains <- 1

out <- abund(formula = ~ abund.cov.1 + abund.cov.2 + abund.cov.3,
             data = data.list, 
             n.batch = n.batch, 
             batch.length = batch.length, 
             inits = inits.list, 
             priors = prior.list, 
             accept.rate = 0.43, 
             n.omp.threads = 1, 
             verbose = TRUE, 
             n.report = 1,
             n.burn = n.burn,
             n.thin = n.thin,
             n.chains = n.chains) 

# Predict at new locations ------------------------------------------------
colnames(X.0) <- c('intercept', 'abund.cov', 'abund.cov.2', 'abund.cov.3')
out.pred <- predict(out, X.0)
mu.0.quants <- apply(out.pred$mu.0.samples, 2, quantile, c(0.025, 0.5, 0.975))
plot(dat$mu[pred.indx], mu.0.quants[2, ], pch = 19, xlab = 'True', 
     ylab = 'Fitted', ylim = c(min(mu.0.quants), max(mu.0.quants)))
segments(dat$mu[pred.indx], mu.0.quants[1, ], dat$mu[pred.indx], mu.0.quants[3, ])
lines(dat$mu[pred.indx], dat$mu[pred.indx])
}
