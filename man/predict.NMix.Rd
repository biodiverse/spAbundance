\name{predict.NMix}
\alias{predict.NMix}
\title{Function for prediction at new locations for single-species N-mixture models}

\description{
  The function \code{predict} collects posterior predictive samples for a set of new locations given an object of class `NMix`. Prediction is possible for both the latent abundance state as well as detection. 
}

\usage{
\method{predict}{NMix}(object, X.0, ignore.RE = FALSE, type = 'abundance', ...)
}

\arguments{

  \item{object}{an object of class NMix}

  \item{X.0}{the design matrix of covariates at the prediction locations. This should include a column of 1s for the intercept if an intercept is included in the model. If random effects are included in the abundance (or detection if \code{type = 'detection'}) portion of the model, the levels of the random effects at the new locations should be included as a column in the design matrix. The ordering of the levels should match the ordering used to fit the data in \code{NMix}. Columns should correspond to the order of how covariates were specified in the corresponding formula argument of \code{NMix}. Column names of all variables must match the names of variables used when fitting the model (for the intercept, use \code{'(Intercept)'}).}

  \item{ignore.RE}{logical value that specifies whether or not to remove random abundance (or detection if \code{type = 'detection'}) effects from the subsequent predictions. If \code{TRUE}, random effects will be included. If \code{FALSE}, random effects will be set to 0 and predictions will only be generated from the fixed effects.}

  \item{type}{a quoted keyword indicating what type of prediction to produce. Valid keywords are 'abundance' to predict latent abundance and expected abundance values (this is the default), or 'detection' to predict detection probability given new values of detection covariates.}

  \item{...}{currently no additional arguments}
}

\note{
  When \code{ignore.RE = FALSE}, both sampled levels and non-sampled levels of random effects are supported for prediction. For sampled levels, the posterior distribution for the random intercept corresponding to that level of the random effect will be used in the prediction. For non-sampled levels, random values are drawn from a normal distribution using the posterior samples of the random effect variance, which results in fully propagated uncertainty in predictions with models that incorporate random effects. 
}

\author{
  Jeffrey W. Doser \email{doserjef@msu.edu}, \cr
}

\value{
  A list object of class \code{predict.NMix}. When \code{type = 'abundance'}, the list consists of: 

  \item{mu.0.samples}{a \code{coda} object of posterior predictive samples for the
    expected abundance values. Note these will be per unit area if an offset
    was used when fitting the model with \code{NMix()}}

  \item{N.0.samples}{a \code{coda} object of posterior predictive samples for the
    latent abundance values. These will be in the same units as \code{mu.0.samples}.}

  When \code{type = 'detection'}, the list consists of: 

  \item{p.0.samples}{a \code{coda} object of posterior predictive samples for the 
    detection probability values.}

  The return object will include additional objects used for standard 
  extractor functions. 
}

\examples{
set.seed(100)
# Simulate Data -----------------------------------------------------------
J.x <- 10
J.y <- 10
J <- J.x * J.y
n.rep <- sample(2:4, J, replace = TRUE)
beta <- c(0.5, 2)
p.abund <- length(beta)
alpha <- c(0, 1)
p.det <- length(alpha)
dat <- simNMix(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, alpha = alpha,
               sp = FALSE)
# Split into fitting and prediction data set
pred.indx <- sample(1:J, round(J * .25), replace = FALSE)
y <- dat$y[-pred.indx, ]
# Abundance covariates
X <- dat$X[-pred.indx, ]
# Prediction covariates
X.0 <- dat$X[pred.indx, ]
# Detection covariates
X.p <- dat$X.p[-pred.indx, , ]

# Package all data into a list
abund.covs <- X[, 2, drop = FALSE]
colnames(abund.covs) <- c('abund.cov')
det.covs <- list(det.cov = X.p[, , 2])
data.list <- list(y = y, 
                  abund.covs = abund.covs,
                  det.covs = det.covs)
# Priors
prior.list <- list(beta.normal = list(mean = rep(0, p.abund),
                                      var = rep(100, p.abund)),
                   alpha.normal = list(mean = rep(0, p.det),
                                       var = rep(2.72, p.det)), 
                   kappa.unif = c(0.001, 10))
# Initial values
inits.list <- list(alpha = rep(0, p.det),
                   beta = rep(0, p.abund),
                   kappa = 0.5,
                   N = apply(y, 1, max, na.rm = TRUE))

n.batch <- 10
batch.length <- 25
n.burn <- 0
n.thin <- 1
n.chains <- 1

out <- NMix(abund.formula = ~ abund.cov, 
            det.formula = ~ det.cov,
            data = data.list, 
            inits = inits.list,
            n.batch = n.batch, 
            batch.length = batch.length,
            family = 'Poisson',
            priors = prior.list,
            n.omp.threads = 1,
            verbose = TRUE,
            n.report = 1, 
            n.burn = n.burn, 
            n.thin = n.thin, 
            n.chains = n.chains)

summary(out)

# Predict at new locations ------------------------------------------------
colnames(X.0) <- c('intercept', 'abund.cov')
out.pred <- predict(out, X.0)
mu.0.quants <- apply(out.pred$mu.0.samples, 2, quantile, c(0.025, 0.5, 0.975))
plot(dat$mu[pred.indx], mu.0.quants[2, ], pch = 19, xlab = 'True', 
     ylab = 'Fitted', ylim = c(min(mu.0.quants), max(mu.0.quants)))
segments(dat$mu[pred.indx], mu.0.quants[1, ], dat$mu[pred.indx], mu.0.quants[3, ])
lines(dat$mu[pred.indx], dat$mu[pred.indx])
}
