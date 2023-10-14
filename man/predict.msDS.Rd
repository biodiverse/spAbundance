\name{predict.msDS}
\alias{predict.msDS}
\title{Function for prediction at new locations for multi-species hierarchical distance sampling models}

\description{
  The function \code{predict} collects posterior predictive samples for a set of new locations given an object of class `msDS`. Prediction is possible for both the latent abundance state as well as detection.
}

\usage{
\method{predict}{msDS}(object, X.0, ignore.RE = FALSE, type = 'abundance', ...)
}

\arguments{

  \item{object}{an object of class msDS}

  \item{X.0}{the design matrix of covariates at the prediction locations. This should include a column of 1s for the intercept if an intercept is included in the model. If random effects are included in the abundance (or detection if \code{type = 'detection'}) portion of the model, the levels of the random effects at the new locations should be included as a column in the design matrix. The ordering of the levels should match the ordering used to fit the data in \code{msDS}. Columns should correspond to the order of how covariates were specified in the corresponding formula argument of \code{msDS}. Column names must match the names of the variables used to fit the model (for the intercept, use \code{'(Intercept)'}).}

  \item{ignore.RE}{a logical value indicating whether to include unstructured random effects for prediction. If TRUE, random effects will be ignored and prediction will only use the fixed effects. If FALSE, random effects will be included in the prediction for both observed and unobserved levels of the random effect.}

  \item{type}{a quoted keyword indicating what type of prediction to produce. Valid keywords are 'abundance' to predict expected abundance and latent abundance values (this is the default), or 'detection' to predict detection probability given new values of detection covariates.}
  
  \item{...}{currently no additional arguments}

}

\note{
  When \code{ignore.RE = FALSE}, both sampled levels and non-sampled levels of random effects are supported for prediction. For sampled levels, the posterior distribution for the random effect corresponding to that level of the random effect will be used in the prediction. For non-sampled levels, random values are drawn from a normal distribution using the posterior samples of the random effect variance, which results in fully propagated uncertainty in predictions with models that incorporate random effects.
}

\author{
  Jeffrey W. Doser \email{doserjef@msu.edu}, \cr
}

\value{
  A list object of class \code{predict.msDS}. When \code{type = 'abundance'}, the list consists of:  

  \item{mu.0.samples}{a three-dimensional array of posterior predictive samples for the
    expected abundance values, or expected abundance values per unit area (i.e., density)
    values when an offset was used when fitting the model with \code{msDS()}.}

  \item{N.0.samples}{a three-dimensional array of posterior predictive samples for the
    latent abundance values. These will be in the same units as \code{mu.0.samples}.}

  When \code{type = 'detection'}, the list consists of: 

  \item{sigma.0.samples}{a three-dimensional array of posterior predictive samples for
    sigma (the parameter controlling detection probability).}

  The return object will include additional objects used for standard 
  extractor functions. 
}

\examples{
set.seed(210)
J.x <- 10
J.y <- 10 
J <- J.x * J.y
# Number of distance bins from which to simulate data. 
n.bins <- 5
# Length of each bin. This should be of length n.bins
bin.width <- c(.10, .10, .20, .3, .1)
# Number of species
n.sp <- 5
# Community-level abundance coefficients
beta.mean <- c(-1, 0.2, 0.3, -0.2)
p.abund <- length(beta.mean)
tau.sq.beta <- c(0.2, 0.3, 0.5, 0.4)
# Detection coefficients
alpha.mean <- c(-1.0, -0.3)
p.det <- length(alpha.mean)
tau.sq.alpha <- c(0.1, 0.2)
# Detection decay function
det.func <- 'halfnormal'
mu.RE <- list()
p.RE <- list()
# Draw species-level effects from community means.
beta <- matrix(NA, nrow = n.sp, ncol = p.abund)
alpha <- matrix(NA, nrow = n.sp, ncol = p.det)
for (i in 1:p.abund) {
  beta[, i] <- rnorm(n.sp, beta.mean[i], sqrt(tau.sq.beta[i]))
}
for (i in 1:p.det) {
  alpha[, i] <- rnorm(n.sp, alpha.mean[i], sqrt(tau.sq.alpha[i]))
}
sp <- FALSE 
family <- 'Poisson'
kappa <- runif(n.sp, 0.3, 3) 
offset <- pi * .8^2
transect <- 'line'
factor.model <- FALSE

dat <- simMsDS(J.x = J.x, J.y = J.y, n.bins = n.bins, bin.width = bin.width,
               n.sp = n.sp, beta = beta, alpha = alpha, det.func = det.func, kappa = kappa, 
               mu.RE = mu.RE, p.RE = p.RE, sp = sp, cov.model = cov.model,
               sigma.sq = sigma.sq, phi = phi, nu = nu, family = family, 
               offset = offset, transect = transect, factor.model = factor.model)

# Split into fitting and prediction data set
pred.indx <- sample(1:J, round(J * .25), replace = FALSE)
y <- dat$y[, -pred.indx, ]
# Occupancy covariates
X <- dat$X[-pred.indx, ]
# Prediction covariates
X.0 <- dat$X[pred.indx, ]
# Detection covariates
X.p <- dat$X.p[-pred.indx, , drop = FALSE]
X.p.0 <- dat$X.p[pred.indx, , drop = FALSE]
coords <- as.matrix(dat$coords[-pred.indx, ])
coords.0 <- as.matrix(dat$coords[pred.indx, ])
dist.breaks <- dat$dist.breaks

covs <- cbind(X, X.p)
colnames(covs) <- c('int.abund', 'abund.cov.1', 'abund.cov.2', 'abund.cov.3', 
                    'int.det', 'det.cov.1')

data.list <- list(y = y, 
                  covs = covs,
                  dist.breaks = dist.breaks, 
                  offset = offset)

# Priors
prior.list <- list(beta.comm.normal = list(mean = 0, var = 10),
		   alpha.comm.normal = list(mean = 0,
			               var = 10), 
                   kappa.unif = list(0, 100), 
                   tau.sq.beta.ig = list(a = 0.1, b = 0.1),
                   tau.sq.alpha.ig = list(a = 0.1, b = 0.1)) 
# Starting values
inits.list <- list(alpha.comm = 0, beta.comm = 0, beta = 0,
		   alpha = 0, kappa = 1)

tuning <- list(beta = 0.1, alpha = 0.1, beta.star = 0.3, alpha.star = 0.1, 
               kappa = 0.8) 

n.batch <- 4 
batch.length <- 25
n.burn <- 0
n.thin <- 1
n.chains <- 1

out <- msDS(abund.formula = ~ abund.cov.1 + abund.cov.2 + abund.cov.3,
            det.formula = ~ det.cov.1,
            data = data.list, 
            n.batch = n.batch, 
            batch.length = batch.length, 
            inits = inits.list, 
            family = 'Poisson',
            det.func = 'halfnormal', 
            transect = transect, 
            tuning = tuning,
            priors = prior.list, 
            accept.rate = 0.43, 
            n.omp.threads = 1, 
            verbose = TRUE, 
            n.report = 10,
            n.burn = n.burn,
            n.thin = n.thin,
            n.chains = n.chains) 
summary(out, level = 'community')

# Predict at new locations ------------------------------------------------
colnames(X.0) <- c('intercept', 'abund.cov.1', 'abund.cov.2', 'abund.cov.3')
out.pred <- predict(out, X.0)
str(out.pred)
}

