\name{waicAbund}
\alias{waicAbund}
\title{Compute Widely Applicable Information Criterion for spAbundance Model Objects}

\usage{
waicAbund(object, N.max, by.species = FALSE, ...)
}

\description{
  Function for computing the Widely Applicable Information Criterion 
  (WAIC; Watanabe 2010) for \code{spAbundance} model objects. 
}

\arguments{
  \item{object}{an object of class \code{NMix}, \code{spNMix}, \code{msNMix}, 
  \code{lfMsNMix}, \code{sfMsNMix}, \code{abund}, \code{spAbund}, 
  \code{msAbund}, \code{lfMsAbund}, \code{sfMsAbund}, \code{DS}, \code{spDS}, 
  \code{msDS}, \code{lfMsDS}, \code{sfMsDS}.}

  \item{N.max}{values indicating the upper limit on the latent abundance 
  values when calculating WAIC for N-mixture models or hierarchical distance sampling models. 
  For single-species models, this can be a single value or a vector of 
  different values for each site. For multi-species models, this can be a 
  single value, a vector of values for each species, or a species by site matrix
  for a separate value for each species/site combination. 
  Defaults to ten plus the largest abundance value for each site/species 
  in the posterior samples \code{object$N.samples}.}

  \item{by.species}{a logical value indicating whether or not WAIC should be
    reported individually for each species (\code{TRUE}) or summed across
    the entire community (\code{FALSE}) for multi-species models. Ignored for 
    single species models.}

  \item{...}{currently no additional arguments}

}

\note{
  When fitting zero-inflated Gaussian models, the WAIC is only calculated for the 
  non-zero values. If fitting a first stage model with 
  \code{spOccupancy} to the binary portion of the 
  zero-inflated model, you can use the \code{spOccupancy::waicOcc} function to 
  calculate WAIC for the binary component. 
}

\references{
  Watanabe, S. (2010). Asymptotic equivalence of Bayes cross
  validation and widely applicable information criterion in 
  singular learning theory. \emph{Journal of Machine Learning Research}, 
  11:3571-3594.

  Gelman, A., J. B. Carlin, H. S. Stern, D. B. Dunson, A. Vehtari, 
  and D. B. Rubin. (2013). Bayesian Data Analysis. 3rd edition. 
  CRC Press, Taylor and Francis Group

  Gelman, A., J. Hwang, and A. Vehtari (2014). Understanding predictive
  information criteria for Bayesian models. \emph{Statistics and Computing}, 
  24:997-1016.
}

\author{
  Jeffrey W. Doser \email{doserjef@msu.edu}, \cr
}

\value{
  Returns a vector with three elements corresponding to 
  estimates of the expected log pointwise predictive density (elpd), the 
  effective number of parameters (pD), and the WAIC. If calculating WAIC
  for a multi-species model and \code{by.species = TRUE}, this will be a 
  data frame with rows corresponding to the different species.
}

\details{
  The effective number of parameters is calculated following the recommendations
  of Gelman et al. (2014).
}

\examples{
set.seed(1010)
J.x <- 15
J.y <- 15
J <- J.x * J.y
n.rep <- sample(3, J, replace = TRUE)
beta <- c(0, -1.5, 0.3, -0.8)
p.abund <- length(beta)
mu.RE <- list(levels = c(30), sigma.sq.mu = c(1.3))
kappa <- 0.5
sp <- FALSE 
family <- 'NB'
dat <- simAbund(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, 
                kappa = kappa, mu.RE = mu.RE, sp = sp, family = 'NB')

y <- dat$y
X <- dat$X
X.re <- dat$X.re

abund.covs <- list(int = X[, , 1], 
                   abund.cov.1 = X[, , 2], 
                   abund.cov.2 = X[, , 3], 
                   abund.cov.3 = X[, , 4],
                   abund.factor.1 = X.re[, , 1])

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

out <- abund(formula = ~ abund.cov.1 + abund.cov.2 + abund.cov.3 + 
                         (1 | abund.factor.1),
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

# Calculate WAIC
waicAbund(out)
}

