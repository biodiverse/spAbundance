\name{ppcAbund}
\alias{ppcAbund}
\title{Function for performing posterior predictive checks}

\usage{
ppcAbund(object, fit.stat, group, type = 'marginal', ...)
}

\description{
  Function for performing posterior predictive checks on \code{spAbundance} model objects.
}

\arguments{
\item{object}{an object of class \code{NMix}, \code{spNMix}, 
	\code{msNMix}, \code{lfMsNMix}, \code{sfMsNMix},
	\code{abund}, \code{spAbund}, \code{msAbund}, 
	\code{lfMsAbund}, \code{sfMsAbund}, \code{DS}, \code{spDS}, 
	\code{msDS}, \code{lfMsDS}, \code{sfMsDS}.}

\item{fit.stat}{a quoted keyword that specifies the fit statistic
to use in the posterior predictive check. Supported fit statistics are
\code{"freeman-tukey"} and \code{"chi-squared"}.}

\item{group}{a positive integer indicating the way to group the
abundance data for the posterior predictive check. Value 
0 will not group the data and use the raw counts, 
1 will group values by row (site), and value 2 will group values by 
column (replicate).} 

\item{type}{a character string indicating whether fitted values should be generated conditional on the estimated latent abundance values (\code{type = 'conditional'}) estimated during the model or based on the marginal expected abundance values (\code{type = 'marginal'}). This is only relevant for N-mixture models.}

\item{...}{currently no additional arguments}
}

\author{
  Jeffrey W. Doser \email{doserjef@msu.edu}, \cr
}

\note{
  \code{ppcAbund} will return an error for Gaussian or zero-inflated Gaussian models. 
  For Gaussian models, standard residual diagnostics can be used to assess model fit.
}

\value{

  An object of class \code{ppcAbund} that is a list comprised of: 

  \item{fit.y}{a numeric vector of posterior samples for the 
  fit statistic calculated on the observed data when \code{object} is 
  of class \code{NMix}, \code{spNMix}, \code{abund}, \code{spAbund}, \code{DS}, or \code{spDS}. 
  When \code{object} is of class \code{msNMix}, \code{lfMsNMix}, \code{sfMsNMix},
  \code{msAbund}, \code{lfMsAbund}, \code{sfMsAbund}, 
  \code{msDS}, \code{lfMsDS}, \code{sfMsDS}, 
  this is a numeric matrix with rows corresponding to posterior samples
  and columns corresponding to species.}

  \item{fit.y.rep}{a numeric vector of posterior samples for the 
  fit statistic calculated on a replicate data set generated from the 
  model when \code{object} is of class \code{NMix}, \code{spNMix}, \code{abund}, 
  \code{spAbund}, \code{DS}, \code{spDS}. When \code{object} is of 
  class \code{msNMix}, \code{lfMsNMix}, \code{sfMsNMix}, \code{msAbund}, 
  \code{lfMsAbund}, \code{sfMsAbund}, \code{msDS}, \code{lfMsDS}, \code{sfMsDS},
  this is a numeric matrix with rows corresponding to posterior samples
  and columns corresponding to species.}

  \item{fit.y.group.quants}{a matrix consisting of posterior quantiles
  for the fit statistic using the observed data for each unique element
  the fit statistic is calculated for (i.e., observations when group = 0, sites when group = 1,
  replicates when group = 2) when \code{object} is of class \code{NMix}, \code{spNMix}, 
  \code{abund}, \code{spAbund}, \code{DS}, or \code{spDS}. 
  When \code{object} is of class \code{msNMix}, \code{lfMsNMix}, \code{sfMsNMix}, 
  \code{msAbund}, \code{lfMsAbund}, \code{sfMsAbund}, \code{msDS}, 
  \code{lfMsDS}, \code{sfMsDS}, this is a three-dimensional
  array with the additional dimension corresponding to species.}

  \item{fit.y.rep.group.quants}{a matrix consisting of posterior quantiles
  for the fit statistic using the model replicated data for each unique element
  the fit statistic is calculated for (i.e., observations when group = 0, sites when group = 1, 
  replicates when group = 2) when \code{object} is of class \code{NMix}, \code{spNMix}, 
  \code{abund}, \code{spAbund}, \code{DS}, \code{spDS}. 
  When \code{object} is of class \code{msNMix}, \code{sfMsNMix}, \code{msAbund}, 
  \code{lfMsAbund}, \code{sfMsAbund}, \code{msDS}, \code{lfMsDS}, 
  \code{sfMsDS}, this is a three-dimensional
  array with the additional dimension corresponding to species.}

  The return object will include additional objects used for 
  standard extractor functions. 

}

\examples{
set.seed(1010)
J.x <- 10
J.y <- 10
J <- J.x * J.y
n.rep <- sample(3, J, replace = TRUE)
beta <- c(0, -1.5)
p.abund <- length(beta)
alpha <- c(0.5, 1.2, -0.5)
p.det <- length(alpha)
mu.RE <- list()
p.RE <- list()
phi <- 3/.6
sigma.sq <- 2
kappa <- 0.3
sp <- FALSE 
cov.model <- 'exponential'
dist <- 'NB'
dat <- simNMix(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, alpha = alpha,
               kappa = kappa, mu.RE = mu.RE, p.RE = p.RE, sp = sp, 
               phi = phi, sigma.sq = sigma.sq, cov.model = cov.model, 
               family = 'NB')

y <- dat$y
X <- dat$X
X.p <- dat$X.p

abund.covs <- X
colnames(abund.covs) <- c('int', 'abund.cov.1')

det.covs <- list(det.cov.1 = X.p[, , 2], det.cov.2 = X.p[, , 3])

data.list <- list(y = y, abund.covs = abund.covs,
                  det.covs = det.covs)

# Priors
prior.list <- list(beta.normal = list(mean = rep(0, p.abund), 
                                      var = rep(100, p.abund)),
                   alpha.normal = list(mean = rep(0, p.det),
                                       var = rep(2.72, p.det)), 
                   kappa.unif = c(0, 10)) 
# Starting values
inits.list <- list(alpha = 0, beta = 0, kappa = kappa, 
                   N = apply(y, 1, max, na.rm = TRUE))


tuning <- 0.5
n.batch <- 4
batch.length <- 25
n.burn <- 50
n.thin <- 1
n.chains <- 1

out <- NMix(abund.formula = ~ abund.cov.1,
            det.formula = ~ det.cov.1 + det.cov.2, 
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

# Posterior predictive check
ppc.out <- ppcAbund(out, fit.stat = 'chi-squared', group = 0)
summary(ppc.out)
}
