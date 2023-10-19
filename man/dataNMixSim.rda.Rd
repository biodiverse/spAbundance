\name{dataNMixSim}

\docType{data}

\alias{dataNMixSim}

\title{Simulated repeated count data of 6 species 
  across 225 sites}

\description{
A simulated data set of repeated count data for 6 species across 225 sites with a maximum
of 3 replicate surveys performed at a given site.
}

\usage{
data(dataNMixSim)
}

\format{
  \code{dataNMixSim} is a list with four elements: 

   \code{y}: a three-dimensional array of count data with 
     dimensions of species (6), sites (225) and replicates (3). 

   \code{abund.covs}: a numeric matrix with 225 rows and two columns consisting
     of a continuous covariate and a categorical variable which may both influence
     abundance of the different species.

   \code{det.covs}: a list of two numeric matrices with 225 rows and 3 columns. 
     Both matrices contain a continuous covariate that may affect detection probability
     of the species

   \code{coords}: a numeric matrix with 225 rows and two columns containing the 
     site coordinates (X and Y). Note the data are generated across 
     a unit square (i.e., the x and y coordinates are both between 0 and 1). 
}

\keyword{datasets}

\examples{
set.seed(6)
J.x <- 15
J.y <- 15
J <- J.x * J.y
n.rep <- sample(3, size = J, replace = TRUE)
# n.rep <- rep(5, J) 
n.sp <- 6
# Community-level covariate effects
# Occurrence
beta.mean <- c(-1, 0.5)
p.abund <- length(beta.mean)
tau.sq.beta <- c(0.4, 1.2)
# Detection
alpha.mean <- c(0, 0.5, 0.8)
tau.sq.alpha <- c(0.2, 1, 1.5)
p.det <- length(alpha.mean)
# Random effects
mu.RE <- list()
mu.RE <- list(levels = c(10),
	       sigma.sq.mu = c(0.5),
               beta.indx = list(1))
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
alpha.true <- alpha
sp <- TRUE 
n.factors <- 3
factor.model <- TRUE
phi <- runif(n.factors, 3/1, 3 / .2)
kappa <- runif(n.sp, 0.1, 1)
family <- 'Poisson'

dat <- simMsNMix(J.x = J.x, J.y = J.y, n.rep = n.rep, n.sp = n.sp, beta = beta, alpha = alpha,
	        mu.RE = mu.RE, p.RE = p.RE, sp = sp, kappa = kappa, family = family, 
                factor.model = factor.model, phi = phi, 
                cov.model = 'exponential', n.factors = n.factors)
table(dat$N)
apply(dat$N, 1, sum)

y <- dat$y
X <- dat$X
X.p <- dat$X.p
X.re <- dat$X.re
X.p.re <- dat$X.p.re
coords <- dat$coords
dimnames(coords)[[2]] <- c('X', 'Y')

# Package all data into a list
abund.covs <- cbind(X, X.re)
colnames(abund.covs) <- c('int', 'abund.cov.1', 'abund.factor.1')
abund.covs <- abund.covs[, -1]
det.covs <- list(det.cov.1 = X.p[, , 2], 
		 det.cov.2 = X.p[, , 3]) 
dataNMixSim <- list(y = y, 
		    abund.covs = abund.covs, 
		    det.covs = det.covs, 
                    coords = coords)
}
