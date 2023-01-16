\name{simMsAbund}
\alias{simMsAbund}
\title{Simulate Multi-Species Count Data without Imperfect Detection}

\description{
  The function \code{simMsAbund} simulates multi-species count data without imperfect detection for simulation studies, power assessments, or function testing. Data can be optionally simulated with a spatial Gaussian Process in the model, as well as an option to allow for species correlations using a factor modeling approach. Non-spatial random intercepts can also be included in the abundance portions of the model.
}

\usage{
simMsAbund(J.x, J.y, n.rep, n.sp, beta, kappa, mu.RE = list(), 
          sp = FALSE, cov.model, 
          sigma.sq, phi, nu, family = 'Poisson',
          factor.model = FALSE, n.factors, ...)
}

\arguments{

\item{J.x}{a single numeric value indicating the number of sites to simulate count data along the horizontal axis. Total number of sites with simulated data is \eqn{J.x \times J.y}{J.x x J.y}. }

\item{J.y}{a single numeric value indicating the number of sites to simulate count data along the vertical axis. Total number of sites with simulated data is \eqn{J.x \times J.y}{J.x x J.y}. }

\item{n.rep}{a numeric vector of length \eqn{J = J.x \times J.y}{J = J.x x J.y} indicating the number of replicate surveys at each of the \eqn{J}{J} sites.}

\item{n.sp}{a single numeric value indicating the number of species to simulate count data.}

\item{beta}{a numeric matrix with \code{n.sp} rows containing the intercept and regression coefficient parameters for the model. Each row corresponds to the regression coefficients for a given species.} 

\item{kappa}{a numeric vector of length \code{n.sp} containing the dispersion parameter for the model for each species. Only relevant when \code{family = 'NB'}.}

\item{mu.RE}{a list used to specify the non-spatial random intercepts included in the model. The list must have two tags: \code{levels} and \code{sigma.sq.mu}. \code{levels} is a vector of length equal to the number of distinct random intercepts to include in the model and contains the number of levels there are in each intercept. \code{sigma.sq.mu} is a vector of length equal to the number of distinct random intercepts to include in the model and contains the variances for each random effect. If not specified, no random effects are included in the model.}

\item{sp}{a logical value indicating whether to simulate a spatially-explicit model with a Gaussian process. By default set to \code{FALSE}.}

\item{cov.model}{a quoted keyword that specifies the covariance function used to model the spatial dependence structure among the abundance values.  Supported covariance model key words are: \code{"exponential"}, \code{"matern"}, \code{"spherical"}, and \code{"gaussian"}.}

\item{sigma.sq}{a numeric vector of length \code{n.sp} containing the spatial variance parameter for each species. Ignored when \code{sp = FALSE} or when \code{factor.model = TRUE}.}

\item{phi}{a numeric vector of length \code{n.sp} containing the spatial decay parameter for each species. Ignored when \code{sp = FALSE}. If \code{factor.model = TRUE}, this should be of length \code{n.factors}.}

\item{nu}{a numeric vector of length \code{n.sp} containing the spatial smoothness parameter for each species. Only used when \code{sp = TRUE} and \code{cov.model = 'matern'}. If \code{factor.model = TRUE}, this should be of length \code{n.factors}.} 

\item{factor.model}{a logical value indicating whether to simulate data following a factor modeling approach that explicitly incoporates species correlations. If \code{sp = TRUE}, the latent factors are simulated from independent spatial processes. If \code{sp = FALSE}, the latent factors are simulated from standard normal distributions.}

\item{n.factors}{a single numeric value specifying the number of latent factors to use to simulate the data if \code{factor.model = TRUE}.}

\item{family}{the distribution to use for the latent abundance process. Currently
  supports \code{'NB'} (negative binomial) and \code{'Poisson'}.}

\item{...}{currently no additional arguments}
}

\author{
  Jeffrey W. Doser \email{doserjef@msu.edu}, \cr
  Andrew O. Finley \email{finleya@msu.edu}
}

\value{
  A list comprised of: 

  \item{X}{a three-dimensional numeric design array of covariates with dimensions corresponding to sites, replicates, and number of covariates (including an intercept) for the model.}
  \item{coords}{a \eqn{J \times 2}{J x 2} numeric matrix of coordinates of each site. Required for spatial models.}
  \item{w}{a \eqn{N \times J}{N x J} matrix of the spatial random effects for each species. Only used to simulate data when \code{sp = TRUE}. If \code{factor.model = TRUE}, the first dimension is \code{n.factors}.}
  \item{mu}{a \code{n.sp x J} matrix of the mean abundances for each species at each site.}
  \item{y}{a \code{n.sp x J x max(n.rep)} array of the raw count data for each species at each site and replicate combination. Sites with fewer than \code{max(n.rep)} replicates will contain \code{NA} values.}
  \item{X.re}{a numeric matrix containing the levels of any abundance random effect included in the model. Only relevant when abundance random effects are specified in \code{mu.RE}.}
  \item{beta.star}{a numeric matrix where each row contains the simulated abundance random effects for each given level of the random effects included in the abundance model. Only relevant when abundance random effects are included in the model.}
}

\examples{
set.seed(408)
J.x <- 8
J.y <- 8
J <- J.x * J.y
n.rep <- sample(3, size = J, replace = TRUE)
n.sp <- 6
# Community-level covariate effects
beta.mean <- c(-2, 0.5)
p.abund <- length(beta.mean)
tau.sq.beta <- c(0.2, 1.2)
# Random effects (two random intercepts)
mu.RE <- list(levels = c(10, 15),
	       sigma.sq.mu = c(0.43, 0.5))
# Draw species-level effects from community means.
beta <- matrix(NA, nrow = n.sp, ncol = p.abund)
for (i in 1:p.abund) {
  beta[, i] <- rnorm(n.sp, beta.mean[i], sqrt(tau.sq.beta[i]))
}
sp <- TRUE
n.factors <- 2
factor.model <- TRUE
phi <- runif(n.factors, 3/1, 3 / .1)
kappa <- runif(n.sp, 0.1, 1)

dat <- simMsAbund(J.x = J.x, J.y = J.y, n.rep = n.rep, n.sp = n.sp, beta = beta,
	          mu.RE = mu.RE, sp = sp, kappa = kappa, family = 'NB',
                  factor.model = factor.model, phi = phi,
                  cov.model = 'spherical', n.factors = n.factors)
}
