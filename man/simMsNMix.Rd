\name{simMsNMix}
\alias{simMsNMix}
\title{Simulate Multi-Species Repeated Count Data with Imperfect Detection}

\description{
  The function \code{simMsNMix} simulates multi-species count data for simulation studies, power assessments, or function testing. Data can be optionally simulated with a spatial Gaussian Process in the abundance portion of the model, as well as an option to allow for species correlations using a factor modeling approach. Non-spatial random effects can also be included in the detection or abundance portions of the model.
}

\usage{
simMsNMix(J.x, J.y, n.rep, n.rep.max, n.sp, beta, alpha, kappa, mu.RE = list(), 
          p.RE = list(), offset = 1, sp = FALSE, cov.model, 
          sigma.sq, phi, nu, family = 'Poisson',
          factor.model = FALSE, n.factors, ...)
}

\arguments{

\item{J.x}{a single numeric value indicating the number of sites to simulate count data along the horizontal axis. Total number of sites with simulated data is \eqn{J.x \times J.y}{J.x x J.y}. }

\item{J.y}{a single numeric value indicating the number of sites to simulate count data along the vertical axis. Total number of sites with simulated data is \eqn{J.x \times J.y}{J.x x J.y}. }

\item{n.rep}{a numeric vector of length \eqn{J = J.x \times J.y}{J = J.x x J.y} indicating the number of repeat visits at each of the \eqn{J}{J} sites.}

\item{n.rep.max}{a single numeric value indicating the maximum number of replicate surveys. This is an optional argument, with its default value set to \code{max(n.rep)}. This can be used to generate data sets with different types of missingness (e.g., simulate data across 20 days (replicate surveys) but sites are only sampled a maximum of ten times each).}

\item{n.sp}{a single numeric value indicating the number of species to simulate count data.}

\item{beta}{a numeric matrix with \code{n.sp} rows containing the intercept and regression coefficient parameters for the abundance portion of the multi-species N-mixture model. Each row corresponds to the regression coefficients for a given species.} 

\item{alpha}{a numeric matrix with \code{n.sp} rows containing the intercept and regression coefficient parameters for the detection portion of the multi-species N-mixture model. Each row corresponds to the regression coefficients for a given species.}

\item{kappa}{a numeric vector of length \code{n.sp} containing the dispersion parameter for the abundance portion of the N-mixture model for each species. Only relevant when \code{family = 'NB'}.}

\item{mu.RE}{a list used to specify the non-spatial random effects included in the abundance portion of the model. The list must have two tags: \code{levels} and \code{sigma.sq.mu}. \code{levels} is a vector of length equal to the number of distinct random effects to include in the model and contains the number of levels there are in each effect. \code{sigma.sq.mu} is a vector of length equal to the number of distinct random effects to include in the model and contains the variances for each random effect. If not specified, no random effects are included in the abundance portion of the model. An optional third tag, \code{beta.indx}, is a list that contains integers denoting the corresponding value of \code{beta} that each random effect corresponds to. This allows specification of random intercepts as well as slopes. By default, all effects are assumed to be random intercepts.}

\item{p.RE}{a list used to specify the non-spatial random effects included in the detection portion of the model. The list must have two tags: \code{levels} and \code{sigma.sq.p}. \code{levels} is a vector of length equal to the number of distinct random effects to include in the model and contains the number of levels there are in each effects. \code{sigma.sq.p} is a vector of length equal to the number of distinct random effects to include in the model and contains the variances for each random effect. If not specified, no random effects are included in the detection portion of the model. An optional third tag, \code{alpha.indx}, is a list that contains integers denoting the corresponding value of \code{alpha} that each random effect corresponds to. This allows specification of random intercepts as well as slopes. By default, all effects are assumed to be random intercepts.}

\item{offset}{either a single numeric value or a vector of length \code{J} that contains the offset for each location in the data set.}

\item{sp}{a logical value indicating whether to simulate a spatially-explicit model with a Gaussian process. By default set to \code{FALSE}.}

\item{cov.model}{a quoted keyword that specifies the covariance function used to model the spatial dependence structure among the latent abundance values.  Supported covariance model key words are: \code{"exponential"}, \code{"matern"}, \code{"spherical"}, and \code{"gaussian"}.}

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
}

\value{
  A list comprised of: 

  \item{X}{a \eqn{J \times p.abund}{J x p.abund} numeric design matrix for the abundance portion of the model.}
  \item{X.p}{a three-dimensional numeric array with dimensions corresponding to sites, repeat visits, and number of detection regression coefficients. This is the design matrix used for the detection portion of the N-mixture model.}
  \item{coords}{a \eqn{J \times 2}{J x 2} numeric matrix of coordinates of each site. Required for spatial models.}
  \item{w}{a \eqn{N \times J}{N x J} matrix of the spatial random effects for each species. Only used to simulate data when \code{sp = TRUE}. If \code{factor.model = TRUE}, the first dimension is \code{n.factors}.}
  \item{mu}{a \code{n.sp x J} matrix of the expected abundances for each species at each site.}
  \item{N}{a \code{n.sp x J} matrix of the latent occurrence states for each species at each site.}
  \item{p}{a \code{n.sp x J x max(n.rep)} array of the detection probabilities for each species at each site and replicate combination. Sites with fewer than \code{max(n.rep)} replicates will contain \code{NA} values.}
  \item{y}{a \code{n.sp x J x max(n.rep)} array of the raw count data for each species at each site and replicate combination. Sites with fewer than \code{max(n.rep)} replicates will contain \code{NA} values.}
  \item{X.p.re}{a three-dimensional numeric array containing the levels of any detection random effect included in the model. Only relevant when detection random effects are specified in \code{p.RE}.}
  \item{X.re}{a numeric matrix containing the levels of any abundance random effect included in the model. Only relevant when abundance random effects are specified in \code{mu.RE}.}
  \item{alpha.star}{a numeric matrix where each row contains the simulated detection random effects for each given level of the random effects included in the detection model. Only relevant when detection random effects are included in the model.}
  \item{beta.star}{a numeric matrix where each row contains the simulated abundance random effects for each given level of the random effects included in the abundance model. Only relevant when abundance random effects are included in the model.}
}

\examples{
J.x <- 8
J.y <- 8
J <- J.x * J.y
n.rep <- sample(2:4, size = J, replace = TRUE)
n.sp <- 10
# Community-level covariate effects
# Abundance 
beta.mean <- c(0.2, -0.15)
p.abund <- length(beta.mean)
tau.sq.beta <- c(0.6, 0.3)
# Detection
alpha.mean <- c(0.5, 0.2)
tau.sq.alpha <- c(0.2, 0.3)
p.det <- length(alpha.mean)
mu.RE <- list(levels = c(10, 12), 
              sigma.sq.mu = c(1.5, 0.3), 
              beta.indx = list(1, 2))
p.RE <- list(levels = c(15, 10), 
             sigma.sq.p = c(0.8, 0.5), 
             alpha.indx = list(1, 2))
# Draw species-level effects from community means.
beta <- matrix(NA, nrow = n.sp, ncol = p.abund)
alpha <- matrix(NA, nrow = n.sp, ncol = p.det)
for (i in 1:p.abund) {
  beta[, i] <- rnorm(n.sp, beta.mean[i], sqrt(tau.sq.beta[i]))
}
for (i in 1:p.det) {
  alpha[, i] <- rnorm(n.sp, alpha.mean[i], sqrt(tau.sq.alpha[i]))
}
factor.model <- TRUE
n.factors <- 3
# Spatial parameters if desired
phi <- runif(n.factors, 3/1, 3/.1)
sp <- TRUE
family <- 'Poisson'

dat <- simMsNMix(J.x = J.x, J.y = J.y, n.rep = n.rep, n.sp = n.sp, beta = beta, 
                alpha = alpha, mu.RE = mu.RE, p.RE = p.RE, sp = TRUE, 
                cov.model = 'exponential', phi = phi, factor.model = factor.model,
                n.factors = n.factors, family = family)
}
