\name{simMsDS}
\alias{simMsDS}
\title{Simulate Multi-Species Distance Sampling Data}

\description{
  The function \code{simMsDS} simulates multi-species distance sampling data for simulation studies, power assessments, or function testing. Data can be optionally simulated with a spatial Gaussian Process in the abundance portion of the model, as well as an option to allow for species correlations using a factor modeling approach. Non-spatial random effects can also be included in the detection or abundance portions of the model.
}

\usage{
simMsDS(J.x, J.y, n.bins, bin.width, n.sp, beta, alpha, 
        det.func, transect = 'line', kappa, mu.RE = list(), 
        p.RE = list(), offset = 1, sp = FALSE, cov.model, 
        sigma.sq, phi, nu, family = 'Poisson',
        factor.model = FALSE, n.factors, ...)
}

\arguments{

\item{J.x}{a single numeric value indicating the number of sites to simulate count data along the horizontal axis. Total number of sites with simulated data is \eqn{J.x \times J.y}{J.x x J.y}. }

\item{J.y}{a single numeric value indicating the number of sites to simulate count data along the vertical axis. Total number of sites with simulated data is \eqn{J.x \times J.y}{J.x x J.y}. }

\item{n.bins}{a single numeric value indicating the number of distance bins from which to generate data.}

\item{bin.width}{a vector of length \code{n.bins} indicating the length of each bin. Lengths can be different for each distance bin or the same across bins.}

\item{n.sp}{a single numeric value indicating the number of species to simulate count data.}

\item{beta}{a numeric matrix with \code{n.sp} rows containing the intercept and regression coefficient parameters for the abundance portion of the multi-species hierarchical distance sampling (HDS) model. Each row corresponds to the regression coefficients for a given species.} 

\item{alpha}{a numeric matrix with \code{n.sp} rows containing the intercept and regression coefficient parameters for the detection portion of the multi-species HDS model. Each row corresponds to the regression coefficients for a given species.}

\item{det.func}{the detection model used to describe how detection probability varies
  with distance. In other software, this is often referred to as the key function. Currently
  supports two functions: half normal (\code{'halfnormal'}) and negative exponential
  (\code{'negexp'}).}

\item{transect}{the type of transect. Currently supports line transects (\code{'line'}) or
  circular transects (i.e., point counts; \code{'point'}).}

\item{kappa}{a numeric vector of length \code{n.sp} containing the dispersion parameter for the abundance portion of the HDS model for each species. Only relevant when \code{family = 'NB'}.}

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
  Jeffrey W. Doser \email{doserjef@msu.edu} \cr
}

\value{
  A list comprised of: 

  \item{X}{a \eqn{J \times p.abund}{J x p.abund} numeric design matrix for the abundance portion of the model.}
  \item{X.p}{a \eqn{J \times p.abund}{J x p.abund} numeric design matrix for the detection portion of the model.}
  \item{coords}{a \eqn{J \times 2}{J x 2} numeric matrix of coordinates of each site. Required for spatial models.}
  \item{w}{a \eqn{N \times J}{N x J} matrix of the spatial random effects for each species. Only used to simulate data when \code{sp = TRUE}. If \code{factor.model = TRUE}, the first dimension is \code{n.factors}.}
  \item{mu}{a \code{n.sp x J} matrix of the expected abundances for each species at each site.}
  \item{N}{a \code{n.sp x J} matrix of the latent occurrence states for each species at each site.}
  \item{p}{a \code{n.sp x J x max(n.rep)} array of the detection probabilities for each species at each site and replicate combination. Sites with fewer than \code{max(n.rep)} replicates will contain \code{NA} values.}
  \item{y}{a \code{n.sp x J x max(n.rep)} array of the raw distance sampling data for each species at each site and and distance bin.}
  \item{X.p.re}{a numeric matrix containing the levels of any detection random effect included in the model. Only relevant when detection random effects are specified in \code{p.RE}.}
  \item{X.re}{a numeric matrix containing the levels of any abundance random effect included in the model. Only relevant when abundance random effects are specified in \code{mu.RE}.}
  \item{alpha.star}{a numeric matrix where each row contains the simulated detection random effects for each given level of the random effects included in the detection model. Only relevant when detection random effects are included in the model.}
  \item{beta.star}{a numeric matrix where each row contains the simulated abundance random effects for each given level of the random effects included in the abundance model. Only relevant when abundance random effects are included in the model.}
}

\examples{
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
family <- 'NB'
kappa <- runif(n.sp, 0.3, 3) 
offset <- pi * .8^2
transect <- 'line'
factor.model <- FALSE

dat <- simMsDS(J.x = J.x, J.y = J.y, n.bins = n.bins, bin.width = bin.width,
               n.sp = n.sp, beta = beta, alpha = alpha, det.func = det.func, kappa = kappa, 
               mu.RE = mu.RE, p.RE = p.RE, sp = sp, cov.model = cov.model,
               sigma.sq = sigma.sq, phi = phi, nu = nu, family = family, 
               offset = offset, transect = transect, factor.model = factor.model)
}
