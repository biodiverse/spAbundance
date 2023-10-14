\name{simDS}
\alias{simDS}
\title{Simulate Single-Species Distance Sampling Data}

\description{
  The function \code{simDS} simulates single-species distance sampling data for simulation studies, power assessments, or function testing. Data can be optionally simulated with a spatial Gaussian Process in the abundance portion of the model. Non-spatial random effects can also be included in the detection or abundance portions of the distance sampling model.
}

\usage{
simDS(J.x, J.y, n.bins, bin.width, beta, alpha, det.func, transect = 'line', 
      kappa, mu.RE = list(), p.RE = list(), offset = 1, 
      sp = FALSE, cov.model, sigma.sq, phi, nu, family = 'Poisson', ...)
}

\arguments{
\item{J.x}{a single numeric value indicating the number of sites to simulate count data along the horizontal axis. Total number of sites with simulated data is \eqn{J.x \times J.y}{J.x x J.y}. }

\item{J.y}{a single numeric value indicating the number of sites to simulate count data along the vertical axis. Total number of sites with simulated data is \eqn{J.x \times J.y}{J.x x J.y}. }

\item{n.bins}{a single numeric value indicating the number of distance bins from which to generate data.}

\item{bin.width}{a vector of length \code{n.bins} indicating the length of each bin. Lengths can be different for each distance bin or the same across bins.}

\item{beta}{a numeric vector containing the intercept and regression coefficient parameters for the abundance portion of the single-species distance sampling model.}

\item{alpha}{a numeric vector containing the intercept and regression coefficient parameters for the detection portion of the single-species distance sampling model.}

\item{det.func}{the detection model used to describe how detection probability varies
  with distance. In other software, this is often referred to as the key function. Currently
  supports two functions: half normal (\code{'halfnormal'}) and negative exponential
  (\code{'negexp'}).}

\item{transect}{the type of transect. Currently supports line transects (\code{'line'}) or
  circular transects (i.e., point counts; \code{'point'}).}

\item{kappa}{a single numeric value containing the dispersion parameter for the abundance portion of the hierarchical distance sampling model. Only relevant when \code{family = 'NB'}.}

\item{mu.RE}{a list used to specify the non-spatial random intercepts included in the abundance portion of the model. The list must have two tags: \code{levels} and \code{sigma.sq.mu}. \code{levels} is a vector of length equal to the number of distinct random intercepts to include in the model and contains the number of levels there are in each intercept. \code{sigma.sq.mu} is a vector of length equal to the number of distinct random intercepts to include in the model and contains the variances for each random effect. If not specified, no random effects are included in the abundance portion of the model.}

\item{p.RE}{a list used to specify the non-spatial random intercepts included in the detection portion of the model. The list must have two tags: \code{levels} and \code{sigma.sq.p}. \code{levels} is a vector of length equal to the number of distinct random intercepts to include in the model and contains the number of levels there are in each intercept. \code{sigma.sq.p} is a vector of length equal to the number of distinct random intercepts to include in the model and contains the variances for each random effect. If not specified, no random effects are included in the detection portion of the model.}

\item{offset}{either a single numeric value or a vector of length \code{J} that contains the offset for each location in the data set.}

\item{sp}{a logical value indicating whether to simulate a spatially-explicit HDS model with a Gaussian process. By default set to \code{FALSE}.}

\item{cov.model}{a quoted keyword that specifies the covariance function used to model the spatial dependence structure among the latent abundance values.  Supported covariance model key words are: \code{"exponential"}, \code{"matern"}, \code{"spherical"}, and \code{"gaussian"}.}

\item{sigma.sq}{a numeric value indicating the spatial variance parameter. Ignored when \code{sp = FALSE}.}

\item{phi}{a numeric value indicating the spatial decay parameter. Ignored when \code{sp = FALSE}. }

\item{nu}{a numeric value indicating the spatial smoothness parameter. Only used when \code{sp = TRUE} and \code{cov.model = "matern"}.} 

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

  \item{X}{a \eqn{J \times p.abund}{J x p.abund} numeric design matrix for the abundance portion of the model.}
  \item{X.p}{a \eqn{J \times p.abund}{J x p.det} numeric design matrix for the detection portion of the model.}
  \item{coords}{a \eqn{J \times 2}{J x 2} numeric matrix of coordinates of each site. Required for spatial models.}
  \item{w}{a \eqn{J \times 1}{J x 1} matrix of the spatial random effects. Only used to simulate data when \code{sp = TRUE}.}
  \item{mu}{a \eqn{J \times 1}{J x 1} matrix of the expected abundance values for each site.}
  \item{N}{a length \eqn{J}{J} vector of the latent abundances at each site.}
  \item{p}{ a length J vector of the detection probabilities at each site.}
  \item{pi.full}{a \code{J x n.bins + 1} vector of the cell-specific detection probabilities for each site, where the last column indicates the probability of not detecting an individual at that site.}
  \item{y}{a \code{J x max(n.bins)} matrix of the raw count data for each site and distance bin.}
  \item{X.p.re}{a numeric matrix containing the levels of any detection random effect included in the model. Only relevant when detection random effects are specified in \code{p.RE}.}
  \item{X.re}{a numeric matrix containing the levels of any abundance random effect included in the model. Only relevant when abundance random effects are specified in \code{mu.RE}.}
  \item{alpha.star}{a numeric vector that contains the simulated detection random effects for each given level of the random effects included in the detection model. Only relevant when detection random effects are included in the model.}
  \item{beta.star}{a numeric vector that contains the simulated abundance random effects for each given level of the random effects included in the HDS model. Only relevant when abundance random effects are included in the model.}
}

\examples{
set.seed(110)
J.x <- 10
J.y <- 10 
J <- J.x * J.y
# Number of distance bins from which to simulate data. 
n.bins <- 5
# Length of each bin. This should be of length n.bins
bin.width <- c(.10, .10, .20, .3, .1)
# Abundance coefficients
beta <- c(1.0, 0.2, 0.3, -0.2)
p.abund <- length(beta)
# Detection coefficients
alpha <- c(-1.0, -0.3)
p.det <- length(alpha)
# Detection decay function
det.func <- 'halfnormal'
mu.RE <- list()
p.RE <- list()
sp <- FALSE 
family <- 'NB'
kappa <- 0.1
offset <- 1.8
transect <- 'point'

dat <- simDS(J.x = J.x, J.y = J.y, n.bins = n.bins, bin.width = bin.width,
             beta = beta, alpha = alpha, det.func = det.func, kappa = kappa, 
             mu.RE = mu.RE, p.RE = p.RE, sp = sp,
             sigma.sq = sigma.sq, phi = phi, nu = nu, family = family, 
             offset = offset, transect = transect)
}
