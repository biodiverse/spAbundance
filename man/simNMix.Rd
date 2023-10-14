\name{simNMix}
\alias{simNMix}
\title{Simulate Single-Species Count Data with Imperfect Detection}

\description{
  The function \code{simNMix} simulates single-species count data for simulation studies, power assessments, or function testing. Data can be optionally simulated with a spatial Gaussian Process in the abundance portion of the model. Non-spatial random intercepts/slopes can also be included in the detection or abundance portions of the N-mixture model.
}

\usage{
simNMix(J.x, J.y, n.rep, n.rep.max, beta, alpha, kappa, mu.RE = list(), 
        p.RE = list(), offset = 1, sp = FALSE, cov.model, sigma.sq, phi, nu, 
        family = 'Poisson', ...)
}

\arguments{
\item{J.x}{a single numeric value indicating the number of sites to simulate count data along the horizontal axis. Total number of sites with simulated data is \eqn{J.x \times J.y}{J.x x J.y}. }

\item{J.y}{a single numeric value indicating the number of sites to simulate count data along the vertical axis. Total number of sites with simulated data is \eqn{J.x \times J.y}{J.x x J.y}. }

\item{n.rep}{a numeric vector of length \eqn{J = J.x \times J.y}{J = J.x x J.y} indicating the number of repeat visits at each of the \eqn{J}{J} sites.}

\item{n.rep.max}{a single numeric value indicating the maximum number of replicate surveys. This is an optional argument, with its default value set to \code{max(n.rep)}. This can be used to generate data sets with different types of missingness (e.g., simulate data across 20 days (replicate surveys) but sites are only sampled a maximum of ten times each).}

\item{beta}{a numeric vector containing the intercept and regression coefficient parameters for the abundance portion of the single-species N-mixture model.}

\item{alpha}{a numeric vector containing the intercept and regression coefficient parameters for the detection portion of the single-species N-mixture model.}

\item{kappa}{a single numeric value containing the dispersion parameter for the abundance portion of the N-mixture model. Only relevant when \code{family = 'NB'}.}

\item{mu.RE}{a list used to specify the non-spatial random effects included in the abundance portion of the model. The list must have two tags: \code{levels} and \code{sigma.sq.mu}. \code{levels} is a vector of length equal to the number of distinct random effects to include in the model and contains the number of levels there are in each effect. \code{sigma.sq.mu} is a vector of length equal to the number of distinct random effects to include in the model and contains the variances for each random effect. If not specified, no random effects are included in the abundance portion of the model. An optional third tag, \code{beta.indx}, is a list that contains integers denoting the corresponding value of \code{beta} that each random effect corresponds to. This allows specification of random intercepts as well as slopes. By default, all effects are assumed to be random intercepts.}

\item{p.RE}{a list used to specify the non-spatial random effects included in the detection portion of the model. The list must have two tags: \code{levels} and \code{sigma.sq.p}. \code{levels} is a vector of length equal to the number of distinct random effects to include in the model and contains the number of levels there are in each effects. \code{sigma.sq.p} is a vector of length equal to the number of distinct random effects to include in the model and contains the variances for each random effect. If not specified, no random effects are included in the detection portion of the model. An optional third tag, \code{alpha.indx}, is a list that contains integers denoting the corresponding value of \code{alpha} that each random effect corresponds to. This allows specification of random intercepts as well as slopes. By default, all effects are assumed to be random intercepts.}

\item{offset}{either a single numeric value or a vector of length \code{J} that contains the offset for each location in the data set.}

\item{sp}{a logical value indicating whether to simulate a spatially-explicit N-mixture model with a Gaussian process. By default set to \code{FALSE}.}

\item{cov.model}{a quoted keyword that specifies the covariance function used to model the spatial dependence structure among the latent abundance values.  Supported covariance model key words are: \code{"exponential"}, \code{"matern"}, \code{"spherical"}, and \code{"gaussian"}.}

\item{sigma.sq}{a numeric value indicating the spatial variance parameter. Ignored when \code{sp = FALSE}.}

\item{phi}{a numeric value indicating the spatial decay parameter. Ignored when \code{sp = FALSE}. }

\item{nu}{a numeric value indicating the spatial smoothness parameter. Only used when \code{sp = TRUE} and \code{cov.model = "matern"}.} 

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
  \item{X.p}{a three-dimensional numeric array with dimensions corresponding to sites, repeat visits, and number of detection regression coefficients. This is the design matrix used for the detection portion of the N-mixture model.}
  \item{coords}{a \eqn{J \times 2}{J x 2} numeric matrix of coordinates of each site. Required for spatial models.}
  \item{w}{a \eqn{J \times 1}{J x 1} matrix of the spatial random effects. Only used to simulate data when \code{sp = TRUE}.}
  \item{mu}{a \eqn{J \times 1}{J x 1} matrix of the expected abundance values for each site.}
  \item{N}{a length \eqn{J}{J} vector of the latent abundances at each site.}
  \item{p}{a \code{J x max(n.rep)} matrix of the detection probabilities for each site and replicate combination. Sites with fewer than \code{max(n.rep)} replicates will contain \code{NA} values.}
  \item{y}{a \code{J x max(n.rep)} matrix of the raw count data for each site and replicate combination.}
  \item{X.p.re}{a three-dimensional numeric array containing the levels of any detection random effect included in the model. Only relevant when detection random effects are specified in \code{p.RE}.}
  \item{X.re}{a numeric matrix containing the levels of any abundance random effect included in the model. Only relevant when abundance random effects are specified in \code{mu.RE}.}
  \item{alpha.star}{a numeric vector that contains the simulated detection random effects for each given level of the random effects included in the detection model. Only relevant when detection random effects are included in the model.}
  \item{beta.star}{a numeric vector that contains the simulated abundance random effects for each given level of the random effects included in the N-mixture model. Only relevant when abundance random effects are included in the model.}
}

\examples{
set.seed(400)
J.x <- 10
J.y <- 10
n.rep <- rep(4, J.x * J.y)
beta <- c(0.5, -0.15)
alpha <- c(0.7, 0.4)
kappa <- 0.5
phi <- 3 / .6
sigma.sq <- 2
mu.RE <- list(levels = 10, sigma.sq.mu = 1.2)
p.RE <- list(levels = 15, sigma.sq.p = 0.8)
dat <- simNMix(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, alpha = alpha,
               kappa = kappa, mu.RE = mu.RE, p.RE = p.RE, sp = TRUE, 
               cov.model = 'spherical', sigma.sq = sigma.sq, phi = phi, 
               family = 'NB')
}
