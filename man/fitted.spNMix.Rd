\name{fitted.spNMix}

\alias{fitted.spNMix}

\title{Extract Model Fitted Values for spNMix Object}

\description{
  Method for extracting model fitted values and detection probabilities from a fitted spatial N-mixture (\code{spNMix}) model. 
}

\usage{
\method{fitted}{spNMix}(object, type = 'marginal', \dots)
}

\arguments{
  \item{object}{object of class \code{spNMix}.}
  \item{type}{a character string indicating whether fitted values should be generated conditional on the estimated latent abundance values (\code{type = 'conditional'}) estimated during the model or based on the marginal expected abundance values (\code{type = 'marginal'}).}
  \item{\dots}{currently no additional arguments}
}

\details{
  A method to the generic \code{\link[stats]{fitted}} function to extract fitted values and detection probabilities for fitted model objects of class \code{spNMix}. 
}

\value{
  A list comprised of: 

  \item{y.rep.samples}{A three-dimensional numeric array of fitted values for use in Goodness of Fit assessments. Array dimensions correspond to MCMC samples, sites, and replicates.}

  \item{p.samples}{A three-dimensional numeric array of detection probability values. Array dimensions correspond to MCMC samples, sites, and replicates.}
}

\keyword{model}

