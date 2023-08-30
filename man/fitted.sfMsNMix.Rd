\name{fitted.sfMsNMix}

\alias{fitted.sfMsNMix}

\title{Extract Model Fitted Values for sfMsNMix Object}

\description{
  Method for extracting model fitted values and detection probability values from a fitted spatial factor multi-species N-mixture (\code{sfMsNMix}) model. 
}

\usage{
\method{fitted}{sfMsNMix}(object, type = 'marginal', \dots)
}

\arguments{
  \item{object}{object of class \code{sfMsNMix}.}
  \item{type}{a character string indicating whether fitted values should be generated conditional on the estimated latent abundance values (\code{type = 'conditional'}) estimated during the model or based on the marginal expected abundance values (\code{type = 'marginal'}).}
  \item{\dots}{currently no additional arguments}
}

\details{
  A method to the generic \code{\link[stats]{fitted}} function to extract fitted values and detection probability values for fitted model objects of class \code{sfMsNMix}. 
}

\value{
  A list comprised of: 

  \item{y.rep.samples}{A four-dimensional numeric array of fitted values for use in Goodness of Fit assessments. Array dimensions correspond to MCMC samples, species, sites, and replicates.}
  \item{p.samples}{A four-dimensional numeric array of detection probability values. Array dimensions correspond to MCMC samples, species, sites, and replicates.}
}

\keyword{model}


