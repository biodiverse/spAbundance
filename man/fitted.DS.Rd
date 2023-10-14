\name{fitted.DS}

\alias{fitted.DS}

\title{Extract Model Fitted Values for DS Object}

\description{
  Method for extracting model fitted values and cell-specific detection probabilities from a hierarchical distance sampling (\code{DS}) model. 
}

\usage{
\method{fitted}{DS}(object, \dots)
}

\arguments{
  \item{object}{object of class \code{DS}.}
  \item{\dots}{currently no additional arguments}
}

\details{
  A method to the generic \code{\link[stats]{fitted}} function to extract fitted values and detection probabilities for fitted model objects of class \code{DS}. 
}

\value{
  A list comprised of: 

  \item{y.rep.samples}{A three-dimensional numeric array of fitted values for use in Goodness of Fit assessments. Array dimensions correspond to MCMC samples, sites, and distance bin.}

  \item{pi.samples}{A three-dimensional numeric array of cell-specific detection probability values. Values correspond to the probability of detecting an individual within a given distance band at a given location. Array dimensions correspond to MCMC samples, sites, and distance band.}
}

\keyword{model}

