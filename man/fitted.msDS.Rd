\name{fitted.msDS}

\alias{fitted.msDS}

\title{Extract Model Fitted Values for msDS Object}

\description{
  Method for extracting model fitted values and cell-specific detection probabilities from a multi-species hierarchical distance sampling (\code{msDS}) model. 
}

\usage{
\method{fitted}{msDS}(object, \dots)
}

\arguments{
  \item{object}{object of class \code{msDS}.}
  \item{\dots}{currently no additional arguments}
}

\details{
  A method to the generic \code{\link[stats]{fitted}} function to extract fitted values and detection probabilities for fitted model objects of class \code{msDS}. 
}

\value{
  A list comprised of: 

  \item{y.rep.samples}{A four-dimensional numeric array of fitted values for use in Goodness of Fit assessments. Array dimensions correspond to MCMC samples, species, sites, and distance bin.}

  \item{pi.samples}{A four-dimensional numeric array of cell-specific detection probability values. Values correspond to the probability of detecting an individual within a given distance band at a given location. Array dimensions correspond to MCMC samples, species, sites, and distance band.}
}

\keyword{model}

