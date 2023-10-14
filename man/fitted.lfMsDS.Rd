\name{fitted.lfMsDS}

\alias{fitted.lfMsDS}

\title{Extract Model Fitted Values for lfMsDS Object}

\description{
  Method for extracting model fitted values and cell-specific detection probabilities from a latent factor multi-species hierarchical distance sampling (\code{lfMsDS}) model. 
}

\usage{
\method{fitted}{lfMsDS}(object, \dots)
}

\arguments{
  \item{object}{object of class \code{lfMsDS}.}
  \item{\dots}{currently no additional arguments}
}

\details{
  A method to the generic \code{\link[stats]{fitted}} function to extract fitted values and detection probabilities for fitted model objects of class \code{lfMsDS}. 
}

\value{
  A list comprised of: 

  \item{y.rep.samples}{A four-dimensional numeric array of fitted values for use in Goodness of Fit assessments. Array dimensions correspond to MCMC samples, species, sites, and distance bin.}

  \item{pi.samples}{A four-dimensional numeric array of cell-specific detection probability values. Values correspond to the probability of detecting an individual within a given distance band at a given location. Array dimensions correspond to MCMC samples, species, sites, and distance band.}
}

\keyword{model}

