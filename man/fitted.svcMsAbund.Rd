\name{fitted.svcMsAbund}

\alias{fitted.svcMsAbund}

\title{Extract Model Fitted Values for svcMsAbund Object}

\description{
  Method for extracting model fitted values from a fitted multivatiate spatially-varying coefficient GLMM (\code{svcMsAbund}). 
}

\usage{
\method{fitted}{svcMsAbund}(object, \dots)
}

\arguments{
  \item{object}{object of class \code{svcMsAbund}.}
  \item{\dots}{currently no additional arguments}
}

\details{
  A method to the generic \code{\link[stats]{fitted}} function to extract fitted values for fitted model objects of class \code{svcMsAbund}. 
}

\value{
  A four-dimensional numeric array of fitted values for use in Goodness of Fit assessments. Array dimensions correspond to MCMC samples, species, sites, and replicates.
}

\keyword{model}

