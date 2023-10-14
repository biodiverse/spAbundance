\name{fitted.svcAbund}

\alias{fitted.svcAbund}

\title{Extract Model Fitted Values for svcAbund Object}

\description{
  Method for extracting model fitted values from a fitted spatially-varying coefficient GLMM (\code{svcAbund}). 
}

\usage{
\method{fitted}{svcAbund}(object, \dots)
}

\arguments{
  \item{object}{object of class \code{svcAbund}.}
  \item{\dots}{currently no additional arguments}
}

\details{
  A method to the generic \code{\link[stats]{fitted}} function to extract fitted values for fitted model objects of class \code{svcAbund}. 
}

\value{
  A three-dimensional numeric array of fitted values for use in Goodness of Fit assessments. Array dimensions correspond to MCMC samples, sites, and replicates
}

\keyword{model}

