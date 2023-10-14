\name{fitted.msAbund}

\alias{fitted.msAbund}

\title{Extract Model Fitted Values for msAbund Object}

\description{
  Method for extracting model fitted values from a fitted multivariate GLMM (\code{msAbund}). 
}

\usage{
\method{fitted}{msAbund}(object, \dots)
}

\arguments{
  \item{object}{object of class \code{msAbund}.}
  \item{\dots}{currently no additional arguments}
}

\details{
  A method to the generic \code{\link[stats]{fitted}} function to extract fitted values for fitted model objects of class \code{msAbund}. 
}

\value{
  A four-dimensional numeric array of fitted values for use in Goodness of Fit assessments. Array dimensions correspond to MCMC samples, species, sites, and replicates.
}

\keyword{model}

