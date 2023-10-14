\name{fitted.sfMsAbund}

\alias{fitted.sfMsAbund}

\title{Extract Model Fitted Values for sfMsAbund Object}

\description{
  Method for extracting model fitted values from a fitted spatial factor multivariate GLMM (\code{sfMsAbund}). 
}

\usage{
\method{fitted}{sfMsAbund}(object, \dots)
}

\arguments{
  \item{object}{object of class \code{sfMsAbund}.}
  \item{\dots}{currently no additional arguments}
}

\details{
  A method to the generic \code{\link[stats]{fitted}} function to extract fitted values for fitted model objects of class \code{sfMsAbund}. 
}

\value{
  A four-dimensional numeric array of fitted values for use in Goodness of Fit assessments. Array dimensions correspond to MCMC samples, species, sites, and replicates.
}

\keyword{model}

