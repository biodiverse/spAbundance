\name{fitted.abund}

\alias{fitted.abund}

\title{Extract Model Fitted Values for abund Object}

\description{
  Method for extracting model fitted values from a fitted GLMM (\code{abund}). 
}

\usage{
\method{fitted}{abund}(object, \dots)
}

\arguments{
  \item{object}{object of class \code{abund}.}
  \item{\dots}{currently no additional arguments}
}

\details{
  A method to the generic \code{\link[stats]{fitted}} function to extract fitted values for fitted model objects of class \code{abund}. 
}

\value{
  A three-dimensional numeric array of fitted values for use in Goodness of Fit assessments. Array dimensions correspond to MCMC samples, sites, and replicates
}

\keyword{model}

