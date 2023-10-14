\name{fitted.spAbund}

\alias{fitted.spAbund}

\title{Extract Model Fitted Values for spAbund Object}

\description{
  Method for extracting model fitted values from a fitted spatial GLMM (\code{spAbund}). 
}

\usage{
\method{fitted}{spAbund}(object, \dots)
}

\arguments{
  \item{object}{object of class \code{spAbund}.}
  \item{\dots}{currently no additional arguments}
}

\details{
  A method to the generic \code{\link[stats]{fitted}} function to extract fitted values for fitted model objects of class \code{spAbund}. 
}

\value{
  A three-dimensional numeric array of fitted values for use in Goodness of Fit assessments. Array dimensions correspond to MCMC samples, sites, and replicates
}

\keyword{model}

