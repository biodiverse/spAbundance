\name{fitted.lfMsAbund}

\alias{fitted.lfMsAbund}

\title{Extract Model Fitted Values for lfMsAbund Object}

\description{
  Method for extracting model fitted values from a fitted latent factor multivariate GLMM (\code{lfMsAbund}). 
}

\usage{
\method{fitted}{lfMsAbund}(object, \dots)
}

\arguments{
  \item{object}{object of class \code{lfMsAbund}.}
  \item{\dots}{currently no additional arguments}
}

\details{
  A method to the generic \code{\link[stats]{fitted}} function to extract fitted values for fitted model objects of class \code{lfMsAbund}. 
}

\value{
  A four-dimensional numeric array of fitted values for use in Goodness of Fit assessments. Array dimensions correspond to MCMC samples, species, sites, and replicates.
}

\keyword{model}

