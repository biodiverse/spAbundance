\name{summary.spNMix}

\alias{summary.spNMix}
\alias{print.spNMix}

\title{Methods for spNMix Object}

\description{
  Methods for extracting information from fitted single-species spatial N-mixture (\code{spNMix}) models.
}

\usage{
\method{summary}{spNMix}(object, quantiles = c(0.025, 0.5, 0.975), 
        digits = max(3L, getOption("digits") - 3L), \dots)
\method{print}{spNMix}(x, \dots)
}

\arguments{
  \item{object, x}{object of class \code{spNMix}.}
  \item{quantiles}{for \code{summary}, posterior distribution quantiles to compute.}
  \item{digits}{for \code{summary}, number of digits to report.}
  \item{\dots}{currently no additional arguments}
}

\details{
  A set of standard extractor functions for fitted model objects of
  class \code{spNMix}, including methods to the generic functions \code{\link[base]{print}} and \code{\link[base]{summary}}. 
}

\value{
  No return value, called to display summary information of a \code{spNMix} object.
}


\keyword{model}
