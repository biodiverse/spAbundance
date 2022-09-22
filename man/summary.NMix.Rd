\name{summary.NMix}

\alias{summary.NMix}
\alias{print.NMix}

\title{Methods for NMix Object}

\description{
  Methods for extracting information from fitted single-species N-mixture (\code{NMix}) model.
}

\usage{
\method{summary}{NMix}(object, quantiles = c(0.025, 0.5, 0.975), 
        digits = max(3L, getOption("digits") - 3L), \dots)
\method{print}{NMix}(x, \dots)
}

\arguments{
  \item{object, x}{object of class \code{NMix}.}
  \item{quantiles}{for \code{summary}, posterior distribution quantiles to compute.}
  \item{digits}{for \code{summary}, number of digits to report.}
  \item{\dots}{currently no additional arguments}
}

\details{
  A set of standard extractor functions for fitted model objects of
  class \code{NMix}, including methods to the generic functions \code{\link[base]{print}} and \code{\link[base]{summary}}. 
}

\value{
  No return value, called to display summary information of a \code{NMix} object.
}


\keyword{model}
