\name{summary.svcAbund}

\alias{summary.svcAbund}
\alias{print.svcAbund}

\title{Methods for svcAbund Object}

\description{
  Methods for extracting information from fitted spatially-varying coefficient abundance GLMs (\code{svcAbund}).
}

\usage{
\method{summary}{svcAbund}(object, quantiles = c(0.025, 0.5, 0.975), 
        digits = max(3L, getOption("digits") - 3L), \dots)
\method{print}{svcAbund}(x, \dots)
}

\arguments{
  \item{object, x}{object of class \code{svcAbund}.}
  \item{quantiles}{for \code{summary}, posterior distribution quantiles to compute.}
  \item{digits}{for \code{summary}, number of digits to report.}
  \item{\dots}{currently no additional arguments}
}

\details{
  A set of standard extractor functions for fitted model objects of
  class \code{svcAbund}, including methods to the generic functions \code{\link[base]{print}} and \code{\link[base]{summary}}. 
}

\value{
  No return value, called to display summary information of a \code{svcAbund} object.
}


\keyword{model}
