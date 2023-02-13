\name{summary.abund}

\alias{summary.abund}
\alias{print.abund}

\title{Methods for abund Object}

\description{
  Methods for extracting information from fitted abundance GLMs (\code{abund}).
}

\usage{
\method{summary}{abund}(object, quantiles = c(0.025, 0.5, 0.975), 
        digits = max(3L, getOption("digits") - 3L), \dots)
\method{print}{abund}(x, \dots)
}

\arguments{
  \item{object, x}{object of class \code{abund}.}
  \item{quantiles}{for \code{summary}, posterior distribution quantiles to compute.}
  \item{digits}{for \code{summary}, number of digits to report.}
  \item{\dots}{currently no additional arguments}
}

\details{
  A set of standard extractor functions for fitted model objects of
  class \code{abund}, including methods to the generic functions \code{\link[base]{print}} and \code{\link[base]{summary}}. 
}

\value{
  No return value, called to display summary information of a \code{abund} object.
}

\keyword{model}
