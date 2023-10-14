\name{summary.spAbund}

\alias{summary.spAbund}
\alias{print.spAbund}
\alias{plot.spAbund}

\title{Methods for spAbund Object}

\description{
  Methods for extracting information from fitted univariate spatially-explicit GLMMs (\code{spAbund}).
}

\usage{
\method{summary}{spAbund}(object, quantiles = c(0.025, 0.5, 0.975), 
        digits = max(3L, getOption("digits") - 3L), \dots)
\method{print}{spAbund}(x, \dots)
\method{plot}{spAbund}(x, param, density = TRUE, \dots)
}

\arguments{
  \item{object, x}{object of class \code{spAbund}.}
  \item{quantiles}{for \code{summary}, posterior distribution quantiles to compute.}
  \item{digits}{for \code{summary}, number of digits to report.}
  \item{param}{parameter name for which to generate a traceplot. Valid names are 
               \code{"beta", "beta.star", "sigma.sq.mu", "tau.sq", "theta"}.}
  \item{density}{logical value indicating whether to also generate a density plot for each parameter
                 in addition to the MCMC traceplot.}
  \item{\dots}{currently no additional arguments}
}

\details{
  A set of standard extractor functions for fitted model objects of
  class \code{spAbund}, including methods to the generic functions \code{\link[base]{print}}, \code{\link[base]{summary}}, \code{\link[base]{plot}}. 
}

\value{
  No return value, called to display summary information of a \code{spAbund} object.
}


\keyword{model}
