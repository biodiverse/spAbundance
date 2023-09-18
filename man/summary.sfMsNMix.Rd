\name{summary.sfMsNMix}

\alias{summary.sfMsNMix}
\alias{print.sfMsNMix}
\alias{plot.sfMsNMix}

\title{Methods for sfMsNMix Object}

\description{
  Methods for extracting information from fitted spatial factor multi-species N-mixture (\code{sfMsNMix}) model.
}

\usage{
\method{summary}{sfMsNMix}(object, level = 'both', quantiles = c(0.025, 0.5, 0.975), 
        digits = max(3L, getOption("digits") - 3L), \dots)
\method{print}{sfMsNMix}(x, \dots)
\method{plot}{sfMsNMix}(x, param, density = TRUE, \dots)
}

\arguments{
  \item{object, x}{object of class \code{sfMsNMix}.}
  \item{level}{a quoted keyword that indicates the level to summarize the 
    model results. Valid key words are: \code{"community"}, 
    \code{"species"}, or \code{"both"}.}
  \item{quantiles}{for \code{summary}, posterior distribution quantiles to compute.}
  \item{digits}{for \code{summary}, number of digits to report.}
  \item{param}{parameter name for which to generate a traceplot. Valid names are 
               \code{"beta", "beta.star", "sigma.sq.mu", "beta.comm", "tau.sq.beta", 
		     "alpha", "alpha.star", "sigma.sq.p", "alpha.comm", "tau.sq.alpha", "lambda", "theta"}.}
  \item{density}{logical value indicating whether to also generate a density plot for each parameter
                 in addition to the MCMC traceplot.}
  \item{\dots}{currently no additional arguments}
}

\details{
  A set of standard extractor functions for fitted model objects of
  class \code{sfMsNMix}, including methods to the generic functions \code{\link[base]{print}}, \code{\link[base]{summary}}, and \code{\link[base]{plot}}. 
}

\value{
  No return value, called to display summary information of a \code{sfMsNMix} object.
}

\keyword{model}
