\name{summary.lfMsNMix}

\alias{summary.lfMsNMix}
\alias{print.lfMsNMix}
\alias{plot.lfMsNMix}

\title{Methods for lfMsNMix Object}

\description{
  Methods for extracting information from fitted latent factor multi-species N-mixture (\code{lfMsNMix}) model.
}

\usage{
\method{summary}{lfMsNMix}(object, level = 'both', quantiles = c(0.025, 0.5, 0.975), 
        digits = max(3L, getOption("digits") - 3L), \dots)
\method{print}{lfMsNMix}(x, \dots)
\method{plot}{lfMsNMix}(x, param, density = TRUE, \dots)
}

\arguments{
  \item{object, x}{object of class \code{lfMsNMix}.}
  \item{level}{a quoted keyword that indicates the level to summarize the 
    model results. Valid key words are: \code{"community"}, 
    \code{"species"}, or \code{"both"}.}
  \item{quantiles}{for \code{summary}, posterior distribution quantiles to compute.}
  \item{digits}{for \code{summary}, number of digits to report.}
  \item{param}{parameter name for which to generate a traceplot. Valid names are 
               \code{"beta", "beta.star", "sigma.sq.mu", "beta.comm", "tau.sq.beta", 
		     "alpha", "alpha.star", "sigma.sq.p", "alpha.comm", "tau.sq.alpha", "lambda"}.}
  \item{density}{logical value indicating whether to also generate a density plot for each parameter
                 in addition to the MCMC traceplot.}
  \item{\dots}{currently no additional arguments}
}

\details{
  A set of standard extractor functions for fitted model objects of
  class \code{lfMsNMix}, including methods to the generic functions \code{\link[base]{print}}, \code{\link[base]{summary}}, and \code{\link[base]{plot}}. 
}

\value{
  No return value, called to display summary information of a \code{lfMsNMix} object.
}

\keyword{model}

