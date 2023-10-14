\name{summary.svcMsAbund}

\alias{summary.svcMsAbund}
\alias{print.svcMsAbund}
\alias{plot.svcMsAbund}

\title{Methods for svcMsAbund Object}

\description{
  Methods for extracting information from fitted multivariate spatially-varying coefficient abundance GLMMs (\code{svcMsAbund}).
}

\usage{
\method{summary}{svcMsAbund}(object, level = 'both', quantiles = c(0.025, 0.5, 0.975), 
        digits = max(3L, getOption("digits") - 3L), \dots)
\method{print}{svcMsAbund}(x, \dots)
\method{plot}{svcMsAbund}(x, param, density = TRUE, \dots)
}

\arguments{
  \item{object, x}{object of class \code{svcMsAbund}.}
  \item{level}{a quoted keyword that indicates the level to summarize the 
    model results. Valid key words are: \code{"community"}, 
    \code{"species"}, or \code{"both"}.}
  \item{quantiles}{for \code{summary}, posterior distribution quantiles to compute.}
  \item{digits}{for \code{summary}, number of digits to report.}
  \item{param}{parameter name for which to generate a traceplot. Valid names are 
               \code{"beta", "beta.star", "sigma.sq.mu", "tau.sq", "beta.comm", 
                     "tau.sq.beta", "lambda", "theta"}.}
  \item{density}{logical value indicating whether to also generate a density plot for each parameter
                 in addition to the MCMC traceplot.}
  \item{\dots}{currently no additional arguments}
}

\details{
  A set of standard extractor functions for fitted model objects of
  class \code{svcMsAbund}, including methods to the generic functions \code{\link[base]{print}}, \code{\link[base]{summary}}, and \code{\link[base]{plot}}. 
}

\value{
  No return value, called to display summary information of a \code{svcMsAbund} object.
}

\keyword{model}
