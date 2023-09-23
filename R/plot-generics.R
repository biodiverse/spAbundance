overallPlot <- function(x, param, density = TRUE, ...) {
  n.post <- x$n.post
  n.chains <- x$n.chains
  curr.samples <- vector(mode = 'list', length = n.chains) 
  indx <- 1:n.post
  if (param == 'beta.comm') {
    if (class(x) %in% c('abund', 'spAbund', 'svcAbund', 'NMix', 
			     'spNMix', 'DS', 'spDS')) {
      stop("beta.comm is not a parameter in the fitted model")
    }
    for (i in 1:n.chains) {
      curr.samples[[i]] <- coda::mcmc(x$beta.comm.samples[indx, , drop = FALSE])
      indx <- (n.post * i + 1):(n.post * (i + 1))
    }
  }
  if (param == 'tau.sq.beta') {
    if (class(x) %in% c('abund', 'spAbund', 'svcAbund', 'NMix', 
			     'spNMix', 'DS', 'spDS')) {
      stop("tau.sq.beta is not a parameter in the fitted model")
    }
    for (i in 1:n.chains) {
      curr.samples[[i]] <- coda::mcmc(x$tau.sq.beta.samples[indx, , drop = FALSE])
      indx <- (n.post * i + 1):(n.post * (i + 1))
    }
  }
  if (param == 'lambda') {
    if (class(x) %in% c('abund', 'spAbund', 'svcAbund', 'NMix', 
			     'spNMix', 'DS', 'spDS', 'msAbund', 'msNMix', 'msDS')) {
      stop("lambda is not a parameter in the fitted model")
    }
    for (i in 1:n.chains) {
      curr.samples[[i]] <- coda::mcmc(x$lambda.samples[indx, , drop = FALSE])
      indx <- (n.post * i + 1):(n.post * (i + 1))
    }
  }
  if (param == 'tau.sq') {
    if (!(class(x) %in% c('abund', 'spAbund', 'svcAbund', 'svcMsAbund', 
			     'msAbund', 'lfMsAbund', 'sfMsAbund'))) {
      stop("tau.sq is not a parameter in the fitted model")
    }
    for (i in 1:n.chains) {
      curr.samples[[i]] <- coda::mcmc(x$tau.sq.samples[indx, , drop = FALSE])
      indx <- (n.post * i + 1):(n.post * (i + 1))
    }
  }
  if (param == 'alpha.comm') {
    if (class(x) %in% c('abund', 'spAbund', 'msAbund', 'lfMsAbund', 'sfMsAbund', 
			     'svcAbund', 'svcMsAbund', 'NMix', 'spNMix', 'DS', 'spDS')) {
      stop("alpha.comm is not a parameter in the fitted model")
    }
    for (i in 1:n.chains) {
      curr.samples[[i]] <- coda::mcmc(x$alpha.comm.samples[indx, , drop = FALSE])
      indx <- (n.post * i + 1):(n.post * (i + 1))
    }
  }
  if (param == 'tau.sq.alpha') {
    if (class(x) %in% c('abund', 'spAbund', 'msAbund', 'lfMsAbund', 'sfMsAbund', 
			     'svcAbund', 'svcMsAbund', 'NMix', 'spNMix', 'DS', 'spDS')) {
      stop("tau.sq.alpha is not a parameter in the fitted model")
    }
    for (i in 1:n.chains) {
      curr.samples[[i]] <- coda::mcmc(x$tau.sq.alpha.samples[indx, , drop = FALSE])
      indx <- (n.post * i + 1):(n.post * (i + 1))
    }
  }
  if (param == 'theta') {
    if (!(class(x) %in% c('spAbund', 'sfMsAbund', 
			     'svcAbund', 'svcMsAbund', 
			     'spNMix', 'sfMsNMix', 
			     'spDS', 'sfMsDS'))) {
      stop("theta is not a parameter in the fitted model")
    }
    for (i in 1:n.chains) {
      curr.samples[[i]] <- coda::mcmc(x$theta.samples[indx, , drop = FALSE])
      indx <- (n.post * i + 1):(n.post * (i + 1))
    }
  }
  if (param == 'beta') {
    for (i in 1:n.chains) {
      curr.samples[[i]] <- coda::mcmc(x$beta.samples[indx, , drop = FALSE])
      indx <- (n.post * i + 1):(n.post * (i + 1))
    }
  }
  if (param == 'alpha') {
    if (class(x) %in% c('abund', 'spAbund', 'msAbund', 'lfMsAbund', 'sfMsAbund', 
			     'svcAbund', 'svcMsAbund')) {
      stop("alpha is not a parameter in the fitted model")
    }
    for (i in 1:n.chains) {
      curr.samples[[i]] <- coda::mcmc(x$alpha.samples[indx, , drop = FALSE])
      indx <- (n.post * i + 1):(n.post * (i + 1))
    }
  }
  if (param == 'sigma.sq.mu') {
    if (!x$muRE) {
      stop("sigma.sq.mu is not a parameter in the fitted model")
    }
    for (i in 1:n.chains) {
      curr.samples[[i]] <- coda::mcmc(x$sigma.sq.mu.samples[indx, , drop = FALSE])
      indx <- (n.post * i + 1):(n.post * (i + 1))
    }
  }
  if (param == 'sigma.sq.p') {
    if (class(x) %in% c('abund', 'spAbund', 'msAbund', 'lfMsAbund', 'sfMsAbund', 
			     'svcAbund', 'svcMsAbund')) {
      stop("N is not a parameter in the fitted model")
    }
    if (!x$pRE) {
      stop("sigma.sq.p is not a parameter in the fitted model")
    }
    for (i in 1:n.chains) {
      curr.samples[[i]] <- coda::mcmc(x$sigma.sq.p.samples[indx, , drop = FALSE])
      indx <- (n.post * i + 1):(n.post * (i + 1))
    }
  }
  # if (param == 'N') {
  #   if (class(x) %in% c('abund', 'spAbund', 'msAbund', 'lfMsAbund', 'sfMsAbund', 
  #       		     'svcAbund', 'svcMsAbund')) {
  #     stop("N is not a parameter in the fitted model")
  #   }
  #   for (i in 1:n.chains) {
  #     curr.samples[[i]] <- coda::mcmc(x$N.samples[indx, , drop = FALSE])
  #     indx <- (n.post * i + 1):(n.post * (i + 1))
  #   }
  # }
  # if (param == 'mu') {
  #   for (i in 1:n.chains) {
  #     curr.samples[[i]] <- coda::mcmc(x$mu.samples[indx, , drop = FALSE])
  #     indx <- (n.post * i + 1):(n.post * (i + 1))
  #   }
  # }
  if (param == 'beta.star') {
    if (!x$muRE) {
      stop("the model was not fit with any abundance random effects (beta.star)")
    }
    for (i in 1:n.chains) {
      curr.samples[[i]] <- coda::mcmc(x$beta.star.samples[indx, , drop = FALSE])
      indx <- (n.post * i + 1):(n.post * (i + 1))
    }
  }
  if (param == 'alpha.star') {
    if (class(x) %in% c('abund', 'spAbund', 'msAbund', 'lfMsAbund', 'sfMsAbund', 
			     'svcAbund', 'svcMsAbund')) {
      stop("alpha.star is not a parameter in the fitted model")
    }
    if (!x$pRE) {
      stop("the model was not fit with any detection random effects (alpha.star)")
    }
    for (i in 1:n.chains) {
      curr.samples[[i]] <- coda::mcmc(x$alpha.star.samples[indx, , drop = FALSE])
      indx <- (n.post * i + 1):(n.post * (i + 1))
    }
  }
  curr.samples <- coda::mcmc.list(curr.samples)
  plot(curr.samples, density = density)
}

plot.NMix <- function(x, param, density = TRUE, ...) {
  overallPlot(x, param, density)
}
plot.spNMix <- function(x, param, density = TRUE, ...) {
  overallPlot(x, param, density)
}
plot.msNMix <- function(x, param, density = TRUE, ...) {
  overallPlot(x, param, density)
}
plot.lfMsNMix <- function(x, param, density = TRUE, ...) {
  overallPlot(x, param, density)
}
plot.sfMsNMix <- function(x, param, density = TRUE, ...) {
  overallPlot(x, param, density)
}
plot.DS <- function(x, param, density = TRUE, ...) {
  overallPlot(x, param, density)
}
plot.spDS <- function(x, param, density = TRUE, ...) {
  overallPlot(x, param, density)
}
plot.msDS <- function(x, param, density = TRUE, ...) {
  overallPlot(x, param, density)
}
plot.lfMsDS <- function(x, param, density = TRUE, ...) {
  overallPlot(x, param, density)
}
plot.sfMsDS <- function(x, param, density = TRUE, ...) {
  overallPlot(x, param, density)
}
plot.abund <- function(x, param, density = TRUE, ...) {
  overallPlot(x, param, density)
}
plot.spAbund <- function(x, param, density = TRUE, ...) {
  overallPlot(x, param, density)
}
plot.msAbund <- function(x, param, density = TRUE, ...) {
  overallPlot(x, param, density)
}
plot.lfMsAbund <- function(x, param, density = TRUE, ...) {
  overallPlot(x, param, density)
}
plot.sfMsAbund <- function(x, param, density = TRUE, ...) {
  overallPlot(x, param, density)
}
plot.svcAbund <- function(x, param, density = TRUE, ...) {
  overallPlot(x, param, density)
}
plot.svcMsAbund <- function(x, param, density = TRUE, ...) {
  overallPlot(x, param, density)
}
plot.svcAbund <- function(x, param, density = TRUE, ...) {
  overallPlot(x, param, density)
}
