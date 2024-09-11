svcMsAbund <- function(formula, data, inits, priors, tuning,
                       svc.cols = 1, cov.model = 'exponential', NNGP = TRUE,
                       n.neighbors = 15, search.type = "cb", n.factors,
                       n.batch, batch.length, accept.rate = 0.43, family = 'Gaussian',
                       n.omp.threads = 1, verbose = TRUE, n.report = 100,
                       n.burn = round(.10 * n.batch * batch.length),
                       n.thin = 1, n.chains = 1, ...){

  ptm <- proc.time()

  if (!(family) %in% c('Poisson', 'NB', 'Gaussian', 'zi-Gaussian')) {
    stop("family must be either 'Poisson', 'NB', 'Gaussian', or 'zi-Gaussian'")
  }
  if (family %in% c('Gaussian', 'zi-Gaussian')) {
    svcMsAbundGaussian(formula, data, inits, priors, tuning, svc.cols, cov.model,
                       NNGP, n.neighbors, search.type, n.factors, n.batch,
                       batch.length, accept.rate, family, n.omp.threads,
                       verbose, n.report, n.burn, n.thin, n.chains)
  } else {
    stop("svcMsAbund is currently only supported for Gaussian and zi-Gaussian")
  }
}
