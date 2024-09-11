simAbund <- function(J.x, J.y, n.rep, n.rep.max, beta, kappa, tau.sq, mu.RE = list(),
		     offset = 1, sp = FALSE, svc.cols = 1, cov.model, sigma.sq, phi,
		     nu, family = 'Poisson', z, trend = FALSE, x.positive = FALSE, ...) {

  # Check for unused arguments ------------------------------------------
  formal.args <- names(formals(sys.function(sys.parent())))
  elip.args <- names(list(...))
  for(i in elip.args){
      if(! i %in% formal.args)
          warning("'",i, "' is not an argument")
  }

  # Check function inputs -------------------------------------------------
  J <- J.x * J.y
  # n.rep -----------------------------
  if (family %in% c('Gaussian', 'zi-Gaussian')) {
    if (missing(n.rep)) {
      n.rep <- rep(1, J)
    }
  }
  if (missing(n.rep)) {
    stop("error: n.rep must be specified.")
  }
  if (length(n.rep) != J) {
    stop(paste("error: n.rep must be a vector of length ", J, sep = ''))
  }
  if (missing(n.rep.max)) {
    n.rep.max <- max(n.rep)
  }
  # beta ------------------------------
  if (missing(beta)) {
    stop("error: beta must be specified.")
  }
  # family ------------------------------
  if (! (family %in% c('NB', 'Poisson', 'Gaussian', 'zi-Gaussian'))) {
    stop("error: family must be one of: NB (negative binomial), Poisson, 'Gaussian', or 'zi-Gaussian'")
  }
  if (family %in% c('Gaussian', 'zi-Gaussian')) {
    if (n.rep.max != 1) {
      stop("n.rep must be one for all sites for Gaussian or zi-Gaussian models")
    }
  }
  # kappa -----------------------------
  if (family == 'NB') {
    if (missing(kappa)) {
      stop("error: kappa (overdispersion parameter) must be specified when family = 'NB'.")
    }
  }
  if (family == 'Poisson' & !missing(kappa)) {
    message("overdispersion parameter (kappa) is ignored when family == 'Poisson'")
  }
  # tau.sq ----------------------------
  if (family %in% c('Gaussian', 'zi-Gaussian')) {
    if (missing(tau.sq)) {
      stop('error: tau.sq (residual variance) must be specified when family is Gaussian or zi-Gaussian')
    }
  }
  # mu.RE ----------------------------
  names(mu.RE) <- tolower(names(mu.RE))
  if (!is.list(mu.RE)) {
    stop("error: if specified, mu.RE must be a list with tags 'levels' and 'sigma.sq.mu'")
  }
  if (length(names(mu.RE)) > 0) {
    if (!'sigma.sq.mu' %in% names(mu.RE)) {
      stop("error: sigma.sq.mu must be a tag in mu.RE with values for the abundance random effect variances")
    }
    if (!'levels' %in% names(mu.RE)) {
      stop("error: levels must be a tag in mu.RE with the number of random effect levels for each abundance random intercept.")
    }
    if (!'beta.indx' %in% names(mu.RE)) {
      mu.RE$beta.indx <- list()
      for (i in 1:length(mu.RE$sigma.sq.mu)) {
        mu.RE$beta.indx[[i]] <- 1
      }
    }
  }
  # Spatial parameters ----------------
  if (length(svc.cols) > 1 & !sp) {
    stop("error: if simulating data with spatially-varying coefficients, set sp = TRUE")
  }
  if (sp) {
    if(missing(sigma.sq)) {
      stop("error: sigma.sq must be specified when sp = TRUE")
    }
    if(missing(phi)) {
      stop("error: phi must be specified when sp = TRUE")
    }
    if(missing(cov.model)) {
      stop("error: cov.model must be specified when sp = TRUE")
    }
    cov.model.names <- c("exponential", "spherical", "matern", "gaussian")
    if(! cov.model %in% cov.model.names){
      stop("error: specified cov.model '",cov.model,"' is not a valid option; choose from ",
           paste(cov.model.names, collapse=", ", sep="") ,".")
    }
    if (cov.model == 'matern' & missing(nu)) {
      stop("error: nu must be specified when cov.model = 'matern'")
    }
    p.svc <- length(svc.cols)
    if (length(phi) != p.svc) {
      stop("error: phi must have the same number of elements as svc.cols")
    }
    if (length(sigma.sq) != p.svc) {
      stop("error: sigma.sq must have the same number of elements as svc.cols")
    }
    if (cov.model == 'matern') {
      if (length(nu) != p.svc) {
        stop("error: nu must have the same number of elements as svc.cols")
      }
    }
  }
  # z values --------------------------
  if (family == 'zi-Gaussian') {
    if (missing(z)) {
      stop('for a zero-inflated Gaussian model, you must supply the z values (binary 0s or 1s)')
    }
  }
  # Subroutines -----------------------------------------------------------
  # MVN
  rmvn <- function(n, mu=0, V = matrix(1)) {
    p <- length(mu)
    if(any(is.na(match(dim(V),p))))
      stop("Dimension problem!")
    D <- chol(V)
    t(matrix(rnorm(n*p), ncol=p)%*%D + rep(mu,rep(n,p)))
  }
  logit <- function(theta, a = 0, b = 1){log((theta-a)/(b-theta))}
  logit.inv <- function(z, a = 0, b = 1){b-(b-a)/(1+exp(z))}

  # Form abundance covariates (if any) ------------------------------------
  n.beta <- length(beta)
  X <- array(NA, dim = c(J, n.rep.max, n.beta))
  # Get index of surveyed replicates for each site.
  rep.indx <- list()
  for (j in 1:J) {
    rep.indx[[j]] <- sample(1:n.rep.max, n.rep[j], replace = FALSE)
  }
  X[, , 1] <- 1
  if (n.beta > 1) {
    if (trend) { # if simulating data with a trend
      for (j in 1:J) {
        X[j, rep.indx[[j]], 2] <- (scale(1:n.rep.max))[rep.indx[[j]]]
        if (n.beta > 2) {
          for (i in 3:n.beta) {
            if (x.positive) {
              X[j, rep.indx[[j]], i] <- runif(n.rep[j], 0, 5)
            } else {
              X[j, rep.indx[[j]], i] <- rnorm(n.rep[j])
            }
          }
        }
      }
    } else {
      for (i in 2:n.beta) {
        for (j in 1:J) {
          if (x.positive) {
            X[j, rep.indx[[j]], i] <- runif(n.rep[j], 0, 5)
          } else {
            X[j, rep.indx[[j]], i] <- rnorm(n.rep[j])
          }
        }
      } # i
    }
  }

  # Simulate spatial random effect ----------------------------------------
  # Matrix of spatial locations
  s.x <- seq(0, 1, length.out = J.x)
  s.y <- seq(0, 1, length.out = J.y)
  coords <- as.matrix(expand.grid(s.x, s.y))
  if (sp) {
    w.mat <- matrix(NA, J, p.svc)
    if (cov.model == 'matern') {
      theta <- cbind(phi, nu)
    } else {
      theta <- as.matrix(phi)
    }
    for (i in 1:p.svc) {
      Sigma <- mkSpCov(coords, as.matrix(sigma.sq[i]), as.matrix(0), theta[i, ], cov.model)
      # Random spatial process
      w.mat[, i] <- rmvn(1, rep(0, J), Sigma)
    }
    X.w <- X[, , svc.cols, drop = FALSE]
  } else {
    w.mat <- NA
    X.w <- NA
  }

  # Random effects --------------------------------------------------------
  if (length(mu.RE) > 0) {
    p.nmix.re <- length(unlist(mu.RE$beta.indx))
    # p.nmix.re <- length(unlist(mu.RE$levels))
    tmp <- sapply(mu.RE$beta.indx, length)
    re.col.indx <- unlist(lapply(1:length(mu.RE$beta.indx), function(a) rep(a, tmp[a])))
    sigma.sq.mu <- mu.RE$sigma.sq.mu
    n.nmix.re.long <- mu.RE$levels[re.col.indx]
    n.nmix.re <- sum(n.nmix.re.long)
    beta.star.indx <- rep(1:p.nmix.re, n.nmix.re.long)
    beta.star <- rep(0, n.nmix.re)
    X.random <- X[, , unlist(mu.RE$beta.indx), drop = FALSE]
    n.random <- dim(X.random)[3]
    X.re <- array(NA, dim = c(J, n.rep.max, length(unique(re.col.indx))))
    for (i in 1:length(unique(re.col.indx))) {
      for (j in 1:J) {
        X.re[j, rep.indx[[j]], i] <- sample(1:mu.RE$levels[i], n.rep[j], replace = TRUE)
      }
    }
    indx.mat <- X.re[, , re.col.indx, drop = FALSE]
    for (i in 1:p.nmix.re) {
      beta.star[which(beta.star.indx == i)] <- rnorm(n.nmix.re.long[i], 0,
						     sqrt(sigma.sq.mu[i]))
    }
    if (length(mu.RE$levels) > 1) {
      for (j in 2:length(mu.RE$levels)) {
        X.re[, , j] <- X.re[, , j] + max(X.re[, , j - 1], na.rm = TRUE)
      }
    }
    if (p.nmix.re > 1) {
      for (j in 2:p.nmix.re) {
        indx.mat[, , j] <- indx.mat[, , j] + max(indx.mat[, , j - 1], na.rm = TRUE)
      }
    }
    beta.star.sites <- matrix(NA, J, n.rep.max)
    for (j in 1:J) {
      for (k in rep.indx[[j]]) {
        beta.star.sites[j, k] <- beta.star[indx.mat[j, k, ]] %*% as.matrix(X.random[j, k, ])
      } # k
    } # j
  } else {
    X.re <- NA
    beta.star <- NA
  }

  # Data formation --------------------------------------------------------
  mu <- matrix(NA, J, n.rep.max)
  y <- matrix(NA, J, n.rep.max)
  # Offset ----------------------------
  # Single value
  if (length(offset) == 1) {
    offset <- matrix(offset, J, n.rep.max)
  } else if (length(dim(offset)) == 1) { # Value for each site
    if (length(offset) != J) {
      stop(paste0("offset must be a single value, vector of length ", J, " or a matrix with ",
	         J, " rows and ", n.rep.max, " columns."))
    }
    offset <- matrix(offset, J, n.rep.max)
  } else if (length(dim(offset)) == 2) { # Value for each site/obs
    if (nrow(offset) != J | ncol(offset) != n.rep.max) {
      stop(paste0("offset must be a single value, vector of length ", J, " or a matrix with ",
                  J, " rows and ", n.rep.max, " columns."))

    }
  }
  for (j in 1:J) {
    for (k in rep.indx[[j]]) {
      if (sp) {
        if (length(mu.RE) > 0) {
          mu[j, k] <- t(as.matrix(X[j, k, ])) %*% as.matrix(beta) +
      			   t(as.matrix(X.w[j, k, ])) %*% as.matrix(w.mat[j, ]) +
      			   beta.star.sites[j, k]
        } else {
          mu[j, k] <- t(as.matrix(X[j, k, ])) %*% as.matrix(beta) +
		      t(as.matrix(X.w[j, k, ])) %*% as.matrix(w.mat[j, ])
        }
      } else {
        if (length(mu.RE) > 0) {
          mu[j, k] <- t(as.matrix(X[j, k, ])) %*% as.matrix(beta) +
      			   beta.star.sites[j, k]
        } else {
          mu[j, k] <- t(as.matrix(X[j, k, ])) %*% as.matrix(beta)
        }
      }
      if (family %in% c('Poisson', 'NB')) {
        mu[j, k] <- exp(mu[j, k])
      }
      if (family == 'NB') {
        y[j, k] <- rnbinom(1, size = kappa, mu = mu[j, k] * offset[j, k])
      }
      if (family == 'Poisson') {
        y[j, k] <- rpois(1, lambda = mu[j, k] * offset[j, k])
      }
      if (family == 'Gaussian') {
        y[j, k] <- rnorm(1, mu[j, k] * offset[j, k], sqrt(tau.sq))
      }
      if (family == 'zi-Gaussian') {
        mu[j, k] <- mu[j, k] * z[j]
        y[j, k] <- rnorm(1, mu[j, k] * offset[j, k], ifelse(z[j] == 1, sqrt(tau.sq), 0))
      }
    } # k
  } # j
  if (family %in% c('Gaussian', 'zi-Gaussian')) {
    y <- y[, 1]
    mu <- mu[, 1]
    X <- X[, 1, ]
    if (length(mu.RE) > 0) {
      X.re <- X.re[, 1, ]
      beta.star.sites <- beta.star.sites[, 1]
    }
  }

  return(
    list(X = X, coords = coords, w = w.mat, mu = mu,
	 y = y, X.re = X.re, beta.star = beta.star)
  )
}
