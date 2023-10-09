simNMix <- function(J.x, J.y, n.rep, n.rep.max, beta, alpha, 
		    kappa, mu.RE = list(), p.RE = list(), 
		    offset = 1, sp = FALSE, cov.model, sigma.sq, 
		    phi, nu, family = 'Poisson', ...) {

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
  # alpha -----------------------------
  if (missing(alpha)) {
    stop("error: alpha must be specified.")
  }
  # family ------------------------------
  if (! (family %in% c('NB', 'Poisson'))) {
    stop("error: family must be either NB (negative binomial) or Poisson")
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
  # p.RE ----------------------------
  names(p.RE) <- tolower(names(p.RE))
  if (!is.list(p.RE)) {
    stop("error: if specified, p.RE must be a list with tags 'levels' and 'sigma.sq.p'")
  }
  if (length(names(p.RE)) > 0) {
    if (!'sigma.sq.p' %in% names(p.RE)) {
      stop("error: sigma.sq.p must be a tag in p.RE with values for the detection random effect variances")
    }
    if (!'levels' %in% names(p.RE)) {
      stop("error: levels must be a tag in p.RE with the number of random effect levels for each detection random intercept.")
    }
    if (!'alpha.indx' %in% names(p.RE)) {
      p.RE$alpha.indx <- list()
      for (i in 1:length(p.RE$sigma.sq.p)) {
        p.RE$alpha.indx[[i]] <- 1
      }
    }
  }
  # Spatial parameters ----------------
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
  X <- matrix(1, nrow = J, ncol = n.beta)
  if (n.beta > 1) {
    for (i in 2:n.beta) {
      X[, i] <- rnorm(J)
    } # i
  }

  # Form detection covariate (if any) -------------------------------------
  n.alpha <- length(alpha)
  X.p <- array(NA, dim = c(J, n.rep.max, n.alpha))
  X.p[, , 1] <- 1
  # Get index of surveyed replicates for each site. 
  rep.indx <- list()
  for (j in 1:J) {
    rep.indx[[j]] <- sample(1:n.rep.max, n.rep[j], replace = FALSE)
  }
  if (n.alpha > 1) {
    for (i in 2:n.alpha) {
      for (j in 1:J) {
        X.p[j, rep.indx[[j]], i] <- rnorm(n.rep[j])
      } # j
    } # i
  }

  # Simulate spatial random effect ----------------------------------------
  # Matrix of spatial locations
  s.x <- seq(0, 1, length.out = J.x)
  s.y <- seq(0, 1, length.out = J.y)
  coords <- as.matrix(expand.grid(s.x, s.y))
  if (sp) {
    if (cov.model == 'matern') {
      theta <- c(phi, nu)
    } else {
      theta <- phi
    }
    Sigma <- mkSpCov(coords, as.matrix(sigma.sq), as.matrix(0), theta, cov.model)
    # Random spatial process
    w <- rmvn(1, rep(0, J), Sigma)
  } else {
    w <- NA
  }

  # Random effects --------------------------------------------------------
  # Abundance -------------------------
  if (length(mu.RE) > 0) {
    p.nmix.re <- length(unlist(mu.RE$beta.indx))
    tmp <- sapply(mu.RE$beta.indx, length)
    re.col.indx <- unlist(lapply(1:length(mu.RE$beta.indx), function(a) rep(a, tmp[a])))
    sigma.sq.mu <- mu.RE$sigma.sq.mu[re.col.indx]
    n.nmix.re.long <- mu.RE$levels[re.col.indx]
    n.nmix.re <- sum(n.nmix.re.long)
    beta.star.indx <- rep(1:p.nmix.re, n.nmix.re.long)
    beta.star <- rep(0, n.nmix.re)
    X.random <- X[, unlist(mu.RE$beta.indx), drop = FALSE]
    n.random <- ncol(X.random)
    X.re <- matrix(NA, J, length(mu.RE$levels))
    for (i in 1:length(mu.RE$levels)) {
      X.re[, i] <- sample(1:mu.RE$levels[i], J, replace = TRUE)  
    }
    indx.mat <- X.re[, re.col.indx, drop = FALSE]
    for (i in 1:p.nmix.re) {
      beta.star[which(beta.star.indx == i)] <- rnorm(n.nmix.re.long[i], 0, 
						     sqrt(sigma.sq.mu[i]))
    }
    if (length(mu.RE$levels) > 1) {
      for (j in 2:length(mu.RE$levels)) {
        X.re[, j] <- X.re[, j] + max(X.re[, j - 1], na.rm = TRUE)
      }
    }
    if (p.nmix.re > 1) {
      for (j in 2:p.nmix.re) {
        indx.mat[, j] <- indx.mat[, j] + max(indx.mat[, j - 1], na.rm = TRUE)
      }
    }
    beta.star.sites <- rep(NA, J)
    for (j in 1:J) {
      beta.star.sites[j] <- beta.star[indx.mat[j, , drop = FALSE]] %*% t(X.random[j, , drop = FALSE])
    }
  } else {
    X.re <- NA
    beta.star <- NA
  }
  # Detection -------------------------
  if (length(p.RE) > 0) {
    p.det.re <- length(unlist(p.RE$alpha.indx))
    tmp <- sapply(p.RE$alpha.indx, length)
    p.re.col.indx <- unlist(lapply(1:length(p.RE$alpha.indx), function(a) rep(a, tmp[a])))
    sigma.sq.p <- p.RE$sigma.sq.p[p.re.col.indx]
    n.det.re.long <- p.RE$levels[p.re.col.indx]
    n.det.re <- sum(n.det.re.long)
    alpha.star.indx <- rep(1:p.det.re, n.det.re.long)
    alpha.star <- rep(0, n.det.re)
    X.p.random <- X.p[, , unlist(p.RE$alpha.indx), drop = FALSE]
    X.p.re <- array(NA, dim = c(J, n.rep.max, length(p.RE$levels)))
    for (i in 1:length(p.RE$levels)) {
      X.p.re[, , i] <- matrix(sample(1:p.RE$levels[i], J * n.rep.max, replace = TRUE), 
		              J, n.rep.max)	      
    }
    for (i in 1:p.det.re) {
      alpha.star[which(alpha.star.indx == i)] <- rnorm(n.det.re.long[i], 0, sqrt(sigma.sq.p[i]))
    }
    for (j in 1:J) {
      X.p.re[j, -rep.indx[[j]], ] <- NA
    }
    indx.mat <- X.p.re[, , p.re.col.indx, drop = FALSE]
    if (length(p.RE$levels) > 1) {
      for (j in 2:length(p.RE$levels)) {
        X.p.re[, , j] <- X.p.re[, , j] + max(X.p.re[, , j - 1], na.rm = TRUE)
      }
    }
    if (p.det.re > 1) {
      for (j in 2:p.det.re) {
        indx.mat[, , j] <- indx.mat[, , j] + max(indx.mat[, , j - 1], na.rm = TRUE) 
      }
    }
    alpha.star.sites <- matrix(NA, J, n.rep.max)
    for (j in 1:J) {
      for (k in rep.indx[[j]]) {
        alpha.star.sites[j, k] <- alpha.star[indx.mat[j, k, ]] %*% X.p.random[j, k, ]
      }
    }
  } else {
    X.p.re <- NA
    alpha.star <- NA
  }

  # Latent abundance process ----------------------------------------------
  if (family == 'NB') {
    if (sp) {
      if (length(mu.RE) > 0) {
        mu <- exp(X %*% as.matrix(beta) + w + beta.star.sites)
      } else {
        mu <- exp(X %*% as.matrix(beta) + w)
      }
    } else {
      if (length(mu.RE) > 0) {
        mu <- exp(X %*% as.matrix(beta) + beta.star.sites)
      } else {
        mu <- exp(X %*% as.matrix(beta))
      }
    }
    # Get mean and overdispersion parameter 
    N <- rnbinom(J, size = kappa, mu = mu * offset)
  } else if (family == 'Poisson') {
    if (sp) {
      if (length(mu.RE) > 0) {
        mu <- exp(X %*% as.matrix(beta) + w + beta.star.sites)
      } else {
        mu <- exp(X %*% as.matrix(beta) + w)
      }
    } else {
      if (length(mu.RE) > 0) {
        mu <- exp(X %*% as.matrix(beta) + beta.star.sites)
      } else {
        mu <- exp(X %*% as.matrix(beta))
      }
    }
    N <- rpois(J, lambda = mu * offset)
  }

  # Data Formation --------------------------------------------------------
  p <- matrix(NA, nrow = J, ncol = n.rep.max)
  y <- matrix(NA, nrow = J, ncol = n.rep.max)
  for (j in 1:J) {
    if (length(p.RE) > 0) {
      p[j, rep.indx[[j]]] <- logit.inv(X.p[j, rep.indx[[j]], ] %*% as.matrix(alpha) + 
				    alpha.star.sites[j, rep.indx[[j]]])
    } else {
      p[j, rep.indx[[j]]] <- logit.inv(X.p[j, rep.indx[[j]], ] %*% as.matrix(alpha))
    }
    y[j, rep.indx[[j]]] <- rbinom(n.rep[j], N[j], p[j, rep.indx[[j]]])
  } # j


  # Data Formation --------------------------------------------------------

  return(
    list(X = X, X.p = X.p, coords = coords, w = w, mu = mu, N = N,
	 y = y, X.re = X.re, X.p.re = X.p.re, beta.star = beta.star,  p = p,
         alpha.star = alpha.star)
  )
}
