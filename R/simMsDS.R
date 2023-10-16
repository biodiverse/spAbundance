simMsDS <- function(J.x, J.y, n.bins, bin.width, n.sp, beta, alpha, 
		    det.func, transect = 'line', kappa, 
		    mu.RE = list(), p.RE = list(), offset = 1, 
                    sp = FALSE, cov.model, sigma.sq, phi, 
		    nu, family = 'Poisson', factor.model = FALSE, 
		    n.factors, ...) {

  # Check for unused arguments ------------------------------------------
  formal.args <- names(formals(sys.function(sys.parent())))
  elip.args <- names(list(...))
  for(i in elip.args){
      if(! i %in% formal.args)
          warning("'",i, "' is not an argument")
  }

  # Check function inputs -------------------------------------------------
  J <- J.x * J.y
  # n.bins -----------------------------
  if (missing(n.bins)) {
    stop("error: n.bins must be specified.")
  }
  # bin.width --------------------------
  if (missing(bin.width)) {
    stop("error: bin.width must be specified.")
  }
  if (length(bin.width) != n.bins) {
    stop(paste("error: bin.width must be of length ", n.bins, ".", sep = ''))
  }
  # transect --------------------------
  if (transect != 'line' & transect != 'point') {
    stop("error: transect must be either 'point' or 'line'")
  }
  # family ------------------------------
  if (! (family %in% c('NB', 'Poisson'))) {
    stop("error: family must be either NB (negative binomial) or Poisson")
  }
  # n.sp ---------------------------------
  if (missing(n.sp)) {
    stop("error: n.sp must be specified")
  }
  if (length(n.sp) != 1) {
    stop("error: n.sp must be a single numeric value.")
  }
  # kappa -----------------------------
  if (family == 'NB') {
    if (missing(kappa)) {
      stop("error: kappa (overdispersion parameter) must be specified when family = 'NB'.")
    }
    if (length(kappa) != n.sp) {
      stop(paste("error: kappa must be a numeric vector with ", n.sp, " values", sep = ''))
    }
  }
  if (family == 'Poisson' & !missing(kappa)) {
    message("overdispersion parameter (kappa) is ignored when family == 'Poisson'")
  }
  # beta ------------------------------
  if (missing(beta)) {
    stop("error: beta must be specified")
  }
  if (!is.matrix(beta)) {
    stop(paste("error: beta must be a numeric matrix with ", n.sp, " rows", sep = ''))
  }
  if (nrow(beta) != n.sp) {
    stop(paste("error: beta must be a numeric matrix with ", n.sp, " rows", sep = ''))
  }
  # alpha -----------------------------
  if (missing(alpha)) {
    stop("error: alpha must be specified.")
  }
  if (!is.matrix(alpha)) {
    stop(paste("error: alpha must be a numeric matrix with ", n.sp, " rows", sep = ''))
  }
  if (nrow(alpha) != n.sp) {
    stop(paste("error: alpha must be a numeric matrix with ", n.sp, " rows", sep = ''))
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
  if (sp & !factor.model) {
    # sigma.sq --------------------------
    if (missing(sigma.sq)) {
      stop("error: sigma.sq must be specified when sp = TRUE")
    }
    if (length(sigma.sq) != n.sp) {
      stop(paste("error: sigma.sq must be a vector of length ", n.sp, sep = ''))
    }
    # phi -------------------------------
    if(missing(phi)) {
      stop("error: phi must be specified when sp = TRUE")
    }
    if (length(phi) != n.sp) {
      stop(paste("error: phi must be a vector of length ", n.sp, sep = ''))
    }
  }
  if (sp) {
    # Covariance model ----------------
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
  if (factor.model) {
    # n.factors -----------------------
    if (missing(n.factors)) {
      stop("error: n.factors must be specified when factor.model = TRUE")
    }
    if (sp) {
      if (!missing(sigma.sq)) {
        message("sigma.sq is specified but will be set to 1 for spatial latent factor model")
      }
      if(missing(phi)) {
        stop("error: phi must be specified when sp = TRUE")
      }
      if (length(phi) != n.factors) {
        stop(paste("error: phi must be a vector of length ", n.factors, sep = ''))
      }
    }
  }
  # det.func -------------------------
  if (missing(det.func)) {
    stop("error: det.func must be specified")
  }
  det.func.names <- c('halfnormal', 'negexp')
  if (! det.func %in% det.func.names) {
    stop("error: specified det.func '", det.func, "' is not a valid option; choose from ", 
	 paste(det.func.names, collapse = ', ', sep = ''), ".")
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
  # Half-normal detection function
  halfNormal <- function(x, sigma, transect) {
    if (transect == 'line') {
      exp(-x^2 / (2 * sigma^2))
    } else {
      exp(-x^2 / (2 * sigma^2)) * x
    }
  }
  # Negative exponential detection function
  negExp <- function(x, sigma, transect) {
    if (transect == 'line') {
      exp(-x / sigma)
    } else {
      exp(-x / sigma) * x
    }
  }

  # Form abundance covariates (if any) ------------------------------------
  p.abund <- ncol(beta)
  X <- matrix(1, nrow = J, ncol = p.abund)
  if (p.abund > 1) {
    for (i in 2:p.abund) {
      X[, i] <- rnorm(J)
    } # i
  }

  # Form detection covariate (if any) -------------------------------------
  p.det <- ncol(alpha)
  X.p <- matrix(1, nrow = J, ncol = p.det)
  if (p.det > 1) {
    for (i in 2:p.det) {
      X.p[, i] <- rnorm(J)
    } # i
  }

  # Simulate latent (spatial) random effect for each species --------------
  # Matrix of spatial locations
  s.x <- seq(0, 1, length.out = J.x)
  s.y <- seq(0, 1, length.out = J.y)
  coords <- as.matrix(expand.grid(s.x, s.y))
  w.star <- matrix(0, nrow = n.sp, ncol = J)
  if (factor.model) {
    lambda <- matrix(rnorm(n.sp * n.factors, 0, 0.5), n.sp, n.factors) 
    # Set diagonals to 1
    diag(lambda) <- 1
    # Set upper tri to 0
    lambda[upper.tri(lambda)] <- 0
    w <- matrix(NA, n.factors, J)
    if (sp) { # sfMsPGOcc
      if (cov.model == 'matern') {
        theta <- cbind(phi, nu)
      } else {
        theta <- as.matrix(phi)
      }
      for (ll in 1:n.factors) {
        Sigma <- mkSpCov(coords, as.matrix(1), as.matrix(0), 
            	     theta[ll, ], cov.model)
        w[ll, ] <- rmvn(1, rep(0, J), Sigma)
      }

    } else { # lsMsPGOcc
      for (ll in 1:n.factors) {
        w[ll, ] <- rnorm(J)
      } # ll  
    }
    for (j in 1:J) {
      w.star[, j] <- lambda %*% w[, j]
    }
  } else {
    if (sp) { # spMsPGOcc
      lambda <- NA
      if (cov.model == 'matern') {
        theta <- cbind(phi, nu)
      } else {
        theta <- as.matrix(phi)
      }
      # Spatial random effects for each species
      for (i in 1:n.sp) {
        Sigma <- mkSpCov(coords, as.matrix(sigma.sq[i]), as.matrix(0), 
            	     theta[i, ], cov.model)
        w.star[i, ] <- rmvn(1, rep(0, J), Sigma)
      }
    }
    # For naming consistency
    w <- w.star
    lambda <- NA
  }

  # Random effects --------------------------------------------------------
  # Abundance -------------------------
  if (length(mu.RE) > 0) {
    p.abund.re <- length(unlist(mu.RE$beta.indx))
    tmp <- sapply(mu.RE$beta.indx, length)
    re.col.indx <- unlist(lapply(1:length(mu.RE$beta.indx), function(a) rep(a, tmp[a])))
    sigma.sq.mu <- mu.RE$sigma.sq.mu[re.col.indx]
    n.abund.re.long <- mu.RE$levels[re.col.indx]
    n.abund.re <- sum(n.abund.re.long)
    beta.star.indx <- rep(1:p.abund.re, n.abund.re.long)
    beta.star <- matrix(0, n.sp, n.abund.re)
    X.random <- X[, unlist(mu.RE$beta.indx), drop = FALSE]
    n.random <- ncol(X.random)
    X.re <- matrix(NA, J, p.abund.re)
    for (l in 1:p.abund.re) {
      X.re[, l] <- sample(1:mu.RE$levels[l], J, replace = TRUE)         
      for (i in 1:n.sp) {
        beta.star[i, which(beta.star.indx == l)] <- rnorm(mu.RE$levels[l], 0, 
        						  sqrt(mu.RE$sigma.sq.mu[l]))
      }
    }
    indx.mat <- X.re[, re.col.indx, drop = FALSE]
    if (p.abund.re > 1) {
      for (j in 2:p.abund.re) {
        X.re[, j] <- X.re[, j] + max(X.re[, j - 1], na.rm = TRUE)
      }
    } 
    if (p.abund.re > 1) {
      for (j in 2:p.abund.re) {
        indx.mat[, j] <- indx.mat[, j] + max(indx.mat[, j - 1], na.rm = TRUE)
      }
    }
    beta.star.sites <- matrix(NA, n.sp, J)
    for (i in 1:n.sp) {
      for (j in 1:J) {
        beta.star.sites[i, j] <- beta.star[i, indx.mat[j, , drop = FALSE]] %*% t(X.random[j, , drop = FALSE])
      }
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
    alpha.star <- matrix(0, n.sp, n.det.re)
    X.p.random <- X.p[, unlist(p.RE$alpha.indx), drop = FALSE]
    n.p.random <- ncol(X.p.random)
    X.p.re <- matrix(NA, J, p.det.re)
    for (l in 1:p.det.re) {
      X.p.re[, l] <- sample(1:p.RE$levels[l], J, replace = TRUE)         
      for (i in 1:n.sp) {
        alpha.star[i, which(alpha.star.indx == l)] <- rnorm(p.RE$levels[l], 0, 
        						  sqrt(p.RE$sigma.sq.p[l]))
      }
    }
    indx.mat <- X.p.re[, p.re.col.indx, drop = FALSE]
    if (p.det.re > 1) {
      for (j in 2:p.det.re) {
        X.p.re[, j] <- X.p.re[, j] + max(X.p.re[, j - 1], na.rm = TRUE)
      }
    } 
    if (p.det.re > 1) {
      for (j in 2:p.det.re) {
        indx.mat[, j] <- indx.mat[, j] + max(indx.mat[, j - 1], na.rm = TRUE)
      }
    }
    alpha.star.sites <- matrix(NA, n.sp, J)
    for (i in 1:n.sp) {
      for (j in 1:J) {
        alpha.star.sites[i, j] <- alpha.star[i, indx.mat[j, , drop = FALSE]] %*% t(X.p.random[j, , drop = FALSE])
      }
    }
  } else {
    X.p.re <- NA
    alpha.star <- NA
  }

  # Latent Abundance Process ----------------------------------------------
  mu <- matrix(NA, nrow = n.sp, ncol = J)
  N <- matrix(NA, nrow = n.sp, ncol = J)
  if (family == 'NB') {
    for (i in 1:n.sp) {
      if (sp | factor.model) {
        if (length(mu.RE) > 0) {
          mu[i, ] <- exp(X %*% as.matrix(beta[i, ]) + w.star[i, ] + beta.star.sites[i, ])
        } else {
          mu[i, ] <- exp(X %*% as.matrix(beta[i, ]) + w.star[i, ])
        }
      } else {
        if (length(mu.RE) > 0) {
          mu[i, ] <- exp(X %*% as.matrix(beta[i, ]) + beta.star.sites[i, ])
        } else {
          mu[i, ] <- exp(X %*% as.matrix(beta[i, ]))
        }
      }
      N[i, ] <- rnbinom(J, size = kappa[i], mu = mu[i, ] * offset)
    }
  } else if (family == 'Poisson') {
    for (i in 1:n.sp) {
      if (sp | factor.model) {
        if (length(mu.RE) > 0) {
          mu[i, ] <- exp(X %*% as.matrix(beta[i, ]) + w.star[i, ] + beta.star.sites[i, ])
        } else {
          mu[i, ] <- exp(X %*% as.matrix(beta[i, ]) + w.star[i, ])
        }
      } else {
        if (length(mu.RE) > 0) {
          mu[i, ] <- exp(X %*% as.matrix(beta[i, ]) + beta.star.sites[i, ])
        } else {
          mu[i, ] <- exp(X %*% as.matrix(beta[i, ]))
        }
      }
      N[i, ] <- rpois(J, lambda = mu[i, ] * offset)
    }
  }

  # Data Formation --------------------------------------------------------
  sigma <- matrix(NA, nrow = n.sp, ncol = J)
  for (i in 1:n.sp) {
    if (length(p.RE) > 0) {
      sigma[i, ] <- exp(X.p %*% as.matrix(alpha[i, ]) + alpha.star.sites[i, ])
    } else {
      sigma[i, ] <- exp(X.p %*% as.matrix(alpha[i, ]))
    }
  }
  # Probability of detecting an individual in a given bin at a given site.
  p <- array(NA, dim = c(n.sp, J, n.bins))
  if (det.func == 'halfnormal') {
    curr.function <- halfNormal
  }
  if (det.func == 'negexp') {
    curr.function <- negExp
  }
  # Create distance bins
  dist.breaks <- rep(0, n.bins + 1)
  tmp <- 0
  for (i in 1:n.bins) {
    tmp <- tmp + bin.width[i]
    dist.breaks[i + 1] <- tmp
  }
  # This is B in the model notation. B = half line transect width or radius of point count
  strip.width <- sum(bin.width)
  for (i in 1:n.sp) {
    for (j in 1:J) {
      for (k in 1:n.bins) {
        p[i, j, k] <- integrate(curr.function, dist.breaks[k], 
          		   dist.breaks[k + 1], sig = sigma[i, j], transect = transect)$value
        if (transect == 'line') {
          p[i, j, k] <- p[i, j, k] / (dist.breaks[k + 1] - dist.breaks[k])
        } else {
          p[i, j, k] <- p[i, j, k] * 2 / (dist.breaks[k + 1]^2 - dist.breaks[k]^2)
        }
      }
    }
  }
  # Probability an individual occurs in each interval.
  psi <- rep(NA, n.bins)
  for (i in 1:n.bins) {
    if (transect == 'line') {
      psi[i] <- bin.width[i] / strip.width
    } else {
      psi[i] <- (dist.breaks[i + 1]^2 - dist.breaks[i]^2) / strip.width^2
    }
  }
  # Multinomial cell probability
  pi.obs <- apply(p,  c(1, 2), function(a) a * psi)
  pi.full <- abind(pi.obs, 1 - apply(pi.obs, c(2, 3), sum), along = 1)
  pi.full <- aperm(pi.full, c(2, 3, 1))

  y <- array(NA, dim = c(n.sp, J, n.bins + 1))
  for (i in 1:n.sp) {
    for (j in 1:J) {
      y[i, j, ] <- rmultinom(1, N[i, j], pi.full[i, j, ])
    }
  }

  # Data Formation --------------------------------------------------------

  return(
    list(X = X, X.p = X.p, coords = coords, w = w, lambda = lambda, mu = mu, N = N,
	 y = y[, , -(n.bins + 1), drop = FALSE], X.re = X.re, 
	 X.p.re = X.p.re, beta.star = beta.star, sigma = sigma,
	 pi.full = pi.full, alpha.star = alpha.star, dist.breaks = dist.breaks, 
         p = p)
  )
}
