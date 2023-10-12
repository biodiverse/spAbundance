simMsGaussian <- function(J.x, J.y, N, beta, mu.RE = list(), 
		          sp = FALSE, svc.cols = 1, joint.cols = 1, cov.model, tau.sq,
		          sigma.sq, phi, nu, factor.model = FALSE, n.factors, 
		          z, ...) {

  # Check for unused arguments ------------------------------------------
  formal.args <- names(formals(sys.function(sys.parent())))
  elip.args <- names(list(...))
  for(i in elip.args){
      if(! i %in% formal.args)
          warning("'",i, "' is not an argument")
  }
  # Check function inputs -------------------------------------------------
  # J.x -------------------------------
  if (missing(J.x)) {
    stop("error: J.x must be specified")
  }
  if (length(J.x) != 1) {
    stop("error: J.x must be a single numeric value.")
  }
  # J.y -------------------------------
  if (missing(J.y)) {
    stop("error: J.y must be specified")
  }
  if (length(J.y) != 1) {
    stop("error: J.y must be a single numeric value.")
  }
  J <- J.x * J.y
  # N ---------------------------------
  if (missing(N)) {
    stop("error: N must be specified")
  }
  if (length(N) != 1) {
    stop("error: N must be a single numeric value.")
  }
  # beta ------------------------------
  if (missing(beta)) {
    stop("error: beta must be specified")
  }
  if (!is.matrix(beta)) {
    stop(paste("error: beta must be a numeric matrix with ", N, " rows", sep = ''))
  }
  if (nrow(beta) != N) {
    stop(paste("error: beta must be a numeric matrix with ", N, " rows", sep = ''))
  }
  # tau.sq ----------------------------
  if (missing(tau.sq)) {
    stop("error: tau.sq must be specified")
  }
  if (length(tau.sq) != N) {
    stop(paste("error: tau.sq must be a vector with ", N, " values", sep = ''))
  }
  # Check spatial stuff ---------------
  if (sp & !factor.model) {
    N.p.svc <- N * length(svc.cols)
    # sigma.sq --------------------------
    if (missing(sigma.sq)) {
      stop("error: sigma.sq must be specified when sp = TRUE")
    }
    if (length(sigma.sq) != N.p.svc) {
      stop(paste("error: sigma.sq must be a vector of length ", N.p.svc, sep = ''))
    }
    # phi -------------------------------
    if(missing(phi)) {
      stop("error: phi must be specified when sp = TRUE")
    }
    if (length(phi) != N.p.svc) {
      stop(paste("error: phi must be a vector of length ", N.p.svc, sep = ''))
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
    if ((length(n.factors) != 1) & (length(n.factors) != length(svc.cols))) {
        stop(paste("error: n.factors must be a vector of length 1 or ", 
		   length(svc.cols), sep = ''))
    }
    if (length(n.factors) == 1) {
      n.factors <- rep(n.factors, length(svc.cols))
    }
    q.p.svc <- sum(n.factors)
    if (sp) {
      if (!missing(sigma.sq)) {
        message("sigma.sq is specified but will be set to 1 for spatial latent factor model")
      }
      if(missing(phi)) {
        stop("error: phi must be specified when sp = TRUE")
      }
      if (length(phi) != q.p.svc) {
        stop(paste("error: phi must be a vector of length ", q.p.svc, sep = ''))
      }
    }
    if (!sp & length(svc.cols) > 1) {
      stop("error: svc.cols > 1 when sp = FALSE. Set sp = TRUE to simulate data with spatially-varying coefficients")
    }
  }
  # mu.RE ----------------------------
  names(mu.RE) <- tolower(names(mu.RE))
  if (!is.list(mu.RE)) {
    stop("error: if specified, mu.RE must be a list with tags 'levels' and 'sigma.sq.mu'")
  }
  if (length(names(mu.RE)) > 0) {
    if (!'sigma.sq.mu' %in% names(mu.RE)) {
      stop("error: sigma.sq.mu must be a tag in mu.RE with values for the occurrence random effect variances")
    }
    if (!'levels' %in% names(mu.RE)) {
      stop("error: levels must be a tag in mu.RE with the number of random effect levels for each occurrence random intercept.")
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

  # Form occupancy covariates (if any) ------------------------------------
  J <- J.x * J.y
  p.occ <- ncol(beta)
  X <- matrix(1, nrow = J, ncol = p.occ) 
  if (p.occ > 1) {
    for (i in 2:p.occ) {
      X[, i] <- rnorm(J)
    } # i
  }

  # Simulate latent (spatial) random effect for each species --------------
  # Matrix of spatial locations
  s.x <- seq(0, 1, length.out = J.x)
  s.y <- seq(0, 1, length.out = J.y)
  p.svc <- length(svc.cols)
  coords <- as.matrix(expand.grid(s.x, s.y))
  w.star <- vector(mode = "list", length = p.svc)
  w <- vector(mode = "list", length = p.svc)
  lambda <- vector(mode = "list", length = p.svc)
  theta.indx <- 1
  # Form spatial process for each spatially-varying covariate
  for (i in 1:p.svc) {
    w.star[[i]] <- matrix(0, nrow = N, ncol = J)
    if (factor.model) {
      # Note: can change this to generate from N(0, 1) instead. 
      if (i %in% joint.cols) {
        lambda[[i]] <- matrix(rnorm(N * n.factors[i], 0, 0.5), N, n.factors[i]) 
      } else {
        lambda[[i]] <- matrix(0, N, n.factors[i])
      }
      # Set diagonals to 1
      diag(lambda[[i]]) <- 1
      # Set upper tri to 0
      lambda[[i]][upper.tri(lambda[[i]])] <- 0
      w[[i]] <- matrix(NA, n.factors[i], J)
      if (sp) { # sfMsPGOcc
        for (ll in 1:n.factors[i]) {
          if (cov.model == 'matern') {
            theta <- cbind(phi[theta.indx], nu[theta.indx])
	  } else {
            theta <- as.matrix(phi[theta.indx])
	  }
          Sigma <- mkSpCov(coords, as.matrix(1), as.matrix(0), 
              	     theta, cov.model)
          w[[i]][ll, ] <- rmvn(1, rep(0, J), Sigma)
	  theta.indx <- theta.indx + 1
        }
      } else { # lsMsPGOcc
        for (ll in 1:n.factors[i]) {
          w[[i]][ll, ] <- rnorm(J)
        } # ll  
      }
      for (j in 1:J) {
        w.star[[i]][, j] <- lambda[[i]] %*% w[[i]][, j]
      }
    } else {
      if (sp) { # spMsPGOcc
        lambda <- NA
        if (cov.model == 'matern') {
          theta <- cbind(phi[((i - 1) * N + 1):(i * N)], 
			 nu[((i - 1) * N + 1):(i * N)])
        } else {
          theta <- as.matrix(phi[((i - 1) * N + 1):(i * N)])
        }
        # Spatial random effects for each species
        for (ll in 1:N) {
          Sigma <- mkSpCov(coords, as.matrix(sigma.sq[(i - 1) * N + ll]), as.matrix(0), 
              	     theta[ll, ], cov.model)
          w.star[[i]][ll, ] <- rmvn(1, rep(0, J), Sigma)
        }
      }
      # For naming consistency
      w <- w.star
      lambda <- NA
    }
  } # i (spatially-varying coefficient)
  # Design matrix for spatially-varying coefficients
  X.w <- X[, svc.cols, drop = FALSE]
  # Create X.tilde, which is a J x J*p.tilde matrix. 
  X.tilde <- matrix(0, J, J * p.svc)
  # Fill in the matrix
  for (j in 1:J) {
    X.tilde[j, ((j - 1) * p.svc + 1):(j * p.svc)] <- X.w[j, ]
  }

  # Random effects --------------------------------------------------------
  if (length(mu.RE) > 0) {
    p.occ.re <- length(mu.RE$levels)
    sigma.sq.mu <- rep(NA, p.occ.re)
    n.occ.re.long <- mu.RE$levels
    n.occ.re <- sum(n.occ.re.long)
    beta.star.indx <- rep(1:p.occ.re, n.occ.re.long)
    beta.star <- matrix(0, N, n.occ.re)
    X.re <- matrix(NA, J, p.occ.re)
    for (l in 1:p.occ.re) {
      X.re[, l] <- sample(1:mu.RE$levels[l], J, replace = TRUE)         
      for (i in 1:N) {
        beta.star[i, which(beta.star.indx == l)] <- rnorm(mu.RE$levels[l], 0, 
							  sqrt(mu.RE$sigma.sq.mu[l]))
      }
    }
    if (p.occ.re > 1) {
      for (j in 2:p.occ.re) {
        X.re[, j] <- X.re[, j] + max(X.re[, j - 1], na.rm = TRUE)
      }
    } 
    beta.star.sites <- matrix(NA, N, J)
    for (i in 1:N) {
      beta.star.sites[i, ] <- apply(X.re, 1, function(a) sum(beta.star[i, a]))
    }
  } else {
    X.re <- NA
    beta.star <- NA
  }

  # Latent Occupancy Process ----------------------------------------------
  mu <- matrix(NA, nrow = N, ncol = J)
  y <- matrix(NA, nrow = N, ncol = J)
  for (i in 1:N) {
    if (sp | factor.model) {
      w.star.curr.mat <- sapply(w.star, function(a) a[i, ])
      w.star.curr <- c(t(w.star.curr.mat))
      if (length(mu.RE) > 0) {
        mu[i, ] <- X %*% as.matrix(beta[i, ]) + 
                   X.tilde %*% w.star.curr + 
                   beta.star.sites[i, ]
      } else {
        mu[i, ] <- X %*% as.matrix(beta[i, ]) +
                   X.tilde %*% w.star.curr
      }
    } else {
      if (length(mu.RE) > 0) {
        mu[i, ] <- X %*% as.matrix(beta[i, ]) + beta.star.sites[i, ]
      } else {
        mu[i, ] <- X %*% as.matrix(beta[i, ])
      }
    }
    y[i, ] <- rnorm(J, mu[i, ], sqrt(tau.sq[i]))
    if (!missing(z)) {
      y[i, ] <- y[i, ] * z[i, ]
    }
  }

  return(
    list(X = X, coords = coords, w = w, mu = mu, y = y, 
	 X.re = X.re, beta.star = beta.star, lambda = lambda, X.w = X.w)
  )
}
