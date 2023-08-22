simMsAbund <- function(J.x, J.y, n.rep, n.sp, beta, kappa, tau.sq, mu.RE = list(), 
		       sp = FALSE, cov.model, sigma.sq, phi, nu, family = 'Poisson',
		       factor.model = FALSE, n.factors, z, ...) {

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
  # n.rep -----------------------------
  if (missing(n.rep)) {
    stop("error: n.rep must be specified.")
  }
  if (length(n.rep) != J) {
    stop(paste("error: n.rep must be a vector of length ", J, sep = ''))
  }
  # family ------------------------------
  if (! (family %in% c('NB', 'Poisson', 'Gaussian', 'zi-Gaussian'))) {
    stop("error: family must be one of: NB (negative binomial), Poisson, 'Gaussian', or 'zi-Gaussian'")
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
  # tau.sq ----------------------------
  if (family %in% c('Gaussian', 'zi-Gaussian')) {
    if (missing(tau.sq)) {
      stop('error: tau.sq (residual variance) must be specified when family is Gaussian or zi-Gaussian')
    }
    if (length(tau.sq) != n.sp) {
      stop(paste("error: tau.sq must be a numeric vector with ", n.sp, " values", sep = ''))
    }
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

  # z values --------------------------
  if (family == 'zi-Gaussian') {
    if (missing(z)) {
      stop('for a zero-inflated Gaussian model, you must supply the z values (binary 0s or 1s)')
    }
    if (!is.matrix(z)) {
      stop(paste0("z must be a matrix with ", n.sp, " rows and ", J.x * J.y, " columns."))
    }
    if (nrow(z) != n.sp | ncol(z) != J.x * J.y) {
      stop(paste0("z must be a matrix with ", n.sp, " rows and ", J.x * J.y, " columns."))
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
  J <- J.x * J.y
  p.abund <- ncol(beta)
  X <- array(NA, dim = c(J, max(n.rep), p.abund))
  X[, , 1] <- 1
  if (p.abund > 1) {
    for (i in 2:p.abund) {
      for (j in 1:J) {
        X[j, 1:n.rep[j], i] <- rnorm(n.rep[j])
      } 
    } # i
  }

  # Simulate latent (spatial) random effect for each species --------------
  # Matrix of spatial locations
  s.x <- seq(0, 1, length.out = J.x)
  s.y <- seq(0, 1, length.out = J.y)
  coords <- as.matrix(expand.grid(s.x, s.y))
  w.star <- matrix(0, nrow = n.sp, ncol = J)
  if (factor.model) {
    lambda <- matrix(rnorm(n.sp * n.factors), n.sp, n.factors) 
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
  if (length(mu.RE) > 0) {
    p.abund.re <- length(unlist(mu.RE$beta.indx))
    tmp <- sapply(mu.RE$beta.indx, length)
    re.col.indx <- unlist(lapply(1:length(mu.RE$beta.indx), function(a) rep(a, tmp[a])))
    sigma.sq.mu <- mu.RE$sigma.sq.mu[re.col.indx]
    n.abund.re.long <- mu.RE$levels[re.col.indx]
    n.abund.re <- sum(n.abund.re.long)
    beta.star.indx <- rep(1:p.abund.re, n.abund.re.long)
    beta.star <- matrix(0, n.sp, n.abund.re)
    X.random <- X[, , unlist(mu.RE$beta.indx), drop = FALSE]
    X.re <- array(NA, dim = c(J, max(n.rep), p.abund.re))
    for (l in 1:p.abund.re) {
      X.re[, , l] <- matrix(sample(1:mu.RE$levels[l], J * max(n.rep), replace = TRUE), 
		              J, max(n.rep))	      
      for (i in 1:n.sp) {
        beta.star[i, which(beta.star.indx == l)] <- rnorm(mu.RE$levels[l], 0, sqrt(mu.RE$sigma.sq[l]))
      }
    }
    for (j in 1:J) {
      X.re[j, -(1:n.rep[j]), ] <- NA
    }
    indx.mat <- X.re[, , re.col.indx, drop = FALSE]
    if (p.abund.re > 1) {
      for (j in 2:p.abund.re) {
        X.re[, , j] <- X.re[, , j] + max(X.re[, , j - 1], na.rm = TRUE) 
      }
    }
    if (p.abund.re > 1) {
      for (j in 2:p.abund.re) {
        indx.mat[, , j] <- indx.mat[, , j] + max(indx.mat[, , j - 1], na.rm = TRUE) 
      }
    }
    beta.star.sites <- array(NA, c(n.sp, J, max(n.rep)))
    for (i in 1:n.sp) {
      for (j in 1:J) {
        for (k in 1:n.rep[j]) {
          beta.star.sites[i, j, k] <- beta.star[i, indx.mat[j, k, ]] %*% X.random[j, k,] 
	}
      }
    }
  } else {
    X.re <- NA
    beta.star <- NA
  }

  # Data formation --------------------------------------------------------
  mu <- array(NA, dim = c(n.sp, J, max(n.rep)))
  y <- array(NA, dim = c(n.sp, J, max(n.rep)))
  for (j in 1:J) {
    for (k in 1:n.rep[j]) {
      for (i in 1:n.sp) {
        if (sp | factor.model) {
          if (length(mu.RE) > 0) {
            mu[i, j, k] <- t(as.matrix(X[j, k, ])) %*% as.matrix(beta[i, ]) + w.star[i, j] + beta.star.sites[i, j, k]
          } else {
            mu[i, j, k] <- t(as.matrix(X[j, k, ])) %*% as.matrix(beta[i, ]) + w.star[i, j]
          }
        } else {
          if (length(mu.RE) > 0) {
            mu[i, j, k] <- t(as.matrix(X[j, k, ])) %*% as.matrix(beta[i, ]) + beta.star.sites[i, j, k]
          } else {
            mu[i, j, k] <- t(as.matrix(X[j, k, ])) %*% as.matrix(beta[i, ])
          }
        }
        if (family %in% c('Poisson', 'NB')) {
          mu[i, j, k] <- exp(mu[i, j, k])
	}
        if (family == 'NB') {
          y[i, j, k] <- rnbinom(1, size = kappa[i], mu = mu[i, j, k])
        }
	if (family == 'Poisson') {
          y[i, j, k] <- rpois(1, lambda = mu[i, j, k])
        }
        if (family == 'Gaussian') {
          y[i, j, k] <- rnorm(1, mu[i, j, k], sqrt(tau.sq[i]))
        }
        if (family == 'zi-Gaussian') {
          mu[i, j, k] <- mu[i, j, k] * z[i, j]
          y[i, j, k] <- rnorm(1, mu[i, j, k], ifelse(z[i, j] == 1, sqrt(tau.sq[i]), 0))
        }
      } # i (species)
    } # k (replicate)
  } # j (site)

  if (family %in% c('Gaussian', 'zi-Gaussian')) {
    y <- y[, , 1]
    mu <- mu[, , 1]
    X <- X[, 1, ]
    if (length(mu.RE) > 0) {
      X.re <- X.re[, 1, ]
      beta.star.sites <- beta.star.sites[, , 1]
    }
  }

  return(
    list(X = X, coords = coords, w = w, lambda = lambda, y = y, 
	 X.re = X.re, beta.star = beta.star, mu = mu)
  )
}
