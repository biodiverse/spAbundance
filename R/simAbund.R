simAbund <- function(J.x, J.y, n.rep, beta, kappa, mu.RE = list(),  
		     sp = FALSE, cov.model, sigma.sq, phi, nu, family = 'NB', ...) {

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
  # beta ------------------------------
  if (missing(beta)) {
    stop("error: beta must be specified.")
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
  X <- array(NA, dim = c(J, max(n.rep), n.beta))
  X[, , 1] <- 1
  if (n.beta > 1) {
    for (i in 2:n.beta) {
      for (j in 1:J) {
        X[j, 1:n.rep[j], i] <- rnorm(n.rep[j])
      } 
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
    X.re <- array(NA, dim = c(J, max(n.rep), length(unique(re.col.indx))))
    for (i in 1:length(unique(re.col.indx))) {
      for (j in 1:J) {
        X.re[j, 1:n.rep[j], i] <- sample(1:mu.RE$levels[i], n.rep[j], replace = TRUE)  
      }
    }
    indx.mat <- X.re[, , re.col.indx, drop = FALSE]
    for (i in 1:p.nmix.re) {
      beta.star[which(beta.star.indx == i)] <- rnorm(n.nmix.re.long[i], 0, 
						     sqrt(sigma.sq.mu[i]))
    }
    if (p.nmix.re > 1) {
      for (j in 2:p.nmix.re) {
        indx.mat[, , j] <- indx.mat[, , j] + max(indx.mat[, , j - 1], na.rm = TRUE)
      }
    }
    beta.star.sites <- matrix(NA, J, max(n.rep))
    for (j in 1:J) {
      for (k in 1:n.rep[j]) {  
        beta.star.sites[j, k] <- beta.star[indx.mat[j, k, ]] %*% as.matrix(X.random[j, k, ])
      } # k
    } # j
  } else {
    X.re <- NA
    beta.star <- NA
  }

  # Data formation --------------------------------------------------------
  mu <- matrix(NA, J, max(n.rep))
  y <- matrix(NA, J, max(n.rep))
  for (j in 1:J) {
    for (k in 1:n.rep[j]) {
      if (sp) {
        if (length(mu.RE) > 0) {
          mu[j, k] <- exp(t(as.matrix(X[j, k, ])) %*% as.matrix(beta) + 
      			   w[j] + 
      			   beta.star.sites[j, k])
        } else {
          mu[j, k] <- exp(t(as.matrix(X[j, k, ])) %*% as.matrix(beta) + w[j])
        }
      } else {
        if (length(mu.RE) > 0) {
          mu[j, k] <- exp(t(as.matrix(X[j, k, ])) %*% as.matrix(beta) + 
      			   beta.star.sites[j, k])
        } else {
          mu[j, k] <- exp(t(as.matrix(X[j, k, ])) %*% as.matrix(beta))
        }
      }
      if (family == 'NB') {
        y[j, k] <- rnbinom(1, size = kappa, mu = mu[j, k])
      } else {
        y[j, k] <- rpois(1, lambda = mu[j, k])
      }
    } # k
  } # j

  return(
    list(X = X, coords = coords, w = w, mu = mu, 
	 y = y, X.re = X.re, beta.star = beta.star)
  )
}
