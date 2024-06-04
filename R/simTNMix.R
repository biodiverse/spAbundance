simTNMix <- function(J.x, J.y, n.time, n.rep, n.rep.max, beta, alpha, sp.only = 0,
                     trend = TRUE, kappa, mu.RE = list(), p.RE = list(),
                     offset = 1, sp = FALSE, cov.model, sigma.sq, phi,
                     nu, family = 'Poisson', ar1 = FALSE, rho, sigma.sq.t, ...) {

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
  # n.time ---------------------------
  if (missing(n.time)) {
    stop("error: n.time must be specified.")
  }
  # n.rep -----------------------------
  if (missing(n.rep)) {
    stop("error: n.rep must be specified.")
  }
  if (!is.matrix(n.rep)) {
    stop(paste("error: n.rep must be a matrix with ", J, " rows and ", max(n.time), " columns", sep = ''))
  }
  if (missing(n.rep.max)) {
    n.rep.max <- max(n.rep, na.rm = TRUE)
  }
  if (nrow(n.rep) != J | ncol(n.rep) != max(n.time)) {
    stop(paste("error: n.rep must be a matrix with ", J, " rows and ", max(n.time), " columns", sep = ''))
  }
  # beta ------------------------------
  if (missing(beta)) {
    stop("error: beta must be specified.")
    if (length(beta) <= 1) {
      stop("error: beta must have at least two elements (intercept and trend)")
    }
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
  # AR1 -------------------------------
  if (ar1) {
    if (missing(rho)) {
      stop("error: rho must be specified when ar1 = TRUE")
    }
    if (missing(sigma.sq.t)) {
      stop("error: sigma.sq.t must be specified when ar1 = TRUE")
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

  # Form abundance covariates if any --------------------------------------
  n.beta <- length(beta)
  n.time.max <- max(n.time, na.rm = TRUE)
  time.indx <- list()
  for (j in 1:J) {
    time.indx[[j]] <- sample(which(!is.na(n.rep[j, ])), n.time[j], replace = FALSE)
  }
  X <- array(NA, dim = c(J, n.time.max, n.beta))
  X[, , 1] <- 1
  if (n.beta > 1) {
    if (trend) { # If simulating data with a trend
      # By default the second simulated covariate is a standardized trend
      X[, , 2] <- scale(c(matrix(rep(1:n.time.max, each = J), nrow = J, ncol = n.time.max)))
      if (n.beta > 2) {
        for (i in 3:n.beta) {
          if (i %in% sp.only) {
            X[, , i] <- rep(rnorm(J), n.time.max)
          } else {
            X[, , i] <- rnorm(J * n.time.max)
          }
        }
      }
    } else { # If not simulating data with a trend
      if (n.beta > 1) {
        for (i in 2:n.beta) {
          if (i %in% sp.only) {
            X[, , i] <- rep(rnorm(J), n.time.max)
          } else {
            X[, , i] <- rnorm(J * n.time.max)
          }
        }
      }
    }
  }
  # Form detection covariates (if any) ------------------------------------
  # Time dependent --------------------
  rep.indx <- list()
  for (j in 1:J) {
    rep.indx[[j]] <- list()
    for (t in time.indx[[j]]) {
      rep.indx[[j]][[t]] <- sample(1:n.rep.max, n.rep[j, t], replace = FALSE)
    }
  }
  n.alpha <- length(alpha)
  X.p <- array(NA, dim = c(J, n.time.max, n.rep.max, n.alpha))
  X.p[, , , 1] <- 1
  if (n.alpha > 1) {
    for (j in 1:J) {
      for (t in time.indx[[j]]) {
        for (k in rep.indx[[j]][[t]]) {
          X.p[j, t, k, 2:n.alpha] <- rnorm(n.alpha - 1)
        } # k
      } # t
    } # j
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
    X.random <- X[, , unlist(mu.RE$beta.indx), drop = FALSE]
    n.random <- dim(X.random)[3]
    X.re <- array(NA, dim = c(J, n.time.max, length(mu.RE$levels)))
    for (i in 1:length(mu.RE$levels)) {
      X.re[, , i] <- sample(1:mu.RE$levels[i], J * n.time.max, replace = TRUE)
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

    beta.star.sites <- matrix(NA, J, n.time.max)
    for (j in 1:J) {
      for (t in 1:n.time.max) {
        beta.star.sites[j, t] <- beta.star[indx.mat[j, t, ]] %*% as.matrix(X.random[j, t, ])
      } # k
    } # j
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
    X.p.random <- X.p[, , , unlist(p.RE$alpha.indx), drop = FALSE]
    X.p.re <- array(NA, dim = c(J, n.time.max, n.rep.max, length(p.RE$levels)))
    for (i in 1:length(p.RE$levels)) {
      X.p.re[, , , i] <- array(sample(1:p.RE$levels[i], J * n.rep.max * n.time.max, replace = TRUE),
                               dim = c(J, n.time.max, n.rep.max))
    }
    for (i in 1:p.det.re) {
      alpha.star[which(alpha.star.indx == i)] <- rnorm(n.det.re.long[i], 0, sqrt(sigma.sq.p[i]))
    }
    for (j in 1:J) {
      for (t in time.indx[[j]]) {
        X.p.re[j, t, -rep.indx[[j]][[t]], ] <- NA
      }
    }
    indx.mat <- X.p.re[, , , p.re.col.indx, drop = FALSE]
    if (length(p.RE$levels) > 1) {
      for (j in 2:length(p.RE$levels)) {
        X.p.re[, , , j] <- X.p.re[, , , j] + max(X.p.re[, , , j - 1], na.rm = TRUE)
      }
    }
    if (p.det.re > 1) {
      for (j in 2:p.det.re) {
        indx.mat[, , , j] <- indx.mat[, , , j] + max(indx.mat[, , , j - 1], na.rm = TRUE)
      }
    }
    alpha.star.sites <- array(NA, dim = c(J, n.time.max, n.rep.max))
    for (j in 1:J) {
      for (t in time.indx[[j]]) {
        for (k in rep.indx[[j]][[t]]) {
          alpha.star.sites[j, t, k] <- alpha.star[indx.mat[j, t, k, ]] %*% X.p.random[j, t, k, ]
        }
      }
    }
  } else {
    X.p.re <- NA
    alpha.star <- NA
  }

  # Simulate temporal (AR1) random effect ---------------------------------
  if (ar1) {
    exponent <- abs(matrix(1:n.time.max - 1, nrow = n.time.max,
                           ncol = n.time.max, byrow = TRUE) - (1:n.time.max - 1))
    Sigma.eta <- sigma.sq.t * rho^exponent
    eta <- rmvn(1, rep(0, n.time.max), Sigma.eta)
  } else {
    eta <- matrix(rep(0, n.time.max))
  }

  # Latent abundance process ----------------------------------------------
  mu <- matrix(NA, J, n.time.max)
  N <- matrix(NA, J, n.time.max)
  for (j in 1:J) {
    for (t in 1:n.time.max) {
      if (sp) {
        if (length(mu.RE) > 0) {
          mu[j, t] <- exp(X[j, t, ] %*% as.matrix(beta) + w[j] +
                          beta.star.sites[j, t] + eta[t])
        } else {
          mu[j, t] <- exp(X[j, t, ] %*% as.matrix(beta) + w[j] + eta[t])
        }
      } else {
        if (length(mu.RE) > 0) {
          mu[j, t] <- exp(X[j, t, ] %*% as.matrix(beta) +
                          beta.star.sites[j, t] + eta[t])
        } else {
          mu[j, t] <- exp(X[j, t, ] %*% as.matrix(beta) + eta[t])
        }
      }
      if (family == 'NB') {
        # Get mean and overdispersion parameter
        N[j, t] <- rnbinom(1, size = kappa, mu = mu[j, t] * offset)
      } else if (family == 'Poisson') {
        N[j, t] <- rpois(1, lambda = mu[j, t] * offset)
      }
    } # t
  } # j

  # Data Formation --------------------------------------------------------
  p <- array(NA, dim = c(J, n.time.max, n.rep.max))
  y <- array(NA, dim = c(J, n.time.max, n.rep.max))
  for (j in 1:J) {
    for (t in time.indx[[j]]) {
      if (length(p.RE) > 0) {
        p[j, t, rep.indx[[j]][[t]]] <- logit.inv(X.p[j, t, rep.indx[[j]][[t]], ]
                                                 %*% as.matrix(alpha) +
                                                 alpha.star.sites[j, t, rep.indx[[j]][[t]]])
      } else {
        p[j, t, rep.indx[[j]][[t]]] <- logit.inv(X.p[j, t, rep.indx[[j]][[t]], ]
                                                 %*% as.matrix(alpha))
      }
      y[j, t, rep.indx[[j]][[t]]] <- rbinom(n.rep[j, t], N[j, t], p[j, t, rep.indx[[j]][[t]]])
    } # t
  } # j

  # Return list -----------------------------------------------------------
  return(
    list(X = X, X.p = X.p, coords = coords, w = w, mu = mu, N = N,
         y = y, X.re = X.re, X.p.re = X.p.re, beta.star = beta.star,  p = p,
         alpha.star = alpha.star, eta = eta)
  )
}
