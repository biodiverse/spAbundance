simTIntAbund <- function(n.data, J.x, J.y, J.obs, n.time, data.seasons, n.rep, n.rep.max, 
                       beta, alpha, kappa, sp.only = 0, trend = TRUE, mu.RE = list(), 
                       p.RE = list(), offset = NULL, sp = FALSE, svc.cols = 1, cov.model, 
                       sigma.sq, phi, nu, family = 'Poisson', ar1 = FALSE, rho, sigma.sq.t, 
                       x.positive = FALSE, det.intercepts = FALSE, ...) {
  
  # Check for unused arguments ------------------------------------------
  formal.args <- names(formals(sys.function(sys.parent())))
  elip.args <- names(list(...))
  for(i in elip.args){
      if(! i %in% formal.args)
          warning("'",i, "' is not an argument")
  }
  # Subroutines -----------------------------------------------------------
  rmvn <- function(n, mu=0, V = matrix(1)){
    p <- length(mu)
    if(any(is.na(match(dim(V),p)))){stop("Dimension problem!")}
    D <- chol(V)
    t(matrix(rnorm(n*p), ncol=p)%*%D + rep(mu,rep(n,p)))
  }
  
  # n.data -------------------------------
  if (missing(n.data)) {
    stop("error: n.data must be specified")
  }
  if (length(n.data) != 1) {
    stop("error: n.data must be a single numeric value.")
  }
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
  # J.obs -----------------------------
  if (missing(J.obs)) {
    stop("error: J.obs must be specified")
  }
  if (length(J.obs) != n.data) {
    stop(paste("error: J.obs must be a vector of length ", n.data, sep = ''))
  }
  # n.time ---------------------------
  if (missing(n.time)) {
    stop("error: n.time must be specified.")
  }
  # n.rep -----------------------------
  if (missing(n.rep)) {
    stop("error: n.rep must be specified.")
  }
  if (!is.list(n.rep)) {
    stop(paste("error: n.rep must be a list of ", n.data, " vectors", sep = ''))
  }
  if (length(n.rep) != n.data) {
    stop(paste("error: n.rep must be a list of ", n.data, " vectors", sep = ''))
  }
  for (i in 1:n.data) {
    if (!is.matrix(n.rep[[i]])) {
      stop(paste("error: n.rep must be a matrix with ", J.obs[i], " rows and ", 
                 max(n.time[[i]]), " columns", sep = ''))
    }
    if (nrow(n.rep[[i]]) != J.obs[i] | ncol(n.rep[[i]]) != max(n.time[[i]])) {
      stop(paste("error: n.rep must be a matrix with ", J.obs[i], " rows and ", 
                 max(n.time[[i]]), " columns", sep = ''))
    }
  }
  if (missing(n.rep.max)) {
    n.rep.max <- sapply(n.rep, max, na.rm = TRUE)
  }
   # family ---------------------------
   if (length(family) == 1) {
     family <- rep(family, n.data)
   }
  for (i in 1:n.data) {
    if (! (family[i] %in% c('NB', 'Poisson'))) {
      stop('family must be one of: NB (negative binomial) or Poisson')
    }
  }
  # kappa -----------------------------
  if (any(family == 'NB')) {
    if (missing(kappa)) {
      stop("kappa (overdispersion parameter) must be specified when family = 'NB'.")
    }
    if (length(kappa) == 1) {
      kappa <- rep(kappa, n.data)
    }
  }

  # beta ------------------------------
  if (missing(beta)) {
    stop("error: beta must be specified.")
  }
  # alpha -----------------------------
  if (missing(alpha)) {
    stop("error: alpha must be specified.")
  }
  if (!is.list(alpha)) {
    stop(paste("error: alpha must be a list with ", n.data, " vectors", sep = ''))
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
  # p.RE ----------------------------
  if (!is.list(p.RE)) {
    stop(paste("error: if species, p.RE must be a list with ", n.data, " lists", sep = ''))
  }
  if (length(p.RE) > 0) {
    for (q in 1:n.data) {
      names(p.RE[[q]]) <- tolower(names(p.RE[[q]]))
      if (!is.list(p.RE[[q]])) {
        stop("error: if specified, p.RE[[", q, "]] must be a list with tags 'levels' and 'sigma.sq.p'")
      }
      if (length(names(p.RE[[q]])) > 0) {
        if (!'sigma.sq.p' %in% names(p.RE[[q]])) {
          stop("error: sigma.sq.p must be a tag in p.RE[[", q, "]] with values for the detection random effect variances")
        }
        if (!'levels' %in% names(p.RE[[q]])) {
          stop("error: levels must be a tag in p.RE[[", q, "]] with the number of random effect levels for each detection random intercept.")
        }
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
  # AR1 -------------------------------
  if (ar1) {
    if (missing(rho)) {
      stop("error: rho must be specified when ar1 = TRUE")
    }
    if (missing(sigma.sq.t)) {
      stop("error: sigma.sq.t must be specified when ar1 = TRUE")
    }
  }
  # data.seasons ------------------------
  if (missing(data.seasons)) {
    stop("error: data.seasons must be specified")
    if (length(data.seasons) != n.data) {
      stop(paste0("data.seasons must be a list with ", n.data, " vectors"))
    }
  }
  n.time.total <- length(unique(unlist(data.seasons)))
  n.time.max <- sapply(data.seasons, length)

  # offset ----------------------------
  if (is.null(offset)) {
    offset <- list()
    for (i in 1:n.data) {
      offset[i] <- 1
    }
  } else if (!is.list(offset)) {
    stop(paste0("offset must be a list with ", n.data, " components."))
  }
  
  # Matrix of spatial locations -------------------------------------------
  s.x <- seq(0, 1, length.out = J.x)
  s.y <- seq(0, 1, length.out = J.y)
  coords <- as.matrix(expand.grid(s.x, s.y))
  
  # Get site ids for each of the data sets --------------------------------
  # Data sources can be obtained at multiple different sites. 
  sites <- list()
  for (i in 1:n.data) {
    sites[[i]] <- sort(sample(1:J, J.obs[i], replace = FALSE))
  }
  
  # Abundurrence ------------------------------------------------------------
  p.abund <- length(beta)
  # A list of list, with each list corresponding a data set, and then the list
  # underneath giving the corresponding years that each site is sampled.
  time.indx <- matrix(0, J, n.time.total) 
  time.indx <- list()
  for (i in 1:n.data) {
    time.indx[[i]] <- list()
    for (j in 1:J.obs[i]) {
      curr.time <- sample(1:n.time.max[i], n.time[[i]][j], replace = FALSE) 
      time.indx[[i]][[j]] <- which(!is.na(n.rep[[i]][j, ]))
    }
  }
  X <- array(NA, dim = c(J, n.time.total, p.abund))
  X[, , 1] <- 1
  if (p.abund > 1) {
    if (trend) { # If simulating data with a trend
      # By default the second simulated covariate is a standardized trend
      X[, , 2] <- scale(c(matrix(rep(1:n.time.total, each = J), nrow = J, ncol = n.time.total)))
      if (p.abund > 2) {
        for (i in 3:p.abund) {
          if (i %in% sp.only) {
            X[, , i] <- rep(rnorm(J), n.time.total)
          } else {
            X[, , i] <- rnorm(J * n.time.total)
          }
        }
      }
    } else { # If not simulating data with a trend
      if (p.abund > 1) {
        for (i in 2:p.abund) {
          if (i %in% sp.only) {
            X[, , i] <- rep(rnorm(J), n.time.total)
          } else {
            X[, , i] <- rnorm(J * n.time.total)
          }
        }
      }
    }
  }
  if (x.positive) {
    if (p.abund > 1) {
      for (i in 2:p.abund) {
        X[, , i] <- runif(J * n.time.total, 0, 1)
      }
    }
  }
  # Form detection covariates (if any) ------------------------------------
  X.p <- list()
  rep.indx <- list()
  for (i in 1:n.data) {
    rep.indx[[i]] <- list()
    n.alpha.curr <- length(alpha[[i]])
    K.curr <- n.rep[[i]]
    J.curr <- J.obs[[i]]
    for (j in 1:J.curr) {
      rep.indx[[i]][[j]] <- list()
      for (t in time.indx[[i]][[j]]) {
        rep.indx[[i]][[j]][[t]] <- sample(1:n.rep.max[i], n.rep[[i]][j, t], 
                                          replace = FALSE) 
      }
    }
    X.p[[i]] <- array(NA, dim = c(J.curr, n.time.max[i], n.rep.max[i], n.alpha.curr))
    if (det.intercepts) {
      X.p[[i]][, , , 1] <- 1
      if (n.alpha.curr > 1) {
        for (q in 2:n.alpha.curr) {
          for (j in 1:J.curr) {
            for (t in time.indx[[i]][[j]]) {
              for (k in rep.indx[[i]][[j]][[t]]) {
                X.p[[i]][j, t, k, 2:n.alpha.curr] <- rnorm(n.alpha.curr - 1)
              } # k
            } # t
          } # j
        } # q
      }
    } else {
      for (q in 1:n.alpha.curr) {
        for (j in 1:J.curr) {
          for (t in time.indx[[i]][[j]]) {
            for (k in rep.indx[[i]][[j]][[t]]) {
              X.p[[i]][j, t, k, 1:n.alpha.curr] <- rnorm(n.alpha.curr)
            } # k
          } # t
        } # j
      } # q
    }
  } # i
  
  # Random effects --------------------------------------------------------
  # Abundance -------------------------
  if (length(mu.RE) > 0) {
    p.abund.re <- length(mu.RE$levels)
    sigma.sq.mu <- rep(NA, p.abund.re)
    n.abund.re.long <- mu.RE$levels
    n.abund.re <- sum(n.abund.re.long)
    beta.star.indx <- rep(1:p.abund.re, n.abund.re.long)
    beta.star <- rep(0, n.abund.re)
    X.re <- array(NA, dim = c(J, n.time.total, p.abund.re))
    for (i in 1:p.abund.re) {
      if (length(mu.RE$site.re) == 0) mu.RE$site.re <- FALSE
      if (mu.RE$site.re == TRUE) {
        if (i == 1) {
          site.vals <- 1:J 
          X.re[, , i] <- site.vals
        } else {
          X.re[, , i] <- sample(1:mu.RE$levels[i], J * n.time.total, replace = TRUE)         
        }
      } else {
        X.re[, , i] <- sample(1:mu.RE$levels[i], J * n.time.total, replace = TRUE)         
      }
      beta.star[which(beta.star.indx == i)] <- rnorm(mu.RE$levels[i], 0, sqrt(mu.RE$sigma.sq.mu[i]))
    }
    if (p.abund.re > 1) {
      for (j in 2:p.abund.re) {
        X.re[, , j] <- X.re[, , j] + max(X.re[, , j - 1])
      }
    } 
    beta.star.sites <- apply(X.re, c(1, 2), function(a) sum(beta.star[a]))
  } else {
    X.re <- NA
    beta.star <- NA
  }
  # Detection -------------------------
  if (length(p.RE) > 0) {
    p.det.re <- list()
    n.det.re.long <- list()
    n.det.re <- list()
    alpha.star.indx <- list()
    alpha.star <- list()
    alpha.star.sites <- list()
    X.p.re <- list()
    for (q in 1:n.data) {
      if (length(p.RE[[q]]) > 0) {
        p.det.re[[q]] <- length(p.RE[[q]]$levels)
        n.det.re.long[[q]] <- p.RE[[q]]$levels
        n.det.re[[q]] <- sum(n.det.re.long[[q]])
        alpha.star.indx[[q]] <- rep(1:p.det.re[[q]], n.det.re.long[[q]])
        alpha.star[[q]] <- rep(0, n.det.re[[q]])
        X.p.re[[q]] <- array(NA, dim = c(J.obs[[q]], n.time.max[q], 
                                         n.rep.max[q], p.det.re[[q]]))
        for (i in 1:p.det.re[[q]]) {
          X.p.re[[q]][, , , i] <- sample(1:p.RE[[q]]$levels[i], 
                                         n.time.max[q] * J.obs[q] * n.rep.max[q], 
                                         replace = TRUE)
          alpha.star[[q]][which(alpha.star.indx[[q]] == i)] <- rnorm(p.RE[[q]]$levels[i], 
                                                                     0, 
                                                                     sqrt(p.RE[[q]]$sigma.sq.p[i]))
        }
        # for (j in 1:J.obs[[q]]) {
        #   for (t in time.indx[[q]][[j]])
        #   X.p.re[[q]][j, t, -rep.indx[[q]][[j]][[t]], ] <- NA
        # }
        if (p.det.re[[q]] > 1) {
          for (j in 2:p.det.re[[q]]) {
            X.p.re[[q]][, , , j] <- X.p.re[[q]][, , , j] + max(X.p.re[[q]][, , , j - 1], 
                                                               na.rm = TRUE) 
          }
        }
          alpha.star.sites[[q]] <- apply(X.p.re[[q]], c(1, 2, 3), function(a) sum(alpha.star[[q]][a]))
      } else {
        X.p.re[[q]] <- NA
        alpha.star[[q]] <- NA
        alpha.star.sites[[q]] <- NA
      }
    } # q (n.data)
  } else {
    X.p.re <- NA
    alpha.star <- NA
  }
  
  # Simulate spatial random effects ---------------------------------------
  if (sp) {
    w.mat <- matrix(NA, J, p.svc)
    if (cov.model == 'matern') {
      theta <- cbind(phi, nu)
    } else {
      theta <- as.matrix(phi)
    }
    for (i in 1:p.svc) {
      Sigma.full <- mkSpCov(coords, as.matrix(sigma.sq[i]), as.matrix(0), theta[i, ], cov.model)
      w.mat[, i] <- rmvn(1, rep(0, nrow(Sigma.full)), Sigma.full)
    }
    X.w <- X[, , svc.cols, drop = FALSE]
    w.sites <- matrix(0, J, n.time.total)
    for (j in 1:J) {
      for (t in 1:n.time.total) {
        w.sites[j, t] <- w.mat[j, ] %*% X.w[j, t, ] 
      }
    }
  } else {
    w.mat <- NA
    w.mat.full <- NA
    w.sites <- matrix(0, J, n.time.total) 
  }
  # Simulate temporal AR random effect ------------------------------------
  if (ar1) {
    exponent <- abs(matrix(1:n.time.total - 1, nrow = n.time.total, 
                           ncol = n.time.total, byrow = TRUE) - (1:n.time.total - 1))
    Sigma.eta <- sigma.sq.t * rho^exponent
    eta <- rmvn(1, rep(0, n.time.total), Sigma.eta)
  } else {
    eta <- matrix(rep(0, n.time.total))
  }
  # Abundance process -----------------------------------------------------
  mu <- matrix(NA, J, n.time.total)
  for (j in 1:J) {
    for (t in 1:n.time.total) {
      if (length(mu.RE) > 0) {
        mu[j, t] <- exp(X[j, t, ] %*% as.matrix(beta) + w.sites[j, t] + beta.star.sites[j, t] + 
                                 eta[t])
      } else {
        mu[j, t] <- exp(X[j, t, ] %*% as.matrix(beta) + w.sites[j, t] + eta[t])
      }
    } # t
  } # j
  # Data formation --------------------------------------------------------
  p <- list()
  y <- list()
  # Offset ----------------------------
  for (i in 1:n.data) {
    K.curr <- n.rep[[i]]
    J.curr <- J.obs[[i]]
    offset.curr <- offset[[i]]
    # TODO: the offset stuff here is not correct. The offset needs to be 
    #       allowed to vary across the same dimensions as y. Just force them 
    #       to specify it as a single value or an entire array.
    # Single value
    if (length(offset.curr) == 1) {
      offset.curr <- matrix(offset.curr, J.curr, n.rep.max[i])
    } else if (length(dim(offset.curr)) == 1) { # Value for each site
      if (length(offset.curr) != J) {
        stop(paste0("offset for each data set must be a single value, vector of length ", J.curr, " or a matrix with ",
	           J.curr, " rows and ", n.rep.max[i], " columns."))
      }
      offset.curr <- matrix(offset.curr, J.curr, n.rep.max[i])
    } else if (length(dim(offset.curr)) == 2) { # Value for each site/obs
      if (nrow(offset.curr) != J.curr | ncol(offset.curr) != n.rep.max[i]) {
        stop(paste0("offset.curr must be a single value, vector of length ", J.curr, " or a matrix with ",
                    J.curr, " rows and ", n.rep.max[i], " columns."))

      }
    }
    p[[i]] <- array(NA, dim = c(J.curr, n.time.max[i], n.rep.max[i]))
    y[[i]] <- array(NA, dim = c(J.curr, n.time.max[i], n.rep.max[i]))
    sites.curr <- sites[[i]]
    seasons.curr <- data.seasons[[i]]
    X.p.curr <- X.p[[i]]
    alpha.curr <- as.matrix(alpha[[i]])
    if (length(p.RE) > 0) {
      alpha.star.sites.curr <- alpha.star.sites[[i]]
    }
    for (j in 1:J.curr) {
      for (t in time.indx[[i]][[j]]) {
        if (length(p.RE) > 0) { # If any detection random effects
          if (length(p.RE[[i]]) > 0) { # If any detection random effects in this data set
            p[[i]][j, t, rep.indx[[i]][[j]][[t]]] <- exp(X.p.curr[j, t, rep.indx[[i]][[j]][[t]], ] %*% 
                                                               alpha.curr + 
                                                               alpha.star.sites.curr[j, t, rep.indx[[i]][[j]][[t]]])
          } else { # Detection random effects, but none in this data set
            p[[i]][j, t, rep.indx[[i]][[j]][[t]]] <- exp(X.p.curr[j, t, rep.indx[[i]][[j]][[t]], ] %*% alpha.curr)
	}  
        }	else { # No detection random effects
          p[[i]][j, t, rep.indx[[i]][[j]][[t]]] <- exp(X.p.curr[j, t, rep.indx[[i]][[j]][[t]], ] %*% alpha.curr)
        }
        if (family[i] == 'NB') {
          y[[i]][j, t, rep.indx[[i]][[j]][[t]]] <- rnbinom(length(rep.indx[[i]][[j]][[t]]), size = kappa[i], mu = p[[i]][j, t, rep.indx[[i]][[j]][[t]]] * mu[sites.curr[j], seasons.curr[t]])
        }
        if (family[i] == 'Poisson') {
          y[[i]][j, t, rep.indx[[i]][[j]][[t]]] <- rpois(length(rep.indx[[i]][[j]][[t]]), lambda = p[[i]][j, t, rep.indx[[i]][[j]][[t]]] * mu[sites.curr[j], seasons.curr[t]])
        }
      }
    }
  }

  # Split into observed and predicted -----------------------------------
  sites.obs <- sort(unique(unlist(sites))) 
  sites.pred <- (1:J)[!(1:J %in% sites.obs)]
  X.obs <- X[sites.obs, , , drop = FALSE]
  X.pred <- X[sites.pred, , , drop = FALSE]
  if (length(mu.RE) > 0) {
    X.re.obs <- X.re[sites.obs, , , drop = FALSE]
    X.re.pred <- X.re[sites.pred, , , drop = FALSE]
  } else {
    X.re.obs <- NA
    X.re.pred <- NA
  }
  coords.obs <- coords[sites.obs,, drop = FALSE]
  coords.pred <- coords[sites.pred,, drop = FALSE]
  if (sp) {
    w.obs <- w.mat[sites.obs, , drop = FALSE]
    w.pred <- w.mat[sites.pred, , drop = FALSE]
  } else {
    w.obs <- NA
    w.pred <- NA
  }
  mu.obs <- mu[sites.obs, , drop = FALSE]
  mu.pred <- mu[sites.pred, , drop = FALSE]
  sites.vec <- unlist(sites)
  # Need to get index of each site value in only the sites.obs
  sites.new <- rep(0, length(sites.vec))
  for (i in 1:length(sites.vec)) {
    sites.new[i] <- which(sites.obs == sites.vec[i])
  }
  sites.return <- list()
  indx <- 1
  for (i in 1:n.data) {
    sites.return[[i]] <- sites.new[indx:(indx + J.obs[i] - 1)]
    indx <- indx + J.obs[i]
  }
  
  
  list(X.obs = X.obs, X.pred = X.pred, X.p = X.p, 
       coords.obs = coords.obs, coords.pred = coords.pred, 
       w.obs = w.obs, w.pred = w.pred, 
       mu.obs = mu.obs, mu.pred = mu.pred, 
       p = p, y = y, sites = sites.return, 
       X.re.obs = X.re.obs, X.re.pred = X.re.pred, beta.star = beta.star, 
       X.p.re = X.p.re, alpha.star = alpha.star, eta = eta 
	)
}
