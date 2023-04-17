simDS <- function(J.x, J.y, n.bins, bin.width, beta, alpha, det.func, transect = 'line',
		  kappa, mu.RE = list(), p.RE = list(), offset = 1, 
		  sp = FALSE, cov.model, sigma.sq, phi, nu, family = 'Poisson', ...) {

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
  # det.func -------------------------
  if (missing(det.func)) {
    stop("error: det.func must be specified")
  }
  det.func.names <- c('halfnormal', 'negexp')
  if (! det.func %in% det.func.names) {
    stop("error: specified det.func '", det.func, "' is not a valid option; choose from ", 
	 paste(det.func.names, collapse = ', ', sep = ''), ".")
  }
  # if (det.func == 'hazard' & missing(b)) {
  #   stop("b (scale parameter) must be specified for a hazard rate detection function")
  # }

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
  # Hazard rate detection function
  hazard <- function(x, sigma, transect, b) {
    if (transect == 'line') {
      1 - exp(-1 * (x / sigma)^(-b)) 
    } else {
      (1 - exp(-1 * (x / sigma)^(-b))) * x
    }
  }

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
  X.p <- matrix(1, nrow = J, ncol = n.alpha)
  if (n.alpha > 1) {
    for (i in 2:n.alpha) {
      X.p[, i] <- rnorm(J)
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
    X.p.random <- X[, unlist(p.RE$alpha.indx), drop = FALSE]
    n.random <- ncol(X.p.random)
    X.p.re <- matrix(NA, J, length(p.RE$levels))
    for (i in 1:length(p.RE$levels)) {
      X.p.re[, i] <- sample(1:p.RE$levels[i], J, replace = TRUE)  
    }
    indx.mat <- X.p.re[, p.re.col.indx, drop = FALSE]
    for (i in 1:p.det.re) {
      alpha.star[which(alpha.star.indx == i)] <- rnorm(n.det.re.long[i], 0, 
						     sqrt(sigma.sq.p[i]))
    }
    if (length(p.RE$levels) > 1) {
      for (j in 2:length(p.RE$levels)) {
        X.p.re[, j] <- X.p.re[, j] + max(X.p.re[, j - 1], na.rm = TRUE)
      }
    }
    if (p.det.re > 1) {
      for (j in 2:p.det.re) {
        indx.mat[, j] <- indx.mat[, j] + max(indx.mat[, j - 1], na.rm = TRUE)
      }
    }
    alpha.star.sites <- rep(NA, J)
    for (j in 1:J) {
      alpha.star.sites[j] <- alpha.star[indx.mat[j, , drop = FALSE]] %*% t(X.p.random[j, , drop = FALSE])
    }
  } else {
    X.p.re <- NA
    alpha.star <- NA
  }

  # Latent abundance process ----------------------------------------------
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
  if (family == 'NB') {
    N <- rnbinom(J, size = kappa, mu = mu * offset)
  } else {
    N <- rpois(J, lambda = mu * offset)
  }

  # Data Formation --------------------------------------------------------
  if (length(p.RE) > 0) {
    sigma <- exp(X.p %*% as.matrix(alpha) + alpha.star.sites)
  } else {
    sigma <- exp(X.p %*% as.matrix(alpha))
  }
  # Probability of detecting an individual in a given bin at a given site.
  p <- matrix(NA, J, n.bins)
  if (det.func == 'halfnormal') {
    curr.function <- halfNormal
  }
  if (det.func == 'negexp') {
    curr.function <- negExp
  }
  # if (det.func == 'hazard') {
  #   curr.function <- hazard
  # }
  # Create distance bins
  dist.breaks <- rep(0, n.bins + 1)
  tmp <- 0
  for (i in 1:n.bins) {
    tmp <- tmp + bin.width[i]
    dist.breaks[i + 1] <- tmp
  }
  # This is B in the model notation. B = half line transect width or radius of point count
  strip.width <- sum(bin.width)
  for (j in 1:J) {
    for (k in 1:n.bins) {
      # if (det.func == 'hazard') {
      #   p[j, k] <- integrate(curr.function, dist.breaks[k], 
      #   		     dist.breaks[k + 1], sig = sigma[j, ], b = b, transect = transect)$value
      # } else {
        p[j, k] <- integrate(curr.function, dist.breaks[k], 
			     dist.breaks[k + 1], sig = sigma[j, ], transect = transect)$value
      # }
      if (transect == 'line') {
        p[j, k] <- p[j, k] / (dist.breaks[k + 1] - dist.breaks[k])
      } else {
        p[j, k] <- p[j, k] * 2 / (dist.breaks[k + 1]^2 - dist.breaks[k]^2)
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
  pi.obs <- t(apply(p,  1, function(a) a * psi))
  pi.full <- cbind(pi.obs, 1 - apply(pi.obs, 1, sum))

  y <- matrix(NA, nrow = J, ncol = n.bins + 1)
  for (j in 1:J) {
    y[j, ] <- rmultinom(1, N[j], pi.full[j, ])
  }

  # Data Formation --------------------------------------------------------

  return(
    list(X = X, X.p = X.p, coords = coords, w = w, mu = mu, N = N,
	 y = y[, -(n.bins + 1), drop = FALSE], X.re = X.re, 
	 X.p.re = X.p.re, beta.star = beta.star, sigma = sigma,
	 pi.full = pi.full, alpha.star = alpha.star, dist.breaks = dist.breaks, 
         p = p)
  )
}
