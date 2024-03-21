DS <- function(abund.formula, det.formula, data, inits, priors, tuning,
               n.batch, batch.length, accept.rate = 0.43, family = 'Poisson',
	       transect = 'line', det.func = 'halfnormal',
               n.omp.threads = 1, verbose = TRUE,
               n.report = 100, n.burn = round(.10 * n.batch * batch.length), n.thin = 1,
               n.chains = 1, ...){

  ptm <- proc.time()

  # Make it look nice
  if (verbose) {
    cat("----------------------------------------\n");
    cat("\tPreparing to run the model\n");
    cat("----------------------------------------\n");
  }

  # Functions ---------------------------------------------------------------
  logit <- function(theta, a = 0, b = 1) {log((theta-a)/(b-theta))}
  logit.inv <- function(z, a = 0, b = 1) {b-(b-a)/(1+exp(z))}
  rigamma <- function(n, a, b){
    1/rgamma(n = n, shape = a, rate = b)
  }
  # Half-normal detection function
  halfNormal <- function(x, sigma) {
    exp(-x^2 / (2 * sigma^2))
  }
  # Negative exponential detection function
  negExp <- function(x, sigma) {
    exp(-x / sigma)
  }

  # Check for unused arguments ------------------------------------------
  formal.args <- names(formals(sys.function(sys.parent())))
  elip.args <- names(list(...))
  for(i in elip.args){
      if(! i %in% formal.args)
          warning("'",i, "' is not an argument")
  }
  # Call ----------------------------------------------------------------
  # Returns a call in which all of the specified arguments are
  # specified by their full names.
  cl <- match.call()

  # Some initial checks -------------------------------------------------
  if (missing(data)) {
    stop("error: data must be specified")
  }
  if (!is.list(data)) {
    stop("error: data must be a list")
  }
  names(data) <- tolower(names(data))
  if (missing(abund.formula)) {
    stop("error: abund.formula must be specified")
  }
  if (missing(det.formula)) {
    stop("error: det.formula must be specified")
  }
  if (!'y' %in% names(data)) {
    stop("error: data y must be specified in data")
  }
  y <- as.matrix(data$y)
  y.mat <- y
  # Offset
  if ('offset' %in% names(data)) {
    offset <- data$offset
    if (length(offset) != nrow(y) & length(offset) != 1) {
      stop(paste("error: data$offset must be of length 1 or ", nrow(y), sep = ''))
    }
    if (length(offset) == 1) {
      offset <- rep(offset, nrow(y))
    }
  } else {
    offset <- rep(1, nrow(y))
  }
  if (!'covs' %in% names(data)) {
    if ((abund.formula == ~ 1) & (det.formula == ~ 1)) {
      if (verbose) {
        message("covariates (covs) not specified in data.\nAssuming intercept only distance sampling model.\n")
      }
      data$covs <- matrix(1, dim(y)[1], 1)
    } else {
      stop("error: covs must be specified in data for a distance sampling model with covariates")
    }
  }
  if (!is.matrix(data$covs) & !is.data.frame(data$covs)) {
    stop("error: covs must be a matrix or data frame")
  }
  if (sum(is.na(data$covs)) > 0) {
    stop("error: missing covariate values in data$covs. Remove these sites from all data or impute non-missing values.")
  }
  if (!'dist.breaks' %in% names(data)) {
    stop("error: distance cut off points (dist.breaks) must be specified in data")
  }
  if (length(data$dist.breaks) != (ncol(y) + 1)) {
    stop(paste('error: dist.breaks must be of length ', ncol(y) + 1, '.', sep = ''))
  }
  dist.breaks <- data$dist.breaks
  if (missing(n.batch)) {
    stop("error: must specify number of MCMC batches")
  }
  if (missing(batch.length)) {
    stop("error: must specify length of each MCMC batch")
  }
  n.samples <- n.batch * batch.length
  if (n.burn > n.samples) {
    stop("error: n.burn must be less than n.samples")
  }
  if (n.thin > n.samples) {
    stop("error: n.thin must be less than n.samples")
  }

  if (!(family) %in% c('Poisson', 'NB')) {
    stop("family must be either 'Poisson' or 'NB'")
  }

  data$covs <- as.data.frame(data$covs)

  # Check whether random effects are sent in as numeric, and
  # return error if they are.
  # Abundance -------------------------
  if (!is.null(findbars(abund.formula))) {
    abund.re.names <- unique(unlist(sapply(findbars(abund.formula), all.vars)))
    for (i in 1:length(abund.re.names)) {
      if (is(data$covs[, abund.re.names[i]], 'factor')) {
        stop(paste("error: random effect variable ", abund.re.names[i], " specified as a factor. Random effect variables must be specified as numeric.", sep = ''))
      }
      if (is(data$covs[, abund.re.names[i]], 'character')) {
        stop(paste("error: random effect variable ", abund.re.names[i], " specified as character. Random effect variables must be specified as numeric.", sep = ''))
      }
    }
  }
  # Detection -----------------------
  if (!is.null(findbars(det.formula))) {
    det.re.names <- unique(unlist(sapply(findbars(det.formula), all.vars)))
    for (i in 1:length(det.re.names)) {
      if (is(data$covs[, det.re.names[i]], 'factor')) {
        stop(paste("error: random effect variable ", det.re.names[i], " specified as a factor. Random effect variables must be specified as numeric.", sep = ''))
      }
      if (is(data$covs[, det.re.names[i]], 'character')) {
        stop(paste("error: random effect variable ", det.re.names[i], " specified as character. Random effect variables must be specified as numeric.", sep = ''))
      }
    }
  }

  # Checking missing values ---------------------------------------------
  # y -------------------------------
  y.na.test <- apply(y, 1, function(a) sum(!is.na(a)))
  if (sum(is.na(data$y)) > 0) {
    stop("error: missing values are not allowed in y for distance sampling models.")
  }
  # covs ------------------------
  if (sum(is.na(data$covs)) != 0) {
    stop("error: missing values in covs. Please remove these sites from all objects in data or somehow replace the NA values with non-missing values (e.g., mean imputation).")
  }

  # Formula -------------------------------------------------------------
  # Abundance -------------------------
  if (is(abund.formula, 'formula')) {
    tmp <- parseFormula(abund.formula, data$covs)
    X <- as.matrix(tmp[[1]])
    X.re <- as.matrix(tmp[[4]])
    x.re.names <- colnames(X.re)
    x.names <- tmp[[2]]
    X.random <- as.matrix(tmp[[5]])
    x.random.names <- colnames(X.random)
  } else {
    stop("error: abund.formula is misspecified")
  }
  # Get RE level names
  re.level.names <- lapply(data$covs[, x.re.names, drop = FALSE],
      		     function (a) sort(unique(a)))
  x.re.names <- x.random.names

  # Detection -----------------------
  if (is(det.formula, 'formula')) {
    tmp <- parseFormula(det.formula, data$covs)
    X.p <- as.matrix(tmp[[1]])
    X.p.re <- as.matrix(tmp[[4]])
    x.p.re.names <- colnames(X.p.re)
    x.p.names <- tmp[[2]]
    X.p.random <- as.matrix(tmp[[5]])
    x.p.random.names <- colnames(X.p.random)
  } else {
    stop("error: det.formula is misspecified")
  }
  p.re.level.names <- lapply(data$covs[, x.p.re.names, drop = FALSE],
                             function (a) sort(unique(a)))
  x.p.re.names <- x.p.random.names

  # Get basic info from inputs ------------------------------------------
  # Number of sites
  J <- nrow(y)
  # Number of abundance parameters
  p.abund <- ncol(X)
  # Number of abundance random effect parameters
  p.abund.re <- ncol(X.re)
  # Number of detection parameters
  p.det <- ncol(X.p)
  # Number of detection random effect parameters
  p.det.re <- ncol(X.p.re)
  # Number of latent abundance random effect values
  n.abund.re <- length(unlist(apply(X.re, 2, unique)))
  n.abund.re.long <- apply(X.re, 2, function(a) length(unique(a)))
  # Number of latent detection random effect values
  n.det.re <- length(unlist(apply(X.p.re, 2, unique)))
  n.det.re.long <- apply(X.p.re, 2, function(a) length(unique(a)))
  if (p.det.re == 0) n.det.re.long <- 0
  # Number of distance bands/bins
  K <- ncol(y)

  # Just to keep things consistent with other functions
  N.long.indx <- rep(1:J, dim(y)[2])
  N.long.indx <- N.long.indx[!is.na(c(y))]
  # Subtract 1 for indices in C
  N.long.indx <- N.long.indx - 1
  # Note that y is ordered by distance bin, then site within bin.
  y <- c(y)
  # Total number of data points
  n.obs <- J * K

  # Get random effect matrices all set ----------------------------------
  X.re <- X.re - 1
  if (p.abund.re > 1) {
    for (j in 2:p.abund.re) {
      X.re[, j] <- X.re[, j] + max(X.re[, j - 1]) + 1
    }
  }
  X.p.re <- X.p.re - 1
  if (p.det.re > 1) {
    for (j in 2:p.det.re) {
      X.p.re[, j] <- X.p.re[, j] + max(X.p.re[, j - 1]) + 1
    }
  }

  # Grab specific distance sampling information ---------------------------
  det.func.names <- c("halfnormal", "negexp")
  if (! det.func %in% det.func.names) {
    stop("error: specified det.func '", det.func, "' is not a valid option; choose from ",
	 paste(det.func.names, collapse = ', ', sep = ''), ".")
  }
  # Obo for det.func lookup on c side
  # halfnormal = 0, negexp = 1
  det.func.indx <- which(det.func == det.func.names) - 1
  if (! transect %in% c('line', 'point')) {
    stop("error: transect must be either 'line', or 'point'")
  }
  # For C side, line = 0, point = 1
  transect.c <- ifelse(transect == 'line', 0, 1)

  # Priors --------------------------------------------------------------
  if (missing(priors)) {
    priors <- list()
  }
  names(priors) <- tolower(names(priors))
  # beta -----------------------
  if ("beta.normal" %in% names(priors)) {
    if (!is.list(priors$beta.normal) | length(priors$beta.normal) != 2) {
      stop("error: beta.normal must be a list of length 2")
    }
    mu.beta <- priors$beta.normal[[1]]
    sigma.beta <- priors$beta.normal[[2]]
    if (length(mu.beta) != p.abund & length(mu.beta) != 1) {
      if (p.abund == 1) {
        stop(paste("error: beta.normal[[1]] must be a vector of length ",
        	     p.abund, " with elements corresponding to betas' mean", sep = ""))
      } else {
        stop(paste("error: beta.normal[[1]] must be a vector of length ",
        	     p.abund, " or 1 with elements corresponding to betas' mean", sep = ""))
      }
    }
    if (length(sigma.beta) != p.abund & length(sigma.beta) != 1) {
      if (p.abund == 1) {
        stop(paste("error: beta.normal[[2]] must be a vector of length ",
      	   p.abund, " with elements corresponding to betas' variance", sep = ""))
      } else {
        stop(paste("error: beta.normal[[2]] must be a vector of length ",
      	   p.abund, " or 1 with elements corresponding to betas' variance", sep = ""))
      }
    }
    if (length(sigma.beta) != p.abund) {
      sigma.beta <- rep(sigma.beta, p.abund)
    }
    if (length(mu.beta) != p.abund) {
      mu.beta <- rep(mu.beta, p.abund)
    }
    Sigma.beta <- sigma.beta * diag(p.abund)
  } else {
    if (verbose) {
      message("No prior specified for beta.normal.\nSetting prior mean to 0 and prior variance to 100\n")
    }
    mu.beta <- rep(0, p.abund)
    sigma.beta <- rep(100, p.abund)
    Sigma.beta <- diag(p.abund) * sigma.beta
  }
  # alpha -----------------------
  if ("alpha.normal" %in% names(priors)) {
    if (!is.list(priors$alpha.normal) | length(priors$alpha.normal) != 2) {
      stop("error: alpha.normal must be a list of length 2")
    }
    mu.alpha <- priors$alpha.normal[[1]]
    sigma.alpha <- priors$alpha.normal[[2]]
    if (length(mu.alpha) != p.det & length(mu.alpha) != 1) {
      if (p.det == 1) {
        stop(paste("error: alpha.normal[[1]] must be a vector of length ",
        	     p.det, " with elements corresponding to alphas' mean", sep = ""))
      } else {
        stop(paste("error: alpha.normal[[1]] must be a vector of length ",
        	     p.det, " or 1 with elements corresponding to alphas' mean", sep = ""))
      }
    }
    if (length(sigma.alpha) != p.det & length(sigma.alpha) != 1) {
      if (p.det == 1) {
        stop(paste("error: alpha.normal[[2]] must be a vector of length ",
      	   p.det, " with elements corresponding to alphas' variance", sep = ""))
      } else {
        stop(paste("error: alpha.normal[[2]] must be a vector of length ",
      	   p.det, " or 1 with elements corresponding to alphas' variance", sep = ""))
      }
    }
    if (length(sigma.alpha) != p.det) {
      sigma.alpha <- rep(sigma.alpha, p.det)
    }
    if (length(mu.alpha) != p.det) {
      mu.alpha <- rep(mu.alpha, p.det)
    }
    Sigma.alpha <- sigma.alpha * diag(p.det)
  } else {
    if (verbose) {
      message("No prior specified for alpha.normal.\nSetting prior mean to 0 and prior variance to 100\n")
    }
    mu.alpha <- rep(0, p.det)
    sigma.alpha <- rep(100, p.det)
    Sigma.alpha <- diag(p.det) * 100
  }
  # sigma.sq.mu --------------------
  if (p.abund.re > 0) {
    if ("sigma.sq.mu.ig" %in% names(priors)) {
      if (!is.list(priors$sigma.sq.mu.ig) | length(priors$sigma.sq.mu.ig) != 2) {
        stop("error: sigma.sq.mu.ig must be a list of length 2")
      }
      sigma.sq.mu.a <- priors$sigma.sq.mu.ig[[1]]
      sigma.sq.mu.b <- priors$sigma.sq.mu.ig[[2]]
      if (length(sigma.sq.mu.a) != p.abund.re & length(sigma.sq.mu.a) != 1) {
        if (p.abund.re == 1) {
        stop(paste("error: sigma.sq.mu.ig[[1]] must be a vector of length ",
        	   p.abund.re, " with elements corresponding to sigma.sq.mus' shape", sep = ""))
        } else {
        stop(paste("error: sigma.sq.mu.ig[[1]] must be a vector of length ",
        	   p.abund.re, " or 1 with elements corresponding to sigma.sq.mus' shape", sep = ""))
        }
      }
      if (length(sigma.sq.mu.b) != p.abund.re & length(sigma.sq.mu.b) != 1) {
        if (p.abund.re == 1) {
          stop(paste("error: sigma.sq.mu.ig[[2]] must be a vector of length ",
        	   p.abund.re, " with elements corresponding to sigma.sq.mus' scale", sep = ""))
        } else {
          stop(paste("error: sigma.sq.mu.ig[[2]] must be a vector of length ",
        	   p.abund.re, " or 1with elements corresponding to sigma.sq.mus' scale", sep = ""))
        }
      }
      if (length(sigma.sq.mu.a) != p.abund.re) {
        sigma.sq.mu.a <- rep(sigma.sq.mu.a, p.abund.re)
      }
      if (length(sigma.sq.mu.b) != p.abund.re) {
        sigma.sq.mu.b <- rep(sigma.sq.mu.b, p.abund.re)
      }
  }   else {
      if (verbose) {
        message("No prior specified for sigma.sq.mu.ig.\nSetting prior shape to 0.1 and prior scale to 0.1\n")
      }
      sigma.sq.mu.a <- rep(0.1, p.abund.re)
      sigma.sq.mu.b <- rep(0.1, p.abund.re)
    }
  } else {
    sigma.sq.mu.a <- 0
    sigma.sq.mu.b <- 0
  }
  # sigma.sq.p --------------------
  if (p.det.re > 0) {
    if ("sigma.sq.p.ig" %in% names(priors)) {
      if (!is.list(priors$sigma.sq.p.ig) | length(priors$sigma.sq.p.ig) != 2) {
        stop("error: sigma.sq.p.ig must be a list of length 2")
      }
      sigma.sq.p.a <- priors$sigma.sq.p.ig[[1]]
      sigma.sq.p.b <- priors$sigma.sq.p.ig[[2]]
      if (length(sigma.sq.p.a) != p.det.re & length(sigma.sq.p.a) != 1) {
        if (p.det.re == 1) {
          stop(paste("error: sigma.sq.p.ig[[1]] must be a vector of length ",
        	   p.det.re, " with elements corresponding to sigma.sq.ps' shape", sep = ""))
        } else {
          stop(paste("error: sigma.sq.p.ig[[1]] must be a vector of length ",
        	   p.det.re, " or 1 with elements corresponding to sigma.sq.ps' shape", sep = ""))
        }
      }
      if (length(sigma.sq.p.b) != p.det.re & length(sigma.sq.p.b) != 1) {
        if (p.det.re == 1) {
          stop(paste("error: sigma.sq.p.ig[[2]] must be a vector of length ",
        	     p.det.re, " with elements corresponding to sigma.sq.ps' scale", sep = ""))
        } else {
          stop(paste("error: sigma.sq.p.ig[[2]] must be a vector of length ",
        	     p.det.re, " or 1 with elements corresponding to sigma.sq.ps' scale", sep = ""))
        }
      }
      if (length(sigma.sq.p.a) != p.det.re) {
        sigma.sq.p.a <- rep(sigma.sq.p.a, p.det.re)
      }
      if (length(sigma.sq.p.b) != p.det.re) {
        sigma.sq.p.b <- rep(sigma.sq.p.b, p.det.re)
      }
  }   else {
      if (verbose) {
        message("No prior specified for sigma.sq.p.ig.\nSetting prior shape to 0.1 and prior scale to 0.1\n")
      }
      sigma.sq.p.a <- rep(0.1, p.det.re)
      sigma.sq.p.b <- rep(0.1, p.det.re)
    }
  } else {
    sigma.sq.p.a <- 0
    sigma.sq.p.b <- 0
  }
  # kappa -----------------------------
  if (family == 'NB') {
    if ("kappa.unif" %in% names(priors)) {
      if (!is.vector(priors$kappa.unif) | !is.atomic(priors$kappa.unif) | length(priors$kappa.unif) != 2) {
        stop("error: kappa.unif must be a vector of length 2 with elements corresponding to kappa's lower and upper bounds")
      }
      kappa.a <- priors$kappa.unif[1]
      kappa.b <- priors$kappa.unif[2]
    } else {
      if (verbose) {
        message("No prior specified for kappa.unif.\nSetting uniform bounds of 0 and 100.\n")
      }
      kappa.a <- 0
      kappa.b <- 100
    }
  } else {
    kappa.a <- 0
    kappa.b <- 0
  }
  # Starting values -----------------------------------------------------
  if (missing(inits)) {
    inits <- list()
  }
  names(inits) <- tolower(names(inits))
  # N -------------------------------
  if ("n" %in% names(inits)) {
    N.inits <- inits$n
    if (!is.vector(N.inits)) {
      stop(paste("error: initial values for N must be a vector of length ",
      	   J, sep = ""))
    }
    if (length(N.inits) != J) {
      stop(paste("error: initial values for N must be a vector of length ",
      	   J, sep = ""))
    }
    N.test <- apply(y.mat, 1, sum, na.rm = TRUE)
    init.test <- sum(N.inits < N.test)
    if (init.test > 0) {
      stop("error: initial values for latent abundance (N) are invalid. Please re-specify inits$N so initial values are greater than or equal to the total number of observed individuals observed at a given site.")
    }
  } else {
    N.inits <- apply(y.mat, 1, sum, na.rm = TRUE)
    if (verbose) {
      message("N is not specified in initial values.\nSetting initial values based on observed data\n")
    }
  }
  # beta -----------------------
  if ("beta" %in% names(inits)) {
    beta.inits <- inits[["beta"]]
    if (length(beta.inits) != p.abund & length(beta.inits) != 1) {
      if (p.abund == 1) {
        stop(paste("error: initial values for beta must be of length ", p.abund,
      	     sep = ""))

      } else {
        stop(paste("error: initial values for beta must be of length ", p.abund, " or 1",
        	     sep = ""))
      }
    }
    if (length(beta.inits) != p.abund) {
      beta.inits <- rep(beta.inits, p.abund)
    }
  } else {
    beta.inits <- rnorm(p.abund)
    if (verbose) {
      message('beta is not specified in initial values.\nSetting initial values to random values from a standard normal distribution\n')
    }
  }
  # alpha -----------------------
  if ("alpha" %in% names(inits)) {
    alpha.input <- TRUE
    alpha.inits <- inits[["alpha"]]
    if (length(alpha.inits) != p.det & length(alpha.inits) != 1) {
      if (p.det == 1) {
      stop(paste("error: initial values for alpha must be of length ", p.det,
      	   sep = ""))
      } else {
        stop(paste("error: initial values for alpha must be of length ", p.det, " or 1",
      	     sep = ""))
      }
    }
    if (length(alpha.inits) != p.det) {
      alpha.inits <- rep(alpha.inits, p.det)
    }
  } else {
    alpha.input <- FALSE
    alpha.inits <- runif(p.det, -10, 10)
    if (verbose) {
      message("alpha is not specified in initial values.\nSetting initial values to random values from a Uniform(-10, 10)\n")
    }
  }

  # sigma.sq.mu -------------------
  if (p.abund.re > 0) {
    if ("sigma.sq.mu" %in% names(inits)) {
      sigma.sq.mu.inits <- inits[["sigma.sq.mu"]]
      if (length(sigma.sq.mu.inits) != p.abund.re & length(sigma.sq.mu.inits) != 1) {
        if (p.abund.re == 1) {
          stop(paste("error: initial values for sigma.sq.mu must be of length ", p.abund.re,
      	       sep = ""))
        } else {
          stop(paste("error: initial values for sigma.sq.mu must be of length ", p.abund.re,
      	       " or 1", sep = ""))
        }
      }
      if (length(sigma.sq.mu.inits) != p.abund.re) {
        sigma.sq.mu.inits <- rep(sigma.sq.mu.inits, p.abund.re)
      }
    } else {
      sigma.sq.mu.inits <- runif(p.abund.re, 0.05, 1)
      if (verbose) {
        message("sigma.sq.mu is not specified in initial values.\nSetting initial values to random values between 0.05 and 1\n")
      }
    }
    beta.star.indx <- rep(0:(p.abund.re - 1), n.abund.re.long)
    beta.star.inits <- rnorm(n.abund.re, sqrt(sigma.sq.mu.inits[beta.star.indx + 1]))
  } else {
    sigma.sq.mu.inits <- 0
    beta.star.indx <- 0
    beta.star.inits <- 0
  }
  # sigma.sq.p ------------------
  if (p.det.re > 0) {
    if ("sigma.sq.p" %in% names(inits)) {
      sigma.sq.p.inits <- inits[["sigma.sq.p"]]
      if (length(sigma.sq.p.inits) != p.det.re & length(sigma.sq.p.inits) != 1) {
        if (p.det.re == 1) {
          stop(paste("error: initial values for sigma.sq.p must be of length ", p.det.re,
      	     sep = ""))
        } else {
          stop(paste("error: initial values for sigma.sq.p must be of length ", p.det.re,
      	       " or 1", sep = ""))

        }
      }
      if (length(sigma.sq.p.inits) != p.det.re) {
        sigma.sq.p.inits <- rep(sigma.sq.p.inits, p.det.re)
      }
    } else {
      sigma.sq.p.inits <- runif(p.det.re, 0.05, 0.5)
      if (verbose) {
        message("sigma.sq.p is not specified in initial values.\nSetting initial values to random values between 0.05 and 0.5\n")
      }
    }
    alpha.star.indx <- rep(0:(p.det.re - 1), n.det.re.long)
    alpha.star.inits <- rep(0, n.det.re)
  } else {
    sigma.sq.p.inits <- 0
    alpha.star.indx <- 0
    alpha.star.inits <- 0
  }
  # kappa ---------------------------
  if (family == 'NB') {
    if ("kappa" %in% names(inits)) {
      kappa.inits <- inits[["kappa"]]
      if (length(kappa.inits) != 1) {
        stop("error: initial values for kappa must be of length 1")
      }
    } else {
      kappa.inits <- runif(1, kappa.a, kappa.b)
      if (verbose) {
        message("kappa is not specified in initial values.\nSetting initial value to random value from the prior distribution\n")
      }
    }
  } else {
    kappa.inits <- 0
  }
  # Should initial values be fixed --
  if ("fix" %in% names(inits)) {
    fix.inits <- inits[["fix"]]
    if ((fix.inits != TRUE) & (fix.inits != FALSE)) {
      stop(paste("error: inits$fix must take value TRUE or FALSE"))
    }
  } else {
    fix.inits <- FALSE
  }
  if (verbose & fix.inits & (n.chains > 1)) {
    message("Fixing initial values across all chains\n")
  }

  # Get tuning values ---------------------------------------------------
  if (missing(tuning)) {
    beta.tuning <- rep(1, p.abund)
    beta.star.tuning <- rep(1, n.abund.re)
    alpha.tuning <- rep(0.5, p.det)
    alpha.star.tuning <- rep(1, n.det.re)
    kappa.tuning <- 1
  } else {
    names(tuning) <- tolower(names(tuning))
    # beta ---------------------------
    if(!"beta" %in% names(tuning)) {
      stop("error: beta must be specified in tuning value list")
    }
    beta.tuning <- tuning$beta
    if (length(beta.tuning) != 1 & length(beta.tuning) != p.abund) {
      stop(paste("error: beta tuning must be a single value or a vector of length ",
        	 p.abund, sep = ''))
    }
    if (length(beta.tuning) == 1) {
      beta.tuning <- rep(beta.tuning, p.abund)
    }
    if (p.abund.re > 0) {
      # beta.star ---------------------------
      if(!"beta.star" %in% names(tuning)) {
        stop("error: beta.star must be specified in tuning value list")
      }
      beta.star.tuning <- tuning$beta.star
      if (length(beta.star.tuning) != 1) {
        stop("error: beta.star tuning must be a single value")
      }
      beta.star.tuning <- rep(beta.star.tuning, n.abund.re)
    } else {
      beta.star.tuning <- NULL
    }
    # alpha ---------------------------
    if(!"alpha" %in% names(tuning)) {
      stop("error: alpha must be specified in tuning value list")
    }
    alpha.tuning <- tuning$alpha
    if (length(alpha.tuning) != 1 & length(alpha.tuning) != p.det) {
      stop(paste("error: alpha tuning must be a single value or a vector of length ",
        	 p.det, sep = ''))
    }
    if (length(alpha.tuning) == 1) {
      alpha.tuning <- rep(alpha.tuning, p.det)
    }
    if (p.det.re > 0) {
      # alpha.star ---------------------------
      if(!"alpha.star" %in% names(tuning)) {
        stop("error: alpha.star must be specified in tuning value list")
      }
      alpha.star.tuning <- tuning$alpha.star
      if (length(alpha.star.tuning) != 1) {
        stop("error: alpha.star tuning must be a single value")
      }
      alpha.star.tuning <- rep(alpha.star.tuning, n.det.re)
    } else {
      alpha.star.tuning <- NULL
    }
    if (family == 'NB') {
      # kappa ---------------------------
      if(!"kappa" %in% names(tuning)) {
        stop("error: kappa must be specified in tuning value list")
      }
      kappa.tuning <- tuning$kappa
      if (length(kappa.tuning) != 1) {
        stop("error: kappa tuning must be a single value")
      }
    } else {
      kappa.tuning <- NULL
    }
  }
  tuning.c <- log(c(beta.tuning, alpha.tuning,
		    beta.star.tuning, alpha.star.tuning,
		    kappa.tuning))
  curr.chain <- 1

  # Get max y values for N update -----------------------------------------
  # Actually a sum, but just keeping as y.max for consistency with NMix()
  y.max <- apply(y.mat, 1, sum, na.rm = TRUE)

  # Other miscellaneous ---------------------------------------------------
  # For prediction with random slopes
  re.cols <- list()
  if (p.abund.re > 0) {
    split.names <- strsplit(x.re.names, "[-]")
    for (j in 1:p.abund.re) {
      re.cols[[j]] <- split.names[[j]][1]
      names(re.cols)[j] <- split.names[[j]][2]
    }
  }
  re.det.cols <- list()
  if (p.det.re > 0) {
    split.names <- strsplit(x.p.re.names, "[-]")
    for (j in 1:p.det.re) {
      re.det.cols[[j]] <- split.names[[j]][1]
      names(re.det.cols)[j] <- split.names[[j]][2]
    }
  }

  # Set storage for all variables ---------------------------------------
  storage.mode(y) <- "double"
  storage.mode(N.inits) <- "double"
  storage.mode(X) <- "double"
  storage.mode(X.p) <- "double"
  storage.mode(y.max) <- "double"
  storage.mode(offset) <- "double"
  consts <- c(J, n.obs, p.abund, p.abund.re, n.abund.re,
	      p.det, p.det.re, n.det.re)
  storage.mode(consts) <- "integer"
  storage.mode(K) <- "integer"
  storage.mode(beta.inits) <- "double"
  storage.mode(alpha.inits) <- "double"
  storage.mode(kappa.inits) <- "double"
  storage.mode(N.long.indx) <- "integer"
  storage.mode(mu.beta) <- "double"
  storage.mode(Sigma.beta) <- "double"
  storage.mode(mu.alpha) <- "double"
  storage.mode(Sigma.alpha) <- "double"
  storage.mode(kappa.a) <- "double"
  storage.mode(kappa.b) <- "double"
  storage.mode(n.batch) <- "integer"
  storage.mode(batch.length) <- "integer"
  storage.mode(accept.rate) <- "double"
  storage.mode(tuning.c) <- "double"
  storage.mode(n.omp.threads) <- "integer"
  storage.mode(verbose) <- "integer"
  storage.mode(n.report) <- "integer"
  chain.info <- c(curr.chain, n.chains)
  storage.mode(chain.info) <- "integer"
  n.post.samples <- length(seq(from = n.burn + 1,
      			 to = n.samples,
      			 by = as.integer(n.thin)))
  storage.mode(n.post.samples) <- "integer"
  samples.info <- c(n.burn, n.thin, n.post.samples)
  storage.mode(samples.info) <- "integer"
  # For detection random effects
  storage.mode(X.p.re) <- "integer"
  storage.mode(X.p.random) <- "double"
  alpha.level.indx <- sort(unique(c(X.p.re)))
  storage.mode(alpha.level.indx) <- "integer"
  storage.mode(n.det.re.long) <- "integer"
  storage.mode(sigma.sq.p.inits) <- "double"
  storage.mode(sigma.sq.p.a) <- "double"
  storage.mode(sigma.sq.p.b) <- "double"
  storage.mode(alpha.star.inits) <- "double"
  storage.mode(alpha.star.indx) <- "integer"
  # For abundance random effects
  storage.mode(X.re) <- "integer"
  storage.mode(X.random) <- "double"
  beta.level.indx <- sort(unique(c(X.re)))
  storage.mode(beta.level.indx) <- "integer"
  storage.mode(sigma.sq.mu.inits) <- "double"
  storage.mode(sigma.sq.mu.a) <- "double"
  storage.mode(sigma.sq.mu.b) <- "double"
  storage.mode(beta.star.inits) <- "double"
  storage.mode(beta.star.indx) <- "integer"
  # NB = 1, Poisson = 0
  family.c <- ifelse(family == 'NB', 1, 0)
  storage.mode(family.c) <- "integer"
  # Distance sampling information
  storage.mode(det.func.indx) <- "integer"
  storage.mode(transect.c) <- 'integer'
  storage.mode(dist.breaks) <- 'double'

  # Fit the model -------------------------------------------------------
  out.tmp <- list()
  out <- list()
  for (i in 1:n.chains) {
    # Change initial values if i > 1
    if ((i > 1) & (!fix.inits)) {
      beta.inits <- rnorm(p.abund)
      alpha.inits <- runif(p.det, -10, 10)
      if (p.abund.re > 0) {
        sigma.sq.mu.inits <- runif(p.abund.re, 0.05, 1)
        beta.star.inits <- rnorm(n.abund.re, sqrt(sigma.sq.mu.inits[beta.star.indx + 1]))
      }
      if (p.det.re > 0) {
        sigma.sq.p.inits <- runif(p.det.re, 0.05, 0.5)
        alpha.star.inits <- runif(n.det.re, -0.5, 0.5)
      }
      if (family == 'NB') {
        kappa.inits <- runif(1, kappa.a, kappa.b)
      }
    }
    storage.mode(chain.info) <- "integer"
    # Check alpha initial values
    tmp <- .Call("checkAlphaDS", y, X.p, X.p.re, X.p.random, y.max,
			 consts, K, n.det.re.long, alpha.inits,
			 sigma.sq.p.inits, alpha.star.inits, N.inits,
			 N.long.indx, alpha.star.indx, alpha.level.indx,
			 mu.alpha, Sigma.alpha, det.func.indx, transect.c, dist.breaks)
    alpha.check <- ifelse(is.nan(tmp$alpha.like.val)[1, 1], TRUE, FALSE)
    if (i == 1 & alpha.input & alpha.check & verbose) {
      message("user-supplied initial values for alpha result in an invalid\nlikelihood. Re-drawing alpha initial values from a Uniform(-10, 10).")
    }
    while(alpha.check) {
      alpha.inits <- runif(p.det, -10, 10)
      tmp <- .Call("checkAlphaDS", y, X.p, X.p.re, X.p.random, y.max,
                   consts, K, n.det.re.long, alpha.inits,
                   sigma.sq.p.inits, alpha.star.inits, N.inits,
                   N.long.indx, alpha.star.indx, alpha.level.indx,
                   mu.alpha, Sigma.alpha, det.func.indx, transect.c, dist.breaks)
      alpha.check <- ifelse(is.nan(tmp$alpha.like.val)[1, 1], TRUE, FALSE)
    }
    # Run the model in C
    out.tmp[[i]] <- .Call("DS", y, X, X.p, X.re, X.p.re, X.random, X.p.random,
        		  y.max, offset, consts, K, n.abund.re.long,
        		  n.det.re.long, beta.inits, alpha.inits, kappa.inits,
        		  sigma.sq.mu.inits, sigma.sq.p.inits, beta.star.inits,
        		  alpha.star.inits, N.inits, N.long.indx, beta.star.indx,
      		          beta.level.indx, alpha.star.indx, alpha.level.indx,
        		  mu.beta, Sigma.beta, mu.alpha, Sigma.alpha,
      		          sigma.sq.mu.a, sigma.sq.mu.b,
        		  sigma.sq.p.a, sigma.sq.p.b, kappa.a, kappa.b,
        		  det.func.indx, transect.c, dist.breaks,
      		          tuning.c, n.batch, batch.length, accept.rate,
      		          n.omp.threads, verbose, n.report, samples.info, chain.info, family.c)
    chain.info[1] <- chain.info[1] + 1
  } # i
  # Calculate R-Hat ---------------
  out <- list()
  out$rhat <- list()
  if (n.chains > 1) {
    # as.vector removes the "Upper CI" when there is only 1 variable.
    out$rhat$beta <- as.vector(gelman.diag(mcmc.list(lapply(out.tmp, function(a)
    					      mcmc(t(a$beta.samples)))),
    			     autoburnin = FALSE, multivariate = FALSE)$psrf[, 2])
    out$rhat$alpha <- as.vector(gelman.diag(mcmc.list(lapply(out.tmp, function(a)
    					      mcmc(t(a$alpha.samples)))),
    			      autoburnin = FALSE, multivariate = FALSE)$psrf[, 2])
    if (p.det.re > 0) {
    out$rhat$sigma.sq.p <- as.vector(gelman.diag(mcmc.list(lapply(out.tmp, function(a)
    					      mcmc(t(a$sigma.sq.p.samples)))),
    			     autoburnin = FALSE, multivariate = FALSE)$psrf[, 2])
    }
    if (p.abund.re > 0) {
    out$rhat$sigma.sq.mu <- as.vector(gelman.diag(mcmc.list(lapply(out.tmp, function(a)
    					      mcmc(t(a$sigma.sq.mu.samples)))),
    			     autoburnin = FALSE, multivariate = FALSE)$psrf[, 2])
    }
    if (family == 'NB') {
        out$rhat$kappa <- as.vector(gelman.diag(mcmc.list(lapply(out.tmp, function(a)
        						       mcmc(t(a$kappa.samples)))),
        					autoburnin = FALSE, multivariate = FALSE)$psrf[, 2])
    }
  } else {
    out$rhat$beta <- rep(NA, p.abund)
    out$rhat$kappa <- NA
    out$rhat$alpha <- rep(NA, p.det)
    if (p.det.re > 0) {
      out$rhat$sigma.sq.p <- rep(NA, p.det.re)
    }
    if (p.abund.re > 0) {
      out$rhat$sigma.sq.mu <- rep(NA, p.abund.re)
    }
  }
  # Put everything into MCMC objects
  out$beta.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$beta.samples))))
  colnames(out$beta.samples) <- x.names
  out$alpha.samples <- mcmc(do.call(rbind,
    				lapply(out.tmp, function(a) t(a$alpha.samples))))
  colnames(out$alpha.samples) <- x.p.names
  if (family == 'NB') {
    out$kappa.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$kappa.samples))))
    colnames(out$kappa.samples) <- c("kappa")
  }
  out$N.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$N.samples))))
  out$mu.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$mu.samples))))
  out$y.rep.samples <- do.call(abind, lapply(out.tmp, function(a) array(a$y.rep.samples,
        								c(K + 1, J, n.post.samples))))
  out$y.rep.samples <- aperm(out$y.rep.samples, c(3, 2, 1))
  out$y.rep.samples <- out$y.rep.samples[, , -c(K + 1)]
  out$pi.samples <- do.call(abind, lapply(out.tmp, function(a) array(a$pi.samples,
        								c(K + 1, J, n.post.samples))))
  out$pi.samples <- aperm(out$pi.samples, c(3, 2, 1))
  out$pi.samples <- out$pi.samples[, , -c(K + 1)]
  if (p.abund.re > 0) {
    out$sigma.sq.mu.samples <- mcmc(
      do.call(rbind, lapply(out.tmp, function(a) t(a$sigma.sq.mu.samples))))
    colnames(out$sigma.sq.mu.samples) <- x.re.names
    out$beta.star.samples <- mcmc(
      do.call(rbind, lapply(out.tmp, function(a) t(a$beta.star.samples))))
    tmp.names <- unlist(re.level.names)
    beta.star.names <- paste(rep(x.re.names, n.abund.re.long), tmp.names, sep = '-')
    colnames(out$beta.star.samples) <- beta.star.names
    out$re.level.names <- re.level.names
  }
  if (p.det.re > 0) {
    out$sigma.sq.p.samples <- mcmc(
      do.call(rbind, lapply(out.tmp, function(a) t(a$sigma.sq.p.samples))))
    colnames(out$sigma.sq.p.samples) <- x.p.re.names
    out$alpha.star.samples <- mcmc(
      do.call(rbind, lapply(out.tmp, function(a) t(a$alpha.star.samples))))
    tmp.names <- unlist(p.re.level.names)
    alpha.star.names <- paste(rep(x.p.re.names, n.det.re.long), tmp.names, sep = '-')
    colnames(out$alpha.star.samples) <- alpha.star.names
    out$p.re.level.names <- p.re.level.names
  }
  # Calculate effective sample sizes
  out$ESS <- list()
  out$ESS$beta <- effectiveSize(out$beta.samples)
  if (family == 'NB') {
    out$ESS$kappa <- effectiveSize(out$kappa.samples)
  }
  out$ESS$alpha <- effectiveSize(out$alpha.samples)
  if (p.det.re > 0) {
    out$ESS$sigma.sq.p <- effectiveSize(out$sigma.sq.p.samples)
  }
  if (p.abund.re > 0) {
    out$ESS$sigma.sq.mu <- effectiveSize(out$sigma.sq.mu.samples)
  }
  out$X <- X
  out$X.p <- X.p
  out$X.re <- X.re
  out$X.p.re <- X.p.re
  out$X.p.random <- X.p.random
  out$y <- y.mat
  out$offset <- offset
  out$n.samples <- n.samples
  out$call <- cl
  out$n.post <- n.post.samples
  out$n.thin <- n.thin
  out$n.burn <- n.burn
  out$n.chains <- n.chains
  out$re.cols <- re.cols
  out$re.det.cols <- re.det.cols
  out$det.func <- det.func
  out$dist.breaks <- dist.breaks
  out$dist <- family
  out$transect <- transect
  if (p.det.re > 0) {
    out$pRE <- TRUE
  } else {
    out$pRE <- FALSE
  }
  if (p.abund.re > 0) {
    out$muRE <- TRUE
  } else {
    out$muRE <- FALSE
  }
  class(out) <- "DS"
  out$run.time <- proc.time() - ptm
  out
}
