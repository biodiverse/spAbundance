sfMsAbundGaussian <- function(formula, data, inits, priors,
                  tuning, cov.model = 'exponential', NNGP = TRUE,
                  n.neighbors = 15, search.type = "cb", n.factors,
                  n.batch, batch.length, accept.rate = 0.43, family = 'Gaussian',
                  n.omp.threads = 1, verbose = TRUE, n.report = 100,
                  n.burn = round(.10 * n.batch * batch.length),
                  n.thin = 1, n.chains = 1, save.fitted = TRUE, ...){

  ptm <- proc.time()

  # Make it look nice
  if (verbose) {
    cat("----------------------------------------\n");
    cat("\tPreparing to run the model\n");
    cat("----------------------------------------\n");
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
  # Only implemented for NNGP
  if (!NNGP) {
    stop("error: sfMsAbundGaussian is currently only implemented for NNGPs, not full Gaussian Processes. Please set NNGP = TRUE.")
  }
  if (missing(data)) {
    stop("error: data must be specified")
  }
  if (!is.list(data)) {
    stop("error: data must be a list")
  }
  names(data) <- tolower(names(data))
  if (!'y' %in% names(data)) {
    stop("error: y must be specified in data")
  }
  if (length(dim(data$y)) != 2) {
    stop("error: data y must be a matrix or data frame with rows corresponding to species and sites corresponding to columns")
  }
  y <- as.matrix(data$y)
  sp.names <- attr(y, 'dimnames')[[1]]
  if (!'covs' %in% names(data)) {
    if (formula == ~ 1) {
      if (verbose) {
        message("covariates (covs) not specified in data.\nAssuming intercept only model.\n")
      }
      data$covs <- list(int = rep(1, dim(y)[2]))
    } else {
      stop("error: covs must be specified in data for a model with covariates")
    }
  }
  if (!is.list(data$covs)) {
    if (is.matrix(data$covs)) {
      data$covs <- data.frame(data$covs)
    } else {
      stop("error: covs must be a list, data frame, or matrix")
    }
  }
  if (!'coords' %in% names(data)) {
    stop("error: coords must be specified in data for a spatial occupancy model.")
  }
  coords <- as.matrix(data$coords)
  if (missing(n.factors)) {
    stop("error: n.factors must be specified for a spatial factor model")
  }
  # Give warning if an offset is specified
  if ('offset' %in% names(data)) {
    message("offsets are not supported with Gaussian or zi-Gaussian GLMMs.\nIgnoring data$offset when fitting the model.\n")
  }

  if (family == 'zi-Gaussian') {
    two.stage <- TRUE
  } else {
    two.stage <- FALSE
  }
  if (two.stage) {
    if (!'z' %in% names(data)) {
      stop("error: z must be specified in data for a two stage model")
    }
    z <- data$z
    if (!is.matrix(z)) {
      stop(paste0("z must be a matrix with ", nrow(y), " rows and ", ncol(y), " columns."))
    }
    if (nrow(z) != nrow(y) | ncol(z) != ncol(y)) {
      stop(paste0("z must be a matrix with ", nrow(y), " rows and ", ncol(y), " columns."))
    }
  } else {
    z <- matrix(1, nrow(y), ncol(y))
  }

  # First subset covariates to only use those that are included in the analysis.
  # Get occurrence covariates in proper format
  # Subset covariates to only use those that are included in the analysis
  data$covs <- data$covs[names(data$covs) %in% all.vars(formula)]
  # Null model support
  if (length(data$covs) == 0) {
    data$covs <- list(int = rep(1, dim(y)[2]))
  }
  # Ordered by rep, then site within rep
  data$covs <- data.frame(lapply(data$covs, function(a) unlist(c(a))))

  # Neighbors and Ordering ----------------------------------------------
  if (NNGP) {
    u.search.type <- 2
    ## Order by x column. Could potentially allow this to be user defined.
    ord <- order(coords[,1])
    # Reorder everything to align with NN ordering
    y <- y[, ord, drop = FALSE]
    z <- z[, ord, drop = FALSE]
    coords <- coords[ord, , drop = FALSE]
    # Covariates
    data$covs <- data$covs[ord, , drop = FALSE]
  }

  # Checking missing values ---------------------------------------------
  # covs ------------------------
  if (sum(is.na(data$covs)) != 0) {
    stop("error: missing values in covs. Please remove these sites from all objects in data or somehow replace the NA values with non-missing values (e.g., mean imputation).")
  }

  # Check whether random effects are sent in as numeric, and
  # return error if they are.
  # Abundance -------------------------
  if (!is.null(findbars(formula))) {
    abund.re.names <- unique(unlist(sapply(findbars(formula), all.vars)))
    for (i in 1:length(abund.re.names)) {
      if (is(data$covs[, abund.re.names[i]], 'factor')) {
        stop(paste("error: random effect variable ", abund.re.names[i], " specified as a factor. Random effect variables must be specified as numeric.", sep = ''))
      }
      if (is(data$covs[, abund.re.names[i]], 'character')) {
        stop(paste("error: random effect variable ", abund.re.names[i], " specified as character. Random effect variables must be specified as numeric.", sep = ''))
      }
    }
  }
  # Check save.fitted ---------------------------------------------------
  if (!(save.fitted %in% c(TRUE, FALSE))) {
    stop("save.fitted must be either TRUE or FALSE")
  }

  # Formula -------------------------------------------------------------
  if (missing(formula)) {
    stop("error: formula must be specified")
  }

  if (is(formula, 'formula')) {
    tmp <- parseFormula(formula, data$covs)
    X <- as.matrix(tmp[[1]])
    X.re <- as.matrix(tmp[[4]])
    x.re.names <- colnames(X.re)
    x.names <- tmp[[2]]
    X.random <- as.matrix(tmp[[5]])
    x.random.names <- colnames(X.random)
  } else {
    stop("error: formula is misspecified")
  }
  # Get RE level names
  re.level.names <- lapply(data$covs[, x.re.names, drop = FALSE],
			   function (a) sort(unique(a)))
  x.re.names <- x.random.names

  # Extract data from inputs --------------------------------------------
  # Number of species
  N <- dim(y)[1]
  # Number of latent factors
  q <- n.factors
  # Number of fixed effects
  p <- ncol(X)
  # Number of random effect parameters
  p.re <- ncol(X.re)
  # Number of latent occupancy random effect values
  n.re <- length(unlist(apply(X.re, 2, unique)))
  n.re.long <- apply(X.re, 2, function(a) length(unique(a)))
  # Number of sites
  J <- nrow(X)
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

  # y is ordered by site, then species within site.
  y.orig <- y
  y <- c(y)

  # Get random effect matrices all set ----------------------------------
  if (p.re > 1) {
    for (j in 2:p.re) {
      X.re[, j] <- X.re[, j] + max(X.re[, j - 1]) + 1
    }
  }
  # Priors --------------------------------------------------------------
  if (missing(priors)) {
    priors <- list()
  }
  names(priors) <- tolower(names(priors))

  # Independent beta parameters -----
  if ('independent.betas' %in% names(priors)) {
    if (priors$independent.betas == TRUE) {
      message("Beta parameters will be estimated independently\n")
      ind.betas <- TRUE
    } else if (priors$independent.betas == FALSE) {
      ind.betas <- FALSE
    }
  } else {
    ind.betas <- FALSE
  }
  # beta.comm -----------------------
  if ("beta.comm.normal" %in% names(priors)) {
    if (!is.list(priors$beta.comm.normal) | length(priors$beta.comm.normal) != 2) {
      stop("error: beta.comm.normal must be a list of length 2")
    }
    mu.beta.comm <- priors$beta.comm.normal[[1]]
    sigma.beta.comm <- priors$beta.comm.normal[[2]]
    if (length(mu.beta.comm) != p & length(mu.beta.comm) != 1) {
      if (p == 1) {
        stop(paste("error: beta.comm.normal[[1]] must be a vector of length ",
        	     p, " with elements corresponding to beta.comms' mean", sep = ""))
      } else {
        stop(paste("error: beta.comm.normal[[1]] must be a vector of length ",
        	     p, " or 1 with elements corresponding to beta.comms' mean", sep = ""))
      }
    }
    if (length(sigma.beta.comm) != p & length(sigma.beta.comm) != 1) {
      if (p == 1) {
        stop(paste("error: beta.comm.normal[[2]] must be a vector of length ",
      	   p, " with elements corresponding to beta.comms' variance", sep = ""))
      } else {
        stop(paste("error: beta.comm.normal[[2]] must be a vector of length ",
      	   p, " or 1 with elements corresponding to beta.comms' variance", sep = ""))
      }
    }
    if (length(sigma.beta.comm) != p) {
      sigma.beta.comm <- rep(sigma.beta.comm, p)
    }
    if (length(mu.beta.comm) != p) {
      mu.beta.comm <- rep(mu.beta.comm, p)
    }
    Sigma.beta.comm <- sigma.beta.comm * diag(p)
  } else {
    if (verbose & !ind.betas) {
      message("No prior specified for beta.comm.normal.\nSetting prior mean to 0 and prior variance to 1000\n")
    }
    mu.beta.comm <- rep(0, p)
    sigma.beta.comm <- rep(1000, p)
    Sigma.beta.comm <- diag(p) * 1000
  }

  # tau.sq.beta -----------------------
  if ("tau.sq.beta.ig" %in% names(priors)) {
    if (!is.list(priors$tau.sq.beta.ig) | length(priors$tau.sq.beta.ig) != 2) {
      stop("error: tau.sq.beta.ig must be a list of length 2")
    }
    tau.sq.beta.a <- priors$tau.sq.beta.ig[[1]]
    tau.sq.beta.b <- priors$tau.sq.beta.ig[[2]]
    if (length(tau.sq.beta.a) != p & length(tau.sq.beta.a) != 1) {
      if (p == 1) {
        stop(paste("error: tau.sq.beta.ig[[1]] must be a vector of length ",
      	   p, " with elements corresponding to tau.sq.betas' shape", sep = ""))
      } else {
        stop(paste("error: tau.sq.beta.ig[[1]] must be a vector of length ",
      	   p, " or 1 with elements corresponding to tau.sq.betas' shape", sep = ""))
      }
    }
    if (length(tau.sq.beta.b) != p & length(tau.sq.beta.b) != 1) {
      if (p == 1) {
        stop(paste("error: tau.sq.beta.ig[[2]] must be a vector of length ",
      	   p, " with elements corresponding to tau.sq.betas' scale", sep = ""))
      } else {
        stop(paste("error: tau.sq.beta.ig[[2]] must be a vector of length ",
      	   p, " or 1 with elements corresponding to tau.sq.betas' scale", sep = ""))
      }
    }
    if (length(tau.sq.beta.a) != p) {
      tau.sq.beta.a <- rep(tau.sq.beta.a, p)
    }
    if (length(tau.sq.beta.b) != p) {
      tau.sq.beta.b <- rep(tau.sq.beta.b, p)
    }
  } else {
    if (verbose & !ind.betas) {
      message("No prior specified for tau.sq.beta.ig.\nSetting prior shape to 0.1 and prior scale to 0.1\n")
    }
    tau.sq.beta.a <- rep(0.1, p)
    tau.sq.beta.b <- rep(0.1, p)
  }

  # tau.sq -----------------------
  if ("tau.sq.ig" %in% names(priors)) {
    if (!is.list(priors$tau.sq.ig) | length(priors$tau.sq.ig) != 2) {
      stop("error: tau.sq.ig must be a list of length 2")
    }
    tau.sq.a <- priors$tau.sq.ig[[1]]
    tau.sq.b <- priors$tau.sq.ig[[2]]
    if (length(tau.sq.a) != N & length(tau.sq.a) != 1) {
      stop(paste("error: tau.sq.ig[[1]] must be a vector of length ",
      	   N, " or 1 with elements corresponding to tau.sqs' shape", sep = ""))
    }
    if (length(tau.sq.b) != N & length(tau.sq.b) != 1) {
      stop(paste("error: tau.sq.ig[[2]] must be a vector of length ",
      	   p, " or 1 with elements corresponding to tau.sqs' scale", sep = ""))
    }
    if (length(tau.sq.a) != N) {
      tau.sq.a <- rep(tau.sq.a, N)
    }
    if (length(tau.sq.b) != N) {
      tau.sq.b <- rep(tau.sq.b, N)
    }
  } else {
    if (verbose) {
      message("No prior specified for tau.sq.ig.\nSetting prior shape to 0.01 and prior scale to 0.01\n")
    }
    tau.sq.a <- rep(0.01, N)
    tau.sq.b <- rep(0.01, N)
  }

    # phi -----------------------------
    # Get distance matrix which is used if priors are not specified
    if ("phi.unif" %in% names(priors)) {
      if (!is.list(priors$phi.unif) | length(priors$phi.unif) != 2) {
        stop("error: phi.unif must be a list of length 2")
      }
      phi.a <- priors$phi.unif[[1]]
      phi.b <- priors$phi.unif[[2]]
      if (length(phi.a) != q & length(phi.a) != 1) {
        stop(paste("error: phi.unif[[1]] must be a vector of length ",
        	   q, " or 1 with elements corresponding to phis' lower bound for each latent factor", sep = ""))
      }
      if (length(phi.b) != q & length(phi.b) != 1) {
        stop(paste("error: phi.unif[[2]] must be a vector of length ",
        	   q, " or 1 with elements corresponding to phis' upper bound for each latent factor", sep = ""))
      }
      if (length(phi.a) != q) {
        phi.a <- rep(phi.a, q)
      }
      if (length(phi.b) != q) {
        phi.b <- rep(phi.b, q)
      }
    } else {
      if (verbose) {
      message("No prior specified for phi.unif.\nSetting uniform bounds based on the range of observed spatial coordinates.\n")
      }
      coords.D <- iDist(coords)
      phi.a <- rep(3 / max(coords.D), q)
      phi.b <- rep(3 / sort(unique(c(coords.D)))[2], q)
    }
    # nu -----------------------------
    if (cov.model == "matern") {
      if (!"nu.unif" %in% names(priors)) {
        stop("error: nu.unif must be specified in priors value list")
      }
      nu.a <- priors$nu.unif[[1]]
      nu.b <- priors$nu.unif[[2]]
      if (!is.list(priors$nu.unif) | length(priors$nu.unif) != 2) {
        stop("error: nu.unif must be a list of length 2")
      }
      if (length(nu.a) != q & length(nu.a) != 1) {
        stop(paste("error: nu.unif[[1]] must be a vector of length ",
        	   q, " or 1 with elements corresponding to nus' lower bound for each latent factor", sep = ""))
      }
      if (length(nu.b) != q & length(nu.b) != 1) {
        stop(paste("error: nu.unif[[2]] must be a vector of length ",
        	   q, " or 1 with elements corresponding to nus' upper bound for each latent factor", sep = ""))
      }
      if (length(nu.a) != q) {
        nu.a <- rep(nu.a, q)
      }
      if (length(nu.b) != q) {
        nu.b <- rep(nu.b, q)
      }
    } else {
      nu.a <- rep(0, q)
      nu.b <- rep(0, q)
    }

  # sigma.sq.mu --------------------
  if (p.re > 0) {
    if ("sigma.sq.mu.ig" %in% names(priors)) {
      if (!is.list(priors$sigma.sq.mu.ig) | length(priors$sigma.sq.mu.ig) != 2) {
        stop("error: sigma.sq.mu.ig must be a list of length 2")
      }
      sigma.sq.mu.a <- priors$sigma.sq.mu.ig[[1]]
      sigma.sq.mu.b <- priors$sigma.sq.mu.ig[[2]]
      if (length(sigma.sq.mu.a) != p.re & length(sigma.sq.mu.a) != 1) {
        if (p.re == 1) {
        stop(paste("error: sigma.sq.mu.ig[[1]] must be a vector of length ",
        	   p.re, " with elements corresponding to sigma.sq.mus' shape", sep = ""))
        } else {
        stop(paste("error: sigma.sq.mu.ig[[1]] must be a vector of length ",
        	   p.re, " or 1 with elements corresponding to sigma.sq.mus' shape", sep = ""))
        }
      }
      if (length(sigma.sq.mu.b) != p.re & length(sigma.sq.mu.b) != 1) {
        if (p.re == 1) {
          stop(paste("error: sigma.sq.mu.ig[[2]] must be a vector of length ",
        	   p.re, " with elements corresponding to sigma.sq.mus' scale", sep = ""))
        } else {
          stop(paste("error: sigma.sq.mu.ig[[2]] must be a vector of length ",
        	   p.re, " or 1with elements corresponding to sigma.sq.mus' scale", sep = ""))
        }
      }
      if (length(sigma.sq.mu.a) != p.re) {
        sigma.sq.mu.a <- rep(sigma.sq.mu.a, p.re)
      }
      if (length(sigma.sq.mu.b) != p.re) {
        sigma.sq.mu.b <- rep(sigma.sq.mu.b, p.re)
      }
  }   else {
      if (verbose) {
        message("No prior specified for sigma.sq.mu.ig.\nSetting prior shape to 0.1 and prior scale to 0.1\n")
      }
      sigma.sq.mu.a <- rep(0.1, p.re)
      sigma.sq.mu.b <- rep(0.1, p.re)
    }
  } else {
    sigma.sq.mu.a <- 0
    sigma.sq.mu.b <- 0
  }

  # Initial values --------------------------------------------------------
  if (missing(inits)) {
    inits <- list()
  }
  names(inits) <- tolower(names(inits))
  # beta.comm -----------------------
  # ORDER: a p vector ordered by the effects in the formula.
  if ("beta.comm" %in% names(inits)) {
    beta.comm.inits <- inits[["beta.comm"]]
    if (length(beta.comm.inits) != p & length(beta.comm.inits) != 1) {
      if (p == 1) {
        stop(paste("error: initial values for beta.comm must be of length ", p,
      	   sep = ""))
      } else {
        stop(paste("error: initial values for beta.comm must be of length ", p,
      	   , " or 1", sep = ""))
      }
    }
    if (length(beta.comm.inits) != p) {
      beta.comm.inits <- rep(beta.comm.inits, p)
    }
  } else {
    beta.comm.inits <- rnorm(p, mu.beta.comm, sqrt(sigma.beta.comm))
    if (verbose) {
      message('beta.comm is not specified in initial values.\nSetting initial values to random values from the prior distribution\n')
    }
  }
  # tau.sq.beta ------------------------
  # ORDER: a p vector ordered by the effects in the occurrence formula
  if ("tau.sq.beta" %in% names(inits)) {
    tau.sq.beta.inits <- inits[["tau.sq.beta"]]
    if (length(tau.sq.beta.inits) != p & length(tau.sq.beta.inits) != 1) {
      if (p == 1) {
        stop(paste("error: initial values for tau.sq.beta must be of length ", p,
      	   sep = ""))
      } else {
        stop(paste("error: initial values for tau.sq.beta must be of length ", p,
      	   " or 1", sep = ""))
      }
    }
    if (length(tau.sq.beta.inits) != p) {
      tau.sq.beta.inits <- rep(tau.sq.beta.inits, p)
    }
  } else {
    tau.sq.beta.inits <- runif(p, 0.5, 10)
    if (verbose) {
      message('tau.sq.beta is not specified in initial values.\nSetting initial values to random values between 0.5 and 10\n')
    }
  }
  # tau.sq ------------------------
  # ORDER: a length N vector
  if ("tau.sq" %in% names(inits)) {
    tau.sq.inits <- inits[["tau.sq"]]
    if (length(tau.sq.inits) != N & length(tau.sq.inits) != 1) {
      stop(paste("error: initial values for tau.sq must be of length ", N,
      	   " or 1", sep = ""))
    }
    if (length(tau.sq.inits) != N) {
      tau.sq.inits <- rep(tau.sq.inits, N)
    }
  } else {
    tau.sq.inits <- runif(N, 0.01, 3)
    if (verbose) {
      message('tau.sq is not specified in initial values.\nSetting initial values to random values between 0.01 and 3\n')
    }
  }
  # beta ----------------------------
  # ORDER: N x p matrix sent in as a column-major vector ordered by
  #        parameter then species within parameter.
  if ("beta" %in% names(inits)) {
    beta.inits <- inits[["beta"]]
    if (is.matrix(beta.inits)) {
      if (ncol(beta.inits) != p | nrow(beta.inits) != N) {
        stop(paste("error: initial values for beta must be a matrix with dimensions ",
        	   N, "x", p, " or a single numeric value", sep = ""))
      }
    }
    if (!is.matrix(beta.inits) & length(beta.inits) != 1) {
      stop(paste("error: initial values for beta must be a matrix with dimensions ",
      	   N, " x ", p, " or a single numeric value", sep = ""))
    }
    if (length(beta.inits) == 1) {
      beta.inits <- matrix(beta.inits, N, p)
    }
  } else {
      beta.inits <- matrix(rnorm(N * p, beta.comm.inits, sqrt(tau.sq.beta.inits)), N, p)
      if (verbose) {
        message('beta is not specified in initial values.\nSetting initial values to random values from the community-level normal distribution\n')
      }
  }
  # Create a N * p x 1 matrix of the species-level regression coefficients.
  # This is ordered by parameter, then species within a parameter.
  beta.inits <- c(beta.inits)
  # phi -----------------------------
  # ORDER: a length N vector ordered by species in the detection-nondetection data.
  if ("phi" %in% names(inits)) {
    phi.inits <- inits[["phi"]]
    if (length(phi.inits) != q & length(phi.inits) != 1) {
      stop(paste("error: initial values for phi must be of length ", q, " or 1",
      	   sep = ""))
    }
    if (length(phi.inits) != q) {
      phi.inits <- rep(phi.inits, q)
    }
  } else {
    phi.inits <- runif(q, phi.a, phi.b)
    if (verbose) {
      message("phi is not specified in initial values.\nSetting initial value to random values from the prior distribution\n")
    }
  }
  # nu ------------------------
  if ("nu" %in% names(inits)) {
    nu.inits <- inits[["nu"]]
    if (length(nu.inits) != q & length(nu.inits) != 1) {
      stop(paste("error: initial values for nu must be of length ", q,  " or 1",
      	   sep = ""))
    }
    if (length(nu.inits) != q) {
      nu.inits <- rep(nu.inits, q)
    }
  } else {
    if (cov.model == 'matern') {
      if (verbose) {
        message("nu is not specified in initial values.\nSetting initial values to random values from the prior distribution\n")
      }
      nu.inits <- runif(q, nu.a, nu.b)
    } else {
      nu.inits <- rep(0, q)
    }
  }
  # lambda ----------------------------
  # ORDER: an n.sp x q matrix sent in as a column-major vector, which is ordered by
  #        factor, then species within factor.
  if ("lambda" %in% names(inits)) {
    lambda.inits <- inits[["lambda"]]
    if (!is.matrix(lambda.inits)) {
      stop(paste("error: initial values for lambda must be a matrix with dimensions ",
        	 N, " x ", q, sep = ""))
    }
    if (nrow(lambda.inits) != N | ncol(lambda.inits) != q) {
      stop(paste("error: initial values for lambda must be a matrix with dimensions ",
        	 N, " x ", q, sep = ""))
    }
    if (!all.equal(diag(lambda.inits), rep(1, q))) {
      stop("error: diagonal of inits$lambda matrix must be all 1s")
    }
    if (sum(lambda.inits[upper.tri(lambda.inits)]) != 0) {
      stop("error: upper triangle of inits$lambda must be all 0s")
    }
  } else {
    lambda.inits <- matrix(0, N, q)
    diag(lambda.inits) <- 1
    lambda.inits[lower.tri(lambda.inits)] <- 0
    if (verbose) {
      message("lambda is not specified in initial values.\nSetting initial values of the lower triangle to 0\n")
    }
    # lambda.inits are organized by factor, then by species. This is necessary for working
    # with dgemv.
    lambda.inits <- c(lambda.inits)
  }
  # w -----------------------------
  if ("w" %in% names(inits)) {
    w.inits <- inits[["w"]]
    if (!is.matrix(w.inits)) {
      stop(paste("error: initial values for w must be a matrix with dimensions ",
      	   q, " x ", J, sep = ""))
    }
    if (nrow(w.inits) != q | ncol(w.inits) != J) {
      stop(paste("error: initial values for w must be a matrix with dimensions ",
      	   q, " x ", J, sep = ""))
    }
    if (NNGP) {
      w.inits <- w.inits[, ord]
    }
  } else {
    w.inits <- matrix(0, q, J)
    if (verbose) {
      message("w is not specified in initial values.\nSetting initial value to 0\n")
    }
  }
  # sigma.sq.mu ------------------
  # ORDER: a length p.re vector ordered by the random effects in the formula.
  if (p.re > 0) {
    if ("sigma.sq.mu" %in% names(inits)) {
      sigma.sq.mu.inits <- inits[["sigma.sq.mu"]]
      if (length(sigma.sq.mu.inits) != p.re & length(sigma.sq.mu.inits) != 1) {
        if (p.re == 1) {
          stop(paste("error: initial values for sigma.sq.mu must be of length ", p.re,
      	     sep = ""))
        } else {
          stop(paste("error: initial values for sigma.sq.mu must be of length ", p.re,
      	     " or 1", sep = ""))
        }
      }
      if (length(sigma.sq.mu.inits) != p.re) {
        sigma.sq.mu.inits <- rep(sigma.sq.mu.inits, p.re)
      }
    } else {
      sigma.sq.mu.inits <- runif(p.re, 0.5, 10)
      if (verbose) {
        message("sigma.sq.mu is not specified in initial values.\nSetting initial values to random values between 0.5 and 10\n")
      }
    }
    beta.star.indx <- rep(0:(p.re - 1), n.re.long)
    beta.star.inits <- rnorm(n.re, sqrt(sigma.sq.mu.inits[beta.star.indx + 1]))
    beta.star.inits <- rep(beta.star.inits, N)
  } else {
    sigma.sq.mu.inits <- 0
    beta.star.indx <- 0
    beta.star.inits <- 0
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
  # Covariance Model ----------------------------------------------------
  # Order must match util.cpp spCor.
  cov.model.names <- c("exponential", "spherical", "matern", "gaussian")
  if(! cov.model %in% cov.model.names){
    stop("error: specified cov.model '",cov.model,"' is not a valid option; choose from ",
         paste(cov.model.names, collapse=", ", sep="") ,".")}
  # Obo for cov model lookup on c side
  cov.model.indx <- which(cov.model == cov.model.names) - 1

  # Get tuning values ---------------------------------------------------
  # Not accessed, but necessary to keep things in line with the underlying functions.
  sigma.sq.tuning <- rep(0, q)
  phi.tuning <- rep(0, q)
  nu.tuning <- rep(0, q)
  if (missing(tuning)) {
    phi.tuning <- rep(1, q)
    if (cov.model == 'matern') {
      nu.tuning <- rep(1, q)
    }
  } else {
    names(tuning) <- tolower(names(tuning))
    # phi ---------------------------
    if(!"phi" %in% names(tuning)) {
      stop("error: phi must be specified in tuning value list")
    }
    phi.tuning <- tuning$phi
    if (length(phi.tuning) == 1) {
      phi.tuning <- rep(tuning$phi, q)
    } else if (length(phi.tuning) != q) {
      stop(paste("error: phi tuning must be either a single value or a vector of length ",
      	   q, sep = ""))
    }
    if (cov.model == 'matern') {
      # nu --------------------------
      if(!"nu" %in% names(tuning)) {
        stop("error: nu must be specified in tuning value list")
      }
      nu.tuning <- tuning$nu
      if (length(nu.tuning) == 1) {
        nu.tuning <- rep(tuning$nu, q)
      } else if (length(nu.tuning) != q) {
        stop(paste("error: nu tuning must be either a single value or a vector of length ",
        	   q, sep = ""))
      }
    }
  }
  tuning.c <- log(c(sigma.sq.tuning, phi.tuning, nu.tuning))
  curr.chain <- 1

  # Other miscellaneous ---------------------------------------------------
  # For prediction with random slopes
  re.cols <- list()
  if (p.re > 0) {
    split.names <- strsplit(x.re.names, "[-]")
    for (j in 1:p.re) {
      re.cols[[j]] <- split.names[[j]][1]
      names(re.cols)[j] <- split.names[[j]][2]
    }
  }

  if (!NNGP) {

    stop("error: sfMsAbund is currently only implemented for NNGPs, not full Gaussian Processes. Please set NNGP = TRUE.")

  } else {

    # Nearest Neighbor Search ---------------------------------------------
    if(verbose){
      cat("----------------------------------------\n");
      cat("\tBuilding the neighbor list\n");
      cat("----------------------------------------\n");
    }

    search.type.names <- c("brute", "cb")

    if(!search.type %in% search.type.names){
      stop("error: specified search.type '",search.type,
	   "' is not a valid option; choose from ",
	   paste(search.type.names, collapse=", ", sep="") ,".")
    }

    ## Indexes
    if(search.type == "brute"){
      indx <- mkNNIndx(coords, n.neighbors, n.omp.threads)
    } else{
      indx <- mkNNIndxCB(coords, n.neighbors, n.omp.threads)
    }

    nn.indx <- indx$nnIndx
    nn.indx.lu <- indx$nnIndxLU
    nn.indx.run.time <- indx$run.time

    if(verbose){
      cat("----------------------------------------\n");
      cat("Building the neighbors of neighbors list\n");
      cat("----------------------------------------\n");
    }

    indx <- mkUIndx(J, n.neighbors, nn.indx, nn.indx.lu, u.search.type)

    u.indx <- indx$u.indx
    u.indx.lu <- indx$u.indx.lu
    ui.indx <- indx$ui.indx
    u.indx.run.time <- indx$run.time

    # Set storage for all variables ---------------------------------------
    storage.mode(y) <- "double"
    storage.mode(X) <- "double"
    storage.mode(z) <- 'double'
    storage.mode(coords) <- "double"
    consts <- c(N, J, p, p.re, n.re, q, ind.betas, save.fitted)
    storage.mode(consts) <- "integer"
    storage.mode(beta.inits) <- "double"
    storage.mode(beta.comm.inits) <- "double"
    storage.mode(tau.sq.inits) <- "double"
    storage.mode(tau.sq.beta.inits) <- "double"
    storage.mode(phi.inits) <- "double"
    storage.mode(lambda.inits) <- "double"
    storage.mode(nu.inits) <- "double"
    storage.mode(w.inits) <- "double"
    storage.mode(mu.beta.comm) <- "double"
    storage.mode(Sigma.beta.comm) <- "double"
    storage.mode(tau.sq.beta.a) <- "double"
    storage.mode(tau.sq.beta.b) <- "double"
    storage.mode(tau.sq.a) <- "double"
    storage.mode(tau.sq.b) <- "double"
    storage.mode(phi.a) <- "double"
    storage.mode(phi.b) <- "double"
    storage.mode(nu.a) <- "double"
    storage.mode(nu.b) <- "double"
    storage.mode(tuning.c) <- "double"
    storage.mode(n.batch) <- "integer"
    storage.mode(batch.length) <- "integer"
    storage.mode(accept.rate) <- "double"
    storage.mode(n.omp.threads) <- "integer"
    storage.mode(verbose) <- "integer"
    storage.mode(n.report) <- "integer"
    storage.mode(nn.indx) <- "integer"
    storage.mode(nn.indx.lu) <- "integer"
    storage.mode(u.indx) <- "integer"
    storage.mode(u.indx.lu) <- "integer"
    storage.mode(ui.indx) <- "integer"
    storage.mode(n.neighbors) <- "integer"
    storage.mode(cov.model.indx) <- "integer"
    # chain.info order: current chain, total number of chains
    chain.info <- c(curr.chain, n.chains)
    storage.mode(chain.info) <- "integer"
    n.post.samples <- length(seq(from = n.burn + 1,
				 to = n.samples,
				 by = as.integer(n.thin)))
    # samples.info order: burn-in, thinning rate, number of posterior samples
    samples.info <- c(n.burn, n.thin, n.post.samples)
    storage.mode(samples.info) <- "integer"
    # For random effects
    storage.mode(X.re) <- "integer"
    storage.mode(X.random) <- "double"
    beta.level.indx <- sort(unique(c(X.re)))
    storage.mode(beta.level.indx) <- "integer"
    storage.mode(sigma.sq.mu.inits) <- "double"
    storage.mode(sigma.sq.mu.a) <- "double"
    storage.mode(sigma.sq.mu.b) <- "double"
    storage.mode(beta.star.inits) <- "double"
    storage.mode(beta.star.indx) <- "integer"
    # Gaussian = 2, zi-Gaussian = 3
    family.c <- ifelse(family == 'Gaussian', 2, 3)
    storage.mode(family.c) <- 'integer'

    # Fit the model -------------------------------------------------------
    out.tmp <- list()
    # Random seed information for each chain of the model
    seeds.list <- list()
    out <- list()
    for (i in 1:n.chains) {
      # Change initial values if i > 1
      if ((i > 1) & (!fix.inits)) {
        if (!ind.betas) {
          beta.comm.inits <- rnorm(p, mu.beta.comm, sqrt(sigma.beta.comm))
          tau.sq.beta.inits <- runif(p, 0.5, 10)
	}
        beta.inits <- matrix(rnorm(N * p, beta.comm.inits,
              		     sqrt(tau.sq.beta.inits)), N, p)
        beta.inits <- c(beta.inits)
        tau.sq.inits <- runif(N, 0.01, 3)
        lambda.inits <- matrix(0, N, q)
        diag(lambda.inits) <- 1
        lambda.inits[lower.tri(lambda.inits)] <- rnorm(sum(lower.tri(lambda.inits)))
        lambda.inits <- c(lambda.inits)
        phi.inits <- runif(q, phi.a, phi.b)
        if (cov.model == 'matern') {
          nu.inits <- runif(q, nu.a, nu.b)
        }
        if (p.re > 0) {
          sigma.sq.mu.inits <- runif(p.re, 0.5, 10)
          beta.star.inits <- rnorm(n.re, sqrt(sigma.sq.mu.inits[beta.star.indx + 1]))
          beta.star.inits <- rep(beta.star.inits, N)
        }
      }

      storage.mode(chain.info) <- "integer"
      # Run the model in C
      out.tmp[[i]] <- .Call("sfMsAbundGaussianNNGP", y, X, coords, X.re,
			    X.random, consts, n.re.long,
        	            n.neighbors, nn.indx, nn.indx.lu, u.indx, u.indx.lu, ui.indx,
        	            beta.inits, beta.comm.inits, tau.sq.beta.inits, tau.sq.inits,
        	            phi.inits, lambda.inits, nu.inits, w.inits, sigma.sq.mu.inits,
        		    beta.star.inits, beta.star.indx, beta.level.indx, mu.beta.comm,
        	            Sigma.beta.comm,
        	            tau.sq.beta.a, tau.sq.beta.b, tau.sq.a, tau.sq.b, phi.a, phi.b,
        	            nu.a, nu.b, sigma.sq.mu.a, sigma.sq.mu.b,
        		    tuning.c, cov.model.indx, n.batch,
        	            batch.length, accept.rate, n.omp.threads, verbose, n.report,
        	            samples.info, chain.info, z, family.c)
      chain.info[1] <- chain.info[1] + 1
    }
    # Calculate R-Hat ---------------
    out$rhat <- list()
    if (n.chains > 1) {
      # as.vector removes the "Upper CI" when there is only 1 variable.
      if (!ind.betas) {
        out$rhat$beta.comm <- as.vector(gelman.diag(mcmc.list(lapply(out.tmp, function(a)
        					      mcmc(t(a$beta.comm.samples)))),
        			     autoburnin = FALSE, multivariate = FALSE)$psrf[, 2])
        out$rhat$tau.sq.beta <- as.vector(gelman.diag(mcmc.list(lapply(out.tmp, function(a)
        					      mcmc(t(a$tau.sq.beta.samples)))),
        			     autoburnin = FALSE, multivariate = FALSE)$psrf[, 2])
      } else {
        out$rhat$beta.comm <- rep(NA, p)
        out$rhat$tau.sq.beta <- rep(NA, p)
      }
      out$rhat$tau.sq <- as.vector(gelman.diag(mcmc.list(lapply(out.tmp, function(a)
      					      mcmc(t(a$tau.sq.samples)))),
      			     autoburnin = FALSE, multivariate = FALSE)$psrf[, 2])
      out$rhat$beta <- as.vector(gelman.diag(mcmc.list(lapply(out.tmp, function(a)
      					         mcmc(t(a$beta.samples)))),
      			     autoburnin = FALSE, multivariate = FALSE)$psrf[, 2])
      out$rhat$theta <- gelman.diag(mcmc.list(lapply(out.tmp, function(a)
      					      mcmc(t(a$theta.samples)))),
      			      autoburnin = FALSE, multivariate = FALSE)$psrf[, 2]
      lambda.mat <- matrix(lambda.inits, N, q)
      out$rhat$lambda.lower.tri <- as.vector(gelman.diag(mcmc.list(lapply(out.tmp, function(a)
        					       mcmc(t(a$lambda.samples[c(lower.tri(lambda.mat)), ])))),
        					       autoburnin = FALSE, multivariate = FALSE)$psrf[, 2])
      if (p.re > 0) {
        out$rhat$sigma.sq.mu <- as.vector(gelman.diag(mcmc.list(lapply(out.tmp, function(a)
        					      mcmc(t(a$sigma.sq.mu.samples)))),
        			     autoburnin = FALSE, multivariate = FALSE)$psrf[, 2])
      }
    } else {
      out$rhat$beta.comm <- rep(NA, p)
      out$rhat$tau.sq.beta <- rep(NA, p)
      out$rhat$tau.sq <- rep(NA, N)
      out$rhat$beta <- rep(NA, p * N)
      out$rhat$theta <- rep(NA, ifelse(cov.model == 'matern', 2 * q, q))
      if (p.re > 0) {
        out$rhat$sigma.sq.mu <- rep(NA, p.re)
      }
    }

    # Put everything into MCMC objects
    out$beta.comm.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$beta.comm.samples))))
    colnames(out$beta.comm.samples) <- x.names
    out$tau.sq.beta.samples <- mcmc(do.call(rbind,
      				lapply(out.tmp, function(a) t(a$tau.sq.beta.samples))))
    colnames(out$tau.sq.beta.samples) <- x.names

    if (is.null(sp.names)) {
      sp.names <- paste('sp', 1:N, sep = '')
    }
    coef.names <- paste(rep(x.names, each = N), sp.names, sep = '-')
    out$beta.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$beta.samples))))
    colnames(out$beta.samples) <- coef.names
    out$tau.sq.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$tau.sq.samples))))
    colnames(out$tau.sq.samples) <- sp.names
    if (p.re > 0) {
      out$sigma.sq.mu.samples <- mcmc(
        do.call(rbind, lapply(out.tmp, function(a) t(a$sigma.sq.mu.samples))))
      colnames(out$sigma.sq.mu.samples) <- x.re.names
      out$beta.star.samples <- mcmc(
        do.call(rbind, lapply(out.tmp, function(a) t(a$beta.star.samples))))
      tmp.names <- unlist(re.level.names)
      beta.star.names <- paste(rep(x.re.names, n.re.long), tmp.names, sep = '-')
      beta.star.names <- paste(beta.star.names, rep(sp.names, each = n.re), sep = '-')
      colnames(out$beta.star.samples) <- beta.star.names
      out$re.level.names <- re.level.names
    }
    loadings.names <- paste(rep(sp.names, times = n.factors), rep(1:n.factors, each = N), sep = '-')
    out$lambda.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$lambda.samples))))
    colnames(out$lambda.samples) <- loadings.names
    out$theta.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$theta.samples))))
    if (cov.model != 'matern') {
      theta.names <- paste(rep(c('phi'), each = q), 1:q, sep = '-')
    } else {
      theta.names <- paste(rep(c('phi', 'nu'), each = q), 1:q, sep = '-')
    }
    colnames(out$theta.samples) <- theta.names

    out$w.samples <- do.call(abind, lapply(out.tmp, function(a) array(a$w.samples,
      								dim = c(q, J, n.post.samples))))
    out$w.samples <- out$w.samples[, order(ord), , drop = FALSE]
    out$w.samples <- aperm(out$w.samples, c(3, 1, 2))
    if (save.fitted) {
      out$mu.samples <- do.call(abind, lapply(out.tmp, function(a) array(a$mu.samples,
        								dim = c(N, J, n.post.samples))))
      out$mu.samples <- out$mu.samples[, order(ord), ]
      out$mu.samples <- aperm(out$mu.samples, c(3, 1, 2))
      out$like.samples <- do.call(abind, lapply(out.tmp, function(a) array(a$like.samples,
        								dim = c(N, J, n.post.samples))))
      out$like.samples <- out$like.samples[, order(ord), ]
      out$like.samples <- aperm(out$like.samples, c(3, 1, 2))
      out$y.rep.samples <- do.call(abind, lapply(out.tmp, function(a) array(a$y.rep.samples,
        								dim = c(N, J, n.post.samples))))
      out$y.rep.samples <- out$y.rep.samples[, order(ord), ]
      out$y.rep.samples <- aperm(out$y.rep.samples, c(3, 1, 2))
    }

    out$X.re <- X.re[order(ord), , drop = FALSE]
    # Calculate effective sample sizes
    out$ESS <- list()
    out$ESS$beta.comm <- effectiveSize(out$beta.comm.samples)
    out$ESS$tau.sq.beta <- effectiveSize(out$tau.sq.beta.samples)
    out$ESS$tau.sq <- effectiveSize(out$tau.sq.samples)
    out$ESS$beta <- effectiveSize(out$beta.samples)
    out$ESS$theta <- effectiveSize(out$theta.samples)
    out$ESS$lambda <- effectiveSize(out$lambda.samples)
    if (p.re > 0) {
      out$ESS$sigma.sq.mu <- effectiveSize(out$sigma.sq.mu.samples)
    }
    out$X <- X[order(ord), , drop = FALSE]
    out$y <- y.orig[, order(ord), drop = FALSE]
    out$call <- cl
    out$n.samples <- n.samples
    out$x.names <- x.names
    out$sp.names <- sp.names
    out$theta.names <- theta.names
    out$type <- "NNGP"
    out$coords <- coords[order(ord), ]
    out$cov.model.indx <- cov.model.indx
    out$n.neighbors <- n.neighbors
    out$q <- q
    out$n.post <- n.post.samples
    out$n.thin <- n.thin
    out$n.burn <- n.burn
    out$n.chains <- n.chains
    out$dist <- family
    out$re.cols <- re.cols
    if (p.re > 0) {
      out$muRE <- TRUE
    } else {
      out$muRE <- FALSE
    }
    class(out) <- "sfMsAbund"
  }

  out$run.time <- proc.time() - ptm
  return(out)
}
