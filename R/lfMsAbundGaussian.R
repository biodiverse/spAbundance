lfMsAbundGaussian <- function(formula, data, inits, priors, tuning, n.factors,
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
    stop("error: coords must be specified in data for a spatial abundance model.")
  }
  if (!is.matrix(data$coords) & !is.data.frame(data$coords)) {
    stop("error: coords must be a matrix or data frame")
  }
  coords <- as.matrix(data$coords)
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
  if (missing(n.factors)) {
    stop("error: n.factors must be specified for a spatial factor GLMM")
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
  } else {
    w.inits <- matrix(0, q, J)
    if (verbose) {
      message("w is not specified in initial values.\nSetting initial value to 0\n")
    }
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
  # Not accessed, just kept to keep order with other abund models
  tuning.c <- 0
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

  # Set storage for all variables ---------------------------------------
  storage.mode(y) <- "double"
  storage.mode(X) <- "double"
  storage.mode(z) <- 'double'
  consts <- c(N, J, p, p.re, n.re, q, save.fitted, ind.betas)
  storage.mode(consts) <- "integer"
  storage.mode(beta.inits) <- "double"
  storage.mode(beta.comm.inits) <- "double"
  storage.mode(tau.sq.inits) <- "double"
  storage.mode(tau.sq.beta.inits) <- "double"
  storage.mode(lambda.inits) <- "double"
  storage.mode(w.inits) <- "double"
  storage.mode(mu.beta.comm) <- "double"
  storage.mode(Sigma.beta.comm) <- "double"
  storage.mode(tau.sq.beta.a) <- "double"
  storage.mode(tau.sq.beta.b) <- "double"
  storage.mode(tau.sq.a) <- "double"
  storage.mode(tau.sq.b) <- "double"
  storage.mode(tuning.c) <- "double"
  storage.mode(n.batch) <- "integer"
  storage.mode(batch.length) <- "integer"
  storage.mode(accept.rate) <- "double"
  storage.mode(n.omp.threads) <- "integer"
  storage.mode(verbose) <- "integer"
  storage.mode(n.report) <- "integer"
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
      lambda.inits <- matrix(0, N, q)
      diag(lambda.inits) <- 1
      lambda.inits[lower.tri(lambda.inits)] <- rnorm(sum(lower.tri(lambda.inits)))
      lambda.inits <- c(lambda.inits)
      tau.sq.inits <- runif(N, 0.01, 3)
      if (p.re > 0) {
        sigma.sq.mu.inits <- runif(p.re, 0.5, 10)
        beta.star.inits <- rnorm(n.re, sqrt(sigma.sq.mu.inits[beta.star.indx + 1]))
        beta.star.inits <- rep(beta.star.inits, N)
      }
    }

    storage.mode(chain.info) <- "integer"
    # Run the model in C
    out.tmp[[i]] <- .Call("lfMsAbundGaussian", y, X, X.re, X.random, consts, n.re.long,
      	            beta.inits, beta.comm.inits, tau.sq.beta.inits, tau.sq.inits,
		    lambda.inits, w.inits,
      	            sigma.sq.mu.inits, beta.star.inits,
		    beta.star.indx, beta.level.indx, mu.beta.comm, Sigma.beta.comm,
      	            tau.sq.beta.a, tau.sq.beta.b, tau.sq.a, tau.sq.b,
		    sigma.sq.mu.a, sigma.sq.mu.b, tuning.c, n.batch,
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
  if (save.fitted) {
    out$mu.samples <- do.call(abind, lapply(out.tmp, function(a) array(a$mu.samples,
      								dim = c(N, J, n.post.samples))))
    out$mu.samples <- aperm(out$mu.samples, c(3, 1, 2))
    out$like.samples <- do.call(abind, lapply(out.tmp, function(a) array(a$like.samples,
      								dim = c(N, J, n.post.samples))))
    out$like.samples <- aperm(out$like.samples, c(3, 1, 2))
    out$y.rep.samples <- do.call(abind, lapply(out.tmp, function(a) array(a$y.rep.samples,
      								dim = c(N, J, n.post.samples))))
    out$y.rep.samples <- aperm(out$y.rep.samples, c(3, 1, 2))
  }
  out$w.samples <- do.call(abind, lapply(out.tmp, function(a) array(a$w.samples,
    								dim = c(q, J, n.post.samples))))
  out$w.samples <- aperm(out$w.samples, c(3, 1, 2))

  out$X.re <- X.re
  # Calculate effective sample sizes
  out$ESS <- list()
  out$ESS$beta.comm <- effectiveSize(out$beta.comm.samples)
  out$ESS$tau.sq.beta <- effectiveSize(out$tau.sq.beta.samples)
  out$ESS$tau.sq <- effectiveSize(out$tau.sq.samples)
  out$ESS$beta <- effectiveSize(out$beta.samples)
  if (p.re > 0) {
    out$ESS$sigma.sq.mu <- effectiveSize(out$sigma.sq.mu.samples)
  }
  out$ESS$lambda <- effectiveSize(out$lambda.samples)
  out$X <- X
  out$y <- y.orig
  out$call <- cl
  out$n.samples <- n.samples
  out$x.names <- x.names
  out$sp.names <- sp.names
  out$n.post <- n.post.samples
  out$n.thin <- n.thin
  out$n.burn <- n.burn
  out$n.chains <- n.chains
  out$dist <- family
  out$coords <- coords
  out$re.cols <- re.cols
  out$q <- q
  if (p.re > 0) {
    out$muRE <- TRUE
  } else {
    out$muRE <- FALSE
  }
  class(out) <- "lfMsAbund"

  out$run.time <- proc.time() - ptm
  return(out)
}
