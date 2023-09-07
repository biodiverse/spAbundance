svcAbund <- function(formula, data, inits, priors, tuning, 
                    svc.cols = 1, cov.model = 'exponential', NNGP = TRUE, 
		    n.neighbors = 15,  search.type = 'cb', n.batch, 
		    batch.length, accept.rate = 0.43, family = 'Gaussian',
		    n.omp.threads = 1, verbose = TRUE, n.report = 100, 
		    n.burn = round(.10 * n.batch * batch.length), 
		    n.thin = 1, n.chains = 1, ...){

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
  if (missing(formula)) {
    stop("error: formula must be specified")
  }
  if (!'y' %in% names(data)) {
    stop("error: detection-nondetection data y must be specified in data")
  }
  y <- c(data$y)
  if (!'covs' %in% names(data)) {
    if (formula == ~ 1) {
      if (verbose) {
        message("covariates (covs) not specified in data.\nAssuming intercept only model.\n")
      }
      data$covs <- matrix(1, length(y), 1)
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
    stop("error: coords must be specified in data for a spatial model.")
  }
  if (!is.matrix(data$coords) & !is.data.frame(data$coords)) {
    stop("error: coords must be a matrix or data frame")
  }
  coords <- as.matrix(data$coords)

  if (!(family) %in% c('Gaussian', 'zi-Gaussian')) {
    stop("svcAbund currently only supports family = 'Gaussian' or 'zi-Gaussian'")
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
  } else {
    z <- rep(1, length(y))
  }

  # First subset covariates to only use those that are included in the analysis. 
  # Get occurrence covariates in proper format
  # Subset covariates to only use those that are included in the analysis
  data$covs <- data$covs[names(data$covs) %in% all.vars(formula)]
  # Null model support
  if (length(data$covs) == 0) {
    data$covs <- list(int = rep(1, length(y)))
  }
  # Ordered by rep, then site within rep
  data$covs <- data.frame(lapply(data$covs, function(a) unlist(c(a))))
  
  # Check first-stage sample ----------------------------------------------
  if (length(z) != length(y)) {
    stop(paste("z must be a vector of length ", length(y), ".", sep = ''))
  }
  # Number of points where z == 1
  J.est <- sum(z == 1)
  # Number of points where z != 1
  J.zero <- sum(z == 0)
  # Index for the sites where z == 1
  z.indx <- which(z == 1)

  # Filter all objects to only use sites with z == 1
  y.orig <- y
  y <- y[z.indx]
  coords <- coords[z.indx, ]
  data$covs <- data$covs[z.indx, , drop = FALSE]

  # Neighbors and Ordering ----------------------------------------------
  if (NNGP) {
    u.search.type <- 2 
    ## Order by x column. Could potentially allow this to be user defined. 
    ord <- order(coords[,1]) 
    # Reorder everything to align with NN ordering
    y <- y[ord, drop = FALSE]
    coords <- coords[ord, , drop = FALSE]
    # Occupancy covariates
    data$covs <- data$covs[ord, , drop = FALSE]
  }

  data$covs <- as.data.frame(data$covs)
  
  # Checking missing values ---------------------------------------------
  # y -------------------------------
  if (sum(is.na(y) > 0)) {
    stop("error: some sites in y have missing values. Remove these sites from all objects in the 'data' argument, then use 'predict' to obtain predictions at these locations if desired.")
  }
  # covs ------------------------
  if (sum(is.na(data$covs)) != 0) {
    stop("error: missing values in covs. Please remove these sites from all objects in data or somehow replace the NA values with non-missing values (e.g., mean imputation).") 
  }

  # Check whether random effects are sent in as numeric, and
  # return error if they are. 
  # Occurrence ----------------------
  if (!is.null(findbars(formula))) {
    re.names <- unique(unlist(sapply(findbars(formula), all.vars)))
    for (i in 1:length(re.names)) {
      if (is(data$covs[, re.names[i]], 'factor')) {
        stop(paste("error: random effect variable ", re.names[i], " specified as a factor. Random effect variables must be specified as numeric.", sep = ''))
      } 
      if (is(data$covs[, re.names[i]], 'character')) {
        stop(paste("error: random effect variable ", re.names[i], " specified as character. Random effect variables must be specified as numeric.", sep = ''))
      }
    }
  }

  # Formula -------------------------------------------------------------
  # Occupancy -----------------------
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

  # Get basic info from inputs ------------------------------------------
  # Number of sites
  J <- nrow(coords)
  # Number of parameters 
  p <- ncol(X)
  # Number of random effect parameters
  p.re <- ncol(X.re)
  # Number of latent random effect values
  n.re <- length(unlist(apply(X.re, 2, unique)))
  n.re.long <- apply(X.re, 2, function(a) length(unique(a)))
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
  n.post.samples <- length(seq(from = n.burn + 1, 
                               to = n.samples, 
                               by = as.integer(n.thin)))

  # Check SVC columns -----------------------------------------------------
  if (is.character(svc.cols)) {
    # Check if all column names in svc are in covs
    if (!all(svc.cols %in% x.names)) {
        missing.cols <- svc.cols[!(svc.cols %in% x.names)]
        stop(paste("error: variable name ", paste(missing.cols, collapse=" and "), " not inurrence covariates", sep=""))
    }
    # Convert desired column names into the numeric column index
    svc.cols <- (1:p)[x.names %in% svc.cols]
    
  } else if (is.numeric(svc.cols)) {
    # Check if all column indices are in 1:p
    if (!all(svc.cols %in% 1:p)) {
        missing.cols <- svc.cols[!(svc.cols %in% (1:p))]
        stop(paste("error: column index ", paste(missing.cols, collapse=" "), " not in design matrix columns", sep=""))
    }
  }
  p.svc <- length(svc.cols)

  # Get random effect matrices all set ----------------------------------
  X.re <- X.re - 1
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
  # beta -----------------------
  if ("beta.normal" %in% names(priors)) {
    if (!is.list(priors$beta.normal) | length(priors$beta.normal) != 2) {
      stop("error: beta.normal must be a list of length 2")
    }
    mu.beta <- priors$beta.normal[[1]]
    sigma.beta <- priors$beta.normal[[2]]
    if (length(mu.beta) != p & length(mu.beta) != 1) {
      if (p == 1) {
        stop(paste("error: beta.normal[[1]] must be a vector of length ",
        	     p, " with elements corresponding to betas' mean", sep = ""))
      } else {
        stop(paste("error: beta.normal[[1]] must be a vector of length ",
        	     p, " or 1 with elements corresponding to betas' mean", sep = ""))
      }
    }
    if (length(sigma.beta) != p & length(sigma.beta) != 1) {
      if (p == 1) {
        stop(paste("error: beta.normal[[2]] must be a vector of length ",
      	   p, " with elements corresponding to betas' variance", sep = ""))
      } else {
        stop(paste("error: beta.normal[[2]] must be a vector of length ",
      	   p, " or 1 with elements corresponding to betas' variance", sep = ""))
      }
    }
    if (length(sigma.beta) != p) {
      sigma.beta <- rep(sigma.beta, p)
    }
    if (length(mu.beta) != p) {
      mu.beta <- rep(mu.beta, p)
    }
    Sigma.beta <- sigma.beta * diag(p)
  } else {
    if (verbose) {
      message("No prior specified for beta.normal.\nSetting prior mean to 0 and prior variance to 100\n")
    }
    mu.beta <- rep(0, p)
    sigma.beta <- rep(100, p)
    Sigma.beta <- diag(p) * 100
  }
  # phi -----------------------------
  # Get distance matrix which is used if priors are not specified
  if ("phi.unif" %in% names(priors)) {
    if (!is.list(priors$phi.unif) | length(priors$phi.unif) != 2) {
      stop("error: phi.unif must be a list of length 2")
    }
    phi.a <- priors$phi.unif[[1]]
    phi.b <- priors$phi.unif[[2]]
    if (length(phi.a) != p.svc & length(phi.a) != 1) {
      stop(paste("error: phi.unif[[1]] must be a vector of length ", 
      	   p.svc, 
           " or 1 with elements corresponding to phis' lower bound for each covariate with spatially-varying coefficients",
           sep = ""))
    }
    if (length(phi.b) != p.svc & length(phi.b) != 1) {
      stop(paste("error: phi.unif[[2]] must be a vector of length ", 
      	   p.svc, 
           " or 1 with elements corresponding to phis' upper bound for each covariate with spatially-varying coefficients", sep = ""))
    }
    if (length(phi.a) != p.svc) {
      phi.a <- rep(phi.a, p.svc)
    }
    if (length(phi.b) != p.svc) {
      phi.b <- rep(phi.b, p.svc)
    }
  } else {
    if (verbose) {
    message("No prior specified for phi.unif.\nSetting uniform bounds based on the range of observed spatial coordinates.\n")
    }
    coords.D <- iDist(coords)
    phi.a <- rep(3 / max(coords.D), p.svc)
    phi.b <- rep(3 / sort(unique(c(coords.D)))[2], p.svc)
  }
  # tau.sq.t ----------------------
  if ("tau.sq.ig" %in% names(priors)) { 
    if (!is.vector(priors$tau.sq.ig) | !is.atomic(priors$tau.sq.ig) | length(priors$tau.sq.ig) != 2) {
      stop("error: tau.sq.ig must be a vector of length 2 with elements corresponding to tau.sq's shape and scale parameters")
    }
    tau.sq.a <- priors$tau.sq.ig[1]
    tau.sq.b <- priors$tau.sq.ig[2]
  } else {
    if (verbose) {
      message("No prior specified for tau.sq.\nUsing an inverse-Gamma prior with the shape parameter set to 2 and scale parameter to 0.5.\n")
    }
    tau.sq.a <- 2
    tau.sq.b <- 0.5
  }
  # sigma.sq -----------------------------
  if ("sigma.sq.ig" %in% names(priors)) {
    if (!is.list(priors$sigma.sq.ig) | length(priors$sigma.sq.ig) != 2) {
      stop("error: sigma.sq.ig must be a list of length 2")
    }
    sigma.sq.a <- priors$sigma.sq.ig[[1]]
    sigma.sq.b <- priors$sigma.sq.ig[[2]]
    if (length(sigma.sq.a) != p.svc & length(sigma.sq.a) != 1) {
      stop(paste("error: sigma.sq.ig[[1]] must be a vector of length ", 
      	   p.svc, " or 1 with elements corresponding to sigma.sqs' shape for each covariate with spatially-varying coefficients", sep = ""))
    }
    if (length(sigma.sq.b) != p.svc & length(sigma.sq.b) != 1) {
      stop(paste("error: sigma.sq.ig[[2]] must be a vector of length ", 
      	   p.svc, " or 1 with elements corresponding to sigma.sqs' scale for each covariate with spatially-varying coefficients", sep = ""))
    }
    if (length(sigma.sq.a) != p.svc) {
      sigma.sq.a <- rep(sigma.sq.a, p.svc)
    }
    if (length(sigma.sq.b) != p.svc) {
      sigma.sq.b <- rep(sigma.sq.b, p.svc)
    }
  } else {
    if (verbose) {
      message("No prior specified for sigma.sq.ig.\nSetting the shape parameter to 2 and scale parameter to 1.\n")
    }
    sigma.sq.a <- rep(2, p.svc)
    sigma.sq.b <- rep(1, p.svc)
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
    if (length(nu.a) != p.svc & length(nu.a) != 1) {
      stop(paste("error: nu.unif[[1]] must be a vector of length ", 
      	   p.svc, " or 1 with elements corresponding to nus' lower bound for each covariate with spatially-varying coefficients", sep = ""))
    }
    if (length(nu.b) != p.svc & length(nu.b) != 1) {
      stop(paste("error: nu.unif[[2]] must be a vector of length ", 
      	   p.svc, " or 1 with elements corresponding to nus' upper bound for each covariate with spatially-varying coefficients", sep = ""))
    }
    if (length(nu.a) != p.svc) {
      nu.a <- rep(nu.a, p.svc)
    }
    if (length(nu.b) != p.svc) {
      nu.b <- rep(nu.b, p.svc)
    }
  } else {
    nu.a <- rep(0, p.svc)
    nu.b <- rep(0, p.svc)
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

  # Starting values -----------------------------------------------------
  if (missing(inits)) {
    inits <- list()
  }
  names(inits) <- tolower(names(inits))
  # beta -----------------------
  if ("beta" %in% names(inits)) {
    beta.inits <- inits[["beta"]]
    if (length(beta.inits) != p & length(beta.inits) != 1) {
      if (p == 1) {
        stop(paste("error: initial values for beta must be of length ", p,
      	     sep = ""))

      } else {
        stop(paste("error: initial values for beta must be of length ", p, " or 1",
        	     sep = ""))
      }
    }
    if (length(beta.inits) != p) {
      beta.inits <- rep(beta.inits, p)
    }
  } else {
    beta.inits <- rnorm(p, mu.beta, sqrt(sigma.beta))
    if (verbose) {
      message('beta is not specified in initial values.\nSetting initial values to random values from the prior distribution\n')
    }
  }
  # phi -----------------------------
  if ("phi" %in% names(inits)) {
    phi.inits <- inits[["phi"]]
    if (length(phi.inits) != p.svc & length(phi.inits) != 1) {
      stop(paste("error: initial values for phi must be of length ", p.svc, " or 1", 
      	   sep = ""))
    }
    if (length(phi.inits) != p.svc) {
      phi.inits <- rep(phi.inits, p.svc)
    }
  } else {
    phi.inits <- runif(p.svc, phi.a, phi.b)
    if (verbose) {
      message("phi is not specified in initial values.\nSetting initial value to random values from the prior distribution\n")
    }
  }
  # tau.sq ------------------------
  if ("tau.sq" %in% names(inits)) {
    tau.sq.inits <- inits[["tau.sq"]]
    if (length(tau.sq.inits) != 1) {
      stop("error: initial values for tau.sq must be of length 1")
    }
  } else {
    tau.sq.inits <- runif(1, 0.5, 10)
    if (verbose) {
      message("tau.sq is not specified in initial values.\nSetting initial value to random value between 0.5 and 10\n")
    }
  }
  # sigma.sq ------------------------
  if ("sigma.sq" %in% names(inits)) {
    sigma.sq.inits <- inits[["sigma.sq"]]
    if (length(sigma.sq.inits) != p.svc & length(sigma.sq.inits) != 1) {
      stop(paste("error: initial values for sigma.sq must be of length ", p.svc,  " or 1",
      	   sep = ""))
    }
    if (length(sigma.sq.inits) != p.svc) {
      sigma.sq.inits <- rep(sigma.sq.inits, p.svc)
    }
  } else {
    sigma.sq.inits <- rigamma(p.svc, sigma.sq.a, sigma.sq.b)
    if (verbose) {
      message("sigma.sq is not specified in initial values.\nSetting initial values to random values from the prior distribution\n")
    }
  }
  # w -----------------------------
  if ("w" %in% names(inits)) {
    w.inits <- inits[["w"]]
    if (!is.matrix(w.inits)) {
      stop(paste("error: initial values for w must be a matrix with dimensions ",
      	   p.svc, " x ", J.est, sep = ""))
    }
    if (nrow(w.inits) != p.svc | ncol(w.inits) != J.est) {
      stop(paste("error: initial values for w must be a matrix with dimensions ",
      	   p.svc, " x ", J.est, sep = ""))
    }
    if (NNGP) {
      w.inits <- w.inits[, ord, drop = FALSE]
    }
  } else {
    w.inits <- matrix(0, p.svc, J.est)
    if (verbose) {
      message("w is not specified in initial values.\nSetting initial value to 0\n")
    }
  }
  # nu ------------------------
  if ("nu" %in% names(inits)) {
    nu.inits <- inits[["nu"]]
    if (length(nu.inits) != p.svc & length(nu.inits) != 1) {
      stop(paste("error: initial values for nu must be of length ", p.svc,  " or 1",
      	   sep = ""))
    }
    if (length(nu.inits) != p.svc) {
      nu.inits <- rep(nu.inits, p.svc)
    }
  } else {
    if (cov.model == 'matern') {
      if (verbose) {
        message("nu is not specified in initial values.\nSetting initial values to random values from the prior distribution\n")
      }
      nu.inits <- runif(p.svc, nu.a, nu.b)
    } else {
      nu.inits <- rep(0, p.svc)
    }
  }
  # sigma.sq.mu -------------------
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
  storage.mode(cov.model.indx) <- "integer"

  # Prep for SVCs ---------------------------------------------------------
  X.w <- X[, svc.cols, drop = FALSE]

  # Get tuning values ---------------------------------------------------
  # Not accessed, but necessary to keep things in line. 
  sigma.sq.tuning <- rep(0, p.svc)
  phi.tuning <- rep(0, p.svc)
  nu.tuning <- rep(0, p.svc)
  if (missing(tuning)) {
    phi.tuning <- rep(1, p.svc)
    if (cov.model == 'matern') {
      nu.tuning <- rep(1, p.svc)
    }
  } else {
    names(tuning) <- tolower(names(tuning))
    # phi ---------------------------
    if(!"phi" %in% names(tuning)) {
      stop("error: phi must be specified in tuning value list")
    }
    phi.tuning <- tuning$phi
    if (length(phi.tuning) == 1) {
      phi.tuning <- rep(tuning$phi, p.svc)
    } else if (length(phi.tuning) != p.svc) {
      stop(paste("error: phi tuning must be either a single value or a vector of length ",
      	   p.svc, sep = ""))
    }
    if (cov.model == 'matern') {
      # nu --------------------------
      if(!"nu" %in% names(tuning)) {
        stop("error: nu must be specified in tuning value list")
      }
      nu.tuning <- tuning$nu
      if (length(nu.tuning) == 1) {
        nu.tuning <- rep(tuning$nu, p.svc)
      } else if (length(nu.tuning) != p.svc) {
        stop(paste("error: nu tuning must be either a single value or a vector of length ",
        	   p.svc, sep = ""))
      }
    }
  }
  tuning.c <- log(c(sigma.sq.tuning, phi.tuning, nu.tuning))
  # Set model.deviance to NA for returning when no cross-validation
  model.deviance <- NA
  curr.chain <- 1

  if (!NNGP) {
    stop("error: svcAbund is currently only implemented for NNGPs, not full Gaussian Processes. Please set NNGP = TRUE.") 

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
    
    storage.mode(n.neighbors) <- "integer"
    storage.mode(n.omp.threads) <- "integer"
    ## Indexes
    if(search.type == "brute"){
      indx <- mkNNIndx(coords, n.neighbors, n.omp.threads)
    } else{
      indx <- mkNNIndxCB(coords, n.neighbors, n.omp.threads)
    }
    
    nn.indx <- indx$nnIndx
    nn.indx.lu <- indx$nnIndxLU
    nn.indx.run.time <- indx$run.time
    
    storage.mode(nn.indx) <- "integer"
    storage.mode(nn.indx.lu) <- "integer"
    storage.mode(u.search.type) <- "integer"
    storage.mode(J.est) <- "integer"

    if(verbose){
      cat("----------------------------------------\n");
      cat("Building the neighbors of neighbors list\n");
      cat("----------------------------------------\n");
    }
    
    indx <- mkUIndx(J.est, n.neighbors, nn.indx, nn.indx.lu, u.search.type)
    
    u.indx <- indx$u.indx
    u.indx.lu <- indx$u.indx.lu
    ui.indx <- indx$ui.indx
    u.indx.run.time <- indx$run.time

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
    storage.mode(X.w) <- "double"
    consts <- c(J.est, p, p.re, n.re, p.svc, J.zero)
    storage.mode(consts) <- "integer"
    storage.mode(coords) <- "double"
    storage.mode(beta.inits) <- "double"
    storage.mode(tau.sq.inits) <- "double"
    storage.mode(phi.inits) <- "double"
    storage.mode(sigma.sq.inits) <- "double"
    storage.mode(nu.inits) <- "double"
    storage.mode(w.inits) <- "double"
    storage.mode(mu.beta) <- "double"
    storage.mode(Sigma.beta) <- "double"
    storage.mode(phi.a) <- "double"
    storage.mode(phi.b) <- "double"
    storage.mode(tau.sq.a) <- "double"
    storage.mode(tau.sq.b) <- "double"
    storage.mode(nu.a) <- "double"
    storage.mode(nu.b) <- "double"
    storage.mode(sigma.sq.a) <- "double"
    storage.mode(sigma.sq.b) <- "double"
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
    chain.info <- c(curr.chain, n.chains)
    storage.mode(chain.info) <- "integer"
    n.post.samples <- length(seq(from = n.burn + 1, 
				 to = n.samples, 
				 by = as.integer(n.thin)))
    storage.mode(n.post.samples) <- "integer"
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
    if (two.stage) {
      storage.mode(z) <- 'double'
    }

    # Fit the model -------------------------------------------------------
    out.tmp <- list()
    # Random seed information for each chain of the model
    seeds.list <- list()
    out <- list()
    for (i in 1:n.chains) {
      # Change initial values if i > 1
      if ((i > 1) & (!fix.inits)) {
        beta.inits <- rnorm(p, mu.beta, sqrt(sigma.beta))
        sigma.sq.inits <- rigamma(p.svc, sigma.sq.a, sigma.sq.b)
	sigma.sq.inits <- runif(p.svc, 0.05, 10)
        phi.inits <- runif(p.svc, phi.a, phi.b)
        if (cov.model == 'matern') {
          nu.inits <- runif(p.svc, nu.a, nu.b)
        }
        if (p.re > 0) {
          sigma.sq.mu.inits <- runif(p.re, 0.5, 10)
          beta.star.inits <- rnorm(n.re, sqrt(sigma.sq.mu.inits[beta.star.indx + 1]))
        }
        tau.sq.inits <- runif(1, 0.1, 10)
      }
      storage.mode(chain.info) <- "integer"
      # Run the model in C    
      out.tmp[[i]] <- .Call("svcAbundNNGP", y, X, 
      		            X.w, coords, 
      		            X.re, X.random, 
      		            consts, n.re.long, 
        	             n.neighbors, nn.indx, nn.indx.lu, u.indx, u.indx.lu, ui.indx, 
                            beta.inits, tau.sq.inits, sigma.sq.mu.inits, beta.star.inits, 
                            w.inits, phi.inits, sigma.sq.inits, nu.inits, 
                            beta.star.indx, beta.level.indx, mu.beta, 
                            Sigma.beta, tau.sq.a, tau.sq.b, phi.a, phi.b, 
                            sigma.sq.a, sigma.sq.b, nu.a, nu.b, 
      		            sigma.sq.mu.a, sigma.sq.mu.b, 
                            tuning.c, cov.model.indx,
                            n.batch, batch.length, 
                            accept.rate, n.omp.threads, verbose, n.report, 
                            samples.info, chain.info)
      chain.info[1] <- chain.info[1] + 1
      seeds.list[[i]] <- .Random.seed
    }
    # Calculate R-Hat ---------------
    out$rhat <- list()
    if (n.chains > 1) {
      # as.vector removes the "Upper CI" when there is only 1 variable. 
      out$rhat$beta <- as.vector(gelman.diag(mcmc.list(lapply(out.tmp, function(a) 
      					         mcmc(t(a$beta.samples)))), 
      			     autoburnin = FALSE)$psrf[, 2])
      out$rhat$tau.sq <- as.vector(gelman.diag(mcmc.list(lapply(out.tmp, function(a) 
      					      mcmc(t(a$tau.sq.samples)))), 
      			     autoburnin = FALSE)$psrf[, 2])
      out$rhat$theta <- gelman.diag(mcmc.list(lapply(out.tmp, function(a) 
      					        mcmc(t(a$theta.samples)))), 
      			      autoburnin = FALSE)$psrf[, 2]
      if (p.re > 0) {
        out$rhat$sigma.sq.mu <- as.vector(gelman.diag(mcmc.list(lapply(out.tmp, function(a) 
        					      mcmc(t(a$sigma.sq.mu.samples)))), 
        			     autoburnin = FALSE)$psrf[, 2])
      }
    } else {
      out$rhat$beta <- rep(NA, p)
      out$rhat$tau.sq <- NA
      out$rhat$theta <- rep(NA, ifelse(cov.model == 'matern', 3 * p.svc, 2 * p.svc))
      if (p.re > 0) {
        out$rhat$sigma.sq.mu <- rep(NA, p.re)
      }
    }
    out$beta.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$beta.samples))))
    colnames(out$beta.samples) <- x.names
    out$theta.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$theta.samples))))
    if (cov.model != 'matern') {
      theta.names <- paste(rep(c('sigma.sq', 'phi'), each = p.svc), x.names[svc.cols], sep = '-')
    } else {
      theta.names <- paste(rep(c('sigma.sq', 'phi', 'nu'), each = p.svc), x.names[svc.cols], sep = '-')
    } 
    colnames(out$theta.samples) <- theta.names
    out$tau.sq.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$tau.sq.samples))))
    # Get everything back in the original order
    out$coords <- coords[order(ord), ]
    if (!two.stage) {
      out$y.rep.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$y.rep.samples))))
      out$y.rep.samples <- mcmc(out$y.rep.samples[, order(ord), drop = FALSE])
      out$like.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$like.samples))))
      out$like.samples <- mcmc(out$like.samples[, order(ord), drop = FALSE])
    } else {
      y.rep.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$y.rep.samples))))
      y.rep.samples <- mcmc(y.rep.samples[, order(ord), drop = FALSE])
      like.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$like.samples))))
      like.samples <- mcmc(like.samples[, order(ord), drop = FALSE])
      y.rep.zero.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$y.rep.zero.samples))))
      out$y.rep.samples <- matrix(NA, n.post.samples * n.chains, J.est + J.zero)
      out$y.rep.samples[, z.indx] <- y.rep.samples
      out$y.rep.samples[, -z.indx] <- y.rep.zero.samples
      out$y.rep.samples <- mcmc(out$y.rep.samples)
      out$like.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$like.samples))))
      out$like.samples <- mcmc(out$like.samples[, order(ord), drop = FALSE])
    }
    out$X <- X[order(ord), , drop = FALSE]
    out$X.re <- X.re[order(ord), , drop = FALSE]
    out$X.w <- X.w[order(ord), , drop = FALSE]
    # Account for case when intercept only spatial model. 
    if (p.svc == 1) {
      tmp <- do.call(rbind, lapply(out.tmp, function(a) t(a$w.samples)))
      tmp <- tmp[, order(ord), drop = FALSE]
      out$w.samples <- array(NA, dim = c(p.svc, J.est, n.post.samples * n.chains))
      out$w.samples[1, , ] <- t(tmp)
    } else {
      out$w.samples <- do.call(abind, lapply(out.tmp, function(a) array(a$w.samples, 
        								dim = c(p.svc, J.est, n.post.samples))))
      out$w.samples <- out$w.samples[, order(ord), ]
    }
    out$w.samples <- aperm(out$w.samples, c(3, 1, 2))
    out$mu.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$mu.samples))))
    out$mu.samples <- mcmc(out$mu.samples[, order(ord), drop = FALSE])
    out$y <- y.orig 
    if (p.re > 0) {
      out$sigma.sq.mu.samples <- mcmc(
        do.call(rbind, lapply(out.tmp, function(a) t(a$sigma.sq.mu.samples))))
      colnames(out$sigma.sq.mu.samples) <- x.re.names
      out$beta.star.samples <- mcmc(
        do.call(rbind, lapply(out.tmp, function(a) t(a$beta.star.samples))))
      tmp.names <- unlist(re.level.names)
      beta.star.names <- paste(rep(x.re.names, n.re.long), tmp.names, sep = '-')
      colnames(out$beta.star.samples) <- beta.star.names
      out$re.level.names <- re.level.names
    }
    # Calculate effective sample sizes
    out$ESS <- list()
    out$ESS$beta <- effectiveSize(out$beta.samples)
    out$ESS$tau.sq <- effectiveSize(out$tau.sq.samples)
    out$ESS$theta <- effectiveSize(out$theta.samples)
    if (p.re > 0) {
      out$ESS$sigma.sq.mu <- effectiveSize(out$sigma.sq.mu.samples)
    }
    out$call <- cl
    out$n.samples <- batch.length * n.batch
    out$n.neighbors <- n.neighbors
    out$cov.model.indx <- cov.model.indx
    out$svc.cols <- svc.cols
    out$type <- "NNGP"
    out$n.post <- n.post.samples
    out$n.thin <- n.thin
    out$n.burn <- n.burn
    out$n.chains <- n.chains
    out$re.cols <- re.cols
    out$dist <- family
    if (p.re > 0) {
      out$muRE <- TRUE
    } else {
      out$muRE <- FALSE
    }
  } # NNGP or GP
  class(out) <- "svcAbund"
  out$run.time <- proc.time() - ptm
  return(out)
}

