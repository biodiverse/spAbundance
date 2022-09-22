sfMsAbund <- function(formula, data, inits, priors,  
		     tuning, cov.model = 'exponential', NNGP = TRUE, 
		     n.neighbors = 15, search.type = 'cb', n.factors,
		     n.batch, batch.length, accept.rate = 0.43, family = 'NB',
		     n.omp.threads = 1, verbose = TRUE, n.report = 100, 
		     n.burn = round(.10 * n.samples), n.thin = 1, n.chains = 1,
		     k.fold, k.fold.threads = 1, k.fold.seed = 100, 
		     k.fold.only = FALSE, ...){

  ptm <- proc.time()

  # Functions -----------------------------------------------------------
  logit <- function(theta, a = 0, b = 1) {log((theta-a)/(b-theta))}
  logit.inv <- function(z, a = 0, b = 1) {b-(b-a)/(1+exp(z))}
  rigamma <- function(n, a, b){
    1/rgamma(n = n, shape = a, rate = b)
  }
 
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
    stop("error: sfMsAbund is currently only implemented for NNGPs, not full Gaussian Processes. Please set NNGP = TRUE.") 
  }
  if (missing(data)) {
    stop("error: data must be specified")
  }
  if (!is.list(data)) {
    stop("error: data must be a list")
  }
  names(data) <- tolower(names(data))
  if (!'y' %in% names(data)) {
    stop("error: count data y must be specified in data")
  }
  if (!(length(dim(data$y)) %in% c(2, 3))) {
    stop("error: count data y must be a two or three-dimensional array with dimensions corresponding to species, sites, and replicates.")
  }
  if (length(dim(data$y)) == 2) {
    data$y <- array(data$y, dim = c(nrow(data$y), ncol(data$y), 1))
  }
  y <- data$y
  sp.names <- attr(y, 'dimnames')[[1]]
  if (!'covs' %in% names(data)) {
    if (formula == ~ 1) {
      if (verbose) {
        message("abundance covariates (covs) not specified in data.\nAssuming intercept only abundance model.\n")
      }
      data$covs <- list(int = array(1, dim = dim(y)))
    } else {
      stop("error: covs must be specified in data for an abundance model with covariates")
    }
  }
  if (!is.list(data$covs)) {
    stop("error: covs must be a list of matrices, data frames, and/or vectors")
  }
  if (!'coords' %in% names(data)) {
    stop("error: coords must be specified in data for a spatial abundance model.")
  }
  if (!is.matrix(data$coords) & !is.data.frame(data$coords)) {
    stop("error: coords must be a matrix or data frame")
  }
  coords <- as.matrix(data$coords)
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
  if (!missing(k.fold)) {
    if (!is.numeric(k.fold) | length(k.fold) != 1 | k.fold < 2) {
      stop("error: k.fold must be a single integer value >= 2")  
    }
  }
  if (missing(n.factors)) {
    stop("error: n.factors must be specified for a spatial factor N-mixture model")
  }

  # Neighbors and Ordering ----------------------------------------------
  if (NNGP) {
    u.search.type <- 2 
    ## Order by x column. Could potentially allow this to be user defined. 
    ord <- order(coords[,1]) 
    # Reorder everything to align with NN ordering
    y <- y[, ord, , drop = FALSE]
    coords <- coords[ord, , drop = FALSE]
    # Covariates
    for (i in 1:length(data$covs)) {
      if (!is.null(dim(data$covs[[i]]))) {
        data$covs[[i]] <- data$covs[[i]][ord, , drop = FALSE]
      } else {
        data$covs[[i]] <- data$covs[[i]][ord]
      }
    }
  }
  y.mat <- y

  # First subset covariates to only use those that are included in the analysis. 
  # Get occurrence covariates in proper format
  # Subset covariates to only use those that are included in the analysis
  data$covs <- data$covs[names(data$covs) %in% all.vars(formula)]
  # Null model support
  if (length(data$covs) == 0) {
    data$covs <- list(int = matrix(1, nrow = dim(y)[1], ncol = dim(y)[2]))
  }
  # Ordered by rep, then site within rep
  data$covs <- data.frame(lapply(data$covs, function(a) unlist(c(a))))
  # Check if only site-level covariates are included
  if (nrow(data$covs) == dim(y)[1]) {
    data$covs <- as.data.frame(mapply(rep, data$covs, dim(y)[2]))
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

  # Checking missing values ---------------------------------------------
  # TODO: I believe these checks will fail if only site-level covariates on 
  #       abundance
  # y -------------------------------
  y.na.test <- apply(y.mat, c(1, 2), function(a) sum(!is.na(a)))
  if (sum(y.na.test == 0) > 0) {
    stop("error: some sites in y have all missing detection histories. Remove these sites from all objects in the 'data' argument, then use 'predict' to obtain predictions at these locations if desired.")
  }
  # covs ------------------------
  for (i in 1:ncol(data$covs)) {
    # Note that this assumes the same detection history for each species.  
    if (sum(is.na(data$covs[, i])) > sum(is.na(y.mat[1, , ]))) {
      stop("error: some elements in covs have missing values where there is an observed data value in y. Please either replace the NA values in covs with non-missing values (e.g., mean imputation) or set the corresponding values in y to NA where the covariate is missing.") 
    }
  }
  # Misalignment between y and covs
  y.missing <- which(is.na(y[1, , ]))
  covs.missing <- lapply(data$covs, function(a) which(is.na(a)))
  for (i in 1:length(covs.missing)) {
    tmp.indx <- !(y.missing %in% covs.missing[[i]])
    if (sum(tmp.indx) > 0) {
      if (i == 1 & verbose) {
        message("There are missing values in data$y with corresponding non-missing values in data$covs.\nRemoving these site/replicate combinations for fitting the model.")
      }
      data$covs[y.missing, i] <- NA
    }
  }
  # Formula -------------------------------------------------------------
  # Abundance -----------------------
  if (missing(formula)) {
    stop("error: formula must be specified")
  }

  if (is(formula, 'formula')) {
    tmp <- parseFormula(formula, data$covs)
    X <- as.matrix(tmp[[1]])
    X.re <- as.matrix(tmp[[4]])
    x.re.names <- colnames(X.re)
    x.names <- tmp[[2]]
  } else {
    stop("error: formula is misspecified")
  }
  # Get RE level names
  re.level.names <- lapply(data$covs[, x.re.names, drop = FALSE],
      		     function (a) sort(unique(a)))


  # Extract data from inputs --------------------------------------------
  # Number of species 
  n.sp <- dim(y)[1]
  # Number of latent factors
  q <- n.factors
  # Number of abundance parameters 
  p.abund <- ncol(X)
  # Number of abundance random effect parameters
  p.abund.re <- ncol(X.re)
  # Number of latent abundance random effect values
  n.abund.re <- length(unlist(apply(X.re, 2, unique)))
  n.abund.re.long <- apply(X.re, 2, function(a) length(unique(a)))
  # Number of sites
  J <- nrow(coords)
  # Number of replicate surveys
  # Note this assumes equivalent detection histories for all species. 
  # May want to change this at some point. 
  n.rep <- apply(y.mat[1, , , drop = FALSE], 2, function(a) sum(!is.na(a)))
  K.max <- max(n.rep)
  # Because I like K better than n.rep
  K <- n.rep

  # Get indices to map N to y -------------------------------------------
  site.indx <- rep(1:J, dim(y.mat)[3]) 
  site.indx <- site.indx[!is.na(c(y.mat[1, , ]))]
  # Subtract 1 for indices in C
  site.indx <- site.indx - 1
  # y is stored in the following order: species, site, visit
  y <- c(y)
  # Assumes the missing data are constant across species, which seems likely, 
  # but may eventually need some updating. 
  names.long <- which(!is.na(c(y.mat[1, , ])))
  # Only need to check this when there are observation level covariates. 
  if (nrow(X) == length(y) / n.sp) {
    X <- X[!is.na(c(y.mat[1, , ])), , drop = FALSE]
  }
  if (nrow(X.re) == length(y) / n.sp & p.abund.re > 0) {
    X.re <- X.re[!is.na(c(y.mat[1, , ])), , drop = FALSE]
  }
  y <- y[!is.na(y)]
  # Number of pseudoreplicates
  n.obs <- nrow(X)

  # Get random effect matrices all set ----------------------------------
  if (p.abund.re > 1) {
    for (j in 2:p.abund.re) {
      X.re[, j] <- X.re[, j] + max(X.re[, j - 1]) + 1
    }
  }

  # Separate out priors -------------------------------------------------
  if (missing(priors)) {
    priors <- list()
  }
  names(priors) <- tolower(names(priors))

  # beta.comm -----------------------
  if ("beta.comm.normal" %in% names(priors)) {
    if (!is.list(priors$beta.comm.normal) | length(priors$beta.comm.normal) != 2) {
      stop("error: beta.comm.normal must be a list of length 2")
    }
    mu.beta.comm <- priors$beta.comm.normal[[1]]
    sigma.beta.comm <- priors$beta.comm.normal[[2]]
    if (length(mu.beta.comm) != p.abund & length(mu.beta.comm) != 1) {
      if (p.abund == 1) {
        stop(paste("error: beta.comm.normal[[1]] must be a vector of length ",
        	     p.abund, " with elements corresponding to beta.comms' mean", sep = ""))
      } else {
        stop(paste("error: beta.comm.normal[[1]] must be a vector of length ",
        	     p.abund, " or 1 with elements corresponding to beta.comms' mean", sep = ""))
      }
    }
    if (length(sigma.beta.comm) != p.abund & length(sigma.beta.comm) != 1) {
      if (p.abund == 1) {
        stop(paste("error: beta.comm.normal[[2]] must be a vector of length ",
      	   p.abund, " with elements corresponding to beta.comms' variance", sep = ""))
      } else {
        stop(paste("error: beta.comm.normal[[2]] must be a vector of length ",
      	   p.abund, " or 1 with elements corresponding to beta.comms' variance", sep = ""))
      }
    }
    if (length(sigma.beta.comm) != p.abund) {
      sigma.beta.comm <- rep(sigma.beta.comm, p.abund)
    }
    if (length(mu.beta.comm) != p.abund) {
      mu.beta.comm <- rep(mu.beta.comm, p.abund)
    }
    Sigma.beta.comm <- sigma.beta.comm * diag(p.abund)
  } else {
    if (verbose) {
      message("No prior specified for beta.comm.normal.\nSetting prior mean to 0 and prior variance to 100\n")
    }
    mu.beta.comm <- rep(0, p.abund)
    sigma.beta.comm <- rep(100, p.abund)
    Sigma.beta.comm <- diag(p.abund) * 100 
  }

  # tau.sq.beta -----------------------
  if ("tau.sq.beta.ig" %in% names(priors)) {
    if (!is.list(priors$tau.sq.beta.ig) | length(priors$tau.sq.beta.ig) != 2) {
      stop("error: tau.sq.beta.ig must be a list of length 2")
    }
    tau.sq.beta.a <- priors$tau.sq.beta.ig[[1]]
    tau.sq.beta.b <- priors$tau.sq.beta.ig[[2]]
    if (length(tau.sq.beta.a) != p.abund & length(tau.sq.beta.a) != 1) {
      if (p.abund == 1) {
        stop(paste("error: tau.sq.beta.ig[[1]] must be a vector of length ", 
      	   p.abund, " with elements corresponding to tau.sq.betas' shape", sep = ""))
      } else {
        stop(paste("error: tau.sq.beta.ig[[1]] must be a vector of length ", 
      	   p.abund, " or 1 with elements corresponding to tau.sq.betas' shape", sep = ""))
      }
    }
    if (length(tau.sq.beta.b) != p.abund & length(tau.sq.beta.b) != 1) {
      if (p.abund == 1) {
        stop(paste("error: tau.sq.beta.ig[[2]] must be a vector of length ", 
      	   p.abund, " with elements corresponding to tau.sq.betas' scale", sep = ""))
      } else {
        stop(paste("error: tau.sq.beta.ig[[2]] must be a vector of length ", 
      	   p.abund, " or 1 with elements corresponding to tau.sq.betas' scale", sep = ""))
      }
    }
    if (length(tau.sq.beta.a) != p.abund) {
      tau.sq.beta.a <- rep(tau.sq.beta.a, p.abund)
    }
    if (length(tau.sq.beta.b) != p.abund) {
      tau.sq.beta.b <- rep(tau.sq.beta.b, p.abund)
    }
  } else {
    if (verbose) {	    
      message("No prior specified for tau.sq.beta.ig.\nSetting prior shape to 0.1 and prior scale to 0.1\n")
    }
    tau.sq.beta.a <- rep(0.1, p.abund)
    tau.sq.beta.b <- rep(0.1, p.abund)
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
  # kappa -----------------------------
  if ("kappa.unif" %in% names(priors)) {
    if (!is.list(priors$kappa.unif) | length(priors$kappa.unif) != 2) {
      stop("error: kappa.unif must be a list of length 2")
    }
    kappa.a <- priors$kappa.unif[[1]]
    kappa.b <- priors$kappa.unif[[2]]
    if (length(kappa.a) != n.sp & length(kappa.a) != 1) {
      stop(paste("error: kappa.unif[[1]] must be a vector of length ", 
      	   n.sp, " or 1 with elements corresponding to kappas' lower bound for each species", sep = ""))
    }
    if (length(kappa.b) != n.sp & length(kappa.b) != 1) {
      stop(paste("error: kappa.unif[[2]] must be a vector of length ", 
      	   n.sp, " or 1 with elements corresponding to kappas' upper bound for each species", sep = ""))
    }
    if (length(kappa.a) != n.sp) {
      kappa.a <- rep(kappa.a, n.sp)
    }
    if (length(kappa.b) != n.sp) {
      kappa.b <- rep(kappa.b, n.sp)
    }
  } else {
    if (verbose) {
    message("No prior specified for kappa.unif.\nSetting uniform bounds of 0 and 10.\n")
    }
    kappa.a <- rep(0, n.sp)
    kappa.b <- rep(10, n.sp)
  }
  # phi -----------------------------
  coords.D <- iDist(coords)
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

  # Starting values -----------------------------------------------------
  if (missing(inits)) {
    inits <- list()
  }
  names(inits) <- tolower(names(inits))
  # beta.comm -----------------------
  if ("beta.comm" %in% names(inits)) {
    beta.comm.inits <- inits[["beta.comm"]]
    if (length(beta.comm.inits) != p.abund & length(beta.comm.inits) != 1) {
      if (p.abund == 1) {
        stop(paste("error: initial values for beta.comm must be of length ", p.abund, 
      	   sep = ""))
      } else {
        stop(paste("error: initial values for beta.comm must be of length ", p.abund, 
      	   , " or 1", sep = ""))
      }
    }
    if (length(beta.comm.inits) != p.abund) {
      beta.comm.inits <- rep(beta.comm.inits, p.abund)
    }
  } else {
    beta.comm.inits <- rnorm(p.abund, mu.beta.comm, sqrt(sigma.beta.comm))
    if (verbose) {
      message('beta.comm is not specified in initial values.\nSetting initial values to random values from the prior distribution\n')
    }
  }
  # tau.sq.beta ------------------------
  if ("tau.sq.beta" %in% names(inits)) {
    tau.sq.beta.inits <- inits[["tau.sq.beta"]]
    if (length(tau.sq.beta.inits) != p.abund & length(tau.sq.beta.inits) != 1) {
      if (p.abund == 1) {
        stop(paste("error: initial values for tau.sq.beta must be of length ", p.abund, 
      	   sep = ""))
      } else {
        stop(paste("error: initial values for tau.sq.beta must be of length ", p.abund, 
      	   " or 1", sep = ""))
      }
    }
    if (length(tau.sq.beta.inits) != p.abund) {
      tau.sq.beta.inits <- rep(tau.sq.beta.inits, p.abund)
    }
  } else {
    tau.sq.beta.inits <- runif(p.abund, 0.5, 10)
    if (verbose) {
      message('tau.sq.beta is not specified in initial values.\nSetting initial values to random values between 0.5 and 10\n')
    }
  }
  # beta ----------------------------
  if ("beta" %in% names(inits)) {
    beta.inits <- inits[["beta"]]
    if (is.matrix(beta.inits)) {
      if (ncol(beta.inits) != p.abund | nrow(beta.inits) != n.sp) {
        stop(paste("error: initial values for beta must be a matrix with dimensions ", 
        	   n.sp, "x", p.abund, " or a single numeric value", sep = ""))
      }
    }
    if (!is.matrix(beta.inits) & length(beta.inits) != 1) {
      stop(paste("error: initial values for beta must be a matrix with dimensions ", 
      	   n.sp, " x ", p.abund, " or a single numeric value", sep = ""))
    }
    if (length(beta.inits) == 1) {
      beta.inits <- matrix(beta.inits, n.sp, p.abund)
    }
  } else {
    beta.inits <- matrix(rnorm(n.sp * p.abund, beta.comm.inits, sqrt(tau.sq.beta.inits)), n.sp, p.abund)
    if (verbose) {
      message('beta is not specified in initial values.\nSetting initial values to random values from the community-level normal distribution\n')
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
      sigma.sq.mu.inits <- runif(p.abund.re, 0.5, 10)
      if (verbose) {
        message("sigma.sq.mu is not specified in initial values.\nSetting initial values to random values between 0.5 and 10\n")
      }
    }
    beta.star.indx <- rep(0:(p.abund.re - 1), n.abund.re.long)
    beta.star.inits <- rnorm(n.abund.re, sqrt(sigma.sq.mu.inits[beta.star.indx + 1]))
    # Starting values for all species 
    beta.star.inits <- rep(beta.star.inits, n.sp)
  } else {
    sigma.sq.mu.inits <- 0
    beta.star.indx <- 0
    beta.star.inits <- 0
  }
  # kappa -----------------------------
  # ORDER: a length n.sp vector ordered by species in the detection-nondetection data.
  if ("kappa" %in% names(inits)) {
    kappa.inits <- inits[["kappa"]]
    if (length(kappa.inits) != n.sp & length(kappa.inits) != 1) {
      stop(paste("error: initial values for kappa must be of length ", n.sp, " or 1", 
      	   sep = ""))
    }
    if (length(kappa.inits) != n.sp) {
      kappa.inits <- rep(kappa.inits, n.sp)
    }
  } else {
    kappa.inits <- runif(n.sp, kappa.a, kappa.b)
    if (verbose) {
      message("kappa is not specified in initial values.\nSetting initial value to random values from the prior distribution\n")
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
		 n.sp, " x ", q, sep = ""))
    }
    if (nrow(lambda.inits) != n.sp | ncol(lambda.inits) != q) {
      stop(paste("error: initial values for lambda must be a matrix with dimensions ",
		 n.sp, " x ", q, sep = ""))
    }
    if (!all.equal(diag(lambda.inits), rep(1, q))) {
      stop("error: diagonal of inits$lambda matrix must be all 1s")
    }
    if (sum(lambda.inits[upper.tri(lambda.inits)]) != 0) {
      stop("error: upper triangle of inits$lambda must be all 0s")
    }
  } else {
    lambda.inits <- matrix(0, n.sp, q)
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
  # Covariance Model ----------------------------------------------------
  # Order must match util.cpp spCor.
  cov.model.names <- c("exponential", "spherical", "matern", "gaussian")
  if(! cov.model %in% cov.model.names){
    stop("error: specified cov.model '",cov.model,"' is not a valid option; choose from ", 
         paste(cov.model.names, collapse=", ", sep="") ,".")}
  # Obo for cov model lookup on c side
  cov.model.indx <- which(cov.model == cov.model.names) - 1

  # Get tuning values ---------------------------------------------------
  sigma.sq.tuning <- rep(0, q)
  phi.tuning <- rep(0, q)
  nu.tuning <- rep(0, q)
  kappa.tuning <- rep(0, n.sp)
  if (missing(tuning)) {
    phi.tuning <- rep(1, q)
    if (cov.model == 'matern') {
      nu.tuning <- rep(1, q)
    }
    kappa.tuning <- rep(1, n.sp)
  } else {
    names(tuning) <- tolower(names(tuning))
    # kappa ---------------------------
    if(!"kappa" %in% names(tuning)) {
      stop("error: kappa must be specified in tuning value list")
    }
    kappa.tuning <- tuning$kappa
    if (length(kappa.tuning) == 1) {
      kappa.tuning <- rep(tuning$kappa, n.sp)
    } else if (length(kappa.tuning) != n.sp) {
      stop(paste("error: kappa tuning must be either a single value or a vector of length ",
      	   n.sp, sep = ""))
    }
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
  tuning.c <- log(c(sigma.sq.tuning, phi.tuning, nu.tuning, kappa.tuning))

  curr.chain <- 1
  if (!NNGP)  {
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
    storage.mode(coords) <- "double"
    consts <- c(n.sp, J, n.obs, p.abund, p.abund.re, n.abund.re, q)
    storage.mode(consts) <- "integer"
    storage.mode(beta.inits) <- "double"
    storage.mode(kappa.inits) <- "double"
    storage.mode(beta.comm.inits) <- "double"
    storage.mode(tau.sq.beta.inits) <- "double"
    storage.mode(lambda.inits) <- "double"
    storage.mode(w.inits) <- "double"
    storage.mode(phi.inits) <- "double"
    storage.mode(nu.inits) <- "double"
    storage.mode(site.indx) <- "integer"
    storage.mode(mu.beta.comm) <- "double"
    storage.mode(Sigma.beta.comm) <- "double"
    storage.mode(kappa.a) <- "double"
    storage.mode(kappa.b) <- "double"
    storage.mode(tau.sq.beta.a) <- "double"
    storage.mode(tau.sq.beta.b) <- "double"
    storage.mode(phi.a) <- "double"
    storage.mode(phi.b) <- "double"
    storage.mode(nu.a) <- "double"
    storage.mode(nu.b) <- "double"
    storage.mode(n.batch) <- "integer"
    storage.mode(batch.length) <- "integer"
    storage.mode(accept.rate) <- "double"
    storage.mode(tuning.c) <- "double"
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
    # samples.info order: burn-in, thinning rate, number of posterior samples
    samples.info <- c(n.burn, n.thin, n.post.samples)
    storage.mode(samples.info) <- "integer"
    # For abundance random effects
    storage.mode(X.re) <- "integer"
    beta.level.indx <- sort(unique(c(X.re)))
    storage.mode(beta.level.indx) <- "integer"
    storage.mode(sigma.sq.mu.inits) <- "double"
    storage.mode(sigma.sq.mu.a) <- "double"
    storage.mode(sigma.sq.mu.b) <- "double"
    storage.mode(beta.star.inits) <- "double"
    storage.mode(beta.star.indx) <- "integer"

    # Fit the model -------------------------------------------------------
    out.tmp <- list()
    for (i in 1:n.chains) {
      # Change initial values if i > 1
      if ((i > 1) & (!fix.inits)) {
        beta.comm.inits <- rnorm(p.abund, mu.beta.comm, sqrt(sigma.beta.comm))
        tau.sq.beta.inits <- runif(p.abund, 0.5, 10)
        beta.inits <- matrix(rnorm(n.sp * p.abund, beta.comm.inits, 
              		     sqrt(tau.sq.beta.inits)), n.sp, p.abund)
        kappa.inits <- runif(n.sp, kappa.a, kappa.b)
        lambda.inits <- matrix(0, n.sp, q)
        diag(lambda.inits) <- 1
        lambda.inits[lower.tri(lambda.inits)] <- rnorm(sum(lower.tri(lambda.inits)))
        lambda.inits <- c(lambda.inits)
        phi.inits <- runif(q, phi.a, phi.b)
        if (cov.model == 'matern') {
          nu.inits <- runif(q, nu.a, nu.b)
        }
        if (p.abund.re > 0) {
          sigma.sq.mu.inits <- runif(p.abund.re, 0.5, 10)
          beta.star.inits <- rnorm(n.abund.re, sqrt(sigma.sq.mu.inits[beta.star.indx + 1]))
          beta.star.inits <- rep(beta.star.inits, n.sp)
        }
      }

      storage.mode(chain.info) <- "integer"
      out.tmp[[i]] <- .Call("sfMsAbundNNGP", y, X, coords, X.re, consts, n.abund.re.long, 
        	            n.neighbors, nn.indx, nn.indx.lu, u.indx, u.indx.lu, ui.indx,
			    beta.inits, kappa.inits, beta.comm.inits, 
        	            tau.sq.beta.inits, 
			    phi.inits, lambda.inits, nu.inits, w.inits,
        		    sigma.sq.mu.inits, beta.star.inits,site.indx, 
        		    beta.star.indx, beta.level.indx,  
        		    mu.beta.comm, Sigma.beta.comm, kappa.a, 
          		    kappa.b, tau.sq.beta.a, tau.sq.beta.b,  
        	            phi.a, phi.b, nu.a, nu.b, 
			    sigma.sq.mu.a, sigma.sq.mu.b,
      	                    tuning.c, cov.model.indx, 
			    n.batch, batch.length, accept.rate, n.omp.threads, 
      	                    verbose, n.report, samples.info, chain.info)
      chain.info[1] <- chain.info[1] + 1
    }
    # Calculate R-Hat ---------------
    out <- list()
    out$rhat <- list()
    if (n.chains > 1) {
      # as.vector removes the "Upper CI" when there is only 1 variable. 
      out$rhat$beta.comm <- as.vector(gelman.diag(mcmc.list(lapply(out.tmp, function(a) 
      					      mcmc(t(a$beta.comm.samples)))), 
      			     autoburnin = FALSE)$psrf[, 2])
      out$rhat$tau.sq.beta <- as.vector(gelman.diag(mcmc.list(lapply(out.tmp, function(a) 
      					      mcmc(t(a$tau.sq.beta.samples)))), 
      			     autoburnin = FALSE)$psrf[, 2])
      out$rhat$beta <- as.vector(gelman.diag(mcmc.list(lapply(out.tmp, function(a) 
      					         mcmc(t(a$beta.samples)))), 
      			     autoburnin = FALSE)$psrf[, 2])
      out$rhat$kappa <- as.vector(gelman.diag(mcmc.list(lapply(out.tmp, function(a) 
      					      mcmc(t(a$kappa.samples)))), 
      			     autoburnin = FALSE)$psrf[, 2])
      out$rhat$theta <- gelman.diag(mcmc.list(lapply(out.tmp, function(a) 
      					      mcmc(t(a$theta.samples)))), 
      			      autoburnin = FALSE)$psrf[, 2]
      lambda.mat <- matrix(lambda.inits, n.sp, q)
      out$rhat$lambda.lower.tri <- as.vector(gelman.diag(mcmc.list(lapply(out.tmp, function(a) 
						       mcmc(t(a$lambda.samples[c(lower.tri(lambda.mat)), ])))), 
						       autoburnin = FALSE)$psrf[, 2])
      if (p.abund.re > 0) {
        out$rhat$sigma.sq.mu <- as.vector(gelman.diag(mcmc.list(lapply(out.tmp, function(a) 
        					      mcmc(t(a$sigma.sq.mu.samples)))), 
        			     autoburnin = FALSE)$psrf[, 2])
      }
    } else {
      out$rhat$beta.comm <- rep(NA, p.abund)
      out$rhat$tau.sq.beta <- rep(NA, p.abund)
      out$rhat$beta <- rep(NA, p.abund * n.sp)
      out$rhat$theta <- rep(NA, ifelse(cov.model == 'matern', 2 * q, q))
      out$rhat$kappa <- rep(NA, n.sp)
      if (p.abund.re > 0) {
        out$rhat$sigma.sq.mu <- rep(NA, p.abund.re)
      }
    }
    # Put everything into MCMC objects
    out$beta.comm.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$beta.comm.samples))))
    colnames(out$beta.comm.samples) <- x.names
    out$tau.sq.beta.samples <- mcmc(do.call(rbind, 
      				lapply(out.tmp, function(a) t(a$tau.sq.beta.samples))))
    colnames(out$tau.sq.beta.samples) <- x.names

    if (is.null(sp.names)) {
      sp.names <- paste('sp', 1:n.sp, sep = '')
    }
    coef.names <- paste(rep(x.names, each = n.sp), sp.names, sep = '-')
    out$beta.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$beta.samples))))
    colnames(out$beta.samples) <- coef.names
    out$kappa.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$kappa.samples))))
    colnames(out$kappa.samples) <- paste('kappa', 1:n.sp, sep = '') 
    loadings.names <- paste(rep(sp.names, times = n.factors), rep(1:n.factors, each = n.sp), sep = '-')
    out$lambda.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$lambda.samples))))
    colnames(out$lambda.samples) <- loadings.names
    out$theta.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$theta.samples))))
    if (cov.model != 'matern') {
      theta.names <- paste(rep(c('phi'), each = q), 1:q, sep = '-')
    } else {
      theta.names <- paste(rep(c('phi', 'nu'), each = q), 1:q, sep = '-')
    } 
    colnames(out$theta.samples) <- theta.names
    y.non.miss.indx <- which(!is.na(y.mat), arr.ind = TRUE)
    out$y.rep.samples <- do.call(abind, lapply(out.tmp, function(a) array(a$y.rep.samples, 
      								dim = c(n.sp, n.obs, n.post.samples))))
    tmp <- array(NA, dim = c(n.post.samples, n.sp, J, K.max))
    for (j in 1:n.obs) {
      curr.indx <- y.non.miss.indx[j, ]
      tmp[, curr.indx[1], curr.indx[2], curr.indx[3]] <- out$y.rep.samples[curr.indx[1], j, ]
    }
    out$y.rep.samples <- tmp[, , order(ord), , drop = FALSE]
    out$mu.samples <- do.call(abind, lapply(out.tmp, function(a) array(a$mu.samples, 
      								dim = c(n.sp, n.obs, n.post.samples))))
    tmp <- array(NA, dim = c(n.post.samples, n.sp, J, K.max))
    for (j in 1:n.obs) {
      curr.indx <- y.non.miss.indx[j, ]
      tmp[, curr.indx[1], curr.indx[2], curr.indx[3]] <- out$mu.samples[curr.indx[1], j, ]
    }
    out$mu.samples <- tmp[, , order(ord), , drop = FALSE]
    out$like.samples <- do.call(abind, lapply(out.tmp, function(a) array(a$like.samples, 
      								dim = c(n.sp, n.obs, n.post.samples))))
    tmp <- array(NA, dim = c(n.post.samples, n.sp, J, K.max))
    for (j in 1:n.obs) {
      curr.indx <- y.non.miss.indx[j, ]
      tmp[, curr.indx[1], curr.indx[2], curr.indx[3]] <- out$like.samples[curr.indx[1], j, ]
    }
    out$like.samples <- tmp[, , order(ord), , drop = FALSE]
    out$w.samples <- do.call(abind, lapply(out.tmp, function(a) array(a$w.samples, 
      								dim = c(q, J, n.post.samples))))
    out$w.samples <- out$w.samples[, order(ord), , drop = FALSE]
    out$w.samples <- aperm(out$w.samples, c(3, 1, 2))
    if (p.abund.re > 0) {
      out$sigma.sq.mu.samples <- mcmc(
        do.call(rbind, lapply(out.tmp, function(a) t(a$sigma.sq.mu.samples))))
      colnames(out$sigma.sq.mu.samples) <- x.re.names
      out$beta.star.samples <- mcmc(
        do.call(rbind, lapply(out.tmp, function(a) t(a$beta.star.samples))))
      tmp.names <- unlist(re.level.names)
      beta.star.names <- paste(rep(x.re.names, n.abund.re.long), tmp.names, sep = '-')
      beta.star.names <- paste(beta.star.names, rep(sp.names, each = n.abund.re), sep = '-')
      colnames(out$beta.star.samples) <- beta.star.names
      out$re.level.names <- re.level.names
    }
    # Calculate effective sample sizes
    out$ESS <- list()
    out$ESS$beta.comm <- effectiveSize(out$beta.comm.samples)
    out$ESS$tau.sq.beta <- effectiveSize(out$tau.sq.beta.samples)
    out$ESS$beta <- effectiveSize(out$beta.samples)
    out$ESS$kappa <- effectiveSize(out$kappa.samples)
    out$ESS$theta <- effectiveSize(out$theta.samples)
    out$ESS$lambda <- effectiveSize(out$lambda.samples)
    if (p.abund.re > 0) {
      out$ESS$sigma.sq.mu <- effectiveSize(out$sigma.sq.mu.samples)
    }
    out$X.re <- X.re[order(ord), , drop = FALSE]
    out$X <- X[order(ord), , drop = FALSE]
    out$y <- y.mat[, order(ord), , drop = FALSE]
    out$call <- cl
    out$n.samples <- n.samples
    out$x.names <- x.names
    out$sp.names <- sp.names
    out$n.post <- n.post.samples
    out$n.thin <- n.thin
    out$n.burn <- n.burn
    out$n.chains <- n.chains
    out$theta.names <- theta.names
    out$type <- "NNGP"
    out$coords <- coords[order(ord), ]
    out$cov.model.indx <- cov.model.indx
    out$n.neighbors <- n.neighbors
    out$dist <- 'NB'
    out$q <- q
    if (p.abund.re > 0) {
      out$muRE <- TRUE
    } else {
      out$muRE <- FALSE
    }
    # TODO: need to update. 
    # K-fold cross-validation -------
    # if (!missing(k.fold)) {
    #   if (verbose) {      
    #     cat("----------------------------------------\n");
    #     cat("\tCross-validation\n");
    #     cat("----------------------------------------\n");
    #     message(paste("Performing ", k.fold, "-fold cross-validation using ", k.fold.threads,
    #   	      " thread(s).", sep = ''))
    #   }
    #   # Currently implemented without parellization. 
    #   set.seed(k.fold.seed)
    #   # Number of sites in each hold out data set. 
    #   sites.random <- sample(1:J)    
    #   sites.k.fold <- split(sites.random, sites.random %% k.fold)
    #   registerDoParallel(k.fold.threads)
    #   model.deviance <- foreach (i = 1:k.fold, .combine = "+") %dopar% {
    #     curr.set <- sort(sites.random[sites.k.fold[[i]]])
    #     if (binom) {
    #       y.indx <- !(1:J %in% curr.set)
    #       y.fit <- y[rep(y.indx, n.sp), drop = FALSE]
    #       y.0 <- y[rep(y.indx, n.sp), drop = FALSE]
    #     } else {
    #       y.indx <- !((site.indx + 1) %in% curr.set)
    #       y.fit <- c(y.mat[, -curr.set, , drop = FALSE])
    #       y.fit <- y.fit[!is.na(y.fit)]
    #       y.0 <- c(y.mat[, curr.set, , drop = FALSE])
    #       y.0 <- y.0[!is.na(y.0)]
    #     }
    #     z.inits.fit <- z.inits[, -curr.set]
    #     y.mat.fit <- y.mat[, -curr.set, , drop = FALSE]
    #     y.mat.0 <- y.mat[, curr.set, , drop = FALSE]
    #     X.p.fit <- X.p[y.indx, , drop = FALSE]
    #     X.p.0 <- X.p[!y.indx, , drop = FALSE]
    #     X.fit <- X[-curr.set, , drop = FALSE]
    #     X.0 <- X[curr.set, , drop = FALSE]
    #     J.fit <- nrow(X.fit)
    #     J.0 <- nrow(X.0)
    #     K.fit <- K[-curr.set]
    #     K.0 <- K[curr.set]
    #     n.obs.fit <- nrow(X.p.fit)
    #     n.obs.0 <- nrow(X.p.0)
    #     # Random detection effects
    #     lambda.p.fit <- lambda.p[y.indx, , drop = FALSE]
    #     lambda.p.0 <- lambda.p[!y.indx, , drop = FALSE]
    #     X.p.re.fit <- X.p.re[y.indx, , drop = FALSE]
    #     X.p.re.0 <- X.p.re[!y.indx, , drop = FALSE]
    #     n.det.re.fit <- length(unique(c(X.p.re.fit)))
    #     n.det.re.long.fit <- apply(X.p.re.fit, 2, function(a) length(unique(a)))
    #     if (p.det.re > 0) {	
    #       alpha.star.indx.fit <- rep(0:(p.det.re - 1), n.det.re.long.fit)
    #       alpha.level.indx.fit <- sort(unique(c(X.p.re.fit)))
    #       alpha.star.inits.fit <- rnorm(n.det.re.fit, 
    #       			      sqrt(sigma.sq.p.inits[alpha.star.indx.fit + 1]))
    #       alpha.star.inits.fit <- rep(alpha.star.inits.fit, n.sp)
    #     } else {
    #       alpha.star.indx.fit <- alpha.star.indx
    #       alpha.level.indx.fit <- alpha.level.indx
    #       alpha.star.inits.fit <- alpha.star.inits
    #     }
    #     # Random abundance effects
    #     X.re.fit <- X.re[-curr.set, , drop = FALSE]
    #     X.re.0 <- X.re[curr.set, , drop = FALSE]
    #     n.abund.re.fit <- length(unique(c(X.re.fit)))
    #     n.abund.re.long.fit <- apply(X.re.fit, 2, function(a) length(unique(a)))
    #     if (p.abund.re > 0) {	
    #       beta.star.indx.fit <- rep(0:(p.abund.re - 1), n.abund.re.long.fit)
    #       beta.level.indx.fit <- sort(unique(c(X.re.fit)))
    #       beta.star.inits.fit <- rnorm(n.abund.re.fit, 
    #       			      sqrt(sigma.sq.mu.inits[beta.star.indx.fit + 1]))
    #       beta.star.inits.fit <- rep(beta.star.inits.fit, n.sp)
    #       re.level.names.fit <- list()
    #       for (t in 1:p.abund.re) {
    #         tmp.indx <- beta.level.indx.fit[beta.star.indx.fit == t - 1]
    #         re.level.names.fit[[t]] <- unlist(re.level.names)[tmp.indx + 1]    
    #       }
    #     } else {
    #       beta.star.indx.fit <- beta.star.indx
    #       beta.level.indx.fit <- beta.level.indx
    #       beta.star.inits.fit <- beta.star.inits
    #       re.level.names.fit <- re.level.names
    #     }
    #     # Gotta be a better way, but will do for now. 
    #     if (binom) {
    #       site.indx.fit <- 0:(J.fit - 1)
    #       z.0.long.indx <- 1:J.0
    #     } else {
    #       site.indx.fit <- matrix(NA, J.fit, max(K.fit))
    #       for (j in 1:J.fit) {
    #         site.indx.fit[j, 1:K.fit[j]] <- j  
    #       }
    #       site.indx.fit <- c(site.indx.fit)
    #       site.indx.fit <- site.indx.fit[!is.na(site.indx.fit)] - 1
    #       z.0.long.indx <- matrix(NA, nrow(X.0), max(K.0))
    #       for (j in 1:nrow(X.0)) {
    #         z.0.long.indx[j, 1:K.0[j]] <- j  
    #       }
    #       z.0.long.indx <- c(z.0.long.indx)
    #       z.0.long.indx <- z.0.long.indx[!is.na(z.0.long.indx)] 
    #     }
    #     verbose.fit <- FALSE
    #     n.omp.threads.fit <- 1

    #     storage.mode(y.fit) <- "double"
    #     storage.mode(z.inits.fit) <- "double"
    #     storage.mode(X.p.fit) <- "double"
    #     storage.mode(X.fit) <- "double"
    #     storage.mode(K.fit) <- "double"
    #     consts.fit <- c(n.sp, J.fit, n.obs.fit, p.abund, p.abund.re, n.abund.re.fit, 
    #                     p.det, p.det.re, n.det.re.fit)
    #     storage.mode(consts.fit) <- "integer"
    #     storage.mode(beta.inits) <- "double"
    #     storage.mode(alpha.inits) <- "double"
    #     storage.mode(site.indx.fit) <- "integer"
    #     storage.mode(n.samples) <- "integer"
    #     storage.mode(n.omp.threads.fit) <- "integer"
    #     storage.mode(verbose.fit) <- "integer"
    #     storage.mode(n.report) <- "integer"
    #     storage.mode(X.p.re.fit) <- "integer"
    #     storage.mode(n.det.re.long.fit) <- "integer"
    #     storage.mode(alpha.star.inits.fit) <- "double"
    #     storage.mode(alpha.star.indx.fit) <- "integer"
    #     storage.mode(alpha.level.indx.fit) <- "integer"
    #     storage.mode(X.re.fit) <- "integer"
    #     storage.mode(n.abund.re.long.fit) <- "integer"
    #     storage.mode(beta.star.inits.fit) <- "double"
    #     storage.mode(beta.star.indx.fit) <- "integer"
    #     storage.mode(beta.level.indx.fit) <- "integer"
    #     chain.info[1] <- 1
    #     storage.mode(chain.info) <- "integer"

    #     out.fit <- .Call("sfMsAbundNNGP", y.fit, X.fit, X.p.fit, X.re.fit, X.p.re.fit, consts.fit, 
    #   	                 K.fit, n.abund.re.long.fit, n.det.re.long.fit,
    #     	         beta.inits, alpha.inits, z.inits.fit, beta.comm.inits, 
    #     	         alpha.comm.inits, tau.sq.beta.inits, tau.sq.alpha.inits, 
    #     		 sigma.sq.mu.inits, sigma.sq.p.inits, 
    #   	                 beta.star.inits.fit, alpha.star.inits.fit, site.indx.fit, 
    #     		 beta.star.indx.fit, beta.level.indx.fit, alpha.star.indx.fit, 
    #     		 alpha.level.indx.fit, mu.beta.comm, mu.alpha.comm, 
    #     		 Sigma.beta.comm, Sigma.alpha.comm, 
    #     	         tau.sq.beta.a, tau.sq.beta.b, tau.sq.alpha.a, 
    #     	         tau.sq.alpha.b, sigma.sq.mu.a, sigma.sq.mu.b, 
    #     		 sigma.sq.p.a, sigma.sq.p.b,
    #   	                 n.samples, n.omp.threads.fit, 
    #   	                 verbose.fit, n.report, samples.info, chain.info)

    #     if (is.null(sp.names)) {
    #       sp.names <- paste('sp', 1:n.sp, sep = '')
    #     }
    #     coef.names <- paste(rep(x.names, each = n.sp), sp.names, sep = '-')
    #     out.fit$beta.samples <- mcmc(t(out.fit$beta.samples))
    #     colnames(out.fit$beta.samples) <- coef.names
    #     out.fit$X <- X.fit
    #     out.fit$y <- y.mat.fit
    #     out.fit$X.p <- X.p.fit
    #     out.fit$call <- cl
    #     out.fit$n.samples <- n.samples
    #     out.fit$n.post <- n.post.samples
    #     out.fit$n.thin <- n.thin
    #     out.fit$n.burn <- n.burn
    #     out.fit$n.chains <- 1
    #     if (p.det.re > 0) {
    #       out.fit$pRE <- TRUE
    #     } else {
    #       out.fit$pRE <- FALSE
    #     }
    #     if (p.abund.re > 0) {
    #       out.fit$sigma.sq.mu.samples <- mcmc(t(out.fit$sigma.sq.mu.samples))
    #       colnames(out.fit$sigma.sq.mu.samples) <- x.re.names
    #       out.fit$beta.star.samples <- mcmc(t(out.fit$beta.star.samples))
    #       tmp.names <- unlist(re.level.names.fit)
    #       beta.star.names <- paste(rep(x.re.names, n.abund.re.long.fit), tmp.names, sep = '-')
    #       beta.star.names <- paste(beta.star.names, rep(sp.names, each = n.abund.re.fit), 
    #     			   sep = '-')
    #       colnames(out.fit$beta.star.samples) <- beta.star.names
    #       out.fit$re.level.names <- re.level.names.fit
    #       out.fit$X.re <- X.re.fit
    #     }
    #     if (p.abund.re > 0) {
    #       out.fit$psiRE <- TRUE
    #     } else {
    #       out.fit$psiRE <- FALSE	
    #     }
    #     class(out.fit) <- "sfMsAbund"

    #     # Predict abundance at new sites. 
    #     if (p.abund.re > 0) {X.0 <- cbind(X.0, X.re.0)}
    #     out.pred <- predict.sfMsAbund(out.fit, X.0)

    #     # Get full random effects if certain levels aren't in the fitted values
    #     if (p.det.re > 0) {
    #       if (n.det.re.fit != n.det.re) {
    #         tmp <- matrix(NA, n.det.re * n.sp, n.post.samples)  
    #         for (l in 1:n.sp) {
    #           tmp[alpha.level.indx.fit + n.det.re * (l - 1) + 1, ] <- out.fit$alpha.star.samples[1:n.det.re.fit + n.det.re.fit * (l - 1), ]
    #         }
    #         out.fit$alpha.star.samples <- tmp
    #       }
    #       # Samples missing NA values
    #       tmp.indx <- which(apply(out.fit$alpha.star.samples, 1, function(a) sum(is.na(a))) == n.post.samples)
    #       curr.indx <- rep(alpha.star.indx, n.sp)
    #       for (l in tmp.indx) {
    #         out.fit$alpha.star.samples[l, ] <- rnorm(n.post.samples, 0, 
    #     					     sqrt(out.fit$sigma.sq.p.samples[curr.indx[l] + 1, ]))
    #       }
    #     }
    #     
    #     # Detection 
    #     sp.indx <- rep(1:n.sp, ncol(X.p.0))
    #     p.0.samples <- array(NA, dim = c(nrow(X.p.0), n.sp, n.post.samples))
    #     if (p.det.re > 0) {
    #       sp.re.indx <- rep(1:n.sp, each = nrow(out.fit$alpha.star.samples) / n.sp)
    #     }
    #     if (binom) {
    #       like.samples <- array(NA, c(n.sp, nrow(X.p.0), dim(y.mat.0)[3]))
    #       for (q in 1:n.sp) {
    #         if (p.det.re > 0) {
    #           p.0.samples[, q, ] <- logit.inv(X.p.0 %*% out.fit$alpha.samples[sp.indx == q, ] +
    #     				      lambda.p.0 %*% out.fit$alpha.star.samples[sp.re.indx == q, ])
    #     
    #         } else {
    #           p.0.samples[, q, ] <- logit.inv(X.p.0 %*% out.fit$alpha.samples[sp.indx == q, ])
    #         }
    #         for (j in 1:nrow(X.p.0)) {
    #           for (k in 1:K.0[j]) {
    #             like.samples[q, j, k] <- mean(dbinom(y.mat.0[q, j, k], 1,
    #             			         p.0.samples[j, q, ] * out.pred$z.0.samples[, q, z.0.long.indx[j]]))
    #           }
    #         }
    #       }
    #     } else {
    #       like.samples <- matrix(NA, n.sp, nrow(X.p.0))
    #       for (q in 1:n.sp) {
    #         if (p.det.re > 0) {
    #           p.0.samples[, q, ] <- logit.inv(X.p.0 %*% out.fit$alpha.samples[sp.indx == q, ] +
    #     				      lambda.p.0 %*% out.fit$alpha.star.samples[sp.re.indx == q, ])
    #         } else {
    #           p.0.samples[, q, ] <- logit.inv(X.p.0 %*% out.fit$alpha.samples[sp.indx == q, ])
    #         }
    #         for (j in 1:nrow(X.p.0)) {
    #           like.samples[q, j] <- mean(dbinom(y.0[n.sp * (j - 1) + q], 1,  
    #             				p.0.samples[j, q, ] * 
    #             			        out.pred$z.0.samples[, q, z.0.long.indx[j]]))
    #         }
    #       }
    #     }
    #     apply(like.samples, 1, function(a) sum(log(a), na.rm = TRUE))
    #   }
    #   model.deviance <- -2 * model.deviance
    #   # Return objects from cross-validation
    #   out$k.fold.deviance <- model.deviance
    #   stopImplicitCluster()
    # }
  } # NNGP
  class(out) <- "sfMsAbund"
  out$run.time <- proc.time() - ptm
  out
}
