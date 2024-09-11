sfMsNMix <- function(abund.formula, det.formula, data, inits, priors,
		   tuning, cov.model = 'exponential', NNGP = TRUE,
		   n.neighbors = 15, search.type = 'cb', n.factors,
		   n.batch, batch.length, accept.rate = 0.43, family = 'Poisson',
		   n.omp.threads = 1, verbose = TRUE, n.report = 100,
		   n.burn = round(.10 * n.batch * batch.length), n.thin = 1,
		   n.chains = 1, ...){

  ptm <- proc.time()

  # Functions -----------------------------------------------------------
  logit <- function(theta, a = 0, b = 1) {log((theta-a)/(b-theta))}
  logit.inv <- function(z, a = 0, b = 1) {b-(b-a)/(1+exp(z))}

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
    stop("sfMsNMix is currently only implemented for NNGPs, not full Gaussian Processes. Please set NNGP = TRUE.")
  }
  if (missing(data)) {
    stop("data must be specified")
  }
  if (!is.list(data)) {
    stop("data must be a list")
  }
  names(data) <- tolower(names(data))
  if (!'y' %in% names(data)) {
    stop("count data y must be specified in data")
  }
  if (length(dim(data$y)) != 3) {
    stop("count data y must be a three-dimensional array with dimensions corresponding to species, sites, and replicates.")
  }
  y <- data$y
  sp.names <- attr(y, 'dimnames')[[1]]
  # Offset
  if ('offset' %in% names(data)) {
    offset <- data$offset
    if (length(offset) != ncol(y) & length(offset) != 1) {
      stop(paste("error: data$offset must be of length 1 or ", ncol(y), sep = ''))
    }
    if (length(offset) == 1) {
      offset <- rep(offset, ncol(y))
    }
  } else {
    offset <- rep(1, ncol(y))
  }
  if (!'abund.covs' %in% names(data)) {
    if (abund.formula == ~ 1) {
      if (verbose) {
        message("abundance covariates (abund.covs) not specified in data.\nAssuming intercept only abundance model.\n")
      }
      data$abund.covs <- matrix(1, dim(y)[2], 1)
    } else {
      stop("abund.covs must be specified in data for an abundance model with covariates")
    }
  }
  if (!'det.covs' %in% names(data)) {
    if (det.formula == ~ 1) {
      if (verbose) {
        message("detection covariates (det.covs) not specified in data.\nAssuming interept only detection model.\n")
      }
      data$det.covs <- list(int = matrix(1, dim(y)[2], dim(y)[3]))
    } else {
      stop("det.covs must be specified in data for a detection model with covariates")
    }
  }
  if (!is.list(data$det.covs)) {
    stop("det.covs must be a list of matrices, data frames, and/or vectors")
  }
  if (!'coords' %in% names(data)) {
    stop("coords must be specified in data for a spatial abundance model.")
  }
  if (!is.matrix(data$coords) & !is.data.frame(data$coords)) {
    stop("coords must be a matrix or data frame")
  }
  coords <- as.matrix(data$coords)
  if (missing(n.batch)) {
    stop("must specify number of MCMC batches")
  }
  if (missing(batch.length)) {
    stop("must specify length of each MCMC batch")
  }
  n.samples <- n.batch * batch.length
  if (n.burn > n.samples) {
    stop("n.burn must be less than n.samples")
  }
  if (n.thin > n.samples) {
    stop("n.thin must be less than n.samples")
  }
  # Check if n.burn, n.thin, and n.samples result in an integer and error if otherwise.
  if (((n.samples - n.burn) / n.thin) %% 1 != 0) {
    stop("the number of posterior samples to save ((n.samples - n.burn) / n.thin) is not a whole number. Please respecify the MCMC criteria such that the number of posterior samples saved is a whole number.")
  }
  if (missing(n.factors)) {
    stop("n.factors must be specified for a spatial factor N-mixture model")
  }

  if (!(family) %in% c('Poisson', 'NB')) {
    stop("family must be either 'Poisson' or 'NB'")
  }

  if (family == 'NB' & verbose) {
    message('**NOTE**: spatial negative binomial models can be difficult to\nestimate as they contain two forms of overdispersion. If experiencing\nvery poor mixing/convergence of MCMC chains (particularly kappa and phi),\nconsider using a spatial Poisson model or more informative\npriors on kappa or phi.\n')
  }

  # Neighbors and Ordering ----------------------------------------------
  if (NNGP) {
    u.search.type <- 2
    ## Order by x column. Could potentially allow this to be user defined.
    ord <- order(coords[,1])
    # Reorder everything to align with NN ordering
    y <- y[, ord, , drop = FALSE]
    coords <- coords[ord, , drop = FALSE]
    # Occupancy covariates
    data$abund.covs <- data$abund.covs[ord, , drop = FALSE]
    offset <- offset[ord]
    for (i in 1:length(data$det.covs)) {
      if (!is.null(dim(data$det.covs[[i]]))) {
        data$det.covs[[i]] <- data$det.covs[[i]][ord, , drop = FALSE]
      } else {
        data$det.covs[[i]] <- data$det.covs[[i]][ord]
      }
    }
  }
  # For later
  y.mat <- y

  # First subset detection covariates to only use those that are included in the analysis.
  data$det.covs <- data$det.covs[names(data$det.covs) %in% all.vars(det.formula)]
  # Null model support
  if (length(data$det.covs) == 0) {
    data$det.covs <- list(int = rep(1, dim(y)[2]))
  }
  # Make both covariates a data frame. Unlist is necessary for when factors
  # are supplied.
  data$det.covs <- data.frame(lapply(data$det.covs, function(a) unlist(c(a))))
  # Indicator of whether all det.covs are site level or not
  site.level.ind <- ifelse(nrow(data$det.covs) == ncol(y), TRUE, FALSE)
  data$abund.covs <- as.data.frame(data$abund.covs)

  # Check whether random effects are sent in as numeric, and
  # return error if they are.
  # Abundance -------------------------
  if (!is.null(findbars(abund.formula))) {
    abund.re.names <- unique(unlist(sapply(findbars(abund.formula), all.vars)))
    for (i in 1:length(abund.re.names)) {
      if (is(data$abund.covs[, abund.re.names[i]], 'factor')) {
        stop(paste("error: random effect variable ", abund.re.names[i], " specified as a factor. Random effect variables must be specified as numeric.", sep = ''))
      }
      if (is(data$abund.covs[, abund.re.names[i]], 'character')) {
        stop(paste("error: random effect variable ", abund.re.names[i], " specified as character. Random effect variables must be specified as numeric.", sep = ''))
      }
    }
  }
  # Detection -----------------------
  if (!is.null(findbars(det.formula))) {
    det.re.names <- unique(unlist(sapply(findbars(det.formula), all.vars)))
    for (i in 1:length(det.re.names)) {
      if (is(data$det.covs[, det.re.names[i]], 'factor')) {
        stop(paste("error: random effect variable ", det.re.names[i], " specified as a factor. Random effect variables must be specified as numeric.", sep = ''))
      }
      if (is(data$det.covs[, det.re.names[i]], 'character')) {
        stop(paste("error: random effect variable ", det.re.names[i], " specified as character. Random effect variables must be specified as numeric.", sep = ''))
      }
    }
  }

  # Checking missing values ---------------------------------------------
  # y -------------------------------
  y.na.test <- apply(y.mat, c(1, 2), function(a) sum(!is.na(a)))
  if (sum(y.na.test == 0) > 0) {
    stop("error: some sites in y have all missing detection histories. Remove these sites from all objects in the 'data' argument, then use 'predict' to obtain predictions at these locations if desired.")
  }
  # abund.covs ------------------------
  if (sum(is.na(data$abund.covs)) != 0) {
    stop("error: missing values in abund.covs. Please remove these sites from all objects in data or somehow replace the NA values with non-missing values (e.g., mean imputation).")
  }
  # det.covs ------------------------
  if (!site.level.ind) {
    for (i in 1:ncol(data$det.covs)) {
      if (sum(is.na(data$det.covs[, i])) > sum(is.na(y.mat[1, , ]))) {
        stop("error: some elements in det.covs have missing values where there is an observed data value in y. Please either replace the NA values in det.covs with non-missing values (e.g., mean imputation) or set the corresponding values in y to NA where the covariate is missing.")
      }
    }
    # Misalignment between y and det.covs
    y.missing <- which(is.na(y[1, , ]))
    det.covs.missing <- lapply(data$det.covs, function(a) which(is.na(a)))
    for (i in 1:length(det.covs.missing)) {
      tmp.indx <- !(y.missing %in% det.covs.missing[[i]])
      if (sum(tmp.indx) > 0) {
        if (i == 1 & verbose) {
          message("There are missing values in data$y with corresponding non-missing values in data$det.covs.\nRemoving these site/replicate combinations for fitting the model.\n")
        }
        data$det.covs[y.missing, i] <- NA
      }
    }
  }
  if (site.level.ind) {
    if (sum(is.na(data$det.covs)) != 0) {
      stop("missing values in site-level det.covs. Please remove these sites from all objects in data or somehow replace the NA values with non-missing values (e.g., mean imputation).")
    }
  }

  # Remove missing values from det.covs in order to ensure formula parsing
  # works when random slopes are provided.
  tmp <- apply(data$det.covs, 1, function (a) sum(is.na(a)))
  data$det.covs <- as.data.frame(data$det.covs[tmp == 0, , drop = FALSE])

  # Formula -------------------------------------------------------------
  # Abundance -----------------------
  if (missing(abund.formula)) {
    stop("error: abund.formula must be specified")
  }

  if (is(abund.formula, 'formula')) {
    tmp <- parseFormula(abund.formula, data$abund.covs)
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
  re.level.names <- lapply(data$abund.covs[, x.re.names, drop = FALSE],
      		     function (a) sort(unique(a)))
  x.re.names <- x.random.names

  # Detection -----------------------
  if (missing(det.formula)) {
    stop("error: det.formula must be specified")
  }

  if (is(det.formula, 'formula')) {
    tmp <- parseFormula(det.formula, data$det.covs)
    X.p <- as.matrix(tmp[[1]])
    X.p.re <- as.matrix(tmp[[4]])
    x.p.re.names <- colnames(X.p.re)
    x.p.names <- tmp[[2]]
    X.p.random <- as.matrix(tmp[[5]])
    x.p.random.names <- colnames(X.p.random)
  } else {
    stop("error: det.formula is misspecified")
  }
  p.re.level.names <- lapply(data$det.covs[, x.p.re.names, drop = FALSE],
      		       function (a) sort(unique(a)))
  x.p.re.names <- x.p.random.names

  # Extract data from inputs --------------------------------------------
  # Number of species
  n.sp <- dim(y)[1]
  # Number of latent factors
  q <- n.factors
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
  # Number of sites
  J <- nrow(X)
  # Number of repeat visits
  # Note this assumes equivalent detection histories for all species.
  # May want to change this at some point.
  n.rep <- apply(y.mat[1, , , drop = FALSE], 2, function(a) sum(!is.na(a)))
  K.max <- dim(y.mat)[3]
  # Because I like K better than n.rep
  K <- n.rep

  # Get indices to map N to y -------------------------------------------
  N.long.indx <- rep(1:J, dim(y.mat)[3])
  N.long.indx <- N.long.indx[!is.na(c(y.mat[1, , ]))]
  # Subtract 1 for indices in C
  N.long.indx <- N.long.indx - 1
  # y is stored in the following order: species, site, visit
  y <- c(y)
  # Assumes the missing data are constant across species, which seems likely,
  # but may eventually need some updating.
  names.long <- which(!is.na(c(y.mat[1, , ])))
  # Make necessary adjustment for site-level covariates only
  if (nrow(X.p) == J) {
    X.p <- X.p[N.long.indx + 1, , drop = FALSE]
    X.p.re <- X.p.re[N.long.indx + 1, , drop = FALSE]
    X.p.random <- X.p.random[N.long.indx + 1, , drop = FALSE]
  }
  # Remove missing observations when the covariate data are available but
  # there are missing detection-nondetection data.
  if (nrow(X.p) == length(y) / n.sp) {
    X.p <- X.p[!is.na(y.mat[1, , ]), , drop = FALSE]
  }
  if (nrow(X.p.re) == length(y) / n.sp & p.det.re > 0) {
    X.p.re <- X.p.re[!is.na(y.mat[1, , ]), , drop = FALSE]
    X.p.random <- X.p.random[!is.na(y.mat[1, , ]), , drop = FALSE]
  }
  y <- y[!is.na(y)]
  # Number of pseudoreplicates
  n.obs <- nrow(X.p)

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

  # alpha.comm -----------------------
  if ("alpha.comm.normal" %in% names(priors)) {
    if (!is.list(priors$alpha.comm.normal) | length(priors$alpha.comm.normal) != 2) {
      stop("error: alpha.comm.normal must be a list of length 2")
    }
    mu.alpha.comm <- priors$alpha.comm.normal[[1]]
    sigma.alpha.comm <- priors$alpha.comm.normal[[2]]
    if (length(mu.alpha.comm) != p.det & length(mu.alpha.comm) != 1) {
      if (p.det == 1) {
        stop(paste("error: alpha.comm.normal[[1]] must be a vector of length ",
        	     p.det, " with elements corresponding to alpha.comms' mean", sep = ""))
      } else {
        stop(paste("error: alpha.comm.normal[[1]] must be a vector of length ",
        	     p.det, " or 1 with elements corresponding to alpha.comms' mean", sep = ""))
      }
    }
    if (length(sigma.alpha.comm) != p.det & length(sigma.alpha.comm) != 1) {
      if (p.det == 1) {
        stop(paste("error: alpha.comm.normal[[2]] must be a vector of length ",
      	   p.det, " with elements corresponding to alpha.comms' variance", sep = ""))
      } else {
        stop(paste("error: alpha.comm.normal[[2]] must be a vector of length ",
      	   p.det, " or 1 with elements corresponding to alpha.comms' variance", sep = ""))
      }
    }
    if (length(sigma.alpha.comm) != p.det) {
      sigma.alpha.comm <- rep(sigma.alpha.comm, p.det)
    }
    if (length(mu.alpha.comm) != p.det) {
      mu.alpha.comm <- rep(mu.alpha.comm, p.det)
    }
    Sigma.alpha.comm <- sigma.alpha.comm * diag(p.det)
  } else {
    if (verbose) {
      message("No prior specified for alpha.comm.normal.\nSetting prior mean to 0 and prior variance to 2.72\n")
    }
    mu.alpha.comm <- rep(0, p.det)
    sigma.alpha.comm <- rep(2.72, p.det)
    Sigma.alpha.comm <- diag(p.det) * 2.72
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

  # tau.sq.alpha -----------------------
  if ("tau.sq.alpha.ig" %in% names(priors)) {
    if (!is.list(priors$tau.sq.alpha.ig) | length(priors$tau.sq.alpha.ig) != 2) {
      stop("error: tau.sq.alpha.ig must be a list of length 2")
    }
    tau.sq.alpha.a <- priors$tau.sq.alpha.ig[[1]]
    tau.sq.alpha.b <- priors$tau.sq.alpha.ig[[2]]
    if (length(tau.sq.alpha.a) != p.det & length(tau.sq.alpha.a) != 1) {
      if (p.det == 1) {
        stop(paste("error: tau.sq.alpha.ig[[1]] must be a vector of length ",
      	   p.det, " with elements corresponding to tau.sq.alphas' shape", sep = ""))
      } else {
        stop(paste("error: tau.sq.alpha.ig[[1]] must be a vector of length ",
      	   p.det, " or 1 with elements corresponding to tau.sq.alphas' shape", sep = ""))
      }
    }
    if (length(tau.sq.alpha.b) != p.det & length(tau.sq.alpha.b) != 1) {
      if (p.det == 1) {
        stop(paste("error: tau.sq.alpha.ig[[2]] must be a vector of length ",
      	   p.det, " with elements corresponding to tau.sq.alphas' scale", sep = ""))
      } else {
        stop(paste("error: tau.sq.alpha.ig[[2]] must be a vector of length ",
      	   p.det, " or 1 with elements corresponding to tau.sq.alphas' scale", sep = ""))
      }
    }
    if (length(tau.sq.alpha.a) != p.det) {
      tau.sq.alpha.a <- rep(tau.sq.alpha.a, p.det)
    }
    if (length(tau.sq.alpha.b) != p.det) {
      tau.sq.alpha.b <- rep(tau.sq.alpha.b, p.det)
    }
  } else {
    if (verbose) {
      message("No prior specified for tau.sq.alpha.ig.\nSetting prior shape to 0.1 and prior scale to 0.1\n")
    }
    tau.sq.alpha.a <- rep(0.1, p.det)
    tau.sq.alpha.b <- rep(0.1, p.det)
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
  } else {
    kappa.a <- rep(0, n.sp)
    kappa.b <- rep(0, n.sp)
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

  # Starting values -----------------------------------------------------
  if (missing(inits)) {
    inits <- list()
  }
  names(inits) <- tolower(names(inits))
  # N -------------------------------
  if ("n" %in% names(inits)) {
    N.inits <- inits$n
    if (!is.matrix(N.inits)) {
      stop(paste("error: initial values for N must be a matrix with dimensions ",
      	   n.sp, " x ", J, sep = ""))
    }
    if (nrow(N.inits) != n.sp | ncol(N.inits) != J) {
      stop(paste("error: initial values for N must be a matrix with dimensions ",
      	   n.sp, " x ", J, sep = ""))
    }
    # Reorder the user supplied inits values for NNGP models
    if (NNGP) {
      N.inits <- N.inits[, ord]
    }
    N.test <- apply(y.mat, c(1, 2), max, na.rm = TRUE)
    init.test <- sum(N.inits < N.test)
    if (init.test > 0) {
      stop("error: initial values for latent abundance (N) are invalid. Please re-specify inits$N so initial values are 1 if the species is observed at that site.")
    }
  } else {
    # In correct order since y is already reordered.
    N.inits <- apply(y.mat, c(1, 2), max, na.rm = TRUE)
    if (verbose) {
      message("N is not specified in initial values.\nSetting initial values based on observed data\n")
    }
  }
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
    beta.comm.inits <- rnorm(p.abund, 0, 1)
    if (verbose) {
      message('beta.comm is not specified in initial values.\nSetting initial values to random values from a standard normal distribution\n')
    }
  }
  # alpha.comm -----------------------
  if ("alpha.comm" %in% names(inits)) {
    alpha.comm.inits <- inits[["alpha.comm"]]
    if (length(alpha.comm.inits) != p.det & length(alpha.comm.inits) != 1) {
      if (p.det == 1) {
        stop(paste("error: initial values for alpha.comm must be of length ", p.det,
      	   sep = ""))
      } else {
        stop(paste("error: initial values for alpha.comm must be of length ", p.det,
      	   , " or 1", sep = ""))
      }
    }
    if (length(alpha.comm.inits) != p.det) {
      alpha.comm.inits <- rep(alpha.comm.inits, p.det)
    }
  } else {
    alpha.comm.inits <- rnorm(p.det, 0, 1)
    if (verbose) {
      message('alpha.comm is not specified in initial values.\nSetting initial values to random values from a standard normal distribution\n')
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
    tau.sq.beta.inits <- runif(p.abund, 0.05, 1)
    if (verbose) {
      message('tau.sq.beta is not specified in initial values.\nSetting initial values to random values between 0.05 and 1\n')
    }
  }
  # tau.sq.alpha -----------------------
  if ("tau.sq.alpha" %in% names(inits)) {
    tau.sq.alpha.inits <- inits[["tau.sq.alpha"]]
    if (length(tau.sq.alpha.inits) != p.det & length(tau.sq.alpha.inits) != 1) {
      if (p.det == 1) {
        stop(paste("error: initial values for tau.sq.alpha must be of length ", p.det,
      	   sep = ""))
      } else {
        stop(paste("error: initial values for tau.sq.alpha must be of length ", p.det,
      	   " or 1", sep = ""))
      }
    }
    if (length(tau.sq.alpha.inits) != p.det) {
      tau.sq.alpha.inits <- rep(tau.sq.alpha.inits, p.det)
    }
  } else {
    tau.sq.alpha.inits <- runif(p.det, 0.05, 1)
    if (verbose) {
      message('tau.sq.alpha is not specified in initial values.\nSetting to initial values to random values between 0.05 and 1\n')
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
  # alpha ----------------------------
  if ("alpha" %in% names(inits)) {
    alpha.inits <- inits[["alpha"]]
    if (is.matrix(alpha.inits)) {
      if (ncol(alpha.inits) != p.det | nrow(alpha.inits) != n.sp) {
        stop(paste("error: initial values for alpha must be a matrix with dimensions ",
        	   n.sp, "x", p.det, " or a single numeric value", sep = ""))
      }
    }
    if (!is.matrix(alpha.inits) & length(alpha.inits) != 1) {
      stop(paste("error: initial values for alpha must be a matrix with dimensions ",
      	   n.sp, " x ", p.det, " or a single numeric value", sep = ""))
    }
    if (length(alpha.inits) == 1) {
      alpha.inits <- matrix(alpha.inits, n.sp, p.det)
    }
  } else {
    alpha.inits <- matrix(rnorm(n.sp * p.det, alpha.comm.inits, sqrt(tau.sq.alpha.inits)), n.sp, p.det)
    if (verbose) {
      message('alpha is not specified in initial values.\nSetting initial values to random values from the community-level normal distribution\n')
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
    beta.star.inits <- rnorm(n.abund.re, 0, sqrt(sigma.sq.mu.inits[beta.star.indx + 1]))
    # Starting values for all species
    beta.star.inits <- rep(beta.star.inits, n.sp)
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
      sigma.sq.p.inits <- runif(p.det.re, 0.05, 1)
      if (verbose) {
        message("sigma.sq.p is not specified in initial values.\nSetting initial values to random values between 0.05 and 1\n")
      }
    }
    alpha.star.indx <- rep(0:(p.det.re - 1), n.det.re.long)
    alpha.star.inits <- rnorm(n.det.re, 0, sqrt(sigma.sq.p.inits[alpha.star.indx + 1]))
    alpha.star.inits <- rep(alpha.star.inits, n.sp)
  } else {
    sigma.sq.p.inits <- 0
    alpha.star.indx <- 0
    alpha.star.inits <- 0
  }
  # kappa -----------------------------
  # ORDER: a length n.sp vector ordered by species in the detection-nondetection data.
  if (family == 'NB') {
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
  } else {
    kappa.inits <- rep(0, n.sp)
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
  # Keep this just for consistency
  sigma.sq.tuning <- rep(0, q)
  if (missing(tuning)) {
    beta.tuning <- rep(1, p.abund * n.sp)
    beta.star.tuning <- rep(1, n.abund.re * n.sp)
    alpha.tuning <- rep(1, p.det * n.sp)
    alpha.star.tuning <- rep(1, n.det.re * n.sp)
    kappa.tuning <- rep(1, n.sp)
    phi.tuning <- rep(1, q)
    if (cov.model == 'matern') {
      nu.tuning <- rep(1, q)
    } else {
      nu.tuning <- NULL
    }
    w.tuning <- rep(1, J * q)
    lambda.tuning <- rep(1, n.sp * q)
  } else {
    names(tuning) <- tolower(names(tuning))
    # beta ---------------------------
    if(!"beta" %in% names(tuning)) {
      stop("error: beta must be specified in tuning value list")
    }
    beta.tuning <- tuning$beta
    if (length(beta.tuning) != 1 & length(beta.tuning) != p.abund * n.sp) {
      stop(paste("error: beta tuning must be a single value or a vector of length ",
        	 p.abund * n.sp, sep = ''))
    }
    if (length(beta.tuning) == 1) {
      beta.tuning <- rep(beta.tuning, p.abund * n.sp)
    }
    # alpha ---------------------------
    if(!"alpha" %in% names(tuning)) {
      stop("error: alpha must be specified in tuning value list")
    }
    alpha.tuning <- tuning$alpha
    if (length(alpha.tuning) != 1 & length(alpha.tuning) != p.det * n.sp) {
      stop(paste("error: alpha tuning must be a single value or a vector of length ",
        	 p.det * n.sp, sep = ''))
    }
    if (length(alpha.tuning) == 1) {
      alpha.tuning <- rep(alpha.tuning, p.det * n.sp)
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
      beta.star.tuning <- rep(beta.star.tuning, n.abund.re * n.sp)
    } else {
      beta.star.tuning <- NULL
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
      alpha.star.tuning <- rep(alpha.star.tuning, n.det.re * n.sp)
    } else {
      alpha.star.tuning <- NULL
    }
    # kappa ---------------------------
    if (family == 'NB') {
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
    } else {
      kappa.tuning <- NULL
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
    } else {
      nu.tuning <- NULL
    }
    # w ---------------------------
    if(!"w" %in% names(tuning)) {
      stop("error: w must be specified in tuning value list")
    }
    w.tuning <- tuning$w
    if (length(w.tuning) != 1 & length(w.tuning) != J * q) {
      stop(paste("error: w tuning must be a single value or a vector of length ",
        	 J * q, sep = ''))
    }
    if (length(w.tuning) == 1) {
      w.tuning <- rep(w.tuning, J * q)
    }
    # lambda ---------------------------
    if(!"lambda" %in% names(tuning)) {
      stop("error: lambda must be specified in tuning value list")
    }
    lambda.tuning <- tuning$lambda
    if (length(lambda.tuning) != 1 & length(lambda.tuning) != n.sp * q) {
      stop(paste("error: lambda tuning must be a single value or a vector of length ",
        	 n.sp * q, sep = ''))
    }
    if (length(lambda.tuning) == 1) {
      lambda.tuning <- rep(lambda.tuning, n.sp * q)
    }
  }
  tuning.c <- log(c(beta.tuning, alpha.tuning, beta.star.tuning,
		    alpha.star.tuning, sigma.sq.tuning, phi.tuning, nu.tuning,
		    lambda.tuning, w.tuning, kappa.tuning))

  # Get max y values for N update -----------------------------------------
  y.max <- apply(y.mat, c(1, 2), max, na.rm = TRUE)

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

  curr.chain <- 1
  if (!NNGP) {

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
    storage.mode(N.inits) <- "double"
    storage.mode(coords) <- "double"
    storage.mode(X.p) <- "double"
    storage.mode(X) <- "double"
    storage.mode(offset) <- "double"
    storage.mode(K) <- "double"
    storage.mode(y.max) <- "double"
    consts <- c(n.sp, J, n.obs, p.abund, p.abund.re, n.abund.re, p.det, p.det.re, n.det.re, q)
    storage.mode(consts) <- "integer"
    storage.mode(beta.inits) <- "double"
    storage.mode(alpha.inits) <- "double"
    storage.mode(kappa.inits) <- "double"
    storage.mode(beta.comm.inits) <- "double"
    storage.mode(alpha.comm.inits) <- "double"
    storage.mode(tau.sq.beta.inits) <- "double"
    storage.mode(tau.sq.alpha.inits) <- "double"
    storage.mode(lambda.inits) <- "double"
    storage.mode(w.inits) <- "double"
    storage.mode(phi.inits) <- "double"
    storage.mode(nu.inits) <- "double"
    storage.mode(N.long.indx) <- "integer"
    storage.mode(mu.beta.comm) <- "double"
    storage.mode(Sigma.beta.comm) <- "double"
    storage.mode(mu.alpha.comm) <- "double"
    storage.mode(Sigma.alpha.comm) <- "double"
    storage.mode(kappa.a) <- "double"
    storage.mode(kappa.b) <- "double"
    storage.mode(tau.sq.beta.a) <- "double"
    storage.mode(tau.sq.beta.b) <- "double"
    storage.mode(tau.sq.alpha.a) <- "double"
    storage.mode(tau.sq.alpha.b) <- "double"
    spatial.priors <- c(phi.a, phi.b, nu.a, nu.b)
    storage.mode(spatial.priors) <- "double"
    storage.mode(n.batch) <- "integer"
    storage.mode(batch.length) <- "integer"
    batch.info <- c(n.batch, batch.length)
    storage.mode(batch.info) <- 'integer'
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

    # Fit the model -------------------------------------------------------
    out.tmp <- list()
    out <- list()
    for (i in 1:n.chains) {
      # Change initial values if i > 1
      if ((i > 1) & (!fix.inits)) {
        beta.comm.inits <- rnorm(p.abund, 0, 1)
        alpha.comm.inits <- rnorm(p.det, 0, 1)
        tau.sq.beta.inits <- runif(p.abund, 0.05, 1)
        tau.sq.alpha.inits <- runif(p.det, 0.05, 1)
        beta.inits <- matrix(rnorm(n.sp * p.abund, beta.comm.inits,
              		     sqrt(tau.sq.beta.inits)), n.sp, p.abund)
        alpha.inits <- matrix(rnorm(n.sp * p.det, alpha.comm.inits,
              		      sqrt(tau.sq.alpha.inits)), n.sp, p.det)
        if (family == 'NB') {
          kappa.inits <- runif(n.sp, kappa.a, kappa.b)
        }
        lambda.inits <- matrix(0, n.sp, q)
        diag(lambda.inits) <- 1
        lambda.inits[lower.tri(lambda.inits)] <- rnorm(sum(lower.tri(lambda.inits)))
        lambda.inits <- c(lambda.inits)
        phi.inits <- runif(q, phi.a, phi.b)
        if (cov.model == 'matern') {
          nu.inits <- runif(q, nu.a, nu.b)
        }
        if (p.abund.re > 0) {
          sigma.sq.mu.inits <- runif(p.abund.re, 0.05, 1)
          beta.star.inits <- rnorm(n.abund.re, 0, sqrt(sigma.sq.mu.inits[beta.star.indx + 1]))
          beta.star.inits <- rep(beta.star.inits, n.sp)
        }
        if (p.det.re > 0) {
          sigma.sq.p.inits <- runif(p.det.re, 0.05, 1)
          alpha.star.inits <- rnorm(n.det.re, 0, sqrt(sigma.sq.p.inits[alpha.star.indx + 1]))
          alpha.star.inits <- rep(alpha.star.inits, n.sp)
        }
      }
      storage.mode(chain.info) <- "integer"
      out.tmp[[i]] <- .Call("sfMsNMixNNGP", y, X, X.p, coords, X.re, X.p.re,
          		    X.random, X.p.random, y.max, consts,
      	                    n.abund.re.long, n.det.re.long,
          	            n.neighbors, nn.indx, nn.indx.lu, u.indx, u.indx.lu, ui.indx,
        		    beta.inits, alpha.inits, kappa.inits, N.inits, beta.comm.inits,
        	            alpha.comm.inits, phi.inits, lambda.inits,
        		    nu.inits, w.inits, tau.sq.beta.inits, tau.sq.alpha.inits,
        		    sigma.sq.mu.inits, sigma.sq.p.inits,
      	                    beta.star.inits, alpha.star.inits, N.long.indx,
        		    beta.star.indx, beta.level.indx, alpha.star.indx,
        		    alpha.level.indx, mu.beta.comm, mu.alpha.comm,
        		    Sigma.beta.comm, Sigma.alpha.comm, kappa.a,
          		    kappa.b, tau.sq.beta.a, tau.sq.beta.b, tau.sq.alpha.a,
        	            tau.sq.alpha.b, spatial.priors,
        		    sigma.sq.mu.a, sigma.sq.mu.b, sigma.sq.p.a, sigma.sq.p.b,
      	                    tuning.c, cov.model.indx, batch.info,
        		    accept.rate, n.omp.threads,
      	                    verbose, n.report, samples.info, chain.info, family.c, offset)
      chain.info[1] <- chain.info[1] + 1
    }
    # Calculate R-Hat ---------------
    out <- list()
    out$rhat <- list()
    if (n.chains > 1) {
      # as.vector removes the "Upper CI" when there is only 1 variable.
      out$rhat$beta.comm <- as.vector(gelman.diag(mcmc.list(lapply(out.tmp, function(a)
      					      mcmc(t(a$beta.comm.samples)))),
      			     autoburnin = FALSE, multivariate = FALSE)$psrf[, 2])
      out$rhat$alpha.comm <- as.vector(gelman.diag(mcmc.list(lapply(out.tmp, function(a)
      					      mcmc(t(a$alpha.comm.samples)))),
      			     autoburnin = FALSE, multivariate = FALSE)$psrf[, 2])
      out$rhat$tau.sq.beta <- as.vector(gelman.diag(mcmc.list(lapply(out.tmp, function(a)
      					      mcmc(t(a$tau.sq.beta.samples)))),
      			     autoburnin = FALSE, multivariate = FALSE)$psrf[, 2])
      out$rhat$tau.sq.alpha <- as.vector(gelman.diag(mcmc.list(lapply(out.tmp, function(a)
      					      mcmc(t(a$tau.sq.alpha.samples)))),
      			     autoburnin = FALSE, multivariate = FALSE)$psrf[, 2])
      out$rhat$beta <- as.vector(gelman.diag(mcmc.list(lapply(out.tmp, function(a)
      					         mcmc(t(a$beta.samples)))),
      			     autoburnin = FALSE, multivariate = FALSE)$psrf[, 2])
      out$rhat$alpha <- as.vector(gelman.diag(mcmc.list(lapply(out.tmp, function(a)
      					      mcmc(t(a$alpha.samples)))),
      			      autoburnin = FALSE, multivariate = FALSE)$psrf[, 2])
      out$rhat$theta <- gelman.diag(mcmc.list(lapply(out.tmp, function(a)
      					      mcmc(t(a$theta.samples)))),
      			      autoburnin = FALSE, multivariate = FALSE)$psrf[, 2]
      lambda.mat <- matrix(lambda.inits, n.sp, q)
      out$rhat$lambda.lower.tri <- as.vector(gelman.diag(mcmc.list(lapply(out.tmp, function(a)
          					       mcmc(t(a$lambda.samples[c(lower.tri(lambda.mat)), ])))),
          					       autoburnin = FALSE, multivariate = FALSE)$psrf[, 2])
      if (family == 'NB') {
        out$rhat$kappa <- as.vector(gelman.diag(mcmc.list(lapply(out.tmp, function(a)
      					      mcmc(t(a$kappa.samples)))),
      			     autoburnin = FALSE, multivariate = FALSE)$psrf[, 2])
      }
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
    } else {
      out$rhat$beta.comm <- rep(NA, p.abund)
      out$rhat$alpha.comm <- rep(NA, p.det)
      out$rhat$tau.sq.beta <- rep(NA, p.abund)
      out$rhat$tau.sq.alpha <- rep(NA, p.det)
      out$rhat$beta <- rep(NA, p.abund * n.sp)
      out$rhat$alpha <- rep(NA, p.det * n.sp)
      out$rhat$theta <- rep(NA, ifelse(cov.model == 'matern', 2 * q, q))
      out$rhat$kappa <- rep(NA, n.sp)
      if (p.det.re > 0) {
        out$rhat$sigma.sq.p <- rep(NA, p.det.re)
      }
      if (p.abund.re > 0) {
        out$rhat$sigma.sq.mu <- rep(NA, p.abund.re)
      }
    }
    # Put everything into MCMC objects
    out$beta.comm.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$beta.comm.samples))))
    colnames(out$beta.comm.samples) <- x.names
    out$alpha.comm.samples <- mcmc(do.call(rbind,
      				lapply(out.tmp, function(a) t(a$alpha.comm.samples))))
    colnames(out$alpha.comm.samples) <- x.p.names
    out$tau.sq.beta.samples <- mcmc(do.call(rbind,
      				lapply(out.tmp, function(a) t(a$tau.sq.beta.samples))))
    colnames(out$tau.sq.beta.samples) <- x.names
    out$tau.sq.alpha.samples <- mcmc(do.call(rbind,
      				lapply(out.tmp, function(a) t(a$tau.sq.alpha.samples))))
    colnames(out$tau.sq.alpha.samples) <- x.p.names

    if (is.null(sp.names)) {
      sp.names <- paste('sp', 1:n.sp, sep = '')
    }
    coef.names <- paste(rep(x.names, each = n.sp), sp.names, sep = '-')
    out$beta.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$beta.samples))))
    colnames(out$beta.samples) <- coef.names
    out$alpha.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$alpha.samples))))
    coef.names.det <- paste(rep(x.p.names, each = n.sp), sp.names, sep = '-')
    colnames(out$alpha.samples) <- coef.names.det
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
    if (family == 'NB') {
      out$kappa.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$kappa.samples))))
      colnames(out$kappa.samples) <- paste('kappa', sp.names, sep = '-')
    }
    out$w.samples <- do.call(abind, lapply(out.tmp, function(a) array(a$w.samples,
      								dim = c(q, J, n.post.samples))))
    out$w.samples <- out$w.samples[, order(ord), , drop = FALSE]
    out$w.samples <- aperm(out$w.samples, c(3, 1, 2))
    out$N.samples <- do.call(abind, lapply(out.tmp, function(a) array(a$N.samples,
      								dim = c(n.sp, J, n.post.samples))))
    out$N.samples <- out$N.samples[, order(ord), , drop = FALSE]
    out$N.samples <- aperm(out$N.samples, c(3, 1, 2))
    out$mu.samples <- do.call(abind, lapply(out.tmp, function(a) array(a$mu.samples,
      								dim = c(n.sp, J, n.post.samples))))
    out$mu.samples <- out$mu.samples[, order(ord), , drop = FALSE]
    out$mu.samples <- aperm(out$mu.samples, c(3, 1, 2))
    if (p.det.re > 0) {
      out$sigma.sq.p.samples <- mcmc(
        do.call(rbind, lapply(out.tmp, function(a) t(a$sigma.sq.p.samples))))
      colnames(out$sigma.sq.p.samples) <- x.p.re.names
      out$alpha.star.samples <- mcmc(
        do.call(rbind, lapply(out.tmp, function(a) t(a$alpha.star.samples))))
      tmp.names <- unlist(p.re.level.names)
      alpha.star.names <- paste(rep(x.p.re.names, n.det.re.long), tmp.names, sep = '-')
      alpha.star.names <- paste(alpha.star.names, rep(sp.names, each = n.det.re), sep = '-')
      colnames(out$alpha.star.samples) <- alpha.star.names
      out$p.re.level.names <- p.re.level.names
    }
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
    out$ESS$alpha.comm <- effectiveSize(out$alpha.comm.samples)
    out$ESS$tau.sq.beta <- effectiveSize(out$tau.sq.beta.samples)
    out$ESS$tau.sq.alpha <- effectiveSize(out$tau.sq.alpha.samples)
    out$ESS$beta <- effectiveSize(out$beta.samples)
    out$ESS$alpha <- effectiveSize(out$alpha.samples)
    out$ESS$theta <- effectiveSize(out$theta.samples)
    out$ESS$lambda <- effectiveSize(out$lambda.samples)
    if (family == 'NB') {
      out$ESS$kappa <- effectiveSize(out$kappa.samples)
    }
    if (p.det.re > 0) {
      out$ESS$sigma.sq.p <- effectiveSize(out$sigma.sq.p.samples)
    }
    if (p.abund.re > 0) {
      out$ESS$sigma.sq.mu <- effectiveSize(out$sigma.sq.mu.samples)
    }
    out$X <- X[order(ord), , drop = FALSE]
    tmp <- matrix(NA, J * K.max, p.det)
    tmp[names.long, ] <- X.p
    tmp <- array(tmp, dim = c(J, K.max, p.det))
    tmp <- tmp[order(ord), , ]
    out$X.p <- matrix(tmp, J * K.max, p.det)
    out$X.p <- out$X.p[apply(out$X.p, 1, function(a) sum(is.na(a))) == 0, , drop = FALSE]
    colnames(out$X.p) <- x.p.names
    tmp <- matrix(NA, J * K.max, p.det.re)
    tmp[names.long, ] <- X.p.re
    tmp <- array(tmp, dim = c(J, K.max, p.det.re))
    tmp <- tmp[order(ord), , ]
    out$X.p.re <- matrix(tmp, J * K.max, p.det.re)
    out$X.p.re <- out$X.p.re[apply(out$X.p.re, 1, function(a) sum(is.na(a))) == 0, , drop = FALSE]
    colnames(out$X.p.re) <- colnames(X.p.re)
    tmp <- matrix(NA, J * K.max, p.det.re)
    tmp[names.long, ] <- X.p.random
    tmp <- array(tmp, dim = c(J, K.max, p.det.re))
    tmp <- tmp[order(ord), , ]
    out$X.p.random <- matrix(tmp, J * K.max, p.det.re)
    out$X.p.random <- out$X.p.random[apply(out$X.p.random, 1, function(a) sum(is.na(a))) == 0, , drop = FALSE]
    colnames(out$X.p.random) <- x.p.random.names
    tmp <- matrix(NA, J * K.max, n.det.re)
    out$X.re <- X.re[order(ord), , drop = FALSE]
    out$X.random <- X.random[order(ord), , drop = FALSE]
    out$y <- y.mat[, order(ord), , drop = FALSE]
    out$offset <- offset[order(ord)]
    out$call <- cl
    out$n.samples <- n.samples
    out$x.names <- x.names
    out$sp.names <- sp.names
    out$x.p.names <- x.p.names
    out$n.post <- n.post.samples
    out$n.thin <- n.thin
    out$n.burn <- n.burn
    out$n.chains <- n.chains
    out$re.cols <- re.cols
    out$re.det.cols <- re.det.cols
    out$theta.names <- theta.names
    out$type <- "NNGP"
    out$coords <- coords[order(ord), ]
    out$cov.model.indx <- cov.model.indx
    out$n.neighbors <- n.neighbors
    out$dist <- family
    out$q <- q
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
  } # NNGP
  class(out) <- "sfMsNMix"
  out$run.time <- proc.time() - ptm
  out
}
