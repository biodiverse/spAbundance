spNMix <- function(abund.formula, det.formula, data, inits, priors, tuning,
                   cov.model = 'exponential', NNGP = TRUE,
                   n.neighbors = 15, search.type = 'cb',
                   n.batch, batch.length, accept.rate = 0.43, family = 'Poisson',
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
    stop("data must be specified")
  }
  if (!is.list(data)) {
    stop("data must be a list")
  }
  names(data) <- tolower(names(data))
  if (missing(abund.formula)) {
    stop("abund.formula must be specified")
  }
  if (missing(det.formula)) {
    stop("error: det.formula must be specified")
  }
  if (!'y' %in% names(data)) {
    stop("error: data y must be specified in data")
  }
  y <- as.matrix(data$y)
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
  if (!'abund.covs' %in% names(data)) {
    if (abund.formula == ~ 1) {
      if (verbose) {
        message("abundance covariates (abund.covs) not specified in data.\nAssuming intercept only abundance model.\n")
      }
      data$abund.covs <- matrix(1, dim(y)[1], 1)
    } else {
      stop("error: abund.covs must be specified in data for an abundance model with covariates")
    }
  }
  if (!is.matrix(data$abund.covs) & !is.data.frame(data$abund.covs)) {
    stop("error: abund.covs must be a matrix or data frame")
  }
  if (sum(is.na(data$abund.covs)) > 0) {
    stop("error: missing covariate values in data$abund.covs. Remove these sites from all data or impute non-missing values.")
  }
  if (!'det.covs' %in% names(data)) {
    if (det.formula == ~ 1) {
      if (verbose) {
        message("detection covariates (det.covs) not specified in data.\nAssuming interept only detection model.\n")
      }
      data$det.covs <- list(int = rep(1, dim(y)[1]))
    } else {
      stop("error: det.covs must be specified in data for a detection model with covariates")
    }
  }
  if (!is.list(data$det.covs)) {
    stop("error: det.covs must be a list of matrices, data frames, and/or vectors")
  }
  if (!'coords' %in% names(data)) {
    stop("error: coords must be specified in data for a spatial N-mixture model.")
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
  # Check if n.burn, n.thin, and n.samples result in an integer and error if otherwise.
  if (((n.samples - n.burn) / n.thin) %% 1 != 0) {
    stop("the number of posterior samples to save ((n.samples - n.burn) / n.thin) is not a whole number. Please respecify the MCMC criteria such that the number of posterior samples saved is a whole number.")
  }

  if (!(family) %in% c('Poisson', 'NB')) {
    stop("family must be either 'Poisson' or 'NB'")
  }

  if (family == 'NB' & verbose) {
    message('**NOTE**: spatial negative binomial models can be difficult to\nestimate as they contain two forms of overdispersion.\nIf experiencing very poor mixing/convergence of MCMC chains (particularly kappa and phi),\nconsider using a spatial Poisson model or more informative priors on kappa or phi.\n')
  }

  # Neighbors and Ordering ----------------------------------------------
  if (NNGP) {
    u.search.type <- 2
    ## Order by x column. Could potentially allow this to be user defined.
    ord <- order(coords[,1])
    # Reorder everything to align with NN ordering
    y <- y[ord, , drop = FALSE]
    coords <- coords[ord, , drop = FALSE]
    # Abundance covariates
    data$abund.covs <- data$abund.covs[ord, , drop = FALSE]
    # Offset
    offset <- offset[ord]
    for (i in 1:length(data$det.covs)) {
      if (!is.null(dim(data$det.covs[[i]]))) {
        data$det.covs[[i]] <- data$det.covs[[i]][ord, , drop = FALSE]
      } else {
        data$det.covs[[i]] <- data$det.covs[[i]][ord]
      }
    }
  }
  y.mat <- y

  # First subset detection covariates to only use those that are included in the analysis.
  data$det.covs <- data$det.covs[names(data$det.covs) %in% all.vars(det.formula)]
  # Null model support
  if (length(data$det.covs) == 0) {
    data$det.covs <- list(int = rep(1, dim(y)[1]))
  }
  # Make both covariates a data frame. Unlist is necessary for when factors
  # are supplied.
  data$det.covs <- data.frame(lapply(data$det.covs, function(a) unlist(c(a))))
  # Indicator of whether all det.covs are site level or not
  site.level.ind <- ifelse(nrow(data$det.covs) == nrow(y), TRUE, FALSE)
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
  y.na.test <- apply(y, 1, function(a) sum(!is.na(a)))
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
      if (sum(is.na(data$det.covs[, i])) > sum(is.na(y))) {
        stop("error: some elements in det.covs have missing values where there is an observed data value in y. Please either replace the NA values in det.covs with non-missing values (e.g., mean imputation) or set the corresponding values in y to NA where the covariate is missing.")
      }
    }
    # Misalignment between y and det.covs
    y.missing <- which(is.na(y))
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
  # Abundance -------------------------
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
  # Number of replicates at each site
  n.rep <- apply(y, 1, function(a) sum(!is.na(a)))
  # Max number of repeat visits
  K.max <- dim(y.mat)[2]
  # Because I like K better than n.rep
  K <- n.rep

  # Get indices to map N to y ---------------------------------------------
  N.long.indx <- rep(1:J, dim(y)[2])
  N.long.indx <- N.long.indx[!is.na(c(y))]
  # Subtract 1 for indices in C
  N.long.indx <- N.long.indx - 1
  y <- c(y)
  names.long <- which(!is.na(y))
  # Make necessary adjustment for site-level covariates only
  if (nrow(X.p) == J) {
    X.p <- X.p[N.long.indx + 1, , drop = FALSE]
    X.p.re <- X.p.re[N.long.indx + 1, , drop = FALSE]
    X.p.random <- X.p.random[N.long.indx + 1, , drop = FALSE]
  }
  # Remove missing observations when the covariate data are available but
  # there are missing detection-nondetection data.
  if (nrow(X.p) == length(y)) {
    X.p <- X.p[!is.na(y), , drop = FALSE]
  }
  if (nrow(X.p.re) == length(y) & p.det.re > 0) {
    X.p.re <- X.p.re[!is.na(y), , drop = FALSE]
    X.p.random <- X.p.random[!is.na(y), , drop = FALSE]
  }
  y <- y[!is.na(y)]
  # Number of data points for the y vector
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
      message("No prior specified for alpha.normal.\nSetting prior mean to 0 and prior variance to 2.72\n")
    }
    mu.alpha <- rep(0, p.det)
    sigma.alpha <- rep(2.72, p.det)
    Sigma.alpha <- diag(p.det) * 2.72
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
  # phi -----------------------------
  if ("phi.unif" %in% names(priors)) {
    if (!is.vector(priors$phi.unif) | !is.atomic(priors$phi.unif) | length(priors$phi.unif) != 2) {
      stop("error: phi.unif must be a vector of length 2 with elements corresponding to phi's lower and upper bounds")
    }
    phi.a <- priors$phi.unif[1]
    phi.b <- priors$phi.unif[2]
  } else {
    if (verbose) {
      message("No prior specified for phi.unif.\nSetting uniform bounds based on the range of observed spatial coordinates.\n")
    }
    # Get distance matrix which is used if priors are not specified
    coords.D <- iDist(coords)
    phi.a <- 3 / max(coords.D)
    phi.b <- 3 / sort(unique(c(coords.D)))[2]
  }
  if (("sigma.sq.ig" %in% names(priors)) & ("sigma.sq.unif" %in% names(priors))) {
    stop("error: cannot specify both an IG and a uniform prior for sigma.sq")
  }
  if ("sigma.sq.ig" %in% names(priors)) { # inverse-gamma prior.
    sigma.sq.ig <- TRUE
    if (!is.vector(priors$sigma.sq.ig) | !is.atomic(priors$sigma.sq.ig) | length(priors$sigma.sq.ig) != 2) {
      stop("error: sigma.sq.ig must be a vector of length 2 with elements corresponding to sigma.sq's shape and scale parameters")
    }
    sigma.sq.a <- priors$sigma.sq.ig[1]
    sigma.sq.b <- priors$sigma.sq.ig[2]
  } else if ('sigma.sq.unif' %in% names(priors)) {
      sigma.sq.ig <- FALSE
      if (!is.vector(priors$sigma.sq.unif) | !is.atomic(priors$sigma.sq.unif) | length(priors$sigma.sq.unif) != 2) {
        stop("error: sigma.sq.unif must be a vector of length 2 with elements corresponding to sigma.sq's lower and upper bounds")
      }
      sigma.sq.a <- priors$sigma.sq.unif[1]
      sigma.sq.b <- priors$sigma.sq.unif[2]
  } else {
    if (verbose) {
      message("No prior specified for sigma.sq.\nUsing an inverse-Gamma prior with the shape parameter set to 2 and scale parameter to 1.\n")
    }
    sigma.sq.ig <- TRUE
    sigma.sq.a <- 2
    sigma.sq.b <- 1
  }
  # nu -----------------------------
  if (cov.model == 'matern') {
    if (!"nu.unif" %in% names(priors)) {
      stop("error: nu.unif must be specified in priors value list")
    }
    if (!is.vector(priors$nu.unif) | !is.atomic(priors$nu.unif) | length(priors$nu.unif) != 2) {
      stop("error: nu.unif must be a vector of length 2 with elements corresponding to nu's lower and upper bounds")
    }
    nu.a <- priors$nu.unif[1]
    nu.b <- priors$nu.unif[2]
  } else {
    nu.a <- 0
    nu.b <- 0
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
    # Reorder the user supplied inits values
    if (NNGP) {
      N.inits <- N.inits[ord]
    }
    N.test <- apply(y.mat, 1, max, na.rm = TRUE)
    init.test <- sum(N.inits < N.test)
    if (init.test > 0) {
      stop("error: initial values for latent abundance (N) are invalid. Please re-specify inits$N so initial values are greater than or equal to the max number of observed individuals at a given site during a given replicate.")
    }
  } else {
    N.inits <- apply(y.mat, 1, max, na.rm = TRUE)
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
    beta.inits <- rnorm(p.abund, 0, 1)
    if (verbose) {
      message('beta is not specified in initial values.\nSetting initial values to random values from a standard normal distribution\n')
    }
  }
  # alpha -----------------------
  if ("alpha" %in% names(inits)) {
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
    alpha.inits <- rnorm(p.det, 0, 1)
    if (verbose) {
      message("alpha is not specified in initial values.\nSetting initial values to random values from a standard normal distribution\n")
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
      sigma.sq.mu.inits <- runif(p.abund.re, 0.05, 2)
      if (verbose) {
        message("sigma.sq.mu is not specified in initial values.\nSetting initial values to random values between 0.05 and 2\n")
      }
    }
    beta.star.indx <- rep(0:(p.abund.re - 1), n.abund.re.long)
    beta.star.inits <- rnorm(n.abund.re, 0, sqrt(sigma.sq.mu.inits[beta.star.indx + 1]))
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
      sigma.sq.p.inits <- runif(p.det.re, 0.05, 2)
      if (verbose) {
        message("sigma.sq.p is not specified in initial values.\nSetting initial values to random values between 0.05 and 2\n")
      }
    }
    alpha.star.indx <- rep(0:(p.det.re - 1), n.det.re.long)
    alpha.star.inits <- rnorm(n.det.re, 0, sqrt(sigma.sq.p.inits[alpha.star.indx + 1]))
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
  # phi -----------------------------
  if ("phi" %in% names(inits)) {
    phi.inits <- inits[["phi"]]
    if (length(phi.inits) != 1) {
      stop("error: initial values for phi must be of length 1")
    }
  } else {
    phi.inits <- runif(1, phi.a, phi.b)
    if (verbose) {
      message("phi is not specified in initial values.\nSetting initial value to random value from the prior distribution\n")
    }
  }
  # sigma.sq ------------------------
  if ("sigma.sq" %in% names(inits)) {
    sigma.sq.inits <- inits[["sigma.sq"]]
    if (length(sigma.sq.inits) != 1) {
      stop("error: initial values for sigma.sq must be of length 1")
    }
  } else {
    if (sigma.sq.ig) {
    sigma.sq.inits <- runif(1, 0.05, 3)
    } else {
      sigma.sq.inits <- runif(1, sigma.sq.a, sigma.sq.b)
    }
    if (verbose) {
      message("sigma.sq is not specified in initial values.\nSetting initial value to random value between 0.05 and 3 or the user-specified bounds if using a uniform prior.\n")
    }
  }
  # w -----------------------------
  if ("w" %in% names(inits)) {
    w.inits <- inits[["w"]]
    if (!is.vector(w.inits)) {
      stop(paste("error: initial values for w must be a vector of length ",
      	   J, sep = ""))
    }
    if (length(w.inits) != J) {
      stop(paste("error: initial values for w must be a vector of length ",
      	   J, sep = ""))
    }
    # Reorder user supplied initial values.
    w.inits <- w.inits[ord]
  } else {
    w.inits <- rep(0, J)
    if (verbose) {
      message("w is not specified in initial values.\nSetting initial value to 0\n")
    }
  }
  # nu ------------------------
  if ("nu" %in% names(inits)) {
    nu.inits <- inits[["nu"]]
    if (length(nu.inits) != 1) {
      stop("error: initial values for nu must be of length 1")
    }
  } else {
    if (cov.model == 'matern') {
      if (verbose) {
        message("nu is not specified in initial values.\nSetting initial value to random value from the prior distribution\n")
      }
      nu.inits <- runif(1, nu.a, nu.b)
    } else {
      nu.inits <- 0
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
  storage.mode(cov.model.indx) <- "integer"

  # Get tuning values ---------------------------------------------------
  sigma.sq.tuning <- 0
  beta.tuning <- 0
  alpha.tuning <- 0
  w.tuning <- 0
  phi.tuning <- 0
  nu.tuning <- 0
  kappa.tuning <- 0
  beta.star.tuning <- 0
  alpha.star.tuning <- 0
  if (missing(tuning)) {
    beta.tuning <- rep(1, p.abund)
    beta.star.tuning <- rep(1, n.abund.re)
    alpha.tuning <- rep(1, p.det)
    alpha.star.tuning <- rep(1, n.det.re)
    phi.tuning <- 1
    kappa.tuning <- 1
    if (cov.model == 'matern') {
      nu.tuning <- 1
    }
    if (!sigma.sq.ig) {
      sigma.sq.tuning <- 1
    }
    w.tuning <- rep(1, J)
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
    # w ---------------------------
    if(!"w" %in% names(tuning)) {
      stop("error: w must be specified in tuning value list")
    }
    w.tuning <- tuning$w
    if (length(w.tuning) != 1 & length(w.tuning) != J) {
      stop(paste("error: w tuning must be a single value or a vector of length ",
        	 J, sep = ''))
    }
    if (length(w.tuning) == 1) {
      w.tuning <- rep(w.tuning, J)
    }
    # phi ---------------------------
    if(!"phi" %in% names(tuning)) {
      stop("error: phi must be specified in tuning value list")
    }
    phi.tuning <- tuning$phi
    if (length(phi.tuning) != 1) {
      stop("error: phi tuning must be a single value")
    }
    if (cov.model == 'matern') {
      # nu --------------------------
      if(!"nu" %in% names(tuning)) {
        stop("error: nu must be specified in tuning value list")
      }
      nu.tuning <- tuning$nu
      if (length(nu.tuning) != 1) {
        stop("error: nu tuning must be a single value")
      }
    } else {
      nu.tuning <- NULL
    }
    # sigma.sq --------------------------
    if (!sigma.sq.ig) {
      if (!"sigma.sq" %in% names(tuning)) {
        stop("error: sigma.sq must be specified in tuning value list")
      }
      sigma.sq.tuning <- tuning$sigma.sq
      if (length(sigma.sq.tuning) != 1) {
        stop("error: sigma.sq tuning must be a single value")
      }
    }
  }
  tuning.c <- log(c(beta.tuning, alpha.tuning, beta.star.tuning,
		    alpha.star.tuning, sigma.sq.tuning, phi.tuning,
		    nu.tuning, w.tuning, kappa.tuning))

  curr.chain <- 1

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

  if (!NNGP) {
    stop("spNMix is currently only implemented with NNGPs. Please set NNGP = TRUE")
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
    storage.mode(J) <- "integer"

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

    # Get max y values for N update -----------------------------------------
    y.max <- apply(y.mat, 1, max, na.rm = TRUE)

    # Set storage for all variables ---------------------------------------
    storage.mode(y) <- "double"
    storage.mode(N.inits) <- "double"
    storage.mode(X) <- "double"
    storage.mode(X.p) <- "double"
    storage.mode(coords) <- "double"
    storage.mode(y.max) <- "double"
    storage.mode(offset) <- "double"
    consts <- c(J, n.obs, p.abund, p.abund.re, n.abund.re,
                p.det, p.det.re, n.det.re)
    storage.mode(consts) <- "integer"
    storage.mode(K) <- "double"
    storage.mode(beta.inits) <- "double"
    storage.mode(alpha.inits) <- "double"
    storage.mode(kappa.inits) <- "double"
    storage.mode(phi.inits) <- "double"
    storage.mode(sigma.sq.inits) <- "double"
    storage.mode(sigma.sq.ig) <- "integer"
    storage.mode(nu.inits) <- "double"
    storage.mode(w.inits) <- "double"
    storage.mode(N.long.indx) <- "integer"
    storage.mode(mu.beta) <- "double"
    storage.mode(Sigma.beta) <- "double"
    storage.mode(mu.alpha) <- "double"
    storage.mode(Sigma.alpha) <- "double"
    storage.mode(kappa.a) <- "double"
    storage.mode(kappa.b) <- "double"
    storage.mode(phi.a) <- "double"
    storage.mode(phi.b) <- "double"
    storage.mode(nu.a) <- "double"
    storage.mode(nu.b) <- "double"
    storage.mode(sigma.sq.a) <- "double"
    storage.mode(sigma.sq.b) <- "double"
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

    # Fit the model -------------------------------------------------------
    out.tmp <- list()
    out <- list()
    for (i in 1:n.chains) {
      # Change initial values if i > 1
      if ((i > 1) & (!fix.inits)) {
        beta.inits <- rnorm(p.abund, 0, 1)
        alpha.inits <- rnorm(p.det, 0, 1)
        if (p.abund.re > 0) {
          sigma.sq.mu.inits <- runif(p.abund.re, 0.05, 2)
          beta.star.inits <- rnorm(n.abund.re, 0, sqrt(sigma.sq.mu.inits[beta.star.indx + 1]))
        }
        if (p.det.re > 0) {
          sigma.sq.p.inits <- runif(p.det.re, 0.5, 2)
          alpha.star.inits <- rnorm(n.det.re, 0, sqrt(sigma.sq.p.inits[alpha.star.indx + 1]))
        }
        kappa.inits <- runif(1, kappa.a, kappa.b)
        if (!sigma.sq.ig) {
          sigma.sq.inits <- runif(1, sigma.sq.a, sigma.sq.b)
        } else {
          sigma.sq.inits <- runif(1, 0.05, 3)
        }
        phi.inits <- runif(1, phi.a, phi.b)
        if (cov.model == 'matern') {
          nu.inits <- runif(1, nu.a, nu.b)
        }
      }
      storage.mode(chain.info) <- "integer"
      # Run the model in C
      out.tmp[[i]] <- .Call("spNMixNNGP", y, X, X.p, coords, X.re, X.p.re, X.random, X.p.random,
        		    y.max, consts,
        		    n.abund.re.long, n.det.re.long,
          	            n.neighbors, nn.indx, nn.indx.lu, u.indx, u.indx.lu, ui.indx,
        		    beta.inits, alpha.inits, kappa.inits,
          		    sigma.sq.mu.inits, sigma.sq.p.inits, beta.star.inits,
          		    alpha.star.inits, N.inits, w.inits, phi.inits, sigma.sq.inits, nu.inits,
        		    N.long.indx, beta.star.indx, beta.level.indx,
        		    alpha.star.indx, alpha.level.indx, mu.beta, Sigma.beta, mu.alpha,
        		    Sigma.alpha, phi.a, phi.b, sigma.sq.a, sigma.sq.b, nu.a, nu.b,
        		    sigma.sq.mu.a, sigma.sq.mu.b,
          		    sigma.sq.p.a, sigma.sq.p.b, kappa.a, kappa.b,
        		    tuning.c, cov.model.indx, n.batch, batch.length, accept.rate,
        		    n.omp.threads, verbose, n.report, samples.info, chain.info,
        		    sigma.sq.ig, family.c, offset)
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
      out$rhat$theta <- gelman.diag(mcmc.list(lapply(out.tmp, function(a)
      					        mcmc(t(a$theta.samples)))),
      			      autoburnin = FALSE, multivariate = FALSE)$psrf[, 2]
    } else {
      out$rhat$beta <- rep(NA, p.abund)
      out$rhat$kappa <- NA
      out$rhat$theta <- rep(NA, ifelse(cov.model == 'matern', 3, 2))
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
    out$theta.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$theta.samples))))
    if (cov.model != 'matern') {
      colnames(out$theta.samples) <- c('sigma.sq', 'phi')
    } else {
      colnames(out$theta.samples) <- c('sigma.sq', 'phi', 'nu')
    }
    out$N.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$N.samples))))
    out$N.samples <- mcmc(out$N.samples[, order(ord), drop = FALSE])
    out$mu.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$mu.samples))))
    out$mu.samples <- mcmc(out$mu.samples[, order(ord), drop = FALSE])
    out$w.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$w.samples))))
    out$w.samples <- mcmc(out$w.samples[, order(ord), drop = FALSE])
    # Get detection covariate stuff in right order.
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
    out$ESS$theta <- effectiveSize(out$theta.samples)
    out$ESS$alpha <- effectiveSize(out$alpha.samples)
    if (p.det.re > 0) {
      out$ESS$sigma.sq.p <- effectiveSize(out$sigma.sq.p.samples)
    }
    if (p.abund.re > 0) {
      out$ESS$sigma.sq.mu <- effectiveSize(out$sigma.sq.mu.samples)
    }
    out$X <- X[order(ord), , drop = FALSE]
    out$X.re <- X.re[order(ord), , drop = FALSE]
    out$y <- y.mat[order(ord), , drop = FALSE]
    out$offset <- offset[order(ord)]
    out$n.samples <- n.samples
    out$call <- cl
    out$n.neighbors <- n.neighbors
    out$coords <- coords[order(ord), ]
    out$cov.model.indx <- cov.model.indx
    out$type <- "NNGP"
    out$n.post <- n.post.samples
    out$n.thin <- n.thin
    out$n.burn <- n.burn
    out$n.chains <- n.chains
    out$re.cols <- re.cols
    out$re.det.cols <- re.det.cols
    out$dist <- family
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
  class(out) <- "spNMix"
  out$run.time <- proc.time() - ptm
  out
}
