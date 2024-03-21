abund <- function(formula, data, inits, priors, tuning,
		  n.batch, batch.length, accept.rate = 0.43, family = 'Poisson',
		  n.omp.threads = 1, verbose = TRUE,
		  n.report = 100, n.burn = round(.10 * n.batch * batch.length), n.thin = 1,
		  n.chains = 1, save.fitted = TRUE, ...) {

  ptm <- proc.time()

  if (!(family) %in% c('Poisson', 'NB', 'Gaussian', 'zi-Gaussian')) {
    stop("family must be either 'Poisson', 'NB', 'Gaussian', or 'zi-Gaussian'")
  }
  if (family %in% c('Gaussian', 'zi-Gaussian')) {
    abundGaussian(formula, data, inits, priors, tuning,
		  n.batch, batch.length, accept.rate, family,
		  n.omp.threads, verbose, n.report, n.burn, n.thin, n.chains, save.fitted)
  } else {

    # Make it look nice
    if (verbose) {
      cat("----------------------------------------\n");
      cat("\tPreparing to run the model\n");
      cat("----------------------------------------\n");
    }

    # Functions ---------------------------------------------------------------
    logit <- function(theta, a = 0, b = 1) {log((theta-a)/(b-theta))}
    logit.inv <- function(z, a = 0, b = 1) {b-(b-a)/(1+exp(z))}

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
      stop("error: data y must be specified in data")
    }
    y <- as.matrix(data$y)
    y.mat <- y
    # Offset
    if ('offset' %in% names(data)) {
      offset <- data$offset
      if (length(offset) == 1) {
        offset <- matrix(offset, nrow(y), ncol(y))
      } else if (length(dim(offset)) == 1) { # Value for each site
        if (length(offset) != nrow(y)) {
          stop(paste0("offset must be a single value, vector of length ", nrow(y), " or a matrix with ",
                     nrow(y), " rows and ", ncol(y), " columns."))
        }
        offset <- matrix(offset, nrow(y), ncol(y))
      } else if (length(dim(offset)) == 2) { # Value for each site/obs
        if (nrow(offset) != nrow(y) | ncol(offset) != ncol(y)) {
          stop(paste0("offset must be a single value, vector of length ", nrow(y), " or a matrix with ",
                      nrow(y), " rows and ", ncol(y), " columns."))

        }
      }
    } else {
      offset <- matrix(1, nrow(y), ncol(y))
    }
    offset.mat <- offset
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
      data$covs <- as.data.frame(lapply(data$covs, rep, dim(y)[2]))
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
    # y -------------------------------
    y.na.test <- apply(y, 1, function(a) sum(!is.na(a)))
    if (sum(y.na.test == 0) > 0) {
      stop("error: some sites in y have all missing detection histories. Remove these sites from all objects in the 'data' argument, then use 'predict' to obtain predictions at these locations if desired.")
    }
    # covs ------------------------------
    for (i in 1:ncol(data$covs)) {
      if (sum(is.na(data$covs[, i])) > sum(is.na(y))) {
        stop("error: some elements in covs have missing values where there is an observed data value in y. Please either replace the NA values in covs with non-missing values (e.g., mean imputation) or set the corresponding values in y to NA where the covariate is missing.")
      }
    }
    # Misalignment between y and covs
    y.missing <- which(is.na(data$y))
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

    # Remove missing values from covs in order to ensure formula parsing
    # works when random slopes are provided.
    tmp <- apply(data$covs, 1, function (a) sum(is.na(a)))
    data$covs <- as.data.frame(data$covs[tmp == 0, , drop = FALSE])

    # Check save.fitted ---------------------------------------------------
    if (!(save.fitted %in% c(TRUE, FALSE))) {
      stop("save.fitted must be either TRUE or FALSE")
    }

    # Formula -------------------------------------------------------------
    # Abundance -------------------------
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
    J <- nrow(y)
    # Number of abundance parameters
    p.abund <- ncol(X)
    # Number of abundance random effect parameters
    p.abund.re <- ncol(X.re)
    # Number of latent abundance random effect values
    n.abund.re <- length(unlist(apply(X.re, 2, unique)))
    n.abund.re.long <- apply(X.re, 2, function(a) length(unique(a)))
    # Number of replicates at each site
    n.rep <- apply(y, 1, function(a) sum(!is.na(a)))
    # Max number of repeat visits
    K.max <- ncol(y)
    # Because I like K better than n.rep
    K <- n.rep

    # Get indices to map N to y ---------------------------------------------
    site.indx <- rep(1:J, dim(y)[2])
    site.indx <- site.indx[!is.na(c(y))]
    # Subtract 1 for indices in C
    site.indx <- site.indx - 1
    y <- c(y)
    offset <- c(offset)
    names.long <- which(!is.na(y))
    # Remove missing observations when the covariate data are available but
    # there are missing abundance data.
    if (nrow(X) == length(y)) {
      X <- X[!is.na(y), , drop = FALSE]
    }
    if (nrow(X.re) == length(y) & p.abund.re > 0) {
      X.re <- X.re[!is.na(y), , drop = FALSE]
    }
    if (nrow(X.random) == length(y) & p.abund.re > 0) {
      X.random <- X.random[!is.na(y), , drop = FALSE]
    }
    y <- y[!is.na(y.mat)]
    offset <- offset[!is.na(y.mat)]
    # Number of data points for the y vector
    n.obs <- nrow(X)

    # Get random effect matrices all set ----------------------------------
    X.re <- X.re - 1
    if (p.abund.re > 1) {
      # Subtract 1 for C
      for (j in 2:p.abund.re) {
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
    tuning.c <- log(c(beta.tuning, beta.star.tuning,
          	    kappa.tuning))
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

    # Set storage for all variables ---------------------------------------
    storage.mode(y) <- "double"
    storage.mode(X) <- "double"
    storage.mode(offset) <- "double"
    consts <- c(J, n.obs, p.abund, p.abund.re, n.abund.re, save.fitted)
    storage.mode(consts) <- "integer"
    storage.mode(beta.inits) <- "double"
    storage.mode(kappa.inits) <- "double"
    storage.mode(site.indx) <- "integer"
    storage.mode(mu.beta) <- "double"
    storage.mode(Sigma.beta) <- "double"
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
        if (family == 'NB') {
          kappa.inits <- runif(1, kappa.a, kappa.b)
        }
        if (p.abund.re > 0) {
          sigma.sq.mu.inits <- runif(p.abund.re, 0.05, 1)
          beta.star.inits <- rnorm(n.abund.re, sqrt(sigma.sq.mu.inits[beta.star.indx + 1]))
        }
      }
      storage.mode(chain.info) <- "integer"
      # Run the model in C
      out.tmp[[i]] <- .Call("abund", y, X, X.re, X.random, consts, n.abund.re.long,
          		   beta.inits, kappa.inits, sigma.sq.mu.inits, beta.star.inits,
        		           site.indx, beta.star.indx,
        		           beta.level.indx, mu.beta, Sigma.beta, kappa.a, kappa.b,
        		           sigma.sq.mu.a, sigma.sq.mu.b,
        		           tuning.c, n.batch, batch.length, accept.rate,
        		           n.omp.threads, verbose, n.report, samples.info, chain.info,
                             family.c, offset)
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
      if (p.abund.re > 0) {
        out$rhat$sigma.sq.mu <- rep(NA, p.abund.re)
      }
    }
    # Put everything into MCMC objects
    out$beta.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$beta.samples))))
    colnames(out$beta.samples) <- x.names
    if (family == 'NB') {
      out$kappa.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$kappa.samples))))
      colnames(out$kappa.samples) <- c("kappa")
    }
    y.non.miss.indx <- which(!is.na(y.mat), arr.ind = TRUE)
    if (save.fitted) {
      out$y.rep.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$y.rep.samples))))
      tmp <- array(NA, dim = c(n.post.samples * n.chains, J, K.max))
      for (j in 1:n.obs) {
        curr.indx <- y.non.miss.indx[j, ]
        tmp[, curr.indx[1], curr.indx[2]] <- out$y.rep.samples[, j]
      }
      out$y.rep.samples <- tmp
      out$mu.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$mu.samples))))
      tmp <- array(NA, dim = c(n.post.samples * n.chains, J, K.max))
      for (j in 1:n.obs) {
        curr.indx <- y.non.miss.indx[j, ]
        tmp[, curr.indx[1], curr.indx[2]] <- out$mu.samples[, j]
      }
      out$mu.samples <- tmp
      out$like.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$like.samples))))
      tmp <- array(NA, dim = c(n.post.samples * n.chains, J, K.max))
      for (j in 1:n.obs) {
        curr.indx <- y.non.miss.indx[j, ]
        tmp[, curr.indx[1], curr.indx[2]] <- out$like.samples[, j]
      }
      out$like.samples <- tmp
    }
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
    # Calculate effective sample sizes
    out$ESS <- list()
    out$ESS$beta <- effectiveSize(out$beta.samples)
    if (family == 'NB') {
      out$ESS$kappa <- effectiveSize(out$kappa.samples)
    }
    if (p.abund.re > 0) {
      out$ESS$sigma.sq.mu <- effectiveSize(out$sigma.sq.mu.samples)
    }
    out$X <- array(NA, dim = c(J, ncol(y.mat), p.abund))
    out$X.re <- array(NA, dim = c(J, ncol(y.mat), p.abund.re))
    for (j in 1:n.obs) {
      curr.indx <- y.non.miss.indx[j, ]
      out$X[curr.indx[1], curr.indx[2], ] <- X[j, ]
      if (p.abund.re > 0) {
        out$X.re[curr.indx[1], curr.indx[2], ] <- X.re[j, ]
      }
    }
    dimnames(out$X)[[3]] <- x.names
    dimnames(out$X.re)[[3]] <- colnames(X.re)
    out$y <- y.mat
    out$offset <- offset.mat
    out$n.samples <- n.samples
    out$call <- cl
    out$n.post <- n.post.samples
    out$n.thin <- n.thin
    out$n.burn <- n.burn
    out$n.chains <- n.chains
    out$dist <- family
    out$re.cols <- re.cols
    if (p.abund.re > 0) {
      out$muRE <- TRUE
    } else {
      out$muRE <- FALSE
    }
    class(out) <- "abund"
    out$run.time <- proc.time() - ptm
    out
  } # Gaussian vs. NB/Poisson
}
