abundGaussian <- function(formula, data, inits, priors, tuning, n.batch,
                          batch.length, accept.rate = 0.43, family = 'Gaussian',
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
      data$covs <- list(int = rep(1, length(y)))
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
  data$covs <- data$covs[z.indx, , drop = FALSE]
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

  # Check save.fitted ---------------------------------------------------
  if (!(save.fitted %in% c(TRUE, FALSE))) {
    stop("save.fitted must be either TRUE or FALSE")
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
  J <- nrow(X)
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
  # Check if n.burn, n.thin, and n.samples result in an integer and error if otherwise.
  if (((n.samples - n.burn) / n.thin) %% 1 != 0) {
    stop("the number of posterior samples to save ((n.samples - n.burn) / n.thin) is not a whole number. Please respecify the MCMC criteria such that the number of posterior samples saved is a whole number.")
  }
  n.post.samples <- length(seq(from = n.burn + 1,
                               to = n.samples,
                               by = as.integer(n.thin)))

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
  # tau.sq.t ----------------------
  if ("tau.sq.ig" %in% names(priors)) {
    if (!is.vector(priors$tau.sq.ig) | !is.atomic(priors$tau.sq.ig) | length(priors$tau.sq.ig) != 2) {
      stop("error: tau.sq.ig must be a vector of length 2 with elements corresponding to tau.sq's shape and scale parameters")
    }
    tau.sq.a <- priors$tau.sq.ig[1]
    tau.sq.b <- priors$tau.sq.ig[2]
  } else {
    if (verbose) {
      message("No prior specified for tau.sq.\nUsing an inverse-Gamma prior with the shape parameter set to 0.01 and scale parameter to 0.01.\n")
    }
    tau.sq.a <- .01
    tau.sq.b <- .01
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
    beta.star.inits <- rnorm(n.re, 0, sqrt(sigma.sq.mu.inits[beta.star.indx + 1]))
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
  consts <- c(J.est, p, p.re, n.re, J.zero, save.fitted)
  storage.mode(consts) <- "integer"
  storage.mode(beta.inits) <- "double"
  storage.mode(tau.sq.inits) <- "double"
  storage.mode(mu.beta) <- "double"
  storage.mode(Sigma.beta) <- "double"
  storage.mode(tau.sq.a) <- "double"
  storage.mode(tau.sq.b) <- "double"
  storage.mode(n.batch) <- "integer"
  storage.mode(batch.length) <- "integer"
  storage.mode(accept.rate) <- "double"
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
      if (p.re > 0) {
        sigma.sq.mu.inits <- runif(p.re, 0.5, 10)
        beta.star.inits <- rnorm(n.re, 0, sqrt(sigma.sq.mu.inits[beta.star.indx + 1]))
      }
      tau.sq.inits <- runif(1, 0.1, 10)
    }
    storage.mode(chain.info) <- "integer"
    # Run the model in C
    out.tmp[[i]] <- .Call("abundGaussian", y, X, X.re, X.random,
    		          consts, n.re.long,
                          beta.inits, tau.sq.inits, sigma.sq.mu.inits, beta.star.inits,
                          beta.star.indx, beta.level.indx, mu.beta,
                          Sigma.beta, tau.sq.a, tau.sq.b,
    		          sigma.sq.mu.a, sigma.sq.mu.b,
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
    			     autoburnin = FALSE, multivariate = FALSE)$psrf[, 2])
    out$rhat$tau.sq <- as.vector(gelman.diag(mcmc.list(lapply(out.tmp, function(a)
    					      mcmc(t(a$tau.sq.samples)))),
    			     autoburnin = FALSE, multivariate = FALSE)$psrf[, 2])
    if (p.re > 0) {
      out$rhat$sigma.sq.mu <- as.vector(gelman.diag(mcmc.list(lapply(out.tmp, function(a)
      					      mcmc(t(a$sigma.sq.mu.samples)))),
      			     autoburnin = FALSE, multivariate = FALSE)$psrf[, 2])
    }
  } else {
    out$rhat$beta <- rep(NA, p)
    out$rhat$tau.sq <- NA
    if (p.re > 0) {
      out$rhat$sigma.sq.mu <- rep(NA, p.re)
    }
  }
  out$beta.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$beta.samples))))
  colnames(out$beta.samples) <- x.names
  out$tau.sq.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$tau.sq.samples))))
  colnames(out$tau.sq.samples) <- 'tau.sq'
  # Get everything back in the original order
  if (save.fitted) {
    if (!two.stage) {
      out$y.rep.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$y.rep.samples))))
      out$like.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$like.samples))))
    } else {
      y.rep.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$y.rep.samples))))
      y.rep.zero.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$y.rep.zero.samples))))
      out$y.rep.samples <- matrix(NA, n.post.samples * n.chains, J.est + J.zero)
      out$y.rep.samples[, z.indx] <- y.rep.samples
      out$y.rep.samples[, -z.indx] <- y.rep.zero.samples
      out$y.rep.samples <- mcmc(out$y.rep.samples)
      out$like.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$like.samples))))
    }
    out$mu.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$mu.samples))))
    out$mu.samples <- mcmc(out$mu.samples)
  }
  out$X <- X
  out$X.re <- X.re
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
  if (p.re > 0) {
    out$ESS$sigma.sq.mu <- effectiveSize(out$sigma.sq.mu.samples)
  }
  out$call <- cl
  out$n.samples <- batch.length * n.batch
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
  class(out) <- "abund"
  out$run.time <- proc.time() - ptm
  return(out)
}

