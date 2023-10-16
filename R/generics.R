# ppcAbund ------------------------------------------------------------------
summary.ppcAbund <- function(object, level = 'both',
			     digits = max(3L, getOption("digits") - 3L), ...) {

  cat("\nCall:", deparse(object$call, width.cutoff = floor(getOption("width") * 0.75)),
      "", sep = "\n")

  n.post <- object$n.post
  n.samples <- object$n.samples
  n.burn <- object$n.burn
  n.thin <- object$n.thin
  n.chains <- object$n.chains

  cat(paste("Samples per Chain: ", n.samples,"\n", sep=""))
  cat(paste("Burn-in: ", n.burn,"\n", sep=""))
  cat(paste("Thinning Rate: ",n.thin,"\n", sep=""))
  cat(paste("Number of Chains: ", n.chains, "\n", sep = ""))
  cat(paste("Total Posterior Samples: ",n.post * n.chains,"\n\n", sep=""))

  if (object$class %in% c('NMix', 'spNMix', 'abund', 'spAbund', 'DS', 'spDS')) {
    cat("Bayesian p-value: ", round(mean(object$fit.y.rep > object$fit.y), digits), "\n")
    cat("Fit statistic: ", object$fit.stat, "\n")
  }

  if (object$class %in% c('msAbund', 'spMsAbund', 'lfMsAbund', 'sfMsAbund', 
			  'msNMix', 'spMsNMix', 'lfMsNMix', 'sfMsNMix', 
			  'msDS', 'lfMsDS', 'sfMsDS')) {

    if (tolower(level) == 'community') {
      cat("----------------------------------------\n");
      cat("\tCommunity Level\n");
      cat("----------------------------------------\n");
      cat("Bayesian p-value: ", round(mean(object$fit.y.rep > object$fit.y), digits), "\n")
      cat("Fit statistic: ", object$fit.stat, "\n")
    }

    if (tolower(level) == 'species') {
      cat("----------------------------------------\n");
      cat("\tSpecies Level\n");
      cat("----------------------------------------\n");
      N <- ncol(object$fit.y)
      for (i in 1:N) {
        cat(paste(object$sp.names[i], " Bayesian p-value: ",
        	  round(mean(object$fit.y.rep[, i] > object$fit.y[, i]), digits), "\n", sep = ''))
      }
      cat("Fit statistic: ", object$fit.stat, "\n")
    }

    if (tolower(level) == 'both') {
      cat("----------------------------------------\n");
      cat("\tCommunity Level\n");
      cat("----------------------------------------\n");
      cat("Bayesian p-value: ", round(mean(object$fit.y.rep > object$fit.y), digits), "\n")
      cat("\n")
      cat("----------------------------------------\n");
      cat("\tSpecies Level\n");
      cat("----------------------------------------\n");
      N <- ncol(object$fit.y)
      for (i in 1:N) {
        cat(paste(object$sp.names[i], " Bayesian p-value: ",
        	  round(mean(object$fit.y.rep[, i] > object$fit.y[, i]), digits), "\n", sep = ''))
      }
      cat("Fit statistic: ", object$fit.stat, "\n")
    }
  }

}

# abund --------------------------------------------------------------------
print.abund <- function(x, ...) {
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.75)),
      "", sep = "\n")
}

fitted.abund <- function(object, ...) {
  return(object$y.rep.samples)
}

summary.abund <- function(object,
			  quantiles = c(0.025, 0.5, 0.975),
			  digits = max(3L, getOption("digits") - 3L), ...) {
  print(object)

  n.post <- object$n.post
  n.samples <- object$n.samples
  n.burn <- object$n.burn
  n.thin <- object$n.thin
  n.chains <- object$n.chains
  run.time <- object$run.time[3] / 60 # minutes

  cat(paste("Samples per Chain: ", n.samples,"\n", sep=""))
  cat(paste("Burn-in: ", n.burn,"\n", sep=""))
  cat(paste("Thinning Rate: ",n.thin,"\n", sep=""))
  cat(paste("Number of Chains: ", n.chains, "\n", sep = ""))
  cat(paste("Total Posterior Samples: ",n.post * n.chains,"\n", sep=""))
  cat(paste("Run Time (min): ", round(run.time, digits), "\n\n", sep = ""))

  # Abundance
  if (object$dist %in% c('Poisson', 'NB')) {
    cat("Abundance (log scale): \n")
  }
  if (object$dist %in% c('Gaussian', 'zi-Gaussian')) {
    cat("Abundance: \n")
  }
  tmp.1 <- t(apply(object$beta.samples, 2,
		   function(x) c(mean(x), sd(x))))
  colnames(tmp.1) <- c("Mean", "SD")
  tmp <- t(apply(object$beta.samples, 2,
		 function(x) quantile(x, prob = quantiles)))
  diags <- matrix(c(object$rhat$beta, round(object$ESS$beta, 0)), ncol = 2)
  colnames(diags) <- c('Rhat', 'ESS')

  print(noquote(round(cbind(tmp.1, tmp, diags), digits)))
  if (object$muRE) {
    cat("\n")
    if (object$dist %in% c('Poisson', 'NB')) {
      cat("Abundance Random Effect Variances (log scale): \n")
    }
    if (object$dist %in% c('Gaussian', 'zi-Gaussian')) {
      cat("Abundance Random Effect Variances: \n")
    }
    tmp.1 <- t(apply(object$sigma.sq.mu.samples, 2,
          	   function(x) c(mean(x), sd(x))))
    colnames(tmp.1) <- c("Mean", "SD")
    tmp <- t(apply(object$sigma.sq.mu.samples, 2,
          	 function(x) quantile(x, prob = quantiles)))
    diags <- matrix(c(object$rhat$sigma.sq.mu, round(object$ESS$sigma.sq.mu, 0)), ncol = 2)
    colnames(diags) <- c('Rhat', 'ESS')

    print(noquote(round(cbind(tmp.1, tmp, diags), digits)))
  }
  if (object$dist == "NB") {
    cat("\n")
    cat("NB overdispersion: \n")
    tmp.1 <- t(apply(object$kappa.samples, 2,
          	   function(x) c(mean(x), sd(x))))
    colnames(tmp.1) <- c("Mean", "SD")
    tmp <- t(apply(object$kappa.samples, 2,
          	 function(x) quantile(x, prob = quantiles)))
    diags <- matrix(c(object$rhat$kappa, round(object$ESS$kappa, 0)), ncol = 2)
    colnames(diags) <- c('Rhat', 'ESS')

    print(noquote(round(cbind(tmp.1, tmp, diags), digits)))
  }
  if (object$dist %in% c('Gaussian', 'zi-Gaussian')) {
    cat("\n")
    cat("Residual variance: \n")
    tmp.1 <- t(apply(object$tau.sq.samples, 2,
          	   function(x) c(mean(x), sd(x))))
    colnames(tmp.1) <- c("Mean", "SD")
    tmp <- t(apply(object$tau.sq.samples, 2,
          	 function(x) quantile(x, prob = quantiles)))
    diags <- matrix(c(object$rhat$tau.sq, round(object$ESS$tau.sq, 0)), ncol = 2)
    colnames(diags) <- c('Rhat', 'ESS')

    print(noquote(round(cbind(tmp.1, tmp, diags), digits)))
  }
}

predict.abund <- function(object, X.0, ignore.RE = FALSE, z.0.samples, ...) {
  # Check for unused arguments ------------------------------------------
  formal.args <- names(formals(sys.function(sys.parent())))
  elip.args <- names(list(...))
  for(i in elip.args){
      if(! i %in% formal.args)
          warning("'",i, "' is not an argument")
  }
  # Call ----------------------------------------------------------------
  cl <- match.call()

  # Functions ---------------------------------------------------------------
  logit <- function(theta, a = 0, b = 1) {log((theta-a)/(b-theta))}
  logit.inv <- function(z, a = 0, b = 1) {b-(b-a)/(1+exp(z))}

  # Some initial checks ---------------------------------------------------
  if (missing(object)) {
    stop("error: predict expects object\n")
  }
  if (!(class(object) %in% c('abund', 'spAbund'))) {
    stop("error: requires an output object of class abund or spAbund\n")
  }

  # Check X.0 -------------------------------------------------------------
  if (missing(X.0)) {
    stop("error: X.0 must be specified\n")
  }
  if (!(length(dim(X.0)) %in% c(2, 3))) {
    stop("error: X.0 must be a matrix with two columns corresponding to site and covariate or a three-dimensional array with dimensions corresponding to site, replicate, and covariate")
  }
  if (length(dim(X.0)) == 2) {
    tmp <- colnames(X.0) 
    X.0 <- array(X.0, dim = c(nrow(X.0), 1, ncol(X.0)))
    dimnames(X.0)[[3]] <- tmp
  }

  # Predictions -----------------------------------------------------------
  p.abund <- dim(object$X)[3]
  # p.design <- p.abund
  # if (object$muRE & !ignore.RE) {
  #   p.design <- p.abund + ncol(object$sigma.sq.mu.samples)
  # }
  # if (dim(X.0)[3] != p.design) {
  #   stop(paste("error: X.0 must have ", p.design, " covariates\n", sep = ''))
  # }
  J.0 <- nrow(X.0)
  K.max.0 <- ncol(X.0)
  n.post <- object$n.post * object$n.chains
  
  # Check z.0.samples -----------------------------------------------------
  # Check stage 1 samples
  if (object$dist == 'zi-Gaussian') {
    if (missing(z.0.samples)) {
      stop("z.0.samples must be supplied for a zi-Gaussian model")
    }
    if (!is.matrix(z.0.samples)) {
      stop(paste("z.0.samples must be a matrix with ", n.post, " rows and ", 
		 J.0, " columns.", sep = ''))
    }
    if (nrow(z.0.samples) != n.post | ncol(z.0.samples) != J.0) {
      stop(paste("z.0.samples must be a matrix with ", n.post, " rows and ", 
		 J.0, " columns.", sep = ''))
    }
  } else {
    if (!missing(z.0.samples)) {
      message("z.0.samples is ignored for the current model family\n")
    }
    z.0.samples <- NA
  }

  # Composition sampling --------------------------------------------------
  re.cols <- object$re.cols
  beta.samples <- as.matrix(object$beta.samples)
  if (object$dist == 'NB') {
    kappa.samples <- as.matrix(object$kappa.samples)
  }
  out <- list()
  if (object$muRE) {
    p.abund.re <- length(object$re.level.names)
  } else {
    p.abund.re <- 0
  }
  if (object$dist %in% c('Gaussian', 'zi-Gaussian')) {
    tau.sq.samples <- as.matrix(object$tau.sq.samples)
  }

  # Get X.0 in long format. 
  tmp.names <- dimnames(X.0)[[3]]
  X.0 <- matrix(X.0, nrow = nrow(X.0) * ncol(X.0), ncol = dim(X.0)[3])
  colnames(X.0) <- tmp.names
  missing.indx <- which(apply(X.0, 1, function(a) sum(is.na(a)) > 0))
  non.missing.indx <- which(apply(X.0, 1, function(a) sum(is.na(a)) == 0))
  X.0 <- X.0[non.missing.indx, , drop = FALSE]

  if (object$muRE & !ignore.RE) {
    beta.star.samples <- object$beta.star.samples
    re.level.names <- object$re.level.names
    # Get columns in design matrix with random effects
    if (object$dist %in% c('Gaussian', 'zi-Gaussian')) {
      x.re.names <- dimnames(object$X.re)[[2]]
    } else {
      x.re.names <- dimnames(object$X.re)[[3]]
    }
    x.0.names <- colnames(X.0)
    # Get the number of times each factor is used. 
    re.long.indx <- sapply(re.cols, length)
    tmp <- sapply(x.re.names, function(a) which(colnames(X.0) %in% a))
    indx <- list()
    for (i in 1:length(tmp)) {
      indx[[i]] <- rep(tmp[i], re.long.indx[i])
    }
    indx <- unlist(indx)
    if (length(indx) == 0) {
      stop("error: column names in X.0 must match variable names in data$abund.covs")
    }
    n.re <- length(indx)
    n.unique.re <- length(unique(indx))
    # Check RE columns
    for (i in 1:n.re) {
      if (is.character(re.cols[[i]])) {
        # Check if all column names in svc are in occ.covs
        if (!all(re.cols[[i]] %in% x.0.names)) {
            missing.cols <- re.cols[[i]][!(re.cols[[i]] %in% x.0.names)]
            stop(paste("error: variable name ", paste(missing.cols, collapse=" and "), " not in occurrence covariates", sep=""))
        }
        # Convert desired column names into the numeric column index
        re.cols[[i]] <- which(x.0.names %in% re.cols[[i]])
        
      } else if (is.numeric(re.cols[[i]])) {
        # Check if all column indices are in 1:p.abund
        if (!all(re.cols %in% 1:p.abund)) {
            missing.cols <- re.cols[[i]][!(re.cols[[i]] %in% (1:p.abund))]
            stop(paste("error: column index ", paste(missing.cols, collapse=" "), " not in design matrix columns", sep=""))
        }
      }
    }
    re.cols <- unlist(re.cols)
    X.re <- as.matrix(X.0[, indx, drop = FALSE])
    X.fix <- as.matrix(X.0[, -indx, drop = FALSE])
    X.random <- as.matrix(X.0[, re.cols, drop = FALSE])
    n.abund.re <- length(unlist(re.level.names))
    X.re.ind <- matrix(NA, nrow(X.re), p.abund.re)

    for (i in 1:p.abund.re) {
      for (j in 1:nrow(X.re)) {
        tmp <- which(re.level.names[[i]] == X.re[j, i])
        if (length(tmp) > 0) {
          X.re.ind[j, i] <- tmp 
        }
      }
    }
    if (p.abund.re > 1) {
      for (j in 2:p.abund.re) {
        X.re.ind[, j] <- X.re.ind[, j] + max(X.re.ind[, j - 1]) 
      }
    }
    # Create the random effects corresponding to each 
    # new location
    # ORDER: ordered by site, then species within site.
    beta.star.sites.0.samples <- matrix(0, n.post,  nrow(X.re))
    for (t in 1:p.abund.re) {
      for (j in 1:nrow(X.re)) {
        if (!is.na(X.re.ind[j, t])) {
          beta.star.sites.0.samples[, j] <- 
            beta.star.samples[, X.re.ind[j, t]] * X.random[j, t] + 
            beta.star.sites.0.samples[, j]
        } else {
          beta.star.sites.0.samples[, j] <- 
            rnorm(n.post, 0, sqrt(object$sigma.sq.mu.samples[, t])) * X.random[j, t] + 
            beta.star.sites.0.samples[, j]
        }
      } # j
    } # t
  } else {
    X.fix <- X.0
    beta.star.sites.0.samples <- matrix(0, n.post, nrow(X.0))
    p.abund.re <- 0
  }

  tmp <- matrix(NA, n.post, length(c(missing.indx, non.missing.indx)))
  if (object$dist %in% c('Poisson', 'NB')) {
  mu.long <- exp(t(X.fix %*% t(beta.samples) + 
	       t(beta.star.sites.0.samples)))
  } 
  if (object$dist %in% c('Gaussian', 'zi-Gaussian')) {
  mu.long <- t(X.fix %*% t(beta.samples) + 
	       t(beta.star.sites.0.samples))
  }
  tmp[, non.missing.indx] <- mu.long
  out$mu.0.samples <- array(tmp, dim = c(n.post, J.0, K.max.0))
  K <- apply(out$mu.0.samples[1, , , drop = FALSE], 2, function(a) sum(!is.na(a)))
  out$y.0.samples <- array(NA, dim(out$mu.0.samples))
  J <- nrow(object$y)
  for (j in 1:J.0) {
    for (k in 1:K.max.0) {
      if (sum(is.na(out$mu.0.samples[, j, k])) == 0) {
        if (object$dist == 'NB') {
          out$y.0.samples[, j, k] <- rnbinom(n.post, kappa.samples, 
          				   mu = out$mu.0.samples[, j, k])
        } else if (object$dist == 'Poisson') {
          out$y.0.samples[, j, k] <- rpois(n.post, out$mu.0.samples[, j, k])
        } else if (object$dist == 'Gaussian') {
          out$y.0.samples[, j, k] <- rnorm(n.post, out$mu.0.samples[, j, k], 
                                           sqrt(tau.sq.samples))
        } else if (object$dist == 'zi-Gaussian') {
          out$y.0.samples[, j, k] <- ifelse(z.0.samples[, j] == 1, 
          				  rnorm(n.post, out$mu.0.samples[, j, k], 
          					sqrt(tau.sq.samples)), 
          				  rnorm(n.post, 0, sqrt(0.0001)))
          out$mu.0.samples[, j, k] <- ifelse(z.0.samples[, j] == 1, 
          				   out$mu.0.samples[, j, k], 0)
        }
      }
    }
  }
  # If Gaussian, collapse to a matrix
  if (object$dist %in% c('Gaussian', 'zi-Gaussian')) {
    out$y.0.samples <- out$y.0.samples[, , 1]  
    out$mu.0.samples <- out$mu.0.samples[, , 1]  
  }

  out$call <- cl

  class(out) <- "predict.abund"
  out
}

# spAbund --------------------------------------------------------------------
print.spAbund <- function(x, ...) {
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.75)),
      "", sep = "\n")
}

fitted.spAbund <- function(object, ...) {
  return(object$y.rep.samples)
}

summary.spAbund <- function(object,
			    quantiles = c(0.025, 0.5, 0.975),
			    digits = max(3L, getOption("digits") - 3L), ...) {
  print(object)

  n.post <- object$n.post
  n.samples <- object$n.samples
  n.burn <- object$n.burn
  n.thin <- object$n.thin
  n.chains <- object$n.chains
  run.time <- object$run.time[3] / 60 # minutes

  cat(paste("Samples per Chain: ", n.samples,"\n", sep=""))
  cat(paste("Burn-in: ", n.burn,"\n", sep=""))
  cat(paste("Thinning Rate: ",n.thin,"\n", sep=""))
  cat(paste("Number of Chains: ", n.chains, "\n", sep = ""))
  cat(paste("Total Posterior Samples: ",n.post * n.chains,"\n", sep=""))
  cat(paste("Run Time (min): ", round(run.time, digits), "\n\n", sep = ""))

  # Abundance
  if (object$dist %in% c('Poisson', 'NB')) {
    cat("Abundance (log scale): \n")
  }
  if (object$dist %in% c('Gaussian', 'zi-Gaussian')) {
    cat("Abundance: \n")
  }
  tmp.1 <- t(apply(object$beta.samples, 2,
		   function(x) c(mean(x), sd(x))))
  colnames(tmp.1) <- c("Mean", "SD")
  tmp <- t(apply(object$beta.samples, 2,
		 function(x) quantile(x, prob = quantiles)))
  diags <- matrix(c(object$rhat$beta, round(object$ESS$beta, 0)), ncol = 2)
  colnames(diags) <- c('Rhat', 'ESS')

  print(noquote(round(cbind(tmp.1, tmp, diags), digits)))
  if (object$muRE) {
    cat("\n")
    if (object$dist %in% c('Poisson', 'NB')) {
      cat("Abundance Random Effect Variances (log scale): \n")
    }
    if (object$dist %in% c('Gaussian', 'zi-Gaussian')) {
      cat("Abundance Random Effect Variances: \n")
    }
    tmp.1 <- t(apply(object$sigma.sq.mu.samples, 2,
          	   function(x) c(mean(x), sd(x))))
    colnames(tmp.1) <- c("Mean", "SD")
    tmp <- t(apply(object$sigma.sq.mu.samples, 2,
          	 function(x) quantile(x, prob = quantiles)))
    diags <- matrix(c(object$rhat$sigma.sq.mu, round(object$ESS$sigma.sq.mu, 0)), ncol = 2)
    colnames(diags) <- c('Rhat', 'ESS')

    print(noquote(round(cbind(tmp.1, tmp, diags), digits)))
  }
  cat("\n")
  # Covariance ------------------------
  cat("Spatial Covariance: \n")
  tmp.1 <- t(apply(object$theta.samples, 2,
		   function(x) c(mean(x), sd(x))))
  colnames(tmp.1) <- c("Mean", "SD")
  tmp <- t(apply(object$theta.samples, 2,
		 function(x) quantile(x, prob = quantiles)))
  diags <- matrix(c(object$rhat$theta, round(object$ESS$theta, 0)), ncol = 2)
  colnames(diags) <- c('Rhat', 'ESS')
  print(noquote(round(cbind(tmp.1, tmp, diags), digits)))

  if (object$dist == "NB") {
    cat("\n")
    cat("NB overdispersion: \n")
    tmp.1 <- t(apply(object$kappa.samples, 2,
          	   function(x) c(mean(x), sd(x))))
    colnames(tmp.1) <- c("Mean", "SD")
    tmp <- t(apply(object$kappa.samples, 2,
          	 function(x) quantile(x, prob = quantiles)))
    diags <- matrix(c(object$rhat$kappa, round(object$ESS$kappa, 0)), ncol = 2)
    colnames(diags) <- c('Rhat', 'ESS')

    print(noquote(round(cbind(tmp.1, tmp, diags), digits)))
  }
  
  if (object$dist %in% c('Gaussian', 'zi-Gaussian')) {
    cat("\n")
    cat("Residual variance: \n")
    tmp.1 <- t(apply(object$tau.sq.samples, 2,
          	   function(x) c(mean(x), sd(x))))
    colnames(tmp.1) <- c("Mean", "SD")
    tmp <- t(apply(object$tau.sq.samples, 2,
          	 function(x) quantile(x, prob = quantiles)))
    diags <- matrix(c(object$rhat$tau.sq, round(object$ESS$tau.sq, 0)), ncol = 2)
    colnames(diags) <- c('Rhat', 'ESS')

    print(noquote(round(cbind(tmp.1, tmp, diags), digits)))
  }
}

predict.spAbund <- function(object, X.0, coords.0,  
			    n.omp.threads = 1, verbose = TRUE, n.report = 100,
			    ignore.RE = FALSE, z.0.samples, include.sp = TRUE, ...) {

  if (object$dist %in% c('Gaussian', 'zi-Gaussian')) {
    predict.svcAbund(object, X.0, coords.0, n.omp.threads, 
		     verbose, n.report, ignore.RE, z.0.samples)
  } else {
    ptm <- proc.time()
    # Check for unused arguments ------------------------------------------
    formal.args <- names(formals(sys.function(sys.parent())))
    elip.args <- names(list(...))
    for(i in elip.args){
        if(! i %in% formal.args)
            warning("'",i, "' is not an argument")
    }
    # Call ----------------------------------------------------------------
    cl <- match.call()

    # Call predict.abund if don't care about spatial effects
    if (!include.sp) {
      out <- predict.abund(object, X.0, ignore.RE, z.0.samples)
    } else {

      # Functions ---------------------------------------------------------------
      logit <- function(theta, a = 0, b = 1) {log((theta-a)/(b-theta))}
      logit.inv <- function(z, a = 0, b = 1) {b-(b-a)/(1+exp(z))}

      # Some initial checks ---------------------------------------------------
      if (missing(object)) {
        stop("error: predict expects object\n")
      }
      if (!(class(object) %in% c('spAbund'))) {
        stop("error: requires an output object of class spAbund\n")
      }

      if (missing(X.0)) {
        stop("error: X.0 must be specified\n")
      }
      if (!(length(dim(X.0)) %in% c(2, 3))) {
        stop("error: X.0 must be a matrix with two columns corresponding to site and covariate or a three-dimensional array with dimensions corresponding to site, replicate, and covariate")
      }
      if (length(dim(X.0)) == 2) {
        tmp <- colnames(X.0) 
        X.0 <- array(X.0, dim = c(nrow(X.0), 1, ncol(X.0)))
        dimnames(X.0)[[3]] <- tmp
      }
      if (missing(coords.0)) {
        stop("error: coords.0 must be specified\n")
      }
      if (!any(is.data.frame(coords.0), is.matrix(coords.0))) {
        stop("error: coords.0 must be a data.frame or matrix\n")
      }
      if (!ncol(coords.0) == 2){
        stop("error: coords.0 must have two columns\n")
      }
      coords.0 <- as.matrix(coords.0)
      sites.0.indx <- rep(0:(nrow(X.0) - 1), times = ncol(X.0))
      K.max.0 <- ncol(X.0)

      # Get X.0 in long format. 
      tmp.names <- dimnames(X.0)[[3]]
      X.0 <- matrix(X.0, nrow = nrow(X.0) * ncol(X.0), ncol = dim(X.0)[3])
      colnames(X.0) <- tmp.names
      missing.indx <- which(apply(X.0, 1, function(a) sum(is.na(a)) > 0))
      non.missing.indx <- which(apply(X.0, 1, function(a) sum(is.na(a)) == 0))
      X.0 <- X.0[non.missing.indx, , drop = FALSE]

      n.post <- object$n.post * object$n.chains
      X <- object$X
      coords <- object$coords
      n.obs <- sum(!is.na(object$y))
      J <- nrow(coords)
      p.abund <- dim(X)[3]
      theta.samples <- object$theta.samples
      beta.samples <- object$beta.samples
      w.samples <- object$w.samples
      n.neighbors <- object$n.neighbors
      cov.model.indx <- object$cov.model.indx
      re.cols <- object$re.cols
      sp.type <- object$type
      if (object$muRE & !ignore.RE) {
        p.abund.re <- length(object$re.level.names)
      } else {
        p.abund.re <- 0
      }
      # Eliminate prediction sites that have already sampled been for now
      match.indx <- match(do.call("paste", as.data.frame(coords.0)), do.call("paste", as.data.frame(coords)))
      coords.0.indx <- which(is.na(match.indx))
      coords.indx <- match.indx[!is.na(match.indx)]
      coords.place.indx <- which(!is.na(match.indx))

      # if (length(coords.indx) == nrow(X.0)) {
      #   stop("error: no new locations to predict at. See object$psi.samples for occurrence probabilities at sampled sites.")
      # }

      if (object$muRE & !ignore.RE) {
        beta.star.samples <- object$beta.star.samples
        re.level.names <- object$re.level.names
        # Get columns in design matrix with random effects
        x.re.names <- dimnames(object$X.re)[[3]]
        x.0.names <- colnames(X.0)
        # Get the number of times each factor is used. 
        re.long.indx <- sapply(re.cols, length)
        # Need sapply to keep the replicate columns if a factor is used
        # in both a random intercept and random slope. 
        tmp <- sapply(x.re.names, function(a) which(colnames(X.0) %in% a))
        # tmp <- which(colnames(X.0) %in% x.re.names)
        indx <- list()
        for (i in 1:length(tmp)) {
          indx[[i]] <- rep(tmp[i], re.long.indx[i])
        }
        indx <- unlist(indx)
        if (length(indx) == 0) {
          stop("error: column names in X.0 must match variable names in data$occ.covs")
        }
        n.re <- length(indx)
        n.unique.re <- length(unique(indx))
        # Check RE columns
        for (i in 1:n.re) {
          if (is.character(re.cols[[i]])) {
            # Check if all column names in svc are in occ.covs
            if (!all(re.cols[[i]] %in% x.0.names)) {
                missing.cols <- re.cols[[i]][!(re.cols[[i]] %in% x.0.names)]
                stop(paste("error: variable name ", paste(missing.cols, collapse=" and "), " not in occurrence covariates", sep=""))
            }
            # Convert desired column names into the numeric column index
            re.cols[[i]] <- which(x.0.names %in% re.cols[[i]])
            
          } else if (is.numeric(re.cols[[i]])) {
            # Check if all column indices are in 1:p.abund
            if (!all(re.cols %in% 1:p.abund)) {
                missing.cols <- re.cols[[i]][!(re.cols[[i]] %in% (1:p.abund))]
                stop(paste("error: column index ", paste(missing.cols, collapse=" "), " not in design matrix columns", sep=""))
            }
          }
        }
        re.cols <- unlist(re.cols)
        X.re <- as.matrix(X.0[, indx, drop = FALSE])
        X.fix <- as.matrix(X.0[, -indx, drop = FALSE])
        X.random <- as.matrix(X.0[, re.cols, drop = FALSE])
        n.abund.re <- length(unlist(re.level.names))
        X.re.ind <- matrix(NA, nrow(X.re), p.abund.re)

        for (i in 1:p.abund.re) {
          for (j in 1:nrow(X.re)) {
            tmp <- which(re.level.names[[i]] == X.re[j, i])
            if (length(tmp) > 0) {
              X.re.ind[j, i] <- tmp 
            }
          }
        }
        if (p.abund.re > 1) {
          for (j in 2:p.abund.re) {
            X.re.ind[, j] <- X.re.ind[, j] + max(X.re.ind[, j - 1]) 
          }
        }
        # Create the random effects corresponding to each 
        # new location
        # ORDER: ordered by site, then species within site.
        beta.star.sites.0.samples <- matrix(0, n.post,  nrow(X.re))
        for (t in 1:p.abund.re) {
          for (j in 1:nrow(X.re)) {
            if (!is.na(X.re.ind[j, t])) {
              beta.star.sites.0.samples[, j] <- 
                beta.star.samples[, X.re.ind[j, t]] * X.random[j, t] + 
                beta.star.sites.0.samples[, j]
            } else {
              beta.star.sites.0.samples[, j] <- 
                rnorm(n.post, 0, sqrt(object$sigma.sq.mu.samples[, t])) * X.random[j, t] + 
                beta.star.sites.0.samples[, j]
            }
          } # j
        } # t
      } else {
        X.fix <- X.0
        beta.star.sites.0.samples <- matrix(0, n.post, nrow(X.0))
        p.re <- 0
      }
      # Sub-sample previous
      beta.samples <- t(beta.samples)
      w.samples <- t(w.samples)
      beta.star.sites.0.samples <- t(beta.star.sites.0.samples)
      if (object$dist == 'NB') {
        kappa.samples <- t(object$kappa.samples)
      } else {
        kappa.samples <- matrix(0, 1, n.post)
      }
      family.c <- ifelse(object$dist == 'NB', 1, 0)
      theta.samples <- t(theta.samples)

      n.obs.0 <- nrow(X.fix)
      J.0 <- length(unique(sites.0.indx))


      # Indicates whether a site has been sampled. 1 = sampled
      sites.0.sampled <- ifelse(!is.na(match.indx), 1, 0)
      sites.link <- rep(NA, J.0)
      sites.link[which(!is.na(match.indx))] <- coords.indx
      # For C
      sites.link <- sites.link - 1

      if (sp.type == 'GP') {
        stop("NNGP = FALSE is not currently supported for spAbund")
      } else {
        # Get nearest neighbors
        # nn2 is a function from RANN.
        nn.indx.0 <- nn2(coords, coords.0, k=n.neighbors)$nn.idx-1

        storage.mode(coords) <- "double"
        storage.mode(J) <- "integer"
        storage.mode(n.obs) <- "integer"
        storage.mode(p.abund) <- "integer"
        storage.mode(n.neighbors) <- "integer"
        storage.mode(X.fix) <- "double"
        storage.mode(coords.0) <- "double"
        storage.mode(J.0) <- "integer"
        storage.mode(n.obs.0) <- "integer"
        storage.mode(sites.link) <- "integer"
        storage.mode(sites.0.sampled) <- "integer"
        storage.mode(sites.0.indx) <- "integer"
        storage.mode(beta.samples) <- "double"
        storage.mode(theta.samples) <- "double"
        storage.mode(kappa.samples) <- "double"
        storage.mode(family.c) <- "integer"
        storage.mode(w.samples) <- "double"
        storage.mode(beta.star.sites.0.samples) <- "double"
        storage.mode(n.post) <- "integer"
        storage.mode(cov.model.indx) <- "integer"
        storage.mode(nn.indx.0) <- "integer"
        storage.mode(n.omp.threads) <- "integer"
        storage.mode(verbose) <- "integer"
        storage.mode(n.report) <- "integer"

        ptm <- proc.time()

        out <- .Call("spAbundNNGPPredict", coords, J, n.obs, p.abund, n.neighbors,
                     X.fix, coords.0, J.0, n.obs.0, sites.link, sites.0.sampled, sites.0.indx, 
            	 nn.indx.0, beta.samples,
                     theta.samples, w.samples, beta.star.sites.0.samples, kappa.samples, 
            	 n.post, cov.model.indx, n.omp.threads, verbose, n.report, family.c)
      }

      tmp <- matrix(NA, n.post, length(c(missing.indx, non.missing.indx)))
      tmp[, non.missing.indx] <- t(out$y.0.samples)
      out$y.0.samples <- array(tmp, dim = c(n.post, J.0, K.max.0))
      tmp <- matrix(NA, n.post, length(c(missing.indx, non.missing.indx)))
      tmp[, non.missing.indx] <- t(out$mu.0.samples)
      out$mu.0.samples <- array(tmp, dim = c(n.post, J.0, K.max.0))
      # This gives the same result, just switched the order
      # tmp <- matrix(NA, length(c(missing.indx, non.missing.indx)), n.post)
      # tmp[non.missing.indx, ] <- out$y.0.samples
      # out$y.0.samples <- array(tmp, dim = c(J.0, K.max.0, n.post))
      # tmp <- matrix(NA, length(c(missing.indx, non.missing.indx)), n.post)
      # tmp[non.missing.indx, ] <- out$mu.0.samples
      # out$mu.0.samples <- array(tmp, dim = c(J.0, K.max.0, n.post))
      out$w.0.samples <- mcmc(t(out$w.0.samples))
      out$run.time <- proc.time() - ptm
    }
    out$call <- cl
    out$object.class <- class(object)
    class(out) <- "predict.spAbund"
    out
  }
}

# msAbund -----------------------------------------------------------------
print.msAbund <- function(x, ...) {
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.75)),
      "", sep = "\n")
}

summary.msAbund <- function(object,
			    level = 'both',
			    quantiles = c(0.025, 0.5, 0.975),
			    digits = max(3L, getOption("digits") - 3L), ...) {

  print(object)

  n.post <- object$n.post
  n.samples <- object$n.samples
  n.burn <- object$n.burn
  n.thin <- object$n.thin
  n.chains <- object$n.chains
  run.time <- object$run.time[3] / 60 # minutes

  cat(paste("Samples per Chain: ", n.samples,"\n", sep=""))
  cat(paste("Burn-in: ", n.burn,"\n", sep=""))
  cat(paste("Thinning Rate: ",n.thin,"\n", sep=""))
  cat(paste("Number of Chains: ", n.chains, "\n", sep = ""))
  cat(paste("Total Posterior Samples: ",n.post * n.chains,"\n", sep=""))
  cat(paste("Run Time (min): ", round(run.time, digits), "\n\n", sep = ""))

  if (tolower(level) %in% c('community', 'both')) {

    cat("----------------------------------------\n");
    cat("\tCommunity Level\n");
    cat("----------------------------------------\n");

    # Abundance 
    # Abundance 
    if (object$dist %in% c('Poisson', 'NB')) {
      cat("Abundance Means (log scale): \n")
    }
    if (object$dist %in% c('Gaussian', 'zi-Gaussian')) {
      cat("Abundance Means: \n")
    }
    tmp.1 <- t(apply(object$beta.comm.samples, 2,
          	   function(x) c(mean(x), sd(x))))
    colnames(tmp.1) <- c("Mean", "SD")
    tmp <- t(apply(object$beta.comm.samples, 2,
          	 function(x) quantile(x, prob = quantiles)))
    diags <- matrix(c(object$rhat$beta.comm, round(object$ESS$beta.comm, 0)), ncol = 2)
    colnames(diags) <- c('Rhat', 'ESS')

    print(noquote(round(cbind(tmp.1, tmp, diags), digits)))

    if (object$dist %in% c('Poisson', 'NB')) {
      cat("\nAbundance Variances (log scale): \n")
    }
    if (object$dist %in% c('Gaussian', 'zi-Gaussian')) {
      cat("\nAbundance Variances: \n")
    }
    tmp.1 <- t(apply(object$tau.sq.beta.samples, 2,
          	   function(x) c(mean(x), sd(x))))
    colnames(tmp.1) <- c("Mean", "SD")
    tmp <- t(apply(object$tau.sq.beta.samples, 2,
          	 function(x) quantile(x, prob = quantiles)))
    diags <- matrix(c(object$rhat$tau.sq.beta, round(object$ESS$tau.sq.beta, 0)), ncol = 2)
    colnames(diags) <- c('Rhat', 'ESS')
    print(noquote(round(cbind(tmp.1, tmp, diags), digits)))
    if (object$muRE) {
      cat("\n")
      if (object$dist %in% c('Poisson', 'NB')) {
        cat("Abundance Random Effect Variances (log scale): \n")
      }
      if (object$dist %in% c('Gaussian', 'zi-Gaussian')) {
        cat("Abundance Random Effect Variances: \n")
      }
      tmp.1 <- t(apply(object$sigma.sq.mu.samples, 2,
            	   function(x) c(mean(x), sd(x))))
      colnames(tmp.1) <- c("Mean", "SD")
      tmp <- t(apply(object$sigma.sq.mu.samples, 2,
            	 function(x) quantile(x, prob = quantiles)))
      diags <- matrix(c(object$rhat$sigma.sq.mu, round(object$ESS$sigma.sq.mu, 0)), ncol = 2)
      colnames(diags) <- c('Rhat', 'ESS')

      print(noquote(round(cbind(tmp.1, tmp, diags), digits)))
    }
  }

  if (tolower(level) %in% c('species', 'both')) {
    if (tolower(level) == 'both') cat("\n")
    cat("----------------------------------------\n");
    cat("\tSpecies Level\n");
    cat("----------------------------------------\n");
    if (object$dist %in% c('Poisson', 'NB')) {
      cat("Abundance (log scale): \n")
    }
    if (object$dist %in% c('Gaussian', 'zi-Gaussian')) {
      cat("Abundance: \n")
    }
    tmp.1 <- t(apply(object$beta.samples, 2,
          	   function(x) c(mean(x), sd(x))))
    colnames(tmp.1) <- c("Mean", "SD")
    tmp <- t(apply(object$beta.samples, 2,
          	 function(x) quantile(x, prob = quantiles)))
    diags <- matrix(c(object$rhat$beta, round(object$ESS$beta, 0)), ncol = 2)
    colnames(diags) <- c('Rhat', 'ESS')

    print(noquote(round(cbind(tmp.1, tmp, diags), digits)))
    
    # NB Overdispersion
    if (object$dist == "NB") {
      cat("\n")
      cat("----------------------------------------\n");
      cat("\tNB overdispersion\n");
      cat("----------------------------------------\n");
      tmp.1 <- t(apply(object$kappa.samples, 2,
            	   function(x) c(mean(x), sd(x))))
      colnames(tmp.1) <- c("Mean", "SD")
      tmp <- t(apply(object$kappa.samples, 2,
            	 function(x) quantile(x, prob = quantiles)))
      diags <- matrix(c(object$rhat$kappa, round(object$ESS$kappa, 0)), ncol = 2)
      colnames(diags) <- c('Rhat', 'ESS')
      print(noquote(round(cbind(tmp.1, tmp, diags), digits)))
    }
  }

  if (object$dist %in% c('Gaussian', 'zi-Gaussian')) {
    cat("\n")
    cat("----------------------------------------\n");
    cat("\tResidual variance\n");
    cat("----------------------------------------\n");
    tmp.1 <- t(apply(object$tau.sq.samples, 2,
          	   function(x) c(mean(x), sd(x))))
    colnames(tmp.1) <- c("Mean", "SD")
    tmp <- t(apply(object$tau.sq.samples, 2,
          	 function(x) quantile(x, prob = quantiles)))
    diags <- matrix(c(object$rhat$tau.sq, round(object$ESS$tau.sq, 0)), ncol = 2)
    colnames(diags) <- c('Rhat', 'ESS')

    print(noquote(round(cbind(tmp.1, tmp, diags), digits)))
  }
}

fitted.msAbund <- function(object, ...) {
  return(object$y.rep.samples)
}

predict.msAbund <- function(object, X.0, ignore.RE = FALSE, z.0.samples, ...) {
  # Check for unused arguments ------------------------------------------
  formal.args <- names(formals(sys.function(sys.parent())))
  elip.args <- names(list(...))
  for(i in elip.args){
      if(! i %in% formal.args)
          warning("'",i, "' is not an argument")
  }
  # Call ----------------------------------------------------------------
  cl <- match.call()

  # Functions ---------------------------------------------------------------
  logit <- function(theta, a = 0, b = 1) {log((theta-a)/(b-theta))}
  logit.inv <- function(z, a = 0, b = 1) {b-(b-a)/(1+exp(z))}

  # Some initial checks ---------------------------------------------------
  if (missing(object)) {
    stop("error: predict expects object\n")
  }
  if (!(class(object) %in% c('msAbund', 'sfMsAbund', 'lfMsAbund'))) {
    stop("error: requires an output object of class msAbund, lfMsAbund, sfMsAbund\n")
  }

  # Check X.0 -------------------------------------------------------------
  if (missing(X.0)) {
    stop("error: X.0 must be specified\n")
  }
  if (!(length(dim(X.0)) %in% c(2, 3))) {
    stop("error: X.0 must be a matrix with two columns corresponding to site and covariate or a three-dimensional array with dimensions corresponding to site, replicate, and covariate")
  }
  if (length(dim(X.0)) == 2) {
    tmp <- colnames(X.0) 
    X.0 <- array(X.0, dim = c(nrow(X.0), 1, ncol(X.0)))
    dimnames(X.0)[[3]] <- tmp
  }
  # Predictions -----------------------------------------------------------
  if (object$dist %in% c('Gaussian', 'zi-Gaussian')) {
    p.abund <- dim(object$X)[2]
  } else {
    p.abund <- dim(object$X)[3]
  }
  # p.design <- p.abund
  # if (object$muRE & !ignore.RE) {
  #   p.design <- p.abund + ncol(object$sigma.sq.mu.samples)
  # }
  # if (dim(X.0)[3] != p.design) {
  #   stop(paste("error: X.0 must have ", p.design, " covariates\n", sep = ''))
  # }
  J.0 <- nrow(X.0)
  K.max.0 <- ncol(X.0)
  n.sp <- dim(object$y)[1]

  # Composition sampling --------------------------------------------------
  re.cols <- object$re.cols
  beta.samples <- as.matrix(object$beta.samples)
  if (object$dist == 'NB') {
    kappa.samples <- as.matrix(object$kappa.samples)
  }
  n.post <- object$n.post * object$n.chains
  out <- list()
  if (object$muRE) {
    p.abund.re <- length(object$re.level.names)
  } else {
    p.abund.re <- 0
  }

  # Get X.0 in long format. 
  tmp.names <- dimnames(X.0)[[3]]
  X.0 <- matrix(X.0, nrow = nrow(X.0) * ncol(X.0), ncol = dim(X.0)[3])
  colnames(X.0) <- tmp.names
  missing.indx <- which(apply(X.0, 1, function(a) sum(is.na(a)) > 0))
  non.missing.indx <- which(apply(X.0, 1, function(a) sum(is.na(a)) == 0))
  X.0 <- X.0[non.missing.indx, , drop = FALSE]

  if (object$muRE & !ignore.RE) {
    beta.star.samples <- object$beta.star.samples
    re.level.names <- object$re.level.names
    # Get columns in design matrix with random effects
    if (object$dist %in% c('Gaussian', 'zi-Gaussian')) {
      x.re.names <- dimnames(object$X.re)[[2]]
    } else {
      x.re.names <- dimnames(object$X.re)[[3]]
    }
    x.0.names <- colnames(X.0)
    # Get the number of times each factor is used. 
    re.long.indx <- sapply(re.cols, length)
    tmp <- sapply(x.re.names, function(a) which(colnames(X.0) %in% a))
    indx <- list()
    for (i in 1:length(tmp)) {
      indx[[i]] <- rep(tmp[i], re.long.indx[i])
    }
    indx <- unlist(indx)
    if (length(indx) == 0) {
      stop("error: column names in X.0 must match variable names in data$occ.covs")
    }
    n.re <- length(indx)
    n.unique.re <- length(unique(indx))
    # Check RE columns
    for (i in 1:n.re) {
      if (is.character(re.cols[[i]])) {
        # Check if all column names in svc are in occ.covs
        if (!all(re.cols[[i]] %in% x.0.names)) {
            missing.cols <- re.cols[[i]][!(re.cols[[i]] %in% x.0.names)]
            stop(paste("error: variable name ", paste(missing.cols, collapse=" and "), " not in occurrence covariates", sep=""))
        }
        # Convert desired column names into the numeric column index
        re.cols[[i]] <- which(x.0.names %in% re.cols[[i]])
        
      } else if (is.numeric(re.cols[[i]])) {
        # Check if all column indices are in 1:p.abund
        if (!all(re.cols %in% 1:p.abund)) {
            missing.cols <- re.cols[[i]][!(re.cols[[i]] %in% (1:p.abund))]
            stop(paste("error: column index ", paste(missing.cols, collapse=" "), " not in design matrix columns", sep=""))
        }
      }
    }
    re.cols <- unlist(re.cols)
    X.re <- as.matrix(X.0[, indx, drop = FALSE])
    X.fix <- as.matrix(X.0[, -indx, drop = FALSE])
    X.random <- as.matrix(X.0[, re.cols, drop = FALSE])
    n.abund.re <- length(unlist(re.level.names))
    X.re.ind <- matrix(NA, nrow(X.re), p.abund.re)

    for (i in 1:p.abund.re) {
      for (j in 1:nrow(X.re)) {
        tmp <- which(re.level.names[[i]] == X.re[j, i])
        if (length(tmp) > 0) {
          X.re.ind[j, i] <- tmp 
        }
      }
    }
    if (p.abund.re > 1) {
      for (j in 2:p.abund.re) {
        X.re.ind[, j] <- X.re.ind[, j] + max(X.re.ind[, j - 1]) 
      }
    }
    # Create the random effects corresponding to each 
    # new location
    beta.star.sites.0.samples <- array(0, dim = c(n.post, n.sp, nrow(X.re)))
    for (i in 1:n.sp) {
      for (t in 1:p.abund.re) {
        for (j in 1:nrow(X.re)) {
          if (!is.na(X.re.ind[j, t])) {
            beta.star.sites.0.samples[, i, j] <- 
              beta.star.samples[, (i - 1) * n.abund.re + X.re.ind[j, t]] * X.random[j, t] + 
              beta.star.sites.0.samples[, i, j]
          } else {
            beta.star.sites.0.samples[, i, j] <- 
              rnorm(n.post, 0, sqrt(object$sigma.sq.mu.samples[, t])) * X.random[j, t] + 
              beta.star.sites.0.samples[, i, j]
          }
        } # j
      } # t
    }
  } else {
    X.fix <- X.0
    beta.star.sites.0.samples <- array(0, dim = c(n.post, n.sp, nrow(X.0)))
    p.abund.re <- 0
  }
  family <- object$dist

  if (family == 'zi-Gaussian') {
    if (missing(z.0.samples)) {
      stop("z.0.samples must be supplied for a zi-Gaussian model")
    }
    if (length(dim(z.0.samples)) != 3) {
      stop(paste("z.0.samples must be a three-dimensional array with dimensions of ", 
		 n.post, ", ", n.sp, ", and ", J.0, sep = ''))
    }
    if (dim(z.0.samples)[1] != n.post | dim(z.0.samples)[2] != n.sp |
	dim(z.0.samples)[3] != J.0) {
      stop(paste("z.0.samples must be a three-dimensional array with dimensions of ", 
		 n.post, ", ", n.sp, ", and ", J.0, sep = ''))
    }
  } else {
    if (!missing(z.0.samples)) {
      message("z.0.samples is ignored for the current model family\n")
    }
    z.0.samples <- array(NA, dim = c(1, 1, 1))
  }

  tmp <- matrix(NA, n.post, length(c(missing.indx, non.missing.indx)))
  sp.indx <- rep(1:n.sp, p.abund)
  out$mu.0.samples <- array(NA, dim = c(n.post, n.sp, J.0, K.max.0))
  mu.long <- array(NA, dim = c(n.post, n.sp, nrow(X.fix)))
  for (i in 1:n.sp) {
    if (object$dist %in% c('Poisson', 'NB')) {
      mu.long[, i, ] <- exp(t(X.fix %*% t(beta.samples[, sp.indx == i]) + 
	                    t(beta.star.sites.0.samples[, i, ])))
    }
    if (object$dist %in% c('Gaussian', 'zi-Gaussian')) {
      mu.long[, i, ] <- t(X.fix %*% t(beta.samples[, sp.indx == i]) + 
	                    t(beta.star.sites.0.samples[, i, ]))
    }
    tmp[, non.missing.indx] <- mu.long[, i, ]
    out$mu.0.samples[, i, , ] <- array(tmp, dim = c(n.post, J.0, K.max.0))
  }
  K <- apply(out$mu.0.samples[1, 1, , , drop = FALSE], 3, function(a) sum(!is.na(a)))
  out$y.0.samples <- array(NA, dim(out$mu.0.samples))
  for (i in 1:n.sp) {
    for (j in 1:J.0) {
      for (k in 1:K.max.0) {
        if (sum(is.na(out$mu.0.samples[, i, j, k])) == 0) {
          if (object$dist == 'NB') {
            out$y.0.samples[, i, j, k] <- rnbinom(n.post, kappa.samples[, i], 
            				        mu = out$mu.0.samples[, i, j, k])
          } 
          if (object$dist == 'Poisson') {
            out$y.0.samples[, i, j, k] <- rpois(n.post, out$mu.0.samples[, i, j, k])
          }
	  if (object$dist == 'Gaussian') {
            out$y.0.samples[, i, j, k] <- rnorm(n.post, out$mu.0.samples[, i, j, k], 
	                                        sqrt(object$tau.sq.samples[, i]))
	  }
	  if (object$dist == 'zi-Gaussian') {
            out$y.0.samples[, i, j, k] <- ifelse(z.0.samples[, i, j] == 1, 
	    				  rnorm(n.post, out$mu.0.samples[, i, j, k], 
	    					sqrt(object$tau.sq.samples[, i])), 
	    				  rnorm(n.post, 0, sqrt(0.0001)))
            out$mu.0.samples[, i, j, k] <- ifelse(z.0.samples[, i, j] == 1, 
	  				   out$mu.0.samples[, i, j, k], 0)
	  }
	}
      }
    }
  }
  out$call <- cl
  # If Gaussian, collapse to a 3D array 
  if (object$dist %in% c('Gaussian', 'zi-Gaussian')) {
    out$y.0.samples <- out$y.0.samples[, , , 1]  
    out$mu.0.samples <- out$mu.0.samples[, , , 1]  
  }

  class(out) <- "predict.msAbund"
  out
}

# sfMsAbund ---------------------------------------------------------------
print.sfMsAbund <- function(x, ...) {
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.75)),
      "", sep = "\n")
}

summary.sfMsAbund <- function(object, level = 'both',
                              quantiles = c(0.025, 0.5, 0.975),
                              digits = max(3L, getOption("digits") - 3L), ...) {

  print(object)

  n.post <- object$n.post
  n.samples <- object$n.samples
  n.burn <- object$n.burn
  n.thin <- object$n.thin
  n.chains <- object$n.chains
  run.time <- object$run.time[3] / 60 # minutes

  cat(paste("Samples per Chain: ", n.samples,"\n", sep=""))
  cat(paste("Burn-in: ", n.burn,"\n", sep=""))
  cat(paste("Thinning Rate: ",n.thin,"\n", sep=""))
  cat(paste("Number of Chains: ", n.chains, "\n", sep = ""))
  cat(paste("Total Posterior Samples: ",n.post * n.chains,"\n", sep=""))
  cat(paste("Run Time (min): ", round(run.time, digits), "\n\n", sep = ""))

  if (tolower(level) %in% c('community', 'both')) {

    cat("----------------------------------------\n");
    cat("\tCommunity Level\n");
    cat("----------------------------------------\n");

    # Abundance 
    if (object$dist %in% c('Poisson', 'NB')) {
      cat("Abundance Means (log scale): \n")
    }
    if (object$dist %in% c('Gaussian', 'zi-Gaussian')) {
      cat("Abundance Means: \n")
    }
    tmp.1 <- t(apply(object$beta.comm.samples, 2,
          	   function(x) c(mean(x), sd(x))))
    colnames(tmp.1) <- c("Mean", "SD")
    tmp <- t(apply(object$beta.comm.samples, 2,
          	 function(x) quantile(x, prob = quantiles)))
    diags <- matrix(c(object$rhat$beta.comm, round(object$ESS$beta.comm, 0)), ncol = 2)
    colnames(diags) <- c('Rhat', 'ESS')

    print(noquote(round(cbind(tmp.1, tmp, diags), digits)))

    if (object$dist %in% c('Poisson', 'NB')) {
      cat("\nAbundance Variances (log scale): \n")
    }
    if (object$dist %in% c('Gaussian', 'zi-Gaussian')) {
      cat("\nAbundance Variances: \n")
    }
    tmp.1 <- t(apply(object$tau.sq.beta.samples, 2,
          	   function(x) c(mean(x), sd(x))))
    colnames(tmp.1) <- c("Mean", "SD")
    tmp <- t(apply(object$tau.sq.beta.samples, 2,
          	 function(x) quantile(x, prob = quantiles)))
    diags <- matrix(c(object$rhat$tau.sq.beta, round(object$ESS$tau.sq.beta, 0)), ncol = 2)
    colnames(diags) <- c('Rhat', 'ESS')
    print(noquote(round(cbind(tmp.1, tmp, diags), digits)))

    if (object$muRE) {
      cat("\n")
      if (object$dist %in% c('Poisson', 'NB')) {
        cat("Abundance Random Effect Variances (log scale): \n")
      }
      if (object$dist %in% c('Gaussian', 'zi-Gaussian')) {
        cat("Abundance Random Effect Variances: \n")
      }
      tmp.1 <- t(apply(object$sigma.sq.mu.samples, 2,
            	   function(x) c(mean(x), sd(x))))
      colnames(tmp.1) <- c("Mean", "SD")
      tmp <- t(apply(object$sigma.sq.mu.samples, 2,
            	 function(x) quantile(x, prob = quantiles)))
      diags <- matrix(c(object$rhat$sigma.sq.mu, round(object$ESS$sigma.sq.mu, 0)), ncol = 2)
      colnames(diags) <- c('Rhat', 'ESS')

      print(noquote(round(cbind(tmp.1, tmp, diags), digits)))
    }
  }

  if (tolower(level) %in% c('species', 'both')) {
    if (tolower(level) == 'both') cat("\n")
    cat("----------------------------------------\n");
    cat("\tSpecies Level\n");
    cat("----------------------------------------\n");
    if (object$dist %in% c('Poisson', 'NB')) {
      cat("Abundance (log scale): \n")
    }
    if (object$dist %in% c('Gaussian', 'zi-Gaussian')) {
      cat("Abundance: \n")
    }
    tmp.1 <- t(apply(object$beta.samples, 2,
          	   function(x) c(mean(x), sd(x))))
    colnames(tmp.1) <- c("Mean", "SD")
    tmp <- t(apply(object$beta.samples, 2,
          	 function(x) quantile(x, prob = quantiles)))
    diags <- matrix(c(object$rhat$beta, round(object$ESS$beta, 0)), ncol = 2)
    colnames(diags) <- c('Rhat', 'ESS')

    print(noquote(round(cbind(tmp.1, tmp, diags), digits)))
  }
  # Covariance
  cat("\n")
  cat("----------------------------------------\n");
  cat("\tSpatial Covariance\n");
  cat("----------------------------------------\n");
  tmp.1 <- t(apply(object$theta.samples, 2,
        	   function(x) c(mean(x), sd(x))))
  colnames(tmp.1) <- c("Mean", "SD")
  tmp <- t(apply(object$theta.samples, 2,
        	 function(x) quantile(x, prob = quantiles)))
  diags <- matrix(c(object$rhat$theta, round(object$ESS$theta, 0)), ncol = 2)
  colnames(diags) <- c('Rhat', 'ESS')
  print(noquote(round(cbind(tmp.1, tmp, diags), digits)))

  # NB Overdispersion
  if (object$dist == "NB" & tolower(level) %in% c('species', 'both')) {
    cat("\n")
    cat("----------------------------------------\n");
    cat("\tNB overdispersion\n");
    cat("----------------------------------------\n");
    tmp.1 <- t(apply(object$kappa.samples, 2,
          	   function(x) c(mean(x), sd(x))))
    colnames(tmp.1) <- c("Mean", "SD")
    tmp <- t(apply(object$kappa.samples, 2,
          	 function(x) quantile(x, prob = quantiles)))
    diags <- matrix(c(object$rhat$kappa, round(object$ESS$kappa, 0)), ncol = 2)
    colnames(diags) <- c('Rhat', 'ESS')
    print(noquote(round(cbind(tmp.1, tmp, diags), digits)))
  }

  if (object$dist %in% c('Gaussian', 'zi-Gaussian')) {
    cat("\n")
    cat("----------------------------------------\n");
    cat("\tResidual variance\n");
    cat("----------------------------------------\n");
    tmp.1 <- t(apply(object$tau.sq.samples, 2,
          	   function(x) c(mean(x), sd(x))))
    colnames(tmp.1) <- c("Mean", "SD")
    tmp <- t(apply(object$tau.sq.samples, 2,
          	 function(x) quantile(x, prob = quantiles)))
    diags <- matrix(c(object$rhat$tau.sq, round(object$ESS$tau.sq, 0)), ncol = 2)
    colnames(diags) <- c('Rhat', 'ESS')

    print(noquote(round(cbind(tmp.1, tmp, diags), digits)))
  }
}

fitted.sfMsAbund <- function(object, ...) {
  return(object$y.rep.samples)
}


predict.sfMsAbund <- function(object, X.0, coords.0, n.omp.threads = 1,
			      verbose = TRUE, n.report = 100,
			      ignore.RE = FALSE, z.0.samples, include.sp = TRUE, ...) {

  if (object$dist %in% c('Gaussian', 'zi-Gaussian')) {
    predict.svcMsAbund(object, X.0, coords.0, n.omp.threads, 
		       verbose, n.report, ignore.RE, z.0.samples)
  } else {

    ptm <- proc.time()
    # Check for unused arguments ------------------------------------------
    formal.args <- names(formals(sys.function(sys.parent())))
    elip.args <- names(list(...))
    for(i in elip.args){
        if(! i %in% formal.args)
            warning("'",i, "' is not an argument")
    }
    # Call ----------------------------------------------------------------
    cl <- match.call()

    # Call predict.msAbund if don't care about spatial effects
    if (!include.sp) {
      out <- predict.msAbund(object, X.0, ignore.RE, z.0.samples)
    } else {

      # Functions ---------------------------------------------------------------
      logit <- function(theta, a = 0, b = 1) {log((theta-a)/(b-theta))}
      logit.inv <- function(z, a = 0, b = 1) {b-(b-a)/(1+exp(z))}

      # Some initial checks ---------------------------------------------------
      if (missing(object)) {
        stop("error: predict expects object\n")
      }
      if (!(class(object) %in% c('sfMsAbund'))) {
        stop("error: requires an output object of class sfMsAbund\n")
      }

      # Check X.0 -------------------------------------------------------------
      if (missing(X.0)) {
        stop("error: X.0 must be specified\n")
      }
      if (!(length(dim(X.0)) %in% c(2, 3))) {
        stop("error: X.0 must be a matrix with two columns corresponding to site and covariate or a three-dimensional array with dimensions corresponding to site, replicate, and covariate")
      }
      if (length(dim(X.0)) == 2) {
        tmp <- colnames(X.0) 
        X.0 <- array(X.0, dim = c(nrow(X.0), 1, ncol(X.0)))
        dimnames(X.0)[[3]] <- tmp
      }
      if (missing(coords.0)) {
        stop("error: coords.0 must be specified\n")
      }
      if (!any(is.data.frame(coords.0), is.matrix(coords.0))) {
        stop("error: coords.0 must be a data.frame or matrix\n")
      }
      if (!ncol(coords.0) == 2){
        stop("error: coords.0 must have two columns\n")
      }
      coords.0 <- as.matrix(coords.0)
      sites.0.indx <- rep(0:(nrow(X.0) - 1), times = ncol(X.0))
      K.max.0 <- ncol(X.0)

      # Get X.0 in long format. 
      tmp.names <- dimnames(X.0)[[3]]
      X.0 <- matrix(X.0, nrow = nrow(X.0) * ncol(X.0), ncol = dim(X.0)[3])
      colnames(X.0) <- tmp.names
      missing.indx <- which(apply(X.0, 1, function(a) sum(is.na(a)) > 0))
      non.missing.indx <- which(apply(X.0, 1, function(a) sum(is.na(a)) == 0))
      X.0 <- X.0[non.missing.indx, , drop = FALSE]

      # Abundance predictions ------------------------------------------------
      n.post <- object$n.post * object$n.chains
      X <- object$X
      y <- object$y
      coords <- object$coords
      J <- nrow(coords)
      n.obs <- sum(!is.na(object$y[1, , ]))
      n.sp <- dim(y)[1]
      q <- object$q
      p.abund <- dim(X)[3]
      theta.samples <- object$theta.samples
      beta.samples <- object$beta.samples
      lambda.samples <- object$lambda.samples
      w.samples <- object$w.samples
      kappa.samples <- object$kappa.samples
      n.neighbors <- object$n.neighbors
      cov.model.indx <- object$cov.model.indx
      family <- object$dist
      family.c <- ifelse(family == 'NB', 1, 0)
      sp.type <- object$type
      if (object$muRE & !ignore.RE) {
        p.abund.re <- length(object$re.level.names)
      } else {
        p.abund.re <- 0
      }
      re.cols <- object$re.cols

      # if (ncol(X.0) != p.abund + p.abund.re){
      #   stop(paste("error: X.0 must have ", p.abund + p.abund.re," columns\n"))
      # }
      X.0 <- as.matrix(X.0)

      if (missing(coords.0)) {
        stop("error: coords.0 must be specified\n")
      }
      if (!any(is.data.frame(coords.0), is.matrix(coords.0))) {
        stop("error: coords.0 must be a data.frame or matrix\n")
      }
      if (!ncol(coords.0) == 2){
        stop("error: coords.0 must have two columns\n")
      }
      coords.0 <- as.matrix(coords.0)

      re.cols <- object$re.cols

      if (object$muRE & !ignore.RE) {
        beta.star.samples <- object$beta.star.samples
        re.level.names <- object$re.level.names
        # Get columns in design matrix with random effects
        x.re.names <- dimnames(object$X.re)[[3]]
        x.0.names <- colnames(X.0)
        # Get the number of times each factor is used. 
        re.long.indx <- sapply(re.cols, length)
        tmp <- sapply(x.re.names, function(a) which(colnames(X.0) %in% a))
        indx <- list()
        for (i in 1:length(tmp)) {
          indx[[i]] <- rep(tmp[i], re.long.indx[i])
        }
        indx <- unlist(indx)
        if (length(indx) == 0) {
          stop("error: column names in X.0 must match variable names in data$occ.covs")
        }
        n.re <- length(indx)
        n.unique.re <- length(unique(indx))
        # Check RE columns
        for (i in 1:n.re) {
          if (is.character(re.cols[[i]])) {
            # Check if all column names in svc are in occ.covs
            if (!all(re.cols[[i]] %in% x.0.names)) {
                missing.cols <- re.cols[[i]][!(re.cols[[i]] %in% x.0.names)]
                stop(paste("error: variable name ", paste(missing.cols, collapse=" and "), " not in occurrence covariates", sep=""))
            }
            # Convert desired column names into the numeric column index
            re.cols[[i]] <- which(x.0.names %in% re.cols[[i]])
            
          } else if (is.numeric(re.cols[[i]])) {
            # Check if all column indices are in 1:p.abund
            if (!all(re.cols %in% 1:p.abund)) {
                missing.cols <- re.cols[[i]][!(re.cols[[i]] %in% (1:p.abund))]
                stop(paste("error: column index ", paste(missing.cols, collapse=" "), " not in design matrix columns", sep=""))
            }
          }
        }
        re.cols <- unlist(re.cols)
        X.re <- as.matrix(X.0[, indx, drop = FALSE])
        X.fix <- as.matrix(X.0[, -indx, drop = FALSE])
        X.random <- as.matrix(X.0[, re.cols, drop = FALSE])
        n.abund.re <- length(unlist(re.level.names))
        X.re.ind <- matrix(NA, nrow(X.re), p.abund.re)

        for (i in 1:p.abund.re) {
          for (j in 1:nrow(X.re)) {
            tmp <- which(re.level.names[[i]] == X.re[j, i])
            if (length(tmp) > 0) {
              X.re.ind[j, i] <- tmp 
            }
          }
        }
        if (p.abund.re > 1) {
          for (j in 2:p.abund.re) {
            X.re.ind[, j] <- X.re.ind[, j] + max(X.re.ind[, j - 1]) 
          }
        }
        # Create the random effects corresponding to each 
        # new location
        beta.star.sites.0.samples <- array(0, dim = c(n.post, n.sp, nrow(X.re)))
        for (i in 1:n.sp) {
          for (t in 1:p.abund.re) {
            for (j in 1:nrow(X.re)) {
              if (!is.na(X.re.ind[j, t])) {
                beta.star.sites.0.samples[, i, j] <- 
                  beta.star.samples[, (i - 1) * n.abund.re + X.re.ind[j, t]] * X.random[j, t] + 
                  beta.star.sites.0.samples[, i, j]
              } else {
                beta.star.sites.0.samples[, i, j] <- 
                  rnorm(n.post, 0, sqrt(object$sigma.sq.mu.samples[, t])) * X.random[j, t] + 
                  beta.star.sites.0.samples[, i, j]
              }
            } # j
          } # t
        } # i
      } else {
        X.fix <- X.0
        beta.star.sites.0.samples <- array(0, dim = c(n.post, n.sp, nrow(X.0)))
        p.abund.re <- 0
      }

      # Sub-sample previous
      theta.samples <- t(theta.samples)
      lambda.samples <- t(lambda.samples)
      beta.samples <- t(beta.samples)
      if (family == 'NB') {
        kappa.samples <- t(kappa.samples)
      } else {
        kappa.samples <- matrix(0, n.sp, n.post)
      }
      w.samples <- aperm(w.samples, c(2, 3, 1))
      beta.star.sites.0.samples <- aperm(beta.star.sites.0.samples, c(2, 3, 1))

      n.obs.0 <- nrow(X.fix)
      J.0 <- length(unique(sites.0.indx))

      match.indx <- match(do.call("paste", as.data.frame(coords.0)), do.call("paste", as.data.frame(coords)))
      coords.0.indx <- which(is.na(match.indx))
      coords.indx <- match.indx[!is.na(match.indx)]
      coords.place.indx <- which(!is.na(match.indx))
      # Indicates whether a site has been sampled. 1 = sampled
      sites.0.sampled <- ifelse(!is.na(match.indx), 1, 0)
      sites.link <- rep(NA, J.0)
      sites.link[which(!is.na(match.indx))] <- coords.indx
      # For C
      sites.link <- sites.link - 1

      if (sp.type == 'GP') {
        # Not currently implemented or accessed.
      } else {
        # Get nearest neighbors
        # nn2 is a function from RANN.
        nn.indx.0 <- nn2(coords, coords.0, k=n.neighbors)$nn.idx-1

        storage.mode(coords) <- "double"
        storage.mode(n.sp) <- "integer"
        storage.mode(J) <- "integer"
        storage.mode(n.obs) <- 'integer'
        storage.mode(p.abund) <- "integer"
        storage.mode(n.neighbors) <- "integer"
        storage.mode(X.fix) <- "double"
        storage.mode(coords.0) <- "double"
        storage.mode(J.0) <- "integer"
        storage.mode(n.obs.0) <- "integer"
        storage.mode(q) <- "integer"
        storage.mode(sites.link) <- "integer"
        storage.mode(sites.0.sampled) <- "integer"
        storage.mode(sites.0.indx) <- "integer"
        storage.mode(beta.samples) <- "double"
        storage.mode(theta.samples) <- "double"
        storage.mode(lambda.samples) <- "double"
        storage.mode(kappa.samples) <- "double"
        storage.mode(beta.star.sites.0.samples) <- "double"
        storage.mode(w.samples) <- "double"
        storage.mode(n.post) <- "integer"
        storage.mode(cov.model.indx) <- "integer"
        storage.mode(nn.indx.0) <- "integer"
        storage.mode(n.omp.threads) <- "integer"
        storage.mode(verbose) <- "integer"
        storage.mode(n.report) <- "integer"
        storage.mode(family.c) <- "integer"

        out <- .Call("sfMsAbundNNGPPredict", coords, J, n.obs, family.c, 
          	         n.sp, q, p.abund, n.neighbors,
                     X.fix, coords.0, J.0, n.obs.0, sites.link, sites.0.sampled, 
            	 sites.0.indx, nn.indx.0, beta.samples,
                     theta.samples, kappa.samples, lambda.samples, w.samples,
            	 beta.star.sites.0.samples, n.post,
                     cov.model.indx, n.omp.threads, verbose, n.report)
      }
      out$y.0.samples <- array(out$y.0.samples, dim = c(n.sp, n.obs.0, n.post))
      out$y.0.samples <- aperm(out$y.0.samples, c(3, 1, 2))
      out$w.0.samples <- array(out$w.0.samples, dim = c(q, J.0, n.post))
      out$w.0.samples <- aperm(out$w.0.samples, c(3, 1, 2))
      out$mu.0.samples <- array(out$mu.0.samples, dim = c(n.sp, n.obs.0, n.post))
      out$mu.0.samples <- aperm(out$mu.0.samples, c(3, 1, 2))

      tmp <- array(NA, dim = c(n.post, n.sp, length(c(missing.indx, non.missing.indx))))
      tmp[, , non.missing.indx] <- out$y.0.samples
      out$y.0.samples <- array(tmp, dim = c(n.post, n.sp, J.0, K.max.0))
      tmp <- array(NA, dim = c(n.post, n.sp, length(c(missing.indx, non.missing.indx))))
      tmp[, , non.missing.indx] <- out$mu.0.samples
      out$mu.0.samples <- array(tmp, dim = c(n.post, n.sp, J.0, K.max.0))
      out$run.time <- proc.time() - ptm
    }
    out$call <- cl
    out$object.class <- class(object)
    class(out) <- "predict.sfMsAbund"
    out
  }
}

# NMix --------------------------------------------------------------------
print.NMix <- function(x, ...) {
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.75)), 
      "", sep = "\n")
}

fitted.NMix <- function(object, type = 'marginal', ...) {
  # Check for unused arguments ------------------------------------------
  formal.args <- names(formals(sys.function(sys.parent())))
  elip.args <- names(list(...))
  for(i in elip.args){
      if(! i %in% formal.args)
          warning("'",i, "' is not an argument")
  }
  # Call ----------------------------------------------------------------
  cl <- match.call()
  # Functions -------------------------------------------------------------
  logit <- function(theta, a = 0, b = 1) {log((theta-a)/(b-theta))}
  logit.inv <- function(z, a = 0, b = 1) {b-(b-a)/(1+exp(z))}

  # Some initial checks -------------------------------------------------
  # Object ----------------------------
  if (missing(object)) {
    stop("error: object must be specified")
  }
  if (!(class(object) %in% c("NMix", "spNMix"))) {
    stop("error: object must be one of class NMix or spNMix\n")
  }
  if (!(type %in% c('marginal', 'conditional'))) {
    stop("type must be either 'marginal', or 'conditional'")
  }
  n.post <- object$n.post * object$n.chains
  X.p <- object$X.p
  y <- object$y
  n.rep <- apply(y, 1, function(a) sum(!is.na(a)))
  K.max <- dim(y)[2]
  J <- nrow(y)
  N.long.indx <- rep(1:J, K.max)
  N.long.indx <- N.long.indx[!is.na(c(y))]
  if (object$pRE) {
    X.p.random <- object$X.p.random
    # Add 1 to get it to R indexing from the C indexing. 
    X.p.re <- object$X.p.re + 1
  }
  y <- c(y)
  y <- y[!is.na(y)]
  if (type == 'conditional') {
    N.samples <- object$N.samples
  }
  if (type == 'marginal') {
    N.samples <- array(NA, dim = dim(object$N.samples))
    for (i in 1:n.post) {
      if (object$dist == 'Poisson') {
        N.samples[i, ] <- rpois(J, object$mu.samples[i, ] * object$offset)
      } else if (object$dist == 'NB') {
        N.samples[i, ] <- rnbinom(J, object$kappa.samples[i], 
				  mu = object$mu.samples[i, ] * object$offset)
      }
    }
  }
  mu.samples <- object$mu.samples
  alpha.samples <- object$alpha.samples
  det.prob.samples <- matrix(NA, n.post, length(y))
  if (object$pRE) {
    for (i in 1:n.post) {
      alpha.star.sum <- apply(apply(X.p.re, 2, function(a) object$alpha.star.samples[i, a]) * X.p.random,
                              1, sum)
      det.prob.samples[i, ] <- logit.inv(X.p %*% alpha.samples[i, ] + 
					 alpha.star.sum)
    }
  } else {
    det.prob.samples <- t(logit.inv(X.p %*% t(alpha.samples)))
  }
  y.rep.samples <- array(NA, dim = dim(det.prob.samples))
  n.obs <- ncol(det.prob.samples)
  for (i in 1:n.post) {
    y.rep.samples[i, ] <- rbinom(n.obs, N.samples[i, N.long.indx], 
				det.prob.samples[i, ]) 
  }
  tmp <- array(NA, dim = c(J * K.max, n.post))
  names.long <- which(!is.na(c(object$y)))
  tmp[names.long, ] <- t(y.rep.samples)
  y.rep.samples <- array(tmp, dim = c(J, K.max, n.post))
  y.rep.samples <- aperm(y.rep.samples, c(3, 1, 2))
  tmp <- array(NA, dim = c(J * K.max, n.post))
  tmp[names.long, ] <- t(det.prob.samples)
  det.prob.samples <- array(tmp, dim = c(J, K.max, n.post))
  det.prob.samples <- aperm(det.prob.samples, c(3, 1, 2))
  out <- list()
  out$y.rep.samples <- y.rep.samples
  out$p.samples <- det.prob.samples
  return(out)
}

predict.NMix <- function(object, X.0, ignore.RE = FALSE, 
			  type = 'abundance', ...) {
  # Check for unused arguments ------------------------------------------
  formal.args <- names(formals(sys.function(sys.parent())))
  elip.args <- names(list(...))
  for(i in elip.args){
      if(! i %in% formal.args)
          warning("'",i, "' is not an argument")
  }
  # Call ----------------------------------------------------------------
  cl <- match.call()

  # Functions ---------------------------------------------------------------
  logit <- function(theta, a = 0, b = 1) {log((theta-a)/(b-theta))}
  logit.inv <- function(z, a = 0, b = 1) {b-(b-a)/(1+exp(z))}

  # Some initial checks ---------------------------------------------------
  if (missing(object)) {
    stop("error: predict expects object\n")
  }
  if (!(tolower(type) %in% c('abundance', 'detection'))) {
    stop("error: prediction type must be either 'abundance' or 'detection'")
  }

  # Check X.0 -------------------------------------------------------------
  if (missing(X.0)) {
    stop("error: X.0 must be specified\n")
  }
  if (!any(is.data.frame(X.0), is.matrix(X.0))) {
    stop("error: X.0 must be a data.frame or matrix\n")
  }
  # Abundance predictions ------------------------------------------------
  if (tolower(type) == 'abundance') {  
    p.abund <- ncol(object$X)
    # p.design <- p.abund
    # if (object$muRE & !ignore.RE) {
    #   p.design <- p.abund + ncol(object$sigma.sq.mu.samples)
    # }
    # if (ncol(X.0) != p.design) {
    #   stop(paste("error: X.0 must have ", p.design, " columns\n", sep = ''))
    # }

    # Composition sampling --------------------------------------------------
    beta.samples <- as.matrix(object$beta.samples)
    if (object$dist == 'NB') {
      kappa.samples <- as.matrix(object$kappa.samples)
    }
    n.post <- object$n.post * object$n.chains
    out <- list()
    if (object$muRE) {
      p.abund.re <- length(object$re.level.names)
    } else {
      p.abund.re <- 0
    }
    re.cols <- object$re.cols

    if (object$muRE & !ignore.RE) {
      beta.star.samples <- object$beta.star.samples
      re.level.names <- object$re.level.names
      # Get columns in design matrix with random effects
      x.re.names <- colnames(object$X.re)
      x.0.names <- colnames(X.0)
      re.long.indx <- sapply(re.cols, length)
      tmp <- sapply(x.re.names, function(a) which(colnames(X.0) %in% a))
      indx <- list()
      for (i in 1:length(tmp)) {
        indx[[i]] <- rep(tmp[i], re.long.indx[i])
      }
      indx <- unlist(indx)
      if (length(indx) == 0) {
        stop("error: column names in X.0 must match variable names in data$occ.covs")
      }
      n.abund.re <- length(indx)
      n.unique.abund.re <- length(unique(indx))
      # Check RE columns
      for (i in 1:n.abund.re) {
        if (is.character(re.cols[[i]])) {
          # Check if all column names in svc are in occ.covs
          if (!all(re.cols[[i]] %in% x.0.names)) {
              missing.cols <- re.cols[[i]][!(re.cols[[i]] %in% x.0.names)]
              stop(paste("error: variable name ", paste(missing.cols, collapse=" and "), " not in abundance covariates", sep=""))
          }
          # Convert desired column names into the numeric column index
          re.cols[[i]] <- which(x.0.names %in% re.cols[[i]])
          
        } else if (is.numeric(re.cols[[i]])) {
          # Check if all column indices are in 1:p.abund
          if (!all(re.cols %in% 1:p.abund)) {
              missing.cols <- re.cols[[i]][!(re.cols[[i]] %in% (1:p.abund))]
              stop(paste("error: column index ", paste(missing.cols, collapse=" "), " not in design matrix columns", sep=""))
          }
        }
      }
      re.cols <- unlist(re.cols)
      X.re <- as.matrix(X.0[, indx, drop = FALSE])
      X.fix <- as.matrix(X.0[, -indx, drop = FALSE])
      X.random <- as.matrix(X.0[, re.cols, drop = FALSE])
      n.abund.re <- length(unlist(re.level.names))
      X.re.ind <- matrix(NA, nrow(X.re), p.abund.re)
      for (i in 1:p.abund.re) {
        for (j in 1:nrow(X.re)) {
          tmp <- which(re.level.names[[i]] == X.re[j, i])
          if (length(tmp) > 0) {
            X.re.ind[j, i] <- tmp 
          }
        }
      }
      if (p.abund.re > 1) {
        for (j in 2:p.abund.re) {
          X.re.ind[, j] <- X.re.ind[, j] + max(X.re.ind[, j - 1]) 
        }
      }
      # Create the random effects corresponding to each 
      # new location
      # ORDER: ordered by site, then species within site.
      beta.star.sites.0.samples <- matrix(0, n.post,  nrow(X.re))
      for (t in 1:p.abund.re) {
        for (j in 1:nrow(X.re)) {
          if (!is.na(X.re.ind[j, t])) {
            beta.star.sites.0.samples[, j] <- 
              beta.star.samples[, X.re.ind[j, t]] * X.random[j, t] + 
              beta.star.sites.0.samples[, j]
          } else {
            beta.star.sites.0.samples[, j] <- 
              rnorm(n.post, 0, sqrt(object$sigma.sq.mu.samples[, t])) * X.random[j, t] + 
              beta.star.sites.0.samples[, j]
          }
        } # j
      } # t
    } else {
      X.fix <- X.0
      beta.star.sites.0.samples <- matrix(0, n.post, nrow(X.0))
      p.abund.re <- 0
    }

    J.str <- nrow(X.0)
    out$mu.0.samples <- mcmc(exp(t(X.fix %*% t(beta.samples) + 
          				t(beta.star.sites.0.samples))))
    if (object$dist == 'NB') {
      out$N.0.samples <- mcmc(matrix(rnbinom(length(out$mu.0.samples), kappa.samples, 
					   mu = c(out$mu.0.samples)), nrow = n.post, 
				   ncol = nrow(X.0)))
    } else {
      out$N.0.samples <- mcmc(matrix(rpois(length(out$mu.0.samples),  
					   lambda = c(out$mu.0.samples)), nrow = n.post, 
				   ncol = nrow(X.0)))
    }
  }
  # Detection predictions -------------------------------------------------
  if (tolower(type) == 'detection') {
    p.det <- ncol(object$X.p)
    # p.design <- p.det
    # if (object$pRE & !ignore.RE) {
    #   p.design <- p.det + ncol(object$sigma.sq.p.samples)
    # }
    # if (ncol(X.0) != p.design) {
    #   stop(paste("error: X.0 must have ", p.design, " columns\n", sep = ''))
    # }
    re.det.cols <- object$re.det.cols

    # Composition sampling --------------------------------------------------
    alpha.samples <- as.matrix(object$alpha.samples)
    n.post <- object$n.post * object$n.chains
    out <- list()
    if (object$pRE) {
      p.det.re <- length(object$p.re.level.names)
    } else {
      p.det.re <- 0
    }

    if (object$pRE & !ignore.RE) {
      alpha.star.samples <- object$alpha.star.samples
      p.re.level.names <- object$p.re.level.names
      # Get columns in design matrix with random effects
      x.p.re.names <- colnames(object$X.p.re)
      x.p.0.names <- colnames(X.0)
      re.long.indx <- sapply(re.det.cols, length)
      tmp <- sapply(x.p.re.names, function(a) which(colnames(X.0) %in% a))
      indx <- list()
      for (i in 1:length(tmp)) {
        indx[[i]] <- rep(tmp[i], re.long.indx[i])
      }
      indx <- unlist(indx)
      if (length(indx) == 0) {
        stop("error: column names in X.0 must match variable names in data$det.covs")
      }
      n.det.re <- length(indx)
      n.unique.det.re <- length(unique(indx))
      # Check RE columns
      for (i in 1:n.det.re) {
        if (is.character(re.det.cols[[i]])) {
          # Check if all column names in svc are in occ.covs
          if (!all(re.det.cols[[i]] %in% x.p.0.names)) {
              missing.cols <- re.det.cols[[i]][!(re.det.cols[[i]] %in% x.p.0.names)]
              stop(paste("error: variable name ", paste(missing.cols, collapse=" and "), " not in detection covariates", sep=""))
          }
          # Convert desired column names into the numeric column index
          re.det.cols[[i]] <- which(x.p.0.names %in% re.det.cols[[i]])
          
        } else if (is.numeric(re.det.cols[[i]])) {
          # Check if all column indices are in 1:p.abund
          if (!all(re.det.cols %in% 1:p.det)) {
              missing.cols <- re.det.cols[[i]][!(re.det.cols[[i]] %in% (1:p.det))]
              stop(paste("error: column index ", paste(missing.cols, collapse=" "), " not in design matrix columns", sep=""))
          }
        }
      }
      re.det.cols <- unlist(re.det.cols)
      X.re <- as.matrix(X.0[, indx, drop = FALSE])
      X.fix <- as.matrix(X.0[, -indx, drop = FALSE])
      X.random <- as.matrix(X.0[, re.det.cols, drop = FALSE])
      n.det.re <- length(unlist(p.re.level.names))
      X.re.ind <- matrix(NA, nrow(X.re), p.det.re)
      for (i in 1:p.det.re) {
        for (j in 1:nrow(X.re)) {
          tmp <- which(p.re.level.names[[i]] == X.re[j, i])
          if (length(tmp) > 0) {
            X.re.ind[j, i] <- tmp 
          }
        }
      }
      if (p.det.re > 1) {
        for (j in 2:p.det.re) {
          X.re.ind[, j] <- X.re.ind[, j] + max(X.re.ind[, j - 1]) 
        }
      }
      # Create the random effects corresponding to each 
      # new location
      # ORDER: ordered by site, then species within site.
      alpha.star.sites.0.samples <- matrix(0, n.post,  nrow(X.re))
      for (t in 1:p.det.re) {
        for (j in 1:nrow(X.re)) {
          if (!is.na(X.re.ind[j, t])) {
            alpha.star.sites.0.samples[, j] <- 
              alpha.star.samples[, X.re.ind[j, t]] * X.random[j, t] + 
              alpha.star.sites.0.samples[, j]
          } else {
            alpha.star.sites.0.samples[, j] <- 
              rnorm(n.post, 0, sqrt(object$sigma.sq.p.samples[, t])) * X.random[j, t] + 
              alpha.star.sites.0.samples[, j]
          }
        } # j
      } # t
    } else {
      X.fix <- X.0
      alpha.star.sites.0.samples <- matrix(0, n.post, nrow(X.0))
      p.det.re <- 0
    }
    out$p.0.samples <- mcmc(logit.inv(t(X.fix %*% t(alpha.samples) + 
          				t(alpha.star.sites.0.samples))))
  }
  out$call <- cl

  class(out) <- "predict.NMix"
  out
}


summary.NMix <- function(object,
			  quantiles = c(0.025, 0.5, 0.975),
			  digits = max(3L, getOption("digits") - 3L), ...) {
  print(object)

  n.post <- object$n.post
  n.samples <- object$n.samples
  n.burn <- object$n.burn
  n.thin <- object$n.thin
  n.chains <- object$n.chains
  run.time <- object$run.time[3] / 60 # minutes

  cat(paste("Samples per Chain: ", n.samples,"\n", sep=""))
  cat(paste("Burn-in: ", n.burn,"\n", sep=""))
  cat(paste("Thinning Rate: ",n.thin,"\n", sep=""))
  cat(paste("Number of Chains: ", n.chains, "\n", sep = ""))
  cat(paste("Total Posterior Samples: ",n.post * n.chains,"\n", sep=""))
  cat(paste("Run Time (min): ", round(run.time, digits), "\n\n", sep = ""))

  # Abundance
  cat("Abundance (log scale): \n")
  tmp.1 <- t(apply(object$beta.samples, 2,
		   function(x) c(mean(x), sd(x))))
  colnames(tmp.1) <- c("Mean", "SD")
  tmp <- t(apply(object$beta.samples, 2,
		 function(x) quantile(x, prob = quantiles)))
  diags <- matrix(c(object$rhat$beta, round(object$ESS$beta, 0)), ncol = 2)
  colnames(diags) <- c('Rhat', 'ESS')

  print(noquote(round(cbind(tmp.1, tmp, diags), digits)))
  if (object$muRE) {
    cat("\n")
    cat("Abundance Random Effect Variances (log scale): \n")
    tmp.1 <- t(apply(object$sigma.sq.mu.samples, 2,
          	   function(x) c(mean(x), sd(x))))
    colnames(tmp.1) <- c("Mean", "SD")
    tmp <- t(apply(object$sigma.sq.mu.samples, 2,
          	 function(x) quantile(x, prob = quantiles)))
    diags <- matrix(c(object$rhat$sigma.sq.mu, round(object$ESS$sigma.sq.mu, 0)), ncol = 2)
    colnames(diags) <- c('Rhat', 'ESS')

    print(noquote(round(cbind(tmp.1, tmp, diags), digits)))
  }
  cat("\n")
  # Detection
  if (is(object, 'NMix')) {
    cat("Detection (logit scale): \n")
  } else if (is(object, 'DS')) {
    cat("Detection (log scale): \n")
  }
  tmp.1 <- t(apply(object$alpha.samples, 2,
		   function(x) c(mean(x), sd(x))))
  colnames(tmp.1) <- c("Mean", "SD")
  tmp <- t(apply(object$alpha.samples, 2,
		 function(x) quantile(x, prob = quantiles)))
  diags <- matrix(c(object$rhat$alpha, round(object$ESS$alpha, 0)), ncol = 2)
  colnames(diags) <- c('Rhat', 'ESS')
  print(noquote(round(cbind(tmp.1, tmp, diags), digits)))
  if (object$pRE) {
    cat("\n")
    if (is(object, 'NMix')) {
      cat("Detection Random Effect Variances (logit scale): \n")
    } else if (is(object, 'DS')) {
      cat("Detection Random Effect Variances (log scale): \n")
    }
    tmp.1 <- t(apply(object$sigma.sq.p.samples, 2,
          	   function(x) c(mean(x), sd(x))))
    colnames(tmp.1) <- c("Mean", "SD")
    tmp <- t(apply(object$sigma.sq.p.samples, 2,
          	 function(x) quantile(x, prob = quantiles)))
    diags <- matrix(c(object$rhat$sigma.sq.p, round(object$ESS$sigma.sq.p, 0)), ncol = 2)
    colnames(diags) <- c('Rhat', 'ESS')

    print(noquote(round(cbind(tmp.1, tmp, diags), digits)))
  }
  if (object$dist == "NB") {
    cat("\n")
    cat("NB overdispersion: \n")
    tmp.1 <- t(apply(object$kappa.samples, 2,
          	   function(x) c(mean(x), sd(x))))
    colnames(tmp.1) <- c("Mean", "SD")
    tmp <- t(apply(object$kappa.samples, 2,
          	 function(x) quantile(x, prob = quantiles)))
    diags <- matrix(c(object$rhat$kappa, round(object$ESS$kappa, 0)), ncol = 2)
    colnames(diags) <- c('Rhat', 'ESS')

    print(noquote(round(cbind(tmp.1, tmp, diags), digits)))

  }
}

# spNMix ------------------------------------------------------------------
print.spNMix <- function(x, ...) {
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.75)),
      "", sep = "\n")
}

summary.spNMix <- function(object,
			   quantiles = c(0.025, 0.5, 0.975),
			   digits = max(3L, getOption("digits") - 3L), ...) {
  print(object)

  n.post <- object$n.post
  n.samples <- object$n.samples
  n.burn <- object$n.burn
  n.thin <- object$n.thin
  n.chains <- object$n.chains
  run.time <- object$run.time[3] / 60 # minutes

  cat(paste("Samples per Chain: ", n.samples,"\n", sep=""))
  cat(paste("Burn-in: ", n.burn,"\n", sep=""))
  cat(paste("Thinning Rate: ",n.thin,"\n", sep=""))
  cat(paste("Number of Chains: ", n.chains, "\n", sep = ""))
  cat(paste("Total Posterior Samples: ",n.post * n.chains,"\n", sep=""))
  cat(paste("Run Time (min): ", round(run.time, digits), "\n\n", sep = ""))

  # Abundance
  cat("Abundance (log scale): \n")
  tmp.1 <- t(apply(object$beta.samples, 2,
		   function(x) c(mean(x), sd(x))))
  colnames(tmp.1) <- c("Mean", "SD")
  tmp <- t(apply(object$beta.samples, 2,
		 function(x) quantile(x, prob = quantiles)))
  diags <- matrix(c(object$rhat$beta, round(object$ESS$beta, 0)), ncol = 2)
  colnames(diags) <- c('Rhat', 'ESS')

  print(noquote(round(cbind(tmp.1, tmp, diags), digits)))
  if (object$muRE) {
    cat("\n")
    cat("Abundance Random Effect Variances (log scale): \n")
    tmp.1 <- t(apply(object$sigma.sq.mu.samples, 2,
          	   function(x) c(mean(x), sd(x))))
    colnames(tmp.1) <- c("Mean", "SD")
    tmp <- t(apply(object$sigma.sq.mu.samples, 2,
          	 function(x) quantile(x, prob = quantiles)))
    diags <- matrix(c(object$rhat$sigma.sq.mu, round(object$ESS$sigma.sq.mu, 0)), ncol = 2)
    colnames(diags) <- c('Rhat', 'ESS')

    print(noquote(round(cbind(tmp.1, tmp, diags), digits)))
  }
  cat("\n")
  # Detection
  if (is(object, 'spNMix')) {
    cat("Detection (logit scale): \n")
  } else if (is(object, 'spDS')) {
    cat("Detection (log scale): \n")
  }
  tmp.1 <- t(apply(object$alpha.samples, 2,
		   function(x) c(mean(x), sd(x))))
  colnames(tmp.1) <- c("Mean", "SD")
  tmp <- t(apply(object$alpha.samples, 2,
		 function(x) quantile(x, prob = quantiles)))
  diags <- matrix(c(object$rhat$alpha, round(object$ESS$alpha, 0)), ncol = 2)
  colnames(diags) <- c('Rhat', 'ESS')
  print(noquote(round(cbind(tmp.1, tmp, diags), digits)))
  if (object$pRE) {
    cat("\n")
    if (is(object, 'spNMix')) {
      cat("Detection Random Effect Variances (logit scale): \n")
    } else if (is(object, 'spDS')) {
      cat("Detection Random Effect Variances (log scale): \n")
    }
    tmp.1 <- t(apply(object$sigma.sq.p.samples, 2,
          	   function(x) c(mean(x), sd(x))))
    colnames(tmp.1) <- c("Mean", "SD")
    tmp <- t(apply(object$sigma.sq.p.samples, 2,
          	 function(x) quantile(x, prob = quantiles)))
    diags <- matrix(c(object$rhat$sigma.sq.p, round(object$ESS$sigma.sq.p, 0)), ncol = 2)
    colnames(diags) <- c('Rhat', 'ESS')

    print(noquote(round(cbind(tmp.1, tmp, diags), digits)))
  }
  cat("\n")
  # Covariance ------------------------
  cat("Spatial Covariance: \n")
  tmp.1 <- t(apply(object$theta.samples, 2, 
		   function(x) c(mean(x), sd(x))))
  colnames(tmp.1) <- c("Mean", "SD")
  tmp <- t(apply(object$theta.samples, 2, 
		 function(x) quantile(x, prob = quantiles)))
  diags <- matrix(c(object$rhat$theta, round(object$ESS$theta, 0)), ncol = 2)
  colnames(diags) <- c('Rhat', 'ESS')
  print(noquote(round(cbind(tmp.1, tmp, diags), digits)))

  if (object$dist == "NB") {
    cat("\n")
    cat("NB overdispersion: \n")
    tmp.1 <- t(apply(object$kappa.samples, 2,
          	   function(x) c(mean(x), sd(x))))
    colnames(tmp.1) <- c("Mean", "SD")
    tmp <- t(apply(object$kappa.samples, 2,
          	 function(x) quantile(x, prob = quantiles)))
    diags <- matrix(c(object$rhat$kappa, round(object$ESS$kappa, 0)), ncol = 2)
    colnames(diags) <- c('Rhat', 'ESS')

    print(noquote(round(cbind(tmp.1, tmp, diags), digits)))
  }
}

fitted.spNMix <- function(object, type = 'marginal', ...) {
  fitted.NMix(object, type = type)
}

predict.spNMix <- function(object, X.0, coords.0, n.omp.threads = 1, 
			   verbose = TRUE, n.report = 100, 
			   ignore.RE = FALSE, 
			   type = 'abundance', include.sp = TRUE, ...) {
  # Check for unused arguments ------------------------------------------
  formal.args <- names(formals(sys.function(sys.parent())))
  elip.args <- names(list(...))
  for(i in elip.args){
      if(! i %in% formal.args)
          warning("'",i, "' is not an argument")
  }
  # Call ----------------------------------------------------------------
  cl <- match.call()

  # Call predict.NMix if don't care about spatial effects
  if (!include.sp) {
    out <- predict.NMix(object, X.0, ignore.RE, type = type)
  } else {

    # Functions ---------------------------------------------------------------
    logit <- function(theta, a = 0, b = 1) {log((theta-a)/(b-theta))}
    logit.inv <- function(z, a = 0, b = 1) {b-(b-a)/(1+exp(z))}

    # Some initial checks ---------------------------------------------------
    if (missing(object)) {
      stop("error: predict expects object\n")
    }
    if (!(tolower(type) %in% c('abundance', 'detection'))) {
      stop("error: prediction type must be either 'abundance' or 'detection'")
    }

    # Check X.0 -------------------------------------------------------------
    if (missing(X.0)) {
      stop("error: X.0 must be specified\n")
    }
    if (!any(is.data.frame(X.0), is.matrix(X.0))) {
      stop("error: X.0 must be a data.frame or matrix\n")
    }

    # Abundance predictions ------------------------------------------------
    if (tolower(type) == 'abundance') {  
      if (missing(coords.0)) {
        stop("error: coords.0 must be specified\n")
      }
      if (!any(is.data.frame(coords.0), is.matrix(coords.0))) {
        stop("error: coords.0 must be a data.frame or matrix\n")
      }
      if (!ncol(coords.0) == 2){
        stop("error: coords.0 must have two columns\n")
      }
      coords.0 <- as.matrix(coords.0)
      
      p.abund <- ncol(object$X)
      p.design <- p.abund
      # if (object$muRE & !ignore.RE) {
      #   p.design <- p.abund + ncol(object$sigma.sq.mu.samples)
      # }
      # if (ncol(X.0) != p.design) {
      #   stop(paste("error: X.0 must have ", p.design, " columns\n", sep = ''))
      # }
      X <- object$X
      coords <- object$coords
      J <- nrow(X)
      beta.samples <- as.matrix(object$beta.samples)
      n.post <- object$n.post * object$n.chains
      family <- object$dist
      family.c <- ifelse(family == 'NB', 1, 0)
      theta.samples <- object$theta.samples
      w.samples <- object$w.samples
      n.neighbors <- object$n.neighbors
      cov.model.indx <- object$cov.model.indx
      re.cols <- object$re.cols
      sp.type <- object$type
      out <- list()
      if (object$muRE) {
        p.abund.re <- length(object$re.level.names)
      } else {
        p.abund.re <- 0
      }

      # Eliminate prediction sites that have already sampled been for now
      match.indx <- match(do.call("paste", as.data.frame(coords.0)), do.call("paste", as.data.frame(coords)))
      coords.0.indx <- which(is.na(match.indx))
      coords.indx <- match.indx[!is.na(match.indx)]
      coords.place.indx <- which(!is.na(match.indx))
      # coords.0.new <- coords.0[coords.0.indx, , drop = FALSE]
      # X.0.new <- X.0[coords.0.indx, , drop = FALSE]

      if (object$muRE & !ignore.RE) {
        beta.star.samples <- object$beta.star.samples
        re.level.names <- object$re.level.names
        # Get columns in design matrix with random effects
        x.re.names <- colnames(object$X.re)
        x.0.names <- colnames(X.0)
        re.long.indx <- sapply(re.cols, length)
        tmp <- sapply(x.re.names, function(a) which(colnames(X.0) %in% a))
        indx <- list()
        for (i in 1:length(tmp)) {
          indx[[i]] <- rep(tmp[i], re.long.indx[i])
        }
        indx <- unlist(indx)
        if (length(indx) == 0) {
          stop("error: column names in X.0 must match variable names in data$abund.covs")
        }
        n.abund.re <- length(indx)
        n.unique.abund.re <- length(unique(indx))
        # Check RE columns
        for (i in 1:n.abund.re) {
          if (is.character(re.cols[[i]])) {
            # Check if all column names in svc are in occ.covs
            if (!all(re.cols[[i]] %in% x.0.names)) {
                missing.cols <- re.cols[[i]][!(re.cols[[i]] %in% x.0.names)]
                stop(paste("error: variable name ", paste(missing.cols, collapse=" and "), " not in abundance covariates", sep=""))
            }
            # Convert desired column names into the numeric column index
            re.cols[[i]] <- which(x.0.names %in% re.cols[[i]])
            
          } else if (is.numeric(re.cols[[i]])) {
            # Check if all column indices are in 1:p.abund
            if (!all(re.cols %in% 1:p.abund)) {
                missing.cols <- re.cols[[i]][!(re.cols[[i]] %in% (1:p.abund))]
                stop(paste("error: column index ", paste(missing.cols, collapse=" "), " not in design matrix columns", sep=""))
            }
          }
        }
        re.cols <- unlist(re.cols)
        X.re <- as.matrix(X.0[, indx, drop = FALSE])
        X.fix <- as.matrix(X.0[, -indx, drop = FALSE])
        X.random <- as.matrix(X.0[, re.cols, drop = FALSE])
        n.abund.re <- length(unlist(re.level.names))
        X.re.ind <- matrix(NA, nrow(X.re), p.abund.re)
        for (i in 1:p.abund.re) {
          for (j in 1:nrow(X.re)) {
            tmp <- which(re.level.names[[i]] == X.re[j, i])
            if (length(tmp) > 0) {
              X.re.ind[j, i] <- tmp 
            }
          }
        }
        if (p.abund.re > 1) {
          for (j in 2:p.abund.re) {
            X.re.ind[, j] <- X.re.ind[, j] + max(X.re.ind[, j - 1]) 
          }
        }
        # Create the random effects corresponding to each 
        # new location
        # ORDER: ordered by site, then species within site.
        beta.star.sites.0.samples <- matrix(0, n.post,  nrow(X.re))
        for (t in 1:p.abund.re) {
          for (j in 1:nrow(X.re)) {
            if (!is.na(X.re.ind[j, t])) {
              beta.star.sites.0.samples[, j] <- 
                beta.star.samples[, X.re.ind[j, t]] * X.random[j, t] + 
                beta.star.sites.0.samples[, j]
            } else {
              beta.star.sites.0.samples[, j] <- 
                rnorm(n.post, 0, sqrt(object$sigma.sq.mu.samples[, t])) * X.random[j, t] + 
                beta.star.sites.0.samples[, j]
            }
          } # j
        } # t
      } else {
        X.fix <- X.0
        beta.star.sites.0.samples <- matrix(0, n.post, nrow(X.0))
        p.abund.re <- 0
      }
      # Get samples in proper format for C++
      beta.samples <- t(beta.samples)
      w.samples <- t(w.samples)
      beta.star.sites.0.samples <- t(beta.star.sites.0.samples)
      if (family == 'NB') {
        kappa.samples <- t(object$kappa.samples)
      } else {
        kappa.samples <- matrix(0, 1, n.post)
      }
      theta.samples <- t(theta.samples)

      sites.0.indx <- 0:(nrow(X.0) - 1)
      J.0 <- length(unique(sites.0.indx))
      sites.0.sampled <- ifelse(!is.na(match.indx), 1, 0)
      sites.link <- rep(NA, J.0)
      sites.link[which(!is.na(match.indx))] <- coords.indx
      # For C
      sites.link <- sites.link - 1

      if (sp.type == 'GP') {
        stop("NNGP = FALSE is not currently supported for spAbund")
      } else {
        # Get nearest neighbors
        # nn2 is a function from RANN.
        nn.indx.0 <- nn2(coords, coords.0, k=n.neighbors)$nn.idx-1

        storage.mode(coords) <- "double"
        storage.mode(J) <- "integer"
        storage.mode(p.abund) <- "integer"
        storage.mode(n.neighbors) <- "integer"
        storage.mode(X.fix) <- "double"
        storage.mode(coords.0) <- "double"
        storage.mode(J.0) <- "integer"
        storage.mode(sites.link) <- "integer"
        storage.mode(sites.0.sampled) <- 'integer'
        storage.mode(sites.0.indx) <- 'integer'
        storage.mode(beta.samples) <- "double"
        storage.mode(theta.samples) <- "double"
        storage.mode(kappa.samples) <- "double"
        storage.mode(family.c) <- "integer"
        storage.mode(w.samples) <- "double"
        storage.mode(beta.star.sites.0.samples) <- "double"
        storage.mode(n.post) <- "integer"
        storage.mode(cov.model.indx) <- "integer"
        storage.mode(nn.indx.0) <- "integer"
        storage.mode(n.omp.threads) <- "integer"
        storage.mode(verbose) <- "integer"
        storage.mode(n.report) <- "integer"

        ptm <- proc.time()

        out <- .Call("spNMixNNGPPredict", coords, J, family.c, p.abund, n.neighbors,
                     X.fix, coords.0, J.0, nn.indx.0, beta.samples,
                     theta.samples, kappa.samples, w.samples, beta.star.sites.0.samples,
          	   sites.link, sites.0.sampled, sites.0.indx,
            	   n.post, cov.model.indx, n.omp.threads, verbose, n.report)
      }
      out$N.0.samples <- mcmc(t(out$N.0.samples))
      out$mu.0.samples <- mcmc(t(out$mu.0.samples))
      out$w.0.samples <- mcmc(t(out$w.0.samples))
    }
    if (tolower(type) == 'detection') {
      out <- predict.NMix(object, X.0, ignore.RE, type, ...)
    }
  }
  out$call <- cl

  class(out) <- "predict.spNMix"
  out

}

# msNMix ------------------------------------------------------------------
print.msNMix <- function(x, ...) {
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.75)),
      "", sep = "\n")
}

summary.msNMix <- function(object,
			   level = 'both',
			   quantiles = c(0.025, 0.5, 0.975),
			   digits = max(3L, getOption("digits") - 3L), ...) {

  print(object)

  n.post <- object$n.post
  n.samples <- object$n.samples
  n.burn <- object$n.burn
  n.thin <- object$n.thin
  n.chains <- object$n.chains
  run.time <- object$run.time[3] / 60 # minutes

  cat(paste("Samples per Chain: ", n.samples,"\n", sep=""))
  cat(paste("Burn-in: ", n.burn,"\n", sep=""))
  cat(paste("Thinning Rate: ",n.thin,"\n", sep=""))
  cat(paste("Number of Chains: ", n.chains, "\n", sep = ""))
  cat(paste("Total Posterior Samples: ",n.post * n.chains,"\n", sep=""))
  cat(paste("Run Time (min): ", round(run.time, digits), "\n\n", sep = ""))

  if (tolower(level) %in% c('community', 'both')) {

    cat("----------------------------------------\n");
    cat("\tCommunity Level\n");
    cat("----------------------------------------\n");

    # Abundance 
    cat("Abundance Means (log scale): \n")
    tmp.1 <- t(apply(object$beta.comm.samples, 2,
          	   function(x) c(mean(x), sd(x))))
    colnames(tmp.1) <- c("Mean", "SD")
    tmp <- t(apply(object$beta.comm.samples, 2,
          	 function(x) quantile(x, prob = quantiles)))
    diags <- matrix(c(object$rhat$beta.comm, round(object$ESS$beta.comm, 0)), ncol = 2)
    colnames(diags) <- c('Rhat', 'ESS')

    print(noquote(round(cbind(tmp.1, tmp, diags), digits)))

    cat("\nAbundance Variances (log scale): \n")
    tmp.1 <- t(apply(object$tau.sq.beta.samples, 2,
          	   function(x) c(mean(x), sd(x))))
    colnames(tmp.1) <- c("Mean", "SD")
    tmp <- t(apply(object$tau.sq.beta.samples, 2,
          	 function(x) quantile(x, prob = quantiles)))
    diags <- matrix(c(object$rhat$tau.sq.beta, round(object$ESS$tau.sq.beta, 0)), ncol = 2)
    colnames(diags) <- c('Rhat', 'ESS')
    print(noquote(round(cbind(tmp.1, tmp, diags), digits)))
    if (object$muRE) {
      cat("\n")
      cat("Abundance Random Effect Variances (log scale): \n")
      tmp.1 <- t(apply(object$sigma.sq.mu.samples, 2,
            	   function(x) c(mean(x), sd(x))))
      colnames(tmp.1) <- c("Mean", "SD")
      tmp <- t(apply(object$sigma.sq.mu.samples, 2,
            	 function(x) quantile(x, prob = quantiles)))
      diags <- matrix(c(object$rhat$sigma.sq.mu, round(object$ESS$sigma.sq.mu, 0)), ncol = 2)
      colnames(diags) <- c('Rhat', 'ESS')

      print(noquote(round(cbind(tmp.1, tmp, diags), digits)))
    }
    cat("\n")
    # Detection
    if (class(object) %in% c('msNMix', 'lfMsNMix', 'sfMsNMix')) {
      cat("Detection Means (logit scale): \n")
    } else if (class(object) %in% c('msDS', 'lfMsDS', 'sfMsDS')) {
      cat("Detection Means (log scale): \n")
    }
    tmp.1 <- t(apply(object$alpha.comm.samples, 2,
          	   function(x) c(mean(x), sd(x))))
    colnames(tmp.1) <- c("Mean", "SD")
    tmp <- t(apply(object$alpha.comm.samples, 2,
          	 function(x) quantile(x, prob = quantiles)))
    diags <- matrix(c(object$rhat$alpha.comm, round(object$ESS$alpha.comm, 0)), ncol = 2)
    colnames(diags) <- c('Rhat', 'ESS')

    print(noquote(round(cbind(tmp.1, tmp, diags), digits)))
    if (class(object) %in% c('msNMix', 'lfMsNMix', 'sfMsNMix')) {
      cat("\nDetection Variances (logit scale): \n")
    } else if (class(object) %in% c('msDS', 'lfMsDS', 'sfMsDS')) {
      cat("\nDetection Variances (log scale): \n")
    }
    tmp.1 <- t(apply(object$tau.sq.alpha.samples, 2,
          	   function(x) c(mean(x), sd(x))))
    colnames(tmp.1) <- c("Mean", "SD")
    tmp <- t(apply(object$tau.sq.alpha.samples, 2,
          	 function(x) quantile(x, prob = quantiles)))
    diags <- matrix(c(object$rhat$tau.sq.alpha, round(object$ESS$tau.sq.alpha, 0)), ncol = 2)
    colnames(diags) <- c('Rhat', 'ESS')

    print(noquote(round(cbind(tmp.1, tmp, diags), digits)))
    if (object$pRE) {
      cat("\n")
    if (class(object) %in% c('msNMix', 'lfMsNMix', 'sfMsNMix')) {
        cat("Detection Random Effect Variances (logit scale): \n")
    } else if (class(object) %in% c('msDS', 'lfMsDS', 'sfMsDS')) {
        cat("Detection Random Effect Variances (log scale): \n")
      }
      tmp.1 <- t(apply(object$sigma.sq.p.samples, 2,
            	   function(x) c(mean(x), sd(x))))
      colnames(tmp.1) <- c("Mean", "SD")
      tmp <- t(apply(object$sigma.sq.p.samples, 2,
            	 function(x) quantile(x, prob = quantiles)))
      diags <- matrix(c(object$rhat$sigma.sq.p, round(object$ESS$sigma.sq.p, 0)), ncol = 2)
      colnames(diags) <- c('Rhat', 'ESS')

      print(noquote(round(cbind(tmp.1, tmp, diags), digits)))
    }
  }

  if (tolower(level) %in% c('species', 'both')) {
    if (tolower(level) == 'both') cat("\n")
    cat("----------------------------------------\n");
    cat("\tSpecies Level\n");
    cat("----------------------------------------\n");
    cat("Abundance (log scale): \n")
    tmp.1 <- t(apply(object$beta.samples, 2,
          	   function(x) c(mean(x), sd(x))))
    colnames(tmp.1) <- c("Mean", "SD")
    tmp <- t(apply(object$beta.samples, 2,
          	 function(x) quantile(x, prob = quantiles)))
    diags <- matrix(c(object$rhat$beta, round(object$ESS$beta, 0)), ncol = 2)
    colnames(diags) <- c('Rhat', 'ESS')

    print(noquote(round(cbind(tmp.1, tmp, diags), digits)))
    cat("\n")
    # Detection
    if (class(object) %in% c('msNMix', 'lfMsNMix', 'sfMsNMix')) {
      cat("Detection (logit scale): \n")
    } else if (class(object) %in% c('msDS', 'lfMsDS', 'sfMsDS')) {
      cat("Detection (log scale): \n")
    }
    tmp.1 <- t(apply(object$alpha.samples, 2,
          	   function(x) c(mean(x), sd(x))))
    colnames(tmp.1) <- c("Mean", "SD")
    tmp <- t(apply(object$alpha.samples, 2,
          	 function(x) quantile(x, prob = quantiles)))
    diags <- matrix(c(object$rhat$alpha, round(object$ESS$alpha, 0)), ncol = 2)
    colnames(diags) <- c('Rhat', 'ESS')
    print(noquote(round(cbind(tmp.1, tmp, diags), digits)))
    
    # NB Overdispersion
    if (object$dist == "NB") {
      cat("\n")
      cat("----------------------------------------\n");
      cat("\tNB overdispersion\n");
      cat("----------------------------------------\n");
      tmp.1 <- t(apply(object$kappa.samples, 2,
            	   function(x) c(mean(x), sd(x))))
      colnames(tmp.1) <- c("Mean", "SD")
      tmp <- t(apply(object$kappa.samples, 2,
            	 function(x) quantile(x, prob = quantiles)))
      diags <- matrix(c(object$rhat$kappa, round(object$ESS$kappa, 0)), ncol = 2)
      colnames(diags) <- c('Rhat', 'ESS')
      print(noquote(round(cbind(tmp.1, tmp, diags), digits)))
    }
  }
}


fitted.msNMix <- function(object, type = 'marginal', ...) {
  # Check for unused arguments ------------------------------------------
  formal.args <- names(formals(sys.function(sys.parent())))
  elip.args <- names(list(...))
  for(i in elip.args){
      if(! i %in% formal.args)
          warning("'",i, "' is not an argument")
  }
  # Call ----------------------------------------------------------------
  cl <- match.call()
  # Functions -------------------------------------------------------------
  logit <- function(theta, a = 0, b = 1) {log((theta-a)/(b-theta))}
  logit.inv <- function(z, a = 0, b = 1) {b-(b-a)/(1+exp(z))}

  # Some initial checks -------------------------------------------------
  # Object ----------------------------
  if (missing(object)) {
    stop("error: object must be specified")
  }
  if (!(class(object) %in% c('msNMix', 'spMsNMix', 'lfMsNMix', 'sfMsNMix'))) {
    stop("error: object must be of class msNMix, spMsNMix, lfMsNMix, or sfMsNMix\n")
  }
  if (!(type %in% c('marginal', 'conditional'))) {
    stop("type must be either 'marginal', or 'conditional'")
  }
  n.post <- object$n.post * object$n.chains
  X.p <- object$X.p
  y <- object$y
  n.rep <- apply(y[1, , , drop = FALSE], 2, function(a) sum(!is.na(a)))
  K.max <- dim(y)[3]
  J <- dim(y)[2]
  n.sp <- dim(y)[1]
  N.long.indx <- rep(1:J, dim(y)[3])
  N.long.indx <- N.long.indx[!is.na(c(y[1, , ]))]
  if (type == 'conditional') {
    N.samples <- object$N.samples
  }
  if (type == 'marginal') {
    N.samples <- array(NA, dim = dim(object$N.samples))
    for (j in 1:n.post) {
      for (i in 1:n.sp) {
        if (object$dist == 'Poisson') {
          N.samples[j, i, ] <- rpois(J, object$mu.samples[j, i, ] * object$offset)
        } else if (object$dist == 'NB') {
          N.samples[j, i, ] <- rnbinom(J, object$kappa.samples[j, i], 
          			  mu = object$mu.samples[j, i, ] * object$offset)
        }
      }
    }
  }
  alpha.samples <- object$alpha.samples
  n.obs <- nrow(X.p)
  det.prob.samples <- array(NA, dim = c(n.obs, n.sp, n.post))
  sp.indx <- rep(1:n.sp, ncol(X.p))
  y <- matrix(y, n.sp, J * K.max)
  y <- y[, apply(y, 2, function(a) !sum(is.na(a)) > 0)]
  for (i in 1:n.sp) {
    if (object$pRE) {
      sp.re.indx <- rep(1:n.sp, each = ncol(object$alpha.star.samples) / n.sp)
      # Add 1 to get it to R indexing. 
      X.p.re <- object$X.p.re + 1
      X.p.random <- object$X.p.random
      # tmp.samples <- matrix(0, n.post, n.obs)
      tmp.alpha.star <- object$alpha.star.samples[, sp.re.indx == i]
      for (j in 1:n.post) {
        alpha.star.sum <- apply(apply(X.p.re, 2, function(a) tmp.alpha.star[j, a]) * X.p.random, 
                               1, sum) 
        det.prob.samples[, i, j] <- logit.inv(X.p %*% alpha.samples[j, sp.indx == i] + alpha.star.sum)
      }
    } else {
      det.prob.samples[, i, ] <- logit.inv(X.p %*% t(alpha.samples[, sp.indx == i]))
    }
  }

  out <- list()
  # Get detection probability
  # Need to be careful here that all arrays line up. 
  det.prob.samples <- aperm(det.prob.samples, c(3, 2, 1))
  tmp <- array(NA, dim = c(n.post, n.sp, J * K.max))
  names.long <- which(!is.na(c(object$y[1, , ])))
  tmp[, , names.long] <- det.prob.samples
  p.samples <- array(tmp, dim = c(n.post, n.sp, J, K.max))
  out$p.samples <- p.samples
  # Get fitted values
  y.rep.samples <- array(NA, dim = dim(det.prob.samples))
  for (i in 1:n.sp) {
    y.rep.samples[, i, ] <- rbinom(n.obs * n.post, N.samples[, i, N.long.indx], 
				   det.prob.samples[, i, ])
  }
  tmp <- array(NA, dim = c(n.post, n.sp, J * K.max))
  names.long <- which(!is.na(c(object$y[1, , ])))
  tmp[, , names.long] <- y.rep.samples
  y.rep.samples <- array(tmp, dim = c(n.post, n.sp, J, K.max))
  out$y.rep.samples <- y.rep.samples
  return(out)
}

predict.msNMix <- function(object, X.0, ignore.RE = FALSE, 
			  type = 'abundance', ...) {
  # Check for unused arguments ------------------------------------------
  formal.args <- names(formals(sys.function(sys.parent())))
  elip.args <- names(list(...))
  for(i in elip.args){
      if(! i %in% formal.args)
          warning("'",i, "' is not an argument")
  }
  # Call ----------------------------------------------------------------
  cl <- match.call()

  # Functions ---------------------------------------------------------------
  logit <- function(theta, a = 0, b = 1) {log((theta-a)/(b-theta))}
  logit.inv <- function(z, a = 0, b = 1) {b-(b-a)/(1+exp(z))}

  # Some initial checks ---------------------------------------------------
  if (missing(object)) {
    stop("error: predict expects object\n")
  }

  if (!(tolower(type) %in% c('abundance', 'detection'))) {
    stop("error: prediction type must be either 'abundance' or 'detection'")
  }

  # Check X.0 -------------------------------------------------------------
  if (missing(X.0)) {
    stop("error: X.0 must be specified\n")
  }
  if (!any(is.data.frame(X.0), is.matrix(X.0))) {
    stop("error: X.0 must be a data.frame or matrix\n")
  }

  # Abundance predictions ------------------------------------------------
  if (tolower(type) == 'abundance') {
    p.abund <- ncol(object$X)
    # Composition sampling --------------------------------------------------
    beta.samples <- as.matrix(object$beta.samples)
    if (object$dist == 'NB') {
      kappa.samples <- as.matrix(object$kappa.samples)
    }
    n.sp <- nrow(object$y)
    sp.indx <- rep(1:n.sp, p.abund)
    n.post <- object$n.post * object$n.chains
    out <- list()
    out$mu.0.samples <- array(NA, dim = c(n.post, n.sp, nrow(X.0)))
    out$N.0.samples <- array(NA, dim = c(n.post, n.sp, nrow(X.0)))
    if (object$muRE) {
      p.abund.re <- length(object$re.level.names)
    } else {
      p.abund.re <- 0
    }
    re.cols <- object$re.cols

    if (object$muRE & !ignore.RE) {
      beta.star.samples <- object$beta.star.samples
      re.level.names <- object$re.level.names
      # Get columns in design matrix with random effects
      x.re.names <- colnames(object$X.re)
      x.0.names <- colnames(X.0)
      re.long.indx <- sapply(re.cols, length)
      tmp <- sapply(x.re.names, function(a) which(colnames(X.0) %in% a))
      indx <- list()
      for (i in 1:length(tmp)) {
        indx[[i]] <- rep(tmp[i], re.long.indx[i])
      }
      indx <- unlist(indx)
      if (length(indx) == 0) {
        stop("error: column names in X.0 must match variable names in data$abund.covs")
      }
      n.abund.re <- length(indx)
      n.unique.abund.re <- length(unique(indx))
      # Check RE columns
      for (i in 1:n.abund.re) {
        if (is.character(re.cols[[i]])) {
          # Check if all column names in svc are in occ.covs
          if (!all(re.cols[[i]] %in% x.0.names)) {
              missing.cols <- re.cols[[i]][!(re.cols[[i]] %in% x.0.names)]
              stop(paste("error: variable name ", paste(missing.cols, collapse=" and "), " not in abundance covariates", sep=""))
          }
          # Convert desired column names into the numeric column index
          re.cols[[i]] <- which(x.0.names %in% re.cols[[i]])
          
        } else if (is.numeric(re.cols[[i]])) {
          # Check if all column indices are in 1:p.abund
          if (!all(re.cols %in% 1:p.abund)) {
              missing.cols <- re.cols[[i]][!(re.cols[[i]] %in% (1:p.abund))]
              stop(paste("error: column index ", paste(missing.cols, collapse=" "), " not in design matrix columns", sep=""))
          }
        }
      }
      re.cols <- unlist(re.cols)
      X.re <- as.matrix(X.0[, indx, drop = FALSE])
      X.fix <- as.matrix(X.0[, -indx, drop = FALSE])
      X.random <- as.matrix(X.0[, re.cols, drop = FALSE])
      n.abund.re <- length(unlist(re.level.names))
      X.re.ind <- matrix(NA, nrow(X.re), p.abund.re)
      for (i in 1:p.abund.re) {
        for (j in 1:nrow(X.re)) {
          tmp <- which(re.level.names[[i]] == X.re[j, i])
          if (length(tmp) > 0) {
            X.re.ind[j, i] <- tmp 
          }
        }
      }
      if (p.abund.re > 1) {
        for (j in 2:p.abund.re) {
          X.re.ind[, j] <- X.re.ind[, j] + max(X.re.ind[, j - 1]) 
        }
      }
      # Create the random effects corresponding to each 
      # new location
      # ORDER: ordered by site, then species within site.
      beta.star.sites.0.samples <- matrix(0, n.post, n.sp * nrow(X.re))
      for (i in 1:n.sp) {
        for (t in 1:p.abund.re) {
          for (j in 1:nrow(X.re)) {
            if (!is.na(X.re.ind[j, t])) {
              beta.star.sites.0.samples[, (j - 1) * n.sp + i] <- 
                beta.star.samples[, (i - 1) * n.abund.re + X.re.ind[j, t]] * X.random[j, t] + 
                beta.star.sites.0.samples[, (j - 1) * n.sp + i]
            } else {
              beta.star.sites.0.samples[, (j - 1) * n.sp + i] <- 
                rnorm(n.post, 0, sqrt(object$sigma.sq.mu.samples[, t])) * X.random[j, t] + 
                beta.star.sites.0.samples[, (j - 1) * n.sp + i]
            }
          } # j
        } # t
      } # i 
    } else {
      X.fix <- X.0
      beta.star.sites.0.samples <- matrix(0, n.post, n.sp * nrow(X.0))
      p.abund.re <- 0
    }
    J.str <- nrow(X.0)
    # Make predictions
    for (i in 1:n.sp) {
      for (j in 1:J.str) {
        out$mu.0.samples[, i, j] <- exp(t(as.matrix(X.fix[j, ])) %*% 
          				     t(beta.samples[, sp.indx == i]) + 
                                               beta.star.sites.0.samples[, (j - 1) * n.sp + i])
        if (object$dist == 'NB') {
          out$N.0.samples[, i, j] <- rnbinom(n.post, kappa.samples[, i], 
					     mu = out$mu.0.samples[, i, j])
	} else {
          out$N.0.samples[, i, j] <- rpois(n.post, out$mu.0.samples[, i, j]) 
	}
      } # j
    } # i
  } # abundance predictions
  # Detection predictions -------------------------------------------------
  if (tolower(type) == 'detection') {
    p.det <- ncol(object$X.p)
    re.det.cols <- object$re.det.cols
    # Composition sampling --------------------------------------------------
    alpha.samples <- as.matrix(object$alpha.samples)
    n.sp <- dim(object$y)[1]
    sp.indx <- rep(1:n.sp, p.det)
    n.post <- object$n.post * object$n.chains
    out <- list()
    out$p.0.samples <- array(NA, dim = c(n.post, n.sp, nrow(X.0)))
    if (object$pRE) {
      p.det.re <- length(object$p.re.level.names)
    } else {
      p.det.re <- 0
    }
    if (object$pRE & !ignore.RE) {
      alpha.star.samples <- object$alpha.star.samples
      p.re.level.names <- object$p.re.level.names
      # Get columns in design matrix with random effects
      x.p.re.names <- colnames(object$X.p.re)
      x.p.0.names <- colnames(X.0)
      re.long.indx <- sapply(re.det.cols, length)
      tmp <- sapply(x.p.re.names, function(a) which(colnames(X.0) %in% a))
      indx <- list()
      for (i in 1:length(tmp)) {
        indx[[i]] <- rep(tmp[i], re.long.indx[i])
      }
      indx <- unlist(indx)
      if (length(indx) == 0) {
        stop("error: column names in X.0 must match variable names in data$det.covs")
      }
      n.det.re <- length(indx)
      n.unique.det.re <- length(unique(indx))
      # Check RE columns
      for (i in 1:n.det.re) {
        if (is.character(re.det.cols[[i]])) {
          # Check if all column names in svc are in occ.covs
          if (!all(re.det.cols[[i]] %in% x.p.0.names)) {
              missing.cols <- re.det.cols[[i]][!(re.det.cols[[i]] %in% x.p.0.names)]
              stop(paste("error: variable name ", paste(missing.cols, collapse=" and "), " not in detection covariates", sep=""))
          }
          # Convert desired column names into the numeric column index
          re.det.cols[[i]] <- which(x.p.0.names %in% re.det.cols[[i]])
          
        } else if (is.numeric(re.det.cols[[i]])) {
          # Check if all column indices are in 1:p.abund
          if (!all(re.det.cols %in% 1:p.det)) {
              missing.cols <- re.det.cols[[i]][!(re.det.cols[[i]] %in% (1:p.det))]
              stop(paste("error: column index ", paste(missing.cols, collapse=" "), " not in design matrix columns", sep=""))
          }
        }
      }
      re.det.cols <- unlist(re.det.cols)
      X.re <- as.matrix(X.0[, indx, drop = FALSE])
      X.fix <- as.matrix(X.0[, -indx, drop = FALSE])
      X.random <- as.matrix(X.0[, re.det.cols, drop = FALSE])
      n.det.re <- length(unlist(p.re.level.names))
      X.re.ind <- matrix(NA, nrow(X.re), p.det.re)
      for (i in 1:p.det.re) {
        for (j in 1:nrow(X.re)) {
          tmp <- which(p.re.level.names[[i]] == X.re[j, i])
          if (length(tmp) > 0) {
            X.re.ind[j, i] <- tmp 
          }
        }
      }
      if (p.det.re > 1) {
        for (j in 2:p.det.re) {
          X.re.ind[, j] <- X.re.ind[, j] + max(X.re.ind[, j - 1]) 
        }
      }
      # Create the random effects corresponding to each 
      # new location
      # ORDER: ordered by site, then species within site.
      alpha.star.sites.0.samples <- matrix(0, n.post, n.sp * nrow(X.re))
      for (i in 1:n.sp) {
        for (t in 1:p.det.re) {
          for (j in 1:nrow(X.re)) {
            if (!is.na(X.re.ind[j, t])) {
              alpha.star.sites.0.samples[, (j - 1) * n.sp + i] <- 
                alpha.star.samples[, (i - 1) * n.det.re + X.re.ind[j, t]] * X.random[j, t] + 
                alpha.star.sites.0.samples[, (j - 1) * n.sp + i]
            } else {
              alpha.star.sites.0.samples[, (j - 1) * n.sp + i] <- 
                rnorm(n.post, 0, sqrt(object$sigma.sq.p.samples[, t])) * X.random[j, t] + 
                alpha.star.sites.0.samples[, (j - 1) * n.sp + i]
            }
          } # j
        } # t
      } # i 
    } else {
      X.fix <- X.0
      alpha.star.sites.0.samples <- matrix(0, n.post, n.sp * nrow(X.0))
      p.det.re <- 0
    }
    J.str <- nrow(X.0)
    # Make predictions
    for (i in 1:n.sp) {
      for (j in 1:J.str) {
        out$p.0.samples[, i, j] <- logit.inv(t(as.matrix(X.fix[j, ])) %*% 
          				     t(alpha.samples[, sp.indx == i]) + 
                                               alpha.star.sites.0.samples[, (j - 1) * n.sp + i])
      } # j
    } # i

  }
  out$call <- cl

  class(out) <- "predict.msNMix"
  out
}
# sfMsNMix ---------------------------------------------------------------
print.sfMsNMix <- function(x, ...) {
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.75)),
      "", sep = "\n")
}

summary.sfMsNMix <- function(object, level = 'both',
                             quantiles = c(0.025, 0.5, 0.975),
                             digits = max(3L, getOption("digits") - 3L), ...) {

  print(object)

  n.post <- object$n.post
  n.samples <- object$n.samples
  n.burn <- object$n.burn
  n.thin <- object$n.thin
  n.chains <- object$n.chains
  run.time <- object$run.time[3] / 60 # minutes

  cat(paste("Samples per Chain: ", n.samples,"\n", sep=""))
  cat(paste("Burn-in: ", n.burn,"\n", sep=""))
  cat(paste("Thinning Rate: ",n.thin,"\n", sep=""))
  cat(paste("Number of Chains: ", n.chains, "\n", sep = ""))
  cat(paste("Total Posterior Samples: ",n.post * n.chains,"\n", sep=""))
  cat(paste("Run Time (min): ", round(run.time, digits), "\n\n", sep = ""))

  if (tolower(level) %in% c('community', 'both')) {

    cat("----------------------------------------\n");
    cat("\tCommunity Level\n");
    cat("----------------------------------------\n");

    # Abundance 
    cat("Abundance Means (log scale): \n")
    tmp.1 <- t(apply(object$beta.comm.samples, 2,
          	   function(x) c(mean(x), sd(x))))
    colnames(tmp.1) <- c("Mean", "SD")
    tmp <- t(apply(object$beta.comm.samples, 2,
          	 function(x) quantile(x, prob = quantiles)))
    diags <- matrix(c(object$rhat$beta.comm, round(object$ESS$beta.comm, 0)), ncol = 2)
    colnames(diags) <- c('Rhat', 'ESS')

    print(noquote(round(cbind(tmp.1, tmp, diags), digits)))

    cat("\nAbundance Variances (log scale): \n")
    tmp.1 <- t(apply(object$tau.sq.beta.samples, 2,
          	   function(x) c(mean(x), sd(x))))
    colnames(tmp.1) <- c("Mean", "SD")
    tmp <- t(apply(object$tau.sq.beta.samples, 2,
          	 function(x) quantile(x, prob = quantiles)))
    diags <- matrix(c(object$rhat$tau.sq.beta, round(object$ESS$tau.sq.beta, 0)), ncol = 2)
    colnames(diags) <- c('Rhat', 'ESS')
    print(noquote(round(cbind(tmp.1, tmp, diags), digits)))

    if (object$muRE) {
      cat("\n")
      cat("Abundance Random Effect Variances (log scale): \n")
      tmp.1 <- t(apply(object$sigma.sq.mu.samples, 2,
            	   function(x) c(mean(x), sd(x))))
      colnames(tmp.1) <- c("Mean", "SD")
      tmp <- t(apply(object$sigma.sq.mu.samples, 2,
            	 function(x) quantile(x, prob = quantiles)))
      diags <- matrix(c(object$rhat$sigma.sq.mu, round(object$ESS$sigma.sq.mu, 0)), ncol = 2)
      colnames(diags) <- c('Rhat', 'ESS')

      print(noquote(round(cbind(tmp.1, tmp, diags), digits)))
    }
    cat("\n")
    # Detection
    if (is(object, 'sfMsNMix')) {
      cat("Detection Means (logit scale): \n")
    } else if (is(object, 'sfMsDS')) {
      cat("Detection Means (log scale): \n")
    }
    tmp.1 <- t(apply(object$alpha.comm.samples, 2,
          	   function(x) c(mean(x), sd(x))))
    colnames(tmp.1) <- c("Mean", "SD")
    tmp <- t(apply(object$alpha.comm.samples, 2,
          	 function(x) quantile(x, prob = quantiles)))
    diags <- matrix(c(object$rhat$alpha.comm, round(object$ESS$alpha.comm, 0)), ncol = 2)
    colnames(diags) <- c('Rhat', 'ESS')

    print(noquote(round(cbind(tmp.1, tmp, diags), digits)))
    if (is(object, 'sfMsNMix')) {
      cat("Detection Variances (logit scale): \n")
    } else if (is(object, 'sfMsDS')) {
      cat("Detection Variances (log scale): \n")
    }
    tmp.1 <- t(apply(object$tau.sq.alpha.samples, 2,
          	   function(x) c(mean(x), sd(x))))
    colnames(tmp.1) <- c("Mean", "SD")
    tmp <- t(apply(object$tau.sq.alpha.samples, 2,
          	 function(x) quantile(x, prob = quantiles)))
    diags <- matrix(c(object$rhat$tau.sq.alpha, round(object$ESS$tau.sq.alpha, 0)), ncol = 2)
    colnames(diags) <- c('Rhat', 'ESS')

    print(noquote(round(cbind(tmp.1, tmp, diags), digits)))
    if (object$pRE) {
      cat("\n")
      if (is(object, 'sfMsNMix')) {
        cat("Detection Random Effect Variances (logit scale): \n")
      } else if (is(object, 'sfMsDS')) {
        cat("Detection Random Effect Variances (log scale): \n")
      }
      tmp.1 <- t(apply(object$sigma.sq.p.samples, 2,
            	   function(x) c(mean(x), sd(x))))
      colnames(tmp.1) <- c("Mean", "SD")
      tmp <- t(apply(object$sigma.sq.p.samples, 2,
            	 function(x) quantile(x, prob = quantiles)))
      diags <- matrix(c(object$rhat$sigma.sq.p, round(object$ESS$sigma.sq.p, 0)), ncol = 2)
      colnames(diags) <- c('Rhat', 'ESS')

      print(noquote(round(cbind(tmp.1, tmp, diags), digits)))
    }
  }

  if (tolower(level) %in% c('species', 'both')) {
    if (tolower(level) == 'both') cat("\n")
    cat("----------------------------------------\n");
    cat("\tSpecies Level\n");
    cat("----------------------------------------\n");
    cat("Abundance (log scale): \n")
    tmp.1 <- t(apply(object$beta.samples, 2,
          	   function(x) c(mean(x), sd(x))))
    colnames(tmp.1) <- c("Mean", "SD")
    tmp <- t(apply(object$beta.samples, 2,
          	 function(x) quantile(x, prob = quantiles)))
    diags <- matrix(c(object$rhat$beta, round(object$ESS$beta, 0)), ncol = 2)
    colnames(diags) <- c('Rhat', 'ESS')

    print(noquote(round(cbind(tmp.1, tmp, diags), digits)))
    cat("\n")
    # Detection
    if (is(object, 'sfMsNMix')) {
      cat("Detection (logit scale): \n")
    } else if (is(object, 'sfMsDS')) {
      cat("Detection (log scale): \n")
    }
    tmp.1 <- t(apply(object$alpha.samples, 2,
          	   function(x) c(mean(x), sd(x))))
    colnames(tmp.1) <- c("Mean", "SD")
    tmp <- t(apply(object$alpha.samples, 2,
          	 function(x) quantile(x, prob = quantiles)))
    diags <- matrix(c(object$rhat$alpha, round(object$ESS$alpha, 0)), ncol = 2)
    colnames(diags) <- c('Rhat', 'ESS')
    print(noquote(round(cbind(tmp.1, tmp, diags), digits)))

  }
    # Covariance
    cat("\n")
    cat("----------------------------------------\n");
    cat("\tSpatial Covariance\n");
    cat("----------------------------------------\n");
    tmp.1 <- t(apply(object$theta.samples, 2,
          	   function(x) c(mean(x), sd(x))))
    colnames(tmp.1) <- c("Mean", "SD")
    tmp <- t(apply(object$theta.samples, 2,
          	 function(x) quantile(x, prob = quantiles)))
    diags <- matrix(c(object$rhat$theta, round(object$ESS$theta, 0)), ncol = 2)
    colnames(diags) <- c('Rhat', 'ESS')
    print(noquote(round(cbind(tmp.1, tmp, diags), digits)))

    # NB Overdispersion
    if (object$dist == "NB" & tolower(level) %in% c('species', 'both')) {
      cat("\n")
      cat("----------------------------------------\n");
      cat("\tNB overdispersion\n");
      cat("----------------------------------------\n");
      tmp.1 <- t(apply(object$kappa.samples, 2,
            	   function(x) c(mean(x), sd(x))))
      colnames(tmp.1) <- c("Mean", "SD")
      tmp <- t(apply(object$kappa.samples, 2,
            	 function(x) quantile(x, prob = quantiles)))
      diags <- matrix(c(object$rhat$kappa, round(object$ESS$kappa, 0)), ncol = 2)
      colnames(diags) <- c('Rhat', 'ESS')
      print(noquote(round(cbind(tmp.1, tmp, diags), digits)))
    }
}
# 
fitted.sfMsNMix <- function(object, type = 'marginal', ...) {
  fitted.msNMix(object, type)
}

predict.sfMsNMix <- function(object, X.0, coords.0, n.omp.threads = 1,
			      verbose = TRUE, n.report = 100,
			      ignore.RE = FALSE, type = 'abundance', 
			      include.sp = TRUE, ...) {

  # Check for unused arguments ------------------------------------------
  formal.args <- names(formals(sys.function(sys.parent())))
  elip.args <- names(list(...))
  for(i in elip.args){
      if(! i %in% formal.args)
          warning("'",i, "' is not an argument")
  }
  # Call ----------------------------------------------------------------
  cl <- match.call()

  # Call predict.msNMix if you don't care about spatial effects
  if (!include.sp) {
    out <- predict.msNMix(object, X.0, ignore.RE, type)
  } else {

    # Functions ---------------------------------------------------------------
    logit <- function(theta, a = 0, b = 1) {log((theta-a)/(b-theta))}
    logit.inv <- function(z, a = 0, b = 1) {b-(b-a)/(1+exp(z))}

    # Some initial checks ---------------------------------------------------
    if (missing(object)) {
      stop("error: predict expects object\n")
    }
    if (!(tolower(type) %in% c('abundance', 'detection'))) {
      stop("error: prediction type must be either 'abundance' or 'detection'")
    }

    # Check X.0 -------------------------------------------------------------
    if (missing(X.0)) {
      stop("error: X.0 must be specified\n")
    }
    if (!any(is.data.frame(X.0), is.matrix(X.0))) {
      stop("error: X.0 must be a data.frame or matrix\n")
    }

    ptm <- proc.time()

    # Abundance predictions ------------------------------------------------
    if (tolower(type == 'abundance')) {
      n.post <- object$n.post * object$n.chains
      X <- object$X
      y <- object$y
      coords <- object$coords
      J <- nrow(X)
      n.sp <- dim(y)[1]
      q <- object$q
      p.abund <- ncol(X)
      theta.samples <- object$theta.samples
      beta.samples <- object$beta.samples
      lambda.samples <- object$lambda.samples
      w.samples <- object$w.samples
      kappa.samples <- object$kappa.samples
      n.neighbors <- object$n.neighbors
      cov.model.indx <- object$cov.model.indx
      family <- object$dist
      family.c <- ifelse(family == 'NB', 1, 0)
      sp.type <- object$type
      if (object$muRE & !ignore.RE) {
        p.abund.re <- length(object$re.level.names)
      } else {
        p.abund.re <- 0
      }
      re.cols <- object$re.cols

      # if (ncol(X.0) != p.abund + p.abund.re){
      #   stop(paste("error: X.0 must have ", p.abund + p.abund.re," columns\n"))
      # }
      X.0 <- as.matrix(X.0)

      if (missing(coords.0)) {
        stop("error: coords.0 must be specified\n")
      }
      if (!any(is.data.frame(coords.0), is.matrix(coords.0))) {
        stop("error: coords.0 must be a data.frame or matrix\n")
      }
      if (!ncol(coords.0) == 2){
        stop("error: coords.0 must have two columns\n")
      }
      coords.0 <- as.matrix(coords.0)

      # Eliminate prediction sites that have already been sampled for now
      match.indx <- match(do.call("paste", as.data.frame(coords.0)), do.call("paste", as.data.frame(coords)))
      coords.0.indx <- which(is.na(match.indx))
      coords.indx <- match.indx[!is.na(match.indx)]
      coords.place.indx <- which(!is.na(match.indx))
      coords.0.new <- coords.0[coords.0.indx, , drop = FALSE]
      X.0.new <- X.0[coords.0.indx, , drop = FALSE]

      if (length(coords.indx) == nrow(X.0)) {
        stop("error: no new locations to predict at. See object$mu.samples for expected abundances at sampled sites.")
      }

      if (object$muRE & !ignore.RE) {
        beta.star.samples <- object$beta.star.samples
        re.level.names <- object$re.level.names
        # Get columns in design matrix with random effects
        x.re.names <- colnames(object$X.re)
        x.0.names <- colnames(X.0.new)
        re.long.indx <- sapply(re.cols, length)
        tmp <- sapply(x.re.names, function(a) which(colnames(X.0.new) %in% a))
        indx <- list()
        for (i in 1:length(tmp)) {
          indx[[i]] <- rep(tmp[i], re.long.indx[i])
        }
        indx <- unlist(indx)
        if (length(indx) == 0) {
          stop("error: column names in X.0 must match variable names in data$abund.covs")
        }
        n.abund.re <- length(indx)
        n.unique.abund.re <- length(unique(indx))
        # Check RE columns
        for (i in 1:n.abund.re) {
          if (is.character(re.cols[[i]])) {
            # Check if all column names in svc are in occ.covs
            if (!all(re.cols[[i]] %in% x.0.names)) {
                missing.cols <- re.cols[[i]][!(re.cols[[i]] %in% x.0.names)]
                stop(paste("error: variable name ", paste(missing.cols, collapse=" and "), " not in abundance covariates", sep=""))
            }
            # Convert desired column names into the numeric column index
            re.cols[[i]] <- which(x.0.names %in% re.cols[[i]])
            
          } else if (is.numeric(re.cols[[i]])) {
            # Check if all column indices are in 1:p.abund
            if (!all(re.cols %in% 1:p.abund)) {
                missing.cols <- re.cols[[i]][!(re.cols[[i]] %in% (1:p.abund))]
                stop(paste("error: column index ", paste(missing.cols, collapse=" "), " not in design matrix columns", sep=""))
            }
          }
        }
        re.cols <- unlist(re.cols)
        X.re <- as.matrix(X.0.new[, indx, drop = FALSE])
        X.fix <- as.matrix(X.0.new[, -indx, drop = FALSE])
        X.random <- as.matrix(X.0.new[, re.cols, drop = FALSE])
        n.abund.re <- length(unlist(re.level.names))
        X.re.ind <- matrix(NA, nrow(X.re), p.abund.re)
        for (i in 1:p.abund.re) {
          for (j in 1:nrow(X.re)) {
            tmp <- which(re.level.names[[i]] == X.re[j, i])
            if (length(tmp) > 0) {
              X.re.ind[j, i] <- tmp 
            }
          }
        }
        if (p.abund.re > 1) {
          for (j in 2:p.abund.re) {
            X.re.ind[, j] <- X.re.ind[, j] + max(X.re.ind[, j - 1]) 
          }
        }
        # Create the random effects corresponding to each 
        # new location
        # ORDER: ordered by site, then species within site.
        beta.star.sites.0.samples <- matrix(0, n.post, n.sp * nrow(X.re))
        for (i in 1:n.sp) {
          for (t in 1:p.abund.re) {
            for (j in 1:nrow(X.re)) {
              if (!is.na(X.re.ind[j, t])) {
                beta.star.sites.0.samples[, (j - 1) * n.sp + i] <- 
                  beta.star.samples[, (i - 1) * n.abund.re + X.re.ind[j, t]] * X.random[j, t] + 
                  beta.star.sites.0.samples[, (j - 1) * n.sp + i]
              } else {
                beta.star.sites.0.samples[, (j - 1) * n.sp + i] <- 
                  rnorm(n.post, 0, sqrt(object$sigma.sq.mu.samples[, t])) * X.random[j, t] + 
                  beta.star.sites.0.samples[, (j - 1) * n.sp + i]
              }
            } # j
          } # t
        } # i 
      } else {
        X.fix <- X.0.new
        beta.star.sites.0.samples <- matrix(0, n.post, n.sp * nrow(X.0.new))
        p.abund.re <- 0
      }

      # Sub-sample previous
      theta.samples <- t(theta.samples)
      lambda.samples <- t(lambda.samples)
      beta.samples <- t(beta.samples)
      if (family == 'NB') {
        kappa.samples <- t(kappa.samples)
      }
      w.samples <- aperm(w.samples, c(2, 3, 1))
      beta.star.sites.0.samples <- t(beta.star.sites.0.samples)

      J.str <- nrow(X.0.new)

      if (sp.type == 'GP') {
        # Not currently implemented or accessed.
      } else {
        # Get nearest neighbors
        # nn2 is a function from RANN.
        nn.indx.0 <- nn2(coords, coords.0.new, k=n.neighbors)$nn.idx-1

        storage.mode(coords) <- "double"
        storage.mode(n.sp) <- "integer"
        storage.mode(J) <- "integer"
        storage.mode(p.abund) <- "integer"
        storage.mode(n.neighbors) <- "integer"
        storage.mode(X.fix) <- "double"
        storage.mode(coords.0.new) <- "double"
        storage.mode(J.str) <- "integer"
        storage.mode(q) <- "integer"
        storage.mode(beta.samples) <- "double"
        storage.mode(theta.samples) <- "double"
        storage.mode(lambda.samples) <- "double"
        storage.mode(kappa.samples) <- "double"
        storage.mode(beta.star.sites.0.samples) <- "double"
        storage.mode(w.samples) <- "double"
        storage.mode(n.post) <- "integer"
        storage.mode(cov.model.indx) <- "integer"
        storage.mode(nn.indx.0) <- "integer"
        storage.mode(n.omp.threads) <- "integer"
        storage.mode(verbose) <- "integer"
        storage.mode(n.report) <- "integer"
        storage.mode(family.c) <- "integer"

        out <- .Call("sfMsNMixNNGPPredict", coords, J, family.c, 
          	   n.sp, q, p.abund, n.neighbors, 
                     X.fix, coords.0.new, J.str, nn.indx.0, beta.samples,
                     theta.samples, kappa.samples, lambda.samples, w.samples,
            	   beta.star.sites.0.samples, n.post,
                     cov.model.indx, n.omp.threads, verbose, n.report)
      }
      out$N.0.samples <- array(out$N.0.samples, dim = c(n.sp, J.str, n.post))
      out$N.0.samples <- aperm(out$N.0.samples, c(3, 1, 2))
      out$w.0.samples <- array(out$w.0.samples, dim = c(q, J.str, n.post))
      out$w.0.samples <- aperm(out$w.0.samples, c(3, 1, 2))
      out$mu.0.samples <- array(out$mu.0.samples, dim = c(n.sp, J.str, n.post))
      out$mu.0.samples <- aperm(out$mu.0.samples, c(3, 1, 2))

      # If some of the sites are sampled
      if (nrow(X.0) != J.str) {
        tmp <- array(NA, dim = c(n.post, n.sp, nrow(X.0)))
        tmp[, , coords.0.indx] <- out$N.0.samples
        tmp[, , coords.place.indx] <- object$N.samples[, , coords.indx]
        out$N.0.samples <- tmp
        tmp <- array(NA, dim = c(n.post, n.sp, nrow(X.0)))
        tmp[, , coords.0.indx] <- out$mu.0.samples
        tmp[, , coords.place.indx] <- object$mu.samples[, , coords.indx]
        out$mu.0.samples <- tmp
        tmp <- array(NA, dim = c(n.post, q, nrow(X.0)))
        tmp[, , coords.0.indx] <- out$w.0.samples
        tmp[, , coords.place.indx] <- object$w.samples[, , coords.indx]
        out$w.0.samples <- tmp
      }
    } # occurrence predictions
    # Detection predictions -------------------------------------------------
    if (tolower(type) == 'detection') {
      out <- predict.msNMix(object, X.0, ignore.RE, type)
    }
  }

  out$call <- cl
  out$object.class <- class(object)

  class(out) <- "predict.sfMsNMix"

  out

}

# lfMsNMix ----------------------------------------------------------------
print.lfMsNMix <- function(x, ...) {
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.75)),
      "", sep = "\n")
}

summary.lfMsNMix <- function(object,
			     level = 'both',
			     quantiles = c(0.025, 0.5, 0.975),
			     digits = max(3L, getOption("digits") - 3L), ...) {
  summary.msNMix(object, level, quantiles, digits)
}

fitted.lfMsNMix <- function(object, type = 'marginal', ...) {
  fitted.msNMix(object, type)
}
predict.lfMsNMix <- function(object, X.0, coords.0, ignore.RE = FALSE, 
			  type = 'abundance', include.w = TRUE, ...) {
  # Check for unused arguments ------------------------------------------
  formal.args <- names(formals(sys.function(sys.parent())))
  elip.args <- names(list(...))
  for(i in elip.args){
      if(! i %in% formal.args)
          warning("'",i, "' is not an argument")
  }
  # Call ----------------------------------------------------------------
  cl <- match.call()

  if (!include.w) {
    out <- predict.msNMix(object, X.0, ignore.RE, type)
  } else {

    # Functions ---------------------------------------------------------------
    logit <- function(theta, a = 0, b = 1) {log((theta-a)/(b-theta))}
    logit.inv <- function(z, a = 0, b = 1) {b-(b-a)/(1+exp(z))}

    # Some initial checks ---------------------------------------------------
    if (missing(object)) {
      stop("error: predict expects object\n")
    }

    if (!(tolower(type) %in% c('abundance', 'detection'))) {
      stop("error: prediction type must be either 'abundance' or 'detection'")
    }

    # Check X.0 -------------------------------------------------------------
    if (missing(X.0)) {
      stop("error: X.0 must be specified\n")
    }
    if (!any(is.data.frame(X.0), is.matrix(X.0))) {
      stop("error: X.0 must be a data.frame or matrix\n")
    }

    # Abundance predictions ------------------------------------------------
    if (tolower(type) == 'abundance') {
      p.abund <- ncol(object$X)
      # Composition sampling --------------------------------------------------
      beta.samples <- as.matrix(object$beta.samples)
      if (object$dist == 'NB') {
        kappa.samples <- as.matrix(object$kappa.samples)
      }
      n.sp <- nrow(object$y)
      J.0 <- nrow(X.0)
      q <- object$q
      sp.indx <- rep(1:n.sp, p.abund)
      n.post <- object$n.post * object$n.chains
      out <- list()
      coords <- out$coords
      out$mu.0.samples <- array(NA, dim = c(n.post, n.sp, nrow(X.0)))
      out$N.0.samples <- array(NA, dim = c(n.post, n.sp, nrow(X.0)))
      lambda.samples <- array(object$lambda.samples, dim = c(n.post, n.sp, q))
      w.samples <- object$w.samples
      if (object$muRE) {
        p.abund.re <- length(object$re.level.names)
      } else {
        p.abund.re <- 0
      }
      re.cols <- object$re.cols

      # Eliminate prediction sites that have already been sampled for now
      if (!missing(coords.0)) {
        match.indx <- match(do.call("paste", as.data.frame(coords.0)), 
              	      do.call("paste", as.data.frame(coords)))
        coords.0.indx <- which(is.na(match.indx))
        coords.indx <- match.indx[!is.na(match.indx)]
        coords.place.indx <- which(!is.na(match.indx))
        coords.0.new <- coords.0[coords.0.indx, , drop = FALSE]
        X.0.new <- X.0[coords.0.indx, , drop = FALSE]
        if (length(coords.indx) == nrow(X.0)) {
          stop("error: no new locations to predict at. See object$psi.samples for occurrence probabilities at sampled sites.")
        }
      } else {
        X.0.new <- X.0
      }


      if (object$muRE & !ignore.RE) {
        beta.star.samples <- object$beta.star.samples
        re.level.names <- object$re.level.names
        # Get columns in design matrix with random effects
        x.re.names <- colnames(object$X.re)
        x.0.names <- colnames(X.0.new)
        re.long.indx <- sapply(re.cols, length)
        tmp <- sapply(x.re.names, function(a) which(colnames(X.0.new) %in% a))
        indx <- list()
        for (i in 1:length(tmp)) {
          indx[[i]] <- rep(tmp[i], re.long.indx[i])
        }
        indx <- unlist(indx)
        if (length(indx) == 0) {
          stop("error: column names in X.0 must match variable names in data$abund.covs")
        }
        n.abund.re <- length(indx)
        n.unique.abund.re <- length(unique(indx))
        # Check RE columns
        for (i in 1:n.abund.re) {
          if (is.character(re.cols[[i]])) {
            # Check if all column names in svc are in occ.covs
            if (!all(re.cols[[i]] %in% x.0.names)) {
                missing.cols <- re.cols[[i]][!(re.cols[[i]] %in% x.0.names)]
                stop(paste("error: variable name ", paste(missing.cols, collapse=" and "), " not in abundance covariates", sep=""))
            }
            # Convert desired column names into the numeric column index
            re.cols[[i]] <- which(x.0.names %in% re.cols[[i]])
            
          } else if (is.numeric(re.cols[[i]])) {
            # Check if all column indices are in 1:p.abund
            if (!all(re.cols %in% 1:p.abund)) {
                missing.cols <- re.cols[[i]][!(re.cols[[i]] %in% (1:p.abund))]
                stop(paste("error: column index ", paste(missing.cols, collapse=" "), " not in design matrix columns", sep=""))
            }
          }
        }
        re.cols <- unlist(re.cols)
        X.re <- as.matrix(X.0.new[, indx, drop = FALSE])
        X.fix <- as.matrix(X.0.new[, -indx, drop = FALSE])
        X.random <- as.matrix(X.0.new[, re.cols, drop = FALSE])
        n.abund.re <- length(unlist(re.level.names))
        X.re.ind <- matrix(NA, nrow(X.re), p.abund.re)
        for (i in 1:p.abund.re) {
          for (j in 1:nrow(X.re)) {
            tmp <- which(re.level.names[[i]] == X.re[j, i])
            if (length(tmp) > 0) {
              X.re.ind[j, i] <- tmp 
            }
          }
        }
        if (p.abund.re > 1) {
          for (j in 2:p.abund.re) {
            X.re.ind[, j] <- X.re.ind[, j] + max(X.re.ind[, j - 1]) 
          }
        }
        # Create the random effects corresponding to each 
        # new location
        # ORDER: ordered by site, then species within site.
        beta.star.sites.0.samples <- matrix(0, n.post, n.sp * nrow(X.re))
        for (i in 1:n.sp) {
          for (t in 1:p.abund.re) {
            for (j in 1:nrow(X.re)) {
              if (!is.na(X.re.ind[j, t])) {
                beta.star.sites.0.samples[, (j - 1) * n.sp + i] <- 
                  beta.star.samples[, (i - 1) * n.abund.re + X.re.ind[j, t]] * X.random[j, t] + 
                  beta.star.sites.0.samples[, (j - 1) * n.sp + i]
              } else {
                beta.star.sites.0.samples[, (j - 1) * n.sp + i] <- 
                  rnorm(n.post, 0, sqrt(object$sigma.sq.mu.samples[, t])) * X.random[j, t] + 
                  beta.star.sites.0.samples[, (j - 1) * n.sp + i]
              }
            } # j
          } # t
        } # i 
      } else {
        X.fix <- X.0.new
        beta.star.sites.0.samples <- matrix(0, n.post, n.sp * nrow(X.0.new))
        p.abund.re <- 0
      }
      J.str <- nrow(X.0.new)
      # Create new random normal latent factors at unobserved sites. 
      w.0.samples <- array(rnorm(n.post * q * J.str), dim = c(n.post, q, J.str))
      w.star.0.samples <- array(NA, dim = c(n.post, n.sp, J.str))

      for (i in 1:n.post) {
        w.star.0.samples[i, , ] <- matrix(lambda.samples[i, , ], n.sp, q) %*%
                                 matrix(w.0.samples[i, , ], q, J.str)
      }
      # Make predictions
      for (i in 1:n.sp) {
        for (j in 1:J.str) {
          out$mu.0.samples[, i, j] <- exp(t(as.matrix(X.fix[j, ])) %*% 
            				     t(beta.samples[, sp.indx == i]) + 
          				     w.star.0.samples[, i, j] + 
                                                 beta.star.sites.0.samples[, (j - 1) * n.sp + i])
          if (object$dist == 'NB') {
            out$N.0.samples[, i, j] <- rnbinom(n.post, kappa.samples[, i], 
          				     mu = out$mu.0.samples[, i, j])
          } else {
            out$N.0.samples[, i, j] <- rpois(n.post, out$mu.0.samples[, i, j]) 
          }
        } # j
      } # i

      # If some of the sites are sampled
      if (nrow(X.0) != J.str) {
        tmp <- array(NA, dim = c(n.post, n.sp, nrow(X.0)))
        tmp[, , coords.0.indx] <- out$N.0.samples
        tmp[, , coords.place.indx] <- object$N.samples[, , coords.indx]
        out$N.0.samples <- tmp
        tmp <- array(NA, dim = c(n.post, n.sp, nrow(X.0)))
        tmp[, , coords.0.indx] <- out$mu.0.samples
        tmp[, , coords.place.indx] <- object$mu.samples[, , coords.indx]
        out$mu.0.samples <- tmp
      }
    } # abundance predictions
    # Detection predictions -------------------------------------------------
    if (tolower(type) == 'detection') {
      out <- predict.msNMix(object, X.0, ignore.RE, type)
    }
  }
  out$call <- cl

  class(out) <- "predict.lfMsNMix"
  out
}

# lfMsAbund ----------------------------------------------------------------
print.lfMsAbund <- function(x, ...) {
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.75)),
      "", sep = "\n")
}

summary.lfMsAbund <- function(object,
			     level = 'both',
			     quantiles = c(0.025, 0.5, 0.975),
			     digits = max(3L, getOption("digits") - 3L), ...) {
  summary.msAbund(object, level, quantiles, digits)
}

fitted.lfMsAbund <- function(object, ...) {
  return(object$y.rep.samples)
}

predict.lfMsAbund <- function(object, X.0, coords.0, ignore.RE = FALSE, 
			      z.0.samples, include.w = TRUE, ...) {
  # Check for unused arguments ------------------------------------------
  formal.args <- names(formals(sys.function(sys.parent())))
  elip.args <- names(list(...))
  for(i in elip.args){
      if(! i %in% formal.args)
          warning("'",i, "' is not an argument")
  }
  # Call ----------------------------------------------------------------
  cl <- match.call()
  
  # Call predict.msAbund if don't care about latent effects
  if (!include.w) {
    out <- predict.msAbund(object, X.0, ignore.RE, z.0.samples)
  } else {

    # Functions ---------------------------------------------------------------
    logit <- function(theta, a = 0, b = 1) {log((theta-a)/(b-theta))}
    logit.inv <- function(z, a = 0, b = 1) {b-(b-a)/(1+exp(z))}

    # Some initial checks ---------------------------------------------------
    if (missing(object)) {
      stop("error: predict expects object\n")
    }
    if (!is(object, "lfMsAbund")) {
      stop("error: requires an output object of class lfMsAbund\n")
    }

    # Check X.0 -------------------------------------------------------------
    if (missing(X.0)) {
      stop("error: X.0 must be specified\n")
    }
    if (!(length(dim(X.0)) %in% c(2, 3))) {
      stop("error: X.0 must be a matrix with two columns corresponding to site and covariate or a three-dimensional array with dimensions corresponding to site, replicate, and covariate")
    }
    if (length(dim(X.0)) == 2) {
      tmp <- colnames(X.0) 
      X.0 <- array(X.0, dim = c(nrow(X.0), 1, ncol(X.0)))
      dimnames(X.0)[[3]] <- tmp
    }
    # Predictions -----------------------------------------------------------
    if (object$dist %in% c('Gaussian', 'zi-Gaussian')) {
      p.abund <- dim(object$X)[2]
    } else {
      p.abund <- dim(object$X)[3]
    }
    J.0 <- nrow(X.0)
    K.max.0 <- ncol(X.0)
    n.sp <- dim(object$y)[1]
    q <- object$q
    n.post <- object$n.post * object$n.chains
    if (missing(coords.0)) {
      stop("error: coords.0 must be specified\n")
    }
    if (!any(is.data.frame(coords.0), is.matrix(coords.0))) {
      stop("error: coords.0 must be a data.frame or matrix\n")
    }
    if (!ncol(coords.0) == 2){
      stop("error: coords.0 must have two columns\n")
    }
    coords.0 <- as.matrix(coords.0)
    sites.0.indx <- rep(1:(nrow(X.0)), times = ncol(X.0))

    # Composition sampling --------------------------------------------------
    re.cols <- object$re.cols
    beta.samples <- as.matrix(object$beta.samples)
    if (object$dist == 'NB') {
      kappa.samples <- as.matrix(object$kappa.samples)
    }
    lambda.samples <- array(object$lambda.samples, dim = c(n.post, n.sp, q))
    w.samples <- object$w.samples
    out <- list()
    if (object$muRE) {
      p.abund.re <- length(object$re.level.names)
    } else {
      p.abund.re <- 0
    }

    # Get X.0 in long format. 
    tmp.names <- dimnames(X.0)[[3]]
    X.0 <- matrix(X.0, nrow = nrow(X.0) * ncol(X.0), ncol = dim(X.0)[3])
    colnames(X.0) <- tmp.names
    missing.indx <- which(apply(X.0, 1, function(a) sum(is.na(a)) > 0))
    non.missing.indx <- which(apply(X.0, 1, function(a) sum(is.na(a)) == 0))
    X.0 <- X.0[non.missing.indx, , drop = FALSE]
    sites.0.indx <- sites.0.indx[non.missing.indx]

    if (object$muRE & !ignore.RE) {
      beta.star.samples <- object$beta.star.samples
      re.level.names <- object$re.level.names
      # Get columns in design matrix with random effects
      # Get columns in design matrix with random effects
      if (object$dist %in% c('Gaussian', 'zi-Gaussian')) {
        x.re.names <- dimnames(object$X.re)[[2]]
      } else {
        x.re.names <- dimnames(object$X.re)[[3]]
      }
      x.0.names <- colnames(X.0)
      # Get the number of times each factor is used. 
      re.long.indx <- sapply(re.cols, length)
      tmp <- sapply(x.re.names, function(a) which(colnames(X.0) %in% a))
      indx <- list()
      for (i in 1:length(tmp)) {
        indx[[i]] <- rep(tmp[i], re.long.indx[i])
      }
      indx <- unlist(indx)
      if (length(indx) == 0) {
        stop("error: column names in X.0 must match variable names in data$occ.covs")
      }
      n.re <- length(indx)
      n.unique.re <- length(unique(indx))
      # Check RE columns
      for (i in 1:n.re) {
        if (is.character(re.cols[[i]])) {
          # Check if all column names in svc are in occ.covs
          if (!all(re.cols[[i]] %in% x.0.names)) {
              missing.cols <- re.cols[[i]][!(re.cols[[i]] %in% x.0.names)]
              stop(paste("error: variable name ", paste(missing.cols, collapse=" and "), " not in occurrence covariates", sep=""))
          }
          # Convert desired column names into the numeric column index
          re.cols[[i]] <- which(x.0.names %in% re.cols[[i]])
          
        } else if (is.numeric(re.cols[[i]])) {
          # Check if all column indices are in 1:p.abund
          if (!all(re.cols %in% 1:p.abund)) {
              missing.cols <- re.cols[[i]][!(re.cols[[i]] %in% (1:p.abund))]
              stop(paste("error: column index ", paste(missing.cols, collapse=" "), " not in design matrix columns", sep=""))
          }
        }
      }
      re.cols <- unlist(re.cols)
      X.re <- as.matrix(X.0[, indx, drop = FALSE])
      X.fix <- as.matrix(X.0[, -indx, drop = FALSE])
      X.random <- as.matrix(X.0[, re.cols, drop = FALSE])
      n.abund.re <- length(unlist(re.level.names))
      X.re.ind <- matrix(NA, nrow(X.re), p.abund.re)

      for (i in 1:p.abund.re) {
        for (j in 1:nrow(X.re)) {
          tmp <- which(re.level.names[[i]] == X.re[j, i])
          if (length(tmp) > 0) {
            X.re.ind[j, i] <- tmp 
          }
        }
      }
      if (p.abund.re > 1) {
        for (j in 2:p.abund.re) {
          X.re.ind[, j] <- X.re.ind[, j] + max(X.re.ind[, j - 1]) 
        }
      }
      # Create the random effects corresponding to each 
      # new location
      beta.star.sites.0.samples <- array(0, dim = c(n.post, n.sp, nrow(X.re)))
      for (i in 1:n.sp) {
        for (t in 1:p.abund.re) {
          for (j in 1:nrow(X.re)) {
            if (!is.na(X.re.ind[j, t])) {
              beta.star.sites.0.samples[, i, j] <- 
                beta.star.samples[, (i - 1) * n.abund.re + X.re.ind[j, t]] * X.random[j, t] + 
                beta.star.sites.0.samples[, i, j]
            } else {
              beta.star.sites.0.samples[, i, j] <- 
                rnorm(n.post, 0, sqrt(object$sigma.sq.mu.samples[, t])) * X.random[j, t] + 
                beta.star.sites.0.samples[, i, j]
            }
          } # j
        } # t
      }
    } else {
      X.fix <- X.0
      beta.star.sites.0.samples <- array(0, dim = c(n.post, n.sp, nrow(X.0)))
      p.abund.re <- 0
    }

    coords <- object$coords
    match.indx <- match(do.call("paste", as.data.frame(coords.0)), do.call("paste", as.data.frame(coords)))
    coords.0.indx <- which(is.na(match.indx))
    coords.indx <- match.indx[!is.na(match.indx)]
    coords.place.indx <- which(!is.na(match.indx))
    # Indicates whether a site has been sampled. 1 = sampled
    sites.0.sampled <- ifelse(!is.na(match.indx), 1, 0)
    sites.link <- rep(NA, J.0)
    sites.link[which(!is.na(match.indx))] <- coords.indx
    family <- object$dist

    if (family == 'zi-Gaussian') {
      if (missing(z.0.samples)) {
        stop("z.0.samples must be supplied for a zi-Gaussian model")
      }
      if (length(dim(z.0.samples)) != 3) {
        stop(paste("z.0.samples must be a three-dimensional array with dimensions of ", 
          	 n.post, ", ", n.sp, ", and ", J.0, sep = ''))
      }
      if (dim(z.0.samples)[1] != n.post | dim(z.0.samples)[2] != n.sp |
          dim(z.0.samples)[3] != J.0) {
        stop(paste("z.0.samples must be a three-dimensional array with dimensions of ", 
          	 n.post, ", ", n.sp, ", and ", J.0, sep = ''))
      }
    } else {
      if (!missing(z.0.samples)) {
        message("z.0.samples is ignored for the current model family\n")
      }
      z.0.samples <- array(NA, dim = c(1, 1, 1))
    }

    # Create new random normal latent factors at all sites.
    w.0.samples <- array(rnorm(n.post * q * J.0), dim = c(n.post, q, J.0))
    # Replace already sampled sites with values from model fit.
    for (j in 1:J.0) {
      if (sites.0.sampled[j] == 1) {
        w.0.samples[, , j] <- w.samples[, , sites.link[j]]
      }
    }
    w.star.0.samples <- array(NA, dim = c(n.post, n.sp, J.0))

    for (i in 1:n.post) {
      w.star.0.samples[i, , ] <- matrix(lambda.samples[i, , ], n.sp, q) %*%
                               matrix(w.0.samples[i, , ], q, J.0)
    }

    tmp <- matrix(NA, n.post, length(c(missing.indx, non.missing.indx)))
    sp.indx <- rep(1:n.sp, p.abund)
    out$mu.0.samples <- array(NA, dim = c(n.post, n.sp, J.0, K.max.0))
    mu.long <- array(NA, dim = c(n.post, n.sp, nrow(X.fix)))
    for (i in 1:n.sp) {
      if (object$dist %in% c('Poisson', 'NB')) {
        mu.long[, i, ] <- exp(t(X.fix %*% t(beta.samples[, sp.indx == i]) + 
                              t(beta.star.sites.0.samples[, i, ]) + 
          		    t(w.star.0.samples[, i, sites.0.indx])))
      }
      if (object$dist %in% c('Gaussian', 'zi-Gaussian')) {
        mu.long[, i, ] <- t(X.fix %*% t(beta.samples[, sp.indx == i]) + 
                              t(beta.star.sites.0.samples[, i, ]) + 
                            t(w.star.0.samples[, i, sites.0.indx]))
      }
      tmp[, non.missing.indx] <- mu.long[, i, ]
      out$mu.0.samples[, i, , ] <- array(tmp, dim = c(n.post, J.0, K.max.0))
    }
    K <- apply(out$mu.0.samples[1, 1, , , drop = FALSE], 3, function(a) sum(!is.na(a)))
    out$y.0.samples <- array(NA, dim(out$mu.0.samples))
    J <- ncol(object$y)
    for (i in 1:n.sp) {
      for (j in 1:J.0) {
        for (k in 1:K.max.0) {
          if (sum(is.na(out$mu.0.samples[, i, j, k])) == 0) {
            if (object$dist == 'NB') {
              out$y.0.samples[, i, j, k] <- rnbinom(n.post, kappa.samples[, i], 
              				        mu = out$mu.0.samples[, i, j, k])
            } 
            if (object$dist == 'Poisson') {
              out$y.0.samples[, i, j, k] <- rpois(n.post, out$mu.0.samples[, i, j, k])
            }
            if (object$dist == 'Gaussian') {
              out$y.0.samples[, i, j, k] <- rnorm(n.post, out$mu.0.samples[, i, j, k], 
                                                  sqrt(object$tau.sq.samples[, i]))
            }
            if (object$dist == 'zi-Gaussian') {
              out$y.0.samples[, i, j, k] <- ifelse(z.0.samples[, i, j] == 1, 
              				  rnorm(n.post, out$mu.0.samples[, i, j, k], 
              					sqrt(object$tau.sq.samples[, i])), 
              				  rnorm(n.post, 0, sqrt(0.0001)))
              out$mu.0.samples[, i, j, k] <- ifelse(z.0.samples[, i, j] == 1, 
            				   out$mu.0.samples[, i, j, k], 0)
            }
          }
        }
      }
    }
    out$w.0.samples <- w.0.samples
  }
  out$call <- cl
  # If Gaussian, collapse to a 3D array 
  if (object$dist %in% c('Gaussian', 'zi-Gaussian')) {
    out$y.0.samples <- out$y.0.samples[, , , 1]  
    out$mu.0.samples <- out$mu.0.samples[, , , 1]  
  }

  class(out) <- "predict.lfMsAbund"
  out
}
# DS ----------------------------------------------------------------------
print.DS <- function(x, ...) {
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.75)),
      "", sep = "\n")
}
summary.DS <- function(object,
                       quantiles = c(0.025, 0.5, 0.975),
                       digits = max(3L, getOption("digits") - 3L), ...) {
  summary.NMix(object, quantiles, digits)
}

predict.DS <- function(object, X.0, ignore.RE = FALSE, 
			  type = 'abundance', ...) {

  # Check for unused arguments ------------------------------------------
  formal.args <- names(formals(sys.function(sys.parent())))
  elip.args <- names(list(...))
  for(i in elip.args){
      if(! i %in% formal.args)
          warning("'",i, "' is not an argument")
  }
  # Call ----------------------------------------------------------------
  cl <- match.call()

  # Functions ---------------------------------------------------------------
  logit <- function(theta, a = 0, b = 1) {log((theta-a)/(b-theta))}
  logit.inv <- function(z, a = 0, b = 1) {b-(b-a)/(1+exp(z))}

  # Some initial checks ---------------------------------------------------
  if (missing(object)) {
    stop("error: predict expects object\n")
  }
  if (!(tolower(type) %in% c('abundance', 'detection'))) {
    stop("error: prediction type must be either 'abundance' or 'detection'")
  }

  # Check X.0 -------------------------------------------------------------
  if (missing(X.0)) {
    stop("error: X.0 must be specified\n")
  }
  if (!any(is.data.frame(X.0), is.matrix(X.0))) {
    stop("error: X.0 must be a data.frame or matrix\n")
  }

  # Abundance predictions ------------------------------------------------- 
  if (tolower(type) == 'abundance') {
    p.abund <- ncol(object$X)

    # Composition sampling --------------------------------------------------
    beta.samples <- as.matrix(object$beta.samples)
    if (object$dist == 'NB') {
      kappa.samples <- as.matrix(object$kappa.samples)
    }
    n.post <- object$n.post * object$n.chains
    out <- list()
    if (object$muRE) {
      p.abund.re <- length(object$re.level.names)
    } else {
      p.abund.re <- 0
    }
    re.cols <- object$re.cols

    if (object$muRE & !ignore.RE) {
      beta.star.samples <- object$beta.star.samples
      re.level.names <- object$re.level.names
      # Get columns in design matrix with random effects
      x.re.names <- colnames(object$X.re)
      x.0.names <- colnames(X.0)
      re.long.indx <- sapply(re.cols, length)
      tmp <- sapply(x.re.names, function(a) which(colnames(X.0) %in% a))
      indx <- list()
      for (i in 1:length(tmp)) {
        indx[[i]] <- rep(tmp[i], re.long.indx[i])
      }
      indx <- unlist(indx)
      if (length(indx) == 0) {
        stop("error: column names in X.0 must match variable names in data$occ.covs")
      }
      n.abund.re <- length(indx)
      n.unique.abund.re <- length(unique(indx))
      # Check RE columns
      for (i in 1:n.abund.re) {
        if (is.character(re.cols[[i]])) {
          # Check if all column names in svc are in occ.covs
          if (!all(re.cols[[i]] %in% x.0.names)) {
              missing.cols <- re.cols[[i]][!(re.cols[[i]] %in% x.0.names)]
              stop(paste("error: variable name ", paste(missing.cols, collapse=" and "), " not in abundance covariates", sep=""))
          }
          # Convert desired column names into the numeric column index
          re.cols[[i]] <- which(x.0.names %in% re.cols[[i]])
          
        } else if (is.numeric(re.cols[[i]])) {
          # Check if all column indices are in 1:p.abund
          if (!all(re.cols %in% 1:p.abund)) {
              missing.cols <- re.cols[[i]][!(re.cols[[i]] %in% (1:p.abund))]
              stop(paste("error: column index ", paste(missing.cols, collapse=" "), " not in design matrix columns", sep=""))
          }
        }
      }
      re.cols <- unlist(re.cols)
      X.re <- as.matrix(X.0[, indx, drop = FALSE])
      X.fix <- as.matrix(X.0[, -indx, drop = FALSE])
      X.random <- as.matrix(X.0[, re.cols, drop = FALSE])
      n.abund.re <- length(unlist(re.level.names))
      X.re.ind <- matrix(NA, nrow(X.re), p.abund.re)
      for (i in 1:p.abund.re) {
        for (j in 1:nrow(X.re)) {
          tmp <- which(re.level.names[[i]] == X.re[j, i])
          if (length(tmp) > 0) {
            X.re.ind[j, i] <- tmp 
          }
        }
      }
      if (p.abund.re > 1) {
        for (j in 2:p.abund.re) {
          X.re.ind[, j] <- X.re.ind[, j] + max(X.re.ind[, j - 1]) 
        }
      }
      # Create the random effects corresponding to each 
      # new location
      # ORDER: ordered by site, then species within site.
      beta.star.sites.0.samples <- matrix(0, n.post,  nrow(X.re))
      for (t in 1:p.abund.re) {
        for (j in 1:nrow(X.re)) {
          if (!is.na(X.re.ind[j, t])) {
            beta.star.sites.0.samples[, j] <- 
              beta.star.samples[, X.re.ind[j, t]] * X.random[j, t] + 
              beta.star.sites.0.samples[, j]
          } else {
            beta.star.sites.0.samples[, j] <- 
              rnorm(n.post, 0, sqrt(object$sigma.sq.mu.samples[, t])) * X.random[j, t] + 
              beta.star.sites.0.samples[, j]
          }
        } # j
      } # t
    } else {
      X.fix <- X.0
      beta.star.sites.0.samples <- matrix(0, n.post, nrow(X.0))
      p.abund.re <- 0
    }

    J.str <- nrow(X.0)
    out$mu.0.samples <- mcmc(exp(t(X.fix %*% t(beta.samples) + 
          				t(beta.star.sites.0.samples))))
    if (object$dist == 'NB') {
      out$N.0.samples <- mcmc(matrix(rnbinom(length(out$mu.0.samples), kappa.samples, 
					   mu = c(t(apply(out$mu.0.samples, 1, 
							  function(a) a)))), 
				     nrow = n.post, 
				     ncol = nrow(X.0)))
    } else {
      out$N.0.samples <- mcmc(matrix(rpois(length(out$mu.0.samples),  
					   lambda = c(t(apply(out$mu.0.samples, 1, 
							      function(a) a)))), 
				     nrow = n.post, 
				     ncol = nrow(X.0)))
    }
  }
  # Detection predictions -------------------------------------------------
  if (tolower(type) == 'detection') {
    p.det <- ncol(object$X.p)
    # p.design <- p.det
    # if (object$pRE & !ignore.RE) {
    #   p.design <- p.det + ncol(object$sigma.sq.p.samples)
    # }
    # if (ncol(X.0) != p.design) {
    #   stop(paste("error: X.0 must have ", p.design, " columns\n", sep = ''))
    # }
    re.det.cols <- object$re.det.cols

    # Composition sampling --------------------------------------------------
    alpha.samples <- as.matrix(object$alpha.samples)
    n.post <- object$n.post * object$n.chains
    out <- list()
    if (object$pRE) {
      p.det.re <- length(object$p.re.level.names)
    } else {
      p.det.re <- 0
    }

    if (object$pRE & !ignore.RE) {
      alpha.star.samples <- object$alpha.star.samples
      p.re.level.names <- object$p.re.level.names
      # Get columns in design matrix with random effects
      x.p.re.names <- colnames(object$X.p.re)
      x.p.0.names <- colnames(X.0)
      re.long.indx <- sapply(re.det.cols, length)
      tmp <- sapply(x.p.re.names, function(a) which(colnames(X.0) %in% a))
      indx <- list()
      for (i in 1:length(tmp)) {
        indx[[i]] <- rep(tmp[i], re.long.indx[i])
      }
      indx <- unlist(indx)
      if (length(indx) == 0) {
        stop("error: column names in X.0 must match variable names in data$det.covs")
      }
      n.det.re <- length(indx)
      n.unique.det.re <- length(unique(indx))
      # Check RE columns
      for (i in 1:n.det.re) {
        if (is.character(re.det.cols[[i]])) {
          # Check if all column names in svc are in occ.covs
          if (!all(re.det.cols[[i]] %in% x.p.0.names)) {
              missing.cols <- re.det.cols[[i]][!(re.det.cols[[i]] %in% x.p.0.names)]
              stop(paste("error: variable name ", paste(missing.cols, collapse=" and "), " not in detection covariates", sep=""))
          }
          # Convert desired column names into the numeric column index
          re.det.cols[[i]] <- which(x.p.0.names %in% re.det.cols[[i]])
          
        } else if (is.numeric(re.det.cols[[i]])) {
          # Check if all column indices are in 1:p.abund
          if (!all(re.det.cols %in% 1:p.det)) {
              missing.cols <- re.det.cols[[i]][!(re.det.cols[[i]] %in% (1:p.det))]
              stop(paste("error: column index ", paste(missing.cols, collapse=" "), " not in design matrix columns", sep=""))
          }
        }
      }
      re.det.cols <- unlist(re.det.cols)
      X.re <- as.matrix(X.0[, indx, drop = FALSE])
      X.fix <- as.matrix(X.0[, -indx, drop = FALSE])
      X.random <- as.matrix(X.0[, re.det.cols, drop = FALSE])
      n.det.re <- length(unlist(p.re.level.names))
      X.re.ind <- matrix(NA, nrow(X.re), p.det.re)
      for (i in 1:p.det.re) {
        for (j in 1:nrow(X.re)) {
          tmp <- which(p.re.level.names[[i]] == X.re[j, i])
          if (length(tmp) > 0) {
            X.re.ind[j, i] <- tmp 
          }
        }
      }
      if (p.det.re > 1) {
        for (j in 2:p.det.re) {
          X.re.ind[, j] <- X.re.ind[, j] + max(X.re.ind[, j - 1]) 
        }
      }
      # Create the random effects corresponding to each 
      # new location
      # ORDER: ordered by site, then species within site.
      alpha.star.sites.0.samples <- matrix(0, n.post,  nrow(X.re))
      for (t in 1:p.det.re) {
        for (j in 1:nrow(X.re)) {
          if (!is.na(X.re.ind[j, t])) {
            alpha.star.sites.0.samples[, j] <- 
              alpha.star.samples[, X.re.ind[j, t]] * X.random[j, t] + 
              alpha.star.sites.0.samples[, j]
          } else {
            alpha.star.sites.0.samples[, j] <- 
              rnorm(n.post, 0, sqrt(object$sigma.sq.p.samples[, t])) * X.random[j, t] + 
              alpha.star.sites.0.samples[, j]
          }
        } # j
      } # t
    } else {
      X.fix <- X.0
      alpha.star.sites.0.samples <- matrix(0, n.post, nrow(X.0))
      p.det.re <- 0
    }
    out$sigma.0.samples <- mcmc(exp(t(X.fix %*% t(alpha.samples) + 
                                t(alpha.star.sites.0.samples))))
  }
  out$call <- cl

  class(out) <- "predict.DS"
  out
}

fitted.DS <- function(object, ...) {
  out <- list()
  out$y.rep.samples <- object$y.rep.samples
  out$pi.samples <- object$pi.samples
  return(out)
}

# spDS --------------------------------------------------------------------
print.spDS <- function(x, ...) {
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.75)),
      "", sep = "\n")
}

summary.spDS <- function(object,
                       quantiles = c(0.025, 0.5, 0.975),
                       digits = max(3L, getOption("digits") - 3L), ...) {
  summary.spNMix(object, quantiles, digits)
}

fitted.spDS <- function(object, ...) {
  out <- fitted.DS(object)
  return(out)
}

predict.spDS <- function(object, X.0, coords.0, n.omp.threads = 1, 
			   verbose = TRUE, n.report = 100, 
			   ignore.RE = FALSE,  
			   type = 'abundance', include.sp = TRUE, ...) {
  if (tolower(type) == 'abundance') {
    out <- predict.spNMix(object, X.0, coords.0, n.omp.threads, 
                          verbose, n.report, 
                          ignore.RE, type, include.sp, ...)  
  } else {
    out <- predict.DS(object, X.0, ignore.RE, type)
  }
  class(out) <- 'spDS'
  return(out)
}

# msDS --------------------------------------------------------------------
print.msDS <- function(x, ...) {
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.75)),
      "", sep = "\n")
}
summary.msDS <- function(object,
			 level = 'both',
			 quantiles = c(0.025, 0.5, 0.975),
			 digits = max(3L, getOption("digits") - 3L), ...) {
  summary.msNMix(object, level, quantiles, digits)
}

fitted.msDS <- function(object, ...) {
  out <- list()
  out$y.rep.samples <- object$y.rep.samples
  out$pi.samples <- object$pi.samples
  return(out)
}

predict.msDS <- function(object, X.0, ignore.RE = FALSE, 
                         type = 'abundance', ...) {

  # Check for unused arguments ------------------------------------------
  formal.args <- names(formals(sys.function(sys.parent())))
  elip.args <- names(list(...))
  for(i in elip.args){
      if(! i %in% formal.args)
          warning("'",i, "' is not an argument")
  }
  # Call ----------------------------------------------------------------
  cl <- match.call()

  # Functions ---------------------------------------------------------------
  logit <- function(theta, a = 0, b = 1) {log((theta-a)/(b-theta))}
  logit.inv <- function(z, a = 0, b = 1) {b-(b-a)/(1+exp(z))}

  # Some initial checks ---------------------------------------------------
  if (missing(object)) {
    stop("error: predict expects object\n")
  }
  if (!(tolower(type) %in% c('abundance', 'detection'))) {
    stop("error: prediction type must be either 'abundance' or 'detection'")
  }

  # Check X.0 -------------------------------------------------------------
  if (missing(X.0)) {
    stop("error: X.0 must be specified\n")
  }
  if (!any(is.data.frame(X.0), is.matrix(X.0))) {
    stop("error: X.0 must be a data.frame or matrix\n")
  }
  # Abundance predictions ------------------------------------------------- 
  if (tolower(type) == 'abundance') {
    out <- predict.msNMix(object, X.0, ignore.RE, type)
  }
  # Detection predictions -------------------------------------------------
  if (tolower(type) == 'detection') {
    p.det <- ncol(object$X.p)
    re.det.cols <- object$re.det.cols
    # Composition sampling --------------------------------------------------
    alpha.samples <- as.matrix(object$alpha.samples)
    n.sp <- dim(object$y)[1]
    sp.indx <- rep(1:n.sp, p.det)
    n.post <- object$n.post * object$n.chains
    out <- list()
    out$sigma.0.samples <- array(NA, dim = c(n.post, n.sp, nrow(X.0)))
    if (object$pRE) {
      p.det.re <- length(object$p.re.level.names)
    } else {
      p.det.re <- 0
    }
    if (object$pRE & !ignore.RE) {
      alpha.star.samples <- object$alpha.star.samples
      p.re.level.names <- object$p.re.level.names
      # Get columns in design matrix with random effects
      x.p.re.names <- colnames(object$X.p.re)
      x.p.0.names <- colnames(X.0)
      re.long.indx <- sapply(re.det.cols, length)
      tmp <- sapply(x.p.re.names, function(a) which(colnames(X.0) %in% a))
      indx <- list()
      for (i in 1:length(tmp)) {
        indx[[i]] <- rep(tmp[i], re.long.indx[i])
      }
      indx <- unlist(indx)
      if (length(indx) == 0) {
        stop("error: column names in X.0 must match variable names in data$det.covs")
      }
      n.det.re <- length(indx)
      n.unique.det.re <- length(unique(indx))
      # Check RE columns
      for (i in 1:n.det.re) {
        if (is.character(re.det.cols[[i]])) {
          # Check if all column names in svc are in occ.covs
          if (!all(re.det.cols[[i]] %in% x.p.0.names)) {
              missing.cols <- re.det.cols[[i]][!(re.det.cols[[i]] %in% x.p.0.names)]
              stop(paste("error: variable name ", paste(missing.cols, collapse=" and "), " not in detection covariates", sep=""))
          }
          # Convert desired column names into the numeric column index
          re.det.cols[[i]] <- which(x.p.0.names %in% re.det.cols[[i]])
          
        } else if (is.numeric(re.det.cols[[i]])) {
          # Check if all column indices are in 1:p.abund
          if (!all(re.det.cols %in% 1:p.det)) {
              missing.cols <- re.det.cols[[i]][!(re.det.cols[[i]] %in% (1:p.det))]
              stop(paste("error: column index ", paste(missing.cols, collapse=" "), " not in design matrix columns", sep=""))
          }
        }
      }
      re.det.cols <- unlist(re.det.cols)
      X.re <- as.matrix(X.0[, indx, drop = FALSE])
      X.fix <- as.matrix(X.0[, -indx, drop = FALSE])
      X.random <- as.matrix(X.0[, re.det.cols, drop = FALSE])
      n.det.re <- length(unlist(p.re.level.names))
      X.re.ind <- matrix(NA, nrow(X.re), p.det.re)
      for (i in 1:p.det.re) {
        for (j in 1:nrow(X.re)) {
          tmp <- which(p.re.level.names[[i]] == X.re[j, i])
          if (length(tmp) > 0) {
            X.re.ind[j, i] <- tmp 
          }
        }
      }
      if (p.det.re > 1) {
        for (j in 2:p.det.re) {
          X.re.ind[, j] <- X.re.ind[, j] + max(X.re.ind[, j - 1]) 
        }
      }
      # Create the random effects corresponding to each 
      # new location
      # ORDER: ordered by site, then species within site.
      alpha.star.sites.0.samples <- matrix(0, n.post, n.sp * nrow(X.re))
      for (i in 1:n.sp) {
        for (t in 1:p.det.re) {
          for (j in 1:nrow(X.re)) {
            if (!is.na(X.re.ind[j, t])) {
              alpha.star.sites.0.samples[, (j - 1) * n.sp + i] <- 
                alpha.star.samples[, (i - 1) * n.det.re + X.re.ind[j, t]] * X.random[j, t] + 
                alpha.star.sites.0.samples[, (j - 1) * n.sp + i]
            } else {
              alpha.star.sites.0.samples[, (j - 1) * n.sp + i] <- 
                rnorm(n.post, 0, sqrt(object$sigma.sq.p.samples[, t])) * X.random[j, t] + 
                alpha.star.sites.0.samples[, (j - 1) * n.sp + i]
            }
          } # j
        } # t
      } # i 
    } else {
      X.fix <- X.0
      alpha.star.sites.0.samples <- matrix(0, n.post, n.sp * nrow(X.0))
      p.det.re <- 0
    }
    J.str <- nrow(X.0)
    # Make predictions
    for (i in 1:n.sp) {
      for (j in 1:J.str) {
        out$sigma.0.samples[, i, j] <- exp(t(as.matrix(X.fix[j, ])) %*% 
                                           t(alpha.samples[, sp.indx == i]) + 
                                           alpha.star.sites.0.samples[, (j - 1) * n.sp + i])
      } # j
    } # i
  }
  class(out) <- 'predict.msDS'
  out
}

# lfMsDS ----------------------------------------------------------------
print.lfMsDS <- function(x, ...) {
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.75)),
      "", sep = "\n")
}

summary.lfMsDS <- function(object,
			     level = 'both',
			     quantiles = c(0.025, 0.5, 0.975),
			     digits = max(3L, getOption("digits") - 3L), ...) {
  summary.msDS(object, level, quantiles, digits)
}

fitted.lfMsDS <- function(object, ...) {
  fitted.msDS(object)
}

predict.lfMsDS <- function(object, X.0, coords.0, ignore.RE = FALSE,
                         type = 'abundance', include.w = TRUE, ...) {

  # Abundance predictions ------------------------------------------------- 
  if (tolower(type) == 'abundance') {
    out <- predict.lfMsNMix(object, X.0, coords.0, ignore.RE, type = 'abundance', include.w)
  }
  # Detection predictions -------------------------------------------------
  if (tolower(type) == 'detection') {
    out <- predict.msDS(object, X.0, ignore.RE, type = 'detection')
  }

  class(out) <- "predict.lfMsDS"
  out
}

# sfMsDS ----------------------------------------------------------------
print.sfMsDS <- function(x, ...) {
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.75)),
      "", sep = "\n")
}

summary.sfMsDS <- function(object,
			     level = 'both',
			     quantiles = c(0.025, 0.5, 0.975),
			     digits = max(3L, getOption("digits") - 3L), ...) {
  summary.sfMsNMix(object, level, quantiles, digits)
}

fitted.sfMsDS <- function(object, ...) {
  fitted.msDS(object)
}

predict.sfMsDS <- function(object, X.0, coords.0, n.omp.threads = 1,
			   verbose = TRUE, n.report = 100,
			   ignore.RE = FALSE, type = 'abundance', 
			   include.sp = TRUE, ...) {
  if (tolower(type) == 'abundance') {
    out <- predict.sfMsNMix(object, X.0, coords.0, n.omp.threads,
                            verbose, n.report,
                            ignore.RE, type = 'abundance', include.sp, ...)
  } else {
    out <- predict.msDS(object, X.0, ignore.RE, type = 'detection')
  }
  class(out) <- 'predict.sfMsDS'
  return(out)
}

# svcAbund --------------------------------------------------------------------
print.svcAbund <- function(x, ...) {
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.75)),
      "", sep = "\n")
}

fitted.svcAbund <- function(object, ...) {
  return(object$y.rep.samples)
}

summary.svcAbund <- function(object,
			    quantiles = c(0.025, 0.5, 0.975),
			    digits = max(3L, getOption("digits") - 3L), ...) {
  summary.spAbund(object, quantiles, digits)

}

predict.svcAbund <- function(object, X.0, coords.0, n.omp.threads = 1, 
                             verbose = TRUE, n.report = 100, 
                             ignore.RE = FALSE, z.0.samples, ...) {
  # Check for unused arguments ------------------------------------------
  formal.args <- names(formals(sys.function(sys.parent())))
  elip.args <- names(list(...))
  for(i in elip.args){
      if(! i %in% formal.args)
          warning("'",i, "' is not an argument")
  }
  # Call ----------------------------------------------------------------
  cl <- match.call()

  # Functions ---------------------------------------------------------------
  logit <- function(theta, a = 0, b = 1) {log((theta-a)/(b-theta))}
  logit.inv <- function(z, a = 0, b = 1) {b-(b-a)/(1+exp(z))}

  # Some initial checks ---------------------------------------------------
  if (missing(object)) {
    stop("error: predict expects object\n")
  }

  # Check X.0 -------------------------------------------------------------
  if (missing(X.0)) {
    stop("error: X.0 must be specified\n")
  }
  if (!any(is.data.frame(X.0), is.matrix(X.0))) {
    stop("error: X.0 must be a data.frame or matrix\n")
  }

  # Abundance predictions ------------------------------------------------
  if (missing(coords.0)) {
    stop("error: coords.0 must be specified\n")
  }
  if (!any(is.data.frame(coords.0), is.matrix(coords.0))) {
    stop("error: coords.0 must be a data.frame or matrix\n")
  }
  if (!ncol(coords.0) == 2){
    stop("error: coords.0 must have two columns\n")
  }
  coords.0 <- as.matrix(coords.0)
  
  p.abund <- ncol(object$X)
  p.design <- p.abund
  X <- object$X
  if (is(object, 'spAbund') & object$dist %in% c('Gaussian', 'zi-Gaussian')) {
    svc.cols <- 1
    p.svc <- 1
    X.w.0 <- matrix(1, nrow = nrow(X.0), ncol = 1)
    X.w <- matrix(1, nrow = nrow(object$X), ncol = 1)
  } else {
    svc.cols <- object$svc.cols
    p.svc <- length(svc.cols)
    X.w.0 <- X.0[, svc.cols, drop = FALSE]
    X.w <- object$X.w
  }
  coords <- object$coords
  J <- nrow(X)
  beta.samples <- as.matrix(object$beta.samples)
  n.post <- object$n.post * object$n.chains
  family <- object$dist
  family.c <- ifelse(family == 'zi-Gaussian', 3, 2)
  theta.samples <- object$theta.samples
  w.samples <- object$w.samples
  n.neighbors <- object$n.neighbors
  cov.model.indx <- object$cov.model.indx
  re.cols <- object$re.cols
  sp.type <- object$type
  out <- list()
  if (object$muRE) {
    p.abund.re <- length(object$re.level.names)
  } else {
    p.abund.re <- 0
  }

  # Eliminate prediction sites that have already sampled been for now
  match.indx <- match(do.call("paste", as.data.frame(coords.0)), do.call("paste", as.data.frame(coords)))
  coords.0.indx <- which(is.na(match.indx))
  coords.indx <- match.indx[!is.na(match.indx)]
  coords.place.indx <- which(!is.na(match.indx))
  # coords.0.new <- coords.0[coords.0.indx, , drop = FALSE]
  # X.0.new <- X.0[coords.0.indx, , drop = FALSE]

  if (object$muRE & !ignore.RE) {
    beta.star.samples <- object$beta.star.samples
    re.level.names <- object$re.level.names
    # Get columns in design matrix with random effects
    x.re.names <- colnames(object$X.re)
    x.0.names <- colnames(X.0)
    re.long.indx <- sapply(re.cols, length)
    tmp <- sapply(x.re.names, function(a) which(colnames(X.0) %in% a))
    indx <- list()
    for (i in 1:length(tmp)) {
      indx[[i]] <- rep(tmp[i], re.long.indx[i])
    }
    indx <- unlist(indx)
    if (length(indx) == 0) {
      stop("error: column names in X.0 must match variable names in data$abund.covs")
    }
    n.abund.re <- length(indx)
    n.unique.abund.re <- length(unique(indx))
    # Check RE columns
    for (i in 1:n.abund.re) {
      if (is.character(re.cols[[i]])) {
        # Check if all column names in svc are in occ.covs
        if (!all(re.cols[[i]] %in% x.0.names)) {
            missing.cols <- re.cols[[i]][!(re.cols[[i]] %in% x.0.names)]
            stop(paste("error: variable name ", paste(missing.cols, collapse=" and "), " not in abundance covariates", sep=""))
        }
        # Convert desired column names into the numeric column index
        re.cols[[i]] <- which(x.0.names %in% re.cols[[i]])
        
      } else if (is.numeric(re.cols[[i]])) {
        # Check if all column indices are in 1:p.abund
        if (!all(re.cols %in% 1:p.abund)) {
            missing.cols <- re.cols[[i]][!(re.cols[[i]] %in% (1:p.abund))]
            stop(paste("error: column index ", paste(missing.cols, collapse=" "), " not in design matrix columns", sep=""))
        }
      }
    }
    re.cols <- unlist(re.cols)
    X.re <- as.matrix(X.0[, indx, drop = FALSE])
    X.fix <- as.matrix(X.0[, -indx, drop = FALSE])
    X.random <- as.matrix(X.0[, re.cols, drop = FALSE])
    n.abund.re <- length(unlist(re.level.names))
    X.re.ind <- matrix(NA, nrow(X.re), p.abund.re)
    for (i in 1:p.abund.re) {
      for (j in 1:nrow(X.re)) {
        tmp <- which(re.level.names[[i]] == X.re[j, i])
        if (length(tmp) > 0) {
          X.re.ind[j, i] <- tmp 
        }
      }
    }
    if (p.abund.re > 1) {
      for (j in 2:p.abund.re) {
        X.re.ind[, j] <- X.re.ind[, j] + max(X.re.ind[, j - 1]) 
      }
    }
    # Create the random effects corresponding to each 
    # new location
    # ORDER: ordered by site, then species within site.
    beta.star.sites.0.samples <- matrix(0, n.post,  nrow(X.re))
    for (t in 1:p.abund.re) {
      for (j in 1:nrow(X.re)) {
        if (!is.na(X.re.ind[j, t])) {
          beta.star.sites.0.samples[, j] <- 
            beta.star.samples[, X.re.ind[j, t]] * X.random[j, t] + 
            beta.star.sites.0.samples[, j]
        } else {
          beta.star.sites.0.samples[, j] <- 
            rnorm(n.post, 0, sqrt(object$sigma.sq.mu.samples[, t])) * X.random[j, t] + 
            beta.star.sites.0.samples[, j]
        }
      } # j
    } # t
  } else {
    X.fix <- X.0
    beta.star.sites.0.samples <- matrix(0, n.post, nrow(X.0))
    p.abund.re <- 0
  }
  # Get samples in proper format for C++
  beta.samples <- t(beta.samples)
  w.samples <- matrix(w.samples, n.post, J * p.svc)
  # Order: iteration, site within iteration, svc within site. 
  # Example: site 1, svc 1, iter 1, site 1, svc 2, iter 1, ..., site 2, svc 1, iter 1
  w.samples <- t(w.samples)
  beta.star.sites.0.samples <- t(beta.star.sites.0.samples)
  if (family %in% c('Gaussian', 'zi-Gaussian')) {
    tau.sq.samples <- t(object$tau.sq.samples)
  } else {
    tau.sq.samples <- matrix(0, 1, n.post)
  }
  theta.samples <- t(theta.samples)

  sites.0.indx <- 0:(nrow(X.0) - 1)
  J.0 <- length(unique(sites.0.indx))
  sites.0.sampled <- ifelse(!is.na(match.indx), 1, 0)
  sites.link <- rep(NA, J.0)
  sites.link[which(!is.na(match.indx))] <- coords.indx
  # For C
  sites.link <- sites.link - 1

  # Check stage 1 samples
  if (family == 'zi-Gaussian') {
    if (missing(z.0.samples)) {
      stop("z.0.samples must be supplied for a zi-Gaussian model")
    }
    if (!is.matrix(z.0.samples)) {
      stop(paste("z.0.samples must be a matrix with ", n.post, " rows and ", 
		 J.0, " columns.", sep = ''))
    }
    if (nrow(z.0.samples) != n.post | ncol(z.0.samples) != J.0) {
      stop(paste("z.0.samples must be a matrix with ", n.post, " rows and ", 
		 J.0, " columns.", sep = ''))
    }
  } else {
    if (!missing(z.0.samples)) {
      message("z.0.samples is ignored for the current model family\n")
    }
    z.0.samples <- NA
  }
  z.0.samples <- t(z.0.samples)

  if (sp.type == 'GP') {
    stop("NNGP = FALSE is not currently supported for svcAbund")
  } else {
    # Get nearest neighbors
    # nn2 is a function from RANN.
    nn.indx.0 <- nn2(coords, coords.0, k=n.neighbors)$nn.idx-1

    storage.mode(coords) <- "double"
    storage.mode(J) <- "integer"
    storage.mode(p.abund) <- "integer"
    storage.mode(p.svc) <- 'integer'
    storage.mode(n.neighbors) <- "integer"
    storage.mode(X.fix) <- "double"
    storage.mode(X.w.0) <- 'double'
    storage.mode(coords.0) <- "double"
    storage.mode(J.0) <- "integer"
    storage.mode(sites.link) <- "integer"
    storage.mode(sites.0.sampled) <- 'integer'
    storage.mode(sites.0.indx) <- 'integer'
    storage.mode(beta.samples) <- "double"
    storage.mode(theta.samples) <- "double"
    storage.mode(tau.sq.samples) <- "double"
    storage.mode(family.c) <- "integer"
    storage.mode(w.samples) <- "double"
    storage.mode(beta.star.sites.0.samples) <- "double"
    storage.mode(n.post) <- "integer"
    storage.mode(cov.model.indx) <- "integer"
    storage.mode(nn.indx.0) <- "integer"
    storage.mode(n.omp.threads) <- "integer"
    storage.mode(verbose) <- "integer"
    storage.mode(n.report) <- "integer"
    storage.mode(z.0.samples) <- "double"

    ptm <- proc.time()

    out <- .Call("svcAbundNNGPPredict", coords, J, family.c, p.abund, p.svc, n.neighbors,
                 X.fix, X.w.0, coords.0, J.0, nn.indx.0, beta.samples,
                 theta.samples, tau.sq.samples, w.samples, beta.star.sites.0.samples,
      	         sites.link, sites.0.sampled, sites.0.indx,
        	 n.post, cov.model.indx, n.omp.threads, verbose, n.report, z.0.samples)
  }
  out$y.0.samples <- mcmc(t(out$y.0.samples))
  out$mu.0.samples <- mcmc(t(out$mu.0.samples))
  out$w.0.samples <- array(out$w.0.samples, dim = c(p.svc, J.0, n.post))
  out$w.0.samples <- aperm(out$w.0.samples, c(3, 1, 2))
  if (is(object, 'spAbund') & object$dist %in% c('Gaussian', 'zi-Gaussian')) {
    out$w.0.samples <- matrix(out$w.0.samples[, 1, ], n.post, J.0)
  }
  out$call <- cl

  class(out) <- "predict.spNMix"
  out

}

# svcMsAbund --------------------------------------------------------------
print.svcMsAbund <- function(x, ...) {
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.75)),
      "", sep = "\n")
}
summary.svcMsAbund <- function(object, level = 'both',
                              quantiles = c(0.025, 0.5, 0.975),
                              digits = max(3L, getOption("digits") - 3L), ...) {
  summary.sfMsAbund(object, level = 'both', 
		    quantiles = c(0.025, 0.5, 0.975), 
		    digits = max(3L, getOption("digits") - 3L))
}
fitted.svcMsAbund <- function(object, ...) {
  return(object$y.rep.samples)
}
predict.svcMsAbund <- function(object, X.0, coords.0, n.omp.threads = 1,
			      verbose = TRUE, n.report = 100,
			      ignore.RE = FALSE, z.0.samples, ...) {

  ptm <- proc.time()
  # Check for unused arguments ------------------------------------------
  formal.args <- names(formals(sys.function(sys.parent())))
  elip.args <- names(list(...))
  for(i in elip.args){
      if(! i %in% formal.args)
          warning("'",i, "' is not an argument")
  }
  # Call ----------------------------------------------------------------
  cl <- match.call()

  # Functions ---------------------------------------------------------------
  logit <- function(theta, a = 0, b = 1) {log((theta-a)/(b-theta))}
  logit.inv <- function(z, a = 0, b = 1) {b-(b-a)/(1+exp(z))}

  # Some initial checks ---------------------------------------------------
  if (missing(object)) {
    stop("error: predict expects object\n")
  }
  if (!(class(object) %in% c('sfMsAbund', 'svcMsAbund'))) {
    stop("error: requires an output object of class sfMsAbund or svcMsAbund\n")
  }

  # Check X.0 -------------------------------------------------------------
  if (missing(X.0)) {
    stop("error: X.0 must be specified\n")
  }
  if (!any(is.data.frame(X.0), is.matrix(X.0))) {
    stop("error: X.0 must be a data.frame or matrix\n")
  }

  if (missing(coords.0)) {
    stop("error: coords.0 must be specified\n")
  }
  if (!any(is.data.frame(coords.0), is.matrix(coords.0))) {
    stop("error: coords.0 must be a data.frame or matrix\n")
  }
  if (!ncol(coords.0) == 2){
    stop("error: coords.0 must have two columns\n")
  }
  coords.0 <- as.matrix(coords.0)

  p.abund <- ncol(object$X)
  p.design <- p.abund
  X <- object$X
  if (is(object, 'sfMsAbund') & object$dist %in% c('Gaussian', 'zi-Gaussian')) {
    svc.cols <- 1
    p.svc <- 1
    X.w.0 <- matrix(1, nrow = nrow(X.0), ncol = 1)
    X.w <- matrix(1, nrow = nrow(object$X), ncol = 1)
  } else {
    svc.cols <- object$svc.cols
    p.svc <- length(svc.cols)
    X.w.0 <- X.0[, svc.cols, drop = FALSE]
    X.w <- object$X.w
  }

  # Abundance predictions ------------------------------------------------
  n.post <- object$n.post * object$n.chains
  X <- object$X
  y <- object$y
  coords <- object$coords
  J <- nrow(coords)
  n.sp <- dim(y)[1]
  q <- object$q
  theta.samples <- object$theta.samples
  beta.samples <- object$beta.samples
  lambda.samples <- object$lambda.samples
  n.neighbors <- object$n.neighbors
  cov.model.indx <- object$cov.model.indx
  family <- object$dist
  family.c <- ifelse(family == 'zi-Gaussian', 3, 2)
  re.cols <- object$re.cols
  sp.type <- object$type
  X.0 <- as.matrix(X.0)
  if (object$muRE & !ignore.RE) {
    p.abund.re <- length(object$re.level.names)
  } else {
    p.abund.re <- 0
  }
  re.cols <- object$re.cols

  if (object$muRE & !ignore.RE) {
    beta.star.samples <- object$beta.star.samples
    re.level.names <- object$re.level.names
    # Get columns in design matrix with random effects
    x.re.names <- dimnames(object$X.re)[[2]]
    x.0.names <- colnames(X.0)
    # Get the number of times each factor is used. 
    re.long.indx <- sapply(re.cols, length)
    tmp <- sapply(x.re.names, function(a) which(colnames(X.0) %in% a))
    indx <- list()
    for (i in 1:length(tmp)) {
      indx[[i]] <- rep(tmp[i], re.long.indx[i])
    }
    indx <- unlist(indx)
    if (length(indx) == 0) {
      stop("error: column names in X.0 must match variable names in data$occ.covs")
    }
    n.re <- length(indx)
    n.unique.re <- length(unique(indx))
    # Check RE columns
    for (i in 1:n.re) {
      if (is.character(re.cols[[i]])) {
        # Check if all column names in svc are in occ.covs
        if (!all(re.cols[[i]] %in% x.0.names)) {
            missing.cols <- re.cols[[i]][!(re.cols[[i]] %in% x.0.names)]
            stop(paste("error: variable name ", paste(missing.cols, collapse=" and "), " not in occurrence covariates", sep=""))
        }
        # Convert desired column names into the numeric column index
        re.cols[[i]] <- which(x.0.names %in% re.cols[[i]])
        
      } else if (is.numeric(re.cols[[i]])) {
        # Check if all column indices are in 1:p.abund
        if (!all(re.cols %in% 1:p.abund)) {
            missing.cols <- re.cols[[i]][!(re.cols[[i]] %in% (1:p.abund))]
            stop(paste("error: column index ", paste(missing.cols, collapse=" "), " not in design matrix columns", sep=""))
        }
      }
    }
    re.cols <- unlist(re.cols)
    X.re <- as.matrix(X.0[, indx, drop = FALSE])
    X.fix <- as.matrix(X.0[, -indx, drop = FALSE])
    X.random <- as.matrix(X.0[, re.cols, drop = FALSE])
    n.abund.re <- length(unlist(re.level.names))
    X.re.ind <- matrix(NA, nrow(X.re), p.abund.re)

    for (i in 1:p.abund.re) {
      for (j in 1:nrow(X.re)) {
        tmp <- which(re.level.names[[i]] == X.re[j, i])
        if (length(tmp) > 0) {
          X.re.ind[j, i] <- tmp 
        }
      }
    }
    if (p.abund.re > 1) {
      for (j in 2:p.abund.re) {
        X.re.ind[, j] <- X.re.ind[, j] + max(X.re.ind[, j - 1]) 
      }
    }
    # Create the random effects corresponding to each 
    # new location
    beta.star.sites.0.samples <- array(0, dim = c(n.post, n.sp, nrow(X.re)))
    for (i in 1:n.sp) {
      for (t in 1:p.abund.re) {
        for (j in 1:nrow(X.re)) {
          if (!is.na(X.re.ind[j, t])) {
            beta.star.sites.0.samples[, i, j] <- 
              beta.star.samples[, (i - 1) * n.abund.re + X.re.ind[j, t]] * X.random[j, t] + 
              beta.star.sites.0.samples[, i, j]
          } else {
            beta.star.sites.0.samples[, i, j] <- 
              rnorm(n.post, 0, sqrt(object$sigma.sq.mu.samples[, t])) * X.random[j, t] + 
              beta.star.sites.0.samples[, i, j]
          }
        } # j
      } # t
    } # i
  } else {
    X.fix <- X.0
    beta.star.sites.0.samples <- array(0, dim = c(n.post, n.sp, nrow(X.0)))
    p.abund.re <- 0
  }

  # Sub-sample previous
  theta.samples <- t(theta.samples)
  lambda.samples <- t(lambda.samples)
  beta.samples <- t(beta.samples)
  if (family == 'NB') {
    kappa.samples <- t(object$kappa.samples)
  } else {
    kappa.samples <- matrix(0, n.sp, n.post)
  }
  if (family %in% c('Gaussian', 'zi-Gaussian')) {
    tau.sq.samples <- t(object$tau.sq.samples)
  } else {
    tau.sq.samples <- matrix(0, 1, n.post)
  }
  w.samples <- object$w.samples
  if (length(dim(w.samples)) == 3) {
    w.samples <- aperm(w.samples, c(2, 3, 1))
  } else {
    w.samples <- aperm(w.samples, c(2, 3, 4, 1))
  }
  beta.star.sites.0.samples <- aperm(beta.star.sites.0.samples, c(2, 3, 1))

  J.0 <- nrow(X.fix)

  match.indx <- match(do.call("paste", as.data.frame(coords.0)), do.call("paste", as.data.frame(coords)))
  coords.0.indx <- which(is.na(match.indx))
  coords.indx <- match.indx[!is.na(match.indx)]
  coords.place.indx <- which(!is.na(match.indx))
  # Indicates whether a site has been sampled. 1 = sampled
  sites.0.sampled <- ifelse(!is.na(match.indx), 1, 0)
  sites.link <- rep(NA, J.0)
  sites.link[which(!is.na(match.indx))] <- coords.indx
  # For C
  sites.link <- sites.link - 1

  if (family == 'zi-Gaussian') {
    if (missing(z.0.samples)) {
      stop("z.0.samples must be supplied for a zi-Gaussian model")
    }
    if (length(dim(z.0.samples)) != 3) {
      stop(paste("z.0.samples must be a three-dimensional array with dimensions of ", 
		 n.post, ", ", n.sp, ", and ", J.0, sep = ''))
    }
    if (dim(z.0.samples)[1] != n.post | dim(z.0.samples)[2] != n.sp |
	dim(z.0.samples)[3] != J.0) {
      stop(paste("z.0.samples must be a three-dimensional array with dimensions of ", 
		 n.post, ", ", n.sp, ", and ", J.0, sep = ''))
    }
  } else {
    if (!missing(z.0.samples)) {
      message("z.0.samples is ignored for the current model family\n")
    }
    z.0.samples <- array(NA, dim = c(1, 1, 1))
  }
  z.0.samples <- aperm(z.0.samples, c(2, 3, 1))

  if (sp.type == 'GP') {
    # Not currently implemented or accessed.
  } else {
    # Get nearest neighbors
    # nn2 is a function from RANN.
    nn.indx.0 <- nn2(coords, coords.0, k=n.neighbors)$nn.idx-1

    storage.mode(coords) <- "double"
    storage.mode(n.sp) <- "integer"
    storage.mode(J) <- "integer"
    storage.mode(p.abund) <- "integer"
    storage.mode(p.svc) <- 'integer'
    storage.mode(n.neighbors) <- "integer"
    storage.mode(X.fix) <- "double"
    storage.mode(X.w.0) <- 'double'
    storage.mode(coords.0) <- "double"
    storage.mode(J.0) <- "integer"
    storage.mode(q) <- "integer"
    storage.mode(sites.link) <- "integer"
    storage.mode(sites.0.sampled) <- "integer"
    storage.mode(beta.samples) <- "double"
    storage.mode(theta.samples) <- "double"
    storage.mode(lambda.samples) <- "double"
    storage.mode(tau.sq.samples) <- "double"
    storage.mode(beta.star.sites.0.samples) <- "double"
    storage.mode(family.c) <- 'integer'
    storage.mode(w.samples) <- "double"
    storage.mode(n.post) <- "integer"
    storage.mode(cov.model.indx) <- "integer"
    storage.mode(nn.indx.0) <- "integer"
    storage.mode(n.omp.threads) <- "integer"
    storage.mode(verbose) <- "integer"
    storage.mode(n.report) <- "integer"
    storage.mode(family.c) <- "integer"
    storage.mode(z.0.samples) <- 'double'

    out <- .Call("svcMsAbundGaussianNNGPPredict", coords, J, family.c, 
      	         n.sp, q, p.abund, p.svc, n.neighbors,
                 X.fix, X.w.0, coords.0, J.0, sites.link, sites.0.sampled, 
		 nn.indx.0, beta.samples,
                 theta.samples, tau.sq.samples, lambda.samples, w.samples,
        	 beta.star.sites.0.samples, n.post,
                 cov.model.indx, n.omp.threads, verbose, n.report, z.0.samples)
  }

  out$y.0.samples <- array(out$y.0.samples, dim = c(n.sp, J.0, n.post))
  out$y.0.samples <- aperm(out$y.0.samples, c(3, 1, 2))
  out$w.0.samples <- array(out$w.0.samples, dim = c(q, J.0, p.svc, n.post))
  out$w.0.samples <- aperm(out$w.0.samples, c(4, 1, 2, 3))
  out$mu.0.samples <- array(out$mu.0.samples, dim = c(n.sp, J.0, n.post))
  out$mu.0.samples <- aperm(out$mu.0.samples, c(3, 1, 2))

  out$run.time <- proc.time() - ptm
  out$call <- cl
  out$object.class <- class(object)
  class(out) <- "predict.svcMsAbund"
  out
}
