# NMix --------------------------------------------------------------------
print.NMix <- function(x, ...) {
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.75)), 
      "", sep = "\n")
}

fitted.NMix <- function(object, ...) {
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
  n.post <- object$n.post * object$n.chains
  X.p <- object$X.p
  y <- object$y
  n.rep <- apply(y, 1, function(a) sum(!is.na(a)))
  K.max <- max(n.rep)
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
  N.samples <- object$N.samples
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
  if (!is(object, "NMix")) {
  # if (class(object) != "NMix") {
    stop("error: requires an output object of class NMix\n")
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
    p.design <- p.abund
    if (object$muRE & !ignore.RE) {
      p.design <- p.abund + ncol(object$sigma.sq.mu.samples)
    }
    if (ncol(X.0) != p.design) {
      stop(paste("error: X.0 must have ", p.design, " columns\n", sep = ''))
    }

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
      tmp <- which(colnames(X.0) %in% x.re.names)
      indx <- list()
      for (i in 1:length(tmp)) {
        indx[[i]] <- rep(tmp[i], re.long.indx[i])
      }
      if (length(indx) == 0) {
        stop("error: column names in X.0 must match variable names in data$occ.covs")
      }
      indx <- unlist(indx)
      n.abund.re <- length(unlist(re.level.names))
      n.unique.abund.re <- length(unique(indx))
      # Check RE columns
      for (i in 1:n.unique.abund.re) {
        if (is.character(re.cols[[i]])) {
          # Check if all column names in svc are in occ.covs
          if (!all(re.cols[[i]] %in% x.0.names)) {
              missing.cols <- re.cols[[i]][!(re.cols[[i]] %in% x.0.names)]
              stop(paste("error: variable name ", paste(missing.cols, collapse=" and "), " not in abundance covariates", sep=""))
          }
          # Convert desired column names into the numeric column index
          re.cols[[i]] <- (1:p.abund)[x.0.names %in% re.cols[[i]]]
          
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
              rnorm(n.post, 0, sqrt(object$sigma.sq.mu.samples[, t])) + 
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
    p.design <- p.det
    if (object$pRE & !ignore.RE) {
      p.design <- p.det + ncol(object$sigma.sq.p.samples)
    }
    if (ncol(X.0) != p.design) {
      stop(paste("error: X.0 must have ", p.design, " columns\n", sep = ''))
    }
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
      tmp <- which(colnames(X.0) %in% x.p.re.names)
      indx <- list()
      for (i in 1:length(tmp)) {
        indx[[i]] <- rep(tmp[i], re.long.indx[i])
      }
      if (length(indx) == 0) {
        stop("error: column names in X.0 must match variable names in data$det.covs")
      }
      indx <- unlist(indx)
      n.det.re <- length(unlist(p.re.level.names))
      n.unique.det.re <- length(unique(indx))
      # Check RE columns
      for (i in 1:n.unique.det.re) {
        if (is.character(re.det.cols[[i]])) {
          # Check if all column names in svc are in occ.covs
          if (!all(re.det.cols[[i]] %in% x.p.0.names)) {
              missing.cols <- re.det.cols[[i]][!(re.det.cols[[i]] %in% x.p.0.names)]
              stop(paste("error: variable name ", paste(missing.cols, collapse=" and "), " not in detection covariates", sep=""))
          }
          # Convert desired column names into the numeric column index
          re.det.cols[[i]] <- (1:p.det)[x.p.0.names %in% re.det.cols[[i]]]
          
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
              rnorm(n.post, 0, sqrt(object$sigma.sq.p.samples[, t])) + 
              alpha.star.sites.0.samples[, j]
          }
        } # j
      } # t
    } else {
      X.fix <- X.0
      alpha.star.sites.0.samples <- matrix(0, n.post, nrow(X.0))
      p.et.re <- 0
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
  cat("Detection (logit scale): \n")
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
    cat("Detection Random Effect Variances (logit scale): \n")
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
  cat("Detection (logit scale): \n")
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
    cat("Detection Random Effect Variances (logit scale): \n")
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

fitted.spNMix <- function(object, ...) {
  fitted.NMix(object)
}


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

  if (object$class %in% c('NMix', 'spNMix', 'abund', 'spAbund')) {
    cat("Bayesian p-value: ", round(mean(object$fit.y.rep > object$fit.y), digits), "\n")
    cat("Fit statistic: ", object$fit.stat, "\n")
  }

  # if (object$class %in% c('msPGOcc', 'spMsPGOcc', 'lfMsPGOcc', 'sfMsNMix')) {

  #   if (tolower(level) == 'community') {
  #     cat("----------------------------------------\n");
  #     cat("\tCommunity Level\n");
  #     cat("----------------------------------------\n");
  #     cat("Bayesian p-value: ", round(mean(object$fit.y.rep > object$fit.y), digits), "\n")
  #     cat("Fit statistic: ", object$fit.stat, "\n")
  #   }

  #   if (tolower(level) == 'species') {
  #     cat("----------------------------------------\n");
  #     cat("\tSpecies Level\n");
  #     cat("----------------------------------------\n");
  #     N <- ncol(object$fit.y)
  #     for (i in 1:N) {
  #       cat(paste(object$sp.names[i], " Bayesian p-value: ",
  #       	  round(mean(object$fit.y.rep[, i] > object$fit.y[, i]), digits), "\n", sep = ''))
  #     }
  #     cat("Fit statistic: ", object$fit.stat, "\n")
  #   }

  #   if (tolower(level) == 'both') {
  #     cat("----------------------------------------\n");
  #     cat("\tCommunity Level\n");
  #     cat("----------------------------------------\n");
  #     cat("Bayesian p-value: ", round(mean(object$fit.y.rep > object$fit.y), digits), "\n")
  #     cat("\n")
  #     cat("----------------------------------------\n");
  #     cat("\tSpecies Level\n");
  #     cat("----------------------------------------\n");
  #     N <- ncol(object$fit.y)
  #     for (i in 1:N) {
  #       cat(paste(object$sp.names[i], " Bayesian p-value: ",
  #       	  round(mean(object$fit.y.rep[, i] > object$fit.y[, i]), digits), "\n", sep = ''))
  #     }
  #     cat("Fit statistic: ", object$fit.stat, "\n")
  #   }
  # }

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
    cat("Detection Means (logit scale): \n")
    tmp.1 <- t(apply(object$alpha.comm.samples, 2,
          	   function(x) c(mean(x), sd(x))))
    colnames(tmp.1) <- c("Mean", "SD")
    tmp <- t(apply(object$alpha.comm.samples, 2,
          	 function(x) quantile(x, prob = quantiles)))
    diags <- matrix(c(object$rhat$alpha.comm, round(object$ESS$alpha.comm, 0)), ncol = 2)
    colnames(diags) <- c('Rhat', 'ESS')

    print(noquote(round(cbind(tmp.1, tmp, diags), digits)))
    cat("\nDetection Variances (logit scale): \n")
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
      cat("Detection Random Effect Variances (logit scale): \n")
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
    cat("Detection (logit scale): \n")
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


fitted.msNMix <- function(object, ...) {
  # TODO: need to update
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
    cat("Detection Means (logit scale): \n")
    tmp.1 <- t(apply(object$alpha.comm.samples, 2,
          	   function(x) c(mean(x), sd(x))))
    colnames(tmp.1) <- c("Mean", "SD")
    tmp <- t(apply(object$alpha.comm.samples, 2,
          	 function(x) quantile(x, prob = quantiles)))
    diags <- matrix(c(object$rhat$alpha.comm, round(object$ESS$alpha.comm, 0)), ncol = 2)
    colnames(diags) <- c('Rhat', 'ESS')

    print(noquote(round(cbind(tmp.1, tmp, diags), digits)))
    cat("\nDetection Variances (logit scale): \n")
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
      cat("Detection Random Effect Variances (logit scale): \n")
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
    cat("Detection (logit scale): \n")
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

fitted.sfMsNMix <- function(object, ...) {
  fitted.msNMix(object)
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

predict.abund <- function(object, X.0, ignore.RE = FALSE, ...) {
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
  if (!is(object, "abund")) {
    stop("error: requires an output object of class abund\n")
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
  p.design <- p.abund
  if (object$muRE & !ignore.RE) {
    p.design <- p.abund + ncol(object$sigma.sq.mu.samples)
  }
  if (dim(X.0)[3] != p.design) {
    stop(paste("error: X.0 must have ", p.design, " covariates\n", sep = ''))
  }
  J.0 <- nrow(X.0)
  K.max.0 <- ncol(X.0)

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

  # TODO: the function currently doesn't handle unbalanced replicates. 
  # Get X.0 in long format. 
  tmp.names <- dimnames(X.0)[[3]]
  X.0 <- matrix(X.0, nrow = nrow(X.0) * ncol(X.0), ncol = dim(X.0)[3])
  colnames(X.0) <- tmp.names

  if (object$muRE & !ignore.RE) {
    beta.star.samples <- object$beta.star.samples
    re.level.names <- object$re.level.names
    # Get columns in design matrix with random effects
    x.re.names <- dimnames(object$X.re)[[3]]
    x.0.names <- colnames(X.0)
    # Get the number of times each factor is used. 
    re.long.indx <- sapply(re.cols, length)
    tmp <- which(colnames(X.0) %in% x.re.names)
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
    for (i in 1:n.unique.re) {
      if (is.character(re.cols[[i]])) {
        # Check if all column names in svc are in occ.covs
        if (!all(re.cols[[i]] %in% x.0.names)) {
            missing.cols <- re.cols[[i]][!(re.cols[[i]] %in% x.0.names)]
            stop(paste("error: variable name ", paste(missing.cols, collapse=" and "), " not in occurrence covariates", sep=""))
        }
        # Convert desired column names into the numeric column index
        re.cols[[i]] <- (1:p.abund)[x.0.names %in% re.cols[[i]]]
        
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
            rnorm(n.post, 0, sqrt(object$sigma.sq.mu.samples[, t])) + 
            beta.star.sites.0.samples[, j]
        }
      } # j
    } # t
  } else {
    X.fix <- X.0
    beta.star.sites.0.samples <- matrix(0, n.post, nrow(X.0))
    p.abund.re <- 0
  }

  out$mu.0.samples <- array(exp(t(X.fix %*% t(beta.samples) + 
        				t(beta.star.sites.0.samples))), 
			    dim = c(n.post, J.0, K.max.0))
  out$y.0.samples <- array(NA, dim(out$mu.0.samples))
  for (j in 1:J.0) {
    for (k in 1:K.max.0) {
      if (object$dist == 'NB') {
        out$y.0.samples[, j, k] <- rnbinom(n.post, kappa.samples, 
					   mu = out$mu.0.samples[, j, k])
      } else {
        out$y.0.samples[, j, k] <- rpois(n.post, out$mu.0.samples[, j, k])
      }
    }
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

predict.spAbund <- function(object, X.0, coords.0,  
			    n.omp.threads = 1, verbose = TRUE, n.report = 100,
			    ignore.RE = FALSE, ...) {

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

  # TODO: the function currently doesn't handle unbalanced replicates. 
  # Get X.0 in long format. 
  tmp.names <- dimnames(X.0)[[3]]
  X.0 <- matrix(X.0, nrow = nrow(X.0) * ncol(X.0), ncol = dim(X.0)[3])
  colnames(X.0) <- tmp.names

  n.post <- object$n.post * object$n.chains
  X <- object$X
  coords <- object$coords
  n.obs <- sum(!is.na(object$y))
  J <- nrow(coords)
  n.years.max <- dim(X.0)[2]
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
  # TODO: also note you got rid of all the "new", so keep that in mind if bugs arise.  
  # coords.0.new <- coords.0[coords.0.indx, , drop = FALSE]
  # X.0.new <- X.0[coords.0.indx, , , drop = FALSE]


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
        re.cols[[i]] <- (1:p.abund)[x.0.names %in% re.cols[[i]]]
        
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
            rnorm(n.post, 0, sqrt(object$sigma.sq.mu.samples[, t])) + 
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

  # TODO: need to make sure the ordering here is right. 
  out$y.0.samples <- array(mcmc(t(out$y.0.samples)), dim = c(n.post, J.0,
							     nrow(out$y.0.samples) / J.0))
  out$mu.0.samples <- array(mcmc(t(out$mu.0.samples)), dim = c(n.post, J.0,
							     nrow(out$mu.0.samples) / J.0))
  out$w.0.samples <- mcmc(t(out$w.0.samples))
  out$run.time <- proc.time() - ptm
  out$call <- cl
  out$object.class <- class(object)
  class(out) <- "predict.spAbund"
  out
}

# spNMix ------------------------------------------------------------------
predict.spNMix <- function(object, X.0, coords.0, n.omp.threads = 1, 
			   verbose = TRUE, ignore.RE = FALSE, type = 'abundance', ...) {
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
  if (!is(object, "spNMix")) {
  # if (class(object) != "NMix") {
    stop("error: requires an output object of class spNMix\n")
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
    if (object$muRE & !ignore.RE) {
      p.design <- p.abund + ncol(object$sigma.sq.mu.samples)
    }
    if (ncol(X.0) != p.design) {
      stop(paste("error: X.0 must have ", p.design, " columns\n", sep = ''))
    }
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
    coords.0.new <- coords.0[coords.0.indx, , drop = FALSE]
    X.0.new <- X.0[coords.0.indx, , drop = FALSE]

    if (length(coords.indx) == nrow(X.0)) {
      stop("error: no new locations to predict at. See object$psi.samples for occurrence probabilities at sampled sites.")
    }

    if (object$muRE & !ignore.RE) {
      beta.star.samples <- object$beta.star.samples
      re.level.names <- object$re.level.names
      # Get columns in design matrix with random effects
      x.re.names <- colnames(object$X.re)
      x.0.names <- colnames(X.0.new)
      re.long.indx <- sapply(re.cols, length)
      tmp <- which(colnames(X.0.new) %in% x.re.names)
      indx <- list()
      for (i in 1:length(tmp)) {
        indx[[i]] <- rep(tmp[i], re.long.indx[i])
      }
      if (length(indx) == 0) {
        stop("error: column names in X.0 must match variable names in data$occ.covs")
      }
      indx <- unlist(indx)
      n.abund.re <- length(unlist(re.level.names))
      n.unique.abund.re <- length(unique(indx))
      # Check RE columns
      for (i in 1:n.unique.abund.re) {
        if (is.character(re.cols[[i]])) {
          # Check if all column names in svc are in occ.covs
          if (!all(re.cols[[i]] %in% x.0.names)) {
              missing.cols <- re.cols[[i]][!(re.cols[[i]] %in% x.0.names)]
              stop(paste("error: variable name ", paste(missing.cols, collapse=" and "), " not in abundance covariates", sep=""))
          }
          # Convert desired column names into the numeric column index
          re.cols[[i]] <- (1:p.abund)[x.0.names %in% re.cols[[i]]]
          
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
      beta.star.sites.0.samples <- matrix(0, n.post,  nrow(X.re))
      for (t in 1:p.abund.re) {
        for (j in 1:nrow(X.re)) {
          if (!is.na(X.re.ind[j, t])) {
            beta.star.sites.0.samples[, j] <- 
              beta.star.samples[, X.re.ind[j, t]] * X.random[j, t] + 
              beta.star.sites.0.samples[, j]
          } else {
            beta.star.sites.0.samples[, j] <- 
              rnorm(n.post, 0, sqrt(object$sigma.sq.mu.samples[, t])) + 
              beta.star.sites.0.samples[, j]
          }
        } # j
      } # t
    } else {
      X.fix <- X.0.new
      beta.star.sites.0.samples <- matrix(0, n.post, nrow(X.0.new))
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
    J.0 <- nrow(X.fix)


    if (sp.type == 'GP') {
      stop("NNGP = FALSE is not currently supported for spAbund")
    } else {
      # Get nearest neighbors
      # nn2 is a function from RANN.
      nn.indx.0 <- nn2(coords, coords.0.new, k=n.neighbors)$nn.idx-1

      storage.mode(coords) <- "double"
      storage.mode(J) <- "integer"
      storage.mode(p.abund) <- "integer"
      storage.mode(n.neighbors) <- "integer"
      storage.mode(X.fix) <- "double"
      storage.mode(coords.0.new) <- "double"
      storage.mode(J.0) <- "integer"
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
          	   n.post, cov.model.indx, n.omp.threads, verbose, n.report)
    }
    if (nrow(X.0) == J.0) {
      out$N.0.samples <- mcmc(t(out$N.0.samples))
      out$mu.0.samples <- mcmc(t(out$mu.0.samples))
      out$w.0.samples <- mcmc(t(out$w.0.samples))
    } else {
      tmp <- matrix(NA, n.post, nrow(X.0))
      tmp[, coords.0.indx] <- t(out$N.0.samples)
      tmp[, coords.place.indx] <- object$N.samples[, coords.indx]
      out$N.0.samples <- mcmc(tmp)
      tmp <- matrix(NA, n.post, nrow(X.0))
      tmp[, coords.0.indx] <- t(out$mu.0.samples)
      tmp[, coords.place.indx] <- object$mu.samples[, coords.indx]
      out$mu.0.samples <- mcmc(tmp)
      tmp <- matrix(NA, n.post, nrow(X.0))
      tmp[, coords.0.indx] <- t(out$w.0.samples)
      tmp[, coords.place.indx] <- object$w.samples[, coords.indx]
      out$w.0.samples <- mcmc(tmp)
    }
  }
  # Detection predictions -------------------------------------------------
  if (tolower(type) == 'detection') {
    p.det <- ncol(object$X.p)
    p.design <- p.det
    if (object$pRE & !ignore.RE) {
      p.design <- p.det + ncol(object$sigma.sq.p.samples)
    }
    if (ncol(X.0) != p.design) {
      stop(paste("error: X.0 must have ", p.design, " columns\n", sep = ''))
    }

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
      indx <- which(colnames(X.0) %in% x.p.re.names)
      if (length(indx) == 0) {
        stop("error: column names in X.0 must match variable names in data$det.covs")
      }
      X.re <- as.matrix(X.0[, indx, drop = FALSE])
      X.fix <- as.matrix(X.0[, -indx, drop = FALSE])
      n.det.re <- length(unlist(p.re.level.names))
      X.re.ind <- matrix(NA, nrow(X.re), p.det.re)
      for (i in 1:p.det.re) {
        for (j in 1:nrow(X.re)) {
          tmp <- which(p.re.level.names[[i]] == X.re[j, i])
          if (length(tmp) > 0) {
            if (i > 1) {
              X.re.ind[j, i] <- tmp + length(p.re.level.names[[i - 1]]) 
            } else {
              X.re.ind[j, i] <- tmp 
            }
          }
        }
      }
      # Create the random effects corresponding to each 
      # new location
      alpha.star.sites.0.samples <- matrix(0, n.post,  nrow(X.re))
      for (t in 1:p.det.re) {
        for (j in 1:nrow(X.re)) {
          if (!is.na(X.re.ind[j, t])) {
            alpha.star.sites.0.samples[, j] <- 
              alpha.star.samples[, X.re.ind[j, t]] + 
              alpha.star.sites.0.samples[, j]
          } else {
            alpha.star.sites.0.samples[, j] <- 
              rnorm(n.post, 0, sqrt(object$sigma.sq.p.samples[, t])) + 
              alpha.star.sites.0.samples[, j]
          }
        } # j
      } # t
    } else {
      X.fix <- X.0
      alpha.star.sites.0.samples <- matrix(0, n.post, nrow(X.0))
      p.det.re <- 0
    }
    J.str <- nrow(X.0)
    out$p.0.samples <- mcmc(logit.inv(t(X.fix %*% t(alpha.samples) + 
          				t(alpha.star.sites.0.samples))))
  }
  out$call <- cl

  class(out) <- "predict.spNMix"
  out

}
