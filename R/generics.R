# NMix --------------------------------------------------------------------
print.NMix <- function(x, ...) {
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.75)), 
      "", sep = "\n")
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
    kappa.samples <- as.matrix(object$kappa.samples)
    n.post <- object$n.post * object$n.chains
    out <- list()
    if (object$muRE) {
      p.abund.re <- length(object$re.level.names)
    } else {
      p.abund.re <- 0
    }

    if (object$muRE & !ignore.RE) {
      beta.star.samples <- object$beta.star.samples
      re.level.names <- object$re.level.names
      # Get columns in design matrix with random effects
      x.re.names <- colnames(object$X.re)
      indx <- which(colnames(X.0) %in% x.re.names)
      if (length(indx) == 0) {
        stop("error: column names in X.0 must match variable names in data$abund.covs")
      }
      X.re <- as.matrix(X.0[, indx, drop = FALSE])
      X.fix <- as.matrix(X.0[, -indx, drop = FALSE])
      n.abund.re <- length(unlist(re.level.names))
      X.re.ind <- matrix(NA, nrow(X.re), p.abund.re)
      for (i in 1:p.abund.re) {
        for (j in 1:nrow(X.re)) {
          tmp <- which(re.level.names[[i]] == X.re[j, i])
          if (length(tmp) > 0) {
            if (i > 1) {
              X.re.ind[j, i] <- tmp + length(re.level.names[[i - 1]]) 
            } else {
              X.re.ind[j, i] <- tmp 
            }
          }
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
              beta.star.samples[, X.re.ind[j, t]] + 
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
    out$mu.0.samples <- exp(t(X.fix %*% t(beta.samples) + 
          				t(beta.star.sites.0.samples)))
    out$mu.0.samples <- mcmc(apply(out$mu.0.samples, 2, 
				       function (a) a * kappa.samples))
    out$N.0.samples <- mcmc(matrix(rnbinom(length(out$mu.0.samples), kappa.samples, 
					   mu = c(out$mu.0.samples)), nrow = n.post, 
				   ncol = nrow(X.0)))
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

  class(out) <- "predict.NMix"
  out
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
  if (nrow(X.p) == nrow(y)) {
    X.p <- do.call(rbind, replicate(ncol(y), X.p, simplify = FALSE))
    X.p <- X.p[!is.na(c(y)), , drop = FALSE]
    if (object$pRE) {
      lambda.p <- do.call(rbind, replicate(ncol(y), object$lambda.p, simplify = FALSE))
      lambda.p <- lambda.p[!is.na(c(y)), , drop = FALSE]
    }
  } else {
    if (object$pRE) {
      lambda.p <- object$lambda.p
    }
  }
  y <- c(y)
  y <- y[!is.na(y)]
  N.samples <- object$N.samples
  alpha.samples <- object$alpha.samples
  if (object$pRE) {
    det.prob.samples <- t(logit.inv(X.p %*% t(alpha.samples) +
      			      lambda.p %*% t(object$alpha.star.samples)))
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
  p.abund <- ncol(object$X)
  p.design <- p.abund
  if (object$muRE & !ignore.RE) {
    p.design <- p.abund + ncol(object$sigma.sq.mu.samples)
  }
  if (dim(X.0)[3] != p.design) {
    stop(paste("error: X.0 must have ", p.design, " covariates\n", sep = ''))
  }

  # Composition sampling --------------------------------------------------
  beta.samples <- as.matrix(object$beta.samples)
  kappa.samples <- as.matrix(object$kappa.samples)
  n.post <- object$n.post * object$n.chains
  out <- list()
  if (object$muRE) {
    p.abund.re <- length(object$re.level.names)
  } else {
    p.abund.re <- 0
  }

  if (object$muRE & !ignore.RE) {
    beta.star.samples <- object$beta.star.samples
    re.level.names <- object$re.level.names
    # Get columns in design matrix with random effects
    x.re.names <- colnames(object$X.re)
    indx <- which(dimnames(X.0)[[3]] %in% x.re.names)
    if (length(indx) == 0) {
      stop("error: column names in X.0 must match variable names in data$abund.covs")
    }
    X.re <- X.0[, , indx, drop = FALSE]
    X.re <- matrix(X.re, nrow = nrow(X.re) * ncol(X.re),
		   ncol = dim(X.re)[3])
    X.fix <- X.0[, , -indx, drop = FALSE]
    X.fix <- matrix(X.fix, nrow = nrow(X.fix) * ncol(X.fix),
		    ncol = dim(X.fix)[3])
    n.abund.re <- length(unlist(re.level.names))
    X.re.ind <- matrix(NA, nrow(X.re), p.abund.re)
    for (i in 1:p.abund.re) {
      for (j in 1:nrow(X.re)) {
        if (!is.na(X.re[j, i])) {
          tmp <- which(re.level.names[[i]] == X.re[j, i])
          if (length(tmp) > 0) {
            if (i > 1) {
              X.re.ind[j, i] <- tmp + length(re.level.names[[i - 1]]) 
            } else {
              X.re.ind[j, i] <- tmp 
            }
          }
	}
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
            beta.star.samples[, X.re.ind[j, t]] + 
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
    X.fix <- matrix(X.fix, nrow = nrow(X.fix) * ncol(X.fix),
		    ncol = dim(X.fix)[3])
    beta.star.sites.0.samples <- matrix(0, n.post, nrow(X.fix))
    p.abund.re <- 0
  }
  J.str <- nrow(X.0)
  K.str.max <- ncol(X.0)
  out$mu.0.samples <- array(exp(t(X.fix %*% t(beta.samples) + 
        				t(beta.star.sites.0.samples))), 
			    dim = c(n.post, J.str, K.str.max))
  out$y.0.samples <- array(NA, dim(out$mu.0.samples))
  for (j in 1:J.str) {
    for (k in 1:K.str.max) {
      if (sum(is.na(X.0[j, k, ])) == 0) {
        out$y.0.samples[, j, k] <- rnbinom(n.post, kappa.samples, 
					   mu = out$mu.0.samples[, j, k])
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

fitted.sfMsAbund <- function(object, ...) {
  return(object$y.rep.samples)
}
