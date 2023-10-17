waicAbund <- function(object, N.max, by.species = FALSE, ...) {

  # Check for unused arguments ------------------------------------------
  formal.args <- names(formals(sys.function(sys.parent())))
  elip.args <- names(list(...))
  for(i in elip.args){
      if(! i %in% formal.args)
          warning("'",i, "' is not an argument")
  }
  # Call ----------------------------------------------------------------
  cl <- match.call()

  # Some initial checks -------------------------------------------------
  # Object ----------------------------
  if (missing(object)) {
    stop("error: object must be specified")
  }
  if (!(class(object) %in% c('NMix', 'spNMix', 'abund', 'spAbund', 
			     'msAbund', 'lfMsAbund', 
			     'sfMsAbund', 'msNMix', 
			     'lfMsNMix', 'sfMsNMix', 'DS', 'spDS', 
			     'msDS', 'lfMsDS', 'sfMsDS', 'svcAbund'))) {
    stop("error: object must be one of the following classes: abund, spAbund, NMix, spNMix, msAbund, lfMsAbund, sfMsAbund, msNMix, lfMsNMix, sfMsNMix, DS, spDS, msDS, lfMsDS, sfMsDS, svcAbund\n")
  }

  if (!(class(object) %in% c('abund', 'spAbund', 'msAbund', 'lfMsAbund', 'sfMsAbund', 'svcAbund'))) {
    if (missing(N.max)) {
      message("N.max not specified. Setting upper index of integration of N to 10 plus\nthe largest estimated abundance value at each site in object$N.samples")
      if (class(object) %in% c('msNMix', 'lfMsNMix', 'sfMsNMix', 
			       'msDS', 'lfMsDS', 'sfMsDS')) {
        N.max <- apply(object$N.samples, c(2, 3), max) + 10
      } else {
        N.max <- apply(object$N.samples, 2, max) + 10 
      }
    }
  }

  if (class(object) %in% c('NMix', 'spNMix', 'DS', 'spDS')) {
    if (!(length(N.max) %in% c(1, nrow(object$y)))) {
      stop(paste("N.max must be of length 1 or ", nrow(object$y), sep = ''))
    }
    if (length(N.max) == 1) {
      N.max <- rep(N.max, nrow(object$y))
    }
  }
  
  if (class(object) %in% c('msNMix', 'lfMsNMix', 'sfMsNMix', 
			   'msDS', 'lfMsDS', 'sfMsDS')) {
    if (is.null(dim(N.max))) {
      if (length(N.max) == 1) {
        N.max <- matrix(N.max, nrow(object$y), ncol(object$y))
      } else if (length(N.max) == nrow(object$y)) {
        N.max <- matrix(N.max, nrow(object$y), ncol(object$y))
      } else {
        stop("N.max is not correctly specified.")
      }
    } else if (length(dim(N.max)) != 2) {
      stop("N.max is not correctly specified.")
    }
  }

  if (!(by.species %in% c(TRUE, FALSE))) {
    stop("by.species must be set to equal TRUE or FALSE")
  }

  n.post <- object$n.post * object$n.chains
  # Functions -------------------------------------------------------------
  logit <- function(theta, a = 0, b = 1) {log((theta-a)/(b-theta))}
  logit.inv <- function(z, a = 0, b = 1) {b-(b-a)/(1+exp(z))}

  # Single-species abundance GLMs -----------------------------------------
  if (class(object) %in% c('abund', 'spAbund', 'svcAbund')) {
    if (object$dist %in% c('Gaussian', 'zi-Gaussian')) {
      if (object$dist == 'zi-Gaussian') {
        message("Calculated WAIC is only for stage 2 of the hurdle model\n")
      }
      elpd <- sum(apply(object$like.samples, c(2), function(a) log(mean(a))), na.rm = TRUE)
      pD <- sum(apply(object$like.samples, c(2), function(a) var(log(a))), na.rm = TRUE)
    }
    if (object$dist %in% c('Poisson', 'NB')) {
      elpd <- sum(apply(object$like.samples, c(2, 3), function(a) log(mean(a))), na.rm = TRUE)
      pD <- sum(apply(object$like.samples, c(2, 3), function(a) var(log(a))), na.rm = TRUE)
    }
    out <- c(elpd, pD, -2 * (elpd - pD))
    names(out) <- c("elpd", "pD", "WAIC")
  }

  # Multi-species abundance GLMs ------------------------------------------
  if (class(object) %in% c('msAbund', 'lfMsAbund', 'sfMsAbund')) {
    if (object$dist %in% c('Gaussian', 'zi-Gaussian')) {
      if (object$dist == 'zi-Gaussian') {
        message("Calculated WAIC is only for stage 2 of the hurdle model\n")
      }
      margin <- c(2, 3)
    }  
    if (object$dist %in% c('Poisson', 'NB')) {
      margin <- c(2, 3, 4)
    }
    if (by.species) {
      elpd <- apply(apply(object$like.samples, margin, function(a) log(mean(a))), 
                    1, sum, na.rm = TRUE) 
      pD <- apply(apply(object$like.samples, margin, function(a) var(log(a))), 
                  1, sum, na.rm = TRUE)
      out <- data.frame(elpd = elpd, 
			pD = pD, 
			WAIC = -2 * (elpd - pD))
    } else {
      elpd <- sum(apply(object$like.samples, margin, function(a) log(mean(a))), na.rm = TRUE)
      pD <- sum(apply(object$like.samples, margin, function(a) var(log(a))), na.rm = TRUE)
      out <- c(elpd, pD, -2 * (elpd - pD))
      names(out) <- c("elpd", "pD", "WAIC")
    }
  }

  # Single-species N-mixture models ---------------------------------------
  if (class(object) %in% c('NMix', 'spNMix')) {
    N.samples <- object$N.samples
    n.samples <- object$n.post * object$n.chains
    kappa.samples <- object$kappa.samples
    y <- object$y
    y.max <- apply(y, 1, max, na.rm = TRUE)
    y.na <- ifelse(is.na(y), 1, 0)
    K.max <- dim(y)[2]
    mu.samples <- t(apply(object$mu.samples, 1, function(a) a * object$offset))
    if (is(object, 'NMix')) {
      p.samples <- fitted.NMix(object)$p.samples
    } else {
      p.samples <- fitted.spNMix(object)$p.samples
    }
    # dist == 1 for NB, 0 for Poisson
    dist <- object$dist
    dist <- ifelse(dist == 'NB', 1, 0)
    J <- nrow(y)
    # Model Type: 0 = single-species N-mixture
    model.type <- 0

    storage.mode(J) <- "integer"
    storage.mode(N.samples) <- "double"
    storage.mode(n.samples) <- "integer"
    storage.mode(kappa.samples) <- "double"
    storage.mode(y) <- "double"
    storage.mode(y.na) <- 'integer'
    storage.mode(mu.samples) <- "double"
    storage.mode(p.samples) <- "double"
    storage.mode(dist) <- "integer"
    storage.mode(model.type) <- "integer"
    storage.mode(N.max) <- "integer"
    storage.mode(K.max) <- "integer"
    storage.mode(y.max) <- "integer"

    tmp <- .Call("waicAbund", J, y.na, dist, model.type, 
		 y, n.samples, N.samples, 
		 kappa.samples, mu.samples, p.samples, 
		 N.max, K.max, y.max)
    elpd <- sum(apply(tmp$like.samples, 1, function(a) log(mean(a))), na.rm = TRUE)
    pD <- sum(apply(tmp$like.samples, 1, function(a) var(log(a))), na.rm = TRUE)
    out <- c(elpd, pD, -2 * (elpd - pD))
    names(out) <- c("elpd", "pD", "WAIC")
  }

  # Multi-species N-mixture models ----------------------------------------
  if (class(object) %in% c('msNMix', 'spMsNMix', 'lfMsNMix', 'sfMsNMix')) {
    n.sp <- nrow(object$y)
    elpd <- rep(NA, n.sp)
    pD <- rep(NA, n.sp)
    p.samples <- fitted.msNMix(object)$p.samples
    for (i in 1:n.sp) {
      message(noquote(paste("Currently on species ", i, " out of ", n.sp, sep = '')))
      N.samples <- object$N.samples[, i, ]
      n.samples <- object$n.post * object$n.chains
      kappa.samples <- object$kappa.samples[, i]
      y <- object$y[i, , ]
      if (length(dim(y)) == 1) {
        y <- as.matrix(y)
      }
      y.max <- apply(y, 1, max, na.rm = TRUE)
      # y.max <- ifelse(y.max == 0, 1, y.max)
      y.na <- ifelse(is.na(y), 1, 0)
      K.max <- dim(y)[2]
      mu.samples <- object$mu.samples[, i, ]
      mu.samples <- t(apply(mu.samples, 1, function(a) a * object$offset))
      p.samples.curr <- p.samples[, i, , ]
      # dist == 1 for NB, 0 for Poisson
      dist <- object$dist
      dist <- ifelse(dist == 'NB', 1, 0)
      J <- nrow(y)
      N.max.curr <- N.max[i, ]
      # Model Type: 0 = single-species N-mixture
      model.type <- 0

      storage.mode(J) <- "integer"
      storage.mode(N.samples) <- "double"
      storage.mode(n.samples) <- "integer"
      storage.mode(kappa.samples) <- "double"
      storage.mode(y) <- "double"
      storage.mode(y.na) <- "integer"
      storage.mode(mu.samples) <- "double"
      storage.mode(p.samples.curr) <- "double"
      storage.mode(dist) <- "integer"
      storage.mode(model.type) <- "integer"
      storage.mode(N.max.curr) <- "integer"
      storage.mode(K.max) <- "integer"
      storage.mode(y.max) <- "integer"

      tmp <- .Call("waicAbund", J, y.na, dist, model.type, 
          	 y, n.samples, N.samples, 
          	 kappa.samples, mu.samples, p.samples.curr, 
          	 N.max.curr, K.max, y.max)
      elpd[i] <- sum(apply(tmp$like.samples, 1, function(a) log(mean(a))), na.rm = TRUE)
      pD[i] <- sum(apply(tmp$like.samples, 1, function(a) var(log(a))), na.rm = TRUE)
    }
    if (by.species == FALSE) {
      out <- c(sum(elpd), sum(pD), -2 * (sum(elpd) - sum(pD)))
      names(out) <- c("elpd", "pD", "WAIC")
    } else {
      out <- data.frame(elpd = elpd, pD = pD, WAIC = -2 * (elpd - pD))
    }
  }

  # Single-species distance sampling --------------------------------------
  if (class(object) %in% c('DS', 'spDS')) {
    N.samples <- object$N.samples
    n.samples <- object$n.post * object$n.chains
    kappa.samples <- object$kappa.samples
    y <- object$y
    J <- nrow(y)
    y.sum <- apply(y, 1, sum, na.rm = TRUE)
    y.na <- ifelse(is.na(y), 1, 0)
    # The number of distance bins
    K.max <- dim(y)[2]
    mu.samples <- t(apply(object$mu.samples, 1, function(a) a * object$offset))
    pi.samples <- object$pi.samples
    pi.samples <- array(NA, dim = c(n.samples, J, K.max + 1))
    pi.samples[, , 1:K.max] <- object$pi.samples
    pi.samples[, , K.max + 1] <- apply(object$pi.samples, c(1, 2),
				       function(a) 1 - sum(a))
    # dist == 1 for NB, 0 for Poisson
    dist <- object$dist
    dist <- ifelse(dist == 'NB', 1, 0)
    # Model Type: 0 = single-species N-mixture, 
    #             1 = single-species distance sampling
    model.type <- 1

    storage.mode(J) <- "integer"
    storage.mode(N.samples) <- "double"
    storage.mode(n.samples) <- "integer"
    storage.mode(kappa.samples) <- "double"
    storage.mode(y) <- "double"
    storage.mode(y.na) <- "integer"
    storage.mode(mu.samples) <- "double"
    storage.mode(pi.samples) <- "double"
    storage.mode(dist) <- "integer"
    storage.mode(model.type) <- "integer"
    storage.mode(N.max) <- "integer"
    storage.mode(K.max) <- "integer"
    storage.mode(y.sum) <- "integer"

    tmp <- .Call("waicAbund", J, y.na, dist, model.type, 
		 y, n.samples, N.samples, 
		 kappa.samples, mu.samples, pi.samples, 
		 N.max, K.max, y.sum)
    elpd <- sum(apply(tmp$like.samples, 1, function(a) log(mean(a))), na.rm = TRUE)
    pD <- sum(apply(tmp$like.samples, 1, function(a) var(log(a))), na.rm = TRUE)
    out <- c(elpd, pD, -2 * (elpd - pD))
    names(out) <- c("elpd", "pD", "WAIC")
  }
 
  # Multi-species distance sampling --------------------------------------- 
  if (class(object) %in% c('msDS', 'sfMsDS', 'lfMsDS')) {
    n.sp <- nrow(object$y)
    elpd <- rep(NA, n.sp)
    pD <- rep(NA, n.sp)
    for (i in 1:n.sp) {
      message(noquote(paste("Currently on species ", i, " out of ", n.sp, sep = '')))
      N.samples <- object$N.samples[, i, ]
      n.samples <- object$n.post * object$n.chains
      kappa.samples <- object$kappa.samples[, i]
      y <- object$y[i, , ]
      if (length(dim(y)) == 1) {
        y <- as.matrix(y)
      }
      y.sum <- apply(y, 1, sum, na.rm = TRUE)
      J <- nrow(y)
      y.na <- ifelse(is.na(y), 1, 0)
      K.max <- dim(y)[2]
      pi.samples <- array(NA, dim = c(n.samples, J, K.max + 1))
      pi.samples[, , 1:K.max] <- object$pi.samples[, i, , ]
      pi.samples[, , K.max + 1] <- apply(pi.samples[, , 1:K.max], c(1, 2),
				       function(a) 1 - sum(a))
      mu.samples <- object$mu.samples[, i, ]
      mu.samples <- t(apply(mu.samples, 1, function(a) a * object$offset))
      # dist == 1 for NB, 0 for Poisson
      dist <- object$dist
      dist <- ifelse(dist == 'NB', 1, 0)
      N.max.curr <- N.max[i, ]
      # Model Type: 0 = single-species N-mixture, 
      #             1 = single-species distance sampling
      model.type <- 1

      storage.mode(J) <- "integer"
      storage.mode(N.samples) <- "double"
      storage.mode(n.samples) <- "integer"
      storage.mode(kappa.samples) <- "double"
      storage.mode(y) <- "double"
      storage.mode(y.na) <- "integer"
      storage.mode(mu.samples) <- "double"
      storage.mode(pi.samples) <- "double"
      storage.mode(dist) <- "integer"
      storage.mode(model.type) <- "integer"
      storage.mode(N.max.curr) <- "integer"
      storage.mode(K.max) <- "integer"
      storage.mode(y.sum) <- "integer"

      tmp <- .Call("waicAbund", J, y.na, dist, model.type, 
          	 y, n.samples, N.samples, 
          	 kappa.samples, mu.samples, pi.samples, 
          	 N.max.curr, K.max, y.sum)
      elpd[i] <- sum(apply(tmp$like.samples, 1, function(a) log(mean(a))), na.rm = TRUE)
      pD[i] <- sum(apply(tmp$like.samples, 1, function(a) var(log(a))), na.rm = TRUE)
    }
    if (by.species == FALSE) {
      out <- c(sum(elpd), sum(pD), -2 * (sum(elpd) - sum(pD)))
      names(out) <- c("elpd", "pD", "WAIC")
    } else {
      out <- data.frame(elpd = elpd, pD = pD, WAIC = -2 * (elpd - pD))
    }
  }
  return(out)

}
