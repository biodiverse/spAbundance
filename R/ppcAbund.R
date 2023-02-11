ppcAbund <- function(object, fit.stat, group, ...) {

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
  # if (class(object) %in% c('lfCount', 'sfCount')) {
  #   stop("error: ppcOcc is not implemented for lfCount and sfCount")
  # }
  if (!(class(object) %in% c('NMix', 'spNMix', 'abund', 'spAbund', 
			     'msAbund', 'sfMsAbund', 'msNMix', 'spNMix', 'lfMsNMix', 'sfMsNMix'))) {
    stop("error: object must be one of the following classes: NMix, spNMix, abund, spAbund, msAbund, sfMsAbund, msNMix")
  }
  # Fit statistic ---------------------
  if (missing(fit.stat)) {
    stop("error: fit.stat must be specified")
  }
  if (!tolower(fit.stat) %in% c('chi-squared', 'freeman-tukey', 'chi-square')) {
    stop("error: fit.stat must be either 'chi-squared' or 'freeman-tukey'")
  }
  fit.stat <- tolower(fit.stat)
  # Group -----------------------------
  if (missing(group)) {
    stop("error: group must be specified")
  }
  if (!(group %in% c(0, 1, 2))) {
    stop("error: group must be 0 (raw data), 1 (sites), or 2 (replicates)")
  }
  if (group != 0 & class(object) %in% c('abund', 'spAbund', 'msAbund')) {
    stop("error: group must be 0 (raw data) for abundance GLM models")
  }

  out <- list()
  # Single-species N-mixture models ---------------------------------------
  if (class(object) %in% c('NMix', 'spNMix')) {
    y <- object$y
    J <- nrow(y)
    if (is(object, 'NMix')) {
      fitted.out <- fitted.NMix(object)
    } else {
      fitted.out <- fitted.spNMix(object)
    }
    N.samples <- object$N.samples
    y.rep.samples <- fitted.out$y.rep.samples
    det.prob <- fitted.out$p.samples
    n.samples <- object$n.post * object$n.chains
    fit.y <- rep(NA, n.samples)
    fit.y.rep <- rep(NA, n.samples)
    K <- apply(y, 1, function(a) sum(!is.na(a)))
    e <- 0.0001
    if (group == 0) {
      if (fit.stat %in% c('chi-squared', 'chi-square')) {
        fit.big.y.rep <- array(NA, dim = c(J, max(K), n.samples))
        fit.big.y <- array(NA, dim = c(J, max(K), n.samples))
        for (i in 1:n.samples) {
	  for (j in 1:J) {
            E <- det.prob[i, j, 1:K[j]] * N.samples[i, j]
	    fit.big.y.rep[j, 1:K[j], i] <- (y.rep.samples[i, j, 1:K[j]] - E)^2 / (E + e)
	    fit.big.y[j, 1:K[j], i] <- (y[j, 1:K[j]] - E)^2 / (E + e)
	  } # j
	  fit.y[i] <- sum(fit.big.y[, , i], na.rm = TRUE)
	  fit.y.rep[i] <- sum(fit.big.y.rep[, , i], na.rm = TRUE)
	} # i
      } else if (fit.stat == 'freeman-tukey') {
        fit.big.y.rep <- array(NA, dim = c(J, max(K), n.samples))
        fit.big.y <- array(NA, dim = c(J, max(K), n.samples))
        for (i in 1:n.samples) {
	  for (j in 1:J) {
            E <- det.prob[i, j, 1:K[j]] * N.samples[i, j]
	    fit.big.y.rep[j, 1:K[j], i] <- (y.rep.samples[i, j, 1:K[j]] - E)^2 / (E + e)
	    fit.big.y.rep[j, 1:K[j], i] <- (sqrt(y.rep.samples[i, j, 1:K[j]]) - sqrt(E))^2
	    fit.big.y[j, 1:K[j], i] <- (sqrt(y[j, 1:K[j]]) - sqrt(E))^2
	  } # j
	  fit.y[i] <- sum(fit.big.y[, , i], na.rm = TRUE)
	  fit.y.rep[i] <- sum(fit.big.y.rep[, , i], na.rm = TRUE)
	} # i
      }
    }
    # Do the stuff 
    if (group == 1) {
      y.grouped <- apply(y, 1, sum, na.rm = TRUE)
      y.rep.grouped <- apply(y.rep.samples, c(1, 2), sum, na.rm = TRUE)
      fit.big.y.rep <- matrix(NA, length(y.grouped), n.samples)
      fit.big.y <- matrix(NA, length(y.grouped), n.samples)
      if (fit.stat %in% c('chi-squared', 'chi-square')) {
        for (i in 1:n.samples) {
          E.grouped <- apply(det.prob[i, , , drop = FALSE] * N.samples[i, ], 2, sum, na.rm = TRUE)
          fit.big.y[, i] <- (y.grouped - E.grouped)^2 / (E.grouped + e)
          fit.y[i] <- sum(fit.big.y[, i])
	  fit.big.y.rep[, i] <- (y.rep.grouped[i,] - E.grouped)^2 / (E.grouped + e)
          fit.y.rep[i] <- sum(fit.big.y.rep[, i])
        }
      } else if (fit.stat == 'freeman-tukey') {
        for (i in 1:n.samples) {
          E.grouped <- apply(det.prob[i, , , drop = FALSE] * N.samples[i, ], 2, sum, na.rm = TRUE)
          fit.big.y[, i] <- (sqrt(y.grouped) - sqrt(E.grouped))^2 
          fit.y[i] <- sum(fit.big.y[, i])
	  fit.big.y.rep[, i] <- (sqrt(y.rep.grouped[i,]) - sqrt(E.grouped))^2 
          fit.y.rep[i] <- sum(fit.big.y.rep[, i])
        }
      }
    } else if (group == 2) {
      y.grouped <- apply(y, 2, sum, na.rm = TRUE)
      y.rep.grouped <- apply(y.rep.samples, c(1, 3), sum, na.rm = TRUE)
      fit.big.y <- matrix(NA, length(y.grouped), n.samples)
      fit.big.y.rep <- matrix(NA, length(y.grouped), n.samples)
      if (fit.stat %in% c('chi-squared', 'chi-square')) {
        for (i in 1:n.samples) {
          E.grouped <- apply(det.prob[i, , , drop = FALSE] * N.samples[i, ], 3, sum, na.rm = TRUE)
          fit.big.y[, i] <- (y.grouped - E.grouped)^2 / (E.grouped + e)
          fit.y[i] <- sum(fit.big.y[, i])
	  fit.big.y.rep[, i] <- (y.rep.grouped[i,] - E.grouped)^2 / (E.grouped + e)
          fit.y.rep[i] <- sum(fit.big.y.rep[, i])
        }
      } else if (fit.stat == 'freeman-tukey') {
        for (i in 1:n.samples) {
          E.grouped <- apply(det.prob[i, , , drop = FALSE] * N.samples[i, ], 3, sum, na.rm = TRUE)
          fit.big.y[, i] <- (sqrt(y.grouped) - sqrt(E.grouped))^2 
          fit.y[i] <- sum(fit.big.y[, i])
	  fit.big.y.rep[, i] <- (sqrt(y.rep.grouped[i,]) - sqrt(E.grouped))^2 
          fit.y.rep[i] <- sum(fit.big.y.rep[, i])
        }
      }
    }
    out$fit.y <- fit.y
    out$fit.y.rep <- fit.y.rep
    if (group != 0) {
      out$fit.y.group.quants <- apply(fit.big.y, 1, quantile, c(0.025, 0.25, 0.5, 0.75, 0.975))
      out$fit.y.rep.group.quants <- apply(fit.big.y.rep, 1, quantile, c(0.025, 0.25, 0.5, 0.75, 0.975))
    } else {
      out$fit.y.group.quants <- apply(fit.big.y, c(1, 2), quantile, c(0.025, 0.25, 0.5, 0.75, 0.975), na.rm = TRUE)
      out$fit.y.rep.group.quants <- apply(fit.big.y.rep, c(1, 2), quantile, c(0.025, 0.25, 0.5, 0.75, 0.975), na.rm = TRUE)
    }
    # For summaries
    out$group <- group
    out$fit.stat <- fit.stat
    out$class <- class(object)
    out$call <- cl
    out$n.samples <- object$n.samples
    out$n.burn <- object$n.burn
    out$n.thin <- object$n.thin
    out$n.post <- object$n.post
    out$n.chains <- object$n.chains
  } 
 
  # Single-species abundance models --------------------------------------- 
  if (class(object) %in% c('abund', 'spAbund')) {
    y <- object$y
    J <- nrow(y)
    if (is(object, 'abund')) {
      y.rep.samples <- fitted.abund(object)
    } else {
      y.rep.samples <- fitted.spAbund(object)
    }
    mu.samples <- object$mu.samples
    n.samples <- object$n.post * object$n.chains
    fit.y <- rep(NA, n.samples)
    fit.y.rep <- rep(NA, n.samples)
    K <- apply(y, 1, function(a) sum(!is.na(a)))
    e <- 0.0001
    if (group == 0) {
      if (fit.stat %in% c('chi-squared', 'chi-square')) {
        fit.big.y.rep <- array(NA, dim = c(J, max(K), n.samples))
        fit.big.y <- array(NA, dim = c(J, max(K), n.samples))
        for (i in 1:n.samples) {
	  for (j in 1:J) {
            E <- mu.samples[i, j, 1:K[j]]
	    fit.big.y.rep[j, 1:K[j], i] <- (y.rep.samples[i, j, 1:K[j]] - E)^2 / (E + e)
	    fit.big.y[j, 1:K[j], i] <- (y[j, 1:K[j]] - E)^2 / (E + e)
	  } # j
	  fit.y[i] <- sum(fit.big.y[, , i], na.rm = TRUE)
	  fit.y.rep[i] <- sum(fit.big.y.rep[, , i], na.rm = TRUE)
	} # i
      } else if (fit.stat == 'freeman-tukey') {
        fit.big.y.rep <- array(NA, dim = c(J, max(K), n.samples))
        fit.big.y <- array(NA, dim = c(J, max(K), n.samples))
        for (i in 1:n.samples) {
	  for (j in 1:J) {
            E <- mu.samples[i, j, 1:K[j]]
	    fit.big.y.rep[j, 1:K[j], i] <- (y.rep.samples[i, j, 1:K[j]] - E)^2 / (E + e)
	    fit.big.y.rep[j, 1:K[j], i] <- (sqrt(y.rep.samples[i, j, 1:K[j]]) - sqrt(E))^2
	    fit.big.y[j, 1:K[j], i] <- (sqrt(y[j, 1:K[j]]) - sqrt(E))^2
	  } # j
	  fit.y[i] <- sum(fit.big.y[, , i], na.rm = TRUE)
	  fit.y.rep[i] <- sum(fit.big.y.rep[, , i], na.rm = TRUE)
	} # i
      }
    }
    out$fit.y <- fit.y
    out$fit.y.rep <- fit.y.rep
    out$fit.y.group.quants <- apply(fit.big.y, c(1, 2), quantile, c(0.025, 0.25, 0.5, 0.75, 0.975), na.rm = TRUE)
    out$fit.y.rep.group.quants <- apply(fit.big.y.rep, c(1, 2), quantile, c(0.025, 0.25, 0.5, 0.75, 0.975), na.rm = TRUE)
    # For summaries
    out$group <- group
    out$fit.stat <- fit.stat
    out$class <- class(object)
    out$call <- cl
    out$n.samples <- object$n.samples
    out$n.burn <- object$n.burn
    out$n.thin <- object$n.thin
    out$n.post <- object$n.post
    out$n.chains <- object$n.chains
  } 

  # Multi-species abundance models ----------------------------------------
  if (class(object) %in% c('msAbund', 'spMsAbund', 'lfMsAbund', 'sfMsAbund')) {
    y <- object$y
    J <- dim(y)[2] 
    n.sp <- dim(y)[1]
    if (is(object, 'msAbund')) {
      y.rep.samples <- fitted.msAbund(object)
    }
    if (is(object, 'spMsAbund')) {
      # y.rep.samples <- fitted.spMsAbund(object)
    }
    if (is(object, 'lfMsAbund')) {
      # y.rep.samples <- fitted.lfMsAbund(object)
    }
    if (is(object, 'sfMsAbund')) {
      y.rep.samples <- fitted.sfMsAbund(object)
    }
    mu.samples <- object$mu.samples
    n.samples <- object$n.post * object$n.chains
    fit.y <- matrix(NA, n.samples, n.sp)
    fit.y.rep <- matrix(NA, n.samples, n.sp)
    K <- apply(y[1, , , drop = FALSE], 2, function(a) sum(!is.na(a)))
    e <- 0.0001
    if (group == 0) {
      for (i in 1:n.sp) {
        message(noquote(paste("Currently on species ", i, " out of ", n.sp, sep = '')))
        if (fit.stat %in% c('chi-squared', 'chi-square')) {
          fit.big.y.rep <- array(NA, dim = c(n.sp, J, max(K), n.samples))
          fit.big.y <- array(NA, dim = c(n.sp, J, max(K), n.samples))
          for (t in 1:n.samples) {
            for (j in 1:J) {
              E <- mu.samples[t, i, j, 1:K[j]]
              fit.big.y.rep[i, j, 1:K[j], t] <- (y.rep.samples[t, i, j, 1:K[j]] - E)^2 / (E + e)
              fit.big.y[i, j, 1:K[j], t] <- (y[i, j, 1:K[j]] - E)^2 / (E + e)
            } # j
            fit.y[t, i] <- sum(fit.big.y[i, , , t], na.rm = TRUE)
            fit.y.rep[t, i] <- sum(fit.big.y.rep[i, , , t], na.rm = TRUE)
          } # t
        } else if (fit.stat == 'freeman-tukey') {
          fit.big.y.rep <- array(NA, dim = c(n.sp, J, max(K), n.samples))
          fit.big.y <- array(NA, dim = c(n.sp, J, max(K), n.samples))
          for (t in 1:n.samples) {
            for (j in 1:J) {
              E <- mu.samples[t, i, j, 1:K[j]]
              fit.big.y.rep[i, j, 1:K[j], t] <- (y.rep.samples[t, i, j, 1:K[j]] - E)^2 / (E + e)
              fit.big.y.rep[i, j, 1:K[j], t] <- (sqrt(y.rep.samples[t, i, j, 1:K[j]]) - sqrt(E))^2
              fit.big.y[i, j, 1:K[j], t] <- (sqrt(y[i, j, 1:K[j]]) - sqrt(E))^2
            } # j
            fit.y[t, i] <- sum(fit.big.y[i, , , t], na.rm = TRUE)
            fit.y.rep[t, i] <- sum(fit.big.y.rep[i, , , t], na.rm = TRUE)
          } # t
        }
      }
    }
    out$fit.y <- fit.y
    out$fit.y.rep <- fit.y.rep
    out$fit.y.group.quants <- apply(fit.big.y, c(1, 2, 3), quantile, c(0.025, 0.25, 0.5, 0.75, 0.975), na.rm = TRUE)
    out$fit.y.rep.group.quants <- apply(fit.big.y.rep, c(1, 2, 3), quantile, c(0.025, 0.25, 0.5, 0.75, 0.975), na.rm = TRUE)
    # For summaries
    out$group <- group
    out$fit.stat <- fit.stat
    out$class <- class(object)
    out$call <- cl
    out$n.samples <- object$n.samples
    out$n.burn <- object$n.burn
    out$n.thin <- object$n.thin
    out$n.post <- object$n.post
    out$n.chains <- object$n.chains
  }

  # Multi-species N-mixture models ----------------------------------------
  if (class(object) %in% c('msNMix', 'spMsNMix', 'lfMsNMix', 'sfMsNMix')) {
    y <- object$y
    J <- dim(y)[2]
    n.sp <- dim(y)[1]
    # Fitted function is the same for all multispecies N-mixtures. 
    fitted.out <- fitted.msNMix(object)
    y.rep.samples <- fitted.out$y.rep.samples
    det.prob <- fitted.out$p.samples
    N.samples <- object$N.samples
    n.samples <- object$n.post * object$n.chains
    fit.y <- matrix(NA, n.samples, n.sp)
    fit.y.rep <- matrix(NA, n.samples, n.sp)
    K <- apply(y[1, , ], 1, function(a) sum(!is.na(a)))
    e <- 0.0001
    # Do the stuff 
    if (group == 0) {
      fit.big.y.rep <- array(NA, dim = dim(y.rep.samples))
      fit.big.y <- array(NA, dim = dim(y.rep.samples))
      for (i in 1:n.sp) {
        message(noquote(paste("Currently on species ", i, " out of ", n.sp, sep = '')))
        if (fit.stat %in% c('chi-squared', 'chi-square')) {
            for (j in 1:n.samples) {
              for (k in 1:J) {
                E.grouped <- det.prob[j, i, k, 1:K[k]] * N.samples[j, i, k]
                fit.big.y[j, i, k, 1:K[k]] <- (y[i, k, 1:K[k]] - E.grouped)^2 / (E.grouped + e)
                fit.y[j, i] <- sum(fit.big.y[j, i, , ], na.rm = TRUE)
                fit.big.y.rep[j, i, k, 1:K[k]] <- (y.rep.samples[j, i, k, 1:K[k]] - E.grouped)^2 / (E.grouped + e)
                fit.y.rep[j, i] <- sum(fit.big.y.rep[j, i, , ], na.rm = TRUE)
	      }
            }
        } else if (fit.stat == 'freeman-tukey') {
          for (j in 1:n.samples) {
            for (k in 1:J) {
              E.grouped <- det.prob[j, i, k, 1:K[k]] * N.samples[j, i, k]
              fit.big.y[j, i, k, 1:K[k]] <- (sqrt(y[i, k, 1:K[k]]) - sqrt(E.grouped))^2 
              fit.y[j, i] <- sum(fit.big.y[j, i, , ], na.rm = TRUE)
              fit.big.y.rep[j, i, k, 1:K[k]] <- (sqrt(y.rep.samples[j, i, k, 1:K[k]]) - sqrt(E.grouped))^2 
              fit.y.rep[j, i] <- sum(fit.big.y.rep[j, i, , ], na.rm = TRUE)
	    }
          }
        }
      }
    }
    if (group == 1) {
      y.grouped <- apply(y, c(1, 2), sum, na.rm = TRUE)
      y.rep.grouped <- apply(y.rep.samples, c(1, 2, 3), sum, na.rm = TRUE)
      fit.big.y.rep <- array(NA, dim = c(n.samples, n.sp, J))
      fit.big.y <- array(NA, dim = c(n.samples, n.sp, J))
      for (i in 1:n.sp) {
        message(noquote(paste("Currently on species ", i, " out of ", n.sp, sep = '')))
        if (fit.stat %in% c('chi-squared', 'chi-square')) {
            for (j in 1:n.samples) {
              E.grouped <- apply(det.prob[j, i, , , drop = FALSE] * N.samples[j, i, ], 3, sum, na.rm = TRUE)
              fit.big.y[j, i, ] <- (y.grouped[i, ] - E.grouped)^2 / (E.grouped + e)
              fit.y[j, i] <- sum(fit.big.y[j, i, ])
              fit.big.y.rep[j, i, ] <- (y.rep.grouped[j, i, ] - E.grouped)^2 / (E.grouped + e)
              fit.y.rep[j, i] <- sum(fit.big.y.rep[j, i, ])
            }
        } else if (fit.stat == 'freeman-tukey') {
          for (j in 1:n.samples) {
            E.grouped <- apply(det.prob[j, i, , , drop = FALSE] * N.samples[j, i, ], 3, sum, na.rm = TRUE)
            fit.big.y[j, i, ] <- (sqrt(y.grouped[i, ]) - sqrt(E.grouped))^2 
            fit.y[j, i] <- sum(fit.big.y[j, i, ])
            fit.big.y.rep[j, i, ] <- (sqrt(y.rep.grouped[j, i, ]) - sqrt(E.grouped))^2 
            fit.y.rep[j, i] <- sum(fit.big.y.rep[j, i, ])
          }
        }
      }
    } else if (group == 2) {
      y.grouped <- apply(y, c(1, 3), sum, na.rm = TRUE)
      y.rep.grouped <- apply(y.rep.samples, c(1, 2, 4), sum, na.rm = TRUE)
      fit.big.y <- array(NA, dim = c(n.samples, n.sp, dim(y)[3]))
      fit.big.y.rep <- array(NA, dim = c(n.samples, n.sp, dim(y)[3]))
      for (i in 1:n.sp) {
        message(noquote(paste("Currently on species ", i, " out of ", n.sp, sep = '')))
        if (fit.stat %in% c('chi-squared', 'chi-square')) {
          for (j in 1:n.samples) {
            E.grouped <- apply(det.prob[j, i, , , drop = FALSE] * N.samples[j, i, ], 4, sum, na.rm = TRUE)
            fit.big.y[j, i, ] <- (y.grouped[i, ] - E.grouped)^2 / (E.grouped + e)
            fit.y[j, i] <- sum(fit.big.y[j, i, ])
            fit.big.y.rep[j, i, ] <- (y.rep.grouped[j, i, ] - E.grouped)^2 / (E.grouped + e)
            fit.y.rep[j, i] <- sum(fit.big.y.rep[j, i, ])
          }
        } else if (fit.stat == 'freeman-tukey') {
          for (j in 1:n.samples) {
            E.grouped <- apply(det.prob[j, i, , , drop = FALSE] * N.samples[j, i, ], 4, sum, na.rm = TRUE)
            fit.big.y[j, i, ] <- (sqrt(y.grouped[i, ]) - sqrt(E.grouped))^2 
            fit.y[j, i] <- sum(fit.big.y[j, i, ])
            fit.big.y.rep[j, i, ] <- (sqrt(y.rep.grouped[j, i, ]) - sqrt(E.grouped))^2 
            fit.y.rep[j, i] <- sum(fit.big.y.rep[j, i, ])
          }
        }
      }
    }
    out$fit.y <- fit.y
    out$fit.y.rep <- fit.y.rep
    if (group != 0) {
      out$fit.y.group.quants <- apply(fit.big.y, c(2, 3), quantile, c(0.025, 0.25, 0.5, 0.75, 0.975))
      out$fit.y.rep.group.quants <- apply(fit.big.y.rep, c(2, 3), quantile, c(0.025, 0.25, 0.5, 0.75, 0.975))
    } else {
      out$fit.y.group.quants <- apply(fit.big.y, c(2, 3, 4), quantile, c(0.025, 0.25, 0.5, 0.75, 0.975), na.rm = TRUE)
      out$fit.y.rep.group.quants <- apply(fit.big.y.rep, c(2, 3, 4), quantile, c(0.025, 0.25, 0.5, 0.75, 0.975), na.rm = TRUE)
    }
    # For summaries
    out$group <- group
    out$fit.stat <- fit.stat
    out$class <- class(object)
    out$call <- cl
    out$n.samples <- object$n.samples
    out$n.burn <- object$n.burn
    out$n.thin <- object$n.thin
    out$n.post <- object$n.post
    out$n.chains <- object$n.chains
    out$sp.names <- object$sp.names
  }

  class(out) <- 'ppcAbund'

  return(out)

}
