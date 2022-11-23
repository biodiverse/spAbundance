waicAbund <- function(object, N.max, ...) {

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
			     'msAbund', 'sfMsAbund'))) {
    stop("error: object must be one of the following classes: abund, spAbund, NMix, spNMix, msAbund, sfMsAbund\n")
  }

  if (!(class(object) %in% c('abund', 'spAbund', 'msAbund', 'sfMsAbund'))) {
    if (missing(N.max)) {
      message("N.max not specified. Setting upper index of integration of\nN to two times the largest estimated abundance value in object$N.samples")
      N.max <- max(object$N.samples) * 2
    }
  }

  n.post <- object$n.post * object$n.chains
  # Functions -------------------------------------------------------------
  logit <- function(theta, a = 0, b = 1) {log((theta-a)/(b-theta))}
  logit.inv <- function(z, a = 0, b = 1) {b-(b-a)/(1+exp(z))}

  # Single-species abundance GLMs -----------------------------------------
  if (class(object) %in% c('abund', 'spAbund')) {
    elpd <- sum(apply(object$like.samples, c(2, 3), function(a) log(mean(a))), na.rm = TRUE)
    pD <- sum(apply(object$like.samples, c(2, 3), function(a) var(log(a))), na.rm = TRUE)
    out <- c(elpd, pD, -2 * (elpd - pD))
    names(out) <- c("elpd", "pD", "WAIC")
  }

  # Multi-species abundance GLMs ------------------------------------------
  if (class(object) %in% c('msAbund', 'sfMsAbund')) {
    elpd <- sum(apply(object$like.samples, c(2, 3, 4), function(a) log(mean(a))), na.rm = TRUE)
    pD <- sum(apply(object$like.samples, c(2, 3, 4), function(a) var(log(a))), na.rm = TRUE)
    out <- c(elpd, pD, -2 * (elpd - pD))
    names(out) <- c("elpd", "pD", "WAIC")
  }

  # Single-species N-mixture models ---------------------------------------
  if (class(object) %in% c('NMix', 'spNMix')) {
    N.samples <- object$N.samples
    n.samples <- object$n.post * object$n.chains
    kappa.samples <- object$kappa.samples
    y <- object$y
    y.max <- apply(y, 1, max, na.rm = TRUE)
    # TODO: 
    # y.max <- ifelse(y.max == 0, 1, y.max)
    K <- apply(y, 1, function(a) sum(!is.na(a)))
    K.max <- max(K)
    mu.samples <- object$mu.samples
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
    storage.mode(K) <- "integer"
    storage.mode(mu.samples) <- "double"
    storage.mode(p.samples) <- "double"
    storage.mode(dist) <- "integer"
    storage.mode(model.type) <- "integer"
    storage.mode(N.max) <- "integer"
    storage.mode(K.max) <- "integer"
    storage.mode(y.max) <- "integer"

    tmp <- .Call("waicAbund", J, K, dist, model.type, 
		 y, n.samples, N.samples, 
		 kappa.samples, mu.samples, p.samples, 
		 N.max, K.max, y.max)
    elpd <- sum(apply(tmp$like.samples, 1, function(a) log(mean(a))), na.rm = TRUE)
    pD <- sum(apply(tmp$like.samples, 1, function(a) var(log(a))), na.rm = TRUE)
    out <- c(elpd, pD, -2 * (elpd - pD))
    names(out) <- c("elpd", "pD", "WAIC")
  }


  return(out)

}
