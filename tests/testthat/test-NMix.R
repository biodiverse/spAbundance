# Test NMix.R -------------------------------------------------------------

skip_on_cran()

# Intercept only ----------------------------------------------------------
set.seed(1000)
J.x <- 15
J.y <- 15
J <- J.x * J.y
n.rep <- sample(2:4, J, replace = TRUE)
beta <- c(0)
p.abund <- length(beta)
alpha <- c(0.5)
p.det <- length(alpha)
mu.RE <- list()
p.RE <- list()
family <- 'NB'
kappa <- 0.5
dat <- simNMix(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, alpha = alpha,
	       mu.RE = mu.RE, p.RE = p.RE, sp = FALSE, family = family, kappa = kappa)

y <- dat$y
X <- dat$X
X.re <- dat$X.re
X.p <- dat$X.p
X.p.re <- dat$X.p.re

abund.covs <- X
colnames(abund.covs) <- c('int')
det.covs <- list(int = X.p[, , 1])
data.list <- list(y = y,
		  abund.covs = abund.covs,
 		  det.covs = det.covs)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 10),
		   alpha.normal = list(mean = 0,
			               var = 2.72))
# Starting values
inits.list <- list(alpha = 0,
		   beta = 0,
		   kappa = 1.3,
		   sigma.sq.mu = 0.5,
		   N = apply(y, 1, max, na.rm = TRUE))
tuning <- list(beta = 0.1, alpha = 0.1, beta.star = 0.1, alpha.star = 0.1,
               kappa = 0.2)
n.batch <- 50
batch.length <- 25
n.burn <- 750
n.thin <- 2
n.chains <- 2
abund.formula <- ~ 1
det.formula <- ~ 1

data = data.list
n.batch = n.batch
batch.length = batch.length
inits = inits.list
priors = prior.list
accept.rate = 0.43
n.omp.threads = 1
verbose = TRUE
family <- 'NB'
n.report = 50
k.fold = 3
k.fold.threads = 3
k.fold.seed <- 100
k.fold.only <- FALSE


out <- NMix(abund.formula = abund.formula,
	    det.formula = det.formula,
	    data = data.list,
	    n.batch = n.batch,
	    batch.length = batch.length,
	    inits = inits.list,
	    family = 'NB',
	    tuning = tuning,
	    priors = prior.list,
	    accept.rate = 0.43,
	    n.omp.threads = 1,
	    verbose = FALSE,
	    n.report = 50,
	    n.burn = n.burn,
	    n.thin = n.thin,
	    n.chains = n.chains,
            k.fold = 2,
            k.fold.threads = 1)

# Test to make sure it worked ---------
test_that("out is of class NMix", {
  expect_s3_class(out, "NMix")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$rmspe), 1)
  expect_type(out$rmspe, "double")
  expect_gt(out$rmspe, 0)
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$pRE, FALSE)
  expect_equal(out$muRE, FALSE)
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, dat$y)
})

# Check Poisson -----------------------
test_that("Poisson works", {
  out <- NMix(abund.formula = abund.formula,
  	    det.formula = det.formula,
  	    data = data.list,
  	    n.batch = n.batch,
  	    batch.length = batch.length,
  	    inits = inits.list,
  	    family = 'Poisson',
  	    tuning = tuning,
  	    priors = prior.list,
  	    accept.rate = 0.43,
  	    n.omp.threads = 1,
  	    verbose = FALSE,
  	    n.report = 50,
  	    n.burn = n.burn,
  	    n.thin = n.thin,
  	    n.chains = 1)
  expect_s3_class(out, "NMix")
})

# Check default priors ----------------
test_that("default priors, multiple chains, inits, burn, thin work", {
  out <- NMix(abund.formula = abund.formula,
	      det.formula = det.formula,
	      data = data.list,
	      n.batch = n.batch,
	      batch.length = batch.length,
	      family = 'NB',
	      accept.rate = 0.43,
	      n.omp.threads = 1,
	      verbose = FALSE,
	      n.report = 50,
	      n.burn = n.burn,
	      n.thin = n.thin,
	      n.chains = n.chains)
  expect_s3_class(out, 'NMix')
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {
  expect_output(NMix(abund.formula = abund.formula,
  	    det.formula = det.formula,
  	    data = data.list,
  	    n.batch = n.batch,
  	    batch.length = batch.length,
  	    inits = inits.list,
  	    family = 'Poisson',
  	    tuning = tuning,
  	    priors = prior.list,
  	    accept.rate = 0.43,
  	    n.omp.threads = 1,
  	    verbose = TRUE,
  	    n.report = 50,
  	    n.burn = n.burn,
  	    n.thin = n.thin,
  	    n.chains = 1))
})

# Check waicAbund -----------------------
test_that("waicAbund works for NMix", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicAbund(out, max(out$N.samples) + 10))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for NMix", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for NMix", {
  X.0 <- dat$X
  pred.out <- predict(out, X.0)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$mu.0.samples), c(out$n.post * out$n.chains, J))
  expect_equal(dim(pred.out$N.0.samples), c(out$n.post * out$n.chains, J))
})
test_that("detection prediction works", {
  J.str <- 100
  X.p.0 <- matrix(1, nrow = J.str, ncol = p.det)
  pred.out <- predict(out, X.p.0, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, J.str))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for NMix", {
  n.post.samples <- out$n.post * out$n.chains
  ppc.out <- ppcAbund(out, 'chi-square', 2)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, max(n.rep)))
  
  ppc.out <- ppcAbund(out, 'chi-square', 1)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J))
  
  ppc.out <- ppcAbund(out, 'freeman-tukey', 1)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J))
  
  ppc.out <- ppcAbund(out, 'freeman-tukey', 2)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, max(n.rep)))

  ppc.out <- ppcAbund(out, 'chi-square', 0)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J, max(n.rep)))
  
  ppc.out <- ppcAbund(out, 'freeman-tukey', 0)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J, max(n.rep)))
})

# Abundance covariate only ------------------------------------------------
set.seed(1000)
J.x <- 15
J.y <- 15
J <- J.x * J.y
n.rep <- sample(2:4, J, replace = TRUE)
beta <- c(0, -0.3)
p.abund <- length(beta)
alpha <- c(0.5)
p.det <- length(alpha)
mu.RE <- list()
p.RE <- list()
family <- 'NB'
kappa <- 0.5
dat <- simNMix(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, alpha = alpha,
	       mu.RE = mu.RE, p.RE = p.RE, sp = FALSE, family = family, kappa = kappa)

y <- dat$y
X <- dat$X
X.re <- dat$X.re
X.p <- dat$X.p
X.p.re <- dat$X.p.re

abund.covs <- X
colnames(abund.covs) <- c('int', 'abund.cov.1')
det.covs <- list(int = X.p[, , 1])
data.list <- list(y = y,
		  abund.covs = abund.covs,
 		  det.covs = det.covs)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 10),
		   alpha.normal = list(mean = 0,
			               var = 2.72))
# Starting values
inits.list <- list(alpha = 0,
		   beta = 0,
		   kappa = 1.3,
		   sigma.sq.mu = 0.5,
		   N = apply(y, 1, max, na.rm = TRUE))
tuning <- list(beta = 0.1, alpha = 0.1, beta.star = 0.1, alpha.star = 0.1,
               kappa = 0.2)
n.batch <- 50
batch.length <- 25
n.burn <- 750
n.thin <- 2
n.chains <- 2
abund.formula <- ~ abund.cov.1
det.formula <- ~ 1

data = data.list
n.batch = n.batch
batch.length = batch.length
inits = inits.list
priors = prior.list
accept.rate = 0.43
n.omp.threads = 1
verbose = TRUE
family <- 'NB'
n.report = 50
k.fold = 3
k.fold.threads = 3
k.fold.seed <- 100
k.fold.only <- FALSE


out <- NMix(abund.formula = abund.formula,
	    det.formula = det.formula,
	    data = data.list,
	    n.batch = n.batch,
	    batch.length = batch.length,
	    inits = inits.list,
	    family = 'NB',
	    tuning = tuning,
	    priors = prior.list,
	    accept.rate = 0.43,
	    n.omp.threads = 1,
	    verbose = FALSE,
	    n.report = 50,
	    n.burn = n.burn,
	    n.thin = n.thin,
	    n.chains = n.chains,
            k.fold = 2,
            k.fold.threads = 1)

# Test to make sure it worked ---------
test_that("out is of class NMix", {
  expect_s3_class(out, "NMix")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$rmspe), 1)
  expect_type(out$rmspe, "double")
  expect_gt(out$rmspe, 0)
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$pRE, FALSE)
  expect_equal(out$muRE, FALSE)
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, dat$y)
})

# Check Poisson -----------------------
test_that("Poisson works", {
  out <- NMix(abund.formula = abund.formula,
  	    det.formula = det.formula,
  	    data = data.list,
  	    n.batch = n.batch,
  	    batch.length = batch.length,
  	    inits = inits.list,
  	    family = 'Poisson',
  	    tuning = tuning,
  	    priors = prior.list,
  	    accept.rate = 0.43,
  	    n.omp.threads = 1,
  	    verbose = FALSE,
  	    n.report = 50,
  	    n.burn = n.burn,
  	    n.thin = n.thin,
  	    n.chains = 1)
  expect_s3_class(out, "NMix")
})

# Check default priors ----------------
test_that("default priors, multiple chains, inits, burn, thin work", {
  out <- NMix(abund.formula = abund.formula,
	      det.formula = det.formula,
	      data = data.list,
	      n.batch = n.batch,
	      batch.length = batch.length,
	      family = 'NB',
	      accept.rate = 0.43,
	      n.omp.threads = 1,
	      verbose = FALSE,
	      n.report = 50,
	      n.burn = n.burn,
	      n.thin = n.thin,
	      n.chains = n.chains)
  expect_s3_class(out, 'NMix')
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {
  expect_output(NMix(abund.formula = abund.formula,
  	    det.formula = det.formula,
  	    data = data.list,
  	    n.batch = n.batch,
  	    batch.length = batch.length,
  	    inits = inits.list,
  	    family = 'Poisson',
  	    tuning = tuning,
  	    priors = prior.list,
  	    accept.rate = 0.43,
  	    n.omp.threads = 1,
  	    verbose = TRUE,
  	    n.report = 50,
  	    n.burn = n.burn,
  	    n.thin = n.thin,
  	    n.chains = 1))
})

# Check waicAbund -----------------------
test_that("waicAbund works for NMix", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicAbund(out, max(out$N.samples) + 10))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for NMix", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for NMix", {
  X.0 <- dat$X
  pred.out <- predict(out, X.0)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$mu.0.samples), c(out$n.post * out$n.chains, J))
  expect_equal(dim(pred.out$N.0.samples), c(out$n.post * out$n.chains, J))
})
test_that("detection prediction works", {
  J.str <- 100
  X.p.0 <- matrix(1, nrow = J.str, ncol = p.det)
  pred.out <- predict(out, X.p.0, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, J.str))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for NMix", {
  n.post.samples <- out$n.post * out$n.chains
  ppc.out <- ppcAbund(out, 'chi-square', 2)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, max(n.rep)))
  
  ppc.out <- ppcAbund(out, 'chi-square', 1)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J))
  
  ppc.out <- ppcAbund(out, 'freeman-tukey', 1)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J))
  
  ppc.out <- ppcAbund(out, 'freeman-tukey', 2)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, max(n.rep)))

  ppc.out <- ppcAbund(out, 'chi-square', 0)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J, max(n.rep)))
  
  ppc.out <- ppcAbund(out, 'freeman-tukey', 0)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J, max(n.rep)))
})

# Detection covariate only ----------------------------------------------------------
set.seed(1000)
J.x <- 15
J.y <- 15
J <- J.x * J.y
n.rep <- sample(2:4, J, replace = TRUE)
beta <- c(0)
p.abund <- length(beta)
alpha <- c(0.5, -0.3)
p.det <- length(alpha)
mu.RE <- list()
p.RE <- list()
family <- 'NB'
kappa <- 0.5
dat <- simNMix(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, alpha = alpha,
	       mu.RE = mu.RE, p.RE = p.RE, sp = FALSE, family = family, kappa = kappa)

y <- dat$y
X <- dat$X
X.re <- dat$X.re
X.p <- dat$X.p
X.p.re <- dat$X.p.re

abund.covs <- X
colnames(abund.covs) <- c('int')
det.covs <- list(int = X.p[, , 1], 
                 det.cov.1 = X.p[, , 2])
data.list <- list(y = y,
		  abund.covs = abund.covs,
 		  det.covs = det.covs)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 10),
		   alpha.normal = list(mean = 0,
			               var = 2.72))
# Starting values
inits.list <- list(alpha = 0,
		   beta = 0,
		   kappa = 1.3,
		   sigma.sq.mu = 0.5,
		   N = apply(y, 1, max, na.rm = TRUE))
tuning <- list(beta = 0.1, alpha = 0.1, beta.star = 0.1, alpha.star = 0.1,
               kappa = 0.2)
n.batch <- 50
batch.length <- 25
n.burn <- 750
n.thin <- 2
n.chains <- 2
abund.formula <- ~ 1
det.formula <- ~ det.cov.1

data = data.list
n.batch = n.batch
batch.length = batch.length
inits = inits.list
priors = prior.list
accept.rate = 0.43
n.omp.threads = 1
verbose = TRUE
family <- 'NB'
n.report = 50
k.fold = 3
k.fold.threads = 3
k.fold.seed <- 100
k.fold.only <- FALSE


out <- NMix(abund.formula = abund.formula,
	    det.formula = det.formula,
	    data = data.list,
	    n.batch = n.batch,
	    batch.length = batch.length,
	    inits = inits.list,
	    family = 'NB',
	    tuning = tuning,
	    priors = prior.list,
	    accept.rate = 0.43,
	    n.omp.threads = 1,
	    verbose = FALSE,
	    n.report = 50,
	    n.burn = n.burn,
	    n.thin = n.thin,
	    n.chains = n.chains,
            k.fold = 2,
            k.fold.threads = 1)

# Test to make sure it worked ---------
test_that("out is of class NMix", {
  expect_s3_class(out, "NMix")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$rmspe), 1)
  expect_type(out$rmspe, "double")
  expect_gt(out$rmspe, 0)
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$pRE, FALSE)
  expect_equal(out$muRE, FALSE)
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, dat$y)
})

# Check Poisson -----------------------
test_that("Poisson works", {
  out <- NMix(abund.formula = abund.formula,
  	    det.formula = det.formula,
  	    data = data.list,
  	    n.batch = n.batch,
  	    batch.length = batch.length,
  	    inits = inits.list,
  	    family = 'Poisson',
  	    tuning = tuning,
  	    priors = prior.list,
  	    accept.rate = 0.43,
  	    n.omp.threads = 1,
  	    verbose = FALSE,
  	    n.report = 50,
  	    n.burn = n.burn,
  	    n.thin = n.thin,
  	    n.chains = 1)
  expect_s3_class(out, "NMix")
})

# Check default priors ----------------
test_that("default priors, multiple chains, inits, burn, thin work", {
  out <- NMix(abund.formula = abund.formula,
	      det.formula = det.formula,
	      data = data.list,
	      n.batch = n.batch,
	      batch.length = batch.length,
	      family = 'NB',
	      accept.rate = 0.43,
	      n.omp.threads = 1,
	      verbose = FALSE,
	      n.report = 50,
	      n.burn = n.burn,
	      n.thin = n.thin,
	      n.chains = n.chains)
  expect_s3_class(out, 'NMix')
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {
  expect_output(NMix(abund.formula = abund.formula,
  	    det.formula = det.formula,
  	    data = data.list,
  	    n.batch = n.batch,
  	    batch.length = batch.length,
  	    inits = inits.list,
  	    family = 'Poisson',
  	    tuning = tuning,
  	    priors = prior.list,
  	    accept.rate = 0.43,
  	    n.omp.threads = 1,
  	    verbose = TRUE,
  	    n.report = 50,
  	    n.burn = n.burn,
  	    n.thin = n.thin,
  	    n.chains = 1))
})

# Check waicAbund -----------------------
test_that("waicAbund works for NMix", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicAbund(out, max(out$N.samples) + 10))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for NMix", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for NMix", {
  X.0 <- dat$X
  pred.out <- predict(out, X.0)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$mu.0.samples), c(out$n.post * out$n.chains, J))
  expect_equal(dim(pred.out$N.0.samples), c(out$n.post * out$n.chains, J))
})
test_that("detection prediction works", {
  J.str <- 100
  X.p.0 <- matrix(1, nrow = J.str, ncol = p.det)
  pred.out <- predict(out, X.p.0, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, J.str))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for NMix", {
  n.post.samples <- out$n.post * out$n.chains
  ppc.out <- ppcAbund(out, 'chi-square', 2)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, max(n.rep)))
  
  ppc.out <- ppcAbund(out, 'chi-square', 1)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J))
  
  ppc.out <- ppcAbund(out, 'freeman-tukey', 1)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J))
  
  ppc.out <- ppcAbund(out, 'freeman-tukey', 2)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, max(n.rep)))

  ppc.out <- ppcAbund(out, 'chi-square', 0)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J, max(n.rep)))
  
  ppc.out <- ppcAbund(out, 'freeman-tukey', 0)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J, max(n.rep)))
})

# Covariates on both ------------------------------------------------------ 
set.seed(1000)
J.x <- 15
J.y <- 15
J <- J.x * J.y
n.rep <- sample(2:4, J, replace = TRUE)
beta <- c(0, 0.3, -0.2)
p.abund <- length(beta)
alpha <- c(0.5, -0.3)
p.det <- length(alpha)
mu.RE <- list()
p.RE <- list()
family <- 'NB'
kappa <- 0.5
dat <- simNMix(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, alpha = alpha,
	       mu.RE = mu.RE, p.RE = p.RE, sp = FALSE, family = family, kappa = kappa)

y <- dat$y
X <- dat$X
X.re <- dat$X.re
X.p <- dat$X.p
X.p.re <- dat$X.p.re

abund.covs <- X
colnames(abund.covs) <- c('int', 'abund.cov.1', 'abund.cov.2')
det.covs <- list(int = X.p[, , 1], 
                 det.cov.1 = X.p[, , 2])
data.list <- list(y = y,
		  abund.covs = abund.covs,
 		  det.covs = det.covs)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 10),
		   alpha.normal = list(mean = 0,
			               var = 2.72))
# Starting values
inits.list <- list(alpha = 0,
		   beta = 0,
		   kappa = 1.3,
		   sigma.sq.mu = 0.5,
		   N = apply(y, 1, max, na.rm = TRUE))
tuning <- list(beta = 0.1, alpha = 0.1, beta.star = 0.1, alpha.star = 0.1,
               kappa = 0.2)
n.batch <- 50
batch.length <- 25
n.burn <- 750
n.thin <- 2
n.chains <- 2
abund.formula <- ~ abund.cov.1 + abund.cov.2
det.formula <- ~ det.cov.1

data = data.list
n.batch = n.batch
batch.length = batch.length
inits = inits.list
priors = prior.list
accept.rate = 0.43
n.omp.threads = 1
verbose = TRUE
family <- 'NB'
n.report = 50
k.fold = 3
k.fold.threads = 3
k.fold.seed <- 100
k.fold.only <- FALSE


out <- NMix(abund.formula = abund.formula,
	    det.formula = det.formula,
	    data = data.list,
	    n.batch = n.batch,
	    batch.length = batch.length,
	    inits = inits.list,
	    family = 'NB',
	    tuning = tuning,
	    priors = prior.list,
	    accept.rate = 0.43,
	    n.omp.threads = 1,
	    verbose = FALSE,
	    n.report = 50,
	    n.burn = n.burn,
	    n.thin = n.thin,
	    n.chains = n.chains,
            k.fold = 2,
            k.fold.threads = 1)

# Test to make sure it worked ---------
test_that("out is of class NMix", {
  expect_s3_class(out, "NMix")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$rmspe), 1)
  expect_type(out$rmspe, "double")
  expect_gt(out$rmspe, 0)
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$pRE, FALSE)
  expect_equal(out$muRE, FALSE)
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, dat$y)
})

# Check Poisson -----------------------
test_that("Poisson works", {
  out <- NMix(abund.formula = abund.formula,
  	    det.formula = det.formula,
  	    data = data.list,
  	    n.batch = n.batch,
  	    batch.length = batch.length,
  	    inits = inits.list,
  	    family = 'Poisson',
  	    tuning = tuning,
  	    priors = prior.list,
  	    accept.rate = 0.43,
  	    n.omp.threads = 1,
  	    verbose = FALSE,
  	    n.report = 50,
  	    n.burn = n.burn,
  	    n.thin = n.thin,
  	    n.chains = 1)
  expect_s3_class(out, "NMix")
})

# Check default priors ----------------
test_that("default priors, multiple chains, inits, burn, thin work", {
  out <- NMix(abund.formula = abund.formula,
	      det.formula = det.formula,
	      data = data.list,
	      n.batch = n.batch,
	      batch.length = batch.length,
	      family = 'NB',
	      accept.rate = 0.43,
	      n.omp.threads = 1,
	      verbose = FALSE,
	      n.report = 50,
	      n.burn = n.burn,
	      n.thin = n.thin,
	      n.chains = n.chains)
  expect_s3_class(out, 'NMix')
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {
  expect_output(NMix(abund.formula = abund.formula,
  	    det.formula = det.formula,
  	    data = data.list,
  	    n.batch = n.batch,
  	    batch.length = batch.length,
  	    inits = inits.list,
  	    family = 'Poisson',
  	    tuning = tuning,
  	    priors = prior.list,
  	    accept.rate = 0.43,
  	    n.omp.threads = 1,
  	    verbose = TRUE,
  	    n.report = 50,
  	    n.burn = n.burn,
  	    n.thin = n.thin,
  	    n.chains = 1))
})

# Check waicAbund -----------------------
test_that("waicAbund works for NMix", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicAbund(out, max(out$N.samples) + 10))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for NMix", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for NMix", {
  X.0 <- dat$X
  pred.out <- predict(out, X.0)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$mu.0.samples), c(out$n.post * out$n.chains, J))
  expect_equal(dim(pred.out$N.0.samples), c(out$n.post * out$n.chains, J))
})
test_that("detection prediction works", {
  X.p.0 <- dat$X.p[, 1, ]
  pred.out <- predict(out, X.p.0, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, J))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for NMix", {
  n.post.samples <- out$n.post * out$n.chains
  ppc.out <- ppcAbund(out, 'chi-square', 2)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, max(n.rep)))
  
  ppc.out <- ppcAbund(out, 'chi-square', 1)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J))
  
  ppc.out <- ppcAbund(out, 'freeman-tukey', 1)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J))
  
  ppc.out <- ppcAbund(out, 'freeman-tukey', 2)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, max(n.rep)))

  ppc.out <- ppcAbund(out, 'chi-square', 0)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J, max(n.rep)))
  
  ppc.out <- ppcAbund(out, 'freeman-tukey', 0)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J, max(n.rep)))
})

# Interactions on both ----------------------------------------------------
set.seed(1000)
J.x <- 15
J.y <- 15
J <- J.x * J.y
n.rep <- sample(2:4, J, replace = TRUE)
beta <- c(0, 0.3, -0.2)
p.abund <- length(beta)
alpha <- c(0.5, -0.3, -0.3)
p.det <- length(alpha)
mu.RE <- list()
p.RE <- list()
family <- 'NB'
kappa <- 0.5
dat <- simNMix(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, alpha = alpha,
	       mu.RE = mu.RE, p.RE = p.RE, sp = FALSE, family = family, kappa = kappa)

y <- dat$y
X <- dat$X
X.re <- dat$X.re
X.p <- dat$X.p
X.p.re <- dat$X.p.re

abund.covs <- X
colnames(abund.covs) <- c('int', 'abund.cov.1', 'abund.cov.2')
det.covs <- list(int = X.p[, , 1], 
                 det.cov.1 = X.p[, , 2], 
                 det.cov.2 = X.p[, , 3])
data.list <- list(y = y,
		  abund.covs = abund.covs,
 		  det.covs = det.covs)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 10),
		   alpha.normal = list(mean = 0,
			               var = 2.72))
# Starting values
inits.list <- list(alpha = 0,
		   beta = 0,
		   kappa = 1.3,
		   sigma.sq.mu = 0.5,
		   N = apply(y, 1, max, na.rm = TRUE))
tuning <- list(beta = 0.1, alpha = 0.1, beta.star = 0.1, alpha.star = 0.1,
               kappa = 0.2)
n.batch <- 50
batch.length <- 25
n.burn <- 750
n.thin <- 2
n.chains <- 2
abund.formula <- ~ abund.cov.1 * abund.cov.2
det.formula <- ~ det.cov.1 * det.cov.2

data = data.list
n.batch = n.batch
batch.length = batch.length
inits = inits.list
priors = prior.list
accept.rate = 0.43
n.omp.threads = 1
verbose = TRUE
family <- 'NB'
n.report = 50
k.fold = 3
k.fold.threads = 3
k.fold.seed <- 100
k.fold.only <- FALSE


out <- NMix(abund.formula = abund.formula,
	    det.formula = det.formula,
	    data = data.list,
	    n.batch = n.batch,
	    batch.length = batch.length,
	    inits = inits.list,
	    family = 'NB',
	    tuning = tuning,
	    priors = prior.list,
	    accept.rate = 0.43,
	    n.omp.threads = 1,
	    verbose = FALSE,
	    n.report = 50,
	    n.burn = n.burn,
	    n.thin = n.thin,
	    n.chains = n.chains,
            k.fold = 2,
            k.fold.threads = 1)

# Test to make sure it worked ---------
test_that("out is of class NMix", {
  expect_s3_class(out, "NMix")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$rmspe), 1)
  expect_type(out$rmspe, "double")
  expect_gt(out$rmspe, 0)
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$pRE, FALSE)
  expect_equal(out$muRE, FALSE)
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, dat$y)
})

# Check Poisson -----------------------
test_that("Poisson works", {
  out <- NMix(abund.formula = abund.formula,
  	    det.formula = det.formula,
  	    data = data.list,
  	    n.batch = n.batch,
  	    batch.length = batch.length,
  	    inits = inits.list,
  	    family = 'Poisson',
  	    tuning = tuning,
  	    priors = prior.list,
  	    accept.rate = 0.43,
  	    n.omp.threads = 1,
  	    verbose = FALSE,
  	    n.report = 50,
  	    n.burn = n.burn,
  	    n.thin = n.thin,
  	    n.chains = 1)
  expect_s3_class(out, "NMix")
})

# Check default priors ----------------
test_that("default priors, multiple chains, inits, burn, thin work", {
  out <- NMix(abund.formula = abund.formula,
	      det.formula = det.formula,
	      data = data.list,
	      n.batch = n.batch,
	      batch.length = batch.length,
	      family = 'NB',
	      accept.rate = 0.43,
	      n.omp.threads = 1,
	      verbose = FALSE,
	      n.report = 50,
	      n.burn = n.burn,
	      n.thin = n.thin,
	      n.chains = n.chains)
  expect_s3_class(out, 'NMix')
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {
  expect_output(NMix(abund.formula = abund.formula,
  	    det.formula = det.formula,
  	    data = data.list,
  	    n.batch = n.batch,
  	    batch.length = batch.length,
  	    inits = inits.list,
  	    family = 'Poisson',
  	    tuning = tuning,
  	    priors = prior.list,
  	    accept.rate = 0.43,
  	    n.omp.threads = 1,
  	    verbose = TRUE,
  	    n.report = 50,
  	    n.burn = n.burn,
  	    n.thin = n.thin,
  	    n.chains = 1))
})

# Check waicAbund -----------------------
test_that("waicAbund works for NMix", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicAbund(out, max(out$N.samples) + 10))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for NMix", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for NMix", {
  X.0 <- dat$X
  X.0 <- cbind(X.0, X.0[, 2] * X.0[, 3])
  pred.out <- predict(out, X.0)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$mu.0.samples), c(out$n.post * out$n.chains, J))
  expect_equal(dim(pred.out$N.0.samples), c(out$n.post * out$n.chains, J))
})
test_that("detection prediction works", {
  X.p.0 <- dat$X.p[, 1, ]
  X.p.0 <- cbind(X.p.0, X.p.0[, 2] * X.p.0[, 3])
  pred.out <- predict(out, X.p.0, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, J))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for NMix", {
  n.post.samples <- out$n.post * out$n.chains
  ppc.out <- ppcAbund(out, 'chi-square', 2)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, max(n.rep)))
  
  ppc.out <- ppcAbund(out, 'chi-square', 1)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J))
  
  ppc.out <- ppcAbund(out, 'freeman-tukey', 1)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J))
  
  ppc.out <- ppcAbund(out, 'freeman-tukey', 2)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, max(n.rep)))

  ppc.out <- ppcAbund(out, 'chi-square', 0)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J, max(n.rep)))
  
  ppc.out <- ppcAbund(out, 'freeman-tukey', 0)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J, max(n.rep)))
})

# Site covariate on detection ----------------------------------------------
set.seed(1000)
J.x <- 15
J.y <- 15
J <- J.x * J.y
n.rep <- sample(2:4, J, replace = TRUE)
beta <- c(0, 0.3, -0.2)
p.abund <- length(beta)
alpha <- c(0.5, -0.3)
p.det <- length(alpha)
mu.RE <- list()
p.RE <- list()
family <- 'NB'
kappa <- 0.5
dat <- simNMix(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, alpha = alpha,
	       mu.RE = mu.RE, p.RE = p.RE, sp = FALSE, family = family, kappa = kappa)

y <- dat$y
X <- dat$X
X.re <- dat$X.re
X.p <- dat$X.p
X.p.re <- dat$X.p.re

abund.covs <- X
colnames(abund.covs) <- c('int', 'abund.cov.1', 'abund.cov.2')
det.covs <- list(det.cov.1 = X[, 2])
data.list <- list(y = y,
		  abund.covs = abund.covs,
 		  det.covs = det.covs)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 10),
		   alpha.normal = list(mean = 0,
			               var = 2.72))
# Starting values
inits.list <- list(alpha = 0,
		   beta = 0,
		   kappa = 1.3,
		   sigma.sq.mu = 0.5,
		   N = apply(y, 1, max, na.rm = TRUE))
tuning <- list(beta = 0.1, alpha = 0.1, beta.star = 0.1, alpha.star = 0.1,
               kappa = 0.2)
n.batch <- 50
batch.length <- 25
n.burn <- 750
n.thin <- 2
n.chains <- 2
abund.formula <- ~ abund.cov.1 + abund.cov.2
det.formula <- ~ det.cov.1

data = data.list
n.batch = n.batch
batch.length = batch.length
inits = inits.list
priors = prior.list
accept.rate = 0.43
n.omp.threads = 1
verbose = TRUE
family <- 'NB'
n.report = 50
k.fold = 3
k.fold.threads = 3
k.fold.seed <- 100
k.fold.only <- FALSE


out <- NMix(abund.formula = abund.formula,
	    det.formula = det.formula,
	    data = data.list,
	    n.batch = n.batch,
	    batch.length = batch.length,
	    inits = inits.list,
	    family = 'NB',
	    tuning = tuning,
	    priors = prior.list,
	    accept.rate = 0.43,
	    n.omp.threads = 1,
	    verbose = FALSE,
	    n.report = 50,
	    n.burn = n.burn,
	    n.thin = n.thin,
	    n.chains = n.chains,
            k.fold = 2,
            k.fold.threads = 1)

# Test to make sure it worked ---------
test_that("out is of class NMix", {
  expect_s3_class(out, "NMix")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$rmspe), 1)
  expect_type(out$rmspe, "double")
  expect_gt(out$rmspe, 0)
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$pRE, FALSE)
  expect_equal(out$muRE, FALSE)
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, dat$y)
})

# Check Poisson -----------------------
test_that("Poisson works", {
  out <- NMix(abund.formula = abund.formula,
  	    det.formula = det.formula,
  	    data = data.list,
  	    n.batch = n.batch,
  	    batch.length = batch.length,
  	    inits = inits.list,
  	    family = 'Poisson',
  	    tuning = tuning,
  	    priors = prior.list,
  	    accept.rate = 0.43,
  	    n.omp.threads = 1,
  	    verbose = FALSE,
  	    n.report = 50,
  	    n.burn = n.burn,
  	    n.thin = n.thin,
  	    n.chains = 1)
  expect_s3_class(out, "NMix")
})

# Check default priors ----------------
test_that("default priors, multiple chains, inits, burn, thin work", {
  out <- NMix(abund.formula = abund.formula,
	      det.formula = det.formula,
	      data = data.list,
	      n.batch = n.batch,
	      batch.length = batch.length,
	      family = 'NB',
	      accept.rate = 0.43,
	      n.omp.threads = 1,
	      verbose = FALSE,
	      n.report = 50,
	      n.burn = n.burn,
	      n.thin = n.thin,
	      n.chains = n.chains)
  expect_s3_class(out, 'NMix')
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {
  expect_output(NMix(abund.formula = abund.formula,
  	    det.formula = det.formula,
  	    data = data.list,
  	    n.batch = n.batch,
  	    batch.length = batch.length,
  	    inits = inits.list,
  	    family = 'Poisson',
  	    tuning = tuning,
  	    priors = prior.list,
  	    accept.rate = 0.43,
  	    n.omp.threads = 1,
  	    verbose = TRUE,
  	    n.report = 50,
  	    n.burn = n.burn,
  	    n.thin = n.thin,
  	    n.chains = 1))
})

# Check waicAbund -----------------------
test_that("waicAbund works for NMix", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicAbund(out, max(out$N.samples) + 10))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for NMix", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for NMix", {
  X.0 <- dat$X
  pred.out <- predict(out, X.0)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$mu.0.samples), c(out$n.post * out$n.chains, J))
  expect_equal(dim(pred.out$N.0.samples), c(out$n.post * out$n.chains, J))
})
test_that("detection prediction works", {
  X.p.0 <- dat$X.p[, 1, ]
  pred.out <- predict(out, X.p.0, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, J))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for NMix", {
  n.post.samples <- out$n.post * out$n.chains
  ppc.out <- ppcAbund(out, 'chi-square', 2)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, max(n.rep)))
  
  ppc.out <- ppcAbund(out, 'chi-square', 1)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J))
  
  ppc.out <- ppcAbund(out, 'freeman-tukey', 1)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J))
  
  ppc.out <- ppcAbund(out, 'freeman-tukey', 2)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, max(n.rep)))

  ppc.out <- ppcAbund(out, 'chi-square', 0)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J, max(n.rep)))
  
  ppc.out <- ppcAbund(out, 'freeman-tukey', 0)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J, max(n.rep)))
})

# Random intercept on abundance -------------------------------------------
set.seed(1000)
J.x <- 15
J.y <- 15
J <- J.x * J.y
n.rep <- sample(2:4, J, replace = TRUE)
beta <- c(0)
p.abund <- length(beta)
alpha <- c(0.5)
p.det <- length(alpha)
mu.RE <- list(levels = c(15),
	       sigma.sq.mu = c(1.3))
p.RE <- list()
family <- 'NB'
kappa <- 0.5
dat <- simNMix(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, alpha = alpha,
	       mu.RE = mu.RE, p.RE = p.RE, sp = FALSE, family = family, kappa = kappa)

y <- dat$y
X <- dat$X
X.re <- dat$X.re
X.p <- dat$X.p
X.p.re <- dat$X.p.re

abund.covs <- cbind(X, X.re)
colnames(abund.covs) <- c('int', 'abund.factor.1')
det.covs <- list(int = X.p[, , 1])
data.list <- list(y = y,
		  abund.covs = abund.covs,
 		  det.covs = det.covs)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 10),
		   alpha.normal = list(mean = 0,
			               var = 2.72))
# Starting values
inits.list <- list(alpha = 0,
		   beta = 0,
		   kappa = 1.3,
		   sigma.sq.mu = 0.5,
		   N = apply(y, 1, max, na.rm = TRUE))
tuning <- list(beta = 0.1, alpha = 0.1, beta.star = 0.1, alpha.star = 0.1,
               kappa = 0.2)
n.batch <- 50
batch.length <- 25
n.burn <- 750
n.thin <- 2
n.chains <- 2
abund.formula <- ~ (1 | abund.factor.1)
det.formula <- ~ 1 

data = data.list
n.batch = n.batch
batch.length = batch.length
inits = inits.list
priors = prior.list
accept.rate = 0.43
n.omp.threads = 1
verbose = TRUE
family <- 'NB'
n.report = 50
k.fold = 3
k.fold.threads = 3
k.fold.seed <- 100
k.fold.only <- FALSE


out <- NMix(abund.formula = abund.formula,
	    det.formula = det.formula,
	    data = data.list,
	    n.batch = n.batch,
	    batch.length = batch.length,
	    inits = inits.list,
	    family = 'NB',
	    tuning = tuning,
	    priors = prior.list,
	    accept.rate = 0.43,
	    n.omp.threads = 1,
	    verbose = FALSE,
	    n.report = 50,
	    n.burn = n.burn,
	    n.thin = n.thin,
	    n.chains = n.chains,
            k.fold = 2,
            k.fold.threads = 1)

# Test to make sure it worked ---------
test_that("out is of class NMix", {
  expect_s3_class(out, "NMix")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$rmspe), 1)
  expect_type(out$rmspe, "double")
  expect_gt(out$rmspe, 0)
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$pRE, FALSE)
  expect_equal(out$muRE, TRUE)
})

test_that("random effect gives error when non-numeric", {
  data.list$abund.covs <- as.data.frame(data.list$abund.covs)
  data.list$abund.covs$abund.factor.1 <- factor(data.list$abund.covs$abund.factor.1)
  expect_error(out <- NMix(abund.formula = abund.formula,
	    det.formula = det.formula,
	    data = data.list,
	    n.batch = n.batch,
	    batch.length = batch.length,
	    inits = inits.list,
	    family = 'NB',
	    tuning = tuning,
	    priors = prior.list,
	    accept.rate = 0.43,
	    n.omp.threads = 1,
	    verbose = FALSE,
	    n.report = 50,
	    n.burn = n.burn,
	    n.thin = n.thin,
	    n.chains = n.chains))
  data.list$abund.covs$abund.factor.1 <- factor(data.list$abund.covs$abund.factor.1)
  expect_error(out <- NMix(abund.formula = abund.formula,
	    det.formula = det.formula,
	    data = data.list,
	    n.batch = n.batch,
	    batch.length = batch.length,
	    inits = inits.list,
	    family = 'NB',
	    tuning = tuning,
	    priors = prior.list,
	    accept.rate = 0.43,
	    n.omp.threads = 1,
	    verbose = FALSE,
	    n.report = 50,
	    n.burn = n.burn,
	    n.thin = n.thin,
	    n.chains = n.chains))
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, dat$y)
})

# Check Poisson -----------------------
test_that("Poisson works", {
  out <- NMix(abund.formula = abund.formula,
  	    det.formula = det.formula,
  	    data = data.list,
  	    n.batch = n.batch,
  	    batch.length = batch.length,
  	    inits = inits.list,
  	    family = 'Poisson',
  	    tuning = tuning,
  	    priors = prior.list,
  	    accept.rate = 0.43,
  	    n.omp.threads = 1,
  	    verbose = FALSE,
  	    n.report = 50,
  	    n.burn = n.burn,
  	    n.thin = n.thin,
  	    n.chains = 1)
  expect_s3_class(out, "NMix")
})

# Check default priors ----------------
test_that("default priors, multiple chains, inits, burn, thin work", {
  out <- NMix(abund.formula = abund.formula,
	      det.formula = det.formula,
	      data = data.list,
	      n.batch = n.batch,
	      batch.length = batch.length,
	      family = 'NB',
	      accept.rate = 0.43,
	      n.omp.threads = 1,
	      verbose = FALSE,
	      n.report = 50,
	      n.burn = n.burn,
	      n.thin = n.thin,
	      n.chains = n.chains)
  expect_s3_class(out, 'NMix')
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {
  expect_output(NMix(abund.formula = abund.formula,
  	    det.formula = det.formula,
  	    data = data.list,
  	    n.batch = n.batch,
  	    batch.length = batch.length,
  	    inits = inits.list,
  	    family = 'Poisson',
  	    tuning = tuning,
  	    priors = prior.list,
  	    accept.rate = 0.43,
  	    n.omp.threads = 1,
  	    verbose = TRUE,
  	    n.report = 50,
  	    n.burn = n.burn,
  	    n.thin = n.thin,
  	    n.chains = 1))
})

# Check waicAbund -----------------------
test_that("waicAbund works for NMix", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicAbund(out, max(out$N.samples) + 10))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for NMix", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for NMix", {
  X.0 <- cbind(X, X.re)
  colnames(X.0) <- c('(Intercept)', 'abund.factor.1')
  pred.out <- predict(out, X.0)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$mu.0.samples), c(out$n.post * out$n.chains, J))
  expect_equal(dim(pred.out$N.0.samples), c(out$n.post * out$n.chains, J))
})
test_that("detection prediction works", {
  X.p.0 <- as.matrix(dat$X.p[, 1, ])
  pred.out <- predict(out, X.p.0, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, J))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for NMix", {
  n.post.samples <- out$n.post * out$n.chains
  ppc.out <- ppcAbund(out, 'chi-square', 2)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, max(n.rep)))
  
  ppc.out <- ppcAbund(out, 'chi-square', 1)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J))
  
  ppc.out <- ppcAbund(out, 'freeman-tukey', 1)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J))
  
  ppc.out <- ppcAbund(out, 'freeman-tukey', 2)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, max(n.rep)))

  ppc.out <- ppcAbund(out, 'chi-square', 0)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J, max(n.rep)))
  
  ppc.out <- ppcAbund(out, 'freeman-tukey', 0)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J, max(n.rep)))
})

# Multiple random intercepts on abundance ---------------------------------
set.seed(1000)
J.x <- 15
J.y <- 15
J <- J.x * J.y
n.rep <- sample(2:4, J, replace = TRUE)
beta <- c(0)
p.abund <- length(beta)
alpha <- c(0.5)
p.det <- length(alpha)
mu.RE <- list(levels = c(15, 10, 12),
	       sigma.sq.mu = c(1.3, 0.2, 0.5))
p.RE <- list()
family <- 'NB'
kappa <- 0.5
dat <- simNMix(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, alpha = alpha,
	       mu.RE = mu.RE, p.RE = p.RE, sp = FALSE, family = family, kappa = kappa)

y <- dat$y
X <- dat$X
X.re <- dat$X.re
X.p <- dat$X.p
X.p.re <- dat$X.p.re

abund.covs <- cbind(X, X.re)
colnames(abund.covs) <- c('int', 'abund.factor.1', 'abund.factor.2', 'abund.factor.3')
det.covs <- list(int = X.p[, , 1])
data.list <- list(y = y,
		  abund.covs = abund.covs,
 		  det.covs = det.covs)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 10),
		   alpha.normal = list(mean = 0,
			               var = 2.72))
# Starting values
inits.list <- list(alpha = 0,
		   beta = 0,
		   kappa = 1.3,
		   sigma.sq.mu = 0.5,
		   N = apply(y, 1, max, na.rm = TRUE))
tuning <- list(beta = 0.1, alpha = 0.1, beta.star = 0.1, alpha.star = 0.1,
               kappa = 0.2)
n.batch <- 50
batch.length <- 25
n.burn <- 750
n.thin <- 2
n.chains <- 2
abund.formula <- ~ (1 | abund.factor.1) + (1 | abund.factor.2) + (1 | abund.factor.3)
det.formula <- ~ 1 

data = data.list
n.batch = n.batch
batch.length = batch.length
inits = inits.list
priors = prior.list
accept.rate = 0.43
n.omp.threads = 1
verbose = TRUE
family <- 'NB'
n.report = 50
k.fold = 3
k.fold.threads = 3
k.fold.seed <- 100
k.fold.only <- FALSE


out <- NMix(abund.formula = abund.formula,
	    det.formula = det.formula,
	    data = data.list,
	    n.batch = n.batch,
	    batch.length = batch.length,
	    inits = inits.list,
	    family = 'NB',
	    tuning = tuning,
	    priors = prior.list,
	    accept.rate = 0.43,
	    n.omp.threads = 1,
	    verbose = FALSE,
	    n.report = 50,
	    n.burn = n.burn,
	    n.thin = n.thin,
	    n.chains = n.chains,
            k.fold = 2,
            k.fold.threads = 1)

# Test to make sure it worked ---------
test_that("out is of class NMix", {
  expect_s3_class(out, "NMix")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$rmspe), 1)
  expect_type(out$rmspe, "double")
  expect_gt(out$rmspe, 0)
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$pRE, FALSE)
  expect_equal(out$muRE, TRUE)
})

test_that("random effect gives error when non-numeric", {
  data.list$abund.covs <- as.data.frame(data.list$abund.covs)
  data.list$abund.covs$abund.factor.1 <- factor(data.list$abund.covs$abund.factor.1)
  expect_error(out <- NMix(abund.formula = abund.formula,
	    det.formula = det.formula,
	    data = data.list,
	    n.batch = n.batch,
	    batch.length = batch.length,
	    inits = inits.list,
	    family = 'NB',
	    tuning = tuning,
	    priors = prior.list,
	    accept.rate = 0.43,
	    n.omp.threads = 1,
	    verbose = FALSE,
	    n.report = 50,
	    n.burn = n.burn,
	    n.thin = n.thin,
	    n.chains = n.chains))
  data.list$abund.covs$abund.factor.1 <- factor(data.list$abund.covs$abund.factor.1)
  expect_error(out <- NMix(abund.formula = abund.formula,
	    det.formula = det.formula,
	    data = data.list,
	    n.batch = n.batch,
	    batch.length = batch.length,
	    inits = inits.list,
	    family = 'NB',
	    tuning = tuning,
	    priors = prior.list,
	    accept.rate = 0.43,
	    n.omp.threads = 1,
	    verbose = FALSE,
	    n.report = 50,
	    n.burn = n.burn,
	    n.thin = n.thin,
	    n.chains = n.chains))
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, dat$y)
})

# Check Poisson -----------------------
test_that("Poisson works", {
  out <- NMix(abund.formula = abund.formula,
  	    det.formula = det.formula,
  	    data = data.list,
  	    n.batch = n.batch,
  	    batch.length = batch.length,
  	    inits = inits.list,
  	    family = 'Poisson',
  	    tuning = tuning,
  	    priors = prior.list,
  	    accept.rate = 0.43,
  	    n.omp.threads = 1,
  	    verbose = FALSE,
  	    n.report = 50,
  	    n.burn = n.burn,
  	    n.thin = n.thin,
  	    n.chains = 1)
  expect_s3_class(out, "NMix")
})

# Check default priors ----------------
test_that("default priors, multiple chains, inits, burn, thin work", {
  out <- NMix(abund.formula = abund.formula,
	      det.formula = det.formula,
	      data = data.list,
	      n.batch = n.batch,
	      batch.length = batch.length,
	      family = 'NB',
	      accept.rate = 0.43,
	      n.omp.threads = 1,
	      verbose = FALSE,
	      n.report = 50,
	      n.burn = n.burn,
	      n.thin = n.thin,
	      n.chains = n.chains)
  expect_s3_class(out, 'NMix')
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {
  expect_output(NMix(abund.formula = abund.formula,
  	    det.formula = det.formula,
  	    data = data.list,
  	    n.batch = n.batch,
  	    batch.length = batch.length,
  	    inits = inits.list,
  	    family = 'Poisson',
  	    tuning = tuning,
  	    priors = prior.list,
  	    accept.rate = 0.43,
  	    n.omp.threads = 1,
  	    verbose = TRUE,
  	    n.report = 50,
  	    n.burn = n.burn,
  	    n.thin = n.thin,
  	    n.chains = 1))
})

# Check waicAbund -----------------------
test_that("waicAbund works for NMix", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicAbund(out, max(out$N.samples) + 10))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for NMix", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for NMix", {
  X.0 <- cbind(X, X.re)
  colnames(X.0) <- c('(Intercept)', 'abund.factor.1', 'abund.factor.2', 'abund.factor.3')
  pred.out <- predict(out, X.0)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$mu.0.samples), c(out$n.post * out$n.chains, J))
  expect_equal(dim(pred.out$N.0.samples), c(out$n.post * out$n.chains, J))
})
test_that("detection prediction works", {
  X.p.0 <- as.matrix(dat$X.p[, 1, ])
  pred.out <- predict(out, X.p.0, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, J))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for NMix", {
  n.post.samples <- out$n.post * out$n.chains
  ppc.out <- ppcAbund(out, 'chi-square', 2)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, max(n.rep)))
  
  ppc.out <- ppcAbund(out, 'chi-square', 1)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J))
  
  ppc.out <- ppcAbund(out, 'freeman-tukey', 1)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J))
  
  ppc.out <- ppcAbund(out, 'freeman-tukey', 2)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, max(n.rep)))

  ppc.out <- ppcAbund(out, 'chi-square', 0)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J, max(n.rep)))
  
  ppc.out <- ppcAbund(out, 'freeman-tukey', 0)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J, max(n.rep)))
})

# Abundance REs + covariates ----------------------------------------------
set.seed(1000)
J.x <- 15
J.y <- 15
J <- J.x * J.y
n.rep <- sample(2:4, J, replace = TRUE)
beta <- c(0, -0.2, 0.3)
p.abund <- length(beta)
alpha <- c(0.5)
p.det <- length(alpha)
mu.RE <- list(levels = c(15, 10, 12),
	       sigma.sq.mu = c(1.3, 0.2, 0.5))
p.RE <- list()
family <- 'NB'
kappa <- 0.5
dat <- simNMix(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, alpha = alpha,
	       mu.RE = mu.RE, p.RE = p.RE, sp = FALSE, family = family, kappa = kappa)

y <- dat$y
X <- dat$X
X.re <- dat$X.re
X.p <- dat$X.p
X.p.re <- dat$X.p.re

abund.covs <- cbind(X, X.re)
colnames(abund.covs) <- c('int', 'abund.cov.1', 'abund.cov.2', 
			  'abund.factor.1', 'abund.factor.2', 'abund.factor.3')
det.covs <- list(int = X.p[, , 1])
data.list <- list(y = y,
		  abund.covs = abund.covs,
 		  det.covs = det.covs)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 10),
		   alpha.normal = list(mean = 0,
			               var = 2.72))
# Starting values
inits.list <- list(alpha = 0,
		   beta = 0,
		   kappa = 1.3,
		   sigma.sq.mu = 0.5,
		   N = apply(y, 1, max, na.rm = TRUE))
tuning <- list(beta = 0.1, alpha = 0.1, beta.star = 0.1, alpha.star = 0.1,
               kappa = 0.2)
n.batch <- 50
batch.length <- 25
n.burn <- 750
n.thin <- 2
n.chains <- 2
abund.formula <- ~ (1 | abund.factor.1) + (1 | abund.factor.2) + (1 | abund.factor.3) + 
                   abund.cov.1 + abund.cov.2
det.formula <- ~ 1 

data = data.list
n.batch = n.batch
batch.length = batch.length
inits = inits.list
priors = prior.list
accept.rate = 0.43
n.omp.threads = 1
verbose = TRUE
family <- 'NB'
n.report = 50
k.fold = 3
k.fold.threads = 3
k.fold.seed <- 100
k.fold.only <- FALSE


out <- NMix(abund.formula = abund.formula,
	    det.formula = det.formula,
	    data = data.list,
	    n.batch = n.batch,
	    batch.length = batch.length,
	    inits = inits.list,
	    family = 'NB',
	    tuning = tuning,
	    priors = prior.list,
	    accept.rate = 0.43,
	    n.omp.threads = 1,
	    verbose = FALSE,
	    n.report = 50,
	    n.burn = n.burn,
	    n.thin = n.thin,
	    n.chains = n.chains,
            k.fold = 2,
            k.fold.threads = 1)

# Test to make sure it worked ---------
test_that("out is of class NMix", {
  expect_s3_class(out, "NMix")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$rmspe), 1)
  expect_type(out$rmspe, "double")
  expect_gt(out$rmspe, 0)
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$pRE, FALSE)
  expect_equal(out$muRE, TRUE)
})

test_that("random effect gives error when non-numeric", {
  data.list$abund.covs <- as.data.frame(data.list$abund.covs)
  data.list$abund.covs$abund.factor.1 <- factor(data.list$abund.covs$abund.factor.1)
  expect_error(out <- NMix(abund.formula = abund.formula,
	    det.formula = det.formula,
	    data = data.list,
	    n.batch = n.batch,
	    batch.length = batch.length,
	    inits = inits.list,
	    family = 'NB',
	    tuning = tuning,
	    priors = prior.list,
	    accept.rate = 0.43,
	    n.omp.threads = 1,
	    verbose = FALSE,
	    n.report = 50,
	    n.burn = n.burn,
	    n.thin = n.thin,
	    n.chains = n.chains))
  data.list$abund.covs$abund.factor.1 <- factor(data.list$abund.covs$abund.factor.1)
  expect_error(out <- NMix(abund.formula = abund.formula,
	    det.formula = det.formula,
	    data = data.list,
	    n.batch = n.batch,
	    batch.length = batch.length,
	    inits = inits.list,
	    family = 'NB',
	    tuning = tuning,
	    priors = prior.list,
	    accept.rate = 0.43,
	    n.omp.threads = 1,
	    verbose = FALSE,
	    n.report = 50,
	    n.burn = n.burn,
	    n.thin = n.thin,
	    n.chains = n.chains))
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, dat$y)
})

# Check Poisson -----------------------
test_that("Poisson works", {
  out <- NMix(abund.formula = abund.formula,
  	    det.formula = det.formula,
  	    data = data.list,
  	    n.batch = n.batch,
  	    batch.length = batch.length,
  	    inits = inits.list,
  	    family = 'Poisson',
  	    tuning = tuning,
  	    priors = prior.list,
  	    accept.rate = 0.43,
  	    n.omp.threads = 1,
  	    verbose = FALSE,
  	    n.report = 50,
  	    n.burn = n.burn,
  	    n.thin = n.thin,
  	    n.chains = 1)
  expect_s3_class(out, "NMix")
})

# Check default priors ----------------
test_that("default priors, multiple chains, inits, burn, thin work", {
  out <- NMix(abund.formula = abund.formula,
	      det.formula = det.formula,
	      data = data.list,
	      n.batch = n.batch,
	      batch.length = batch.length,
	      family = 'NB',
	      accept.rate = 0.43,
	      n.omp.threads = 1,
	      verbose = FALSE,
	      n.report = 50,
	      n.burn = n.burn,
	      n.thin = n.thin,
	      n.chains = n.chains)
  expect_s3_class(out, 'NMix')
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {
  expect_output(NMix(abund.formula = abund.formula,
  	    det.formula = det.formula,
  	    data = data.list,
  	    n.batch = n.batch,
  	    batch.length = batch.length,
  	    inits = inits.list,
  	    family = 'Poisson',
  	    tuning = tuning,
  	    priors = prior.list,
  	    accept.rate = 0.43,
  	    n.omp.threads = 1,
  	    verbose = TRUE,
  	    n.report = 50,
  	    n.burn = n.burn,
  	    n.thin = n.thin,
  	    n.chains = 1))
})

# Check waicAbund -----------------------
test_that("waicAbund works for NMix", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicAbund(out, max(out$N.samples) + 10))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for NMix", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for NMix", {
  X.0 <- cbind(X, X.re)
  colnames(X.0) <- c('(Intercept)', 'abund.cov.1', 'abund.cov.2',
		     'abund.factor.1', 'abund.factor.2', 'abund.factor.3')
  pred.out <- predict(out, X.0)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$mu.0.samples), c(out$n.post * out$n.chains, J))
  expect_equal(dim(pred.out$N.0.samples), c(out$n.post * out$n.chains, J))
})
test_that("detection prediction works", {
  X.p.0 <- as.matrix(dat$X.p[, 1, ])
  pred.out <- predict(out, X.p.0, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, J))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for NMix", {
  n.post.samples <- out$n.post * out$n.chains
  ppc.out <- ppcAbund(out, 'chi-square', 2)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, max(n.rep)))
  
  ppc.out <- ppcAbund(out, 'chi-square', 1)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J))
  
  ppc.out <- ppcAbund(out, 'freeman-tukey', 1)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J))
  
  ppc.out <- ppcAbund(out, 'freeman-tukey', 2)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, max(n.rep)))

  ppc.out <- ppcAbund(out, 'chi-square', 0)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J, max(n.rep)))
  
  ppc.out <- ppcAbund(out, 'freeman-tukey', 0)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J, max(n.rep)))
})

# Abundance REs + covariates in all ---------------------------------------
set.seed(1000)
J.x <- 15
J.y <- 15
J <- J.x * J.y
n.rep <- sample(2:4, J, replace = TRUE)
beta <- c(0, -0.2, 0.3)
p.abund <- length(beta)
alpha <- c(0.5, -0.2)
p.det <- length(alpha)
mu.RE <- list(levels = c(15, 10, 12),
	       sigma.sq.mu = c(1.3, 0.2, 0.5))
p.RE <- list()
family <- 'NB'
kappa <- 0.5
dat <- simNMix(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, alpha = alpha,
	       mu.RE = mu.RE, p.RE = p.RE, sp = FALSE, family = family, kappa = kappa)

y <- dat$y
X <- dat$X
X.re <- dat$X.re
X.p <- dat$X.p
X.p.re <- dat$X.p.re

abund.covs <- cbind(X, X.re)
colnames(abund.covs) <- c('int', 'abund.cov.1', 'abund.cov.2', 
			  'abund.factor.1', 'abund.factor.2', 'abund.factor.3')
det.covs <- list(int = X.p[, , 1], 
                 det.cov.1 = X.p[, , 2])
data.list <- list(y = y,
		  abund.covs = abund.covs,
 		  det.covs = det.covs)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 10),
		   alpha.normal = list(mean = 0,
			               var = 2.72))
# Starting values
inits.list <- list(alpha = 0,
		   beta = 0,
		   kappa = 1.3,
		   sigma.sq.mu = 0.5,
		   N = apply(y, 1, max, na.rm = TRUE))
tuning <- list(beta = 0.1, alpha = 0.1, beta.star = 0.1, alpha.star = 0.1,
               kappa = 0.2)
n.batch <- 50
batch.length <- 25
n.burn <- 750
n.thin <- 2
n.chains <- 2
abund.formula <- ~ (1 | abund.factor.1) + (1 | abund.factor.2) + (1 | abund.factor.3) + 
                   abund.cov.1 + abund.cov.2
det.formula <- ~ det.cov.1

data = data.list
n.batch = n.batch
batch.length = batch.length
inits = inits.list
priors = prior.list
accept.rate = 0.43
n.omp.threads = 1
verbose = TRUE
family <- 'NB'
n.report = 50
k.fold = 3
k.fold.threads = 3
k.fold.seed <- 100
k.fold.only <- FALSE


out <- NMix(abund.formula = abund.formula,
	    det.formula = det.formula,
	    data = data.list,
	    n.batch = n.batch,
	    batch.length = batch.length,
	    inits = inits.list,
	    family = 'NB',
	    tuning = tuning,
	    priors = prior.list,
	    accept.rate = 0.43,
	    n.omp.threads = 1,
	    verbose = FALSE,
	    n.report = 50,
	    n.burn = n.burn,
	    n.thin = n.thin,
	    n.chains = n.chains,
            k.fold = 2,
            k.fold.threads = 1)

# Test to make sure it worked ---------
test_that("out is of class NMix", {
  expect_s3_class(out, "NMix")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$rmspe), 1)
  expect_type(out$rmspe, "double")
  expect_gt(out$rmspe, 0)
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$pRE, FALSE)
  expect_equal(out$muRE, TRUE)
})

test_that("random effect gives error when non-numeric", {
  data.list$abund.covs <- as.data.frame(data.list$abund.covs)
  data.list$abund.covs$abund.factor.1 <- factor(data.list$abund.covs$abund.factor.1)
  expect_error(out <- NMix(abund.formula = abund.formula,
	    det.formula = det.formula,
	    data = data.list,
	    n.batch = n.batch,
	    batch.length = batch.length,
	    inits = inits.list,
	    family = 'NB',
	    tuning = tuning,
	    priors = prior.list,
	    accept.rate = 0.43,
	    n.omp.threads = 1,
	    verbose = FALSE,
	    n.report = 50,
	    n.burn = n.burn,
	    n.thin = n.thin,
	    n.chains = n.chains))
  data.list$abund.covs$abund.factor.1 <- factor(data.list$abund.covs$abund.factor.1)
  expect_error(out <- NMix(abund.formula = abund.formula,
	    det.formula = det.formula,
	    data = data.list,
	    n.batch = n.batch,
	    batch.length = batch.length,
	    inits = inits.list,
	    family = 'NB',
	    tuning = tuning,
	    priors = prior.list,
	    accept.rate = 0.43,
	    n.omp.threads = 1,
	    verbose = FALSE,
	    n.report = 50,
	    n.burn = n.burn,
	    n.thin = n.thin,
	    n.chains = n.chains))
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, dat$y)
})

# Check Poisson -----------------------
test_that("Poisson works", {
  out <- NMix(abund.formula = abund.formula,
  	    det.formula = det.formula,
  	    data = data.list,
  	    n.batch = n.batch,
  	    batch.length = batch.length,
  	    inits = inits.list,
  	    family = 'Poisson',
  	    tuning = tuning,
  	    priors = prior.list,
  	    accept.rate = 0.43,
  	    n.omp.threads = 1,
  	    verbose = FALSE,
  	    n.report = 50,
  	    n.burn = n.burn,
  	    n.thin = n.thin,
  	    n.chains = 1)
  expect_s3_class(out, "NMix")
})

# Check default priors ----------------
test_that("default priors, multiple chains, inits, burn, thin work", {
  out <- NMix(abund.formula = abund.formula,
	      det.formula = det.formula,
	      data = data.list,
	      n.batch = n.batch,
	      batch.length = batch.length,
	      family = 'NB',
	      accept.rate = 0.43,
	      n.omp.threads = 1,
	      verbose = FALSE,
	      n.report = 50,
	      n.burn = n.burn,
	      n.thin = n.thin,
	      n.chains = n.chains)
  expect_s3_class(out, 'NMix')
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {
  expect_output(NMix(abund.formula = abund.formula,
  	    det.formula = det.formula,
  	    data = data.list,
  	    n.batch = n.batch,
  	    batch.length = batch.length,
  	    inits = inits.list,
  	    family = 'Poisson',
  	    tuning = tuning,
  	    priors = prior.list,
  	    accept.rate = 0.43,
  	    n.omp.threads = 1,
  	    verbose = TRUE,
  	    n.report = 50,
  	    n.burn = n.burn,
  	    n.thin = n.thin,
  	    n.chains = 1))
})

# Check waicAbund -----------------------
test_that("waicAbund works for NMix", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicAbund(out, max(out$N.samples) + 10))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for NMix", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for NMix", {
  X.0 <- cbind(X, X.re)
  colnames(X.0) <- c('(Intercept)', 'abund.cov.1', 'abund.cov.2',
		     'abund.factor.1', 'abund.factor.2', 'abund.factor.3')
  pred.out <- predict(out, X.0)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$mu.0.samples), c(out$n.post * out$n.chains, J))
  expect_equal(dim(pred.out$N.0.samples), c(out$n.post * out$n.chains, J))
})
test_that("detection prediction works", {
  X.p.0 <- as.matrix(dat$X.p[, 1, ])
  pred.out <- predict(out, X.p.0, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, J))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for NMix", {
  n.post.samples <- out$n.post * out$n.chains
  ppc.out <- ppcAbund(out, 'chi-square', 2)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, max(n.rep)))
  
  ppc.out <- ppcAbund(out, 'chi-square', 1)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J))
  
  ppc.out <- ppcAbund(out, 'freeman-tukey', 1)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J))
  
  ppc.out <- ppcAbund(out, 'freeman-tukey', 2)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, max(n.rep)))

  ppc.out <- ppcAbund(out, 'chi-square', 0)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J, max(n.rep)))
  
  ppc.out <- ppcAbund(out, 'freeman-tukey', 0)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J, max(n.rep)))
})

# Random intercepts on detection ------------------------------------------
set.seed(1000)
J.x <- 15
J.y <- 15
J <- J.x * J.y
n.rep <- sample(2:4, J, replace = TRUE)
beta <- c(0)
p.abund <- length(beta)
alpha <- c(0.5)
p.det <- length(alpha)
mu.RE <- list()
p.RE <- list(levels = c(25), 
             sigma.sq.p = c(0.5))
family <- 'NB'
kappa <- 0.5
dat <- simNMix(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, alpha = alpha,
	       mu.RE = mu.RE, p.RE = p.RE, sp = FALSE, family = family, kappa = kappa)

y <- dat$y
X <- dat$X
X.re <- dat$X.re
X.p <- dat$X.p
X.p.re <- dat$X.p.re

abund.covs <- cbind(X)
colnames(abund.covs) <- c('int')
det.covs <- list(int = X.p[, , 1], 
                 det.factor.1 = X.p.re[, , 1])
data.list <- list(y = y,
		  abund.covs = abund.covs,
 		  det.covs = det.covs)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 10),
		   alpha.normal = list(mean = 0,
			               var = 2.72))
# Starting values
inits.list <- list(alpha = 0,
		   beta = 0,
		   kappa = 1.3,
		   sigma.sq.mu = 0.5,
		   N = apply(y, 1, max, na.rm = TRUE))
tuning <- list(beta = 0.1, alpha = 0.1, beta.star = 0.1, alpha.star = 0.1,
               kappa = 0.2)
n.batch <- 50
batch.length <- 25
n.burn <- 750
n.thin <- 2
n.chains <- 2
abund.formula <- ~ 1
det.formula <- ~ (1 | det.factor.1) 

data = data.list
n.batch = n.batch
batch.length = batch.length
inits = inits.list
priors = prior.list
accept.rate = 0.43
n.omp.threads = 1
verbose = TRUE
family <- 'NB'
n.report = 50
k.fold = 3
k.fold.threads = 3
k.fold.seed <- 100
k.fold.only <- FALSE


out <- NMix(abund.formula = abund.formula,
	    det.formula = det.formula,
	    data = data.list,
	    n.batch = n.batch,
	    batch.length = batch.length,
	    inits = inits.list,
	    family = 'NB',
	    tuning = tuning,
	    priors = prior.list,
	    accept.rate = 0.43,
	    n.omp.threads = 1,
	    verbose = FALSE,
	    n.report = 50,
	    n.burn = n.burn,
	    n.thin = n.thin,
	    n.chains = n.chains,
            k.fold = 2,
            k.fold.threads = 1)

# Test to make sure it worked ---------
test_that("out is of class NMix", {
  expect_s3_class(out, "NMix")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$rmspe), 1)
  expect_type(out$rmspe, "double")
  expect_gt(out$rmspe, 0)
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$pRE, TRUE)
  expect_equal(out$muRE, FALSE)
})

test_that("random effect gives error when non-numeric", {
  data.list$det.covs$det.factor.1 <- factor(data.list$det.covs$det.factor.1)
  expect_error(out <- NMix(abund.formula = abund.formula,
	    det.formula = det.formula,
	    data = data.list,
	    n.batch = n.batch,
	    batch.length = batch.length,
	    inits = inits.list,
	    family = 'NB',
	    tuning = tuning,
	    priors = prior.list,
	    accept.rate = 0.43,
	    n.omp.threads = 1,
	    verbose = FALSE,
	    n.report = 50,
	    n.burn = n.burn,
	    n.thin = n.thin,
	    n.chains = n.chains))
  data.list$det.covs$det.factor.1 <- as.character(factor(data.list$det.covs$det.factor.1))
  expect_error(out <- NMix(abund.formula = abund.formula,
	    det.formula = det.formula,
	    data = data.list,
	    n.batch = n.batch,
	    batch.length = batch.length,
	    inits = inits.list,
	    family = 'NB',
	    tuning = tuning,
	    priors = prior.list,
	    accept.rate = 0.43,
	    n.omp.threads = 1,
	    verbose = FALSE,
	    n.report = 50,
	    n.burn = n.burn,
	    n.thin = n.thin,
	    n.chains = n.chains))
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, dat$y)
})

# Check Poisson -----------------------
test_that("Poisson works", {
  out <- NMix(abund.formula = abund.formula,
  	    det.formula = det.formula,
  	    data = data.list,
  	    n.batch = n.batch,
  	    batch.length = batch.length,
  	    inits = inits.list,
  	    family = 'Poisson',
  	    tuning = tuning,
  	    priors = prior.list,
  	    accept.rate = 0.43,
  	    n.omp.threads = 1,
  	    verbose = FALSE,
  	    n.report = 50,
  	    n.burn = n.burn,
  	    n.thin = n.thin,
  	    n.chains = 1)
  expect_s3_class(out, "NMix")
})

# Check default priors ----------------
test_that("default priors, multiple chains, inits, burn, thin work", {
  out <- NMix(abund.formula = abund.formula,
	      det.formula = det.formula,
	      data = data.list,
	      n.batch = n.batch,
	      batch.length = batch.length,
	      family = 'NB',
	      accept.rate = 0.43,
	      n.omp.threads = 1,
	      verbose = FALSE,
	      n.report = 50,
	      n.burn = n.burn,
	      n.thin = n.thin,
	      n.chains = n.chains)
  expect_s3_class(out, 'NMix')
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {
  expect_output(NMix(abund.formula = abund.formula,
  	    det.formula = det.formula,
  	    data = data.list,
  	    n.batch = n.batch,
  	    batch.length = batch.length,
  	    inits = inits.list,
  	    family = 'Poisson',
  	    tuning = tuning,
  	    priors = prior.list,
  	    accept.rate = 0.43,
  	    n.omp.threads = 1,
  	    verbose = TRUE,
  	    n.report = 50,
  	    n.burn = n.burn,
  	    n.thin = n.thin,
  	    n.chains = 1))
})

# Check waicAbund -----------------------
test_that("waicAbund works for NMix", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicAbund(out, max(out$N.samples) + 10))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for NMix", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for NMix", {
  X.0 <- cbind(X)
  pred.out <- predict(out, X.0)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$mu.0.samples), c(out$n.post * out$n.chains, J))
  expect_equal(dim(pred.out$N.0.samples), c(out$n.post * out$n.chains, J))
})
test_that("detection prediction works", {
  X.p.0 <- cbind(dat$X.p[, 1, ], dat$X.p.re[, 1, 1])
  colnames(X.p.0) <- c('(Intercept)', 'det.factor.1')
  pred.out <- predict(out, X.p.0, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, J))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for NMix", {
  n.post.samples <- out$n.post * out$n.chains
  ppc.out <- ppcAbund(out, 'chi-square', 2)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, max(n.rep)))
  
  ppc.out <- ppcAbund(out, 'chi-square', 1)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J))
  
  ppc.out <- ppcAbund(out, 'freeman-tukey', 1)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J))
  
  ppc.out <- ppcAbund(out, 'freeman-tukey', 2)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, max(n.rep)))

  ppc.out <- ppcAbund(out, 'chi-square', 0)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J, max(n.rep)))
  
  ppc.out <- ppcAbund(out, 'freeman-tukey', 0)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J, max(n.rep)))
})

# Multiple random intercepts on detection ---------------------------------
set.seed(1000)
J.x <- 15
J.y <- 15
J <- J.x * J.y
n.rep <- sample(2:4, J, replace = TRUE)
beta <- c(0)
p.abund <- length(beta)
alpha <- c(0.5)
p.det <- length(alpha)
mu.RE <- list()
p.RE <- list(levels = c(25, 20), 
             sigma.sq.p = c(0.5, 0.75))
family <- 'NB'
kappa <- 0.5
dat <- simNMix(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, alpha = alpha,
	       mu.RE = mu.RE, p.RE = p.RE, sp = FALSE, family = family, kappa = kappa)

y <- dat$y
X <- dat$X
X.re <- dat$X.re
X.p <- dat$X.p
X.p.re <- dat$X.p.re

abund.covs <- cbind(X)
colnames(abund.covs) <- c('int')
det.covs <- list(int = X.p[, , 1], 
                 det.factor.1 = X.p.re[, , 1], 
                 det.factor.2 = X.p.re[, , 2])
data.list <- list(y = y,
		  abund.covs = abund.covs,
 		  det.covs = det.covs)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 10),
		   alpha.normal = list(mean = 0,
			               var = 2.72))
# Starting values
inits.list <- list(alpha = 0,
		   beta = 0,
		   kappa = 1.3,
		   sigma.sq.mu = 0.5,
		   N = apply(y, 1, max, na.rm = TRUE))
tuning <- list(beta = 0.1, alpha = 0.1, beta.star = 0.1, alpha.star = 0.1,
               kappa = 0.2)
n.batch <- 50
batch.length <- 25
n.burn <- 750
n.thin <- 2
n.chains <- 2
abund.formula <- ~ 1
det.formula <- ~ (1 | det.factor.1) + (1 | det.factor.2)

data = data.list
n.batch = n.batch
batch.length = batch.length
inits = inits.list
priors = prior.list
accept.rate = 0.43
n.omp.threads = 1
verbose = TRUE
family <- 'NB'
n.report = 50
k.fold = 3
k.fold.threads = 3
k.fold.seed <- 100
k.fold.only <- FALSE


out <- NMix(abund.formula = abund.formula,
	    det.formula = det.formula,
	    data = data.list,
	    n.batch = n.batch,
	    batch.length = batch.length,
	    inits = inits.list,
	    family = 'NB',
	    tuning = tuning,
	    priors = prior.list,
	    accept.rate = 0.43,
	    n.omp.threads = 1,
	    verbose = FALSE,
	    n.report = 50,
	    n.burn = n.burn,
	    n.thin = n.thin,
	    n.chains = n.chains,
            k.fold = 2,
            k.fold.threads = 1)

# Test to make sure it worked ---------
test_that("out is of class NMix", {
  expect_s3_class(out, "NMix")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$rmspe), 1)
  expect_type(out$rmspe, "double")
  expect_gt(out$rmspe, 0)
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$pRE, TRUE)
  expect_equal(out$muRE, FALSE)
})

test_that("random effect gives error when non-numeric", {
  data.list$det.covs$det.factor.1 <- factor(data.list$det.covs$det.factor.1)
  expect_error(out <- NMix(abund.formula = abund.formula,
	    det.formula = det.formula,
	    data = data.list,
	    n.batch = n.batch,
	    batch.length = batch.length,
	    inits = inits.list,
	    family = 'NB',
	    tuning = tuning,
	    priors = prior.list,
	    accept.rate = 0.43,
	    n.omp.threads = 1,
	    verbose = FALSE,
	    n.report = 50,
	    n.burn = n.burn,
	    n.thin = n.thin,
	    n.chains = n.chains))
  data.list$det.covs$det.factor.1 <- as.character(factor(data.list$det.covs$det.factor.1))
  expect_error(out <- NMix(abund.formula = abund.formula,
	    det.formula = det.formula,
	    data = data.list,
	    n.batch = n.batch,
	    batch.length = batch.length,
	    inits = inits.list,
	    family = 'NB',
	    tuning = tuning,
	    priors = prior.list,
	    accept.rate = 0.43,
	    n.omp.threads = 1,
	    verbose = FALSE,
	    n.report = 50,
	    n.burn = n.burn,
	    n.thin = n.thin,
	    n.chains = n.chains))
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, dat$y)
})

# Check Poisson -----------------------
test_that("Poisson works", {
  out <- NMix(abund.formula = abund.formula,
  	    det.formula = det.formula,
  	    data = data.list,
  	    n.batch = n.batch,
  	    batch.length = batch.length,
  	    inits = inits.list,
  	    family = 'Poisson',
  	    tuning = tuning,
  	    priors = prior.list,
  	    accept.rate = 0.43,
  	    n.omp.threads = 1,
  	    verbose = FALSE,
  	    n.report = 50,
  	    n.burn = n.burn,
  	    n.thin = n.thin,
  	    n.chains = 1)
  expect_s3_class(out, "NMix")
})

# Check default priors ----------------
test_that("default priors, multiple chains, inits, burn, thin work", {
  out <- NMix(abund.formula = abund.formula,
	      det.formula = det.formula,
	      data = data.list,
	      n.batch = n.batch,
	      batch.length = batch.length,
	      family = 'NB',
	      accept.rate = 0.43,
	      n.omp.threads = 1,
	      verbose = FALSE,
	      n.report = 50,
	      n.burn = n.burn,
	      n.thin = n.thin,
	      n.chains = n.chains)
  expect_s3_class(out, 'NMix')
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {
  expect_output(NMix(abund.formula = abund.formula,
  	    det.formula = det.formula,
  	    data = data.list,
  	    n.batch = n.batch,
  	    batch.length = batch.length,
  	    inits = inits.list,
  	    family = 'Poisson',
  	    tuning = tuning,
  	    priors = prior.list,
  	    accept.rate = 0.43,
  	    n.omp.threads = 1,
  	    verbose = TRUE,
  	    n.report = 50,
  	    n.burn = n.burn,
  	    n.thin = n.thin,
  	    n.chains = 1))
})

# Check waicAbund -----------------------
test_that("waicAbund works for NMix", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicAbund(out, max(out$N.samples) + 10))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for NMix", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for NMix", {
  X.0 <- cbind(X)
  pred.out <- predict(out, X.0)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$mu.0.samples), c(out$n.post * out$n.chains, J))
  expect_equal(dim(pred.out$N.0.samples), c(out$n.post * out$n.chains, J))
})
test_that("detection prediction works", {
  X.p.0 <- cbind(dat$X.p[, 1, ], dat$X.p.re[, 1, ])
  colnames(X.p.0) <- c('(Intercept)', 'det.factor.1', 'det.factor.2')
  pred.out <- predict(out, X.p.0, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, J))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for NMix", {
  n.post.samples <- out$n.post * out$n.chains
  ppc.out <- ppcAbund(out, 'chi-square', 2)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, max(n.rep)))
  
  ppc.out <- ppcAbund(out, 'chi-square', 1)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J))
  
  ppc.out <- ppcAbund(out, 'freeman-tukey', 1)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J))
  
  ppc.out <- ppcAbund(out, 'freeman-tukey', 2)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, max(n.rep)))

  ppc.out <- ppcAbund(out, 'chi-square', 0)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J, max(n.rep)))
  
  ppc.out <- ppcAbund(out, 'freeman-tukey', 0)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J, max(n.rep)))
})

# Multiple random intercepts with covariates ------------------------------
set.seed(1000)
J.x <- 15
J.y <- 15
J <- J.x * J.y
n.rep <- sample(2:4, J, replace = TRUE)
beta <- c(0)
p.abund <- length(beta)
alpha <- c(0.5, -0.2)
p.det <- length(alpha)
mu.RE <- list()
p.RE <- list(levels = c(25, 20), 
             sigma.sq.p = c(0.5, 0.75))
family <- 'NB'
kappa <- 0.5
dat <- simNMix(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, alpha = alpha,
	       mu.RE = mu.RE, p.RE = p.RE, sp = FALSE, family = family, kappa = kappa)

y <- dat$y
X <- dat$X
X.re <- dat$X.re
X.p <- dat$X.p
X.p.re <- dat$X.p.re

abund.covs <- cbind(X)
colnames(abund.covs) <- c('int')
det.covs <- list(int = X.p[, , 1], 
		 det.cov.1 = X.p[, , 2],
                 det.factor.1 = X.p.re[, , 1], 
                 det.factor.2 = X.p.re[, , 2])
data.list <- list(y = y,
		  abund.covs = abund.covs,
 		  det.covs = det.covs)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 10),
		   alpha.normal = list(mean = 0,
			               var = 2.72))
# Starting values
inits.list <- list(alpha = 0,
		   beta = 0,
		   kappa = 1.3,
		   sigma.sq.mu = 0.5,
		   N = apply(y, 1, max, na.rm = TRUE))
tuning <- list(beta = 0.1, alpha = 0.1, beta.star = 0.1, alpha.star = 0.1,
               kappa = 0.2)
n.batch <- 50
batch.length <- 25
n.burn <- 750
n.thin <- 2
n.chains <- 2
abund.formula <- ~ 1
det.formula <- ~ det.cov.1 + (1 | det.factor.1) + (1 | det.factor.2)

data = data.list
n.batch = n.batch
batch.length = batch.length
inits = inits.list
priors = prior.list
accept.rate = 0.43
n.omp.threads = 1
verbose = TRUE
family <- 'NB'
n.report = 50
k.fold = 3
k.fold.threads = 3
k.fold.seed <- 100
k.fold.only <- FALSE


out <- NMix(abund.formula = abund.formula,
	    det.formula = det.formula,
	    data = data.list,
	    n.batch = n.batch,
	    batch.length = batch.length,
	    inits = inits.list,
	    family = 'NB',
	    tuning = tuning,
	    priors = prior.list,
	    accept.rate = 0.43,
	    n.omp.threads = 1,
	    verbose = FALSE,
	    n.report = 50,
	    n.burn = n.burn,
	    n.thin = n.thin,
	    n.chains = n.chains,
            k.fold = 2,
            k.fold.threads = 1)

# Test to make sure it worked ---------
test_that("out is of class NMix", {
  expect_s3_class(out, "NMix")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$rmspe), 1)
  expect_type(out$rmspe, "double")
  expect_gt(out$rmspe, 0)
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$pRE, TRUE)
  expect_equal(out$muRE, FALSE)
})

test_that("random effect gives error when non-numeric", {
  data.list$det.covs$det.factor.1 <- factor(data.list$det.covs$det.factor.1)
  expect_error(out <- NMix(abund.formula = abund.formula,
	    det.formula = det.formula,
	    data = data.list,
	    n.batch = n.batch,
	    batch.length = batch.length,
	    inits = inits.list,
	    family = 'NB',
	    tuning = tuning,
	    priors = prior.list,
	    accept.rate = 0.43,
	    n.omp.threads = 1,
	    verbose = FALSE,
	    n.report = 50,
	    n.burn = n.burn,
	    n.thin = n.thin,
	    n.chains = n.chains))
  data.list$det.covs$det.factor.1 <- as.character(factor(data.list$det.covs$det.factor.1))
  expect_error(out <- NMix(abund.formula = abund.formula,
	    det.formula = det.formula,
	    data = data.list,
	    n.batch = n.batch,
	    batch.length = batch.length,
	    inits = inits.list,
	    family = 'NB',
	    tuning = tuning,
	    priors = prior.list,
	    accept.rate = 0.43,
	    n.omp.threads = 1,
	    verbose = FALSE,
	    n.report = 50,
	    n.burn = n.burn,
	    n.thin = n.thin,
	    n.chains = n.chains))
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, dat$y)
})

# Check Poisson -----------------------
test_that("Poisson works", {
  out <- NMix(abund.formula = abund.formula,
  	    det.formula = det.formula,
  	    data = data.list,
  	    n.batch = n.batch,
  	    batch.length = batch.length,
  	    inits = inits.list,
  	    family = 'Poisson',
  	    tuning = tuning,
  	    priors = prior.list,
  	    accept.rate = 0.43,
  	    n.omp.threads = 1,
  	    verbose = FALSE,
  	    n.report = 50,
  	    n.burn = n.burn,
  	    n.thin = n.thin,
  	    n.chains = 1)
  expect_s3_class(out, "NMix")
})

# Check default priors ----------------
test_that("default priors, multiple chains, inits, burn, thin work", {
  out <- NMix(abund.formula = abund.formula,
	      det.formula = det.formula,
	      data = data.list,
	      n.batch = n.batch,
	      batch.length = batch.length,
	      family = 'NB',
	      accept.rate = 0.43,
	      n.omp.threads = 1,
	      verbose = FALSE,
	      n.report = 50,
	      n.burn = n.burn,
	      n.thin = n.thin,
	      n.chains = n.chains)
  expect_s3_class(out, 'NMix')
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {
  expect_output(NMix(abund.formula = abund.formula,
  	    det.formula = det.formula,
  	    data = data.list,
  	    n.batch = n.batch,
  	    batch.length = batch.length,
  	    inits = inits.list,
  	    family = 'Poisson',
  	    tuning = tuning,
  	    priors = prior.list,
  	    accept.rate = 0.43,
  	    n.omp.threads = 1,
  	    verbose = TRUE,
  	    n.report = 50,
  	    n.burn = n.burn,
  	    n.thin = n.thin,
  	    n.chains = 1))
})

# Check waicAbund -----------------------
test_that("waicAbund works for NMix", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicAbund(out, max(out$N.samples) + 10))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for NMix", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for NMix", {
  X.0 <- cbind(X)
  pred.out <- predict(out, X.0)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$mu.0.samples), c(out$n.post * out$n.chains, J))
  expect_equal(dim(pred.out$N.0.samples), c(out$n.post * out$n.chains, J))
})
test_that("detection prediction works", {
  X.p.0 <- cbind(dat$X.p[, 1, ], dat$X.p.re[, 1, ])
  colnames(X.p.0) <- c('(Intercept)', 'det.cov.1', 'det.factor.1', 'det.factor.2')
  pred.out <- predict(out, X.p.0, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, J))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for NMix", {
  n.post.samples <- out$n.post * out$n.chains
  ppc.out <- ppcAbund(out, 'chi-square', 2)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, max(n.rep)))
  
  ppc.out <- ppcAbund(out, 'chi-square', 1)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J))
  
  ppc.out <- ppcAbund(out, 'freeman-tukey', 1)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J))
  
  ppc.out <- ppcAbund(out, 'freeman-tukey', 2)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, max(n.rep)))

  ppc.out <- ppcAbund(out, 'chi-square', 0)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J, max(n.rep)))
  
  ppc.out <- ppcAbund(out, 'freeman-tukey', 0)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J, max(n.rep)))
})

# Multiple random intercepts with covariates on both ----------------------
set.seed(1000)
J.x <- 15
J.y <- 15
J <- J.x * J.y
n.rep <- sample(2:4, J, replace = TRUE)
beta <- c(0, -0.2, 0.3)
p.abund <- length(beta)
alpha <- c(0.5, -0.2)
p.det <- length(alpha)
mu.RE <- list()
p.RE <- list(levels = c(25, 20), 
             sigma.sq.p = c(0.5, 0.75))
family <- 'NB'
kappa <- 0.5
dat <- simNMix(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, alpha = alpha,
	       mu.RE = mu.RE, p.RE = p.RE, sp = FALSE, family = family, kappa = kappa)

y <- dat$y
X <- dat$X
X.re <- dat$X.re
X.p <- dat$X.p
X.p.re <- dat$X.p.re

abund.covs <- cbind(X)
colnames(abund.covs) <- c('int', 'abund.cov.1', 'abund.cov.2')
det.covs <- list(int = X.p[, , 1], 
		 det.cov.1 = X.p[, , 2],
                 det.factor.1 = X.p.re[, , 1], 
                 det.factor.2 = X.p.re[, , 2])
data.list <- list(y = y,
		  abund.covs = abund.covs,
 		  det.covs = det.covs)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 10),
		   alpha.normal = list(mean = 0,
			               var = 2.72))
# Starting values
inits.list <- list(alpha = 0,
		   beta = 0,
		   kappa = 1.3,
		   sigma.sq.mu = 0.5,
		   N = apply(y, 1, max, na.rm = TRUE))
tuning <- list(beta = 0.1, alpha = 0.1, beta.star = 0.1, alpha.star = 0.1,
               kappa = 0.2)
n.batch <- 50
batch.length <- 25
n.burn <- 750
n.thin <- 2
n.chains <- 2
abund.formula <- ~ abund.cov.1 + abund.cov.2
det.formula <- ~ det.cov.1 + (1 | det.factor.1) + (1 | det.factor.2)

data = data.list
n.batch = n.batch
batch.length = batch.length
inits = inits.list
priors = prior.list
accept.rate = 0.43
n.omp.threads = 1
verbose = TRUE
family <- 'NB'
n.report = 50
k.fold = 3
k.fold.threads = 3
k.fold.seed <- 100
k.fold.only <- FALSE


out <- NMix(abund.formula = abund.formula,
	    det.formula = det.formula,
	    data = data.list,
	    n.batch = n.batch,
	    batch.length = batch.length,
	    inits = inits.list,
	    family = 'NB',
	    tuning = tuning,
	    priors = prior.list,
	    accept.rate = 0.43,
	    n.omp.threads = 1,
	    verbose = FALSE,
	    n.report = 50,
	    n.burn = n.burn,
	    n.thin = n.thin,
	    n.chains = n.chains,
            k.fold = 2,
            k.fold.threads = 1)

# Test to make sure it worked ---------
test_that("out is of class NMix", {
  expect_s3_class(out, "NMix")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$rmspe), 1)
  expect_type(out$rmspe, "double")
  expect_gt(out$rmspe, 0)
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$pRE, TRUE)
  expect_equal(out$muRE, FALSE)
})

test_that("random effect gives error when non-numeric", {
  data.list$det.covs$det.factor.1 <- factor(data.list$det.covs$det.factor.1)
  expect_error(out <- NMix(abund.formula = abund.formula,
	    det.formula = det.formula,
	    data = data.list,
	    n.batch = n.batch,
	    batch.length = batch.length,
	    inits = inits.list,
	    family = 'NB',
	    tuning = tuning,
	    priors = prior.list,
	    accept.rate = 0.43,
	    n.omp.threads = 1,
	    verbose = FALSE,
	    n.report = 50,
	    n.burn = n.burn,
	    n.thin = n.thin,
	    n.chains = n.chains))
  data.list$det.covs$det.factor.1 <- as.character(factor(data.list$det.covs$det.factor.1))
  expect_error(out <- NMix(abund.formula = abund.formula,
	    det.formula = det.formula,
	    data = data.list,
	    n.batch = n.batch,
	    batch.length = batch.length,
	    inits = inits.list,
	    family = 'NB',
	    tuning = tuning,
	    priors = prior.list,
	    accept.rate = 0.43,
	    n.omp.threads = 1,
	    verbose = FALSE,
	    n.report = 50,
	    n.burn = n.burn,
	    n.thin = n.thin,
	    n.chains = n.chains))
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, dat$y)
})

# Check Poisson -----------------------
test_that("Poisson works", {
  out <- NMix(abund.formula = abund.formula,
  	    det.formula = det.formula,
  	    data = data.list,
  	    n.batch = n.batch,
  	    batch.length = batch.length,
  	    inits = inits.list,
  	    family = 'Poisson',
  	    tuning = tuning,
  	    priors = prior.list,
  	    accept.rate = 0.43,
  	    n.omp.threads = 1,
  	    verbose = FALSE,
  	    n.report = 50,
  	    n.burn = n.burn,
  	    n.thin = n.thin,
  	    n.chains = 1)
  expect_s3_class(out, "NMix")
})

# Check default priors ----------------
test_that("default priors, multiple chains, inits, burn, thin work", {
  out <- NMix(abund.formula = abund.formula,
	      det.formula = det.formula,
	      data = data.list,
	      n.batch = n.batch,
	      batch.length = batch.length,
	      family = 'NB',
	      accept.rate = 0.43,
	      n.omp.threads = 1,
	      verbose = FALSE,
	      n.report = 50,
	      n.burn = n.burn,
	      n.thin = n.thin,
	      n.chains = n.chains)
  expect_s3_class(out, 'NMix')
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {
  expect_output(NMix(abund.formula = abund.formula,
  	    det.formula = det.formula,
  	    data = data.list,
  	    n.batch = n.batch,
  	    batch.length = batch.length,
  	    inits = inits.list,
  	    family = 'Poisson',
  	    tuning = tuning,
  	    priors = prior.list,
  	    accept.rate = 0.43,
  	    n.omp.threads = 1,
  	    verbose = TRUE,
  	    n.report = 50,
  	    n.burn = n.burn,
  	    n.thin = n.thin,
  	    n.chains = 1))
})

# Check waicAbund -----------------------
test_that("waicAbund works for NMix", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicAbund(out, max(out$N.samples) + 10))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for NMix", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for NMix", {
  X.0 <- cbind(X)
  colnames(X.0) <- c('(Intercept)', 'abund.cov.1', 'abund.cov.2')
  pred.out <- predict(out, X.0)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$mu.0.samples), c(out$n.post * out$n.chains, J))
  expect_equal(dim(pred.out$N.0.samples), c(out$n.post * out$n.chains, J))
})
test_that("detection prediction works", {
  X.p.0 <- cbind(dat$X.p[, 1, ], dat$X.p.re[, 1, ])
  colnames(X.p.0) <- c('(Intercept)', 'det.cov.1', 'det.factor.1', 'det.factor.2')
  pred.out <- predict(out, X.p.0, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, J))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for NMix", {
  n.post.samples <- out$n.post * out$n.chains
  ppc.out <- ppcAbund(out, 'chi-square', 2)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, max(n.rep)))
  
  ppc.out <- ppcAbund(out, 'chi-square', 1)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J))
  
  ppc.out <- ppcAbund(out, 'freeman-tukey', 1)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J))
  
  ppc.out <- ppcAbund(out, 'freeman-tukey', 2)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, max(n.rep)))

  ppc.out <- ppcAbund(out, 'chi-square', 0)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J, max(n.rep)))
  
  ppc.out <- ppcAbund(out, 'freeman-tukey', 0)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J, max(n.rep)))
})

# Multiple random intercepts and covariates -------------------------------
set.seed(1000)
J.x <- 15
J.y <- 15
J <- J.x * J.y
n.rep <- sample(2:4, J, replace = TRUE)
beta <- c(0, -0.2, 0.3)
p.abund <- length(beta)
alpha <- c(0.5, -0.2)
p.det <- length(alpha)
mu.RE <- list(levels = c(10, 15), 
              sigma.sq.mu = c(0.3, 0.9))
p.RE <- list(levels = c(25, 20), 
             sigma.sq.p = c(0.5, 0.75))
family <- 'NB'
kappa <- 0.5
dat <- simNMix(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, alpha = alpha,
	       mu.RE = mu.RE, p.RE = p.RE, sp = FALSE, family = family, kappa = kappa)

y <- dat$y
X <- dat$X
X.re <- dat$X.re
X.p <- dat$X.p
X.p.re <- dat$X.p.re

abund.covs <- cbind(X, X.re)
colnames(abund.covs) <- c('int', 'abund.cov.1', 'abund.cov.2', 'abund.factor.1', 
                          'abund.factor.2')
det.covs <- list(int = X.p[, , 1], 
		 det.cov.1 = X.p[, , 2],
                 det.factor.1 = X.p.re[, , 1], 
                 det.factor.2 = X.p.re[, , 2])
data.list <- list(y = y,
		  abund.covs = abund.covs,
 		  det.covs = det.covs)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 10),
		   alpha.normal = list(mean = 0,
			               var = 2.72))
# Starting values
inits.list <- list(alpha = 0,
		   beta = 0,
		   kappa = 1.3,
		   sigma.sq.mu = 0.5,
		   N = apply(y, 1, max, na.rm = TRUE))
tuning <- list(beta = 0.1, alpha = 0.1, beta.star = 0.1, alpha.star = 0.1,
               kappa = 0.2)
n.batch <- 50
batch.length <- 25
n.burn <- 750
n.thin <- 2
n.chains <- 2
abund.formula <- ~ abund.cov.1 + abund.cov.2 + (1 | abund.factor.1) + (1 | abund.factor.2)
det.formula <- ~ det.cov.1 + (1 | det.factor.1) + (1 | det.factor.2)

data = data.list
n.batch = n.batch
batch.length = batch.length
inits = inits.list
priors = prior.list
accept.rate = 0.43
n.omp.threads = 1
verbose = TRUE
family <- 'NB'
n.report = 50
k.fold = 3
k.fold.threads = 3
k.fold.seed <- 100
k.fold.only <- FALSE


out <- NMix(abund.formula = abund.formula,
	    det.formula = det.formula,
	    data = data.list,
	    n.batch = n.batch,
	    batch.length = batch.length,
	    inits = inits.list,
	    family = 'NB',
	    tuning = tuning,
	    priors = prior.list,
	    accept.rate = 0.43,
	    n.omp.threads = 1,
	    verbose = FALSE,
	    n.report = 50,
	    n.burn = n.burn,
	    n.thin = n.thin,
	    n.chains = n.chains,
            k.fold = 2,
            k.fold.threads = 1)

# Test to make sure it worked ---------
test_that("out is of class NMix", {
  expect_s3_class(out, "NMix")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$rmspe), 1)
  expect_type(out$rmspe, "double")
  expect_gt(out$rmspe, 0)
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$pRE, TRUE)
  expect_equal(out$muRE, TRUE)
})

test_that("random effect gives error when non-numeric", {
  data.list$det.covs$det.factor.1 <- factor(data.list$det.covs$det.factor.1)
  expect_error(out <- NMix(abund.formula = abund.formula,
	    det.formula = det.formula,
	    data = data.list,
	    n.batch = n.batch,
	    batch.length = batch.length,
	    inits = inits.list,
	    family = 'NB',
	    tuning = tuning,
	    priors = prior.list,
	    accept.rate = 0.43,
	    n.omp.threads = 1,
	    verbose = FALSE,
	    n.report = 50,
	    n.burn = n.burn,
	    n.thin = n.thin,
	    n.chains = n.chains))
  data.list$det.covs$det.factor.1 <- as.character(factor(data.list$det.covs$det.factor.1))
  expect_error(out <- NMix(abund.formula = abund.formula,
	    det.formula = det.formula,
	    data = data.list,
	    n.batch = n.batch,
	    batch.length = batch.length,
	    inits = inits.list,
	    family = 'NB',
	    tuning = tuning,
	    priors = prior.list,
	    accept.rate = 0.43,
	    n.omp.threads = 1,
	    verbose = FALSE,
	    n.report = 50,
	    n.burn = n.burn,
	    n.thin = n.thin,
	    n.chains = n.chains))
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, dat$y)
})

# Check Poisson -----------------------
test_that("Poisson works", {
  out <- NMix(abund.formula = abund.formula,
  	    det.formula = det.formula,
  	    data = data.list,
  	    n.batch = n.batch,
  	    batch.length = batch.length,
  	    inits = inits.list,
  	    family = 'Poisson',
  	    tuning = tuning,
  	    priors = prior.list,
  	    accept.rate = 0.43,
  	    n.omp.threads = 1,
  	    verbose = FALSE,
  	    n.report = 50,
  	    n.burn = n.burn,
  	    n.thin = n.thin,
  	    n.chains = 1)
  expect_s3_class(out, "NMix")
})

# Check default priors ----------------
test_that("default priors, multiple chains, inits, burn, thin work", {
  out <- NMix(abund.formula = abund.formula,
	      det.formula = det.formula,
	      data = data.list,
	      n.batch = n.batch,
	      batch.length = batch.length,
	      family = 'NB',
	      accept.rate = 0.43,
	      n.omp.threads = 1,
	      verbose = FALSE,
	      n.report = 50,
	      n.burn = n.burn,
	      n.thin = n.thin,
	      n.chains = n.chains)
  expect_s3_class(out, 'NMix')
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {
  expect_output(NMix(abund.formula = abund.formula,
  	    det.formula = det.formula,
  	    data = data.list,
  	    n.batch = n.batch,
  	    batch.length = batch.length,
  	    inits = inits.list,
  	    family = 'Poisson',
  	    tuning = tuning,
  	    priors = prior.list,
  	    accept.rate = 0.43,
  	    n.omp.threads = 1,
  	    verbose = TRUE,
  	    n.report = 50,
  	    n.burn = n.burn,
  	    n.thin = n.thin,
  	    n.chains = 1))
})

# Check waicAbund -----------------------
test_that("waicAbund works for NMix", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicAbund(out, max(out$N.samples) + 10))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for NMix", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for NMix", {
  X.0 <- cbind(X, X.re)
  colnames(X.0) <- c('(Intercept)', 'abund.cov.1', 'abund.cov.2', 'abund.factor.1', 'abund.factor.2')
  pred.out <- predict(out, X.0)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$mu.0.samples), c(out$n.post * out$n.chains, J))
  expect_equal(dim(pred.out$N.0.samples), c(out$n.post * out$n.chains, J))
})
test_that("detection prediction works", {
  X.p.0 <- cbind(dat$X.p[, 1, ], dat$X.p.re[, 1, ])
  colnames(X.p.0) <- c('(Intercept)', 'det.cov.1', 'det.factor.1', 'det.factor.2')
  pred.out <- predict(out, X.p.0, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, J))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for NMix", {
  n.post.samples <- out$n.post * out$n.chains
  ppc.out <- ppcAbund(out, 'chi-square', 2)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, max(n.rep)))
  
  ppc.out <- ppcAbund(out, 'chi-square', 1)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J))
  
  ppc.out <- ppcAbund(out, 'freeman-tukey', 1)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J))
  
  ppc.out <- ppcAbund(out, 'freeman-tukey', 2)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, max(n.rep)))

  ppc.out <- ppcAbund(out, 'chi-square', 0)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J, max(n.rep)))
  
  ppc.out <- ppcAbund(out, 'freeman-tukey', 0)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J, max(n.rep)))
})

# Multiple random intercepts, covariates, and random slopes ---------------
set.seed(1000)
J.x <- 15
J.y <- 15
J <- J.x * J.y
n.rep <- sample(2:4, J, replace = TRUE)
beta <- c(0, -0.2, 0.3)
p.abund <- length(beta)
alpha <- c(0.5, -0.2)
p.det <- length(alpha)
mu.RE <- list(levels = c(10, 15), 
              sigma.sq.mu = c(0.3, 0.9), 
              beta.indx = c(1, 2))
p.RE <- list(levels = c(25, 20), 
             sigma.sq.p = c(0.5, 0.75), 
             alpha.indx = c(1, 2))
family <- 'NB'
kappa <- 0.5
dat <- simNMix(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, alpha = alpha,
	       mu.RE = mu.RE, p.RE = p.RE, sp = FALSE, family = family, kappa = kappa)

y <- dat$y
X <- dat$X
X.re <- dat$X.re
X.p <- dat$X.p
X.p.re <- dat$X.p.re

abund.covs <- cbind(X, X.re)
colnames(abund.covs) <- c('int', 'abund.cov.1', 'abund.cov.2', 'abund.factor.1', 
                          'abund.factor.2')
det.covs <- list(int = X.p[, , 1], 
		 det.cov.1 = X.p[, , 2],
                 det.factor.1 = X.p.re[, , 1], 
                 det.factor.2 = X.p.re[, , 2])
data.list <- list(y = y,
		  abund.covs = abund.covs,
 		  det.covs = det.covs)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 10),
		   alpha.normal = list(mean = 0,
			               var = 2.72))
# Starting values
inits.list <- list(alpha = 0,
		   beta = 0,
		   kappa = 1.3,
		   sigma.sq.mu = 0.5,
		   N = apply(y, 1, max, na.rm = TRUE))
tuning <- list(beta = 0.1, alpha = 0.1, beta.star = 0.1, alpha.star = 0.1,
               kappa = 0.2)
n.batch <- 50
batch.length <- 25
n.burn <- 750
n.thin <- 2
n.chains <- 2
abund.formula <- ~ abund.cov.1 + abund.cov.2 + (1 | abund.factor.1) + (abund.cov.1 | abund.factor.2)
det.formula <- ~ det.cov.1 + (1 | det.factor.1) + (det.cov.1 | det.factor.2)

data = data.list
n.batch = n.batch
batch.length = batch.length
inits = inits.list
priors = prior.list
accept.rate = 0.43
n.omp.threads = 1
verbose = TRUE
family <- 'NB'
n.report = 50
k.fold = 3
k.fold.threads = 3
k.fold.seed <- 100
k.fold.only <- FALSE


out <- NMix(abund.formula = abund.formula,
	    det.formula = det.formula,
	    data = data.list,
	    n.batch = n.batch,
	    batch.length = batch.length,
	    inits = inits.list,
	    family = 'NB',
	    tuning = tuning,
	    priors = prior.list,
	    accept.rate = 0.43,
	    n.omp.threads = 1,
	    verbose = FALSE,
	    n.report = 50,
	    n.burn = n.burn,
	    n.thin = n.thin,
	    n.chains = n.chains,
            k.fold = 2,
            k.fold.threads = 1)

# Test to make sure it worked ---------
test_that("out is of class NMix", {
  expect_s3_class(out, "NMix")
})

# Check cross-validation --------------
test_that("cross-validation works", {
  expect_equal(length(out$rmspe), 1)
  expect_type(out$rmspe, "double")
  expect_gt(out$rmspe, 0)
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$pRE, TRUE)
  expect_equal(out$muRE, TRUE)
})

test_that("random effect gives error when non-numeric", {
  data.list$det.covs$det.factor.1 <- factor(data.list$det.covs$det.factor.1)
  expect_error(out <- NMix(abund.formula = abund.formula,
	    det.formula = det.formula,
	    data = data.list,
	    n.batch = n.batch,
	    batch.length = batch.length,
	    inits = inits.list,
	    family = 'NB',
	    tuning = tuning,
	    priors = prior.list,
	    accept.rate = 0.43,
	    n.omp.threads = 1,
	    verbose = FALSE,
	    n.report = 50,
	    n.burn = n.burn,
	    n.thin = n.thin,
	    n.chains = n.chains))
  data.list$det.covs$det.factor.1 <- as.character(factor(data.list$det.covs$det.factor.1))
  expect_error(out <- NMix(abund.formula = abund.formula,
	    det.formula = det.formula,
	    data = data.list,
	    n.batch = n.batch,
	    batch.length = batch.length,
	    inits = inits.list,
	    family = 'NB',
	    tuning = tuning,
	    priors = prior.list,
	    accept.rate = 0.43,
	    n.omp.threads = 1,
	    verbose = FALSE,
	    n.report = 50,
	    n.burn = n.burn,
	    n.thin = n.thin,
	    n.chains = n.chains))
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, dat$y)
})

# Check Poisson -----------------------
test_that("Poisson works", {
  out <- NMix(abund.formula = abund.formula,
  	    det.formula = det.formula,
  	    data = data.list,
  	    n.batch = n.batch,
  	    batch.length = batch.length,
  	    inits = inits.list,
  	    family = 'Poisson',
  	    tuning = tuning,
  	    priors = prior.list,
  	    accept.rate = 0.43,
  	    n.omp.threads = 1,
  	    verbose = FALSE,
  	    n.report = 50,
  	    n.burn = n.burn,
  	    n.thin = n.thin,
  	    n.chains = 1)
  expect_s3_class(out, "NMix")
})

# Check default priors ----------------
test_that("default priors, multiple chains, inits, burn, thin work", {
  out <- NMix(abund.formula = abund.formula,
	      det.formula = det.formula,
	      data = data.list,
	      n.batch = n.batch,
	      batch.length = batch.length,
	      family = 'NB',
	      accept.rate = 0.43,
	      n.omp.threads = 1,
	      verbose = FALSE,
	      n.report = 50,
	      n.burn = n.burn,
	      n.thin = n.thin,
	      n.chains = n.chains)
  expect_s3_class(out, 'NMix')
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {
  expect_output(NMix(abund.formula = abund.formula,
  	    det.formula = det.formula,
  	    data = data.list,
  	    n.batch = n.batch,
  	    batch.length = batch.length,
  	    inits = inits.list,
  	    family = 'Poisson',
  	    tuning = tuning,
  	    priors = prior.list,
  	    accept.rate = 0.43,
  	    n.omp.threads = 1,
  	    verbose = TRUE,
  	    n.report = 50,
  	    n.burn = n.burn,
  	    n.thin = n.thin,
  	    n.chains = 1))
})

# Check waicAbund -----------------------
test_that("waicAbund works for NMix", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicAbund(out, max(out$N.samples) + 10))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for NMix", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for NMix", {
  X.0 <- cbind(X, X.re)
  colnames(X.0) <- c('(Intercept)', 'abund.cov.1', 'abund.cov.2', 'abund.factor.1', 'abund.factor.2')
  pred.out <- predict(out, X.0)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$mu.0.samples), c(out$n.post * out$n.chains, J))
  expect_equal(dim(pred.out$N.0.samples), c(out$n.post * out$n.chains, J))
})
test_that("detection prediction works", {
  X.p.0 <- cbind(dat$X.p[, 1, ], dat$X.p.re[, 1, ])
  colnames(X.p.0) <- c('(Intercept)', 'det.cov.1', 'det.factor.1', 'det.factor.2')
  pred.out <- predict(out, X.p.0, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, J))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for NMix", {
  n.post.samples <- out$n.post * out$n.chains
  ppc.out <- ppcAbund(out, 'chi-square', 2)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, max(n.rep)))
  
  ppc.out <- ppcAbund(out, 'chi-square', 1)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J))
  
  ppc.out <- ppcAbund(out, 'freeman-tukey', 1)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J))
  
  ppc.out <- ppcAbund(out, 'freeman-tukey', 2)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, max(n.rep)))

  ppc.out <- ppcAbund(out, 'chi-square', 0)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J, max(n.rep)))
  
  ppc.out <- ppcAbund(out, 'freeman-tukey', 0)
  expect_type(ppc.out, "list")
  expect_equal(length(ppc.out$fit.y), n.post.samples)
  expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J, max(n.rep)))
})
