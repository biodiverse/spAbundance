# Test lfMsNMix.R -------------------------------------------------------------

skip_on_cran()

# Intercept only ----------------------------------------------------------
J.x <- 10
J.y <- 10
J <- J.x * J.y
n.rep <- sample(3, size = J, replace = TRUE)
n.sp <- 6
# Community-level covariate effects
# Occurrence
beta.mean <- c(0)
p.abund <- length(beta.mean)
tau.sq.beta <- c(0.3)
# Detection
alpha.mean <- c(0)
tau.sq.alpha <- c(0.2)
p.det <- length(alpha.mean)
# Random effects
mu.RE <- list()
p.RE <- list()
# Draw species-level effects from community means.
beta <- matrix(NA, nrow = n.sp, ncol = p.abund)
alpha <- matrix(NA, nrow = n.sp, ncol = p.det)
for (i in 1:p.abund) {
  beta[, i] <- rnorm(n.sp, beta.mean[i], sqrt(tau.sq.beta[i]))
}
for (i in 1:p.det) {
  alpha[, i] <- rnorm(n.sp, alpha.mean[i], sqrt(tau.sq.alpha[i]))
}
kappa <- runif(n.sp, 0.3, 3)
n.factors <- 2
factor.model <- TRUE
family <- 'NB'

dat <- simMsNMix(J.x = J.x, J.y = J.y, n.rep = n.rep, n.sp = n.sp, beta = beta, alpha = alpha,
	        mu.RE = mu.RE, p.RE = p.RE, sp = FALSE, kappa = kappa, family = family, 
                n.factors = n.factors, factor.model = factor.model)

# Just to keep WAIC relatively fast
y <- ifelse(dat$y > 50, sample(1:50, 100, replace = TRUE), dat$y)
X <- dat$X
X.p <- dat$X.p

abund.covs <- X
colnames(abund.covs) <- c('int')
det.covs <- list(int = X.p[, , 1])
data.list <- list(y = y,
		  abund.covs = abund.covs,
		  coords = dat$coords,
 		  det.covs = det.covs)
# Priors
prior.list <- list(beta.comm.normal = list(mean = 0, var = 10),
		   alpha.comm.normal = list(mean = 0,
			                    var = 2.72), 
                   kappa.unif = list(a = 0, b = 100))
# Starting values
lambda.inits <- matrix(rnorm(n.sp * n.factors), n.sp, n.factors)
diag(lambda.inits) <- 1
lambda.inits[upper.tri(lambda.inits)] <- 0
inits.list <- list(alpha = 0, alpha.comm = 0,
		   beta = 0, beta.comm = 0,
		   sigma.sq.mu = 0.5, sigma.sq.p = 0.5,
		   kappa = 1.3, lambda = lambda.inits, 
		   N = apply(y, c(1, 2), max, na.rm = TRUE))
tuning <- list(beta = 0.1, alpha = 0.1, beta.star = 0.1, alpha.star = 0.1,
               kappa = 0.2, w = 0.5, lambda = 0.5)
n.batch <- 50
batch.length <- 25
n.burn <- 750
n.thin <- 2
n.chains <- 2
abund.formula <- ~ 1
det.formula <- ~ 1


out <- lfMsNMix(abund.formula = abund.formula,
	    det.formula = det.formula,
	    data = data.list,
	    n.batch = n.batch,
	    batch.length = batch.length,
	    inits = inits.list,
	    family = 'NB',
	    n.factors = n.factors,
	    tuning = tuning,
	    priors = prior.list,
	    accept.rate = 0.43,
	    n.omp.threads = 1,
	    verbose = FALSE,
	    n.report = 50,
	    n.burn = n.burn,
	    n.thin = n.thin,
	    n.chains = n.chains)

# Test to make sure it worked ---------
test_that("out is of class lfMsNMix", {
  expect_s3_class(out, "lfMsNMix")
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$pRE, FALSE)
  expect_equal(out$muRE, FALSE)
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, y)
})

# Check Poisson -----------------------
test_that("Poisson works", {
  out <- lfMsNMix(abund.formula = abund.formula,
  	    det.formula = det.formula,
  	    data = data.list,
  	    n.batch = n.batch,
  	    batch.length = batch.length,
  	    inits = inits.list,
  	    family = 'Poisson',
  	    tuning = tuning,
	    n.factors = n.factors,
  	    priors = prior.list,
  	    accept.rate = 0.43,
  	    n.omp.threads = 1,
  	    verbose = FALSE,
  	    n.report = 50,
  	    n.burn = n.burn,
  	    n.thin = n.thin,
  	    n.chains = 1)
  expect_s3_class(out, "lfMsNMix")
})

# Check default priors ----------------
test_that("default priors, multiple chains, inits, burn, thin work", {
  out <- lfMsNMix(abund.formula = abund.formula,
	      det.formula = det.formula,
	      data = data.list,
	      n.batch = n.batch,
	      batch.length = batch.length,
	      family = 'NB',
	      n.factors = n.factors,
	      accept.rate = 0.43,
	      n.omp.threads = 1,
	      verbose = FALSE,
	      n.report = 50,
	      n.burn = n.burn,
	      n.thin = n.thin,
	      n.chains = n.chains)
  expect_s3_class(out, 'lfMsNMix')
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {
  expect_output(lfMsNMix(abund.formula = abund.formula,
  	    det.formula = det.formula,
  	    data = data.list,
  	    n.batch = n.batch,
  	    batch.length = batch.length,
  	    inits = inits.list,
  	    family = 'Poisson',
  	    tuning = tuning,
	    n.factors = n.factors,
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
test_that("waicAbund works for lfMsNMix", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicAbund(out, max(out$N.samples) + 10))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for lfMsNMix", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for lfMsNMix", {
  X.0 <- dat$X
  pred.out <- predict(out, X.0, dat$coords)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$mu.0.samples), c(out$n.post * out$n.chains, n.sp, J))
  expect_equal(dim(pred.out$N.0.samples), c(out$n.post * out$n.chains, n.sp, J))
})
test_that("detection prediction works", {
  J.str <- 100
  X.p.0 <- matrix(1, nrow = J.str, ncol = p.det)
  pred.out <- predict(out, X.p.0, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, n.sp, J.str))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for lfMsNMix", {
  n.post.samples <- out$n.post * out$n.chains
  ppc.out <- ppcAbund(out, 'chi-square', 2)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, n.sp, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, n.sp, max(n.rep)))
  
  ppc.out <- ppcAbund(out, 'chi-square', 1)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, n.sp, J))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, n.sp, J))
  
  ppc.out <- ppcAbund(out, 'freeman-tukey', 1)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, n.sp, J))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, n.sp, J))
  
  ppc.out <- ppcAbund(out, 'freeman-tukey', 2)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, n.sp, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, n.sp, max(n.rep)))

  ppc.out <- ppcAbund(out, 'chi-square', 0)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, n.sp, J, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, n.sp, J, max(n.rep)))
  
  ppc.out <- ppcAbund(out, 'freeman-tukey', 0)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, n.sp, J, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, n.sp, J, max(n.rep)))
})

# Abundance covariate only ------------------------------------------------
J.x <- 10
J.y <- 10
J <- J.x * J.y
n.rep <- sample(3, size = J, replace = TRUE)
n.sp <- 6
# Community-level covariate effects
# Occurrence
beta.mean <- c(0, 0.3)
p.abund <- length(beta.mean)
tau.sq.beta <- c(0.3, 0.5)
# Detection
alpha.mean <- c(0)
tau.sq.alpha <- c(0.2)
p.det <- length(alpha.mean)
# Random effects
mu.RE <- list()
p.RE <- list()
# Draw species-level effects from community means.
beta <- matrix(NA, nrow = n.sp, ncol = p.abund)
alpha <- matrix(NA, nrow = n.sp, ncol = p.det)
for (i in 1:p.abund) {
  beta[, i] <- rnorm(n.sp, beta.mean[i], sqrt(tau.sq.beta[i]))
}
for (i in 1:p.det) {
  alpha[, i] <- rnorm(n.sp, alpha.mean[i], sqrt(tau.sq.alpha[i]))
}
kappa <- runif(n.sp, 0.3, 3)
n.factors <- 2
factor.model <- TRUE
family <- 'NB'

dat <- simMsNMix(J.x = J.x, J.y = J.y, n.rep = n.rep, n.sp = n.sp, beta = beta, alpha = alpha,
	        mu.RE = mu.RE, p.RE = p.RE, sp = FALSE, kappa = kappa, family = family, 
                n.factors = n.factors, factor.model = factor.model)

# Just to keep WAIC relatively fast
y <- ifelse(dat$y > 50, sample(1:50, 100, replace = TRUE), dat$y)

X <- dat$X
X.re <- dat$X.re
X.p <- dat$X.p
X.p.re <- dat$X.p.re

abund.covs <- X
colnames(abund.covs) <- c('int', 'abund.cov.1')
det.covs <- list(int = X.p[, , 1])
data.list <- list(y = y,
		  abund.covs = abund.covs,
		  coords = dat$coords,
 		  det.covs = det.covs)
# Priors
prior.list <- list(beta.comm.normal = list(mean = 0, var = 10),
		   alpha.comm.normal = list(mean = 0,
			                    var = 2.72), 
                   kappa.unif = list(a = 0, b = 100))
# Starting values
inits.list <- list(alpha = 0, alpha.comm = 0,
		   beta = 0, beta.comm = 0,
		   sigma.sq.mu = 0.5, sigma.sq.p = 0.5,
		   kappa = 1.3,
		   N = apply(y, c(1, 2), max, na.rm = TRUE))
tuning <- list(beta = 0.1, alpha = 0.1, beta.star = 0.1, alpha.star = 0.1,
               kappa = 0.2, w = 0.5, lambda = 0.5)
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


out <- lfMsNMix(abund.formula = abund.formula,
	    det.formula = det.formula,
	    data = data.list,
	    n.batch = n.batch,
	    batch.length = batch.length,
	    inits = inits.list,
	    family = 'NB',
	    tuning = tuning,
	    n.factors = n.factors,
	    priors = prior.list,
	    accept.rate = 0.43,
	    n.omp.threads = 1,
	    verbose = FALSE,
	    n.report = 50,
	    n.burn = n.burn,
	    n.thin = n.thin,
	    n.chains = n.chains)

# Test to make sure it worked ---------
test_that("out is of class lfMsNMix", {
  expect_s3_class(out, "lfMsNMix")
})


# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$pRE, FALSE)
  expect_equal(out$muRE, FALSE)
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, y)
})

# Check Poisson -----------------------
test_that("Poisson works", {
  out <- lfMsNMix(abund.formula = abund.formula,
  	    det.formula = det.formula,
  	    data = data.list,
  	    n.batch = n.batch,
  	    batch.length = batch.length,
  	    inits = inits.list,
  	    family = 'Poisson',
  	    tuning = tuning,
	    n.factors = n.factors,
  	    priors = prior.list,
  	    accept.rate = 0.43,
  	    n.omp.threads = 1,
  	    verbose = FALSE,
  	    n.report = 50,
  	    n.burn = n.burn,
  	    n.thin = n.thin,
  	    n.chains = 1)
  expect_s3_class(out, "lfMsNMix")
})

# Check default priors ----------------
test_that("default priors, multiple chains, inits, burn, thin work", {
  out <- lfMsNMix(abund.formula = abund.formula,
	      det.formula = det.formula,
	      data = data.list,
	      n.batch = n.batch,
	      batch.length = batch.length,
	      n.factors = n.factors,
	      family = 'NB',
	      accept.rate = 0.43,
	      n.omp.threads = 1,
	      verbose = FALSE,
	      n.report = 50,
	      n.burn = n.burn,
	      n.thin = n.thin,
	      n.chains = n.chains)
  expect_s3_class(out, 'lfMsNMix')
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {
  expect_output(lfMsNMix(abund.formula = abund.formula,
  	    det.formula = det.formula,
  	    data = data.list,
  	    n.batch = n.batch,
  	    batch.length = batch.length,
  	    inits = inits.list,
  	    family = 'Poisson',
  	    tuning = tuning,
	    n.factors = n.factors,
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
test_that("waicAbund works for lfMsNMix", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicAbund(out, max(out$N.samples) + 10))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for lfMsNMix", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for lfMsNMix", {
  X.0 <- dat$X
  pred.out <- predict(out, X.0, dat$coords)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$mu.0.samples), c(out$n.post * out$n.chains, n.sp, J))
  expect_equal(dim(pred.out$N.0.samples), c(out$n.post * out$n.chains, n.sp, J))
})
test_that("detection prediction works", {
  J.str <- 100
  X.p.0 <- matrix(1, nrow = J.str, ncol = p.det)
  pred.out <- predict(out, X.p.0, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, n.sp, J.str))
})

# Detection covariate only ----------------------------------------------------------
J.x <- 10
J.y <- 10
J <- J.x * J.y
n.rep <- sample(3, size = J, replace = TRUE)
n.sp <- 6
# Community-level covariate effects
# Occurrence
beta.mean <- c(0)
p.abund <- length(beta.mean)
tau.sq.beta <- c(0.3)
# Detection
alpha.mean <- c(0, -0.2)
tau.sq.alpha <- c(0.2, 0.4)
p.det <- length(alpha.mean)
# Random effects
mu.RE <- list()
p.RE <- list()
# Draw species-level effects from community means.
beta <- matrix(NA, nrow = n.sp, ncol = p.abund)
alpha <- matrix(NA, nrow = n.sp, ncol = p.det)
for (i in 1:p.abund) {
  beta[, i] <- rnorm(n.sp, beta.mean[i], sqrt(tau.sq.beta[i]))
}
for (i in 1:p.det) {
  alpha[, i] <- rnorm(n.sp, alpha.mean[i], sqrt(tau.sq.alpha[i]))
}
kappa <- runif(n.sp, 0.3, 3)
n.factors <- 2
factor.model <- TRUE
family <- 'NB'

dat <- simMsNMix(J.x = J.x, J.y = J.y, n.rep = n.rep, n.sp = n.sp, beta = beta, alpha = alpha,
	        mu.RE = mu.RE, p.RE = p.RE, sp = FALSE, kappa = kappa, family = family, 
                n.factors = n.factors, factor.model = factor.model)

# Just to keep WAIC relatively fast
y <- ifelse(dat$y > 50, sample(1:50, 100, replace = TRUE), dat$y)
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
		  coords = dat$coords,
 		  det.covs = det.covs)
# Priors
prior.list <- list(beta.comm.normal = list(mean = 0, var = 10),
		   alpha.comm.normal = list(mean = 0,
			                    var = 2.72), 
                   kappa.unif = list(a = 0, b = 100))
# Starting values
inits.list <- list(alpha = 0, alpha.comm = 0,
		   beta = 0, beta.comm = 0,
		   sigma.sq.mu = 0.5, sigma.sq.p = 0.5,
		   kappa = 1.3,
		   N = apply(y, c(1, 2), max, na.rm = TRUE))
tuning <- list(beta = 0.1, alpha = 0.1, beta.star = 0.1, alpha.star = 0.1,
               kappa = 0.2, w = 0.5, lambda = 0.5)
n.batch <- 50
batch.length <- 25
n.burn <- 750
n.thin <- 2
n.chains <- 2
abund.formula <- ~ 1
det.formula <- ~ det.cov.1

out <- lfMsNMix(abund.formula = abund.formula,
	    det.formula = det.formula,
	    data = data.list,
	    n.batch = n.batch,
	    batch.length = batch.length,
	    inits = inits.list,
	    family = 'NB',
	    tuning = tuning,
	    n.factors = n.factors,
	    priors = prior.list,
	    accept.rate = 0.43,
	    n.omp.threads = 1,
	    verbose = FALSE,
	    n.report = 50,
	    n.burn = n.burn,
	    n.thin = n.thin,
	    n.chains = n.chains)

# Test to make sure it worked ---------
test_that("out is of class lfMsNMix", {
  expect_s3_class(out, "lfMsNMix")
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$pRE, FALSE)
  expect_equal(out$muRE, FALSE)
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, y)
})

# Check Poisson -----------------------
test_that("Poisson works", {
  out <- lfMsNMix(abund.formula = abund.formula,
  	    det.formula = det.formula,
  	    data = data.list,
  	    n.batch = n.batch,
  	    batch.length = batch.length,
  	    inits = inits.list,
  	    family = 'Poisson',
  	    tuning = tuning,
	    n.factors = n.factors,
  	    priors = prior.list,
  	    accept.rate = 0.43,
  	    n.omp.threads = 1,
  	    verbose = FALSE,
  	    n.report = 50,
  	    n.burn = n.burn,
  	    n.thin = n.thin,
  	    n.chains = 1)
  expect_s3_class(out, "lfMsNMix")
})

# Check default priors ----------------
test_that("default priors, multiple chains, inits, burn, thin work", {
  out <- lfMsNMix(abund.formula = abund.formula,
	      det.formula = det.formula,
	      data = data.list,
	      n.batch = n.batch,
	      n.factors = n.factors,
	      batch.length = batch.length,
	      family = 'NB',
	      accept.rate = 0.43,
	      n.omp.threads = 1,
	      verbose = FALSE,
	      n.report = 50,
	      n.burn = n.burn,
	      n.thin = n.thin,
	      n.chains = n.chains)
  expect_s3_class(out, 'lfMsNMix')
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {
  expect_output(lfMsNMix(abund.formula = abund.formula,
  	    det.formula = det.formula,
  	    data = data.list,
  	    n.batch = n.batch,
  	    batch.length = batch.length,
  	    inits = inits.list,
  	    family = 'Poisson',
  	    tuning = tuning,
	    n.factors = n.factors,
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
test_that("waicAbund works for lfMsNMix", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicAbund(out, max(out$N.samples) + 10))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for lfMsNMix", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for lfMsNMix", {
  X.0 <- dat$X
  pred.out <- predict(out, X.0, dat$coords)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$mu.0.samples), c(out$n.post * out$n.chains, n.sp, J))
  expect_equal(dim(pred.out$N.0.samples), c(out$n.post * out$n.chains, n.sp, J))
})
test_that("detection prediction works", {
  J.str <- 100
  X.p.0 <- matrix(1, nrow = J.str, ncol = p.det)
  pred.out <- predict(out, X.p.0, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, n.sp, J.str))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for lfMsNMix", {
  n.post.samples <- out$n.post * out$n.chains
  ppc.out <- ppcAbund(out, 'chi-square', 2)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, n.sp, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, n.sp, max(n.rep)))
  
  ppc.out <- ppcAbund(out, 'chi-square', 1)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, n.sp, J))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, n.sp, J))
  
  ppc.out <- ppcAbund(out, 'freeman-tukey', 1)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, n.sp, J))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, n.sp, J))
  
  ppc.out <- ppcAbund(out, 'freeman-tukey', 2)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, n.sp, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, n.sp, max(n.rep)))

  ppc.out <- ppcAbund(out, 'chi-square', 0)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, n.sp, J, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, n.sp, J, max(n.rep)))
  
  ppc.out <- ppcAbund(out, 'freeman-tukey', 0)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, n.sp, J, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, n.sp, J, max(n.rep)))
})

# Covariates on both ------------------------------------------------------ 
J.x <- 10
J.y <- 10
J <- J.x * J.y
n.rep <- sample(3, size = J, replace = TRUE)
n.sp <- 6
# Community-level covariate effects
# Occurrence
beta.mean <- c(0, 0.2, -0.3)
p.abund <- length(beta.mean)
tau.sq.beta <- c(0.3, 0.3, 0.5)
# Detection
alpha.mean <- c(0, -0.5)
tau.sq.alpha <- c(0.2, 0.6)
p.det <- length(alpha.mean)
# Random effects
mu.RE <- list()
p.RE <- list()
# Draw species-level effects from community means.
beta <- matrix(NA, nrow = n.sp, ncol = p.abund)
alpha <- matrix(NA, nrow = n.sp, ncol = p.det)
for (i in 1:p.abund) {
  beta[, i] <- rnorm(n.sp, beta.mean[i], sqrt(tau.sq.beta[i]))
}
for (i in 1:p.det) {
  alpha[, i] <- rnorm(n.sp, alpha.mean[i], sqrt(tau.sq.alpha[i]))
}
kappa <- runif(n.sp, 0.3, 3)
n.factors <- 2
factor.model <- TRUE
family <- 'NB'

dat <- simMsNMix(J.x = J.x, J.y = J.y, n.rep = n.rep, n.sp = n.sp, beta = beta, alpha = alpha,
	        mu.RE = mu.RE, p.RE = p.RE, sp = FALSE, kappa = kappa, family = family, 
                n.factors = n.factors, factor.model = factor.model)

# Just to keep WAIC relatively fast
y <- ifelse(dat$y > 50, sample(1:50, 100, replace = TRUE), dat$y)
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
		  coords = dat$coords,
 		  det.covs = det.covs)
# Priors
prior.list <- list(beta.comm.normal = list(mean = 0, var = 10),
		   alpha.comm.normal = list(mean = 0,
			                    var = 2.72), 
		   sigma.sq.mu.ig = list(a = 2, b = 0.1),
		   sigma.sq.p.ig = list(a = 2, b = 0.1),
                   kappa.unif = list(a = 0, b = 100))
# Starting values
inits.list <- list(alpha = 0, alpha.comm = 0,
		   beta = 0, beta.comm = 0,
		   sigma.sq.mu = 0.5, sigma.sq.p = 0.5,
		   kappa = 1.3,
		   N = apply(y, c(1, 2), max, na.rm = TRUE))
tuning <- list(beta = 0.1, alpha = 0.1, beta.star = 0.1, alpha.star = 0.1,
               kappa = 0.2, w = 0.5, lambda = 0.5)
n.batch <- 50
batch.length <- 25
n.burn <- 750
n.thin <- 2
n.chains <- 2
abund.formula <- ~ abund.cov.1 + abund.cov.2
det.formula <- ~ det.cov.1


out <- lfMsNMix(abund.formula = abund.formula,
	    det.formula = det.formula,
	    data = data.list,
	    n.batch = n.batch,
	    batch.length = batch.length,
	    inits = inits.list,
	    family = 'NB',
	    tuning = tuning,
	    n.factors = n.factors,
	    priors = prior.list,
	    accept.rate = 0.43,
	    n.omp.threads = 1,
	    verbose = FALSE,
	    n.report = 50,
	    n.burn = n.burn,
	    n.thin = n.thin,
	    n.chains = n.chains)

# Test to make sure it worked ---------
test_that("out is of class lfMsNMix", {
  expect_s3_class(out, "lfMsNMix")
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$pRE, FALSE)
  expect_equal(out$muRE, FALSE)
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, y)
})

# Check Poisson -----------------------
test_that("Poisson works", {
  out <- lfMsNMix(abund.formula = abund.formula,
  	    det.formula = det.formula,
  	    data = data.list,
  	    n.batch = n.batch,
  	    batch.length = batch.length,
  	    inits = inits.list,
  	    family = 'Poisson',
  	    tuning = tuning,
	    n.factors = n.factors,
  	    priors = prior.list,
  	    accept.rate = 0.43,
  	    n.omp.threads = 1,
  	    verbose = FALSE,
  	    n.report = 50,
  	    n.burn = n.burn,
  	    n.thin = n.thin,
  	    n.chains = 1)
  expect_s3_class(out, "lfMsNMix")
})

# Check default priors ----------------
test_that("default priors, multiple chains, inits, burn, thin work", {
  out <- lfMsNMix(abund.formula = abund.formula,
	      det.formula = det.formula,
	      data = data.list,
	      n.batch = n.batch,
	      n.factors = n.factors,
	      batch.length = batch.length,
	      family = 'NB',
	      accept.rate = 0.43,
	      n.omp.threads = 1,
	      verbose = FALSE,
	      n.report = 50,
	      n.burn = n.burn,
	      n.thin = n.thin,
	      n.chains = n.chains)
  expect_s3_class(out, 'lfMsNMix')
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {
  expect_output(lfMsNMix(abund.formula = abund.formula,
  	    det.formula = det.formula,
  	    data = data.list,
  	    n.batch = n.batch,
  	    batch.length = batch.length,
  	    inits = inits.list,
  	    family = 'Poisson',
  	    tuning = tuning,
	    n.factors = n.factors,
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
test_that("waicAbund works for lfMsNMix", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicAbund(out, max(out$N.samples) + 10))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for lfMsNMix", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for lfMsNMix", {
  X.0 <- dat$X
  pred.out <- predict(out, X.0, dat$coords)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$mu.0.samples), c(out$n.post * out$n.chains, n.sp, J))
  expect_equal(dim(pred.out$N.0.samples), c(out$n.post * out$n.chains, n.sp, J))
})
test_that("detection prediction works", {
  X.p.0 <- dat$X.p[, 1, ]
  pred.out <- predict(out, X.p.0, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, n.sp, J))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for lfMsNMix", {
  n.post.samples <- out$n.post * out$n.chains
  ppc.out <- ppcAbund(out, 'chi-square', 2)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, n.sp, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, n.sp, max(n.rep)))
  
  ppc.out <- ppcAbund(out, 'chi-square', 1)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, n.sp, J))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, n.sp, J))
  
  ppc.out <- ppcAbund(out, 'freeman-tukey', 1)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, n.sp, J))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, n.sp, J))
  
  ppc.out <- ppcAbund(out, 'freeman-tukey', 2)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, n.sp, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, n.sp, max(n.rep)))

  ppc.out <- ppcAbund(out, 'chi-square', 0)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, n.sp, J, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, n.sp, J, max(n.rep)))
  
  ppc.out <- ppcAbund(out, 'freeman-tukey', 0)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, n.sp, J, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, n.sp, J, max(n.rep)))
})

# Interactions on both ----------------------------------------------------
J.x <- 10
J.y <- 10
J <- J.x * J.y
n.rep <- sample(3, size = J, replace = TRUE)
n.sp <- 6
# Community-level covariate effects
# Occurrence
beta.mean <- c(0, 0.4, 0.3)
p.abund <- length(beta.mean)
tau.sq.beta <- c(0.3, 0.4, 0.4)
# Detection
alpha.mean <- c(0, 0.2, 0.2)
tau.sq.alpha <- c(0.2, 0.3, 0.6)
p.det <- length(alpha.mean)
# Random effects
mu.RE <- list()
p.RE <- list()
# Draw species-level effects from community means.
beta <- matrix(NA, nrow = n.sp, ncol = p.abund)
alpha <- matrix(NA, nrow = n.sp, ncol = p.det)
for (i in 1:p.abund) {
  beta[, i] <- rnorm(n.sp, beta.mean[i], sqrt(tau.sq.beta[i]))
}
for (i in 1:p.det) {
  alpha[, i] <- rnorm(n.sp, alpha.mean[i], sqrt(tau.sq.alpha[i]))
}
kappa <- runif(n.sp, 0.3, 3)
n.factors <- 2
factor.model <- TRUE
family <- 'NB'

dat <- simMsNMix(J.x = J.x, J.y = J.y, n.rep = n.rep, n.sp = n.sp, beta = beta, alpha = alpha,
	        mu.RE = mu.RE, p.RE = p.RE, sp = FALSE, kappa = kappa, family = family, 
                n.factors = n.factors, factor.model = factor.model)

# Just to keep WAIC relatively fast
y <- ifelse(dat$y > 50, sample(1:50, 100, replace = TRUE), dat$y)
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
		  coords = dat$coords,
 		  det.covs = det.covs)
# Priors
prior.list <- list(beta.comm.normal = list(mean = 0, var = 10),
		   alpha.comm.normal = list(mean = 0,
			                    var = 2.72), 
		   sigma.sq.mu.ig = list(a = 2, b = 0.1),
		   sigma.sq.p.ig = list(a = 2, b = 0.1),
                   kappa.unif = list(a = 0, b = 100))
# Starting values
inits.list <- list(alpha = 0, alpha.comm = 0,
		   beta = 0, beta.comm = 0,
		   sigma.sq.mu = 0.5, sigma.sq.p = 0.5,
		   kappa = 1.3,
		   N = apply(y, c(1, 2), max, na.rm = TRUE))
tuning <- list(beta = 0.1, alpha = 0.1, beta.star = 0.1, alpha.star = 0.1,
               kappa = 0.2, w = 0.5, lambda = 0.5)
n.batch <- 50
batch.length <- 25
n.burn <- 750
n.thin <- 2
n.chains <- 2
abund.formula <- ~ abund.cov.1 * abund.cov.2
det.formula <- ~ det.cov.1 * det.cov.2


out <- lfMsNMix(abund.formula = abund.formula,
	    det.formula = det.formula,
	    data = data.list,
	    n.batch = n.batch,
	    batch.length = batch.length,
	    inits = inits.list,
	    family = 'NB',
	    tuning = tuning,
	    n.factors = n.factors,
	    priors = prior.list,
	    accept.rate = 0.43,
	    n.omp.threads = 1,
	    verbose = FALSE,
	    n.report = 50,
	    n.burn = n.burn,
	    n.thin = n.thin,
	    n.chains = n.chains)

# Test to make sure it worked ---------
test_that("out is of class lfMsNMix", {
  expect_s3_class(out, "lfMsNMix")
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$pRE, FALSE)
  expect_equal(out$muRE, FALSE)
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, y)
})

# Check Poisson -----------------------
test_that("Poisson works", {
  out <- lfMsNMix(abund.formula = abund.formula,
  	    det.formula = det.formula,
  	    data = data.list,
  	    n.batch = n.batch,
  	    batch.length = batch.length,
  	    inits = inits.list,
  	    family = 'Poisson',
  	    tuning = tuning,
	    n.factors = n.factors,
  	    priors = prior.list,
  	    accept.rate = 0.43,
  	    n.omp.threads = 1,
  	    verbose = FALSE,
  	    n.report = 50,
  	    n.burn = n.burn,
  	    n.thin = n.thin,
  	    n.chains = 1)
  expect_s3_class(out, "lfMsNMix")
})

# Check default priors ----------------
test_that("default priors, multiple chains, inits, burn, thin work", {
  out <- lfMsNMix(abund.formula = abund.formula,
	      det.formula = det.formula,
	      data = data.list,
	      n.batch = n.batch,
	      n.factors = n.factors,
	      batch.length = batch.length,
	      family = 'NB',
	      accept.rate = 0.43,
	      n.omp.threads = 1,
	      verbose = FALSE,
	      n.report = 50,
	      n.burn = n.burn,
	      n.thin = n.thin,
	      n.chains = n.chains)
  expect_s3_class(out, 'lfMsNMix')
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {
  expect_output(lfMsNMix(abund.formula = abund.formula,
  	    det.formula = det.formula,
  	    data = data.list,
  	    n.batch = n.batch,
  	    batch.length = batch.length,
  	    inits = inits.list,
  	    family = 'Poisson',
  	    tuning = tuning,
	    n.factors = n.factors,
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
test_that("waicAbund works for lfMsNMix", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicAbund(out, max(out$N.samples) + 10))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for lfMsNMix", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for lfMsNMix", {
  X.0 <- dat$X
  X.0 <- cbind(X.0, X.0[, 2] * X.0[, 3])
  pred.out <- predict(out, X.0, dat$coords)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$mu.0.samples), c(out$n.post * out$n.chains, n.sp, J))
  expect_equal(dim(pred.out$N.0.samples), c(out$n.post * out$n.chains, n.sp, J))
})
test_that("detection prediction works", {
  X.p.0 <- dat$X.p[, 1, ]
  X.p.0 <- cbind(X.p.0, X.p.0[, 2] * X.p.0[, 3])
  pred.out <- predict(out, X.p.0, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, n.sp, J))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for lfMsNMix", {
  n.post.samples <- out$n.post * out$n.chains
  ppc.out <- ppcAbund(out, 'chi-square', 2)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, n.sp, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, n.sp, max(n.rep)))
  
  ppc.out <- ppcAbund(out, 'chi-square', 1)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, n.sp, J))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, n.sp, J))
  
  ppc.out <- ppcAbund(out, 'freeman-tukey', 1)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, n.sp, J))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, n.sp, J))
  
  ppc.out <- ppcAbund(out, 'freeman-tukey', 2)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, n.sp, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, n.sp, max(n.rep)))

  ppc.out <- ppcAbund(out, 'chi-square', 0)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, n.sp, J, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, n.sp, J, max(n.rep)))
  
  ppc.out <- ppcAbund(out, 'freeman-tukey', 0)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, n.sp, J, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, n.sp, J, max(n.rep)))
})

# Site covariate on detection ----------------------------------------------
J.x <- 10
J.y <- 10
J <- J.x * J.y
n.rep <- sample(3, size = J, replace = TRUE)
n.sp <- 6
# Community-level covariate effects
# Occurrence
beta.mean <- c(0, 0.5, -0.5)
p.abund <- length(beta.mean)
tau.sq.beta <- c(0.3, 0.8, 0.3)
# Detection
alpha.mean <- c(0)
tau.sq.alpha <- c(0.2)
p.det <- length(alpha.mean)
# Random effects
mu.RE <- list()
p.RE <- list()
# Draw species-level effects from community means.
beta <- matrix(NA, nrow = n.sp, ncol = p.abund)
alpha <- matrix(NA, nrow = n.sp, ncol = p.det)
for (i in 1:p.abund) {
  beta[, i] <- rnorm(n.sp, beta.mean[i], sqrt(tau.sq.beta[i]))
}
for (i in 1:p.det) {
  alpha[, i] <- rnorm(n.sp, alpha.mean[i], sqrt(tau.sq.alpha[i]))
}
kappa <- runif(n.sp, 0.3, 3)
n.factors <- 2
factor.model <- TRUE
family <- 'NB'

dat <- simMsNMix(J.x = J.x, J.y = J.y, n.rep = n.rep, n.sp = n.sp, beta = beta, alpha = alpha,
	        mu.RE = mu.RE, p.RE = p.RE, sp = FALSE, kappa = kappa, family = family, 
                n.factors = n.factors, factor.model = factor.model)

# Just to keep WAIC relatively fast
y <- ifelse(dat$y > 50, sample(1:50, 100, replace = TRUE), dat$y)
X <- dat$X
X.re <- dat$X.re
X.p <- dat$X.p
X.p.re <- dat$X.p.re

abund.covs <- X
colnames(abund.covs) <- c('int', 'abund.cov.1', 'abund.cov.2')
det.covs <- list(det.cov.1 = X[, 2])
data.list <- list(y = y,
		  abund.covs = abund.covs,
		  coords = dat$coords,
 		  det.covs = det.covs)
# Priors
prior.list <- list(beta.comm.normal = list(mean = 0, var = 10),
		   alpha.comm.normal = list(mean = 0,
			                    var = 2.72), 
		   sigma.sq.mu.ig = list(a = 2, b = 0.1),
		   sigma.sq.p.ig = list(a = 2, b = 0.1),
                   kappa.unif = list(a = 0, b = 100))
# Starting values
inits.list <- list(alpha = 0, alpha.comm = 0,
		   beta = 0, beta.comm = 0,
		   sigma.sq.mu = 0.5, sigma.sq.p = 0.5,
		   kappa = 1.3,
		   N = apply(y, c(1, 2), max, na.rm = TRUE))
tuning <- list(beta = 0.1, alpha = 0.1, beta.star = 0.1, alpha.star = 0.1,
               kappa = 0.2, w = 0.5, lambda = 0.5)
n.batch <- 50
batch.length <- 25
n.burn <- 750
n.thin <- 2
n.chains <- 2
abund.formula <- ~ abund.cov.1 + abund.cov.2
det.formula <- ~ det.cov.1


out <- lfMsNMix(abund.formula = abund.formula,
	    det.formula = det.formula,
	    data = data.list,
	    n.batch = n.batch,
	    batch.length = batch.length,
	    inits = inits.list,
	    family = 'NB',
	    tuning = tuning,
	    n.factors = n.factors,
	    priors = prior.list,
	    accept.rate = 0.43,
	    n.omp.threads = 1,
	    verbose = FALSE,
	    n.report = 50,
	    n.burn = n.burn,
	    n.thin = n.thin,
	    n.chains = n.chains)

# Test to make sure it worked ---------
test_that("out is of class lfMsNMix", {
  expect_s3_class(out, "lfMsNMix")
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$pRE, FALSE)
  expect_equal(out$muRE, FALSE)
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y,y)
})

# Check Poisson -----------------------
test_that("Poisson works", {
  out <- lfMsNMix(abund.formula = abund.formula,
  	    det.formula = det.formula,
  	    data = data.list,
  	    n.batch = n.batch,
  	    batch.length = batch.length,
  	    inits = inits.list,
  	    family = 'Poisson',
  	    tuning = tuning,
	    n.factors = n.factors,
  	    priors = prior.list,
  	    accept.rate = 0.43,
  	    n.omp.threads = 1,
  	    verbose = FALSE,
  	    n.report = 50,
  	    n.burn = n.burn,
  	    n.thin = n.thin,
  	    n.chains = 1)
  expect_s3_class(out, "lfMsNMix")
})

# Check default priors ----------------
test_that("default priors, multiple chains, inits, burn, thin work", {
  out <- lfMsNMix(abund.formula = abund.formula,
	      det.formula = det.formula,
	      data = data.list,
	      n.batch = n.batch,
	      n.factors = n.factors,
	      batch.length = batch.length,
	      family = 'NB',
	      accept.rate = 0.43,
	      n.omp.threads = 1,
	      verbose = FALSE,
	      n.report = 50,
	      n.burn = n.burn,
	      n.thin = n.thin,
	      n.chains = n.chains)
  expect_s3_class(out, 'lfMsNMix')
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {
  expect_output(lfMsNMix(abund.formula = abund.formula,
  	    det.formula = det.formula,
  	    data = data.list,
  	    n.batch = n.batch,
  	    batch.length = batch.length,
  	    inits = inits.list,
  	    family = 'Poisson',
  	    tuning = tuning,
	    n.factors = n.factors,
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
test_that("waicAbund works for lfMsNMix", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicAbund(out, max(out$N.samples) + 10))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for lfMsNMix", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for lfMsNMix", {
  X.0 <- dat$X
  pred.out <- predict(out, X.0, dat$coords)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$mu.0.samples), c(out$n.post * out$n.chains, n.sp, J))
  expect_equal(dim(pred.out$N.0.samples), c(out$n.post * out$n.chains, n.sp, J))
})
test_that("detection prediction works", {
  X.p.0 <- cbind(as.matrix(dat$X.p[, 1, ]), dat$X[, 2])
  pred.out <- predict(out, X.p.0, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, n.sp, J))
})

# Abundance REs + covariates ----------------------------------------------
J.x <- 10
J.y <- 10
J <- J.x * J.y
n.rep <- sample(3, size = J, replace = TRUE)
n.sp <- 6
# Community-level covariate effects
# Occurrence
beta.mean <- c(0, -0.3, 0.3)
p.abund <- length(beta.mean)
tau.sq.beta <- c(0.3, 0.4, 0.4)
# Detection
alpha.mean <- c(0)
tau.sq.alpha <- c(0.2)
p.det <- length(alpha.mean)
# Random effects
mu.RE <- list(levels = c(10, 15, 12),
	       sigma.sq.mu = c(0.3, 0.5, 0.4),
               beta.indx = list(1, 1, 1))
p.RE <- list()
# Draw species-level effects from community means.
beta <- matrix(NA, nrow = n.sp, ncol = p.abund)
alpha <- matrix(NA, nrow = n.sp, ncol = p.det)
for (i in 1:p.abund) {
  beta[, i] <- rnorm(n.sp, beta.mean[i], sqrt(tau.sq.beta[i]))
}
for (i in 1:p.det) {
  alpha[, i] <- rnorm(n.sp, alpha.mean[i], sqrt(tau.sq.alpha[i]))
}
kappa <- runif(n.sp, 0.3, 3)
n.factors <- 2
factor.model <- TRUE
family <- 'NB'

dat <- simMsNMix(J.x = J.x, J.y = J.y, n.rep = n.rep, n.sp = n.sp, beta = beta, alpha = alpha,
	        mu.RE = mu.RE, p.RE = p.RE, sp = FALSE, kappa = kappa, family = family, 
                n.factors = n.factors, factor.model = factor.model)

# Just to keep WAIC relatively fast
y <- ifelse(dat$y > 50, sample(1:50, 100, replace = TRUE), dat$y)
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
		  coords = dat$coords,
 		  det.covs = det.covs)
# Priors
prior.list <- list(beta.comm.normal = list(mean = 0, var = 10),
		   alpha.comm.normal = list(mean = 0,
			                    var = 2.72), 
		   sigma.sq.mu.ig = list(a = 2, b = 0.1),
		   sigma.sq.p.ig = list(a = 2, b = 0.1),
                   kappa.unif = list(a = 0, b = 100))
# Starting values
inits.list <- list(alpha = 0, alpha.comm = 0,
		   beta = 0, beta.comm = 0,
		   sigma.sq.mu = 0.5, sigma.sq.p = 0.5,
		   kappa = 1.3,
		   N = apply(y, c(1, 2), max, na.rm = TRUE))
tuning <- list(beta = 0.1, alpha = 0.1, beta.star = 0.1, alpha.star = 0.1,
               kappa = 0.2, w = 0.5, lambda = 0.5)
n.batch <- 50
batch.length <- 25
n.burn <- 750
n.thin <- 2
n.chains <- 2
abund.formula <- ~ (1 | abund.factor.1) + (1 | abund.factor.2) + (1 | abund.factor.3) + 
                   abund.cov.1 + abund.cov.2
det.formula <- ~ 1 

out <- lfMsNMix(abund.formula = abund.formula,
	    det.formula = det.formula,
	    data = data.list,
	    n.batch = n.batch,
	    batch.length = batch.length,
	    inits = inits.list,
	    family = 'NB',
	    tuning = tuning,
	    n.factors = n.factors,
	    priors = prior.list,
	    accept.rate = 0.43,
	    n.omp.threads = 1,
	    verbose = FALSE,
	    n.report = 50,
	    n.burn = n.burn,
	    n.thin = n.thin,
	    n.chains = n.chains)

# Test to make sure it worked ---------
test_that("out is of class lfMsNMix", {
  expect_s3_class(out, "lfMsNMix")
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$pRE, FALSE)
  expect_equal(out$muRE, TRUE)
})

test_that("random effect gives error when non-numeric", {
  data.list$abund.covs <- as.data.frame(data.list$abund.covs)
  data.list$abund.covs$abund.factor.1 <- factor(data.list$abund.covs$abund.factor.1)
  expect_error(out <- lfMsNMix(abund.formula = abund.formula,
	    det.formula = det.formula,
	    data = data.list,
	    n.batch = n.batch,
	    batch.length = batch.length,
	    inits = inits.list,
	    family = 'NB',
	    tuning = tuning,
	    n.factors = n.factors,
	    priors = prior.list,
	    accept.rate = 0.43,
	    n.omp.threads = 1,
	    verbose = FALSE,
	    n.report = 50,
	    n.burn = n.burn,
	    n.thin = n.thin,
	    n.chains = n.chains))
  data.list$abund.covs$abund.factor.1 <- factor(data.list$abund.covs$abund.factor.1)
  expect_error(out <- lfMsNMix(abund.formula = abund.formula,
	    det.formula = det.formula,
	    data = data.list,
	    n.batch = n.batch,
	    batch.length = batch.length,
	    inits = inits.list,
	    family = 'NB',
	    tuning = tuning,
	    n.factors = n.factors,
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
  expect_equal(out$y, y)
})

# Check Poisson -----------------------
test_that("Poisson works", {
  out <- lfMsNMix(abund.formula = abund.formula,
  	    det.formula = det.formula,
  	    data = data.list,
  	    n.batch = n.batch,
  	    batch.length = batch.length,
  	    inits = inits.list,
  	    family = 'Poisson',
  	    tuning = tuning,
	    n.factors = n.factors,
  	    priors = prior.list,
  	    accept.rate = 0.43,
  	    n.omp.threads = 1,
  	    verbose = FALSE,
  	    n.report = 50,
  	    n.burn = n.burn,
  	    n.thin = n.thin,
  	    n.chains = 1)
  expect_s3_class(out, "lfMsNMix")
})

# Check default priors ----------------
test_that("default priors, multiple chains, inits, burn, thin work", {
  out <- lfMsNMix(abund.formula = abund.formula,
	      det.formula = det.formula,
	      data = data.list,
	      n.factors = n.factors,
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
  expect_s3_class(out, 'lfMsNMix')
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {
  expect_output(lfMsNMix(abund.formula = abund.formula,
  	    det.formula = det.formula,
  	    data = data.list,
  	    n.batch = n.batch,
  	    batch.length = batch.length,
  	    inits = inits.list,
  	    family = 'Poisson',
  	    tuning = tuning,
	    n.factors = n.factors,
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
test_that("waicAbund works for lfMsNMix", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicAbund(out, max(out$N.samples) + 10))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for lfMsNMix", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for lfMsNMix", {
  X.0 <- cbind(X, X.re)
  colnames(X.0) <- c('(Intercept)', 'abund.cov.1', 'abund.cov.2',
		     'abund.factor.1', 'abund.factor.2', 'abund.factor.3')
  pred.out <- predict(out, X.0, dat$coords)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$mu.0.samples), c(out$n.post * out$n.chains, n.sp, J))
  expect_equal(dim(pred.out$N.0.samples), c(out$n.post * out$n.chains, n.sp, J))
})
test_that("detection prediction works", {
  X.p.0 <- as.matrix(dat$X.p[, 1, ])
  pred.out <- predict(out, X.p.0, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, n.sp, J))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for lfMsNMix", {
  n.post.samples <- out$n.post * out$n.chains
  ppc.out <- ppcAbund(out, 'chi-square', 2)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, n.sp, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, n.sp, max(n.rep)))
  
  ppc.out <- ppcAbund(out, 'chi-square', 1)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, n.sp, J))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, n.sp, J))
  
  ppc.out <- ppcAbund(out, 'freeman-tukey', 1)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, n.sp, J))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, n.sp, J))
  
  ppc.out <- ppcAbund(out, 'freeman-tukey', 2)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, n.sp, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, n.sp, max(n.rep)))

  ppc.out <- ppcAbund(out, 'chi-square', 0)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, n.sp, J, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, n.sp, J, max(n.rep)))
  
  ppc.out <- ppcAbund(out, 'freeman-tukey', 0)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, n.sp, J, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, n.sp, J, max(n.rep)))
})

# Multiple random intercepts with covariates ------------------------------
J.x <- 10
J.y <- 10
J <- J.x * J.y
n.rep <- sample(3, size = J, replace = TRUE)
n.sp <- 6
# Community-level covariate effects
# Occurrence
beta.mean <- c(0)
p.abund <- length(beta.mean)
tau.sq.beta <- c(0.3)
# Detection
alpha.mean <- c(0, 0.3)
tau.sq.alpha <- c(0.2, 0.4)
p.det <- length(alpha.mean)
# Random effects
mu.RE <- list()
p.RE <- list(levels = c(10, 12),
	     sigma.sq.p = c(0.5, 0.8),
             alpha.indx = list(1, 2))
# Draw species-level effects from community means.
beta <- matrix(NA, nrow = n.sp, ncol = p.abund)
alpha <- matrix(NA, nrow = n.sp, ncol = p.det)
for (i in 1:p.abund) {
  beta[, i] <- rnorm(n.sp, beta.mean[i], sqrt(tau.sq.beta[i]))
}
for (i in 1:p.det) {
  alpha[, i] <- rnorm(n.sp, alpha.mean[i], sqrt(tau.sq.alpha[i]))
}
kappa <- runif(n.sp, 0.3, 3)
n.factors <- 2
factor.model <- TRUE
family <- 'NB'

dat <- simMsNMix(J.x = J.x, J.y = J.y, n.rep = n.rep, n.sp = n.sp, beta = beta, alpha = alpha,
	        mu.RE = mu.RE, p.RE = p.RE, sp = FALSE, kappa = kappa, family = family, 
                n.factors = n.factors, factor.model = factor.model)

# Just to keep WAIC relatively fast
y <- ifelse(dat$y > 50, sample(1:50, 100, replace = TRUE), dat$y)
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
		  coords = dat$coords,
 		  det.covs = det.covs)
# Priors
prior.list <- list(beta.comm.normal = list(mean = 0, var = 10),
		   alpha.comm.normal = list(mean = 0,
			                    var = 2.72), 
		   sigma.sq.mu.ig = list(a = 2, b = 0.1),
		   sigma.sq.p.ig = list(a = 2, b = 0.1),
                   kappa.unif = list(a = 0, b = 100))
# Starting values
inits.list <- list(alpha = 0, alpha.comm = 0,
		   beta = 0, beta.comm = 0,
		   sigma.sq.mu = 0.5, sigma.sq.p = 0.5,
		   kappa = 1.3,
		   N = apply(y, c(1, 2), max, na.rm = TRUE))
tuning <- list(beta = 0.1, alpha = 0.1, beta.star = 0.1, alpha.star = 0.1,
               kappa = 0.2, w = 0.5, lambda = 0.5)
n.batch <- 50
batch.length <- 25
n.burn <- 750
n.thin <- 2
n.chains <- 2
abund.formula <- ~ 1
det.formula <- ~ det.cov.1 + (1 | det.factor.1) + (1 | det.factor.2)


out <- lfMsNMix(abund.formula = abund.formula,
	    det.formula = det.formula,
	    data = data.list,
	    n.batch = n.batch,
	    batch.length = batch.length,
	    inits = inits.list,
	    family = 'NB',
	    tuning = tuning,
	    n.factors = n.factors,
	    priors = prior.list,
	    accept.rate = 0.43,
	    n.omp.threads = 1,
	    verbose = FALSE,
	    n.report = 50,
	    n.burn = n.burn,
	    n.thin = n.thin,
	    n.chains = n.chains)

# Test to make sure it worked ---------
test_that("out is of class lfMsNMix", {
  expect_s3_class(out, "lfMsNMix")
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$pRE, TRUE)
  expect_equal(out$muRE, FALSE)
})

test_that("random effect gives error when non-numeric", {
  data.list$det.covs$det.factor.1 <- factor(data.list$det.covs$det.factor.1)
  expect_error(out <- lfMsNMix(abund.formula = abund.formula,
	    det.formula = det.formula,
	    data = data.list,
	    n.batch = n.batch,
	    batch.length = batch.length,
	    inits = inits.list,
	    family = 'NB',
	    tuning = tuning,
	    n.factors = n.factors,
	    priors = prior.list,
	    accept.rate = 0.43,
	    n.omp.threads = 1,
	    verbose = FALSE,
	    n.report = 50,
	    n.burn = n.burn,
	    n.thin = n.thin,
	    n.chains = n.chains))
  data.list$det.covs$det.factor.1 <- as.character(factor(data.list$det.covs$det.factor.1))
  expect_error(out <- lfMsNMix(abund.formula = abund.formula,
	    det.formula = det.formula,
	    data = data.list,
	    n.batch = n.batch,
	    batch.length = batch.length,
	    inits = inits.list,
	    family = 'NB',
	    tuning = tuning,
	    n.factors = n.factors,
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
  expect_equal(out$y, y)
})

# Check Poisson -----------------------
test_that("Poisson works", {
  out <- lfMsNMix(abund.formula = abund.formula,
  	    det.formula = det.formula,
  	    data = data.list,
  	    n.batch = n.batch,
  	    batch.length = batch.length,
  	    inits = inits.list,
  	    family = 'Poisson',
  	    tuning = tuning,
	    n.factors = n.factors,
  	    priors = prior.list,
  	    accept.rate = 0.43,
  	    n.omp.threads = 1,
  	    verbose = FALSE,
  	    n.report = 50,
  	    n.burn = n.burn,
  	    n.thin = n.thin,
  	    n.chains = 1)
  expect_s3_class(out, "lfMsNMix")
})

# Check default priors ----------------
test_that("default priors, multiple chains, inits, burn, thin work", {
  out <- lfMsNMix(abund.formula = abund.formula,
	      det.formula = det.formula,
	      data = data.list,
	      n.batch = n.batch,
	      n.factors = n.factors,
	      batch.length = batch.length,
	      family = 'NB',
	      accept.rate = 0.43,
	      n.omp.threads = 1,
	      verbose = FALSE,
	      n.report = 50,
	      n.burn = n.burn,
	      n.thin = n.thin,
	      n.chains = n.chains)
  expect_s3_class(out, 'lfMsNMix')
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {
  expect_output(lfMsNMix(abund.formula = abund.formula,
  	    det.formula = det.formula,
  	    data = data.list,
  	    n.batch = n.batch,
  	    batch.length = batch.length,
  	    inits = inits.list,
  	    family = 'Poisson',
  	    tuning = tuning,
	    n.factors = n.factors,
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
test_that("waicAbund works for lfMsNMix", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicAbund(out, max(out$N.samples) + 10))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for lfMsNMix", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for lfMsNMix", {
  X.0 <- cbind(X)
  pred.out <- predict(out, X.0, dat$coords)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$mu.0.samples), c(out$n.post * out$n.chains, n.sp, J))
  expect_equal(dim(pred.out$N.0.samples), c(out$n.post * out$n.chains, n.sp, J))
})
test_that("detection prediction works", {
  X.p.0 <- cbind(dat$X.p[, 1, ], dat$X.p.re[, 1, ])
  colnames(X.p.0) <- c('(Intercept)', 'det.cov.1', 'det.factor.1', 'det.factor.2')
  pred.out <- predict(out, X.p.0, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, n.sp, J))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for lfMsNMix", {
  n.post.samples <- out$n.post * out$n.chains
  ppc.out <- ppcAbund(out, 'chi-square', 2)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, n.sp, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, n.sp, max(n.rep)))
  
  ppc.out <- ppcAbund(out, 'chi-square', 1)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, n.sp, J))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, n.sp, J))
  
  ppc.out <- ppcAbund(out, 'freeman-tukey', 1)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, n.sp, J))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, n.sp, J))
  
  ppc.out <- ppcAbund(out, 'freeman-tukey', 2)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, n.sp, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, n.sp, max(n.rep)))

  ppc.out <- ppcAbund(out, 'chi-square', 0)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, n.sp, J, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, n.sp, J, max(n.rep)))
  
  ppc.out <- ppcAbund(out, 'freeman-tukey', 0)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, n.sp, J, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, n.sp, J, max(n.rep)))
})

# Multiple random intercepts, covariates, and random slopes ---------------
J.x <- 10
J.y <- 10
J <- J.x * J.y
n.rep <- sample(3, size = J, replace = TRUE)
n.sp <- 6
# Community-level covariate effects
# Occurrence
beta.mean <- c(0, 0.3, -0.2)
p.abund <- length(beta.mean)
tau.sq.beta <- c(0.3, 0.4, 0.5)
# Detection
alpha.mean <- c(0, -0.3)
tau.sq.alpha <- c(0.2, 0.5)
p.det <- length(alpha.mean)
# Random effects
mu.RE <- list(levels = c(10, 15),
	       sigma.sq.mu = c(0.3, 0.5),
               beta.indx = list(1, 2))
p.RE <- list(levels = c(10, 12), 
	     sigma.sq.p = c(0.5, 0.8), 
             alpha.indx = list(1, 2))
# Draw species-level effects from community means.
beta <- matrix(NA, nrow = n.sp, ncol = p.abund)
alpha <- matrix(NA, nrow = n.sp, ncol = p.det)
for (i in 1:p.abund) {
  beta[, i] <- rnorm(n.sp, beta.mean[i], sqrt(tau.sq.beta[i]))
}
for (i in 1:p.det) {
  alpha[, i] <- rnorm(n.sp, alpha.mean[i], sqrt(tau.sq.alpha[i]))
}
kappa <- runif(n.sp, 0.3, 3)
n.factors <- 2
factor.model <- TRUE
family <- 'NB'

dat <- simMsNMix(J.x = J.x, J.y = J.y, n.rep = n.rep, n.sp = n.sp, beta = beta, alpha = alpha,
	        mu.RE = mu.RE, p.RE = p.RE, sp = FALSE, kappa = kappa, family = family, 
                n.factors = n.factors, factor.model = factor.model)

# Just to keep WAIC relatively fast
y <- ifelse(dat$y > 50, sample(1:50, 100, replace = TRUE), dat$y)
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
		  coords = dat$coords,
 		  det.covs = det.covs)
# Priors
prior.list <- list(beta.comm.normal = list(mean = 0, var = 10),
		   alpha.comm.normal = list(mean = 0,
			                    var = 2.72), 
		   sigma.sq.mu.ig = list(a = 2, b = 0.1),
		   sigma.sq.p.ig = list(a = 2, b = 0.1),
                   kappa.unif = list(a = 0, b = 100))
# Starting values
inits.list <- list(alpha = 0, alpha.comm = 0,
		   beta = 0, beta.comm = 0,
		   sigma.sq.mu = 0.5, sigma.sq.p = 0.5,
		   kappa = 1.3,
		   N = apply(y, c(1, 2), max, na.rm = TRUE))
tuning <- list(beta = 0.1, alpha = 0.1, beta.star = 0.1, alpha.star = 0.1,
               kappa = 0.2, w = 0.5, lambda = 0.5)
n.batch <- 50
batch.length <- 25
n.burn <- 750
n.thin <- 2
n.chains <- 2
abund.formula <- ~ abund.cov.1 + abund.cov.2 + (1 | abund.factor.1) + (abund.cov.1 | abund.factor.2)
det.formula <- ~ det.cov.1 + (1 | det.factor.1) + (det.cov.1 | det.factor.2)


out <- lfMsNMix(abund.formula = abund.formula,
	    det.formula = det.formula,
	    data = data.list,
	    n.batch = n.batch,
	    batch.length = batch.length,
	    inits = inits.list,
	    family = 'NB',
	    tuning = tuning,
	    n.factors = n.factors,
	    priors = prior.list,
	    accept.rate = 0.43,
	    n.omp.threads = 1,
	    verbose = FALSE,
	    n.report = 50,
	    n.burn = n.burn,
	    n.thin = n.thin,
	    n.chains = n.chains)

# Test to make sure it worked ---------
test_that("out is of class lfMsNMix", {
  expect_s3_class(out, "lfMsNMix")
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$pRE, TRUE)
  expect_equal(out$muRE, TRUE)
})

test_that("random effect gives error when non-numeric", {
  data.list$det.covs$det.factor.1 <- factor(data.list$det.covs$det.factor.1)
  expect_error(out <- lfMsNMix(abund.formula = abund.formula,
	    det.formula = det.formula,
	    data = data.list,
	    n.batch = n.batch,
	    batch.length = batch.length,
	    inits = inits.list,
	    family = 'NB',
	    tuning = tuning,
	    n.factors = n.factors,
	    priors = prior.list,
	    accept.rate = 0.43,
	    n.omp.threads = 1,
	    verbose = FALSE,
	    n.report = 50,
	    n.burn = n.burn,
	    n.thin = n.thin,
	    n.chains = n.chains))
  data.list$det.covs$det.factor.1 <- as.character(factor(data.list$det.covs$det.factor.1))
  expect_error(out <- lfMsNMix(abund.formula = abund.formula,
	    det.formula = det.formula,
	    data = data.list,
	    n.batch = n.batch,
	    batch.length = batch.length,
	    inits = inits.list,
	    family = 'NB',
	    tuning = tuning,
	    n.factors = n.factors,
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
  expect_equal(out$y, y)
})

# Check Poisson -----------------------
test_that("Poisson works", {
  out <- lfMsNMix(abund.formula = abund.formula,
  	    det.formula = det.formula,
  	    data = data.list,
  	    n.batch = n.batch,
  	    batch.length = batch.length,
  	    inits = inits.list,
  	    family = 'Poisson',
  	    tuning = tuning,
	    n.factors = n.factors,
  	    priors = prior.list,
  	    accept.rate = 0.43,
  	    n.omp.threads = 1,
  	    verbose = FALSE,
  	    n.report = 50,
  	    n.burn = n.burn,
  	    n.thin = n.thin,
  	    n.chains = 1)
  expect_s3_class(out, "lfMsNMix")
})

# Check default priors ----------------
test_that("default priors, multiple chains, inits, burn, thin work", {
  out <- lfMsNMix(abund.formula = abund.formula,
	      det.formula = det.formula,
	      data = data.list,
	      n.batch = n.batch,
	      n.factors = n.factors,
	      batch.length = batch.length,
	      family = 'NB',
	      accept.rate = 0.43,
	      n.omp.threads = 1,
	      verbose = FALSE,
	      n.report = 50,
	      n.burn = n.burn,
	      n.thin = n.thin,
	      n.chains = n.chains)
  expect_s3_class(out, 'lfMsNMix')
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {
  expect_output(lfMsNMix(abund.formula = abund.formula,
  	    det.formula = det.formula,
  	    data = data.list,
  	    n.batch = n.batch,
  	    batch.length = batch.length,
  	    inits = inits.list,
  	    family = 'Poisson',
  	    tuning = tuning,
	    n.factors = n.factors,
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
test_that("waicAbund works for lfMsNMix", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicAbund(out, max(out$N.samples) + 10))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for lfMsNMix", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for lfMsNMix", {
  X.0 <- cbind(X, X.re)
  colnames(X.0) <- c('(Intercept)', 'abund.cov.1', 'abund.cov.2', 'abund.factor.1', 'abund.factor.2')
  pred.out <- predict(out, X.0, dat$coords)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$mu.0.samples), c(out$n.post * out$n.chains, n.sp, J))
  expect_equal(dim(pred.out$N.0.samples), c(out$n.post * out$n.chains, n.sp, J))
})
test_that("detection prediction works", {
  X.p.0 <- cbind(dat$X.p[, 1, ], dat$X.p.re[, 1, ])
  colnames(X.p.0) <- c('(Intercept)', 'det.cov.1', 'det.factor.1', 'det.factor.2')
  pred.out <- predict(out, X.p.0, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$p.0.samples), c(out$n.post * out$n.chains, n.sp, J))
})
# Check PPCs --------------------------
test_that("posterior predictive checks work for lfMsNMix", {
  n.post.samples <- out$n.post * out$n.chains
  ppc.out <- ppcAbund(out, 'chi-square', 2)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, n.sp, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, n.sp, max(n.rep)))
  
  ppc.out <- ppcAbund(out, 'chi-square', 1)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, n.sp, J))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, n.sp, J))
  
  ppc.out <- ppcAbund(out, 'freeman-tukey', 1)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, n.sp, J))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, n.sp, J))
  
  ppc.out <- ppcAbund(out, 'freeman-tukey', 2)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, n.sp, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, n.sp, max(n.rep)))

  ppc.out <- ppcAbund(out, 'chi-square', 0)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, n.sp, J, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, n.sp, J, max(n.rep)))
  
  ppc.out <- ppcAbund(out, 'freeman-tukey', 0)
  expect_type(ppc.out, "list")
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, n.sp, J, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, n.sp, J, max(n.rep)))
})

