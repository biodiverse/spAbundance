# Test msAbund.R ----------------------------------------------------------
skip_on_cran()

# Intercept only ----------------------------------------------------------
set.seed(10101)
J.x <- 10
J.y <- 10
J <- J.x * J.y
n.rep <- sample(3, size = J, replace = TRUE)
n.sp <- 6
# Community-level covariate effects
beta.mean <- c(0.5)
p.abund <- length(beta.mean)
tau.sq.beta <- c(0.3)
# Random effects
mu.RE <- list()
# Draw species-level effects from community means.
beta <- matrix(NA, nrow = n.sp, ncol = p.abund)
for (i in 1:p.abund) {
  beta[, i] <- rnorm(n.sp, beta.mean[i], sqrt(tau.sq.beta[i]))
}
sp <- FALSE
factor.model <- FALSE
kappa <- runif(n.sp, 0.1, 1)
family <- 'NB'

dat <- simMsAbund(J.x = J.x, J.y = J.y, n.rep = n.rep, n.sp = n.sp, beta = beta,
	          mu.RE = mu.RE, sp = sp, kappa = kappa, family = family,
                  factor.model = factor.model, n.factors = n.factors)

y <- dat$y
X <- dat$X
coords <- dat$coords

data.list <- list(y = y)

prior.list <- list(beta.comm.normal = list(mean = 0, var = 100),
                   kappa.unif = list(a = 0, b = 10),
                   tau.sq.beta.ig = list(a = .1, b = .1))

inits.list <- list(beta.comm = 0,
		   beta = 0,
		   tau.sq.beta = 1)
tuning.list <- list(beta = 0.1, beta.star = 0.1, kappa = 0.5)

n.batch <- 40
batch.length <- 25
n.burn <- 400
n.thin <- 2
n.chains <- 1
formula <- ~ 1

out <- msAbund(formula = formula,
	       data = data.list,
	       n.batch = n.batch,
	       batch.length = batch.length,
	       inits = inits.list,
	       family = 'NB',
	       tuning = tuning.list,
	       priors = prior.list,
	       accept.rate = 0.43,
	       n.omp.threads = 1,
	       verbose = FALSE,
	       n.report = 10,
	       n.burn = n.burn,
	       n.thin = n.thin,
	       n.chains = 1)

# Test to make sure it worked ---------
test_that("out is of class msAbund", {
  expect_s3_class(out, "msAbund")
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$muRE, FALSE)
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, dat$y)
})

# Check Poisson -----------------------
test_that("Poisson works", {
  out <- msAbund(formula = formula,
  	     data = data.list, 
  	     n.batch = n.batch, 
  	     batch.length = batch.length, 
  	     inits = inits.list, 
  	     family = 'Poisson',
  	     tuning = tuning.list,
  	     priors = prior.list, 
  	     accept.rate = 0.43, 
  	     n.omp.threads = 1, 
  	     verbose = FALSE, 
  	     n.report = 10,
  	     n.burn = n.burn,
  	     n.thin = n.thin,
  	     n.chains = 1) 
  expect_s3_class(out, "msAbund")
})

# Check default priors ----------------
test_that("default priors, multiple chains, inits, burn, thin work", {
  out <- msAbund(formula = formula,
  	     data = data.list, 
  	     n.batch = n.batch, 
  	     batch.length = batch.length, 
  	     family = 'NB',
  	     accept.rate = 0.43, 
  	     n.omp.threads = 1, 
  	     verbose = FALSE, 
  	     n.report = 10,
  	     n.burn = n.burn,
  	     n.thin = n.thin,
  	     n.chains = 2) 
  expect_s3_class(out, "msAbund")
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {
  expect_output(msAbund(formula = formula,
  	     data = data.list, 
  	     n.batch = n.batch, 
  	     batch.length = batch.length, 
  	     family = 'NB',
  	     accept.rate = 0.43, 
  	     n.omp.threads = 1, 
  	     verbose = TRUE, 
  	     n.report = 10,
  	     n.burn = n.burn,
  	     n.thin = n.thin,
  	     n.chains = 1))
})

# Check waicAbund -----------------------
test_that("waicAbund works for msAbund", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicAbund(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for msAbund", {
  fitted.out <- fitted(out)
  expect_is(fitted.out, "array")
})

# Check predictions -------------------
test_that("predict works for msAbund", {
  X.0 <- dat$X
  pred.out <- predict(out, X.0)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$mu.0.samples), c(out$n.post * out$n.chains, n.sp, J, max(n.rep)))
  expect_equal(dim(pred.out$y.0.samples), c(out$n.post * out$n.chains, n.sp, J, max(n.rep)))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for msAbund", {
  n.post.samples <- out$n.post * out$n.chains
  ppc.out <- ppcAbund(out, 'chi-square', 0)
  expect_type(ppc.out, "list")
  expect_equal(nrow(ppc.out$fit.y), n.post.samples)
  expect_equal(nrow(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, n.sp, J, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, n.sp, J, max(n.rep)))
  
  ppc.out <- ppcAbund(out, 'freeman-tukey', 0)
  expect_type(ppc.out, "list")
  expect_equal(nrow(ppc.out$fit.y), n.post.samples)
  expect_equal(nrow(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, n.sp, J, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, n.sp, J, max(n.rep)))
})

# Covariates --------------------------------------------------------------
J.x <- 10
J.y <- 10
J <- J.x * J.y
n.rep <- sample(3, size = J, replace = TRUE)
n.sp <- 6
# Community-level covariate effects
beta.mean <- c(0, 0.9, -0.3)
p.abund <- length(beta.mean)
tau.sq.beta <- c(0.3, 0.2, 0.5)
# Random effects
mu.RE <- list()
# Draw species-level effects from community means.
beta <- matrix(NA, nrow = n.sp, ncol = p.abund)
for (i in 1:p.abund) {
  beta[, i] <- rnorm(n.sp, beta.mean[i], sqrt(tau.sq.beta[i]))
}
sp <- FALSE
factor.model <- FALSE
kappa <- runif(n.sp, 0.1, 1)
family <- 'NB'

dat <- simMsAbund(J.x = J.x, J.y = J.y, n.rep = n.rep, n.sp = n.sp, beta = beta,
	          mu.RE = mu.RE, sp = sp, kappa = kappa, family = family,
                  factor.model = factor.model, n.factors = n.factors)

y <- dat$y
X <- dat$X
coords <- dat$coords

covs <- list(cov.1 = X[, , 2], 
	     cov.2 = X[, , 3])

data.list <- list(y = y, 
                  covs = covs)

prior.list <- list(beta.comm.normal = list(mean = 0, var = 100),
                   kappa.unif = list(a = 0, b = 10),
                   tau.sq.beta.ig = list(a = .1, b = .1))

inits.list <- list(beta.comm = 0,
		   beta = 0,
		   tau.sq.beta = 1)
tuning.list <- list(beta = 0.1, beta.star = 0.1, kappa = 0.5)

n.batch <- 40
batch.length <- 25
n.burn <- 400
n.thin <- 2
n.chains <- 1
formula <- ~ cov.1 + cov.2

out <- msAbund(formula = formula,
	       data = data.list,
	       n.batch = n.batch,
	       batch.length = batch.length,
	       inits = inits.list,
	       family = 'NB',
	       tuning = tuning.list,
	       priors = prior.list,
	       accept.rate = 0.43,
	       n.omp.threads = 1,
	       verbose = FALSE,
	       n.report = 10,
	       n.burn = n.burn,
	       n.thin = n.thin,
	       n.chains = 1)

# Test to make sure it worked ---------
test_that("out is of class msAbund", {
  expect_s3_class(out, "msAbund")
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$muRE, FALSE)
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, dat$y)
})

# Check Poisson -----------------------
test_that("Poisson works", {
  out <- msAbund(formula = formula,
  	     data = data.list, 
  	     n.batch = n.batch, 
  	     batch.length = batch.length, 
  	     inits = inits.list, 
  	     family = 'Poisson',
  	     tuning = tuning.list,
  	     priors = prior.list, 
  	     accept.rate = 0.43, 
  	     n.omp.threads = 1, 
  	     verbose = FALSE, 
  	     n.report = 10,
  	     n.burn = n.burn,
  	     n.thin = n.thin,
  	     n.chains = 1) 
  expect_s3_class(out, "msAbund")
})

# Check default priors ----------------
test_that("default priors, multiple chains, inits, burn, thin work", {
  out <- msAbund(formula = formula,
  	     data = data.list, 
  	     n.batch = n.batch, 
  	     batch.length = batch.length, 
  	     family = 'NB',
  	     accept.rate = 0.43, 
  	     n.omp.threads = 1, 
  	     verbose = FALSE, 
  	     n.report = 10,
  	     n.burn = n.burn,
  	     n.thin = n.thin,
  	     n.chains = 2) 
  expect_s3_class(out, "msAbund")
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {
  expect_output(msAbund(formula = formula,
  	     data = data.list, 
  	     n.batch = n.batch, 
  	     batch.length = batch.length, 
  	     family = 'NB',
  	     accept.rate = 0.43, 
  	     n.omp.threads = 1, 
  	     verbose = TRUE, 
  	     n.report = 10,
  	     n.burn = n.burn,
  	     n.thin = n.thin,
  	     n.chains = 1))
})

# Check waicAbund -----------------------
test_that("waicAbund works for msAbund", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicAbund(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for msAbund", {
  fitted.out <- fitted(out)
  expect_is(fitted.out, "array")
})

# Check predictions -------------------
test_that("predict works for msAbund", {
  X.0 <- dat$X
  pred.out <- predict(out, X.0)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$mu.0.samples), c(out$n.post * out$n.chains, n.sp, J, max(n.rep)))
  expect_equal(dim(pred.out$y.0.samples), c(out$n.post * out$n.chains, n.sp, J, max(n.rep)))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for msAbund", {
  n.post.samples <- out$n.post * out$n.chains
  ppc.out <- ppcAbund(out, 'chi-square', 0)
  expect_type(ppc.out, "list")
  expect_equal(nrow(ppc.out$fit.y), n.post.samples)
  expect_equal(nrow(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, n.sp, J, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, n.sp, J, max(n.rep)))
  
  ppc.out <- ppcAbund(out, 'freeman-tukey', 0)
  expect_type(ppc.out, "list")
  expect_equal(nrow(ppc.out$fit.y), n.post.samples)
  expect_equal(nrow(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, n.sp, J, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, n.sp, J, max(n.rep)))
})

# Interactions ------------------------------------------------------------
J.x <- 10
J.y <- 10
J <- J.x * J.y
n.rep <- rep(3, J)
n.sp <- 6
# Community-level covariate effects
beta.mean <- c(0, 0.9, -0.3)
p.abund <- length(beta.mean)
tau.sq.beta <- c(0.3, 0.2, 0.5)
# Random effects
mu.RE <- list()
# Draw species-level effects from community means.
beta <- matrix(NA, nrow = n.sp, ncol = p.abund)
for (i in 1:p.abund) {
  beta[, i] <- rnorm(n.sp, beta.mean[i], sqrt(tau.sq.beta[i]))
}
sp <- FALSE
factor.model <- FALSE
kappa <- runif(n.sp, 0.1, 1)
family <- 'NB'

dat <- simMsAbund(J.x = J.x, J.y = J.y, n.rep = n.rep, n.sp = n.sp, beta = beta,
	          mu.RE = mu.RE, sp = sp, kappa = kappa, family = family,
                  factor.model = factor.model, n.factors = n.factors)

y <- dat$y
X <- dat$X
coords <- dat$coords

covs <- list(cov.1 = X[, , 2], 
	     cov.2 = X[, , 3])

data.list <- list(y = y, 
                  covs = covs)

prior.list <- list(beta.comm.normal = list(mean = 0, var = 100),
                   kappa.unif = list(a = 0, b = 10),
                   tau.sq.beta.ig = list(a = .1, b = .1))

inits.list <- list(beta.comm = 0,
		   beta = 0,
		   tau.sq.beta = 1)
tuning.list <- list(beta = 0.1, beta.star = 0.1, kappa = 0.5)

n.batch <- 40
batch.length <- 25
n.burn <- 400
n.thin <- 2
n.chains <- 1
formula <- ~ cov.1 * cov.2

out <- msAbund(formula = formula,
	       data = data.list,
	       n.batch = n.batch,
	       batch.length = batch.length,
	       inits = inits.list,
	       family = 'NB',
	       tuning = tuning.list,
	       priors = prior.list,
	       accept.rate = 0.43,
	       n.omp.threads = 1,
	       verbose = FALSE,
	       n.report = 10,
	       n.burn = n.burn,
	       n.thin = n.thin,
	       n.chains = 1)

# Test to make sure it worked ---------
test_that("out is of class msAbund", {
  expect_s3_class(out, "msAbund")
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$muRE, FALSE)
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, dat$y)
})

# Check Poisson -----------------------
test_that("Poisson works", {
  out <- msAbund(formula = formula,
  	     data = data.list, 
  	     n.batch = n.batch, 
  	     batch.length = batch.length, 
  	     inits = inits.list, 
  	     family = 'Poisson',
  	     tuning = tuning.list,
  	     priors = prior.list, 
  	     accept.rate = 0.43, 
  	     n.omp.threads = 1, 
  	     verbose = FALSE, 
  	     n.report = 10,
  	     n.burn = n.burn,
  	     n.thin = n.thin,
  	     n.chains = 1) 
  expect_s3_class(out, "msAbund")
})

# Check default priors ----------------
test_that("default priors, multiple chains, inits, burn, thin work", {
  out <- msAbund(formula = formula,
  	     data = data.list, 
  	     n.batch = n.batch, 
  	     batch.length = batch.length, 
  	     family = 'NB',
  	     accept.rate = 0.43, 
  	     n.omp.threads = 1, 
  	     verbose = FALSE, 
  	     n.report = 10,
  	     n.burn = n.burn,
  	     n.thin = n.thin,
  	     n.chains = 2) 
  expect_s3_class(out, "msAbund")
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {
  expect_output(msAbund(formula = formula,
  	     data = data.list, 
  	     n.batch = n.batch, 
  	     batch.length = batch.length, 
  	     family = 'NB',
  	     accept.rate = 0.43, 
  	     n.omp.threads = 1, 
  	     verbose = TRUE, 
  	     n.report = 10,
  	     n.burn = n.burn,
  	     n.thin = n.thin,
  	     n.chains = 1))
})

# Check waicAbund -----------------------
test_that("waicAbund works for msAbund", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicAbund(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for msAbund", {
  fitted.out <- fitted(out)
  expect_is(fitted.out, "array")
})

# Check predictions -------------------
test_that("predict works for msAbund", {
  X.0 <- abind(dat$X, dat$X[, , 2] * dat$X[, , 3], along = 3)
  pred.out <- predict(out, X.0)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$mu.0.samples), c(out$n.post * out$n.chains, n.sp, J, max(n.rep)))
  expect_equal(dim(pred.out$y.0.samples), c(out$n.post * out$n.chains, n.sp, J, max(n.rep)))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for msAbund", {
  n.post.samples <- out$n.post * out$n.chains
  ppc.out <- ppcAbund(out, 'chi-square', 0)
  expect_type(ppc.out, "list")
  expect_equal(nrow(ppc.out$fit.y), n.post.samples)
  expect_equal(nrow(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, n.sp, J, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, n.sp, J, max(n.rep)))
  
  ppc.out <- ppcAbund(out, 'freeman-tukey', 0)
  expect_type(ppc.out, "list")
  expect_equal(nrow(ppc.out$fit.y), n.post.samples)
  expect_equal(nrow(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, n.sp, J, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, n.sp, J, max(n.rep)))
})

# Site-level covariate on abundance ---------------------------------------
J.x <- 10
J.y <- 10
J <- J.x * J.y
n.rep <- sample(3, size = J, replace = TRUE)
# n.rep <- rep(5, J)
n.sp <- 6
# Community-level covariate effects
beta.mean <- c(0, 0.9, -0.3)
p.abund <- length(beta.mean)
tau.sq.beta <- c(0.3, 0.2, 0.5)
# Random effects
mu.RE <- list()
# Draw species-level effects from community means.
beta <- matrix(NA, nrow = n.sp, ncol = p.abund)
for (i in 1:p.abund) {
  beta[, i] <- rnorm(n.sp, beta.mean[i], sqrt(tau.sq.beta[i]))
}
sp <- FALSE
factor.model <- FALSE
kappa <- runif(n.sp, 0.1, 1)
family <- 'NB'

dat <- simMsAbund(J.x = J.x, J.y = J.y, n.rep = n.rep, n.sp = n.sp, beta = beta,
	          mu.RE = mu.RE, sp = sp, kappa = kappa, family = family,
                  factor.model = factor.model, n.factors = n.factors)

y <- dat$y
X <- dat$X
coords <- dat$coords

covs <- list(cov.1 = X[, 1, 2], 
	     cov.2 = X[, 1, 3])

data.list <- list(y = y, 
                  covs = covs)

prior.list <- list(beta.comm.normal = list(mean = 0, var = 100),
                   kappa.unif = list(a = 0, b = 10),
                   tau.sq.beta.ig = list(a = .1, b = .1))

inits.list <- list(beta.comm = 0,
		   beta = 0,
		   tau.sq.beta = 1)
tuning.list <- list(beta = 0.1, beta.star = 0.1, kappa = 0.5)

n.batch <- 40
batch.length <- 25
n.burn <- 400
n.thin <- 2
n.chains <- 1
formula <- ~ cov.1 + cov.2


out <- msAbund(formula = formula,
	       data = data.list,
	       n.batch = n.batch,
	       batch.length = batch.length,
	       inits = inits.list,
	       family = 'NB',
	       tuning = tuning.list,
	       priors = prior.list,
	       accept.rate = 0.43,
	       n.omp.threads = 1,
	       verbose = FALSE,
	       n.report = 10,
	       n.burn = n.burn,
	       n.thin = n.thin,
	       n.chains = 1)

# Test to make sure it worked ---------
test_that("out is of class msAbund", {
  expect_s3_class(out, "msAbund")
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$muRE, FALSE)
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, dat$y)
})

# Check Poisson -----------------------
test_that("Poisson works", {
  out <- msAbund(formula = formula,
  	     data = data.list, 
  	     n.batch = n.batch, 
  	     batch.length = batch.length, 
  	     inits = inits.list, 
  	     family = 'Poisson',
  	     tuning = tuning.list,
  	     priors = prior.list, 
  	     accept.rate = 0.43, 
  	     n.omp.threads = 1, 
  	     verbose = FALSE, 
  	     n.report = 10,
  	     n.burn = n.burn,
  	     n.thin = n.thin,
  	     n.chains = 1) 
  expect_s3_class(out, "msAbund")
})

# Check default priors ----------------
test_that("default priors, multiple chains, inits, burn, thin work", {
  out <- msAbund(formula = formula,
  	     data = data.list, 
  	     n.batch = n.batch, 
  	     batch.length = batch.length, 
  	     family = 'NB',
  	     accept.rate = 0.43, 
  	     n.omp.threads = 1, 
  	     verbose = FALSE, 
  	     n.report = 10,
  	     n.burn = n.burn,
  	     n.thin = n.thin,
  	     n.chains = 2) 
  expect_s3_class(out, "msAbund")
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {
  expect_output(msAbund(formula = formula,
  	     data = data.list, 
  	     n.batch = n.batch, 
  	     batch.length = batch.length, 
  	     family = 'NB',
  	     accept.rate = 0.43, 
  	     n.omp.threads = 1, 
  	     verbose = TRUE, 
  	     n.report = 10,
  	     n.burn = n.burn,
  	     n.thin = n.thin,
  	     n.chains = 1))
})

# Check waicAbund -----------------------
test_that("waicAbund works for msAbund", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicAbund(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for msAbund", {
  fitted.out <- fitted(out)
  expect_is(fitted.out, "array")
})

# Check predictions -------------------
test_that("predict works for msAbund", {
  X.0 <- dat$X
  pred.out <- predict(out, X.0)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$mu.0.samples), c(out$n.post * out$n.chains, n.sp, J, max(n.rep)))
  expect_equal(dim(pred.out$y.0.samples), c(out$n.post * out$n.chains, n.sp, J, max(n.rep)))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for msAbund", {
  n.post.samples <- out$n.post * out$n.chains
  ppc.out <- ppcAbund(out, 'chi-square', 0)
  expect_type(ppc.out, "list")
  expect_equal(nrow(ppc.out$fit.y), n.post.samples)
  expect_equal(nrow(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, n.sp, J, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, n.sp, J, max(n.rep)))
  
  ppc.out <- ppcAbund(out, 'freeman-tukey', 0)
  expect_type(ppc.out, "list")
  expect_equal(nrow(ppc.out$fit.y), n.post.samples)
  expect_equal(nrow(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, n.sp, J, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, n.sp, J, max(n.rep)))
})

# Random intercept --------------------------------------------------------
J.x <- 10
J.y <- 10
J <- J.x * J.y
# n.rep <- sample(3, size = J, replace = TRUE)
n.rep <- rep(1, J)
n.sp <- 6
# Community-level covariate effects
beta.mean <- c(0)
p.abund <- length(beta.mean)
tau.sq.beta <- c(0.3)
# Random effects
mu.RE <- list(levels = c(20),
              sigma.sq.mu = c(1))
# Draw species-level effects from community means.
beta <- matrix(NA, nrow = n.sp, ncol = p.abund)
for (i in 1:p.abund) {
  beta[, i] <- rnorm(n.sp, beta.mean[i], sqrt(tau.sq.beta[i]))
}
sp <- FALSE
factor.model <- FALSE
kappa <- runif(n.sp, 0.1, 1)
family <- 'NB'

dat <- simMsAbund(J.x = J.x, J.y = J.y, n.rep = n.rep, n.sp = n.sp, beta = beta,
	          mu.RE = mu.RE, sp = sp, kappa = kappa, family = family,
                  factor.model = factor.model, n.factors = n.factors)

y <- dat$y
X <- dat$X
X.re <- dat$X.re
coords <- dat$coords

covs <- list(factor.1 = X.re[, , 1])

data.list <- list(y = y, 
                  covs = covs)

prior.list <- list(beta.comm.normal = list(mean = 0, var = 100),
                   kappa.unif = list(a = 0, b = 10),
                   tau.sq.beta.ig = list(a = .1, b = .1))

inits.list <- list(beta.comm = 0,
		   beta = 0,
		   tau.sq.beta = 1)
tuning.list <- list(beta = 0.1, beta.star = 0.1, kappa = 0.5)

n.batch <- 40
batch.length <- 25
n.burn <- 400
n.thin <- 2
n.chains <- 1
formula <- ~ (1 | factor.1)


out <- msAbund(formula = formula,
	       data = data.list,
	       n.batch = n.batch,
	       batch.length = batch.length,
	       inits = inits.list,
	       family = 'NB',
	       tuning = tuning.list,
	       priors = prior.list,
	       accept.rate = 0.43,
	       n.omp.threads = 1,
	       verbose = FALSE,
	       n.report = 10,
	       n.burn = n.burn,
	       n.thin = n.thin,
	       n.chains = 1)

# Test to make sure it worked ---------
test_that("out is of class msAbund", {
  expect_s3_class(out, "msAbund")
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$muRE, TRUE)
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, dat$y)
})

# Check RE error ----------------------
test_that("random effect gives error when non-numeric", {
  data.list$covs$factor.1 <- factor(data.list$covs$factor.1)
  expect_error(out <- abund(formula = formula,
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
	       n.report = 10,
	       n.burn = n.burn,
	       n.thin = n.thin,
	       n.chains = 1))
  data.list$covs$factor.1 <- as.character(factor(data.list$covs$factor.1))
  expect_error(out <- abund(formula = formula,
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
	       n.report = 10,
	       n.burn = n.burn,
	       n.thin = n.thin,
	       n.chains = 1))
})

# Check Poisson -----------------------
test_that("Poisson works", {
  out <- msAbund(formula = formula,
  	     data = data.list, 
  	     n.batch = n.batch, 
  	     batch.length = batch.length, 
  	     inits = inits.list, 
  	     family = 'Poisson',
  	     tuning = tuning.list,
  	     priors = prior.list, 
  	     accept.rate = 0.43, 
  	     n.omp.threads = 1, 
  	     verbose = FALSE, 
  	     n.report = 10,
  	     n.burn = n.burn,
  	     n.thin = n.thin,
  	     n.chains = 1) 
  expect_s3_class(out, "msAbund")
})

# Check default priors ----------------
test_that("default priors, multiple chains, inits, burn, thin work", {
  out <- msAbund(formula = formula,
  	     data = data.list, 
  	     n.batch = n.batch, 
  	     batch.length = batch.length, 
  	     family = 'NB',
  	     accept.rate = 0.43, 
  	     n.omp.threads = 1, 
  	     verbose = FALSE, 
  	     n.report = 10,
  	     n.burn = n.burn,
  	     n.thin = n.thin,
  	     n.chains = 2) 
  expect_s3_class(out, "msAbund")
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {
  expect_output(msAbund(formula = formula,
  	     data = data.list, 
  	     n.batch = n.batch, 
  	     batch.length = batch.length, 
  	     family = 'NB',
  	     accept.rate = 0.43, 
  	     n.omp.threads = 1, 
  	     verbose = TRUE, 
  	     n.report = 10,
  	     n.burn = n.burn,
  	     n.thin = n.thin,
  	     n.chains = 1))
})

# Check waicAbund -----------------------
test_that("waicAbund works for msAbund", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicAbund(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for msAbund", {
  fitted.out <- fitted(out)
  expect_is(fitted.out, "array")
})

# Check predictions -------------------
test_that("predict works for msAbund", {
  X.0 <- abind::abind(dat$X, dat$X.re)
  dimnames(X.0)[[3]] <- c('(Intercept)', 'factor.1')
  pred.out <- predict(out, X.0)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$mu.0.samples), c(out$n.post * out$n.chains, n.sp, J, max(n.rep)))
  expect_equal(dim(pred.out$y.0.samples), c(out$n.post * out$n.chains, n.sp, J, max(n.rep)))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for msAbund", {
  n.post.samples <- out$n.post * out$n.chains
  ppc.out <- ppcAbund(out, 'chi-square', 0)
  expect_type(ppc.out, "list")
  expect_equal(nrow(ppc.out$fit.y), n.post.samples)
  expect_equal(nrow(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, n.sp, J, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, n.sp, J, max(n.rep)))
  
  ppc.out <- ppcAbund(out, 'freeman-tukey', 0)
  expect_type(ppc.out, "list")
  expect_equal(nrow(ppc.out$fit.y), n.post.samples)
  expect_equal(nrow(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, n.sp, J, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, n.sp, J, max(n.rep)))
})

# Multiple random intercepts ----------------------------------------------
J.x <- 10
J.y <- 10
J <- J.x * J.y
# n.rep <- sample(3, size = J, replace = TRUE)
n.rep <- rep(1, J)
n.sp <- 6
# Community-level covariate effects
beta.mean <- c(0)
p.abund <- length(beta.mean)
tau.sq.beta <- c(0.3)
# Random effects
mu.RE <- list(levels = c(20, 15),
              sigma.sq.mu = c(1, 0.5))
# Draw species-level effects from community means.
beta <- matrix(NA, nrow = n.sp, ncol = p.abund)
for (i in 1:p.abund) {
  beta[, i] <- rnorm(n.sp, beta.mean[i], sqrt(tau.sq.beta[i]))
}
sp <- FALSE
factor.model <- FALSE
kappa <- runif(n.sp, 0.1, 1)
family <- 'NB'

dat <- simMsAbund(J.x = J.x, J.y = J.y, n.rep = n.rep, n.sp = n.sp, beta = beta,
	          mu.RE = mu.RE, sp = sp, kappa = kappa, family = family,
                  factor.model = factor.model, n.factors = n.factors)

y <- dat$y
X <- dat$X
X.re <- dat$X.re
coords <- dat$coords

covs <- list(factor.1 = X.re[, , 1], 
             factor.2 = X.re[, , 2])

data.list <- list(y = y, 
                  covs = covs)

prior.list <- list(beta.comm.normal = list(mean = 0, var = 100),
                   kappa.unif = list(a = 0, b = 10),
                   tau.sq.beta.ig = list(a = .1, b = .1))

inits.list <- list(beta.comm = 0,
		   beta = 0,
		   tau.sq.beta = 1)
tuning.list <- list(beta = 0.1, beta.star = 0.1, kappa = 0.5)

n.batch <- 40
batch.length <- 25
n.burn <- 400
n.thin <- 2
n.chains <- 1
formula <- ~ (1 | factor.1) + (1 | factor.2)


out <- msAbund(formula = formula,
	       data = data.list,
	       n.batch = n.batch,
	       batch.length = batch.length,
	       inits = inits.list,
	       family = 'NB',
	       tuning = tuning.list,
	       priors = prior.list,
	       accept.rate = 0.43,
	       n.omp.threads = 1,
	       verbose = FALSE,
	       n.report = 10,
	       n.burn = n.burn,
	       n.thin = n.thin,
	       n.chains = 1)

# Test to make sure it worked ---------
test_that("out is of class msAbund", {
  expect_s3_class(out, "msAbund")
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$muRE, TRUE)
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, dat$y)
})

# Check RE error ----------------------
test_that("random effect gives error when non-numeric", {
  data.list$covs$factor.1 <- factor(data.list$covs$factor.1)
  expect_error(out <- abund(formula = formula,
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
	       n.report = 10,
	       n.burn = n.burn,
	       n.thin = n.thin,
	       n.chains = 1))
  data.list$covs$factor.1 <- as.character(factor(data.list$covs$factor.1))
  expect_error(out <- abund(formula = formula,
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
	       n.report = 10,
	       n.burn = n.burn,
	       n.thin = n.thin,
	       n.chains = 1))
})

# Check Poisson -----------------------
test_that("Poisson works", {
  out <- msAbund(formula = formula,
  	     data = data.list, 
  	     n.batch = n.batch, 
  	     batch.length = batch.length, 
  	     inits = inits.list, 
  	     family = 'Poisson',
  	     tuning = tuning.list,
  	     priors = prior.list, 
  	     accept.rate = 0.43, 
  	     n.omp.threads = 1, 
  	     verbose = FALSE, 
  	     n.report = 10,
  	     n.burn = n.burn,
  	     n.thin = n.thin,
  	     n.chains = 1) 
  expect_s3_class(out, "msAbund")
})

# Check default priors ----------------
test_that("default priors, multiple chains, inits, burn, thin work", {
  out <- msAbund(formula = formula,
  	     data = data.list, 
  	     n.batch = n.batch, 
  	     batch.length = batch.length, 
  	     family = 'NB',
  	     accept.rate = 0.43, 
  	     n.omp.threads = 1, 
  	     verbose = FALSE, 
  	     n.report = 10,
  	     n.burn = n.burn,
  	     n.thin = n.thin,
  	     n.chains = 2) 
  expect_s3_class(out, "msAbund")
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {
  expect_output(msAbund(formula = formula,
  	     data = data.list, 
  	     n.batch = n.batch, 
  	     batch.length = batch.length, 
  	     family = 'NB',
  	     accept.rate = 0.43, 
  	     n.omp.threads = 1, 
  	     verbose = TRUE, 
  	     n.report = 10,
  	     n.burn = n.burn,
  	     n.thin = n.thin,
  	     n.chains = 1))
})

# Check waicAbund -----------------------
test_that("waicAbund works for msAbund", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicAbund(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for msAbund", {
  fitted.out <- fitted(out)
  expect_is(fitted.out, "array")
})

# Check predictions -------------------
test_that("predict works for msAbund", {
  X.0 <- abind::abind(dat$X, dat$X.re)
  dimnames(X.0)[[3]] <- c('(Intercept)', 'factor.1', 'factor.2')
  pred.out <- predict(out, X.0)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$mu.0.samples), c(out$n.post * out$n.chains, n.sp, J, max(n.rep)))
  expect_equal(dim(pred.out$y.0.samples), c(out$n.post * out$n.chains, n.sp, J, max(n.rep)))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for msAbund", {
  n.post.samples <- out$n.post * out$n.chains
  ppc.out <- ppcAbund(out, 'chi-square', 0)
  expect_type(ppc.out, "list")
  expect_equal(nrow(ppc.out$fit.y), n.post.samples)
  expect_equal(nrow(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, n.sp, J, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, n.sp, J, max(n.rep)))
  
  ppc.out <- ppcAbund(out, 'freeman-tukey', 0)
  expect_type(ppc.out, "list")
  expect_equal(nrow(ppc.out$fit.y), n.post.samples)
  expect_equal(nrow(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, n.sp, J, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, n.sp, J, max(n.rep)))
})

# Random intercept and slope ----------------------------------------------
J.x <- 10
J.y <- 10
J <- J.x * J.y
n.rep <- sample(3, size = J, replace = TRUE)
# n.rep <- rep(1, J)
n.sp <- 6
# Community-level covariate effects
beta.mean <- c(0, -0.3)
p.abund <- length(beta.mean)
tau.sq.beta <- c(0.3, 0.7)
# Random effects
mu.RE <- list(levels = c(20, 15),
              sigma.sq.mu = c(1, 0.5), 
              beta.indx = c(1, 2))
# Draw species-level effects from community means.
beta <- matrix(NA, nrow = n.sp, ncol = p.abund)
for (i in 1:p.abund) {
  beta[, i] <- rnorm(n.sp, beta.mean[i], sqrt(tau.sq.beta[i]))
}
sp <- FALSE
factor.model <- FALSE
kappa <- runif(n.sp, 0.1, 1)
family <- 'NB'

dat <- simMsAbund(J.x = J.x, J.y = J.y, n.rep = n.rep, n.sp = n.sp, beta = beta,
	          mu.RE = mu.RE, sp = sp, kappa = kappa, family = family,
                  factor.model = factor.model, n.factors = n.factors)

y <- dat$y
X <- dat$X
X.re <- dat$X.re
coords <- dat$coords

covs <- list(factor.1 = X.re[, , 1], 
             factor.2 = X.re[, , 2], 
             cov.1 = X[, , 2])

data.list <- list(y = y, 
                  covs = covs)

prior.list <- list(beta.comm.normal = list(mean = 0, var = 100),
                   kappa.unif = list(a = 0, b = 10),
                   tau.sq.beta.ig = list(a = .1, b = .1))

inits.list <- list(beta.comm = 0,
		   beta = 0,
		   tau.sq.beta = 1)
tuning.list <- list(beta = 0.1, beta.star = 0.1, kappa = 0.5)

n.batch <- 40
batch.length <- 25
n.burn <- 400
n.thin <- 2
n.chains <- 1
formula <- ~ cov.1 + (1 | factor.1) + (cov.1 | factor.2)


out <- msAbund(formula = formula,
	       data = data.list,
	       n.batch = n.batch,
	       batch.length = batch.length,
	       inits = inits.list,
	       family = 'NB',
	       tuning = tuning.list,
	       priors = prior.list,
	       accept.rate = 0.43,
	       n.omp.threads = 1,
	       verbose = FALSE,
	       n.report = 10,
	       n.burn = n.burn,
	       n.thin = n.thin,
	       n.chains = 1)

# Test to make sure it worked ---------
test_that("out is of class msAbund", {
  expect_s3_class(out, "msAbund")
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$muRE, TRUE)
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, dat$y)
})

# Check RE error ----------------------
test_that("random effect gives error when non-numeric", {
  data.list$covs$factor.1 <- factor(data.list$covs$factor.1)
  expect_error(out <- abund(formula = formula,
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
	       n.report = 10,
	       n.burn = n.burn,
	       n.thin = n.thin,
	       n.chains = 1))
  data.list$covs$factor.1 <- as.character(factor(data.list$covs$factor.1))
  expect_error(out <- abund(formula = formula,
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
	       n.report = 10,
	       n.burn = n.burn,
	       n.thin = n.thin,
	       n.chains = 1))
})

# Check Poisson -----------------------
test_that("Poisson works", {
  out <- msAbund(formula = formula,
  	     data = data.list, 
  	     n.batch = n.batch, 
  	     batch.length = batch.length, 
  	     inits = inits.list, 
  	     family = 'Poisson',
  	     tuning = tuning.list,
  	     priors = prior.list, 
  	     accept.rate = 0.43, 
  	     n.omp.threads = 1, 
  	     verbose = FALSE, 
  	     n.report = 10,
  	     n.burn = n.burn,
  	     n.thin = n.thin,
  	     n.chains = 1) 
  expect_s3_class(out, "msAbund")
})

# Check default priors ----------------
test_that("default priors, multiple chains, inits, burn, thin work", {
  out <- msAbund(formula = formula,
  	     data = data.list, 
  	     n.batch = n.batch, 
  	     batch.length = batch.length, 
  	     family = 'NB',
  	     accept.rate = 0.43, 
  	     n.omp.threads = 1, 
  	     verbose = FALSE, 
  	     n.report = 10,
  	     n.burn = n.burn,
  	     n.thin = n.thin,
  	     n.chains = 2) 
  expect_s3_class(out, "msAbund")
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {
  expect_output(msAbund(formula = formula,
  	     data = data.list, 
  	     n.batch = n.batch, 
  	     batch.length = batch.length, 
  	     family = 'NB',
  	     accept.rate = 0.43, 
  	     n.omp.threads = 1, 
  	     verbose = TRUE, 
  	     n.report = 10,
  	     n.burn = n.burn,
  	     n.thin = n.thin,
  	     n.chains = 1))
})

# Check waicAbund -----------------------
test_that("waicAbund works for msAbund", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicAbund(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for msAbund", {
  fitted.out <- fitted(out)
  expect_is(fitted.out, "array")
})

# Check predictions -------------------
test_that("predict works for msAbund", {
  X.0 <- abind::abind(dat$X, dat$X.re)
  dimnames(X.0)[[3]] <- c('(Intercept)', 'cov.1', 'factor.1', 'factor.2')
  pred.out <- predict(out, X.0)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$mu.0.samples), c(out$n.post * out$n.chains, n.sp, J, max(n.rep)))
  expect_equal(dim(pred.out$y.0.samples), c(out$n.post * out$n.chains, n.sp, J, max(n.rep)))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for msAbund", {
  n.post.samples <- out$n.post * out$n.chains
  ppc.out <- ppcAbund(out, 'chi-square', 0)
  expect_type(ppc.out, "list")
  expect_equal(nrow(ppc.out$fit.y), n.post.samples)
  expect_equal(nrow(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, n.sp, J, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, n.sp, J, max(n.rep)))
  
  ppc.out <- ppcAbund(out, 'freeman-tukey', 0)
  expect_type(ppc.out, "list")
  expect_equal(nrow(ppc.out$fit.y), n.post.samples)
  expect_equal(nrow(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, n.sp, J, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, n.sp, J, max(n.rep)))
})

# Random intercept and slope with the same factor -------------------------
J.x <- 10
J.y <- 10
J <- J.x * J.y
n.rep <- sample(3, size = J, replace = TRUE)
# n.rep <- rep(1, J)
n.sp <- 6
# Community-level covariate effects
beta.mean <- c(0, -0.3)
p.abund <- length(beta.mean)
tau.sq.beta <- c(0.3, 0.7)
# Random effects
mu.RE <- list(levels = c(20, 15),
              sigma.sq.mu = c(1, 0.5), 
              beta.indx = c(1, 2))
# Draw species-level effects from community means.
beta <- matrix(NA, nrow = n.sp, ncol = p.abund)
for (i in 1:p.abund) {
  beta[, i] <- rnorm(n.sp, beta.mean[i], sqrt(tau.sq.beta[i]))
}
sp <- FALSE
factor.model <- FALSE
kappa <- runif(n.sp, 0.1, 1)
family <- 'NB'

dat <- simMsAbund(J.x = J.x, J.y = J.y, n.rep = n.rep, n.sp = n.sp, beta = beta,
	          mu.RE = mu.RE, sp = sp, kappa = kappa, family = family,
                  factor.model = factor.model, n.factors = n.factors)

y <- dat$y
X <- dat$X
X.re <- dat$X.re
coords <- dat$coords

covs <- list(factor.1 = X.re[, , 1], 
             factor.2 = X.re[, , 2], 
             cov.1 = X[, , 2])

data.list <- list(y = y, 
                  covs = covs)

prior.list <- list(beta.comm.normal = list(mean = 0, var = 100),
                   kappa.unif = list(a = 0, b = 10),
                   tau.sq.beta.ig = list(a = .1, b = .1))

inits.list <- list(beta.comm = 0,
		   beta = 0,
		   tau.sq.beta = 1)
tuning.list <- list(beta = 0.1, beta.star = 0.1, kappa = 0.5)

n.batch <- 40
batch.length <- 25
n.burn <- 400
n.thin <- 2
n.chains <- 1
formula <- ~ cov.1 + (1 | factor.1) + (cov.1 | factor.2) + (cov.1 | factor.1)

out <- msAbund(formula = formula,
	       data = data.list,
	       n.batch = n.batch,
	       batch.length = batch.length,
	       inits = inits.list,
	       family = 'NB',
	       tuning = tuning.list,
	       priors = prior.list,
	       accept.rate = 0.43,
	       n.omp.threads = 1,
	       verbose = FALSE,
	       n.report = 10,
	       n.burn = n.burn,
	       n.thin = n.thin,
	       n.chains = 1)

# Test to make sure it worked ---------
test_that("out is of class msAbund", {
  expect_s3_class(out, "msAbund")
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$muRE, TRUE)
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, dat$y)
})

# Check RE error ----------------------
test_that("random effect gives error when non-numeric", {
  data.list$covs$factor.1 <- factor(data.list$covs$factor.1)
  expect_error(out <- abund(formula = formula,
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
	       n.report = 10,
	       n.burn = n.burn,
	       n.thin = n.thin,
	       n.chains = 1))
  data.list$covs$factor.1 <- as.character(factor(data.list$covs$factor.1))
  expect_error(out <- abund(formula = formula,
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
	       n.report = 10,
	       n.burn = n.burn,
	       n.thin = n.thin,
	       n.chains = 1))
})

# Check Poisson -----------------------
test_that("Poisson works", {
  out <- msAbund(formula = formula,
  	     data = data.list, 
  	     n.batch = n.batch, 
  	     batch.length = batch.length, 
  	     inits = inits.list, 
  	     family = 'Poisson',
  	     tuning = tuning.list,
  	     priors = prior.list, 
  	     accept.rate = 0.43, 
  	     n.omp.threads = 1, 
  	     verbose = FALSE, 
  	     n.report = 10,
  	     n.burn = n.burn,
  	     n.thin = n.thin,
  	     n.chains = 1) 
  expect_s3_class(out, "msAbund")
})

# Check default priors ----------------
test_that("default priors, multiple chains, inits, burn, thin work", {
  out <- msAbund(formula = formula,
  	     data = data.list, 
  	     n.batch = n.batch, 
  	     batch.length = batch.length, 
  	     family = 'NB',
  	     accept.rate = 0.43, 
  	     n.omp.threads = 1, 
  	     verbose = FALSE, 
  	     n.report = 10,
  	     n.burn = n.burn,
  	     n.thin = n.thin,
  	     n.chains = 2) 
  expect_s3_class(out, "msAbund")
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {
  expect_output(msAbund(formula = formula,
  	     data = data.list, 
  	     n.batch = n.batch, 
  	     batch.length = batch.length, 
  	     family = 'NB',
  	     accept.rate = 0.43, 
  	     n.omp.threads = 1, 
  	     verbose = TRUE, 
  	     n.report = 10,
  	     n.burn = n.burn,
  	     n.thin = n.thin,
  	     n.chains = 1))
})

# Check waicAbund -----------------------
test_that("waicAbund works for msAbund", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicAbund(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for msAbund", {
  fitted.out <- fitted(out)
  expect_is(fitted.out, "array")
})

# Check predictions -------------------
test_that("predict works for msAbund", {
  X.0 <- abind::abind(dat$X, dat$X.re)
  dimnames(X.0)[[3]] <- c('(Intercept)', 'cov.1', 'factor.1', 'factor.2')
  pred.out <- predict(out, X.0)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$mu.0.samples), c(out$n.post * out$n.chains, n.sp, J, max(n.rep)))
  expect_equal(dim(pred.out$y.0.samples), c(out$n.post * out$n.chains, n.sp, J, max(n.rep)))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for msAbund", {
  n.post.samples <- out$n.post * out$n.chains
  ppc.out <- ppcAbund(out, 'chi-square', 0)
  expect_type(ppc.out, "list")
  expect_equal(nrow(ppc.out$fit.y), n.post.samples)
  expect_equal(nrow(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, n.sp, J, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, n.sp, J, max(n.rep)))
  
  ppc.out <- ppcAbund(out, 'freeman-tukey', 0)
  expect_type(ppc.out, "list")
  expect_equal(nrow(ppc.out$fit.y), n.post.samples)
  expect_equal(nrow(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, n.sp, J, max(n.rep)))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, n.sp, J, max(n.rep)))
})
