# Test sfMsAbund.R ----------------------------------------------------------
skip_on_cran()

# Intercept only ----------------------------------------------------------
set.seed(10101)
J.x <- 10
J.y <- 10
J <- J.x * J.y
n.rep <- rep(1, J)
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
sp <- TRUE
tau.sq <- runif(n.sp, 0.1, 3)
n.factors <- 3
factor.model <- TRUE
phi <- runif(n.factors, 3/1, 3 / .3)
nu <- rep(1, n.factors)
cov.model <- 'matern'
family <- 'zi-Gaussian'
z <- matrix(rbinom(J * n.sp, 1, 0.5), n.sp, J)

dat <- simMsAbund(J.x = J.x, J.y = J.y, n.rep = n.rep, n.sp = n.sp, beta = beta,
	          mu.RE = mu.RE, sp = sp, tau.sq = tau.sq, family = family,
                  factor.model = factor.model, n.factors = n.factors, 
                  cov.model = cov.model, phi = phi, nu = nu, z = z)

pred.indx <- sample(1:J, round(J * .25), replace = FALSE)
y <- dat$y[, -pred.indx, drop = FALSE]
y.0 <- dat$y[, pred.indx, drop = FALSE]
X <- as.matrix(dat$X[-pred.indx])
coords <- as.matrix(dat$coords[-pred.indx, , drop = FALSE])
X.0 <- as.matrix(dat$X[pred.indx])
coords.0 <- as.matrix(dat$coords[pred.indx, , drop = FALSE])
z.0 <- z[, pred.indx]
z.fit <- z[, -pred.indx]

data.list <- list(y = y, coords = coords, z = z.fit)

prior.list <- list(beta.comm.normal = list(mean = 0, var = 100),
                   tau.sq.unif = list(a = 2, b = 1),
		   phi.unif = list(3 / 1, 3 / .1),
		   nu.unif = list(0.1, 3),
                   tau.sq.beta.ig = list(a = .1, b = .1))

inits.list <- list(beta.comm = 0,
		   beta = 0, tau.sq = 1,
		   phi = 3 / .5, nu = 2,
		   tau.sq.beta = 1)

n.batch <- 40
batch.length <- 25
n.burn <- 400
n.thin <- 2
n.chains <- 1
formula <- ~ 1

out <- sfMsAbund(formula = formula,
	       data = data.list,
	       n.batch = n.batch,
	       batch.length = batch.length,
	       inits = inits.list,
	       family = 'Gaussian',
	       NNGP = TRUE, 
	       n.neighbors = 5,
	       n.factors = n.factors,
	       priors = prior.list,
	       accept.rate = 0.43,
	       n.omp.threads = 1,
	       verbose = FALSE,
	       n.report = 10,
	       n.burn = n.burn,
	       n.thin = n.thin,
	       n.chains = 1)

# Test to make sure it worked ---------
test_that("out is of class sfMsAbund", {
  expect_s3_class(out, "sfMsAbund")
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$muRE, FALSE)
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, y)
})

# Check zi Gaussian -------------------
test_that("zi-Gaussian works", {
  out <- sfMsAbund(formula = formula,
  	     data = data.list, 
  	     n.batch = n.batch, 
  	     batch.length = batch.length, 
  	     inits = inits.list, 
  	     family = 'zi-Gaussian',
	     NNGP = TRUE, 
	     n.neighbors = 5, 
	     cov.model = 'exponential',
	     n.factors = n.factors, 
  	     priors = prior.list, 
  	     accept.rate = 0.43, 
  	     n.omp.threads = 1, 
  	     verbose = FALSE, 
  	     n.report = 10,
  	     n.burn = n.burn,
  	     n.thin = n.thin,
  	     n.chains = 1) 
  expect_s3_class(out, "sfMsAbund")
})

# Check default priors ----------------
test_that("default priors, multiple chains, inits, burn, thin work", {
  out <- sfMsAbund(formula = formula,
  	     data = data.list, 
  	     n.batch = n.batch, 
  	     batch.length = batch.length, 
  	     family = 'Gaussian',
	     NNGP = TRUE, 
	     cov.model = 'spherical',
	     n.neighbors = 5, 
	     n.factors = n.factors, 
  	     accept.rate = 0.43, 
  	     n.omp.threads = 1, 
  	     verbose = FALSE, 
  	     n.report = 10,
  	     n.burn = n.burn,
  	     n.thin = n.thin,
  	     n.chains = 2) 
  expect_s3_class(out, "sfMsAbund")
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {
  expect_output(sfMsAbund(formula = formula,
  	     data = data.list, 
  	     n.batch = n.batch, 
  	     batch.length = batch.length, 
  	     family = 'Gaussian',
	     cov.model = 'gaussian',
	     NNGP = TRUE, 
	     n.neighbors = 5, 
	     n.factors = n.factors, 
  	     accept.rate = 0.43, 
  	     n.omp.threads = 1, 
  	     verbose = TRUE, 
  	     n.report = 10,
  	     n.burn = n.burn,
  	     n.thin = n.thin,
  	     n.chains = 1))
})

# Check waicAbund -----------------------
test_that("waicAbund works for sfMsAbund", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicAbund(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for sfMsAbund", {
  fitted.out <- fitted(out)
  expect_is(fitted.out, "array")
})

# Check predictions -------------------
test_that("predict works for sfMsAbund", {
  pred.out <- predict(out, X.0, coords.0, verbose = FALSE)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$mu.0.samples), c(out$n.post * out$n.chains, n.sp, nrow(coords.0)))
  expect_equal(dim(pred.out$y.0.samples), c(out$n.post * out$n.chains, n.sp, nrow(coords.0)))
})

# Covariates --------------------------------------------------------------
J.x <- 10
J.y <- 10
J <- J.x * J.y
n.rep <- rep(1, J)
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
sp <- TRUE
tau.sq <- runif(n.sp, 0.1, 3)
n.factors <- 3
phi <- runif(n.factors, 3/1, 3 / .3)
nu <- rep(1, n.factors)
cov.model <- 'matern'
family <- 'zi-Gaussian'
z <- matrix(rbinom(J * n.sp, 1, 0.5), n.sp, J)

dat <- simMsAbund(J.x = J.x, J.y = J.y, n.rep = n.rep, n.sp = n.sp, beta = beta,
	          mu.RE = mu.RE, sp = sp, tau.sq = tau.sq, family = family,
                  factor.model = factor.model, n.factors = n.factors, 
                  cov.model = cov.model, phi = phi, nu = nu, z = z)

y <- dat$y
X <- dat$X
coords <- dat$coords

covs <- list(cov.1 = X[, 2], 
	     cov.2 = X[, 3])

data.list <- list(y = y, 
                  covs = covs, coords = coords, z = z)

prior.list <- list(beta.comm.normal = list(mean = 0, var = 100),
                   tau.sq.ig = list(2, 1),
		   nu.unif = list(a = 0.1, b = 3),
                   tau.sq.beta.ig = list(a = .1, b = .1))

inits.list <- list(beta.comm = 0,
		   beta = 0, tau.sq = tau.sq,
		   tau.sq.beta = 1)

n.batch <- 40
batch.length <- 25
n.burn <- 400
n.thin <- 2
n.chains <- 1
formula <- ~ cov.1 + cov.2

out <- sfMsAbund(formula = formula,
	       data = data.list,
	       n.batch = n.batch,
	       batch.length = batch.length,
	       inits = inits.list,
	       family = 'Gaussian',
	       NNGP = TRUE, 
	       n.neighbors = 5, 
	       cov.model = 'matern',
	       n.factors = n.factors, 
	       priors = prior.list,
	       accept.rate = 0.43,
	       n.omp.threads = 1,
	       verbose = FALSE,
	       n.report = 10,
	       n.burn = n.burn,
	       n.thin = n.thin,
	       n.chains = 1)

# Test to make sure it worked ---------
test_that("out is of class sfMsAbund", {
  expect_s3_class(out, "sfMsAbund")
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
  out <- sfMsAbund(formula = formula,
  	     data = data.list, 
  	     n.batch = n.batch, 
  	     batch.length = batch.length, 
  	     inits = inits.list, 
  	     family = 'Gaussian',
	     NNGP = TRUE, 
	     cov.model = 'spherical',
	     n.neighbors = 5, 
	     n.factors = n.factors, 
  	     priors = prior.list, 
  	     accept.rate = 0.43, 
  	     n.omp.threads = 1, 
  	     verbose = FALSE, 
  	     n.report = 10,
  	     n.burn = n.burn,
  	     n.thin = n.thin,
  	     n.chains = 1) 
  expect_s3_class(out, "sfMsAbund")
})

# Check default priors ----------------
test_that("default priors, multiple chains, inits, burn, thin work", {
  out <- sfMsAbund(formula = formula,
  	     data = data.list, 
  	     n.batch = n.batch, 
  	     batch.length = batch.length, 
  	     family = 'Gaussian',
	     NNGP = TRUE, 
	     cov.model = 'gaussian',
	     n.neighbors = 5, 
	     n.factors = n.factors, 
  	     accept.rate = 0.43, 
  	     n.omp.threads = 1, 
  	     verbose = FALSE, 
  	     n.report = 10,
  	     n.burn = n.burn,
  	     n.thin = n.thin,
  	     n.chains = 2) 
  expect_s3_class(out, "sfMsAbund")
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {
  expect_output(sfMsAbund(formula = formula,
  	     data = data.list, 
  	     n.batch = n.batch, 
  	     batch.length = batch.length, 
  	     family = 'Gaussian',
	     NNGP = TRUE, 
	     cov.model = 'exponential',
	     n.neighbors = 5, 
	     n.factors = n.factors, 
  	     accept.rate = 0.43, 
  	     n.omp.threads = 1, 
  	     verbose = TRUE, 
  	     n.report = 10,
  	     n.burn = n.burn,
  	     n.thin = n.thin,
  	     n.chains = 1))
})

# Check waicAbund -----------------------
test_that("waicAbund works for sfMsAbund", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicAbund(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for sfMsAbund", {
  fitted.out <- fitted(out)
  expect_is(fitted.out, "array")
})

# Check predictions -------------------
test_that("predict works for sfMsAbund", {
  pred.out <- predict(out, X.0, coords.0 = coords.0, verbose = FALSE)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$mu.0.samples), c(out$n.post * out$n.chains, dim(y.0)))
  expect_equal(dim(pred.out$y.0.samples), c(out$n.post * out$n.chains, dim(y.0)))
})

# Random intercept and slope ----------------------------------------------
J.x <- 10
J.y <- 10
J <- J.x * J.y
n.rep <- rep(1, J)
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
sp <- TRUE
tau.sq <- runif(n.sp, 0.1, 3)
n.factors <- 3
phi <- runif(n.factors, 3/1, 3 / .3)
nu <- rep(1, n.factors)
cov.model <- 'matern'
family <- 'zi-Gaussian'
z <- matrix(rbinom(J * n.sp, 1, 0.5), n.sp, J)

dat <- simMsAbund(J.x = J.x, J.y = J.y, n.rep = n.rep, n.sp = n.sp, beta = beta,
	          mu.RE = mu.RE, sp = sp, tau.sq = tau.sq, family = family,
                  factor.model = factor.model, n.factors = n.factors, 
                  cov.model = cov.model, phi = phi, nu = nu, z = z)

pred.indx <- sample(1:J, round(J * .25), replace = FALSE)
y <- dat$y[, -pred.indx, drop = FALSE]
y.0 <- dat$y[, pred.indx, drop = FALSE]
X <- dat$X[-pred.indx, , drop = FALSE]
X.re <- dat$X.re[-pred.indx, , drop = FALSE]
coords <- as.matrix(dat$coords[-pred.indx, , drop = FALSE])
X.0 <- dat$X[pred.indx, , drop = FALSE]
X.re.0 <- dat$X.re[pred.indx, , drop = FALSE]
coords.0 <- as.matrix(dat$coords[pred.indx, , drop = FALSE])
z.0 <- z[, pred.indx]
z.fit <- z[, -pred.indx]

covs <- list(factor.1 = X.re[, 1], 
             factor.2 = X.re[, 2], 
             cov.1 = X[, 2])

data.list <- list(y = y, 
                  covs = covs, coords = coords, z = z.fit)

prior.list <- list(beta.comm.normal = list(mean = 0, var = 100),
                   tau.sq.ig = list(a = 2, b = 1),
		   nu.unif = list(a = 0.1, b = 5),
                   tau.sq.beta.ig = list(a = .1, b = .1))

inits.list <- list(beta.comm = 0,
		   beta = 0,
		   tau.sq.beta = 1)

n.batch <- 40
batch.length <- 25
n.burn <- 400
n.thin <- 2
n.chains <- 1
formula <- ~ cov.1 + (1 | factor.1) + (cov.1 | factor.2)


out <- sfMsAbund(formula = formula,
	       data = data.list,
	       n.batch = n.batch,
	       batch.length = batch.length,
	       inits = inits.list,
	       family = 'zi-Gaussian',
	       cov.model = 'matern',
	     NNGP = TRUE, 
	     n.neighbors = 5, 
	     n.factors = n.factors, 
	       priors = prior.list,
	       accept.rate = 0.43,
	       n.omp.threads = 1,
	       verbose = FALSE,
	       n.report = 10,
	       n.burn = n.burn,
	       n.thin = n.thin,
	       n.chains = 1)

# Test to make sure it worked ---------
test_that("out is of class sfMsAbund", {
  expect_s3_class(out, "sfMsAbund")
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$muRE, TRUE)
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, y)
})

# Check RE error ----------------------
test_that("random effect gives error when non-numeric", {
  data.list$covs$factor.1 <- factor(data.list$covs$factor.1)
  expect_error(out <- sfMsAbund(formula = formula,
	       data = data.list,
	       n.batch = n.batch,
	       batch.length = batch.length,
	       inits = inits.list,
	       family = 'zi-Gaussian',
	     NNGP = TRUE, 
	     n.neighbors = 5, 
	     n.factors = n.factors, 
	       priors = prior.list,
	       accept.rate = 0.43,
	       n.omp.threads = 1,
	       verbose = FALSE,
	       n.report = 10,
	       n.burn = n.burn,
	       n.thin = n.thin,
	       n.chains = 1))
  data.list$covs$factor.1 <- as.character(factor(data.list$covs$factor.1))
  expect_error(out <- sfMsAbund(formula = formula,
	       data = data.list,
	       n.batch = n.batch,
	       batch.length = batch.length,
	       inits = inits.list,
	       family = 'zi-Gaussian',
	     NNGP = TRUE, 
	     n.neighbors = 5, 
	     n.factors = n.factors, 
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
test_that("Gaussian works", {
  out <- sfMsAbund(formula = formula,
  	     data = data.list, 
  	     n.batch = n.batch, 
  	     batch.length = batch.length, 
  	     inits = inits.list, 
  	     family = 'Gaussian',
	     NNGP = TRUE, 
	     cov.model = 'spherical',
	     n.neighbors = 5, 
	     n.factors = n.factors, 
  	     priors = prior.list, 
  	     accept.rate = 0.43, 
  	     n.omp.threads = 1, 
  	     verbose = FALSE, 
  	     n.report = 10,
  	     n.burn = n.burn,
  	     n.thin = n.thin,
  	     n.chains = 1) 
  expect_s3_class(out, "sfMsAbund")
})

# Check default priors ----------------
test_that("default priors, multiple chains, inits, burn, thin work", {
  out <- sfMsAbund(formula = formula,
  	     data = data.list, 
  	     n.batch = n.batch, 
  	     batch.length = batch.length, 
  	     family = 'zi-Gaussian',
	     cov.model = 'gaussian',
	     NNGP = TRUE, 
	     n.neighbors = 5, 
	     n.factors = n.factors, 
  	     accept.rate = 0.43, 
  	     n.omp.threads = 1, 
  	     verbose = FALSE, 
  	     n.report = 10,
  	     n.burn = n.burn,
  	     n.thin = n.thin,
  	     n.chains = 2) 
  expect_s3_class(out, "sfMsAbund")
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {
  expect_output(sfMsAbund(formula = formula,
  	     data = data.list, 
  	     n.batch = n.batch, 
  	     batch.length = batch.length, 
  	     family = 'zi-Gaussian',
	     NNGP = TRUE, 
	     cov.model = 'exponential',
	     n.neighbors = 5, 
	     n.factors = n.factors, 
  	     accept.rate = 0.43, 
  	     n.omp.threads = 1, 
  	     verbose = TRUE, 
  	     n.report = 10,
  	     n.burn = n.burn,
  	     n.thin = n.thin,
  	     n.chains = 1))
})

# Check waicAbund -----------------------
test_that("waicAbund works for sfMsAbund", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicAbund(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for sfMsAbund", {
  fitted.out <- fitted(out)
  expect_is(fitted.out, "array")
})

# Check predictions -------------------
test_that("predict works for sfMsAbund", {
  X.0 <- cbind(X.0, X.re.0)
  colnames(X.0) <- c('(Intercept)', 'cov.1', 'factor.1', 'factor.2')
  z.0.samples <- array(NA, dim = c(nrow(out$beta.samples), n.sp, ncol(z.0)))
  for (i in 1:nrow(z.0.samples)) {
    z.0.samples[i, , ] <- z.0
  }
  pred.out <- predict(out, X.0, coords.0 = coords.0, z.0.samples = z.0.samples, verbose = FALSE)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$mu.0.samples), c(out$n.post * out$n.chains, dim(y.0)))
  expect_equal(dim(pred.out$y.0.samples), c(out$n.post * out$n.chains, dim(y.0)))
})
