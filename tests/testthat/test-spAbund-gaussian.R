# Test spAbund-Gaussian.R -------------------------------------------------

skip_on_cran()

# Intercept Only ----------------------------------------------------------
J.x <- 15
J.y <- 15
J <- J.x * J.y
n.rep <- rep(1, J)
beta <- c(-1)
p.abund <- length(beta)
mu.RE <- list()
sp <- TRUE 
cov.model <- 'matern'
family <- 'zi-Gaussian' 
phi <- 3 / .6
nu <- 3
sigma.sq <- 0.4
tau.sq <- 0.2
z <- rbinom(J, 1, 0.5)
dat <- simAbund(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta,
	        mu.RE = mu.RE, sp = sp, phi = phi, tau.sq = tau.sq, nu = nu,
		sigma.sq = sigma.sq, cov.model = cov.model, family = family, z = z)

pred.indx <- sample(1:J, round(J * .25), replace = FALSE)
y <- dat$y[-pred.indx]
y.0 <- dat$y[pred.indx]
X <- as.matrix(dat$X[-pred.indx])
coords <- as.matrix(dat$coords[-pred.indx, , drop = FALSE])
X.0 <- as.matrix(dat$X[pred.indx])
coords.0 <- as.matrix(dat$coords[pred.indx, , drop = FALSE])
z.fit <- z[-pred.indx]
z.0 <- z[pred.indx]

data.list <- list(y = y, coords = coords, z = z.fit) 

# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 100), 
                   tau.sq.ig = c(2, 1), 
                   sigma.sq.ig = c(2, 1),
                   phi.unif = c(3 / 1, 3 / .1), 
                   nu.unif = c(0.1, 5)) 
# Starting values
inits.list <- list(beta = 0,
                   sigma.sq.mu = 0.3,
                   tau.sq = 1, phi = 3 / .5, sigma.sq = 0.5, nu = 2)

# Small
n.batch <- 10
batch.length <- 25
n.burn <- 50 
n.thin <- 2
n.chains <- 1
formula <- ~ 1

out <- spAbund(formula = formula,
	       data = data.list, 
	       n.batch = n.batch, 
	       batch.length = batch.length, 
	       inits = inits.list, 
	       family = 'Gaussian',
	       priors = prior.list, 
	       cov.model = 'matern',
	       NNGP = TRUE,
	       n.neighbors = 5,
	       accept.rate = 0.43, 
	       n.omp.threads = 1, 
	       verbose = FALSE, 
	       n.report = 10,
	       n.burn = n.burn,
	       n.thin = n.thin,
	       n.chains = 1) 

# Test to make sure it worked ---------
test_that("out is of class spAbund", {
  expect_s3_class(out, "spAbund")
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$muRE, FALSE)
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, y)
})

# Check zi-Gaussian -----------------------
test_that("zi-Gaussian works", {
  out <- spAbund(formula = formula,
  	         data = data.list, 
  	         n.batch = n.batch, 
  	         batch.length = batch.length, 
  	         inits = inits.list, 
  	         family = 'zi-Gaussian',
  	         priors = prior.list, 
	         NNGP = TRUE,
	         n.neighbors = 5,
		 cov.model = 'spherical',
  	         accept.rate = 0.43, 
  	         n.omp.threads = 1, 
  	         verbose = FALSE, 
  	         n.report = 10,
  	         n.burn = n.burn,
  	         n.thin = n.thin,
  	         n.chains = 1) 
  expect_s3_class(out, "spAbund")
})

# Check default priors ----------------
test_that("default priors, multiple chains, inits, burn, thin work", {
  out <- spAbund(formula = formula,
  	     data = data.list, 
  	     n.batch = n.batch, 
  	     batch.length = batch.length, 
  	     family = 'Gaussian',
  	     accept.rate = 0.43, 
  	     n.omp.threads = 1, 
	     NNGP = TRUE,
	     n.neighbors = 5,
	     cov.model = 'exponential',
  	     verbose = FALSE, 
  	     n.report = 10,
  	     n.burn = n.burn,
  	     n.thin = n.thin,
  	     n.chains = 2) 
  expect_s3_class(out, "spAbund")
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {
  expect_output(spAbund(formula = formula,
  	     data = data.list, 
  	     n.batch = n.batch, 
  	     batch.length = batch.length, 
  	     family = 'Gaussian',
	     NNGP = TRUE,
	     n.neighbors = 5,
	     cov.model = 'gaussian',
  	     accept.rate = 0.43, 
  	     n.omp.threads = 1, 
  	     verbose = TRUE, 
  	     n.report = 10,
  	     n.burn = n.burn,
  	     n.thin = n.thin,
  	     n.chains = 1))
})

# Check waicAbund -----------------------
test_that("waicAbund works for spAbund", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicAbund(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for spAbund", {
  fitted.out <- fitted(out)
  expect_is(fitted.out, "mcmc")
})

# Check predictions -------------------
test_that("predict works for spAbund", {
  pred.out <- predict(out, X.0, coords.0, verbose = FALSE)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$mu.0.samples), c(out$n.post * out$n.chains, length(y.0)))
  expect_equal(dim(pred.out$y.0.samples), c(out$n.post * out$n.chains, length(y.0)))
})

# Covariates --------------------------------------------------------------
J.x <- 15
J.y <- 15
J <- J.x * J.y
n.rep <- rep(1, J)
beta <- c(-1, 0.3, 0.4)
p.abund <- length(beta)
mu.RE <- list()
sp <- TRUE 
cov.model <- 'matern'
family <- 'zi-Gaussian' 
phi <- 3 / .6
nu <- 3
sigma.sq <- 0.4
tau.sq <- 0.2
z <- rbinom(J, 1, 0.5)
dat <- simAbund(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta,
	        mu.RE = mu.RE, sp = sp, phi = phi, tau.sq = tau.sq, nu = nu, z = z,
		sigma.sq = sigma.sq, cov.model = cov.model, family = family)

pred.indx <- sample(1:J, round(J * .25), replace = FALSE)
y <- dat$y[-pred.indx]
y.0 <- dat$y[pred.indx]
X <- dat$X[-pred.indx, , drop = FALSE]
coords <- as.matrix(dat$coords[-pred.indx, , drop = FALSE])
X.0 <- dat$X[pred.indx, , drop = FALSE]
coords.0 <- as.matrix(dat$coords[pred.indx, , drop = FALSE])
z.fit <- z[-pred.indx]
z.0 <- z[pred.indx]

covs <- list(cov.1 = X[, 2],
	     cov.2 = X[, 3])
data.list <- list(y = y, covs = covs, coords = coords, z = z.fit)

# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 100), 
                   tau.sq.ig = c(2, 1), 
                   sigma.sq.ig = c(2, 1),
                   phi.unif = c(3 / 1, 3 / .1), 
                   nu.unif = c(0.1, 5)) 
# Starting values
inits.list <- list(beta = 0,
                   sigma.sq.mu = 0.3,
                   tau.sq = 1, phi = 3 / .5, sigma.sq = 0.5, nu = 2)

# Small
n.batch <- 10
batch.length <- 25
n.burn <- 50 
n.thin <- 2
n.chains <- 1
formula <- ~ cov.1 + cov.2

out <- spAbund(formula = formula,
	       data = data.list, 
	       n.batch = n.batch, 
	       batch.length = batch.length, 
	       inits = inits.list, 
	       family = 'Gaussian',
	       priors = prior.list, 
	       cov.model = 'matern',
	       NNGP = TRUE,
	       n.neighbors = 5,
	       accept.rate = 0.43, 
	       n.omp.threads = 1, 
	       verbose = FALSE, 
	       n.report = 10,
	       n.burn = n.burn,
	       n.thin = n.thin,
	       n.chains = 1) 

# Test to make sure it worked ---------
test_that("out is of class spAbund", {
  expect_s3_class(out, "spAbund")
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$muRE, FALSE)
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, y)
})

# Check zi-Gaussian -----------------------
test_that("zi-Gaussian works", {
  out <- spAbund(formula = formula,
  	         data = data.list, 
  	         n.batch = n.batch, 
  	         batch.length = batch.length, 
  	         inits = inits.list, 
  	         family = 'zi-Gaussian',
  	         priors = prior.list, 
	         NNGP = TRUE,
	         n.neighbors = 5,
		 cov.model = 'spherical',
  	         accept.rate = 0.43, 
  	         n.omp.threads = 1, 
  	         verbose = FALSE, 
  	         n.report = 10,
  	         n.burn = n.burn,
  	         n.thin = n.thin,
  	         n.chains = 1) 
  expect_s3_class(out, "spAbund")
})

# Check default priors ----------------
test_that("default priors, multiple chains, inits, burn, thin work", {
  out <- spAbund(formula = formula,
  	     data = data.list, 
  	     n.batch = n.batch, 
  	     batch.length = batch.length, 
  	     family = 'Gaussian',
  	     accept.rate = 0.43, 
  	     n.omp.threads = 1, 
	     NNGP = TRUE,
	     n.neighbors = 5,
	     cov.model = 'exponential',
  	     verbose = FALSE, 
  	     n.report = 10,
  	     n.burn = n.burn,
  	     n.thin = n.thin,
  	     n.chains = 2) 
  expect_s3_class(out, "spAbund")
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {
  expect_output(spAbund(formula = formula,
  	     data = data.list, 
  	     n.batch = n.batch, 
  	     batch.length = batch.length, 
  	     family = 'Gaussian',
	     NNGP = TRUE,
	     n.neighbors = 5,
	     cov.model = 'gaussian',
  	     accept.rate = 0.43, 
  	     n.omp.threads = 1, 
  	     verbose = TRUE, 
  	     n.report = 10,
  	     n.burn = n.burn,
  	     n.thin = n.thin,
  	     n.chains = 1))
})

# Check waicAbund -----------------------
test_that("waicAbund works for spAbund", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicAbund(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for spAbund", {
  fitted.out <- fitted(out)
  expect_is(fitted.out, "mcmc")
})

# Check predictions -------------------
test_that("predict works for spAbund", {
  pred.out <- predict(out, X.0, coords.0, verbose = FALSE)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$mu.0.samples), c(out$n.post * out$n.chains, length(y.0)))
  expect_equal(dim(pred.out$y.0.samples), c(out$n.post * out$n.chains, length(y.0)))
})


# Random intercept and slope ----------------------------------------------
J.x <- 15
J.y <- 15
J <- J.x * J.y
n.rep <- rep(1, J)
beta <- c(-1, 0.2, -0.3)
p.abund <- length(beta)
mu.RE <- list(levels = c(20, 15),
              sigma.sq.mu = c(1, 0.3), 
              beta.indx = list(1, 2))
sp <- TRUE 
cov.model <- 'matern'
family <- 'zi-Gaussian' 
phi <- 3 / .6
nu <- 3
sigma.sq <- 0.4
tau.sq <- 0.2
z <- rbinom(J, 1, 0.5)
dat <- simAbund(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta,
	        mu.RE = mu.RE, sp = sp, phi = phi, tau.sq = tau.sq, nu = nu,
		sigma.sq = sigma.sq, cov.model = cov.model, family = family, z = z)

pred.indx <- sample(1:J, round(J * .25), replace = FALSE)
y <- dat$y[-pred.indx]
y.0 <- dat$y[pred.indx]
X <- dat$X[-pred.indx, , drop = FALSE]
X.re <- dat$X.re[-pred.indx, , drop = FALSE]
coords <- as.matrix(dat$coords[-pred.indx, , drop = FALSE])
X.0 <- dat$X[pred.indx, , drop = FALSE]
X.re.0 <- dat$X.re[pred.indx, , drop = FALSE]
coords.0 <- as.matrix(dat$coords[pred.indx, , drop = FALSE])
z.fit <- z[-pred.indx]
z.0 <- z[pred.indx]

covs <- list(cov.1 = X[, 2], 
	     cov.2 = X[, 3],
	     abund.factor.1 = X.re[, 1], 
             abund.factor.2 = X.re[, 2])
data.list <- list(y = y, covs = covs, coords = coords, z = z.fit)

# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 100), 
                   tau.sq.ig = c(2, 1), 
                   sigma.sq.ig = c(2, 1),
                   phi.unif = c(3 / 1, 3 / .1), 
                   nu.unif = c(0.1, 5)) 
# Starting values
inits.list <- list(beta = 0,
                   sigma.sq.mu = 0.3,
                   tau.sq.ig = 1, phi = 3 / .5, sigma.sq = 0.5, nu = 2)

# Small
n.batch <- 10
batch.length <- 25
n.burn <- 50 
n.thin <- 2
n.chains <- 1
formula <- ~ (1 | abund.factor.1) + (cov.1 | abund.factor.2) + cov.1 + cov.2

out <- spAbund(formula = formula,
	       data = data.list, 
	       n.batch = n.batch, 
	       batch.length = batch.length, 
	       inits = inits.list, 
	       family = 'zi-Gaussian',
	       priors = prior.list, 
	       cov.model = 'matern',
	       NNGP = TRUE,
	       n.neighbors = 5,
	       accept.rate = 0.43, 
	       n.omp.threads = 1, 
	       verbose = FALSE, 
	       n.report = 10,
	       n.burn = n.burn,
	       n.thin = n.thin,
	       n.chains = 1) 

# Test to make sure it worked ---------
test_that("out is of class spAbund", {
  expect_s3_class(out, "spAbund")
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$muRE, TRUE)
})

# Check RE error ----------------------
test_that("random effect gives error when non-numeric", {
  data.list$covs$abund.factor.1 <- factor(data.list$covs$abund.factor.1)
  expect_error(out <- spAbund(formula = formula,
	       data = data.list, 
	       n.batch = n.batch, 
	       batch.length = batch.length, 
	       inits = inits.list, 
	       family = 'Gaussian',
	       tuning = tuning,
	       priors = prior.list, 
	       accept.rate = 0.43, 
	       n.omp.threads = 1, 
	       verbose = FALSE, 
	       n.report = 10,
	       n.burn = n.burn,
	       n.thin = n.thin,
	       n.chains = 1)) 
  data.list$covs$abund.factor.1 <- as.character(factor(data.list$covs$abund.factor.1))
  expect_error(out <- spAbund(formula = formula,
	       data = data.list, 
	       n.batch = n.batch, 
	       batch.length = batch.length, 
	       inits = inits.list, 
	       family = 'Gaussian',
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

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, y)
})

# Check zi-Gaussian -----------------------
test_that("Gaussian works", {
  out <- spAbund(formula = formula,
  	         data = data.list, 
  	         n.batch = n.batch, 
  	         batch.length = batch.length, 
  	         inits = inits.list, 
  	         family = 'Gaussian',
  	         priors = prior.list, 
	         NNGP = TRUE,
	         n.neighbors = 5,
		 cov.model = 'spherical',
  	         accept.rate = 0.43, 
  	         n.omp.threads = 1, 
  	         verbose = FALSE, 
  	         n.report = 10,
  	         n.burn = n.burn,
  	         n.thin = n.thin,
  	         n.chains = 1) 
  expect_s3_class(out, "spAbund")
})

# Check default priors ----------------
test_that("default priors, multiple chains, inits, burn, thin work", {
  out <- spAbund(formula = formula,
  	     data = data.list, 
  	     n.batch = n.batch, 
  	     batch.length = batch.length, 
  	     family = 'zi-Gaussian',
  	     accept.rate = 0.43, 
  	     n.omp.threads = 1, 
	     NNGP = TRUE,
	     n.neighbors = 5,
	     cov.model = 'exponential',
  	     verbose = FALSE, 
  	     n.report = 10,
  	     n.burn = n.burn,
  	     n.thin = n.thin,
  	     n.chains = 2) 
  expect_s3_class(out, "spAbund")
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {
  expect_output(spAbund(formula = formula,
  	     data = data.list, 
  	     n.batch = n.batch, 
  	     batch.length = batch.length, 
  	     family = 'Gaussian',
	     NNGP = TRUE,
	     n.neighbors = 5,
	     cov.model = 'gaussian',
  	     accept.rate = 0.43, 
  	     n.omp.threads = 1, 
  	     verbose = TRUE, 
  	     n.report = 10,
  	     n.burn = n.burn,
  	     n.thin = n.thin,
  	     n.chains = 1))
})

# Check waicAbund -----------------------
test_that("waicAbund works for spAbund", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicAbund(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for spAbund", {
  fitted.out <- fitted(out)
  expect_is(fitted.out, "mcmc")
})

# Check predictions -------------------
test_that("predict works for spAbund", {
  X.0.full <- cbind(X.0, X.re.0)
  colnames(X.0.full) <- c('(Intercept)', 'cov.1', 'cov.2', 
			  'abund.factor.1', 'abund.factor.2')
  z.0.samples <- matrix(z.0, out$n.post * out$n.chains, length(y.0), byrow = TRUE)
  pred.out <- predict(out, X.0.full, coords.0, z.0.samples = z.0.samples, verbose = FALSE)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$mu.0.samples), c(out$n.post * out$n.chains, length(y.0)))
  expect_equal(dim(pred.out$y.0.samples), c(out$n.post * out$n.chains, length(y.0)))
})
