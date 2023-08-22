# Test abundGaussian.R ------------------------------------------------------------

skip_on_cran()

# Intercept Only ----------------------------------------------------------
set.seed(200)
J.x <- 15
J.y <- 15
J <- J.x * J.y
n.rep <- rep(1, J)
beta <- c(3)
p.abund <- length(beta)
mu.RE <- list()
sp <- FALSE 
family <- 'zi-Gaussian' 
tau.sq <- 1
z <- rbinom(J, 1, 0.5)
dat <- simAbund(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta,
	        mu.RE = mu.RE, sp = sp, phi = phi, tau.sq = tau.sq,
		sigma.sq = sigma.sq, cov.model = cov.model, family = family, 
                z = z)

y <- dat$y
X <- dat$X
coords <- dat$coords

data.list <- list(y = y, z = z) 

# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 1000), 
                   tau.sq.ig = c(2, 1)) 
# Starting values
inits.list <- list(beta = 0,
                   tau.sq = 1)

# Small
n.batch <- 40
batch.length <- 25
n.burn <- 400 
n.thin <- 2
n.chains <- 1
formula <- ~ 1

out <- abund(formula = formula,
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
	     n.chains = 1) 

# Test to make sure it worked ---------
test_that("out is of class abund", {
  expect_s3_class(out, "abund")
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$muRE, FALSE)
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, dat$y)
})

# Check zi-Gaussian -------------------
test_that("Poisson works", {
  out <- abund(formula = formula,
  	     data = data.list, 
  	     n.batch = n.batch, 
  	     batch.length = batch.length, 
  	     inits = inits.list, 
  	     family = 'zi-Gaussian',
  	     tuning = tuning,
  	     priors = prior.list, 
  	     accept.rate = 0.43, 
  	     n.omp.threads = 1, 
  	     verbose = FALSE, 
  	     n.report = 10,
  	     n.burn = n.burn,
  	     n.thin = n.thin,
  	     n.chains = 1) 
  expect_s3_class(out, "abund")
})

# Check default priors ----------------
test_that("default priors, multiple chains, inits, burn, thin work", {
  out <- abund(formula = formula,
  	     data = data.list, 
  	     n.batch = n.batch, 
  	     batch.length = batch.length, 
  	     family = 'Gaussian',
  	     accept.rate = 0.43, 
  	     n.omp.threads = 1, 
  	     verbose = FALSE, 
  	     n.report = 10,
  	     n.burn = n.burn,
  	     n.thin = n.thin,
  	     n.chains = 2) 
  expect_s3_class(out, "abund")
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {
  expect_output(abund(formula = formula,
  	     data = data.list, 
  	     n.batch = n.batch, 
  	     batch.length = batch.length, 
  	     family = 'Gaussian',
  	     accept.rate = 0.43, 
  	     n.omp.threads = 1, 
  	     verbose = TRUE, 
  	     n.report = 10,
  	     n.burn = n.burn,
  	     n.thin = n.thin,
  	     n.chains = 1))
})

# Check waicAbund -----------------------
test_that("waicAbund works for abund", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicAbund(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for abund", {
  fitted.out <- fitted(out)
  expect_is(fitted.out, "mcmc")
})

# Check predictions -------------------
test_that("predict works for abund", {
  X.0 <- as.matrix(dat$X)
  pred.out <- predict(out, X.0)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$mu.0.samples), c(out$n.post * out$n.chains, J))
  expect_equal(dim(pred.out$y.0.samples), c(out$n.post * out$n.chains, J))
})

# Covariates --------------------------------------------------------------
J.x <- 15
J.y <- 15
J <- J.x * J.y
n.rep <- rep(1, J)
beta <- c(10, 0.5, -0.3)
p.abund <- length(beta)
mu.RE <- list()
sp <- FALSE 
family <- 'zi-Gaussian' 
tau.sq <- 3
z <- rbinom(J, 1, 0.5)
dat <- simAbund(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta,
	        mu.RE = mu.RE, sp = sp, phi = phi, tau.sq = tau.sq,
		sigma.sq = sigma.sq, cov.model = cov.model, family = family, 
                z = z)

y <- dat$y
X <- dat$X
coords <- dat$coords

covs <- list(cov.1 = X[, 2], 
	     cov.2 = X[, 3])
data.list <- list(y = y, covs = covs, z = z) 

# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 100), 
                   tau.sq = c(0.01, .01)) 
# Starting values
inits.list <- list(beta = 0,
                   sigma.sq.mu = 0.5,
                   kappa = 1)

# Small
n.batch <- 40
batch.length <- 25
n.burn <- 400 
n.thin <- 2
n.chains <- 1
formula <- ~ cov.1 + cov.2

out <- abund(formula = formula,
	     data = data.list, 
	     n.batch = n.batch, 
	     batch.length = batch.length, 
	     inits = inits.list, 
	     family = 'Gaussian',
	     priors = prior.list, 
	     accept.rate = 0.43, 
	     n.omp.threads = 1, 
	     verbose = FALSE, 
	     n.report = 10,
	     n.burn = n.burn,
	     n.thin = n.thin,
	     n.chains = 2) 

# Test to make sure it worked ---------
test_that("out is of class abund", {
  expect_s3_class(out, "abund")
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$muRE, FALSE)
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, dat$y)
})

# Check ZI Gaussian -----------------------
test_that("ZI Gaussian works", {
  out <- abund(formula = formula,
  	     data = data.list, 
  	     n.batch = n.batch, 
  	     batch.length = batch.length, 
  	     inits = inits.list, 
  	     family = 'zi-Gaussian',
  	     tuning = tuning,
  	     priors = prior.list, 
  	     accept.rate = 0.43, 
  	     n.omp.threads = 1, 
  	     verbose = FALSE, 
  	     n.report = 10,
  	     n.burn = n.burn,
  	     n.thin = n.thin,
  	     n.chains = 1) 
  expect_s3_class(out, "abund")
})

# Check default priors ----------------
test_that("default priors, multiple chains, inits, burn, thin work", {
  out <- abund(formula = formula,
  	     data = data.list, 
  	     n.batch = n.batch, 
  	     batch.length = batch.length, 
  	     family = 'zi-Gaussian',
  	     accept.rate = 0.43, 
  	     n.omp.threads = 1, 
  	     verbose = FALSE, 
  	     n.report = 10,
  	     n.burn = n.burn,
  	     n.thin = n.thin,
  	     n.chains = 2) 
  expect_s3_class(out, "abund")
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {
  expect_output(abund(formula = formula,
  	     data = data.list, 
  	     n.batch = n.batch, 
  	     batch.length = batch.length, 
  	     family = 'Gaussian',
  	     accept.rate = 0.43, 
  	     n.omp.threads = 1, 
  	     verbose = TRUE, 
  	     n.report = 10,
  	     n.burn = n.burn,
  	     n.thin = n.thin,
  	     n.chains = 1))
})

# Check waicAbund ---------------------
test_that("waicAbund works for abund", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicAbund(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for abund", {
  fitted.out <- fitted(out)
  expect_is(fitted.out, "mcmc")
})

# Check predictions -------------------
test_that("predict works for abund", {
  X.0 <- dat$X
  pred.out <- predict(out, X.0)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$mu.0.samples), c(out$n.post * out$n.chains, J))
  expect_equal(dim(pred.out$y.0.samples), c(out$n.post * out$n.chains, J))
})

# Interactions ------------------------------------------------------------
J.x <- 15
J.y <- 15
J <- J.x * J.y
n.rep <- sample(1, J, replace = TRUE)
beta <- c(-1, 0.5, -0.3)
p.abund <- length(beta)
mu.RE <- list()
sp <- FALSE 
family <- 'zi-Gaussian' 
tau.sq <- 1 
z <- rbinom(J, 1, 0.5)
dat <- simAbund(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta,
	        mu.RE = mu.RE, sp = sp, phi = phi, tau.sq = tau.sq,
		sigma.sq = sigma.sq, cov.model = cov.model, family = family, 
                z = z)

y <- dat$y
X <- dat$X
coords <- dat$coords

covs <- list(cov.1 = X[, 2], 
	     cov.2 = X[, 3])
data.list <- list(y = y, covs = covs, z = z) 

# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 1000), 
                   tau.sq.ig = c(2, 1)) 
# Starting values
inits.list <- list(beta = 0,
                   tau.sq = 1)

# Small
n.batch <- 40
batch.length <- 25
n.burn <- 400 
n.thin <- 2
n.chains <- 1
formula <- ~ cov.1 * cov.2

out <- abund(formula = formula,
	     data = data.list, 
	     n.batch = n.batch, 
	     batch.length = batch.length, 
	     inits = inits.list, 
	     family = 'Gaussian',
	     priors = prior.list, 
	     accept.rate = 0.43, 
	     n.omp.threads = 1, 
	     verbose = FALSE, 
	     n.report = 10,
	     n.burn = n.burn,
	     n.thin = n.thin,
	     n.chains = 2) 

# Test to make sure it worked ---------
test_that("out is of class abund", {
  expect_s3_class(out, "abund")
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$muRE, FALSE)
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, dat$y)
})

# Check ZI Gaussian -----------------------
test_that("ZI Gaussian works", {
  out <- abund(formula = formula,
  	     data = data.list, 
  	     n.batch = n.batch, 
  	     batch.length = batch.length, 
  	     inits = inits.list, 
  	     family = 'zi-Gaussian',
  	     tuning = tuning,
  	     priors = prior.list, 
  	     accept.rate = 0.43, 
  	     n.omp.threads = 1, 
  	     verbose = FALSE, 
  	     n.report = 10,
  	     n.burn = n.burn,
  	     n.thin = n.thin,
  	     n.chains = 1) 
  expect_s3_class(out, "abund")
})

# Check default priors ----------------
test_that("default priors, multiple chains, inits, burn, thin work", {
  out <- abund(formula = formula,
  	     data = data.list, 
  	     n.batch = n.batch, 
  	     batch.length = batch.length, 
  	     family = 'Gaussian',
  	     accept.rate = 0.43, 
  	     n.omp.threads = 1, 
  	     verbose = FALSE, 
  	     n.report = 10,
  	     n.burn = n.burn,
  	     n.thin = n.thin,
  	     n.chains = 2) 
  expect_s3_class(out, "abund")
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {
  expect_output(abund(formula = formula,
  	     data = data.list, 
  	     n.batch = n.batch, 
  	     batch.length = batch.length, 
  	     family = 'Gaussian',
  	     accept.rate = 0.43, 
  	     n.omp.threads = 1, 
  	     verbose = TRUE, 
  	     n.report = 10,
  	     n.burn = n.burn,
  	     n.thin = n.thin,
  	     n.chains = 1))
})

# Check waicAbund ---------------------
test_that("waicAbund works for abund", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicAbund(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for abund", {
  fitted.out <- fitted(out)
  expect_is(fitted.out, "mcmc")
})

# Check predictions -------------------
test_that("predict works for abund", {
  X.0 <- cbind(dat$X, dat$X[, 2] * dat$X[, 3])
  pred.out <- predict(out, X.0)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$mu.0.samples), c(out$n.post * out$n.chains, J))
  expect_equal(dim(pred.out$y.0.samples), c(out$n.post * out$n.chains, J))
})

# Random intercept --------------------------------------------------------
J.x <- 15
J.y <- 15
J <- J.x * J.y
n.rep <- rep(1, J)
beta <- c(-1)
p.abund <- length(beta)
mu.RE <- list(levels = c(20),
              sigma.sq.mu = c(1))
sp <- FALSE 
family <- 'zi-Gaussian' 
tau.sq <- 1
z <- rbinom(J, 1, 0.5)
dat <- simAbund(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta,
	        mu.RE = mu.RE, sp = sp, phi = phi, tau.sq = tau.sq,
		sigma.sq = sigma.sq, cov.model = cov.model, family = family, 
                z = z)

y <- dat$y
X <- dat$X
X.re <- dat$X.re
coords <- dat$coords

covs <- list(factor.1 = X.re)
data.list <- list(y = y, covs = covs, z = z) 

# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 100), 
                   tau.sq.ig = c(0.001, .001), 
                   sigma.sq.mu.ig = list(0.1, 0.1)) 
# Starting values
inits.list <- list(beta = 0,
                   sigma.sq.mu = 0.3,
                   tau.sq = 1)

# Small
n.batch <- 40
batch.length <- 25
n.burn <- 400 
n.thin <- 2
n.chains <- 1
formula <- ~ (1 | factor.1)

out <- abund(formula = formula,
	     data = data.list, 
	     n.batch = n.batch, 
	     batch.length = batch.length, 
	     inits = inits.list, 
	     family = 'zi-Gaussian',
	     tuning = tuning,
	     priors = prior.list, 
	     accept.rate = 0.43, 
	     n.omp.threads = 1, 
	     verbose = FALSE, 
	     n.report = 10,
	     n.burn = n.burn,
	     n.thin = n.thin,
	     n.chains = 3) 

# Test to make sure it worked ---------
test_that("out is of class abund", {
  expect_s3_class(out, "abund")
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$muRE, TRUE)
})

# Check RE error ----------------------
test_that("random effect gives error when non-numeric", {
  data.list$covs$factor.1 <- factor(data.list$covs$factor.1)
  expect_error(out <- abund(formula = formula,
	       data = data.list, 
	       n.batch = n.batch, 
	       batch.length = batch.length, 
	       inits = inits.list, 
	       family = 'zi-Gaussian',
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
	       family = 'zi-Gaussian',
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
  expect_equal(out$y, dat$y)
})

# Check Gaussian -----------------------
test_that("Gaussian works", {
  out <- abund(formula = formula,
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
  	     n.chains = 1) 
  expect_s3_class(out, "abund")
})

# Check default priors ----------------
test_that("default priors, multiple chains, inits, burn, thin work", {
  out <- abund(formula = formula,
  	     data = data.list, 
  	     n.batch = n.batch, 
  	     batch.length = batch.length, 
  	     family = 'zi-Gaussian',
  	     accept.rate = 0.43, 
  	     n.omp.threads = 1, 
  	     verbose = FALSE, 
  	     n.report = 10,
  	     n.burn = n.burn,
  	     n.thin = n.thin,
  	     n.chains = 2) 
  expect_s3_class(out, "abund")
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {
  expect_output(abund(formula = formula,
  	     data = data.list, 
  	     n.batch = n.batch, 
  	     batch.length = batch.length, 
  	     family = 'zi-Gaussian',
  	     accept.rate = 0.43, 
  	     n.omp.threads = 1, 
  	     verbose = TRUE, 
  	     n.report = 10,
  	     n.burn = n.burn,
  	     n.thin = n.thin,
  	     n.chains = 1))
})

# Check waicAbund -----------------------
test_that("waicAbund works for abund", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicAbund(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for abund", {
  fitted.out <- fitted(out)
  expect_is(fitted.out, "mcmc")
})

# Check predictions -------------------
test_that("predict works for abund", {
  X.0 <- cbind(X, X.re)
  colnames(X.0) <- c('(Intercept)', 'factor.1')
  z.0.samples <- matrix(z, out$n.post * out$n.chains, J, byrow = TRUE)
  pred.out <- predict(out, X.0, z.0.samples = z.0.samples)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$mu.0.samples), c(out$n.post * out$n.chains, J))
  expect_equal(dim(pred.out$y.0.samples), c(out$n.post * out$n.chains, J))
})

# Random intercept and slope ----------------------------------------------
J.x <- 15
J.y <- 15
J <- J.x * J.y
n.rep <- rep(1, J)
beta <- c(-1, 0.5)
p.abund <- length(beta)
mu.RE <- list(levels = c(20, 15),
              sigma.sq.mu = c(1, 0.5), 
              beta.indx = c(1, 2))
sp <- FALSE 
family <- 'Gaussian' 
tau.sq <- 2 
dat <- simAbund(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta,
	        mu.RE = mu.RE, sp = sp, phi = phi, tau.sq = tau.sq,
		sigma.sq = sigma.sq, cov.model = cov.model, family = family)

y <- dat$y
X <- dat$X
X.re <- dat$X.re
coords <- dat$coords

covs <- list(cov.1 = X[, 2], 
	     factor.1 = X.re[, 1], 
             factor.2 = X.re[, 2])
data.list <- list(y = y, covs = covs) 

# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 100), 
                   tau.sq.ig = c(2, 1), 
                   sigma.sq.mu.ig = list(0.1, 0.1)) 
# Starting values
inits.list <- list(beta = 0,
                   sigma.sq.mu = 0.3,
                   tau.sq = 1)

# Small
n.batch <- 40
batch.length <- 25
n.burn <- 400 
n.thin <- 2
n.chains <- 1
formula <- ~ cov.1 + (1 | factor.1) + (cov.1 | factor.2)

out <- abund(formula = formula,
	     data = data.list, 
	     n.batch = n.batch, 
	     batch.length = batch.length, 
	     inits = inits.list, 
	     family = 'Gaussian',
	     priors = prior.list, 
	     accept.rate = 0.43, 
	     n.omp.threads = 1, 
	     verbose = FALSE, 
	     n.report = 10,
	     n.burn = n.burn,
	     n.thin = n.thin,
	     n.chains = 1) 

# Test to make sure it worked ---------
test_that("out is of class abund", {
  expect_s3_class(out, "abund")
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$muRE, TRUE)
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, dat$y)
})

# Check default priors ----------------
test_that("default priors, multiple chains, inits, burn, thin work", {
  out <- abund(formula = formula,
  	     data = data.list, 
  	     n.batch = n.batch, 
  	     batch.length = batch.length, 
  	     family = 'Gaussian',
  	     accept.rate = 0.43, 
  	     n.omp.threads = 1, 
  	     verbose = FALSE, 
  	     n.report = 10,
  	     n.burn = n.burn,
  	     n.thin = n.thin,
  	     n.chains = 2) 
  expect_s3_class(out, "abund")
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {
  expect_output(abund(formula = formula,
  	     data = data.list, 
  	     n.batch = n.batch, 
  	     batch.length = batch.length, 
  	     family = 'Gaussian',
  	     accept.rate = 0.43, 
  	     n.omp.threads = 1, 
  	     verbose = TRUE, 
  	     n.report = 10,
  	     n.burn = n.burn,
  	     n.thin = n.thin,
  	     n.chains = 1))
})

# Check waicAbund -----------------------
test_that("waicAbund works for abund", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicAbund(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for abund", {
  fitted.out <- fitted(out)
  expect_is(fitted.out, "mcmc")
})

# Check predictions -------------------
test_that("predict works for abund", {
  X.0 <- cbind(dat$X, dat$X.re)
  colnames(X.0) <- c('(Intercept)', 'cov.1', 'factor.1', 'factor.2')
  pred.out <- predict(out, X.0)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$mu.0.samples), c(out$n.post * out$n.chains, J))
  expect_equal(dim(pred.out$y.0.samples), c(out$n.post * out$n.chains, J))
})
