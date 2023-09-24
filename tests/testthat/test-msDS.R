# Test msDS.R ---------------------------------------------------------------

skip_on_cran()

# Intercept only ----------------------------------------------------------
J.x <- 10
J.y <- 10
J <- J.x * J.y
# Number of distance bins from which to simulate data.
n.bins <- 5
# Length of each bin. This should be of length n.bins
bin.width <- c(.10, .10, .20, .3, .1)
# Number of species
n.sp <- 5
# Community-level abundance coefficients
beta.mean <- c(-1)
p.abund <- length(beta.mean)
tau.sq.beta <- c(0.2)
# Detection coefficients
alpha.mean <- c(-1.0)
p.det <- length(alpha.mean)
tau.sq.alpha <- c(0.1)
# Detection decay function
det.func <- 'halfnormal'
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
sp <- FALSE
family <- 'NB'
kappa <- runif(n.sp, 0.3, 3)
offset <- pi * .8^2
transect <- 'point'
factor.model <- FALSE

dat <- simMsDS(J.x = J.x, J.y = J.y, n.bins = n.bins, bin.width = bin.width,
	     beta = beta, alpha = alpha, det.func = det.func, kappa = kappa,
	     mu.RE = mu.RE, p.RE = p.RE, sp = sp, cov.model = cov.model,
	     sigma.sq = sigma.sq, phi = phi, nu = nu, family = family,
             offset = offset, transect = transect, n.sp = n.sp)

y <- dat$y
X <- dat$X
X.p <- dat$X.p
dist.breaks <- dat$dist.breaks

covs <- cbind(X, X.p)
colnames(covs) <- c('int.abund', 'int.det')

data.list <- list(y = y,
		  covs = covs,
 		  dist.breaks = dist.breaks,
                  offset = offset)

# Priors
prior.list <- list(beta.comm.normal = list(mean = 0, var = 10),
		   alpha.comm.normal = list(mean = 0,
			                    var = 10),
		   tau.sq.beta.ig = list(a = 0.1, b = 0.1),
		   tau.sq.alpha.ig = list(a = 0.1, b = 0.1),
                   kappa.unif = list(0, 100))
# Starting values
inits.list <- list(alpha.comm = 0, alpha = 0,
		   beta.comm = 0, beta = 0, sigma.sq.p = 1, sigma.sq.mu = 1,
		   beta.star = 0, alpha.star = 0,
		   N = apply(y, c(1, 2), sum, na.rm = TRUE),
		   kappa = 1)

tuning <- list(beta = 0.1, alpha = 0.1, beta.star = 0.3, alpha.star = 0.1,
               kappa = 0.2)

n.batch <- 50
batch.length <- 25
n.burn <- 250
n.thin <- 2
n.chains <- 2
abund.formula <- ~ 1
det.formula <- ~ 1

out <- msDS(abund.formula = abund.formula,
	    det.formula = det.formula,
	    data = data.list,
	    n.batch = n.batch,
	    batch.length = batch.length,
	    inits = inits.list,
	    family = 'NB',
	    det.func = 'halfnormal',
	    transect = transect,
	    tuning = tuning,
	    priors = prior.list,
	    accept.rate = 0.43,
	    n.omp.threads = 1,
	    verbose = FALSE,
	    n.report = 100,
	    n.burn = n.burn,
	    n.thin = n.thin,
	    n.chains = n.chains)

# Test to make sure it worked ---------
test_that("out is of class msDS", {
  expect_s3_class(out, "msDS")
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
test_that("Poisson works and negExp works", {
  out <- msDS(abund.formula = abund.formula,
  	    det.formula = det.formula,
  	    data = data.list,
  	    n.batch = n.batch,
  	    batch.length = batch.length,
  	    inits = inits.list,
  	    family = 'Poisson',
	    det.func = 'negexp',
  	    tuning = tuning,
  	    priors = prior.list,
  	    accept.rate = 0.43,
  	    n.omp.threads = 1,
  	    verbose = FALSE,
  	    n.report = 50,
  	    n.burn = n.burn,
  	    n.thin = n.thin,
  	    n.chains = 1)
  expect_s3_class(out, "msDS")
})

# Check default priors ----------------
test_that("default priors, inits, burn, thin work", {
  out <- msDS(abund.formula = ~ 1,
	    det.formula = ~ 1,
	    data = data.list,
	    n.batch = n.batch,
	    batch.length = batch.length,
	    family = 'NB',
	    det.func = 'halfnormal',
	    transect = transect,
	    tuning = tuning,
	    accept.rate = 0.43,
	    n.omp.threads = 1,
	    n.report = 100,
	    verbose = FALSE,
	    n.chains = 1)
  expect_s3_class(out, 'msDS')
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {
  expect_output(msDS(abund.formula = abund.formula,
  	    det.formula = det.formula,
  	    data = data.list,
  	    n.batch = n.batch,
  	    batch.length = batch.length,
  	    inits = inits.list,
  	    family = 'Poisson',
	    det.func = 'halfnormal',
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
test_that("waicAbund works for msDS", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicAbund(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for msDS", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for msDS", {
  X.0 <- dat$X
  pred.out <- predict(out, X.0)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$mu.0.samples), c(out$n.post * out$n.chains, n.sp, J))
  expect_equal(dim(pred.out$N.0.samples), c(out$n.post * out$n.chains, n.sp, J))
})
test_that("detection prediction works", {
  X.p.0 <- dat$X.p
  pred.out <- predict(out, X.p.0, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$sigma.0.samples), c(out$n.post * out$n.chains, n.sp, J))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for msDS", {
  n.post.samples <- out$n.post * out$n.chains
  ppc.out <- ppcAbund(out, 'chi-square', 0)
  expect_type(ppc.out, "list")
  expect_equal(nrow(ppc.out$fit.y), n.post.samples)
  expect_equal(nrow(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, n.sp, J, n.bins))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, n.sp, J, n.bins))

  ppc.out <- ppcAbund(out, 'chi-square', 1)
  expect_type(ppc.out, "list")
  expect_equal(nrow(ppc.out$fit.y), n.post.samples)
  expect_equal(nrow(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, n.sp, J))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, n.sp, J))

  ppc.out <- ppcAbund(out, 'freeman-tukey', 0)
  expect_type(ppc.out, "list")
  expect_equal(nrow(ppc.out$fit.y), n.post.samples)
  expect_equal(nrow(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, n.sp, J, n.bins))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, n.sp, J, n.bins))

  ppc.out <- ppcAbund(out, 'freeman-tukey', 1)
  expect_type(ppc.out, "list")
  expect_equal(nrow(ppc.out$fit.y), n.post.samples)
  expect_equal(nrow(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, n.sp, J))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, n.sp, J))
})

# Abundance covariate only ------------------------------------------------
J.x <- 10
J.y <- 10
J <- J.x * J.y
# Number of distance bins from which to simulate data.
n.bins <- 5
# Length of each bin. This should be of length n.bins
bin.width <- c(.10, .10, .20, .3, .1)
# Number of species
n.sp <- 5
# Community-level abundance coefficients
beta.mean <- c(-1, 0.2)
p.abund <- length(beta.mean)
tau.sq.beta <- c(0.2, 1.0)
# Detection coefficients
alpha.mean <- c(-1.0)
p.det <- length(alpha.mean)
tau.sq.alpha <- c(0.1)
# Detection decay function
det.func <- 'halfnormal'
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
sp <- FALSE
family <- 'NB'
kappa <- runif(n.sp, 0.3, 3)
offset <- pi * .8^2
transect <- 'point'
factor.model <- FALSE

dat <- simMsDS(J.x = J.x, J.y = J.y, n.bins = n.bins, bin.width = bin.width,
	     beta = beta, alpha = alpha, det.func = det.func, kappa = kappa,
	     mu.RE = mu.RE, p.RE = p.RE, sp = sp, cov.model = cov.model,
	     sigma.sq = sigma.sq, phi = phi, nu = nu, family = family,
             offset = offset, transect = transect, n.sp = n.sp)

y <- dat$y
X <- dat$X
X.re <- dat$X.re
X.p <- dat$X.p
X.p.re <- dat$X.p.re
dist.breaks <- dat$dist.breaks

covs <- cbind(X, X.p)
colnames(covs) <- c('int.abund', 'abund.cov.1', 'int.det')

data.list <- list(y = y,
		  covs = covs,
 		  dist.breaks = dist.breaks,
                  offset = offset)

# Priors
prior.list <- list(beta.comm.normal = list(mean = 0, var = 10),
		   alpha.comm.normal = list(mean = 0,
			                    var = 10),
		   tau.sq.beta.ig = list(a = 0.1, b = 0.1),
		   tau.sq.alpha.ig = list(a = 0.1, b = 0.1),
                   kappa.unif = list(0, 100))
# Starting values
inits.list <- list(alpha.comm = 0, alpha = 0,
		   beta.comm = 0, beta = 0, sigma.sq.p = 1, sigma.sq.mu = 1,
		   beta.star = 0, alpha.star = 0,
		   N = apply(y, c(1, 2), sum, na.rm = TRUE),
		   kappa = 1)

tuning <- list(beta = 0.1, alpha = 0.1, beta.star = 0.3, alpha.star = 0.1,
               kappa = 0.2)

n.batch <- 50
batch.length <- 25
n.burn <- 250
n.thin <- 2
n.chains <- 2
abund.formula <- ~ abund.cov.1
det.formula <- ~ 1

out <- msDS(abund.formula = abund.formula,
	    det.formula = det.formula,
	    data = data.list,
	    n.batch = n.batch,
	    batch.length = batch.length,
	    inits = inits.list,
	    family = 'NB',
	    det.func = 'halfnormal',
	    transect = transect,
	    tuning = tuning,
	    priors = prior.list,
	    accept.rate = 0.43,
	    n.omp.threads = 1,
	    verbose = FALSE,
	    n.report = 100,
	    n.burn = n.burn,
	    n.thin = n.thin,
	    n.chains = n.chains)

# Test to make sure it worked ---------
test_that("out is of class msDS", {
  expect_s3_class(out, "msDS")
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
test_that("Poisson works and negExp works", {
  out <- msDS(abund.formula = abund.formula,
  	    det.formula = det.formula,
  	    data = data.list,
  	    n.batch = n.batch,
  	    batch.length = batch.length,
  	    inits = inits.list,
  	    family = 'Poisson',
	    det.func = 'negexp',
  	    tuning = tuning,
  	    priors = prior.list,
  	    accept.rate = 0.43,
  	    n.omp.threads = 1,
  	    verbose = FALSE,
  	    n.report = 50,
  	    n.burn = n.burn,
  	    n.thin = n.thin,
  	    n.chains = 1)
  expect_s3_class(out, "msDS")
})

# Check default priors ----------------
test_that("default priors, inits, burn, thin work", {
  out <- msDS(abund.formula = ~ 1,
	    det.formula = ~ 1,
	    data = data.list,
	    n.batch = n.batch,
	    batch.length = batch.length,
	    family = 'NB',
	    det.func = 'halfnormal',
	    transect = transect,
	    tuning = tuning,
	    accept.rate = 0.43,
	    n.omp.threads = 1,
	    n.report = 100,
	    verbose = FALSE,
	    n.chains = 1)
  expect_s3_class(out, 'msDS')
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {
  expect_output(msDS(abund.formula = abund.formula,
  	    det.formula = det.formula,
  	    data = data.list,
  	    n.batch = n.batch,
  	    batch.length = batch.length,
  	    inits = inits.list,
  	    family = 'Poisson',
	    det.func = 'halfnormal',
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
test_that("waicAbund works for msDS", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicAbund(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for msDS", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for msDS", {
  X.0 <- dat$X
  pred.out <- predict(out, X.0)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$mu.0.samples), c(out$n.post * out$n.chains, n.sp, J))
  expect_equal(dim(pred.out$N.0.samples), c(out$n.post * out$n.chains, n.sp, J))
})
test_that("detection prediction works", {
  X.p.0 <- dat$X.p
  pred.out <- predict(out, X.p.0, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$sigma.0.samples), c(out$n.post * out$n.chains, n.sp, J))
})

# Detection covariate only ------------------------------------------------
J.x <- 10
J.y <- 10
J <- J.x * J.y
# Number of distance bins from which to simulate data.
n.bins <- 5
# Length of each bin. This should be of length n.bins
bin.width <- c(.10, .10, .20, .3, .1)
# Number of species
n.sp <- 5
# Community-level abundance coefficients
beta.mean <- c(-1)
p.abund <- length(beta.mean)
tau.sq.beta <- c(0.2)
# Detection coefficients
alpha.mean <- c(-1.0, 0.2)
p.det <- length(alpha.mean)
tau.sq.alpha <- c(0.1, 0.5)
# Detection decay function
det.func <- 'halfnormal'
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
sp <- FALSE
family <- 'NB'
kappa <- runif(n.sp, 0.3, 3)
offset <- pi * .8^2
transect <- 'point'
factor.model <- FALSE

dat <- simMsDS(J.x = J.x, J.y = J.y, n.bins = n.bins, bin.width = bin.width,
	     beta = beta, alpha = alpha, det.func = det.func, kappa = kappa,
	     mu.RE = mu.RE, p.RE = p.RE, sp = sp, cov.model = cov.model,
	     sigma.sq = sigma.sq, phi = phi, nu = nu, family = family,
             offset = offset, transect = transect, n.sp = n.sp)

y <- dat$y
X <- dat$X
X.re <- dat$X.re
X.p <- dat$X.p
X.p.re <- dat$X.p.re
dist.breaks <- dat$dist.breaks

covs <- cbind(X, X.p)
colnames(covs) <- c('int.abund', 'int.det', 'det.cov.1')

data.list <- list(y = y,
		  covs = covs,
 		  dist.breaks = dist.breaks,
                  offset = offset)

# Priors
prior.list <- list(beta.comm.normal = list(mean = 0, var = 10),
		   alpha.comm.normal = list(mean = 0,
			                    var = 10),
		   tau.sq.beta.ig = list(a = 0.1, b = 0.1),
		   tau.sq.alpha.ig = list(a = 0.1, b = 0.1),
                   kappa.unif = list(0, 100))
# Starting values
inits.list <- list(alpha.comm = 0, alpha = 0,
		   beta.comm = 0, beta = 0, sigma.sq.p = 1, sigma.sq.mu = 1,
		   beta.star = 0, alpha.star = 0,
		   N = apply(y, c(1, 2), sum, na.rm = TRUE),
		   kappa = 1)

tuning <- list(beta = 0.1, alpha = 0.1, beta.star = 0.3, alpha.star = 0.1,
               kappa = 0.2)

n.batch <- 50
batch.length <- 25
n.burn <- 250
n.thin <- 2
n.chains <- 2
abund.formula <- ~ 1 
det.formula <- ~ det.cov.1

out <- msDS(abund.formula = abund.formula,
	    det.formula = det.formula,
	    data = data.list,
	    n.batch = n.batch,
	    batch.length = batch.length,
	    inits = inits.list,
	    family = 'NB',
	    det.func = 'halfnormal',
	    transect = transect,
	    tuning = tuning,
	    priors = prior.list,
	    accept.rate = 0.43,
	    n.omp.threads = 1,
	    verbose = FALSE,
	    n.report = 100,
	    n.burn = n.burn,
	    n.thin = n.thin,
	    n.chains = n.chains)

# Test to make sure it worked ---------
test_that("out is of class msDS", {
  expect_s3_class(out, "msDS")
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
test_that("Poisson works and negExp works", {
  out <- msDS(abund.formula = abund.formula,
  	    det.formula = det.formula,
  	    data = data.list,
  	    n.batch = n.batch,
  	    batch.length = batch.length,
  	    inits = inits.list,
  	    family = 'Poisson',
	    det.func = 'negexp',
  	    tuning = tuning,
  	    priors = prior.list,
  	    accept.rate = 0.43,
  	    n.omp.threads = 1,
  	    verbose = FALSE,
  	    n.report = 50,
  	    n.burn = n.burn,
  	    n.thin = n.thin,
  	    n.chains = 1)
  expect_s3_class(out, "msDS")
})

# Check default priors ----------------
test_that("default priors, inits, burn, thin work", {
  out <- msDS(abund.formula = ~ 1,
	    det.formula = ~ 1,
	    data = data.list,
	    n.batch = n.batch,
	    batch.length = batch.length,
	    family = 'NB',
	    det.func = 'halfnormal',
	    transect = transect,
	    tuning = tuning,
	    accept.rate = 0.43,
	    n.omp.threads = 1,
	    n.report = 100,
	    verbose = FALSE,
	    n.chains = 1)
  expect_s3_class(out, 'msDS')
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {
  expect_output(msDS(abund.formula = abund.formula,
  	    det.formula = det.formula,
  	    data = data.list,
  	    n.batch = n.batch,
  	    batch.length = batch.length,
  	    inits = inits.list,
  	    family = 'Poisson',
	    det.func = 'halfnormal',
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
test_that("waicAbund works for msDS", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicAbund(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for msDS", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for msDS", {
  X.0 <- dat$X
  pred.out <- predict(out, X.0)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$mu.0.samples), c(out$n.post * out$n.chains, n.sp, J))
  expect_equal(dim(pred.out$N.0.samples), c(out$n.post * out$n.chains, n.sp, J))
})
test_that("detection prediction works", {
  X.p.0 <- dat$X.p
  pred.out <- predict(out, X.p.0, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$sigma.0.samples), c(out$n.post * out$n.chains, n.sp, J))
})

# Covariates on both ------------------------------------------------------
J.x <- 10
J.y <- 10
J <- J.x * J.y
# Number of distance bins from which to simulate data.
n.bins <- 5
# Length of each bin. This should be of length n.bins
bin.width <- c(.10, .10, .20, .3, .1)
# Number of species
n.sp <- 5
# Community-level abundance coefficients
beta.mean <- c(-1, 0.2, 0.3)
p.abund <- length(beta.mean)
tau.sq.beta <- c(0.2, 0.4, 0.8)
# Detection coefficients
alpha.mean <- c(-1.0, 0.2)
p.det <- length(alpha.mean)
tau.sq.alpha <- c(0.1, 0.5)
# Detection decay function
det.func <- 'halfnormal'
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
sp <- FALSE
family <- 'NB'
kappa <- runif(n.sp, 0.3, 3)
offset <- pi * .8^2
transect <- 'point'
factor.model <- FALSE

dat <- simMsDS(J.x = J.x, J.y = J.y, n.bins = n.bins, bin.width = bin.width,
	     beta = beta, alpha = alpha, det.func = det.func, kappa = kappa,
	     mu.RE = mu.RE, p.RE = p.RE, sp = sp, cov.model = cov.model,
	     sigma.sq = sigma.sq, phi = phi, nu = nu, family = family,
             offset = offset, transect = transect, n.sp = n.sp)

y <- dat$y
X <- dat$X
X.re <- dat$X.re
X.p <- dat$X.p
X.p.re <- dat$X.p.re
dist.breaks <- dat$dist.breaks

covs <- cbind(X, X.p)
colnames(covs) <- c('int.abund', 'abund.cov.1', 'abund.cov.2', 
		    'int.det', 'det.cov.1')

data.list <- list(y = y,
		  covs = covs,
 		  dist.breaks = dist.breaks,
                  offset = offset)

# Priors
prior.list <- list(beta.comm.normal = list(mean = 0, var = 10),
		   alpha.comm.normal = list(mean = 0,
			                    var = 10),
		   tau.sq.beta.ig = list(a = 0.1, b = 0.1),
		   tau.sq.alpha.ig = list(a = 0.1, b = 0.1),
                   kappa.unif = list(0, 100))
# Starting values
inits.list <- list(alpha.comm = 0, alpha = 0,
		   beta.comm = 0, beta = 0, sigma.sq.p = 1, sigma.sq.mu = 1,
		   beta.star = 0, alpha.star = 0,
		   N = apply(y, c(1, 2), sum, na.rm = TRUE),
		   kappa = 1)

tuning <- list(beta = 0.1, alpha = 0.1, beta.star = 0.3, alpha.star = 0.1,
               kappa = 0.2)

n.batch <- 50
batch.length <- 25
n.burn <- 250
n.thin <- 2
n.chains <- 2
abund.formula <- ~ abund.cov.1 + abund.cov.2 
det.formula <- ~ det.cov.1

out <- msDS(abund.formula = abund.formula,
	    det.formula = det.formula,
	    data = data.list,
	    n.batch = n.batch,
	    batch.length = batch.length,
	    inits = inits.list,
	    family = 'NB',
	    det.func = 'halfnormal',
	    transect = transect,
	    tuning = tuning,
	    priors = prior.list,
	    accept.rate = 0.43,
	    n.omp.threads = 1,
	    verbose = FALSE,
	    n.report = 100,
	    n.burn = n.burn,
	    n.thin = n.thin,
	    n.chains = n.chains)

# Test to make sure it worked ---------
test_that("out is of class msDS", {
  expect_s3_class(out, "msDS")
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
test_that("Poisson works and negExp works", {
  out <- msDS(abund.formula = abund.formula,
  	    det.formula = det.formula,
  	    data = data.list,
  	    n.batch = n.batch,
  	    batch.length = batch.length,
  	    inits = inits.list,
  	    family = 'Poisson',
	    det.func = 'negexp',
  	    tuning = tuning,
  	    priors = prior.list,
  	    accept.rate = 0.43,
  	    n.omp.threads = 1,
  	    verbose = FALSE,
  	    n.report = 50,
  	    n.burn = n.burn,
  	    n.thin = n.thin,
  	    n.chains = 1)
  expect_s3_class(out, "msDS")
})

# Check default priors ----------------
test_that("default priors, inits, burn, thin work", {
  out <- msDS(abund.formula = ~ 1,
	    det.formula = ~ 1,
	    data = data.list,
	    n.batch = n.batch,
	    batch.length = batch.length,
	    family = 'NB',
	    det.func = 'halfnormal',
	    transect = transect,
	    tuning = tuning,
	    accept.rate = 0.43,
	    n.omp.threads = 1,
	    n.report = 100,
	    verbose = FALSE,
	    n.chains = 1)
  expect_s3_class(out, 'msDS')
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {
  expect_output(msDS(abund.formula = abund.formula,
  	    det.formula = det.formula,
  	    data = data.list,
  	    n.batch = n.batch,
  	    batch.length = batch.length,
  	    inits = inits.list,
  	    family = 'Poisson',
	    det.func = 'halfnormal',
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
test_that("waicAbund works for msDS", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicAbund(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for msDS", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for msDS", {
  X.0 <- dat$X
  pred.out <- predict(out, X.0)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$mu.0.samples), c(out$n.post * out$n.chains, n.sp, J))
  expect_equal(dim(pred.out$N.0.samples), c(out$n.post * out$n.chains, n.sp, J))
})
test_that("detection prediction works", {
  X.p.0 <- dat$X.p
  pred.out <- predict(out, X.p.0, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$sigma.0.samples), c(out$n.post * out$n.chains, n.sp, J))
})

# Interactions on both ----------------------------------------------------
J.x <- 10
J.y <- 10
J <- J.x * J.y
# Number of distance bins from which to simulate data.
n.bins <- 5
# Length of each bin. This should be of length n.bins
bin.width <- c(.10, .10, .20, .3, .1)
# Number of species
n.sp <- 5
# Community-level abundance coefficients
beta.mean <- c(-1, 0.2, 0.3)
p.abund <- length(beta.mean)
tau.sq.beta <- c(0.2, 0.5, 0.8)
# Detection coefficients
alpha.mean <- c(-1.0, 0.2, -0.2, 0.3)
p.det <- length(alpha.mean)
tau.sq.alpha <- c(0.1, 0.5, 0.4, 0.8)
# Detection decay function
det.func <- 'halfnormal'
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
sp <- FALSE
family <- 'NB'
kappa <- runif(n.sp, 0.3, 3)
offset <- pi * .8^2
transect <- 'point'
factor.model <- FALSE

dat <- simMsDS(J.x = J.x, J.y = J.y, n.bins = n.bins, bin.width = bin.width,
	     beta = beta, alpha = alpha, det.func = det.func, kappa = kappa,
	     mu.RE = mu.RE, p.RE = p.RE, sp = sp, cov.model = cov.model,
	     sigma.sq = sigma.sq, phi = phi, nu = nu, family = family,
             offset = offset, transect = transect, n.sp = n.sp)

y <- dat$y
X <- dat$X
X.re <- dat$X.re
X.p <- dat$X.p
X.p.re <- dat$X.p.re
dist.breaks <- dat$dist.breaks

covs <- cbind(X, X.p)
colnames(covs) <- c('int.abund', 'abund.cov.1', 'abund.cov.2', 
		    'int.det', 'det.cov.1', 'det.cov.2', 'det.cov.3')

data.list <- list(y = y,
		  covs = covs,
 		  dist.breaks = dist.breaks,
                  offset = offset)

# Priors
prior.list <- list(beta.comm.normal = list(mean = 0, var = 10),
		   alpha.comm.normal = list(mean = 0,
			                    var = 10),
		   tau.sq.beta.ig = list(a = 0.1, b = 0.1),
		   tau.sq.alpha.ig = list(a = 0.1, b = 0.1),
                   kappa.unif = list(0, 100))
# Starting values
inits.list <- list(alpha.comm = 0, alpha = 0,
		   beta.comm = 0, beta = 0, sigma.sq.p = 1, sigma.sq.mu = 1,
		   beta.star = 0, alpha.star = 0,
		   N = apply(y, c(1, 2), sum, na.rm = TRUE),
		   kappa = 1)

tuning <- list(beta = 0.1, alpha = 0.1, beta.star = 0.3, alpha.star = 0.1,
               kappa = 0.2)

n.batch <- 50
batch.length <- 25
n.burn <- 250
n.thin <- 2
n.chains <- 2
abund.formula <- ~ abund.cov.1 * abund.cov.2 
det.formula <- ~ det.cov.1 + det.cov.2 * det.cov.3

out <- msDS(abund.formula = abund.formula,
	    det.formula = det.formula,
	    data = data.list,
	    n.batch = n.batch,
	    batch.length = batch.length,
	    inits = inits.list,
	    family = 'NB',
	    det.func = 'halfnormal',
	    transect = transect,
	    tuning = tuning,
	    priors = prior.list,
	    accept.rate = 0.43,
	    n.omp.threads = 1,
	    verbose = FALSE,
	    n.report = 100,
	    n.burn = n.burn,
	    n.thin = n.thin,
	    n.chains = n.chains)

# Test to make sure it worked ---------
test_that("out is of class msDS", {
  expect_s3_class(out, "msDS")
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
test_that("Poisson works and negExp works", {
  out <- msDS(abund.formula = abund.formula,
  	    det.formula = det.formula,
  	    data = data.list,
  	    n.batch = n.batch,
  	    batch.length = batch.length,
  	    inits = inits.list,
  	    family = 'Poisson',
	    det.func = 'negexp',
  	    tuning = tuning,
  	    priors = prior.list,
  	    accept.rate = 0.43,
  	    n.omp.threads = 1,
  	    verbose = FALSE,
  	    n.report = 50,
  	    n.burn = n.burn,
  	    n.thin = n.thin,
  	    n.chains = 1)
  expect_s3_class(out, "msDS")
})

# Check default priors ----------------
test_that("default priors, inits, burn, thin work", {
  out <- msDS(abund.formula = abund.formula,
	    det.formula = det.formula,
	    data = data.list,
	    n.batch = n.batch,
	    batch.length = batch.length,
	    family = 'NB',
	    det.func = 'halfnormal',
	    transect = transect,
	    tuning = tuning,
	    accept.rate = 0.43,
	    n.omp.threads = 1,
	    n.report = 100,
	    verbose = FALSE,
	    n.chains = 1)
  expect_s3_class(out, 'msDS')
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {
  expect_output(msDS(abund.formula = abund.formula,
  	    det.formula = det.formula,
  	    data = data.list,
  	    n.batch = n.batch,
  	    batch.length = batch.length,
  	    inits = inits.list,
  	    family = 'Poisson',
	    det.func = 'halfnormal',
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
test_that("waicAbund works for msDS", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicAbund(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for msDS", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for msDS", {
  X.0 <- dat$X
  X.0 <- cbind(X.0, X.0[, 2] * X.0[, 3])
  pred.out <- predict(out, X.0)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$mu.0.samples), c(out$n.post * out$n.chains, n.sp, J))
  expect_equal(dim(pred.out$N.0.samples), c(out$n.post * out$n.chains, n.sp, J))
})
test_that("detection prediction works", {
  X.p.0 <- dat$X.p
  X.p.0 <- cbind(X.p.0, X.p.0[, 4] * X.p.0[, 3])
  pred.out <- predict(out, X.p.0, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$sigma.0.samples), c(out$n.post * out$n.chains, n.sp, J))
})

# Random intercept on Abundance -------------------------------------------
J.x <- 10
J.y <- 10
J <- J.x * J.y
# Number of distance bins from which to simulate data.
n.bins <- 5
# Length of each bin. This should be of length n.bins
bin.width <- c(.10, .10, .20, .3, .1)
# Number of species
n.sp <- 5
# Community-level abundance coefficients
beta.mean <- c(-1)
p.abund <- length(beta.mean)
tau.sq.beta <- c(0.2)
# Detection coefficients
alpha.mean <- c(-1.0)
p.det <- length(alpha.mean)
tau.sq.alpha <- c(0.1)
# Detection decay function
det.func <- 'halfnormal'
mu.RE <- list(levels = c(15), 
	      sigma.sq.mu = c(0.4))
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
sp <- FALSE
family <- 'NB'
kappa <- runif(n.sp, 0.3, 3)
offset <- pi * .8^2
transect <- 'point'
factor.model <- FALSE

dat <- simMsDS(J.x = J.x, J.y = J.y, n.bins = n.bins, bin.width = bin.width,
	     beta = beta, alpha = alpha, det.func = det.func, kappa = kappa,
	     mu.RE = mu.RE, p.RE = p.RE, sp = sp, cov.model = cov.model,
	     sigma.sq = sigma.sq, phi = phi, nu = nu, family = family,
             offset = offset, transect = transect, n.sp = n.sp)

y <- dat$y
X <- dat$X
X.re <- dat$X.re
X.p <- dat$X.p
X.p.re <- dat$X.p.re
dist.breaks <- dat$dist.breaks

covs <- cbind(X, X.re, X.p)
colnames(covs) <- c('int.abund', 'abund.factor.1', 'int.det')

data.list <- list(y = y,
		  covs = covs,
 		  dist.breaks = dist.breaks,
                  offset = offset)

# Priors
prior.list <- list(beta.comm.normal = list(mean = 0, var = 10),
		   alpha.comm.normal = list(mean = 0,
			                    var = 10),
		   tau.sq.beta.ig = list(a = 0.1, b = 0.1),
		   sigma.sq.mu.ig = list(a = 0.1, b = 0.1),
		   sigma.sq.p.ig = list(a = 0.1, b = 0.1),
		   tau.sq.alpha.ig = list(a = 0.1, b = 0.1),
                   kappa.unif = list(0, 100))
# Starting values
inits.list <- list(alpha.comm = 0, alpha = 0,
		   beta.comm = 0, beta = 0, sigma.sq.p = 1, sigma.sq.mu = 1,
		   beta.star = 0, alpha.star = 0,
		   N = apply(y, c(1, 2), sum, na.rm = TRUE),
		   kappa = 1)

tuning <- list(beta = 0.1, alpha = 0.1, beta.star = 0.3, alpha.star = 0.1,
               kappa = 0.2)

n.batch <- 50
batch.length <- 25
n.burn <- 250
n.thin <- 2
n.chains <- 2
abund.formula <- ~ (1 | abund.factor.1)
det.formula <- ~ 1

out <- msDS(abund.formula = abund.formula,
	    det.formula = det.formula,
	    data = data.list,
	    n.batch = n.batch,
	    batch.length = batch.length,
	    inits = inits.list,
	    family = 'NB',
	    det.func = 'halfnormal',
	    transect = transect,
	    tuning = tuning,
	    priors = prior.list,
	    accept.rate = 0.43,
	    n.omp.threads = 1,
	    verbose = FALSE,
	    n.report = 100,
	    n.burn = n.burn,
	    n.thin = n.thin,
	    n.chains = n.chains)

# Test to make sure it worked ---------
test_that("out is of class msDS", {
  expect_s3_class(out, "msDS")
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$pRE, FALSE)
  expect_equal(out$muRE, TRUE)
})

test_that("random effect gives error when non-numeric", {
  data.list$covs <- as.data.frame(data.list$covs)
  data.list$covs$abund.factor.1 <- factor(data.list$covs$abund.factor.1)
  expect_error(out <- msDS(abund.formula = abund.formula,
	    det.formula = det.formula,
	    data = data.list,
	    n.batch = n.batch,
	    batch.length = batch.length,
	    inits = inits.list,
	    family = 'NB',
	    tuning = tuning,
	    priors = prior.list,
	    det.func = 'halfnormal',
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
test_that("Poisson works and negExp works", {
  out <- msDS(abund.formula = abund.formula,
  	    det.formula = det.formula,
  	    data = data.list,
  	    n.batch = n.batch,
  	    batch.length = batch.length,
  	    inits = inits.list,
  	    family = 'Poisson',
	    det.func = 'negexp',
  	    tuning = tuning,
  	    priors = prior.list,
  	    accept.rate = 0.43,
  	    n.omp.threads = 1,
  	    verbose = FALSE,
  	    n.report = 50,
  	    n.burn = n.burn,
  	    n.thin = n.thin,
  	    n.chains = 1)
  expect_s3_class(out, "msDS")
})

# Check default priors ----------------
test_that("default priors, inits, burn, thin work", {
  out <- msDS(abund.formula = ~ 1,
	    det.formula = ~ 1,
	    data = data.list,
	    n.batch = n.batch,
	    batch.length = batch.length,
	    family = 'NB',
	    det.func = 'halfnormal',
	    transect = transect,
	    tuning = tuning,
	    accept.rate = 0.43,
	    n.omp.threads = 1,
	    n.report = 100,
	    verbose = FALSE,
	    n.chains = 1)
  expect_s3_class(out, 'msDS')
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {
  expect_output(msDS(abund.formula = abund.formula,
  	    det.formula = det.formula,
  	    data = data.list,
  	    n.batch = n.batch,
  	    batch.length = batch.length,
  	    inits = inits.list,
  	    family = 'Poisson',
	    det.func = 'halfnormal',
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
test_that("waicAbund works for msDS", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicAbund(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for msDS", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for msDS", {
  X.0 <- cbind(dat$X, dat$X.re)
  colnames(X.0) <- c('(Intercept)', 'abund.factor.1')
  pred.out <- predict(out, X.0)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$mu.0.samples), c(out$n.post * out$n.chains, n.sp, J))
  expect_equal(dim(pred.out$N.0.samples), c(out$n.post * out$n.chains, n.sp, J))
})
test_that("detection prediction works", {
  X.p.0 <- dat$X.p
  pred.out <- predict(out, X.p.0, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$sigma.0.samples), c(out$n.post * out$n.chains, n.sp, J))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for msDS", {
  n.post.samples <- out$n.post * out$n.chains
  ppc.out <- ppcAbund(out, 'chi-square', 0)
  expect_type(ppc.out, "list")
  expect_equal(nrow(ppc.out$fit.y), n.post.samples)
  expect_equal(nrow(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, n.sp, J, n.bins))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, n.sp, J, n.bins))

  ppc.out <- ppcAbund(out, 'chi-square', 1)
  expect_type(ppc.out, "list")
  expect_equal(nrow(ppc.out$fit.y), n.post.samples)
  expect_equal(nrow(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, n.sp, J))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, n.sp, J))

  ppc.out <- ppcAbund(out, 'freeman-tukey', 0)
  expect_type(ppc.out, "list")
  expect_equal(nrow(ppc.out$fit.y), n.post.samples)
  expect_equal(nrow(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, n.sp, J, n.bins))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, n.sp, J, n.bins))

  ppc.out <- ppcAbund(out, 'freeman-tukey', 1)
  expect_type(ppc.out, "list")
  expect_equal(nrow(ppc.out$fit.y), n.post.samples)
  expect_equal(nrow(ppc.out$fit.y.rep), n.post.samples)
  expect_equal(dim(ppc.out$fit.y.group.quants), c(5, n.sp, J))
  expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, n.sp, J))
})

# Abundance REs + covariates in all ---------------------------------------
J.x <- 10
J.y <- 10
J <- J.x * J.y
# Number of distance bins from which to simulate data.
n.bins <- 5
# Length of each bin. This should be of length n.bins
bin.width <- c(.10, .10, .20, .3, .1)
# Number of species
n.sp <- 5
# Community-level abundance coefficients
beta.mean <- c(-1, 0.3, 0.2)
p.abund <- length(beta.mean)
tau.sq.beta <- c(0.2, 0.4, 0.8)
# Detection coefficients
alpha.mean <- c(-1.0, -0.2)
p.det <- length(alpha.mean)
tau.sq.alpha <- c(0.1, 0.5)
# Detection decay function
det.func <- 'halfnormal'
mu.RE <- list(levels = c(15, 10),
	      sigma.sq.mu = c(0.4, 0.9))
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
sp <- FALSE
family <- 'NB'
kappa <- runif(n.sp, 0.3, 3)
offset <- pi * .8^2
transect <- 'point'
factor.model <- FALSE

dat <- simMsDS(J.x = J.x, J.y = J.y, n.bins = n.bins, bin.width = bin.width,
	     beta = beta, alpha = alpha, det.func = det.func, kappa = kappa,
	     mu.RE = mu.RE, p.RE = p.RE, sp = sp, cov.model = cov.model,
	     sigma.sq = sigma.sq, phi = phi, nu = nu, family = family,
             offset = offset, transect = transect, n.sp = n.sp)

y <- dat$y
X <- dat$X
X.re <- dat$X.re
X.p <- dat$X.p
X.p.re <- dat$X.p.re
dist.breaks <- dat$dist.breaks

covs <- cbind(X, X.re, X.p)
colnames(covs) <- c('int.abund', 'abund.cov.1', 'abund.cov.2', 
		    'abund.factor.1', 'abund.factor.2', 'int.det', 'det.cov.1')

data.list <- list(y = y,
		  covs = covs,
 		  dist.breaks = dist.breaks,
                  offset = offset)

# Priors
prior.list <- list(beta.comm.normal = list(mean = 0, var = 10),
		   alpha.comm.normal = list(mean = 0,
			                    var = 10),
		   tau.sq.beta.ig = list(a = 0.1, b = 0.1),
		   sigma.sq.mu.ig = list(a = 0.1, b = 0.1),
		   sigma.sq.p.ig = list(a = 0.1, b = 0.1),
		   tau.sq.alpha.ig = list(a = 0.1, b = 0.1),
                   kappa.unif = list(0, 100))
# Starting values
inits.list <- list(alpha.comm = 0, alpha = 0,
		   beta.comm = 0, beta = 0, sigma.sq.p = 1, sigma.sq.mu = 1,
		   beta.star = 0, alpha.star = 0,
		   N = apply(y, c(1, 2), sum, na.rm = TRUE),
		   kappa = 1)

tuning <- list(beta = 0.1, alpha = 0.1, beta.star = 0.3, alpha.star = 0.1,
               kappa = 0.2)

n.batch <- 50
batch.length <- 25
n.burn <- 250
n.thin <- 2
n.chains <- 2
abund.formula <- ~ abund.cov.1 + (1 | abund.factor.1) + (1 | abund.factor.2) + abund.cov.2
det.formula <- ~ det.cov.1

out <- msDS(abund.formula = abund.formula,
	    det.formula = det.formula,
	    data = data.list,
	    n.batch = n.batch,
	    batch.length = batch.length,
	    inits = inits.list,
	    family = 'NB',
	    det.func = 'halfnormal',
	    transect = transect,
	    tuning = tuning,
	    priors = prior.list,
	    accept.rate = 0.43,
	    n.omp.threads = 1,
	    verbose = FALSE,
	    n.report = 100,
	    n.burn = n.burn,
	    n.thin = n.thin,
	    n.chains = n.chains)

# Test to make sure it worked ---------
test_that("out is of class msDS", {
  expect_s3_class(out, "msDS")
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$pRE, FALSE)
  expect_equal(out$muRE, TRUE)
})

test_that("random effect gives error when non-numeric", {
  data.list$covs <- as.data.frame(data.list$covs)
  data.list$covs$abund.factor.1 <- factor(data.list$covs$abund.factor.1)
  expect_error(out <- msDS(abund.formula = abund.formula,
	    det.formula = det.formula,
	    data = data.list,
	    n.batch = n.batch,
	    batch.length = batch.length,
	    inits = inits.list,
	    family = 'NB',
	    tuning = tuning,
	    priors = prior.list,
	    det.func = 'halfnormal',
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
test_that("Poisson works and negExp works", {
  out <- msDS(abund.formula = abund.formula,
  	    det.formula = det.formula,
  	    data = data.list,
  	    n.batch = n.batch,
  	    batch.length = batch.length,
  	    inits = inits.list,
  	    family = 'Poisson',
	    det.func = 'negexp',
  	    tuning = tuning,
  	    priors = prior.list,
  	    accept.rate = 0.43,
  	    n.omp.threads = 1,
  	    verbose = FALSE,
  	    n.report = 50,
  	    n.burn = n.burn,
  	    n.thin = n.thin,
  	    n.chains = 1)
  expect_s3_class(out, "msDS")
})

# Check default priors ----------------
test_that("default priors, inits, burn, thin work", {
  out <- msDS(abund.formula = ~ 1,
	    det.formula = ~ 1,
	    data = data.list,
	    n.batch = n.batch,
	    batch.length = batch.length,
	    family = 'NB',
	    det.func = 'halfnormal',
	    transect = transect,
	    tuning = tuning,
	    accept.rate = 0.43,
	    n.omp.threads = 1,
	    n.report = 100,
	    verbose = FALSE,
	    n.chains = 1)
  expect_s3_class(out, 'msDS')
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {
  expect_output(msDS(abund.formula = abund.formula,
  	    det.formula = det.formula,
  	    data = data.list,
  	    n.batch = n.batch,
  	    batch.length = batch.length,
  	    inits = inits.list,
  	    family = 'Poisson',
	    det.func = 'halfnormal',
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
test_that("waicAbund works for msDS", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicAbund(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for msDS", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for msDS", {
  X.0 <- cbind(dat$X, dat$X.re)
  colnames(X.0) <- c('(Intercept)', 'abund.cov.1', 'abund.cov.2', 
		     'abund.factor.1', 'abund.factor.2')
  pred.out <- predict(out, X.0)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$mu.0.samples), c(out$n.post * out$n.chains, n.sp, J))
  expect_equal(dim(pred.out$N.0.samples), c(out$n.post * out$n.chains, n.sp, J))
})
test_that("detection prediction works", {
  X.p.0 <- dat$X.p
  pred.out <- predict(out, X.p.0, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$sigma.0.samples), c(out$n.post * out$n.chains, n.sp, J))
})

# Multiple random intercepts on detection ---------------------------------
J.x <- 10
J.y <- 10
J <- J.x * J.y
# Number of distance bins from which to simulate data.
n.bins <- 5
# Length of each bin. This should be of length n.bins
bin.width <- c(.10, .10, .20, .3, .1)
# Number of species
n.sp <- 5
# Community-level abundance coefficients
beta.mean <- c(-1)
p.abund <- length(beta.mean)
tau.sq.beta <- c(0.2)
# Detection coefficients
alpha.mean <- c(-1.0)
p.det <- length(alpha.mean)
tau.sq.alpha <- c(0.1)
# Detection decay function
det.func <- 'halfnormal'
mu.RE <- list()
p.RE <- list(levels = c(15, 20), 
             sigma.sq.p = c(0.5, 0.8))
# Draw species-level effects from community means.
beta <- matrix(NA, nrow = n.sp, ncol = p.abund)
alpha <- matrix(NA, nrow = n.sp, ncol = p.det)
for (i in 1:p.abund) {
  beta[, i] <- rnorm(n.sp, beta.mean[i], sqrt(tau.sq.beta[i]))
}
for (i in 1:p.det) {
  alpha[, i] <- rnorm(n.sp, alpha.mean[i], sqrt(tau.sq.alpha[i]))
}
sp <- FALSE
family <- 'NB'
kappa <- runif(n.sp, 0.3, 3)
offset <- pi * .8^2
transect <- 'point'
factor.model <- FALSE

dat <- simMsDS(J.x = J.x, J.y = J.y, n.bins = n.bins, bin.width = bin.width,
	     beta = beta, alpha = alpha, det.func = det.func, kappa = kappa,
	     mu.RE = mu.RE, p.RE = p.RE, sp = sp, cov.model = cov.model,
	     sigma.sq = sigma.sq, phi = phi, nu = nu, family = family,
             offset = offset, transect = transect, n.sp = n.sp)

y <- dat$y
X <- dat$X
X.re <- dat$X.re
X.p <- dat$X.p
X.p.re <- dat$X.p.re
dist.breaks <- dat$dist.breaks

covs <- cbind(X, X.p, X.p.re)
colnames(covs) <- c('int.abund',  
		    'int.det', 'det.factor.1', 'det.factor.2')

data.list <- list(y = y,
		  covs = covs,
 		  dist.breaks = dist.breaks,
                  offset = offset)

# Priors
prior.list <- list(beta.comm.normal = list(mean = 0, var = 10),
		   alpha.comm.normal = list(mean = 0,
			                    var = 10),
		   tau.sq.beta.ig = list(a = 0.1, b = 0.1),
		   sigma.sq.mu.ig = list(a = 0.1, b = 0.1),
		   sigma.sq.p.ig = list(a = 0.1, b = 0.1),
		   tau.sq.alpha.ig = list(a = 0.1, b = 0.1),
                   kappa.unif = list(0, 100))
# Starting values
inits.list <- list(alpha.comm = 0, alpha = 0,
		   beta.comm = 0, beta = 0, sigma.sq.p = 1, sigma.sq.mu = 1,
		   beta.star = 0, alpha.star = 0,
		   N = apply(y, c(1, 2), sum, na.rm = TRUE),
		   kappa = 1)

tuning <- list(beta = 0.1, alpha = 0.1, beta.star = 0.3, alpha.star = 0.1,
               kappa = 0.2)

n.batch <- 50
batch.length <- 25
n.burn <- 250
n.thin <- 2
n.chains <- 2
abund.formula <- ~ 1
det.formula <- ~ (1 | det.factor.1) + (1 | det.factor.2)

out <- msDS(abund.formula = abund.formula,
	    det.formula = det.formula,
	    data = data.list,
	    n.batch = n.batch,
	    batch.length = batch.length,
	    inits = inits.list,
	    family = 'NB',
	    det.func = 'halfnormal',
	    transect = transect,
	    tuning = tuning,
	    priors = prior.list,
	    accept.rate = 0.43,
	    n.omp.threads = 1,
	    verbose = FALSE,
	    n.report = 100,
	    n.burn = n.burn,
	    n.thin = n.thin,
	    n.chains = n.chains)

# Test to make sure it worked ---------
test_that("out is of class msDS", {
  expect_s3_class(out, "msDS")
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$pRE, TRUE)
  expect_equal(out$muRE, FALSE)
})

test_that("random effect gives error when non-numeric", {
  data.list$covs <- as.data.frame(data.list$covs)
  data.list$covs$det.factor.1 <- factor(data.list$covs$det.factor.1)
  expect_error(out <- msDS(abund.formula = abund.formula,
	    det.formula = det.formula,
	    data = data.list,
	    n.batch = n.batch,
	    batch.length = batch.length,
	    inits = inits.list,
	    family = 'NB',
	    tuning = tuning,
	    priors = prior.list,
	    det.func = 'halfnormal',
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
test_that("Poisson works and negExp works", {
  out <- msDS(abund.formula = abund.formula,
  	    det.formula = det.formula,
  	    data = data.list,
  	    n.batch = n.batch,
  	    batch.length = batch.length,
  	    inits = inits.list,
  	    family = 'Poisson',
	    det.func = 'negexp',
  	    tuning = tuning,
  	    priors = prior.list,
  	    accept.rate = 0.43,
  	    n.omp.threads = 1,
  	    verbose = FALSE,
  	    n.report = 50,
  	    n.burn = n.burn,
  	    n.thin = n.thin,
  	    n.chains = 1)
  expect_s3_class(out, "msDS")
})

# Check default priors ----------------
test_that("default priors, inits, burn, thin work", {
  out <- msDS(abund.formula = ~ 1,
	    det.formula = ~ 1,
	    data = data.list,
	    n.batch = n.batch,
	    batch.length = batch.length,
	    family = 'NB',
	    det.func = 'halfnormal',
	    transect = transect,
	    tuning = tuning,
	    accept.rate = 0.43,
	    n.omp.threads = 1,
	    n.report = 100,
	    verbose = FALSE,
	    n.chains = 1)
  expect_s3_class(out, 'msDS')
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {
  expect_output(msDS(abund.formula = abund.formula,
  	    det.formula = det.formula,
  	    data = data.list,
  	    n.batch = n.batch,
  	    batch.length = batch.length,
  	    inits = inits.list,
  	    family = 'Poisson',
	    det.func = 'halfnormal',
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
test_that("waicAbund works for msDS", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicAbund(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for msDS", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for msDS", {
  X.0 <- cbind(dat$X)
  colnames(X.0) <- c('(Intercept)')
  pred.out <- predict(out, X.0)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$mu.0.samples), c(out$n.post * out$n.chains, n.sp, J))
  expect_equal(dim(pred.out$N.0.samples), c(out$n.post * out$n.chains, n.sp, J))
})
test_that("detection prediction works", {
  X.p.0 <- cbind(dat$X.p, dat$X.p.re)
  colnames(X.p.0) <- c('(Intercept)', 'det.factor.1', 'det.factor.2')
  pred.out <- predict(out, X.p.0, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$sigma.0.samples), c(out$n.post * out$n.chains, n.sp, J))
})

# Multiple random intercepts, covariates, and random slopes ---------------
J.x <- 10
J.y <- 10
J <- J.x * J.y
# Number of distance bins from which to simulate data.
n.bins <- 5
# Length of each bin. This should be of length n.bins
bin.width <- c(.10, .10, .20, .3, .1)
# Number of species
n.sp <- 5
# Community-level abundance coefficients
beta.mean <- c(-1, 0.2)
p.abund <- length(beta.mean)
tau.sq.beta <- c(0.2, 0.5)
# Detection coefficients
alpha.mean <- c(-1.0, 0.5, -0.3)
p.det <- length(alpha.mean)
tau.sq.alpha <- c(0.1, 0.4, 0.8)
# Detection decay function
det.func <- 'halfnormal'
mu.RE <- list(levels = c(15, 10),
	      sigma.sq.mu = c(0.4, 0.3), 
              beta.indx = c(1, 2))
p.RE <- list(levels = c(15, 20),
             sigma.sq.p = c(0.5, 0.8), 
             alpha.indx = c(1, 2))
# Draw species-level effects from community means.
beta <- matrix(NA, nrow = n.sp, ncol = p.abund)
alpha <- matrix(NA, nrow = n.sp, ncol = p.det)
for (i in 1:p.abund) {
  beta[, i] <- rnorm(n.sp, beta.mean[i], sqrt(tau.sq.beta[i]))
}
for (i in 1:p.det) {
  alpha[, i] <- rnorm(n.sp, alpha.mean[i], sqrt(tau.sq.alpha[i]))
}
sp <- FALSE
family <- 'NB'
kappa <- runif(n.sp, 0.3, 3)
offset <- pi * .8^2
transect <- 'point'
factor.model <- FALSE

dat <- simMsDS(J.x = J.x, J.y = J.y, n.bins = n.bins, bin.width = bin.width,
	     beta = beta, alpha = alpha, det.func = det.func, kappa = kappa,
	     mu.RE = mu.RE, p.RE = p.RE, sp = sp, cov.model = cov.model,
	     sigma.sq = sigma.sq, phi = phi, nu = nu, family = family,
             offset = offset, transect = transect, n.sp = n.sp)

y <- dat$y
X <- dat$X
X.re <- dat$X.re
X.p <- dat$X.p
X.p.re <- dat$X.p.re
dist.breaks <- dat$dist.breaks

covs <- cbind(X, X.re, X.p, X.p.re)
colnames(covs) <- c('int.abund', 'abund.cov.1', 'abund.factor.1', 'abund.factor.2',
		    'int.det', 'det.cov.1', 'det.cov.2', 'det.factor.1', 'det.factor.2')

data.list <- list(y = y,
		  covs = covs,
 		  dist.breaks = dist.breaks,
                  offset = offset)

# Priors
prior.list <- list(beta.comm.normal = list(mean = 0, var = 10),
		   alpha.comm.normal = list(mean = 0,
			                    var = 10),
		   tau.sq.beta.ig = list(a = 0.1, b = 0.1),
		   sigma.sq.mu.ig = list(a = 0.1, b = 0.1),
		   sigma.sq.p.ig = list(a = 0.1, b = 0.1),
		   tau.sq.alpha.ig = list(a = 0.1, b = 0.1),
                   kappa.unif = list(0, 100))
# Starting values
inits.list <- list(alpha.comm = 0, alpha = 0,
		   beta.comm = 0, beta = 0, sigma.sq.p = 1, sigma.sq.mu = 1,
		   beta.star = 0, alpha.star = 0,
		   N = apply(y, c(1, 2), sum, na.rm = TRUE),
		   kappa = 1)

tuning <- list(beta = 0.1, alpha = 0.1, beta.star = 0.3, alpha.star = 0.1,
               kappa = 0.2)

n.batch <- 50
batch.length <- 25
n.burn <- 250
n.thin <- 2
n.chains <- 2
abund.formula <- ~ abund.cov.1 + (1 | abund.factor.1) + (abund.cov.1 | abund.factor.2)
det.formula <- ~ det.cov.1 + det.cov.2 + (1 | det.factor.1) + (det.cov.1 | det.factor.2)

out <- msDS(abund.formula = abund.formula,
	    det.formula = det.formula,
	    data = data.list,
	    n.batch = n.batch,
	    batch.length = batch.length,
	    inits = inits.list,
	    family = 'NB',
	    det.func = 'halfnormal',
	    transect = transect,
	    tuning = tuning,
	    priors = prior.list,
	    accept.rate = 0.43,
	    n.omp.threads = 1,
	    verbose = FALSE,
	    n.report = 100,
	    n.burn = n.burn,
	    n.thin = n.thin,
	    n.chains = n.chains)

# Test to make sure it worked ---------
test_that("out is of class msDS", {
  expect_s3_class(out, "msDS")
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$pRE, TRUE)
  expect_equal(out$muRE, TRUE)
})

test_that("random effect gives error when non-numeric", {
  data.list$covs <- as.data.frame(data.list$covs)
  data.list$covs$det.factor.1 <- factor(data.list$covs$det.factor.1)
  expect_error(out <- msDS(abund.formula = abund.formula,
	    det.formula = det.formula,
	    data = data.list,
	    n.batch = n.batch,
	    batch.length = batch.length,
	    inits = inits.list,
	    family = 'NB',
	    tuning = tuning,
	    priors = prior.list,
	    det.func = 'halfnormal',
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
test_that("Poisson works and negExp works", {
  out <- msDS(abund.formula = abund.formula,
  	    det.formula = det.formula,
  	    data = data.list,
  	    n.batch = n.batch,
  	    batch.length = batch.length,
  	    inits = inits.list,
  	    family = 'Poisson',
	    det.func = 'negexp',
  	    tuning = tuning,
  	    priors = prior.list,
  	    accept.rate = 0.43,
  	    n.omp.threads = 1,
  	    verbose = FALSE,
  	    n.report = 50,
  	    n.burn = n.burn,
  	    n.thin = n.thin,
  	    n.chains = 1)
  expect_s3_class(out, "msDS")
})

# Check default priors ----------------
test_that("default priors, inits, burn, thin work", {
  out <- msDS(abund.formula = ~ 1,
	    det.formula = ~ 1,
	    data = data.list,
	    n.batch = n.batch,
	    batch.length = batch.length,
	    family = 'NB',
	    det.func = 'halfnormal',
	    transect = transect,
	    tuning = tuning,
	    accept.rate = 0.43,
	    n.omp.threads = 1,
	    n.report = 100,
	    verbose = FALSE,
	    n.chains = 1)
  expect_s3_class(out, 'msDS')
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {
  expect_output(msDS(abund.formula = abund.formula,
  	    det.formula = det.formula,
  	    data = data.list,
  	    n.batch = n.batch,
  	    batch.length = batch.length,
  	    inits = inits.list,
  	    family = 'Poisson',
	    det.func = 'halfnormal',
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
test_that("waicAbund works for msDS", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicAbund(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for msDS", {
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 2)
})

# Check predictions -------------------
test_that("predict works for msDS", {
  X.0 <- cbind(dat$X, dat$X.re)
  colnames(X.0) <- c('(Intercept)', 'abund.cov.1', 'abund.factor.1', 'abund.factor.2')
  pred.out <- predict(out, X.0)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$mu.0.samples), c(out$n.post * out$n.chains, n.sp, J))
  expect_equal(dim(pred.out$N.0.samples), c(out$n.post * out$n.chains, n.sp, J))
})
test_that("detection prediction works", {
  X.p.0 <- cbind(dat$X.p, dat$X.p.re)
  colnames(X.p.0) <- c('(Intercept)', 'det.cov.1' ,'det.cov.2', 'det.factor.1', 'det.factor.2')
  pred.out <- predict(out, X.p.0, type = 'detection')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$sigma.0.samples), c(out$n.post * out$n.chains, n.sp, J))
})

