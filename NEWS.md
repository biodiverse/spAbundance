# spAbundance 0.1.3

+ Added in the `independent.betas` tag into the `priors` list for certain multi-species model types to allow for specifying an independent prior on the species-specific random effects as opposed to treating species-specific effects as random effects. This can be useful under certain circumstances when the distribution of effects across species may not be adequately represented by a Gaussian distribution. This tag is available for the following functions: `lfMsAbund` (Gaussian only), `sfMsAbund` (Gaussian only), and `svcMsAbund`.
+ Fixed a bug in zero-inflated Gaussian latent factor abundance models (`lfMsAbund`). 
+ Fixed a bug in `waicAbund()` that prevented it from working with `svcMsAbund` models.
+ Fixed minor issues in C++ to pass CRAN additional checks.

# spAbundance 0.1.1

+ Updated the `neonDWP` data set after NEON announced an error in some of the bird point count data (more information [here](https://www.neonscience.org/impact/observatory-blog/bird-point-ids-within-grids-were-transposed-resulting-inaccurate-point). Also updated the associated vignette that uses these data.
+ Updated the `msAbund()` function to make it compatible with `spOccupancy::updateMCMC()`, which allows for picking up a model run where it left off instead of having to restart the MCMC from scratch. I am planning to allow this functionality for all model fitting functions in `spAbundance` in a future version.
+ Updated all model fitting functions to fix a potential error that could arise in the calculation of Rhat that could lead to the function failing under certain circumstances. This was particularly the case in multi-species models with large amounts of rare species. 

# spAbundance 0.1.0

+ This is the first release of `spAbundance`
