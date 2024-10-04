# spAbundance 0.2.1

+ Fixed a `C++` memory issue in `predict.spAbund()` that could result in crashes under certain situations.

# spAbundance 0.2.0

+ Updated `svcAbund()` to now work with Poisson and negative binomial families. Note that the function now defaults to use `family = 'Poisson'`, which differs from the previous implementation when only `family = 'Gaussian'` was supported. This switch was done to maintain consistency with other spAbundance model-fitting functions. 
+ Added functionality for independent priors on the species-specific effects to allow species-effects to be treated as fixed effects as opposed to random effects from a common community-level distribution for the following model types: `msNMix()`.  
+ Added in a check at the top of all model fitting functions to return an error when the number of posterior samples saved based on the MCMC criteria (`n.batch`, `batch.length`, `n.samples`, `n.burn`, `n.thin`, `n.chains`) are specified in a way that leads to a non-integer value. In such situations, models would previously run and return without an error, but sometimes the last posterior sample in any given chain could have widely inaccurate values, or values that prevented subsequent functions from working. Thanks to Wendy Leuenberger for bringing this to my attention. 
+ Fixed some typos in the documentation.
+ Updated C++ code to adhere to the new lack of re-mapping of functions in Rinternals.h and R_ext/Error.h when building packages on CRAN. 
+ Fixed a typo in the generation of initial values for latent unstructured random effects in all model functions. The typo had no major ramifications, if anything it just led to slower convergence, as it resulted in very large (or very small) initial values for the latent random effects that are not really viable on the log scale.
+ Moved GitHub development page to the "biodiverse" group. The GitHub link for the development repository is now [https://github.com/biodiverse/spAbundance](https://github.com/biodiverse/spAbundance). 


# spAbundance 0.1.3

+ Added in the `independent.betas` tag into the `priors` list for certain multi-species model types to allow for specifying an independent prior on the species-specific effects as opposed to treating species-specific effects as random effects. This can be useful under certain circumstances when the distribution of effects across species may not be adequately represented by a Gaussian distribution. This tag is available for the following functions: `lfMsAbund` (Gaussian only), `sfMsAbund` (Gaussian only), and `svcMsAbund`.
+ Fixed a bug in zero-inflated Gaussian latent factor abundance models (`lfMsAbund`). 
+ Fixed a bug in `waicAbund()` that prevented it from working with `svcMsAbund` models.
+ Fixed minor issues in C++ to pass CRAN additional checks.

# spAbundance 0.1.1

+ Updated the `neonDWP` data set after NEON announced an error in some of the bird point count data (more information [here](https://www.neonscience.org/impact/observatory-blog/bird-point-ids-within-grids-were-transposed-resulting-inaccurate-point). Also updated the associated vignette that uses these data.
+ Updated the `msAbund()` function to make it compatible with `spOccupancy::updateMCMC()`, which allows for picking up a model run where it left off instead of having to restart the MCMC from scratch. I am planning to allow this functionality for all model fitting functions in `spAbundance` in a future version.
+ Updated all model fitting functions to fix a potential error that could arise in the calculation of Rhat that could lead to the function failing under certain circumstances. This was particularly the case in multi-species models with large amounts of rare species. 

# spAbundance 0.1.0

+ This is the first release of `spAbundance`
