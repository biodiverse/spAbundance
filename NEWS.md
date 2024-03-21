# spAbundance 0.1.1

+ Updated the `neonDWP` data set after NEON announced an error in some of the bird point count data (more information [here](https://www.neonscience.org/impact/observatory-blog/bird-point-ids-within-grids-were-transposed-resulting-inaccurate-point). Also updated the associated vignette that uses these data.
+ Updated the `msAbund()` function to make it compatible with `spOccupancy::updateMCMC()`, which allows for picking up a model run where it left off instead of having to restart the MCMC from scratch. I am planning to allow this functionality for all model fitting functions in `spAbundance` in a future version.
+ Updated all model fitting functions to fix a potential error that could arise in the calculation of Rhat that could lead to the function failing under certain circumstances. This was particularly the case in multi-species models with large amounts of rare species. 

# spAbundance 0.1.0

+ This is the first release of `spAbundance`
