
<!-- README.md is generated from README.Rmd. Please edit that file -->

# spAbundance

`spAbundance` fits single-species and multi-species spatial abundance
models using Markov Chain Monte Carlo (MCMC). The package is in initial
stages of development and testing of existing functionality is not
complete so use of the package prior to it’s first release on CRAN is at
your own risk. `spAbundance` uses analogous syntax to its “sister
package” [spOccupancy](https://www.jeffdoser.com/files/spoccupancy-web/)
(Doser et al. 2022).

## Installation

`spAbundance` is in the initial stages of active development. The
functionality currently implemented is not yet stable, and should be
used with caution. To download the development version of the package,
you can use `devtools` as follows:

``` r
devtools::install_github("doserjef/spAbundance")
```

Note that because we implement the MCMC in C++, you will need a C++
compiler on your computer to install the package from GitHub. To compile
C++ on Windows, you can install
[`RTools`](https://cran.r-project.org/bin/windows/Rtools/). To compile
C++ on a Mac, you can install `XCode` from the Mac app store.

## Functionality

Below is a running list of currently implemented functionality, along
with each function’s [lifecycle
stage](https://lifecycle.r-lib.org/articles/stages.html). This will be
updated as new functionality is included

| `spAbundance` Function | Description                                                       |
| ---------------------- | ----------------------------------------------------------------- |
| `NMix()`               | *stable* Single-species N-mixture model                           |
| `spNMix()`             | *stable* Single-species spatial N-mixture model                   |
| `DS()`                 | *experimental* Single-species distance sampling model             |
| `abund()`              | *stable* Single-species abundance GLM                             |
| `spAbund()`            | *stable* Single-species spatial abundance GLM                     |
| `msAbund()`            | *stable* Multi-species abundance GLM                              |
| `sfMsAbund()`          | *experimental* Multi-species spatial abundance GLM                |
| `ppcAbund()`           | *stable* Posterior predictive check using Bayesian p-values       |
| `waicAbund()`          | *stable* Calculate Widely Applicable Information Criterion (WAIC) |
| `simNMix()`            | *stable* Simulate single-species repeated count data              |
| `simAbund()`           | *stable* Simulate single-species count data                       |
| `simDS()`              | *experimental* Simulate single-species distance sampling data     |

All model fitting functions allow for Poisson and negative binomial
distributions for the abundance portion of the model.

## Functionality to be included prior to first release on CRAN

  - Single-species and multi-species spatial and non-spatial N-mixture
    models with Poisson and Negative Binomial distributions for latent
    abundance.
  - Single-species and multi-species spatial and non-spatial abundance
    GLMS without imperfect detection with Poisson and Negative Binomial
    likelihoods.
  - Single-species and multi-species spatial and non-spatial distance
    sampling models with Poisson and Negative Binomial distributions for
    latent abundance.
  - Posterior predcitve checks for all model types.
  - Simulation functions for all model types.
  - WAIC and k-fold cross-validation for all model types.

## Citing `spAbundance`

Please cite `spAbundance` as:

Doser, J. (2023). spAbundance: Single-species and multi-species
spatially-explicit abundance models". R package version 0.1.0.
<https://github.com/doserjef/spAbundance>

## References

Doser, J. W., Finley, A. O., Kéry, M., and Zipkin, E. F. (2022a).
spOccupancy: An R package for single-species, multi-species, and
integrated spatial occupancy models. Methods in Ecology and Evolution.
<https://doi.org/10.1111/2041-210X.13897>.
