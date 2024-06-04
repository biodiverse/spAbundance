#define USE_FC_LEN_T
#include <string>
#include "util.h"
#include "rpg.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#define R_NO_REMAP
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Linpack.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#ifndef FCONE
# define FCONE
#endif

extern "C" {
  SEXP abundGaussian(SEXP y_r, SEXP X_r, SEXP XRE_r, SEXP XRandom_r,
               SEXP consts_r, SEXP nRELong_r,
               SEXP betaStarting_r, SEXP tauSqStarting_r, SEXP sigmaSqMuStarting_r,
               SEXP betaStarStarting_r,
               SEXP betaStarIndx_r, SEXP betaLevelIndx_r,
               SEXP muBeta_r, SEXP SigmaBeta_r,
               SEXP tauSqA_r, SEXP tauSqB_r,
               SEXP sigmaSqMuA_r, SEXP sigmaSqMuB_r, SEXP nBatch_r,
               SEXP batchLength_r, SEXP acceptRate_r, SEXP nThreads_r, SEXP verbose_r,
               SEXP nReport_r, SEXP samplesInfo_r, SEXP chainInfo_r){

    /**********************************************************************
     * Initial constants
     * *******************************************************************/
    int i, j, l, ll, s, r, q, info, nProtect=0;
    int status = 0; // For AMCMC.
    const int inc = 1;
    const double one = 1.0;
    const double zero = 0.0;
    char const *lower = "L";
    char const *ntran = "N";
    char const *ytran = "T";

    /**********************************************************************
     * Get Inputs
     * *******************************************************************/
    double *y = REAL(y_r);
    double *X = REAL(X_r);
    int *XRE = INTEGER(XRE_r);
    double *XRandom = REAL(XRandom_r);
    // Load constants
    int J = INTEGER(consts_r)[0];
    int p = INTEGER(consts_r)[1];
    int pRE = INTEGER(consts_r)[2];
    int nRE = INTEGER(consts_r)[3];
    int JZero = INTEGER(consts_r)[4];
    int saveFitted = INTEGER(consts_r)[5];
    int pp = p * p;
    int JpRE = J * pRE;
    // Priors
    double *muBeta = (double *) R_alloc(p, sizeof(double));
    F77_NAME(dcopy)(&p, REAL(muBeta_r), &inc, muBeta, &inc);
    double *SigmaBetaInv = (double *) R_alloc(pp, sizeof(double));
    F77_NAME(dcopy)(&pp, REAL(SigmaBeta_r), &inc, SigmaBetaInv, &inc);
    double tauSqA = REAL(tauSqA_r)[0];
    double tauSqB = REAL(tauSqB_r)[0];
    double *sigmaSqMuA = REAL(sigmaSqMuA_r);
    double *sigmaSqMuB = REAL(sigmaSqMuB_r);
    int *nRELong = INTEGER(nRELong_r);
    int *betaStarIndx = INTEGER(betaStarIndx_r);
    int *betaLevelIndx = INTEGER(betaLevelIndx_r);
    int nBatch = INTEGER(nBatch_r)[0];
    int batchLength = INTEGER(batchLength_r)[0];
    int nSamples = nBatch * batchLength;
    int nBurn = INTEGER(samplesInfo_r)[0];
    int nThin = INTEGER(samplesInfo_r)[1];
    int nPost = INTEGER(samplesInfo_r)[2];
    int currChain = INTEGER(chainInfo_r)[0];
    int nChain = INTEGER(chainInfo_r)[1];
    int nThreads = INTEGER(nThreads_r)[0];
    int verbose = INTEGER(verbose_r)[0];
    int nReport = INTEGER(nReport_r)[0];
    int thinIndx = 0;
    int sPost = 0;

#ifdef _OPENMP
    omp_set_num_threads(nThreads);
#else
    if(nThreads > 1){
      Rf_warning("n.omp.threads > %i, but source not compiled with OpenMP support.", nThreads);
      nThreads = 1;
    }
#endif

    /**********************************************************************
     * Print Information
     * *******************************************************************/
    if(verbose){
      if (currChain == 1) {
        Rprintf("----------------------------------------\n");
        Rprintf("\tModel description\n");
        Rprintf("----------------------------------------\n");
	if (JZero > 0) {
          Rprintf("Zero-inflated Gaussian model with %i non-zero sites.\n\n", J);
	} else {
          Rprintf("Gaussian model with %i sites.\n\n", J);
	}
        Rprintf("Samples per chain: %i (%i batches of length %i)\n", nSamples, nBatch, batchLength);
        Rprintf("Burn-in: %i \n", nBurn);
        Rprintf("Thinning Rate: %i \n", nThin);
        Rprintf("Number of Chains: %i \n", nChain);
        Rprintf("Total Posterior Samples: %i \n\n", nPost * nChain);
#ifdef _OPENMP
        Rprintf("Source compiled with OpenMP support and model fit using %i thread(s).\n\n", nThreads);
#else
        Rprintf("Source not compiled with OpenMP support.\n\n");
#endif
      }
      Rprintf("----------------------------------------\n");
      Rprintf("\tChain %i\n", currChain);
      Rprintf("----------------------------------------\n");
      Rprintf("Sampling ... \n");
      #ifdef Win32
        R_FlushConsole();
      #endif
    }

    /**********************************************************************
     * Parameters
     * *******************************************************************/
    double *beta = (double *) R_alloc(p, sizeof(double));
    F77_NAME(dcopy)(&p, REAL(betaStarting_r), &inc, beta, &inc);
    // Occupancy random effect variances
    double *sigmaSqMu = (double *) R_alloc(pRE, sizeof(double));
    F77_NAME(dcopy)(&pRE, REAL(sigmaSqMuStarting_r), &inc, sigmaSqMu, &inc);
    // Latent occupancy random effects
    double *betaStar = (double *) R_alloc(nRE, sizeof(double));
    F77_NAME(dcopy)(&nRE, REAL(betaStarStarting_r), &inc, betaStar, &inc);
    // Nugget
    double tauSq = REAL(tauSqStarting_r)[0];

    /**********************************************************************
     * Return Stuff
     * *******************************************************************/
    SEXP betaSamples_r;
    PROTECT(betaSamples_r = Rf_allocMatrix(REALSXP, p, nPost)); nProtect++;
    zeros(REAL(betaSamples_r), p * nPost);
    SEXP yRepSamples_r;
    SEXP yRepZeroSamples_r;
    SEXP muSamples_r;
    SEXP likeSamples_r;
    if (saveFitted == 1) {
      PROTECT(yRepSamples_r = Rf_allocMatrix(REALSXP, J, nPost)); nProtect++;
      zeros(REAL(yRepSamples_r), J * nPost);
      PROTECT(yRepZeroSamples_r = Rf_allocMatrix(REALSXP, JZero, nPost)); nProtect++;
      zeros(REAL(yRepZeroSamples_r), JZero * nPost);
      PROTECT(muSamples_r = Rf_allocMatrix(REALSXP, J, nPost)); nProtect++;
      zeros(REAL(muSamples_r), J * nPost);
      PROTECT(likeSamples_r = Rf_allocMatrix(REALSXP, J, nPost)); nProtect++;
      zeros(REAL(likeSamples_r), J * nPost);
    }
    // Occurrence random effects
    SEXP sigmaSqMuSamples_r;
    SEXP betaStarSamples_r;
    if (pRE > 0) {
      PROTECT(sigmaSqMuSamples_r = Rf_allocMatrix(REALSXP, pRE, nPost)); nProtect++;
      zeros(REAL(sigmaSqMuSamples_r), pRE * nPost);
      PROTECT(betaStarSamples_r = Rf_allocMatrix(REALSXP, nRE, nPost)); nProtect++;
      zeros(REAL(betaStarSamples_r), nRE * nPost);
    }
    SEXP tauSqSamples_r;
    PROTECT(tauSqSamples_r = Rf_allocMatrix(REALSXP, inc, nPost)); nProtect++;
    zeros(REAL(tauSqSamples_r), nPost);
    // Likelihood samples for WAIC.

    /**********************************************************************
     * Other initial starting stuff
     * *******************************************************************/
    int Jp = J * p;
    double tmp_0, tmp_02;
    double *tmp_pp = (double *) R_alloc(pp, sizeof(double));
    double *tmp_p = (double *) R_alloc(p, sizeof(double));
    double *tmp_p2 = (double *) R_alloc(p, sizeof(double));
    double *tmp_one = (double *) R_alloc(1, sizeof(double));
    int *tmp_J = (int *) R_alloc(J, sizeof(int));
    for (j = 0; j < J; j++) {
      tmp_J[j] = zero;
    }
    double *tmp_Jp = (double *) R_alloc(Jp, sizeof(double));
    double *tmp_J1 = (double *) R_alloc(J, sizeof(double));

    // For latent occupancy
    double *mu = (double *) R_alloc(J, sizeof(double));
    zeros(mu, J);
    double *like = (double *) R_alloc(J, sizeof(double)); zeros(like, J);
    double *yRep = (double *) R_alloc(J, sizeof(double)); zeros(yRep, J);
    double *yRepZero = (double *) R_alloc(JZero, sizeof(double)); zeros(yRepZero, JZero);

    // For normal priors
    // Occupancy regression coefficient priors.
    F77_NAME(dpotrf)(lower, &p, SigmaBetaInv, &p, &info FCONE);
    if(info != 0){Rf_error("c++ error: dpotrf SigmaBetaInv failed\n");}
    F77_NAME(dpotri)(lower, &p, SigmaBetaInv, &p, &info FCONE);
    if(info != 0){Rf_error("c++ error: dpotri SigmaBetaInv failed\n");}
    double *SigmaBetaInvMuBeta = (double *) R_alloc(p, sizeof(double));
    F77_NAME(dsymv)(lower, &p, &one, SigmaBetaInv, &p, muBeta, &inc, &zero,
        	    SigmaBetaInvMuBeta, &inc FCONE);

    /**********************************************************************
     * Prep for random effects
     * *******************************************************************/
    // Site-level sums of the occurrence random effects
    double *betaStarSites = (double *) R_alloc(J, sizeof(double));
    int *betaStarLongIndx = (int *) R_alloc(JpRE, sizeof(int));
    zeros(betaStarSites, J);
    // Initial sums
    for (j = 0; j < J; j++) {
      for (l = 0; l < pRE; l++) {
        betaStarLongIndx[l * J + j] = which(XRE[l * J + j], betaLevelIndx, nRE);
        betaStarSites[j] += betaStar[betaStarLongIndx[l * J + j]] * XRandom[l * J + j];
      }
    }
    // Starting index for occurrence random effects
    int *betaStarStart = (int *) R_alloc(pRE, sizeof(int));
    for (l = 0; l < pRE; l++) {
      betaStarStart[l] = which(l, betaStarIndx, nRE);
    }

    GetRNGstate();

    /**********************************************************************
     * Begin Sampler
     * *******************************************************************/
    for (s = 0, q = 0; s < nBatch; s++) {
      for (r = 0; r < batchLength; r++, q++) {

        /********************************************************************
         *Update Regression Coefficients
         *******************************************************************/
        for (j = 0; j < J; j++) {
          tmp_J1[j] = (y[j] - betaStarSites[j]) / tauSq;
        } // j
        /********************************
         * Compute b.beta
         *******************************/
        F77_NAME(dgemv)(ytran, &J, &p, &one, X, &J, tmp_J1, &inc, &zero, tmp_p, &inc FCONE);
        for (j = 0; j < p; j++) {
          tmp_p[j] += SigmaBetaInvMuBeta[j];
        } // j

        /********************************
         * Compute A.beta
         * *****************************/
        for(j = 0; j < J; j++){
          for(i = 0; i < p; i++){
            tmp_Jp[i*J+j] = X[i*J+j] / tauSq;
          }
        }

        F77_NAME(dgemm)(ytran, ntran, &p, &p, &J, &one, X, &J, tmp_Jp, &J, &zero, tmp_pp, &p FCONE FCONE);
        for (j = 0; j < pp; j++) {
          tmp_pp[j] += SigmaBetaInv[j];
        } // j

        F77_NAME(dpotrf)(lower, &p, tmp_pp, &p, &info FCONE);
        if(info != 0){Rf_error("c++ error: dpotrf here failed\n");}
        F77_NAME(dpotri)(lower, &p, tmp_pp, &p, &info FCONE);
        if(info != 0){Rf_error("c++ error: dpotri here failed\n");}
        F77_NAME(dsymv)(lower, &p, &one, tmp_pp, &p, tmp_p, &inc, &zero, tmp_p2, &inc FCONE);
        F77_NAME(dpotrf)(lower, &p, tmp_pp, &p, &info FCONE);
	if(info != 0){Rf_error("c++ error: dpotrf here failed\n");}
        mvrnorm(beta, tmp_p2, tmp_pp, p);

        /********************************************************************
         *Update random effects variance
         *******************************************************************/
        for (l = 0; l < pRE; l++) {
          tmp_0 = F77_NAME(ddot)(&nRELong[l], &betaStar[betaStarStart[l]], &inc, &betaStar[betaStarStart[l]], &inc);
          tmp_0 *= 0.5;
          sigmaSqMu[l] = rigamma(sigmaSqMuA[l] + nRELong[l] / 2.0, sigmaSqMuB[l] + tmp_0);
        }

        /********************************************************************
         *Update random effects
         *******************************************************************/
        if (pRE > 0) {
          // Update each individual random effect one by one.
          for (l = 0; l < nRE; l++) {
            /********************************
             * Compute b.beta.star
             *******************************/
            zeros(tmp_one, inc);
            tmp_0 = 0.0;
            // Only allow information to come from when XRE == betaLevelIndx[l].
            // aka information only comes from the sites with any given level
            // of a random effect.
            for (j = 0; j < J; j++) {
              if (XRE[betaStarIndx[l] * J + j] == betaLevelIndx[l]) {
                tmp_02 = 0.0;
                for (ll = 0; ll < pRE; ll++) {
                  tmp_02 += betaStar[betaStarLongIndx[ll * J + j]] * XRandom[ll * J + j];
                }
                tmp_one[0] += XRandom[betaStarIndx[l] * J + j] * (y[j] - F77_NAME(ddot)(&p, &X[j], &J, beta, &inc) -
          		    tmp_02 + (betaStar[l] * XRandom[betaStarIndx[l] * J + j])) / tauSq;
                tmp_0 += XRandom[betaStarIndx[l] * J + j] * XRandom[betaStarIndx[l] * J + j] / tauSq;
              }
            }
            /********************************
             * Compute A.beta.star
             *******************************/
            tmp_0 += 1.0 / sigmaSqMu[betaStarIndx[l]];
            tmp_0 = 1.0 / tmp_0;
            betaStar[l] = rnorm(tmp_0 * tmp_one[0], sqrt(tmp_0));
          }

          // Update the RE sums for the current species
          zeros(betaStarSites, J);
          for (j = 0; j < J; j++) {
            for (l = 0; l < pRE; l++) {
              betaStarSites[j] += betaStar[betaStarLongIndx[l * J + j]] * XRandom[l * J + j];
            }
          }
        }

        /********************************************************************
         *Update tau.sq
         *******************************************************************/
        for (j = 0; j < J; j++) {
	  tmp_J1[j] = y[j] -
	      	F77_NAME(ddot)(&p, &X[j], &J, beta, &inc) - betaStarSites[j];
        }
        tauSq = rigamma(tauSqA + J / 2.0, tauSqB + 0.5 *
			   F77_NAME(ddot)(&J, tmp_J1, &inc, tmp_J1, &inc));

        /********************************************************************
         *Get fitted values and likelihood for WAIC
         *******************************************************************/
	if (saveFitted == 1) {
         for (j = 0; j < J; j++) {
           mu[j] = F77_NAME(ddot)(&p, &X[j], &J, beta, &inc) + betaStarSites[j];
           yRep[j] = rnorm(mu[j], sqrt(tauSq));
           like[j] = dnorm(y[j], mu[j], sqrt(tauSq), 0);
	 } // j
	}
        /********************************************************************
         *Get fitted values and likelihood for WAIC for the zero values
         *******************************************************************/
	if (saveFitted == 1) {
          for (j = 0; j < JZero; j++) {
            yRepZero[j] = rnorm(0.0, sqrt(0.0001));
	  } // j
        }
        /********************************************************************
         *Save samples
         *******************************************************************/
	if (q >= nBurn) {
          thinIndx++;
	  if (thinIndx == nThin) {
            F77_NAME(dcopy)(&p, beta, &inc, &REAL(betaSamples_r)[sPost*p], &inc);
	    REAL(tauSqSamples_r)[sPost] = tauSq;
	    if (saveFitted == 1) {
	      F77_NAME(dcopy)(&J, yRep, &inc, &REAL(yRepSamples_r)[sPost*J], &inc);
              F77_NAME(dcopy)(&J, mu, &inc, &REAL(muSamples_r)[sPost*J], &inc);
	      F77_NAME(dcopy)(&JZero, yRepZero, &inc, &REAL(yRepZeroSamples_r)[sPost*JZero], &inc);
              F77_NAME(dcopy)(&J, like, &inc,
        		      &REAL(likeSamples_r)[sPost*J], &inc);
	    }
            if (pRE > 0) {
              F77_NAME(dcopy)(&pRE, sigmaSqMu, &inc,
                  	    &REAL(sigmaSqMuSamples_r)[sPost*pRE], &inc);
              F77_NAME(dcopy)(&nRE, betaStar, &inc,
                  	    &REAL(betaStarSamples_r)[sPost*nRE], &inc);
            }
	    sPost++;
	    thinIndx = 0;
	  }
	}

        R_CheckUserInterrupt();
      } // r (end batch)

      /********************************************************************
       *Report
       *******************************************************************/
      if (verbose) {
	if (status == nReport) {
	  Rprintf("Batch: %i of %i, %3.2f%%\n", s, nBatch, 100.0*s/nBatch);
	  Rprintf("-------------------------------------------------\n");
          #ifdef Win32
	  R_FlushConsole();
          #endif
	  status = 0;
	}
      }
      status++;
    } // s (sample loop)
    if (verbose) {
      Rprintf("Batch: %i of %i, %3.2f%%\n", s, nBatch, 100.0*s/nBatch);
    }


    // This is necessary when generating random numbers in C.
    PutRNGstate();

    //make return object (which is a list)
    SEXP result_r, resultName_r;
    int nResultListObjs = 6;
    if (pRE > 0) {
      nResultListObjs += 2;
    }

    PROTECT(result_r = Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;
    PROTECT(resultName_r = Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;

    // Setting the components of the output list.
    SET_VECTOR_ELT(result_r, 0, betaSamples_r);
    SET_VECTOR_ELT(result_r, 1, tauSqSamples_r);
    if (saveFitted == 1) {
      SET_VECTOR_ELT(result_r, 2, yRepSamples_r);
      SET_VECTOR_ELT(result_r, 3, muSamples_r);
      SET_VECTOR_ELT(result_r, 4, likeSamples_r);
      SET_VECTOR_ELT(result_r, 5, yRepZeroSamples_r);
    }
    if (pRE > 0) {
      SET_VECTOR_ELT(result_r, 6, sigmaSqMuSamples_r);
      SET_VECTOR_ELT(result_r, 7, betaStarSamples_r);
    }

    // Rf_mkChar turns a C string into a CHARSXP
    SET_VECTOR_ELT(resultName_r, 0, Rf_mkChar("beta.samples"));
    SET_VECTOR_ELT(resultName_r, 1, Rf_mkChar("tau.sq.samples"));
    if (saveFitted == 1) {
      SET_VECTOR_ELT(resultName_r, 2, Rf_mkChar("y.rep.samples"));
      SET_VECTOR_ELT(resultName_r, 3, Rf_mkChar("mu.samples"));
      SET_VECTOR_ELT(resultName_r, 4, Rf_mkChar("like.samples"));
      SET_VECTOR_ELT(resultName_r, 5, Rf_mkChar("y.rep.zero.samples"));
    }
    if (pRE > 0) {
      SET_VECTOR_ELT(resultName_r, 6, Rf_mkChar("sigma.sq.mu.samples"));
      SET_VECTOR_ELT(resultName_r, 7, Rf_mkChar("beta.star.samples"));
    }

    // Set the names of the output list.
    Rf_namesgets(result_r, resultName_r);

    //unprotect
    UNPROTECT(nProtect);

    return(result_r);
  }
}


