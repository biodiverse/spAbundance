#define USE_FC_LEN_T
#include <string>
#include "util.h"
#include "rpg.h"

#ifdef _OPENMP
#include <omp.h>
#endif

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
  SEXP abund(SEXP y_r, SEXP X_r, SEXP XRE_r, SEXP XRandom_r,
             SEXP consts_r,SEXP nAbundRELong_r, 
             SEXP betaStarting_r, SEXP kappaStarting_r, 
	     SEXP sigmaSqMuStarting_r, SEXP betaStarStarting_r, 
             SEXP siteIndx_r, SEXP betaStarIndx_r, SEXP betaLevelIndx_r, 
	     SEXP muBeta_r, SEXP SigmaBeta_r, 
	     SEXP kappaA_r, SEXP kappaB_r, SEXP sigmaSqMuA_r, 
	     SEXP sigmaSqMuB_r, SEXP tuning_r,
	     SEXP nBatch_r, SEXP batchLength_r, SEXP acceptRate_r, SEXP nThreads_r, 
             SEXP verbose_r, SEXP nReport_r, SEXP samplesInfo_r,
	     SEXP chainInfo_r, SEXP family_r, SEXP offset_r){
   
    /**********************************************************************
     * Initial constants
     * *******************************************************************/
    int g, t, j, s, l, k, ll, nProtect=0;
    const int inc = 1;
    
    /**********************************************************************
     * Get Inputs
     * *******************************************************************/
    double *y = REAL(y_r);
    double *X = REAL(X_r);
    double *offset = REAL(offset_r);
    int *XRE = INTEGER(XRE_r); 
    // Load constants
    int J = INTEGER(consts_r)[0];
    double *XRandom = REAL(XRandom_r);
    int nObs = INTEGER(consts_r)[1];
    int pAbund = INTEGER(consts_r)[2];
    int pAbundRE = INTEGER(consts_r)[3];
    int nAbundRE = INTEGER(consts_r)[4];
    int saveFitted = INTEGER(consts_r)[5];
    int ppAbund = pAbund * pAbund; 
    double *muBeta = (double *) R_alloc(pAbund, sizeof(double));   
    F77_NAME(dcopy)(&pAbund, REAL(muBeta_r), &inc, muBeta, &inc);
    double *SigmaBeta = (double *) R_alloc(ppAbund, sizeof(double));
    F77_NAME(dcopy)(&ppAbund, REAL(SigmaBeta_r), &inc, SigmaBeta, &inc);
    double kappaA = REAL(kappaA_r)[0];
    double kappaB = REAL(kappaB_r)[0];
    double *sigmaSqMuA = REAL(sigmaSqMuA_r); 
    double *sigmaSqMuB = REAL(sigmaSqMuB_r); 
    int *nAbundRELong = INTEGER(nAbundRELong_r); 
    int *betaStarIndx = INTEGER(betaStarIndx_r); 
    int *betaLevelIndx = INTEGER(betaLevelIndx_r);
    int nBatch = INTEGER(nBatch_r)[0]; 
    int batchLength = INTEGER(batchLength_r)[0]; 
    int nSamples = nBatch * batchLength; 
    int nBurn = INTEGER(samplesInfo_r)[0]; 
    int nThin = INTEGER(samplesInfo_r)[1];
    int nPost = INTEGER(samplesInfo_r)[2]; 
    int currChain = INTEGER(chainInfo_r)[0];
    double acceptRate = REAL(acceptRate_r)[0];
    double* tuning = REAL(tuning_r);
    int nChain = INTEGER(chainInfo_r)[1];
    int nThreads = INTEGER(nThreads_r)[0];
    int verbose = INTEGER(verbose_r)[0];
    int nReport = INTEGER(nReport_r)[0];
    int status = 0; 
    int thinIndx = 0;
    int sPost = 0;  
    // NB = 1, Poisson = 0;
    int family = INTEGER(family_r)[0];

#ifdef _OPENMP
    omp_set_num_threads(nThreads);
#else
    if(nThreads > 1){
      warning("n.omp.threads > 1, but source not compiled with OpenMP support.");
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
	if (family == 1) {
          Rprintf("Negative binomial abundance model fit with %i sites.\n\n", J);
	} else {
          Rprintf("Poisson abundance model fit with %i sites.\n\n", J);
	}
        Rprintf("Samples per Chain: %i (%i batches of length %i)\n", nSamples, nBatch, batchLength);
        Rprintf("Burn-in: %i \n", nBurn); 
        Rprintf("Thinning Rate: %i \n", nThin); 
        Rprintf("Number of Chains: %i \n", nChain);
        Rprintf("Total Posterior Samples: %i \n\n", nPost * nChain); 
#ifdef _OPENMP
        Rprintf("Source compiled with OpenMP support and model fit using %i thread(s).\n\n", nThreads);
#else
        Rprintf("Source not compiled with OpenMP support.\n\n");
#endif
        Rprintf("Adaptive Metropolis with target acceptance rate: %.1f\n", 100*acceptRate);
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
    // Abundance Covariates 
    double *beta = (double *) R_alloc(pAbund, sizeof(double));   
    F77_NAME(dcopy)(&pAbund, REAL(betaStarting_r), &inc, beta, &inc);
    // Abundance random effect variances
    double *sigmaSqMu = (double *) R_alloc(pAbundRE, sizeof(double)); 
    F77_NAME(dcopy)(&pAbundRE, REAL(sigmaSqMuStarting_r), &inc, sigmaSqMu, &inc); 
    // Overdispersion parameter for NB;
    double kappa = REAL(kappaStarting_r)[0];
    // Latent random effects
    double *betaStar = (double *) R_alloc(nAbundRE, sizeof(double)); 
    F77_NAME(dcopy)(&nAbundRE, REAL(betaStarStarting_r), &inc, betaStar, &inc); 
    // Abundance fitted values
    double *yRep = (double *) R_alloc(nObs, sizeof(double)); zeros(yRep, nObs);

    /**********************************************************************
     * Return Stuff
     * *******************************************************************/
    SEXP betaSamples_r;
    PROTECT(betaSamples_r = allocMatrix(REALSXP, pAbund, nPost)); nProtect++;
    zeros(REAL(betaSamples_r), pAbund * nPost);
    SEXP yRepSamples_r; 
    SEXP muSamples_r; 
    SEXP likeSamples_r;
    if (saveFitted == 1) {
      PROTECT(yRepSamples_r = allocMatrix(REALSXP, nObs, nPost)); nProtect++; 
      zeros(REAL(yRepSamples_r), nObs * nPost);
      PROTECT(muSamples_r = allocMatrix(REALSXP, nObs, nPost)); nProtect++; 
      zeros(REAL(muSamples_r), nObs * nPost);
      PROTECT(likeSamples_r = allocMatrix(REALSXP, nObs, nPost)); nProtect++;
      zeros(REAL(likeSamples_r), nObs * nPost);
    }
    SEXP kappaSamples_r;
    if (family == 1) {
      PROTECT(kappaSamples_r = allocMatrix(REALSXP, inc, nPost)); nProtect++;
      zeros(REAL(kappaSamples_r), nPost);
    }
    // Abundance random effects
    SEXP sigmaSqMuSamples_r; 
    SEXP betaStarSamples_r; 
    if (pAbundRE > 0) {
      PROTECT(sigmaSqMuSamples_r = allocMatrix(REALSXP, pAbundRE, nPost)); nProtect++;
      zeros(REAL(sigmaSqMuSamples_r), pAbundRE * nPost);
      PROTECT(betaStarSamples_r = allocMatrix(REALSXP, nAbundRE, nPost)); nProtect++;
      zeros(REAL(betaStarSamples_r), nAbundRE * nPost);
    }
    
    /********************************************************************
      Some constants and temporary variables to be used later
    ********************************************************************/
    int nObspAbundRE = nObs * pAbundRE;
    double tmp_0; 
    double *tmp_nObs = (double *) R_alloc(nObs, sizeof(double)); 
   
    // For latent abundance and WAIC
    double *like = (double *) R_alloc(nObs, sizeof(double)); zeros(like, nObs);
    double *mu = (double *) R_alloc(nObs, sizeof(double)); 
    zeros(mu, nObs); 

    /********************************************************************
      Set up MH stuff
    ********************************************************************/
    double logPostBetaCurr = 0.0, logPostBetaCand = 0.0;
    double logPostKappaCurr = 0.0, logPostKappaCand = 0.0;
    double *logPostBetaStarCand = (double *) R_alloc(nAbundRE, sizeof(double));
    double *logPostBetaStarCurr = (double *) R_alloc(nAbundRE, sizeof(double));
    for (j = 0; j < nAbundRE; j++) {
      logPostBetaStarCurr[j] = R_NegInf;
      logPostBetaStarCand[j] = logPostBetaStarCurr[j];
    }
    int nAMCMC = 0;
    if (pAbundRE > 0) {
      nAMCMC = pAbund + nAbundRE;
    } else {
      nAMCMC = pAbund;
    }
    if (family == 1) {
      nAMCMC++;
    }
    int betaAMCMCIndx = 0;
    int betaStarAMCMCIndx = betaAMCMCIndx + pAbund;
    int kappaAMCMCIndx = betaStarAMCMCIndx + nAbundRE; // nAbundRE will be 0 if no REs. 
    double *accept = (double *) R_alloc(nAMCMC, sizeof(double)); zeros(accept, nAMCMC);
    // Set the initial candidate values for everything to the inital values. 
    double *betaCand = (double *) R_alloc(pAbund, sizeof(double)); 
    for (j = 0; j < pAbund; j++) {
      betaCand[j] = beta[j];
    } 
    double *betaStarCand = (double *) R_alloc(nAbundRE, sizeof(double));
    for (j = 0; j < nAbundRE; j++) {
      betaStarCand[j] = betaStar[j];
    }
    double kappaCand = 0.0;
    kappaCand = kappa;
    SEXP acceptSamples_r; 
    PROTECT(acceptSamples_r = allocMatrix(REALSXP, nAMCMC, nBatch)); nProtect++; 
    zeros(REAL(acceptSamples_r), nAMCMC * nBatch);
    SEXP tuningSamples_r; 
    PROTECT(tuningSamples_r = allocMatrix(REALSXP, nAMCMC, nBatch)); nProtect++; 
    zeros(REAL(tuningSamples_r), nAMCMC * nBatch);

    /**********************************************************************
     * Prep for random effects
     * *******************************************************************/
    // Site-level sums of the abundance random effects
    double *betaStarSites = (double *) R_alloc(nObs, sizeof(double)); 
    zeros(betaStarSites, nObs); 
    double *betaStarSitesCand = (double *) R_alloc(nObs, sizeof(double)); 
    int *betaStarLongIndx = (int *) R_alloc(nObspAbundRE, sizeof(int));
    // Initial sums
    for (j = 0; j < nObs; j++) {
      for (l = 0; l < pAbundRE; l++) {
        betaStarLongIndx[l * nObs + j] = which(XRE[l * nObs + j], betaLevelIndx, nAbundRE);
        betaStarSites[j] += betaStar[betaStarLongIndx[l * nObs + j]] * XRandom[l * nObs + j];
      }
      betaStarSitesCand[j] = betaStarSites[j];
    }
    // Starting index for abundance random effects
    int *betaStarStart = (int *) R_alloc(pAbundRE, sizeof(int)); 
    for (l = 0; l < pAbundRE; l++) {
      betaStarStart[l] = which(l, betaStarIndx, nAbundRE); 
    }

    logPostBetaCurr = R_NegInf;

    GetRNGstate(); 

    for (s = 0, g = 0; s < nBatch; s++) {
      for (t = 0; t < batchLength; t++, g++) {
        /********************************************************************
         *Update Abundance Regression Coefficients
         *******************************************************************/
        // Proposal
        for (k = 0; k < pAbund; k++) {
          logPostBetaCand = 0.0;
	  logPostBetaCurr = 0.0;
          betaCand[k] = rnorm(beta[k], exp(tuning[betaAMCMCIndx + k]));
          logPostBetaCand += dnorm(betaCand[k], muBeta[k], sqrt(SigmaBeta[k * pAbund + k]), 1);
	  logPostBetaCurr += dnorm(beta[k], muBeta[k], sqrt(SigmaBeta[k * pAbund + k]), 1);
          for (j = 0; j < nObs; j++) {
            tmp_nObs[j] = exp(F77_NAME(ddot)(&pAbund, &X[j], &nObs, betaCand, &inc) + betaStarSites[j]);
	    if (family == 1) {
              logPostBetaCand += dnbinom_mu(y[j], kappa, tmp_nObs[j] * offset[j], 1);
	    } else {
              logPostBetaCand += dpois(y[j], tmp_nObs[j] * offset[j], 1);
	    }
            tmp_nObs[j] = exp(F77_NAME(ddot)(&pAbund, &X[j], &nObs, beta, &inc) + betaStarSites[j]);
	    if (family == 1) {
              logPostBetaCurr += dnbinom_mu(y[j], kappa, tmp_nObs[j] * offset[j], 1);
	    } else {
              logPostBetaCurr += dpois(y[j], tmp_nObs[j] * offset[j], 1);
	    }
          }
          if (runif(0.0, 1.0) <= exp(logPostBetaCand - logPostBetaCurr)) {
            beta[k] = betaCand[k];
            accept[betaAMCMCIndx + k]++;
          } else {
            betaCand[k] = beta[k];
	  }
        }
        
        /********************************************************************
         *Update abundance random effects variance
         *******************************************************************/
        for (l = 0; l < pAbundRE; l++) {
          tmp_0 = F77_NAME(ddot)(&nAbundRELong[l], &betaStar[betaStarStart[l]], &inc, &betaStar[betaStarStart[l]], &inc); 
          tmp_0 *= 0.5; 
          sigmaSqMu[l] = rigamma(sigmaSqMuA[l] + nAbundRELong[l] / 2.0, sigmaSqMuB[l] + tmp_0); 
        }

        /********************************************************************
         *Update abundance random effects
         *******************************************************************/
        if (pAbundRE > 0) {
          for (l = 0; l < nAbundRE; l++) {
	    betaStarCand[l] = rnorm(betaStar[l], exp(tuning[betaStarAMCMCIndx + l]));
            logPostBetaStarCand[l] = dnorm(betaStarCand[l], 0.0, 
			                   sqrt(sigmaSqMu[betaStarIndx[l]]), 1);
            logPostBetaStarCurr[l] = dnorm(betaStar[l], 0.0, 
			                   sqrt(sigmaSqMu[betaStarIndx[l]]), 1);
	    for (j = 0; j < nObs; j++) {
              if (XRE[betaStarIndx[l] * nObs + j] == betaLevelIndx[l]) {
                // Candidate
                betaStarSitesCand[j] = 0.0;
                for (ll = 0; ll < pAbundRE; ll++) {
                  betaStarSitesCand[j] += betaStarCand[betaStarLongIndx[ll * nObs + j]] * 
	                              XRandom[ll * nObs + j];
                }
                tmp_nObs[j] = exp(F77_NAME(ddot)(&pAbund, &X[j], &nObs, beta, &inc) + 
				  betaStarSitesCand[j]);
                if (family == 1) {
		  logPostBetaStarCand[l] += dnbinom_mu(y[j], kappa, tmp_nObs[j] * offset[j], 1);
		} else {
		  logPostBetaStarCand[l] += dpois(y[j], tmp_nObs[j] * offset[j], 1);
		}
		// Current
                betaStarSites[j] = 0.0;
                for (ll = 0; ll < pAbundRE; ll++) {
                  betaStarSites[j] += betaStar[betaStarLongIndx[ll * nObs + j]] * 
	                              XRandom[ll * nObs + j];
                }
                tmp_nObs[j] = exp(F77_NAME(ddot)(&pAbund, &X[j], &nObs, beta, &inc) + 
				  betaStarSites[j]);
                if (family == 1) {
		  logPostBetaStarCurr[l] += dnbinom_mu(y[j], kappa, tmp_nObs[j] * offset[j], 1);
		} else {
		  logPostBetaStarCurr[l] += dpois(y[j], tmp_nObs[j] * offset[j], 1);
		}
	      }
	    }
	    if (runif (0.0, 1.0) <= exp(logPostBetaStarCand[l] - logPostBetaStarCurr[l])) {
              betaStar[l] = betaStarCand[l];
	      F77_NAME(dcopy)(&nObs, betaStarSitesCand, &inc, betaStarSites, &inc);
	      accept[betaStarAMCMCIndx + l]++;
	    } else {
              betaStarCand[l] = betaStar[l];
	      F77_NAME(dcopy)(&nObs, betaStarSites, &inc, betaStarSitesCand, &inc);
	    }
	  }
	}

        /********************************************************************
         *Update kappa (the NB size parameter)
         *******************************************************************/
        if (family == 1) {
          kappaCand = logitInv(rnorm(logit(kappa, kappaA, kappaB), exp(tuning[kappaAMCMCIndx])), 
			       kappaA, kappaB);
	  logPostKappaCurr = 0.0;
	  logPostKappaCand = 0.0;
	  for (j = 0; j < nObs; j++) {
            mu[j] = exp(F77_NAME(ddot)(&pAbund, &X[j], &nObs, beta, &inc) + 
                        betaStarSites[j]);
            logPostKappaCurr += dnbinom_mu(y[j], kappa, mu[j] * offset[j], 1);
	    logPostKappaCand += dnbinom_mu(y[j], kappaCand, mu[j] * offset[j], 1);
	  }
	  // Jacobian adjustment
	  logPostKappaCurr += log(kappa - kappaA) + log(kappaB - kappa);
	  logPostKappaCand += log(kappaCand - kappaA) + log(kappaB - kappaCand);
          if (runif(0.0, 1.0) <= exp(logPostKappaCand - logPostKappaCurr)) {
            kappa = kappaCand;
	    accept[kappaAMCMCIndx]++;
	  }
	}

        /********************************************************************
         *Get fitted values
         *******************************************************************/
	if (saveFitted == 1) {
          for (j = 0; j < nObs; j++) {
            // Only calculate if Poisson since it's already calculated in kappa update
	    if (family == 0) {
              mu[j] = exp(F77_NAME(ddot)(&pAbund, &X[j], &nObs, beta, &inc) + 
                          betaStarSites[j]);
              yRep[j] = rpois(mu[j] * offset[j]);
              like[j] = dpois(y[j], mu[j] * offset[j], 0);
	    } else {
              yRep[j] = rnbinom_mu(kappa, mu[j] * offset[j]);
              like[j] = dnbinom_mu(y[j], kappa, mu[j] * offset[j], 0);
	    }
          }
	}

        /********************************************************************
         *Save samples
         *******************************************************************/
        if (g >= nBurn) {
          thinIndx++; 
          if (thinIndx == nThin) {
            F77_NAME(dcopy)(&pAbund, beta, &inc, &REAL(betaSamples_r)[sPost*pAbund], &inc);
	    if (saveFitted == 1) {
              F77_NAME(dcopy)(&nObs, yRep, &inc, &REAL(yRepSamples_r)[sPost*nObs], &inc); 
              F77_NAME(dcopy)(&nObs, mu, &inc, &REAL(muSamples_r)[sPost*nObs], &inc); 
              F77_NAME(dcopy)(&nObs, like, &inc, 
          		      &REAL(likeSamples_r)[sPost*nObs], &inc);
	    }
	    if (family == 1) {
	      REAL(kappaSamples_r)[sPost] = kappa;
	    }
            if (pAbundRE > 0) {
              F77_NAME(dcopy)(&pAbundRE, sigmaSqMu, &inc, 
          		    &REAL(sigmaSqMuSamples_r)[sPost*pAbundRE], &inc);
              F77_NAME(dcopy)(&nAbundRE, betaStar, &inc, 
          		    &REAL(betaStarSamples_r)[sPost*nAbundRE], &inc);
            }
            sPost++; 
            thinIndx = 0; 
          }
        }
        R_CheckUserInterrupt();
      } // t (end batch)
      /********************************************************************
       *Adjust tuning 
       *******************************************************************/
      for (j = 0; j < nAMCMC; j++) {
        REAL(acceptSamples_r)[s * nAMCMC + j] = accept[j]/batchLength; 
        REAL(tuningSamples_r)[s * nAMCMC + j] = tuning[j]; 
        if (accept[j] / batchLength > acceptRate) {
          tuning[j] += std::min(0.01, 1.0/sqrt(static_cast<double>(s)));
        } else{
            tuning[j] -= std::min(0.01, 1.0/sqrt(static_cast<double>(s)));
          }
        accept[j] = 0;
      }

      /********************************************************************
       *Report 
       *******************************************************************/
      if (verbose) {
        if (status == nReport) {
          Rprintf("Batch: %i of %i, %3.2f%%\n", s, nBatch, 100.0*s/nBatch);
          Rprintf("\tParameter\tAcceptance\tTuning\n");	  
          for (j = 0; j < pAbund; j++) {
            Rprintf("\tbeta[%i]\t\t%3.1f\t\t%1.5f\n", j + 1, 100.0*REAL(acceptSamples_r)[s * nAMCMC + betaAMCMCIndx + j], exp(tuning[betaAMCMCIndx + j]));
          }
	  if (family == 1) {
            Rprintf("\tkappa\t\t%3.1f\t\t%1.5f\n", 100.0*REAL(acceptSamples_r)[s * nAMCMC + kappaAMCMCIndx], exp(tuning[kappaAMCMCIndx]));
	  }
          Rprintf("-------------------------------------------------\n");
          #ifdef Win32
          R_FlushConsole();
          #endif
          status = 0;
        }
      } 
      status++;        

    } // all batches
    if (verbose) {
      Rprintf("Batch: %i of %i, %3.2f%%\n", s, nBatch, 100.0*s/nBatch);
    }
    PutRNGstate();

    SEXP result_r, resultName_r;
    int nResultListObjs = 4;
    if (pAbundRE > 0) {
      nResultListObjs += 2;
    }
    if (family == 1) {
      nResultListObjs += 1;
    }

    PROTECT(result_r = allocVector(VECSXP, nResultListObjs)); nProtect++;
    PROTECT(resultName_r = allocVector(VECSXP, nResultListObjs)); nProtect++;

    // Setting the components of the output list.
    SET_VECTOR_ELT(result_r, 0, betaSamples_r);
    if (saveFitted == 1) {
      SET_VECTOR_ELT(result_r, 1, yRepSamples_r); 
      SET_VECTOR_ELT(result_r, 2, muSamples_r);
      SET_VECTOR_ELT(result_r, 3, likeSamples_r);
    }
    if (pAbundRE > 0) {
      SET_VECTOR_ELT(result_r, 4, sigmaSqMuSamples_r);
      SET_VECTOR_ELT(result_r, 5, betaStarSamples_r);
    }
    if (family == 1) {
      if (pAbundRE > 0) {
        tmp_0 = 6;
      } else {
        tmp_0 = 4;
      }
      SET_VECTOR_ELT(result_r, tmp_0, kappaSamples_r);
    }

    SET_VECTOR_ELT(resultName_r, 0, mkChar("beta.samples")); 
    if (saveFitted == 1) {
      SET_VECTOR_ELT(resultName_r, 1, mkChar("y.rep.samples")); 
      SET_VECTOR_ELT(resultName_r, 2, mkChar("mu.samples"));
      SET_VECTOR_ELT(resultName_r, 3, mkChar("like.samples"));
    }
    if (pAbundRE > 0) {
      SET_VECTOR_ELT(resultName_r, 4, mkChar("sigma.sq.mu.samples")); 
      SET_VECTOR_ELT(resultName_r, 5, mkChar("beta.star.samples")); 
    }
    if (family == 1) {
      SET_VECTOR_ELT(resultName_r, tmp_0, mkChar("kappa.samples")); 
    }
   
    namesgets(result_r, resultName_r);
    
    UNPROTECT(nProtect);
    
    return(result_r);
  }
}

