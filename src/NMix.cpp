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
  SEXP NMix(SEXP y_r, SEXP X_r, SEXP Xp_r, SEXP XRE_r, SEXP XpRE_r, 
            SEXP XRandom_r, SEXP XpRandom_r, SEXP yMax_r,
            SEXP consts_r, SEXP K_r, SEXP nAbundRELong_r, SEXP nDetRELong_r,
            SEXP betaStarting_r, SEXP alphaStarting_r, SEXP kappaStarting_r,
            SEXP sigmaSqMuStarting_r, SEXP sigmaSqPStarting_r,
            SEXP betaStarStarting_r, SEXP alphaStarStarting_r, SEXP NStarting_r,
            SEXP NLongIndx_r, SEXP betaStarIndx_r, SEXP betaLevelIndx_r, 
	    SEXP alphaStarIndx_r, SEXP alphaLevelIndx_r,
	    SEXP muBeta_r, SEXP SigmaBeta_r, 
	    SEXP muAlpha_r, SEXP SigmaAlpha_r, 
            SEXP sigmaSqMuA_r, SEXP sigmaSqMuB_r, 
            SEXP sigmaSqPA_r, SEXP sigmaSqPB_r, 
	    SEXP kappaA_r, SEXP kappaB_r, SEXP tuning_r,
	    SEXP nBatch_r, SEXP batchLength_r, SEXP acceptRate_r, SEXP nThreads_r, 
            SEXP verbose_r, SEXP nReport_r, SEXP samplesInfo_r,
	    SEXP chainInfo_r, SEXP family_r, SEXP offset_r){
   
    /**********************************************************************
     * Initial constants
     * *******************************************************************/
    int i, g, t, j, s, l, k, ll, nProtect=0;
    const int inc = 1;
    const double one = 1.0;
    const double zero = 0.0;
    
    /**********************************************************************
     * Get Inputs
     * *******************************************************************/
    double *y = REAL(y_r);
    double *X = REAL(X_r);
    double *Xp = REAL(Xp_r);
    double *yMax = REAL(yMax_r);
    double *offset = REAL(offset_r);
    int *XRE = INTEGER(XRE_r); 
    int *XpRE = INTEGER(XpRE_r);
    double *XRandom = REAL(XRandom_r);
    double *XpRandom = REAL(XpRandom_r);
    // Load constants
    int J = INTEGER(consts_r)[0];
    int nObs = INTEGER(consts_r)[1];
    int pAbund = INTEGER(consts_r)[2];
    int pAbundRE = INTEGER(consts_r)[3];
    int nAbundRE = INTEGER(consts_r)[4];
    int pDet = INTEGER(consts_r)[5];
    int pDetRE = INTEGER(consts_r)[6];
    int nDetRE = INTEGER(consts_r)[7];
    int ppDet = pDet * pDet;
    int ppAbund = pAbund * pAbund; 
    double *muBeta = (double *) R_alloc(pAbund, sizeof(double));   
    F77_NAME(dcopy)(&pAbund, REAL(muBeta_r), &inc, muBeta, &inc);
    double *muAlpha = (double *) R_alloc(pDet, sizeof(double));   
    F77_NAME(dcopy)(&pDet, REAL(muAlpha_r), &inc, muAlpha, &inc);
    double *SigmaBeta = (double *) R_alloc(ppAbund, sizeof(double));   
    F77_NAME(dcopy)(&ppAbund, REAL(SigmaBeta_r), &inc, SigmaBeta, &inc);
    double *SigmaAlpha = (double *) R_alloc(ppDet, sizeof(double));   
    F77_NAME(dcopy)(&ppDet, REAL(SigmaAlpha_r), &inc, SigmaAlpha, &inc);
    double *sigmaSqMuA = REAL(sigmaSqMuA_r); 
    double *sigmaSqMuB = REAL(sigmaSqMuB_r); 
    double *sigmaSqPA = REAL(sigmaSqPA_r); 
    double *sigmaSqPB = REAL(sigmaSqPB_r); 
    double kappaA = REAL(kappaA_r)[0];
    double kappaB = REAL(kappaB_r)[0];
    int *nAbundRELong = INTEGER(nAbundRELong_r); 
    int *nDetRELong = INTEGER(nDetRELong_r); 
    int *NLongIndx = INTEGER(NLongIndx_r); 
    int *alphaStarIndx = INTEGER(alphaStarIndx_r); 
    int *alphaLevelIndx = INTEGER(alphaLevelIndx_r);
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
    double *tuning = REAL(tuning_r);
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
          Rprintf("Negative Binomial N-mixture model with %i sites.\n\n", J);
	} else {
          Rprintf("Poisson N-mixture model with %i sites.\n\n", J);
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
    // Latent random effects
    double *betaStar = (double *) R_alloc(nAbundRE, sizeof(double)); 
    F77_NAME(dcopy)(&nAbundRE, REAL(betaStarStarting_r), &inc, betaStar, &inc); 
    // Detection covariates
    double *alpha = (double *) R_alloc(pDet, sizeof(double));   
    F77_NAME(dcopy)(&pDet, REAL(alphaStarting_r), &inc, alpha, &inc);
    // Detection random effect variances
    double *sigmaSqP = (double *) R_alloc(pDetRE, sizeof(double)); 
    F77_NAME(dcopy)(&pDetRE, REAL(sigmaSqPStarting_r), &inc, sigmaSqP, &inc); 
    // Latent detection random effects
    double *alphaStar = (double *) R_alloc(nDetRE, sizeof(double)); 
    F77_NAME(dcopy)(&nDetRE, REAL(alphaStarStarting_r), &inc, alphaStar, &inc); 
    // Size parameter
    double kappa = REAL(kappaStarting_r)[0];
    // Latent Abundance
    double *N = (double *) R_alloc(J, sizeof(double));   
    F77_NAME(dcopy)(&J, REAL(NStarting_r), &inc, N, &inc);

    /**********************************************************************
     * Return Stuff
     * *******************************************************************/
    SEXP betaSamples_r;
    PROTECT(betaSamples_r = allocMatrix(REALSXP, pAbund, nPost)); nProtect++;
    zeros(REAL(betaSamples_r), pAbund * nPost);
    SEXP alphaSamples_r; 
    PROTECT(alphaSamples_r = allocMatrix(REALSXP, pDet, nPost)); nProtect++;
    zeros(REAL(alphaSamples_r), pDet * nPost);
    SEXP NSamples_r; 
    PROTECT(NSamples_r = allocMatrix(REALSXP, J, nPost)); nProtect++; 
    zeros(REAL(NSamples_r), J * nPost);
    SEXP kappaSamples_r;
    if (family == 1) {
      PROTECT(kappaSamples_r = allocMatrix(REALSXP, inc, nPost)); nProtect++;
      zeros(REAL(kappaSamples_r), nPost);
    }
    SEXP muSamples_r; 
    PROTECT(muSamples_r = allocMatrix(REALSXP, J, nPost)); nProtect++; 
    zeros(REAL(muSamples_r), J * nPost);
    // Detection random effects
    SEXP sigmaSqPSamples_r; 
    SEXP alphaStarSamples_r; 
    if (pDetRE > 0) {
      PROTECT(sigmaSqPSamples_r = allocMatrix(REALSXP, pDetRE, nPost)); nProtect++;
      zeros(REAL(sigmaSqPSamples_r), pDetRE * nPost);
      PROTECT(alphaStarSamples_r = allocMatrix(REALSXP, nDetRE, nPost)); nProtect++;
      zeros(REAL(alphaStarSamples_r), nDetRE * nPost);
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
    double tmp_0, tmp_02; 
    double *tmp_nObs = (double *) R_alloc(nObs, sizeof(double)); 
    double *tmp_J = (double *) R_alloc(J, sizeof(double));
   
    // For latent abundance
    double *detProb = (double *) R_alloc(nObs, sizeof(double)); zeros(detProb, nObs);
    double *mu = (double *) R_alloc(J, sizeof(double)); 
    zeros(mu, J); 

    /********************************************************************
      Set up MH stuff
    ********************************************************************/
    double logPostBetaCurr = 0.0, logPostBetaCand = 0.0;
    double logPostAlphaCurr = 0.0, logPostAlphaCand = 0.0;
    double logPostKappaCurr = 0.0, logPostKappaCand = 0.0;
    double *logPostBetaStarCand = (double *) R_alloc(nAbundRE, sizeof(double));
    double *logPostBetaStarCurr = (double *) R_alloc(nAbundRE, sizeof(double));
    for (j = 0; j < nAbundRE; j++) {
      logPostBetaStarCurr[j] = R_NegInf;
      logPostBetaStarCand[j] = logPostBetaStarCurr[j];
    }
    double *logPostAlphaStarCand = (double *) R_alloc(nDetRE, sizeof(double));
    double *logPostAlphaStarCurr = (double *) R_alloc(nDetRE, sizeof(double));
    for (j = 0; j < nDetRE; j++) {
      logPostAlphaStarCurr[j] = R_NegInf;
      logPostAlphaStarCand[j] = logPostAlphaStarCurr[j];
    }
    int nAMCMC = pAbund + pDet;
    if (pAbundRE > 0) {
      nAMCMC += nAbundRE;
    }
    if (pDetRE > 0) {
      nAMCMC += nDetRE;
    }
    if (family == 1) {
      nAMCMC++;
    }
    int betaAMCMCIndx = 0;
    int alphaAMCMCIndx = betaAMCMCIndx + pAbund;
    int betaStarAMCMCIndx = alphaAMCMCIndx + pDet;
    int alphaStarAMCMCIndx = betaStarAMCMCIndx + nAbundRE;
    int kappaAMCMCIndx = alphaStarAMCMCIndx + nDetRE;
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
    double *alphaCand = (double *) R_alloc(pDet, sizeof(double)); 
    for (j = 0; j < pDet; j++) {
      alphaCand[j] = alpha[j];
    } 
    double *alphaStarCand = (double *) R_alloc(nDetRE, sizeof(double));
    for (j = 0; j < nDetRE; j++) {
      alphaStarCand[j] = alphaStar[j];
    }
    double kappaCand = 0.0;
    kappaCand = kappa;
    double *logPostCurrN = (double *) R_alloc(J, sizeof(double));
    zeros(logPostCurrN, J);
    double *logPostCandN = (double *) R_alloc(J, sizeof(double));
    zeros(logPostCandN, J);
    double epsilonN = 1;
    double *NCand = (double *) R_alloc(J, sizeof(double));
    for (j = 0; j < J; j++) {
      NCand[j] = N[j];
    }
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
    double *betaStarSites = (double *) R_alloc(J, sizeof(double)); 
    zeros(betaStarSites, J); 
    double *betaStarSitesCand = (double *) R_alloc(J, sizeof(double)); 
    // Initial sums
    for (j = 0; j < J; j++) {
      for (l = 0; l < pAbundRE; l++) {
        betaStarSites[j] += betaStar[which(XRE[l * J + j], betaLevelIndx, nAbundRE)] * 
                            XRandom[l * J + j];
      }
      betaStarSitesCand[j] = betaStarSites[j];
    }
    // Observation-level sums of the detection random effects
    double *alphaStarObs = (double *) R_alloc(nObs, sizeof(double)); 
    zeros(alphaStarObs, nObs); 
    double *alphaStarObsCand = (double *) R_alloc(nObs, sizeof(double));
    // Get sums of the current REs for each site/visit combo
    for (i = 0; i < nObs; i++) {
      for (l = 0; l < pDetRE; l++) {
        alphaStarObs[i] += alphaStar[which(XpRE[l * nObs + i], alphaLevelIndx, nDetRE)] * 
                           XpRandom[l * nObs + i];
      }
      alphaStarObsCand[i] = alphaStarObs[i];
    }
    // Starting index for abundance random effects
    int *betaStarStart = (int *) R_alloc(pAbundRE, sizeof(int)); 
    for (l = 0; l < pAbundRE; l++) {
      betaStarStart[l] = which(l, betaStarIndx, nAbundRE); 
    }
    // Starting index for detection random effects
    int *alphaStarStart = (int *) R_alloc(pDetRE, sizeof(int)); 
    for (l = 0; l < pDetRE; l++) {
      alphaStarStart[l] = which(l, alphaStarIndx, nDetRE); 
    }

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
          for (i = 0; i < pAbund; i++) {
            logPostBetaCand += dnorm(betaCand[i], muBeta[i], sqrt(SigmaBeta[i * pAbund + i]), 1);
	    logPostBetaCurr += dnorm(beta[i], muBeta[i], sqrt(SigmaBeta[i * pAbund + i]), 1);
          }
          for (j = 0; j < J; j++) {
            tmp_J[j] = exp(F77_NAME(ddot)(&pAbund, &X[j], &J, betaCand, &inc) + betaStarSites[j]);
	    if (family == 1) {
              logPostBetaCand += nb_logpost(kappa, N[j], tmp_J[j], offset[j]);
	    } else {
              logPostBetaCand += poisson_logpost(N[j], tmp_J[j], offset[j]);
	    }
            tmp_J[j] = exp(F77_NAME(ddot)(&pAbund, &X[j], &J, beta, &inc) + betaStarSites[j]);
	    if (family == 1) {
              logPostBetaCurr += nb_logpost(kappa, N[j], tmp_J[j], offset[j]);
	    } else {
              logPostBetaCurr += poisson_logpost(N[j], tmp_J[j], offset[j]);
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
         *Update Detection Regression Coefficients
         *******************************************************************/
        // Proposal
        for (k = 0; k < pDet; k++) {
          logPostAlphaCand = 0.0;
	  logPostAlphaCurr = 0.0;
          alphaCand[k] = rnorm(alpha[k], exp(tuning[alphaAMCMCIndx + k]));
          for (i = 0; i < pDet; i++) {
            logPostAlphaCand += dnorm(alphaCand[i], muAlpha[i], sqrt(SigmaAlpha[i * pDet + i]), 1);
	    logPostAlphaCurr += dnorm(alpha[i], muAlpha[i], sqrt(SigmaAlpha[i * pDet + i]), 1);
          }
          for (i = 0; i < nObs; i++) {
            if (N[NLongIndx[i]] > 0.0) {
              tmp_nObs[i] = logitInv(F77_NAME(ddot)(&pDet, &Xp[i], &nObs, alphaCand, &inc) + 
			             alphaStarObs[i], zero, one);
              logPostAlphaCand += dbinom(y[i], N[NLongIndx[i]], tmp_nObs[i], 1);	     
              tmp_nObs[i] = logitInv(F77_NAME(ddot)(&pDet, &Xp[i], &nObs, alpha, &inc) + 
			             alphaStarObs[i], zero, one);
              logPostAlphaCurr += dbinom(y[i], N[NLongIndx[i]], tmp_nObs[i], 1);	     
	    }
          }
          if (runif(0.0, 1.0) <= exp(logPostAlphaCand - logPostAlphaCurr)) {
            alpha[k] = alphaCand[k];
            accept[alphaAMCMCIndx + k]++;
          } else {
            alphaCand[k] = alpha[k];
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
         *Update detection  random effects variance
         *******************************************************************/
        for (l = 0; l < pDetRE; l++) {
          tmp_0 = F77_NAME(ddot)(&nDetRELong[l], &alphaStar[alphaStarStart[l]], &inc, &alphaStar[alphaStarStart[l]], &inc); 
          tmp_0 *= 0.5; 
          sigmaSqP[l] = rigamma(sigmaSqPA[l] + nDetRELong[l] / 2.0, sigmaSqPB[l] + tmp_0); 
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
	    for (j = 0; j < J; j++) {
              if (XRE[betaStarIndx[l] * J + j] == betaLevelIndx[l]) {
                // Candidate
                betaStarSitesCand[j] = 0.0;
                for (ll = 0; ll < pAbundRE; ll++) {
                  betaStarSitesCand[j] += betaStarCand[which(XRE[ll * J + j], 
				                         betaLevelIndx, nAbundRE)] * 
	                              XRandom[ll * J + j];
                }
                tmp_J[j] = exp(F77_NAME(ddot)(&pAbund, &X[j], &J, beta, &inc) + 
				  betaStarSitesCand[j]);
                if (family == 1) {
		  logPostBetaStarCand[l] += nb_logpost(kappa, N[j], tmp_J[j], offset[j]);
		} else {
		  logPostBetaStarCand[l] += poisson_logpost(N[j], tmp_J[j], offset[j]);
		}
		// Current
                betaStarSites[j] = 0.0;
                for (ll = 0; ll < pAbundRE; ll++) {
                  betaStarSites[j] += betaStar[which(XRE[ll * J + j], 
				               betaLevelIndx, nAbundRE)] * 
	                              XRandom[ll * J + j];
                }
                tmp_J[j] = exp(F77_NAME(ddot)(&pAbund, &X[j], &J, beta, &inc) + 
				  betaStarSites[j]);
                if (family == 1) {
		  logPostBetaStarCurr[l] += nb_logpost(kappa, N[j], tmp_J[j], offset[j]);
		} else {
		  logPostBetaStarCurr[l] += poisson_logpost(N[j], tmp_J[j], offset[j]);
		}
	      }
	    }
	    if (runif (0.0, 1.0) <= exp(logPostBetaStarCand[l] - logPostBetaStarCurr[l])) {
              betaStar[l] = betaStarCand[l];
	      F77_NAME(dcopy)(&J, betaStarSitesCand, &inc, betaStarSites, &inc);
	      accept[betaStarAMCMCIndx + l]++;
	    } else {
              betaStarCand[l] = betaStar[l];
	      F77_NAME(dcopy)(&J, betaStarSites, &inc, betaStarSitesCand, &inc);
	    }
	  }
	}

        /********************************************************************
         *Update detection random effects
         *******************************************************************/
        if (pDetRE > 0) {
          for (l = 0; l < nDetRE; l++) {
	    alphaStarCand[l] = rnorm(alphaStar[l], exp(tuning[alphaStarAMCMCIndx + l]));
            logPostAlphaStarCand[l] = dnorm(alphaStarCand[l], 0.0, 
			                   sqrt(sigmaSqP[alphaStarIndx[l]]), 1);
            logPostAlphaStarCurr[l] = dnorm(alphaStar[l], 0.0, 
			                   sqrt(sigmaSqP[alphaStarIndx[l]]), 1);
	    for (i = 0; i < nObs; i++) {
              if ((N[NLongIndx[i]] > 0) && (XpRE[alphaStarIndx[l] * nObs + i] == alphaLevelIndx[l])) {
                // Candidate
                alphaStarObsCand[i] = 0.0;
                for (ll = 0; ll < pDetRE; ll++) {
                  alphaStarObsCand[i] += alphaStarCand[which(XpRE[ll * nObs + i], 
				                         alphaLevelIndx, nDetRE)] * 
	                              XpRandom[ll * nObs + i];
                }
                tmp_nObs[i] = logitInv(F77_NAME(ddot)(&pDet, &Xp[i], &nObs, alpha, &inc) + 
			               alphaStarObsCand[i], zero, one);
                logPostAlphaStarCand[l] += dbinom(y[i], N[NLongIndx[i]], tmp_nObs[i], 1);	     
		// Current
                alphaStarObs[i] = 0.0;
                for (ll = 0; ll < pDetRE; ll++) {
                  alphaStarObs[i] += alphaStar[which(XpRE[ll * nObs + i], 
				               alphaLevelIndx, nDetRE)] * 
	                              XpRandom[ll * nObs + i];
                }
                tmp_nObs[i] = logitInv(F77_NAME(ddot)(&pDet, &Xp[i], &nObs, alpha, &inc) + 
			               alphaStarObs[i], zero, one);
                logPostAlphaStarCurr[l] += dbinom(y[i], N[NLongIndx[i]], tmp_nObs[i], 1);	     
	      }
	    }
	    if (runif (0.0, 1.0) <= exp(logPostAlphaStarCand[l] - logPostAlphaStarCurr[l])) {
              alphaStar[l] = alphaStarCand[l];
	      F77_NAME(dcopy)(&nObs, alphaStarObsCand, &inc, alphaStarObs, &inc);
	      accept[alphaStarAMCMCIndx + l]++;
	    } else {
              alphaStarCand[l] = alphaStar[l];
	      F77_NAME(dcopy)(&nObs, alphaStarObs, &inc, alphaStarObsCand, &inc);
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
	  for (j = 0; j < J; j++) {
            mu[j] = exp(F77_NAME(ddot)(&pAbund, &X[j], &J, beta, &inc) + 
                        betaStarSites[j]);
            logPostKappaCurr += nb_logpost(kappa, N[j], mu[j], offset[j]);
	    logPostKappaCand += nb_logpost(kappaCand, N[j], mu[j], offset[j]);
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
         *Update Latent Abundance 
         *******************************************************************/
	zeros(logPostCurrN, J);
	zeros(logPostCandN, J);
	// Proposal
	for (j = 0; j < J; j++) {
          NCand[j] = rpois(N[j] + epsilonN);
	  // Only calculate if Poisson since its already calculated in kappa update
	  if (family == 0) {
            mu[j] = exp(F77_NAME(ddot)(&pAbund, &X[j], &J, beta, &inc) + 
                        betaStarSites[j]);
	  }
	}
	// Likelihood contribution to Metropolis ratios
	for (i = 0; i < nObs; i++) {
          detProb[i] = logitInv(F77_NAME(ddot)(&pDet, &Xp[i], &nObs, alpha, &inc) + 
	        	        alphaStarObs[i], zero, one);
          logPostCurrN[NLongIndx[i]] += dbinom(y[i], N[NLongIndx[i]], detProb[i], 1);
	  logPostCandN[NLongIndx[i]] += dbinom(y[i], NCand[NLongIndx[i]], detProb[i], 1);
	}
        for (j = 0; j < J; j++) {
	  if (NCand[j] >= yMax[j]) {
            /********************************
             * Current
             *******************************/
	    // Contribution from NB or Poisson
	    if (family == 0) {
              logPostCurrN[j] += poisson_logpost(N[j], mu[j], offset[j]);
	    } else {
              logPostCurrN[j] += nb_logpost(kappa, N[j], mu[j], offset[j]);
	    }
	    // MH contribution for assymetric proposal distribution.
	    logPostCurrN[j] += poisson_logpost(NCand[j], N[j] + epsilonN, 1.0);
            /********************************
             * Candidate
             *******************************/
	    if (family == 0) {
              logPostCandN[j] += poisson_logpost(NCand[j], mu[j], offset[j]);
	    } else {
              logPostCandN[j] += nb_logpost(kappa, NCand[j], mu[j], offset[j]);
	    }
	    // MH contribution for assymetric proposal distribution.
	    logPostCandN[j] += poisson_logpost(N[j], NCand[j] + epsilonN, 1.0);
            if (runif(0.0,1.0) <= exp(logPostCandN[j] - logPostCurrN[j])) {
              N[j] = NCand[j];
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
            F77_NAME(dcopy)(&pDet, alpha, &inc, &REAL(alphaSamples_r)[sPost*pDet], &inc);
	    if (family == 1) {
	      REAL(kappaSamples_r)[sPost] = kappa;
	    }
            F77_NAME(dcopy)(&J, mu, &inc, &REAL(muSamples_r)[sPost*J], &inc); 
            F77_NAME(dcopy)(&J, N, &inc, &REAL(NSamples_r)[sPost*J], &inc); 
            if (pAbundRE > 0) {
              F77_NAME(dcopy)(&pAbundRE, sigmaSqMu, &inc, 
          		    &REAL(sigmaSqMuSamples_r)[sPost*pAbundRE], &inc);
              F77_NAME(dcopy)(&nAbundRE, betaStar, &inc, 
          		    &REAL(betaStarSamples_r)[sPost*nAbundRE], &inc);
            }
            if (pDetRE > 0) {
              F77_NAME(dcopy)(&pDetRE, sigmaSqP, &inc, 
                  	      &REAL(sigmaSqPSamples_r)[sPost*pDetRE], &inc);
              F77_NAME(dcopy)(&nDetRE, alphaStar, &inc, 
                  	      &REAL(alphaStarSamples_r)[sPost*nDetRE], &inc);
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
          for (j = 0; j < pDet; j++) {
            Rprintf("\talpha[%i]\t%3.1f\t\t%1.5f\n", j + 1, 100.0*REAL(acceptSamples_r)[s * nAMCMC + alphaAMCMCIndx + j], exp(tuning[alphaAMCMCIndx + j]));
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
    int nResultListObjs = 5;
    if (pDetRE > 0) {
      nResultListObjs += 2; 
    }
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
    SET_VECTOR_ELT(result_r, 1, alphaSamples_r);
    SET_VECTOR_ELT(result_r, 2, NSamples_r); 
    SET_VECTOR_ELT(result_r, 3, muSamples_r);
    SET_VECTOR_ELT(result_r, 4, tuningSamples_r);
    if (pDetRE > 0) {
      tmp_0 = 5; // Needed to make tracking kappa easier. 
      SET_VECTOR_ELT(result_r, 5, sigmaSqPSamples_r);
      SET_VECTOR_ELT(result_r, 6, alphaStarSamples_r);
    }
    if (pAbundRE > 0) {
      if (pDetRE > 0) {
        tmp_0 = 7;
      } else {
        tmp_0 = 5;
      }
      SET_VECTOR_ELT(result_r, tmp_0, sigmaSqMuSamples_r);
      SET_VECTOR_ELT(result_r, tmp_0 + 1, betaStarSamples_r);
    }
    if (family == 1) {
      if ((pDetRE > 0) || (pAbundRE > 0)) {
        tmp_02 = tmp_0 + 2;
      } else {
        tmp_02 = 5;
      }  
      SET_VECTOR_ELT(result_r, tmp_02, kappaSamples_r);
    }

    SET_VECTOR_ELT(resultName_r, 0, mkChar("beta.samples")); 
    SET_VECTOR_ELT(resultName_r, 1, mkChar("alpha.samples")); 
    SET_VECTOR_ELT(resultName_r, 2, mkChar("N.samples")); 
    SET_VECTOR_ELT(resultName_r, 3, mkChar("mu.samples"));
    SET_VECTOR_ELT(resultName_r, 4, mkChar("tuning.samples"));
    if (pDetRE > 0) {
      SET_VECTOR_ELT(resultName_r, 5, mkChar("sigma.sq.p.samples")); 
      SET_VECTOR_ELT(resultName_r, 6, mkChar("alpha.star.samples")); 
    }
    if (pAbundRE > 0) {
      SET_VECTOR_ELT(resultName_r, tmp_0, mkChar("sigma.sq.mu.samples")); 
      SET_VECTOR_ELT(resultName_r, tmp_0 + 1, mkChar("beta.star.samples")); 
    }
    if (family == 1) {
      SET_VECTOR_ELT(resultName_r, tmp_02, mkChar("kappa.samples")); 
    }
   
    namesgets(result_r, resultName_r);
    
    UNPROTECT(nProtect);
    
    return(result_r);
  }
}

