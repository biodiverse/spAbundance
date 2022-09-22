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
  SEXP NMix(SEXP y_r, SEXP X_r, SEXP Xp_r, SEXP XRE_r, SEXP XpRE_r, SEXP yMax_r,
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
	    SEXP chainInfo_r){
   
    /**********************************************************************
     * Initial constants
     * *******************************************************************/
    int i, g, t, j, s, r, l, info, nProtect=0;
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
    double *Xp = REAL(Xp_r);
    double *yMax = REAL(yMax_r);
    int *XRE = INTEGER(XRE_r); 
    int *XpRE = INTEGER(XpRE_r);
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
    double *SigmaBetaInv = (double *) R_alloc(ppAbund, sizeof(double));   
    F77_NAME(dcopy)(&ppAbund, REAL(SigmaBeta_r), &inc, SigmaBetaInv, &inc);
    double *SigmaAlphaInv = (double *) R_alloc(ppDet, sizeof(double));   
    F77_NAME(dcopy)(&ppDet, REAL(SigmaAlpha_r), &inc, SigmaAlphaInv, &inc);
    double *sigmaSqMuA = REAL(sigmaSqMuA_r); 
    double *sigmaSqMuB = REAL(sigmaSqMuB_r); 
    double *sigmaSqPA = REAL(sigmaSqPA_r); 
    double *sigmaSqPB = REAL(sigmaSqPB_r); 
    double kappaA = REAL(kappaA_r)[0];
    double kappaB = REAL(kappaB_r)[0];
    int *nAbundRELong = INTEGER(nAbundRELong_r); 
    int *nDetRELong = INTEGER(nDetRELong_r); 
    double *K = REAL(K_r); 
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
    double tuning = REAL(tuning_r)[0];
    int nChain = INTEGER(chainInfo_r)[1];
    int nThreads = INTEGER(nThreads_r)[0];
    int verbose = INTEGER(verbose_r)[0];
    int nReport = INTEGER(nReport_r)[0];
    int status = 0; 
    int thinIndx = 0;
    int sPost = 0;  

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
        Rprintf("Negative Binomial N-mixture model with Polya-Gamma latent\nvariable fit with %i sites.\n\n", J);
        Rprintf("Samples per Chain: %i (%i batches of length %i)\n", nSamples, nBatch, batchLength);
        Rprintf("Burn-in: %i \n", nBurn); 
        Rprintf("Thinning Rate: %i \n", nThin); 
        Rprintf("Number of Chains: %i \n", nChain);
        Rprintf("Total Posterior Samples: %i \n\n", nPost * nChain); 
#ifdef _OPENMP
        Rprintf("\nSource compiled with OpenMP support and model fit using %i thread(s).\n\n", nThreads);
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
    // Auxiliary variables
    double *omegaAbund = (double *) R_alloc(J, sizeof(double)); zeros(omegaAbund, J);
    double *omegaDet = (double *) R_alloc(nObs, sizeof(double)); zeros(omegaDet, nObs);
    double *NStar = (double *) R_alloc(J, sizeof(double)); zeros(NStar, J);
    double *yStar = (double *) R_alloc(nObs, sizeof(double)); zeros(yStar, nObs);

    /**********************************************************************
     * Return Stuff
     * *******************************************************************/
    SEXP betaSamples_r;
    PROTECT(betaSamples_r = allocMatrix(REALSXP, pAbund, nPost)); nProtect++;
    SEXP alphaSamples_r; 
    PROTECT(alphaSamples_r = allocMatrix(REALSXP, pDet, nPost)); nProtect++;
    SEXP NSamples_r; 
    PROTECT(NSamples_r = allocMatrix(REALSXP, J, nPost)); nProtect++; 
    SEXP kappaSamples_r;
    PROTECT(kappaSamples_r = allocMatrix(REALSXP, inc, nPost)); nProtect++;
    SEXP muSamples_r; 
    PROTECT(muSamples_r = allocMatrix(REALSXP, J, nPost)); nProtect++; 
    // Detection random effects
    SEXP sigmaSqPSamples_r; 
    SEXP alphaStarSamples_r; 
    if (pDetRE > 0) {
      PROTECT(sigmaSqPSamples_r = allocMatrix(REALSXP, pDetRE, nPost)); nProtect++;
      PROTECT(alphaStarSamples_r = allocMatrix(REALSXP, nDetRE, nPost)); nProtect++;
    }
    // Abundance random effects
    SEXP sigmaSqMuSamples_r; 
    SEXP betaStarSamples_r; 
    if (pAbundRE > 0) {
      PROTECT(sigmaSqMuSamples_r = allocMatrix(REALSXP, pAbundRE, nPost)); nProtect++;
      PROTECT(betaStarSamples_r = allocMatrix(REALSXP, nAbundRE, nPost)); nProtect++;
    }
    
    /********************************************************************
      Some constants and temporary variables to be used later
    ********************************************************************/
    int JpAbund = J * pAbund; 
    int nObspDet = nObs * pDet;
    double tmp_0; 
    double *tmp_one = (double *) R_alloc(inc, sizeof(double)); 
    double *tmp_ppDet = (double *) R_alloc(ppDet, sizeof(double));
    double *tmp_ppAbund = (double *) R_alloc(ppAbund, sizeof(double)); 
    double *tmp_pDet = (double *) R_alloc(pDet, sizeof(double));
    double *tmp_pAbund = (double *) R_alloc(pAbund, sizeof(double));
    double *tmp_pDet2 = (double *) R_alloc(pDet, sizeof(double));
    double *tmp_pAbund2 = (double *) R_alloc(pAbund, sizeof(double));
    int *tmp_J = (int *) R_alloc(J, sizeof(int));
    for (j = 0; j < J; j++) {
      tmp_J[j] = 0; 
    }
    double *tmp_nObs = (double *) R_alloc(nObs, sizeof(double)); 
    double *tmp_JpAbund = (double *) R_alloc(JpAbund, sizeof(double));
    double *tmp_nObspDet = (double *) R_alloc(nObspDet, sizeof(double));
    double *tmp_J1 = (double *) R_alloc(J, sizeof(double));
   
    // For latent abundance and WAIC
    double muNum; 
    double *detProb = (double *) R_alloc(nObs, sizeof(double)); zeros(detProb, nObs);
    double *psi = (double *) R_alloc(J, sizeof(double)); 
    zeros(psi, J); 
    double *mu = (double *) R_alloc(J, sizeof(double)); 
    zeros(mu, J); 
    double *piProd = (double *) R_alloc(J, sizeof(double)); 
    ones(piProd, J); 
    double *piProdWAIC = (double *) R_alloc(J, sizeof(double)); 
    ones(piProdWAIC, J); 
    double *ySum = (double *) R_alloc(J, sizeof(double)); zeros(ySum, J);

    // For normal priors
    // Abundupancy regression coefficient priors. 
    F77_NAME(dpotrf)(lower, &pAbund, SigmaBetaInv, &pAbund, &info FCONE); 
    if(info != 0){error("c++ error: dpotrf SigmaBetaInv failed\n");}
    F77_NAME(dpotri)(lower, &pAbund, SigmaBetaInv, &pAbund, &info FCONE); 
    if(info != 0){error("c++ error: dpotri SigmaBetaInv failed\n");}
    double *SigmaBetaInvMuBeta = (double *) R_alloc(pAbund, sizeof(double)); 
    F77_NAME(dsymv)(lower, &pAbund, &one, SigmaBetaInv, &pAbund, muBeta, &inc, &zero, 
        	    SigmaBetaInvMuBeta, &inc FCONE);
    // Detection regression coefficient priors. 
    F77_NAME(dpotrf)(lower, &pDet, SigmaAlphaInv, &pDet, &info FCONE); 
    if(info != 0){error("c++ error: dpotrf SigmaAlphaInv failed\n");}
    F77_NAME(dpotri)(lower, &pDet, SigmaAlphaInv, &pDet, &info FCONE); 
    if(info != 0){error("c++ error: dpotri SigmaAlphaInv failed\n");}
    double *SigmaAlphaInvMuAlpha = (double *) R_alloc(pDet, sizeof(double)); 
    F77_NAME(dsymv)(lower, &pDet, &one, SigmaAlphaInv, &pDet, muAlpha, &inc, &zero, 
                   SigmaAlphaInvMuAlpha, &inc FCONE);

    // For NB PG sampling (recommendations from Polson et al. 2013). 
    int trunc = 200;
    double *tmp_trunc = (double *) R_alloc(trunc, sizeof(double));
    
    /********************************************************************
      Set up MH stuff
    ********************************************************************/
    double accept = 0.0;
    double logPostCurr = 0.0, logPostCand = 0.0;
    double *logPostCurrN = (double *) R_alloc(J, sizeof(double));
    zeros(logPostCurrN, J);
    double *logPostCandN = (double *) R_alloc(J, sizeof(double));
    zeros(logPostCandN, J);
    double epsilon = 1.0;
    double kappaCand = 0.0;
    double *NCand = (double *) R_alloc(J, sizeof(double));
    zeros(NCand, J);
    SEXP acceptSamples_r; 
    PROTECT(acceptSamples_r = allocMatrix(REALSXP, inc, nBatch)); nProtect++; 
    SEXP tuningSamples_r; 
    PROTECT(tuningSamples_r = allocMatrix(REALSXP, inc, nBatch)); nProtect++; 

    /**********************************************************************
     * Prep for random effects
     * *******************************************************************/
    // Site-level sums of the abundance random effects
    double *betaStarSites = (double *) R_alloc(J, sizeof(double)); 
    zeros(betaStarSites, J); 
    // Initial sums
    for (j = 0; j < J; j++) {
      for (l = 0; l < pAbundRE; l++) {
        betaStarSites[j] += betaStar[which(XRE[l * J + j], betaLevelIndx, nAbundRE)];
      }
    }
    // Observation-level sums of the detection random effects
    double *alphaStarObs = (double *) R_alloc(nObs, sizeof(double)); 
    zeros(alphaStarObs, nObs); 
    // Get sums of the current REs for each site/visit combo
    for (i = 0; i < nObs; i++) {
      for (l = 0; l < pDetRE; l++) {
        alphaStarObs[i] += alphaStar[which(XpRE[l * nObs + i], alphaLevelIndx, nDetRE)];
      }
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
         *Update Abundance Auxiliary Variables 
         *******************************************************************/
        for (j = 0; j < J; j++) {
          omegaAbund[j] = rpgGamma(kappa + N[j], F77_NAME(ddot)(&pAbund, &X[j], &J, beta, &inc) + betaStarSites[j], trunc, tmp_trunc);
        } // j
      
        /********************************************************************
        *Update Detection Auxiliary Variables 
        *******************************************************************/
	zeros(omegaDet, nObs);
        for (i = 0; i < nObs; i++) {
          if (N[NLongIndx[i]] > 0.0) {
            omegaDet[i] = rpg(N[NLongIndx[i]], 
			      F77_NAME(ddot)(&pDet, &Xp[i], &nObs, alpha, &inc) + 
        		      alphaStarObs[i]);
          }
        } // i

        /********************************************************************
         *Update Abundance Regression Coefficients
         *******************************************************************/
        for (j = 0; j < J; j++) {
          NStar[j] = (N[j] - kappa) / (2.0 * omegaAbund[j]);
          tmp_J1[j] = (NStar[j] - betaStarSites[j]) * omegaAbund[j];
        } // j

        /********************************
         * Compute b.beta
         *******************************/
        F77_NAME(dgemv)(ytran, &J, &pAbund, &one, X, &J, tmp_J1, &inc, &zero, tmp_pAbund, &inc FCONE); 	 
        for (j = 0; j < pAbund; j++) {
          tmp_pAbund[j] += SigmaBetaInvMuBeta[j]; 
        } // j 

        /********************************
         * Compute A.beta
         * *****************************/
        // tmp_Jp is X %*% omega. 
        for(j = 0; j < J; j++){
          for(i = 0; i < pAbund; i++){
            tmp_JpAbund[i*J+j] = X[i*J+j]*omegaAbund[j];
          }
        }

        // This finishes off A.beta
        // 1 * X * tmp_Jp + 0 * tmp_pp = tmp_pp
        F77_NAME(dgemm)(ytran, ntran, &pAbund, &pAbund, &J, &one, X, &J, 
			tmp_JpAbund, &J, &zero, tmp_ppAbund, &pAbund FCONE FCONE);
        for (j = 0; j < ppAbund; j++) {
          tmp_ppAbund[j] += SigmaBetaInv[j]; 
        } // j

        F77_NAME(dpotrf)(lower, &pAbund, tmp_ppAbund, &pAbund, &info FCONE); 
        if(info != 0){error("c++ error: dpotrf here failed\n");}
        F77_NAME(dpotri)(lower, &pAbund, tmp_ppAbund, &pAbund, &info FCONE); 
        if(info != 0){error("c++ error: dpotri here failed\n");}
        F77_NAME(dsymv)(lower, &pAbund, &one, tmp_ppAbund, &pAbund, tmp_pAbund,
		       	&inc, &zero, tmp_pAbund2, &inc FCONE);
        F77_NAME(dpotrf)(lower, &pAbund, tmp_ppAbund, &pAbund, &info FCONE); 
        if(info != 0){error("c++ error: dpotrf here failed\n");}
        mvrnorm(beta, tmp_pAbund2, tmp_ppAbund, pAbund);
        
        /********************************************************************
         *Update Detection Regression Coefficients
         *******************************************************************/
        /********************************
         * Compute b.alpha
         *******************************/
        zeros(tmp_nObs, nObs);
	zeros(yStar, nObs);
	zeros(tmp_nObspDet, nObspDet);
        for (i = 0; i < nObs; i++) {
          if (N[NLongIndx[i]] > 0.0) {
            yStar[i] = (y[i] - N[NLongIndx[i]]/2.0) / omegaDet[i];
            tmp_nObs[i] = (yStar[i] - alphaStarObs[i]) * omegaDet[i]; 
	  }
        } // i
        
        F77_NAME(dgemv)(ytran, &nObs, &pDet, &one, Xp, &nObs, tmp_nObs, &inc, &zero, tmp_pDet, &inc FCONE); 	  
        for (j = 0; j < pDet; j++) {
          tmp_pDet[j] += SigmaAlphaInvMuAlpha[j]; 
        } // j

        /********************************
         * Compute A.alpha
         * *****************************/
        for (j = 0; j < nObs; j++) {
          if (N[NLongIndx[j]] > 0.0) {
            for (i = 0; i < pDet; i++) {
              tmp_nObspDet[i*nObs + j] = Xp[i * nObs + j] * omegaDet[j];
            } // i
	  }
        } // j

        F77_NAME(dgemm)(ytran, ntran, &pDet, &pDet, &nObs, &one, Xp, &nObs, tmp_nObspDet, &nObs, &zero, tmp_ppDet, &pDet FCONE FCONE);

        for (j = 0; j < ppDet; j++) {
          tmp_ppDet[j] += SigmaAlphaInv[j]; 
        } // j

        F77_NAME(dpotrf)(lower, &pDet, tmp_ppDet, &pDet, &info FCONE); 
        if(info != 0){error("c++ error: dpotrf A.alpha failed\n");}
        F77_NAME(dpotri)(lower, &pDet, tmp_ppDet, &pDet, &info FCONE); 
        if(info != 0){error("c++ error: dpotri A.alpha failed\n");}
        F77_NAME(dsymv)(lower, &pDet, &one, tmp_ppDet, &pDet, tmp_pDet, &inc, &zero, tmp_pDet2, &inc FCONE);
        F77_NAME(dpotrf)(lower, &pDet, tmp_ppDet, &pDet, &info FCONE); 
        if(info != 0){error("c++ error: dpotrf here failed\n");}
        mvrnorm(alpha, tmp_pDet2, tmp_ppDet, pDet);
        
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
          // Update each individual random effect one by one. 
          for (l = 0; l < nAbundRE; l++) {
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
                tmp_one[0] += (NStar[j] - F77_NAME(ddot)(&pAbund, &X[j], &J, beta, &inc) - 
          		    betaStarSites[j] + betaStar[l]) * omegaAbund[j];
                tmp_0 += omegaAbund[j];
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
            for (l = 0; l < pAbundRE; l++) {
              betaStarSites[j] += betaStar[which(XRE[l * J + j], betaLevelIndx, nAbundRE)];
            }
          }
        }
        
	/********************************************************************
         *Update detection random effects
         *******************************************************************/
        if (pDetRE > 0) {
          // Update each individual random effect one by one. 
          for (l = 0; l < nDetRE; l++) {
            /********************************
             * Compute b.alpha.star
             *******************************/
            // Only allow information to come from when z[r] == 1 and XpRE == alphaLevelIndx[l]
            zeros(tmp_one, inc);
            tmp_0 = 0.0;
            for (r = 0; r < nObs; r++) {
              if ((N[NLongIndx[r]] > 0.0) && (XpRE[alphaStarIndx[l] * nObs + r] == alphaLevelIndx[l])) {
                tmp_one[0] += (yStar[r] - (F77_NAME(ddot)(&pDet, &Xp[r], &nObs, alpha, &inc) + alphaStarObs[r] - alphaStar[l])) * omegaDet[r];
        	      tmp_0 += omegaDet[r];
              }
            }
            /********************************
             * Compute A.alpha.star
             *******************************/
            tmp_0 += 1.0 / sigmaSqP[alphaStarIndx[l]]; 
            tmp_0 = 1.0 / tmp_0; 
            alphaStar[l] = rnorm(tmp_0 * tmp_one[0], sqrt(tmp_0)); 
          }
          zeros(alphaStarObs, nObs); 
          // Update the RE sums for the current species
          for (r = 0; r < nObs; r++) {
            for (l = 0; l < pDetRE; l++) {
              alphaStarObs[r] += alphaStar[which(XpRE[l * nObs + r], alphaLevelIndx, nDetRE)]; 
            }
          }
        }
     
        /********************************************************************
         *Update kappa (the NB size parameter)
         *******************************************************************/
        for (j = 0; j < J; j++) {
          tmp_J1[j] = F77_NAME(ddot)(&pAbund, &X[j], &J, beta, &inc) + betaStarSites[j];
          psi[j] = logitInv(tmp_J1[j], zero, one);
        }
        /********************************
         * Current
         *******************************/
        // Log gamma function is lgammafn
        // Likelihood contribution
        logPostCurr = 0.0;
        for (j = 0; j < J; j++) {
          logPostCurr += lgammafn(N[j] + kappa) - lgammafn(kappa) + kappa * log(1 - psi[j]);
        }
        logPostCurr += log(kappa - kappaA) + log(kappaB - kappa);
        /********************************
         * Candidate
         *******************************/
        kappaCand = logitInv(rnorm(logit(kappa, kappaA, kappaB), exp(tuning)), kappaA, kappaB);
        logPostCand = 0.0;
        for (j = 0; j < J; j++) {
          logPostCand += lgammafn(N[j] + kappaCand) - lgammafn(kappaCand) + kappaCand * log(1 - psi[j]);
        }
        logPostCand += log(kappaCand - kappaA) + log(kappaB - kappaCand);

        if (runif(0.0,1.0) <= exp(logPostCand - logPostCurr)) {
          kappa = kappaCand;
          accept++;
        }

        /********************************************************************
         *Update Latent Abundance 
         *******************************************************************/
	zeros(logPostCurrN, J);
	zeros(logPostCandN, J);
	// Proposal
	for (j = 0; j < J; j++) {
          NCand[j] = rpois(N[j] + epsilon);
	  mu[j] = exp(tmp_J1[j]) * kappa;
	  // Rprintf("NCand[%i]: %f\n", j, NCand[j]); 
	}
	// Likelihood contribution to Metropolis ratios
	for (i = 0; i < nObs; i++) {
          detProb[i] = logitInv(F77_NAME(ddot)(&pDet, &Xp[i], &nObs, alpha, &inc) + 
			        alphaStarObs[i], zero, one);
          logPostCurrN[NLongIndx[i]] += dbinom(y[i], N[NLongIndx[i]], detProb[i], 1);
	  logPostCandN[NLongIndx[i]] += dbinom(y[i], NCand[NLongIndx[i]], detProb[i], 1);
	}
        /********************************
         * Current
         *******************************/
        // Log gamma function is lgammafn
        // Likelihood contribution
        for (j = 0; j < J; j++) {
	  if (NCand[j] >= yMax[j]) {
	    // Contribution from NB(N)
	    logPostCurrN[j] += dnbinom_mu(N[j], kappa, mu[j], 1);
	    // MH contribution for assymetric proposal distribution.
	    logPostCurrN[j] += dpois(NCand[j], N[j] + epsilon, 1);
            /********************************
             * Candidate
             *******************************/
            logPostCandN[j] += dnbinom_mu(NCand[j], kappa, mu[j], 1);
	    logPostCandN[j] += dpois(N[j], NCand[j] + epsilon, 1);
	    // Rprintf("logPostCurr[%i]: %f\n", j, logPostCurrN[j]);
	    // Rprintf("logPostCand[%i]: %f\n", j, logPostCandN[j]);
	    // Rprintf("N[%i]: %f\n", j, N[j]);
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
	    REAL(kappaSamples_r)[sPost] = kappa;
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
      REAL(acceptSamples_r)[s] = accept/batchLength; 
      REAL(tuningSamples_r)[s] = tuning; 
      if (accept / batchLength > acceptRate) {
        tuning += std::min(0.01, 1.0/sqrt(static_cast<double>(s)));
      } else{
          tuning -= std::min(0.01, 1.0/sqrt(static_cast<double>(s)));
        }
      accept = 0.0;

      /********************************************************************
       *Report 
       *******************************************************************/
      if (verbose) {
        if (status == nReport) {
          Rprintf("Batch: %i of %i, %3.2f%%\n", s, nBatch, 100.0*s/nBatch);
          Rprintf("\tParameter\tAcceptance\tTuning\n");	  
          Rprintf("\tkappa\t\t%3.1f\t\t%1.5f\n", 100.0*REAL(acceptSamples_r)[s], exp(tuning));
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
    int nResultListObjs = 6;
    if (pDetRE > 0) {
      nResultListObjs += 2; 
    }
    if (pAbundRE > 0) {
      nResultListObjs += 2;
    }

    PROTECT(result_r = allocVector(VECSXP, nResultListObjs)); nProtect++;
    PROTECT(resultName_r = allocVector(VECSXP, nResultListObjs)); nProtect++;

    // Setting the components of the output list.
    SET_VECTOR_ELT(result_r, 0, betaSamples_r);
    SET_VECTOR_ELT(result_r, 1, alphaSamples_r);
    SET_VECTOR_ELT(result_r, 2, kappaSamples_r);
    SET_VECTOR_ELT(result_r, 3, NSamples_r); 
    SET_VECTOR_ELT(result_r, 4, muSamples_r);
    if (pDetRE > 0) {
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
    SET_VECTOR_ELT(resultName_r, 0, mkChar("beta.samples")); 
    SET_VECTOR_ELT(resultName_r, 1, mkChar("alpha.samples")); 
    SET_VECTOR_ELT(resultName_r, 2, mkChar("kappa.samples")); 
    SET_VECTOR_ELT(resultName_r, 3, mkChar("N.samples")); 
    SET_VECTOR_ELT(resultName_r, 4, mkChar("mu.samples"));
    if (pDetRE > 0) {
      SET_VECTOR_ELT(resultName_r, 5, mkChar("sigma.sq.p.samples")); 
      SET_VECTOR_ELT(resultName_r, 6, mkChar("alpha.star.samples")); 
    }
    if (pAbundRE > 0) {
      SET_VECTOR_ELT(resultName_r, tmp_0, mkChar("sigma.sq.mu.samples")); 
      SET_VECTOR_ELT(resultName_r, tmp_0 + 1, mkChar("beta.star.samples")); 
    }
   
    namesgets(result_r, resultName_r);
    
    UNPROTECT(nProtect);
    
    return(result_r);
  }
}

