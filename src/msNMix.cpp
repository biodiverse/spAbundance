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
  SEXP msNMix(SEXP y_r, SEXP X_r, SEXP Xp_r, SEXP XRE_r, SEXP XpRE_r, SEXP yMax_r, 
	      SEXP consts_r, SEXP nAbundRELong_r, SEXP nDetRELong_r, 
	      SEXP betaStarting_r, SEXP alphaStarting_r, SEXP kappaStarting_r, SEXP NStarting_r, 
	      SEXP betaCommStarting_r, SEXP alphaCommStarting_r, 
	      SEXP tauSqBetaStarting_r, SEXP tauSqAlphaStarting_r, 
	      SEXP sigmaSqMuStarting_r, SEXP sigmaSqPStarting_r, 
	      SEXP betaStarStarting_r, SEXP alphaStarStarting_r, 
	      SEXP NLongIndx_r, SEXP betaStarIndx_r, SEXP betaLevelIndx_r, 
	      SEXP alphaStarIndx_r, SEXP alphaLevelIndx_r, 
	      SEXP muBetaComm_r, SEXP muAlphaComm_r, 
	      SEXP SigmaBetaComm_r, SEXP SigmaAlphaComm_r, SEXP kappaA_r, 
	      SEXP kappaB_r, SEXP tauSqBetaA_r, 
	      SEXP tauSqBetaB_r, SEXP tauSqAlphaA_r, SEXP tauSqAlphaB_r, 
	      SEXP sigmaSqMuA_r, SEXP sigmaSqMuB_r, 
	      SEXP sigmaSqPA_r, SEXP sigmaSqPB_r, 
	      SEXP tuning_r, SEXP nBatch_r, SEXP batchLength_r, SEXP acceptRate_r,
	      SEXP nThreads_r, SEXP verbose_r, SEXP nReport_r, 
	      SEXP samplesInfo_r, SEXP chainInfo_r){
   
    /**********************************************************************
     * Initial constants
     * *******************************************************************/
    int i, j, s, t, g, l, q, r, info, nProtect=0;
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
    int *XRE = INTEGER(XRE_r); 
    int *XpRE = INTEGER(XpRE_r); 
    double *yMax = REAL(yMax_r);
    // Load constants
    int nSp = INTEGER(consts_r)[0]; 
    int J = INTEGER(consts_r)[1];
    int nObs = INTEGER(consts_r)[2]; 
    int pAbund = INTEGER(consts_r)[3];
    int pAbundRE = INTEGER(consts_r)[4];
    int nAbundRE = INTEGER(consts_r)[5];
    int pDet = INTEGER(consts_r)[6];
    int pDetRE = INTEGER(consts_r)[7];
    int nDetRE = INTEGER(consts_r)[8];
    int ppDet = pDet * pDet;
    int ppAbund = pAbund * pAbund; 
    double *muBetaComm = REAL(muBetaComm_r); 
    double *muAlphaComm = REAL(muAlphaComm_r); 
    double *SigmaBetaCommInv = (double *) R_alloc(ppAbund, sizeof(double));   
    F77_NAME(dcopy)(&ppAbund, REAL(SigmaBetaComm_r), &inc, SigmaBetaCommInv, &inc);
    double *SigmaAlphaCommInv = (double *) R_alloc(ppDet, sizeof(double));   
    F77_NAME(dcopy)(&ppDet, REAL(SigmaAlphaComm_r), &inc, SigmaAlphaCommInv, &inc);
    double *tauSqBetaA = REAL(tauSqBetaA_r); 
    double *tauSqBetaB = REAL(tauSqBetaB_r); 
    double *tauSqAlphaA = REAL(tauSqAlphaA_r); 
    double *tauSqAlphaB = REAL(tauSqAlphaB_r); 
    double *sigmaSqMuA = REAL(sigmaSqMuA_r); 
    double *sigmaSqMuB = REAL(sigmaSqMuB_r); 
    double *sigmaSqPA = REAL(sigmaSqPA_r); 
    double *sigmaSqPB = REAL(sigmaSqPB_r); 
    double *kappaA = REAL(kappaA_r); 
    double *kappaB = REAL(kappaB_r); 
    int *nAbundRELong = INTEGER(nAbundRELong_r); 
    int *nDetRELong = INTEGER(nDetRELong_r); 
    int *NLongIndx = INTEGER(NLongIndx_r); 
    int *betaStarIndx = INTEGER(betaStarIndx_r);
    int *betaLevelIndx = INTEGER(betaLevelIndx_r);
    int *alphaStarIndx = INTEGER(alphaStarIndx_r); 
    int *alphaLevelIndx = INTEGER(alphaLevelIndx_r);
    int nBatch = INTEGER(nBatch_r)[0]; 
    int batchLength = INTEGER(batchLength_r)[0]; 
    double acceptRate = REAL(acceptRate_r)[0];
    int nSamples = nBatch * batchLength; 
    double *tuning = REAL(tuning_r);
    int nThreads = INTEGER(nThreads_r)[0];
    int verbose = INTEGER(verbose_r)[0];
    int nReport = INTEGER(nReport_r)[0];
    int nBurn = INTEGER(samplesInfo_r)[0]; 
    int nThin = INTEGER(samplesInfo_r)[1];
    int nPost = INTEGER(samplesInfo_r)[2]; 
    int currChain = INTEGER(chainInfo_r)[0];
    int nChain = INTEGER(chainInfo_r)[1];
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
        Rprintf("Multispecies Negative Binomial N-Mixture Model with Polya-Gamma latent\nvariable fit with %i sites and %i species.\n\n", J, nSp);
        Rprintf("Samples per Chain: %i \n", nSamples);
        Rprintf("Burn-in: %i \n", nBurn); 
        Rprintf("Thinning Rate: %i \n", nThin); 
	Rprintf("Number of Chains: %i \n", nChain);
        Rprintf("Total Posterior Samples: %i \n\n", nPost * nChain); 
#ifdef _OPENMP
        Rprintf("\nSource compiled with OpenMP support and model fit using %i thread(s).\n\n", nThreads);
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
       Some constants and temporary variables to be used later
     * *******************************************************************/
    int pAbundnSp = pAbund * nSp; 
    int pDetnSp = pDet * nSp; 
    int nObsnSp = nObs * nSp; 
    int nAbundREnSp = nAbundRE * nSp; 
    int nDetREnSp = nDetRE * nSp; 
    int JnSp = J * nSp;
    int JpAbund = J * pAbund; 
    int nObspDet = nObs * pDet;
    double tmp_0; 
    double *tmp_one = (double *) R_alloc(inc, sizeof(double)); 
    double *tmp_ppDet = (double *) R_alloc(ppDet, sizeof(double));
    double *tmp_ppAbund = (double *) R_alloc(ppAbund, sizeof(double)); 
    double *tmp_pDet = (double *) R_alloc(pDet, sizeof(double));
    double *tmp_pAbund = (double *) R_alloc(pAbund, sizeof(double));
    double *tmp_beta = (double *) R_alloc(pAbund, sizeof(double));
    double *tmp_alpha = (double *) R_alloc(pDet, sizeof(double));
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

    /**********************************************************************
     * Parameters
     * *******************************************************************/
    double *betaComm = (double *) R_alloc(pAbund, sizeof(double)); 
    F77_NAME(dcopy)(&pAbund, REAL(betaCommStarting_r), &inc, betaComm, &inc);
    double *tauSqBeta = (double *) R_alloc(pAbund, sizeof(double)); 
    F77_NAME(dcopy)(&pAbund, REAL(tauSqBetaStarting_r), &inc, tauSqBeta, &inc);
    double *alphaComm = (double *) R_alloc(pDet, sizeof(double));   
    F77_NAME(dcopy)(&pDet, REAL(alphaCommStarting_r), &inc, alphaComm, &inc);
    double *tauSqAlpha = (double *) R_alloc(pDet, sizeof(double)); 
    F77_NAME(dcopy)(&pDet, REAL(tauSqAlphaStarting_r), &inc, tauSqAlpha, &inc);
    double *beta = (double *) R_alloc(pAbundnSp, sizeof(double));   
    F77_NAME(dcopy)(&pAbundnSp, REAL(betaStarting_r), &inc, beta, &inc);
    double *alpha = (double *) R_alloc(pDetnSp, sizeof(double));   
    F77_NAME(dcopy)(&pDetnSp, REAL(alphaStarting_r), &inc, alpha, &inc);
    // Abundance random effect variances
    double *sigmaSqMu = (double *) R_alloc(pAbundRE, sizeof(double)); 
    F77_NAME(dcopy)(&pAbundRE, REAL(sigmaSqMuStarting_r), &inc, sigmaSqMu, &inc); 
    // Latent occupancy random effects
    double *betaStar = (double *) R_alloc(nAbundREnSp, sizeof(double)); 
    F77_NAME(dcopy)(&nAbundREnSp, REAL(betaStarStarting_r), &inc, betaStar, &inc); 
    // Detection random effect variances
    double *sigmaSqP = (double *) R_alloc(pDetRE, sizeof(double)); 
    F77_NAME(dcopy)(&pDetRE, REAL(sigmaSqPStarting_r), &inc, sigmaSqP, &inc); 
    // Latent detection random effects
    double *alphaStar = (double *) R_alloc(nDetREnSp, sizeof(double)); 
    F77_NAME(dcopy)(&nDetREnSp, REAL(alphaStarStarting_r), &inc, alphaStar, &inc); 
    // Overdispersion parameter
    double *kappa = (double *) R_alloc(nSp, sizeof(double)); 
    F77_NAME(dcopy)(&nSp, REAL(kappaStarting_r), &inc, kappa, &inc); 
    // Latent Abundance
    double *N = (double *) R_alloc(JnSp, sizeof(double));   
    F77_NAME(dcopy)(&JnSp, REAL(NStarting_r), &inc, N, &inc);
    // Auxiliary variables
    double *omegaDet = (double *) R_alloc(nObs, sizeof(double)); zeros(omegaDet, nObs);
    double *omegaAbund = (double *) R_alloc(JnSp, sizeof(double)); zeros(omegaAbund, JnSp);
    double *NStar = (double *) R_alloc(JnSp, sizeof(double)); zeros(NStar, JnSp);
    double *yStar = (double *) R_alloc(nObs, sizeof(double)); zeros(yStar, nObs);

    /**********************************************************************
     * Return Stuff
     * *******************************************************************/
    // Community level
    SEXP betaCommSamples_r; 
    PROTECT(betaCommSamples_r = allocMatrix(REALSXP, pAbund, nPost)); nProtect++;
    SEXP alphaCommSamples_r;
    PROTECT(alphaCommSamples_r = allocMatrix(REALSXP, pDet, nPost)); nProtect++;
    SEXP tauSqBetaSamples_r; 
    PROTECT(tauSqBetaSamples_r = allocMatrix(REALSXP, pAbund, nPost)); nProtect++; 
    SEXP tauSqAlphaSamples_r; 
    PROTECT(tauSqAlphaSamples_r = allocMatrix(REALSXP, pDet, nPost)); nProtect++; 
    // Species level
    SEXP betaSamples_r;
    PROTECT(betaSamples_r = allocMatrix(REALSXP, pAbundnSp, nPost)); nProtect++;
    SEXP alphaSamples_r; 
    PROTECT(alphaSamples_r = allocMatrix(REALSXP, pDetnSp, nPost)); nProtect++;
    SEXP NSamples_r; 
    PROTECT(NSamples_r = allocMatrix(REALSXP, JnSp, nPost)); nProtect++; 
    SEXP muSamples_r; 
    PROTECT(muSamples_r = allocMatrix(REALSXP, JnSp, nPost)); nProtect++; 
    // Detection random effects
    SEXP sigmaSqPSamples_r; 
    SEXP alphaStarSamples_r; 
    if (pDetRE > 0) {
      PROTECT(sigmaSqPSamples_r = allocMatrix(REALSXP, pDetRE, nPost)); nProtect++;
      PROTECT(alphaStarSamples_r = allocMatrix(REALSXP, nDetREnSp, nPost)); nProtect++;
    }
    // Abundance random effects
    SEXP sigmaSqMuSamples_r; 
    SEXP betaStarSamples_r; 
    if (pAbundRE > 0) {
      PROTECT(sigmaSqMuSamples_r = allocMatrix(REALSXP, pAbundRE, nPost)); nProtect++;
      PROTECT(betaStarSamples_r = allocMatrix(REALSXP, nAbundREnSp, nPost)); nProtect++;
    }
    SEXP kappaSamples_r;
    PROTECT(kappaSamples_r = allocMatrix(REALSXP, nSp, nPost)); nProtect++;
    // Likelihood samples for WAIC. 
    
    /**********************************************************************
     * Additional Sampler Prep
     * *******************************************************************/
    // For latent occupancy
    double muNum; 
    double *detProb = (double *) R_alloc(nObsnSp, sizeof(double)); 
    double *psi = (double *) R_alloc(J, sizeof(double)); 
    zeros(psi, J); 
    double *mu = (double *) R_alloc(JnSp, sizeof(double)); 
    zeros(mu, JnSp); 
    double *piProd = (double *) R_alloc(JnSp, sizeof(double)); 
    ones(piProd, JnSp); 
    double *piProdWAIC = (double *) R_alloc(JnSp, sizeof(double)); 
    ones(piProdWAIC, JnSp); 
    double *ySum = (double *) R_alloc(JnSp, sizeof(double)); zeros(ySum, nSp);

    // For normal priors
    F77_NAME(dpotrf)(lower, &pAbund, SigmaBetaCommInv, &pAbund, &info FCONE); 
    if(info != 0){error("c++ error: dpotrf SigmaBetaCommInv failed\n");}
    F77_NAME(dpotri)(lower, &pAbund, SigmaBetaCommInv, &pAbund, &info FCONE); 
    if(info != 0){error("c++ error: dpotri SigmaBetaCommInv failed\n");}
    double *SigmaBetaCommInvMuBeta = (double *) R_alloc(pAbund, sizeof(double)); 
    F77_NAME(dsymv)(lower, &pAbund, &one, SigmaBetaCommInv, &pAbund, muBetaComm, &inc, &zero, 
        	    SigmaBetaCommInvMuBeta, &inc FCONE);
    // Detection regression coefficient priors. 
    F77_NAME(dpotrf)(lower, &pDet, SigmaAlphaCommInv, &pDet, &info FCONE); 
    if(info != 0){error("c++ error: dpotrf SigmaAlphaCommInv failed\n");}
    F77_NAME(dpotri)(lower, &pDet, SigmaAlphaCommInv, &pDet, &info FCONE); 
    if(info != 0){error("c++ error: dpotri SigmaAlphaCommInv failed\n");}
    double *SigmaAlphaCommInvMuAlpha = (double *) R_alloc(pDet, sizeof(double)); 
    F77_NAME(dsymv)(lower, &pDet, &one, SigmaAlphaCommInv, &pDet, muAlphaComm, &inc, &zero, 
                   SigmaAlphaCommInvMuAlpha, &inc FCONE);
    // Put community level variances in a pAbund x PAbund matrix.
    double *TauBetaInv = (double *) R_alloc(ppAbund, sizeof(double)); zeros(TauBetaInv, ppAbund); 
    for (i = 0; i < pAbund; i++) {
      TauBetaInv[i * pAbund + i] = tauSqBeta[i]; 
    } // i
    F77_NAME(dpotrf)(lower, &pAbund, TauBetaInv, &pAbund, &info FCONE); 
    if(info != 0){error("c++ error: dpotrf TauBetaInv failed\n");}
    F77_NAME(dpotri)(lower, &pAbund, TauBetaInv, &pAbund, &info FCONE); 
    if(info != 0){error("c++ error: dpotri TauBetaInv failed\n");}
    // Put community level variances in a pDet x pDet matrix. 
    double *TauAlphaInv = (double *) R_alloc(ppDet, sizeof(double)); zeros(TauAlphaInv, ppDet); 
    for (i = 0; i < pDet; i++) {
      TauAlphaInv[i * pDet + i] = tauSqAlpha[i]; 
    } // i
    F77_NAME(dpotrf)(lower, &pDet, TauAlphaInv, &pDet, &info FCONE); 
    if(info != 0){error("c++ error: dpotrf TauAlphaInv failed\n");}
    F77_NAME(dpotri)(lower, &pDet, TauAlphaInv, &pDet, &info FCONE); 
    if(info != 0){error("c++ error: dpotri TauAlphaInv failed\n");}

    // For NB PG sampling (recommendations from Polson et al. 2013). 
    int trunc = 200;
    double *tmp_trunc = (double *) R_alloc(trunc, sizeof(double));

    /**********************************************************************
     * Prep for random effects
     * *******************************************************************/
    // Site-level sums of the occurrence random effects
    double *betaStarSites = (double *) R_alloc(JnSp, sizeof(double)); 
    zeros(betaStarSites, JnSp); 
    // Initial sums (initiate with the first species)
    for (i = 0; i < nSp; i++) {
      for (j = 0; j < J; j++) {
        for (l = 0; l < pAbundRE; l++) {
          betaStarSites[i * J + j] += betaStar[i * nAbundRE + which(XRE[l * J + j], betaLevelIndx, nAbundRE)];
        }
      }
    }
    // Observation-level sums of the detection random effects
    double *alphaStarObs = (double *) R_alloc(nObsnSp, sizeof(double)); 
    zeros(alphaStarObs, nObsnSp); 
    // Get sums of the current REs for each site/visit combo for all species
    for (i = 0; i < nSp; i++) {
      for (r = 0; r < nObs; r++) {
        for (l = 0; l < pDetRE; l++) {
          alphaStarObs[i * nObs + r] += alphaStar[i * nDetRE + which(XpRE[l * nObs + r], alphaLevelIndx, nDetRE)];
        }
      }
    }
    // Starting index for occurrence random effects
    int *betaStarStart = (int *) R_alloc(pAbundRE, sizeof(int)); 
    for (l = 0; l < pAbundRE; l++) {
      betaStarStart[l] = which(l, betaStarIndx, nAbundRE); 
    }
    // Starting index for detection random effects
    int *alphaStarStart = (int *) R_alloc(pDetRE, sizeof(int)); 
    for (l = 0; l < pDetRE; l++) {
      alphaStarStart[l] = which(l, alphaStarIndx, nDetRE); 
    }

    /********************************************************************
      Set up MH stuff
    ********************************************************************/
    double logPostCurr = 0.0, logPostCand = 0.0, logDet;
    double kappaCand = 0.0;  
    double *logPostCurrN = (double *) R_alloc(J, sizeof(double));
    zeros(logPostCurrN, J);
    double *logPostCandN = (double *) R_alloc(J, sizeof(double));
    zeros(logPostCandN, J);
    double epsilon = 1.0;
    double *NCand = (double *) R_alloc(J, sizeof(double));
    zeros(NCand, J);
    int nAMCMC = nSp;
    int kappaAMCMCIndx = 0;
    double *accept = (double *) R_alloc(nAMCMC, sizeof(double)); zeros(accept, nAMCMC); 
    SEXP acceptSamples_r; 
    PROTECT(acceptSamples_r = allocMatrix(REALSXP, nAMCMC, nBatch)); nProtect++; 
    SEXP tuningSamples_r; 
    PROTECT(tuningSamples_r = allocMatrix(REALSXP, nAMCMC, nBatch)); nProtect++; 
    
    GetRNGstate();

    for (s = 0, g = 0; s < nBatch; s++) {
      for (t = 0; t < batchLength; t++, g++) {
        /********************************************************************
         Update Community level Abundance Coefficients
         *******************************************************************/
        /********************************
         Compute b.beta.comm
         *******************************/
        zeros(tmp_pAbund, pAbund); 
        for (i = 0; i < nSp; i++) {
          F77_NAME(dgemv)(ytran, &pAbund, &pAbund, &one, TauBetaInv, &pAbund, &beta[i], &nSp, &one, tmp_pAbund, &inc FCONE); 
        } // i
        for (q = 0; q < pAbund; q++) {
          tmp_pAbund[q] += SigmaBetaCommInvMuBeta[q];  
        } // j

        /********************************
         Compute A.beta.comm
         *******************************/
        for (q = 0; q < ppAbund; q++) {
          tmp_ppAbund[q] = SigmaBetaCommInv[q] + nSp * TauBetaInv[q]; 
        }
        F77_NAME(dpotrf)(lower, &pAbund, tmp_ppAbund, &pAbund, &info FCONE); 
        if(info != 0){error("c++ error: dpotrf ABetaComm failed\n");}
        F77_NAME(dpotri)(lower, &pAbund, tmp_ppAbund, &pAbund, &info FCONE); 
        if(info != 0){error("c++ error: dpotri ABetaComm failed\n");}
        F77_NAME(dsymv)(lower, &pAbund, &one, tmp_ppAbund, &pAbund, tmp_pAbund, &inc, &zero, tmp_pAbund2, &inc FCONE);
        F77_NAME(dpotrf)(lower, &pAbund, tmp_ppAbund, &pAbund, &info FCONE); 
        if(info != 0){error("c++ error: dpotrf ABetaComm failed\n");}
        mvrnorm(betaComm, tmp_pAbund2, tmp_ppAbund, pAbund);
        /********************************************************************
         Update Community level Detection Coefficients
         *******************************************************************/
        /********************************
         * Compute b.alpha.comm
         *******************************/
         zeros(tmp_pDet, pDet); 
         for (i = 0; i < nSp; i++) {
           F77_NAME(dgemv)(ytran, &pDet, &pDet, &one, TauAlphaInv, &pDet, &alpha[i], &nSp, &one, tmp_pDet, &inc FCONE); 
         } // i
         for (q = 0; q < pDet; q++) {
           tmp_pDet[q] += SigmaAlphaCommInvMuAlpha[q];  
         } // j
        /********************************
         * Compute A.alpha.comm
         *******************************/
        for (q = 0; q < ppDet; q++) {
          tmp_ppDet[q] = SigmaAlphaCommInv[q] + nSp * TauAlphaInv[q]; 
        }
        F77_NAME(dpotrf)(lower, &pDet, tmp_ppDet, &pDet, &info FCONE); 
        if(info != 0){error("c++ error: dpotrf AAlphaComm failed\n");}
        F77_NAME(dpotri)(lower, &pDet, tmp_ppDet, &pDet, &info FCONE); 
        if(info != 0){error("c++ error: dpotri AAlphaComm failed\n");}
        F77_NAME(dsymv)(lower, &pDet, &one, tmp_ppDet, &pDet, tmp_pDet, &inc, &zero, tmp_pDet2, &inc FCONE);
        F77_NAME(dpotrf)(lower, &pDet, tmp_ppDet, &pDet, &info FCONE); 
        if(info != 0){error("c++ error: dpotrf AAlphaComm failed\n");}
        mvrnorm(alphaComm, tmp_pDet2, tmp_ppDet, pDet);

        /********************************************************************
         Update Community Abundance Variance Parameter
        ********************************************************************/
        for (q = 0; q < pAbund; q++) {
          tmp_0 = 0.0;  
          for (i = 0; i < nSp; i++) {
            tmp_0 += (beta[q * nSp + i] - betaComm[q]) * (beta[q * nSp + i] - betaComm[q]);
          } // i
          tmp_0 *= 0.5;
          tauSqBeta[q] = rigamma(tauSqBetaA[q] + nSp / 2.0, tauSqBetaB[q] + tmp_0); 
        } // q
        for (q = 0; q < pAbund; q++) {
          TauBetaInv[q * pAbund + q] = tauSqBeta[q]; 
        } // q
        F77_NAME(dpotrf)(lower, &pAbund, TauBetaInv, &pAbund, &info FCONE); 
        if(info != 0){error("c++ error: dpotrf TauBetaInv failed\n");}
        F77_NAME(dpotri)(lower, &pAbund, TauBetaInv, &pAbund, &info FCONE); 
        if(info != 0){error("c++ error: dpotri TauBetaInv failed\n");}
        /********************************************************************
         Update Community Detection Variance Parameter
        ********************************************************************/
        for (q = 0; q < pDet; q++) {
          tmp_0 = 0.0;  
          for (i = 0; i < nSp; i++) {
            tmp_0 += (alpha[q * nSp + i] - alphaComm[q]) * (alpha[q * nSp + i] - alphaComm[q]);
          } // i
          tmp_0 *= 0.5;
          tauSqAlpha[q] = rigamma(tauSqAlphaA[q] + nSp / 2.0, tauSqAlphaB[q] + tmp_0); 
        } // q
        for (q = 0; q < pDet; q++) {
          TauAlphaInv[q * pDet + q] = tauSqAlpha[q]; 
        } // q
        F77_NAME(dpotrf)(lower, &pDet, TauAlphaInv, &pDet, &info FCONE); 
        if(info != 0){error("c++ error: dpotrf TauAlphaInv failed\n");}
        F77_NAME(dpotri)(lower, &pDet, TauAlphaInv, &pDet, &info FCONE); 
        if(info != 0){error("c++ error: dpotri TauAlphaInv failed\n");}

        /********************************************************************
         *Update Abundance random effects variance
         *******************************************************************/
        for (l = 0; l < pAbundRE; l++) {
          tmp_0 = 0.0; 
          for (i = 0; i < nSp; i++) {
            tmp_0 += F77_NAME(ddot)(&nAbundRELong[l], &betaStar[i*nAbundRE + betaStarStart[l]], &inc, &betaStar[i*nAbundRE + betaStarStart[l]], &inc); 
          }
          tmp_0 *= 0.5; 
          sigmaSqMu[l] = rigamma(sigmaSqMuA[l] + nAbundRELong[l] * nSp / 2.0, sigmaSqMuB[l] + tmp_0);
        }

        /********************************************************************
         *Update Detection random effects variance
         *******************************************************************/
        for (l = 0; l < pDetRE; l++) {
          tmp_0 = 0.0; 
          for (i = 0; i < nSp; i++) {
            tmp_0 += F77_NAME(ddot)(&nDetRELong[l], &alphaStar[i*nDetRE + alphaStarStart[l]], &inc, &alphaStar[i*nDetRE + alphaStarStart[l]], &inc); 
          }
          tmp_0 *= 0.5; 
          sigmaSqP[l] = rigamma(sigmaSqPA[l] + nDetRELong[l] * nSp / 2.0, sigmaSqPB[l] + tmp_0);
        }

         
        /**********************************************************************
         *Species specific variables  
         *********************************************************************/
        for (i = 0; i < nSp; i++) {  
          /********************************************************************
           *Update Abundance Auxiliary Variables 
           *******************************************************************/
          for (j = 0; j < J; j++) {
            omegaAbund[j * nSp + i] = rpgGamma(kappa[i] + N[j * nSp + i], F77_NAME(ddot)(&pAbund, &X[j], &J, &beta[i], &nSp) + betaStarSites[i * J + j], trunc, tmp_trunc);
          } // j
          /********************************************************************
           *Update Detection Auxiliary Variables 
           *******************************************************************/
	  zeros(omegaDet, nObs);
          for (r = 0; r < nObs; r++) {
            if (N[NLongIndx[r] * nSp + i] > 0.0) {
              omegaDet[r] = rpg(N[NLongIndx[r] * nSp + i], F77_NAME(ddot)(&pDet, &Xp[r], &nObs, &alpha[i], &nSp) + alphaStarObs[i * nObs + r]);
            }
          } // r
             
          /********************************************************************
           *Update Abundance Regression Coefficients
           *******************************************************************/
          for (j = 0; j < J; j++) {
            NStar[j * nSp + i] = (N[j * nSp + i] - kappa[i]) / (2.0 * omegaAbund[j * nSp + i]);
            tmp_J1[j] = (NStar[j * nSp + i] - betaStarSites[i * J + j]) * omegaAbund[j * nSp + i]; 
          } // j
          /********************************
           * Compute b.beta
           *******************************/
          F77_NAME(dgemv)(ytran, &J, &pAbund, &one, X, &J, tmp_J1, &inc, &zero, tmp_pAbund, &inc FCONE);           // TauBetaInv %*% betaComm + tmp_pAbund = tmp_pAbund
          F77_NAME(dgemv)(ntran, &pAbund, &pAbund, &one, TauBetaInv, &pAbund, betaComm, &inc, &one, tmp_pAbund, &inc FCONE); 

          /********************************
           * Compute A.beta
           * *****************************/
          for(j = 0; j < J; j++){
            for(q = 0; q < pAbund; q++){
              tmp_JpAbund[q*J+j] = X[q*J+j] * omegaAbund[j * nSp + i];
            }
          }
          F77_NAME(dgemm)(ytran, ntran, &pAbund, &pAbund, &J, &one, X, &J, tmp_JpAbund, &J, &zero, tmp_ppAbund, &pAbund FCONE FCONE);
          for (q = 0; q < ppAbund; q++) {
            tmp_ppAbund[q] += TauBetaInv[q]; 
          } // q
          F77_NAME(dpotrf)(lower, &pAbund, tmp_ppAbund, &pAbund, &info FCONE); 
          if(info != 0){error("c++ error: dpotrf ABeta failed\n");}
          F77_NAME(dpotri)(lower, &pAbund, tmp_ppAbund, &pAbund, &info FCONE); 
          if(info != 0){error("c++ error: dpotri ABeta failed\n");}
          F77_NAME(dsymv)(lower, &pAbund, &one, tmp_ppAbund, &pAbund, tmp_pAbund, &inc, &zero, tmp_pAbund2, &inc FCONE);
          F77_NAME(dpotrf)(lower, &pAbund, tmp_ppAbund, &pAbund, &info FCONE); 
          if(info != 0){error("c++ error: dpotrf here failed\n");}
          mvrnorm(tmp_beta, tmp_pAbund2, tmp_ppAbund, pAbund);
          for (q = 0; q < pAbund; q++) {
            beta[q * nSp + i] = tmp_beta[q]; 
	    // Rprintf("Species: %i, beta: %f\n", i, beta[q * nSp + i]);
          }
        
          /********************************************************************
           *Update Detection Regression Coefficients
           *******************************************************************/
          /********************************
           * Compute b.alpha
           *******************************/
	  zeros(tmp_nObs, nObs);
	  zeros(yStar, nObs);
	  zeros(tmp_nObspDet, nObspDet);
          for (r = 0; r < nObs; r++) {
            if (N[NLongIndx[r] * nSp + i] > 0.0) {
              yStar[r] = (y[r * nSp + i] - N[NLongIndx[r] * nSp + i] / 2.0) / omegaDet[r];
              tmp_nObs[r] = (yStar[r] - alphaStarObs[i * nObs + r]) * omegaDet[r];
	    }
          } // r
          
          F77_NAME(dgemv)(ytran, &nObs, &pDet, &one, Xp, &nObs, tmp_nObs, &inc, &zero, tmp_pDet, &inc FCONE); 	  
          F77_NAME(dgemv)(ntran, &pDet, &pDet, &one, TauAlphaInv, &pDet, alphaComm, &inc, &one, tmp_pDet, &inc FCONE); 
          /********************************
           * Compute A.alpha
           * *****************************/
          for (r = 0; r < nObs; r++) {
            if (N[NLongIndx[r] * nSp + i] > 0.0) {
              for (q = 0; q < pDet; q++) {
                tmp_nObspDet[q*nObs + r] = Xp[q * nObs + r] * omegaDet[r];
              } // i
	    }
          } // j

          F77_NAME(dgemm)(ytran, ntran, &pDet, &pDet, &nObs, &one, Xp, &nObs, tmp_nObspDet, &nObs, &zero, tmp_ppDet, &pDet FCONE FCONE);

          for (q = 0; q < ppDet; q++) {
            tmp_ppDet[q] += TauAlphaInv[q]; 
          } // q
          F77_NAME(dpotrf)(lower, &pDet, tmp_ppDet, &pDet, &info FCONE); 
          if(info != 0){error("c++ error: dpotrf A.alpha failed\n");}
          F77_NAME(dpotri)(lower, &pDet, tmp_ppDet, &pDet, &info FCONE); 
          if(info != 0){error("c++ error: dpotri A.alpha failed\n");}
          F77_NAME(dsymv)(lower, &pDet, &one, tmp_ppDet, &pDet, tmp_pDet, &inc, &zero, tmp_pDet2, &inc FCONE);
          F77_NAME(dpotrf)(lower, &pDet, tmp_ppDet, &pDet, &info FCONE); 
          if(info != 0){error("c++ error: dpotrf here failed\n");}
          mvrnorm(tmp_alpha, tmp_pDet2, tmp_ppDet, pDet);
          for (q = 0; q < pDet; q++) {
            alpha[q * nSp + i] = tmp_alpha[q];
          }

          /********************************************************************
           *Update Abundance random effects
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
                  tmp_one[0] += (NStar[j * nSp + i] - F77_NAME(ddot)(&pAbund, &X[j], &J, &beta[i], &nSp) - betaStarSites[i * J + j] + betaStar[i * nAbundRE + l]) * omegaAbund[j * nSp + i];
                  tmp_0 += omegaAbund[j * nSp + i];
                }
              }
              /********************************
               * Compute A.beta.star
               *******************************/
              tmp_0 += 1.0 / sigmaSqMu[betaStarIndx[l]]; 
              tmp_0 = 1.0 / tmp_0; 
              betaStar[i * nAbundRE + l] = rnorm(tmp_0 * tmp_one[0], sqrt(tmp_0)); 
            }

            // Update the RE sums for the current species
            zeros(&betaStarSites[i * J], J);
            for (j = 0; j < J; j++) {
              for (l = 0; l < pAbundRE; l++) {
                betaStarSites[i * J + j] += betaStar[i * nAbundRE + which(XRE[l * J + j], betaLevelIndx, nAbundRE)];
              }
            }
          }

          /********************************************************************
           *Update Detection random effects
           *******************************************************************/
          if (pDetRE > 0) {
            // Update each individual random effect one by one. 
            for (l = 0; l < nDetRE; l++) {
              /********************************
               * Compute b.alpha.star
               *******************************/
              zeros(tmp_one, inc);
              tmp_0 = 0.0;
              for (r = 0; r < nObs; r++) {
                if ((N[NLongIndx[r] * nSp + i] > 0.0) && (XpRE[alphaStarIndx[l] * nObs + r] == alphaLevelIndx[l])) {
                  tmp_one[0] += (yStar[r] - F77_NAME(ddot)(&pDet, &Xp[r], &nObs, &alpha[i], &nSp) - alphaStarObs[i * nObs + r] + alphaStar[i * nDetRE + l]) * omegaDet[r];
          	  tmp_0 += omegaDet[r];
                }
              }
              /********************************
               * Compute A.alpha.star
               *******************************/
              tmp_0 += 1.0 / sigmaSqP[alphaStarIndx[l]]; 
              tmp_0 = 1.0 / tmp_0; 
              alphaStar[i * nDetRE + l] = rnorm(tmp_0 * tmp_one[0], sqrt(tmp_0)); 
            }
            zeros(&alphaStarObs[i * nObs], nObs); 
            // Update the RE sums for the current species
            for (r = 0; r < nObs; r++) {
              for (l = 0; l < pDetRE; l++) {
                alphaStarObs[i * nObs + r] += alphaStar[i * nDetRE + which(XpRE[l * nObs + r], alphaLevelIndx, nDetRE)]; 
              }
            }
          }
          /********************************************************************
           *Update kappa (the NB size parameter)
           *******************************************************************/
          for (j = 0; j < J; j++) {
            tmp_J1[j] = F77_NAME(ddot)(&pAbund, &X[j], &J, &beta[i], &nSp) + betaStarSites[i * J + j];
            psi[j] = logitInv(tmp_J1[j], zero, one);
          }
          /********************************
           * Current
           *******************************/
          // Log gamma function is lgammafn
          // Likelihood contribution
          logPostCurr = 0.0;
          for (j = 0; j < J; j++) {
            logPostCurr += lgammafn(N[j * nSp + i] + kappa[i]) - lgammafn(kappa[i]) + kappa[i] * log(1 - psi[j]);
          }
          logPostCurr += log(kappa[i] - kappaA[i]) + log(kappaB[i] - kappa[i]);
          /********************************
           * Candidate
           *******************************/
          kappaCand = logitInv(rnorm(logit(kappa[i], kappaA[i], kappaB[i]), 
                                     exp(tuning[kappaAMCMCIndx * nSp + i])), kappaA[i], kappaB[i]);
          logPostCand = 0.0;
          for (j = 0; j < J; j++) {
            logPostCand += lgammafn(N[j * nSp + i] + kappaCand) - lgammafn(kappaCand) + kappaCand * log(1 - psi[j]);
          }
          logPostCand += log(kappaCand - kappaA[i]) + log(kappaB[i] - kappaCand);

          if (runif(0.0,1.0) <= exp(logPostCand - logPostCurr)) {
            kappa[i] = kappaCand;
            accept[kappaAMCMCIndx * nSp + i]++;
          }
          /********************************************************************
           *Update Latent Abundance
           *******************************************************************/
	  zeros(logPostCurrN, J);
	  zeros(logPostCandN, J);
	  // Proposal
	  for (j = 0; j < J; j++) {
            NCand[j] = rpois(N[j * nSp + i] + epsilon);
	    mu[j * nSp + i] = exp(tmp_J1[j]) * kappa[i];
	  }
	  // Likelihood contribution to Metropolis ratios
	  for (r = 0; r < nObs; r++) {
            detProb[r * nSp + i] = logitInv(F77_NAME(ddot)(&pDet, &Xp[r], &nObs, &alpha[i], &nSp) + 
	  		        alphaStarObs[i * nObs + r], zero, one);
            logPostCurrN[NLongIndx[r]] += dbinom(y[r * nSp + i], N[NLongIndx[r] * nSp + i], 
	        	                         detProb[r * nSp + i], 1);
	    logPostCandN[NLongIndx[r]] += dbinom(y[r * nSp + i], NCand[NLongIndx[r]], 
	        	                  detProb[r * nSp + i], 1);
	  }
          /********************************
           * Current
           *******************************/
          // Log gamma function is lgammafn
          // Likelihood contribution
          for (j = 0; j < J; j++) {
	    if (NCand[j] >= yMax[j * nSp + i]) {
	      // Contribution from NB(N)
	      logPostCurrN[j] += dnbinom_mu(N[j * nSp + i], kappa[i], mu[j * nSp + i], 1);
	      // MH contribution for assymetric proposal distribution.
	      logPostCurrN[j] += dpois(NCand[j], N[j * nSp + i] + epsilon, 1);
              /********************************
               * Candidate
               *******************************/
              logPostCandN[j] += dnbinom_mu(NCand[j], kappa[i], mu[j * nSp + i], 1);
	      logPostCandN[j] += dpois(N[j * nSp + i], NCand[j] + epsilon, 1);
              if (runif(0.0,1.0) <= exp(logPostCandN[j] - logPostCurrN[j])) {
                N[j * nSp + i] = NCand[j];
              }
	    }
	  }
	} // i (species)

        /********************************************************************
         *Save samples
        *******************************************************************/
        if (g >= nBurn) {
          thinIndx++; 
          if (thinIndx == nThin) {
            F77_NAME(dcopy)(&pAbund, betaComm, &inc, &REAL(betaCommSamples_r)[sPost*pAbund], &inc);
            F77_NAME(dcopy)(&pDet, alphaComm, &inc, &REAL(alphaCommSamples_r)[sPost*pDet], &inc);
            F77_NAME(dcopy)(&pAbund, tauSqBeta, &inc, &REAL(tauSqBetaSamples_r)[sPost*pAbund], &inc);
            F77_NAME(dcopy)(&pDet, tauSqAlpha, &inc, &REAL(tauSqAlphaSamples_r)[sPost*pDet], &inc);
            F77_NAME(dcopy)(&pAbundnSp, beta, &inc, &REAL(betaSamples_r)[sPost*pAbundnSp], &inc); 
            F77_NAME(dcopy)(&pDetnSp, alpha, &inc, &REAL(alphaSamples_r)[sPost*pDetnSp], &inc); 
            F77_NAME(dcopy)(&nSp, kappa, &inc, &REAL(kappaSamples_r)[sPost*nSp], &inc); 
            F77_NAME(dcopy)(&JnSp, N, &inc, &REAL(NSamples_r)[sPost*JnSp], &inc); 
            F77_NAME(dcopy)(&JnSp, mu, &inc, &REAL(muSamples_r)[sPost*JnSp], &inc); 
            if (pDetRE > 0) {
              F77_NAME(dcopy)(&pDetRE, sigmaSqP, &inc, &REAL(sigmaSqPSamples_r)[sPost*pDetRE], &inc);
              F77_NAME(dcopy)(&nDetREnSp, alphaStar, &inc, &REAL(alphaStarSamples_r)[sPost*nDetREnSp], &inc);
            }
            if (pAbundRE > 0) {
              F77_NAME(dcopy)(&pAbundRE, sigmaSqMu, &inc, &REAL(sigmaSqMuSamples_r)[sPost*pAbundRE], &inc);
              F77_NAME(dcopy)(&nAbundREnSp, betaStar, &inc, &REAL(betaStarSamples_r)[sPost*nAbundREnSp], &inc);
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
	  Rprintf("\tSpecies\t\tParameter\tAcceptance\tTuning\n");	  
	  for (i = 0; i < nSp; i++) {
	    Rprintf("\t%i\t\tkappa\t\t%3.1f\t\t%1.5f\n", i + 1, 100.0*REAL(acceptSamples_r)[s * nAMCMC + kappaAMCMCIndx * nSp + i], exp(tuning[kappaAMCMCIndx * nSp + i]));
	  } // i
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
    int nResultListObjs = 9;
    if (pDetRE > 0) {
      nResultListObjs += 2; 
    }
    if (pAbundRE > 0) {
      nResultListObjs += 2;
    }

    PROTECT(result_r = allocVector(VECSXP, nResultListObjs)); nProtect++;
    PROTECT(resultName_r = allocVector(VECSXP, nResultListObjs)); nProtect++;

    SET_VECTOR_ELT(result_r, 0, betaCommSamples_r);
    SET_VECTOR_ELT(result_r, 1, alphaCommSamples_r);
    SET_VECTOR_ELT(result_r, 2, tauSqBetaSamples_r);
    SET_VECTOR_ELT(result_r, 3, tauSqAlphaSamples_r);
    SET_VECTOR_ELT(result_r, 4, betaSamples_r);
    SET_VECTOR_ELT(result_r, 5, alphaSamples_r);
    SET_VECTOR_ELT(result_r, 6, NSamples_r);
    SET_VECTOR_ELT(result_r, 7, muSamples_r);
    SET_VECTOR_ELT(result_r, 8, kappaSamples_r);
    if (pDetRE > 0) {
      SET_VECTOR_ELT(result_r, 9, sigmaSqPSamples_r);
      SET_VECTOR_ELT(result_r, 10, alphaStarSamples_r);
    }
    if (pAbundRE > 0) {
      if (pDetRE > 0) {
        tmp_0 = 11;
      } else {
        tmp_0 = 9;
      }
      SET_VECTOR_ELT(result_r, tmp_0, sigmaSqMuSamples_r);
      SET_VECTOR_ELT(result_r, tmp_0 + 1, betaStarSamples_r);
    }
    SET_VECTOR_ELT(resultName_r, 0, mkChar("beta.comm.samples")); 
    SET_VECTOR_ELT(resultName_r, 1, mkChar("alpha.comm.samples")); 
    SET_VECTOR_ELT(resultName_r, 2, mkChar("tau.sq.beta.samples")); 
    SET_VECTOR_ELT(resultName_r, 3, mkChar("tau.sq.alpha.samples")); 
    SET_VECTOR_ELT(resultName_r, 4, mkChar("beta.samples")); 
    SET_VECTOR_ELT(resultName_r, 5, mkChar("alpha.samples")); 
    SET_VECTOR_ELT(resultName_r, 6, mkChar("N.samples")); 
    SET_VECTOR_ELT(resultName_r, 7, mkChar("mu.samples")); 
    SET_VECTOR_ELT(resultName_r, 8, mkChar("kappa.samples")); 
    if (pDetRE > 0) {
      SET_VECTOR_ELT(resultName_r, 9, mkChar("sigma.sq.p.samples")); 
      SET_VECTOR_ELT(resultName_r, 10, mkChar("alpha.star.samples")); 
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


