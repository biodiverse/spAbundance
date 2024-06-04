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
  SEXP msNMix(SEXP y_r, SEXP X_r, SEXP Xp_r, SEXP XRE_r, SEXP XpRE_r,
              SEXP XRandom_r, SEXP XpRandom_r, SEXP yMax_r,
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
              SEXP samplesInfo_r, SEXP chainInfo_r, SEXP family_r, SEXP offset_r){

    /**********************************************************************
     * Initial constants
     * *******************************************************************/
    int i, j, k, s, t, ll, g, l, q, r, info, nProtect=0;
    const int inc = 1;
    const double one = 1.0;
    const double zero = 0.0;
    char const *lower = "L";
    char const *ytran = "T";

    /**********************************************************************
     * Get Inputs
     * *******************************************************************/
    double *y = REAL(y_r);
    double *X = REAL(X_r);
    double *Xp = REAL(Xp_r);
    int *XRE = INTEGER(XRE_r);
    double *XRandom = REAL(XRandom_r);
    int *XpRE = INTEGER(XpRE_r);
    double *XpRandom = REAL(XpRandom_r);
    double *yMax = REAL(yMax_r);
    double *offset = REAL(offset_r);
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
    int indBetas = INTEGER(consts_r)[9];
    int indAlphas = INTEGER(consts_r)[10];
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
    // NB = 1, Poisson = 0;
    int family = INTEGER(family_r)[0];

#ifdef _OPENMP
    omp_set_num_threads(nThreads);
#else
    if(nThreads > 1){
      Rf_warning("n.omp.threads > 1, but source not compiled with OpenMP support.");
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
          Rprintf("Multi-species Negative Binomial N-Mixture model fit with %i sites and %i species.\n\n", J, nSp);
	} else {
          Rprintf("Multi-species Poisson N-Mixture model fit with %i sites and %i species.\n\n", J, nSp);
	}
        Rprintf("Samples per Chain: %i \n", nSamples);
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
       Some constants and temporary variables to be used later
     * *******************************************************************/
    int pAbundnSp = pAbund * nSp;
    int pDetnSp = pDet * nSp;
    int nObsnSp = nObs * nSp;
    int nAbundREnSp = nAbundRE * nSp;
    int nDetREnSp = nDetRE * nSp;
    int JnSp = J * nSp;
    int JpAbundRE = J * pAbundRE;
    int nObspDetRE = nObs * pDetRE;
    double tmp_0, tmp_02;
    double *tmp_ppDet = (double *) R_alloc(ppDet, sizeof(double));
    double *tmp_ppAbund = (double *) R_alloc(ppAbund, sizeof(double));
    double *tmp_pDet = (double *) R_alloc(pDet, sizeof(double));
    double *tmp_pAbund = (double *) R_alloc(pAbund, sizeof(double));
    double *tmp_pDet2 = (double *) R_alloc(pDet, sizeof(double));
    double *tmp_pAbund2 = (double *) R_alloc(pAbund, sizeof(double));
    int *tmp_JInt = (int *) R_alloc(J, sizeof(int));
    for (j = 0; j < J; j++) {
      tmp_JInt[j] = 0;
    }
    double *tmp_nObs = (double *) R_alloc(nObs, sizeof(double));
    double *tmp_J = (double *) R_alloc(J, sizeof(double));
    zeros(tmp_J, J);

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

    /**********************************************************************
     * Return Stuff
     * *******************************************************************/
    // Community level
    SEXP betaCommSamples_r;
    PROTECT(betaCommSamples_r = Rf_allocMatrix(REALSXP, pAbund, nPost)); nProtect++;
    zeros(REAL(betaCommSamples_r), pAbund * nPost);
    SEXP alphaCommSamples_r;
    PROTECT(alphaCommSamples_r = Rf_allocMatrix(REALSXP, pDet, nPost)); nProtect++;
    zeros(REAL(alphaCommSamples_r), pDet * nPost);
    SEXP tauSqBetaSamples_r;
    PROTECT(tauSqBetaSamples_r = Rf_allocMatrix(REALSXP, pAbund, nPost)); nProtect++;
    zeros(REAL(tauSqBetaSamples_r), pAbund * nPost);
    SEXP tauSqAlphaSamples_r;
    PROTECT(tauSqAlphaSamples_r = Rf_allocMatrix(REALSXP, pDet, nPost)); nProtect++;
    zeros(REAL(tauSqAlphaSamples_r), pDet * nPost);
    // Species level
    SEXP betaSamples_r;
    PROTECT(betaSamples_r = Rf_allocMatrix(REALSXP, pAbundnSp, nPost)); nProtect++;
    zeros(REAL(betaSamples_r), pAbundnSp * nPost);
    SEXP alphaSamples_r;
    PROTECT(alphaSamples_r = Rf_allocMatrix(REALSXP, pDetnSp, nPost)); nProtect++;
    zeros(REAL(alphaSamples_r), pDetnSp * nPost);
    SEXP NSamples_r;
    PROTECT(NSamples_r = Rf_allocMatrix(REALSXP, JnSp, nPost)); nProtect++;
    zeros(REAL(NSamples_r), JnSp * nPost);
    SEXP muSamples_r;
    PROTECT(muSamples_r = Rf_allocMatrix(REALSXP, JnSp, nPost)); nProtect++;
    zeros(REAL(muSamples_r), JnSp * nPost);
    // Detection random effects
    SEXP sigmaSqPSamples_r;
    SEXP alphaStarSamples_r;
    if (pDetRE > 0) {
      PROTECT(sigmaSqPSamples_r = Rf_allocMatrix(REALSXP, pDetRE, nPost)); nProtect++;
      zeros(REAL(sigmaSqPSamples_r), pDetRE * nPost);
      PROTECT(alphaStarSamples_r = Rf_allocMatrix(REALSXP, nDetREnSp, nPost)); nProtect++;
      zeros(REAL(alphaStarSamples_r), nDetREnSp * nPost);
    }
    // Abundance random effects
    SEXP sigmaSqMuSamples_r;
    SEXP betaStarSamples_r;
    if (pAbundRE > 0) {
      PROTECT(sigmaSqMuSamples_r = Rf_allocMatrix(REALSXP, pAbundRE, nPost)); nProtect++;
      zeros(REAL(sigmaSqMuSamples_r), pAbundRE * nPost);
      PROTECT(betaStarSamples_r = Rf_allocMatrix(REALSXP, nAbundREnSp, nPost)); nProtect++;
      zeros(REAL(betaStarSamples_r), nAbundREnSp * nPost);
    }
    SEXP kappaSamples_r;
    if (family == 1) {
      PROTECT(kappaSamples_r = Rf_allocMatrix(REALSXP, nSp, nPost)); nProtect++;
      zeros(REAL(kappaSamples_r), nSp * nPost);
    }
    // Likelihood samples for WAIC.

    /**********************************************************************
     * Additional Sampler Prep
     * *******************************************************************/
    double *detProb = (double *) R_alloc(nObsnSp, sizeof(double));
    double *mu = (double *) R_alloc(JnSp, sizeof(double));
    zeros(mu, JnSp);

    // For normal priors
    F77_NAME(dpotrf)(lower, &pAbund, SigmaBetaCommInv, &pAbund, &info FCONE);
    if(info != 0){Rf_error("c++ error: dpotrf SigmaBetaCommInv failed\n");}
    F77_NAME(dpotri)(lower, &pAbund, SigmaBetaCommInv, &pAbund, &info FCONE);
    if(info != 0){Rf_error("c++ error: dpotri SigmaBetaCommInv failed\n");}
    double *SigmaBetaCommInvMuBeta = (double *) R_alloc(pAbund, sizeof(double));
    F77_NAME(dsymv)(lower, &pAbund, &one, SigmaBetaCommInv, &pAbund, muBetaComm,
		    &inc, &zero, SigmaBetaCommInvMuBeta, &inc FCONE);
    // Detection regression coefficient priors.
    F77_NAME(dpotrf)(lower, &pDet, SigmaAlphaCommInv, &pDet, &info FCONE);
    if(info != 0){Rf_error("c++ error: dpotrf SigmaAlphaCommInv failed\n");}
    F77_NAME(dpotri)(lower, &pDet, SigmaAlphaCommInv, &pDet, &info FCONE);
    if(info != 0){Rf_error("c++ error: dpotri SigmaAlphaCommInv failed\n");}
    double *SigmaAlphaCommInvMuAlpha = (double *) R_alloc(pDet, sizeof(double));
    F77_NAME(dsymv)(lower, &pDet, &one, SigmaAlphaCommInv, &pDet, muAlphaComm, &inc, &zero,
                   SigmaAlphaCommInvMuAlpha, &inc FCONE);
    // Put community level variances in a pAbund x PAbund matrix.
    double *TauBetaInv = (double *) R_alloc(ppAbund, sizeof(double)); zeros(TauBetaInv, ppAbund);
    for (i = 0; i < pAbund; i++) {
      TauBetaInv[i * pAbund + i] = tauSqBeta[i];
    } // i
    F77_NAME(dpotrf)(lower, &pAbund, TauBetaInv, &pAbund, &info FCONE);
    if(info != 0){Rf_error("c++ error: dpotrf TauBetaInv failed\n");}
    F77_NAME(dpotri)(lower, &pAbund, TauBetaInv, &pAbund, &info FCONE);
    if(info != 0){Rf_error("c++ error: dpotri TauBetaInv failed\n");}
    // Put community level variances in a pDet x pDet matrix.
    double *TauAlphaInv = (double *) R_alloc(ppDet, sizeof(double)); zeros(TauAlphaInv, ppDet);
    for (i = 0; i < pDet; i++) {
      TauAlphaInv[i * pDet + i] = tauSqAlpha[i];
    } // i
    F77_NAME(dpotrf)(lower, &pDet, TauAlphaInv, &pDet, &info FCONE);
    if(info != 0){Rf_error("c++ error: dpotrf TauAlphaInv failed\n");}
    F77_NAME(dpotri)(lower, &pDet, TauAlphaInv, &pDet, &info FCONE);
    if(info != 0){Rf_error("c++ error: dpotri TauAlphaInv failed\n");}

    /**********************************************************************
     * Prep for random effects
     * *******************************************************************/
    // Site-level sums of the abundance random effects
    double *betaStarSites = (double *) R_alloc(JnSp, sizeof(double));
    double *betaStarSitesCand = (double *) R_alloc(JnSp, sizeof(double));
    zeros(betaStarSites, JnSp);
    int *betaStarLongIndx = (int *) R_alloc(JpAbundRE, sizeof(int));
    // Initial sums
    for (j = 0; j < J; j++) {
      for (l = 0; l < pAbundRE; l++) {
        betaStarLongIndx[l * J + j] = which(XRE[l * J + j], betaLevelIndx, nAbundRE);
        for (i = 0; i < nSp; i++) {
          betaStarSites[i * J + j] += betaStar[i * nAbundRE + betaStarLongIndx[l * J + j]] * XRandom[l * J + j];
	  betaStarSitesCand[i * J + j] = betaStarSites[i * J + j];
        }
      }
    }
    // Observation-level sums of the detection random effects
    double *alphaStarObs = (double *) R_alloc(nObsnSp, sizeof(double));
    double *alphaStarObsCand = (double *) R_alloc(nObsnSp, sizeof(double));
    zeros(alphaStarObs, nObsnSp);
    int *alphaStarLongIndx = (int *) R_alloc(nObspDetRE, sizeof(int));
    // Get sums of the current REs for each site/visit combo for all species
    for (r = 0; r < nObs; r++) {
      for (l = 0; l < pDetRE; l++) {
        alphaStarLongIndx[l * nObs + r] = which(XpRE[l * nObs + r], alphaLevelIndx, nDetRE);
        for (i = 0; i < nSp; i++) {
          alphaStarObs[i * nObs + r] += alphaStar[i * nDetRE + alphaStarLongIndx[l * nObs + r]] * XpRandom[l * nObs + r];
	  alphaStarObsCand[i * nObs + r] = alphaStarObs[i * nObs + r];
        }
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
    int nAMCMC = 0;
    if (pAbundRE > 0) {
      nAMCMC = (pAbund + nAbundRE) * nSp;
    } else {
      nAMCMC = pAbund * nSp;
    }
    if (pDetRE > 0) {
      nAMCMC += (pDet + nDetRE) * nSp;
    } else {
      nAMCMC += pDet * nSp;
    }
    if (family == 1) {
      nAMCMC += nSp;
    }
    int betaAMCMCIndx = 0;
    int alphaAMCMCIndx = betaAMCMCIndx + pAbund * nSp;
    int betaStarAMCMCIndx = alphaAMCMCIndx + pDet * nSp;
    int alphaStarAMCMCIndx = betaStarAMCMCIndx + nAbundRE * nSp;
    int kappaAMCMCIndx = alphaStarAMCMCIndx + nDetRE * nSp;
    double *accept = (double *) R_alloc(nAMCMC, sizeof(double)); zeros(accept, nAMCMC);
    // Set the initial candidate values for everything to the inital values.
    double *betaCand = (double *) R_alloc(pAbund * nSp, sizeof(double));
    // beta is sorted by parameter, then species within parameter.
    for (k = 0; k < pAbund; k++) {
      for (i = 0; i < nSp; i++) {
        betaCand[k * nSp + i] = beta[k * nSp + i];
      }
    }
    double *betaStarCand = (double *) R_alloc(nAbundREnSp, sizeof(double));
    // betaStar is sorted by species, then parameter within species.
    for (j = 0; j < nAbundRE; j++) {
      for (i = 0; i < nSp; i++) {
        betaStarCand[i * nAbundRE + j] = betaStar[i * nAbundRE + j];
      }
    }
    // Set the initial candidate values for everything to the inital values.
    double *alphaCand = (double *) R_alloc(pDet * nSp, sizeof(double));
    // alpha is sorted by parameter, then species within parameter.
    for (k = 0; k < pDet; k++) {
      for (i = 0; i < nSp; i++) {
        alphaCand[k * nSp + i] = alpha[k * nSp + i];
      }
    }
    double *alphaStarCand = (double *) R_alloc(nDetREnSp, sizeof(double));
    // betaStar is sorted by species, then parameter within species.
    for (j = 0; j < nDetRE; j++) {
      for (i = 0; i < nSp; i++) {
        alphaStarCand[i * nDetRE + j] = alphaStar[i * nDetRE + j];
      }
    }
    double kappaCand = 0.0;
    double *logPostCurrN = (double *) R_alloc(J, sizeof(double));
    zeros(logPostCurrN, J);
    double *logPostCandN = (double *) R_alloc(J, sizeof(double));
    zeros(logPostCandN, J);
    double epsilonN = 1;
    kappaCand = kappa[0];
    double *NCand = (double *) R_alloc(J, sizeof(double));
    for (j = 0; j < J; j++) {
      NCand[j] = N[j];
    }
    SEXP acceptSamples_r;
    PROTECT(acceptSamples_r = Rf_allocMatrix(REALSXP, nAMCMC, nBatch)); nProtect++;
    zeros(REAL(acceptSamples_r), nAMCMC * nBatch);
    SEXP tuningSamples_r;
    PROTECT(tuningSamples_r = Rf_allocMatrix(REALSXP, nAMCMC, nBatch)); nProtect++;
    zeros(REAL(tuningSamples_r), nAMCMC * nBatch);

    GetRNGstate();

    for (s = 0, g = 0; s < nBatch; s++) {
      for (t = 0; t < batchLength; t++, g++) {
        /********************************************************************
         Update Community level Abundance Coefficients
         *******************************************************************/
        if (indBetas == 0) {
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
          if(info != 0){Rf_error("c++ error: dpotrf ABetaComm failed\n");}
          F77_NAME(dpotri)(lower, &pAbund, tmp_ppAbund, &pAbund, &info FCONE);
          if(info != 0){Rf_error("c++ error: dpotri ABetaComm failed\n");}
          F77_NAME(dsymv)(lower, &pAbund, &one, tmp_ppAbund, &pAbund, tmp_pAbund, &inc, &zero, tmp_pAbund2, &inc FCONE);
          F77_NAME(dpotrf)(lower, &pAbund, tmp_ppAbund, &pAbund, &info FCONE);
          if(info != 0){Rf_error("c++ error: dpotrf ABetaComm failed\n");}
          mvrnorm(betaComm, tmp_pAbund2, tmp_ppAbund, pAbund);
        }
        /********************************************************************
         Update Community level Detection Coefficients
         *******************************************************************/
        if (indAlphas == 0) {
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
          if(info != 0){Rf_error("c++ error: dpotrf AAlphaComm failed\n");}
          F77_NAME(dpotri)(lower, &pDet, tmp_ppDet, &pDet, &info FCONE);
          if(info != 0){Rf_error("c++ error: dpotri AAlphaComm failed\n");}
          F77_NAME(dsymv)(lower, &pDet, &one, tmp_ppDet, &pDet, tmp_pDet, &inc, &zero, tmp_pDet2, &inc FCONE);
          F77_NAME(dpotrf)(lower, &pDet, tmp_ppDet, &pDet, &info FCONE);
          if(info != 0){Rf_error("c++ error: dpotrf AAlphaComm failed\n");}
          mvrnorm(alphaComm, tmp_pDet2, tmp_ppDet, pDet);
        }

        /********************************************************************
         Update Community Abundance Variance Parameter
        ********************************************************************/
        if (indBetas == 0) {
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
          if(info != 0){Rf_error("c++ error: dpotrf TauBetaInv failed\n");}
          F77_NAME(dpotri)(lower, &pAbund, TauBetaInv, &pAbund, &info FCONE);
          if(info != 0){Rf_error("c++ error: dpotri TauBetaInv failed\n");}
        }
        /********************************************************************
         Update Community Detection Variance Parameter
        ********************************************************************/
        if (indAlphas == 0) {
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
          if(info != 0){Rf_error("c++ error: dpotrf TauAlphaInv failed\n");}
          F77_NAME(dpotri)(lower, &pDet, TauAlphaInv, &pDet, &info FCONE);
          if(info != 0){Rf_error("c++ error: dpotri TauAlphaInv failed\n");}
        }

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
          // beta is ordered by parameter, then species.
          /********************************************************************
           *Update Abundance Regression Coefficients
           *******************************************************************/
          // Proposal
          for (k = 0; k < pAbund; k++) {
            logPostBetaCand = 0.0;
	    logPostBetaCurr = 0.0;
            betaCand[k * nSp + i] = rnorm(beta[k * nSp + i], exp(tuning[betaAMCMCIndx + k * nSp + i]));
            logPostBetaCand += dnorm(betaCand[k * nSp + i], betaComm[k], sqrt(tauSqBeta[k]), 1);
	    logPostBetaCurr += dnorm(beta[k * nSp + i], betaComm[k], sqrt(tauSqBeta[k]), 1);
            for (j = 0; j < J; j++) {
              tmp_J[j] = exp(F77_NAME(ddot)(&pAbund, &X[j], &J, &betaCand[i], &nSp) +
                                betaStarSites[i * J + j]);
	      if (family == 1) {
                logPostBetaCand += nb_logpost(kappa[i], N[j * nSp + i], tmp_J[j], offset[j]);
	      } else {
                logPostBetaCand += poisson_logpost(N[j * nSp + i], tmp_J[j], offset[j]);
	      }
              tmp_J[j] = exp(F77_NAME(ddot)(&pAbund, &X[j], &J, &beta[i], &nSp) +
                                betaStarSites[i * J + j]);
	      if (family == 1) {
                logPostBetaCurr += nb_logpost(kappa[i], N[j * nSp + i], tmp_J[j], offset[j]);
	      } else {
                logPostBetaCurr += poisson_logpost(N[j * nSp + i], tmp_J[j], offset[j]);
	      }
            }
            if (runif(0.0, 1.0) <= exp(logPostBetaCand - logPostBetaCurr)) {
              beta[k * nSp + i] = betaCand[k * nSp + i];
              accept[betaAMCMCIndx + k * nSp + i]++;
            } else {
              betaCand[k * nSp + i] = beta[k * nSp + i];
	    }
          }

	  /********************************************************************
           *Update Detection Regression Coefficients
           *******************************************************************/
          // Proposal
          zeros(tmp_nObs, nObs);
          for (k = 0; k < pDet; k++) {
            logPostAlphaCand = 0.0;
	    logPostAlphaCurr = 0.0;
            alphaCand[k * nSp + i] = rnorm(alpha[k * nSp + i],
			                   exp(tuning[alphaAMCMCIndx + k * nSp + i]));
            logPostAlphaCand += dnorm(alphaCand[k * nSp + i], alphaComm[k],
			              sqrt(tauSqAlpha[k]), 1);
	    logPostAlphaCurr += dnorm(alpha[k * nSp + i], alphaComm[k],
                                      sqrt(tauSqAlpha[k]), 1);
            for (r = 0; r < nObs; r++) {
              if (N[NLongIndx[r] * nSp + i] > 0.0) {
                tmp_nObs[r] = logitInv(F77_NAME(ddot)(&pDet, &Xp[r], &nObs, &alphaCand[i], &nSp) +
	  		             alphaStarObs[i * nObs + r], zero, one);
                logPostAlphaCand += dbinom(y[r * nSp + i], N[NLongIndx[r] * nSp + i], tmp_nObs[r], 1);
                tmp_nObs[r] = logitInv(F77_NAME(ddot)(&pDet, &Xp[r], &nObs, &alpha[i], &nSp) +
	  		             alphaStarObs[i * nObs + r], zero, one);
                logPostAlphaCurr += dbinom(y[r * nSp + i], N[NLongIndx[r] * nSp + i], tmp_nObs[r], 1);
	      }
            }
            if (runif(0.0, 1.0) <= exp(logPostAlphaCand - logPostAlphaCurr)) {
              alpha[k * nSp + i] = alphaCand[k * nSp + i];
              accept[alphaAMCMCIndx + k * nSp + i]++;
            } else {
              alphaCand[k * nSp + i] = alpha[k * nSp + i];
	    }
          }

          /********************************************************************
           *Update abundance random effects
           *******************************************************************/
          if (pAbundRE > 0) {
            for (l = 0; l < nAbundRE; l++) {
	      betaStarCand[i * nAbundRE + l] = rnorm(betaStar[i * nAbundRE + l],
                                                     exp(tuning[betaStarAMCMCIndx + i * nAbundRE + l]));
              logPostBetaStarCand[l] = dnorm(betaStarCand[i * nAbundRE + l], 0.0,
	  		                   sqrt(sigmaSqMu[betaStarIndx[l]]), 1);
              logPostBetaStarCurr[l] = dnorm(betaStar[i * nAbundRE + l], 0.0,
	  		                   sqrt(sigmaSqMu[betaStarIndx[l]]), 1);
	      for (j = 0; j < J; j++) {
                if (XRE[betaStarIndx[l] * J + j] == betaLevelIndx[l]) {
                  // Candidate
                  betaStarSitesCand[i * J + j] = 0.0;
                  for (ll = 0; ll < pAbundRE; ll++) {
                    betaStarSitesCand[i * J + j] += betaStarCand[i * nAbundRE + betaStarLongIndx[ll * J + j]] *
	                                XRandom[ll * J + j];
                  }
                  tmp_J[j] = exp(F77_NAME(ddot)(&pAbund, &X[j], &J, &beta[i], &nSp) +
	  			  betaStarSitesCand[i * J + j]);
                  if (family == 1) {
	  	    logPostBetaStarCand[l] += nb_logpost(kappa[i], N[j * nSp + i], tmp_J[j], offset[j]);
	  	  } else {
	  	    logPostBetaStarCand[l] += poisson_logpost(N[j * nSp + i], tmp_J[j], offset[j]);
	  	  }
	  	  // Current
                    betaStarSites[i * J + j] = 0.0;
                    for (ll = 0; ll < pAbundRE; ll++) {
                      betaStarSites[i * J + j] += betaStar[i * nAbundRE + betaStarLongIndx[ll * J + j]] *
	                                  XRandom[ll * J + j];
                    }
                    tmp_J[j] = exp(F77_NAME(ddot)(&pAbund, &X[j], &J, &beta[i], &nSp) +
	  	  		  betaStarSites[i * J + j]);
                    if (family == 1) {
	  	    logPostBetaStarCurr[l] += nb_logpost(kappa[i], N[j * nSp + i], tmp_J[j], offset[j]);
	  	  } else {
	  	    logPostBetaStarCurr[l] += poisson_logpost(N[j * nSp + i], tmp_J[j], offset[j]);
	  	  }
	        }
	      }
	      if (runif (0.0, 1.0) <= exp(logPostBetaStarCand[l] - logPostBetaStarCurr[l])) {
                betaStar[i * nAbundRE + l] = betaStarCand[i * nAbundRE + l];
	        F77_NAME(dcopy)(&JnSp, betaStarSitesCand, &inc, betaStarSites, &inc);
	        accept[betaStarAMCMCIndx + i * nAbundRE + l]++;
	      } else {
                betaStarCand[i * nAbundRE + l] = betaStar[i * nAbundRE + l];
	        F77_NAME(dcopy)(&JnSp, betaStarSites, &inc, betaStarSitesCand, &inc);
	      }
	    }
	  }

          /********************************************************************
           *Update detection random effects
           *******************************************************************/
          if (pDetRE > 0) {
            for (l = 0; l < nDetRE; l++) {
	      alphaStarCand[i * nDetRE + l] = rnorm(alphaStar[i * nDetRE + l],
                                                    exp(tuning[alphaStarAMCMCIndx + i * nDetRE + l]));
              logPostAlphaStarCand[l] = dnorm(alphaStarCand[i * nDetRE + l], 0.0,
	  		                   sqrt(sigmaSqP[alphaStarIndx[l]]), 1);
              logPostAlphaStarCurr[l] = dnorm(alphaStar[i * nDetRE + l], 0.0,
	  		                   sqrt(sigmaSqP[alphaStarIndx[l]]), 1);
              zeros(tmp_nObs, nObs);
	      for (r = 0; r < nObs; r++) {
                if ((N[NLongIndx[r] * nSp + i] > 0) && (XpRE[alphaStarIndx[l] * nObs + r] == alphaLevelIndx[l])) {
                  // Candidate
                  alphaStarObsCand[i * nObs + r] = 0.0;
                  for (ll = 0; ll < pDetRE; ll++) {
                    alphaStarObsCand[i * nObs + r] += alphaStarCand[i * nDetRE + alphaStarLongIndx[ll * nObs + r]] *
	                                XpRandom[ll * nObs + r];
                  }
                  tmp_nObs[r] = logitInv(F77_NAME(ddot)(&pDet, &Xp[r], &nObs, &alpha[i], &nSp) +
	  		               alphaStarObsCand[i * nObs + r], zero, one);
                  logPostAlphaStarCand[l] += dbinom(y[r * nSp + i], N[NLongIndx[r] * nSp + i], tmp_nObs[r], 1);
	  	  // Current
                  alphaStarObs[i * nObs + r] = 0.0;
                  for (ll = 0; ll < pDetRE; ll++) {
                    alphaStarObs[i * nObs + r] += alphaStar[i * nDetRE + alphaStarLongIndx[ll * nObs + r]] *
	                                XpRandom[ll * nObs + r];
                  }
                  tmp_nObs[r] = logitInv(F77_NAME(ddot)(&pDet, &Xp[r], &nObs, &alpha[i], &nSp) +
	  		               alphaStarObs[i * nObs + r], zero, one);
                  logPostAlphaStarCurr[l] += dbinom(y[r * nSp + i], N[NLongIndx[r] * nSp + i], tmp_nObs[r], 1);
	        }
	      }
	      if (runif (0.0, 1.0) <= exp(logPostAlphaStarCand[l] - logPostAlphaStarCurr[l])) {
                alphaStar[i * nDetRE + l] = alphaStarCand[i * nDetRE + l];
	        F77_NAME(dcopy)(&nObsnSp, alphaStarObsCand, &inc, alphaStarObs, &inc);
	        accept[alphaStarAMCMCIndx + i * nDetRE + l]++;
	      } else {
                alphaStarCand[i * nDetRE + l] = alphaStar[i * nDetRE + l];
	        F77_NAME(dcopy)(&nObsnSp, alphaStarObs, &inc, alphaStarObsCand, &inc);
	      }
	    }
	  }
          /********************************************************************
           *Update kappa (the NB size parameter)
           *******************************************************************/
	  if (family == 1) {
            kappaCand = logitInv(rnorm(logit(kappa[i], kappaA[i], kappaB[i]),
				       exp(tuning[kappaAMCMCIndx + i])),
	  		         kappaA[i], kappaB[i]);
	    logPostKappaCurr = 0.0;
	    logPostKappaCand = 0.0;
	    for (j = 0; j < J; j++) {
              mu[j * nSp + i] = exp(F77_NAME(ddot)(&pAbund, &X[j], &J, &beta[i], &nSp) +
                          betaStarSites[i * J + j]);
              logPostKappaCurr += nb_logpost(kappa[i], N[j * nSp + i], mu[j * nSp + i], offset[j]);
	      logPostKappaCand += nb_logpost(kappaCand, N[j * nSp + i], mu[j * nSp + i], offset[j]);
	    }
	    // Jacobian adjustment
	    logPostKappaCurr += log(kappa[i] - kappaA[i]) + log(kappaB[i] - kappa[i]);
	    logPostKappaCand += log(kappaCand - kappaA[i]) + log(kappaB[i] - kappaCand);
            if (runif(0.0, 1.0) <= exp(logPostKappaCand - logPostKappaCurr)) {
              kappa[i] = kappaCand;
	      accept[kappaAMCMCIndx + i]++;
	    }
	  }

          /********************************************************************
           *Update Latent Abundance
           *******************************************************************/
	  zeros(logPostCurrN, J);
	  zeros(logPostCandN, J);
	  // Proposal
	  for (j = 0; j < J; j++) {
            NCand[j] = rpois(N[j * nSp + i] + epsilonN);
	    // Only calculate if Poisson since its already calculated in kappa update
	    if (family == 0) {
              mu[j * nSp + i] = exp(F77_NAME(ddot)(&pAbund, &X[j], &J, &beta[i], &nSp) +
                          betaStarSites[i * J + j]);
	    }
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
          for (j = 0; j < J; j++) {
	    if (NCand[j] >= yMax[j * nSp + i]) {
              /********************************
               * Current
               *******************************/
	      // Contribution from NB or Poisson
	      if (family == 0) {
	        logPostCurrN[j] += poisson_logpost(N[j * nSp + i], mu[j * nSp + i], offset[j]);
	      } else {
	        logPostCurrN[j] += nb_logpost(kappa[i], N[j * nSp + i], mu[j * nSp + i], offset[j]);
	      }
	      // MH contribution for assymetric proposal distribution.
	      logPostCurrN[j] += poisson_logpost(NCand[j], N[j * nSp + i] + epsilonN, 1.0);
              /********************************
               * Candidate
               *******************************/
	      if (family == 0) {
                logPostCandN[j] += poisson_logpost(NCand[j], mu[j * nSp + i], offset[j]);
	      } else {
                logPostCandN[j] += nb_logpost(kappa[i], NCand[j], mu[j * nSp + i], offset[j]);
	      }
	      // MH contribution for assymetric proposal distribution
	      logPostCandN[j] += poisson_logpost(N[j * nSp + i], NCand[j] + epsilonN, 1.0);
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
	    if (family == 1) {
              F77_NAME(dcopy)(&nSp, kappa, &inc, &REAL(kappaSamples_r)[sPost*nSp], &inc);
	    }
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
          Rprintf("\tSpecies\tParameter\tAcceptance\tTuning\n");
          for (ll = 0; ll < nSp; ll++) {
            Rprintf("\t%i\tbeta[1]\t\t%3.1f\t\t%1.5f\n", ll + 1, 100.0*REAL(acceptSamples_r)[s * nAMCMC + betaAMCMCIndx + ll], exp(tuning[betaAMCMCIndx + ll]));
            Rprintf("\t%i\talpha[1]\t%3.1f\t\t%1.5f\n", ll + 1, 100.0*REAL(acceptSamples_r)[s * nAMCMC + alphaAMCMCIndx + ll], exp(tuning[alphaAMCMCIndx + ll]));
          } // ll
	  if (family == 1) {
	    for (i = 0; i < nSp; i++) {
	      Rprintf("\t%i\tkappa\t\t%3.1f\t\t%1.5f\n", i + 1, 100.0*REAL(acceptSamples_r)[s * nAMCMC + kappaAMCMCIndx + i], exp(tuning[kappaAMCMCIndx + i]));
	    } // i
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
    int nResultListObjs = 8;
    if (pDetRE > 0) {
      nResultListObjs += 2;
    }
    if (pAbundRE > 0) {
      nResultListObjs += 2;
    }
    if (family == 1) {
      nResultListObjs += 2;
    }

    PROTECT(result_r = Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;
    PROTECT(resultName_r = Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;

    SET_VECTOR_ELT(result_r, 0, betaCommSamples_r);
    SET_VECTOR_ELT(result_r, 1, alphaCommSamples_r);
    SET_VECTOR_ELT(result_r, 2, tauSqBetaSamples_r);
    SET_VECTOR_ELT(result_r, 3, tauSqAlphaSamples_r);
    SET_VECTOR_ELT(result_r, 4, betaSamples_r);
    SET_VECTOR_ELT(result_r, 5, alphaSamples_r);
    SET_VECTOR_ELT(result_r, 6, NSamples_r);
    SET_VECTOR_ELT(result_r, 7, muSamples_r);
    if (pDetRE > 0) {
      tmp_0 = 8; // Needed to make tracking kappa easier.
      SET_VECTOR_ELT(result_r, 8, sigmaSqPSamples_r);
      SET_VECTOR_ELT(result_r, 9, alphaStarSamples_r);
    }
    if (pAbundRE > 0) {
      if (pDetRE > 0) {
        tmp_0 = 10;
      } else {
        tmp_0 = 8;
      }
      SET_VECTOR_ELT(result_r, tmp_0, sigmaSqMuSamples_r);
      SET_VECTOR_ELT(result_r, tmp_0 + 1, betaStarSamples_r);
    }
    if (family == 1) {
      if ((pDetRE > 0) || (pAbundRE > 0)) {
        tmp_02 = tmp_0 + 2;
      } else {
        tmp_02 = 8;
      }
      SET_VECTOR_ELT(result_r, tmp_02, kappaSamples_r);
    }
    SET_VECTOR_ELT(resultName_r, 0, Rf_mkChar("beta.comm.samples"));
    SET_VECTOR_ELT(resultName_r, 1, Rf_mkChar("alpha.comm.samples"));
    SET_VECTOR_ELT(resultName_r, 2, Rf_mkChar("tau.sq.beta.samples"));
    SET_VECTOR_ELT(resultName_r, 3, Rf_mkChar("tau.sq.alpha.samples"));
    SET_VECTOR_ELT(resultName_r, 4, Rf_mkChar("beta.samples"));
    SET_VECTOR_ELT(resultName_r, 5, Rf_mkChar("alpha.samples"));
    SET_VECTOR_ELT(resultName_r, 6, Rf_mkChar("N.samples"));
    SET_VECTOR_ELT(resultName_r, 7, Rf_mkChar("mu.samples"));
    if (pDetRE > 0) {
      SET_VECTOR_ELT(resultName_r, 8, Rf_mkChar("sigma.sq.p.samples"));
      SET_VECTOR_ELT(resultName_r, 9, Rf_mkChar("alpha.star.samples"));
    }
    if (pAbundRE > 0) {
      SET_VECTOR_ELT(resultName_r, tmp_0, Rf_mkChar("sigma.sq.mu.samples"));
      SET_VECTOR_ELT(resultName_r, tmp_0 + 1, Rf_mkChar("beta.star.samples"));
    }
    if (family == 1) {
      SET_VECTOR_ELT(resultName_r, tmp_02, Rf_mkChar("kappa.samples"));
    }

    Rf_namesgets(result_r, resultName_r);

    UNPROTECT(nProtect);

    return(result_r);
  }
}


