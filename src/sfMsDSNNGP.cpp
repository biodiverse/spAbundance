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

void updateBFSFDS(double *B, double *F, double *c, double *C, double *coords, int *nnIndx, int *nnIndxLU, int n, int m, double sigmaSq, double phi, double nu, int covModel, double *bk, double nuUnifb){

  int i, k, l;
  int info = 0;
  int inc = 1;
  double one = 1.0;
  double zero = 0.0;
  char lower = 'L';

  //bk must be 1+(int)floor(alpha) * nthread
  int nb = 1+static_cast<int>(floor(nuUnifb));
  int threadID = 0;
  double e;
  int mm = m*m;

#ifdef _OPENMP
#pragma omp parallel for private(k, l, info, threadID, e)
#endif
    for(i = 0; i < n; i++){
#ifdef _OPENMP
      threadID = omp_get_thread_num();
#endif
      if(i > 0){
	for(k = 0; k < nnIndxLU[n+i]; k++){
	  e = dist2(coords[i], coords[n+i], coords[nnIndx[nnIndxLU[i]+k]], coords[n+nnIndx[nnIndxLU[i]+k]]);
	  c[m*threadID+k] = sigmaSq*spCor(e, phi, nu, covModel, &bk[threadID*nb]);
	  for(l = 0; l <= k; l++){
	    e = dist2(coords[nnIndx[nnIndxLU[i]+k]], coords[n+nnIndx[nnIndxLU[i]+k]], coords[nnIndx[nnIndxLU[i]+l]], coords[n+nnIndx[nnIndxLU[i]+l]]);
	    C[mm*threadID+l*nnIndxLU[n+i]+k] = sigmaSq*spCor(e, phi, nu, covModel, &bk[threadID*nb]);
	  }
	}
	F77_NAME(dpotrf)(&lower, &nnIndxLU[n+i], &C[mm*threadID], &nnIndxLU[n+i], &info FCONE); if(info != 0){error("c++ error: dpotrf failed\n");}
	F77_NAME(dpotri)(&lower, &nnIndxLU[n+i], &C[mm*threadID], &nnIndxLU[n+i], &info FCONE); if(info != 0){error("c++ error: dpotri failed\n");}
	F77_NAME(dsymv)(&lower, &nnIndxLU[n+i], &one, &C[mm*threadID], &nnIndxLU[n+i], &c[m*threadID], &inc, &zero, &B[nnIndxLU[i]], &inc FCONE);
	F[i] = sigmaSq - F77_NAME(ddot)(&nnIndxLU[n+i], &B[nnIndxLU[i]], &inc, &c[m*threadID], &inc);
      }else{
	B[i] = 0;
	F[i] = sigmaSq;
      }
    }

}

extern "C" {
  SEXP sfMsDSNNGP(SEXP y_r, SEXP X_r, SEXP Xp_r, SEXP coords_r, SEXP XRE_r, SEXP XpRE_r, 
                  SEXP XRandom_r, SEXP XpRandom_r, SEXP yMax_r, SEXP offset_r,
                  SEXP consts_r, SEXP nAbundRELong_r, SEXP nDetRELong_r,
                  SEXP nnIndx_r, SEXP nnIndxLU_r, SEXP uIndx_r, SEXP uIndxLU_r, SEXP uiIndx_r, 
                  SEXP betaStarting_r, SEXP alphaStarting_r, SEXP kappaStarting_r,
                  SEXP NStarting_r, SEXP betaCommStarting_r, SEXP alphaCommStarting_r,
	          SEXP phiStarting_r, SEXP lambdaStarting_r, SEXP nuStarting_r, SEXP wStarting_r,
	          SEXP tauSqBetaStarting_r, SEXP tauSqAlphaStarting_r,
                  SEXP sigmaSqMuStarting_r, SEXP sigmaSqPStarting_r,
                  SEXP betaStarStarting_r, SEXP alphaStarStarting_r, 
                  SEXP NLongIndx_r, SEXP betaStarIndx_r, SEXP betaLevelIndx_r, 
                  SEXP alphaStarIndx_r, SEXP alphaLevelIndx_r,
                  SEXP muBetaComm_r, SEXP SigmaBetaComm_r, 
                  SEXP muAlphaComm_r, SEXP SigmaAlphaComm_r, 
                  SEXP sigmaSqMuA_r, SEXP sigmaSqMuB_r, 
                  SEXP sigmaSqPA_r, SEXP sigmaSqPB_r, 
                  SEXP kappaA_r, SEXP kappaB_r, 
	          SEXP tauSqBetaA_r, SEXP tauSqBetaB_r, 
	          SEXP tauSqAlphaA_r, SEXP tauSqAlphaB_r, SEXP spatialPriors_r, 
	          SEXP transect_r, SEXP distBreaks_r, SEXP tuning_r, 
                  SEXP nBatch_r, SEXP batchLength_r, SEXP acceptRate_r, SEXP nThreads_r, 
                  SEXP verbose_r, SEXP nReport_r, SEXP samplesInfo_r,
                  SEXP chainInfo_r){
   
    /**********************************************************************
     * Initial constants
     * *******************************************************************/
    int i, g, t, j, s, r, l, k, ll, ii, jj, info, nProtect=0;
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
    double *coords = REAL(coords_r); 
    double *yMax = REAL(yMax_r);
    double *offset = REAL(offset_r);
    int *XRE = INTEGER(XRE_r); 
    int *XpRE = INTEGER(XpRE_r);
    double *XRandom = REAL(XRandom_r);
    double *XpRandom = REAL(XpRandom_r);
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
    int q = INTEGER(consts_r)[9]; 
    int K = INTEGER(consts_r)[10];
    int m = INTEGER(consts_r)[11];
    // NB = 1, Poisson = 0;
    int family = INTEGER(consts_r)[12];
    int covModel = INTEGER(consts_r)[13];
    std::string corName = getCorName(covModel);
    int detModel = INTEGER(consts_r)[14];
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
    double *phiA = REAL(spatialPriors_r); 
    double *phiB = &REAL(spatialPriors_r)[q]; 
    double *nuA = &REAL(spatialPriors_r)[2*q]; 
    double *nuB = &REAL(spatialPriors_r)[3*q]; 
    int *nnIndx = INTEGER(nnIndx_r);
    int *nnIndxLU = INTEGER(nnIndxLU_r);
    int *uIndx = INTEGER(uIndx_r);
    int *uIndxLU = INTEGER(uIndxLU_r);
    int *nAbundRELong = INTEGER(nAbundRELong_r); 
    int *nDetRELong = INTEGER(nDetRELong_r); 
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
    // DS stuff
    int transect = INTEGER(transect_r)[0];
    double *distBreaks = REAL(distBreaks_r);

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
          Rprintf("Spatial Factor NNGP Multi-species Negative Binomial HDS model with %i sites and %i species.\n\n", J, nSp);
	} else {
          Rprintf("Spatial Factor Multi-species Poisson HDS model with %i sites and %i species.\n\n", J, nSp);
	}
        Rprintf("Samples per Chain: %i (%i batches of length %i)\n", nSamples, nBatch, batchLength);
        Rprintf("Burn-in: %i \n", nBurn); 
        Rprintf("Thinning Rate: %i \n", nThin); 
        Rprintf("Number of Chains: %i \n", nChain);
        Rprintf("Total Posterior Samples: %i \n\n", nPost * nChain); 
        Rprintf("Using the %s spatial correlation model.\n\n", corName.c_str());
        Rprintf("Using %i latent spatial factors.\n", q);
        Rprintf("Using %i nearest neighbors.\n\n", m);
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

    /********************************************************************
      Some constants and temporary variables to be used later
    ********************************************************************/
    int pAbundnSp = pAbund * nSp; 
    int pDetnSp = pDet * nSp; 
    int nAbundREnSp = nAbundRE * nSp; 
    int nDetREnSp = nDetRE * nSp; 
    int JnSp = J * nSp;
    int JpAbundRE = J * pAbundRE;
    int JpDetRE = J * pDetRE;
    int KFull = K + 1;
    int nObsFull = KFull * J;
    int nObsFullnSp = nObsFull * nSp;
    int Jq = J * q;
    int nSpq = nSp * q;
    double tmp_0, tmp_02; 
    double *tmp_ppDet = (double *) R_alloc(ppDet, sizeof(double));
    double *tmp_ppAbund = (double *) R_alloc(ppAbund, sizeof(double)); 
    double *tmp_pDet = (double *) R_alloc(pDet, sizeof(double));
    double *tmp_pAbund = (double *) R_alloc(pAbund, sizeof(double));
    double *tmp_pDet2 = (double *) R_alloc(pDet, sizeof(double));
    double *tmp_pAbund2 = (double *) R_alloc(pAbund, sizeof(double));
    double *tmp_J = (double *) R_alloc(J, sizeof(double));
    double *tmp_KFull = (double *) R_alloc(KFull, sizeof(double));
    int *tmp_KFullInt = (int *) R_alloc(KFull, sizeof(int));
   
    // For latent abundance
    double *mu = (double *) R_alloc(JnSp, sizeof(double)); 
    zeros(mu, JnSp); 

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
    // Latent factors
    double *w = (double *) R_alloc(Jq, sizeof(double));
    F77_NAME(dcopy)(&Jq, REAL(wStarting_r), &inc, w, &inc);
    // Latent factor loadings
    double *lambda = (double *) R_alloc(nSpq, sizeof(double));
    F77_NAME(dcopy)(&nSpq, REAL(lambdaStarting_r), &inc, lambda, &inc);
    // Spatial range parameter
    double *phi = (double *) R_alloc(q, sizeof(double)); 
    F77_NAME(dcopy)(&q, REAL(phiStarting_r), &inc, phi, &inc); 
    // Spatial smoothing parameter for Matern
    double *nu = (double *) R_alloc(q, sizeof(double)); 
    F77_NAME(dcopy)(&q, REAL(nuStarting_r), &inc, nu, &inc); 
    // Abundance random effect variances
    double *sigmaSqMu = (double *) R_alloc(pAbundRE, sizeof(double)); 
    F77_NAME(dcopy)(&pAbundRE, REAL(sigmaSqMuStarting_r), &inc, sigmaSqMu, &inc); 
    // Latent random effects
    double *betaStar = (double *) R_alloc(nAbundREnSp, sizeof(double)); 
    F77_NAME(dcopy)(&nAbundREnSp, REAL(betaStarStarting_r), &inc, betaStar, &inc); 
    // Detection random effect variances
    double *sigmaSqP = (double *) R_alloc(pDetRE, sizeof(double)); 
    F77_NAME(dcopy)(&pDetRE, REAL(sigmaSqPStarting_r), &inc, sigmaSqP, &inc); 
    // Latent detection random effects
    double *alphaStar = (double *) R_alloc(nDetREnSp, sizeof(double)); 
    F77_NAME(dcopy)(&nDetREnSp, REAL(alphaStarStarting_r), &inc, alphaStar, &inc); 
    // Size parameter
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
    PROTECT(betaCommSamples_r = allocMatrix(REALSXP, pAbund, nPost)); nProtect++;
    zeros(REAL(betaCommSamples_r), pAbund * nPost);
    SEXP alphaCommSamples_r;
    PROTECT(alphaCommSamples_r = allocMatrix(REALSXP, pDet, nPost)); nProtect++;
    zeros(REAL(alphaCommSamples_r), pDet * nPost);
    SEXP tauSqBetaSamples_r; 
    PROTECT(tauSqBetaSamples_r = allocMatrix(REALSXP, pAbund, nPost)); nProtect++; 
    zeros(REAL(tauSqBetaSamples_r), pAbund * nPost);
    SEXP tauSqAlphaSamples_r; 
    PROTECT(tauSqAlphaSamples_r = allocMatrix(REALSXP, pDet, nPost)); nProtect++; 
    zeros(REAL(tauSqAlphaSamples_r), pDet * nPost);
    // Species level
    SEXP betaSamples_r;
    PROTECT(betaSamples_r = allocMatrix(REALSXP, pAbundnSp, nPost)); nProtect++;
    zeros(REAL(betaSamples_r), pAbundnSp * nPost);
    SEXP alphaSamples_r; 
    PROTECT(alphaSamples_r = allocMatrix(REALSXP, pDetnSp, nPost)); nProtect++;
    zeros(REAL(alphaSamples_r), pDetnSp * nPost);
    SEXP NSamples_r; 
    PROTECT(NSamples_r = allocMatrix(REALSXP, JnSp, nPost)); nProtect++; 
    zeros(REAL(NSamples_r), JnSp * nPost);
    SEXP muSamples_r; 
    PROTECT(muSamples_r = allocMatrix(REALSXP, JnSp, nPost)); nProtect++; 
    zeros(REAL(muSamples_r), JnSp * nPost);
    // Latent factor parameters
    SEXP lambdaSamples_r; 
    PROTECT(lambdaSamples_r = allocMatrix(REALSXP, nSpq, nPost)); nProtect++;
    zeros(REAL(lambdaSamples_r), nSpq * nPost);
    SEXP wSamples_r; 
    PROTECT(wSamples_r = allocMatrix(REALSXP, Jq, nPost)); nProtect++; 
    zeros(REAL(wSamples_r), Jq * nPost);
    // Detection random effects
    SEXP sigmaSqPSamples_r; 
    SEXP alphaStarSamples_r; 
    if (pDetRE > 0) {
      PROTECT(sigmaSqPSamples_r = allocMatrix(REALSXP, pDetRE, nPost)); nProtect++;
      zeros(REAL(sigmaSqPSamples_r), pDetRE * nPost);
      PROTECT(alphaStarSamples_r = allocMatrix(REALSXP, nDetREnSp, nPost)); nProtect++;
      zeros(REAL(alphaStarSamples_r), nDetREnSp * nPost);
    }
    // Abundance random effects
    SEXP sigmaSqMuSamples_r; 
    SEXP betaStarSamples_r; 
    if (pAbundRE > 0) {
      PROTECT(sigmaSqMuSamples_r = allocMatrix(REALSXP, pAbundRE, nPost)); nProtect++;
      zeros(REAL(sigmaSqMuSamples_r), pAbundRE * nPost);
      PROTECT(betaStarSamples_r = allocMatrix(REALSXP, nAbundREnSp, nPost)); nProtect++;
      zeros(REAL(betaStarSamples_r), nAbundREnSp * nPost);
    }
    SEXP kappaSamples_r;
    if (family == 1) {
      PROTECT(kappaSamples_r = allocMatrix(REALSXP, nSp, nPost)); nProtect++;
      zeros(REAL(kappaSamples_r), nSp * nPost);
    }
    SEXP yRepSamples_r; 
    PROTECT(yRepSamples_r = allocMatrix(INTSXP, nObsFullnSp, nPost)); nProtect++; 
    zerosInt(INTEGER(yRepSamples_r), nObsFullnSp * nPost);
    SEXP piFullSamples_r; 
    PROTECT(piFullSamples_r = allocMatrix(REALSXP, nObsFullnSp, nPost)); nProtect++; 
    zeros(REAL(piFullSamples_r), nObsFullnSp * nPost);

    /**********************************************************************
     * Additional Sampler Prep
     * *******************************************************************/
    // For normal priors
    F77_NAME(dpotrf)(lower, &pAbund, SigmaBetaCommInv, &pAbund, &info FCONE); 
    if(info != 0){error("c++ error: dpotrf SigmaBetaCommInv failed\n");}
    F77_NAME(dpotri)(lower, &pAbund, SigmaBetaCommInv, &pAbund, &info FCONE); 
    if(info != 0){error("c++ error: dpotri SigmaBetaCommInv failed\n");}
    double *SigmaBetaCommInvMuBeta = (double *) R_alloc(pAbund, sizeof(double)); 
    F77_NAME(dsymv)(lower, &pAbund, &one, SigmaBetaCommInv, &pAbund, muBetaComm, 
		    &inc, &zero, SigmaBetaCommInvMuBeta, &inc FCONE);
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

    /**********************************************************************
     Set up spatial stuff
     * *******************************************************************/
    // Note that even though sigmaSq is fixed at 1 for spatial factor models, 
    // still maintaining the same approach to leverage the same correlation
    // functions underneath. 
    int nTheta, sigmaSqIndx, phiIndx, nuIndx;
    if (corName != "matern") {
      nTheta = 2; // sigma^2, phi 
      sigmaSqIndx = 0; phiIndx = 1; 
    } else {
      nTheta = 3; // sigma^2, phi, nu 
      sigmaSqIndx = 0; phiIndx = 1; nuIndx = 2; 
    }
    int nThetaq = nTheta * q; 
    int nThetaqSave = (nTheta - 1) * q;
    double *theta = (double *) R_alloc(nThetaq, sizeof(double));
    for (ll = 0; ll < q; ll++) {
      theta[phiIndx * q + ll] = phi[ll];
      // sigmaSq by default is 1 for spatial factor models. 
      theta[sigmaSqIndx * q + ll] = 1.0;
      if (corName == "matern") {
        theta[nuIndx * q + ll] = nu[ll]; 
      } 
    } // ll
    SEXP thetaSamples_r; 
    PROTECT(thetaSamples_r = allocMatrix(REALSXP, nThetaqSave, nPost)); nProtect++; 
    zeros(REAL(thetaSamples_r), nThetaqSave * nPost);
    // Species-level spatial random effects
    double *wStar = (double *) R_alloc(JnSp, sizeof(double)); zeros(wStar, JnSp);
    // Multiply Lambda %*% w[j] to get wStar. 
    for (j = 0; j < J; j++) {
      F77_NAME(dgemv)(ntran, &nSp, &q, &one, lambda, &nSp, &w[j*q], &inc, &zero, &wStar[j * nSp], &inc FCONE);
    }
    // For NNGP
    double b, e, aa; 
    double *a = (double *) R_alloc(q, sizeof(double));
    double *v = (double *) R_alloc(q, sizeof(double));

    // Allocate for the U index vector that keep track of which locations have 
    // the i-th location as a neighbor
    int nIndx = static_cast<int>(static_cast<double>(1+m)/2*m+(J-m-1)*m);

    // For NNGP. Create a copy of these for each species. Increases storage 
    // space that is needed, but reduces amount of computations. 
    int mm = m*m;
    double *B = (double *) R_alloc(nIndx * q, sizeof(double)); 
    double *F = (double *) R_alloc(J * q, sizeof(double));
    // Only need one of these. 
    double *BCand = (double *) R_alloc(nIndx, sizeof(double));
    double *FCand = (double *) R_alloc(J, sizeof(double));
    double *c =(double *) R_alloc(m*nThreads*q, sizeof(double));
    double *C = (double *) R_alloc(mm*nThreads*q, sizeof(double));
    int sizeBK = nThreads*(1.0+static_cast<int>(floor(nuB[0])));
    double *bk = (double *) R_alloc(q*sizeBK, sizeof(double));

    // Initiate B and F for each species
    for (ll = 0; ll < q; ll++) {
      updateBFSFDS(&B[ll * nIndx], &F[ll*J], &c[ll * m*nThreads], &C[ll * mm * nThreads], coords, nnIndx, nnIndxLU, J, m, theta[sigmaSqIndx * q + ll], theta[phiIndx * q + ll], nu[ll], covModel, &bk[ll * sizeBK], nuB[0]);
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
    double *logPostWCand = (double *) R_alloc(Jq, sizeof(double));
    double *logPostWCurr = (double *) R_alloc(Jq, sizeof(double));
    for (j = 0; j < Jq; j++) {
      logPostWCand[j] = R_NegInf;
      logPostWCurr[j] = R_NegInf;
    }
    double *logPostLambdaCand = (double *) R_alloc(nSpq, sizeof(double));
    double *logPostLambdaCurr = (double *) R_alloc(nSpq, sizeof(double));
    for (j = 0; j < nSpq; j++) {
      logPostLambdaCand[j] = R_NegInf;
      logPostLambdaCurr[j] = R_NegInf;
    }
    double *logPostThetaCurr = (double *) R_alloc(q, sizeof(double));
    double *logPostThetaCand = (double *) R_alloc(q, sizeof(double));
    for (j = 0; j < q; j++) {
      logPostThetaCurr[j] = R_NegInf;
      logPostThetaCand[j] = R_NegInf;
    }
    int nAMCMC = 0;
    if (pAbundRE > 0) {
      nAMCMC = (pAbund + nAbundRE) * nSp + nSpq + Jq + nThetaq;
    } else {
      nAMCMC = pAbund * nSp + nSpq + Jq + nThetaq;
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
    int sigmaSqAMCMCIndx = alphaStarAMCMCIndx + nDetRE * nSp;
    int phiAMCMCIndx = sigmaSqAMCMCIndx + q;
    int nuAMCMCIndx;
    if (corName == "matern") {
      nuAMCMCIndx = phiAMCMCIndx + q;
    } else {
      nuAMCMCIndx = phiAMCMCIndx;
    }
    int lambdaAMCMCIndx = nuAMCMCIndx + q;
    int wAMCMCIndx = lambdaAMCMCIndx + nSpq;
    int kappaAMCMCIndx = wAMCMCIndx + Jq;
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
    double *wCand = (double *) R_alloc(Jq, sizeof(double));
    for (j = 0; j < J; j++) {
      for (ll = 0; ll < q; ll++) {
        wCand[j * q + ll] = w[j * q + ll];
      }
    } // j
    double *wStarCand = (double *) R_alloc(JnSp, sizeof(double));
    for (j = 0; j < J; j++) {
      for (i = 0; i < nSp; i++) {
        wStarCand[j * nSp + i] = wStar[j * nSp + i];
      } // i
    } // j
    // lambda is ordered by factor, then species within factor.
    double *lambdaCand = (double *) R_alloc(nSpq, sizeof(double));
    for (i = 0; i < nSp; i++) {
      for (ll = 0; ll < q; ll++) {
        lambdaCand[ll * nSp + i] = lambda[ll * nSp + i];
      } // ll 
    } // i
    // phi is ordered by factor
    double logDet;
    double *phiCand = (double *) R_alloc(q, sizeof(double));
    double *nuCand = (double *) R_alloc(q, sizeof(double)); zeros(nuCand, q);
    for (ll = 0; ll < q; ll++) {
      phiCand[ll] = phi[ll];
      nuCand[ll] = nu[ll];
    }
    double *logPostCurrN = (double *) R_alloc(J, sizeof(double));
    zeros(logPostCurrN, J);
    double *logPostCandN = (double *) R_alloc(J, sizeof(double));
    zeros(logPostCandN, J);
    double epsilon = 1;
    kappaCand = kappa[0];
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
    double *alphaStarSites = (double *) R_alloc(JnSp, sizeof(double)); 
    double *alphaStarSitesCand = (double *) R_alloc(JnSp, sizeof(double)); 
    zeros(alphaStarSites, JnSp); 
    int *alphaStarLongIndx = (int *) R_alloc(JpDetRE, sizeof(int));
    // Get sums of the current REs for each site/visit combo for all species
    for (r = 0; r < J; r++) {
      for (l = 0; l < pDetRE; l++) {
        alphaStarLongIndx[l * J + r] = which(XpRE[l * J + r], alphaLevelIndx, nDetRE);
        for (i = 0; i < nSp; i++) {
          alphaStarSites[i * J + r] += alphaStar[i * nDetRE + alphaStarLongIndx[l * J + r]] * XpRandom[l * J + r];
	  alphaStarSitesCand[i * J + r] = alphaStarSites[i * J + r];
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
    
    /**********************************************************************
     * DS Prep 
     * *******************************************************************/
    double *sigma = (double *) R_alloc(J, sizeof(double));
    double *p = (double *) R_alloc(nObs, sizeof(double)); zeros(p, nObs);
    // Number of break points for integration
    int nInt = 5;
    double *binWidth = (double *) R_alloc(K, sizeof(double));
    double stripWidth = 0.0;
    double *psi = (double *) R_alloc(K, sizeof(double));
    for (k = 0; k < K; k++) {
      binWidth[k] = distBreaks[k + 1] - distBreaks[k];
      stripWidth += binWidth[k];
    }
    for (k = 0; k < K; k++) {
      if (transect == 0) {
        psi[k] = binWidth[k] / stripWidth;
      } else {
        psi[k] = (pow(distBreaks[k + 1], 2) - pow(distBreaks[k], 2)) / pow(stripWidth, 2);
      }
    }
    double *piFull = (double *) R_alloc(nObsFullnSp, sizeof(double)); 
    zeros(piFull, nObsFullnSp);
    double *piFullCand = (double *) R_alloc(nObsFullnSp, sizeof(double)); 
    zeros(piFullCand, nObsFullnSp);
    double likeVal = 0.0;

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
        for (ll = 0; ll < pAbund; ll++) {
          tmp_pAbund[ll] += SigmaBetaCommInvMuBeta[ll];  
        } // j

        /********************************
         Compute A.beta.comm
         *******************************/
        for (ll = 0; ll < ppAbund; ll++) {
          tmp_ppAbund[ll] = SigmaBetaCommInv[ll] + nSp * TauBetaInv[ll]; 
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
         for (ll = 0; ll < pDet; ll++) {
           tmp_pDet[ll] += SigmaAlphaCommInvMuAlpha[ll];  
         } // j
        /********************************
         * Compute A.alpha.comm
         *******************************/
        for (ll = 0; ll < ppDet; ll++) {
          tmp_ppDet[ll] = SigmaAlphaCommInv[ll] + nSp * TauAlphaInv[ll]; 
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
        for (ll = 0; ll < pAbund; ll++) {
          tmp_0 = 0.0;  
          for (i = 0; i < nSp; i++) {
            tmp_0 += (beta[ll * nSp + i] - betaComm[ll]) * (beta[ll * nSp + i] - betaComm[ll]);
          } // i
          tmp_0 *= 0.5;
          tauSqBeta[ll] = rigamma(tauSqBetaA[ll] + nSp / 2.0, tauSqBetaB[ll] + tmp_0); 
        } // ll 
        for (ll = 0; ll < pAbund; ll++) {
          TauBetaInv[ll * pAbund + ll] = tauSqBeta[ll]; 
        } // ll 
        F77_NAME(dpotrf)(lower, &pAbund, TauBetaInv, &pAbund, &info FCONE); 
        if(info != 0){error("c++ error: dpotrf TauBetaInv failed\n");}
        F77_NAME(dpotri)(lower, &pAbund, TauBetaInv, &pAbund, &info FCONE); 
        if(info != 0){error("c++ error: dpotri TauBetaInv failed\n");}
        /********************************************************************
         Update Community Detection Variance Parameter
        ********************************************************************/
        for (ll = 0; ll < pDet; ll++) {
          tmp_0 = 0.0;  
          for (i = 0; i < nSp; i++) {
            tmp_0 += (alpha[ll * nSp + i] - alphaComm[ll]) * (alpha[ll * nSp + i] - alphaComm[ll]);
          } // i
          tmp_0 *= 0.5;
          tauSqAlpha[ll] = rigamma(tauSqAlphaA[ll] + nSp / 2.0, tauSqAlphaB[ll] + tmp_0); 
        } // ll 
        for (ll = 0; ll < pDet; ll++) {
          TauAlphaInv[ll * pDet + ll] = tauSqAlpha[ll]; 
        } // ll 
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
                                betaStarSites[i * J + j] + wStar[j * nSp + i]);
	      if (family == 1) {
                logPostBetaCand += nb_logpost(kappa[i], N[j * nSp + i], 
                                              tmp_J[j], offset[j]);
	      } else {
                logPostBetaCand += poisson_logpost(N[j * nSp + i], tmp_J[j], offset[j]);
	      }
              tmp_J[j] = exp(F77_NAME(ddot)(&pAbund, &X[j], &J, &beta[i], &nSp) + 
                                betaStarSites[i * J + j] + wStar[j * nSp + i]);
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
	  // alpha is ordered by parameter, then species
          // Proposal
	  for (l = 0; l < pDet; l++) {
            logPostAlphaCand = 0.0;
	    logPostAlphaCurr = 0.0;
            alphaCand[l * nSp + i] = rnorm(alpha[l * nSp + i], 
                                 exp(tuning[alphaAMCMCIndx + l * nSp + i]));
            logPostAlphaCand += dnorm(alphaCand[l * nSp + i], alphaComm[l], 
			              sqrt(tauSqAlpha[l]), 1);
	    logPostAlphaCurr += dnorm(alpha[l * nSp + i], alphaComm[l], 
                                      sqrt(tauSqAlpha[l]), 1);
	    for (j = 0; j < J; j++) {
              /********************************
               * Candidate 
               *******************************/
              tmp_0 = 0.0; 
	      likeVal = 0.0;
              sigma[j] = exp(F77_NAME(ddot)(&pDet, &Xp[j], &J, &alphaCand[i], &nSp) +
                             alphaStarSites[i * J + j]);
              for (k = 0; k < K; k++) {
                p[k * J + j] = integrate(detModel, distBreaks[k], distBreaks[k + 1], sigma[j], 
                                         nInt, transect); 
	        if (transect == 0) {
                  p[k * J + j] /= (distBreaks[k + 1] - distBreaks[k]);
	        } else {
                  p[k * J + j] = p[k * J + j] * 2.0 / (pow(distBreaks[k + 1], 2) - pow(distBreaks[k], 2));
	        }
	        piFullCand[k * JnSp + j * nSp + i] = p[k * J + j] * psi[k];
	        tmp_0 += piFullCand[k * JnSp + j * nSp + i];
	        likeVal += y[k * JnSp + j * nSp + i] * log(piFullCand[k * JnSp + j * nSp + i]);
	      } // k (bins)
	      piFullCand[k * JnSp + j * nSp + i] = 1.0 - tmp_0;
	      likeVal += (N[j * nSp + i] - yMax[j * nSp + i]) * log(piFullCand[K * JnSp + j * nSp + i]);
	      logPostAlphaCand += likeVal;
              /********************************
               * Current 
               *******************************/
              tmp_0 = 0.0; 
	      likeVal = 0.0;
              sigma[j] = exp(F77_NAME(ddot)(&pDet, &Xp[j], &J, &alpha[i], &nSp) + 
                             alphaStarSites[i * J + j]);
              for (k = 0; k < K; k++) {
                p[k * J + j] = integrate(detModel, distBreaks[k], distBreaks[k + 1], sigma[j], 
                                         nInt, transect); 
	        if (transect == 0) {
                  p[k * J + j] /= (distBreaks[k + 1] - distBreaks[k]);
	        } else {
                  p[k * J + j] = p[k * J + j] * 2.0 / (pow(distBreaks[k + 1], 2) - pow(distBreaks[k], 2));
	        }
	        piFull[k * JnSp + j * nSp + i] = p[k * J + j] * psi[k];
	        tmp_0 += piFull[k * JnSp + j * nSp + i];
	        likeVal += y[k * JnSp + j * nSp + i] * log(piFull[k * JnSp + j * nSp + i]);
	      } // k (bins)
	      piFull[k * JnSp + j * nSp + i] = 1.0 - tmp_0;
	      likeVal += (N[j * nSp + i] - yMax[j * nSp + i]) * log(piFull[K * JnSp + j * nSp + i]);
	      logPostAlphaCurr += likeVal;
	    } // j (sites)
            if (runif(0.0, 1.0) <= exp(logPostAlphaCand - logPostAlphaCurr)) {
              alpha[l * nSp + i] = alphaCand[l * nSp + i];
              accept[alphaAMCMCIndx + l * nSp + i]++;
	      F77_NAME(dcopy)(&nObsFullnSp, piFullCand, &inc, piFull, &inc);
            } else {
              alphaCand[l * nSp + i] = alpha[l * nSp + i];
	      F77_NAME(dcopy)(&nObsFullnSp, piFull, &inc, piFullCand, &inc);
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
	  			  betaStarSitesCand[i * J + j] + wStar[j * nSp + i]);
                  if (family == 1) {
	            logPostBetaStarCand[l] += nb_logpost(kappa[i], N[j * nSp + i], 
                                                         tmp_J[j], offset[j]); 	         
	  	  } else {
	            logPostBetaStarCand[l] += poisson_logpost(N[j * nSp + i], tmp_J[j],
                                                              offset[j]); 
	  	  }
	  	  // Current
                    betaStarSites[i * J + j] = 0.0;
                    for (ll = 0; ll < pAbundRE; ll++) {
                      betaStarSites[i * J + j] += betaStar[i * nAbundRE + betaStarLongIndx[ll * J + j]] * 
	                                  XRandom[ll * J + j];
                    }
                    tmp_J[j] = exp(F77_NAME(ddot)(&pAbund, &X[j], &J, &beta[i], &nSp) + 
	  	  		  betaStarSites[i * J + j] + wStar[j * nSp + i]);
                  if (family == 1) {
	            logPostBetaStarCurr[l] += nb_logpost(kappa[i], N[j * nSp + i], 
                                                         tmp_J[j], offset[j]); 	         
	  	  } else {
	            logPostBetaStarCurr[l] += poisson_logpost(N[j * nSp + i], tmp_J[j], 
				                              offset[j]); 
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
	      for (j = 0; j < J; j++) {
                likeVal = 0.0;
                tmp_0 = 0.0; 
                if ((XpRE[alphaStarIndx[l] * J + j] == alphaLevelIndx[l])) {
                  /********************************
                   * Candidate 
                   *******************************/
                  alphaStarSitesCand[i * J + j] = 0.0;
                  for (ll = 0; ll < pDetRE; ll++) {
                    alphaStarSitesCand[i * J + j] += alphaStarCand[i * nDetRE + alphaStarLongIndx[ll * J + j]] * 
	                                XpRandom[ll * J + j];
                  }
                  sigma[j] = exp(F77_NAME(ddot)(&pDet, &Xp[j], &J, &alpha[i], &nSp) + 
                                 alphaStarSitesCand[i * J + j]);
                  for (k = 0; k < K; k++) {
                    p[k * J + j] = integrate(detModel, distBreaks[k], distBreaks[k + 1], sigma[j], 
                                             nInt, transect); 
	            if (transect == 0) {
                      p[k * J + j] /= (distBreaks[k + 1] - distBreaks[k]);
	            } else {
                      p[k * J + j] = p[k * J + j] * 2.0 / (pow(distBreaks[k + 1], 2) - pow(distBreaks[k], 2));
	            }
	            piFullCand[k * JnSp + j * nSp + i] = p[k * J + j] * psi[k];
	            tmp_0 += piFullCand[k * JnSp + j * nSp + i];
	            likeVal += y[k * JnSp + j * nSp + i] * log(piFullCand[k * JnSp + j * nSp + i]);
	          } // k (bins)
	          piFullCand[k * JnSp + j * nSp + i] = 1.0 - tmp_0;
	          likeVal += (N[j * nSp + i] - yMax[j * nSp + i]) * log(piFullCand[K * JnSp + j * nSp + i]);
	          logPostAlphaStarCand[l] += likeVal;
                  /********************************
                   * Current
                   *******************************/
	  	  likeVal = 0.0;
	  	  tmp_0 = 0.0;
                  alphaStarSites[i * J + j] = 0.0;
                  for (ll = 0; ll < pDetRE; ll++) {
                    alphaStarSites[i * J + j] += alphaStar[i * nDetRE + alphaStarLongIndx[ll * J + j]] * 
	                                XpRandom[ll * J + j];
                  }
                  sigma[j] = exp(F77_NAME(ddot)(&pDet, &Xp[j], &J, &alpha[i], &nSp) + 
                                 alphaStarSites[i * J + j]);
                  for (k = 0; k < K; k++) {
                    p[k * J + j] = integrate(detModel, distBreaks[k], distBreaks[k + 1], sigma[j], 
                                             nInt, transect); 
	            if (transect == 0) {
                      p[k * J + j] /= (distBreaks[k + 1] - distBreaks[k]);
	            } else {
                      p[k * J + j] = p[k * J + j] * 2.0 / (pow(distBreaks[k + 1], 2) - pow(distBreaks[k], 2));
	            }
	            piFull[k * JnSp + j * nSp + i] = p[k * J + j] * psi[k];
	            tmp_0 += piFull[k * JnSp + j * nSp + i];
	            likeVal += y[k * JnSp + j * nSp + i] * log(piFull[k * JnSp + j * nSp + i]);
	          } // k (bins)
	          piFull[k * JnSp + j * nSp + i] = 1.0 - tmp_0;
	          likeVal += (N[j * nSp + i] - yMax[j * nSp + i]) * log(piFull[K * JnSp + j * nSp + i]);
	  	  logPostAlphaStarCurr[l] += likeVal;
	        }
	      }
	      if (runif (0.0, 1.0) <= exp(logPostAlphaStarCand[l] - logPostAlphaStarCurr[l])) {
                alphaStar[i * nDetRE + l] = alphaStarCand[i * nDetRE + l];
	        F77_NAME(dcopy)(&JnSp, alphaStarSitesCand, &inc, alphaStarSites, &inc);
	        accept[alphaStarAMCMCIndx + i * nDetRE + l]++;
	        F77_NAME(dcopy)(&nObsFullnSp, piFullCand, &inc, piFull, &inc);
	      } else {
                alphaStarCand[i * nDetRE + l] = alphaStar[i * nDetRE + l];
	        F77_NAME(dcopy)(&JnSp, alphaStarSites, &inc, alphaStarSitesCand, &inc);
	        F77_NAME(dcopy)(&nObsFullnSp, piFull, &inc, piFullCand, &inc);
	      }
	    }
	  }
	} // i (species)

        /********************************************************************
         *Update Spatial Random Effects (w)
         *******************************************************************/
	// Update w for each factor and each location. 
	for (ll = 0; ll < q; ll++) {
	  for (j = 0; j < J; j++) {
            // Propose new value
	    logPostWCand[j * q + ll] = 0.0;
	    wCand[j * q + ll] = rnorm(w[j * q + ll], exp(tuning[wAMCMCIndx + j * q + ll]));
            /*****************************
	     *Candidate
             *****************************/
            a[ll] = 0; 
	    v[ll] = 0; 
	    // MVN for any neighbors of j
	    if (uIndxLU[J + j] > 0){ // is j is a neighbor for anybody
	      for (r = 0; r < uIndxLU[J+j]; r++){ // how many locations have j as a neighbor
	        // now the neighbors for the jth location who has ii as a neighbor
	        jj = uIndx[uIndxLU[j]+r]; // jj is the index of the rth location who has j as a neighbor
                e = 0;
	        for (k = 0; k < nnIndxLU[J + jj]; k++) {
                  e += B[ll * nIndx + nnIndxLU[jj] + k] * wCand[nnIndx[nnIndxLU[jj] + k] * q + ll];
	        }
	        b = wCand[jj * q + ll] - e;
	        a[ll] += b * b / F[ll * J + jj];
	      } // r
	    }
	    // MVN for j
	    if (nnIndxLU[J + j] > 0) { // if j has any neighbors.
              e = 0; 
	      for (r = 0; r < nnIndxLU[J + j]; r++) {
                e += B[ll * nIndx + nnIndxLU[j] + r] * wCand[nnIndx[nnIndxLU[j] + r] * q + ll];
	      } // r 
              b = wCand[j * q + ll] - e;
	    } else {
              b = wCand[j * q + ll];
	    }
	    a[ll] += b * b / F[ll * J + j];
	    logPostWCand[j * q + ll] = -0.5 * a[ll];

            /*****************************
	     *Current
             *****************************/
            a[ll] = 0; 
	    v[ll] = 0; 
	    // MVN for any neighbors of j
	    if (uIndxLU[J + j] > 0){ // is j is a neighbor for anybody
	      for (r = 0; r < uIndxLU[J+j]; r++){ // how many locations have j as a neighbor
	        // now the neighbors for the jth location who has ii as a neighbor
	        jj = uIndx[uIndxLU[j]+r]; // jj is the index of the rth location who has j as a neighbor
                e = 0;
	        for (k = 0; k < nnIndxLU[J + jj]; k++) { // these are the neighbors of jj
                  e += B[ll * nIndx + nnIndxLU[jj] + k] * w[nnIndx[nnIndxLU[jj] + k] * q + ll];
	        }
	        b = w[jj * q + ll] - e;
	        a[ll] += b * b / F[ll * J + jj];
	      } // r
	    }
	    // MVN for j
	    if (nnIndxLU[J + j] > 0) { // if j has any neighbors.
              e = 0; 
	      for (r = 0; r < nnIndxLU[J + j]; r++) {
                e += B[ll * nIndx + nnIndxLU[j] + r] * w[nnIndx[nnIndxLU[j] + r] * q + ll];
	      } // r 
              b = w[j * q + ll] - e;
	    } else {
              b = w[j * q + ll];
	    }
	    a[ll] += b * b / F[ll * J + j];
	    logPostWCurr[j * q + ll] = -0.5 * a[ll];

	    // Likelihood for proposal
	    for (i = 0; i < nSp; i++) {
	      wStarCand[j * nSp + i] = 0.0;
	      for (ii = 0; ii < q; ii++) {
                wStarCand[j * nSp + i] += lambda[ii * nSp + i] * wCand[j * q + ii];
	      } // ii 
              /*****************************
	       *Candidate
               *****************************/
              tmp_J[j] = exp(F77_NAME(ddot)(&pAbund, &X[j], &J, &beta[i], &nSp) + 
                                betaStarSites[i * J + j] + wStarCand[j * nSp + i]);
	      if (family == 1) {
		logPostWCand[j * q + ll] += nb_logpost(kappa[i], N[j * nSp + i], tmp_J[j], 
				                       offset[j]);
	      } else {
		logPostWCand[j * q + ll] += poisson_logpost(N[j * nSp + i], tmp_J[j], 
                                                            offset[j]);
	      }
              /*****************************
	       *Current
               *****************************/
              tmp_J[j] = exp(F77_NAME(ddot)(&pAbund, &X[j], &J, &beta[i], &nSp) + 
                                betaStarSites[i * J + j] + wStar[j * nSp + i]);
	      if (family == 1) {
		logPostWCurr[j * q + ll] += nb_logpost(kappa[i], N[j * nSp + i], tmp_J[j], 
				                       offset[j]);
	      } else {
		logPostWCurr[j * q + ll] += poisson_logpost(N[j * nSp + i], tmp_J[j], 
                                                            offset[j]);
	      }
	    } // i
	    if (runif(0.0, 1.0) <= exp(logPostWCand[j * q + ll] - logPostWCurr[j * q + ll])) {
              w[j * q + ll] = wCand[j * q + ll];
	      for (i = 0; i < nSp; i++) {
	        wStar[j * nSp + i] = wStarCand[j * nSp + i];
	      }
	      accept[wAMCMCIndx + j * q + ll]++;
	    } else {
              // Reset everything back to what it was.
              for (i = 0; i < nSp; i++) {
                wStarCand[j * nSp + i] = wStar[j * nSp + i];
	      }
	      wCand[j * q + ll] = w[j * q + ll];
	    }
	  } // j
        } // ll 

        /********************************************************************
         *Update spatial factors (lambda)
         *******************************************************************/
	for (i = 0; i < nSp; i++) {
          for (ll = 0; ll < q; ll++) {
            if (ll < i) { // only update lower triangle
	      lambdaCand[ll * nSp + i] = rnorm(lambda[ll * nSp + i],
                                               exp(tuning[lambdaAMCMCIndx + ll * nSp + i]));
	      logPostLambdaCand[ll * nSp + i] = dnorm(lambdaCand[ll * nSp + i], 0.0, 1.0, 1);
              logPostLambdaCurr[ll * nSp + i] = dnorm(lambda[ll * nSp + i], 0.0, 1.0, 1);
	      for (j = 0; j < J; j++) {
                wStarCand[j * nSp + i] = 0.0;
		for (ii = 0; ii < q; ii++) {
                  wStarCand[j * nSp + i] += lambdaCand[ii * nSp + i] * w[j * q + ii];
		} // ii
                // Candidate
                tmp_J[j] = exp(F77_NAME(ddot)(&pAbund, &X[j], &J, &beta[i], &nSp) + 
                                  betaStarSites[i * J + j] + wStarCand[j * nSp + i]);
	        if (family == 1) {
                  logPostLambdaCand[ll * nSp + i] += nb_logpost(kappa[i], N[j * nSp + i],  
	      		                                        tmp_J[j], offset[j]);
	        } else {
                  logPostLambdaCand[ll * nSp + i] += poisson_logpost(N[j * nSp + i], tmp_J[j], 
                                                                     offset[j]);
	        }
		// Current
                tmp_J[j] = exp(F77_NAME(ddot)(&pAbund, &X[j], &J, &beta[i], &nSp) + 
                               betaStarSites[i * J + j] + wStar[j * nSp + i]);
	        if (family == 1) {
                  logPostLambdaCurr[ll * nSp + i] += nb_logpost(kappa[i], N[j * nSp + i], 
	      		                                        tmp_J[j], offset[j]);
	        } else {
                  logPostLambdaCurr[ll * nSp + i] += poisson_logpost(N[j * nSp + i], tmp_J[j], 
                                                                     offset[j]);
	        }
	      } // j
              if (runif(0.0, 1.0) <= exp(logPostLambdaCand[ll * nSp + i] - 
				         logPostLambdaCurr[ll * nSp + i])) {
                lambda[ll * nSp + i] = lambdaCand[ll * nSp + i];
		for (j = 0; j < J; j++) {
                  wStar[j * nSp + i] = wStarCand[j * nSp + i];
		}
		accept[lambdaAMCMCIndx + ll * nSp + i]++;
	      } else {
                // Reset everything back to what it was
		for (j = 0; j < J; j++) {
                  wStarCand[j * nSp + i] = wStar[j * nSp + i];
		}
                lambdaCand[ll * nSp + i] = lambda[ll * nSp + i];		
	      }
	    }
	  } // ll (factor) 
	} // i (species)

        /********************************************************************
         *Update phi (and nu if matern)
         *******************************************************************/
	for (ll = 0; ll < q; ll++) {
          // Current
          if (corName == "matern"){ 
	    nu[ll] = theta[nuIndx * q + ll];
       	  }
          updateBFSFDS(&B[ll * nIndx], &F[ll*J], &c[ll * m*nThreads], &C[ll * mm * nThreads], coords, nnIndx, nnIndxLU, J, m, theta[sigmaSqIndx * q + ll], theta[phiIndx * q + ll], nu[ll], covModel, &bk[ll * sizeBK], nuB[ll]);
          
	  aa = 0;
          logDet = 0;

#ifdef _OPENMP
#pragma omp parallel for private (e, ii, b) reduction(+:aa, logDet)
#endif
          for (j = 0; j < J; j++){
            if (nnIndxLU[J+j] > 0){
              e = 0;
              for (ii = 0; ii < nnIndxLU[J+j]; ii++){
                e += B[ll * nIndx + nnIndxLU[j]+ii]*w[nnIndx[nnIndxLU[j]+ii] * q + ll];
              }
              b = w[j * q + ll] - e;
            } else{
              b = w[j * q + ll];
            }	
            aa += b*b/F[ll * J + j];
            logDet += log(F[ll * J + j]);
          }
      
          logPostThetaCurr[ll] = -0.5 * logDet - 0.5 * aa;
          logPostThetaCurr[ll] += log(theta[phiIndx * q + ll] - phiA[ll]) + log(phiB[ll] - theta[phiIndx * q + ll]); 
          if(corName == "matern"){
       	    logPostThetaCurr[ll] += log(theta[nuIndx * q + ll] - nuA[ll]) + log(nuB[ll] - theta[nuIndx * q + ll]); 
          }
          
          // Candidate
          phiCand[ll] = logitInv(rnorm(logit(theta[phiIndx * q + ll], phiA[ll], phiB[ll]), exp(tuning[phiAMCMCIndx + ll])), phiA[ll], phiB[ll]);
          if (corName == "matern"){
      	    nuCand[ll] = logitInv(rnorm(logit(theta[nuIndx * q + ll], nuA[ll], nuB[ll]), exp(tuning[nuIndx * q + ll])), nuA[ll], nuB[ll]);
          }
      
          updateBFSFDS(BCand, FCand, &c[ll * m*nThreads], &C[ll * mm * nThreads], coords, nnIndx, nnIndxLU, J, m, theta[sigmaSqIndx * q + ll], phiCand[ll], nuCand[ll], covModel, &bk[ll * sizeBK], nuB[ll]);
      
          aa = 0;
          logDet = 0;
      
#ifdef _OPENMP
#pragma omp parallel for private (e, ii, b) reduction(+:aa, logDet)
#endif
          for (j = 0; j < J; j++){
            if (nnIndxLU[J+j] > 0){
              e = 0;
              for (ii = 0; ii < nnIndxLU[J+j]; ii++){
                e += BCand[nnIndxLU[j]+ii]*w[nnIndx[nnIndxLU[j]+ii] * q + ll];
              }
              b = w[j * q + ll] - e;
            } else{
              b = w[j * q + ll];
              }	
              aa += b*b/FCand[j];
              logDet += log(FCand[j]);
          }
          
          logPostThetaCand[ll] = -0.5*logDet - 0.5*aa;      
          logPostThetaCand[ll] += log(phiCand[ll] - phiA[ll]) + log(phiB[ll] - phiCand[ll]); 
          if (corName == "matern"){
            logPostThetaCand[ll] += log(nuCand[ll] - nuA[ll]) + log(nuB[ll] - nuCand[ll]); 
          }

          if (runif(0.0,1.0) <= exp(logPostThetaCand[ll] - logPostThetaCurr[ll])) {

            F77_NAME(dcopy)(&nIndx, BCand, &inc, &B[ll * nIndx], &inc);
            F77_NAME(dcopy)(&J, FCand, &inc, &F[ll * J], &inc);
            
	    theta[phiIndx * q + ll] = phiCand[ll];
            accept[phiAMCMCIndx + ll]++;
            if (corName == "matern") {
              nu[ll] = nuCand[ll]; 
	      theta[nuIndx * q + ll] = nu[ll]; 
              accept[nuAMCMCIndx + ll]++; 
            }
          }
	} // ll

	for (i = 0; i < nSp; i++) {
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
                          betaStarSites[i * J + j] + wStar[j * nSp + i]);
	      logPostKappaCurr += nb_logpost(kappa[i], N[j * nSp + i], mu[j * nSp + i],
			                     offset[j]);
	      logPostKappaCand += nb_logpost(kappaCand, N[j * nSp + i], mu[j * nSp + i], 
			                     offset[j]);
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
            NCand[j] = rpois(N[j * nSp + i] + epsilon);
	    // Only calculate if Poisson since its already calculated in kappa update
	    if (family == 0) {
              mu[j * nSp + i] = exp(F77_NAME(ddot)(&pAbund, &X[j], &J, &beta[i], &nSp) + 
                                    betaStarSites[i * J + j] + wStar[j * nSp + i]);
	    }
	    // Likelihood contribution
	    if (NCand[j] >= yMax[j * nSp + i]) {
	      likeVal = 0.0;
	      tmp_0 = 0.0;
	      for (k = 0; k < K; k++) {
                likeVal += y[k * JnSp + j * nSp + i] * log(piFull[k * JnSp + j * nSp + i]);  
	      }
              /********************************
               * Current
               *******************************/
	      logPostCurrN[j] += likeVal + 
                                 ((N[j * nSp + i] - yMax[j * nSp + i]) * log(piFull[K * JnSp + j * nSp + i])) + 
                                 lgammafn(N[j * nSp + i] + 1.0) -
	  		       lgammafn(N[j * nSp + i] - yMax[j * nSp + i] + 1.0);
	      // Contribution from NB or Poisson
	      if (family == 0) {
                logPostCurrN[j] += poisson_logpost(N[j * nSp + i], mu[j * nSp + i], offset[j]);
	      } else {
	        logPostCurrN[j] += nb_logpost(kappa[i], N[j * nSp + i], mu[j * nSp + i], offset[j]);
	      }
	      // MH contribution for assymetric proposal distribution.
	      logPostCurrN[j] += poisson_logpost(NCand[j], N[j * nSp + i] + epsilon, 1.0);
              /********************************
               * Candidate
               *******************************/
	      logPostCandN[j] += likeVal + 
                                 ((NCand[j] - yMax[j * nSp + i]) * log(piFull[K * JnSp + j * nSp + i])) + 
                                 lgammafn(NCand[j] + 1.0) -
	  		       lgammafn(NCand[j] - yMax[j * nSp + i] + 1.0);
	      if (family == 0) {
                logPostCandN[j] += poisson_logpost(NCand[j], mu[j * nSp + i], offset[j]);
	      } else {
	        logPostCandN[j] += nb_logpost(kappa[i], NCand[j], mu[j * nSp + i], offset[j]);
	      }
	      // MH contribution for assymetric proposal distribution.
	      logPostCandN[j] += poisson_logpost(N[j * nSp + i], NCand[j] + epsilon, 1.0);
              if (runif(0.0,1.0) <= exp(logPostCandN[j] - logPostCurrN[j])) {
                N[j * nSp + i] = NCand[j];
              }
	    }
	  }
	}

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
            F77_NAME(dcopy)(&Jq, w, &inc, &REAL(wSamples_r)[sPost*Jq], &inc); 
            F77_NAME(dcopy)(&nSpq, lambda, &inc, &REAL(lambdaSamples_r)[sPost*nSpq], &inc); 
	    if (family == 1) {
              F77_NAME(dcopy)(&nSp, kappa, &inc, &REAL(kappaSamples_r)[sPost*nSp], &inc); 
	    }
            F77_NAME(dcopy)(&JnSp, N, &inc, &REAL(NSamples_r)[sPost*JnSp], &inc); 
            F77_NAME(dcopy)(&JnSp, mu, &inc, &REAL(muSamples_r)[sPost*JnSp], &inc); 
            F77_NAME(dcopy)(&nThetaqSave, &theta[phiIndx * q], &inc, 
			    &REAL(thetaSamples_r)[sPost*nThetaqSave], &inc); 
	    // Generate and save replicate values. 
	    for (i = 0; i < nSp; i++) {
	      for (j = 0; j < J; j++) { 
                for (k = 0; k < KFull; k++) {
                  tmp_KFull[k] = piFull[k * JnSp + j * nSp + i];
	          REAL(piFullSamples_r)[sPost * nObsFullnSp + k * JnSp + j * nSp + i] = piFull[k * JnSp + j * nSp + i];
	        } 
	        if (family == 0) {
                  tmp_0 = rpois(mu[j * nSp + i] * offset[j]);
	        } else {
                  tmp_0 = rnbinom_mu(mu[j * nSp + i] * offset[j], kappa[i]);
	        }
	        rmultinom(tmp_0, tmp_KFull, KFull, tmp_KFullInt);
		for (k = 0; k < KFull; k++) {
                  INTEGER(yRepSamples_r)[sPost * nObsFullnSp + k * JnSp + j * nSp + i] = tmp_KFullInt[k];
		}
	      }
	    }
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
          Rprintf("\tNumber\tParameter\tAcceptance\tTuning\n");	  
          for (ll = 0; ll < nSp; ll++) {
            Rprintf("\t%i\tbeta[1]\t\t%3.1f\t\t%1.5f\n", ll + 1, 100.0*REAL(acceptSamples_r)[s * nAMCMC + betaAMCMCIndx + ll], exp(tuning[betaAMCMCIndx + ll]));
            Rprintf("\t%i\talpha[1]\t%3.1f\t\t%1.5f\n", ll + 1, 100.0*REAL(acceptSamples_r)[s * nAMCMC + alphaAMCMCIndx + ll], exp(tuning[alphaAMCMCIndx + ll]));
          } // ll
          for (ll = 0; ll < q; ll++) {
            Rprintf("\t%i\tphi\t\t%3.1f\t\t%1.5f\n", ll + 1, 100.0*REAL(acceptSamples_r)[s * nAMCMC + phiAMCMCIndx + ll], exp(tuning[phiAMCMCIndx + ll]));
	    if (corName == "matern") {
            Rprintf("\t%i\tnu\t\t%3.1f\t\t%1.5f\n", ll + 1, 100.0*REAL(acceptSamples_r)[s * nAMCMC + nuAMCMCIndx + ll], exp(tuning[nuAMCMCIndx + ll]));
	    }
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
    int nResultListObjs = 14;
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
    SET_VECTOR_ELT(result_r, 0, betaCommSamples_r);
    SET_VECTOR_ELT(result_r, 1, alphaCommSamples_r);
    SET_VECTOR_ELT(result_r, 2, tauSqBetaSamples_r);
    SET_VECTOR_ELT(result_r, 3, tauSqAlphaSamples_r);
    SET_VECTOR_ELT(result_r, 4, betaSamples_r);
    SET_VECTOR_ELT(result_r, 5, alphaSamples_r);
    SET_VECTOR_ELT(result_r, 6, NSamples_r); 
    SET_VECTOR_ELT(result_r, 7, muSamples_r);
    SET_VECTOR_ELT(result_r, 8, tuningSamples_r);
    SET_VECTOR_ELT(result_r, 9, yRepSamples_r);
    SET_VECTOR_ELT(result_r, 10, piFullSamples_r);
    SET_VECTOR_ELT(result_r, 11, lambdaSamples_r);
    SET_VECTOR_ELT(result_r, 12, wSamples_r); 
    SET_VECTOR_ELT(result_r, 13, thetaSamples_r);
    if (pDetRE > 0) {
      tmp_0 = 14; // Needed to make tracking kappa easier. 
      SET_VECTOR_ELT(result_r, 14, sigmaSqPSamples_r);
      SET_VECTOR_ELT(result_r, 15, alphaStarSamples_r);
    }
    if (pAbundRE > 0) {
      if (pDetRE > 0) {
        tmp_0 = 16;
      } else {
        tmp_0 = 14;
      }
      SET_VECTOR_ELT(result_r, tmp_0, sigmaSqMuSamples_r);
      SET_VECTOR_ELT(result_r, tmp_0 + 1, betaStarSamples_r);
    }
    if (family == 1) {
      if ((pDetRE > 0) || (pAbundRE > 0)) {
        tmp_02 = tmp_0 + 2;
      } else {
        tmp_02 = 14;
      }  
      SET_VECTOR_ELT(result_r, tmp_02, kappaSamples_r);
    }

    SET_VECTOR_ELT(resultName_r, 0, mkChar("beta.comm.samples")); 
    SET_VECTOR_ELT(resultName_r, 1, mkChar("alpha.comm.samples")); 
    SET_VECTOR_ELT(resultName_r, 2, mkChar("tau.sq.beta.samples")); 
    SET_VECTOR_ELT(resultName_r, 3, mkChar("tau.sq.alpha.samples")); 
    SET_VECTOR_ELT(resultName_r, 4, mkChar("beta.samples")); 
    SET_VECTOR_ELT(resultName_r, 5, mkChar("alpha.samples")); 
    SET_VECTOR_ELT(resultName_r, 6, mkChar("N.samples")); 
    SET_VECTOR_ELT(resultName_r, 7, mkChar("mu.samples"));
    SET_VECTOR_ELT(resultName_r, 8, mkChar("tuning.samples"));
    SET_VECTOR_ELT(resultName_r, 9, mkChar("y.rep.samples"));
    SET_VECTOR_ELT(resultName_r, 10, mkChar("pi.samples"));
    SET_VECTOR_ELT(resultName_r, 11, mkChar("lambda.samples")); 
    SET_VECTOR_ELT(resultName_r, 12, mkChar("w.samples")); 
    SET_VECTOR_ELT(resultName_r, 13, mkChar("theta.samples"));
    if (pDetRE > 0) {
      SET_VECTOR_ELT(resultName_r, 14, mkChar("sigma.sq.p.samples")); 
      SET_VECTOR_ELT(resultName_r, 15, mkChar("alpha.star.samples")); 
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

