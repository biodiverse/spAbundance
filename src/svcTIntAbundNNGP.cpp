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

void updateBFsvcTIntAbund(double *B, double *F, double *c, double *C, double *coords, int *nnIndx, int *nnIndxLU, int n, int m, double sigmaSq, double phi, double nu, int covModel, double *bk, double nuUnifb){

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
	F77_NAME(dpotrf)(&lower, &nnIndxLU[n+i], &C[mm*threadID], &nnIndxLU[n+i], &info FCONE); if(info != 0){Rf_error("c++ error: dpotrf failed\n");}
	F77_NAME(dpotri)(&lower, &nnIndxLU[n+i], &C[mm*threadID], &nnIndxLU[n+i], &info FCONE); if(info != 0){Rf_error("c++ error: dpotri failed\n");}
	F77_NAME(dsymv)(&lower, &nnIndxLU[n+i], &one, &C[mm*threadID], &nnIndxLU[n+i], &c[m*threadID], &inc, &zero, &B[nnIndxLU[i]], &inc FCONE);
	F[i] = sigmaSq - F77_NAME(ddot)(&nnIndxLU[n+i], &B[nnIndxLU[i]], &inc, &c[m*threadID], &inc);
      }else{
	B[i] = 0;
	F[i] = sigmaSq;
      }
    }

}

extern "C" {
  SEXP svcTIntAbundNNGP(SEXP y_r, SEXP X_r, SEXP Xw_r, SEXP Xp_r, SEXP coords_r, SEXP XRE_r,
                        SEXP XpRE_r, SEXP consts_r, SEXP pDetLong_r, 
                        SEXP JLong_r, SEXP nObsLong_r, SEXP nAbundRELong_r, SEXP nDetRELong_r,
                        SEXP nnIndx_r, SEXP nnIndxLU_r, SEXP uIndx_r, SEXP uIndxLU_r, 
                        SEXP betaStarting_r, SEXP alphaStarting_r, 
                        SEXP sigmaSqMuStarting_r, SEXP sigmaSqPStarting_r, 
                        SEXP betaStarStarting_r, SEXP alphaStarStarting_r,  
                        SEXP phiStarting_r, SEXP sigmaSqStarting_r, SEXP nuStarting_r,
                        SEXP wStarting_r, SEXP kappaStarting_r,  
                        SEXP zLongIndx_r, SEXP dataIndx_r, SEXP alphaIndx_r, SEXP zYearIndx_r, 
                        SEXP zDatIndx_r, SEXP zSiteIndx_r, SEXP siteIndx_r, 
                        SEXP betaStarIndx_r, SEXP betaLevelIndx_r,
                        SEXP alphaStarIndx_r, SEXP alphaLevelIndx_r, SEXP alphaNREIndx_r,
                        SEXP alphaColIndx_r, SEXP muBeta_r, SEXP SigmaBeta_r,
                        SEXP muAlpha_r, SEXP sigmaAlpha_r,
                        SEXP phiA_r, SEXP phiB_r,
                        SEXP sigmaSqA_r, SEXP sigmaSqB_r, SEXP nuA_r, SEXP nuB_r, 
                        SEXP sigmaSqMuA_r, SEXP sigmaSqMuB_r, SEXP sigmaSqPA_r, 
                        SEXP sigmaSqPB_r, SEXP kappaA_r, SEXP kappaB_r, 
                        SEXP tuning_r, SEXP acceptRate_r, SEXP chainInfo_r, 
                        SEXP waicNObsIndx_r, SEXP waicCellIndx_r, SEXP offset_r){

    /**********************************************************************
     * Initial constants
     * *******************************************************************/
    int i, g, t, j, k, jj, s, r, l, ll, rr, nProtect=0;
    const int inc = 1;

    /**********************************************************************
     * Get Inputs
     * *******************************************************************/
    double *y = REAL(y_r);
    double *X = REAL(X_r);
    double *Xw = REAL(Xw_r);
    double *Xp = REAL(Xp_r);
    double *offset = REAL(offset_r);
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
    int nYearsMax = INTEGER(consts_r)[8];
    int nData = INTEGER(consts_r)[9];
    int nWAIC = INTEGER(consts_r)[10];
    int sigmaSqIG = INTEGER(consts_r)[11];
    int covModel = INTEGER(consts_r)[12];
    std::string corName = getCorName(covModel);
    int m = INTEGER(consts_r)[13]; 
    int nBatch = INTEGER(consts_r)[14]; 
    int batchLength = INTEGER(consts_r)[15]; 
    int nSamples = nBatch * batchLength; 
    int nThreads = INTEGER(consts_r)[16];
    int verbose = INTEGER(consts_r)[17];
    int nReport = INTEGER(consts_r)[18];
    int nThin = INTEGER(consts_r)[19]; 
    int nBurn = INTEGER(consts_r)[20]; 
    int nPost = INTEGER(consts_r)[21]; 
    int pTilde = INTEGER(consts_r)[22]; 
    int *pDetLong = INTEGER(pDetLong_r); 
    int saveFitted = INTEGER(consts_r)[23];
    int anyFamily = INTEGER(consts_r)[24];
    int *family = (int *) R_alloc(nData, sizeof(int));
    for (j = 0; j < nData; j++) {
      family[j] = INTEGER(consts_r)[25 + j];
    }
    int JpTilde = J * pTilde;
    int ppAbund = pAbund * pAbund;
    double *muBeta = (double *) R_alloc(pAbund, sizeof(double));
    F77_NAME(dcopy)(&pAbund, REAL(muBeta_r), &inc, muBeta, &inc);
    double *SigmaBeta = (double *) R_alloc(ppAbund, sizeof(double));
    F77_NAME(dcopy)(&ppAbund, REAL(SigmaBeta_r), &inc, SigmaBeta, &inc);
    double *muAlpha = (double *) R_alloc(pDet, sizeof(double));   
    F77_NAME(dcopy)(&pDet, REAL(muAlpha_r), &inc, muAlpha, &inc);
    double *sigmaAlpha = (double *) R_alloc(pDet, sizeof(double)); 
    F77_NAME(dcopy)(&pDet, REAL(sigmaAlpha_r), &inc, sigmaAlpha, &inc);
    double *kappaA = REAL(kappaA_r);
    double *kappaB = REAL(kappaB_r);
    double *sigmaSqMuA = REAL(sigmaSqMuA_r);
    double *sigmaSqMuB = REAL(sigmaSqMuB_r);
    double *sigmaSqPA = REAL(sigmaSqPA_r); 
    double *sigmaSqPB = REAL(sigmaSqPB_r); 
    double *phiA = REAL(phiA_r);
    double *phiB = REAL(phiB_r);
    double *nuA = REAL(nuA_r);
    double *nuB = REAL(nuB_r);
    double *sigmaSqA = REAL(sigmaSqA_r);
    double *sigmaSqB = REAL(sigmaSqB_r);
    int *JLong = INTEGER(JLong_r);
    // Total number of sites across all data sets, including
    // sites that are sampled by multiple data sources
    int JSum = 0; 
    for (i = 0; i < nData; i++) {
      JSum += JLong[i];
    }
    int *nObsLong = INTEGER(nObsLong_r); 
    double *tuning = REAL(tuning_r);
    double *coords = REAL(coords_r);
    int *nDetRELong = INTEGER(nDetRELong_r); 
    int *nAbundRELong = INTEGER(nAbundRELong_r);
    int *betaStarIndx = INTEGER(betaStarIndx_r);
    int *betaLevelIndx = INTEGER(betaLevelIndx_r);
    int *zLongIndx = INTEGER(zLongIndx_r); 
    int *dataIndx = INTEGER(dataIndx_r); 
    int *alphaIndx = INTEGER(alphaIndx_r); 
    int *zYearIndx = INTEGER(zYearIndx_r); 
    int *zDatIndx = INTEGER(zDatIndx_r); 
    int *siteIndx = INTEGER(siteIndx_r);
    int *zSiteIndx = INTEGER(zSiteIndx_r);
    int *alphaStarIndx = INTEGER(alphaStarIndx_r); 
    int *alphaLevelIndx = INTEGER(alphaLevelIndx_r);
    int *alphaNREIndx = INTEGER(alphaNREIndx_r);
    int *alphaColIndx = INTEGER(alphaColIndx_r);
    int *nnIndx = INTEGER(nnIndx_r);
    int *nnIndxLU = INTEGER(nnIndxLU_r);
    int *uIndx = INTEGER(uIndx_r);
    int *uIndxLU = INTEGER(uIndxLU_r);
    int currChain = INTEGER(chainInfo_r)[0];
    double acceptRate = REAL(acceptRate_r)[0];
    int nChain = INTEGER(chainInfo_r)[1];
    int status = 0;
    int thinIndx = 0; 
    int sPost = 0; 
    // For looping through data sets
    int stNObs = 0; 
    // For getting likelihood values for WAIC
    int *waicCellIndx = INTEGER(waicCellIndx_r);
    int *waicNObsIndx = INTEGER(waicNObsIndx_r);

    // More constant
    int JnYears = J * nYearsMax;

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
        Rprintf("SVC Integrated Multi-season Abundance model fit with %i sites and %i primary time periods.\n\n", J, nYearsMax);
        Rprintf("Integrating %i data sets.\n\n", nData); 
        Rprintf("Samples per Chain: %i (%i batches of length %i)\n", nSamples, nBatch, batchLength);
        Rprintf("Burn-in: %i \n", nBurn);
        Rprintf("Thinning Rate: %i \n", nThin);
        Rprintf("Number of Chains: %i \n", nChain);
        Rprintf("Total Posterior Samples: %i \n\n", nPost * nChain);
        Rprintf("Number of spatially-varying coefficients: %i \n", pTilde);
        Rprintf("Using the %s spatial correlation model.\n\n", corName.c_str());
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

    /**********************************************************************
     * Parameters
     * *******************************************************************/
    // Abundance Covariates
    double *beta = (double *) R_alloc(pAbund, sizeof(double));
    F77_NAME(dcopy)(&pAbund, REAL(betaStarting_r), &inc, beta, &inc);
    // Abundance random effect variances
    double *sigmaSqMu = (double *) R_alloc(pAbundRE, sizeof(double));
    F77_NAME(dcopy)(&pAbundRE, REAL(sigmaSqMuStarting_r), &inc, sigmaSqMu, &inc);
    // Latent abundance random effects
    double *betaStar = (double *) R_alloc(nAbundRE, sizeof(double));
    F77_NAME(dcopy)(&nAbundRE, REAL(betaStarStarting_r), &inc, betaStar, &inc);
    // Overdispersion parameter for NB for each data set;
    double *kappa = (double *) R_alloc(nData, sizeof(double));
    F77_NAME(dcopy)(&nData, REAL(kappaStarting_r), &inc, kappa, &inc);
    // Detection fixed effects
    double *alpha = (double *) R_alloc(pDet, sizeof(double));   
    F77_NAME(dcopy)(&pDet, REAL(alphaStarting_r), &inc, alpha, &inc);
    // Detection random effect variances
    double *sigmaSqP = (double *) R_alloc(pDetRE, sizeof(double)); 
    F77_NAME(dcopy)(&pDetRE, REAL(sigmaSqPStarting_r), &inc, sigmaSqP, &inc); 
    // Latent detection random effects
    double *alphaStar = (double *) R_alloc(nDetRE, sizeof(double)); 
    F77_NAME(dcopy)(&nDetRE, REAL(alphaStarStarting_r), &inc, alphaStar, &inc); 
    // Spatial parameters
    double *w = (double *) R_alloc(JpTilde, sizeof(double));
    F77_NAME(dcopy)(&JpTilde, REAL(wStarting_r), &inc, w, &inc);
    // Latent Abundance
    double *yRep = (double *) R_alloc(nObs, sizeof(double)); zeros(yRep, nObs);

    /**********************************************************************
     * Return Stuff
     * *******************************************************************/
    SEXP betaSamples_r;
    PROTECT(betaSamples_r = Rf_allocMatrix(REALSXP, pAbund, nPost)); nProtect++;
    zeros(REAL(betaSamples_r), pAbund * nPost);
    SEXP alphaSamples_r; 
    PROTECT(alphaSamples_r = Rf_allocMatrix(REALSXP, pDet, nPost)); nProtect++;
    zeros(REAL(alphaSamples_r), pDet * nPost);
    SEXP yRepSamples_r;
    SEXP lambdaSamples_r;
    SEXP likeSamples_r;
    if (saveFitted == 1) {
      PROTECT(yRepSamples_r = Rf_allocMatrix(REALSXP, nObs, nPost)); nProtect++;
      zeros(REAL(yRepSamples_r), nObs * nPost);
      PROTECT(lambdaSamples_r = Rf_allocMatrix(REALSXP, nObs, nPost)); nProtect++;
      zeros(REAL(lambdaSamples_r), nObs * nPost);
      PROTECT(likeSamples_r = Rf_allocMatrix(REALSXP, nObs, nPost)); nProtect++;
      zeros(REAL(likeSamples_r), nObs * nPost);
    }
    SEXP wSamples_r;
    PROTECT(wSamples_r = Rf_allocMatrix(REALSXP, JpTilde, nPost)); nProtect++;
    zeros(REAL(wSamples_r), JpTilde * nPost);
    SEXP kappaSamples_r;
    if (anyFamily == 1) {
      PROTECT(kappaSamples_r = Rf_allocMatrix(REALSXP, nData, nPost)); nProtect++;
      zeros(REAL(kappaSamples_r), nPost);
    }
    // Abundance random effects
    SEXP sigmaSqMuSamples_r;
    SEXP betaStarSamples_r;
    if (pAbundRE > 0) {
      PROTECT(sigmaSqMuSamples_r = Rf_allocMatrix(REALSXP, pAbundRE, nPost)); nProtect++;
      zeros(REAL(sigmaSqMuSamples_r), pAbundRE * nPost);
      PROTECT(betaStarSamples_r = Rf_allocMatrix(REALSXP, nAbundRE, nPost)); nProtect++;
      zeros(REAL(betaStarSamples_r), nAbundRE * nPost);
    }
    // Detection random effects
    SEXP sigmaSqPSamples_r; 
    SEXP alphaStarSamples_r; 
    if (pDetRE > 0) {
      PROTECT(sigmaSqPSamples_r = Rf_allocMatrix(REALSXP, pDetRE, nPost)); nProtect++;
      zeros(REAL(sigmaSqPSamples_r), pDetRE * nPost);
      PROTECT(alphaStarSamples_r = Rf_allocMatrix(REALSXP, nDetRE, nPost)); nProtect++;
      zeros(REAL(alphaStarSamples_r), nDetRE * nPost);
    }
    // Only the process values that disregard the observational effects
    SEXP muSamples_r; 
    PROTECT(muSamples_r = Rf_allocMatrix(REALSXP, JnYears, nPost)); nProtect++; 
    zeros(REAL(muSamples_r), JnYears * nPost);

    /********************************************************************
      Some constants and temporary variables to be used later
    ********************************************************************/
    double tmp_0, tmp_02;
    double *tmp_nObs = (double *) R_alloc(nObs, sizeof(double));

    // For latent abundance and WAIC
    double *like = (double *) R_alloc(nObs, sizeof(double)); zeros(like, nObs);
    double *mu = (double *) R_alloc(JnYears, sizeof(double));
    zeros(mu, JnYears);
    double *lambda = (double *) R_alloc(nObs, sizeof(double));
    zeros(lambda, nObs);
    int *stAlpha = (int *) R_alloc(nObs, sizeof(double));
    for (i = 0; i < nObs; i++) {
      stAlpha[i] = which(dataIndx[i], alphaIndx, pDet);
    }

    /**********************************************************************
     * Set up spatial stuff
     * *******************************************************************/
    int nTheta, sigmaSqIndx, phiIndx, nuIndx;
    if (corName != "matern") {
      nTheta = 2; // sigma^2, phi
      sigmaSqIndx = 0; phiIndx = 1;
    } else {
      nTheta = 3; // sigma^2, phi, nu
      sigmaSqIndx = 0; phiIndx = 1; nuIndx = 2;
    }
    int nThetapTilde = nTheta * pTilde;
    // theta is ordered by parameter, then SVC within parameter
    double *theta = (double *) R_alloc(nThetapTilde, sizeof(double));
    double *nu = (double *) R_alloc(pTilde, sizeof(double));
    SEXP thetaSamples_r;
    PROTECT(thetaSamples_r = Rf_allocMatrix(REALSXP, nThetapTilde, nPost)); nProtect++;
    zeros(REAL(thetaSamples_r), nThetapTilde * nPost);
    double a, b, e;
    // Initiate spatial values
    for (i = 0; i < pTilde; i++) {
      theta[sigmaSqIndx * pTilde + i] = REAL(sigmaSqStarting_r)[i];
      theta[phiIndx * pTilde + i] = REAL(phiStarting_r)[i];
      if (corName == "matern") {
        theta[nuIndx * pTilde + i] = REAL(nuStarting_r)[i];
        nu[i] = theta[nuIndx * pTilde + i];
      }
    } // i
    // Allocate for the U index vector that keep track of which locations have
    // the i-th location as a neighbor
    int nIndx = static_cast<int>(static_cast<double>(1+m)/2*m+(J-m-1)*m);

    // For NNGP
    int mm = m*m;
    double *B = (double *) R_alloc(nIndx * pTilde, sizeof(double));
    double *F = (double *) R_alloc(J * pTilde, sizeof(double));
    double *BCand = (double *) R_alloc(nIndx, sizeof(double));
    double *FCand = (double *) R_alloc(J, sizeof(double));
    double *c =(double *) R_alloc(m*nThreads * pTilde, sizeof(double));
    double *C = (double *) R_alloc(mm*nThreads * pTilde, sizeof(double));
    int sizeBK = nThreads*(1.0+static_cast<int>(floor(nuB[0])));
    double *bk = (double *) R_alloc(nThreads* pTilde * sizeBK, sizeof(double));

    for (i = 0; i < pTilde; i++) {
      updateBFsvcTIntAbund(&B[i * nIndx], &F[i * J], &c[i * m * nThreads],
                           &C[i * mm * nThreads], coords, nnIndx, nnIndxLU, J, m,
                           theta[sigmaSqIndx * pTilde + i], theta[phiIndx * pTilde + i], nu[i],
                           covModel, &bk[i * sizeBK], nuB[i]);
    }

    // Spatial process sums for each observation. 
    double *wObs = (double *) R_alloc(nObs, sizeof(double));
    double *wObsCand = (double *) R_alloc(nObs, sizeof(double));
    double wSite = 0.0;
    // For each location, multiply w x Xw
    for (i = 0; i < nObs; i++) {
      wObs[i] = 0.0;
      for (ll = 0; ll < pTilde; ll++) {
        wObs[i] += w[siteIndx[i] * pTilde + ll] * Xw[ll * JnYears + zLongIndx[i]];
      }
      wObsCand[i] = wObs[i];
    }

    /********************************************************************
      Set up MH stuff
    ********************************************************************/
    double logPostBetaCurr = 0.0, logPostBetaCand = 0.0;
    double logPostAlphaCurr = 0.0, logPostAlphaCand = 0.0;
    double logPostKappaCurr = 0.0, logPostKappaCand = 0.0;
    double logPostThetaCurr = 0.0, logPostThetaCand = 0.0;
    double *logPostWCand = (double *) R_alloc(JpTilde, sizeof(double));
    double *logPostWCurr = (double *) R_alloc(JpTilde, sizeof(double));
    for (j = 0; j < JpTilde; j++) {
      logPostWCurr[j] = R_NegInf;
      logPostWCand[j] = logPostWCurr[j];
    }
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
    double logDet;
    double *phiCand = (double *) R_alloc(pTilde, sizeof(double));
    double *sigmaSqCand = (double *) R_alloc(pTilde, sizeof(double));
    double *nuCand = (double *) R_alloc(pTilde, sizeof(double));
    for (i = 0; i < pTilde; i++) {
      phiCand[i] = theta[phiIndx * pTilde + i];
      sigmaSqCand[i] = theta[sigmaSqIndx * pTilde + i];
      if (corName == "matern") {
        nuCand[i] = theta[nuIndx * pTilde + i];
      }
    }
    double *betaCand = (double *) R_alloc(pAbund, sizeof(double));
    for (j = 0; j < pAbund; j++) {
      betaCand[j] = beta[j];
    }
    double *alphaCand = (double *) R_alloc(pDet, sizeof(double));
    for (j = 0; j < pDet; j++) {
      alphaCand[j] = alpha[j];
    }
    // w is ordered by site, then SVC within site.
    double *wCand = (double *) R_alloc(JpTilde, sizeof(double));
    for (i = 0; i < pTilde; i++) {
      for (j = 0; j < J; j++) {
        wCand[j * pTilde + i] = w[j * pTilde + i];
      }
    }
    double *betaStarCand = (double *) R_alloc(nAbundRE, sizeof(double));
    for (j = 0; j < nAbundRE; j++) {
      betaStarCand[j] = betaStar[j];
    }
    double *alphaStarCand = (double *) R_alloc(nDetRE, sizeof(double));
    for (j = 0; j < nDetRE; j++) {
      alphaStarCand[j] = alphaStar[j];
    }
    double *kappaCand = (double *) R_alloc(nData, sizeof(double));
    for (j = 0; j < nData; j++) {
      kappaCand[j] = kappa[j];
    }
    // theta, beta, alpha, w
    int nAMCMC = nThetapTilde + pAbund + pDet + JpTilde;
    if (pAbundRE > 0) {
      nAMCMC += nAbundRE;
    }
    if (pDetRE > 0) {
      nAMCMC += nDetRE;
    }
    if (anyFamily == 1) {
      nAMCMC += nData;
    }
    int betaAMCMCIndx = 0;
    int alphaAMCMCIndx = betaAMCMCIndx + pAbund;
    int sigmaSqAMCMCIndx = alphaAMCMCIndx + pDet;
    int phiAMCMCIndx = sigmaSqAMCMCIndx + pTilde;
    int nuAMCMCIndx;
    if (corName == "matern") {
      nuAMCMCIndx = phiAMCMCIndx + pTilde;
    } else {
      nuAMCMCIndx = phiAMCMCIndx;
    }
    int wAMCMCIndx = nuAMCMCIndx + pTilde;
    int betaStarAMCMCIndx = wAMCMCIndx + JpTilde;
    int alphaStarAMCMCIndx = betaStarAMCMCIndx + nAbundRE;
    int kappaAMCMCIndx = alphaStarAMCMCIndx + nDetRE;
    double *accept = (double *) R_alloc(nAMCMC, sizeof(double)); zeros(accept, nAMCMC);
    SEXP acceptSamples_r;
    PROTECT(acceptSamples_r = Rf_allocMatrix(REALSXP, nAMCMC, nBatch)); nProtect++;
    zeros(REAL(acceptSamples_r), nAMCMC * nBatch);
    SEXP tuningSamples_r;
    PROTECT(tuningSamples_r = Rf_allocMatrix(REALSXP, nAMCMC, nBatch)); nProtect++;
    zeros(REAL(tuningSamples_r), nAMCMC * nBatch);

    /**********************************************************************
     * Prep for random effects
     * *******************************************************************/
    // Site-level sums of the abundance random effects
    double *betaStarObs = (double *) R_alloc(nObs, sizeof(double));
    double *betaStarObsCand = (double *) R_alloc(nObs, sizeof(double));
    double betaStarSite = 0.0;
    zeros(betaStarObs, nObs);
    // Initial sums
    for (j = 0; j < nObs; j++) {
      for (l = 0; l < pAbundRE; l++) {
        betaStarObs[j] += betaStar[which(XRE[l * JnYears + zLongIndx[j]], betaLevelIndx, nAbundRE)];
      }
      betaStarObsCand[j] = betaStarObs[j];
    }
    // Starting index for abundance random effects
    int *betaStarStart = (int *) R_alloc(pAbundRE, sizeof(int));
    for (l = 0; l < pAbundRE; l++) {
      betaStarStart[l] = which(l, betaStarIndx, nAbundRE);
    }
    // Observation-level sums of the detection random effects
    double *alphaStarObs = (double *) R_alloc(nObs, sizeof(double)); 
    double *alphaStarObsCand = (double *) R_alloc(nObs, sizeof(double)); 
    zeros(alphaStarObs, nObs); 
    // Get sums of the current REs for each site/visit combo
    for (i = 0; i < nObs; i++) {
      for (l = 0; l < alphaNREIndx[i]; l++) {
        alphaStarObs[i] += alphaStar[which(XpRE[l * nObs + i], alphaLevelIndx, nDetRE)];
      }
      alphaStarObsCand[i] = alphaStarObs[i];
    }
    // Starting index for detection random effects
    int *alphaStarStart = (int *) R_alloc(pDetRE, sizeof(int)); 
    for (l = 0; l < pDetRE; l++) {
      alphaStarStart[l] = which(l, alphaStarIndx, nDetRE); 
    }

    logPostBetaCurr = R_NegInf;
    logPostThetaCurr = R_NegInf;
    GetRNGstate();

    for (s = 0, g = 0; s < nBatch; s++) {
      for (t = 0; t < batchLength; t++, g++) {
        /********************************************************************
         *Update Abundance Regression Coefficients
         *******************************************************************/
        for (k = 0; k < pAbund; k++) {
          logPostBetaCand = 0.0;
          logPostBetaCurr = 0.0;
          betaCand[k] = rnorm(beta[k], exp(tuning[betaAMCMCIndx + k]));
          logPostBetaCand += dnorm(betaCand[k], muBeta[k], sqrt(SigmaBeta[k * pAbund + k]), 1);
          logPostBetaCurr += dnorm(beta[k], muBeta[k], sqrt(SigmaBeta[k * pAbund + k]), 1);
          for (j = 0; j < nObs; j++) {
            tmp_nObs[j] = exp(F77_NAME(ddot)(&pAbund, &X[zLongIndx[j]], &JnYears, betaCand, &inc) + 
                              betaStarObs[j] +
                              wObs[j] + 
                              F77_NAME(ddot)(&pDetLong[dataIndx[j]], &Xp[j], &nObs, 
                                             &alpha[stAlpha[j]], &inc) +  
                              alphaStarObs[j]);
            if (family[dataIndx[j]] == 1) {
              logPostBetaCand += dnbinom_mu(y[j], kappa[dataIndx[j]], tmp_nObs[j] * offset[j], 1);
            } else {
              logPostBetaCand += dpois(y[j], tmp_nObs[j] * offset[j], 1);
            }
            tmp_nObs[j] = exp(F77_NAME(ddot)(&pAbund, &X[zLongIndx[j]], &JnYears, beta, &inc) + 
                              betaStarObs[j] +
                              wObs[j] + 
                              F77_NAME(ddot)(&pDetLong[dataIndx[j]], &Xp[j], &nObs, 
                                             &alpha[stAlpha[j]], &inc) + 
                              alphaStarObs[j]);
            if (family[dataIndx[j]] == 1) {
              logPostBetaCurr += dnbinom_mu(y[j], kappa[dataIndx[j]], tmp_nObs[j] * offset[j], 1);
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
         *Update Detection Regression Coefficients
         *******************************************************************/
        for (k = 0; k < pDet; k++) {
          logPostAlphaCand = 0.0;
          logPostAlphaCurr = 0.0;
          alphaCand[k] = rnorm(alpha[k], exp(tuning[alphaAMCMCIndx + k]));
          logPostAlphaCand += dnorm(alphaCand[k], muAlpha[k], sqrt(sigmaAlpha[k]), 1);
          logPostAlphaCurr += dnorm(alpha[k], muAlpha[k], sqrt(sigmaAlpha[k]), 1);
          for (j = 0; j < nObs; j++) {
            if (dataIndx[j] == alphaIndx[k]) { 
              tmp_nObs[j] = exp(F77_NAME(ddot)(&pAbund, &X[zLongIndx[j]], &JnYears, beta, &inc) + 
                                betaStarObs[j] +
                                wObs[j] + 
                                F77_NAME(ddot)(&pDetLong[dataIndx[j]], &Xp[j], &nObs, 
                                               &alphaCand[stAlpha[j]], &inc) + 
                                alphaStarObs[j]);
              if (family[dataIndx[j]] == 1) {
                logPostAlphaCand += dnbinom_mu(y[j], kappa[dataIndx[j]], tmp_nObs[j] * offset[j], 1);
              } else {
                logPostAlphaCand += dpois(y[j], tmp_nObs[j] * offset[j], 1);
              }
              tmp_nObs[j] = exp(F77_NAME(ddot)(&pAbund, &X[zLongIndx[j]], &JnYears, beta, &inc) + 
                                betaStarObs[j] +
                                wObs[j] + 
                                F77_NAME(ddot)(&pDetLong[dataIndx[j]], &Xp[j], &nObs, 
                                               &alpha[stAlpha[j]], &inc) + 
                                alphaStarObs[j]);
              if (family[dataIndx[j]] == 1) {
                logPostAlphaCurr += dnbinom_mu(y[j], kappa[dataIndx[j]], tmp_nObs[j] * offset[j], 1);
              } else {
                logPostAlphaCurr += dpois(y[j], tmp_nObs[j] * offset[j], 1);
              }
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
          tmp_0 = F77_NAME(ddot)(&nAbundRELong[l], &betaStar[betaStarStart[l]],
                                 &inc, &betaStar[betaStarStart[l]], &inc);
          tmp_0 *= 0.5;
          sigmaSqMu[l] = rigamma(sigmaSqMuA[l] + nAbundRELong[l] / 2.0, sigmaSqMuB[l] + tmp_0);
        }

        /********************************************************************
         *Update Detection random effects variance
         *******************************************************************/
        for (l = 0; l < pDetRE; l++) {
          tmp_0 = F77_NAME(ddot)(&nDetRELong[l], &alphaStar[alphaStarStart[l]], 
                                 &inc, &alphaStar[alphaStarStart[l]], &inc); 
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
            for (j = 0; j < nObs; j++) {
              if (XRE[betaStarIndx[l] * JnYears + zLongIndx[j]] == betaLevelIndx[l]) {
                // Candidate
                betaStarObsCand[j] = 0.0;
                for (ll = 0; ll < pAbundRE; ll++) {
                  betaStarObsCand[j] += betaStarCand[which(XRE[ll * JnYears + zLongIndx[j]],
                                                             betaLevelIndx, nAbundRE)];
                }
                tmp_nObs[j] = exp(F77_NAME(ddot)(&pAbund, &X[zLongIndx[j]], &JnYears, beta, &inc) + 
                                  betaStarObsCand[j] +
                                  wObs[j] + 
                                  F77_NAME(ddot)(&pDetLong[dataIndx[j]], &Xp[j], &nObs, 
                                                 &alpha[stAlpha[j]], &inc) + 
                                  alphaStarObs[j]);
                if (family[dataIndx[j]] == 1) {
                  logPostBetaStarCand[l] += dnbinom_mu(y[j], kappa[dataIndx[j]], tmp_nObs[j] * offset[j], 1);
                } else {
                  logPostBetaStarCand[l] += dpois(y[j], tmp_nObs[j] * offset[j], 1);
                }
                // Current
                betaStarObs[j] = 0.0;
                for (ll = 0; ll < pAbundRE; ll++) {
                  betaStarObs[j] += betaStar[which(XRE[ll * JnYears + zLongIndx[j]],
                                                     betaLevelIndx, nAbundRE)];
                }
                tmp_nObs[j] = exp(F77_NAME(ddot)(&pAbund, &X[zLongIndx[j]], &JnYears, beta, &inc) + 
                                  betaStarObs[j] +
                                  wObs[j] + 
                                  F77_NAME(ddot)(&pDetLong[dataIndx[j]], &Xp[j], &nObs, 
                                                 &alpha[stAlpha[j]], &inc) + 
                                  alphaStarObs[j]);
                if (family[dataIndx[j]] == 1) {
                  logPostBetaStarCurr[l] += dnbinom_mu(y[j], kappa[dataIndx[j]], 
                                                       tmp_nObs[j] * offset[j], 1);
                } else {
                  logPostBetaStarCurr[l] += dpois(y[j], tmp_nObs[j] * offset[j], 1);
                }
              }
            }
            if (runif (0.0, 1.0) <= exp(logPostBetaStarCand[l] - logPostBetaStarCurr[l])) {
              betaStar[l] = betaStarCand[l];
              F77_NAME(dcopy)(&nObs, betaStarObsCand, &inc, betaStarObs, &inc);
              accept[betaStarAMCMCIndx + l]++;
            } else {
              betaStarCand[l] = betaStar[l];
              F77_NAME(dcopy)(&nObs, betaStarObs, &inc, betaStarObsCand, &inc);
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
            for (j = 0; j < nObs; j++) {
              if (alphaNREIndx[j] > 0 && XpRE[alphaColIndx[l] * nObs + j] == alphaLevelIndx[l]) {
                // Candidate
                alphaStarObsCand[j] = 0.0;
                for (ll = 0; ll < alphaNREIndx[j]; ll++) {
                  alphaStarObsCand[j] += alphaStarCand[which(XpRE[ll * nObs + j],
                                                             alphaLevelIndx, nDetRE)];
                }
                tmp_nObs[j] = exp(F77_NAME(ddot)(&pAbund, &X[zLongIndx[j]], &JnYears, beta, &inc) + 
                                  betaStarObs[j] +
                                  wObs[j] + 
                                  F77_NAME(ddot)(&pDetLong[dataIndx[j]], &Xp[j], &nObs, 
                                                 &alpha[stAlpha[j]], &inc) + 
                                  alphaStarObsCand[j]);
                if (family[dataIndx[j]] == 1) {
                  logPostAlphaStarCand[l] += dnbinom_mu(y[j], kappa[dataIndx[j]], 
                                                        tmp_nObs[j] * offset[j], 1);
                } else {
                  logPostAlphaStarCand[l] += dpois(y[j], tmp_nObs[j] * offset[j], 1);
                }
                // Current
                alphaStarObs[j] = 0.0;
                for (ll = 0; ll < alphaNREIndx[j]; ll++) {
                  alphaStarObs[j] += alphaStar[which(XpRE[ll * nObs + j],
                                                     alphaLevelIndx, nDetRE)];
                }
                tmp_nObs[j] = exp(F77_NAME(ddot)(&pAbund, &X[zLongIndx[j]], &JnYears, beta, &inc) + 
                                  betaStarObs[j] +
                                  wObs[j] + 
                                  F77_NAME(ddot)(&pDetLong[dataIndx[j]], &Xp[j], &nObs, 
                                                 &alpha[stAlpha[j]], &inc) + 
                                  alphaStarObs[j]);
                if (family[dataIndx[j]] == 1) {
                  logPostAlphaStarCurr[l] += dnbinom_mu(y[j], kappa[dataIndx[j]], 
                                                       tmp_nObs[j] * offset[j], 1);
                } else {
                  logPostAlphaStarCurr[l] += dpois(y[j], tmp_nObs[j] * offset[j], 1);
                }
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
         * Update all spatial parameters one SVC at a time
         *******************************************************************/
        for (ll = 0; ll < pTilde; ll++) {
          /********************************************************************
           *Update w (spatial random effects)
           *******************************************************************/
          for (j = 0; j < J; j++) {
            // Proposal
            a = 0.0;
            // Propose new value
            logPostWCand[j * pTilde + ll] = 0.0;
            wCand[j * pTilde + ll] = rnorm(w[j * pTilde + ll],
                                           exp(tuning[wAMCMCIndx + j * pTilde + ll]));
	          // MVN for any neighbors of j
	          if (uIndxLU[J + j] > 0) { // if current location j is a neighbor for anybody
	            for (r = 0; r < uIndxLU[J + j]; r++) { // how many locations have j as a neighbor
	              jj = uIndx[uIndxLU[j] + r]; // jj is the index of the rth location who has j as a neighbor
	              e = 0;
	              for (i = 0; i < nnIndxLU[J+jj]; i++){ // neighbors of the jjth location
	                e += B[ll * nIndx + nnIndxLU[jj]+i]*wCand[(nnIndx[nnIndxLU[jj]+i]) * pTilde + ll];
                }
                b = wCand[jj * pTilde + ll] - e;
                a += b*b/F[ll * J + jj];
              }
	          }
	          // MVN for j
	          if (nnIndxLU[J + j] > 0) { // if j has any neighbors
              e = 0;
              for(i = 0; i < nnIndxLU[J+j]; i++){
                e += B[ll * nIndx + nnIndxLU[j]+i]*wCand[(nnIndx[nnIndxLU[j]+i]) * pTilde + ll];
              }
              b = wCand[j * pTilde + ll] - e;
            } else{
              b = wCand[j * pTilde + ll];
            }
            a += b*b/F[ll * J + j];
            logPostWCand[j * pTilde + ll] = -0.5*a;
            for (i = 0; i < nObs; i++) {
              if (siteIndx[i] == j) {
                wObsCand[i] = 0.0;
                for (rr = 0; rr < pTilde; rr++) {
                  wObsCand[i] += wCand[j * pTilde + rr] * Xw[rr * JnYears + zLongIndx[i]];
                }
                tmp_nObs[i] = exp(F77_NAME(ddot)(&pAbund, &X[zLongIndx[i]], &JnYears, beta, &inc) + 
                                  betaStarObs[i] +
                                  wObsCand[i] + 
                                  F77_NAME(ddot)(&pDetLong[dataIndx[i]], &Xp[i], &nObs, 
                                                 &alpha[stAlpha[i]], &inc) + 
                                  alphaStarObs[i]);
                if (family[dataIndx[i]] == 1) {
                  logPostWCand[j * pTilde + ll] += dnbinom_mu(y[i], kappa[dataIndx[i]], 
                                                              tmp_nObs[i] * offset[i], 1);
                } else {
                  logPostWCand[j * pTilde + ll] += dpois(y[i], tmp_nObs[i] * offset[i], 1);
                }
              }
            }
            a = 0.0;
            // MVN for any neighbors of j
            if (uIndxLU[J + j] > 0) { // if current location j is a neighbor for anybody
              for (r = 0; r < uIndxLU[J + j]; r++) { // how many locations have j as a neighbor
                jj = uIndx[uIndxLU[j] + r]; // jj is the index of the rth location who has j as a neighbor
                e = 0;
                for (i = 0; i < nnIndxLU[J+jj]; i++){ // neighbors of the jjth location
                  e += B[ll * nIndx + nnIndxLU[jj]+i]*w[(nnIndx[nnIndxLU[jj]+i]) * pTilde + ll];
                }
                b = w[jj * pTilde + ll] - e;
                a += b*b/F[ll * J + jj];
              }
            }
            // MVN for j
            if(nnIndxLU[J+j] > 0){ // if j has any neighbors
              e = 0;
              for(i = 0; i < nnIndxLU[J+j]; i++){
                e += B[ll * nIndx + nnIndxLU[j]+i]*w[(nnIndx[nnIndxLU[j]+i]) * pTilde + ll];
              }
              b = w[j * pTilde + ll] - e;
            } else{
              b = w[j * pTilde + ll];
            }
            a += b*b/F[ll * J + j];
            logPostWCurr[j * pTilde + ll] = -0.5*a;
            for (i = 0; i < nObs; i++) {
              if (siteIndx[i] == j) {
                tmp_nObs[i] = exp(F77_NAME(ddot)(&pAbund, &X[zLongIndx[i]], &JnYears, beta, &inc) + 
                                  betaStarObs[i] +
                                  wObs[i] + 
                                  F77_NAME(ddot)(&pDetLong[dataIndx[i]], &Xp[i], &nObs, 
                                                 &alpha[stAlpha[i]], &inc) + 
                                  alphaStarObs[i]);
                if (family[dataIndx[i]] == 1) {
                  logPostWCurr[j * pTilde + ll] += dnbinom_mu(y[i], kappa[dataIndx[i]], 
                                                              tmp_nObs[i] * offset[i], 1);
                } else {
                  logPostWCurr[j * pTilde + ll] += dpois(y[i], tmp_nObs[i] * offset[i], 1);
                }
              }
            }
            if (runif(0.0, 1.0) <= exp(logPostWCand[j * pTilde + ll] - logPostWCurr[j * pTilde + ll])) {
              w[j * pTilde + ll] = wCand[j * pTilde + ll];
              for (i = 0; i < nObs; i++) {
                wObs[i] = wObsCand[i];
              }
              accept[wAMCMCIndx + j * pTilde + ll]++;
            } else {
              wCand[j * pTilde + ll] = w[j * pTilde + ll];
              for (i = 0; i < nObs; i++) {
                wObsCand[i] = wObs[i];
              }
            }
          } // j

          /********************************************************************
           *Update sigmaSq
           *******************************************************************/
	        a = 0;
#ifdef _OPENMP
#pragma omp parallel for private (e, i, b) reduction(+:a)
#endif
          for (j = 0; j < J; j++){
            if (nnIndxLU[J+j] > 0){
              e = 0;
              for(i = 0; i < nnIndxLU[J+j]; i++){
                e += B[ll * nIndx + nnIndxLU[j]+i]*w[nnIndx[nnIndxLU[j]+i] * pTilde + ll];
              }
              b = w[j * pTilde + ll] - e;
            } else{
              b = w[j * pTilde + ll];
            }
            a += b*b/F[ll * J + j];
          }

	        theta[sigmaSqIndx * pTilde + ll] = rigamma(sigmaSqA[ll] + J / 2.0,
                                                     sigmaSqB[ll] + 0.5 * a * theta[sigmaSqIndx * pTilde + ll]);

          /********************************************************************
           *Update phi (and nu if matern)
           *******************************************************************/
          if (corName == "matern"){
            nu[ll] = theta[nuIndx * pTilde + ll];
          }
          updateBFsvcTIntAbund(&B[ll * nIndx], &F[ll * J], &c[ll * m * nThreads],
                           &C[ll * mm * nThreads], coords, nnIndx, nnIndxLU, J, m,
                           theta[sigmaSqIndx * pTilde + ll], theta[phiIndx * pTilde + ll], nu[ll],
                           covModel, &bk[ll * sizeBK], nuB[ll]);
          a = 0;
          logDet = 0;

#ifdef _OPENMP
#pragma omp parallel for private (e, i, b) reduction(+:a, logDet)
#endif
          for (j = 0; j < J; j++){
            if (nnIndxLU[J+j] > 0){
              e = 0;
              for (i = 0; i < nnIndxLU[J+j]; i++){
                e += B[ll * nIndx + nnIndxLU[j]+i]*w[(nnIndx[nnIndxLU[j]+i]) * pTilde + ll];
              }
              b = w[j * pTilde + ll] - e;
            } else{
              b = w[j * pTilde + ll];
            }
            a += b*b/F[ll * J + j];
            logDet += log(F[ll * J + j]);
          }

          logPostThetaCurr = -0.5*logDet - 0.5*a;
          logPostThetaCurr += log(theta[phiIndx * pTilde + ll] - phiA[ll]) +
                              log(phiB[ll] - theta[phiIndx * pTilde + ll]);
          if (corName == "matern"){
           logPostThetaCurr += log(theta[nuIndx * pTilde + ll] - nuA[ll]) +
                                log(nuB[ll] - theta[nuIndx * pTilde + ll]);
          }

          // Candidate
          phiCand[ll] = logitInv(rnorm(logit(theta[phiIndx * pTilde + ll], phiA[ll], phiB[ll]),
                                       exp(tuning[phiAMCMCIndx + ll])), phiA[ll], phiB[ll]);
          if (corName == "matern"){
            nuCand[ll] = logitInv(rnorm(logit(theta[nuIndx * pTilde + ll], nuA[ll], nuB[ll]),
                                        exp(tuning[nuAMCMCIndx + ll])), nuA[ll], nuB[ll]);
          }

          updateBFsvcTIntAbund(BCand, FCand, &c[ll * m * nThreads], &C[ll * mm * nThreads], 
                           coords, nnIndx, nnIndxLU, J, m,
                           theta[sigmaSqIndx * pTilde + ll], phiCand[ll], nuCand[ll], covModel, 
                           &bk[ll * sizeBK], nuB[ll]);

          a = 0;
          logDet = 0;

#ifdef _OPENMP
#pragma omp parallel for private (e, i, b) reduction(+:a, logDet)
#endif
          for (j = 0; j < J; j++){
            if (nnIndxLU[J+j] > 0){
              e = 0;
              for (i = 0; i < nnIndxLU[J+j]; i++){
                e += BCand[nnIndxLU[j]+i]*w[(nnIndx[nnIndxLU[j]+i]) * pTilde + ll];
              }
              b = w[j * pTilde + ll] - e;
            } else{
              b = w[j * pTilde + ll];
            }
            a += b*b/FCand[j];
            logDet += log(FCand[j]);
          }

          logPostThetaCand = -0.5*logDet - 0.5*a;
          logPostThetaCand += log(phiCand[ll] - phiA[ll]) + log(phiB[ll] - phiCand[ll]);
          if (corName == "matern"){
            logPostThetaCand += log(nuCand[ll] - nuA[ll]) + log(nuB[ll] - nuCand[ll]);
          }
          
          if (runif(0.0,1.0) <= exp(logPostThetaCand - logPostThetaCurr)) {

            F77_NAME(dcopy)(&nIndx, BCand, &inc, &B[ll * nIndx], &inc);
            F77_NAME(dcopy)(&J, FCand, &inc, &F[ll * J], &inc);

            theta[phiIndx * pTilde + ll] = phiCand[ll];
            accept[phiAMCMCIndx + ll]++;
            if(corName == "matern"){
              theta[nuIndx * pTilde + ll] = nuCand[ll];
              accept[nuAMCMCIndx + ll]++;
            }
          }
        } // svc

        /********************************************************************
         *Update kappa (the NB size parameter)
         *******************************************************************/
        for (k = 0; k < nData; k++) {
          if (family[k] == 1) {
            kappaCand[k] = logitInv(rnorm(logit(kappa[k], kappaA[k], kappaB[k]), 
                                    exp(tuning[kappaAMCMCIndx + k])),
                                    kappaA[k], kappaB[k]);
            logPostKappaCurr = 0.0;
            logPostKappaCand = 0.0;
            for (j = 0; j < nObs; j++) {
              if (dataIndx[j] == k) {
                tmp_nObs[j] = exp(F77_NAME(ddot)(&pAbund, &X[zLongIndx[j]], &JnYears, beta, &inc) + 
                                  betaStarObs[j] +
                                  wObs[j] + 
                                  F77_NAME(ddot)(&pDetLong[dataIndx[j]], &Xp[j], &nObs, 
                                                 &alpha[stAlpha[j]], &inc) + 
                                  alphaStarObs[j]);
                logPostKappaCurr += dnbinom_mu(y[j], kappa[k], tmp_nObs[j] * offset[j], 1);
                logPostKappaCand += dnbinom_mu(y[j], kappaCand[k], tmp_nObs[j] * offset[j], 1);
              }
            }
            // Jacobian adjustment
            logPostKappaCurr += log(kappa[k] - kappaA[k]) + log(kappaB[k] - kappa[k]);
            logPostKappaCand += log(kappaCand[k] - kappaA[k]) + log(kappaB[k] - kappaCand[k]);
            if (runif(0.0, 1.0) <= exp(logPostKappaCand - logPostKappaCurr)) {
              kappa[k] = kappaCand[k];
              accept[kappaAMCMCIndx + k]++;
            }
          }
        }

        /********************************************************************
         *Get fitted values
         *******************************************************************/
        if (saveFitted == 1) {
          // Process-level means
          for (j = 0; j < JnYears; j++) {
            wSite = 0.0;
            betaStarSite = 0.0;
            for (ll = 0; ll < pTilde; ll++) {
              wSite += w[zSiteIndx[j] * pTilde + ll] * Xw[ll * JnYears + j];
            }
            for (ll = 0; ll < pAbundRE; ll++) {
              betaStarSite += betaStar[which(XRE[ll * JnYears + j], betaLevelIndx, nAbundRE)];
            }
            mu[j] = exp(F77_NAME(ddot)(&pAbund, &X[j], &JnYears, beta, &inc) + 
                        betaStarSite + 
                        wSite);
          }
          // Observation-level means and replicate data
          for (j = 0; j < nObs; j++) {
            lambda[j] = exp(F77_NAME(ddot)(&pAbund, &X[zLongIndx[j]], &JnYears, beta, &inc) + 
                            betaStarObs[j] +
                            wObs[j] + 
                            F77_NAME(ddot)(&pDetLong[dataIndx[j]], &Xp[j], &nObs, 
                                           &alpha[stAlpha[j]], &inc) + 
                            alphaStarObs[j]);
            if (family[dataIndx[j]] == 0) {
              yRep[j] = rpois(lambda[j] * offset[j]);
              like[j] = dpois(y[j], lambda[j] * offset[j], 0);
            } else {
              yRep[j] = rnbinom_mu(kappa[dataIndx[j]], lambda[j] * offset[j]);
              like[j] = dnbinom_mu(y[j], kappa[dataIndx[j]], lambda[j] * offset[j], 0);
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
            if (saveFitted == 1) {
              F77_NAME(dcopy)(&nObs, lambda, &inc, &REAL(lambdaSamples_r)[sPost*nObs], &inc);
              F77_NAME(dcopy)(&nObs, yRep, &inc, &REAL(yRepSamples_r)[sPost*nObs], &inc);
              F77_NAME(dcopy)(&nObs, like, &inc,
                              &REAL(likeSamples_r)[sPost*nObs], &inc);
            }
            F77_NAME(dcopy)(&JnYears, mu, &inc, &REAL(muSamples_r)[sPost*JnYears], &inc);
            F77_NAME(dcopy)(&nThetapTilde, theta, &inc, &REAL(thetaSamples_r)[sPost*nThetapTilde], &inc);
            F77_NAME(dcopy)(&JpTilde, w, &inc, &REAL(wSamples_r)[sPost*JpTilde], &inc);
            if (anyFamily == 1) {
              F77_NAME(dcopy)(&nData, kappa, &inc, &REAL(kappaSamples_r)[sPost*nData], &inc);
            }
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
          for (ll = 0; ll < pTilde; ll++) {
	          Rprintf("\tphi[%i]\t\t%3.1f\t\t%1.5f\n", ll + 1, 100.0*REAL(acceptSamples_r)[s * nAMCMC + phiAMCMCIndx + ll], exp(tuning[phiAMCMCIndx + ll]));
	          if (corName == "matern") {
	            Rprintf("\tnu[%i]\t\t%3.1f\t\t%1.5f\n", ll + 1, 100.0*REAL(acceptSamples_r)[s * nAMCMC + nuAMCMCIndx + ll], exp(tuning[nuAMCMCIndx + ll]));
	          }
          }
          for (ll = 0; ll < nData; ll++) {
            if (family[ll] == 1) {
              Rprintf("\tkappa[%i]\t%3.1f\t\t%1.5f\n", ll + 1, 100.0*REAL(acceptSamples_r)[s * nAMCMC + kappaAMCMCIndx + ll], exp(tuning[kappaAMCMCIndx + ll]));
            }
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
    if (pAbundRE > 0) {
      nResultListObjs += 2;
    }
    if (pDetRE > 0) {
      nResultListObjs += 2;
    }
    if (anyFamily == 1) {
      nResultListObjs += nData;
    }

    PROTECT(result_r = Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;
    PROTECT(resultName_r = Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;

    // Setting the components of the output list.
    SET_VECTOR_ELT(result_r, 0, betaSamples_r);
    if (saveFitted == 1) {
      SET_VECTOR_ELT(result_r, 1, yRepSamples_r);
      SET_VECTOR_ELT(result_r, 2, lambdaSamples_r);
      SET_VECTOR_ELT(result_r, 3, likeSamples_r);
    }
    SET_VECTOR_ELT(result_r, 4, muSamples_r);
    SET_VECTOR_ELT(result_r, 5, alphaSamples_r);
    SET_VECTOR_ELT(result_r, 6, wSamples_r);
    SET_VECTOR_ELT(result_r, 7, thetaSamples_r);
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
    if (anyFamily == 1) {
      if ((pDetRE > 0) || (pAbundRE > 0)) {
        tmp_02 = tmp_0 + 2;
      } else {
        tmp_02 = 8;
      }
      SET_VECTOR_ELT(result_r, tmp_02, kappaSamples_r);
    }

    SET_VECTOR_ELT(resultName_r, 0, Rf_mkChar("beta.samples"));
    if (saveFitted == 1) {
      SET_VECTOR_ELT(resultName_r, 1, Rf_mkChar("y.rep.samples"));
      SET_VECTOR_ELT(resultName_r, 2, Rf_mkChar("lambda.samples"));
      SET_VECTOR_ELT(resultName_r, 3, Rf_mkChar("like.samples"));
    }
    SET_VECTOR_ELT(resultName_r, 4, Rf_mkChar("mu.samples"));
    SET_VECTOR_ELT(resultName_r, 5, Rf_mkChar("alpha.samples"));
    SET_VECTOR_ELT(resultName_r, 6, Rf_mkChar("w.samples"));
    SET_VECTOR_ELT(resultName_r, 7, Rf_mkChar("theta.samples"));
    if (pDetRE > 0) {
      SET_VECTOR_ELT(resultName_r, 8, Rf_mkChar("sigma.sq.p.samples"));
      SET_VECTOR_ELT(resultName_r, 9, Rf_mkChar("alpha.star.samples"));
    }
    if (pAbundRE > 0) {
      SET_VECTOR_ELT(resultName_r, tmp_0, Rf_mkChar("sigma.sq.mu.samples"));
      SET_VECTOR_ELT(resultName_r, tmp_0 + 1, Rf_mkChar("beta.star.samples"));
    }
    if (anyFamily == 1) {
      SET_VECTOR_ELT(resultName_r, tmp_02, Rf_mkChar("kappa.samples"));
    }

    Rf_namesgets(result_r, resultName_r);

    UNPROTECT(nProtect);

    return(result_r);
  }
}
