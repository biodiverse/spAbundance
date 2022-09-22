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

void updateBFSFAbund(double *B, double *F, double *c, double *C, double *coords, int *nnIndx, int *nnIndxLU, int n, int m, double sigmaSq, double phi, double nu, int covModel, double *bk, double nuUnifb){

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
  SEXP sfMsAbundNNGP(SEXP y_r, SEXP X_r, SEXP coords_r, SEXP XRE_r,  
	            SEXP consts_r, SEXP nAbundRELong_r, 
                    SEXP m_r, SEXP nnIndx_r, SEXP nnIndxLU_r, SEXP uIndx_r, SEXP uIndxLU_r, SEXP uiIndx_r,
	            SEXP betaStarting_r, SEXP kappaStarting_r,  
	            SEXP betaCommStarting_r,  
	            SEXP tauSqBetaStarting_r, 
		    SEXP phiStarting_r, SEXP lambdaStarting_r, SEXP nuStarting_r, SEXP wStarting_r,
	            SEXP sigmaSqMuStarting_r,  
	            SEXP betaStarStarting_r,  
	            SEXP siteIndx_r, SEXP betaStarIndx_r, SEXP betaLevelIndx_r, 
	            SEXP muBetaComm_r,  
	            SEXP SigmaBetaComm_r, SEXP kappaA_r, 
	            SEXP kappaB_r, SEXP tauSqBetaA_r, 
	            SEXP tauSqBetaB_r, 
		    SEXP phiA_r, SEXP phiB_r, SEXP nuA_r, SEXP nuB_r,
	            SEXP sigmaSqMuA_r, SEXP sigmaSqMuB_r, 
	            SEXP tuning_r, SEXP covModel_r, 
		    SEXP nBatch_r, SEXP batchLength_r, SEXP acceptRate_r,
	            SEXP nThreads_r, SEXP verbose_r, SEXP nReport_r, 
	            SEXP samplesInfo_r, SEXP chainInfo_r){
   
    /**********************************************************************
     * Initial constants
     * *******************************************************************/
    int i, j, s, t, g, l, h, r, ll, ii, jj, kk, k, info, nProtect=0;
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
    double *coords = REAL(coords_r); 
    int *XRE = INTEGER(XRE_r); 
    int m = INTEGER(m_r)[0]; 
    // Load constants
    int nSp = INTEGER(consts_r)[0]; 
    int J = INTEGER(consts_r)[1];
    int nObs = INTEGER(consts_r)[2]; 
    int pAbund = INTEGER(consts_r)[3];
    int pAbundRE = INTEGER(consts_r)[4];
    int nAbundRE = INTEGER(consts_r)[5];
    int q = INTEGER(consts_r)[6]; 
    int ppAbund = pAbund * pAbund; 
    double *muBetaComm = REAL(muBetaComm_r); 
    double *SigmaBetaCommInv = (double *) R_alloc(ppAbund, sizeof(double));   
    F77_NAME(dcopy)(&ppAbund, REAL(SigmaBetaComm_r), &inc, SigmaBetaCommInv, &inc);
    double *tauSqBetaA = REAL(tauSqBetaA_r); 
    double *tauSqBetaB = REAL(tauSqBetaB_r); 
    double *sigmaSqMuA = REAL(sigmaSqMuA_r); 
    double *sigmaSqMuB = REAL(sigmaSqMuB_r); 
    double *kappaA = REAL(kappaA_r); 
    double *kappaB = REAL(kappaB_r); 
    double *phiA = REAL(phiA_r); 
    double *phiB = REAL(phiB_r); 
    double *nuA = REAL(nuA_r); 
    double *nuB = REAL(nuB_r); 
    int *nnIndx = INTEGER(nnIndx_r);
    int *nnIndxLU = INTEGER(nnIndxLU_r);
    int *uIndx = INTEGER(uIndx_r);
    int *uIndxLU = INTEGER(uIndxLU_r);
    int *uiIndx = INTEGER(uiIndx_r);
    int covModel = INTEGER(covModel_r)[0];
    std::string corName = getCorName(covModel);
    int *nAbundRELong = INTEGER(nAbundRELong_r); 
    int *siteIndx = INTEGER(siteIndx_r); 
    int *betaStarIndx = INTEGER(betaStarIndx_r);
    int *betaLevelIndx = INTEGER(betaLevelIndx_r);
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
    int currDim = 0;

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
        Rprintf("Spatial Factor NNGP Multispecies Negative Binomial Abundance model with Polya-Gamma latent\nvariable fit with %i sites and %i species.\n\n", J, nSp);
        Rprintf("Samples per Chain: %i \n", nSamples);
        Rprintf("Burn-in: %i \n", nBurn); 
        Rprintf("Thinning Rate: %i \n", nThin); 
	Rprintf("Number of Chains: %i \n", nChain);
        Rprintf("Total Posterior Samples: %i \n\n", nPost * nChain); 
        Rprintf("Using the %s spatial correlation model.\n\n", corName.c_str());
        Rprintf("Using %i latent spatial factors.\n", q);
        Rprintf("Using %i nearest neighbors.\n\n", m);
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
    int nObsnSp = nObs * nSp; 
    int nAbundREnSp = nAbundRE * nSp; 
    int JnSp = J * nSp;
    int JpAbund = J * pAbund; 
    int nObspAbund = nObs * pAbund;
    int Jq = J * q;
    int nObsq = nObs * q;
    int qq = q * q;
    int nSpq = nSp * q;
    double tmp_0; 
    double *tmp_one = (double *) R_alloc(inc, sizeof(double)); 
    double *tmp_ppAbund = (double *) R_alloc(ppAbund, sizeof(double)); 
    double *tmp_pAbund = (double *) R_alloc(pAbund, sizeof(double));
    double *tmp_beta = (double *) R_alloc(pAbund, sizeof(double));
    double *tmp_pAbund2 = (double *) R_alloc(pAbund, sizeof(double));
    int *tmp_JInt = (int *) R_alloc(J, sizeof(int));
    for (j = 0; j < J; j++) {
      tmp_JInt[j] = 0; 
    }
    double *tmp_nObs = (double *) R_alloc(nObs, sizeof(double)); 
    double *tmp_JpAbund = (double *) R_alloc(JpAbund, sizeof(double));
    double *tmp_nObspAbund = (double *) R_alloc(nObspAbund, sizeof(double));
    double *tmp_J1 = (double *) R_alloc(J, sizeof(double));
    double *tmp_qq = (double *) R_alloc(qq, sizeof(double));
    double *tmp_q = (double *) R_alloc(q, sizeof(double));
    double *tmp_q2 = (double *) R_alloc(q, sizeof(double));
    double *tmp_qq2 = (double *) R_alloc(qq, sizeof(double));
    double *tmp_Jq = (double *) R_alloc(Jq, sizeof(double));
    double *tmp_nObsq = (double *) R_alloc(nObsq, sizeof(double));
    double *tmp_nSpq = (double *) R_alloc(nSpq, sizeof(double));
    double *tmp_nSp = (double *) R_alloc(nSp, sizeof(double));

    /**********************************************************************
     * Parameters
     * *******************************************************************/
    double *betaComm = (double *) R_alloc(pAbund, sizeof(double)); 
    F77_NAME(dcopy)(&pAbund, REAL(betaCommStarting_r), &inc, betaComm, &inc);
    double *tauSqBeta = (double *) R_alloc(pAbund, sizeof(double)); 
    F77_NAME(dcopy)(&pAbund, REAL(tauSqBetaStarting_r), &inc, tauSqBeta, &inc);
    double *beta = (double *) R_alloc(pAbundnSp, sizeof(double));   
    F77_NAME(dcopy)(&pAbundnSp, REAL(betaStarting_r), &inc, beta, &inc);
    // Spatial random effects
    double *w = (double *) R_alloc(Jq, sizeof(double));
    F77_NAME(dcopy)(&Jq, REAL(wStarting_r), &inc, w, &inc);
    // Latent spatial factors
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
    // Latent abundance random effects
    double *betaStar = (double *) R_alloc(nAbundREnSp, sizeof(double)); 
    F77_NAME(dcopy)(&nAbundREnSp, REAL(betaStarStarting_r), &inc, betaStar, &inc); 
    // Overdispersion parameter
    double *kappa = (double *) R_alloc(nSp, sizeof(double)); 
    F77_NAME(dcopy)(&nSp, REAL(kappaStarting_r), &inc, kappa, &inc); 
    // Latent Abundance
    double *yRep = (double *) R_alloc(nObsnSp, sizeof(double)); zeros(yRep, nObsnSp); 
    // Auxiliary variables
    double *omegaAbund = (double *) R_alloc(nObsnSp, sizeof(double)); zeros(omegaAbund, nObsnSp);
    double *yStar = (double *) R_alloc(nObsnSp, sizeof(double)); zeros(yStar, nObsnSp);

    /**********************************************************************
     * Return Stuff
     * *******************************************************************/
    // Community level
    SEXP betaCommSamples_r; 
    PROTECT(betaCommSamples_r = allocMatrix(REALSXP, pAbund, nPost)); nProtect++;
    SEXP tauSqBetaSamples_r; 
    PROTECT(tauSqBetaSamples_r = allocMatrix(REALSXP, pAbund, nPost)); nProtect++; 
    // Species level
    SEXP betaSamples_r;
    PROTECT(betaSamples_r = allocMatrix(REALSXP, pAbundnSp, nPost)); nProtect++;
    SEXP yRepSamples_r; 
    PROTECT(yRepSamples_r = allocMatrix(REALSXP, nObsnSp, nPost)); nProtect++; 
    SEXP muSamples_r; 
    PROTECT(muSamples_r = allocMatrix(REALSXP, nObsnSp, nPost)); nProtect++; 
    // Spatial parameters
    SEXP lambdaSamples_r; 
    PROTECT(lambdaSamples_r = allocMatrix(REALSXP, nSpq, nPost)); nProtect++;
    SEXP wSamples_r; 
    PROTECT(wSamples_r = allocMatrix(REALSXP, Jq, nPost)); nProtect++; 
    // Abundance random effects
    SEXP sigmaSqMuSamples_r; 
    SEXP betaStarSamples_r; 
    if (pAbundRE > 0) {
      PROTECT(sigmaSqMuSamples_r = allocMatrix(REALSXP, pAbundRE, nPost)); nProtect++;
      PROTECT(betaStarSamples_r = allocMatrix(REALSXP, nAbundREnSp, nPost)); nProtect++;
    }
    // Overdispersion
    SEXP kappaSamples_r;
    PROTECT(kappaSamples_r = allocMatrix(REALSXP, nSp, nPost)); nProtect++;
    // Likelihood samples for WAIC. 
    SEXP likeSamples_r;
    PROTECT(likeSamples_r = allocMatrix(REALSXP, nObsnSp, nPost)); nProtect++;
    
    /**********************************************************************
     * Additional Sampler Prep
     * *******************************************************************/
    double *like = (double *) R_alloc(nObsnSp, sizeof(double)); 
    zeros(like, nObsnSp);
    double *psi = (double *) R_alloc(nObs, sizeof(double)); 
    zeros(psi, nObs); 
    double *mu = (double *) R_alloc(nObsnSp, sizeof(double)); 
    zeros(mu, nObsnSp); 

    // For normal priors
    F77_NAME(dpotrf)(lower, &pAbund, SigmaBetaCommInv, &pAbund, &info FCONE); 
    if(info != 0){error("c++ error: dpotrf SigmaBetaCommInv failed\n");}
    F77_NAME(dpotri)(lower, &pAbund, SigmaBetaCommInv, &pAbund, &info FCONE); 
    if(info != 0){error("c++ error: dpotri SigmaBetaCommInv failed\n");}
    double *SigmaBetaCommInvMuBeta = (double *) R_alloc(pAbund, sizeof(double)); 
    F77_NAME(dsymv)(lower, &pAbund, &one, SigmaBetaCommInv, &pAbund, muBetaComm, &inc, &zero, 
        	    SigmaBetaCommInvMuBeta, &inc FCONE);
    // Put community level variances in a pAbund x PAbund matrix.
    double *TauBetaInv = (double *) R_alloc(ppAbund, sizeof(double)); zeros(TauBetaInv, ppAbund); 
    for (i = 0; i < pAbund; i++) {
      TauBetaInv[i * pAbund + i] = tauSqBeta[i]; 
    } // i
    F77_NAME(dpotrf)(lower, &pAbund, TauBetaInv, &pAbund, &info FCONE); 
    if(info != 0){error("c++ error: dpotrf TauBetaInv failed\n");}
    F77_NAME(dpotri)(lower, &pAbund, TauBetaInv, &pAbund, &info FCONE); 
    if(info != 0){error("c++ error: dpotri TauBetaInv failed\n");}

    // For NB PG sampling (recommendations from Polson et al. 2013). 
    int trunc = 200;
    double *tmp_trunc = (double *) R_alloc(trunc, sizeof(double));

    /**********************************************************************
     * Prep for random effects
     * *******************************************************************/
    // Site-level sums of the occurrence random effects
    double *betaStarSites = (double *) R_alloc(nObsnSp, sizeof(double)); 
    zeros(betaStarSites, nObsnSp); 
    // Initial sums (initiate with the first species)
    for (i = 0; i < nSp; i++) {
      for (j = 0; j < nObs; j++) {
        for (l = 0; l < pAbundRE; l++) {
          betaStarSites[i * nObs + j] += betaStar[i * nAbundRE + which(XRE[l * nObs + j], betaLevelIndx, nAbundRE)];
        }
      }
    }
    // Starting index for occurrence random effects
    int *betaStarStart = (int *) R_alloc(pAbundRE, sizeof(int)); 
    for (l = 0; l < pAbundRE; l++) {
      betaStarStart[l] = which(l, betaStarIndx, nAbundRE); 
    }

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
    // Species-level spatial random effects
    double *wStar = (double *) R_alloc(JnSp, sizeof(double)); zeros(wStar, JnSp);
    // Multiply Lambda %*% w[j] to get wStar. 
    for (j = 0; j < J; j++) {
      F77_NAME(dgemv)(ntran, &nSp, &q, &one, lambda, &nSp, &w[j*q], &inc, &zero, &wStar[j * nSp], &inc FCONE);
    }
    // For NNGP
    double b, e, aij, aa; 
    double *a = (double *) R_alloc(q, sizeof(double));
    double *v = (double *) R_alloc(q, sizeof(double));
    double *muNNGP = (double *) R_alloc(q, sizeof(double));
    double *var = (double *) R_alloc(qq, sizeof(double));
    double *ff = (double *) R_alloc(q, sizeof(double));
    double *gg = (double *) R_alloc(q, sizeof(double));

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
      updateBFSFAbund(&B[ll * nIndx], &F[ll*J], &c[ll * m*nThreads], &C[ll * mm * nThreads], coords, nnIndx, nnIndxLU, J, m, theta[sigmaSqIndx * q + ll], theta[phiIndx * q + ll], nu[ll], covModel, &bk[ll * sizeBK], nuB[0]);
    }

    /********************************************************************
      Set up MH stuff
    ********************************************************************/
    double logPostCurr = 0.0, logPostCand = 0.0, logDet;
    double kappaCand = 0.0, phiCand = 0.0, nuCand = 0.0; 
    // Also tuning kappa (the overdispersion parameter), which has n.sp values. 
    // Note I also include nu regardless, which differs from previous stuff.  
    int nAMCMC = 3 * q + nSp;
    int sigmaSqAMCMCIndx = 0, phiAMCMCIndx = 1, nuAMCMCIndx = 2, kappaAMCMCIndx = 3;
    double *accept = (double *) R_alloc(nAMCMC, sizeof(double)); zeros(accept, nAMCMC); 
    SEXP acceptSamples_r; 
    PROTECT(acceptSamples_r = allocMatrix(REALSXP, nAMCMC, nBatch)); nProtect++; 
    SEXP tuningSamples_r; 
    PROTECT(tuningSamples_r = allocMatrix(REALSXP, nAMCMC, nBatch)); nProtect++; 
    
    GetRNGstate();

    /**********************************************************************
     Start sampling
     ********************************************************************/
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
        for (h = 0; h < pAbund; h++) {
          tmp_pAbund[h] += SigmaBetaCommInvMuBeta[h];  
        } // j

        /********************************
         Compute A.beta.comm
         *******************************/
        for (h = 0; h < ppAbund; h++) {
          tmp_ppAbund[h] = SigmaBetaCommInv[h] + nSp * TauBetaInv[h]; 
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
         Update Community Abundance Variance Parameter
        ********************************************************************/
        for (h = 0; h < pAbund; h++) {
          tmp_0 = 0.0;  
          for (i = 0; i < nSp; i++) {
            tmp_0 += (beta[h * nSp + i] - betaComm[h]) * (beta[h * nSp + i] - betaComm[h]);
          } // i
          tmp_0 *= 0.5;
          tauSqBeta[h] = rigamma(tauSqBetaA[h] + nSp / 2.0, tauSqBetaB[h] + tmp_0); 
        } // h
        for (h = 0; h < pAbund; h++) {
          TauBetaInv[h * pAbund + h] = tauSqBeta[h]; 
        } // h
        F77_NAME(dpotrf)(lower, &pAbund, TauBetaInv, &pAbund, &info FCONE); 
        if(info != 0){error("c++ error: dpotrf TauBetaInv failed\n");}
        F77_NAME(dpotri)(lower, &pAbund, TauBetaInv, &pAbund, &info FCONE); 
        if(info != 0){error("c++ error: dpotri TauBetaInv failed\n");}

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

        /**********************************************************************
         *Species specific variables  
         *********************************************************************/
        for (i = 0; i < nSp; i++) {  
          /********************************************************************
           *Update Abundance Auxiliary Variables 
           *******************************************************************/
          for (j = 0; j < nObs; j++) {
            omegaAbund[j * nSp + i] = rpgGamma(kappa[i] + y[j * nSp + i], F77_NAME(ddot)(&pAbund, &X[j], &nObs, &beta[i], &nSp) + betaStarSites[i * nObs + j] + wStar[siteIndx[j] * nSp + i], trunc, tmp_trunc);
          } // j
          /********************************************************************
           *Update Abundance Regression Coefficients
           *******************************************************************/
          for (j = 0; j < nObs; j++) {
            yStar[j * nSp + i] = (y[j * nSp + i] - kappa[i]) / (2.0 * omegaAbund[j * nSp + i]);
            tmp_nObs[j] = (yStar[j * nSp + i] - betaStarSites[i * nObs + j] - 
	        	 wStar[siteIndx[j] * nSp + i]) * omegaAbund[j * nSp + i]; 
          } // j
          /********************************
           * Compute b.beta
           *******************************/
          F77_NAME(dgemv)(ytran, &nObs, &pAbund, &one, X, &nObs, tmp_nObs, &inc, &zero, tmp_pAbund, &inc FCONE);           // TauBetaInv %*% betaComm + tmp_pAbund = tmp_pAbund
          F77_NAME(dgemv)(ntran, &pAbund, &pAbund, &one, TauBetaInv, &pAbund, betaComm, &inc, &one, tmp_pAbund, &inc FCONE); 

          /********************************
           * Compute A.beta
           * *****************************/
          for(j = 0; j < nObs; j++){
            for(h = 0; h < pAbund; h++){
              tmp_nObspAbund[h*nObs+j] = X[h*nObs+j] * omegaAbund[j * nSp + i];
            }
          }
          F77_NAME(dgemm)(ytran, ntran, &pAbund, &pAbund, &nObs, &one, X, &nObs, tmp_nObspAbund, &nObs, &zero, tmp_ppAbund, &pAbund FCONE FCONE);
          for (h = 0; h < ppAbund; h++) {
            tmp_ppAbund[h] += TauBetaInv[h]; 
          } // h
          F77_NAME(dpotrf)(lower, &pAbund, tmp_ppAbund, &pAbund, &info FCONE); 
          if(info != 0){error("c++ error: dpotrf ABeta failed\n");}
          F77_NAME(dpotri)(lower, &pAbund, tmp_ppAbund, &pAbund, &info FCONE); 
          if(info != 0){error("c++ error: dpotri ABeta failed\n");}
          F77_NAME(dsymv)(lower, &pAbund, &one, tmp_ppAbund, &pAbund, tmp_pAbund, &inc, &zero, tmp_pAbund2, &inc FCONE);
          F77_NAME(dpotrf)(lower, &pAbund, tmp_ppAbund, &pAbund, &info FCONE); 
          if(info != 0){error("c++ error: dpotrf here failed\n");}
          mvrnorm(tmp_beta, tmp_pAbund2, tmp_ppAbund, pAbund);
          for (h = 0; h < pAbund; h++) {
            beta[h * nSp + i] = tmp_beta[h]; 
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
              for (j = 0; j < nObs; j++) {
                if (XRE[betaStarIndx[l] * nObs + j] == betaLevelIndx[l]) {
                  tmp_one[0] += (yStar[j * nSp + i] - F77_NAME(ddot)(&pAbund, &X[j], &nObs, &beta[i], &nSp) - wStar[siteIndx[j] * nSp + i] - betaStarSites[i * nObs + j] + betaStar[i * nAbundRE + l]) * omegaAbund[j * nSp + i];
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
            zeros(&betaStarSites[i * nObs], nObs);
            for (j = 0; j < nObs; j++) {
              for (l = 0; l < pAbundRE; l++) {
                betaStarSites[i * nObs + j] += betaStar[i * nAbundRE + which(XRE[l * nObs + j], betaLevelIndx, nAbundRE)];
              }
            }
          }
	} // i (species)
        /********************************************************************
         *Update Spatial Random Effects (w)
         *******************************************************************/
	// Update B and F
        for (ll = 0; ll < q; ll++) {
          updateBFSFAbund(&B[ll * nIndx], &F[ll*J], &c[ll * m*nThreads], &C[ll * mm * nThreads], coords, nnIndx, nnIndxLU, J, m, theta[sigmaSqIndx * q + ll], theta[phiIndx * q + ll], nu[ll], covModel, &bk[ll * sizeBK], nuB[0]);
        }

	for (ii = 0; ii < J; ii++) {
          zeros(tmp_qq, qq);
	  for (j = 0; j < nObs; j++) {
            if(siteIndx[j] == ii) {
              zeros(tmp_nSpq, nSpq);
	      for (i = 0; i < nSp; i++) {
                for (ll = 0; ll < q; ll++) {
                  tmp_nSpq[ll * nSp + i] += lambda[ll * nSp + i] * omegaAbund[j * nSp + i];
                } // ll
              } // i
	      F77_NAME(dgemm)(ytran, ntran, &q, &q, &nSp, &one, tmp_nSpq, &nSp, lambda, 
			      &nSp, &one, tmp_qq, &q FCONE FCONE);
	    }
	  } // j

	  for (ll = 0; ll < q; ll++) {

            a[ll] = 0; 
	    v[ll] = 0; 

	    if (uIndxLU[J + ii] > 0){ // is ii a neighbor for anybody
	      for (j = 0; j < uIndxLU[J+ii]; j++){ // how many locations have ii as a neighbor
	        b = 0;
	        // now the neighbors for the jth location who has ii as a neighbor
	        jj = uIndx[uIndxLU[ii]+j]; // jj is the index of the jth location who has ii as a neighbor
	        for(k = 0; k < nnIndxLU[J+jj]; k++){ // these are the neighbors of the jjth location
	          kk = nnIndx[nnIndxLU[jj]+k]; // kk is the index for the jth locations neighbors
	          if(kk != ii){ //if the neighbor of jj is not ii
	    	    b += B[ll*nIndx + nnIndxLU[jj]+k]*w[kk * q + ll]; //covariance between jj and kk and the random effect of kk
	          }
	        } // k
	        aij = w[jj * q + ll] - b;
	        a[ll] += B[ll*nIndx + nnIndxLU[jj]+uiIndx[uIndxLU[ii]+j]]*aij/F[ll*J + jj];
	        v[ll] += pow(B[ll * nIndx + nnIndxLU[jj]+uiIndx[uIndxLU[ii]+j]],2)/F[ll * J + jj];
	      } // j
	    }
	    
	    e = 0;
	    for(j = 0; j < nnIndxLU[J+ii]; j++){
	      e += B[ll * nIndx + nnIndxLU[ii]+j]*w[nnIndx[nnIndxLU[ii]+j] * q + ll];
	    }

	    ff[ll] = 1.0 / F[ll * J + ii];
	    gg[ll] = e / F[ll * J + ii];
	  } // ll

	  // var
	  F77_NAME(dcopy)(&qq, tmp_qq, &inc, var, &inc);
	  for (k = 0; k < q; k++) {
            var[k * q + k] += ff[k] + v[k]; 
          } // k
	  F77_NAME(dpotrf)(lower, &q, var, &q, &info FCONE);
          if(info != 0){error("c++ error: dpotrf var failed\n");}
	  F77_NAME(dpotri)(lower, &q, var, &q, &info FCONE);
          if(info != 0){error("c++ error: dpotri var failed\n");}

	  // muNNGP
	  zeros(muNNGP, q);
	  for (j = 0; j < nObs; j++) {
            if (siteIndx[j] == ii) {
              zeros(tmp_nSp, nSp);
	      for (k = 0; k < nSp; k++) {
                tmp_nSp[k] = (yStar[j * nSp + k] - F77_NAME(ddot)(&pAbund, &X[j], &nObs, &beta[k], &nSp) -
			      betaStarSites[k * nObs + j]) * omegaAbund[j * nSp + k];
              } // k (species)
	      F77_NAME(dgemv)(ytran, &nSp, &q, &one, lambda, &nSp, tmp_nSp, 
			      &inc, &one, muNNGP, &inc FCONE);
	    } // if current observation is from the current site. 
	  } // j (observation)

	  for (k = 0; k < q; k++) {
            muNNGP[k] += gg[k] + a[k];
	  } // k

	  F77_NAME(dsymv)(lower, &q, &one, var, &q, muNNGP, &inc, &zero, tmp_nSp, &inc FCONE);

	  F77_NAME(dpotrf)(lower, &q, var, &q, &info FCONE); 
          if(info != 0){error("c++ error: dpotrf var 2 failed\n");}

	  mvrnorm(&w[ii * q], tmp_nSp, var, q);

        } // ii

        /********************************************************************
         *Update spatial factors (lambda)
         *******************************************************************/
        for (i = 1; i < nSp; i++) {
          zeros(tmp_qq, qq);
          zeros(tmp_q, q);
	  zeros(tmp_qq2, qq);
	  // W' %*% S_beta %*% W
          for (k = 0; k < q; k++) {
            for (l = 0; l < q; l++) {
              for (j = 0; j < J; j++) {
                for (ll = 0; ll < nObs; ll++) {
                  if (siteIndx[ll] == j) {
                    tmp_qq[k * q + l] += w[j * q + k] * w[j * q + l] * omegaAbund[ll * nSp + i];
		  }
		} // ll
              } // j
            } // l
          } // k

          // currDim gives the mean dimension of mu and var. 
	  if (i < q) {
            currDim = i;  
          } else {
            currDim = q;
          }
          /*****************************
	   *mu
           *****************************/
	  // NStar - X %*% beta
	  for (j = 0; j < J; j++) {
            for (ll = 0; ll < nObs; ll++) {
              if (siteIndx[ll] == j) {
                tmp_nObs[ll] = yStar[ll * nSp + i] - 
	                   F77_NAME(ddot)(&pAbund, &X[ll], &nObs, &beta[i], &nSp) - 
	                   betaStarSites[i * nObs + ll];
	        if (i < q) {
                  tmp_nObs[ll] -= w[j * q + i];
                }
	      }
	    } // ll (nObs)
          } // j

	  // S_beta %*% W' = tmp_Jq
	  // aka multiply W[j, ] by omegaOcc[j] of the current species you're on. 
	  for (j = 0; j < J; j++) {
            for (ll = 0; ll < nObs; ll++) {
              if (siteIndx[ll] == j) {
                for (l = 0; l < q; l++) {
                  tmp_nObsq[ll * q + l] = omegaAbund[ll * nSp + i] * w[j * q + l];  
                }
	      }
	    }
          }

	  // tmp_Jq %*% tmp_J
	  for (k = 0; k < currDim; k++) {
            for (j = 0; j < nObs; j++) {
              tmp_q[k] += tmp_nObsq[j * q + k] * tmp_nObs[j];
            } // j
          } // k

          /*****************************
	   *var
           *****************************/
	  // Only get relevant columns of t(W) %*% W
	  for (k = 0, l = 0; k < currDim; k++) {
            for (j = 0; j < currDim; j++, l++) {
              tmp_qq2[l] = tmp_qq[k * q + j];
	    } // j
          } // k

	  // Add 1
	  for (j = 0; j < currDim; j++) {
            tmp_qq2[j * currDim + j] += 1.0;  
          } // j

          F77_NAME(dpotrf)(lower, &currDim, tmp_qq2, &currDim, &info FCONE); 
	  if(info != 0){error("c++ error: dpotrf for spatial factors failed\n");}
          F77_NAME(dpotri)(lower, &currDim, tmp_qq2, &currDim, &info FCONE); 
	  if(info != 0){error("c++ error: dpotri for spatial factors failed\n");}

          F77_NAME(dsymv)(lower, &currDim, &one, tmp_qq2, &currDim, tmp_q, &inc, &zero, tmp_q2, &inc FCONE);

          F77_NAME(dpotrf)(lower, &currDim, tmp_qq2, &currDim, &info FCONE); 
	  if(info != 0){error("c++ error: dpotrf for spatial factors 2 failed\n");}
          
          mvrnorm(tmp_q, tmp_q2, tmp_qq2, currDim);
          F77_NAME(dcopy)(&currDim, tmp_q, &inc, &lambda[i], &nSp);
        } // i

        // Multiply Lambda %*% w[j] to get wStar. 
        for (j = 0; j < J; j++) {
          F77_NAME(dgemv)(ntran, &nSp, &q, &one, lambda, &nSp, &w[j*q], &inc, &zero, &wStar[j * nSp], &inc FCONE);
        } // j

        /********************************************************************
         *Update phi (and nu if matern)
         *******************************************************************/
	for (ll = 0; ll < q; ll++) {
          // Current
          if (corName == "matern"){ 
	    nu[ll] = theta[nuIndx * q + ll];
       	  }
          updateBFSFAbund(&B[ll * nIndx], &F[ll*J], &c[ll * m*nThreads], &C[ll * mm * nThreads], coords, nnIndx, nnIndxLU, J, m, theta[sigmaSqIndx * q + ll], theta[phiIndx * q + ll], nu[ll], covModel, &bk[ll * sizeBK], nuB[ll]);
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
      
          logPostCurr = -0.5 * logDet - 0.5 * aa;
          logPostCurr += log(theta[phiIndx * q + ll] - phiA[ll]) + log(phiB[ll] - theta[phiIndx * q + ll]); 
          if(corName == "matern"){
       	    logPostCurr += log(theta[nuIndx * q + ll] - nuA[ll]) + log(nuB[ll] - theta[nuIndx * q + ll]); 
          }
          
          // Candidate
          phiCand = logitInv(rnorm(logit(theta[phiIndx * q + ll], phiA[ll], phiB[ll]), exp(tuning[phiIndx * q + ll])), phiA[ll], phiB[ll]);
          if (corName == "matern"){
      	    nuCand = logitInv(rnorm(logit(theta[nuIndx * q + ll], nuA[ll], nuB[ll]), exp(tuning[nuIndx * q + ll])), nuA[ll], nuB[ll]);
          }
      
          updateBFSFAbund(BCand, FCand, &c[ll * m*nThreads], &C[ll * mm * nThreads], coords, nnIndx, nnIndxLU, J, m, theta[sigmaSqIndx * q + ll], phiCand, nuCand, covModel, &bk[ll * sizeBK], nuB[ll]);
      
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
          
          logPostCand = -0.5*logDet - 0.5*aa;      
          logPostCand += log(phiCand - phiA[ll]) + log(phiB[ll] - phiCand); 
          if (corName == "matern"){
            logPostCand += log(nuCand - nuA[ll]) + log(nuB[ll] - nuCand); 
          }

          if (runif(0.0,1.0) <= exp(logPostCand - logPostCurr)) {

            F77_NAME(dcopy)(&nIndx, BCand, &inc, &B[ll * nIndx], &inc);
            F77_NAME(dcopy)(&J, FCand, &inc, &F[ll * J], &inc);
            
	    theta[phiIndx * q + ll] = phiCand;
            accept[phiIndx * q + ll]++;
            if (corName == "matern") {
              nu[ll] = nuCand; 
	      theta[nuIndx * q + ll] = nu[ll]; 
              accept[nuIndx * q + ll]++; 
            }
          }
	} // ll

        for (i = 0; i < nSp; i++) {	
          /********************************************************************
           *Update kappa (the NB size parameter)
           *******************************************************************/
          for (j = 0; j < nObs; j++) {
            tmp_nObs[j] = F77_NAME(ddot)(&pAbund, &X[j], &nObs, &beta[i], &nSp) + 
	                  betaStarSites[i * nObs + j] +
	                  wStar[siteIndx[j] * nSp + i];
            psi[j] = logitInv(tmp_nObs[j], zero, one);
          }
          /********************************
           * Current
           *******************************/
          // Log gamma function is lgammafn
          // Likelihood contribution
          logPostCurr = 0.0;
          for (j = 0; j < nObs; j++) {
            logPostCurr += lgammafn(y[j * nSp + i] + kappa[i]) - lgammafn(kappa[i]) + kappa[i] * log(1 - psi[j]);
          }
          logPostCurr += log(kappa[i] - kappaA[i]) + log(kappaB[i] - kappa[i]);
          /********************************
           * Candidate
           *******************************/
          kappaCand = logitInv(rnorm(logit(kappa[i], kappaA[i], kappaB[i]), 
                                     exp(tuning[kappaAMCMCIndx * q + i])), kappaA[i], kappaB[i]);
          logPostCand = 0.0;
          for (j = 0; j < nObs; j++) {
            logPostCand += lgammafn(y[j * nSp + i] + kappaCand) - lgammafn(kappaCand) + kappaCand * log(1 - psi[j]);
          }
          logPostCand += log(kappaCand - kappaA[i]) + log(kappaB[i] - kappaCand);

          if (runif(0.0,1.0) <= exp(logPostCand - logPostCurr)) {
            kappa[i] = kappaCand;
            accept[kappaAMCMCIndx * q + i]++;
          }

          /********************************************************************
           *Get fitted values
           *******************************************************************/
          for (j = 0; j < nObs; j++) {
            mu[j * nSp + i] = exp(tmp_nObs[j]) * kappa[i];
            yRep[j * nSp + i] = rnbinom_mu(kappa[i], mu[j * nSp + i]);
	    like[j * nSp + i] = dnbinom_mu(y[j * nSp + i], kappa[i], mu[j * nSp + i], 0);
          }
	} // i (species)

        /********************************************************************
         *Save samples
        *******************************************************************/
        if (g >= nBurn) {
          thinIndx++; 
          if (thinIndx == nThin) {
            F77_NAME(dcopy)(&pAbund, betaComm, &inc, &REAL(betaCommSamples_r)[sPost*pAbund], &inc);
            F77_NAME(dcopy)(&pAbund, tauSqBeta, &inc, &REAL(tauSqBetaSamples_r)[sPost*pAbund], &inc);
            F77_NAME(dcopy)(&pAbundnSp, beta, &inc, &REAL(betaSamples_r)[sPost*pAbundnSp], &inc); 
            F77_NAME(dcopy)(&nSp, kappa, &inc, &REAL(kappaSamples_r)[sPost*nSp], &inc); 
            F77_NAME(dcopy)(&nObsnSp, yRep, &inc, &REAL(yRepSamples_r)[sPost*nObsnSp], &inc); 
            F77_NAME(dcopy)(&nObsnSp, mu, &inc, &REAL(muSamples_r)[sPost*nObsnSp], &inc); 
            F77_NAME(dcopy)(&Jq, w, &inc, &REAL(wSamples_r)[sPost*Jq], &inc); 
            F77_NAME(dcopy)(&nSpq, lambda, &inc, &REAL(lambdaSamples_r)[sPost*nSpq], &inc); 
            F77_NAME(dcopy)(&nThetaqSave, &theta[phiIndx * q], &inc, 
			    &REAL(thetaSamples_r)[sPost*nThetaqSave], &inc); 
            F77_NAME(dcopy)(&nObsnSp, like, &inc, &REAL(likeSamples_r)[sPost*nObsnSp], &inc); 
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
      for (ll = 0; ll < q; ll++) {
        for (k = 0; k < nTheta; k++) {
          REAL(acceptSamples_r)[s * nAMCMC + k * q + ll] = accept[k * q + ll]/batchLength; 
          if (accept[k * q + ll] / batchLength > acceptRate) {
            tuning[k * q + ll] += std::min(0.01, 1.0/sqrt(static_cast<double>(s)));
          } else{
            tuning[k * q + ll] -= std::min(0.01, 1.0/sqrt(static_cast<double>(s)));
            }
          accept[k * q + ll] = 0.0;
        } // k
      } // i
      for (i = 0; i < nSp; i++) {
        REAL(acceptSamples_r)[s * nAMCMC + kappaAMCMCIndx * q + i] = accept[kappaAMCMCIndx * q + i]/batchLength; 
        if (accept[kappaAMCMCIndx * q + i] / batchLength > acceptRate) {
          tuning[kappaAMCMCIndx * q + i] += std::min(0.01, 1.0/sqrt(static_cast<double>(s)));
        } else{
            tuning[kappaAMCMCIndx * q + i] -= std::min(0.01, 1.0/sqrt(static_cast<double>(s)));
          }
        accept[kappaAMCMCIndx * q + i] = 0.0;
      }
      /********************************************************************
       *Report 
       *******************************************************************/
      if (verbose) {
	if (status == nReport) {
	  Rprintf("Batch: %i of %i, %3.2f%%\n", s, nBatch, 100.0*s/nBatch);
          Rprintf("\tNumber\tParameter\tAcceptance\tTuning\n");	  
          for (ll = 0; ll < q; ll++) {
            Rprintf("\t%i\tphi\t\t%3.1f\t\t%1.5f\n", ll + 1, 100.0*REAL(acceptSamples_r)[s * nAMCMC + phiIndx * q + ll], exp(tuning[phiIndx * q + ll]));
	    if (corName == "matern") {
            Rprintf("\t%i\tnu\t\t%3.1f\t\t%1.5f\n", ll + 1, 100.0*REAL(acceptSamples_r)[s * nAMCMC + nuIndx * q + ll], exp(tuning[nuIndx * q + ll]));
	    }
          } // ll
	  for (i = 0; i < nSp; i++) {
	    Rprintf("\t%i\tkappa\t\t%3.1f\t\t%1.5f\n", i + 1, 100.0*REAL(acceptSamples_r)[s * nAMCMC + kappaAMCMCIndx * q + i], exp(tuning[kappaAMCMCIndx * q + i]));
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
    int nResultListObjs = 12;
    if (pAbundRE > 0) {
      nResultListObjs += 2;
    }

    PROTECT(result_r = allocVector(VECSXP, nResultListObjs)); nProtect++;
    PROTECT(resultName_r = allocVector(VECSXP, nResultListObjs)); nProtect++;

    SET_VECTOR_ELT(result_r, 0, betaCommSamples_r);
    SET_VECTOR_ELT(result_r, 1, tauSqBetaSamples_r);
    SET_VECTOR_ELT(result_r, 2, betaSamples_r);
    SET_VECTOR_ELT(result_r, 3, yRepSamples_r);
    SET_VECTOR_ELT(result_r, 4, muSamples_r);
    SET_VECTOR_ELT(result_r, 5, kappaSamples_r);
    SET_VECTOR_ELT(result_r, 6, lambdaSamples_r);
    SET_VECTOR_ELT(result_r, 7, wSamples_r); 
    SET_VECTOR_ELT(result_r, 8, thetaSamples_r); 
    SET_VECTOR_ELT(result_r, 9, tuningSamples_r); 
    SET_VECTOR_ELT(result_r, 10, acceptSamples_r); 
    SET_VECTOR_ELT(result_r, 11, likeSamples_r); 
    if (pAbundRE > 0) {
      SET_VECTOR_ELT(result_r, 12, sigmaSqMuSamples_r);
      SET_VECTOR_ELT(result_r, 13, betaStarSamples_r);
    }

    SET_VECTOR_ELT(resultName_r, 0, mkChar("beta.comm.samples")); 
    SET_VECTOR_ELT(resultName_r, 1, mkChar("tau.sq.beta.samples")); 
    SET_VECTOR_ELT(resultName_r, 2, mkChar("beta.samples")); 
    SET_VECTOR_ELT(resultName_r, 3, mkChar("y.rep.samples")); 
    SET_VECTOR_ELT(resultName_r, 4, mkChar("mu.samples")); 
    SET_VECTOR_ELT(resultName_r, 5, mkChar("kappa.samples")); 
    SET_VECTOR_ELT(resultName_r, 6, mkChar("lambda.samples")); 
    SET_VECTOR_ELT(resultName_r, 7, mkChar("w.samples")); 
    SET_VECTOR_ELT(resultName_r, 8, mkChar("theta.samples")); 
    SET_VECTOR_ELT(resultName_r, 9, mkChar("tune")); 
    SET_VECTOR_ELT(resultName_r, 10, mkChar("accept")); 
    SET_VECTOR_ELT(resultName_r, 11, mkChar("like.samples")); 
    if (pAbundRE > 0) {
      SET_VECTOR_ELT(resultName_r, 12, mkChar("sigma.sq.mu.samples")); 
      SET_VECTOR_ELT(resultName_r, 13, mkChar("beta.star.samples")); 
    }
   
    namesgets(result_r, resultName_r);
    
    UNPROTECT(nProtect);
    
    return(result_r);
  }
}


