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

void updateBFsfGaus(double *B, double *F, double *c, double *C, double *coords, int *nnIndx, int *nnIndxLU, int n, int m, double sigmaSq, double phi, double nu, int covModel, double *bk, double nuUnifb){
    
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
  SEXP sfMsAbundGaussianNNGP(SEXP y_r, SEXP X_r, SEXP coords_r, SEXP XRE_r,
                             SEXP XRandom_r, SEXP consts_r, SEXP nRELong_r, 
			     SEXP m_r, SEXP nnIndx_r, 
                             SEXP nnIndxLU_r, SEXP uIndx_r, SEXP uIndxLU_r, SEXP uiIndx_r, 
                             SEXP betaStarting_r, SEXP betaCommStarting_r,  
                             SEXP tauSqBetaStarting_r, SEXP tauSqStarting_r,
                             SEXP phiStarting_r, SEXP lambdaStarting_r, SEXP nuStarting_r, 
                             SEXP wStarting_r, SEXP sigmaSqMuStarting_r, 
                             SEXP betaStarStarting_r, SEXP betaStarIndx_r, 
			     SEXP betaLevelIndx_r, 
                             SEXP muBetaComm_r, SEXP SigmaBetaComm_r, SEXP tauSqBetaA_r, 
                             SEXP tauSqBetaB_r, SEXP tauSqA_r, SEXP tauSqB_r, SEXP phiA_r, 
                             SEXP phiB_r, SEXP nuA_r, SEXP nuB_r, 
                             SEXP sigmaSqMuA_r, SEXP sigmaSqMuB_r, 
                             SEXP tuning_r, SEXP covModel_r,
			     SEXP nBatch_r, SEXP batchLength_r, 
                             SEXP acceptRate_r, SEXP nThreads_r, 
			     SEXP verbose_r, SEXP nReport_r, 
                             SEXP samplesInfo_r, SEXP chainInfo_r,
			     SEXP z_r, SEXP family_r){
   
    /**********************************************************************
     * Initial constants
     * *******************************************************************/
    // Indexing key (this is general, might deviate in a couple spots): 
    // i = species 
    // j = site
    // ll = latent factor
    // h = coefficients
    // rr = spatially-varying coefficient. 
    int i, j, k, s, g, t, h, l, info, nProtect=0, ii, ll;    

    const int inc = 1;
    const double one = 1.0;
    const double zero = 0.0;
    char const *lower = "L";
    char const *ntran = "N";
    char const *ytran = "T";
    
    /**********************************************************************
     * Get Inputs
     * *******************************************************************/
    // Sorted by site, then by species. 
    double *y = REAL(y_r);
    double *X = REAL(X_r);
    double *XRandom = REAL(XRandom_r);
    // Order: covariate, site
    double *coords = REAL(coords_r); 
    int *XRE = INTEGER(XRE_r);
    int m = INTEGER(m_r)[0]; 
    // Load constants
    int N = INTEGER(consts_r)[0]; 
    int J = INTEGER(consts_r)[1];
    int p = INTEGER(consts_r)[2];
    int pRE = INTEGER(consts_r)[3];
    int nRE = INTEGER(consts_r)[4];
    int q = INTEGER(consts_r)[5]; 
    int indBetas = INTEGER(consts_r)[6];
    int saveFitted = INTEGER(consts_r)[7];
    int pp = p * p; 
    double *muBetaComm = REAL(muBetaComm_r); 
    double *SigmaBetaCommInv = (double *) R_alloc(pp, sizeof(double));   
    F77_NAME(dcopy)(&pp, REAL(SigmaBetaComm_r), &inc, SigmaBetaCommInv, &inc);
    double *tauSqBetaA = REAL(tauSqBetaA_r); 
    double *tauSqBetaB = REAL(tauSqBetaB_r); 
    double *tauSqA = REAL(tauSqA_r); 
    double *tauSqB = REAL(tauSqB_r); 
    double *phiA = REAL(phiA_r); 
    double *phiB = REAL(phiB_r); 
    double *nuA = REAL(nuA_r); 
    double *nuB = REAL(nuB_r); 
    double *sigmaSqMuA = REAL(sigmaSqMuA_r); 
    double *sigmaSqMuB = REAL(sigmaSqMuB_r); 
    double *tuning = REAL(tuning_r); 
    int *nnIndx = INTEGER(nnIndx_r);
    int *nnIndxLU = INTEGER(nnIndxLU_r);
    int *uIndx = INTEGER(uIndx_r);
    int *uIndxLU = INTEGER(uIndxLU_r);
    int *uiIndx = INTEGER(uiIndx_r);
    int covModel = INTEGER(covModel_r)[0];
    std::string corName = getCorName(covModel);
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
    double acceptRate = REAL(acceptRate_r)[0];
    int nThreads = INTEGER(nThreads_r)[0];
    int verbose = INTEGER(verbose_r)[0];
    int nReport = INTEGER(nReport_r)[0];
    int status = 0; 
    int thinIndx = 0; 
    int sPost = 0; 
    // First stage samples
    double *z = REAL(z_r);
    int family = INTEGER(family_r)[0];

#ifdef _OPENMP
    omp_set_num_threads(nThreads);
#else
    if(nThreads > 1){
      warning("n.omp.threads > %i, but source not compiled with OpenMP support.", nThreads);
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
	if (family == 2) {
          Rprintf("Spatial Factor NNGP Multi-species Gaussian Model\nwith %i sites and %i species.\n\n", J, N);
	} else {
          Rprintf("Spatial Factor NNGP Multi-species Zero-Inflated Gaussian Model\nwith %i sites and %i species.\n\n", J, N);
	}
        Rprintf("Samples per chain: %i (%i batches of length %i)\n", nSamples, nBatch, batchLength);
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
    }

    /**********************************************************************
       Some constants and temporary variables to be used later
     * *******************************************************************/
    int pN = p * N; 
    int nREN = nRE * N; 
    int Jq = J * q;
    int qq = q * q;
    int JN = J * N;
    int Nq = N * q;
    int Jp = J * p; 
    int JpRE = J * pRE;
    int jj, kk;
    double tmp_0, tmp_02; 
    double *tmp_one = (double *) R_alloc(inc, sizeof(double)); 
    double *tmp_pp = (double *) R_alloc(pp, sizeof(double)); 
    double *tmp_p = (double *) R_alloc(p, sizeof(double));
    double *tmp_beta = (double *) R_alloc(p, sizeof(double));
    double *tmp_p2 = (double *) R_alloc(p, sizeof(double));
    int *tmp_JInt = (int *) R_alloc(J, sizeof(int));
    for (j = 0; j < J; j++) {
      tmp_JInt[j] = 0; 
    }
    double *tmp_J = (double *) R_alloc(J, sizeof(double));
    zeros(tmp_J, J);
    double *tmp_J1 = (double *) R_alloc(J, sizeof(double));
    double *tmp_Jp = (double *) R_alloc(Jp, sizeof(double));
    double *tmp_qq = (double *) R_alloc(qq, sizeof(double)); zeros(tmp_qq, qq);
    double *tmp_q = (double *) R_alloc(q, sizeof(double)); zeros(tmp_q, q);
    double *tmp_q2 = (double *) R_alloc(q, sizeof(double)); zeros(tmp_q2, q);
    double *tmp_qq2 = (double *) R_alloc(qq, sizeof(double)); zeros(tmp_qq2, qq);
    double *tmp_Jq = (double *) R_alloc(Jq, sizeof(double)); zeros(tmp_Jq, Jq);
    double *tmp_Jq2 = (double *) R_alloc(Jq, sizeof(double)); zeros(tmp_Jq2, Jq);
    double *tmp_Nq = (double *) R_alloc(Nq, sizeof(double)); zeros(tmp_Nq, Nq);
    double *tmp_Nq2 = (double *) R_alloc(Nq, sizeof(double)); zeros(tmp_Nq2, Nq);
    double *tmp_N = (double *) R_alloc(N, sizeof(double)); zeros(tmp_N, N);
    double *tmp_N2 = (double *) R_alloc(N, sizeof(double)); zeros(tmp_N2, N);
    int currDim = 0;

    /**********************************************************************
     * Parameters
     * *******************************************************************/
    // Community level
    double *betaComm = (double *) R_alloc(p, sizeof(double)); 
    F77_NAME(dcopy)(&p, REAL(betaCommStarting_r), &inc, betaComm, &inc);
    double *tauSqBeta = (double *) R_alloc(p, sizeof(double)); 
    F77_NAME(dcopy)(&p, REAL(tauSqBetaStarting_r), &inc, tauSqBeta, &inc);
    // Species level
    double *beta = (double *) R_alloc(pN, sizeof(double));   
    F77_NAME(dcopy)(&pN, REAL(betaStarting_r), &inc, beta, &inc);
    // Nuggets
    double *tauSq = (double *) R_alloc(N, sizeof(double)); 
    F77_NAME(dcopy)(&N, REAL(tauSqStarting_r), &inc, tauSq, &inc);
    // Occurrence random effect variances
    double *sigmaSqMu = (double *) R_alloc(pRE, sizeof(double)); 
    F77_NAME(dcopy)(&pRE, REAL(sigmaSqMuStarting_r), &inc, sigmaSqMu, &inc); 
    // Spatial random effects
    double *w = (double *) R_alloc(Jq, sizeof(double));
    F77_NAME(dcopy)(&Jq, REAL(wStarting_r), &inc, w, &inc); 
    // Latent spatial factor loadings 
    double *lambda = (double *) R_alloc(Nq, sizeof(double));
    F77_NAME(dcopy)(&Nq, REAL(lambdaStarting_r), &inc, lambda, &inc);
    // Latent occurrence random effects
    double *betaStar = (double *) R_alloc(nREN, sizeof(double)); 
    F77_NAME(dcopy)(&nREN, REAL(betaStarStarting_r), &inc, betaStar, &inc); 
    // Spatial range parameter
    double *phi = (double *) R_alloc(q, sizeof(double)); 
    F77_NAME(dcopy)(&q, REAL(phiStarting_r), &inc, phi, &inc); 
    // Spatial smoothing parameter for Matern
    double *nu = (double *) R_alloc(q, sizeof(double)); 
    F77_NAME(dcopy)(&q, REAL(nuStarting_r), &inc, nu, &inc); 
    // Fitted vals
    double *yRep = (double *) R_alloc(JN, sizeof(double)); zeros(yRep, JN);

    /**********************************************************************
     * Return Stuff
     * *******************************************************************/
    // Community level
    SEXP betaCommSamples_r; 
    PROTECT(betaCommSamples_r = allocMatrix(REALSXP, p, nPost)); nProtect++;
    zeros(REAL(betaCommSamples_r), p * nPost);
    SEXP tauSqBetaSamples_r; 
    PROTECT(tauSqBetaSamples_r = allocMatrix(REALSXP, p, nPost)); nProtect++; 
    zeros(REAL(tauSqBetaSamples_r), p * nPost);
    // Species level
    SEXP betaSamples_r;
    PROTECT(betaSamples_r = allocMatrix(REALSXP, pN, nPost)); nProtect++;
    zeros(REAL(betaSamples_r), pN * nPost);
    SEXP tauSqSamples_r;
    PROTECT(tauSqSamples_r = allocMatrix(REALSXP, N, nPost)); nProtect++;
    zeros(REAL(tauSqSamples_r), N * nPost);
    SEXP yRepSamples_r; 
    SEXP muSamples_r; 
    SEXP likeSamples_r;
    if (saveFitted == 1) {
      PROTECT(yRepSamples_r = allocMatrix(REALSXP, JN, nPost)); nProtect++; 
      zeros(REAL(yRepSamples_r), JN * nPost);
      PROTECT(muSamples_r = allocMatrix(REALSXP, JN, nPost)); nProtect++; 
      zeros(REAL(muSamples_r), JN * nPost);
      PROTECT(likeSamples_r = allocMatrix(REALSXP, JN, nPost)); nProtect++;
      zeros(REAL(likeSamples_r), JN * nPost);
    }
    // Spatial parameters
    SEXP lambdaSamples_r; 
    PROTECT(lambdaSamples_r = allocMatrix(REALSXP, Nq, nPost)); nProtect++;
    zeros(REAL(lambdaSamples_r), Nq * nPost);
    SEXP wSamples_r; 
    PROTECT(wSamples_r = allocMatrix(REALSXP, Jq, nPost)); nProtect++; 
    zeros(REAL(wSamples_r), Jq * nPost);
    // Random effects
    SEXP sigmaSqMuSamples_r; 
    SEXP betaStarSamples_r; 
    if (pRE > 0) {
      PROTECT(sigmaSqMuSamples_r = allocMatrix(REALSXP, pRE, nPost)); nProtect++;
      zeros(REAL(sigmaSqMuSamples_r), pRE * nPost);
      PROTECT(betaStarSamples_r = allocMatrix(REALSXP, nREN, nPost)); nProtect++;
      zeros(REAL(betaStarSamples_r), nREN * nPost);
    }
    
    /**********************************************************************
     * Additional Sampler Prep
     * *******************************************************************/
    double *mu = (double *) R_alloc(JN, sizeof(double)); 
    zeros(mu, JN); 
    double *like = (double *) R_alloc(JN, sizeof(double)); zeros(like, JN); 

    // For normal community-level priors
    // Occurrence coefficients
    F77_NAME(dpotrf)(lower, &p, SigmaBetaCommInv, &p, &info FCONE); 
    if(info != 0){error("c++ error: dpotrf SigmaBetaCommInv failed\n");}
    F77_NAME(dpotri)(lower, &p, SigmaBetaCommInv, &p, &info FCONE); 
    if(info != 0){error("c++ error: dpotri SigmaBetaCommInv failed\n");}
    double *SigmaBetaCommInvMuBeta = (double *) R_alloc(p, sizeof(double)); 
    F77_NAME(dsymv)(lower, &p, &one, SigmaBetaCommInv, &p, muBetaComm, &inc, &zero, 
        	    SigmaBetaCommInvMuBeta, &inc FCONE);
    // Put community level occurrence variances in a p x p matrix.
    double *TauBetaInv = (double *) R_alloc(pp, sizeof(double)); zeros(TauBetaInv, pp); 
    for (i = 0; i < p; i++) {
      TauBetaInv[i * p + i] = tauSqBeta[i]; 
    } // i
    F77_NAME(dpotrf)(lower, &p, TauBetaInv, &p, &info FCONE); 
    if(info != 0){error("c++ error: dpotrf TauBetaInv failed\n");}
    F77_NAME(dpotri)(lower, &p, TauBetaInv, &p, &info FCONE); 
    if(info != 0){error("c++ error: dpotri TauBetaInv failed\n");}

    /**********************************************************************
     * Prep for random effects (if they exist)
     * *******************************************************************/
    // Site-level sums of the occurrence random effects
    double *betaStarSites = (double *) R_alloc(JN, sizeof(double)); 
    zeros(betaStarSites, JN); 
    int *betaStarLongIndx = (int *) R_alloc(JpRE, sizeof(int));
    // Initial sums (initiate with the first species)
    for (j = 0; j < J; j++) {
      for (l = 0; l < pRE; l++) {
        betaStarLongIndx[l * J + j] = which(XRE[l * J + j], betaLevelIndx, nRE);
        for (i = 0; i < N; i++) {
          betaStarSites[i * J + j] += betaStar[i * nRE + betaStarLongIndx[l * J + j]] * XRandom[l * J + j];
        }
      }
    }

    // Starting index for occurrence random effects
    int *betaStarStart = (int *) R_alloc(pRE, sizeof(int)); 
    for (l = 0; l < pRE; l++) {
      betaStarStart[l] = which(l, betaStarIndx, nRE); 
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
    zeros(REAL(thetaSamples_r), nThetaqSave * nPost);
    // Species-level spatial random effects
    double *wStar = (double *) R_alloc(JN, sizeof(double)); zeros(wStar, JN);
    // Multiply Lambda %*% w[j] to get wStar. 
    for (j = 0; j < J; j++) {
      F77_NAME(dgemv)(ntran, &N, &q, &one, lambda, &N, &w[j*q], &inc, &zero, &wStar[j * N], &inc FCONE);
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
      updateBFsfGaus(&B[ll * nIndx], &F[ll*J], &c[ll * m*nThreads], &C[ll * mm * nThreads], coords, nnIndx, nnIndxLU, J, m, theta[sigmaSqIndx * q + ll], theta[phiIndx * q + ll], nu[ll], covModel, &bk[ll * sizeBK], nuB[0]);
    }

    /**********************************************************************
     Set up stuff for Adaptive MH and other misc
     * *******************************************************************/
    double logPostCurr = 0.0, logPostCand = 0.0; 
    logPostCurr = R_NegInf; 
    double *accept = (double *) R_alloc(nThetaq, sizeof(double)); zeros(accept, nThetaq); 
    double phiCand = 0.0, nuCand = 0.0; 
    double logDet; 
    // MCMC info if desired
    double *accept2 = (double *) R_alloc(nThetaq, sizeof(double)); zeros(accept2, nThetaq);
    SEXP tuningSamples_r; 
    PROTECT(tuningSamples_r = allocMatrix(REALSXP, nThetaq, nBatch)); nProtect++; 
    zeros(REAL(tuningSamples_r), nThetaq * nBatch);
    // For current number of nonzero z values
    double *currJ = (double *) R_alloc(N, sizeof(double)); zeros(currJ, N);
    for (i = 0; i < N; i++) {
      for (j = 0; j < J; j++) {
        currJ[i] += z[j * N + i];
      }
    }

    GetRNGstate();

    /**********************************************************************
     Start sampling
     * *******************************************************************/
    for (s = 0, g = 0; s < nBatch; s++) {
      for (t = 0; t < batchLength; t++, g++) {

        /********************************************************************
         Update Community level Coefficients
         *******************************************************************/
        if (indBetas == 0) {
          /********************************
           Compute b.beta.comm
           *******************************/
          zeros(tmp_p, p); 
          for (i = 0; i < N; i++) {
            F77_NAME(dgemv)(ytran, &p, &p, &one, TauBetaInv, &p, &beta[i], &N, &one, tmp_p, &inc FCONE); 
          } // i
          for (h = 0; h < p; h++) {
            tmp_p[h] += SigmaBetaCommInvMuBeta[h];  
          } // j

          /********************************
           Compute A.beta.comm
           *******************************/
          for (h = 0; h < pp; h++) {
            tmp_pp[h] = SigmaBetaCommInv[h] + N * TauBetaInv[h]; 
          }
          F77_NAME(dpotrf)(lower, &p, tmp_pp, &p, &info FCONE); 
          if(info != 0){error("c++ error: dpotrf ABetaComm failed\n");}
          F77_NAME(dpotri)(lower, &p, tmp_pp, &p, &info FCONE); 
          if(info != 0){error("c++ error: dpotri ABetaComm failed\n");}
          // A.beta.inv %*% b.beta
          // 1 * tmp_pp * tmp_p + 0 * tmp_p2  = tmp_p2
          F77_NAME(dsymv)(lower, &p, &one, tmp_pp, &p, tmp_p, &inc, &zero, tmp_p2, &inc FCONE);
          // Computes cholesky of tmp_pp again stored back in tmp_pp. This chol(A.beta.inv)
          F77_NAME(dpotrf)(lower, &p, tmp_pp, &p, &info FCONE); 
          if(info != 0){error("c++ error: dpotrf ABetaComm failed\n");}
          // Args: destination, mu, cholesky of the inverse covariance matrix, dimension
          mvrnorm(betaComm, tmp_p2, tmp_pp, p);
	}
        /********************************************************************
         Update Community Occupancy Variance Parameter
        ********************************************************************/
        if (indBetas == 0) {
          for (h = 0; h < p; h++) {
            tmp_0 = 0.0;  
            for (i = 0; i < N; i++) {
              tmp_0 += (beta[h * N + i] - betaComm[h]) * (beta[h * N + i] - betaComm[h]);
            } // i
            tmp_0 *= 0.5;
            tauSqBeta[h] = rigamma(tauSqBetaA[h] + N / 2.0, tauSqBetaB[h] + tmp_0); 
          } // h
          for (h = 0; h < p; h++) {
            TauBetaInv[h * p + h] = tauSqBeta[h]; 
          } // i
          F77_NAME(dpotrf)(lower, &p, TauBetaInv, &p, &info FCONE); 
          if(info != 0){error("c++ error: dpotrf TauBetaInv failed\n");}
          F77_NAME(dpotri)(lower, &p, TauBetaInv, &p, &info FCONE); 
          if(info != 0){error("c++ error: dpotri TauBetaInv failed\n");}
	}
        /********************************************************************
         *Update Occupancy random effects variance
         *******************************************************************/
        for (l = 0; l < pRE; l++) {
          tmp_0 = 0.0; 
          for (i = 0; i < N; i++) {
            tmp_0 += F77_NAME(ddot)(&nRELong[l], &betaStar[i*nRE + betaStarStart[l]], &inc, &betaStar[i*nRE + betaStarStart[l]], &inc); 
          }
          tmp_0 *= 0.5; 
          sigmaSqMu[l] = rigamma(sigmaSqMuA[l] + nRELong[l] * N / 2.0, sigmaSqMuB[l] + tmp_0);
        }
        /********************************************************************
         *Update Species-Specific Regression Parameters
         *******************************************************************/
        for (i = 0; i < N; i++) {  
          /********************************************************************
           *Update Occupancy Regression Coefficients
           *******************************************************************/
          zeros(tmp_J1, J);
          for (j = 0; j < J; j++) {
            if (z[j * N + i] == 1.0) {
	      tmp_J1[j] = (y[j * N + i] - wStar[j * N + i] - 
	        	   betaStarSites[i * J + j]) / tauSq[i];
	    }
          } // j
          /********************************
           * Compute b.beta
           *******************************/
          // t(X) * tmp_J1 + 0 * tmp_p = tmp_p. 
          F77_NAME(dgemv)(ytran, &J, &p, &one, X, &J, tmp_J1, &inc, &zero, tmp_p, &inc FCONE); 	 
          // TauBetaInv %*% betaComm + tmp_p = tmp_p
          F77_NAME(dgemv)(ntran, &p, &p, &one, TauBetaInv, &p, betaComm, &inc, &one, tmp_p, &inc FCONE); 

          /********************************
           * Compute A.beta
           * *****************************/
          // t(X) %*% diag(omegaOcc)
	  zeros(tmp_Jp, Jp);
          for(j = 0; j < J; j++){
            if (z[j * N + i] == 1.0) {
              for(h = 0; h < p; h++){
                tmp_Jp[h*J+j] = X[h*J+j] / tauSq[i];
              }
	    }
          }
          // This finishes off A.beta
          // 1 * X * tmp_Jp + 0 * tmp_pp = tmp_pp
          F77_NAME(dgemm)(ytran, ntran, &p, &p, &J, &one, X, &J, tmp_Jp, &J, &zero, tmp_pp, &p FCONE FCONE);
          for (h = 0; h < pp; h++) {
            tmp_pp[h] += TauBetaInv[h]; 
          } // j
          F77_NAME(dpotrf)(lower, &p, tmp_pp, &p, &info FCONE); 
          if(info != 0){error("c++ error: dpotrf ABeta failed\n");}
          F77_NAME(dpotri)(lower, &p, tmp_pp, &p, &info FCONE); 
          if(info != 0){error("c++ error: dpotri ABeta failed\n");}
          // A.beta.inv %*% b.beta
          F77_NAME(dsymv)(lower, &p, &one, tmp_pp, &p, tmp_p, &inc, &zero, tmp_p2, &inc FCONE);
          F77_NAME(dpotrf)(lower, &p, tmp_pp, &p, &info FCONE); 
	  if(info != 0){error("c++ error: dpotrf A.beta 2 failed\n");}
          // Args: destination, mu, cholesky of the covariance matrix, dimension
          mvrnorm(tmp_beta, tmp_p2, tmp_pp, p);
          // Can eventually get rid of this and change order of beta. 
          for (h = 0; h < p; h++) {
            beta[h * N + i] = tmp_beta[h]; 
          }
          /********************************************************************
           *Update tau.sq
           *******************************************************************/
	   zeros(tmp_J1, J);
           for(j = 0; j < J; j++){
             if (z[j * N + i] == 1.0) {
	       tmp_J1[j] = y[j * N + i] - wStar[j * N + i] - 
	         	  F77_NAME(ddot)(&p, &X[j], &J, &beta[i], &N) - 
	         	  betaStarSites[i * J + j];
	     }
           }
           tauSq[i] = rigamma(tauSqA[i] + currJ[i] / 2.0, tauSqB[i] + 0.5 * 
                              F77_NAME(ddot)(&J, tmp_J1, &inc, tmp_J1, &inc));
          /********************************************************************
           *Update Occupancy random effects
           *******************************************************************/
	  if (pRE > 0) {
            // Reset betaStar to 0 since some levels may not be sampled for a given species
	    // zeros(betaStar, nREN);
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
                if ((z[j * N + i] == 1.0) && (XRE[betaStarIndx[l] * J + j] == betaLevelIndx[l])) {
                  tmp_02 = 0.0;
		  for (ll = 0; ll < pRE; ll++) {
                    tmp_02 += betaStar[i * nRE + betaStarLongIndx[ll * J + j]] * XRandom[ll * J + j];
		  }
                  tmp_one[0] += XRandom[betaStarIndx[l] * J + j] * (y[j * N + i] - F77_NAME(ddot)(&p, &X[j], &J, &beta[i], &N) - wStar[j * N + i] + (betaStar[i * nRE + l] * XRandom[betaStarIndx[l] * J + j]) - tmp_02) / tauSq[i];
                tmp_0 += XRandom[betaStarIndx[l] * J + j] * XRandom[betaStarIndx[l] * J + j] / tauSq[i];
	        }
              }
              /********************************
               * Compute A.beta.star
               *******************************/
              tmp_0 += 1.0 / sigmaSqMu[betaStarIndx[l]]; 
              tmp_0 = 1.0 / tmp_0; 
              betaStar[i * nRE + l] = rnorm(tmp_0 * tmp_one[0], sqrt(tmp_0)); 
            }

            // Update the RE sums for the current species
            zeros(&betaStarSites[i * J], J);
            for (j = 0; j < J; j++) {
              if (z[j * N + i] == 1.0) {
                for (l = 0; l < pRE; l++) {
                  betaStarSites[i * J + j] += betaStar[i * nRE + betaStarLongIndx[l * J + j]] * XRandom[l * J + j];
                }
	      }
            }
	  }
	} // i
        /********************************************************************
         *Update Spatial Random Effects (w)
         *******************************************************************/
	// Update B and F
        for (ll = 0; ll < q; ll++) {
          updateBFsfGaus(&B[ll * nIndx], &F[ll*J], &c[ll * m*nThreads], &C[ll * mm * nThreads], coords, nnIndx, nnIndxLU, J, m, theta[sigmaSqIndx * q + ll], theta[phiIndx * q + ll], nu[ll], covModel, &bk[ll * sizeBK], nuB[0]);
        }

	for (ii = 0; ii < J; ii++) {
          zeros(tmp_Nq, Nq);
          // tmp_qq = lambda' S_beta lambda 
	  for (i = 0; i < N; i++) {
            if (z[ii * N + i] == 1.0) {
              for (ll = 0; ll < q; ll++) {
                tmp_Nq[ll * N + i] = lambda[ll * N + i] / tauSq[i];
              } // ll
	    }
          } // i
	  F77_NAME(dgemm)(ytran, ntran, &q, &q, &N, &one, tmp_Nq, &N, lambda, &N, &zero, tmp_qq, &q FCONE FCONE);

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

	  // mu
	  zeros(tmp_N, N);
	  for (k = 0; k < N; k++) {
            if (z[ii * N + k] == 1.0) {
              tmp_N[k] = (y[ii * N + k] - F77_NAME(ddot)(&p, &X[ii], &J, &beta[k], &N) - betaStarSites[k * J + ii]) / tauSq[k];
	    }
          } // k

	  F77_NAME(dgemv)(ytran, &N, &q, &one, lambda, &N, tmp_N, &inc, &zero, muNNGP, &inc FCONE);

	  for (k = 0; k < q; k++) {
            muNNGP[k] += gg[k] + a[k];
	  } // k

	  F77_NAME(dsymv)(lower, &q, &one, var, &q, muNNGP, &inc, &zero, tmp_N, &inc FCONE);

	  F77_NAME(dpotrf)(lower, &q, var, &q, &info FCONE); 
          if(info != 0){error("c++ error: dpotrf var 2 failed\n");}

	  mvrnorm(&w[ii * q], tmp_N, var, q);
        } // ii

        /********************************************************************
         *Update spatial factors (lambda)
         *******************************************************************/
        for (i = 1; i < N; i++) {
          zeros(tmp_qq, qq);
          zeros(tmp_q, q);
	  zeros(tmp_qq2, qq);
	  // W' %*% S_beta %*% W
          for (k = 0; k < q; k++) {
            for (l = 0; l < q; l++) {
              for (j = 0; j < J; j++) {
                if (z[j * N + i] == 1.0) {
                  tmp_qq[k * q + l] += w[j * q + k] * w[j * q + l] / tauSq[i]; 
		}
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
	  // zStar - X %*% beta
	  zeros(tmp_J, J);
	  zeros(tmp_Jq, Jq);
	  for (j = 0; j < J; j++) {
            if (z[j * N + i] == 1.0) {
              tmp_J[j] = y[j * N + i] - F77_NAME(ddot)(&p, &X[j], &J, &beta[i], &N) - 
	                 betaStarSites[i * J + j];
	      if (i < q) {
                tmp_J[j] -= w[j * q + i];
              }
	    }
          } // j

	  // S_beta %*% W' = tmp_Jq
	  // aka multiply W[j, ] by omegaOcc[j] of the current species you're on. 
	  for (j = 0; j < J; j++) {
            if (z[j * N + i] == 1.0) {
              for (ll = 0; ll < q; ll++) {
                tmp_Jq[j * q + ll] = w[j * q + ll] / tauSq[i];  
              }
	    }
          }

	  // tmp_Jq %*% tmp_J
	  for (k = 0; k < currDim; k++) {
            for (j = 0; j < J; j++) {
              if (z[j * N + i] == 1.0) {
                tmp_q[k] += tmp_Jq[j * q + k] * tmp_J[j];
	      }
            } // j
          } // k

          /*****************************
	   *var
           *****************************/
	  // Only get relevant columns of t(W) %*% W
	  for (k = 0; k < currDim; k++) {
            for (j = 0; j < currDim; j++) {
              tmp_qq2[k * currDim + j] = tmp_qq[k * q + j];
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

          F77_NAME(dsymv)(lower, &currDim, &one, tmp_qq2, &currDim, tmp_q, 
			  &inc, &zero, tmp_q2, &inc FCONE);

          F77_NAME(dpotrf)(lower, &currDim, tmp_qq2, &currDim, &info FCONE); 
	  if(info != 0){error("c++ error: dpotrf for spatial factors 2 failed\n");}
          
          mvrnorm(tmp_q, tmp_q2, tmp_qq2, currDim);
          F77_NAME(dcopy)(&currDim, tmp_q, &inc, &lambda[i], &N);
        } // i

        // Multiply Lambda %*% w[j] to get wStar. 
        for (j = 0; j < J; j++) {
          F77_NAME(dgemv)(ntran, &N, &q, &one, lambda, &N, &w[j*q], &inc, &zero, &wStar[j * N], &inc FCONE);
        } // j

        /********************************************************************
         *Update phi (and nu if matern)
         *******************************************************************/
	for (ll = 0; ll < q; ll++) {
          // Current
          if (corName == "matern"){ 
	    nu[ll] = theta[nuIndx * q + ll];
       	  }
          updateBFsfGaus(&B[ll * nIndx], &F[ll*J], &c[ll * m*nThreads], &C[ll * mm * nThreads], coords, nnIndx, nnIndxLU, J, m, theta[sigmaSqIndx * q + ll], theta[phiIndx * q + ll], nu[ll], covModel, &bk[ll * sizeBK], nuB[ll]);
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
      
          updateBFsfGaus(BCand, FCand, &c[ll * m*nThreads], &C[ll * mm * nThreads], coords, nnIndx, nnIndxLU, J, m, theta[sigmaSqIndx * q + ll], phiCand, nuCand, covModel, &bk[ll * sizeBK], nuB[ll]);
      
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

        /********************************************************************
         *Get fitted values and likelihood for WAIC
         *******************************************************************/
	if (saveFitted == 1) {
	  for (i = 0; i < N; i++) {
            for (j = 0; j < J; j++) {
              if (z[j * N + i] == 1.0) {
                mu[j * N + i] = F77_NAME(ddot)(&p, &X[j], &J, &beta[i], &N) + wStar[j * N + i] + betaStarSites[i * J + j];
                yRep[j * N + i] = rnorm(mu[j * N + i], sqrt(tauSq[i]));
                like[j * N + i] = dnorm(y[j * N + i], mu[j * N + i], sqrt(tauSq[i]), 0);
	      } else {
                mu[j * N + i] = 0.0;
	        yRep[j * N + i] = rnorm(mu[j * N + i], sqrt(0.0001));
	        like[j * N + i] = 1.0;
	      }
	    } // j
	  } // i
        }

        /********************************************************************
         *Save samples
         *******************************************************************/
        if (g >= nBurn) {
          thinIndx++;
          if (thinIndx == nThin) {
            F77_NAME(dcopy)(&p, betaComm, &inc, &REAL(betaCommSamples_r)[sPost*p], &inc);
            F77_NAME(dcopy)(&p, tauSqBeta, &inc, &REAL(tauSqBetaSamples_r)[sPost*p], &inc);
            F77_NAME(dcopy)(&N, tauSq, &inc, &REAL(tauSqSamples_r)[sPost*N], &inc);
            F77_NAME(dcopy)(&pN, beta, &inc, &REAL(betaSamples_r)[sPost*pN], &inc); 
            F77_NAME(dcopy)(&Nq, lambda, &inc, &REAL(lambdaSamples_r)[sPost*Nq], &inc); 
	    if (saveFitted == 1) {
              F77_NAME(dcopy)(&JN, mu, &inc, &REAL(muSamples_r)[sPost*JN], &inc); 
              F77_NAME(dcopy)(&JN, yRep, &inc, &REAL(yRepSamples_r)[sPost*JN], &inc); 
              F77_NAME(dcopy)(&JN, like, &inc, &REAL(likeSamples_r)[sPost*JN], &inc); 
	    }
            F77_NAME(dcopy)(&Jq, w, &inc, &REAL(wSamples_r)[sPost*Jq], &inc); 
            F77_NAME(dcopy)(&nThetaqSave, &theta[phiIndx * q], &inc, &REAL(thetaSamples_r)[sPost*nThetaqSave], &inc); 
	    if (pRE > 0) {
              F77_NAME(dcopy)(&pRE, sigmaSqMu, &inc, &REAL(sigmaSqMuSamples_r)[sPost*pRE], &inc);
              F77_NAME(dcopy)(&nREN, betaStar, &inc, &REAL(betaStarSamples_r)[sPost*nREN], &inc);
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
	  accept2[k * q + ll] = accept[k * q + ll] / batchLength;
          if (accept[k * q + ll] / batchLength > acceptRate) {
            tuning[k * q + ll] += std::min(0.01, 1.0/sqrt(static_cast<double>(s)));
          } else{
              tuning[k * q + ll] -= std::min(0.01, 1.0/sqrt(static_cast<double>(s)));
            }
          accept[k * q + ll] = 0.0;
          if (s == (nBatch - 1)) {
          REAL(tuningSamples_r)[k * q + ll] = tuning[k * q + ll]; 
	  }
        } // k
      } // i

      /********************************************************************
       *Report 
       *******************************************************************/
      if (verbose) {
        if (status == nReport) {
          Rprintf("Batch: %i of %i, %3.2f%%\n", s, nBatch, 100.0*s/nBatch);
          Rprintf("\tLatent Factor\tParameter\tAcceptance\tTuning\n");	  
          for (ll = 0; ll < q; ll++) {
            Rprintf("\t%i\t\tphi\t\t%3.1f\t\t%1.5f\n", ll + 1, 100.0*accept2[phiIndx * q + ll], exp(tuning[phiIndx * q + ll]));
	    if (corName == "matern") {
            Rprintf("\t%i\t\tnu\t\t%3.1f\t\t%1.5f\n", ll + 1, 100.0*accept2[nuIndx * q + ll], exp(tuning[nuIndx * q + ll]));
	    }
          } // ll
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
  
    // This is necessary when generating random numbers in C.     
    PutRNGstate();

    // make return object (which is a list)
    SEXP result_r, resultName_r;
    int nResultListObjs = 10;
    if (pRE > 0) {
      nResultListObjs += 2;
    }

    PROTECT(result_r = allocVector(VECSXP, nResultListObjs)); nProtect++;
    PROTECT(resultName_r = allocVector(VECSXP, nResultListObjs)); nProtect++;

    // Setting the components of the output list.
    SET_VECTOR_ELT(result_r, 0, betaCommSamples_r);
    SET_VECTOR_ELT(result_r, 1, tauSqBetaSamples_r);
    SET_VECTOR_ELT(result_r, 2, betaSamples_r);
    if (saveFitted == 1) {
      SET_VECTOR_ELT(result_r, 3, yRepSamples_r);
      SET_VECTOR_ELT(result_r, 4, muSamples_r);
      SET_VECTOR_ELT(result_r, 9, likeSamples_r); 
    }
    SET_VECTOR_ELT(result_r, 5, lambdaSamples_r);
    SET_VECTOR_ELT(result_r, 6, wSamples_r); 
    SET_VECTOR_ELT(result_r, 7, thetaSamples_r); 
    SET_VECTOR_ELT(result_r, 8, tauSqSamples_r); 
    if (pRE > 0) {
      SET_VECTOR_ELT(result_r, 10, sigmaSqMuSamples_r);
      SET_VECTOR_ELT(result_r, 11, betaStarSamples_r);
    }

    // mkChar turns a C string into a CHARSXP
    SET_VECTOR_ELT(resultName_r, 0, mkChar("beta.comm.samples")); 
    SET_VECTOR_ELT(resultName_r, 1, mkChar("tau.sq.beta.samples")); 
    SET_VECTOR_ELT(resultName_r, 2, mkChar("beta.samples")); 
    if (saveFitted == 1) {
      SET_VECTOR_ELT(resultName_r, 3, mkChar("y.rep.samples")); 
      SET_VECTOR_ELT(resultName_r, 4, mkChar("mu.samples")); 
      SET_VECTOR_ELT(resultName_r, 9, mkChar("like.samples")); 
    }
    SET_VECTOR_ELT(resultName_r, 5, mkChar("lambda.samples")); 
    SET_VECTOR_ELT(resultName_r, 6, mkChar("w.samples")); 
    SET_VECTOR_ELT(resultName_r, 7, mkChar("theta.samples")); 
    SET_VECTOR_ELT(resultName_r, 8, mkChar("tau.sq.samples")); 
    if (pRE > 0) {
      SET_VECTOR_ELT(resultName_r, 10, mkChar("sigma.sq.mu.samples")); 
      SET_VECTOR_ELT(resultName_r, 11, mkChar("beta.star.samples")); 
    }
   
    // Set the names of the output list.  
    namesgets(result_r, resultName_r);
    
    //unprotect
    UNPROTECT(nProtect);
    
    return(result_r);
  }
}


