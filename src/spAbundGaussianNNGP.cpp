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

void updateBFGaus(double *B, double *F, double *c, double *C, double *coords, int *nnIndx, int *nnIndxLU, int n, int m, double sigmaSq, double phi, double nu, int covModel, double *bk, double nuUnifb){

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
  SEXP spAbundGaussianNNGP(SEXP y_r, SEXP X_r, SEXP coords_r, SEXP XRE_r, SEXP XRandom_r,
                           SEXP consts_r, SEXP nRELong_r, SEXP m_r, SEXP nnIndx_r, 
                           SEXP nnIndxLU_r, SEXP uIndx_r, SEXP uIndxLU_r, SEXP uiIndx_r,
                           SEXP betaStarting_r, SEXP tauSqStarting_r, SEXP sigmaSqMuStarting_r,
                           SEXP betaStarStarting_r,  
                           SEXP wStarting_r, SEXP phiStarting_r, 
                           SEXP sigmaSqStarting_r, SEXP nuStarting_r, 
                           SEXP betaStarIndx_r, SEXP betaLevelIndx_r, 
                           SEXP muBeta_r, SEXP SigmaBeta_r, 
	                   SEXP tauSqA_r, SEXP tauSqB_r, SEXP phiA_r, SEXP phiB_r, 
                           SEXP sigmaSqA_r, SEXP sigmaSqB_r, SEXP nuA_r, SEXP nuB_r, 
                           SEXP sigmaSqMuA_r, SEXP sigmaSqMuB_r, 
                           SEXP tuning_r, SEXP covModel_r, SEXP nBatch_r, 
                           SEXP batchLength_r, SEXP acceptRate_r, 
			   SEXP nThreads_r, SEXP verbose_r, 
                           SEXP nReport_r, SEXP samplesInfo_r, SEXP chainInfo_r, 
			   SEXP sigmaSqIG_r){
   
    /**********************************************************************
     * Initial constants
     * *******************************************************************/
    int i, j, l, ll, k, s, r, q, info, nProtect=0;
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
    int m = INTEGER(m_r)[0]; 
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
    double phiA = REAL(phiA_r)[0];
    double phiB = REAL(phiB_r)[0]; 
    double nuA = REAL(nuA_r)[0]; 
    double nuB = REAL(nuB_r)[0]; 
    double sigmaSqA = REAL(sigmaSqA_r)[0]; 
    double sigmaSqB = REAL(sigmaSqB_r)[0]; 
    double tauSqA = REAL(tauSqA_r)[0]; 
    double tauSqB = REAL(tauSqB_r)[0]; 
    double *sigmaSqMuA = REAL(sigmaSqMuA_r); 
    double *sigmaSqMuB = REAL(sigmaSqMuB_r); 
    double *tuning = REAL(tuning_r); 
    double *coords = REAL(coords_r);
    int *nRELong = INTEGER(nRELong_r); 
    int *nnIndx = INTEGER(nnIndx_r);
    int *nnIndxLU = INTEGER(nnIndxLU_r);
    int *uIndx = INTEGER(uIndx_r);
    int *uIndxLU = INTEGER(uIndxLU_r);
    int *uiIndx = INTEGER(uiIndx_r);
    int covModel = INTEGER(covModel_r)[0];
    std::string corName = getCorName(covModel);
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
    int sigmaSqIG = INTEGER(sigmaSqIG_r)[0];
    int thinIndx = 0; 
    int sPost = 0; 

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
	if (JZero > 0) {
          Rprintf("Spatial NNGP zero-inflated Gaussian model with %i non-zero sites.\n\n", J);
	} else {
          Rprintf("Spatial NNGP Gaussian model with %i sites.\n\n", J);
	}
        Rprintf("Samples per chain: %i (%i batches of length %i)\n", nSamples, nBatch, batchLength);
        Rprintf("Burn-in: %i \n", nBurn); 
        Rprintf("Thinning Rate: %i \n", nThin); 
        Rprintf("Number of Chains: %i \n", nChain);
        Rprintf("Total Posterior Samples: %i \n\n", nPost * nChain); 
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
    double *beta = (double *) R_alloc(p, sizeof(double));   
    F77_NAME(dcopy)(&p, REAL(betaStarting_r), &inc, beta, &inc);
    // Occupancy random effect variances
    double *sigmaSqMu = (double *) R_alloc(pRE, sizeof(double)); 
    F77_NAME(dcopy)(&pRE, REAL(sigmaSqMuStarting_r), &inc, sigmaSqMu, &inc); 
    // Latent occupancy random effects
    double *betaStar = (double *) R_alloc(nRE, sizeof(double)); 
    F77_NAME(dcopy)(&nRE, REAL(betaStarStarting_r), &inc, betaStar, &inc); 
    // Spatial processes
    double *w = (double *) R_alloc(J, sizeof(double));   
    F77_NAME(dcopy)(&J, REAL(wStarting_r), &inc, w, &inc);
    // Spatial variance
    double sigmaSq = REAL(sigmaSqStarting_r)[0];
    // Nugget
    double tauSq = REAL(tauSqStarting_r)[0];
    // Spatial range parameter
    double phi = REAL(phiStarting_r)[0];
    // Spatial smoothing parameter for Matern
    double nu = REAL(nuStarting_r)[0];
    
    /**********************************************************************
     * Return Stuff
     * *******************************************************************/
    SEXP betaSamples_r;
    PROTECT(betaSamples_r = allocMatrix(REALSXP, p, nPost)); nProtect++;
    zeros(REAL(betaSamples_r), p * nPost);
    SEXP yRepSamples_r; 
    SEXP yRepZeroSamples_r;
    SEXP muSamples_r; 
    SEXP likeSamples_r;
    if (saveFitted == 1) {
      PROTECT(yRepSamples_r = allocMatrix(REALSXP, J, nPost)); nProtect++; 
      zeros(REAL(yRepSamples_r), J * nPost);
      PROTECT(yRepZeroSamples_r = allocMatrix(REALSXP, JZero, nPost)); nProtect++; 
      zeros(REAL(yRepZeroSamples_r), JZero * nPost);
      PROTECT(muSamples_r = allocMatrix(REALSXP, J, nPost)); nProtect++; 
      zeros(REAL(muSamples_r), J * nPost);
      PROTECT(likeSamples_r = allocMatrix(REALSXP, J, nPost)); nProtect++;
      zeros(REAL(likeSamples_r), J * nPost);
    }
    SEXP wSamples_r; 
    PROTECT(wSamples_r = allocMatrix(REALSXP, J, nPost)); nProtect++; 
    zeros(REAL(wSamples_r), J * nPost);
    // Occurrence random effects
    SEXP sigmaSqMuSamples_r; 
    SEXP betaStarSamples_r; 
    if (pRE > 0) {
      PROTECT(sigmaSqMuSamples_r = allocMatrix(REALSXP, pRE, nPost)); nProtect++;
      zeros(REAL(sigmaSqMuSamples_r), pRE * nPost);
      PROTECT(betaStarSamples_r = allocMatrix(REALSXP, nRE, nPost)); nProtect++;
      zeros(REAL(betaStarSamples_r), nRE * nPost);
    }
    SEXP tauSqSamples_r;
    PROTECT(tauSqSamples_r = allocMatrix(REALSXP, inc, nPost)); nProtect++;
    zeros(REAL(tauSqSamples_r), nPost);
    
    /**********************************************************************
     * Other initial starting stuff
     * *******************************************************************/
    int Jp = J * p; 
    int jj, kk;
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
    if(info != 0){error("c++ error: dpotrf SigmaBetaInv failed\n");}
    F77_NAME(dpotri)(lower, &p, SigmaBetaInv, &p, &info FCONE); 
    if(info != 0){error("c++ error: dpotri SigmaBetaInv failed\n");}
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

    /**********************************************************************
     * Set up spatial stuff and MH stuff
     * *******************************************************************/
    int nTheta, sigmaSqIndx, phiIndx, nuIndx;
    if (corName != "matern") {
      nTheta = 2; // sigma^2, phi 
      sigmaSqIndx = 0; phiIndx = 1; 
    } else {
      nTheta = 3; // sigma^2, phi, nu 
      sigmaSqIndx = 0; phiIndx = 1; nuIndx = 2; 
    }  
    double *accept = (double *) R_alloc(nTheta, sizeof(double)); zeros(accept, nTheta); 
    double *theta = (double *) R_alloc(nTheta, sizeof(double));
    double logPostThetaCurr = 0.0, logPostThetaCand = 0.0;
    double logDet;  
    double phiCand = 0.0, nuCand = 0.0, sigmaSqCand = 0.0;  
    SEXP acceptSamples_r; 
    PROTECT(acceptSamples_r = allocMatrix(REALSXP, nTheta, nBatch)); nProtect++; 
    zeros(REAL(acceptSamples_r), nTheta * nBatch);
    SEXP tuningSamples_r; 
    PROTECT(tuningSamples_r = allocMatrix(REALSXP, nTheta, nBatch)); nProtect++; 
    zeros(REAL(tuningSamples_r), nTheta * nBatch);
    SEXP thetaSamples_r; 
    PROTECT(thetaSamples_r = allocMatrix(REALSXP, nTheta, nPost)); nProtect++; 
    zeros(REAL(thetaSamples_r), nTheta * nPost);
    double a, v, b, e, muNNGP, var, aij; 
    theta[sigmaSqIndx] = sigmaSq;
    theta[phiIndx] = phi;
    if (corName == "matern") {
      theta[nuIndx] = nu;
    }
    // Allocate for the U index vector that keep track of which locations have 
    // the i-th location as a neighbor
    int nIndx = static_cast<int>(static_cast<double>(1+m)/2*m+(J-m-1)*m);

    // For NNGP
    int mm = m*m;
    double *B = (double *) R_alloc(nIndx, sizeof(double));
    double *F = (double *) R_alloc(J, sizeof(double));
    double *BCand = (double *) R_alloc(nIndx, sizeof(double));
    double *FCand = (double *) R_alloc(J, sizeof(double));
    double *c =(double *) R_alloc(m*nThreads, sizeof(double));
    double *C = (double *) R_alloc(mm*nThreads, sizeof(double));
    int sizeBK = nThreads*(1.0+static_cast<int>(floor(nuB)));
    double *bk = (double *) R_alloc(sizeBK, sizeof(double));

    // Initiate B and F for each SVC
    updateBFGaus(B, F, c, C, coords, nnIndx, nnIndxLU, J, m, 
                 theta[sigmaSqIndx], theta[phiIndx], nu, covModel, bk, nuB);

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
          tmp_J1[j] = (y[j] - w[j] - betaStarSites[j]) / tauSq;
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
        if(info != 0){error("c++ error: dpotrf here failed\n");}
        F77_NAME(dpotri)(lower, &p, tmp_pp, &p, &info FCONE); 
        if(info != 0){error("c++ error: dpotri here failed\n");}
        F77_NAME(dsymv)(lower, &p, &one, tmp_pp, &p, tmp_p, &inc, &zero, tmp_p2, &inc FCONE);
        F77_NAME(dpotrf)(lower, &p, tmp_pp, &p, &info FCONE); 
	if(info != 0){error("c++ error: dpotrf here failed\n");}
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
          		    tmp_02 + (betaStar[l] * XRandom[betaStarIndx[l] * J + j]) - 
			    w[j]) / tauSq;
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
        for(j = 0; j < J; j++){
	  tmp_J1[j] = y[j] - w[j] - 
	      	F77_NAME(ddot)(&p, &X[j], &J, beta, &inc) - betaStarSites[j];
        }
        tauSq = rigamma(tauSqA + J / 2.0, tauSqB + 0.5 * 
			   F77_NAME(ddot)(&J, tmp_J1, &inc, tmp_J1, &inc));

        /********************************************************************
         *Update w (spatial random effects)
         *******************************************************************/
	for (i = 0; i < J; i++ ) {
          a = 0;
	  v = 0;
	  if (uIndxLU[J + i] > 0){ // is i a neighbor for anybody
	    for (j = 0; j < uIndxLU[J+i]; j++){ // how many locations have i as a neighbor
	      b = 0;
	      // now the neighbors for the jth location who has i as a neighbor
	      jj = uIndx[uIndxLU[i]+j]; // jj is the index of the jth location who has i as a neighbor
	      for(k = 0; k < nnIndxLU[J+jj]; k++){ // these are the neighbors of the jjth location
	        kk = nnIndx[nnIndxLU[jj]+k]; // kk is the index for the jth locations neighbors
	        if(kk != i){ //if the neighbor of jj is not i
	  	b += B[nnIndxLU[jj]+k]*w[kk]; //covariance between jj and kk and the random effect of kk
	        }
	      }
	      aij = w[jj] - b;
	      a += B[nnIndxLU[jj]+uiIndx[uIndxLU[i]+j]]*aij/F[jj];
	      v += pow(B[nnIndxLU[jj]+uiIndx[uIndxLU[i]+j]],2)/F[jj];
	    }
	  }
	  
	  e = 0;
	  for(j = 0; j < nnIndxLU[J+i]; j++){
	    e += B[nnIndxLU[i]+j]*w[nnIndx[nnIndxLU[i]+j]];
	  }
	  
	  muNNGP = (y[i] - F77_NAME(ddot)(&p, &X[i], &J, beta, &inc) - betaStarSites[i]) / tauSq + e/F[i] + a;
	  
	  var = 1.0/(1.0/tauSq + 1.0/F[i] + v);
	  
	  w[i] = rnorm(muNNGP * var, sqrt(var));

        } // i 


        /********************************************************************
         *Update sigmaSq
         *******************************************************************/
	if (sigmaSqIG == 1) {
	  a = 0;
	  logDet = 0;
#ifdef _OPENMP
#pragma omp parallel for private (e, i, b) reduction(+:a, logDet)
#endif
          for (j = 0; j < J; j++){
            if(nnIndxLU[J+j] > 0){
              e = 0;
              for(i = 0; i < nnIndxLU[J+j]; i++){
                e += B[nnIndxLU[j]+i]*w[nnIndx[nnIndxLU[j]+i]];
              }
              b = w[j] - e;
            }else{
              b = w[j];
            }
            a += b*b/F[j];
          }

	  theta[sigmaSqIndx] = rigamma(sigmaSqA + J / 2.0, 
			               sigmaSqB + 0.5 * a * theta[sigmaSqIndx]);
	}

        /********************************************************************
         *Update phi (and nu if matern)
         *******************************************************************/
        if (corName == "matern"){ nu = theta[nuIndx]; }
        updateBFGaus(B, F, c, C, coords, nnIndx, nnIndxLU, J, m, theta[sigmaSqIndx], 
		   theta[phiIndx], nu, covModel, bk, nuB);
        
        a = 0;
        logDet = 0;

#ifdef _OPENMP
#pragma omp parallel for private (e, i, b) reduction(+:a, logDet)
#endif
        for (j = 0; j < J; j++){
          if (nnIndxLU[J+j] > 0){
            e = 0;
            for (i = 0; i < nnIndxLU[J+j]; i++){
              e += B[nnIndxLU[j]+i]*w[nnIndx[nnIndxLU[j]+i]];
            }
            b = w[j] - e;
          } else{
            b = w[j];
          }	
          a += b*b/F[j];
          logDet += log(F[j]);
        }
      
        logPostThetaCurr = -0.5*logDet - 0.5*a;
        logPostThetaCurr += log(theta[phiIndx] - phiA) + log(phiB - theta[phiIndx]); 
        if(corName == "matern"){
        	logPostThetaCurr += log(theta[nuIndx] - nuA) + log(nuB - theta[nuIndx]); 
        }
	if (sigmaSqIG == 0) {
          logPostThetaCurr += log(theta[sigmaSqIndx] - sigmaSqA) + log(sigmaSqB - theta[sigmaSqIndx]);
	}
        
        // Candidate
        phiCand = logitInv(rnorm(logit(theta[phiIndx], phiA, phiB), 
				exp(tuning[phiIndx])), phiA, phiB);
	if (sigmaSqIG == 0) {
	  sigmaSqCand = logitInv(rnorm(logit(theta[sigmaSqIndx], sigmaSqA, sigmaSqB), 
				 exp(tuning[sigmaSqIndx])), sigmaSqA, sigmaSqB); 
	}
        if (corName == "matern"){
      	  nuCand = logitInv(rnorm(logit(theta[nuIndx], nuA, nuB), exp(tuning[nuIndx])), nuA, nuB);
        }

        if (sigmaSqIG) {
          updateBFGaus(BCand, FCand, c, C, coords, nnIndx, nnIndxLU, J, m, theta[sigmaSqIndx], phiCand, nuCand, covModel, bk, nuB);
	  } else {
            updateBFGaus(BCand, FCand, c, C, coords, nnIndx, nnIndxLU, J, m, sigmaSqCand, phiCand, nuCand, covModel, bk, nuB);
	}
      
        a = 0;
        logDet = 0;
      
#ifdef _OPENMP
#pragma omp parallel for private (e, i, b) reduction(+:a, logDet)
#endif
        for (j = 0; j < J; j++){
          if (nnIndxLU[J+j] > 0){
            e = 0;
            for (i = 0; i < nnIndxLU[J+j]; i++){
              e += BCand[nnIndxLU[j]+i]*w[nnIndx[nnIndxLU[j]+i]];
            }
            b = w[j] - e;
          } else{
            b = w[j];
            }	
            a += b*b/FCand[j];
            logDet += log(FCand[j]);
        }
        
        logPostThetaCand = -0.5*logDet - 0.5*a;      
        logPostThetaCand += log(phiCand - phiA) + log(phiB - phiCand); 
        if (corName == "matern"){
          logPostThetaCand += log(nuCand - nuA) + log(nuB - nuCand); 
        }
	  if (sigmaSqIG == 0) {
            logPostThetaCand += log(sigmaSqCand - sigmaSqA) + log(sigmaSqB - sigmaSqCand);
	  }

        if (runif(0.0,1.0) <= exp(logPostThetaCand - logPostThetaCurr)) {

          std::swap(BCand, B);
          std::swap(FCand, F);
          
          theta[phiIndx] = phiCand;
          accept[phiIndx]++;
          if(corName == "matern"){
            theta[nuIndx] = nuCand; 
            accept[nuIndx]++; 
          }
	  if (sigmaSqIG == 0) {
            theta[sigmaSqIndx] = sigmaSqCand;
	    accept[sigmaSqIndx]++;
	  }
        }

        /********************************************************************
         *Get fitted values and likelihood for WAIC
         *******************************************************************/
	if (saveFitted == 1) {
          for (j = 0; j < J; j++) {
            mu[j] = F77_NAME(ddot)(&p, &X[j], &J, beta, &inc) + w[j] + betaStarSites[j];
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
            F77_NAME(dcopy)(&J, w, &inc, &REAL(wSamples_r)[sPost*J], &inc); 
	    REAL(tauSqSamples_r)[sPost] = tauSq;
	    F77_NAME(dcopy)(&nTheta, theta, &inc, 
			    &REAL(thetaSamples_r)[sPost*nTheta], &inc); 
	    if (saveFitted == 1) {
	      F77_NAME(dcopy)(&J, yRep, &inc, &REAL(yRepSamples_r)[sPost*J], &inc); 
	      F77_NAME(dcopy)(&JZero, yRepZero, &inc, &REAL(yRepZeroSamples_r)[sPost*JZero], &inc); 
              F77_NAME(dcopy)(&J, mu, &inc, &REAL(muSamples_r)[sPost*J], &inc); 
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
       *Adjust tuning 
       *******************************************************************/
      for (j = 0; j < nTheta; j++) {
        REAL(acceptSamples_r)[s * nTheta + j] = accept[j]/batchLength; 
        REAL(tuningSamples_r)[s * nTheta + j] = tuning[j]; 
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
	  Rprintf("\tphi\t\t%3.1f\t\t%1.5f\n", 100.0*REAL(acceptSamples_r)[s * nTheta + phiIndx], exp(tuning[phiIndx]));
	  if (corName == "matern") {
	    Rprintf("\tnu\t\t%3.1f\t\t%1.5f\n", 100.0*REAL(acceptSamples_r)[s * nTheta + nuIndx], exp(tuning[nuIndx]));
	  }
	  if (sigmaSqIG == 0) {
	    Rprintf("\tsigmaSq\t\t%3.1f\t\t%1.5f\n", 100.0*REAL(acceptSamples_r)[s * nTheta + sigmaSqIndx], exp(tuning[sigmaSqIndx]));
	  }
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
    int nResultListObjs = 10;
    if (pRE > 0) {
      nResultListObjs += 2;
    }

    PROTECT(result_r = allocVector(VECSXP, nResultListObjs)); nProtect++;
    PROTECT(resultName_r = allocVector(VECSXP, nResultListObjs)); nProtect++;

    // Setting the components of the output list.
    SET_VECTOR_ELT(result_r, 0, betaSamples_r);
    SET_VECTOR_ELT(result_r, 1, tauSqSamples_r);
    if (saveFitted == 1) {
      SET_VECTOR_ELT(result_r, 2, yRepSamples_r); 
      SET_VECTOR_ELT(result_r, 3, muSamples_r);
      SET_VECTOR_ELT(result_r, 8, likeSamples_r); 
      SET_VECTOR_ELT(result_r, 9, yRepZeroSamples_r); 
    }
    SET_VECTOR_ELT(result_r, 4, thetaSamples_r); 
    SET_VECTOR_ELT(result_r, 5, wSamples_r); 
    SET_VECTOR_ELT(result_r, 6, tuningSamples_r); 
    SET_VECTOR_ELT(result_r, 7, acceptSamples_r); 
    if (pRE > 0) {
      SET_VECTOR_ELT(result_r, 10, sigmaSqMuSamples_r);
      SET_VECTOR_ELT(result_r, 11, betaStarSamples_r);
    }

    // mkChar turns a C string into a CHARSXP
    SET_VECTOR_ELT(resultName_r, 0, mkChar("beta.samples")); 
    SET_VECTOR_ELT(resultName_r, 1, mkChar("tau.sq.samples")); 
    if (saveFitted == 1) {
      SET_VECTOR_ELT(resultName_r, 2, mkChar("y.rep.samples")); 
      SET_VECTOR_ELT(resultName_r, 3, mkChar("mu.samples"));
      SET_VECTOR_ELT(resultName_r, 8, mkChar("like.samples")); 
      SET_VECTOR_ELT(resultName_r, 9, mkChar("y.rep.zero.samples")); 
    }
    SET_VECTOR_ELT(resultName_r, 4, mkChar("theta.samples")); 
    SET_VECTOR_ELT(resultName_r, 5, mkChar("w.samples")); 
    SET_VECTOR_ELT(resultName_r, 6, mkChar("tune")); 
    SET_VECTOR_ELT(resultName_r, 7, mkChar("accept")); 
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


