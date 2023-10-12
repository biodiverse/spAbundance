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
  SEXP msAbundGaussian(SEXP y_r, SEXP X_r, SEXP XRE_r,
                       SEXP XRandom_r, SEXP consts_r, SEXP nRELong_r, 
                       SEXP betaStarting_r, SEXP betaCommStarting_r,  
                       SEXP tauSqBetaStarting_r, SEXP tauSqStarting_r,
                       SEXP sigmaSqMuStarting_r, SEXP betaStarStarting_r, 
		       SEXP betaStarIndx_r, SEXP betaLevelIndx_r, 
                       SEXP muBetaComm_r, SEXP SigmaBetaComm_r, SEXP tauSqBetaA_r, 
                       SEXP tauSqBetaB_r, SEXP tauSqA_r, SEXP tauSqB_r,
		       SEXP sigmaSqMuA_r, SEXP sigmaSqMuB_r, 
                       SEXP tuning_r, SEXP nBatch_r, SEXP batchLength_r, 
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
    int i, j, s, g, t, h, l, info, nProtect=0, ll;    

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
    int *XRE = INTEGER(XRE_r);
    // Load constants
    int N = INTEGER(consts_r)[0]; 
    int J = INTEGER(consts_r)[1];
    int p = INTEGER(consts_r)[2];
    int pRE = INTEGER(consts_r)[3];
    int nRE = INTEGER(consts_r)[4];
    int saveFitted = INTEGER(consts_r)[5];
    int pp = p * p; 
    double *muBetaComm = REAL(muBetaComm_r); 
    double *SigmaBetaCommInv = (double *) R_alloc(pp, sizeof(double));   
    F77_NAME(dcopy)(&pp, REAL(SigmaBetaComm_r), &inc, SigmaBetaCommInv, &inc);
    double *tauSqBetaA = REAL(tauSqBetaA_r); 
    double *tauSqBetaB = REAL(tauSqBetaB_r); 
    double *tauSqA = REAL(tauSqA_r); 
    double *tauSqB = REAL(tauSqB_r); 
    double *sigmaSqMuA = REAL(sigmaSqMuA_r); 
    double *sigmaSqMuB = REAL(sigmaSqMuB_r); 
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
          Rprintf("Multi-species Gaussian Model with %i sites and %i species.\n\n", J, N);
	} else {
          Rprintf("Multi-species Zero-Inflated Gaussian Model with %i sites and %i species.\n\n", J, N);
	}
        Rprintf("Samples per chain: %i (%i batches of length %i)\n", nSamples, nBatch, batchLength);
        Rprintf("Burn-in: %i \n", nBurn); 
        Rprintf("Thinning Rate: %i \n", nThin); 
        Rprintf("Number of Chains: %i \n", nChain);
        Rprintf("Total Posterior Samples: %i \n\n", nPost * nChain); 
#ifdef _OPENMP
        Rprintf("Source compiled with OpenMP support and model fit using %i thread(s).\n\n", nThreads);
#else
        Rprintf("Source not compiled with OpenMP support.\n\n");
#endif
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
    int JN = J * N;
    int Jp = J * p; 
    int JpRE = J * p;
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
    double *tmp_N = (double *) R_alloc(N, sizeof(double)); zeros(tmp_N, N);
    double *tmp_N2 = (double *) R_alloc(N, sizeof(double)); zeros(tmp_N2, N);

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
    // Latent occurrence random effects
    double *betaStar = (double *) R_alloc(nREN, sizeof(double)); 
    F77_NAME(dcopy)(&nREN, REAL(betaStarStarting_r), &inc, betaStar, &inc); 
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
    // Random effects
    SEXP sigmaSqMuSamples_r; 
    SEXP betaStarSamples_r; 
    if (pRE > 0) {
      PROTECT(sigmaSqMuSamples_r = allocMatrix(REALSXP, pRE, nPost)); nProtect++;
      zeros(REAL(sigmaSqMuSamples_r), pRE * nPost);
      PROTECT(betaStarSamples_r = allocMatrix(REALSXP, nREN, nPost)); nProtect++;
      zeros(REAL(betaStarSamples_r), nREN * nPost);
    }
    // Likelihood samples for WAIC. 
    
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
    // For current number of nonzero z values
    double *currJ = (double *) R_alloc(N, sizeof(double)); zeros(currJ, N);
    for (i = 0; i < N; i++) {
      for (j = 0; j < J; j++) {
        currJ[i] += z[j * N + i];
      }
    }

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

    GetRNGstate();

    /**********************************************************************
     Start sampling
     * *******************************************************************/
    for (s = 0, g = 0; s < nBatch; s++) {
      for (t = 0; t < batchLength; t++, g++) {

        /********************************************************************
         Update Community level Coefficients
         *******************************************************************/
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
        /********************************************************************
         Update Community Occupancy Variance Parameter
        ********************************************************************/
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
	      tmp_J1[j] = (y[j * N + i] - betaStarSites[i * J + j]) / tauSq[i];
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
	       tmp_J1[j] = y[j * N + i] - 
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
                  tmp_one[0] += XRandom[betaStarIndx[l] * J + j] * (y[j * N + i] - F77_NAME(ddot)(&p, &X[j], &J, &beta[i], &N) + (betaStar[i * nRE + l] * XRandom[betaStarIndx[l] * J + j]) - tmp_02) / tauSq[i];
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
         *Get fitted values and likelihood for WAIC
         *******************************************************************/
	if (saveFitted == 1) {
	  for (i = 0; i < N; i++) {
            for (j = 0; j < J; j++) {
              if (z[j * N + i] == 1.0) {
                mu[j * N + i] = F77_NAME(ddot)(&p, &X[j], &J, &beta[i], &N) + betaStarSites[i * J + j];
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
	    if (saveFitted == 1) {
              F77_NAME(dcopy)(&JN, mu, &inc, &REAL(muSamples_r)[sPost*JN], &inc); 
              F77_NAME(dcopy)(&JN, yRep, &inc, &REAL(yRepSamples_r)[sPost*JN], &inc); 
              F77_NAME(dcopy)(&JN, like, &inc, &REAL(likeSamples_r)[sPost*JN], &inc); 
	    }
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
       *Report 
       *******************************************************************/
      if (verbose) {
	if (status == nReport) {
	  Rprintf("Batch: %i of %i, %3.2f%%\n", s, nBatch, 100.0*s/nBatch);
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
    int nResultListObjs = 7;
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
      SET_VECTOR_ELT(result_r, 6, likeSamples_r); 
    }
    SET_VECTOR_ELT(result_r, 5, tauSqSamples_r); 
    if (pRE > 0) {
      SET_VECTOR_ELT(result_r, 7, sigmaSqMuSamples_r);
      SET_VECTOR_ELT(result_r, 8, betaStarSamples_r);
    }

    // mkChar turns a C string into a CHARSXP
    SET_VECTOR_ELT(resultName_r, 0, mkChar("beta.comm.samples")); 
    SET_VECTOR_ELT(resultName_r, 1, mkChar("tau.sq.beta.samples")); 
    SET_VECTOR_ELT(resultName_r, 2, mkChar("beta.samples")); 
    if (saveFitted == 1) {
      SET_VECTOR_ELT(resultName_r, 3, mkChar("y.rep.samples")); 
      SET_VECTOR_ELT(resultName_r, 4, mkChar("mu.samples")); 
      SET_VECTOR_ELT(resultName_r, 6, mkChar("like.samples")); 
    }
    SET_VECTOR_ELT(resultName_r, 5, mkChar("tau.sq.samples")); 
    if (pRE > 0) {
      SET_VECTOR_ELT(resultName_r, 7, mkChar("sigma.sq.mu.samples")); 
      SET_VECTOR_ELT(resultName_r, 8, mkChar("beta.star.samples")); 
    }
   
    // Set the names of the output list.  
    namesgets(result_r, resultName_r);
    
    //unprotect
    UNPROTECT(nProtect);
    
    return(result_r);
  }
}
