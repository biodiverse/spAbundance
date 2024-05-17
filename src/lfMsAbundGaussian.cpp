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
  SEXP lfMsAbundGaussian(SEXP y_r, SEXP X_r, SEXP XRE_r,
                       SEXP XRandom_r, SEXP consts_r, SEXP nRELong_r,
                       SEXP betaStarting_r, SEXP betaCommStarting_r,
                       SEXP tauSqBetaStarting_r, SEXP tauSqStarting_r,
                       SEXP lambdaStarting_r, SEXP wStarting_r,
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
    int *XRE = INTEGER(XRE_r);
    // Load constants
    int N = INTEGER(consts_r)[0];
    int J = INTEGER(consts_r)[1];
    int p = INTEGER(consts_r)[2];
    int pRE = INTEGER(consts_r)[3];
    int nRE = INTEGER(consts_r)[4];
    int q = INTEGER(consts_r)[5];
    int saveFitted = INTEGER(consts_r)[6];
    int indBetas = INTEGER(consts_r)[7];
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
      Rf_warning("n.omp.threads > %i, but source not compiled with OpenMP support.", nThreads);
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
          Rprintf("Latent Factor Multi-species Gaussian Model\nwith %i sites and %i species.\n\n", J, N);
	} else {
          Rprintf("Latent Factor Multi-species Zero-Inflated Gaussian Model\nwith %i sites and %i species.\n\n", J, N);
	}
        Rprintf("Samples per chain: %i (%i batches of length %i)\n", nSamples, nBatch, batchLength);
        Rprintf("Burn-in: %i \n", nBurn);
        Rprintf("Thinning Rate: %i \n", nThin);
        Rprintf("Number of Chains: %i \n", nChain);
        Rprintf("Total Posterior Samples: %i \n\n", nPost * nChain);
        Rprintf("Using %i latent factors.\n\n", q);
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
    int Jq = J * q;
    int qq = q * q;
    int Nq = N * q;
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
    double *tmp_qq = (double *) R_alloc(qq, sizeof(double)); zeros(tmp_qq, qq);
    double *tmp_q = (double *) R_alloc(q, sizeof(double)); zeros(tmp_q, q);
    double *tmp_q2 = (double *) R_alloc(q, sizeof(double)); zeros(tmp_q2, q);
    double *tmp_qq2 = (double *) R_alloc(qq, sizeof(double)); zeros(tmp_qq2, qq);
    double *tmp_Jq = (double *) R_alloc(Jq, sizeof(double)); zeros(tmp_Jq, Jq);
    double *tmp_Nq = (double *) R_alloc(Nq, sizeof(double)); zeros(tmp_Nq, Nq);
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
    // Latent factor stuff
    double *w = (double *) R_alloc(Jq, sizeof(double));
    F77_NAME(dcopy)(&Jq, REAL(wStarting_r), &inc, w, &inc);
    // Latent factor loadings
    double *lambda = (double *) R_alloc(Nq, sizeof(double));
    F77_NAME(dcopy)(&Nq, REAL(lambdaStarting_r), &inc, lambda, &inc);
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
    PROTECT(betaCommSamples_r = Rf_allocMatrix(REALSXP, p, nPost)); nProtect++;
    zeros(REAL(betaCommSamples_r), p * nPost);
    SEXP tauSqBetaSamples_r;
    PROTECT(tauSqBetaSamples_r = Rf_allocMatrix(REALSXP, p, nPost)); nProtect++;
    zeros(REAL(tauSqBetaSamples_r), p * nPost);
    // Species level
    SEXP betaSamples_r;
    PROTECT(betaSamples_r = Rf_allocMatrix(REALSXP, pN, nPost)); nProtect++;
    zeros(REAL(betaSamples_r), pN * nPost);
    SEXP tauSqSamples_r;
    PROTECT(tauSqSamples_r = Rf_allocMatrix(REALSXP, N, nPost)); nProtect++;
    zeros(REAL(tauSqSamples_r), N * nPost);
    SEXP yRepSamples_r;
    SEXP muSamples_r;
    SEXP likeSamples_r;
    if (saveFitted == 1) {
      PROTECT(yRepSamples_r = Rf_allocMatrix(REALSXP, JN, nPost)); nProtect++;
      zeros(REAL(yRepSamples_r), JN * nPost);
      PROTECT(muSamples_r = Rf_allocMatrix(REALSXP, JN, nPost)); nProtect++;
      zeros(REAL(muSamples_r), JN * nPost);
      PROTECT(likeSamples_r = Rf_allocMatrix(REALSXP, JN, nPost)); nProtect++;
      zeros(REAL(likeSamples_r), JN * nPost);
    }
    SEXP lambdaSamples_r;
    PROTECT(lambdaSamples_r = Rf_allocMatrix(REALSXP, Nq, nPost)); nProtect++;
    zeros(REAL(lambdaSamples_r), Nq * nPost);
    SEXP wSamples_r;
    PROTECT(wSamples_r = Rf_allocMatrix(REALSXP, Jq, nPost)); nProtect++;
    zeros(REAL(wSamples_r), Jq * nPost);
    // Random effects
    SEXP sigmaSqMuSamples_r;
    SEXP betaStarSamples_r;
    if (pRE > 0) {
      PROTECT(sigmaSqMuSamples_r = Rf_allocMatrix(REALSXP, pRE, nPost)); nProtect++;
      zeros(REAL(sigmaSqMuSamples_r), pRE * nPost);
      PROTECT(betaStarSamples_r = Rf_allocMatrix(REALSXP, nREN, nPost)); nProtect++;
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
    if(info != 0){Rf_error("c++ error: dpotrf SigmaBetaCommInv failed\n");}
    F77_NAME(dpotri)(lower, &p, SigmaBetaCommInv, &p, &info FCONE);
    if(info != 0){Rf_error("c++ error: dpotri SigmaBetaCommInv failed\n");}
    double *SigmaBetaCommInvMuBeta = (double *) R_alloc(p, sizeof(double));
    F77_NAME(dsymv)(lower, &p, &one, SigmaBetaCommInv, &p, muBetaComm, &inc, &zero,
        	    SigmaBetaCommInvMuBeta, &inc FCONE);
    // Put community level occurrence variances in a p x p matrix.
    double *TauBetaInv = (double *) R_alloc(pp, sizeof(double)); zeros(TauBetaInv, pp);
    for (i = 0; i < p; i++) {
      TauBetaInv[i * p + i] = tauSqBeta[i];
    } // i
    F77_NAME(dpotrf)(lower, &p, TauBetaInv, &p, &info FCONE);
    if(info != 0){Rf_error("c++ error: dpotrf TauBetaInv failed\n");}
    F77_NAME(dpotri)(lower, &p, TauBetaInv, &p, &info FCONE);
    if(info != 0){Rf_error("c++ error: dpotri TauBetaInv failed\n");}
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

    /**********************************************************************
     Latent factor stuff
     * *******************************************************************/
    // Species-level latent factor effects
    double *wStar = (double *) R_alloc(JN, sizeof(double)); zeros(wStar, JN);
    // Multiply Lambda %*% w[j] to get wStar.
    for (j = 0; j < J; j++) {
      F77_NAME(dgemv)(ntran, &N, &q, &one, lambda, &N, &w[j*q],
                      &inc, &zero, &wStar[j * N], &inc FCONE);
    }
    double *var = (double *) R_alloc(qq, sizeof(double)); zeros(var, qq);

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
          if(info != 0){Rf_error("c++ error: dpotrf ABetaComm failed\n");}
          F77_NAME(dpotri)(lower, &p, tmp_pp, &p, &info FCONE);
          if(info != 0){Rf_error("c++ error: dpotri ABetaComm failed\n");}
          // A.beta.inv %*% b.beta
          // 1 * tmp_pp * tmp_p + 0 * tmp_p2  = tmp_p2
          F77_NAME(dsymv)(lower, &p, &one, tmp_pp, &p, tmp_p, &inc, &zero, tmp_p2, &inc FCONE);
          // Computes cholesky of tmp_pp again stored back in tmp_pp. This chol(A.beta.inv)
          F77_NAME(dpotrf)(lower, &p, tmp_pp, &p, &info FCONE);
          if(info != 0){Rf_error("c++ error: dpotrf ABetaComm failed\n");}
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
          if(info != 0){Rf_error("c++ error: dpotrf TauBetaInv failed\n");}
          F77_NAME(dpotri)(lower, &p, TauBetaInv, &p, &info FCONE);
          if(info != 0){Rf_error("c++ error: dpotri TauBetaInv failed\n");}
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
	      tmp_J1[j] = (y[j * N + i] - betaStarSites[i * J + j] - wStar[j * N + i]) / tauSq[i];
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
          if(info != 0){Rf_error("c++ error: dpotrf ABeta failed\n");}
          F77_NAME(dpotri)(lower, &p, tmp_pp, &p, &info FCONE);
          if(info != 0){Rf_error("c++ error: dpotri ABeta failed\n");}
          // A.beta.inv %*% b.beta
          F77_NAME(dsymv)(lower, &p, &one, tmp_pp, &p, tmp_p, &inc, &zero, tmp_p2, &inc FCONE);
          F77_NAME(dpotrf)(lower, &p, tmp_pp, &p, &info FCONE);
	  if(info != 0){Rf_error("c++ error: dpotrf A.beta 2 failed\n");}
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
         *Update latent effects (w)
         *******************************************************************/
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
          F77_NAME(dgemm)(ytran, ntran, &q, &q, &N, &one, tmp_Nq, &N, lambda, &N,
			  &zero, var, &q FCONE FCONE);

          // var
          for (k = 0; k < q; k++) {
            var[k * q + k] += 1.0;
          } // k
          F77_NAME(dpotrf)(lower, &q, var, &q, &info FCONE);
          if (info != 0){Rf_error("c++ error: dpotrf var failed\n");}
          F77_NAME(dpotri)(lower, &q, var, &q, &info FCONE);
          if (info != 0){Rf_error("c++ error: dpotri var failed\n");}

          // mu
	  zeros(tmp_N, N);
          for (k = 0; k < N; k++) {
            if (z[ii * N + k] == 1.0) {
              tmp_N[k] = (y[ii * N + k] - F77_NAME(ddot)(&p, &X[ii], &J, &beta[k], &N) - betaStarSites[k * J + ii]) / tauSq[k];
	    }
          } // k

          F77_NAME(dgemv)(ytran, &N, &q, &one, lambda, &N, tmp_N, &inc, &zero, mu, &inc FCONE);

          F77_NAME(dsymv)(lower, &q, &one, var, &q, mu, &inc, &zero, tmp_N, &inc FCONE);

          F77_NAME(dpotrf)(lower, &q, var, &q, &info FCONE);
          if(info != 0){Rf_error("c++ error: dpotrf var 2 failed\n");}

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
          if(info != 0){Rf_error("c++ error: dpotrf for spatial factors failed\n");}
          F77_NAME(dpotri)(lower, &currDim, tmp_qq2, &currDim, &info FCONE);
          if(info != 0){Rf_error("c++ error: dpotri for spatial factors failed\n");}

          F77_NAME(dsymv)(lower, &currDim, &one, tmp_qq2, &currDim, tmp_q, &inc, &zero, tmp_q2, &inc FCONE);

          F77_NAME(dpotrf)(lower, &currDim, tmp_qq2, &currDim, &info FCONE);
          if(info != 0){Rf_error("c++ error: dpotrf for spatial factors 2 failed\n");}

          mvrnorm(tmp_q, tmp_q2, tmp_qq2, currDim);
          F77_NAME(dcopy)(&currDim, tmp_q, &inc, &lambda[i], &N);
        } // i

        // Multiply Lambda %*% w[j] to get wStar.
        for (j = 0; j < J; j++) {
          F77_NAME(dgemv)(ntran, &N, &q, &one, lambda, &N, &w[j*q], &inc, &zero, &wStar[j * N], &inc FCONE);
        } // j

        /********************************************************************
         *Get fitted values and likelihood for WAIC
         *******************************************************************/
        if (saveFitted == 1) {
	  for (i = 0; i < N; i++) {
            for (j = 0; j < J; j++) {
              if (z[j * N + i] == 1.0) {
                mu[j * N + i] = F77_NAME(ddot)(&p, &X[j], &J, &beta[i], &N) + betaStarSites[i * J + j] +
	  	              wStar[j * N + i];
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
            F77_NAME(dcopy)(&Jq, w, &inc, &REAL(wSamples_r)[sPost*Jq], &inc);
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
    int nResultListObjs = 9;
    if (pRE > 0) {
      nResultListObjs += 2;
    }

    PROTECT(result_r = Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;
    PROTECT(resultName_r = Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;

    // Setting the components of the output list.
    SET_VECTOR_ELT(result_r, 0, betaCommSamples_r);
    SET_VECTOR_ELT(result_r, 1, tauSqBetaSamples_r);
    SET_VECTOR_ELT(result_r, 2, betaSamples_r);
    if (saveFitted == 1) {
      SET_VECTOR_ELT(result_r, 3, yRepSamples_r);
      SET_VECTOR_ELT(result_r, 6, likeSamples_r);
      SET_VECTOR_ELT(result_r, 4, muSamples_r);
    }
    SET_VECTOR_ELT(result_r, 5, tauSqSamples_r);
    SET_VECTOR_ELT(result_r, 7, lambdaSamples_r);
    SET_VECTOR_ELT(result_r, 8, wSamples_r);
    if (pRE > 0) {
      SET_VECTOR_ELT(result_r, 9, sigmaSqMuSamples_r);
      SET_VECTOR_ELT(result_r, 10, betaStarSamples_r);
    }

    // Rf_mkChar turns a C string into a CHARSXP
    SET_VECTOR_ELT(resultName_r, 0, Rf_mkChar("beta.comm.samples"));
    SET_VECTOR_ELT(resultName_r, 1, Rf_mkChar("tau.sq.beta.samples"));
    SET_VECTOR_ELT(resultName_r, 2, Rf_mkChar("beta.samples"));
    if (saveFitted == 1) {
      SET_VECTOR_ELT(resultName_r, 3, Rf_mkChar("y.rep.samples"));
      SET_VECTOR_ELT(resultName_r, 4, Rf_mkChar("mu.samples"));
      SET_VECTOR_ELT(resultName_r, 6, Rf_mkChar("like.samples"));
    }
    SET_VECTOR_ELT(resultName_r, 5, Rf_mkChar("tau.sq.samples"));
    SET_VECTOR_ELT(resultName_r, 7, Rf_mkChar("lambda.samples"));
    SET_VECTOR_ELT(resultName_r, 8, Rf_mkChar("w.samples"));
    if (pRE > 0) {
      SET_VECTOR_ELT(resultName_r, 9, Rf_mkChar("sigma.sq.mu.samples"));
      SET_VECTOR_ELT(resultName_r, 10, Rf_mkChar("beta.star.samples"));
    }

    // Set the names of the output list.
    Rf_namesgets(result_r, resultName_r);

    //unprotect
    UNPROTECT(nProtect);

    return(result_r);
  }
}


