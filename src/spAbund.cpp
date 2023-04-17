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
  SEXP spAbund(SEXP y_r, SEXP X_r, SEXP coordsD_r, SEXP XRE_r, 
               SEXP XRandom_r, SEXP consts_r, SEXP nAbundRELong_r, 
               SEXP m_r, SEXP nnIndx_r, 
               SEXP nnIndxLU_r, SEXP uIndx_r, SEXP uIndxLU_r, SEXP uiIndx_r,
               SEXP betaStarting_r, SEXP kappaStarting_r,
               SEXP sigmaSqMuStarting_r, SEXP betaStarStarting_r, 
               SEXP wStarting_r, SEXP phiStarting_r, SEXP sigmaSqStarting_r, 
               SEXP nuStarting_r,
               SEXP siteIndx_r, SEXP betaStarIndx_r, SEXP betaLevelIndx_r, 
               SEXP muBeta_r, SEXP SigmaBeta_r, 
	       SEXP kappaA_r, SEXP kappaB_r,
	       SEXP sigmaSqMuA_r, SEXP sigmaSqMuB_r, 
               SEXP phiA_r, SEXP phiB_r, 
               SEXP sigmaSqA_r, SEXP sigmaSqB_r, SEXP nuA_r, SEXP nuB_r, SEXP tuning_r,
               SEXP covModel_r,
               SEXP nBatch_r, SEXP batchLength_r, SEXP acceptRate_r, SEXP nThreads_r, 
               SEXP verbose_r, SEXP nReport_r, SEXP samplesInfo_r,
               SEXP chainInfo_r, SEXP sigmaSqIG_r, SEXP family_r){
   
    /**********************************************************************
     * Initial constants
     * *******************************************************************/
    int i, ii, g, t, j, kk, k, jj, s, r, l, ll, info, nProtect=0;
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
    // Load constants
    int J = INTEGER(consts_r)[0];
    int nObs = INTEGER(consts_r)[1];
    int pAbund = INTEGER(consts_r)[2];
    int pAbundRE = INTEGER(consts_r)[3];
    int nAbundRE = INTEGER(consts_r)[4];
    int ppAbund = pAbund * pAbund; 
    double *muBeta = (double *) R_alloc(pAbund, sizeof(double));   
    F77_NAME(dcopy)(&pAbund, REAL(muBeta_r), &inc, muBeta, &inc);
    double *SigmaBeta = (double *) R_alloc(ppAbund, sizeof(double));
    F77_NAME(dcopy)(&ppAbund, REAL(SigmaBeta_r), &inc, SigmaBeta, &inc);
    double kappaA = REAL(kappaA_r)[0];
    double kappaB = REAL(kappaB_r)[0];
    double *sigmaSqMuA = REAL(sigmaSqMuA_r); 
    double *sigmaSqMuB = REAL(sigmaSqMuB_r); 
    double phiA = REAL(phiA_r)[0];
    double phiB = REAL(phiB_r)[0]; 
    double nuA = REAL(nuA_r)[0]; 
    double nuB = REAL(nuB_r)[0]; 
    double sigmaSqA = REAL(sigmaSqA_r)[0]; 
    double sigmaSqB = REAL(sigmaSqB_r)[0]; 
    int *nAbundRELong = INTEGER(nAbundRELong_r); 
    int *siteIndx = INTEGER(siteIndx_r); 
    int *betaStarIndx = INTEGER(betaStarIndx_r); 
    int *betaLevelIndx = INTEGER(betaLevelIndx_r);
    int nBatch = INTEGER(nBatch_r)[0]; 
    int batchLength = INTEGER(batchLength_r)[0]; 
    int nSamples = nBatch * batchLength; 
    int nBurn = INTEGER(samplesInfo_r)[0]; 
    int nThin = INTEGER(samplesInfo_r)[1];
    int nPost = INTEGER(samplesInfo_r)[2]; 
    int m = INTEGER(m_r)[0];
    int *nnIndx = INTEGER(nnIndx_r);
    int *nnIndxLU = INTEGER(nnIndxLU_r);
    int *uIndx = INTEGER(uIndx_r);
    int *uIndxLU = INTEGER(uIndxLU_r);
    int *uiIndx = INTEGER(uiIndx_r);
    int covModel = INTEGER(covModel_r)[0];
    std::string corName = getCorName(covModel);
    int currChain = INTEGER(chainInfo_r)[0];
    double acceptRate = REAL(acceptRate_r)[0];
    double *tuning = REAL(tuning_r);
    double *coordsD = REAL(coordsD_r);
    int nChain = INTEGER(chainInfo_r)[1];
    int nThreads = INTEGER(nThreads_r)[0];
    int verbose = INTEGER(verbose_r)[0];
    int nReport = INTEGER(nReport_r)[0];
    int sigmaSqIG = INTEGER(sigmaSqIG_r)[0];
    int status = 0; 
    int thinIndx = 0;
    int sPost = 0;  
    // NB = 1, Poisson = 0;
    int family = INTEGER(family_r)[0];

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
          Rprintf("Spatial Negative Binomial Abundance model fit with %i sites.\n\n", J);
	} else {
          Rprintf("Spatial Poisson Abundance model fit with %i sites.\n\n", J);
	}
        Rprintf("Samples per Chain: %i (%i batches of length %i)\n", nSamples, nBatch, batchLength);
        Rprintf("Burn-in: %i \n", nBurn); 
        Rprintf("Thinning Rate: %i \n", nThin); 
        Rprintf("Number of Chains: %i \n", nChain);
        Rprintf("Total Posterior Samples: %i \n\n", nPost * nChain); 
        Rprintf("Using the %s spatial correlation model.\n\n", corName.c_str());
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
    // Overdispersion parameter for NB;
    double kappa = REAL(kappaStarting_r)[0];
    double* epsilon = (double *) R_alloc(nObs, sizeof(double)); ones(epsilon, nObs);
    // Latent random effects
    double *betaStar = (double *) R_alloc(nAbundRE, sizeof(double)); 
    F77_NAME(dcopy)(&nAbundRE, REAL(betaStarStarting_r), &inc, betaStar, &inc); 
    // Spatial parameters
    double *w = (double *) R_alloc(J, sizeof(double));   
    F77_NAME(dcopy)(&J, REAL(wStarting_r), &inc, w, &inc);
    // Latent Abundance
    double *yRep = (double *) R_alloc(nObs, sizeof(double)); zeros(yRep, nObs);

    /**********************************************************************
     * Return Stuff
     * *******************************************************************/
    SEXP betaSamples_r;
    PROTECT(betaSamples_r = allocMatrix(REALSXP, pAbund, nPost)); nProtect++;
    SEXP yRepSamples_r; 
    PROTECT(yRepSamples_r = allocMatrix(REALSXP, nObs, nPost)); nProtect++; 
    SEXP muSamples_r; 
    PROTECT(muSamples_r = allocMatrix(REALSXP, nObs, nPost)); nProtect++; 
    SEXP wSamples_r;
    PROTECT(wSamples_r = allocMatrix(REALSXP, J, nPost)); nProtect++;
    SEXP kappaSamples_r;
    if (family == 1) {
      PROTECT(kappaSamples_r = allocMatrix(REALSXP, inc, nPost)); nProtect++;
    }
    // Abundance random effects
    SEXP sigmaSqMuSamples_r; 
    SEXP betaStarSamples_r; 
    if (pAbundRE > 0) {
      PROTECT(sigmaSqMuSamples_r = allocMatrix(REALSXP, pAbundRE, nPost)); nProtect++;
      PROTECT(betaStarSamples_r = allocMatrix(REALSXP, nAbundRE, nPost)); nProtect++;
    }
    // Likelihood samples for WAIC. 
    SEXP likeSamples_r;
    PROTECT(likeSamples_r = allocMatrix(REALSXP, nObs, nPost)); nProtect++;
    
    /********************************************************************
      Some constants and temporary variables to be used later
    ********************************************************************/
    int JpAbund = J * pAbund; 
    int nObspAbund = nObs * pAbund;
    int JJ = J * J;
    double tmp_0; 
    double *tmp_one = (double *) R_alloc(inc, sizeof(double)); 
    double *tmp_ppAbund = (double *) R_alloc(ppAbund, sizeof(double)); 
    double *tmp_pAbund = (double *) R_alloc(pAbund, sizeof(double));
    double *tmp_pAbund2 = (double *) R_alloc(pAbund, sizeof(double));
    double *tmp_nObs = (double *) R_alloc(nObs, sizeof(double)); 
    double *tmp_nObspAbund = (double *) R_alloc(nObspAbund, sizeof(double)); 
    double *tmp_JpAbund = (double *) R_alloc(JpAbund, sizeof(double));
    double *tmp_J1 = (double *) R_alloc(J, sizeof(double));
   
    // For latent abundance and WAIC
    double *like = (double *) R_alloc(nObs, sizeof(double)); zeros(like, nObs);
    double *psi = (double *) R_alloc(nObs, sizeof(double)); 
    zeros(psi, nObs); 
    double *mu = (double *) R_alloc(nObs, sizeof(double)); 
    zeros(mu, nObs); 

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
    double *theta = (double *) R_alloc(nTheta, sizeof(double));
    SEXP thetaSamples_r; 
    PROTECT(thetaSamples_r = allocMatrix(REALSXP, nTheta, nPost)); nProtect++; 
    double a, v, b, e, muNNGP, var, aij; 
    // Initiate spatial values
    theta[sigmaSqIndx] = REAL(sigmaSqStarting_r)[0]; 
    theta[phiIndx] = REAL(phiStarting_r)[0]; 
    double phi = theta[phiIndx];
    double sigmaSq = theta[sigmaSqIndx];
    double nu = REAL(nuStarting_r)[0]; 
    if (corName == "matern") {
      theta[nuIndx] = nu; 
    } 

    double *C = (double *) R_alloc(JJ, sizeof(double));
    double *CCand = (double *) R_alloc(JJ, sizeof(double));
    double *tmp_JD = (double *) R_alloc(J, sizeof(double));
    double *tmp_JD2 = (double *) R_alloc(J, sizeof(double));
    double *R = (double *) R_alloc(JJ, sizeof(double)); 
    if (sigmaSqIG) {
      spCorLT(coordsD, J, theta, corName, R); 
    }
    spCovLT(coordsD, J, theta, corName, C); 
    F77_NAME(dpotrf)(lower, &J, C, &J, &info FCONE); 
    if(info != 0){error("c++ error: Cholesky failed in initial covariance matrix\n");}
    F77_NAME(dpotri)(lower, &J, C, &J, &info FCONE); 
    if(info != 0){error("c++ error: Cholesky inverse failed in initial covariance matrix\n");}
    // For sigmaSq sampler
    double aSigmaSqPost = 0.5 * J + sigmaSqA; 
    double bSigmaSqPost = 0.0; 
    double *wTRInv = (double *) R_alloc(J, sizeof(double)); 

    /********************************************************************
      Set up MH stuff
    ********************************************************************/
    double logPostBetaCurr = 0.0, logPostBetaCand = 0.0;
    double logPostKappaCurr = 0.0, logPostKappaCand = 0.0;
    double logPostThetaCurr = 0.0, logPostThetaCand = 0.0;
    double *logPostWCand = (double *) R_alloc(J, sizeof(double));
    double *logPostWCurr = (double *) R_alloc(J, sizeof(double));
    for (j = 0; j < J; j++) {
      logPostWCurr[j] = R_NegInf;
      logPostWCand[j] = logPostWCurr[j];
    }
    double *logPostBetaStarCand = (double *) R_alloc(nAbundRE, sizeof(double));
    double *logPostBetaStarCurr = (double *) R_alloc(nAbundRE, sizeof(double));
    for (j = 0; j < nAbundRE; j++) {
      logPostBetaStarCurr[j] = R_NegInf;
      logPostBetaStarCand[j] = logPostBetaStarCurr[j];
    }
    double logDet, detCand, detCurr; 
    double phiCand = 0.0, nuCand = 0.0, sigmaSqCand = 0.0;  
    double *betaCand = (double *) R_alloc(pAbund, sizeof(double)); 
    for (j = 0; j < pAbund; j++) {
      betaCand[j] = beta[j];
    } 
    double *wCand = (double *) R_alloc(J, sizeof(double));
    for (j = 0; j < J; j++) {
      wCand[j] = w[j];
    }
    double *betaStarCand = (double *) R_alloc(nAbundRE, sizeof(double));
    for (j = 0; j < nAbundRE; j++) {
      betaStarCand[j] = betaStar[j];
    }
    double kappaCand = 0.0;
    kappaCand = kappa;
    // theta, beta, and w
    int nAMCMC = 0;
    if (pAbundRE > 0) {
      nAMCMC = nTheta + pAbund + J + nAbundRE;
    } else {
      nAMCMC = nTheta + pAbund + J;
    }
    if (family == 1) {
      nAMCMC++;
    }
    int betaAMCMCIndx = 0;
    int sigmaSqAMCMCIndx = betaAMCMCIndx + pAbund; 
    int phiAMCMCIndx = sigmaSqAMCMCIndx + 1; 
    int nuAMCMCIndx; 
    if (corName == "matern") {
      nuAMCMCIndx = phiAMCMCIndx + 1;
    } else {
      nuAMCMCIndx = phiAMCMCIndx;
    }
    int wAMCMCIndx = nuAMCMCIndx + 1;
    int betaStarAMCMCIndx = wAMCMCIndx + J;
    int kappaAMCMCIndx = betaStarAMCMCIndx + nAbundRE; 
    double *accept = (double *) R_alloc(nAMCMC, sizeof(double)); zeros(accept, nAMCMC); 
    // TODO: will probably want to cut this back eventually. 
    SEXP acceptSamples_r; 
    PROTECT(acceptSamples_r = allocMatrix(REALSXP, nAMCMC, nBatch)); nProtect++; 
    SEXP tuningSamples_r; 
    PROTECT(tuningSamples_r = allocMatrix(REALSXP, nAMCMC, nBatch)); nProtect++; 

    /**********************************************************************
     * Prep for random effects
     * *******************************************************************/
    // Site-level sums of the abundance random effects
    double *betaStarSites = (double *) R_alloc(nObs, sizeof(double)); 
    double *betaStarSitesCand = (double *) R_alloc(nObs, sizeof(double)); 
    zeros(betaStarSites, nObs); 
    // Initial sums
    for (j = 0; j < nObs; j++) {
      for (l = 0; l < pAbundRE; l++) {
        betaStarSites[j] += betaStar[which(XRE[l * nObs + j], betaLevelIndx, nAbundRE)] * 
		            XRandom[l * nObs + j];
      }
      betaStarSitesCand[j] = betaStarSites[j];
    }
    // Starting index for abundance random effects
    int *betaStarStart = (int *) R_alloc(pAbundRE, sizeof(int)); 
    for (l = 0; l < pAbundRE; l++) {
      betaStarStart[l] = which(l, betaStarIndx, nAbundRE); 
    }

    logPostBetaCurr = R_NegInf;
    logPostThetaCurr = R_NegInf;
    GetRNGstate(); 

    for (s = 0, g = 0; s < nBatch; s++) {
      for (t = 0; t < batchLength; t++, g++) {
        /********************************************************************
         *Update Abundance Regression Coefficients
         *******************************************************************/
        // Proposal
        for (k = 0; k < pAbund; k++) {
          logPostBetaCand = 0.0;
	  logPostBetaCurr = 0.0;
          betaCand[k] = rnorm(beta[k], exp(tuning[betaAMCMCIndx + k]));
          for (i = 0; i < pAbund; i++) {
            logPostBetaCand += dnorm(betaCand[i], muBeta[i], sqrt(SigmaBeta[i * pAbund + i]), 1);
	    logPostBetaCurr += dnorm(beta[i], muBeta[i], sqrt(SigmaBeta[i * pAbund + i]), 1);
          }
          for (j = 0; j < nObs; j++) {
            tmp_nObs[j] = exp(F77_NAME(ddot)(&pAbund, &X[j], &nObs, betaCand, &inc) + betaStarSites[j] + 
			      w[siteIndx[j]]);
	    if (family == 1) {
              logPostBetaCand += dnbinom_mu(y[j], kappa, tmp_nObs[j], 1);
	    } else {
              logPostBetaCand += dpois(y[j], tmp_nObs[j], 1);
	    }
            tmp_nObs[j] = exp(F77_NAME(ddot)(&pAbund, &X[j], &nObs, beta, &inc) + betaStarSites[j] + 
			      w[siteIndx[j]]);
	    if (family == 1) {
              logPostBetaCurr += dnbinom_mu(y[j], kappa, tmp_nObs[j], 1);
	    } else {
              logPostBetaCurr += dpois(y[j], tmp_nObs[j], 1);
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
         *Update abundance random effects variance
         *******************************************************************/
        for (l = 0; l < pAbundRE; l++) {
          tmp_0 = F77_NAME(ddot)(&nAbundRELong[l], &betaStar[betaStarStart[l]], &inc, &betaStar[betaStarStart[l]], &inc); 
          tmp_0 *= 0.5; 
          sigmaSqMu[l] = rigamma(sigmaSqMuA[l] + nAbundRELong[l] / 2.0, sigmaSqMuB[l] + tmp_0); 
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
              if (XRE[betaStarIndx[l] * nObs + j] == betaLevelIndx[l]) {
                // Candidate
                betaStarSitesCand[j] = 0.0;
                for (ll = 0; ll < pAbundRE; ll++) {
                  betaStarSitesCand[j] += betaStarCand[which(XRE[ll * nObs + j], 
				                         betaLevelIndx, nAbundRE)] * 
	                              XRandom[ll * nObs + j];
                }
                tmp_nObs[j] = exp(F77_NAME(ddot)(&pAbund, &X[j], &nObs, beta, &inc) + 
				  betaStarSitesCand[j] + w[siteIndx[j]]);
                if (family == 1) {
		  logPostBetaStarCand[l] += dnbinom_mu(y[j], kappa, tmp_nObs[j], 1);
		} else {
		  logPostBetaStarCand[l] += dpois(y[j], tmp_nObs[j], 1);
		}
		// Current
                betaStarSites[j] = 0.0;
                for (ll = 0; ll < pAbundRE; ll++) {
                  betaStarSites[j] += betaStar[which(XRE[ll * nObs + j], 
				               betaLevelIndx, nAbundRE)] * 
	                              XRandom[ll * nObs + j];
                }
                tmp_nObs[j] = exp(F77_NAME(ddot)(&pAbund, &X[j], &nObs, beta, &inc) + 
				  betaStarSites[j] + w[siteIndx[j]]);
                if (family == 1) {
		  logPostBetaStarCurr[l] += dnbinom_mu(y[j], kappa, tmp_nObs[j], 1);
		} else {
		  logPostBetaStarCurr[l] += dpois(y[j], tmp_nObs[j], 1);
		}
	      }
	    }
	    if (runif (0.0, 1.0) <= exp(logPostBetaStarCand[l] - logPostBetaStarCurr[l])) {
              betaStar[l] = betaStarCand[l];
	      F77_NAME(dcopy)(&nObs, betaStarSitesCand, &inc, betaStarSites, &inc);
	      accept[betaStarAMCMCIndx + l]++;
	    } else {
              betaStarCand[l] = betaStar[l];
	      F77_NAME(dcopy)(&nObs, betaStarSites, &inc, betaStarSitesCand, &inc);
	    }
	  }
	}

	/********************************************************************
         *Update sigmaSq
         *******************************************************************/
        if (sigmaSqIG) {
	  // Get inverse correlation matrix in reverse from inverse covariance matrix
	  // Remember: C currently contains the inverse of covariance matrix. 
	  fillUTri(C, J); 
	  for (j = 0; j < JJ; j++) {
            R[j] = theta[sigmaSqIndx] * C[j]; 
	  } // j
	  // Compute t(w) %*% R^-1 %*% w / 
	  // t(w) %*% R^-1
	  // Def a better way to do this operation. 
	  for (j = 0; j < J; j++) {
            wTRInv[j] = F77_NAME(ddot)(&J, &R[j], &J, w, &inc);  
          } // j
	  bSigmaSqPost = F77_NAME(ddot)(&J, wTRInv, &inc, w, &inc); 
	  bSigmaSqPost /= 2.0; 
	  bSigmaSqPost += sigmaSqB; 
	  theta[sigmaSqIndx] = rigamma(aSigmaSqPost, bSigmaSqPost); 
	}

        /********************************************************************
         *Update phi (and nu if matern and sigmaSq if uniform prior)
         *******************************************************************/
	if (corName == "matern") {
          nu = theta[nuIndx]; 
	  nuCand = logitInv(rnorm(logit(theta[nuIndx], nuA, nuB), exp(tuning[nuAMCMCIndx])), nuA, nuB); 
          theta[nuIndx] = nuCand; 
        }
	phi = theta[phiIndx]; 
	phiCand = logitInv(rnorm(logit(phi, phiA, phiB), exp(tuning[phiAMCMCIndx])), phiA, phiB); 
	theta[phiIndx] = phiCand; 
	if (sigmaSqIG == 0) {
	  sigmaSq = theta[sigmaSqIndx]; 
	  sigmaSqCand = logitInv(rnorm(logit(sigmaSq, sigmaSqA, sigmaSqB), 
				 exp(tuning[sigmaSqAMCMCIndx])), sigmaSqA, sigmaSqB); 
	  theta[sigmaSqIndx] = sigmaSqCand; 
	}

	// Construct covariance matrix (stored in C). 
	spCovLT(coordsD, J, theta, corName, CCand); 

        /********************************
         * Proposal
         *******************************/
	// Invert CCand and log det cov. 
        detCand = 0.0;
	F77_NAME(dpotrf)(lower, &J, CCand, &J, &info FCONE); 
	if(info != 0){error("c++ error: Cholesky failed in covariance matrix\n");}
	// Get log of the determinant of the covariance matrix. 
	for (k = 0; k < J; k++) {
	  detCand += 2.0 * log(CCand[k*J+k]);
	} // k
	F77_NAME(dpotri)(lower, &J, CCand, &J, &info FCONE); 
	if(info != 0){error("c++ error: Cholesky inverse failed in covariance matrix\n");}
        logPostThetaCand = 0.0; 
	// Jacobian and Uniform prior. 
	logPostThetaCand += log(phiCand - phiA) + log(phiB - phiCand); 
	// (-1/2) * tmp_JD` *  C^-1 * tmp_JD
	F77_NAME(dsymv)(lower, &J, &one,  CCand, &J, w, &inc, &zero, tmp_JD, &inc FCONE);
	logPostThetaCand += -0.5*detCand-0.5*F77_NAME(ddot)(&J, w, &inc, tmp_JD, &inc);
        if (corName == "matern"){
          logPostThetaCand += log(nuCand - nuA) + log(nuB - nuCand); 
        }
	if (sigmaSqIG == 0) {
          logPostThetaCand += log(sigmaSqCand - sigmaSqA) + log(sigmaSqB - sigmaSqCand);
	}

        /********************************
         * Current
         *******************************/
	if (corName == "matern") {
	  theta[nuIndx] = nu; 
	}
	theta[phiIndx] = phi; 
	if (sigmaSqIG == 0) {
          theta[sigmaSqIndx] = sigmaSq;
	}
	spCovLT(coordsD, J, theta, corName, C); 
        detCurr = 0.0;
	F77_NAME(dpotrf)(lower, &J, C, &J, &info FCONE); 
	if(info != 0){error("c++ error: Cholesky failed in covariance matrix\n");}
	for (k = 0; k < J; k++) {
	  detCurr += 2.0 * log(C[k*J+k]);
	} // k
	F77_NAME(dpotri)(lower, &J, C, &J, &info FCONE); 
	if(info != 0){error("c++ error: Cholesky inverse failed in covariance matrix\n");}
        logPostThetaCurr = 0.0; 
	logPostThetaCurr += log(phi - phiA) + log(phiB - phi); 
	// (-1/2) * tmp_JD` *  C^-1 * tmp_JD
	F77_NAME(dsymv)(lower, &J, &one, C, &J, w, &inc, &zero, tmp_JD, &inc FCONE);
	logPostThetaCurr += -0.5*detCurr-0.5*F77_NAME(ddot)(&J, w, &inc, tmp_JD, &inc);
        if (corName == "matern"){
          logPostThetaCurr += log(nu - nuA) + log(nuB - nu); 
        }
	if (sigmaSqIG == 0) {
          logPostThetaCurr += log(sigmaSq - sigmaSqA) + log(sigmaSqB - sigmaSq);
	}

	// MH Accept/Reject
	if (runif(0.0, 1.0) <= exp(logPostThetaCand - logPostThetaCurr)) {
          theta[phiIndx] = phiCand;
          accept[phiAMCMCIndx]++;
          if (corName == "matern") {
            theta[nuIndx] = nuCand; 
            accept[nuAMCMCIndx]++; 
          }
	  if (sigmaSqIG == 0) {
            theta[sigmaSqIndx] = sigmaSqCand;
	    accept[sigmaSqAMCMCIndx]++;
	  }
	  F77_NAME(dcopy)(&JJ, CCand, &inc, C, &inc); 
        }

        /********************************************************************
         *Update w (spatial random effects)
         *******************************************************************/
	for (j = 0; j < J; j++) {
          // Proposal
          a = 0.0;
          // Propose new value
	  logPostWCand[j] = 0.0;
	  wCand[j] = rnorm(w[j], exp(tuning[wAMCMCIndx + j]));
          /********************************
           * Proposal
           *******************************/
	  // (-1/2) * tmp_JD` *  C^-1 * tmp_JD
	  F77_NAME(dsymv)(lower, &J, &one, C, &J, wCand, &inc, &zero, tmp_JD, &inc FCONE);
	  logPostWCand[j] += -0.5*F77_NAME(ddot)(&J, wCand, &inc, tmp_JD, &inc);
          for (i = 0; i < nObs; i++) {
            if (siteIndx[i] == j) { 
              tmp_nObs[i] = exp(F77_NAME(ddot)(&pAbund, &X[i], &nObs, beta, &inc) + 
                                betaStarSites[i] + wCand[j]);
              if (family == 1) {
                logPostWCand[j] += dnbinom_mu(y[i], kappa, tmp_nObs[i], 1);
	      } else {
                logPostWCand[j] += dpois(y[i], tmp_nObs[i], 1);
	      }
	    }
	  }
          /********************************
           * Current 
           *******************************/
	  logPostWCurr[j] = 0.0;
	  // (-1/2) * tmp_JD` *  C^-1 * tmp_JD
	  F77_NAME(dsymv)(lower, &J, &one, C, &J, w, &inc, &zero, tmp_JD, &inc FCONE);
	  logPostWCurr[j] += -0.5*F77_NAME(ddot)(&J, w, &inc, tmp_JD, &inc);
          for (i = 0; i < nObs; i++) {
            if (siteIndx[i] == j) { 
              tmp_nObs[i] = exp(F77_NAME(ddot)(&pAbund, &X[i], &nObs, beta, &inc) + 
                                betaStarSites[i] + w[j]);
              if (family == 1) {
                logPostWCurr[j] += dnbinom_mu(y[i], kappa, tmp_nObs[i], 1);
	      } else {
                logPostWCurr[j] += dpois(y[i], tmp_nObs[i], 1);
	      }
	    }
	  }
	  if (runif(0.0, 1.0) <= exp(logPostWCand[j] - logPostWCurr[j])) {
	    w[j] = wCand[j];
	    accept[wAMCMCIndx + j]++;
	  } else {
            wCand[j] = w[j];
	  }
        }

        /********************************************************************
         *Update kappa (the NB size parameter)
         *******************************************************************/
	if (family == 1) {
          kappaCand = logitInv(rnorm(logit(kappa, kappaA, kappaB), exp(tuning[kappaAMCMCIndx])), 
			       kappaA, kappaB);
	  logPostKappaCurr = 0.0;
	  logPostKappaCand = 0.0;
	  for (j = 0; j < nObs; j++) {
            mu[j] = exp(F77_NAME(ddot)(&pAbund, &X[j], &nObs, beta, &inc) + 
                        betaStarSites[j] + w[siteIndx[j]]);
            logPostKappaCurr += dnbinom_mu(y[j], kappa, mu[j], 1);
	    logPostKappaCand += dnbinom_mu(y[j], kappaCand, mu[j], 1);
	  }
	  // Jacobian adjustment
	  logPostKappaCurr += log(kappa - kappaA) + log(kappaB - kappa);
	  logPostKappaCand += log(kappaCand - kappaA) + log(kappaB - kappaCand);
          if (runif(0.0, 1.0) <= exp(logPostKappaCand - logPostKappaCurr)) {
            kappa = kappaCand;
	    accept[kappaAMCMCIndx]++;
	  }
	}

        /********************************************************************
         *Get fitted values
         *******************************************************************/
        for (j = 0; j < nObs; j++) {
          // Only calculate if Poisson since it's already calculated in kappa update
          if (family == 0) {
            mu[j] = exp(F77_NAME(ddot)(&pAbund, &X[j], &nObs, beta, &inc) + 
                        betaStarSites[j] + w[siteIndx[j]]);
            yRep[j] = rpois(mu[j]);
            like[j] = dpois(y[j], mu[j], 0);
	  } else {
            yRep[j] = rnbinom_mu(kappa, mu[j]);
            like[j] = dnbinom_mu(y[j], kappa, mu[j], 0);
	  }
        }

        /********************************************************************
         *Save samples
         *******************************************************************/
        if (g >= nBurn) {
          thinIndx++; 
          if (thinIndx == nThin) {
            F77_NAME(dcopy)(&pAbund, beta, &inc, &REAL(betaSamples_r)[sPost*pAbund], &inc);
            F77_NAME(dcopy)(&nObs, mu, &inc, &REAL(muSamples_r)[sPost*nObs], &inc); 
            F77_NAME(dcopy)(&nObs, yRep, &inc, &REAL(yRepSamples_r)[sPost*nObs], &inc); 
	    F77_NAME(dcopy)(&nTheta, theta, &inc, &REAL(thetaSamples_r)[sPost*nTheta], &inc); 
            F77_NAME(dcopy)(&J, w, &inc, &REAL(wSamples_r)[sPost*J], &inc); 
	    if (family == 1) {
	      REAL(kappaSamples_r)[sPost] = kappa;
	    }
            if (pAbundRE > 0) {
              F77_NAME(dcopy)(&pAbundRE, sigmaSqMu, &inc, 
          		    &REAL(sigmaSqMuSamples_r)[sPost*pAbundRE], &inc);
              F77_NAME(dcopy)(&nAbundRE, betaStar, &inc, 
          		    &REAL(betaStarSamples_r)[sPost*nAbundRE], &inc);
            }
            F77_NAME(dcopy)(&nObs, like, &inc, 
          		  &REAL(likeSamples_r)[sPost*nObs], &inc);
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
	  Rprintf("\tphi\t\t%3.1f\t\t%1.5f\n", 100.0*REAL(acceptSamples_r)[s * nAMCMC + phiAMCMCIndx], exp(tuning[phiAMCMCIndx]));
	  if (corName == "matern") {
	    Rprintf("\tnu\t\t%3.1f\t\t%1.5f\n", 100.0*REAL(acceptSamples_r)[s * nAMCMC + nuAMCMCIndx], exp(tuning[nuAMCMCIndx]));
	  }
	  if (sigmaSqIG == 0) {
	    Rprintf("\tsigmaSq\t\t%3.1f\t\t%1.5f\n", 100.0*REAL(acceptSamples_r)[s * nAMCMC + sigmaSqAMCMCIndx], exp(tuning[sigmaSqAMCMCIndx]));
	  }
	  if (family == 1) {
            Rprintf("\tkappa\t\t%3.1f\t\t%1.5f\n", 100.0*REAL(acceptSamples_r)[s * nAMCMC + kappaAMCMCIndx], exp(tuning[kappaAMCMCIndx]));
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
    int nResultListObjs = 6;
    if (pAbundRE > 0) {
      nResultListObjs += 2;
    }
    if (family == 1) {
      nResultListObjs += 1;
    }

    PROTECT(result_r = allocVector(VECSXP, nResultListObjs)); nProtect++;
    PROTECT(resultName_r = allocVector(VECSXP, nResultListObjs)); nProtect++;

    // Setting the components of the output list.
    SET_VECTOR_ELT(result_r, 0, betaSamples_r);
    SET_VECTOR_ELT(result_r, 1, yRepSamples_r); 
    SET_VECTOR_ELT(result_r, 2, muSamples_r);
    SET_VECTOR_ELT(result_r, 3, likeSamples_r);
    SET_VECTOR_ELT(result_r, 4, wSamples_r);
    SET_VECTOR_ELT(result_r, 5, thetaSamples_r);
    if (pAbundRE > 0) {
      SET_VECTOR_ELT(result_r, 6, sigmaSqMuSamples_r);
      SET_VECTOR_ELT(result_r, 7, betaStarSamples_r);
    }
    if (family == 1) {
      if (pAbundRE > 0) {
        tmp_0 = 8;
      } else {
        tmp_0 = 6;
      }
      SET_VECTOR_ELT(result_r, tmp_0, kappaSamples_r);
    }

    SET_VECTOR_ELT(resultName_r, 0, mkChar("beta.samples")); 
    SET_VECTOR_ELT(resultName_r, 1, mkChar("y.rep.samples")); 
    SET_VECTOR_ELT(resultName_r, 2, mkChar("mu.samples"));
    SET_VECTOR_ELT(resultName_r, 3, mkChar("like.samples"));
    SET_VECTOR_ELT(resultName_r, 4, mkChar("w.samples"));
    SET_VECTOR_ELT(resultName_r, 5, mkChar("theta.samples"));
    if (pAbundRE > 0) {
      SET_VECTOR_ELT(resultName_r, 6, mkChar("sigma.sq.mu.samples")); 
      SET_VECTOR_ELT(resultName_r, 7, mkChar("beta.star.samples")); 
    }
    if (family == 1) {
      SET_VECTOR_ELT(resultName_r, tmp_0, mkChar("kappa.samples")); 
    }
   
    namesgets(result_r, resultName_r);
    
    UNPROTECT(nProtect);
    
    return(result_r);
  }
}

