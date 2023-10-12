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
  SEXP lfMsAbund(SEXP y_r, SEXP X_r, SEXP XRE_r, SEXP XRandom_r,
	         SEXP consts_r, SEXP nAbundRELong_r, 
	         SEXP betaStarting_r, SEXP kappaStarting_r, SEXP betaCommStarting_r, 
		 SEXP tauSqBetaStarting_r, SEXP lambdaStarting_r, SEXP wStarting_r,
	         SEXP sigmaSqMuStarting_r, SEXP betaStarStarting_r,  
	         SEXP siteIndx_r, SEXP betaStarIndx_r, SEXP betaLevelIndx_r, 
	         SEXP muBetaComm_r, SEXP SigmaBetaComm_r, SEXP kappaA_r, 
	         SEXP kappaB_r, SEXP tauSqBetaA_r, SEXP tauSqBetaB_r, 
	         SEXP sigmaSqMuA_r, SEXP sigmaSqMuB_r, SEXP tuning_r,  
		 SEXP nBatch_r, SEXP batchLength_r, SEXP acceptRate_r,
	         SEXP nThreads_r, SEXP verbose_r, SEXP nReport_r, 
	         SEXP samplesInfo_r, SEXP chainInfo_r, SEXP family_r, SEXP offset_r){
   
    /**********************************************************************
     * Initial constants
     * *******************************************************************/
    int i, j, s, t, g, l, h, r, ll, ii, k, info, nProtect=0;
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
    double *offset = REAL(offset_r);
    int *XRE = INTEGER(XRE_r); 
    double *XRandom = REAL(XRandom_r);
    // Load constants
    int nSp = INTEGER(consts_r)[0]; 
    int J = INTEGER(consts_r)[1];
    int nObs = INTEGER(consts_r)[2]; 
    int pAbund = INTEGER(consts_r)[3];
    int pAbundRE = INTEGER(consts_r)[4];
    int nAbundRE = INTEGER(consts_r)[5];
    int q = INTEGER(consts_r)[6]; 
    int saveFitted = INTEGER(consts_r)[7];
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
          Rprintf("Latent Factor Multi-species Negative Binomial\nAbundance model fit with %i sites and %i species.\n\n", J, nSp);
	} else {
          Rprintf("Latent Factor Multi-species Poisson Abundance\nmodel fit with %i sites and %i species.\n\n", J, nSp);
	}
        Rprintf("Samples per Chain: %i \n", nSamples);
        Rprintf("Burn-in: %i \n", nBurn); 
        Rprintf("Thinning Rate: %i \n", nThin); 
	Rprintf("Number of Chains: %i \n", nChain);
        Rprintf("Total Posterior Samples: %i \n\n", nPost * nChain); 
        Rprintf("Using %i latent factors.\n", q);
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
       Some constants and temporary variables to be used later
     * *******************************************************************/
    int pAbundnSp = pAbund * nSp; 
    int nObsnSp = nObs * nSp; 
    int nAbundREnSp = nAbundRE * nSp; 
    int JnSp = J * nSp;
    int nObspAbundRE = nObs * pAbundRE;
    int Jq = J * q;
    int nSpq = nSp * q;
    double tmp_0; 
    double *tmp_ppAbund = (double *) R_alloc(ppAbund, sizeof(double)); 
    double *tmp_pAbund = (double *) R_alloc(pAbund, sizeof(double));
    double *tmp_pAbund2 = (double *) R_alloc(pAbund, sizeof(double));
    int *tmp_JInt = (int *) R_alloc(J, sizeof(int));
    for (j = 0; j < J; j++) {
      tmp_JInt[j] = 0; 
    }
    double *tmp_nObs = (double *) R_alloc(nObs, sizeof(double)); 

    /**********************************************************************
     * Parameters
     * *******************************************************************/
    double *betaComm = (double *) R_alloc(pAbund, sizeof(double)); 
    F77_NAME(dcopy)(&pAbund, REAL(betaCommStarting_r), &inc, betaComm, &inc);
    double *tauSqBeta = (double *) R_alloc(pAbund, sizeof(double)); 
    F77_NAME(dcopy)(&pAbund, REAL(tauSqBetaStarting_r), &inc, tauSqBeta, &inc);
    double *beta = (double *) R_alloc(pAbundnSp, sizeof(double));   
    F77_NAME(dcopy)(&pAbundnSp, REAL(betaStarting_r), &inc, beta, &inc);
    // Latent factor stuff 
    double *w = (double *) R_alloc(Jq, sizeof(double));
    F77_NAME(dcopy)(&Jq, REAL(wStarting_r), &inc, w, &inc);
    // Latent factor loadings
    double *lambda = (double *) R_alloc(nSpq, sizeof(double));
    F77_NAME(dcopy)(&nSpq, REAL(lambdaStarting_r), &inc, lambda, &inc);
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

    /**********************************************************************
     * Return Stuff
     * *******************************************************************/
    // Community level
    SEXP betaCommSamples_r; 
    PROTECT(betaCommSamples_r = allocMatrix(REALSXP, pAbund, nPost)); nProtect++;
    zeros(REAL(betaCommSamples_r), pAbund * nPost);
    SEXP tauSqBetaSamples_r; 
    PROTECT(tauSqBetaSamples_r = allocMatrix(REALSXP, pAbund, nPost)); nProtect++; 
    zeros(REAL(tauSqBetaSamples_r), pAbund * nPost);
    // Species level
    SEXP betaSamples_r;
    PROTECT(betaSamples_r = allocMatrix(REALSXP, pAbundnSp, nPost)); nProtect++;
    zeros(REAL(betaSamples_r), pAbundnSp * nPost);
    SEXP yRepSamples_r; 
    SEXP muSamples_r; 
    SEXP likeSamples_r;
    if (saveFitted == 1) {
      PROTECT(yRepSamples_r = allocMatrix(REALSXP, nObsnSp, nPost)); nProtect++; 
      zeros(REAL(yRepSamples_r), nObsnSp * nPost);
      PROTECT(muSamples_r = allocMatrix(REALSXP, nObsnSp, nPost)); nProtect++; 
      zeros(REAL(muSamples_r), nObsnSp * nPost);
      PROTECT(likeSamples_r = allocMatrix(REALSXP, nObsnSp, nPost)); nProtect++;
      zeros(REAL(likeSamples_r), nObsnSp * nPost);
    }
    // Spatial parameters
    SEXP lambdaSamples_r; 
    PROTECT(lambdaSamples_r = allocMatrix(REALSXP, nSpq, nPost)); nProtect++;
    zeros(REAL(lambdaSamples_r), nSpq * nPost);
    SEXP wSamples_r; 
    PROTECT(wSamples_r = allocMatrix(REALSXP, Jq, nPost)); nProtect++; 
    zeros(REAL(wSamples_r), Jq * nPost);
    // Abundance random effects
    SEXP sigmaSqMuSamples_r; 
    SEXP betaStarSamples_r; 
    if (pAbundRE > 0) {
      PROTECT(sigmaSqMuSamples_r = allocMatrix(REALSXP, pAbundRE, nPost)); nProtect++;
      zeros(REAL(sigmaSqMuSamples_r), pAbundRE * nPost);
      PROTECT(betaStarSamples_r = allocMatrix(REALSXP, nAbundREnSp, nPost)); nProtect++;
      zeros(REAL(betaStarSamples_r), nAbundREnSp * nPost);
    }
    // Overdispersion
    SEXP kappaSamples_r;
    if (family == 1) {
      PROTECT(kappaSamples_r = allocMatrix(REALSXP, nSp, nPost)); nProtect++;
      zeros(REAL(kappaSamples_r), nSp * nPost);
    }
    
    /**********************************************************************
     * Additional Sampler Prep
     * *******************************************************************/
    double *like = (double *) R_alloc(nObsnSp, sizeof(double)); 
    zeros(like, nObsnSp);
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

    /**********************************************************************
     * Prep for random effects
     * *******************************************************************/
    // Site-level sums of the occurrence random effects
    double *betaStarSites = (double *) R_alloc(nObsnSp, sizeof(double)); 
    zeros(betaStarSites, nObsnSp); 
    double *betaStarSitesCand = (double *) R_alloc(nObsnSp, sizeof(double)); 
    int *betaStarLongIndx = (int *) R_alloc(nObspAbundRE, sizeof(int));
    // Initial sums (initiate with the first species)
    for (j = 0; j < nObs; j++) {
      for (l = 0; l < pAbundRE; l++) {
        betaStarLongIndx[l * nObs + j] = which(XRE[l * nObs + j], betaLevelIndx, nAbundRE);
        for (i = 0; i < nSp; i++) {
          betaStarSites[i * nObs + j] += betaStar[i * nAbundRE + betaStarLongIndx[l * nObs + j]] * XRandom[l * nObs + j];
	  betaStarSitesCand[i * nObs + j] = betaStarSites[i * nObs + j];
        }
      }
    }
    // Starting index for abundance random effects
    int *betaStarStart = (int *) R_alloc(pAbundRE, sizeof(int)); 
    for (l = 0; l < pAbundRE; l++) {
      betaStarStart[l] = which(l, betaStarIndx, nAbundRE); 
    }

    /**********************************************************************
     Latent factor stuff 
     * *******************************************************************/
    // Species-level latent factor effects 
    double *wStar = (double *) R_alloc(JnSp, sizeof(double)); zeros(wStar, JnSp);
    // Multiply Lambda %*% w[j] to get wStar. 
    for (j = 0; j < J; j++) {
      F77_NAME(dgemv)(ntran, &nSp, &q, &one, lambda, &nSp, &w[j*q], &inc, &zero, &wStar[j * nSp], &inc FCONE);
    }

    /********************************************************************
      Set up MH stuff
    ********************************************************************/
    double logPostBetaCurr = 0.0, logPostBetaCand = 0.0;
    double logPostKappaCurr = 0.0, logPostKappaCand = 0.0;
    double *logPostBetaStarCand = (double *) R_alloc(nAbundRE, sizeof(double));
    double *logPostBetaStarCurr = (double *) R_alloc(nAbundRE, sizeof(double));
    for (j = 0; j < nAbundRE; j++) {
      logPostBetaStarCurr[j] = R_NegInf;
      logPostBetaStarCand[j] = logPostBetaStarCurr[j];
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
    int nAMCMC = 0;
    if (pAbundRE > 0) {
      nAMCMC = (pAbund + nAbundRE) * nSp + nSpq + Jq;
    } else {
      nAMCMC = pAbund * nSp + nSpq + Jq;
    }
    if (family == 1) {
      nAMCMC += nSp;
    }
    int betaAMCMCIndx = 0;
    int lambdaAMCMCIndx = betaAMCMCIndx + pAbundnSp;
    int wAMCMCIndx = lambdaAMCMCIndx + nSpq;
    int betaStarAMCMCIndx = wAMCMCIndx + Jq;
    int kappaAMCMCIndx = betaStarAMCMCIndx + nAbundRE * nSp;

    double *accept = (double *) R_alloc(nAMCMC, sizeof(double)); zeros(accept, nAMCMC); 
    SEXP acceptSamples_r; 
    PROTECT(acceptSamples_r = allocMatrix(REALSXP, nAMCMC, nBatch)); nProtect++; 
    zeros(REAL(acceptSamples_r), nAMCMC * nBatch);
    SEXP tuningSamples_r; 
    PROTECT(tuningSamples_r = allocMatrix(REALSXP, nAMCMC, nBatch)); nProtect++; 
    zeros(REAL(tuningSamples_r), nAMCMC * nBatch);
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
    double kappaCand = 0.0;
    kappaCand = kappa[0];
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
            for (j = 0; j < nObs; j++) {
              tmp_nObs[j] = exp(F77_NAME(ddot)(&pAbund, &X[j], &nObs, &betaCand[i], &nSp) + 
                                betaStarSites[i * nObs + j] + wStar[siteIndx[j] * nSp + i]);
	      if (family == 1) {
                logPostBetaCand += dnbinom_mu(y[j * nSp + i], kappa[i], tmp_nObs[j] * offset[j], 1);
	      } else {
                logPostBetaCand += dpois(y[j * nSp + i], tmp_nObs[j] * offset[j], 1);
	      }
              tmp_nObs[j] = exp(F77_NAME(ddot)(&pAbund, &X[j], &nObs, &beta[i], &nSp) + 
                                betaStarSites[i * nObs + j] + wStar[siteIndx[j] * nSp + i]);
	      if (family == 1) {
                logPostBetaCurr += dnbinom_mu(y[j * nSp + i], kappa[i], tmp_nObs[j] * offset[j], 1);
	      } else {
                logPostBetaCurr += dpois(y[j * nSp + i], tmp_nObs[j] * offset[j], 1);
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
	      for (j = 0; j < nObs; j++) {
                if (XRE[betaStarIndx[l] * nObs + j] == betaLevelIndx[l]) {
                  // Candidate
                  betaStarSitesCand[i * nObs + j] = 0.0;
                  for (ll = 0; ll < pAbundRE; ll++) {
                    betaStarSitesCand[i * nObs + j] += betaStarCand[i * nAbundRE + betaStarLongIndx[ll * nObs + j]] * 
	                                XRandom[ll * nObs + j];
                  }
                  tmp_nObs[j] = exp(F77_NAME(ddot)(&pAbund, &X[j], &nObs, &beta[i], &nSp) + 
	  			  betaStarSitesCand[i * nObs + j] + wStar[siteIndx[j] * nSp + i]);
                  if (family == 1) {
	  	    logPostBetaStarCand[l] += dnbinom_mu(y[j * nSp + i], kappa[i], 
                                                         tmp_nObs[j] * offset[j], 1);
	  	  } else {
	  	    logPostBetaStarCand[l] += dpois(y[j * nSp + i], tmp_nObs[j] * offset[j], 1);
	  	  }
	  	  // Current
                    betaStarSites[i * nObs + j] = 0.0;
                    for (ll = 0; ll < pAbundRE; ll++) {
                      betaStarSites[i * nObs + j] += betaStar[i * nAbundRE + betaStarLongIndx[ll * nObs + j]] * 
	                                  XRandom[ll * nObs + j];
                    }
                    tmp_nObs[j] = exp(F77_NAME(ddot)(&pAbund, &X[j], &nObs, &beta[i], &nSp) + 
	  	  		  betaStarSites[i * nObs + j] + wStar[siteIndx[j] * nSp + i]);
                    if (family == 1) {
	  	    logPostBetaStarCurr[l] += dnbinom_mu(y[j * nSp + i], kappa[i], 
                                                         tmp_nObs[j] * offset[j], 1);
	  	  } else {
	  	    logPostBetaStarCurr[l] += dpois(y[j * nSp + i], tmp_nObs[j] * offset[j], 1);
	  	  }
	        }
	      }
	      if (runif (0.0, 1.0) <= exp(logPostBetaStarCand[l] - logPostBetaStarCurr[l])) {
                betaStar[i * nAbundRE + l] = betaStarCand[i * nAbundRE + l];
	        F77_NAME(dcopy)(&nObsnSp, betaStarSitesCand, &inc, betaStarSites, &inc);
	        accept[betaStarAMCMCIndx + i * nAbundRE + l]++;
	      } else {
                betaStarCand[i * nAbundRE + l] = betaStar[i * nAbundRE + l];
	        F77_NAME(dcopy)(&nObsnSp, betaStarSites, &inc, betaStarSitesCand, &inc);
	      }
	    }
	  }
	} // species loop

        /********************************************************************
         *Update Latent Factors (w) 
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
	    // N(0, 1) prior for w
	    logPostWCand[j * q + ll] = dnorm(wCand[j * q + ll], 0.0, 1.0, 1);
            /*****************************
	     *Current
             *****************************/
            logPostWCurr[j * q + ll] = dnorm(w[j * q + ll], 0.0, 1.0, 1);
	    // Likelihood for proposal
	    for (i = 0; i < nSp; i++) {
	      wStarCand[j * nSp + i] = 0.0;
	      for (ii = 0; ii < q; ii++) {
                wStarCand[j * nSp + i] += lambda[ii * nSp + i] * wCand[j * q + ii];
	      } // ii 
              for (r = 0; r < nObs; r++) {
                if (siteIndx[r] == j) {
                  /*****************************
	           *Candidate
                   *****************************/
                  tmp_nObs[r] = exp(F77_NAME(ddot)(&pAbund, &X[r], &nObs, &beta[i], &nSp) + 
                                    betaStarSites[i * nObs + r] + wStarCand[j * nSp + i]);
	          if (family == 1) {
                    logPostWCand[j * q + ll] += dnbinom_mu(y[r * nSp + i], kappa[i], 
                                                           tmp_nObs[r] * offset[r], 1);
	          } else {
                    logPostWCand[j * q + ll] += dpois(y[r * nSp + i], tmp_nObs[r] * offset[r], 1);
	          }
                  /*****************************
	           *Current
                   *****************************/
                  tmp_nObs[r] = exp(F77_NAME(ddot)(&pAbund, &X[r], &nObs, &beta[i], &nSp) + 
                                    betaStarSites[i * nObs + r] + wStar[j * nSp + i]);
	          if (family == 1) {
                    logPostWCurr[j * q + ll] += dnbinom_mu(y[r * nSp + i], kappa[i], 
                                                           tmp_nObs[r] * offset[r], 1);
	          } else {
                    logPostWCurr[j * q + ll] += dpois(y[r * nSp + i], tmp_nObs[r] * offset[r], 1);
	          }
	        }
	      } // r
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
         *Update factor loadings (lambda) 
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
	      } // j
              for (r = 0; r < nObs; r++) {
                // Candidate
                tmp_nObs[r] = exp(F77_NAME(ddot)(&pAbund, &X[r], &nObs, &beta[i], &nSp) + 
                                  betaStarSites[i * nObs + r] + wStarCand[siteIndx[r] * nSp + i]);
	        if (family == 1) {
                  logPostLambdaCand[ll * nSp + i] += dnbinom_mu(y[r * nSp + i], kappa[i], 
                                                                tmp_nObs[r] * offset[r], 1);
	        } else {
                  logPostLambdaCand[ll * nSp + i] += dpois(y[r * nSp + i], tmp_nObs[r] * offset[r], 1);
	        }
		// Current
                tmp_nObs[r] = exp(F77_NAME(ddot)(&pAbund, &X[r], &nObs, &beta[i], &nSp) + 
                                  betaStarSites[i * nObs + r] + wStar[siteIndx[r] * nSp + i]);
	        if (family == 1) {
                  logPostLambdaCurr[ll * nSp + i] += dnbinom_mu(y[r * nSp + i], kappa[i], 
                                                                tmp_nObs[r] * offset[r], 1);
	        } else {
                  logPostLambdaCurr[ll * nSp + i] += dpois(y[r * nSp + i], tmp_nObs[r] * offset[r], 1);
	        }
	      } // r
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
	    for (j = 0; j < nObs; j++) {
              mu[j * nSp + i] = exp(F77_NAME(ddot)(&pAbund, &X[j], &nObs, &beta[i], &nSp) + 
                          betaStarSites[i * nObs + j] + wStar[siteIndx[j] * nSp + i]);
              logPostKappaCurr += dnbinom_mu(y[j * nSp + i], kappa[i], mu[j * nSp + i] * offset[j], 1);
	      logPostKappaCand += dnbinom_mu(y[j * nSp + i], kappaCand, mu[j * nSp + i] * offset[j], 1);
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
           *Get fitted values
           *******************************************************************/
	  if (saveFitted == 1) {
            for (r = 0; r < nObs; r++) {
              // Only calculate mu if Poisson since it's already calculated in kappa update
	      if (family == 0) {
                mu[r * nSp + i] = exp(F77_NAME(ddot)(&pAbund, &X[r], &nObs, &beta[i], &nSp) + 
                                      betaStarSites[i * nObs + r] + wStar[siteIndx[r] * nSp + i]);
                yRep[r * nSp + i] = rpois(mu[r * nSp + i] * offset[r]);
                like[r * nSp + i] = dpois(y[r * nSp + i], mu[r * nSp + i] * offset[r], 0);
	      } else {
                yRep[r * nSp + i] = rnbinom_mu(kappa[i], mu[r * nSp + i] * offset[r]);
                like[r * nSp + i] = dnbinom_mu(y[r * nSp + i], kappa[i], mu[r * nSp + i] * offset[r], 0);
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
            F77_NAME(dcopy)(&pAbund, tauSqBeta, &inc, &REAL(tauSqBetaSamples_r)[sPost*pAbund], &inc);
            F77_NAME(dcopy)(&pAbundnSp, beta, &inc, &REAL(betaSamples_r)[sPost*pAbundnSp], &inc); 
	    if (family == 1) {
              F77_NAME(dcopy)(&nSp, kappa, &inc, &REAL(kappaSamples_r)[sPost*nSp], &inc); 
	    }
	    if (saveFitted == 1) {
              F77_NAME(dcopy)(&nObsnSp, yRep, &inc, &REAL(yRepSamples_r)[sPost*nObsnSp], &inc); 
              F77_NAME(dcopy)(&nObsnSp, mu, &inc, &REAL(muSamples_r)[sPost*nObsnSp], &inc); 
              F77_NAME(dcopy)(&nObsnSp, like, &inc, &REAL(likeSamples_r)[sPost*nObsnSp], &inc); 
	    }
            F77_NAME(dcopy)(&Jq, w, &inc, &REAL(wSamples_r)[sPost*Jq], &inc); 
            F77_NAME(dcopy)(&nSpq, lambda, &inc, &REAL(lambdaSamples_r)[sPost*nSpq], &inc); 
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
    int nResultListObjs = 11;
    if (pAbundRE > 0) {
      nResultListObjs += 2;
    }
    if (family == 1) {
      nResultListObjs += 1;
    }

    PROTECT(result_r = allocVector(VECSXP, nResultListObjs)); nProtect++;
    PROTECT(resultName_r = allocVector(VECSXP, nResultListObjs)); nProtect++;

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
    SET_VECTOR_ELT(result_r, 7, tuningSamples_r); 
    SET_VECTOR_ELT(result_r, 8, acceptSamples_r); 
    if (pAbundRE > 0) {
      SET_VECTOR_ELT(result_r, 10, sigmaSqMuSamples_r);
      SET_VECTOR_ELT(result_r, 11, betaStarSamples_r);
    }
    if (family == 1) {
      if (pAbundRE > 0) {
        tmp_0 = 12;
      } else {
        tmp_0 = 10;
      }
      SET_VECTOR_ELT(result_r, tmp_0, kappaSamples_r);
    }

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
    SET_VECTOR_ELT(resultName_r, 7, mkChar("tune")); 
    SET_VECTOR_ELT(resultName_r, 8, mkChar("accept")); 
    if (pAbundRE > 0) {
      SET_VECTOR_ELT(resultName_r, 10, mkChar("sigma.sq.mu.samples")); 
      SET_VECTOR_ELT(resultName_r, 11, mkChar("beta.star.samples")); 
    }
    if (family == 1) {
      SET_VECTOR_ELT(resultName_r, tmp_0, mkChar("kappa.samples")); 
    }
   
    namesgets(result_r, resultName_r);
    
    UNPROTECT(nProtect);
    
    return(result_r);
  }
}


