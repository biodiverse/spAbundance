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
  SEXP abund(SEXP y_r, SEXP X_r, SEXP XRE_r, 
            SEXP consts_r, SEXP K_r, SEXP nAbundRELong_r, 
            SEXP betaStarting_r, SEXP kappaStarting_r,
            SEXP sigmaSqMuStarting_r, SEXP betaStarStarting_r, 
            SEXP siteIndx_r, SEXP betaStarIndx_r, SEXP betaLevelIndx_r, 
	    SEXP muBeta_r, SEXP SigmaBeta_r, SEXP sigmaSqMuA_r, SEXP sigmaSqMuB_r, 
	    SEXP kappaA_r, SEXP kappaB_r, SEXP tuning_r,
	    SEXP nBatch_r, SEXP batchLength_r, SEXP acceptRate_r, SEXP nThreads_r, 
            SEXP verbose_r, SEXP nReport_r, SEXP samplesInfo_r,
	    SEXP chainInfo_r){
   
    /**********************************************************************
     * Initial constants
     * *******************************************************************/
    int i, g, t, j, s, r, l, info, nProtect=0;
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
    // Load constants
    int J = INTEGER(consts_r)[0];
    int nObs = INTEGER(consts_r)[1];
    int pAbund = INTEGER(consts_r)[2];
    int pAbundRE = INTEGER(consts_r)[3];
    int nAbundRE = INTEGER(consts_r)[4];
    int ppAbund = pAbund * pAbund; 
    double *muBeta = (double *) R_alloc(pAbund, sizeof(double));   
    F77_NAME(dcopy)(&pAbund, REAL(muBeta_r), &inc, muBeta, &inc);
    double *SigmaBetaInv = (double *) R_alloc(ppAbund, sizeof(double));   
    F77_NAME(dcopy)(&ppAbund, REAL(SigmaBeta_r), &inc, SigmaBetaInv, &inc);
    double *sigmaSqMuA = REAL(sigmaSqMuA_r); 
    double *sigmaSqMuB = REAL(sigmaSqMuB_r); 
    double kappaA = REAL(kappaA_r)[0];
    double kappaB = REAL(kappaB_r)[0];
    int *nAbundRELong = INTEGER(nAbundRELong_r); 
    double *K = REAL(K_r); 
    int *siteIndx = INTEGER(siteIndx_r); 
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
    double tuning = REAL(tuning_r)[0];
    int nChain = INTEGER(chainInfo_r)[1];
    int nThreads = INTEGER(nThreads_r)[0];
    int verbose = INTEGER(verbose_r)[0];
    int nReport = INTEGER(nReport_r)[0];
    int status = 0; 
    int thinIndx = 0;
    int sPost = 0;  

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
        Rprintf("Negative Binomial Abundance model with Polya-Gamma latent\nvariable fit with %i sites.\n\n", J);
        Rprintf("Samples per Chain: %i (%i batches of length %i)\n", nSamples, nBatch, batchLength);
        Rprintf("Burn-in: %i \n", nBurn); 
        Rprintf("Thinning Rate: %i \n", nThin); 
        Rprintf("Number of Chains: %i \n", nChain);
        Rprintf("Total Posterior Samples: %i \n\n", nPost * nChain); 
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
    // Latent random effects
    double *betaStar = (double *) R_alloc(nAbundRE, sizeof(double)); 
    F77_NAME(dcopy)(&nAbundRE, REAL(betaStarStarting_r), &inc, betaStar, &inc); 
    // Size parameter
    double kappa = REAL(kappaStarting_r)[0];
    // Latent Abundance
    double *yRep = (double *) R_alloc(nObs, sizeof(double)); zeros(yRep, nObs);
    // Auxiliary variables
    double *omegaAbund = (double *) R_alloc(nObs, sizeof(double)); zeros(omegaAbund, nObs);
    double *yStar = (double *) R_alloc(nObs, sizeof(double)); zeros(yStar, nObs);

    /**********************************************************************
     * Return Stuff
     * *******************************************************************/
    SEXP betaSamples_r;
    PROTECT(betaSamples_r = allocMatrix(REALSXP, pAbund, nPost)); nProtect++;
    SEXP yRepSamples_r; 
    PROTECT(yRepSamples_r = allocMatrix(REALSXP, nObs, nPost)); nProtect++; 
    SEXP kappaSamples_r;
    PROTECT(kappaSamples_r = allocMatrix(REALSXP, inc, nPost)); nProtect++;
    SEXP muSamples_r; 
    PROTECT(muSamples_r = allocMatrix(REALSXP, nObs, nPost)); nProtect++; 
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

    // For normal priors
    // Abundupancy regression coefficient priors. 
    F77_NAME(dpotrf)(lower, &pAbund, SigmaBetaInv, &pAbund, &info FCONE); 
    if(info != 0){error("c++ error: dpotrf SigmaBetaInv failed\n");}
    F77_NAME(dpotri)(lower, &pAbund, SigmaBetaInv, &pAbund, &info FCONE); 
    if(info != 0){error("c++ error: dpotri SigmaBetaInv failed\n");}
    double *SigmaBetaInvMuBeta = (double *) R_alloc(pAbund, sizeof(double)); 
    F77_NAME(dsymv)(lower, &pAbund, &one, SigmaBetaInv, &pAbund, muBeta, &inc, &zero, 
        	    SigmaBetaInvMuBeta, &inc FCONE);

    // For NB PG sampling (recommendations from Polson et al. 2013). 
    int trunc = 200;
    double *tmp_trunc = (double *) R_alloc(trunc, sizeof(double));
    
    /********************************************************************
      Set up MH stuff
    ********************************************************************/
    double accept = 0.0;
    double logPostCurr = 0.0, logPostCand = 0.0;
    double epsilon = 1.0;
    double kappaCand = 0.0;
    SEXP acceptSamples_r; 
    PROTECT(acceptSamples_r = allocMatrix(REALSXP, inc, nBatch)); nProtect++; 
    SEXP tuningSamples_r; 
    PROTECT(tuningSamples_r = allocMatrix(REALSXP, inc, nBatch)); nProtect++; 

    /**********************************************************************
     * Prep for random effects
     * *******************************************************************/
    // Site-level sums of the abundance random effects
    double *betaStarSites = (double *) R_alloc(nObs, sizeof(double)); 
    zeros(betaStarSites, nObs); 
    // Initial sums
    for (j = 0; j < nObs; j++) {
      for (l = 0; l < pAbundRE; l++) {
        betaStarSites[j] += betaStar[which(XRE[l * nObs + j], betaLevelIndx, nAbundRE)];
      }
    }
    // Starting index for abundance random effects
    int *betaStarStart = (int *) R_alloc(pAbundRE, sizeof(int)); 
    for (l = 0; l < pAbundRE; l++) {
      betaStarStart[l] = which(l, betaStarIndx, nAbundRE); 
    }

    GetRNGstate(); 

    for (s = 0, g = 0; s < nBatch; s++) {
      for (t = 0; t < batchLength; t++, g++) {
        /********************************************************************
         *Update Abundance Auxiliary Variables 
         *******************************************************************/
        for (j = 0; j < nObs; j++) {
          omegaAbund[j] = rpgGamma(kappa + y[j], F77_NAME(ddot)(&pAbund, &X[j], &nObs, beta, &inc) + betaStarSites[j], trunc, tmp_trunc);
        } // j
      
        /********************************************************************
         *Update Abundance Regression Coefficients
         *******************************************************************/
        for (j = 0; j < nObs; j++) {
          yStar[j] = (y[j] - kappa) / (2.0 * omegaAbund[j]);
          tmp_nObs[j] = (yStar[j] - betaStarSites[j]) * omegaAbund[j];
        } // j

        /********************************
         * Compute b.beta
         *******************************/
        F77_NAME(dgemv)(ytran, &nObs, &pAbund, &one, X, &nObs, tmp_nObs, &inc, &zero, tmp_pAbund, &inc FCONE); 	 
        for (j = 0; j < pAbund; j++) {
          tmp_pAbund[j] += SigmaBetaInvMuBeta[j]; 
        } // j 

        /********************************
         * Compute A.beta
         * *****************************/
        // tmp_Jp is X %*% omega. 
        for(j = 0; j < nObs; j++){
          for(i = 0; i < pAbund; i++){
            tmp_nObspAbund[i*nObs+j] = X[i*nObs+j]*omegaAbund[j];
          }
        }

        // This finishes off A.beta
        // 1 * X * tmp_Jp + 0 * tmp_pp = tmp_pp
        F77_NAME(dgemm)(ytran, ntran, &pAbund, &pAbund, &nObs, &one, X, &nObs, 
			tmp_nObspAbund, &nObs, &zero, tmp_ppAbund, &pAbund FCONE FCONE);
        for (j = 0; j < ppAbund; j++) {
          tmp_ppAbund[j] += SigmaBetaInv[j]; 
        } // j

        F77_NAME(dpotrf)(lower, &pAbund, tmp_ppAbund, &pAbund, &info FCONE); 
        if(info != 0){error("c++ error: dpotrf here failed\n");}
        F77_NAME(dpotri)(lower, &pAbund, tmp_ppAbund, &pAbund, &info FCONE); 
        if(info != 0){error("c++ error: dpotri here failed\n");}
        F77_NAME(dsymv)(lower, &pAbund, &one, tmp_ppAbund, &pAbund, tmp_pAbund,
		       	&inc, &zero, tmp_pAbund2, &inc FCONE);
        F77_NAME(dpotrf)(lower, &pAbund, tmp_ppAbund, &pAbund, &info FCONE); 
        if(info != 0){error("c++ error: dpotrf here failed\n");}
        mvrnorm(beta, tmp_pAbund2, tmp_ppAbund, pAbund);
        
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
                tmp_one[0] += (yStar[j] - F77_NAME(ddot)(&pAbund, &X[j], &nObs, beta, &inc) - 
          		    betaStarSites[j] + betaStar[l]) * omegaAbund[j];
                tmp_0 += omegaAbund[j];
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
          zeros(betaStarSites, nObs);
          for (j = 0; j < nObs; j++) {
            for (l = 0; l < pAbundRE; l++) {
              betaStarSites[j] += betaStar[which(XRE[l * nObs + j], betaLevelIndx, nAbundRE)];
            }
          }
        }
        
        /********************************************************************
         *Update kappa (the NB size parameter)
         *******************************************************************/
        for (j = 0; j < nObs; j++) {
          tmp_nObs[j] = F77_NAME(ddot)(&pAbund, &X[j], &nObs, beta, &inc) + betaStarSites[j];
          psi[j] = logitInv(tmp_nObs[j], zero, one);
        }
        /********************************
         * Current
         *******************************/
        // Log gamma function is lgammafn
        // Likelihood contribution
        logPostCurr = 0.0;
        for (j = 0; j < nObs; j++) {
          logPostCurr += lgammafn(y[j] + kappa) - lgammafn(kappa) + kappa * log(1 - psi[j]);
        }
        logPostCurr += log(kappa - kappaA) + log(kappaB - kappa);
        /********************************
         * Candidate
         *******************************/
        kappaCand = logitInv(rnorm(logit(kappa, kappaA, kappaB), exp(tuning)), kappaA, kappaB);
        logPostCand = 0.0;
        for (j = 0; j < nObs; j++) {
          logPostCand += lgammafn(y[j] + kappaCand) - lgammafn(kappaCand) + kappaCand * log(1 - psi[j]);
        }
        logPostCand += log(kappaCand - kappaA) + log(kappaB - kappaCand);

        if (runif(0.0,1.0) <= exp(logPostCand - logPostCurr)) {
          kappa = kappaCand;
          accept++;
        }
        /********************************************************************
         *Get fitted values
         *******************************************************************/
        for (j = 0; j < nObs; j++) {
          mu[j] = exp(tmp_nObs[j]) * kappa;
          yRep[j] = rnbinom_mu(kappa, mu[j]);
	  like[j] = dnbinom_mu(y[j], kappa, mu[j], 0);
        }

        /********************************************************************
         *Save samples
         *******************************************************************/
        if (g >= nBurn) {
          thinIndx++; 
          if (thinIndx == nThin) {
            F77_NAME(dcopy)(&pAbund, beta, &inc, &REAL(betaSamples_r)[sPost*pAbund], &inc);
	    REAL(kappaSamples_r)[sPost] = kappa;
            F77_NAME(dcopy)(&nObs, mu, &inc, &REAL(muSamples_r)[sPost*nObs], &inc); 
            F77_NAME(dcopy)(&nObs, yRep, &inc, &REAL(yRepSamples_r)[sPost*nObs], &inc); 
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
      REAL(acceptSamples_r)[s] = accept/batchLength; 
      REAL(tuningSamples_r)[s] = tuning; 
      if (accept / batchLength > acceptRate) {
        tuning += std::min(0.01, 1.0/sqrt(static_cast<double>(s)));
      } else{
          tuning -= std::min(0.01, 1.0/sqrt(static_cast<double>(s)));
        }
      accept = 0.0;

      /********************************************************************
       *Report 
       *******************************************************************/
      if (verbose) {
        if (status == nReport) {
          Rprintf("Batch: %i of %i, %3.2f%%\n", s, nBatch, 100.0*s/nBatch);
          Rprintf("\tParameter\tAcceptance\tTuning\n");	  
          Rprintf("\tkappa\t\t%3.1f\t\t%1.5f\n", 100.0*REAL(acceptSamples_r)[s], exp(tuning));
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
    int nResultListObjs = 5;
    if (pAbundRE > 0) {
      nResultListObjs += 2;
    }

    PROTECT(result_r = allocVector(VECSXP, nResultListObjs)); nProtect++;
    PROTECT(resultName_r = allocVector(VECSXP, nResultListObjs)); nProtect++;

    // Setting the components of the output list.
    SET_VECTOR_ELT(result_r, 0, betaSamples_r);
    SET_VECTOR_ELT(result_r, 1, kappaSamples_r);
    SET_VECTOR_ELT(result_r, 2, yRepSamples_r); 
    SET_VECTOR_ELT(result_r, 3, muSamples_r);
    SET_VECTOR_ELT(result_r, 4, likeSamples_r);
    if (pAbundRE > 0) {
      SET_VECTOR_ELT(result_r, 5, sigmaSqMuSamples_r);
      SET_VECTOR_ELT(result_r, 6, betaStarSamples_r);
    }

    SET_VECTOR_ELT(resultName_r, 0, mkChar("beta.samples")); 
    SET_VECTOR_ELT(resultName_r, 1, mkChar("kappa.samples")); 
    SET_VECTOR_ELT(resultName_r, 2, mkChar("y.rep.samples")); 
    SET_VECTOR_ELT(resultName_r, 3, mkChar("mu.samples"));
    SET_VECTOR_ELT(resultName_r, 4, mkChar("like.samples"));
    if (pAbundRE > 0) {
      SET_VECTOR_ELT(resultName_r, 5, mkChar("sigma.sq.mu.samples")); 
      SET_VECTOR_ELT(resultName_r, 6, mkChar("beta.star.samples")); 
    }
   
    namesgets(result_r, resultName_r);
    
    UNPROTECT(nProtect);
    
    return(result_r);
  }
}

