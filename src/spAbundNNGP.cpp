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

void updateBFPNNGP(double *B, double *F, double *c, double *C, double *coords, int *nnIndx, int *nnIndxLU, int n, int m, double sigmaSq, double phi, double nu, int covModel, double *bk, double nuUnifb){

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
  SEXP spAbundNNGP(SEXP y_r, SEXP X_r, SEXP coords_r, SEXP XRE_r,
                   SEXP XRandom_r, SEXP consts_r, SEXP nAbundRELong_r,
                   SEXP m_r, SEXP nnIndx_r,
                   SEXP nnIndxLU_r, SEXP uIndx_r, SEXP uIndxLU_r,
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
                   SEXP chainInfo_r, SEXP sigmaSqIG_r, SEXP family_r, SEXP offset_r){

    /**********************************************************************
     * Initial constants
     * *******************************************************************/
    int i, g, t, j, k, jj, s, r, l, ll, nProtect=0;
    const int inc = 1;

    /**********************************************************************
     * Get Inputs
     * *******************************************************************/
    double *y = REAL(y_r);
    double *X = REAL(X_r);
    double *offset = REAL(offset_r);
    int *XRE = INTEGER(XRE_r);
    double *XRandom = REAL(XRandom_r);
    // Load constants
    int J = INTEGER(consts_r)[0];
    int nObs = INTEGER(consts_r)[1];
    int pAbund = INTEGER(consts_r)[2];
    int pAbundRE = INTEGER(consts_r)[3];
    int nAbundRE = INTEGER(consts_r)[4];
    int saveFitted = INTEGER(consts_r)[5];
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
    int covModel = INTEGER(covModel_r)[0];
    std::string corName = getCorName(covModel);
    int currChain = INTEGER(chainInfo_r)[0];
    double acceptRate = REAL(acceptRate_r)[0];
    double *tuning = REAL(tuning_r);
    double *coords = REAL(coords_r);
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
	if (family == 1) {
          Rprintf("Spatial NNGP Negative Binomial Abundance model fit with %i sites.\n\n", J);
	} else {
          Rprintf("Spatial NNGP Poisson Abundance model fit with %i sites.\n\n", J);
	}
        Rprintf("Samples per Chain: %i (%i batches of length %i)\n", nSamples, nBatch, batchLength);
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
    // Abundance Covariates
    double *beta = (double *) R_alloc(pAbund, sizeof(double));
    F77_NAME(dcopy)(&pAbund, REAL(betaStarting_r), &inc, beta, &inc);
    // Abundance random effect variances
    double *sigmaSqMu = (double *) R_alloc(pAbundRE, sizeof(double));
    F77_NAME(dcopy)(&pAbundRE, REAL(sigmaSqMuStarting_r), &inc, sigmaSqMu, &inc);
    // Overdispersion parameter for NB;
    double kappa = REAL(kappaStarting_r)[0];
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
    PROTECT(betaSamples_r = Rf_allocMatrix(REALSXP, pAbund, nPost)); nProtect++;
    zeros(REAL(betaSamples_r), pAbund * nPost);
    SEXP yRepSamples_r;
    SEXP muSamples_r;
    SEXP likeSamples_r;
    if (saveFitted == 1) {
      PROTECT(yRepSamples_r = Rf_allocMatrix(REALSXP, nObs, nPost)); nProtect++;
      zeros(REAL(yRepSamples_r), nObs * nPost);
      PROTECT(muSamples_r = Rf_allocMatrix(REALSXP, nObs, nPost)); nProtect++;
      zeros(REAL(muSamples_r), nObs * nPost);
      PROTECT(likeSamples_r = Rf_allocMatrix(REALSXP, nObs, nPost)); nProtect++;
      zeros(REAL(likeSamples_r), nObs * nPost);
    }
    SEXP wSamples_r;
    PROTECT(wSamples_r = Rf_allocMatrix(REALSXP, J, nPost)); nProtect++;
    zeros(REAL(wSamples_r), J * nPost);
    SEXP kappaSamples_r;
    if (family == 1) {
      PROTECT(kappaSamples_r = Rf_allocMatrix(REALSXP, inc, nPost)); nProtect++;
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

    /********************************************************************
      Some constants and temporary variables to be used later
    ********************************************************************/
    double tmp_0;
    double *tmp_nObs = (double *) R_alloc(nObs, sizeof(double));

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
    PROTECT(thetaSamples_r = Rf_allocMatrix(REALSXP, nTheta, nPost)); nProtect++;
    zeros(REAL(thetaSamples_r), nTheta * nPost);
    double a, b, e;
    // Initiate spatial values
    theta[sigmaSqIndx] = REAL(sigmaSqStarting_r)[0];
    theta[phiIndx] = REAL(phiStarting_r)[0];
    double nu = REAL(nuStarting_r)[0];
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

    double *bk = (double *) R_alloc(nThreads*(1.0+static_cast<int>(floor(nuB))), sizeof(double));

    if (corName == "matern") {
      nu = theta[nuIndx];
    }
    updateBFPNNGP(B, F, c, C, coords, nnIndx, nnIndxLU, J, m,
               theta[sigmaSqIndx], theta[phiIndx], nu, covModel, bk, nuB);

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
    double logDet;
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
          logPostBetaCand += dnorm(betaCand[k], muBeta[k], sqrt(SigmaBeta[k * pAbund + k]), 1);
	  logPostBetaCurr += dnorm(beta[k], muBeta[k], sqrt(SigmaBeta[k * pAbund + k]), 1);
          for (j = 0; j < nObs; j++) {
            tmp_nObs[j] = exp(F77_NAME(ddot)(&pAbund, &X[j], &nObs, betaCand, &inc) + betaStarSites[j] +
			      w[siteIndx[j]]);
	    if (family == 1) {
              logPostBetaCand += dnbinom_mu(y[j], kappa, tmp_nObs[j] * offset[j], 1);
	    } else {
              logPostBetaCand += dpois(y[j], tmp_nObs[j] * offset[j], 1);
	    }
            tmp_nObs[j] = exp(F77_NAME(ddot)(&pAbund, &X[j], &nObs, beta, &inc) + betaStarSites[j] +
			      w[siteIndx[j]]);
	    if (family == 1) {
              logPostBetaCurr += dnbinom_mu(y[j], kappa, tmp_nObs[j] * offset[j], 1);
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
		  logPostBetaStarCand[l] += dnbinom_mu(y[j], kappa, tmp_nObs[j] * offset[j], 1);
		} else {
		  logPostBetaStarCand[l] += dpois(y[j], tmp_nObs[j] * offset[j], 1);
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
		  logPostBetaStarCurr[l] += dnbinom_mu(y[j], kappa, tmp_nObs[j] * offset[j], 1);
		} else {
		  logPostBetaStarCurr[l] += dpois(y[j], tmp_nObs[j] * offset[j], 1);
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
         *Update w (spatial random effects)
         *******************************************************************/
	for (j = 0; j < J; j++) {
          // Proposal
          a = 0.0;
          // Propose new value
	  logPostWCand[j] = 0.0;
	  wCand[j] = rnorm(w[j], exp(tuning[wAMCMCIndx + j]));
	  // MVN for any neighbors of j
	  if (uIndxLU[J + j] > 0) { // if current location j is a neighbor for anybody
            for (r = 0; r < uIndxLU[J + j]; r++) { // how many locations have j as a neighbor
              jj = uIndx[uIndxLU[j] + r]; // jj is the index of the rth location who has j as a neighbor
              e = 0;
              for (i = 0; i < nnIndxLU[J+jj]; i++){ // neighbors of the jjth location
                e += B[nnIndxLU[jj]+i]*wCand[nnIndx[nnIndxLU[jj]+i]];
              }
              b = wCand[jj] - e;
              a += b*b/F[jj];
            }
	  }
	  // MVN for j
          if(nnIndxLU[J+j] > 0){ // if j has any neighbors
            e = 0;
            for(i = 0; i < nnIndxLU[J+j]; i++){
              e += B[nnIndxLU[j]+i]*wCand[nnIndx[nnIndxLU[j]+i]];
            }
            b = wCand[j] - e;
          }else{
            b = wCand[j];
          }
          a += b*b/F[j];
          logPostWCand[j] = -0.5*a;
          for (i = 0; i < nObs; i++) {
            if (siteIndx[i] == j) {
              tmp_nObs[i] = exp(F77_NAME(ddot)(&pAbund, &X[i], &nObs, beta, &inc) +
                                betaStarSites[i] + wCand[j]);
              if (family == 1) {
                logPostWCand[j] += dnbinom_mu(y[i], kappa, tmp_nObs[i] * offset[i], 1);
	      } else {
                logPostWCand[j] += dpois(y[i], tmp_nObs[i] * offset[i], 1);
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
                e += B[nnIndxLU[jj]+i]*w[nnIndx[nnIndxLU[jj]+i]];
              }
              b = w[jj] - e;
              a += b*b/F[jj];
            }
	  }
	  // MVN for j
          if(nnIndxLU[J+j] > 0){ // if j has any neighbors
            e = 0;
            for(i = 0; i < nnIndxLU[J+j]; i++){
              e += B[nnIndxLU[j]+i]*w[nnIndx[nnIndxLU[j]+i]];
            }
            b = w[j] - e;
          }else{
            b = w[j];
          }
          a += b*b/F[j];
          logPostWCurr[j] = -0.5*a;
          for (i = 0; i < nObs; i++) {
            if (siteIndx[i] == j) {
              tmp_nObs[i] = exp(F77_NAME(ddot)(&pAbund, &X[i], &nObs, beta, &inc) +
                                betaStarSites[i] + w[j]);
	      if (family == 1) {
                logPostWCurr[j] += dnbinom_mu(y[i], kappa, tmp_nObs[i] * offset[i], 1);
	      } else {
                logPostWCurr[j] += dpois(y[i], tmp_nObs[i] * offset[i], 1);
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
         *Update sigmaSq
         *******************************************************************/
	if (sigmaSqIG == 1) {
	  a = 0;
#ifdef _OPENMP
#pragma omp parallel for private (e, i, b) reduction(+:a)
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
        updateBFPNNGP(B, F, c, C, coords, nnIndx, nnIndxLU, J, m, theta[sigmaSqIndx],
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
				exp(tuning[phiAMCMCIndx])), phiA, phiB);
	if (sigmaSqIG == 0) {
	  sigmaSqCand = logitInv(rnorm(logit(theta[sigmaSqIndx], sigmaSqA, sigmaSqB),
				 exp(tuning[sigmaSqAMCMCIndx])), sigmaSqA, sigmaSqB);
	}
        if (corName == "matern"){
      	  nuCand = logitInv(rnorm(logit(theta[nuIndx], nuA, nuB), exp(tuning[nuAMCMCIndx])), nuA, nuB);
        }

        if (sigmaSqIG) {
          updateBFPNNGP(BCand, FCand, c, C, coords, nnIndx, nnIndxLU, J, m, theta[sigmaSqIndx], phiCand, nuCand, covModel, bk, nuB);
	  } else {
            updateBFPNNGP(BCand, FCand, c, C, coords, nnIndx, nnIndxLU, J, m, sigmaSqCand, phiCand, nuCand, covModel, bk, nuB);
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
          accept[phiAMCMCIndx]++;
          if(corName == "matern"){
            theta[nuIndx] = nuCand;
            accept[nuAMCMCIndx]++;
          }
	  if (sigmaSqIG == 0) {
            theta[sigmaSqIndx] = sigmaSqCand;
	    accept[sigmaSqAMCMCIndx]++;
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
            logPostKappaCurr += dnbinom_mu(y[j], kappa, mu[j] * offset[j], 1);
	    logPostKappaCand += dnbinom_mu(y[j], kappaCand, mu[j] * offset[j], 1);
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
	if (saveFitted == 1) {
          for (j = 0; j < nObs; j++) {
            // Only calculate if Poisson since it's already calculated in kappa update
            if (family == 0) {
              mu[j] = exp(F77_NAME(ddot)(&pAbund, &X[j], &nObs, beta, &inc) +
                          betaStarSites[j] + w[siteIndx[j]]);
              yRep[j] = rpois(mu[j] * offset[j]);
              like[j] = dpois(y[j], mu[j] * offset[j], 0);
	    } else {
              yRep[j] = rnbinom_mu(kappa, mu[j] * offset[j]);
              like[j] = dnbinom_mu(y[j], kappa, mu[j] * offset[j], 0);
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
	    if (saveFitted == 1) {
              F77_NAME(dcopy)(&nObs, mu, &inc, &REAL(muSamples_r)[sPost*nObs], &inc);
              F77_NAME(dcopy)(&nObs, yRep, &inc, &REAL(yRepSamples_r)[sPost*nObs], &inc);
              F77_NAME(dcopy)(&nObs, like, &inc,
          		      &REAL(likeSamples_r)[sPost*nObs], &inc);
	    }
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

    PROTECT(result_r = Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;
    PROTECT(resultName_r = Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;

    // Setting the components of the output list.
    SET_VECTOR_ELT(result_r, 0, betaSamples_r);
    if (saveFitted == 1) {
      SET_VECTOR_ELT(result_r, 1, yRepSamples_r);
      SET_VECTOR_ELT(result_r, 2, muSamples_r);
      SET_VECTOR_ELT(result_r, 3, likeSamples_r);
    }
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

    SET_VECTOR_ELT(resultName_r, 0, Rf_mkChar("beta.samples"));
    if (saveFitted == 1) {
      SET_VECTOR_ELT(resultName_r, 1, Rf_mkChar("y.rep.samples"));
      SET_VECTOR_ELT(resultName_r, 2, Rf_mkChar("mu.samples"));
      SET_VECTOR_ELT(resultName_r, 3, Rf_mkChar("like.samples"));
    }
    SET_VECTOR_ELT(resultName_r, 4, Rf_mkChar("w.samples"));
    SET_VECTOR_ELT(resultName_r, 5, Rf_mkChar("theta.samples"));
    if (pAbundRE > 0) {
      SET_VECTOR_ELT(resultName_r, 6, Rf_mkChar("sigma.sq.mu.samples"));
      SET_VECTOR_ELT(resultName_r, 7, Rf_mkChar("beta.star.samples"));
    }
    if (family == 1) {
      SET_VECTOR_ELT(resultName_r, tmp_0, Rf_mkChar("kappa.samples"));
    }

    Rf_namesgets(result_r, resultName_r);

    UNPROTECT(nProtect);

    return(result_r);
  }
}


