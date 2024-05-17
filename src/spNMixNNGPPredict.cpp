#define USE_FC_LEN_T
#include <string>
#include "util.h"

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

  SEXP spNMixNNGPPredict(SEXP coords_r, SEXP J_r, SEXP family_r,
                         SEXP pAbund_r, SEXP m_r, SEXP X0_r, SEXP coords0_r,
                         SEXP J0_r, SEXP nnIndx0_r, SEXP betaSamples_r,
                         SEXP thetaSamples_r, SEXP kappaSamples_r, SEXP wSamples_r,
                         SEXP betaStarSiteSamples_r,
                         SEXP sitesLink_r, SEXP sites0Sampled_r, SEXP sites0_r,
                         SEXP nSamples_r, SEXP covModel_r, SEXP nThreads_r, SEXP verbose_r,
                         SEXP nReport_r){

    int i, k, l, s, info, nProtect=0;
    const int inc = 1;
    const double one = 1.0;
    const double zero = 0.0;
    char const *lower = "L";

    double *coords = REAL(coords_r);
    int J = INTEGER(J_r)[0];
    int pAbund = INTEGER(pAbund_r)[0];
    int family = INTEGER(family_r)[0];

    double *X0 = REAL(X0_r);
    double *coords0 = REAL(coords0_r);
    int J0 = INTEGER(J0_r)[0];
    int m = INTEGER(m_r)[0];
    int mm = m * m;
    int *sitesLink = INTEGER(sitesLink_r);
    int *sites0Sampled = INTEGER(sites0Sampled_r);

    int *nnIndx0 = INTEGER(nnIndx0_r);
    double *beta = REAL(betaSamples_r);
    double *theta = REAL(thetaSamples_r);
    double *kappa = REAL(kappaSamples_r);
    double *w = REAL(wSamples_r);
    double *betaStarSite = REAL(betaStarSiteSamples_r);

    int nSamples = INTEGER(nSamples_r)[0];
    int covModel = INTEGER(covModel_r)[0];
    std::string corName = getCorName(covModel);
    int nThreads = INTEGER(nThreads_r)[0];
    int verbose = INTEGER(verbose_r)[0];
    int nReport = INTEGER(nReport_r)[0];

#ifdef _OPENMP
    omp_set_num_threads(nThreads);
#else
    if(nThreads > 1){
      Rf_warning("n.omp.threads > %i, but source not compiled with OpenMP support.", nThreads);
      nThreads = 1;
    }
#endif

    if(verbose){
      Rprintf("----------------------------------------\n");
      Rprintf("\tPrediction description\n");
      Rprintf("----------------------------------------\n");
      Rprintf("NNGP spatial abundance model with %i observations.\n\n", J);
      Rprintf("Number of covariates %i (including intercept if specified).\n\n", pAbund);
      Rprintf("Using the %s spatial correlation model.\n\n", corName.c_str());
      Rprintf("Using %i nearest neighbors.\n\n", m);
      Rprintf("Number of MCMC samples %i.\n\n", nSamples);
      Rprintf("Predicting at %i locations.\n", J0);
#ifdef _OPENMP
      Rprintf("\nSource compiled with OpenMP support and model fit using %i threads.\n", nThreads);
#else
      Rprintf("\nSource not compiled with OpenMP support.\n");
#endif
    }

    // parameters
    int nTheta, sigmaSqIndx,  phiIndx, nuIndx;

    if (corName != "matern") {
      nTheta = 2; //sigma^2, phi
      sigmaSqIndx = 0; phiIndx = 1;
      } else{
	nTheta = 3; //sigma^2, phi, nu
	sigmaSqIndx = 0; phiIndx = 1; nuIndx = 2;
      }

    // get max nu
    double nuMax = 0;
    int nb = 0;

    if(corName == "matern"){
      for(i = 0; i < nSamples; i++){
	if(theta[i*nTheta+nuIndx] > nuMax){
	  nuMax = theta[i*nTheta+nuIndx];
	}
      }

      nb = 1+static_cast<int>(floor(nuMax));
    }

    double *bk = (double *) R_alloc(nThreads*nb, sizeof(double));

    double *C = (double *) R_alloc(nThreads*mm, sizeof(double)); zeros(C, nThreads*mm);
    double *c = (double *) R_alloc(nThreads*m, sizeof(double)); zeros(c, nThreads*m);
    double *tmp_m  = (double *) R_alloc(nThreads*m, sizeof(double));
    double phi = 0, nu = 0, sigmaSq = 0, d;
    int threadID = 0, status = 0;

    SEXP N0_r, w0_r, mu0_r;
    PROTECT(N0_r = Rf_allocMatrix(REALSXP, J0, nSamples)); nProtect++;
    PROTECT(mu0_r = Rf_allocMatrix(REALSXP, J0, nSamples)); nProtect++;
    PROTECT(w0_r = Rf_allocMatrix(REALSXP, J0, nSamples)); nProtect++;
    double *N0 = REAL(N0_r);
    double *mu0 = REAL(mu0_r);
    double *w0 = REAL(w0_r);

    if (verbose) {
      Rprintf("-------------------------------------------------\n");
      Rprintf("\t\tPredicting\n");
      Rprintf("-------------------------------------------------\n");
      #ifdef Win32
        R_FlushConsole();
      #endif
    }

    int vIndx = -1;
    double *wV = (double *) R_alloc(J0*nSamples, sizeof(double));

    GetRNGstate();

    for(i = 0; i < J0*nSamples; i++){
      wV[i] = rnorm(0.0,1.0);
    }

    for(i = 0; i < J0; i++){
#ifdef _OPENMP
#pragma omp parallel for private(threadID, phi, nu, sigmaSq, k, l, d, info)
#endif
      for(s = 0; s < nSamples; s++){
#ifdef _OPENMP
	threadID = omp_get_thread_num();
#endif
        if (sites0Sampled[i] == 1) {
          w0[s * J0 + i] = w[s * J + sitesLink[i]];
        } else {
	  phi = theta[s*nTheta+phiIndx];
	  if(corName == "matern"){
	    nu = theta[s*nTheta+nuIndx];
	  }
	  sigmaSq = theta[s*nTheta+sigmaSqIndx];

	  for(k = 0; k < m; k++){
	    d = dist2(coords[nnIndx0[i+J0*k]], coords[J+nnIndx0[i+J0*k]], coords0[i], coords0[J0+i]);
	    c[threadID*m+k] = sigmaSq*spCor(d, phi, nu, covModel, &bk[threadID*nb]);
	    for(l = 0; l < m; l++){
	      d = dist2(coords[nnIndx0[i+J0*k]], coords[J+nnIndx0[i+J0*k]], coords[nnIndx0[i+J0*l]], coords[J+nnIndx0[i+J0*l]]);
	      C[threadID*mm+l*m+k] = sigmaSq*spCor(d, phi, nu, covModel, &bk[threadID*nb]);
	    }
	  }

	  F77_NAME(dpotrf)(lower, &m, &C[threadID*mm], &m, &info FCONE);
	  if(info != 0){Rf_error("c++ error: dpotrf failed\n");}
	  F77_NAME(dpotri)(lower, &m, &C[threadID*mm], &m, &info FCONE);
	  if(info != 0){Rf_error("c++ error: dpotri failed\n");}

	  F77_NAME(dsymv)(lower, &m, &one, &C[threadID*mm], &m, &c[threadID*m], &inc, &zero, &tmp_m[threadID*m], &inc FCONE);

	  d = 0;
	  for(k = 0; k < m; k++){
	    d += tmp_m[threadID*m+k]*w[s*J+nnIndx0[i+J0*k]];
	  }

	  #ifdef _OPENMP
          #pragma omp atomic
          #endif
	  vIndx++;

	  w0[s*J0+i] = sqrt(sigmaSq - F77_NAME(ddot)(&m, &tmp_m[threadID*m], &inc, &c[threadID*m], &inc))*wV[vIndx] + d;
	}
	mu0[s*J0+i] = exp(F77_NAME(ddot)(&pAbund, &X0[i], &J0, &beta[s*pAbund], &inc) +
                          w0[s*J0+i] +
                          betaStarSite[s * J0 + i]);
      }

      if(verbose){
	if(status == nReport){
	  Rprintf("Location: %i of %i, %3.2f%%\n", i, J0, 100.0*i/J0);
          #ifdef Win32
	  R_FlushConsole();
          #endif
	  status = 0;
	}
      }
      status++;
      R_CheckUserInterrupt();
    }

    if(verbose){
      Rprintf("Location: %i of %i, %3.2f%%\n", i, J0, 100.0*i/J0);
      #ifdef Win32
      R_FlushConsole();
      #endif
    }

    // Generate latent occurrence state after the fact.
    // Temporary fix. Will embed this in the above loop at some point.
    if (verbose) {
      Rprintf("Generating latent abundance state\n");
    }
    for(i = 0; i < J0; i++){
      for(s = 0; s < nSamples; s++){
        if (family == 1) {
          N0[s * J0 + i] = rnbinom_mu(kappa[s], mu0[s * J0 + i]);
	} else {
          N0[s * J0 + i] = rpois(mu0[s * J0 + i]);
	}
      } // s
    } // i

    PutRNGstate();


    //make return object
    SEXP result_r, resultName_r;
    int nResultListObjs = 3;

    PROTECT(result_r = Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;
    PROTECT(resultName_r = Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;

    SET_VECTOR_ELT(result_r, 0, N0_r);
    SET_VECTOR_ELT(resultName_r, 0, Rf_mkChar("N.0.samples"));

    SET_VECTOR_ELT(result_r, 1, w0_r);
    SET_VECTOR_ELT(resultName_r, 1, Rf_mkChar("w.0.samples"));

    SET_VECTOR_ELT(result_r, 2, mu0_r);
    SET_VECTOR_ELT(resultName_r, 2, Rf_mkChar("mu.0.samples"));

    Rf_namesgets(result_r, resultName_r);

    //unprotect
    UNPROTECT(nProtect);

    return(result_r);

  }
}


