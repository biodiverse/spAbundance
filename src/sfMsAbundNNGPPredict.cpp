#define USE_FC_LEN_T
#include <string>
#include "util.h"

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

  SEXP sfMsAbundNNGPPredict(SEXP coords_r, SEXP J_r, SEXP nObs_r, SEXP family_r, SEXP nSp_r, 
		            SEXP q_r, SEXP pAbund_r, SEXP m_r, SEXP X0_r, SEXP coords0_r, 
			    SEXP J0_r, SEXP nObs0_r, SEXP sitesLink_r, SEXP sites0Sampled_r, 
			    SEXP sites0_r, SEXP nnIndx0_r, SEXP betaSamples_r, 
			    SEXP thetaSamples_r, SEXP kappaSamples_r, SEXP lambdaSamples_r, 
			    SEXP wSamples_r, SEXP betaStarSiteSamples_r, 
			    SEXP nSamples_r, SEXP covModel_r, SEXP nThreads_r, SEXP verbose_r, 
			    SEXP nReport_r){

    int i, j, k, l, s, info, nProtect=0, ll;
    const int inc = 1;
    const double one = 1.0;
    char const *ntran = "N";
    const double zero = 0.0;
    char const *lower = "L";
    
    double *coords = REAL(coords_r);
    int J = INTEGER(J_r)[0];
    int nSp = INTEGER(nSp_r)[0];
    int q = INTEGER(q_r)[0]; 
    int pAbund = INTEGER(pAbund_r)[0];
    int pAbundnSp = pAbund * nSp;
    int Jq = J * q;
    int nSpq = nSp * q;
    int family = INTEGER(family_r)[0];

    double *X0 = REAL(X0_r);
    double *coords0 = REAL(coords0_r);
    int J0 = INTEGER(J0_r)[0];
    int J0nSp = J0 * nSp;
    int J0q = J0 * q;
    int nObs0 = INTEGER(nObs0_r)[0];
    int nObs0nSp = nObs0 * nSp;
    int m = INTEGER(m_r)[0]; 
    int mm = m * m; 
    int *sitesLink = INTEGER(sitesLink_r);
    int *sites0Sampled = INTEGER(sites0Sampled_r);
    int *sites0 = INTEGER(sites0_r);

    int *nnIndx0 = INTEGER(nnIndx0_r);        
    double *beta = REAL(betaSamples_r);
    double *kappa = REAL(kappaSamples_r);
    double *theta = REAL(thetaSamples_r);
    double *lambda = REAL(lambdaSamples_r);
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
      warning("n.omp.threads > 1, but source not compiled with OpenMP support.");
      nThreads = 1;
    }
#endif
    
    if(verbose){
      Rprintf("----------------------------------------\n");
      Rprintf("\tPrediction description\n");
      Rprintf("----------------------------------------\n");
      Rprintf("Spatial Factor NNGP GLMM with %i observations.\n\n", J);
      Rprintf("Number of covariates %i (including intercept if specified).\n\n", pAbund);
      Rprintf("Using the %s spatial correlation model.\n\n", corName.c_str());
      Rprintf("Using %i nearest neighbors.\n", m);
      Rprintf("Using %i latent spatial factors.\n\n", q);
      Rprintf("Number of MCMC samples %i.\n\n", nSamples);
      Rprintf("Predicting at %i non-sampled locations.\n\n", J0);  
#ifdef _OPENMP
      Rprintf("\nSource compiled with OpenMP support and model fit using %i threads.\n", nThreads);
#else
      Rprintf("\n\nSource not compiled with OpenMP support.\n");
#endif
    } 
    
    // parameters
    int nTheta, phiIndx, nuIndx;

    // NOTE: this differs from other non-factor modeling functions
    if (corName != "matern") {
      nTheta = 1; //phi
      phiIndx = 0;
    } else{
	nTheta = 2; //phi, nu
	phiIndx = 0; nuIndx = 1;
      }
   
    int nThetaq = nTheta * q; 
    // get max nu
    double *nuMax = (double *) R_alloc(q, sizeof(double));
    int *nb = (int *) R_alloc(q, sizeof(nb)); 
    // Fill in with zeros
    for (ll = 0; ll < q; ll++) {
      nuMax[ll] = 0.0; 
      nb[ll] = 0; 
    }
    
    if(corName == "matern"){
      for (ll = 0; ll < q; ll++) {
        for(s = 0; s < nSamples; s++){
          if(theta[s*nThetaq + nuIndx * q + ll] > nuMax[ll]){
            nuMax[ll] = theta[s*nThetaq + nuIndx * q + ll];
          }
        }
        nb[ll] = 1+static_cast<int>(floor(nuMax[ll]));
      }
    }

    int nbMax = 0; 
    for (ll = 0; ll < q; ll++) {
      if (nb[ll] > nbMax) {
        nbMax = nb[ll]; 
      }
    }

    double *bk = (double *) R_alloc(nThreads*nbMax, sizeof(double));
    double *C = (double *) R_alloc(nThreads*mm, sizeof(double)); zeros(C, nThreads*mm);
    double *c = (double *) R_alloc(nThreads*m, sizeof(double)); zeros(c, nThreads*m);
    double *tmp_m  = (double *) R_alloc(nThreads*m, sizeof(double));
    double phi = 0, nu = 0, sigmaSq = 1.0, d;
    int threadID = 0, status = 0;

    SEXP y0_r, w0_r, mu0_r;
    PROTECT(y0_r = allocMatrix(REALSXP, nObs0nSp, nSamples)); nProtect++; 
    PROTECT(mu0_r = allocMatrix(REALSXP, nObs0nSp, nSamples)); nProtect++; 
    PROTECT(w0_r = allocMatrix(REALSXP, J0q, nSamples)); nProtect++;
    double *y0 = REAL(y0_r);
    double *mu0 = REAL(mu0_r); 
    double *w0 = REAL(w0_r);
    double *w0Star = (double *) R_alloc(nSamples * J0nSp, sizeof(double));
    if (verbose) {
      Rprintf("-------------------------------------------------\n");
      Rprintf("\t\tPredicting\n");
      Rprintf("-------------------------------------------------\n");
      #ifdef Win32
        R_FlushConsole();
      #endif
    }

    int vIndx = -1;
    double *wV = (double *) R_alloc(J0q*nSamples, sizeof(double));

    GetRNGstate();
    
    for(i = 0; i < J0q*nSamples; i++){
      wV[i] = rnorm(0.0,1.0);
    }
    
    for(j = 0; j < J0; j++){
      for (ll = 0; ll < q; ll++) {
#ifdef _OPENMP
#pragma omp parallel for private(threadID, phi, nu, sigmaSq, k, l, d, info)
#endif     
        for(s = 0; s < nSamples; s++){
#ifdef _OPENMP
	  threadID = omp_get_thread_num();
#endif 	
	  if (sites0Sampled[j] == 1) {
            w0[s * J0q + j * q + ll] = w[s * Jq + sitesLink[j] * q + ll];
	  } else {
	    phi = theta[s * nThetaq + phiIndx * q + ll];
	    if(corName == "matern"){
	      nu = theta[s * nThetaq + nuIndx * q + ll];
	    }
	    sigmaSq = 1.0;

	    for(k = 0; k < m; k++){
	      d = dist2(coords[nnIndx0[j+J0*k]], coords[J+nnIndx0[j+J0*k]], coords0[j], coords0[J0+j]);
	      c[threadID*m+k] = sigmaSq*spCor(d, phi, nu, covModel, &bk[threadID*nb[ll]]);
	      for(l = 0; l < m; l++){
	        d = dist2(coords[nnIndx0[j+J0*k]], coords[J+nnIndx0[j+J0*k]], coords[nnIndx0[j+J0*l]], coords[J+nnIndx0[j+J0*l]]);
	        C[threadID*mm+l*m+k] = sigmaSq*spCor(d, phi, nu, covModel, &bk[threadID*nb[ll]]);
	      }
	    }

	    F77_NAME(dpotrf)(lower, &m, &C[threadID*mm], &m, &info FCONE); 
	    if(info != 0){error("c++ error: dpotrf failed\n");}
	    F77_NAME(dpotri)(lower, &m, &C[threadID*mm], &m, &info FCONE); 
	    if(info != 0){error("c++ error: dpotri failed\n");}

	    F77_NAME(dsymv)(lower, &m, &one, &C[threadID*mm], &m, &c[threadID*m], &inc, &zero, &tmp_m[threadID*m], &inc FCONE);

	    d = 0;
	    for(k = 0; k < m; k++){
	      d += tmp_m[threadID*m+k]*w[s*Jq+nnIndx0[j+J0*k] * q + ll];
	    }

	    #ifdef _OPENMP
            #pragma omp atomic
            #endif   
	    vIndx++;

	    w0[s * J0q + j * q + ll] = sqrt(sigmaSq - F77_NAME(ddot)(&m, &tmp_m[threadID*m], &inc, &c[threadID*m], &inc))*wV[vIndx] + d;
	  }

        } // sample
      } // factor

      if(verbose){
	if(status == nReport){
	  Rprintf("Location: %i of %i, %3.2f%%\n", j, J0, 100.0*j/J0);
          #ifdef Win32
	  R_FlushConsole();
          #endif
	  status = 0;
	}
      }
      status++;
      R_CheckUserInterrupt();
    } // location

    if(verbose){
      Rprintf("Location: %i of %i, %3.2f%%\n", j, J0, 100.0*j/J0);
      #ifdef Win32
      R_FlushConsole();
      #endif
    }

    // Generate abundance after the fact.
    if (verbose) {
      Rprintf("Generating abundance predictions\n");
    }
    for (j = 0; j < nObs0; j++) {
      for (s = 0; s < nSamples; s++) {
          F77_NAME(dgemv)(ntran, &nSp, &q, &one, &lambda[s * nSpq], &nSp, &w0[s * J0q + sites0[j]*q], &inc, &zero, &w0Star[s * J0nSp + sites0[j] * nSp], &inc FCONE);
	  for (i = 0; i < nSp; i++) {
	    mu0[s * nObs0nSp + j * nSp + i] = exp(F77_NAME(ddot)(&pAbund, &X0[j], &nObs0, &beta[s*pAbundnSp + i], &nSp) + w0Star[s * J0nSp + sites0[j] * nSp + i] + betaStarSite[s*nObs0nSp + j * nSp + i]);
	    if (family == 1) {
	      y0[s * nObs0nSp + j * nSp + i] = rnbinom_mu(kappa[s * nSp + i], mu0[s * nObs0nSp + j * nSp + i]);
	    } else {
	      y0[s * nObs0nSp + j * nSp + i] = rpois(mu0[s * nObs0nSp + j * nSp + i]);
	    }
	  } // i
      } // s
      R_CheckUserInterrupt();
    } // j

    PutRNGstate();
    
    //make return object
    SEXP result_r, resultName_r;
    int nResultListObjs = 3;

    PROTECT(result_r = allocVector(VECSXP, nResultListObjs)); nProtect++;
    PROTECT(resultName_r = allocVector(VECSXP, nResultListObjs)); nProtect++;

    SET_VECTOR_ELT(result_r, 0, y0_r);
    SET_VECTOR_ELT(resultName_r, 0, mkChar("y.0.samples")); 
    
    SET_VECTOR_ELT(result_r, 1, w0_r);
    SET_VECTOR_ELT(resultName_r, 1, mkChar("w.0.samples"));

    SET_VECTOR_ELT(result_r, 2, mu0_r);
    SET_VECTOR_ELT(resultName_r, 2, mkChar("mu.0.samples")); 

    namesgets(result_r, resultName_r);
    
    //unprotect
    UNPROTECT(nProtect);
    
    return(result_r);
  
  }
}

    
