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
  SEXP checkAlphaDS(SEXP y_r, SEXP Xp_r, SEXP XpRE_r,
                    SEXP XpRandom_r, SEXP yMax_r,
                    SEXP consts_r, SEXP K_r, SEXP nDetRELong_r,
                    SEXP alphaStarting_r,
                    SEXP sigmaSqPStarting_r,
                    SEXP alphaStarStarting_r, SEXP NStarting_r,
                    SEXP NLongIndx_r,
                    SEXP alphaStarIndx_r, SEXP alphaLevelIndx_r,
                    SEXP muAlpha_r, SEXP SigmaAlpha_r,
                    SEXP detModel_r,
                    SEXP transect_r, SEXP distBreaks_r){

    /**********************************************************************
     * Initial constants
     * *******************************************************************/
    int i, j, l, k, nProtect=0;
    const int inc = 1;

    /**********************************************************************
     * Get Inputs
     * *******************************************************************/
    double *y = REAL(y_r);
    double *Xp = REAL(Xp_r);
    double *yMax = REAL(yMax_r);
    int *XpRE = INTEGER(XpRE_r);
    double *XpRandom = REAL(XpRandom_r);
    // Load constants
    int J = INTEGER(consts_r)[0];
    int nObs = INTEGER(consts_r)[1];
    int pDet = INTEGER(consts_r)[5];
    int pDetRE = INTEGER(consts_r)[6];
    int nDetRE = INTEGER(consts_r)[7];
    int ppDet = pDet * pDet;
    double *muAlpha = (double *) R_alloc(pDet, sizeof(double));
    F77_NAME(dcopy)(&pDet, REAL(muAlpha_r), &inc, muAlpha, &inc);
    double *SigmaAlpha = (double *) R_alloc(ppDet, sizeof(double));
    F77_NAME(dcopy)(&ppDet, REAL(SigmaAlpha_r), &inc, SigmaAlpha, &inc);
    int K = INTEGER(K_r)[0];
    int *alphaStarIndx = INTEGER(alphaStarIndx_r);
    int *alphaLevelIndx = INTEGER(alphaLevelIndx_r);
    // DS stuff
    int detModel = INTEGER(detModel_r)[0];
    int transect = INTEGER(transect_r)[0];
    double *distBreaks = REAL(distBreaks_r);

    /**********************************************************************
     * Parameters
     * *******************************************************************/
    // Detection covariates
    double *alpha = (double *) R_alloc(pDet, sizeof(double));
    F77_NAME(dcopy)(&pDet, REAL(alphaStarting_r), &inc, alpha, &inc);
    // Detection random effect variances
    double *sigmaSqP = (double *) R_alloc(pDetRE, sizeof(double));
    F77_NAME(dcopy)(&pDetRE, REAL(sigmaSqPStarting_r), &inc, sigmaSqP, &inc);
    // Latent detection random effects
    double *alphaStar = (double *) R_alloc(nDetRE, sizeof(double));
    F77_NAME(dcopy)(&nDetRE, REAL(alphaStarStarting_r), &inc, alphaStar, &inc);
    // Latent Abundance
    double *N = (double *) R_alloc(J, sizeof(double));
    F77_NAME(dcopy)(&J, REAL(NStarting_r), &inc, N, &inc);

    /**********************************************************************
     * Return Stuff
     * *******************************************************************/
    SEXP alphaSuccess_r;
    PROTECT(alphaSuccess_r = Rf_allocMatrix(REALSXP, inc, inc)); nProtect++;

    /********************************************************************
      Some constants and temporary variables to be used later
    ********************************************************************/
    int JpDetRE = J * pDetRE;
    double tmp_0;

    /********************************************************************
      Set up MH stuff
    ********************************************************************/
    double logPostAlphaCurr = 0.0;

    /**********************************************************************
     * Prep for random effects
     * *******************************************************************/
    // Observation-level sums of the detection random effects
    double *alphaStarSites = (double *) R_alloc(J, sizeof(double));
    zeros(alphaStarSites, J);
    double *alphaStarSitesCand = (double *) R_alloc(J, sizeof(double));
    int *alphaStarLongIndx = (int *) R_alloc(JpDetRE, sizeof(int));
    // Get sums of the current REs for each site/visit combo
    for (j = 0; j < J; j++) {
      for (l = 0; l < pDetRE; l++) {
        alphaStarLongIndx[l * J + j] = which(XpRE[l * J + j], alphaLevelIndx, nDetRE);
        alphaStarSites[j] += alphaStar[alphaStarLongIndx[l * J + j]] *
                           XpRandom[l * J + j];
      }
      alphaStarSitesCand[j] = alphaStarSites[j];
    }
    // Starting index for detection random effects
    int *alphaStarStart = (int *) R_alloc(pDetRE, sizeof(int));
    for (l = 0; l < pDetRE; l++) {
      alphaStarStart[l] = which(l, alphaStarIndx, nDetRE);
    }

    /**********************************************************************
     * DS Prep
     * *******************************************************************/
    double *sigma = (double *) R_alloc(J, sizeof(double));
    double *p = (double *) R_alloc(nObs, sizeof(double)); zeros(p, nObs);
    // Number of break points for integration
    int nInt = 5;
    double *binWidth = (double *) R_alloc(K, sizeof(double));
    double stripWidth = 0.0;
    double *psi = (double *) R_alloc(K, sizeof(double));
    for (k = 0; k < K; k++) {
      binWidth[k] = distBreaks[k + 1] - distBreaks[k];
      stripWidth += binWidth[k];
    }
    for (k = 0; k < K; k++) {
      if (transect == 0) {
        psi[k] = binWidth[k] / stripWidth;
      } else {
        psi[k] = (pow(distBreaks[k + 1], 2) - pow(distBreaks[k], 2)) / pow(stripWidth, 2);
      }
    }
    int KFull = K + 1;
    int nObsFull = KFull * J;
    double *piFull = (double *) R_alloc(nObsFull, sizeof(double));
    zeros(piFull, nObsFull);
    double likeVal = 0.0;

    GetRNGstate();

    /********************************************************************
     *Update Detection Regression Coefficients
     *******************************************************************/
    for (l = 0; l < pDet; l++) {
      logPostAlphaCurr = 0.0;
      for (i = 0; i < pDet; i++ ) {
        logPostAlphaCurr += dnorm(alpha[i], muAlpha[i], sqrt(SigmaAlpha[i * pDet + i]), 1);
      }
      for (j = 0; j < J; j++) {
        /********************************
         * Current
         *******************************/
        tmp_0 = 0.0;
        likeVal = 0.0;
        sigma[j] = exp(F77_NAME(ddot)(&pDet, &Xp[j], &J, alpha, &inc) +
                       alphaStarSites[j]);
        for (k = 0; k < K; k++) {
          p[k * J + j] = integrate(detModel, distBreaks[k], distBreaks[k + 1], sigma[j],
                                   nInt, transect);
          if (transect == 0) {
            p[k * J + j] /= (distBreaks[k + 1] - distBreaks[k]);
          } else {
            p[k * J + j] = p[k * J + j] * 2.0 / (pow(distBreaks[k + 1], 2) - pow(distBreaks[k], 2));
          }
          piFull[k * J + j] = p[k * J + j] * psi[k];
          tmp_0 += piFull[k * J + j];
          likeVal += y[k * J + j] * log(piFull[k * J + j]);
        } // k (bins)
        piFull[K * J + j] = 1.0 - tmp_0;
        likeVal += (N[j] - yMax[j]) * log(piFull[K * J + j]);
        logPostAlphaCurr += likeVal;
      } // j (sites)
    }
    REAL(alphaSuccess_r)[0] = logPostAlphaCurr;
    PutRNGstate();

    SEXP result_r, resultName_r;
    int nResultListObjs = 1;

    PROTECT(result_r = Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;
    PROTECT(resultName_r = Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;

    // Setting the components of the output list.
    SET_VECTOR_ELT(result_r, 0, alphaSuccess_r);
    SET_VECTOR_ELT(resultName_r, 0, Rf_mkChar("alpha.like.val"));

    Rf_namesgets(result_r, resultName_r);

    UNPROTECT(nProtect);

    return(result_r);
  }
}

