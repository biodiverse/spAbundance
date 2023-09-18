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
  SEXP checkMSAlphaDS(SEXP y_r, SEXP Xp_r,  SEXP XpRE_r, 
                      SEXP XpRandom_r, SEXP yMax_r, 
                      SEXP consts_r, SEXP K_r, 
                      SEXP alphaStarting_r, 
                      SEXP alphaCommStarting_r,
	              SEXP tauSqAlphaStarting_r,
                      SEXP sigmaSqPStarting_r,
                      SEXP alphaStarStarting_r, SEXP NStarting_r,
                      SEXP alphaStarIndx_r, SEXP alphaLevelIndx_r,
                      SEXP muAlphaComm_r, SEXP SigmaAlphaComm_r, 
	              SEXP detModel_r, 
	              SEXP transect_r, SEXP distBreaks_r){
   
    /**********************************************************************
     * Initial constants
     * *******************************************************************/
    int i, j, r, l, k, info, nProtect=0;
    const int inc = 1;
    const double one = 1.0;
    const double zero = 0.0;
    char const *lower = "L";
    
    /**********************************************************************
     * Get Inputs
     * *******************************************************************/
    double *y = REAL(y_r);
    double *Xp = REAL(Xp_r);
    double *yMax = REAL(yMax_r);
    int *XpRE = INTEGER(XpRE_r);
    double *XpRandom = REAL(XpRandom_r);
    // Load constants
    int nSp = INTEGER(consts_r)[0]; 
    int J = INTEGER(consts_r)[1];
    int nObs = INTEGER(consts_r)[2];
    int pDet = INTEGER(consts_r)[6];
    int pDetRE = INTEGER(consts_r)[7];
    int nDetRE = INTEGER(consts_r)[8];
    int ppDet = pDet * pDet;
    double *muAlphaComm = REAL(muAlphaComm_r); 
    double *SigmaAlphaCommInv = (double *) R_alloc(ppDet, sizeof(double));   
    F77_NAME(dcopy)(&ppDet, REAL(SigmaAlphaComm_r), &inc, SigmaAlphaCommInv, &inc);
    int K = INTEGER(K_r)[0]; 
    int *alphaStarIndx = INTEGER(alphaStarIndx_r); 
    int *alphaLevelIndx = INTEGER(alphaLevelIndx_r);
    // DS stuff
    int detModel = INTEGER(detModel_r)[0];
    int transect = INTEGER(transect_r)[0];
    double *distBreaks = REAL(distBreaks_r);

    /********************************************************************
      Some constants and temporary variables to be used later
    ********************************************************************/
    int pDetnSp = pDet * nSp; 
    int nDetREnSp = nDetRE * nSp; 
    int JpDetRE = J * pDetRE;
    int JnSp = J * nSp;
    int KFull = K + 1;
    int nObsFull = KFull * J;
    int nObsFullnSp = nObsFull * nSp;
    double tmp_0; 
   
    /**********************************************************************
     * Parameters
     * *******************************************************************/
    double *alphaComm = (double *) R_alloc(pDet, sizeof(double));   
    F77_NAME(dcopy)(&pDet, REAL(alphaCommStarting_r), &inc, alphaComm, &inc);
    double *tauSqAlpha = (double *) R_alloc(pDet, sizeof(double)); 
    F77_NAME(dcopy)(&pDet, REAL(tauSqAlphaStarting_r), &inc, tauSqAlpha, &inc);
    double *alpha = (double *) R_alloc(pDetnSp, sizeof(double));   
    F77_NAME(dcopy)(&pDetnSp, REAL(alphaStarting_r), &inc, alpha, &inc);
    // Detection random effect variances
    double *sigmaSqP = (double *) R_alloc(pDetRE, sizeof(double)); 
    F77_NAME(dcopy)(&pDetRE, REAL(sigmaSqPStarting_r), &inc, sigmaSqP, &inc); 
    // Latent detection random effects
    double *alphaStar = (double *) R_alloc(nDetREnSp, sizeof(double)); 
    F77_NAME(dcopy)(&nDetREnSp, REAL(alphaStarStarting_r), &inc, alphaStar, &inc); 
    // Latent Abundance
    double *N = (double *) R_alloc(JnSp, sizeof(double));   
    F77_NAME(dcopy)(&JnSp, REAL(NStarting_r), &inc, N, &inc);

    /**********************************************************************
     * Return Stuff
     * *******************************************************************/
    SEXP alphaSuccess_r;
    PROTECT(alphaSuccess_r = allocMatrix(REALSXP, nSp, inc)); nProtect++;

    /**********************************************************************
     * Additional Sampler Prep
     * *******************************************************************/
    // Detection regression coefficient priors. 
    F77_NAME(dpotrf)(lower, &pDet, SigmaAlphaCommInv, &pDet, &info FCONE); 
    if(info != 0){error("c++ error: dpotrf SigmaAlphaCommInv failed\n");}
    F77_NAME(dpotri)(lower, &pDet, SigmaAlphaCommInv, &pDet, &info FCONE); 
    if(info != 0){error("c++ error: dpotri SigmaAlphaCommInv failed\n");}
    double *SigmaAlphaCommInvMuAlpha = (double *) R_alloc(pDet, sizeof(double)); 
    F77_NAME(dsymv)(lower, &pDet, &one, SigmaAlphaCommInv, &pDet, muAlphaComm, &inc, &zero, 
                   SigmaAlphaCommInvMuAlpha, &inc FCONE);
    // Put community level variances in a pDet x pDet matrix. 
    double *TauAlphaInv = (double *) R_alloc(ppDet, sizeof(double)); zeros(TauAlphaInv, ppDet); 
    for (i = 0; i < pDet; i++) {
      TauAlphaInv[i * pDet + i] = tauSqAlpha[i]; 
    } // i
    F77_NAME(dpotrf)(lower, &pDet, TauAlphaInv, &pDet, &info FCONE); 
    if(info != 0){error("c++ error: dpotrf TauAlphaInv failed\n");}
    F77_NAME(dpotri)(lower, &pDet, TauAlphaInv, &pDet, &info FCONE); 
    if(info != 0){error("c++ error: dpotri TauAlphaInv failed\n");}

    /********************************************************************
      Set up MH stuff
    ********************************************************************/
    double logPostAlphaCurr = 0.0;

    /**********************************************************************
     * Prep for random effects
     * *******************************************************************/
    // Observation-level sums of the detection random effects
    double *alphaStarSites = (double *) R_alloc(JnSp, sizeof(double)); 
    double *alphaStarSitesCand = (double *) R_alloc(JnSp, sizeof(double)); 
    zeros(alphaStarSites, JnSp); 
    int *alphaStarLongIndx = (int *) R_alloc(JpDetRE, sizeof(int));
    // Get sums of the current REs for each site/visit combo for all species
    for (r = 0; r < J; r++) {
      for (l = 0; l < pDetRE; l++) {
        alphaStarLongIndx[l * J + r] = which(XpRE[l * J + r], alphaLevelIndx, nDetRE);
        for (i = 0; i < nSp; i++) {
          alphaStarSites[i * J + r] += alphaStar[i * nDetRE + alphaStarLongIndx[l * J + r]] * XpRandom[l * J + r];
	  alphaStarSitesCand[i * J + r] = alphaStarSites[i * J + r];
        }
      }
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
    double *piFull = (double *) R_alloc(nObsFullnSp, sizeof(double)); 
    zeros(piFull, nObsFullnSp);
    double likeVal = 0.0;

    GetRNGstate(); 

    /**********************************************************************
     *Species specific variables  
     *********************************************************************/
    for (i = 0; i < nSp; i++) { 
      /********************************************************************
       *Update Detection Regression Coefficients
       *******************************************************************/
      // alpha is ordered by parameter, then species
      // Proposal
      for (l = 0; l < pDet; l++) {
        logPostAlphaCurr = 0.0;
        logPostAlphaCurr += dnorm(alpha[l * nSp + i], alphaComm[l], 
                                  sqrt(tauSqAlpha[l]), 1);
        for (j = 0; j < J; j++) {
          /********************************
           * Current 
           *******************************/
          tmp_0 = 0.0; 
          likeVal = 0.0;
          sigma[j] = exp(F77_NAME(ddot)(&pDet, &Xp[j], &J, &alpha[i], &nSp) + 
                         alphaStarSites[i * J + j]);
          for (k = 0; k < K; k++) {
            p[k * J + j] = integrate(detModel, distBreaks[k], distBreaks[k + 1], sigma[j], 
                                     nInt, transect); 
            if (transect == 0) {
              p[k * J + j] /= (distBreaks[k + 1] - distBreaks[k]);
            } else {
              p[k * J + j] = p[k * J + j] * 2.0 / (pow(distBreaks[k + 1], 2) - pow(distBreaks[k], 2));
            }
            piFull[k * JnSp + j * nSp + i] = p[k * J + j] * psi[k];
            tmp_0 += piFull[k * JnSp + j * nSp + i];
            likeVal += y[k * JnSp + j * nSp + i] * log(piFull[k * JnSp + j * nSp + i]);
          } // k (bins)
          piFull[k * JnSp + j * nSp + i] = 1.0 - tmp_0;
          likeVal += (N[j * nSp + i] - yMax[j * nSp + i]) * log(piFull[K * JnSp + j * nSp + i]);
          logPostAlphaCurr += likeVal;
        } // j (sites)
      }
      REAL(alphaSuccess_r)[i] = logPostAlphaCurr;
    }
    PutRNGstate();

    SEXP result_r, resultName_r;
    int nResultListObjs = 1;

    PROTECT(result_r = allocVector(VECSXP, nResultListObjs)); nProtect++;
    PROTECT(resultName_r = allocVector(VECSXP, nResultListObjs)); nProtect++;

    // Setting the components of the output list.
    SET_VECTOR_ELT(result_r, 0, alphaSuccess_r);

    SET_VECTOR_ELT(resultName_r, 0, mkChar("alpha.like.val")); 
   
    namesgets(result_r, resultName_r);
    
    UNPROTECT(nProtect);
    
    return(result_r);
  }
}

