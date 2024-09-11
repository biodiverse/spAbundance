#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include "spAbundance.h"

static const R_CallMethodDef CallEntries[] = {
    {"abund", (DL_FUNC) &abund, 30},
    {"abundGaussian", (DL_FUNC) &abundGaussian, 26},
    {"spAbundNNGP", (DL_FUNC) &spAbundNNGP, 48},
    {"spAbundGaussianNNGP", (DL_FUNC) &spAbundGaussianNNGP, 46},
    {"svcAbundNNGPPredict", (DL_FUNC) &svcAbundNNGPPredict, 26},
    {"msAbund", (DL_FUNC) &msAbund, 34},
    {"msAbundGaussian", (DL_FUNC) &msAbundGaussian, 33},
    {"lfMsAbund", (DL_FUNC) &lfMsAbund, 36},
    {"lfMsAbundGaussian", (DL_FUNC) &lfMsAbundGaussian, 35},
    {"sfMsAbundNNGP", (DL_FUNC) &sfMsAbundNNGP, 50},
    {"sfMsAbundGaussianNNGP", (DL_FUNC) &sfMsAbundGaussianNNGP, 49},
    {"sfMsAbundNNGPPredict", (DL_FUNC) &sfMsAbundNNGPPredict, 27},
    {"NMix", (DL_FUNC) &NMix, 46},
    {"spNMixNNGP", (DL_FUNC) &spNMixNNGP, 64},
    {"spNMixNNGPPredict", (DL_FUNC) &spNMixNNGPPredict, 22},
    {"msNMix", (DL_FUNC) &msNMix, 53},
    {"lfMsNMix", (DL_FUNC) &lfMsNMix, 55},
    {"sfMsNMixNNGP", (DL_FUNC) &sfMsNMixNNGP, 65},
    {"sfMsNMixNNGPPredict", (DL_FUNC) &sfMsNMixNNGPPredict, 24},
    {"DS", (DL_FUNC) &DS, 49},
    {"spDSNNGP", (DL_FUNC) &spDSNNGP, 63},
    {"msDS", (DL_FUNC) &msDS, 57},
    {"lfMsDS", (DL_FUNC) &lfMsDS, 59},
    {"sfMsDSNNGP", (DL_FUNC) &sfMsDSNNGP, 65},
    {"waicAbund", (DL_FUNC) &waicAbund, 13},
    {"checkAlphaDS", (DL_FUNC) &checkAlphaDS, 20},
    {"checkMSAlphaDS", (DL_FUNC) &checkMSAlphaDS, 20},
    {"svcAbundNNGP", (DL_FUNC) &svcAbundNNGP, 48},
    {"svcAbundGaussianNNGP", (DL_FUNC) &svcAbundNNGP, 46},
    {"svcAbundGaussianNNGPPredict", (DL_FUNC) &svcAbundGaussianNNGPPredict, 25},
    {"svcMsAbundGaussianNNGP", (DL_FUNC) &svcMsAbundGaussianNNGP, 50},
    {"svcMsAbundGaussianNNGPPredict", (DL_FUNC) &svcMsAbundGaussianNNGPPredict, 27},
    {NULL, NULL, 0}
};

void R_init_spAbundance(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
