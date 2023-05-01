#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include "spAbundance.h"

static const R_CallMethodDef CallEntries[] = {
    {"abund", (DL_FUNC) &abund, 30},
    {"spAbundNNGP", (DL_FUNC) &spAbundNNGP, 48},
    {"spAbundNNGPPredict", (DL_FUNC) &spAbundNNGPPredict, 24},
    {"msAbund", (DL_FUNC) &msAbund, 33},
    {"lfMsAbund", (DL_FUNC) &lfMsAbund, 35},
    {"sfMsAbundNNGP", (DL_FUNC) &sfMsAbundNNGP, 49},
    {"sfMsAbundNNGPPredict", (DL_FUNC) &sfMsAbundNNGPPredict, 27},
    {"NMix", (DL_FUNC) &NMix, 45},
    {"spNMixNNGP", (DL_FUNC) &spNMixNNGP, 63},
    {"spNMixNNGPPredict", (DL_FUNC) &spNMixNNGPPredict, 23},
    {"msNMix", (DL_FUNC) &msNMix, 52},
    {"lfMsNMix", (DL_FUNC) &lfMsNMix, 54},
    {"sfMsNMixNNGP", (DL_FUNC) &sfMsNMixNNGP, 65},
    {"sfMsNMixNNGPPredict", (DL_FUNC) &sfMsNMixNNGPPredict, 22},
    {"DS", (DL_FUNC) &DS, 49},
    {"spDSNNGP", (DL_FUNC) &spDSNNGP, 63},
    {"waicAbund", (DL_FUNC) &waicAbund, 13},
    {"svcAbundNNGP", (DL_FUNC) &svcAbundNNGP, 46},
    {"abundGaussian", (DL_FUNC) &abundGaussian, 26},
    {"svcAbundNNGPPredict", (DL_FUNC) &svcAbundNNGPPredict, 25},
    {"spAbundGaussianNNGP", (DL_FUNC) &spAbundGaussianNNGP, 46},
    {NULL, NULL, 0}
};

void R_init_spAbundance(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
