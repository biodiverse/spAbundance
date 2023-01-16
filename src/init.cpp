#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include "spAbundance.h"

static const R_CallMethodDef CallEntries[] = {
    {"abund", (DL_FUNC) &abund, 30},
    {"DS", (DL_FUNC) &DS, 49},
    {"msAbund", (DL_FUNC) &msAbund, 33},
    {"NMix", (DL_FUNC) &NMix, 45},
    {"sfMsAbundNNGP", (DL_FUNC) &sfMsAbundNNGP, 49},
    {"spAbundNNGP", (DL_FUNC) &spAbundNNGP, 48},
    {"spAbundNNGPPredict", (DL_FUNC) &spAbundNNGPPredict, 24},
    {"spNMixNNGP", (DL_FUNC) &spNMixNNGP, 63},
    {"spNMixNNGPPredict", (DL_FUNC) &spNMixNNGPPredict, 19},
    {"waicAbund", (DL_FUNC) &waicAbund, 13},
    {NULL, NULL, 0}
};

void R_init_spAbundance(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
