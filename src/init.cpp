#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include "spAbundance.h"

static const R_CallMethodDef CallEntries[] = {
    {"NMix", (DL_FUNC) &NMix, 42},
    {"spNMixNNGP", (DL_FUNC) &spNMixNNGP, 59},
    {"msNMix", (DL_FUNC) &msNMix, 49},
    {"sfMsNMixNNGP", (DL_FUNC) &sfMsNMixNNGP, 65},
    {"abund", (DL_FUNC) &abund, 28},
    {"spAbundNNGP", (DL_FUNC) &spAbundNNGP, 46},
    {"sfMsAbundNNGP", (DL_FUNC) &sfMsAbundNNGP, 47},
    {"waicAbund", (DL_FUNC) &waicAbund, 13},
    {NULL, NULL, 0}
};

void R_init_spAbundance(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
