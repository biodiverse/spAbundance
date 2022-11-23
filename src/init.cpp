#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include "spAbundance.h"

static const R_CallMethodDef CallEntries[] = {
    {"NMix", (DL_FUNC) &NMix, 42},
    {"spNMixNNGP", (DL_FUNC) &spNMixNNGP, 59},
    {"abund", (DL_FUNC) &abund, 30},
    {"spAbundNNGP", (DL_FUNC) &spAbundNNGP, 48},
    {"sfMsAbundNNGP", (DL_FUNC) &sfMsAbundNNGP, 49},
    {"waicAbund", (DL_FUNC) &waicAbund, 13},
    {NULL, NULL, 0}
};

void R_init_spAbundance(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
