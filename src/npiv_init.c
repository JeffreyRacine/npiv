#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include "mgcv.h"

/* .C calls */
extern int gsl_bspline(double *, int *, int *, int *, double *, double *, double *, int *, double *);
extern int gsl_bspline_deriv(double *, int *, int *, int *, int *, int *, double *, double *, double *, int *, double *);

static const R_CMethodDef CEntries[] = {
    {"gsl_bspline",       (DL_FUNC) &gsl_bspline,        9},
    {"gsl_bspline_deriv", (DL_FUNC) &gsl_bspline_deriv, 11},
    {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
    {"glp_model_tmm",     (DL_FUNC) &glp_model_tmm,     7}, 
    {"mgcv_tmm",          (DL_FUNC) &mgcv_tmm,          5}, 
    {NULL, NULL, 0}
};

void R_init_npiv(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
