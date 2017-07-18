#include <R.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>

extern SEXP ouch_weights (SEXP object, SEXP lambda, SEXP S, SEXP beta);
extern SEXP ouch_covar (SEXP object, SEXP lambda, SEXP S, SEXP sigmasq);

static const R_CallMethodDef callMethods[] = {
  {"ouch_weights", (DL_FUNC) &ouch_weights, 4},
  {"ouch_covar", (DL_FUNC) &ouch_covar, 4},
  {NULL, NULL, 0}
};

void R_init_ouch (DllInfo *info) {
  // Register routines
  R_registerRoutines(info,NULL,callMethods,NULL,NULL);
  R_useDynamicSymbols(info,TRUE);
}
