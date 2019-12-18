#ifndef R_R_H
# include <R.h>
#endif

#ifndef R_EXT_DYNLOAD_H_
# include <R_ext/Rdynload.h>
#endif

#include <Rinternals.h>
#include <stdlib.h> // for NULL

/* register native routines ------------------------------------------------ */

/* .Fortran calls */

void F77_NAME(inithiv)(void *);
void F77_NAME(derivshiv)(int *, double *, double *, double *, double *, int *);    
 
R_FortranMethodDef FEntries[] = {
    {"inithiv",       (DL_FUNC) &F77_SUB(inithiv),      1},
    {"derivshiv",     (DL_FUNC) &F77_SUB(derivshiv),    6},
    {NULL, NULL, 0}
};

/* Initialization ---------------------------------------------------------- */
void R_init_FME(DllInfo *dll) {

  R_registerRoutines(dll, NULL, NULL, FEntries, NULL);

  // the following two lines protect against accidentially finding entry points

  R_useDynamicSymbols(dll, FALSE); // disable dynamic searching
}
