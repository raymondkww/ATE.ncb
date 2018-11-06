#include <R.h>
#include <Rinternals.h>
#include <R_ext/Lapack.h>

// wrapper of DLAED4

SEXP dpr1_eigval(SEXP RI, SEXP RD, SEXP RZ, SEXP RRHO){
  int N;
  SEXP RDELTA, RDLAM, RINFO, Rout;
  
  N = length(RD);

  RDELTA = PROTECT(allocVector(REALSXP, N));
  RDLAM = PROTECT(allocVector(REALSXP, 1));
  RINFO = PROTECT(allocVector(INTSXP, 1));

  dlaed4_(&N, INTEGER(RI), REAL(RD), REAL(RZ), REAL(RDELTA), REAL(RRHO),
      REAL(RDLAM), INTEGER(RINFO));

  Rout = PROTECT(allocVector(VECSXP, 3));
  SET_VECTOR_ELT(Rout, 0, RDELTA);
  SET_VECTOR_ELT(Rout, 1, RDLAM);
  SET_VECTOR_ELT(Rout, 2, RINFO);

  UNPROTECT(4);

  return Rout;
}
