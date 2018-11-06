#include <R.h>
#include <Rinternals.h>
#include <R_ext/Lapack.h>

// wrapper of DLAED9

SEXP dpr1_eigen(SEXP KSTART, SEXP KSTOP, SEXP RD, SEXP RZ, SEXP RRHO){
  int N;
  SEXP RQ, RS, RDLAM, RINFO, Rout;
  
  N = length(RD);

  RQ = PROTECT(allocMatrix(REALSXP, N, N));
  RS = PROTECT(allocMatrix(REALSXP, N, N));
  RDLAM = PROTECT(allocVector(REALSXP, INTEGER(KSTOP)[0]-INTEGER(KSTART)[0]+1));
  RINFO = PROTECT(allocVector(INTSXP, 1));

  dlaed9_(INTEGER(KSTOP), INTEGER(KSTART), INTEGER(KSTOP), &N, REAL(RDLAM), REAL(RQ),
      &N, REAL(RRHO), REAL(RD), REAL(RZ), REAL(RS), &N, INTEGER(RINFO));

  Rout = PROTECT(allocVector(VECSXP, 3));
  SET_VECTOR_ELT(Rout, 0, RS);
  SET_VECTOR_ELT(Rout, 1, RDLAM);
  SET_VECTOR_ELT(Rout, 2, RINFO);

  UNPROTECT(5);

  return Rout;
}
