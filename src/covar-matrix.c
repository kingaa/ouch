#include "ouch.h"

static void ouch_covar_matrix (int *nchar, 
			       double *bt, 
			       double *lambda, 
			       double *S, 
			       double *sigmasq, 
			       int *nterm, 
			       double *V) {
  double *U, *W;
  double sij, ti, tj;
  double *elti, *eltj;
  double tmp;
  int n = *nchar, nt = *nterm;
  int i, j, k, l, r, s;
  U = R_Calloc(n*n,double);
  W = R_Calloc(n*n,double);
  elti = R_Calloc(n,double);
  eltj = R_Calloc(n,double);
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      U[i+j*n] = 0;
      for (k = 0; k < n; k++) {
	for (l = 0; l < n; l++) {
	  U[i+j*n] += S[k+i*n]*sigmasq[k+l*n]*S[l+j*n];
	}
      }
    }
  }
  for (i = 0; i < nt; i++) {
    for (j = 0; j <= i; j++) {
      ti = bt[i+i*nt];
      sij = bt[i+j*nt];
      tj = bt[j+j*nt];
      for (k = 0; k < n; k++) {
	elti[k] = exp(-lambda[k]*(ti-sij));
	eltj[k] = exp(-lambda[k]*(tj-sij));
      }
      for (k = 0; k < n; k++) {
	for (l = 0; l < n; l++) {
	  V[i+nt*(k+n*(j+nt*l))] = 0; 
	  V[j+nt*(k+n*(i+nt*l))] = 0; 
	  W[k+l*n] = elti[k]*U[k+l*n]*eltj[l]/(lambda[k]+lambda[l]);
	}
      }
      for (k = 0; k < n; k++) {
	for (l = 0; l < n; l++) {
	  for (r = 0; r < n; r++) {
	    for (s = 0; s < n; s++) {
	      tmp = S[k+r*n]*W[r+s*n]*S[l+s*n];
	      V[i+nt*(k+n*(j+nt*l))] += tmp;
	      if (j != i) 
		V[j+nt*(l+n*(i+nt*k))] += tmp;
	    }
	  }
	}
      }
    }
  }
  R_Free(U);
  R_Free(W);
  R_Free(elti);
  R_Free(eltj);
}

SEXP ouch_covar (SEXP object, SEXP lambda, SEXP S, SEXP sigmasq) {
  int nprotect = 0;
  SEXP V, nm, nterm, bt;
  int nchar, nt, vdim[2];
  PROTECT(nm = NEW_CHARACTER(1)); nprotect++;
  SET_STRING_ELT(nm,0,mkChar("nterm"));
  PROTECT(nterm = GET_SLOT(object,nm)); nprotect++;
  nt = INTEGER(nterm)[0];
  SET_STRING_ELT(nm,0,mkChar("branch.times"));
  PROTECT(bt = GET_SLOT(object,nm)); nprotect++;
  nchar = GET_LENGTH(lambda);
  vdim[0] = nt*nchar; vdim[1] = vdim[0];
  PROTECT(V = makearray(2,vdim)); nprotect++;
  ouch_covar_matrix(&nchar,REAL(bt),REAL(lambda),REAL(S),REAL(sigmasq),&nt,REAL(V));
  UNPROTECT(nprotect);
  return V;
}


