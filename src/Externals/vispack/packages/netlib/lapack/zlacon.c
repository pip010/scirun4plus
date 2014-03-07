#include "f2c.h"

/* Subroutine */ int zlacon_(integer *n, doublecomplex *v, doublecomplex *x, 
	doublereal *est, integer *kase)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    ZLACON estimates the 1-norm of a square, complex matrix A.   
    Reverse communication is used for evaluating matrix-vector products. 
  

    Arguments   
    =========   

    N      (input) INTEGER   
           The order of the matrix.  N >= 1.   

    V      (workspace) COMPLEX*16 array, dimension (N)   
           On the final return, V = A*W,  where  EST = norm(V)/norm(W)   
           (W is not returned).   

    X      (input/output) COMPLEX*16 array, dimension (N)   
           On an intermediate return, X should be overwritten by   
                 A * X,   if KASE=1,   
                 A' * X,  if KASE=2,   
           where A' is the conjugate transpose of A, and ZLACON must be   
           re-called with all the other parameters unchanged.   

    EST    (output) DOUBLE PRECISION   
           An estimate (a lower bound) for norm(A).   

    KASE   (input/output) INTEGER   
           On the initial call to ZLACON, KASE should be 0.   
           On an intermediate return, KASE will be 1 or 2, indicating   
           whether X should be overwritten by A * X  or A' * X.   
           On the final return from ZLACON, KASE will again be 0.   

    Further Details   
    ======= =======   

    Contributed by Nick Higham, University of Manchester.   
    Originally named CONEST, dated March 16, 1988.   

    Reference: N.J. Higham, "FORTRAN codes for estimating the one-norm of 
  
    a real or complex matrix, with applications to condition estimation", 
  
    ACM Trans. Math. Soft., vol. 14, no. 4, pp. 381-396, December 1988.   

    ===================================================================== 
  


    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static integer c__1 = 1;
    
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;
    doublecomplex z__1, z__2;
    /* Builtin functions */
    double z_abs(doublecomplex *);
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *);
    /* Local variables */
    static integer iter;
    static doublereal temp;
    static integer jump, i, j, jlast;
    extern /* Subroutine */ int zcopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    extern integer izmax1_(integer *, doublecomplex *, integer *);
    extern doublereal dzsum1_(integer *, doublecomplex *, integer *), dlamch_(
	    char *);
    static doublereal safmin, altsgn, estold;



#define X(I) x[(I)-1]
#define V(I) v[(I)-1]


    safmin = dlamch_("Safe minimum");
    if (*kase == 0) {
	i__1 = *n;
	for (i = 1; i <= *n; ++i) {
	    i__2 = i;
	    d__1 = 1. / (doublereal) (*n);
	    z__1.r = d__1, z__1.i = 0.;
	    X(i).r = z__1.r, X(i).i = z__1.i;
/* L10: */
	}
	*kase = 1;
	jump = 1;
	return 0;
    }

    switch (jump) {
	case 1:  goto L20;
	case 2:  goto L40;
	case 3:  goto L70;
	case 4:  goto L90;
	case 5:  goto L120;
    }

/*     ................ ENTRY   (JUMP = 1)   
       FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY A*X. */

L20:
    if (*n == 1) {
	V(1).r = X(1).r, V(1).i = X(1).i;
	*est = z_abs(&V(1));
/*        ... QUIT */
	goto L130;
    }
    *est = dzsum1_(n, &X(1), &c__1);

    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
	if (z_abs(&X(i)) > safmin) {
	    i__2 = i;
	    d__1 = z_abs(&X(i));
	    z__2.r = d__1, z__2.i = 0.;
	    z_div(&z__1, &X(i), &z__2);
	    X(i).r = z__1.r, X(i).i = z__1.i;
	} else {
	    i__2 = i;
	    X(i).r = 1., X(i).i = 0.;
	}
/* L30: */
    }
    *kase = 2;
    jump = 2;
    return 0;

/*     ................ ENTRY   (JUMP = 2)   
       FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY ZTRANS(A)*X. */

L40:
    j = izmax1_(n, &X(1), &c__1);
    iter = 2;

/*     MAIN LOOP - ITERATIONS 2,3,...,ITMAX. */

L50:
    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
	i__2 = i;
	X(i).r = 0., X(i).i = 0.;
/* L60: */
    }
    i__1 = j;
    X(j).r = 1., X(j).i = 0.;
    *kase = 1;
    jump = 3;
    return 0;

/*     ................ ENTRY   (JUMP = 3)   
       X HAS BEEN OVERWRITTEN BY A*X. */

L70:
    zcopy_(n, &X(1), &c__1, &V(1), &c__1);
    estold = *est;
    *est = dzsum1_(n, &V(1), &c__1);

/*     TEST FOR CYCLING. */
    if (*est <= estold) {
	goto L100;
    }

    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
	if (z_abs(&X(i)) > safmin) {
	    i__2 = i;
	    d__1 = z_abs(&X(i));
	    z__2.r = d__1, z__2.i = 0.;
	    z_div(&z__1, &X(i), &z__2);
	    X(i).r = z__1.r, X(i).i = z__1.i;
	} else {
	    i__2 = i;
	    X(i).r = 1., X(i).i = 0.;
	}
/* L80: */
    }
    *kase = 2;
    jump = 4;
    return 0;

/*     ................ ENTRY   (JUMP = 4)   
       X HAS BEEN OVERWRITTEN BY ZTRANS(A)*X. */

L90:
    jlast = j;
    j = izmax1_(n, &X(1), &c__1);
    i__1 = jlast;
    i__2 = j;
    if (X(jlast).r != (d__1 = X(j).r, abs(d__1)) && iter < 5) {
	++iter;
	goto L50;
    }

/*     ITERATION COMPLETE.  FINAL STAGE. */

L100:
    altsgn = 1.;
    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
	i__2 = i;
	d__1 = altsgn * ((doublereal) (i - 1) / (doublereal) (*n - 1) + 1.);
	z__1.r = d__1, z__1.i = 0.;
	X(i).r = z__1.r, X(i).i = z__1.i;
	altsgn = -altsgn;
/* L110: */
    }
    *kase = 1;
    jump = 5;
    return 0;

/*     ................ ENTRY   (JUMP = 5)   
       X HAS BEEN OVERWRITTEN BY A*X. */

L120:
    temp = dzsum1_(n, &X(1), &c__1) / (doublereal) (*n * 3) * 2.;
    if (temp > *est) {
	zcopy_(n, &X(1), &c__1, &V(1), &c__1);
	*est = temp;
    }

L130:
    *kase = 0;
    return 0;

/*     End of ZLACON */

} /* zlacon_ */

