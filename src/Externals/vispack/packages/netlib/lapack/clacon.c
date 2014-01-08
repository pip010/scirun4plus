#include "f2c.h"

/* Subroutine */ int clacon_(integer *n, complex *v, complex *x, real *est, 
	integer *kase)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    CLACON estimates the 1-norm of a square, complex matrix A.   
    Reverse communication is used for evaluating matrix-vector products. 
  

    Arguments   
    =========   

    N      (input) INTEGER   
           The order of the matrix.  N >= 1.   

    V      (workspace) COMPLEX array, dimension (N)   
           On the final return, V = A*W,  where  EST = norm(V)/norm(W)   
           (W is not returned).   

    X      (input/output) COMPLEX array, dimension (N)   
           On an intermediate return, X should be overwritten by   
                 A * X,   if KASE=1,   
                 A' * X,  if KASE=2,   
           where A' is the conjugate transpose of A, and CLACON must be   
           re-called with all the other parameters unchanged.   

    EST    (output) REAL   
           An estimate (a lower bound) for norm(A).   

    KASE   (input/output) INTEGER   
           On the initial call to CLACON, KASE should be 0.   
           On an intermediate return, KASE will be 1 or 2, indicating   
           whether X should be overwritten by A * X  or A' * X.   
           On the final return from CLACON, KASE will again be 0.   

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
    real r__1;
    doublereal d__1;
    complex q__1, q__2;
    /* Builtin functions */
    double c_abs(complex *);
    void c_div(complex *, complex *, complex *);
    /* Local variables */
    static integer iter;
    static real temp;
    static integer jump, i, j, jlast;
    extern /* Subroutine */ int ccopy_(integer *, complex *, integer *, 
	    complex *, integer *);
    extern integer icmax1_(integer *, complex *, integer *);
    extern doublereal scsum1_(integer *, complex *, integer *), slamch_(char *
	    );
    static real safmin, altsgn, estold;



#define X(I) x[(I)-1]
#define V(I) v[(I)-1]


    safmin = slamch_("Safe minimum");
    if (*kase == 0) {
	i__1 = *n;
	for (i = 1; i <= *n; ++i) {
	    i__2 = i;
	    d__1 = 1.f / (real) (*n);
	    q__1.r = d__1, q__1.i = 0.f;
	    X(i).r = q__1.r, X(i).i = q__1.i;
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
	*est = c_abs(&V(1));
/*        ... QUIT */
	goto L130;
    }
    *est = scsum1_(n, &X(1), &c__1);

    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
	if (c_abs(&X(i)) > safmin) {
	    i__2 = i;
	    d__1 = c_abs(&X(i));
	    q__2.r = d__1, q__2.i = 0.f;
	    c_div(&q__1, &X(i), &q__2);
	    X(i).r = q__1.r, X(i).i = q__1.i;
	} else {
	    i__2 = i;
	    X(i).r = 1.f, X(i).i = 0.f;
	}
/* L30: */
    }
    *kase = 2;
    jump = 2;
    return 0;

/*     ................ ENTRY   (JUMP = 2)   
       FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY CTRANS(A)*X. */

L40:
    j = icmax1_(n, &X(1), &c__1);
    iter = 2;

/*     MAIN LOOP - ITERATIONS 2,3,...,ITMAX. */

L50:
    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
	i__2 = i;
	X(i).r = 0.f, X(i).i = 0.f;
/* L60: */
    }
    i__1 = j;
    X(j).r = 1.f, X(j).i = 0.f;
    *kase = 1;
    jump = 3;
    return 0;

/*     ................ ENTRY   (JUMP = 3)   
       X HAS BEEN OVERWRITTEN BY A*X. */

L70:
    ccopy_(n, &X(1), &c__1, &V(1), &c__1);
    estold = *est;
    *est = scsum1_(n, &V(1), &c__1);

/*     TEST FOR CYCLING. */
    if (*est <= estold) {
	goto L100;
    }

    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
	if (c_abs(&X(i)) > safmin) {
	    i__2 = i;
	    d__1 = c_abs(&X(i));
	    q__2.r = d__1, q__2.i = 0.f;
	    c_div(&q__1, &X(i), &q__2);
	    X(i).r = q__1.r, X(i).i = q__1.i;
	} else {
	    i__2 = i;
	    X(i).r = 1.f, X(i).i = 0.f;
	}
/* L80: */
    }
    *kase = 2;
    jump = 4;
    return 0;

/*     ................ ENTRY   (JUMP = 4)   
       X HAS BEEN OVERWRITTEN BY CTRANS(A)*X. */

L90:
    jlast = j;
    j = icmax1_(n, &X(1), &c__1);
    i__1 = jlast;
    i__2 = j;
    if (X(jlast).r != (r__1 = X(j).r, dabs(r__1)) && iter < 5) {
	++iter;
	goto L50;
    }

/*     ITERATION COMPLETE.  FINAL STAGE. */

L100:
    altsgn = 1.f;
    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
	i__2 = i;
	d__1 = altsgn * ((real) (i - 1) / (real) (*n - 1) + 1.f);
	q__1.r = d__1, q__1.i = 0.f;
	X(i).r = q__1.r, X(i).i = q__1.i;
	altsgn = -(doublereal)altsgn;
/* L110: */
    }
    *kase = 1;
    jump = 5;
    return 0;

/*     ................ ENTRY   (JUMP = 5)   
       X HAS BEEN OVERWRITTEN BY A*X. */

L120:
    temp = scsum1_(n, &X(1), &c__1) / (real) (*n * 3) * 2.f;
    if (temp > *est) {
	ccopy_(n, &X(1), &c__1, &V(1), &c__1);
	*est = temp;
    }

L130:
    *kase = 0;
    return 0;

/*     End of CLACON */

} /* clacon_ */

