#include "f2c.h"

/* Subroutine */ int dpttrf_(integer *n, doublereal *d, doublereal *e, 
	integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       March 31, 1993   


    Purpose   
    =======   

    DPTTRF computes the factorization of a real symmetric positive   
    definite tridiagonal matrix A.   

    If the subdiagonal elements of A are supplied in the array E, the   
    factorization has the form A = L*D*L**T, where D is diagonal and L   
    is unit lower bidiagonal; if the superdiagonal elements of A are   
    supplied, it has the form A = U**T*D*U, where U is unit upper   
    bidiagonal.  (The two forms are equivalent if A is real.)   

    Arguments   
    =========   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    D       (input/output) DOUBLE PRECISION array, dimension (N)   
            On entry, the n diagonal elements of the tridiagonal matrix   
            A.  On exit, the n diagonal elements of the diagonal matrix   
            D from the L*D*L**T factorization of A.   

    E       (input/output) DOUBLE PRECISION array, dimension (N-1)   
            On entry, the (n-1) off-diagonal elements of the tridiagonal 
  
            matrix A.   
            On exit, the (n-1) off-diagonal elements of the unit   
            bidiagonal factor L or U from the factorization of A.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   
            > 0:  if INFO = i, the leading minor of order i is not   
                  positive definite; if i < N, the factorization could   
                  not be completed, while if i = N, the factorization was 
  
                  completed, but D(N) = 0.   

    ===================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    integer i__1;
    /* Local variables */
    static integer i;
    static doublereal di, ei;
    extern /* Subroutine */ int xerbla_(char *, integer *);


#define E(I) e[(I)-1]
#define D(I) d[(I)-1]


    *info = 0;
    if (*n < 0) {
	*info = -1;
	i__1 = -(*info);
	xerbla_("DPTTRF", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return 0;
    }

/*     Compute the L*D*L' (or U'*D*U) factorization of A. */

    i__1 = *n - 1;
    for (i = 1; i <= *n-1; ++i) {

/*        Drop out of the loop if d(i) <= 0: the matrix is not positiv
e   
          definite. */

	di = D(i);
	if (di <= 0.) {
	    goto L20;
	}

/*        Solve for e(i) and d(i+1). */

	ei = E(i);
	E(i) = ei / di;
	D(i + 1) -= E(i) * ei;
/* L10: */
    }

/*     Check d(n) for positive definiteness. */

    i = *n;
    if (D(i) > 0.) {
	goto L30;
    }

L20:
    *info = i;

L30:
    return 0;

/*     End of DPTTRF */

} /* dpttrf_ */

