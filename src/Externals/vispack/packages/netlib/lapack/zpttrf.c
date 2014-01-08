#include "f2c.h"

/* Subroutine */ int zpttrf_(integer *n, doublereal *d, doublecomplex *e, 
	integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       March 31, 1993   


    Purpose   
    =======   

    ZPTTRF computes the factorization of a complex Hermitian positive   
    definite tridiagonal matrix A.   

    If the subdiagonal elements of A are supplied in the array E, the   
    factorization has the form A = L*D*L**H, where D is diagonal and L   
    is unit lower bidiagonal; if the superdiagonal elements of A are   
    supplied, it has the form A = U**H*D*U, where U is unit upper   
    bidiagonal.   

    Arguments   
    =========   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    D       (input/output) DOUBLE PRECISION array, dimension (N)   
            On entry, the n diagonal elements of the tridiagonal matrix   
            A.  On exit, the n diagonal elements of the diagonal matrix   
            D from the L*D*L**H factorization of A.   

    E       (input/output) COMPLEX*16 array, dimension (N-1)   
            On entry, the (n-1) off-diagonal elements of the tridiagonal 
  
            matrix A.  On exit, the (n-1) off-diagonal elements of the   
            unit bidiagonal factor L or U from the factorization of A.   

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
    integer i__1, i__2;
    doublecomplex z__1;
    /* Builtin functions */
    double d_imag(doublecomplex *);
    /* Local variables */
    static doublereal f, g;
    static integer i;
    static doublereal di;
    extern /* Subroutine */ int xerbla_(char *, integer *);
    static doublereal eii, eir;


#define E(I) e[(I)-1]
#define D(I) d[(I)-1]


    *info = 0;
    if (*n < 0) {
	*info = -1;
	i__1 = -(*info);
	xerbla_("ZPTTRF", &i__1);
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

	i__2 = i;
	eir = E(i).r;
	eii = d_imag(&E(i));
	f = eir / di;
	g = eii / di;
	i__2 = i;
	z__1.r = f, z__1.i = g;
	E(i).r = z__1.r, E(i).i = z__1.i;
	D(i + 1) = D(i + 1) - f * eir - g * eii;
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

/*     End of ZPTTRF */

} /* zpttrf_ */

