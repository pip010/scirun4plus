#include "f2c.h"

/* Subroutine */ int dlagtf_(integer *n, doublereal *a, doublereal *lambda, 
	doublereal *b, doublereal *c, doublereal *tol, doublereal *d, integer 
	*in, integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    DLAGTF factorizes the matrix (T - lambda*I), where T is an n by n   
    tridiagonal matrix and lambda is a scalar, as   

       T - lambda*I = PLU,   

    where P is a permutation matrix, L is a unit lower tridiagonal matrix 
  
    with at most one non-zero sub-diagonal elements per column and U is   
    an upper triangular matrix with at most two non-zero super-diagonal   
    elements per column.   

    The factorization is obtained by Gaussian elimination with partial   
    pivoting and implicit row scaling.   

    The parameter LAMBDA is included in the routine so that DLAGTF may   
    be used, in conjunction with DLAGTS, to obtain eigenvectors of T by   
    inverse iteration.   

    Arguments   
    =========   

    N       (input) INTEGER   
            The order of the matrix T.   

    A       (input/output) DOUBLE PRECISION array, dimension (N)   
            On entry, A must contain the diagonal elements of T.   

            On exit, A is overwritten by the n diagonal elements of the   
            upper triangular matrix U of the factorization of T.   

    LAMBDA  (input) DOUBLE PRECISION   
            On entry, the scalar lambda.   

    B       (input/output) DOUBLE PRECISION array, dimension (N-1)   
            On entry, B must contain the (n-1) super-diagonal elements of 
  
            T.   

            On exit, B is overwritten by the (n-1) super-diagonal   
            elements of the matrix U of the factorization of T.   

    C       (input/output) DOUBLE PRECISION array, dimension (N-1)   
            On entry, C must contain the (n-1) sub-diagonal elements of   
            T.   

            On exit, C is overwritten by the (n-1) sub-diagonal elements 
  
            of the matrix L of the factorization of T.   

    TOL     (input) DOUBLE PRECISION   
            On entry, a relative tolerance used to indicate whether or   
            not the matrix (T - lambda*I) is nearly singular. TOL should 
  
            normally be chose as approximately the largest relative error 
  
            in the elements of T. For example, if the elements of T are   
            correct to about 4 significant figures, then TOL should be   
            set to about 5*10**(-4). If TOL is supplied as less than eps, 
  
            where eps is the relative machine precision, then the value   
            eps is used in place of TOL.   

    D       (output) DOUBLE PRECISION array, dimension (N-2)   
            On exit, D is overwritten by the (n-2) second super-diagonal 
  
            elements of the matrix U of the factorization of T.   

    IN      (output) INTEGER array, dimension (N)   
            On exit, IN contains details of the permutation matrix P. If 
  
            an interchange occurred at the kth step of the elimination,   
            then IN(k) = 1, otherwise IN(k) = 0. The element IN(n)   
            returns the smallest positive integer j such that   

               abs( u(j,j) ).le. norm( (T - lambda*I)(j) )*TOL,   

            where norm( A(j) ) denotes the sum of the absolute values of 
  
            the jth row of the matrix A. If no such j exists then IN(n)   
            is returned as zero. If IN(n) is returned as positive, then a 
  
            diagonal element of U is small, indicating that   
            (T - lambda*I) is singular or nearly singular,   

    INFO    (output)   
            = 0   : successful exit   
            .lt. 0: if INFO = -k, the kth argument had an illegal value   

   ===================================================================== 
  


    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;
    /* Local variables */
    static doublereal temp, mult;
    static integer k;
    static doublereal scale1, scale2;
    extern doublereal dlamch_(char *);
    static doublereal tl;
    extern /* Subroutine */ int xerbla_(char *, integer *);
    static doublereal eps, piv1, piv2;


#define IN(I) in[(I)-1]
#define D(I) d[(I)-1]
#define C(I) c[(I)-1]
#define B(I) b[(I)-1]
#define A(I) a[(I)-1]


    *info = 0;
    if (*n < 0) {
	*info = -1;
	i__1 = -(*info);
	xerbla_("DLAGTF", &i__1);
	return 0;
    }

    if (*n == 0) {
	return 0;
    }

    A(1) -= *lambda;
    IN(*n) = 0;
    if (*n == 1) {
	if (A(1) == 0.) {
	    IN(1) = 1;
	}
	return 0;
    }

    eps = dlamch_("Epsilon");

    tl = max(*tol,eps);
    scale1 = abs(A(1)) + abs(B(1));
    i__1 = *n - 1;
    for (k = 1; k <= *n-1; ++k) {
	A(k + 1) -= *lambda;
	scale2 = (d__1 = C(k), abs(d__1)) + (d__2 = A(k + 1), abs(d__2));
	if (k < *n - 1) {
	    scale2 += (d__1 = B(k + 1), abs(d__1));
	}
	if (A(k) == 0.) {
	    piv1 = 0.;
	} else {
	    piv1 = (d__1 = A(k), abs(d__1)) / scale1;
	}
	if (C(k) == 0.) {
	    IN(k) = 0;
	    piv2 = 0.;
	    scale1 = scale2;
	    if (k < *n - 1) {
		D(k) = 0.;
	    }
	} else {
	    piv2 = (d__1 = C(k), abs(d__1)) / scale2;
	    if (piv2 <= piv1) {
		IN(k) = 0;
		scale1 = scale2;
		C(k) /= A(k);
		A(k + 1) -= C(k) * B(k);
		if (k < *n - 1) {
		    D(k) = 0.;
		}
	    } else {
		IN(k) = 1;
		mult = A(k) / C(k);
		A(k) = C(k);
		temp = A(k + 1);
		A(k + 1) = B(k) - mult * temp;
		if (k < *n - 1) {
		    D(k) = B(k + 1);
		    B(k + 1) = -mult * D(k);
		}
		B(k) = temp;
		C(k) = mult;
	    }
	}
	if (max(piv1,piv2) <= tl && IN(*n) == 0) {
	    IN(*n) = k;
	}
/* L10: */
    }
    if ((d__1 = A(*n), abs(d__1)) <= scale1 * tl && IN(*n) == 0) {
	IN(*n) = *n;
    }

    return 0;

/*     End of DLAGTF */

} /* dlagtf_ */

