#include "f2c.h"

/* Subroutine */ int dlapll_(integer *n, doublereal *x, integer *incx, 
	doublereal *y, integer *incy, doublereal *ssmin)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       March 31, 1993   


    Purpose   
    =======   

    Given two column vectors X and Y, let   

                         A = ( X Y ).   

    The subroutine first computes the QR factorization of A = Q*R,   
    and then computes the SVD of the 2-by-2 upper triangular matrix R.   
    The smaller singular value of R is returned in SSMIN, which is used   
    as the measurement of the linear dependency of the vectors X and Y.   

    Arguments   
    =========   

    N       (input) INTEGER   
            The length of the vectors X and Y.   

    X       (input/output) DOUBLE PRECISION array,   
                           dimension (1+(N-1)*INCX)   
            On entry, X contains the N-vector X.   
            On exit, X is overwritten.   

    INCX    (input) INTEGER   
            The increment between successive elements of X. INCX > 0.   

    Y       (input/output) DOUBLE PRECISION array,   
                           dimension (1+(N-1)*INCY)   
            On entry, Y contains the N-vector Y.   
            On exit, Y is overwritten.   

    INCY    (input) INTEGER   
            The increment between successive elements of Y. INCY > 0.   

    SSMIN   (output) DOUBLE PRECISION   
            The smallest singular value of the N-by-2 matrix A = ( X Y ). 
  

    ===================================================================== 
  


       Quick return if possible   

    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    integer i__1;
    /* Local variables */
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    extern /* Subroutine */ int dlas2_(doublereal *, doublereal *, doublereal 
	    *, doublereal *, doublereal *);
    static doublereal c;
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);
    static doublereal ssmax, a11, a12, a22;
    extern /* Subroutine */ int dlarfg_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *);
    static doublereal tau;


#define Y(I) y[(I)-1]
#define X(I) x[(I)-1]


    if (*n <= 1) {
	*ssmin = 0.;
	return 0;
    }

/*     Compute the QR factorization of the N-by-2 matrix ( X Y ) */

    dlarfg_(n, &X(1), &X(*incx + 1), incx, &tau);
    a11 = X(1);
    X(1) = 1.;

    c = -tau * ddot_(n, &X(1), incx, &Y(1), incy);
    daxpy_(n, &c, &X(1), incx, &Y(1), incy);

    i__1 = *n - 1;
    dlarfg_(&i__1, &Y(*incy + 1), &Y((*incy << 1) + 1), incy, &tau);

    a12 = Y(1);
    a22 = Y(*incy + 1);

/*     Compute the SVD of 2-by-2 Upper triangular matrix. */

    dlas2_(&a11, &a12, &a22, ssmin, &ssmax);

    return 0;

/*     End of DLAPLL */

} /* dlapll_ */

