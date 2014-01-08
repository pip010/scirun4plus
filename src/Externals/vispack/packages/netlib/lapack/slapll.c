#include "f2c.h"

/* Subroutine */ int slapll_(integer *n, real *x, integer *incx, real *y, 
	integer *incy, real *ssmin)
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

    X       (input/output) REAL array,   
                           dimension (1+(N-1)*INCX)   
            On entry, X contains the N-vector X.   
            On exit, X is overwritten.   

    INCX    (input) INTEGER   
            The increment between successive elements of X. INCX > 0.   

    Y       (input/output) REAL array,   
                           dimension (1+(N-1)*INCY)   
            On entry, Y contains the N-vector Y.   
            On exit, Y is overwritten.   

    INCY    (input) INTEGER   
            The increment between successive elements of Y. INCY > 0.   

    SSMIN   (output) REAL   
            The smallest singular value of the N-by-2 matrix A = ( X Y ). 
  

    ===================================================================== 
  


       Quick return if possible   

    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    integer i__1;
    /* Local variables */
    extern doublereal sdot_(integer *, real *, integer *, real *, integer *);
    extern /* Subroutine */ int slas2_(real *, real *, real *, real *, real *)
	    ;
    static real c, ssmax;
    extern /* Subroutine */ int saxpy_(integer *, real *, real *, integer *, 
	    real *, integer *);
    static real a11, a12, a22;
    extern /* Subroutine */ int slarfg_(integer *, real *, real *, integer *, 
	    real *);
    static real tau;


#define Y(I) y[(I)-1]
#define X(I) x[(I)-1]


    if (*n <= 1) {
	*ssmin = 0.f;
	return 0;
    }

/*     Compute the QR factorization of the N-by-2 matrix ( X Y ) */

    slarfg_(n, &X(1), &X(*incx + 1), incx, &tau);
    a11 = X(1);
    X(1) = 1.f;

    c = -(doublereal)tau * sdot_(n, &X(1), incx, &Y(1), incy);
    saxpy_(n, &c, &X(1), incx, &Y(1), incy);

    i__1 = *n - 1;
    slarfg_(&i__1, &Y(*incy + 1), &Y((*incy << 1) + 1), incy, &tau);

    a12 = Y(1);
    a22 = Y(*incy + 1);

/*     Compute the SVD of 2-by-2 Upper triangular matrix. */

    slas2_(&a11, &a12, &a22, ssmin, &ssmax);

    return 0;

/*     End of SLAPLL */

} /* slapll_ */

