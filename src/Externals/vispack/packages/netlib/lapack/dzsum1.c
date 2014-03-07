#include "f2c.h"

doublereal dzsum1_(integer *n, doublecomplex *cx, integer *incx)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    DZSUM1 takes the sum of the absolute values of a complex   
    vector and returns a double precision result.   

    Based on DZASUM from the Level 1 BLAS.   
    The change is to use the 'genuine' absolute value.   

    Contributed by Nick Higham for use with ZLACON.   

    Arguments   
    =========   

    N       (input) INTEGER   
            The number of elements in the vector CX.   

    CX      (input) COMPLEX*16 array, dimension (N)   
            The vector whose elements will be summed.   

    INCX    (input) INTEGER   
            The spacing between successive values of CX.  INCX > 0.   

    ===================================================================== 
  


    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    integer i__1, i__2;
    doublereal ret_val;
    /* Builtin functions */
    double z_abs(doublecomplex *);
    /* Local variables */
    static integer i, nincx;
    static doublereal stemp;


#define CX(I) cx[(I)-1]


    ret_val = 0.;
    stemp = 0.;
    if (*n <= 0) {
	return ret_val;
    }
    if (*incx == 1) {
	goto L20;
    }

/*     CODE FOR INCREMENT NOT EQUAL TO 1 */

    nincx = *n * *incx;
    i__1 = nincx;
    i__2 = *incx;
    for (i = 1; *incx < 0 ? i >= nincx : i <= nincx; i += *incx) {

/*        NEXT LINE MODIFIED. */

	stemp += z_abs(&CX(i));
/* L10: */
    }
    ret_val = stemp;
    return ret_val;

/*     CODE FOR INCREMENT EQUAL TO 1 */

L20:
    i__2 = *n;
    for (i = 1; i <= *n; ++i) {

/*        NEXT LINE MODIFIED. */

	stemp += z_abs(&CX(i));
/* L30: */
    }
    ret_val = stemp;
    return ret_val;

/*     End of DZSUM1 */

} /* dzsum1_ */

