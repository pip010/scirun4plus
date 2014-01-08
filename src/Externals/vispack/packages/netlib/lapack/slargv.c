#include "f2c.h"

/* Subroutine */ int slargv_(integer *n, real *x, integer *incx, real *y, 
	integer *incy, real *c, integer *incc)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   


    Purpose   
    =======   

    SLARGV generates a vector of real plane rotations, determined by   
    elements of the real vectors x and y. For i = 1,2,...,n   

       (  c(i)  s(i) ) ( x(i) ) = ( a(i) )   
       ( -s(i)  c(i) ) ( y(i) ) = (   0  )   

    Arguments   
    =========   

    N       (input) INTEGER   
            The number of plane rotations to be generated.   

    X       (input/output) REAL array,   
                           dimension (1+(N-1)*INCX)   
            On entry, the vector x.   
            On exit, x(i) is overwritten by a(i), for i = 1,...,n.   

    INCX    (input) INTEGER   
            The increment between elements of X. INCX > 0.   

    Y       (input/output) REAL array,   
                           dimension (1+(N-1)*INCY)   
            On entry, the vector y.   
            On exit, the sines of the plane rotations.   

    INCY    (input) INTEGER   
            The increment between elements of Y. INCY > 0.   

    C       (output) REAL array, dimension (1+(N-1)*INCC)   
            The cosines of the plane rotations.   

    INCC    (input) INTEGER   
            The increment between elements of C. INCC > 0.   

    ===================================================================== 
  


    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    integer i__1;
    real r__1, r__2;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    static integer i;
    static real w;
    static integer ic, ix, iy;
    static real xi, yi, tt;


#define C(I) c[(I)-1]
#define Y(I) y[(I)-1]
#define X(I) x[(I)-1]


    ix = 1;
    iy = 1;
    ic = 1;
    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
	xi = X(ix);
	yi = Y(iy);
	if (xi == 0.f) {
	    C(ic) = 0.f;
	    Y(iy) = 1.f;
	    X(ix) = yi;
	} else {
/* Computing MAX */
	    r__1 = dabs(xi), r__2 = dabs(yi);
	    w = dmax(r__1,r__2);
	    xi /= w;
	    yi /= w;
	    tt = sqrt(xi * xi + yi * yi);
	    C(ic) = xi / tt;
	    Y(iy) = yi / tt;
	    X(ix) = w * tt;
	}
	ix += *incx;
	iy += *incy;
	ic += *incc;
/* L10: */
    }
    return 0;

/*     End of SLARGV */

} /* slargv_ */

