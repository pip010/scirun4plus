#include "f2c.h"

/* Subroutine */ int clargv_(integer *n, complex *x, integer *incx, complex *
	y, integer *incy, real *c, integer *incc)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   


    Purpose   
    =======   

    CLARGV generates a vector of complex plane rotations with real   
    cosines, determined by elements of the complex vectors x and y.   
    For i = 1,2,...,n   

       (        c(i)   s(i) ) ( x(i) ) = ( a(i) )   
       ( -conjg(s(i))  c(i) ) ( y(i) ) = (   0  )   

    Arguments   
    =========   

    N       (input) INTEGER   
            The number of plane rotations to be generated.   

    X       (input/output) COMPLEX array, dimension (1+(N-1)*INCX)   
            On entry, the vector x.   
            On exit, x(i) is overwritten by a(i), for i = 1,...,n.   

    INCX    (input) INTEGER   
            The increment between elements of X. INCX > 0.   

    Y       (input/output) COMPLEX array, dimension (1+(N-1)*INCY)   
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
    integer i__1, i__2;
    doublereal d__1;
    complex q__1, q__2, q__3;
    /* Builtin functions */
    double c_abs(complex *), sqrt(doublereal);
    void r_cnjg(complex *, complex *);
    /* Local variables */
    static real absx, absy;
    static integer i;
    static complex t;
    static real w;
    static integer ic, ix, iy;
    static complex xi, yi;
    static real tt;


#define C(I) c[(I)-1]
#define Y(I) y[(I)-1]
#define X(I) x[(I)-1]


    ix = 1;
    iy = 1;
    ic = 1;
    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
	i__2 = ix;
	xi.r = X(ix).r, xi.i = X(ix).i;
	i__2 = iy;
	yi.r = Y(iy).r, yi.i = Y(iy).i;
	absx = c_abs(&xi);
	if (absx == 0.f) {
	    C(ic) = 0.f;
	    i__2 = iy;
	    Y(iy).r = 1.f, Y(iy).i = 0.f;
	    i__2 = ix;
	    X(ix).r = yi.r, X(ix).i = yi.i;
	} else {
	    absy = c_abs(&yi);
	    w = dmax(absx,absy);
	    q__1.r = xi.r / absx, q__1.i = xi.i / absx;
	    t.r = q__1.r, t.i = q__1.i;
	    absx /= w;
	    absy /= w;
	    tt = sqrt(absx * absx + absy * absy);
	    C(ic) = absx / tt;
	    i__2 = iy;
	    r_cnjg(&q__3, &yi);
	    q__2.r = t.r * q__3.r - t.i * q__3.i, q__2.i = t.r * q__3.i + t.i 
		    * q__3.r;
	    d__1 = w * tt;
	    q__1.r = q__2.r / d__1, q__1.i = q__2.i / d__1;
	    Y(iy).r = q__1.r, Y(iy).i = q__1.i;
	    i__2 = ix;
	    d__1 = w * tt;
	    q__1.r = d__1 * t.r, q__1.i = d__1 * t.i;
	    X(ix).r = q__1.r, X(ix).i = q__1.i;
	}
	ix += *incx;
	iy += *incy;
	ic += *incc;
/* L10: */
    }
    return 0;

/*     End of CLARGV */

} /* clargv_ */

