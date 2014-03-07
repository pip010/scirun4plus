#include "f2c.h"

/* Subroutine */ int zlartv_(integer *n, doublecomplex *x, integer *incx, 
	doublecomplex *y, integer *incy, doublereal *c, doublecomplex *s, 
	integer *incc)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   


    Purpose   
    =======   

    ZLARTV applies a vector of complex plane rotations with real cosines 
  
    to elements of the complex vectors x and y. For i = 1,2,...,n   

       ( x(i) ) := (        c(i)   s(i) ) ( x(i) )   
       ( y(i) )    ( -conjg(s(i))  c(i) ) ( y(i) )   

    Arguments   
    =========   

    N       (input) INTEGER   
            The number of plane rotations to be applied.   

    X       (input/output) COMPLEX*16 array, dimension (1+(N-1)*INCX)   
            The vector x.   

    INCX    (input) INTEGER   
            The increment between elements of X. INCX > 0.   

    Y       (input/output) COMPLEX*16 array, dimension (1+(N-1)*INCY)   
            The vector y.   

    INCY    (input) INTEGER   
            The increment between elements of Y. INCY > 0.   

    C       (input) DOUBLE PRECISION array, dimension (1+(N-1)*INCC)   
            The cosines of the plane rotations.   

    S       (input) COMPLEX*16 array, dimension (1+(N-1)*INCC)   
            The sines of the plane rotations.   

    INCC    (input) INTEGER   
            The increment between elements of C and S. INCC > 0.   

    ===================================================================== 
  


    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    doublecomplex z__1, z__2, z__3, z__4;
    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);
    /* Local variables */
    static integer i, ic, ix, iy;
    static doublecomplex xi, yi;


#define S(I) s[(I)-1]
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
	i__2 = ix;
	i__3 = ic;
	z__2.r = C(ic) * xi.r, z__2.i = C(ic) * xi.i;
	i__4 = ic;
	z__3.r = S(ic).r * yi.r - S(ic).i * yi.i, z__3.i = S(ic).r * 
		yi.i + S(ic).i * yi.r;
	z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
	X(ix).r = z__1.r, X(ix).i = z__1.i;
	i__2 = iy;
	i__3 = ic;
	z__2.r = C(ic) * yi.r, z__2.i = C(ic) * yi.i;
	d_cnjg(&z__4, &S(ic));
	z__3.r = z__4.r * xi.r - z__4.i * xi.i, z__3.i = z__4.r * xi.i + 
		z__4.i * xi.r;
	z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
	Y(iy).r = z__1.r, Y(iy).i = z__1.i;
	ix += *incx;
	iy += *incy;
	ic += *incc;
/* L10: */
    }
    return 0;

/*     End of ZLARTV */

} /* zlartv_ */

