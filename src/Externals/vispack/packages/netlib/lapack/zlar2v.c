#include "f2c.h"

/* Subroutine */ int zlar2v_(integer *n, doublecomplex *x, doublecomplex *y, 
	doublecomplex *z, integer *incx, doublereal *c, doublecomplex *s, 
	integer *incc)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   


    Purpose   
    =======   

    ZLAR2V applies a vector of complex plane rotations with real cosines 
  
    from both sides to a sequence of 2-by-2 complex Hermitian matrices,   
    defined by the elements of the vectors x, y and z. For i = 1,2,...,n 
  

       (       x(i)  z(i) ) :=   
       ( conjg(z(i)) y(i) )   

         (  c(i) conjg(s(i)) ) (       x(i)  z(i) ) ( c(i) -conjg(s(i)) ) 
  
         ( -s(i)       c(i)  ) ( conjg(z(i)) y(i) ) ( s(i)        c(i)  ) 
  

    Arguments   
    =========   

    N       (input) INTEGER   
            The number of plane rotations to be applied.   

    X       (input/output) COMPLEX*16 array, dimension (1+(N-1)*INCX)   
            The vector x; the elements of x are assumed to be real.   

    Y       (input/output) COMPLEX*16 array, dimension (1+(N-1)*INCX)   
            The vector y; the elements of y are assumed to be real.   

    Z       (input/output) COMPLEX*16 array, dimension (1+(N-1)*INCX)   
            The vector z.   

    INCX    (input) INTEGER   
            The increment between elements of X, Y and Z. INCX > 0.   

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
    integer i__1, i__2;
    doublereal d__1;
    doublecomplex z__1, z__2, z__3, z__4, z__5;
    /* Builtin functions */
    double d_imag(doublecomplex *);
    void d_cnjg(doublecomplex *, doublecomplex *);
    /* Local variables */
    static integer i;
    static doublecomplex t2, t3, t4;
    static doublereal t5, t6;
    static integer ic;
    static doublereal ci;
    static doublecomplex si;
    static integer ix;
    static doublereal xi, yi;
    static doublecomplex zi;
    static doublereal t1i, t1r, sii, zii, sir, zir;


#define S(I) s[(I)-1]
#define C(I) c[(I)-1]
#define Z(I) z[(I)-1]
#define Y(I) y[(I)-1]
#define X(I) x[(I)-1]


    ix = 1;
    ic = 1;
    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
	i__2 = ix;
	xi = X(ix).r;
	i__2 = ix;
	yi = Y(ix).r;
	i__2 = ix;
	zi.r = Z(ix).r, zi.i = Z(ix).i;
	zir = zi.r;
	zii = d_imag(&zi);
	ci = C(ic);
	i__2 = ic;
	si.r = S(ic).r, si.i = S(ic).i;
	sir = si.r;
	sii = d_imag(&si);
	t1r = sir * zir - sii * zii;
	t1i = sir * zii + sii * zir;
	z__1.r = ci * zi.r, z__1.i = ci * zi.i;
	t2.r = z__1.r, t2.i = z__1.i;
	d_cnjg(&z__3, &si);
	z__2.r = xi * z__3.r, z__2.i = xi * z__3.i;
	z__1.r = t2.r - z__2.r, z__1.i = t2.i - z__2.i;
	t3.r = z__1.r, t3.i = z__1.i;
	d_cnjg(&z__2, &t2);
	z__3.r = yi * si.r, z__3.i = yi * si.i;
	z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
	t4.r = z__1.r, t4.i = z__1.i;
	t5 = ci * xi + t1r;
	t6 = ci * yi - t1r;
	i__2 = ix;
	d__1 = ci * t5 + (sir * t4.r + sii * d_imag(&t4));
	X(ix).r = d__1, X(ix).i = 0.;
	i__2 = ix;
	d__1 = ci * t6 - (sir * t3.r - sii * d_imag(&t3));
	Y(ix).r = d__1, Y(ix).i = 0.;
	i__2 = ix;
	z__2.r = ci * t3.r, z__2.i = ci * t3.i;
	d_cnjg(&z__4, &si);
	z__5.r = t6, z__5.i = t1i;
	z__3.r = z__4.r * z__5.r - z__4.i * z__5.i, z__3.i = z__4.r * z__5.i 
		+ z__4.i * z__5.r;
	z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
	Z(ix).r = z__1.r, Z(ix).i = z__1.i;
	ix += *incx;
	ic += *incc;
/* L10: */
    }
    return 0;

/*     End of ZLAR2V */

} /* zlar2v_ */

