#include "f2c.h"

/* Subroutine */ int zrot_(integer *n, doublecomplex *cx, integer *incx, 
	doublecomplex *cy, integer *incy, doublereal *c, doublecomplex *s)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    ZROT   applies a plane rotation, where the cos (C) is real and the   
    sin (S) is complex, and the vectors CX and CY are complex.   

    Arguments   
    =========   

    N       (input) INTEGER   
            The number of elements in the vectors CX and CY.   

    CX      (input/output) COMPLEX*16 array, dimension (N)   
            On input, the vector X.   
            On output, CX is overwritten with C*X + S*Y.   

    INCX    (input) INTEGER   
            The increment between successive values of CY.  INCX <> 0.   

    CY      (input/output) COMPLEX*16 array, dimension (N)   
            On input, the vector Y.   
            On output, CY is overwritten with -CONJG(S)*X + C*Y.   

    INCY    (input) INTEGER   
            The increment between successive values of CY.  INCX <> 0.   

    C       (input) DOUBLE PRECISION   
    S       (input) COMPLEX*16   
            C and S define a rotation   
               [  C          S  ]   
               [ -conjg(S)   C  ]   
            where C*C + S*CONJG(S) = 1.0.   

   ===================================================================== 
  


    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    doublecomplex z__1, z__2, z__3, z__4;
    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);
    /* Local variables */
    static integer i;
    static doublecomplex stemp;
    static integer ix, iy;


#define CY(I) cy[(I)-1]
#define CX(I) cx[(I)-1]


    if (*n <= 0) {
	return 0;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*     Code for unequal increments or equal increments not equal to 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
	i__2 = ix;
	z__2.r = *c * CX(ix).r, z__2.i = *c * CX(ix).i;
	i__3 = iy;
	z__3.r = s->r * CY(iy).r - s->i * CY(iy).i, z__3.i = s->r * CY(
		iy).i + s->i * CY(iy).r;
	z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
	stemp.r = z__1.r, stemp.i = z__1.i;
	i__2 = iy;
	i__3 = iy;
	z__2.r = *c * CY(iy).r, z__2.i = *c * CY(iy).i;
	d_cnjg(&z__4, s);
	i__4 = ix;
	z__3.r = z__4.r * CX(ix).r - z__4.i * CX(ix).i, z__3.i = z__4.r * 
		CX(ix).i + z__4.i * CX(ix).r;
	z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
	CY(iy).r = z__1.r, CY(iy).i = z__1.i;
	i__2 = ix;
	CX(ix).r = stemp.r, CX(ix).i = stemp.i;
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    return 0;

/*     Code for both increments equal to 1 */

L20:
    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
	i__2 = i;
	z__2.r = *c * CX(i).r, z__2.i = *c * CX(i).i;
	i__3 = i;
	z__3.r = s->r * CY(i).r - s->i * CY(i).i, z__3.i = s->r * CY(
		i).i + s->i * CY(i).r;
	z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
	stemp.r = z__1.r, stemp.i = z__1.i;
	i__2 = i;
	i__3 = i;
	z__2.r = *c * CY(i).r, z__2.i = *c * CY(i).i;
	d_cnjg(&z__4, s);
	i__4 = i;
	z__3.r = z__4.r * CX(i).r - z__4.i * CX(i).i, z__3.i = z__4.r * 
		CX(i).i + z__4.i * CX(i).r;
	z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
	CY(i).r = z__1.r, CY(i).i = z__1.i;
	i__2 = i;
	CX(i).r = stemp.r, CX(i).i = stemp.i;
/* L30: */
    }
    return 0;
} /* zrot_ */

