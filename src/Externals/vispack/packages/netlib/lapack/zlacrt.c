#include "f2c.h"

/* Subroutine */ int zlacrt_(integer *n, doublecomplex *cx, integer *incx, 
	doublecomplex *cy, integer *incy, doublecomplex *c, doublecomplex *s)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    ZLACRT applies a plane rotation, where the cos and sin (C and S) are 
  
    complex and the vectors CX and CY are complex.   

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
            On output, CY is overwritten with -S*X + C*Y.   

    INCY    (input) INTEGER   
            The increment between successive values of CY.  INCX <> 0.   

    C       (input) COMPLEX*16   
    S       (input) COMPLEX*16   
            C and S define a complex rotation   
               [  C   S  ]   
               [ -S   C  ]   
            where C*C + S*S = 1.0.   

   ===================================================================== 
  


    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    doublecomplex z__1, z__2, z__3;
    /* Local variables */
    static integer i;
    static doublecomplex ctemp;
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
	z__2.r = c->r * CX(ix).r - c->i * CX(ix).i, z__2.i = c->r * CX(
		ix).i + c->i * CX(ix).r;
	i__3 = iy;
	z__3.r = s->r * CY(iy).r - s->i * CY(iy).i, z__3.i = s->r * CY(
		iy).i + s->i * CY(iy).r;
	z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
	ctemp.r = z__1.r, ctemp.i = z__1.i;
	i__2 = iy;
	i__3 = iy;
	z__2.r = c->r * CY(iy).r - c->i * CY(iy).i, z__2.i = c->r * CY(
		iy).i + c->i * CY(iy).r;
	i__4 = ix;
	z__3.r = s->r * CX(ix).r - s->i * CX(ix).i, z__3.i = s->r * CX(
		ix).i + s->i * CX(ix).r;
	z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
	CY(iy).r = z__1.r, CY(iy).i = z__1.i;
	i__2 = ix;
	CX(ix).r = ctemp.r, CX(ix).i = ctemp.i;
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
	z__2.r = c->r * CX(i).r - c->i * CX(i).i, z__2.i = c->r * CX(
		i).i + c->i * CX(i).r;
	i__3 = i;
	z__3.r = s->r * CY(i).r - s->i * CY(i).i, z__3.i = s->r * CY(
		i).i + s->i * CY(i).r;
	z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
	ctemp.r = z__1.r, ctemp.i = z__1.i;
	i__2 = i;
	i__3 = i;
	z__2.r = c->r * CY(i).r - c->i * CY(i).i, z__2.i = c->r * CY(
		i).i + c->i * CY(i).r;
	i__4 = i;
	z__3.r = s->r * CX(i).r - s->i * CX(i).i, z__3.i = s->r * CX(
		i).i + s->i * CX(i).r;
	z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
	CY(i).r = z__1.r, CY(i).i = z__1.i;
	i__2 = i;
	CX(i).r = ctemp.r, CX(i).i = ctemp.i;
/* L30: */
    }
    return 0;
} /* zlacrt_ */

