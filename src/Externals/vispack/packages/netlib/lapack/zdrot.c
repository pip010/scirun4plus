#include "f2c.h"

/* Subroutine */ int zdrot_(integer *n, doublecomplex *cx, integer *incx, 
	doublecomplex *cy, integer *incy, doublereal *c, doublereal *s)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    doublecomplex z__1, z__2, z__3;
    /* Local variables */
    static integer i;
    static doublecomplex ctemp;
    static integer ix, iy;
/*     applies a plane rotation, where the cos and sin (c and s) are real 
  
       and the vectors cx and cy are complex.   
       jack dongarra, linpack, 3/11/78.   
   ===================================================================== 
  
    
   Parameter adjustments   
       Function Body */
#define CY(I) cy[(I)-1]
#define CX(I) cx[(I)-1]
    if (*n <= 0) {
	return 0;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }
/*        code for unequal increments or equal increments not equal   
            to 1 */
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
	z__3.r = *s * CY(iy).r, z__3.i = *s * CY(iy).i;
	z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
	ctemp.r = z__1.r, ctemp.i = z__1.i;
	i__2 = iy;
	i__3 = iy;
	z__2.r = *c * CY(iy).r, z__2.i = *c * CY(iy).i;
	i__4 = ix;
	z__3.r = *s * CX(ix).r, z__3.i = *s * CX(ix).i;
	z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
	CY(iy).r = z__1.r, CY(iy).i = z__1.i;
	i__2 = ix;
	CX(ix).r = ctemp.r, CX(ix).i = ctemp.i;
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    return 0;
/*        code for both increments equal to 1 */
L20:
    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
	i__2 = i;
	z__2.r = *c * CX(i).r, z__2.i = *c * CX(i).i;
	i__3 = i;
	z__3.r = *s * CY(i).r, z__3.i = *s * CY(i).i;
	z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
	ctemp.r = z__1.r, ctemp.i = z__1.i;
	i__2 = i;
	i__3 = i;
	z__2.r = *c * CY(i).r, z__2.i = *c * CY(i).i;
	i__4 = i;
	z__3.r = *s * CX(i).r, z__3.i = *s * CX(i).i;
	z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
	CY(i).r = z__1.r, CY(i).i = z__1.i;
	i__2 = i;
	CX(i).r = ctemp.r, CX(i).i = ctemp.i;
/* L30: */
    }
    return 0;
} /* zdrot_ */


