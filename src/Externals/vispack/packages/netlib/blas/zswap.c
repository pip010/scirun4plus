
/*  -- translated by f2c (version 19940927).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Subroutine */ int zswap_(integer *n, doublecomplex *zx, integer *incx, 
	doublecomplex *zy, integer *incy)
{


    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer i;
    static doublecomplex ztemp;
    static integer ix, iy;


/*     interchanges two vectors.   
       jack dongarra, 3/11/78.   
       modified 12/3/93, array(1) declarations changed to array(*)   


    
   Parameter adjustments   
       Function Body */
#define ZY(I) zy[(I)-1]
#define ZX(I) zx[(I)-1]


    if (*n <= 0) {
	return 0;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*       code for unequal increments or equal increments not equal   
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
	ztemp.r = ZX(ix).r, ztemp.i = ZX(ix).i;
	i__2 = ix;
	i__3 = iy;
	ZX(ix).r = ZY(iy).r, ZX(ix).i = ZY(iy).i;
	i__2 = iy;
	ZY(iy).r = ztemp.r, ZY(iy).i = ztemp.i;
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    return 0;

/*       code for both increments equal to 1 */
L20:
    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
	i__2 = i;
	ztemp.r = ZX(i).r, ztemp.i = ZX(i).i;
	i__2 = i;
	i__3 = i;
	ZX(i).r = ZY(i).r, ZX(i).i = ZY(i).i;
	i__2 = i;
	ZY(i).r = ztemp.r, ZY(i).i = ztemp.i;
/* L30: */
    }
    return 0;
} /* zswap_ */

