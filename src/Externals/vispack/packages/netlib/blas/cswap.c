
/*  -- translated by f2c (version 19940927).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Subroutine */ int cswap_(integer *n, complex *cx, integer *incx, complex *
	cy, integer *incy)
{


    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer i;
    static complex ctemp;
    static integer ix, iy;


/*     interchanges two vectors.   
       jack dongarra, linpack, 3/11/78.   
       modified 12/3/93, array(1) declarations changed to array(*)   


    
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
	ctemp.r = CX(ix).r, ctemp.i = CX(ix).i;
	i__2 = ix;
	i__3 = iy;
	CX(ix).r = CY(iy).r, CX(ix).i = CY(iy).i;
	i__2 = iy;
	CY(iy).r = ctemp.r, CY(iy).i = ctemp.i;
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
	ctemp.r = CX(i).r, ctemp.i = CX(i).i;
	i__2 = i;
	i__3 = i;
	CX(i).r = CY(i).r, CX(i).i = CY(i).i;
	i__2 = i;
	CY(i).r = ctemp.r, CY(i).i = ctemp.i;
/* L30: */
    }
    return 0;
} /* cswap_ */

