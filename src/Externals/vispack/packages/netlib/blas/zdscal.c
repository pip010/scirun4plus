
/*  -- translated by f2c (version 19940927).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Subroutine */ int zdscal_(integer *n, doublereal *da, doublecomplex *zx, 
	integer *incx)
{


    /* System generated locals */
    integer i__1, i__2, i__3;
    doublecomplex z__1, z__2;

    /* Local variables */
    static integer i, ix;


/*     scales a vector by a constant.   
       jack dongarra, 3/11/78.   
       modified 3/93 to return if incx .le. 0.   
       modified 12/3/93, array(1) declarations changed to array(*)   


    
   Parameter adjustments   
       Function Body */
#define ZX(I) zx[(I)-1]


    if (*n <= 0 || *incx <= 0) {
	return 0;
    }
    if (*incx == 1) {
	goto L20;
    }

/*        code for increment not equal to 1 */

    ix = 1;
    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
	i__2 = ix;
	z__2.r = *da, z__2.i = 0.;
	i__3 = ix;
	z__1.r = z__2.r * ZX(ix).r - z__2.i * ZX(ix).i, z__1.i = z__2.r * 
		ZX(ix).i + z__2.i * ZX(ix).r;
	ZX(ix).r = z__1.r, ZX(ix).i = z__1.i;
	ix += *incx;
/* L10: */
    }
    return 0;

/*        code for increment equal to 1 */

L20:
    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
	i__2 = i;
	z__2.r = *da, z__2.i = 0.;
	i__3 = i;
	z__1.r = z__2.r * ZX(i).r - z__2.i * ZX(i).i, z__1.i = z__2.r * 
		ZX(i).i + z__2.i * ZX(i).r;
	ZX(i).r = z__1.r, ZX(i).i = z__1.i;
/* L30: */
    }
    return 0;
} /* zdscal_ */

