/*  -- translated by f2c (version 19940927).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static real c_b4 = 1.f;

/* Subroutine */ int srotg_(real *sa, real *sb, real *c, real *s)
{


    /* System generated locals */
    real r__1, r__2;

    /* Builtin functions */
    double sqrt(doublereal), r_sign(real *, real *);

    /* Local variables */
    static real r, scale, z, roe;


/*     construct givens plane rotation.   
       jack dongarra, linpack, 3/11/78. */


    roe = *sb;
    if (dabs(*sa) > dabs(*sb)) {
	roe = *sa;
    }
    scale = dabs(*sa) + dabs(*sb);
    if (scale != 0.f) {
	goto L10;
    }
    *c = 1.f;
    *s = 0.f;
    r = 0.f;
    z = 0.f;
    goto L20;
L10:
/* Computing 2nd power */
    r__1 = *sa / scale;
/* Computing 2nd power */
    r__2 = *sb / scale;
    r = scale * sqrt(r__1 * r__1 + r__2 * r__2);
    r = r_sign(&c_b4, &roe) * r;
    *c = *sa / r;
    *s = *sb / r;
    z = 1.f;
    if (dabs(*sa) > dabs(*sb)) {
	z = *s;
    }
    if (dabs(*sb) >= dabs(*sa) && *c != 0.f) {
	z = 1.f / *c;
    }
L20:
    *sa = r;
    *sb = z;
    return 0;
} /* srotg_ */

