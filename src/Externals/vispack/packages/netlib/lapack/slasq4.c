#include "f2c.h"

/* Subroutine */ int slasq4_(integer *n, real *q, real *e, real *tau, real *
	sup)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


       Purpose   
       =======   

       SLASQ4 estimates TAU, the smallest eigenvalue of a matrix. This   
       routine improves the input value of SUP which is an upper bound   
       for the smallest eigenvalue for this matrix .   

       Arguments   
       =========   

    N       (input) INTEGER   
            On entry, N specifies the number of rows and columns   
            in the matrix. N must be at least 0.   

    Q       (input) REAL array, dimension (N)   
            Q array   

    E       (input) REAL array, dimension (N)   
            E array   

    TAU     (output) REAL   
            Estimate of the shift   

    SUP     (input/output) REAL   
            Upper bound for the smallest singular value   

    ===================================================================== 
  

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static real c_b4 = .7f;
    
    /* System generated locals */
    integer i__1;
    real r__1, r__2;
    /* Builtin functions */
    double pow_ri(real *, integer *);
    /* Local variables */
    static real xinf, d;
    static integer i;
    static real dm;
    static integer ifl;



#define E(I) e[(I)-1]
#define Q(I) q[(I)-1]


    ifl = 1;
/* Computing MIN */
    r__1 = min(*sup,Q(1)), r__1 = min(r__1,Q(2)), r__1 = min(r__1,Q(3)), r__2 
	    = Q(*n), r__1 = min(r__1,r__2), r__2 = Q(*n - 1), r__1 = min(r__1,
	    r__2), r__2 = Q(*n - 2);
    *sup = dmin(r__1,r__2);
    *tau = *sup * .9999f;
    xinf = 0.f;
L10:
    if (ifl == 5) {
	*tau = xinf;
	return 0;
    }
    d = Q(1) - *tau;
    dm = d;
    i__1 = *n - 2;
    for (i = 1; i <= *n-2; ++i) {
	d = d / (d + E(i)) * Q(i + 1) - *tau;
	if (dm > d) {
	    dm = d;
	}
	if (d < 0.f) {
	    *sup = *tau;
/* Computing MAX */
	    r__1 = *sup * pow_ri(&c_b4, &ifl), r__2 = d + *tau;
	    *tau = dmax(r__1,r__2);
	    ++ifl;
	    goto L10;
	}
/* L20: */
    }
    d = d / (d + E(*n - 1)) * Q(*n) - *tau;
    if (dm > d) {
	dm = d;
    }
    if (d < 0.f) {
	*sup = *tau;
/* Computing MAX */
	r__1 = xinf, r__2 = d + *tau;
	xinf = dmax(r__1,r__2);
	if (*sup * pow_ri(&c_b4, &ifl) <= xinf) {
	    *tau = xinf;
	} else {
	    *tau = *sup * pow_ri(&c_b4, &ifl);
	    ++ifl;
	    goto L10;
	}
    } else {
/* Computing MIN */
	r__1 = *sup, r__2 = dm + *tau;
	*sup = dmin(r__1,r__2);
    }
    return 0;

/*     End of SLASQ4 */

} /* slasq4_ */

