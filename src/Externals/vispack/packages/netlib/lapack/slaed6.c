#include "f2c.h"

/* Subroutine */ int slaed6_(integer *kniter, logical *orgati, real *rho, 
	real *d, real *z, real *finit, real *tau, integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Oak Ridge National Lab, Argonne National Lab, 
  
       Courant Institute, NAG Ltd., and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    SLAED6 computes the positive or negative root (closest to the origin) 
  
    of   
                     z(1)        z(2)        z(3)   
    f(x) =   rho + --------- + ---------- + ---------   
                    d(1)-x      d(2)-x      d(3)-x   

    It is assumed that   

          if ORGATI = .true. the root is between d(2) and d(3);   
          otherwise it is between d(1) and d(2)   

    This routine will be called by SLAED4 when necessary. In most cases, 
  
    the root sought is the smallest in magnitude, though it might not be 
  
    in some extremely rare situations.   

    Arguments   
    =========   

    KNITER       (input) INTEGER   
                 Refer to SLAED4 for its significance.   

    ORGATI       (input) LOGICAL   
                 If ORGATI is true, the needed root is between d(2) and   
                 d(3); otherwise it is between d(1) and d(2).  See   
                 SLAED4 for further details.   

    RHO          (input) REAL   
                 Refer to the equation f(x) above.   

    D            (input) REAL array, dimension (3)   
                 D satisfies d(1) < d(2) < d(3).   

    Z            (input) REAL array, dimension (3)   
                 Each of the elements in z must be positive.   

    FINIT        (input) REAL   
                 The value of f at 0. It is more accurate than the one   
                 evaluated inside this routine (if someone wants to do   
                 so).   

    TAU          (output) REAL   
                 The root of the equation f(x).   

    INFO         (output) INTEGER   
                 = 0: successful exit   
                 > 0: if INFO = 1, failure to converge   

    ===================================================================== 
  

    
   Parameter adjustments   
       Function Body */
    /* Initialized data */
    static logical first = TRUE_;
    /* System generated locals */
    integer i__1;
    real r__1, r__2, r__3, r__4;
    /* Builtin functions */
    double sqrt(doublereal), log(doublereal), pow_ri(real *, integer *);
    /* Local variables */
    static real base;
    static integer iter;
    static real temp, temp1, temp2, temp3, a, b, c, f;
    static integer i;
    static logical scale;
    static integer niter;
    static real small1, small2, fc, df, sminv1, sminv2, dscale[3], sclfac;
    extern doublereal slamch_(char *);
    static real zscale[3], erretm, sclinv, pretau, ddf, eta, eps;


#define Z(I) z[(I)-1]
#define D(I) d[(I)-1]



    *info = 0;

    niter = 1;
    *tau = 0.f;
    if (*kniter == 2) {
	if (*orgati) {
	    temp = (D(3) - D(2)) / 2.f;
	    c = *rho + Z(1) / (D(1) - D(2) - temp);
	    a = c * (D(2) + D(3)) + Z(2) + Z(3);
	    b = c * D(2) * D(3) + Z(2) * D(3) + Z(3) * D(2);
	} else {
	    temp = (D(1) - D(2)) / 2.f;
	    c = *rho + Z(3) / (D(3) - D(2) - temp);
	    a = c * (D(1) + D(2)) + Z(1) + Z(2);
	    b = c * D(1) * D(2) + Z(1) * D(2) + Z(2) * D(1);
	}
	if (c == 0.f) {
	    *tau = b / a;
	} else if (a <= 0.f) {
	    *tau = (a - sqrt((r__1 = a * a - b * 4.f * c, dabs(r__1)))) / (c *
		     2.f);
	} else {
	    *tau = b * 2.f / (a + sqrt((r__1 = a * a - b * 4.f * c, dabs(r__1)
		    )));
	}
	temp = *rho + Z(1) / (D(1) - *tau) + Z(2) / (D(2) - *tau) + Z(3) / (D(
		3) - *tau);
	if (dabs(*finit) <= dabs(temp)) {
	    *tau = 0.f;
	}
    }

/*     On first call to routine, get machine parameters for   
       possible scaling to avoid overflow */

    if (first) {
	eps = slamch_("Epsilon");
	base = slamch_("Base");
	i__1 = (integer) (log(slamch_("SafMin")) / log(base) / 3.f);
	small1 = pow_ri(&base, &i__1);
	sminv1 = 1.f / small1;
	small2 = small1 * small1;
	sminv2 = sminv1 * sminv1;
	first = FALSE_;
    }

/*     Determine if scaling of inputs necessary to avoid overflow   
       when computing 1/TEMP**3 */

    if (*orgati) {
/* Computing MIN */
	r__3 = (r__1 = D(2) - *tau, dabs(r__1)), r__4 = (r__2 = D(3) - *tau, 
		dabs(r__2));
	temp = dmin(r__3,r__4);
    } else {
/* Computing MIN */
	r__3 = (r__1 = D(1) - *tau, dabs(r__1)), r__4 = (r__2 = D(2) - *tau, 
		dabs(r__2));
	temp = dmin(r__3,r__4);
    }
    scale = FALSE_;
    if (temp <= small1) {
	scale = TRUE_;
	if (temp <= small2) {

/*        Scale up by power of radix nearest 1/SAFMIN**(2/3) */

	    sclfac = sminv2;
	    sclinv = small2;
	} else {

/*        Scale up by power of radix nearest 1/SAFMIN**(1/3) */

	    sclfac = sminv1;
	    sclinv = small1;
	}

/*        Scaling up safe because D, Z, TAU scaled elsewhere to be O(1
) */

	for (i = 1; i <= 3; ++i) {
	    dscale[i - 1] = D(i) * sclfac;
	    zscale[i - 1] = Z(i) * sclfac;
/* L10: */
	}
	*tau *= sclfac;
    } else {

/*        Copy D and Z to DSCALE and ZSCALE */

	for (i = 1; i <= 3; ++i) {
	    dscale[i - 1] = D(i);
	    zscale[i - 1] = Z(i);
/* L20: */
	}
    }

    fc = 0.f;
    df = 0.f;
    ddf = 0.f;
    for (i = 1; i <= 3; ++i) {
	temp = 1.f / (dscale[i - 1] - *tau);
	temp1 = zscale[i - 1] * temp;
	temp2 = temp1 * temp;
	temp3 = temp2 * temp;
	fc += temp1 / dscale[i - 1];
	df += temp2;
	ddf += temp3;
/* L30: */
    }
    f = *finit + *tau * fc;

    if (dabs(f) <= 0.f) {
	goto L60;
    }

/*        Iteration begins   

       It is not hard to see that   

             1) Iterations will go up monotonically   
                if FINIT < 0;   

             2) Iterations will go down monotonically   
                if FINIT > 0. */

    iter = niter + 1;

    for (niter = iter; niter <= 20; ++niter) {

	if (*orgati) {
	    temp1 = dscale[1] - *tau;
	    temp2 = dscale[2] - *tau;
	} else {
	    temp1 = dscale[0] - *tau;
	    temp2 = dscale[1] - *tau;
	}
	a = (temp1 + temp2) * f - temp1 * temp2 * df;
	b = temp1 * temp2 * f;
	c = f - (temp1 + temp2) * df + temp1 * temp2 * ddf;
	if (c == 0.f) {
	    eta = b / a;
	} else if (a <= 0.f) {
	    eta = (a - sqrt((r__1 = a * a - b * 4.f * c, dabs(r__1)))) / (c * 
		    2.f);
	} else {
	    eta = b * 2.f / (a + sqrt((r__1 = a * a - b * 4.f * c, dabs(r__1))
		    ));
	}
	if (f * eta >= 0.f) {
	    eta = -(doublereal)f / df;
	}

	temp = eta + *tau;
	if (*orgati) {
	    if (eta > 0.f && temp >= dscale[2]) {
		eta = (dscale[2] - *tau) / 2.f;
	    }
	    if (eta < 0.f && temp <= dscale[1]) {
		eta = (dscale[1] - *tau) / 2.f;
	    }
	} else {
	    if (eta > 0.f && temp >= dscale[1]) {
		eta = (dscale[1] - *tau) / 2.f;
	    }
	    if (eta < 0.f && temp <= dscale[0]) {
		eta = (dscale[0] - *tau) / 2.f;
	    }
	}
	pretau = *tau;
	*tau += eta;

	fc = 0.f;
	erretm = dabs(*rho);
	df = 0.f;
	ddf = 0.f;
	for (i = 1; i <= 3; ++i) {
	    temp = 1.f / (dscale[i - 1] - *tau);
	    temp1 = zscale[i - 1] * temp;
	    temp2 = temp1 * temp;
	    temp3 = temp2 * temp;
	    fc += temp1 / (dscale[i - 1] - pretau);
	    erretm += dabs(temp1);
	    df += temp2;
	    ddf += temp3;
/* L40: */
	}
	f += eta * fc;
	erretm = erretm * 8.f + dabs(*tau) * df;
	if (dabs(f) <= eps * erretm) {
	    goto L60;
	}
/* L50: */
    }
    *info = 1;
L60:

/*     Undo scaling */

    if (scale) {
	*tau *= sclinv;
    }
    return 0;

/*     End of SLAED6 */

} /* slaed6_ */

