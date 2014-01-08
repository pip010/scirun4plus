#include "f2c.h"

/* Subroutine */ int dlanv2_(doublereal *a, doublereal *b, doublereal *c, 
	doublereal *d, doublereal *rt1r, doublereal *rt1i, doublereal *rt2r, 
	doublereal *rt2i, doublereal *cs, doublereal *sn)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DLANV2 computes the Schur factorization of a real 2-by-2 nonsymmetric 
  
    matrix in standard form:   

         [ A  B ] = [ CS -SN ] [ AA  BB ] [ CS  SN ]   
         [ C  D ]   [ SN  CS ] [ CC  DD ] [-SN  CS ]   

    where either   
    1) CC = 0 so that AA and DD are real eigenvalues of the matrix, or   
    2) AA = DD and BB*CC < 0, so that AA + or - sqrt(BB*CC) are complex   
    conjugate eigenvalues.   

    Arguments   
    =========   

    A       (input/output) DOUBLE PRECISION   
    B       (input/output) DOUBLE PRECISION   
    C       (input/output) DOUBLE PRECISION   
    D       (input/output) DOUBLE PRECISION   
            On entry, the elements of the input matrix.   
            On exit, they are overwritten by the elements of the   
            standardised Schur form.   

    RT1R    (output) DOUBLE PRECISION   
    RT1I    (output) DOUBLE PRECISION   
    RT2R    (output) DOUBLE PRECISION   
    RT2I    (output) DOUBLE PRECISION   
            The real and imaginary parts of the eigenvalues. If the   
            eigenvalues are both real, abs(RT1R) >= abs(RT2R); if the   
            eigenvalues are a complex conjugate pair, RT1I > 0.   

    CS      (output) DOUBLE PRECISION   
    SN      (output) DOUBLE PRECISION   
            Parameters of the rotation matrix.   

    ===================================================================== 
  


       Initialize CS and SN */
    /* Table of constant values */
    static doublereal c_b3 = 1.;
    
    /* System generated locals */
    doublereal d__1;
    /* Builtin functions */
    double d_sign(doublereal *, doublereal *), sqrt(doublereal);
    /* Local variables */
    static doublereal temp, p, sigma;
    extern doublereal dlapy2_(doublereal *, doublereal *);
    static doublereal aa, bb, cc, dd, cs1, sn1, sab, sac, tau;



    *cs = 1.;
    *sn = 0.;

    if (*c == 0.) {
	goto L10;

    } else if (*b == 0.) {

/*        Swap rows and columns */

	*cs = 0.;
	*sn = 1.;
	temp = *d;
	*d = *a;
	*a = temp;
	*b = -(*c);
	*c = 0.;
	goto L10;
    } else if (*a - *d == 0. && d_sign(&c_b3, b) != d_sign(&c_b3, c)) {
	goto L10;
    } else {

/*        Make diagonal elements equal */

	temp = *a - *d;
	p = temp * .5;
	sigma = *b + *c;
	tau = dlapy2_(&sigma, &temp);
	cs1 = sqrt((abs(sigma) / tau + 1.) * .5);
	sn1 = -(p / (tau * cs1)) * d_sign(&c_b3, &sigma);

/*        Compute [ AA  BB ] = [ A  B ] [ CS1 -SN1 ]   
                  [ CC  DD ]   [ C  D ] [ SN1  CS1 ] */

	aa = *a * cs1 + *b * sn1;
	bb = -(*a) * sn1 + *b * cs1;
	cc = *c * cs1 + *d * sn1;
	dd = -(*c) * sn1 + *d * cs1;

/*        Compute [ A  B ] = [ CS1  SN1 ] [ AA  BB ]   
                  [ C  D ]   [-SN1  CS1 ] [ CC  DD ] */

	*a = aa * cs1 + cc * sn1;
	*b = bb * cs1 + dd * sn1;
	*c = -aa * sn1 + cc * cs1;
	*d = -bb * sn1 + dd * cs1;

/*        Accumulate transformation */

	temp = *cs * cs1 - *sn * sn1;
	*sn = *cs * sn1 + *sn * cs1;
	*cs = temp;

	temp = (*a + *d) * .5;
	*a = temp;
	*d = temp;

	if (*c != 0.) {
	    if (*b != 0.) {
		if (d_sign(&c_b3, b) == d_sign(&c_b3, c)) {

/*                 Real eigenvalues: reduce to upper trian
gular form */

		    sab = sqrt((abs(*b)));
		    sac = sqrt((abs(*c)));
		    d__1 = sab * sac;
		    p = d_sign(&d__1, c);
		    tau = 1. / sqrt((d__1 = *b + *c, abs(d__1)));
		    *a = temp + p;
		    *d = temp - p;
		    *b -= *c;
		    *c = 0.;
		    cs1 = sab * tau;
		    sn1 = sac * tau;
		    temp = *cs * cs1 - *sn * sn1;
		    *sn = *cs * sn1 + *sn * cs1;
		    *cs = temp;
		}
	    } else {
		*b = -(*c);
		*c = 0.;
		temp = *cs;
		*cs = -(*sn);
		*sn = temp;
	    }
	}
    }

L10:

/*     Store eigenvalues in (RT1R,RT1I) and (RT2R,RT2I). */

    *rt1r = *a;
    *rt2r = *d;
    if (*c == 0.) {
	*rt1i = 0.;
	*rt2i = 0.;
    } else {
	*rt1i = sqrt((abs(*b))) * sqrt((abs(*c)));
	*rt2i = -(*rt1i);
    }
    return 0;

/*     End of DLANV2 */

} /* dlanv2_ */

