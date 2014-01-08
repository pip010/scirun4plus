#include "f2c.h"

/* Subroutine */ int slanv2_(real *a, real *b, real *c, real *d, real *rt1r, 
	real *rt1i, real *rt2r, real *rt2i, real *cs, real *sn)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    SLANV2 computes the Schur factorization of a real 2-by-2 nonsymmetric 
  
    matrix in standard form:   

         [ A  B ] = [ CS -SN ] [ AA  BB ] [ CS  SN ]   
         [ C  D ]   [ SN  CS ] [ CC  DD ] [-SN  CS ]   

    where either   
    1) CC = 0 so that AA and DD are real eigenvalues of the matrix, or   
    2) AA = DD and BB*CC < 0, so that AA + or - sqrt(BB*CC) are complex   
    conjugate eigenvalues.   

    Arguments   
    =========   

    A       (input/output) REAL   
    B       (input/output) REAL   
    C       (input/output) REAL   
    D       (input/output) REAL   
            On entry, the elements of the input matrix.   
            On exit, they are overwritten by the elements of the   
            standardised Schur form.   

    RT1R    (output) REAL   
    RT1I    (output) REAL   
    RT2R    (output) REAL   
    RT2I    (output) REAL   
            The real and imaginary parts of the eigenvalues. If the   
            eigenvalues are both real, abs(RT1R) >= abs(RT2R); if the   
            eigenvalues are a complex conjugate pair, RT1I > 0.   

    CS      (output) REAL   
    SN      (output) REAL   
            Parameters of the rotation matrix.   

    ===================================================================== 
  


       Initialize CS and SN */
    /* Table of constant values */
    static real c_b3 = 1.f;
    
    /* System generated locals */
    real r__1;
    /* Builtin functions */
    double r_sign(real *, real *), sqrt(doublereal);
    /* Local variables */
    static real temp, p, sigma, aa, bb, cc, dd;
    extern doublereal slapy2_(real *, real *);
    static real cs1, sn1, sab, sac, tau;



    *cs = 1.f;
    *sn = 0.f;

    if (*c == 0.f) {
	goto L10;

    } else if (*b == 0.f) {

/*        Swap rows and columns */

	*cs = 0.f;
	*sn = 1.f;
	temp = *d;
	*d = *a;
	*a = temp;
	*b = -(doublereal)(*c);
	*c = 0.f;
	goto L10;
    } else if (*a - *d == 0.f && r_sign(&c_b3, b) != r_sign(&c_b3, c)) {
	goto L10;
    } else {

/*        Make diagonal elements equal */

	temp = *a - *d;
	p = temp * .5f;
	sigma = *b + *c;
	tau = slapy2_(&sigma, &temp);
	cs1 = sqrt((dabs(sigma) / tau + 1.f) * .5f);
	sn1 = -(doublereal)(p / (tau * cs1)) * r_sign(&c_b3, &sigma);

/*        Compute [ AA  BB ] = [ A  B ] [ CS1 -SN1 ]   
                  [ CC  DD ]   [ C  D ] [ SN1  CS1 ] */

	aa = *a * cs1 + *b * sn1;
	bb = -(doublereal)(*a) * sn1 + *b * cs1;
	cc = *c * cs1 + *d * sn1;
	dd = -(doublereal)(*c) * sn1 + *d * cs1;

/*        Compute [ A  B ] = [ CS1  SN1 ] [ AA  BB ]   
                  [ C  D ]   [-SN1  CS1 ] [ CC  DD ] */

	*a = aa * cs1 + cc * sn1;
	*b = bb * cs1 + dd * sn1;
	*c = -(doublereal)aa * sn1 + cc * cs1;
	*d = -(doublereal)bb * sn1 + dd * cs1;

/*        Accumulate transformation */

	temp = *cs * cs1 - *sn * sn1;
	*sn = *cs * sn1 + *sn * cs1;
	*cs = temp;

	temp = (*a + *d) * .5f;
	*a = temp;
	*d = temp;

	if (*c != 0.f) {
	    if (*b != 0.f) {
		if (r_sign(&c_b3, b) == r_sign(&c_b3, c)) {

/*                 Real eigenvalues: reduce to upper trian
gular form */

		    sab = sqrt((dabs(*b)));
		    sac = sqrt((dabs(*c)));
		    r__1 = sab * sac;
		    p = r_sign(&r__1, c);
		    tau = 1.f / sqrt((r__1 = *b + *c, dabs(r__1)));
		    *a = temp + p;
		    *d = temp - p;
		    *b -= *c;
		    *c = 0.f;
		    cs1 = sab * tau;
		    sn1 = sac * tau;
		    temp = *cs * cs1 - *sn * sn1;
		    *sn = *cs * sn1 + *sn * cs1;
		    *cs = temp;
		}
	    } else {
		*b = -(doublereal)(*c);
		*c = 0.f;
		temp = *cs;
		*cs = -(doublereal)(*sn);
		*sn = temp;
	    }
	}
    }

L10:

/*     Store eigenvalues in (RT1R,RT1I) and (RT2R,RT2I). */

    *rt1r = *a;
    *rt2r = *d;
    if (*c == 0.f) {
	*rt1i = 0.f;
	*rt2i = 0.f;
    } else {
	*rt1i = sqrt((dabs(*b))) * sqrt((dabs(*c)));
	*rt2i = -(doublereal)(*rt1i);
    }
    return 0;

/*     End of SLANV2 */

} /* slanv2_ */

