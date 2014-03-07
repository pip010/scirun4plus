#include "f2c.h"

/* Subroutine */ int slaein_(logical *rightv, logical *noinit, integer *n, 
	real *h, integer *ldh, real *wr, real *wi, real *vr, real *vi, real *
	b, integer *ldb, real *work, real *eps3, real *smlnum, real *bignum, 
	integer *info)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    SLAEIN uses inverse iteration to find a right or left eigenvector   
    corresponding to the eigenvalue (WR,WI) of a real upper Hessenberg   
    matrix H.   

    Arguments   
    =========   

    RIGHTV   (input) LOGICAL   
            = .TRUE. : compute right eigenvector;   
            = .FALSE.: compute left eigenvector.   

    NOINIT   (input) LOGICAL   
            = .TRUE. : no initial vector supplied in (VR,VI).   
            = .FALSE.: initial vector supplied in (VR,VI).   

    N       (input) INTEGER   
            The order of the matrix H.  N >= 0.   

    H       (input) REAL array, dimension (LDH,N)   
            The upper Hessenberg matrix H.   

    LDH     (input) INTEGER   
            The leading dimension of the array H.  LDH >= max(1,N).   

    WR      (input) REAL   
    WI      (input) REAL   
            The real and imaginary parts of the eigenvalue of H whose   
            corresponding right or left eigenvector is to be computed.   

    VR      (input/output) REAL array, dimension (N)   
    VI      (input/output) REAL array, dimension (N)   
            On entry, if NOINIT = .FALSE. and WI = 0.0, VR must contain   
            a real starting vector for inverse iteration using the real   
            eigenvalue WR; if NOINIT = .FALSE. and WI.ne.0.0, VR and VI   
            must contain the real and imaginary parts of a complex   
            starting vector for inverse iteration using the complex   
            eigenvalue (WR,WI); otherwise VR and VI need not be set.   
            On exit, if WI = 0.0 (real eigenvalue), VR contains the   
            computed real eigenvector; if WI.ne.0.0 (complex eigenvalue), 
  
            VR and VI contain the real and imaginary parts of the   
            computed complex eigenvector. The eigenvector is normalized   
            so that the component of largest magnitude has magnitude 1;   
            here the magnitude of a complex number (x,y) is taken to be   
            |x| + |y|.   
            VI is not referenced if WI = 0.0.   

    B       (workspace) REAL array, dimension (LDB,N)   

    LDB     (input) INTEGER   
            The leading dimension of the array B.  LDB >= N+1.   

    WORK   (workspace) REAL array, dimension (N)   

    EPS3    (input) REAL   
            A small machine-dependent value which is used to perturb   
            close eigenvalues, and to replace zero pivots.   

    SMLNUM  (input) REAL   
            A machine-dependent value close to the underflow threshold.   

    BIGNUM  (input) REAL   
            A machine-dependent value close to the overflow threshold.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            = 1:  inverse iteration did not converge; VR is set to the   
                  last iterate, and so is VI if WI.ne.0.0.   

    ===================================================================== 
  


    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static integer c__1 = 1;
    
    /* System generated locals */
    integer b_dim1, b_offset, h_dim1, h_offset, i__1, i__2, i__3, i__4;
    real r__1, r__2, r__3, r__4;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    static integer ierr;
    static real temp, norm, vmax;
    extern doublereal snrm2_(integer *, real *, integer *);
    static integer i, j;
    static real scale, w, x, y;
    extern /* Subroutine */ int sscal_(integer *, real *, real *, integer *);
    static char trans[1];
    static real vcrit;
    extern doublereal sasum_(integer *, real *, integer *);
    static integer i1, i2, i3;
    static real rootn, vnorm, w1;
    extern doublereal slapy2_(real *, real *);
    static real ei, ej, absbii, absbjj, xi, xr;
    extern integer isamax_(integer *, real *, integer *);
    extern /* Subroutine */ int sladiv_(real *, real *, real *, real *, real *
	    , real *);
    static char normin[1];
    static real nrmsml;
    extern /* Subroutine */ int slatrs_(char *, char *, char *, char *, 
	    integer *, real *, integer *, real *, real *, real *, integer *);
    static real growto, rec;
    static integer its;



#define VR(I) vr[(I)-1]
#define VI(I) vi[(I)-1]
#define WORK(I) work[(I)-1]

#define H(I,J) h[(I)-1 + ((J)-1)* ( *ldh)]
#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]

    *info = 0;

/*     GROWTO is the threshold used in the acceptance test for an   
       eigenvector. */

    rootn = sqrt((real) (*n));
    growto = .1f / rootn;
/* Computing MAX */
    r__1 = 1.f, r__2 = *eps3 * rootn;
    nrmsml = dmax(r__1,r__2) * *smlnum;

/*     Form B = H - (WR,WI)*I (except that the subdiagonal elements and   
       the imaginary parts of the diagonal elements are not stored). */

    i__1 = *n;
    for (j = 1; j <= *n; ++j) {
	i__2 = j - 1;
	for (i = 1; i <= j-1; ++i) {
	    B(i,j) = H(i,j);
/* L10: */
	}
	B(j,j) = H(j,j) - *wr;
/* L20: */
    }

    if (*wi == 0.f) {

/*        Real eigenvalue. */

	if (*noinit) {

/*           Set initial vector. */

	    i__1 = *n;
	    for (i = 1; i <= *n; ++i) {
		VR(i) = *eps3;
/* L30: */
	    }
	} else {

/*           Scale supplied initial vector. */

	    vnorm = snrm2_(n, &VR(1), &c__1);
	    r__1 = *eps3 * rootn / dmax(vnorm,nrmsml);
	    sscal_(n, &r__1, &VR(1), &c__1);
	}

	if (*rightv) {

/*           LU decomposition with partial pivoting of B, replacin
g zero   
             pivots by EPS3. */

	    i__1 = *n - 1;
	    for (i = 1; i <= *n-1; ++i) {
		ei = H(i+1,i);
		if ((r__1 = B(i,i), dabs(r__1)) < dabs(ei)) {

/*                 Interchange rows and eliminate. */

		    x = B(i,i) / ei;
		    B(i,i) = ei;
		    i__2 = *n;
		    for (j = i + 1; j <= *n; ++j) {
			temp = B(i+1,j);
			B(i+1,j) = B(i,j) - x * temp;
			B(i,j) = temp;
/* L40: */
		    }
		} else {

/*                 Eliminate without interchange. */

		    if (B(i,i) == 0.f) {
			B(i,i) = *eps3;
		    }
		    x = ei / B(i,i);
		    if (x != 0.f) {
			i__2 = *n;
			for (j = i + 1; j <= *n; ++j) {
			    B(i+1,j) -= x * B(i,j);
/* L50: */
			}
		    }
		}
/* L60: */
	    }
	    if (B(*n,*n) == 0.f) {
		B(*n,*n) = *eps3;
	    }

	    *(unsigned char *)trans = 'N';

	} else {

/*           UL decomposition with partial pivoting of B, replacin
g zero   
             pivots by EPS3. */

	    for (j = *n; j >= 2; --j) {
		ej = H(j,j-1);
		if ((r__1 = B(j,j), dabs(r__1)) < dabs(ej)) {

/*                 Interchange columns and eliminate. */

		    x = B(j,j) / ej;
		    B(j,j) = ej;
		    i__1 = j - 1;
		    for (i = 1; i <= j-1; ++i) {
			temp = B(i,j-1);
			B(i,j-1) = B(i,j) - x * 
				temp;
			B(i,j) = temp;
/* L70: */
		    }
		} else {

/*                 Eliminate without interchange. */

		    if (B(j,j) == 0.f) {
			B(j,j) = *eps3;
		    }
		    x = ej / B(j,j);
		    if (x != 0.f) {
			i__1 = j - 1;
			for (i = 1; i <= j-1; ++i) {
			    B(i,j-1) -= x * B(i,j);
/* L80: */
			}
		    }
		}
/* L90: */
	    }
	    if (B(1,1) == 0.f) {
		B(1,1) = *eps3;
	    }

	    *(unsigned char *)trans = 'T';

	}

	*(unsigned char *)normin = 'N';
	i__1 = *n;
	for (its = 1; its <= *n; ++its) {

/*           Solve U*x = scale*v for a right eigenvector   
               or U'*x = scale*v for a left eigenvector,   
             overwriting x on v. */

	    slatrs_("Upper", trans, "Nonunit", normin, n, &B(1,1), ldb, &
		    VR(1), &scale, &WORK(1), &ierr);
	    *(unsigned char *)normin = 'Y';

/*           Test for sufficient growth in the norm of v. */

	    vnorm = sasum_(n, &VR(1), &c__1);
	    if (vnorm >= growto * scale) {
		goto L120;
	    }

/*           Choose new orthogonal starting vector and try again. 
*/

	    temp = *eps3 / (rootn + 1.f);
	    VR(1) = *eps3;
	    i__2 = *n;
	    for (i = 2; i <= *n; ++i) {
		VR(i) = temp;
/* L100: */
	    }
	    VR(*n - its + 1) -= *eps3 * rootn;
/* L110: */
	}

/*        Failure to find eigenvector in N iterations. */

	*info = 1;

L120:

/*        Normalize eigenvector. */

	i = isamax_(n, &VR(1), &c__1);
	r__2 = 1.f / (r__1 = VR(i), dabs(r__1));
	sscal_(n, &r__2, &VR(1), &c__1);
    } else {

/*        Complex eigenvalue. */

	if (*noinit) {

/*           Set initial vector. */

	    i__1 = *n;
	    for (i = 1; i <= *n; ++i) {
		VR(i) = *eps3;
		VI(i) = 0.f;
/* L130: */
	    }
	} else {

/*           Scale supplied initial vector. */

	    r__1 = snrm2_(n, &VR(1), &c__1);
	    r__2 = snrm2_(n, &VI(1), &c__1);
	    norm = slapy2_(&r__1, &r__2);
	    rec = *eps3 * rootn / dmax(norm,nrmsml);
	    sscal_(n, &rec, &VR(1), &c__1);
	    sscal_(n, &rec, &VI(1), &c__1);
	}

	if (*rightv) {

/*           LU decomposition with partial pivoting of B, replacin
g zero   
             pivots by EPS3.   

             The imaginary part of the (i,j)-th element of U is st
ored in   
             B(j+1,i). */

	    B(2,1) = -(doublereal)(*wi);
	    i__1 = *n;
	    for (i = 2; i <= *n; ++i) {
		B(i+1,1) = 0.f;
/* L140: */
	    }

	    i__1 = *n - 1;
	    for (i = 1; i <= *n-1; ++i) {
		absbii = slapy2_(&B(i,i), &B(i+1,i));
		ei = H(i+1,i);
		if (absbii < dabs(ei)) {

/*                 Interchange rows and eliminate. */

		    xr = B(i,i) / ei;
		    xi = B(i+1,i) / ei;
		    B(i,i) = ei;
		    B(i+1,i) = 0.f;
		    i__2 = *n;
		    for (j = i + 1; j <= *n; ++j) {
			temp = B(i+1,j);
			B(i+1,j) = B(i,j) - xr * temp;
			B(j+1,i+1) = B(j+1,i) - 
				xi * temp;
			B(i,j) = temp;
			B(j+1,i) = 0.f;
/* L150: */
		    }
		    B(i+2,i) = -(doublereal)(*wi);
		    B(i+1,i+1) -= xi * *wi;
		    B(i+2,i+1) += xr * *wi;
		} else {

/*                 Eliminate without interchanging rows. 
*/

		    if (absbii == 0.f) {
			B(i,i) = *eps3;
			B(i+1,i) = 0.f;
			absbii = *eps3;
		    }
		    ei = ei / absbii / absbii;
		    xr = B(i,i) * ei;
		    xi = -(doublereal)B(i+1,i) * ei;
		    i__2 = *n;
		    for (j = i + 1; j <= *n; ++j) {
			B(i+1,j) = B(i+1,j) - xr * 
				B(i,j) + xi * B(j+1,i)
				;
			B(j+1,i+1) = -(doublereal)xr * B(j+1,i) - xi * B(i,j);
/* L160: */
		    }
		    B(i+2,i+1) -= *wi;
		}

/*              Compute 1-norm of offdiagonal elements of i-th
 row. */

		i__2 = *n - i;
		i__3 = *n - i;
		WORK(i) = sasum_(&i__2, &B(i,i+1), ldb) + 
			sasum_(&i__3, &B(i+2,i), &c__1);
/* L170: */
	    }
	    if (B(*n,*n) == 0.f && B(*n+1,*n) == 0.f) 
		    {
		B(*n,*n) = *eps3;
	    }
	    WORK(*n) = 0.f;

	    i1 = *n;
	    i2 = 1;
	    i3 = -1;
	} else {

/*           UL decomposition with partial pivoting of conjg(B), 
  
             replacing zero pivots by EPS3.   

             The imaginary part of the (i,j)-th element of U is st
ored in   
             B(j+1,i). */

	    B(*n+1,*n) = *wi;
	    i__1 = *n - 1;
	    for (j = 1; j <= *n-1; ++j) {
		B(*n+1,j) = 0.f;
/* L180: */
	    }

	    for (j = *n; j >= 2; --j) {
		ej = H(j,j-1);
		absbjj = slapy2_(&B(j,j), &B(j+1,j));
		if (absbjj < dabs(ej)) {

/*                 Interchange columns and eliminate */

		    xr = B(j,j) / ej;
		    xi = B(j+1,j) / ej;
		    B(j,j) = ej;
		    B(j+1,j) = 0.f;
		    i__1 = j - 1;
		    for (i = 1; i <= j-1; ++i) {
			temp = B(i,j-1);
			B(i,j-1) = B(i,j) - xr * 
				temp;
			B(j,i) = B(j+1,i) - xi * temp;
			B(i,j) = temp;
			B(j+1,i) = 0.f;
/* L190: */
		    }
		    B(j+1,j-1) = *wi;
		    B(j-1,j-1) += xi * *wi;
		    B(j,j-1) -= xr * *wi;
		} else {

/*                 Eliminate without interchange. */

		    if (absbjj == 0.f) {
			B(j,j) = *eps3;
			B(j+1,j) = 0.f;
			absbjj = *eps3;
		    }
		    ej = ej / absbjj / absbjj;
		    xr = B(j,j) * ej;
		    xi = -(doublereal)B(j+1,j) * ej;
		    i__1 = j - 1;
		    for (i = 1; i <= j-1; ++i) {
			B(i,j-1) = B(i,j-1) - 
				xr * B(i,j) + xi * B(j+1,i);
			B(j,i) = -(doublereal)xr * B(j+1,i) - xi * B(i,j);
/* L200: */
		    }
		    B(j,j-1) += *wi;
		}

/*              Compute 1-norm of offdiagonal elements of j-th
 column. */

		i__1 = j - 1;
		i__2 = j - 1;
		WORK(j) = sasum_(&i__1, &B(1,j), &c__1) + sasum_(&
			i__2, &B(j+1,1), ldb);
/* L210: */
	    }
	    if (B(1,1) == 0.f && B(2,1) == 0.f) {
		B(1,1) = *eps3;
	    }
	    WORK(1) = 0.f;

	    i1 = 1;
	    i2 = *n;
	    i3 = 1;
	}

	i__1 = *n;
	for (its = 1; its <= *n; ++its) {
	    scale = 1.f;
	    vmax = 1.f;
	    vcrit = *bignum;

/*           Solve U*(xr,xi) = scale*(vr,vi) for a right eigenvect
or,   
               or U'*(xr,xi) = scale*(vr,vi) for a left eigenvecto
r,   
             overwriting (xr,xi) on (vr,vi). */

	    i__2 = i2;
	    i__3 = i3;
	    for (i = i1; i3 < 0 ? i >= i2 : i <= i2; i += i3) {

		if (WORK(i) > vcrit) {
		    rec = 1.f / vmax;
		    sscal_(n, &rec, &VR(1), &c__1);
		    sscal_(n, &rec, &VI(1), &c__1);
		    scale *= rec;
		    vmax = 1.f;
		    vcrit = *bignum;
		}

		xr = VR(i);
		xi = VI(i);
		if (*rightv) {
		    i__4 = *n;
		    for (j = i + 1; j <= *n; ++j) {
			xr = xr - B(i,j) * VR(j) + B(j+1,i) * VI(j);
			xi = xi - B(i,j) * VI(j) - B(j+1,i) * VR(j);
/* L220: */
		    }
		} else {
		    i__4 = i - 1;
		    for (j = 1; j <= i-1; ++j) {
			xr = xr - B(j,i) * VR(j) + B(i+1,j) * VI(j);
			xi = xi - B(j,i) * VI(j) - B(i+1,j) * VR(j);
/* L230: */
		    }
		}

		w = (r__1 = B(i,i), dabs(r__1)) + (r__2 = B(i+1,i), dabs(r__2));
		if (w > *smlnum) {
		    if (w < 1.f) {
			w1 = dabs(xr) + dabs(xi);
			if (w1 > w * *bignum) {
			    rec = 1.f / w1;
			    sscal_(n, &rec, &VR(1), &c__1);
			    sscal_(n, &rec, &VI(1), &c__1);
			    xr = VR(i);
			    xi = VI(i);
			    scale *= rec;
			    vmax *= rec;
			}
		    }

/*                 Divide by diagonal element of B. */

		    sladiv_(&xr, &xi, &B(i,i), &B(i+1,i), &VR(i), &VI(i));
/* Computing MAX */
		    r__3 = (r__1 = VR(i), dabs(r__1)) + (r__2 = VI(i), dabs(
			    r__2));
		    vmax = dmax(r__3,vmax);
		    vcrit = *bignum / vmax;
		} else {
		    i__4 = *n;
		    for (j = 1; j <= *n; ++j) {
			VR(j) = 0.f;
			VI(j) = 0.f;
/* L240: */
		    }
		    VR(i) = 1.f;
		    VI(i) = 1.f;
		    scale = 0.f;
		    vmax = 1.f;
		    vcrit = *bignum;
		}
/* L250: */
	    }

/*           Test for sufficient growth in the norm of (VR,VI). */

	    vnorm = sasum_(n, &VR(1), &c__1) + sasum_(n, &VI(1), &c__1);
	    if (vnorm >= growto * scale) {
		goto L280;
	    }

/*           Choose a new orthogonal starting vector and try again
. */

	    y = *eps3 / (rootn + 1.f);
	    VR(1) = *eps3;
	    VI(1) = 0.f;

	    i__3 = *n;
	    for (i = 2; i <= *n; ++i) {
		VR(i) = y;
		VI(i) = 0.f;
/* L260: */
	    }
	    VR(*n - its + 1) -= *eps3 * rootn;
/* L270: */
	}

/*        Failure to find eigenvector in N iterations */

	*info = 1;

L280:

/*        Normalize eigenvector. */

	vnorm = 0.f;
	i__1 = *n;
	for (i = 1; i <= *n; ++i) {
/* Computing MAX */
	    r__3 = vnorm, r__4 = (r__1 = VR(i), dabs(r__1)) + (r__2 = VI(i), 
		    dabs(r__2));
	    vnorm = dmax(r__3,r__4);
/* L290: */
	}
	r__1 = 1.f / vnorm;
	sscal_(n, &r__1, &VR(1), &c__1);
	r__1 = 1.f / vnorm;
	sscal_(n, &r__1, &VI(1), &c__1);

    }

    return 0;

/*     End of SLAEIN */

} /* slaein_ */

