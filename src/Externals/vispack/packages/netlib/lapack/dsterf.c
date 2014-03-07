#include "f2c.h"

/* Subroutine */ int dsterf_(integer *n, doublereal *d, doublereal *e, 
	integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DSTERF computes all eigenvalues of a symmetric tridiagonal matrix   
    using the Pal-Walker-Kahan variant of the QL or QR algorithm.   

    Arguments   
    =========   

    N       (input) INTEGER   
            The order of the matrix.  N >= 0.   

    D       (input/output) DOUBLE PRECISION array, dimension (N)   
            On entry, the n diagonal elements of the tridiagonal matrix. 
  
            On exit, if INFO = 0, the eigenvalues in ascending order.   

    E       (input/output) DOUBLE PRECISION array, dimension (N-1)   
            On entry, the (n-1) subdiagonal elements of the tridiagonal   
            matrix.   
            On exit, E has been destroyed.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   
            > 0:  the algorithm failed to find all of the eigenvalues in 
  
                  a total of 30*N iterations; if INFO = i, then i   
                  elements of E have not converged to zero.   

    ===================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static integer c__0 = 0;
    static integer c__1 = 1;
    static doublereal c_b32 = 1.;
    
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;
    /* Builtin functions */
    double sqrt(doublereal), d_sign(doublereal *, doublereal *);
    /* Local variables */
    static doublereal oldc;
    static integer lend, jtot;
    extern /* Subroutine */ int dlae2_(doublereal *, doublereal *, doublereal 
	    *, doublereal *, doublereal *);
    static doublereal c;
    static integer i, l, m;
    static doublereal p, gamma, r, s, alpha, sigma, anorm;
    static integer l1, lendm1, lendp1;
    extern doublereal dlapy2_(doublereal *, doublereal *);
    static doublereal bb;
    extern doublereal dlamch_(char *);
    static integer iscale;
    extern /* Subroutine */ int dlascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *);
    static doublereal oldgam, safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *);
    static doublereal safmax;
    extern doublereal dlanst_(char *, integer *, doublereal *, doublereal *);
    extern /* Subroutine */ int dlasrt_(char *, integer *, doublereal *, 
	    integer *);
    static integer lendsv;
    static doublereal ssfmin;
    static integer nmaxit;
    static doublereal ssfmax;
    static integer lm1, mm1, nm1;
    static doublereal rt1, rt2, eps, rte;
    static integer lsv;
    static doublereal tst, eps2;



#define E(I) e[(I)-1]
#define D(I) d[(I)-1]


    *info = 0;

/*     Quick return if possible */

    if (*n < 0) {
	*info = -1;
	i__1 = -(*info);
	xerbla_("DSTERF", &i__1);
	return 0;
    }
    if (*n <= 1) {
	return 0;
    }

/*     Determine the unit roundoff for this environment. */

    eps = dlamch_("E");
/* Computing 2nd power */
    d__1 = eps;
    eps2 = d__1 * d__1;
    safmin = dlamch_("S");
    safmax = 1. / safmin;
    ssfmax = sqrt(safmax) / 3.;
    ssfmin = sqrt(safmin) / eps2;

/*     Compute the eigenvalues of the tridiagonal matrix. */

    nmaxit = *n * 30;
    sigma = 0.;
    jtot = 0;

/*     Determine where the matrix splits and choose QL or QR iteration   
       for each block, according to whether top or bottom diagonal   
       element is smaller. */

    l1 = 1;
    nm1 = *n - 1;

L10:
    if (l1 > *n) {
	goto L170;
    }
    if (l1 > 1) {
	E(l1 - 1) = 0.;
    }
    if (l1 <= nm1) {
	i__1 = nm1;
	for (m = l1; m <= nm1; ++m) {
	    tst = (d__1 = E(m), abs(d__1));
	    if (tst == 0.) {
		goto L30;
	    }
	    if (tst <= sqrt((d__1 = D(m), abs(d__1))) * sqrt((d__2 = D(m + 1),
		     abs(d__2))) * eps) {
		E(m) = 0.;
		goto L30;
	    }
/* L20: */
	}
    }
    m = *n;

L30:
    l = l1;
    lsv = l;
    lend = m;
    lendsv = lend;
    l1 = m + 1;
    if (lend == l) {
	goto L10;
    }

/*     Scale submatrix in rows and columns L to LEND */

    i__1 = lend - l + 1;
    anorm = dlanst_("I", &i__1, &D(l), &E(l));
    iscale = 0;
    if (anorm > ssfmax) {
	iscale = 1;
	i__1 = lend - l + 1;
	dlascl_("G", &c__0, &c__0, &anorm, &ssfmax, &i__1, &c__1, &D(l), n, 
		info);
	i__1 = lend - l;
	dlascl_("G", &c__0, &c__0, &anorm, &ssfmax, &i__1, &c__1, &E(l), n, 
		info);
    } else if (anorm < ssfmin) {
	iscale = 2;
	i__1 = lend - l + 1;
	dlascl_("G", &c__0, &c__0, &anorm, &ssfmin, &i__1, &c__1, &D(l), n, 
		info);
	i__1 = lend - l;
	dlascl_("G", &c__0, &c__0, &anorm, &ssfmin, &i__1, &c__1, &E(l), n, 
		info);
    }

    i__1 = lend - 1;
    for (i = l; i <= lend-1; ++i) {
/* Computing 2nd power */
	d__1 = E(i);
	E(i) = d__1 * d__1;
/* L40: */
    }

/*     Choose between QL and QR iteration */

    if ((d__1 = D(lend), abs(d__1)) < (d__2 = D(l), abs(d__2))) {
	lend = lsv;
	l = lendsv;
    }

    if (lend >= l) {

/*        QL Iteration   

          Look for small subdiagonal element. */

L50:
	if (l != lend) {
	    lendm1 = lend - 1;
	    i__1 = lendm1;
	    for (m = l; m <= lendm1; ++m) {
		tst = (d__1 = E(m), abs(d__1));
		if (tst <= eps2 * (d__1 = D(m) * D(m + 1), abs(d__1))) {
		    goto L70;
		}
/* L60: */
	    }
	}

	m = lend;

L70:
	if (m < lend) {
	    E(m) = 0.;
	}
	p = D(l);
	if (m == l) {
	    goto L90;
	}

/*        If remaining matrix is 2 by 2, use DLAE2 to compute its   
          eigenvalues. */

	if (m == l + 1) {
	    rte = sqrt(E(l));
	    dlae2_(&D(l), &rte, &D(l + 1), &rt1, &rt2);
	    D(l) = rt1;
	    D(l + 1) = rt2;
	    E(l) = 0.;
	    l += 2;
	    if (l <= lend) {
		goto L50;
	    }
	    goto L150;
	}

	if (jtot == nmaxit) {
	    goto L150;
	}
	++jtot;

/*        Form shift. */

	rte = sqrt(E(l));
	sigma = (D(l + 1) - p) / (rte * 2.);
	r = dlapy2_(&sigma, &c_b32);
	sigma = p - rte / (sigma + d_sign(&r, &sigma));

	c = 1.;
	s = 0.;
	gamma = D(m) - sigma;
	p = gamma * gamma;

/*        Inner loop */

	mm1 = m - 1;
	i__1 = l;
	for (i = mm1; i >= l; --i) {
	    bb = E(i);
	    r = p + bb;
	    if (i != m - 1) {
		E(i + 1) = s * r;
	    }
	    oldc = c;
	    c = p / r;
	    s = bb / r;
	    oldgam = gamma;
	    alpha = D(i);
	    gamma = c * (alpha - sigma) - s * oldgam;
	    D(i + 1) = oldgam + (alpha - gamma);
	    if (c != 0.) {
		p = gamma * gamma / c;
	    } else {
		p = oldc * bb;
	    }
/* L80: */
	}

	E(l) = s * p;
	D(l) = sigma + gamma;
	goto L50;

/*        Eigenvalue found. */

L90:
	D(l) = p;

	++l;
	if (l <= lend) {
	    goto L50;
	}
	goto L150;

    } else {

/*        QR Iteration   

          Look for small superdiagonal element. */

L100:
	if (l != lend) {
	    lendp1 = lend + 1;
	    i__1 = lendp1;
	    for (m = l; m >= lendp1; --m) {
		tst = (d__1 = E(m - 1), abs(d__1));
		if (tst <= eps2 * (d__1 = D(m) * D(m - 1), abs(d__1))) {
		    goto L120;
		}
/* L110: */
	    }
	}

	m = lend;

L120:
	if (m > lend) {
	    E(m - 1) = 0.;
	}
	p = D(l);
	if (m == l) {
	    goto L140;
	}

/*        If remaining matrix is 2 by 2, use DLAE2 to compute its   
          eigenvalues. */

	if (m == l - 1) {
	    rte = sqrt(E(l - 1));
	    dlae2_(&D(l), &rte, &D(l - 1), &rt1, &rt2);
	    D(l) = rt1;
	    D(l - 1) = rt2;
	    E(l - 1) = 0.;
	    l += -2;
	    if (l >= lend) {
		goto L100;
	    }
	    goto L150;
	}

	if (jtot == nmaxit) {
	    goto L150;
	}
	++jtot;

/*        Form shift. */

	rte = sqrt(E(l - 1));
	sigma = (D(l - 1) - p) / (rte * 2.);
	r = dlapy2_(&sigma, &c_b32);
	sigma = p - rte / (sigma + d_sign(&r, &sigma));

	c = 1.;
	s = 0.;
	gamma = D(m) - sigma;
	p = gamma * gamma;

/*        Inner loop */

	lm1 = l - 1;
	i__1 = lm1;
	for (i = m; i <= lm1; ++i) {
	    bb = E(i);
	    r = p + bb;
	    if (i != m) {
		E(i - 1) = s * r;
	    }
	    oldc = c;
	    c = p / r;
	    s = bb / r;
	    oldgam = gamma;
	    alpha = D(i + 1);
	    gamma = c * (alpha - sigma) - s * oldgam;
	    D(i) = oldgam + (alpha - gamma);
	    if (c != 0.) {
		p = gamma * gamma / c;
	    } else {
		p = oldc * bb;
	    }
/* L130: */
	}

	E(lm1) = s * p;
	D(l) = sigma + gamma;
	goto L100;

/*        Eigenvalue found. */

L140:
	D(l) = p;

	--l;
	if (l >= lend) {
	    goto L100;
	}
	goto L150;

    }

/*     Undo scaling if necessary */

L150:
    if (iscale == 1) {
	i__1 = lendsv - lsv + 1;
	dlascl_("G", &c__0, &c__0, &ssfmax, &anorm, &i__1, &c__1, &D(lsv), n, 
		info);
    }
    if (iscale == 2) {
	i__1 = lendsv - lsv + 1;
	dlascl_("G", &c__0, &c__0, &ssfmin, &anorm, &i__1, &c__1, &D(lsv), n, 
		info);
    }

/*     Check for no convergence to an eigenvalue after a total   
       of N*MAXIT iterations. */

    if (jtot == nmaxit) {
	i__1 = *n - 1;
	for (i = 1; i <= *n-1; ++i) {
	    if (E(i) != 0.) {
		++(*info);
	    }
/* L160: */
	}
	return 0;
    }
    goto L10;

/*     Sort eigenvalues in increasing order. */

L170:
    dlasrt_("I", n, &D(1), info);

    return 0;

/*     End of DSTERF */

} /* dsterf_ */

