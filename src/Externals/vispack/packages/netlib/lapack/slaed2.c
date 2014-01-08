#include "f2c.h"

/* Subroutine */ int slaed2_(integer *k, integer *n, real *d, real *q, 
	integer *ldq, integer *indxq, real *rho, integer *cutpnt, real *z, 
	real *dlamda, real *q2, integer *ldq2, integer *indxc, real *w, 
	integer *indxp, integer *indx, integer *coltyp, integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Oak Ridge National Lab, Argonne National Lab, 
  
       Courant Institute, NAG Ltd., and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    SLAED2 merges the two sets of eigenvalues together into a single   
    sorted set.  Then it tries to deflate the size of the problem.   
    There are two ways in which deflation can occur:  when two or more   
    eigenvalues are close together or if there is a tiny entry in the   
    Z vector.  For each such occurrence the order of the related secular 
  
    equation problem is reduced by one.   

    Arguments   
    =========   

    K      (output) INTEGER   
           The number of non-deflated eigenvalues, and the order of the   
           related secular equation. 0 <= K <=N.   

    N      (input) INTEGER   
           The dimension of the symmetric tridiagonal matrix.  N >= 0.   

    D      (input/output) REAL array, dimension (N)   
           On entry, D contains the eigenvalues of the two submatrices to 
  
           be combined.   
           On exit, D contains the trailing (N-K) updated eigenvalues   
           (those which were deflated) sorted into increasing order.   

    Q      (input/output) REAL array, dimension (LDQ, N)   
           On entry, Q contains the eigenvectors of two submatrices in   
           the two square blocks with corners at (1,1), (CUTPNT,CUTPNT)   
           and (CUTPNT+1, CUTPNT+1), (N,N).   
           On exit, Q contains the trailing (N-K) updated eigenvectors   
           (those which were deflated) in its last N-K columns.   

    LDQ    (input) INTEGER   
           The leading dimension of the array Q.  LDQ >= max(1,N).   

    INDXQ  (input/output) INTEGER array, dimension (N)   
           The permutation which separately sorts the two sub-problems   
           in D into ascending order.  Note that elements in the second   
           half of this permutation must first have CUTPNT added to their 
  
           values. Destroyed on exit.   

    RHO    (input/output) REAL   
           On entry, the off-diagonal element associated with the rank-1 
  
           cut which originally split the two submatrices which are now   
           being recombined.   
           On exit, RHO has been modified to the value required by   
           SLAED3.   

    CUTPNT (input) INTEGER   
           The location of the last eigenvalue in the leading sub-matrix. 
  
           min(1,N) <= CUTPNT <= N.   

    Z      (input) REAL array, dimension (N)   
           On entry, Z contains the updating vector (the last   
           row of the first sub-eigenvector matrix and the first row of   
           the second sub-eigenvector matrix).   
           On exit, the contents of Z have been destroyed by the updating 
  
           process.   

    DLAMDA (output) REAL array, dimension (N)   
           A copy of the first K eigenvalues which will be used by   
           SLAED3 to form the secular equation.   

    Q2     (output) REAL array, dimension (LDQ2, N)   
           A copy of the first K eigenvectors which will be used by   
           SLAED3 in a matrix multiply (SGEMM) to solve for the new   
           eigenvectors.   Q2 is arranged into three blocks.  The   
           first block contains non-zero elements only at and above   
           CUTPNT, the second contains non-zero elements only below   
           CUTPNT, and the third is dense.   

    LDQ2   (input) INTEGER   
           The leading dimension of the array Q2.  LDQ2 >= max(1,N).   

    INDXC  (output) INTEGER array, dimension (N)   
           The permutation used to arrange the columns of the deflated   
           Q matrix into three groups:  the first group contains non-zero 
  
           elements only at and above CUTPNT, the second contains   
           non-zero elements only below CUTPNT, and the third is dense.   

    W      (output) REAL array, dimension (N)   
           The first k values of the final deflation-altered z-vector   
           which will be passed to SLAED3.   

    INDXP  (workspace) INTEGER array, dimension (N)   
           The permutation used to place deflated values of D at the end 
  
           of the array.  INDXP(1:K) points to the nondeflated D-values   
           and INDXP(K+1:N) points to the deflated eigenvalues.   

    INDX   (workspace) INTEGER array, dimension (N)   
           The permutation used to sort the contents of D into ascending 
  
           order.   

    COLTYP (workspace/output) INTEGER array, dimension (N)   
           During execution, a label which will indicate which of the   
           following types a column in the Q2 matrix is:   
           1 : non-zero in the upper half only;   
           2 : non-zero in the lower half only;   
           3 : dense;   
           4 : deflated.   
           On exit, COLTYP(i) is the number of columns of type i,   
           for i=1 to 4 only.   

    INFO   (output) INTEGER   
            = 0:  successful exit.   
            < 0:  if INFO = -i, the i-th argument had an illegal value.   

    ===================================================================== 
  


       Test the input parameters.   

       Parameter adjustments */
    /* Table of constant values */
    static real c_b3 = -1.f;
    static integer c__1 = 1;
    
    /* System generated locals */
    integer q_dim1, q_offset, q2_dim1, q2_offset, i__1;
    real r__1, r__2, r__3, r__4;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    static integer jlam, imax, jmax, ctot[4];
    extern /* Subroutine */ int srot_(integer *, real *, integer *, real *, 
	    integer *, real *, real *);
    static real c;
    static integer i, j;
    static real s, t;
    extern /* Subroutine */ int sscal_(integer *, real *, real *, integer *);
    static integer k2;
    extern /* Subroutine */ int scopy_(integer *, real *, integer *, real *, 
	    integer *);
    static integer n1, n2;
    extern doublereal slapy2_(real *, real *);
    static integer ct, jp;
    extern doublereal slamch_(char *);
    extern /* Subroutine */ int xerbla_(char *, integer *);
    extern integer isamax_(integer *, real *, integer *);
    extern /* Subroutine */ int slamrg_(integer *, integer *, real *, integer 
	    *, integer *, integer *), slacpy_(char *, integer *, integer *, 
	    real *, integer *, real *, integer *);
    static integer n1p1;
    static real eps, tau, tol;
    static integer psm[4];


    --d;
    q_dim1 = *ldq;
    q_offset = q_dim1 + 1;
    q -= q_offset;
    --indxq;
    --z;
    --dlamda;
    q2_dim1 = *ldq2;
    q2_offset = q2_dim1 + 1;
    q2 -= q2_offset;
    --indxc;
    --w;
    --indxp;
    --indx;
    --coltyp;

    /* Function Body */
    *info = 0;

    if (*n < 0) {
	*info = -2;
    } else if (*ldq < max(1,*n)) {
	*info = -5;
    } else if (min(1,*n) > *cutpnt || *n < *cutpnt) {
	*info = -8;
    } else if (*ldq2 < max(1,*n)) {
	*info = -12;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("SLAED2", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return 0;
    }

    n1 = *cutpnt;
    n2 = *n - n1;
    n1p1 = n1 + 1;

    if (*rho < 0.f) {
	sscal_(&n2, &c_b3, &z[n1p1], &c__1);
    }

/*     Normalize z so that norm(z) = 1.  Since z is the concatenation of 
  
       two normalized vectors, norm2(z) = sqrt(2). */

    t = 1.f / sqrt(2.f);
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	indx[j] = j;
/* L10: */
    }
    sscal_(n, &t, &z[1], &c__1);

/*     RHO = ABS( norm(z)**2 * RHO ) */

    *rho = (r__1 = *rho * 2.f, dabs(r__1));

    i__1 = *cutpnt;
    for (i = 1; i <= i__1; ++i) {
	coltyp[i] = 1;
/* L20: */
    }
    i__1 = *n;
    for (i = *cutpnt + 1; i <= i__1; ++i) {
	coltyp[i] = 2;
/* L30: */
    }

/*     Sort the eigenvalues into increasing order */

    i__1 = *n;
    for (i = *cutpnt + 1; i <= i__1; ++i) {
	indxq[i] += *cutpnt;
/* L40: */
    }

/*     re-integrate the deflated parts from the last pass */

    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	dlamda[i] = d[indxq[i]];
	w[i] = z[indxq[i]];
	indxc[i] = coltyp[indxq[i]];
/* L50: */
    }
    slamrg_(&n1, &n2, &dlamda[1], &c__1, &c__1, &indx[1]);
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	d[i] = dlamda[indx[i]];
	z[i] = w[indx[i]];
	coltyp[i] = indxc[indx[i]];
/* L60: */
    }

/*     Calculate the allowable deflation tolerance */

    imax = isamax_(n, &z[1], &c__1);
    jmax = isamax_(n, &d[1], &c__1);
    eps = slamch_("Epsilon");
/* Computing MAX */
    r__3 = (r__1 = d[jmax], dabs(r__1)), r__4 = (r__2 = z[imax], dabs(r__2));
    tol = eps * 8.f * dmax(r__3,r__4);

/*     If the rank-1 modifier is small enough, no more needs to be done   
       except to reorganize Q so that its columns correspond with the   
       elements in D. */

    if (*rho * (r__1 = z[imax], dabs(r__1)) <= tol) {
	*k = 0;
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    scopy_(n, &q[indxq[indx[j]] * q_dim1 + 1], &c__1, &q2[j * q2_dim1 
		    + 1], &c__1);
/* L70: */
	}
	slacpy_("A", n, n, &q2[q2_offset], ldq2, &q[q_offset], ldq);
	goto L180;
    }

/*     If there are multiple eigenvalues then the problem deflates.  Here 
  
       the number of equal eigenvalues are found.  As each equal   
       eigenvalue is found, an elementary reflector is computed to rotate 
  
       the corresponding eigensubspace so that the corresponding   
       components of Z are zero in this new basis. */

    *k = 0;
    k2 = *n + 1;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	if (*rho * (r__1 = z[j], dabs(r__1)) <= tol) {

/*           Deflate due to small z component. */

	    --k2;
	    indxp[k2] = j;
	    coltyp[j] = 4;
	    if (j == *n) {
		goto L120;
	    }
	} else {
	    jlam = j;
	    goto L90;
	}
/* L80: */
    }
L90:
    ++j;
    if (j > *n) {
	goto L110;
    }
    if (*rho * (r__1 = z[j], dabs(r__1)) <= tol) {

/*        Deflate due to small z component. */

	--k2;
	indxp[k2] = j;
	coltyp[j] = 4;
    } else {

/*        Check if eigenvalues are close enough to allow deflation. */

	s = z[jlam];
	c = z[j];

/*        Find sqrt(a**2+b**2) without overflow or   
          destructive underflow. */

	tau = slapy2_(&c, &s);
	t = d[j] - d[jlam];
	c /= tau;
	s = -(doublereal)s / tau;
	if ((r__1 = t * c * s, dabs(r__1)) <= tol) {

/*           Deflation is possible. */

	    z[j] = tau;
	    z[jlam] = 0.f;
	    if (coltyp[j] != coltyp[jlam]) {
		coltyp[j] = 3;
	    }
	    coltyp[jlam] = 4;
	    srot_(n, &q[indxq[indx[jlam]] * q_dim1 + 1], &c__1, &q[indxq[indx[
		    j]] * q_dim1 + 1], &c__1, &c, &s);
/* Computing 2nd power */
	    r__1 = c;
/* Computing 2nd power */
	    r__2 = s;
	    t = d[jlam] * (r__1 * r__1) + d[j] * (r__2 * r__2);
/* Computing 2nd power */
	    r__1 = s;
/* Computing 2nd power */
	    r__2 = c;
	    d[j] = d[jlam] * (r__1 * r__1) + d[j] * (r__2 * r__2);
	    d[jlam] = t;
	    --k2;
	    i = 1;
L100:
	    if (k2 + i <= *n) {
		if (d[jlam] < d[indxp[k2 + i]]) {
		    indxp[k2 + i - 1] = indxp[k2 + i];
		    indxp[k2 + i] = jlam;
		    ++i;
		    goto L100;
		} else {
		    indxp[k2 + i - 1] = jlam;
		}
	    } else {
		indxp[k2 + i - 1] = jlam;
	    }
	    jlam = j;
	} else {
	    ++(*k);
	    w[*k] = z[jlam];
	    dlamda[*k] = d[jlam];
	    indxp[*k] = jlam;
	    jlam = j;
	}
    }
    goto L90;
L110:

/*     Record the last eigenvalue. */

    ++(*k);
    w[*k] = z[jlam];
    dlamda[*k] = d[jlam];
    indxp[*k] = jlam;

L120:

/*     Count up the total number of the various types of columns, then   
       form a permutation which positions the four column types into   
       four uniform groups (although one or more of these groups may be   
       empty). */

    for (j = 1; j <= 4; ++j) {
	ctot[j - 1] = 0;
/* L130: */
    }
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	ct = coltyp[j];
	++ctot[ct - 1];
/* L140: */
    }

/*     PSM(*) = Position in SubVISMatrix (of types 1 through 4) */

    psm[0] = 1;
    psm[1] = ctot[0] + 1;
    psm[2] = psm[1] + ctot[1];
    psm[3] = psm[2] + ctot[2];

/*     Fill out the INDXC array so that the permutation which it induces 
  
       will place all type-1 columns first, all type-2 columns next,   
       then all type-3's, and finally all type-4's. */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	jp = indxp[j];
	ct = coltyp[jp];
	indxc[psm[ct - 1]] = j;
	++psm[ct - 1];
/* L150: */
    }

/*     Sort the eigenvalues and corresponding eigenvectors into DLAMDA   
       and Q2 respectively.  The eigenvalues/vectors which were not   
       deflated go into the first K slots of DLAMDA and Q2 respectively, 
  
       while those which were deflated go into the last N - K slots. */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	jp = indxp[j];
	dlamda[j] = d[jp];
	scopy_(n, &q[indxq[indx[indxp[indxc[j]]]] * q_dim1 + 1], &c__1, &q2[j 
		* q2_dim1 + 1], &c__1);
/* L160: */
    }

/*     The deflated eigenvalues and their corresponding vectors go back   
       into the last N - K slots of D and Q respectively. */

    i__1 = *n - *k;
    scopy_(&i__1, &dlamda[*k + 1], &c__1, &d[*k + 1], &c__1);
    i__1 = *n - *k;
    slacpy_("A", n, &i__1, &q2[(*k + 1) * q2_dim1 + 1], ldq2, &q[(*k + 1) * 
	    q_dim1 + 1], ldq);

/*     Copy CTOT into COLTYP for referencing in SLAED3. */

    for (j = 1; j <= 4; ++j) {
	coltyp[j] = ctot[j - 1];
/* L170: */
    }

L180:
    return 0;

/*     End of SLAED2 */

} /* slaed2_ */

