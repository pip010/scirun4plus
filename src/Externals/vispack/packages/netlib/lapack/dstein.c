#include "f2c.h"

/* Subroutine */ int dstein_(integer *n, doublereal *d, doublereal *e, 
	integer *m, doublereal *w, integer *iblock, integer *isplit, 
	doublereal *z, integer *ldz, doublereal *work, integer *iwork, 
	integer *ifail, integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DSTEIN computes the eigenvectors of a real symmetric tridiagonal   
    matrix T corresponding to specified eigenvalues, using inverse   
    iteration.   

    The maximum number of iterations allowed for each eigenvector is   
    specified by an internal parameter MAXITS (currently set to 5).   

    Arguments   
    =========   

    N       (input) INTEGER   
            The order of the matrix.  N >= 0.   

    D       (input) DOUBLE PRECISION array, dimension (N)   
            The n diagonal elements of the tridiagonal matrix T.   

    E       (input) DOUBLE PRECISION array, dimension (N)   
            The (n-1) subdiagonal elements of the tridiagonal matrix   
            T, in elements 1 to N-1.  E(N) need not be set.   

    M       (input) INTEGER   
            The number of eigenvectors to be found.  0 <= M <= N.   

    W       (input) DOUBLE PRECISION array, dimension (N)   
            The first M elements of W contain the eigenvalues for   
            which eigenvectors are to be computed.  The eigenvalues   
            should be grouped by split-off block and ordered from   
            smallest to largest within the block.  ( The output array   
            W from DSTEBZ with ORDER = 'B' is expected here. )   

    IBLOCK  (input) INTEGER array, dimension (N)   
            The submatrix indices associated with the corresponding   
            eigenvalues in W; IBLOCK(i)=1 if eigenvalue W(i) belongs to   
            the first submatrix from the top, =2 if W(i) belongs to   
            the second submatrix, etc.  ( The output array IBLOCK   
            from DSTEBZ is expected here. )   

    ISPLIT  (input) INTEGER array, dimension (N)   
            The splitting points, at which T breaks up into submatrices. 
  
            The first submatrix consists of rows/columns 1 to   
            ISPLIT( 1 ), the second of rows/columns ISPLIT( 1 )+1   
            through ISPLIT( 2 ), etc.   
            ( The output array ISPLIT from DSTEBZ is expected here. )   

    Z       (output) DOUBLE PRECISION array, dimension (LDZ, M)   
            The computed eigenvectors.  The eigenvector associated   
            with the eigenvalue W(i) is stored in the i-th column of   
            Z.  Any vector which fails to converge is set to its current 
  
            iterate after MAXITS iterations.   

    LDZ     (input) INTEGER   
            The leading dimension of the array Z.  LDZ >= max(1,N).   

    WORK    (workspace) DOUBLE PRECISION array, dimension (5*N)   

    IWORK   (workspace) INTEGER array, dimension (N)   

    IFAIL   (output) INTEGER array, dimension (M)   
            On normal exit, all elements of IFAIL are zero.   
            If one or more eigenvectors fail to converge after   
            MAXITS iterations, then their indices are stored in   
            array IFAIL.   

    INFO    (output) INTEGER   
            = 0: successful exit.   
            < 0: if INFO = -i, the i-th argument had an illegal value   
            > 0: if INFO = i, then i eigenvectors failed to converge   
                 in MAXITS iterations.  Their indices are stored in   
                 array IFAIL.   

    Internal Parameters   
    ===================   

    MAXITS  INTEGER, default = 5   
            The maximum number of iterations performed.   

    EXTRA   INTEGER, default = 2   
            The number of iterations performed after norm growth   
            criterion is satisfied, should be at least 1.   

    ===================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static integer c__2 = 2;
    static integer c__1 = 1;
    static integer c_n1 = -1;
    
    /* System generated locals */
    integer z_dim1, z_offset, i__1, i__2, i__3;
    doublereal d__1, d__2, d__3, d__4, d__5;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    static integer jblk, nblk;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static integer jmax;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    static integer i, j;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static integer iseed[4], gpind, iinfo;
    extern doublereal dasum_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static integer b1;
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);
    static integer j1;
    static doublereal ortol;
    static integer indrv1, indrv2, indrv3, indrv4, indrv5, bn;
    extern doublereal dlamch_(char *);
    extern /* Subroutine */ int dlagtf_(integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *, integer *
	    , integer *);
    static doublereal xj;
    extern integer idamax_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int xerbla_(char *, integer *), dlagts_(
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *);
    static integer nrmchk;
    extern /* Subroutine */ int dlarnv_(integer *, integer *, integer *, 
	    doublereal *);
    static integer blksiz;
    static doublereal onenrm, dtpcrt, pertol, scl, eps, sep, nrm, tol;
    static integer its;
    static doublereal xjm, ztr, eps1;



#define ISEED(I) iseed[(I)]
#define D(I) d[(I)-1]
#define E(I) e[(I)-1]
#define W(I) w[(I)-1]
#define IBLOCK(I) iblock[(I)-1]
#define ISPLIT(I) isplit[(I)-1]
#define WORK(I) work[(I)-1]
#define IWORK(I) iwork[(I)-1]
#define IFAIL(I) ifail[(I)-1]

#define Z(I,J) z[(I)-1 + ((J)-1)* ( *ldz)]

    *info = 0;
    i__1 = *m;
    for (i = 1; i <= *m; ++i) {
	IFAIL(i) = 0;
/* L10: */
    }

    if (*n < 0) {
	*info = -1;
    } else if (*m < 0 || *m > *n) {
	*info = -4;
    } else if (*ldz < max(1,*n)) {
	*info = -9;
    } else {
	i__1 = *m;
	for (j = 2; j <= *m; ++j) {
	    if (IBLOCK(j) < IBLOCK(j - 1)) {
		*info = -6;
		goto L30;
	    }
	    if (IBLOCK(j) == IBLOCK(j - 1) && W(j) < W(j - 1)) {
		*info = -5;
		goto L30;
	    }
/* L20: */
	}
L30:
	;
    }

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DSTEIN", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0 || *m == 0) {
	return 0;
    } else if (*n == 1) {
	Z(1,1) = 1.;
	return 0;
    }

/*     Get machine constants. */

    eps = dlamch_("Precision");

/*     Initialize seed for random number generator DLARNV. */

    for (i = 1; i <= 4; ++i) {
	ISEED(i - 1) = 1;
/* L40: */
    }

/*     Initialize pointers. */

    indrv1 = 0;
    indrv2 = indrv1 + *n;
    indrv3 = indrv2 + *n;
    indrv4 = indrv3 + *n;
    indrv5 = indrv4 + *n;

/*     Compute eigenvectors of matrix blocks. */

    j1 = 1;
    i__1 = IBLOCK(*m);
    for (nblk = 1; nblk <= IBLOCK(*m); ++nblk) {

/*        Find starting and ending indices of block nblk. */

	if (nblk == 1) {
	    b1 = 1;
	} else {
	    b1 = ISPLIT(nblk - 1) + 1;
	}
	bn = ISPLIT(nblk);
	blksiz = bn - b1 + 1;
	if (blksiz == 1) {
	    goto L60;
	}
	gpind = b1;

/*        Compute reorthogonalization criterion and stopping criterion
. */

	onenrm = (d__1 = D(b1), abs(d__1)) + (d__2 = E(b1), abs(d__2));
/* Computing MAX */
	d__3 = onenrm, d__4 = (d__1 = D(bn), abs(d__1)) + (d__2 = E(bn - 1), 
		abs(d__2));
	onenrm = max(d__3,d__4);
	i__2 = bn - 1;
	for (i = b1 + 1; i <= bn-1; ++i) {
/* Computing MAX */
	    d__4 = onenrm, d__5 = (d__1 = D(i), abs(d__1)) + (d__2 = E(i - 1),
		     abs(d__2)) + (d__3 = E(i), abs(d__3));
	    onenrm = max(d__4,d__5);
/* L50: */
	}
	ortol = onenrm * .001;

	dtpcrt = sqrt(.1 / blksiz);

/*        Loop through eigenvalues of block nblk. */

L60:
	jblk = 0;
	i__2 = *m;
	for (j = j1; j <= *m; ++j) {
	    if (IBLOCK(j) != nblk) {
		j1 = j;
		goto L160;
	    }
	    ++jblk;
	    xj = W(j);

/*           Skip all the work if the block size is one. */

	    if (blksiz == 1) {
		WORK(indrv1 + 1) = 1.;
		goto L120;
	    }

/*           If eigenvalues j and j-1 are too close, add a relativ
ely   
             small perturbation. */

	    if (jblk > 1) {
		eps1 = (d__1 = eps * xj, abs(d__1));
		pertol = eps1 * 10.;
		sep = xj - xjm;
		if (sep < pertol) {
		    xj = xjm + pertol;
		}
	    }

	    its = 0;
	    nrmchk = 0;

/*           Get random starting vector. */

	    dlarnv_(&c__2, iseed, &blksiz, &WORK(indrv1 + 1));

/*           Copy the matrix T so it won't be destroyed in factori
zation. */

	    dcopy_(&blksiz, &D(b1), &c__1, &WORK(indrv4 + 1), &c__1);
	    i__3 = blksiz - 1;
	    dcopy_(&i__3, &E(b1), &c__1, &WORK(indrv2 + 2), &c__1);
	    i__3 = blksiz - 1;
	    dcopy_(&i__3, &E(b1), &c__1, &WORK(indrv3 + 1), &c__1);

/*           Compute LU factors with partial pivoting  ( PT = LU )
 */

	    tol = 0.;
	    dlagtf_(&blksiz, &WORK(indrv4 + 1), &xj, &WORK(indrv2 + 2), &WORK(
		    indrv3 + 1), &tol, &WORK(indrv5 + 1), &IWORK(1), &iinfo);

/*           Update iteration count. */

L70:
	    ++its;
	    if (its > 5) {
		goto L100;
	    }

/*           Normalize and scale the righthand side vector Pb.   

   Computing MAX */
	    d__2 = eps, d__3 = (d__1 = WORK(indrv4 + blksiz), abs(d__1));
	    scl = blksiz * onenrm * max(d__2,d__3) / dasum_(&blksiz, &WORK(
		    indrv1 + 1), &c__1);
	    dscal_(&blksiz, &scl, &WORK(indrv1 + 1), &c__1);

/*           Solve the system LU = Pb. */

	    dlagts_(&c_n1, &blksiz, &WORK(indrv4 + 1), &WORK(indrv2 + 2), &
		    WORK(indrv3 + 1), &WORK(indrv5 + 1), &IWORK(1), &WORK(
		    indrv1 + 1), &tol, &iinfo);

/*           Reorthogonalize by modified Gram-Schmidt if eigenvalu
es are   
             close enough. */

	    if (jblk == 1) {
		goto L90;
	    }
	    if ((d__1 = xj - xjm, abs(d__1)) > ortol) {
		gpind = j;
	    }
	    if (gpind != j) {
		i__3 = j - 1;
		for (i = gpind; i <= j-1; ++i) {
		    ztr = -ddot_(&blksiz, &WORK(indrv1 + 1), &c__1, &Z(b1,i), &c__1);
		    daxpy_(&blksiz, &ztr, &Z(b1,i), &c__1, &WORK(
			    indrv1 + 1), &c__1);
/* L80: */
		}
	    }

/*           Check the infinity norm of the iterate. */

L90:
	    jmax = idamax_(&blksiz, &WORK(indrv1 + 1), &c__1);
	    nrm = (d__1 = WORK(indrv1 + jmax), abs(d__1));

/*           Continue for additional iterations after norm reaches
   
             stopping criterion. */

	    if (nrm < dtpcrt) {
		goto L70;
	    }
	    ++nrmchk;
	    if (nrmchk < 3) {
		goto L70;
	    }

	    goto L110;

/*           If stopping criterion was not satisfied, update info 
and   
             store eigenvector number in array ifail. */

L100:
	    ++(*info);
	    IFAIL(*info) = j;

/*           Accept iterate as jth eigenvector. */

L110:
	    scl = 1. / dnrm2_(&blksiz, &WORK(indrv1 + 1), &c__1);
	    jmax = idamax_(&blksiz, &WORK(indrv1 + 1), &c__1);
	    if (WORK(indrv1 + jmax) < 0.) {
		scl = -scl;
	    }
	    dscal_(&blksiz, &scl, &WORK(indrv1 + 1), &c__1);
L120:
	    i__3 = *n;
	    for (i = 1; i <= *n; ++i) {
		Z(i,j) = 0.;
/* L130: */
	    }
	    i__3 = blksiz;
	    for (i = 1; i <= blksiz; ++i) {
		Z(b1+i-1,j) = WORK(indrv1 + i);
/* L140: */
	    }

/*           Save the shift to check eigenvalue spacing at next   
             iteration. */

	    xjm = xj;

/* L150: */
	}
L160:
	;
    }

    return 0;

/*     End of DSTEIN */

} /* dstein_ */

