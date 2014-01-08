#include "f2c.h"

/* Subroutine */ int cstein_(integer *n, real *d, real *e, integer *m, real *
	w, integer *iblock, integer *isplit, complex *z, integer *ldz, real *
	work, integer *iwork, integer *ifail, integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    CSTEIN computes the eigenvectors of a real symmetric tridiagonal   
    matrix T corresponding to specified eigenvalues, using inverse   
    iteration.   

    The maximum number of iterations allowed for each eigenvector is   
    specified by an internal parameter MAXITS (currently set to 5).   

    Although the eigenvectors are real, they are stored in a complex   
    array, which may be passed to CUNMTR or CUPMTR for back   
    transformation to the eigenvectors of a complex Hermitian matrix   
    which was reduced to tridiagonal form.   


    Arguments   
    =========   

    N       (input) INTEGER   
            The order of the matrix.  N >= 0.   

    D       (input) REAL array, dimension (N)   
            The n diagonal elements of the tridiagonal matrix T.   

    E       (input) REAL array, dimension (N)   
            The (n-1) subdiagonal elements of the tridiagonal matrix   
            T, stored in elements 1 to N-1; E(N) need not be set.   

    M       (input) INTEGER   
            The number of eigenvectors to be found.  0 <= M <= N.   

    W       (input) REAL array, dimension (N)   
            The first M elements of W contain the eigenvalues for   
            which eigenvectors are to be computed.  The eigenvalues   
            should be grouped by split-off block and ordered from   
            smallest to largest within the block.  ( The output array   
            W from SSTEBZ with ORDER = 'B' is expected here. )   

    IBLOCK  (input) INTEGER array, dimension (N)   
            The submatrix indices associated with the corresponding   
            eigenvalues in W; IBLOCK(i)=1 if eigenvalue W(i) belongs to   
            the first submatrix from the top, =2 if W(i) belongs to   
            the second submatrix, etc.  ( The output array IBLOCK   
            from SSTEBZ is expected here. )   

    ISPLIT  (input) INTEGER array, dimension (N)   
            The splitting points, at which T breaks up into submatrices. 
  
            The first submatrix consists of rows/columns 1 to   
            ISPLIT( 1 ), the second of rows/columns ISPLIT( 1 )+1   
            through ISPLIT( 2 ), etc.   
            ( The output array ISPLIT from SSTEBZ is expected here. )   

    Z       (output) COMPLEX array, dimension (LDZ, M)   
            The computed eigenvectors.  The eigenvector associated   
            with the eigenvalue W(i) is stored in the i-th column of   
            Z.  Any vector which fails to converge is set to its current 
  
            iterate after MAXITS iterations.   
            The imaginary parts of the eigenvectors are set to zero.   

    LDZ     (input) INTEGER   
            The leading dimension of the array Z.  LDZ >= max(1,N).   

    WORK    (workspace) REAL array, dimension (5*N)   

    IWORK   (workspace) INTEGER array, dimension (N)   

    IFAIL   (output) INTEGER array, dimension (M)   
            On normal exit, all elements of IFAIL are zero.   
            If one or more eigenvectors fail to converge after   
            MAXITS iterations, then their indices are stored in   
            array IFAIL.   

    INFO    (output) INTEGER   
            = 0: successful exit   
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
    integer z_dim1, z_offset, i__1, i__2, i__3, i__4, i__5;
    real r__1, r__2, r__3, r__4, r__5;
    complex q__1;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    static integer jblk, nblk, jmax;
    extern doublereal snrm2_(integer *, real *, integer *);
    static integer i, j, iseed[4], gpind, iinfo;
    extern /* Subroutine */ int sscal_(integer *, real *, real *, integer *);
    static integer b1;
    extern doublereal sasum_(integer *, real *, integer *);
    static integer j1;
    extern /* Subroutine */ int scopy_(integer *, real *, integer *, real *, 
	    integer *);
    static real ortol;
    static integer indrv1, indrv2, indrv3, indrv4, indrv5, bn, jr;
    static real xj;
    extern doublereal slamch_(char *);
    extern /* Subroutine */ int xerbla_(char *, integer *), slagtf_(
	    integer *, real *, real *, real *, real *, real *, real *, 
	    integer *, integer *);
    static integer nrmchk;
    extern integer isamax_(integer *, real *, integer *);
    extern /* Subroutine */ int slagts_(integer *, integer *, real *, real *, 
	    real *, real *, integer *, real *, real *, integer *);
    static integer blksiz;
    static real onenrm, pertol;
    extern /* Subroutine */ int slarnv_(integer *, integer *, integer *, real 
	    *);
    static real stpcrt, scl, eps, ctr, sep, nrm, tol;
    static integer its;
    static real xjm, eps1;



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
	xerbla_("CSTEIN", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0 || *m == 0) {
	return 0;
    } else if (*n == 1) {
	i__1 = z_dim1 + 1;
	Z(1,1).r = 1.f, Z(1,1).i = 0.f;
	return 0;
    }

/*     Get machine constants. */

    eps = slamch_("Precision");

/*     Initialize seed for random number generator SLARNV. */

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

	onenrm = (r__1 = D(b1), dabs(r__1)) + (r__2 = E(b1), dabs(r__2));
/* Computing MAX */
	r__3 = onenrm, r__4 = (r__1 = D(bn), dabs(r__1)) + (r__2 = E(bn - 1), 
		dabs(r__2));
	onenrm = dmax(r__3,r__4);
	i__2 = bn - 1;
	for (i = b1 + 1; i <= bn-1; ++i) {
/* Computing MAX */
	    r__4 = onenrm, r__5 = (r__1 = D(i), dabs(r__1)) + (r__2 = E(i - 1)
		    , dabs(r__2)) + (r__3 = E(i), dabs(r__3));
	    onenrm = dmax(r__4,r__5);
/* L50: */
	}
	ortol = onenrm * .001f;

	stpcrt = sqrt(.1f / blksiz);

/*        Loop through eigenvalues of block nblk. */

L60:
	jblk = 0;
	i__2 = *m;
	for (j = j1; j <= *m; ++j) {
	    if (IBLOCK(j) != nblk) {
		j1 = j;
		goto L180;
	    }
	    ++jblk;
	    xj = W(j);

/*           Skip all the work if the block size is one. */

	    if (blksiz == 1) {
		WORK(indrv1 + 1) = 1.f;
		goto L140;
	    }

/*           If eigenvalues j and j-1 are too close, add a relativ
ely   
             small perturbation. */

	    if (jblk > 1) {
		eps1 = (r__1 = eps * xj, dabs(r__1));
		pertol = eps1 * 10.f;
		sep = xj - xjm;
		if (sep < pertol) {
		    xj = xjm + pertol;
		}
	    }

	    its = 0;
	    nrmchk = 0;

/*           Get random starting vector. */

	    slarnv_(&c__2, iseed, &blksiz, &WORK(indrv1 + 1));

/*           Copy the matrix T so it won't be destroyed in factori
zation. */

	    scopy_(&blksiz, &D(b1), &c__1, &WORK(indrv4 + 1), &c__1);
	    i__3 = blksiz - 1;
	    scopy_(&i__3, &E(b1), &c__1, &WORK(indrv2 + 2), &c__1);
	    i__3 = blksiz - 1;
	    scopy_(&i__3, &E(b1), &c__1, &WORK(indrv3 + 1), &c__1);

/*           Compute LU factors with partial pivoting  ( PT = LU )
 */

	    tol = 0.f;
	    slagtf_(&blksiz, &WORK(indrv4 + 1), &xj, &WORK(indrv2 + 2), &WORK(
		    indrv3 + 1), &tol, &WORK(indrv5 + 1), &IWORK(1), &iinfo);

/*           Update iteration count. */

L70:
	    ++its;
	    if (its > 5) {
		goto L120;
	    }

/*           Normalize and scale the righthand side vector Pb.   

   Computing MAX */
	    r__2 = eps, r__3 = (r__1 = WORK(indrv4 + blksiz), dabs(r__1));
	    scl = blksiz * onenrm * dmax(r__2,r__3) / sasum_(&blksiz, &WORK(
		    indrv1 + 1), &c__1);
	    sscal_(&blksiz, &scl, &WORK(indrv1 + 1), &c__1);

/*           Solve the system LU = Pb. */

	    slagts_(&c_n1, &blksiz, &WORK(indrv4 + 1), &WORK(indrv2 + 2), &
		    WORK(indrv3 + 1), &WORK(indrv5 + 1), &IWORK(1), &WORK(
		    indrv1 + 1), &tol, &iinfo);

/*           Reorthogonalize by modified Gram-Schmidt if eigenvalu
es are   
             close enough. */

	    if (jblk == 1) {
		goto L110;
	    }
	    if ((r__1 = xj - xjm, dabs(r__1)) > ortol) {
		gpind = j;
	    }
	    if (gpind != j) {
		i__3 = j - 1;
		for (i = gpind; i <= j-1; ++i) {
		    ctr = 0.f;
		    i__4 = blksiz;
		    for (jr = 1; jr <= blksiz; ++jr) {
			i__5 = b1 - 1 + jr + i * z_dim1;
			ctr += WORK(indrv1 + jr) * Z(b1-1+jr,i).r;
/* L80: */
		    }
		    i__4 = blksiz;
		    for (jr = 1; jr <= blksiz; ++jr) {
			i__5 = b1 - 1 + jr + i * z_dim1;
			WORK(indrv1 + jr) -= ctr * Z(b1-1+jr,i).r;
/* L90: */
		    }
/* L100: */
		}
	    }

/*           Check the infinity norm of the iterate. */

L110:
	    jmax = isamax_(&blksiz, &WORK(indrv1 + 1), &c__1);
	    nrm = (r__1 = WORK(indrv1 + jmax), dabs(r__1));

/*           Continue for additional iterations after norm reaches
   
             stopping criterion. */

	    if (nrm < stpcrt) {
		goto L70;
	    }
	    ++nrmchk;
	    if (nrmchk < 3) {
		goto L70;
	    }

	    goto L130;

/*           If stopping criterion was not satisfied, update info 
and   
             store eigenvector number in array ifail. */

L120:
	    ++(*info);
	    IFAIL(*info) = j;

/*           Accept iterate as jth eigenvector. */

L130:
	    scl = 1.f / snrm2_(&blksiz, &WORK(indrv1 + 1), &c__1);
	    jmax = isamax_(&blksiz, &WORK(indrv1 + 1), &c__1);
	    if (WORK(indrv1 + jmax) < 0.f) {
		scl = -(doublereal)scl;
	    }
	    sscal_(&blksiz, &scl, &WORK(indrv1 + 1), &c__1);
L140:
	    i__3 = *n;
	    for (i = 1; i <= *n; ++i) {
		i__4 = i + j * z_dim1;
		Z(i,j).r = 0.f, Z(i,j).i = 0.f;
/* L150: */
	    }
	    i__3 = blksiz;
	    for (i = 1; i <= blksiz; ++i) {
		i__4 = b1 + i - 1 + j * z_dim1;
		i__5 = indrv1 + i;
		q__1.r = WORK(indrv1+i), q__1.i = 0.f;
		Z(b1+i-1,j).r = q__1.r, Z(b1+i-1,j).i = q__1.i;
/* L160: */
	    }

/*           Save the shift to check eigenvalue spacing at next   
             iteration. */

	    xjm = xj;

/* L170: */
	}
L180:
	;
    }

    return 0;

/*     End of CSTEIN */

} /* cstein_ */

