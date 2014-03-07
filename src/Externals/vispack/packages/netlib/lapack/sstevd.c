#include "f2c.h"

/* Subroutine */ int sstevd_(char *jobz, integer *n, real *d, real *e, real *
	z, integer *ldz, real *work, integer *lwork, integer *iwork, integer *
	liwork, integer *info)
{
/*  -- LAPACK driver routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    SSTEVD computes all eigenvalues and, optionally, eigenvectors of a   
    real symmetric tridiagonal matrix. If eigenvectors are desired, it   
    uses a divide and conquer algorithm.   

    The divide and conquer algorithm makes very mild assumptions about   
    floating point arithmetic. It will work on machines with a guard   
    digit in add/subtract, or on those binary machines without guard   
    digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or   
    Cray-2. It could conceivably fail on hexadecimal or decimal machines 
  
    without guard digits, but we know of none.   

    Arguments   
    =========   

    JOBZ    (input) CHARACTER*1   
            = 'N':  Compute eigenvalues only;   
            = 'V':  Compute eigenvalues and eigenvectors.   

    N       (input) INTEGER   
            The order of the matrix.  N >= 0.   

    D       (input/output) REAL array, dimension (N)   
            On entry, the n diagonal elements of the tridiagonal matrix   
            A.   
            On exit, if INFO = 0, the eigenvalues in ascending order.   

    E       (input/output) REAL array, dimension (N)   
            On entry, the (n-1) subdiagonal elements of the tridiagonal   
            matrix A, stored in elements 1 to N-1 of E; E(N) need not   
            be set, but is used by the routine.   
            On exit, the contents of E are destroyed.   

    Z       (output) REAL array, dimension (LDZ, N)   
            If JOBZ = 'V', then if INFO = 0, Z contains the orthonormal   
            eigenvectors of the matrix A, with the i-th column of Z   
            holding the eigenvector associated with D(i).   
            If JOBZ = 'N', then Z is not referenced.   

    LDZ     (input) INTEGER   
            The leading dimension of the array Z.  LDZ >= 1, and if   
            JOBZ = 'V', LDZ >= max(1,N).   

    WORK    (workspace/output) REAL array,   
                                           dimension (LWORK)   
            On exit, if LWORK > 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INTEGER   
            The dimension of the array WORK.   
            If JOBZ  = 'N' or N <= 1 then LWORK must be at least 1.   
            If JOBZ  = 'V' and N > 1 then LWORK must be at least   
                           ( 1 + 3*N + 2*N*lg N + 2*N**2 ),   
                           where lg( N ) = smallest integer k such   
                                           that 2**k >= N.   

    IWORK   (workspace/output) INTEGER array, dimension (LIWORK)   
            On exit, if LIWORK > 0, IWORK(1) returns the optimal LIWORK. 
  

    LIWORK  (input) INTEGER   
            The dimension of the array IWORK.   
            If JOBZ  = 'N' or N <= 1 then LIWORK must be at least 1.   
            If JOBZ  = 'V' and N > 1 then LIWORK must be at least 2+5*N. 
  

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   
            > 0:  if INFO = i, the algorithm failed to converge; i   
                  off-diagonal elements of E did not converge to zero.   

    ===================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static integer c__2 = 2;
    static integer c__1 = 1;
    
    /* System generated locals */
    integer z_dim1, z_offset, i__1;
    real r__1;
    /* Builtin functions */
    double log(doublereal);
    integer pow_ii(integer *, integer *);
    double sqrt(doublereal);
    /* Local variables */
    static real rmin, rmax, tnrm, sigma;
    extern logical lsame_(char *, char *);
    extern /* Subroutine */ int sscal_(integer *, real *, real *, integer *);
    static integer lwmin;
    static logical wantz;
    static integer iscale;
    extern doublereal slamch_(char *);
    static real safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *);
    static real bignum;
    extern /* Subroutine */ int sstedc_(char *, integer *, real *, real *, 
	    real *, integer *, real *, integer *, integer *, integer *, 
	    integer *);
    static integer liwmin;
    extern doublereal slanst_(char *, integer *, real *, real *);
    extern /* Subroutine */ int ssterf_(integer *, real *, real *, integer *);
    static real smlnum;
    static integer lgn;
    static real eps;



#define D(I) d[(I)-1]
#define E(I) e[(I)-1]
#define WORK(I) work[(I)-1]
#define IWORK(I) iwork[(I)-1]

#define Z(I,J) z[(I)-1 + ((J)-1)* ( *ldz)]

    wantz = lsame_(jobz, "V");

    *info = 0;
    liwmin = 1;
    lwmin = 1;
    if (! (wantz || lsame_(jobz, "N"))) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*ldz < 1 || wantz && *ldz < *n) {
	*info = -6;
    } else if (*n > 1 && wantz) {
	lgn = (integer) (log((real) (*n)) / log(2.f));
	if (pow_ii(&c__2, &lgn) < *n) {
	    ++lgn;
	}
	if (pow_ii(&c__2, &lgn) < *n) {
	    ++lgn;
	}
/* Computing 2nd power */
	i__1 = *n;
	lwmin = *n * 3 + 1 + (*n << 1) * lgn + (i__1 * i__1 << 1);
	liwmin = *n * 5 + 2;
	if (*lwork < lwmin) {
	    *info = -8;
	} else if (*liwork < liwmin) {
	    *info = -10;
	}
    } else if (*lwork < 1) {
	*info = -8;
    } else if (*liwork < 1) {
	*info = -10;
    }

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("SSTEVD", &i__1);
	goto L10;
    }

/*     Quick return if possible */

    if (*n == 0) {
	goto L10;
    }

    if (*n == 1) {
	if (wantz) {
	    Z(1,1) = 1.f;
	}
	goto L10;
    }

/*     Get machine constants. */

    safmin = slamch_("Safe minimum");
    eps = slamch_("Precision");
    smlnum = safmin / eps;
    bignum = 1.f / smlnum;
    rmin = sqrt(smlnum);
    rmax = sqrt(bignum);

/*     Scale matrix to allowable range, if necessary. */

    iscale = 0;
    tnrm = slanst_("M", n, &D(1), &E(1));
    if (tnrm > 0.f && tnrm < rmin) {
	iscale = 1;
	sigma = rmin / tnrm;
    } else if (tnrm > rmax) {
	iscale = 1;
	sigma = rmax / tnrm;
    }
    if (iscale == 1) {
	sscal_(n, &sigma, &D(1), &c__1);
	i__1 = *n - 1;
	sscal_(&i__1, &sigma, &E(1), &c__1);
    }

/*     For eigenvalues only, call SSTERF.  For eigenvalues and   
       eigenvectors, call SSTEDC. */

    if (! wantz) {
	ssterf_(n, &D(1), &E(1), info);
    } else {
	sstedc_("I", n, &D(1), &E(1), &Z(1,1), ldz, &WORK(1), lwork, &
		IWORK(1), liwork, info);
    }

/*     If matrix was scaled, then rescale eigenvalues appropriately. */

    if (iscale == 1) {
	r__1 = 1.f / sigma;
	sscal_(n, &r__1, &D(1), &c__1);
    }

L10:
    if (*lwork > 0) {
	WORK(1) = (real) lwmin;
    }
    if (*liwork > 0) {
	IWORK(1) = liwmin;
    }
    return 0;

/*     End of SSTEVD */

} /* sstevd_ */

