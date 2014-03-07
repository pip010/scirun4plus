#include "f2c.h"

/* Subroutine */ int dsyevd_(char *jobz, char *uplo, integer *n, doublereal *
	a, integer *lda, doublereal *w, doublereal *work, integer *lwork, 
	integer *iwork, integer *liwork, integer *info)
{
/*  -- LAPACK driver routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DSYEVD computes all eigenvalues and, optionally, eigenvectors of a   
    real symmetric matrix A. If eigenvectors are desired, it uses a   
    divide and conquer algorithm.   

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

    UPLO    (input) CHARACTER*1   
            = 'U':  Upper triangle of A is stored;   
            = 'L':  Lower triangle of A is stored.   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    A       (input/output) DOUBLE PRECISION array, dimension (LDA, N)   
            On entry, the symmetric matrix A.  If UPLO = 'U', the   
            leading N-by-N upper triangular part of A contains the   
            upper triangular part of the matrix A.  If UPLO = 'L',   
            the leading N-by-N lower triangular part of A contains   
            the lower triangular part of the matrix A.   
            On exit, if JOBZ = 'V', then if INFO = 0, A contains the   
            orthonormal eigenvectors of the matrix A.   
            If JOBZ = 'N', then on exit the lower triangle (if UPLO='L') 
  
            or the upper triangle (if UPLO='U') of A, including the   
            diagonal, is destroyed.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(1,N).   

    W       (output) DOUBLE PRECISION array, dimension (N)   
            If INFO = 0, the eigenvalues in ascending order.   

    WORK    (workspace/output) DOUBLE PRECISION array,   
                                           dimension (LWORK)   
            On exit, if LWORK > 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INTEGER   
            The dimension of the array WORK.   
            If N <= 1,               LWORK must be at least 1.   
            If JOBZ = 'N' and N > 1, LWORK must be at least 2*N+1.   
            If JOBZ = 'V' and N > 1, LWORK must be at least   
                           1 + 5*N + 2*N*lg N + 3*N**2,   
                           where lg( N ) = smallest integer k such   
                                           that 2**k >= N.   

    IWORK   (workspace/output) INTEGER array, dimension (LIWORK)   
            On exit, if LIWORK > 0, IWORK(1) returns the optimal LIWORK. 
  

    LIWORK  (input) INTEGER   
            The dimension of the array IWORK.   
            If N <= 1,                LIWORK must be at least 1.   
            If JOBZ  = 'N' and N > 1, LIWORK must be at least 1.   
            If JOBZ  = 'V' and N > 1, LIWORK must be at least 2 + 5*N.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   
            > 0:  if INFO = i, the algorithm failed to converge; i   
                  off-diagonal elements of an intermediate tridiagonal   
                  form did not converge to zero.   

    ===================================================================== 
  



       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static integer c__2 = 2;
    static integer c__0 = 0;
    static doublereal c_b15 = 1.;
    static integer c__1 = 1;
    
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublereal d__1;
    /* Builtin functions */
    double log(doublereal);
    integer pow_ii(integer *, integer *);
    double sqrt(doublereal);
    /* Local variables */
    static integer inde;
    static doublereal anrm, rmin, rmax;
    static integer lopt;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static doublereal sigma;
    extern logical lsame_(char *, char *);
    static integer iinfo, lwmin, liopt;
    static logical lower, wantz;
    static integer indwk2, llwrk2;
    extern doublereal dlamch_(char *);
    static integer iscale;
    extern /* Subroutine */ int dlascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *), dstedc_(char *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     integer *, integer *, integer *, integer *), dlacpy_(
	    char *, integer *, integer *, doublereal *, integer *, doublereal 
	    *, integer *);
    static doublereal safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *);
    static doublereal bignum;
    static integer indtau;
    extern /* Subroutine */ int dsterf_(integer *, doublereal *, doublereal *,
	     integer *);
    extern doublereal dlansy_(char *, char *, integer *, doublereal *, 
	    integer *, doublereal *);
    static integer indwrk, liwmin;
    extern /* Subroutine */ int dormtr_(char *, char *, char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *), dsytrd_(char *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *);
    static integer llwork;
    static doublereal smlnum;
    static integer lgn;
    static doublereal eps;



#define W(I) w[(I)-1]
#define WORK(I) work[(I)-1]
#define IWORK(I) iwork[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    wantz = lsame_(jobz, "V");
    lower = lsame_(uplo, "L");

    *info = 0;
    if (*n <= 1) {
	lgn = 0;
	liwmin = 1;
	lwmin = 1;
	lopt = lwmin;
	liopt = liwmin;
    } else {
	lgn = (integer) (log((doublereal) (*n)) / log(2.));
	if (pow_ii(&c__2, &lgn) < *n) {
	    ++lgn;
	}
	if (pow_ii(&c__2, &lgn) < *n) {
	    ++lgn;
	}
	if (wantz) {
	    liwmin = *n * 5 + 2;
/* Computing 2nd power */
	    i__1 = *n;
	    lwmin = *n * 5 + 1 + (*n << 1) * lgn + i__1 * i__1 * 3;
	} else {
	    liwmin = 1;
	    lwmin = (*n << 1) + 1;
	}
	lopt = lwmin;
	liopt = liwmin;
    }
    if (! (wantz || lsame_(jobz, "N"))) {
	*info = -1;
    } else if (! (lower || lsame_(uplo, "U"))) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*lda < max(1,*n)) {
	*info = -5;
    } else if (*lwork < lwmin) {
	*info = -8;
    } else if (*liwork < liwmin) {
	*info = -10;
    }

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DSYEVD ", &i__1);
	goto L10;
    }

/*     Quick return if possible */

    if (*n == 0) {
	goto L10;
    }

    if (*n == 1) {
	W(1) = A(1,1);
	if (wantz) {
	    A(1,1) = 1.;
	}
	goto L10;
    }

/*     Get machine constants. */

    safmin = dlamch_("Safe minimum");
    eps = dlamch_("Precision");
    smlnum = safmin / eps;
    bignum = 1. / smlnum;
    rmin = sqrt(smlnum);
    rmax = sqrt(bignum);

/*     Scale matrix to allowable range, if necessary. */

    anrm = dlansy_("M", uplo, n, &A(1,1), lda, &WORK(1));
    iscale = 0;
    if (anrm > 0. && anrm < rmin) {
	iscale = 1;
	sigma = rmin / anrm;
    } else if (anrm > rmax) {
	iscale = 1;
	sigma = rmax / anrm;
    }
    if (iscale == 1) {
	dlascl_(uplo, &c__0, &c__0, &c_b15, &sigma, n, n, &A(1,1), lda, 
		info);
    }

/*     Call DSYTRD to reduce symmetric matrix to tridiagonal form. */

    inde = 1;
    indtau = inde + *n;
    indwrk = indtau + *n;
    llwork = *lwork - indwrk + 1;
    indwk2 = indwrk + *n * *n;
    llwrk2 = *lwork - indwk2 + 1;

    dsytrd_(uplo, n, &A(1,1), lda, &W(1), &WORK(inde), &WORK(indtau), &
	    WORK(indwrk), &llwork, &iinfo);
    lopt = (integer) ((*n << 1) + WORK(indwrk));

/*     For eigenvalues only, call DSTERF.  For eigenvectors, first call   
       DSTEDC to generate the eigenvector matrix, WORK(INDWRK), of the   
       tridiagonal matrix, then call DORMTR to multiply it by the   
       Householder transformations stored in A. */

    if (! wantz) {
	dsterf_(n, &W(1), &WORK(inde), info);
    } else {
	dstedc_("I", n, &W(1), &WORK(inde), &WORK(indwrk), n, &WORK(indwk2), &
		llwrk2, &IWORK(1), liwork, info);
	dormtr_("L", uplo, "N", n, n, &A(1,1), lda, &WORK(indtau), &WORK(
		indwrk), n, &WORK(indwk2), &llwrk2, &iinfo);
	dlacpy_("A", n, n, &WORK(indwrk), n, &A(1,1), lda);
/* Computing MAX   
   Computing 2nd power */
	i__3 = *n;
	i__1 = lopt, i__2 = *n * 5 + 1 + (*n << 1) * lgn + i__3 * i__3 * 3;
	lopt = max(i__1,i__2);
    }

/*     If matrix was scaled, then rescale eigenvalues appropriately. */

    if (iscale == 1) {
	d__1 = 1. / sigma;
	dscal_(n, &d__1, &W(1), &c__1);
    }
L10:
    if (*lwork > 0) {
	WORK(1) = (doublereal) lopt;
    }
    if (*liwork > 0) {
	IWORK(1) = liopt;
    }
    return 0;

/*     End of DSYEVD */

} /* dsyevd_ */

