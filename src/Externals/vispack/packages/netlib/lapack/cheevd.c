#include "f2c.h"

/* Subroutine */ int cheevd_(char *jobz, char *uplo, integer *n, complex *a, 
	integer *lda, real *w, complex *work, integer *lwork, real *rwork, 
	integer *lrwork, integer *iwork, integer *liwork, integer *info)
{
/*  -- LAPACK driver routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    CHEEVD computes all eigenvalues and, optionally, eigenvectors of a   
    complex Hermitian matrix A.  If eigenvectors are desired, it uses a   
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

    A       (input/output) COMPLEX array, dimension (LDA, N)   
            On entry, the Hermitian matrix A.  If UPLO = 'U', the   
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

    W       (output) REAL array, dimension (N)   
            If INFO = 0, the eigenvalues in ascending order.   

    WORK    (workspace/output) COMPLEX array, dimension (LWORK)   
            On exit, if LWORK > 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INTEGER   
            The length of the array WORK.   
            If N <= 1,                LWORK must be at least 1.   
            If JOBZ  = 'N' and N > 1, LWORK must be at least N + 1.   
            If JOBZ  = 'V' and N > 1, LWORK must be at least 2*N + N**2. 
  

    RWORK   (workspace/output) REAL array,   
                                           dimension (LRWORK)   
            On exit, if LRWORK > 0, RWORK(1) returns the optimal LRWORK. 
  

    LRWORK  (input) INTEGER   
            The dimension of the array RWORK.   
            If N <= 1,                LRWORK must be at least 1.   
            If JOBZ  = 'N' and N > 1, LRWORK must be at least N.   
            If JOBZ  = 'V' and N > 1, LRWORK must be at least   
                           1 + 4*N + 2*N*lg N + 3*N**2 ,   
                           where lg( N ) = smallest integer k such   
                           that 2**k >= N .   

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
    static real c_b16 = 1.f;
    static integer c__1 = 1;
    
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;
    real r__1, r__2;
    /* Builtin functions */
    double log(doublereal);
    integer pow_ii(integer *, integer *);
    double sqrt(doublereal);
    /* Local variables */
    static integer inde;
    static real anrm;
    static integer imax;
    static real rmin, rmax;
    static integer lopt;
    static real sigma;
    extern logical lsame_(char *, char *);
    static integer iinfo;
    extern /* Subroutine */ int sscal_(integer *, real *, real *, integer *);
    static integer lwmin, liopt;
    static logical lower;
    static integer llrwk, lropt;
    static logical wantz;
    static integer indwk2, llwrk2;
    extern doublereal clanhe_(char *, char *, integer *, complex *, integer *,
	     real *);
    static integer iscale;
    extern /* Subroutine */ int clascl_(char *, integer *, integer *, real *, 
	    real *, integer *, integer *, complex *, integer *, integer *), cstedc_(char *, integer *, real *, real *, complex *, 
	    integer *, complex *, integer *, real *, integer *, integer *, 
	    integer *, integer *);
    extern doublereal slamch_(char *);
    extern /* Subroutine */ int chetrd_(char *, integer *, complex *, integer 
	    *, real *, real *, complex *, complex *, integer *, integer *), clacpy_(char *, integer *, integer *, complex *, integer 
	    *, complex *, integer *);
    static real safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *);
    static real bignum;
    static integer indtau, indrwk, indwrk, liwmin;
    extern /* Subroutine */ int ssterf_(integer *, real *, real *, integer *);
    static integer lrwmin;
    extern /* Subroutine */ int cunmtr_(char *, char *, char *, integer *, 
	    integer *, complex *, integer *, complex *, complex *, integer *, 
	    complex *, integer *, integer *);
    static integer llwork;
    static real smlnum;
    static integer lgn;
    static real eps;



#define W(I) w[(I)-1]
#define WORK(I) work[(I)-1]
#define RWORK(I) rwork[(I)-1]
#define IWORK(I) iwork[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    wantz = lsame_(jobz, "V");
    lower = lsame_(uplo, "L");

    *info = 0;
    if (*n <= 1) {
	lgn = 0;
	lwmin = 1;
	lrwmin = 1;
	liwmin = 1;
	lopt = lwmin;
	lropt = lrwmin;
	liopt = liwmin;
    } else {
	lgn = (integer) (log((real) (*n)) / log(2.f));
	if (pow_ii(&c__2, &lgn) < *n) {
	    ++lgn;
	}
	if (pow_ii(&c__2, &lgn) < *n) {
	    ++lgn;
	}
	if (wantz) {
	    lwmin = (*n << 1) + *n * *n;
/* Computing 2nd power */
	    i__1 = *n;
	    lrwmin = (*n << 2) + 1 + (*n << 1) * lgn + i__1 * i__1 * 3;
	    liwmin = *n * 5 + 2;
	} else {
	    lwmin = *n + 1;
	    lrwmin = *n;
	    liwmin = 1;
	}
	lopt = lwmin;
	lropt = lrwmin;
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
    } else if (*lrwork < lrwmin) {
	*info = -10;
    } else if (*liwork < liwmin) {
	*info = -12;
    }

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("CHEEVD ", &i__1);
	goto L10;
    }

/*     Quick return if possible */

    if (*n == 0) {
	goto L10;
    }

    if (*n == 1) {
	i__1 = a_dim1 + 1;
	W(1) = A(1,1).r;
	if (wantz) {
	    i__1 = a_dim1 + 1;
	    A(1,1).r = 1.f, A(1,1).i = 0.f;
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

    anrm = clanhe_("M", uplo, n, &A(1,1), lda, &RWORK(1));
    iscale = 0;
    if (anrm > 0.f && anrm < rmin) {
	iscale = 1;
	sigma = rmin / anrm;
    } else if (anrm > rmax) {
	iscale = 1;
	sigma = rmax / anrm;
    }
    if (iscale == 1) {
	clascl_(uplo, &c__0, &c__0, &c_b16, &sigma, n, n, &A(1,1), lda, 
		info);
    }

/*     Call CHETRD to reduce Hermitian matrix to tridiagonal form. */

    inde = 1;
    indtau = 1;
    indwrk = indtau + *n;
    indrwk = inde + *n;
    indwk2 = indwrk + *n * *n;
    llwork = *lwork - indwrk + 1;
    llwrk2 = *lwork - indwk2 + 1;
    llrwk = *lrwork - indrwk + 1;
    chetrd_(uplo, n, &A(1,1), lda, &W(1), &RWORK(inde), &WORK(indtau), &
	    WORK(indwrk), &llwork, &iinfo);
/* Computing MAX */
    i__1 = indwrk;
    r__1 = (real) lopt, r__2 = (real) (*n) + WORK(indwrk).r;
    lopt = dmax(r__1,r__2);

/*     For eigenvalues only, call SSTERF.  For eigenvectors, first call   
       CSTEDC to generate the eigenvector matrix, WORK(INDWRK), of the   
       tridiagonal matrix, then call CUNMTR to multiply it to the   
       Householder transformations represented as Householder vectors in 
  
       A. */

    if (! wantz) {
	ssterf_(n, &W(1), &RWORK(inde), info);
    } else {
	cstedc_("I", n, &W(1), &RWORK(inde), &WORK(indwrk), n, &WORK(indwk2), 
		&llwrk2, &RWORK(indrwk), &llrwk, &IWORK(1), liwork, info);
	cunmtr_("L", uplo, "N", n, n, &A(1,1), lda, &WORK(indtau), &WORK(
		indwrk), n, &WORK(indwk2), &llwrk2, &iinfo);
	clacpy_("A", n, n, &WORK(indwrk), n, &A(1,1), lda);
/* Computing MAX   
   Computing 2nd power */
	i__3 = *n;
	i__4 = indwk2;
	i__1 = lopt, i__2 = *n + i__3 * i__3 + (integer) WORK(indwk2).r;
	lopt = max(i__1,i__2);
    }

/*     If matrix was scaled, then rescale eigenvalues appropriately. */

    if (iscale == 1) {
	if (*info == 0) {
	    imax = *n;
	} else {
	    imax = *info - 1;
	}
	r__1 = 1.f / sigma;
	sscal_(&imax, &r__1, &W(1), &c__1);
    }

L10:
    if (*lwork > 0) {
	WORK(1).r = (real) lopt, WORK(1).i = 0.f;
    }
    if (*lrwork > 0) {
	RWORK(1) = (real) lropt;
    }
    if (*liwork > 0) {
	IWORK(1) = liopt;
    }
    return 0;

/*     End of CHEEVD */

} /* cheevd_ */

