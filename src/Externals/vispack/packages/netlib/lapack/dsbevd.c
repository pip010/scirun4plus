#include "f2c.h"

/* Subroutine */ int dsbevd_(char *jobz, char *uplo, integer *n, integer *kd, 
	doublereal *ab, integer *ldab, doublereal *w, doublereal *z, integer *
	ldz, doublereal *work, integer *lwork, integer *iwork, integer *
	liwork, integer *info)
{
/*  -- LAPACK driver routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DSBEVD computes all the eigenvalues and, optionally, eigenvectors of 
  
    a real symmetric band matrix A. If eigenvectors are desired, it uses 
  
    a divide and conquer algorithm.   

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

    KD      (input) INTEGER   
            The number of superdiagonals of the matrix A if UPLO = 'U',   
            or the number of subdiagonals if UPLO = 'L'.  KD >= 0.   

    AB      (input/output) DOUBLE PRECISION array, dimension (LDAB, N)   
            On entry, the upper or lower triangle of the symmetric band   
            matrix A, stored in the first KD+1 rows of the array.  The   
            j-th column of A is stored in the j-th column of the array AB 
  
            as follows:   
            if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j; 
  
            if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd). 
  

            On exit, AB is overwritten by values generated during the   
            reduction to tridiagonal form.  If UPLO = 'U', the first   
            superdiagonal and the diagonal of the tridiagonal matrix T   
            are returned in rows KD and KD+1 of AB, and if UPLO = 'L',   
            the diagonal and first subdiagonal of T are returned in the   
            first two rows of AB.   

    LDAB    (input) INTEGER   
            The leading dimension of the array AB.  LDAB >= KD + 1.   

    W       (output) DOUBLE PRECISION array, dimension (N)   
            If INFO = 0, the eigenvalues in ascending order.   

    Z       (output) DOUBLE PRECISION array, dimension (LDZ, N)   
            If JOBZ = 'V', then if INFO = 0, Z contains the orthonormal   
            eigenvectors of the matrix A, with the i-th column of Z   
            holding the eigenvector associated with W(i).   
            If JOBZ = 'N', then Z is not referenced.   

    LDZ     (input) INTEGER   
            The leading dimension of the array Z.  LDZ >= 1, and if   
            JOBZ = 'V', LDZ >= max(1,N).   

    WORK    (workspace/output) DOUBLE PRECISION array,   
                                           dimension (LWORK)   
            On exit, if LWORK > 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INTEGER   
            The dimension of the array WORK.   
            IF N <= 1,                LWORK must be at least 1.   
            If JOBZ  = 'N' and N > 2, LWORK must be at least 2*N.   
            If JOBZ  = 'V' and N > 2, LWORK must be at least   
                           ( 1 + 4*N + 2*N*lg N + 3*N**2 ),   
                           where lg( N ) = smallest integer k such   
                                           that 2**k >= N.   

    IWORK   (workspace/output) INTEGER array, dimension (LIWORK)   
            On exit, if LIWORK > 0, IWORK(1) returns the optimal LIWORK. 
  

    LIWORK  (input) INTEGER   
            The dimension of the array LIWORK.   
            If JOBZ  = 'N' or N <= 1, LIWORK must be at least 1.   
            If JOBZ  = 'V' and N > 2, LIWORK must be at least 2 + 5*N.   

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
    static doublereal c_b14 = 1.;
    static doublereal c_b21 = 0.;
    static integer c__1 = 1;
    
    /* System generated locals */
    integer ab_dim1, ab_offset, z_dim1, z_offset, i__1;
    doublereal d__1;
    /* Builtin functions */
    double log(doublereal);
    integer pow_ii(integer *, integer *);
    double sqrt(doublereal);
    /* Local variables */
    static integer inde;
    static doublereal anrm, rmin, rmax;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), dgemm_(char *, char *, integer *, integer *, integer *
	    , doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *);
    static doublereal sigma;
    extern logical lsame_(char *, char *);
    static integer iinfo, lwmin;
    static logical lower, wantz;
    static integer indwk2, llwrk2;
    extern doublereal dlamch_(char *);
    static integer iscale;
    extern /* Subroutine */ int dlascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *);
    extern doublereal dlansb_(char *, char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *);
    extern /* Subroutine */ int dstedc_(char *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, integer *, integer *), dlacpy_(char *, integer 
	    *, integer *, doublereal *, integer *, doublereal *, integer *);
    static doublereal safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *);
    static doublereal bignum;
    extern /* Subroutine */ int dsbtrd_(char *, char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *), dsterf_(
	    integer *, doublereal *, doublereal *, integer *);
    static integer indwrk, liwmin;
    static doublereal smlnum;
    static integer lgn;
    static doublereal eps;



#define W(I) w[(I)-1]
#define WORK(I) work[(I)-1]
#define IWORK(I) iwork[(I)-1]

#define AB(I,J) ab[(I)-1 + ((J)-1)* ( *ldab)]
#define Z(I,J) z[(I)-1 + ((J)-1)* ( *ldz)]

    wantz = lsame_(jobz, "V");
    lower = lsame_(uplo, "L");

    *info = 0;
    if (*n <= 1) {
	lgn = 0;
	liwmin = 1;
	lwmin = 1;
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
	    lwmin = (*n << 2) + 1 + (*n << 1) * lgn + i__1 * i__1 * 3;
	} else {
	    liwmin = 1;
	    lwmin = *n << 1;
	}
    }
    if (! (wantz || lsame_(jobz, "N"))) {
	*info = -1;
    } else if (! (lower || lsame_(uplo, "U"))) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*kd < 0) {
	*info = -4;
    } else if (*ldab < *kd + 1) {
	*info = -6;
    } else if (*ldz < 1 || wantz && *ldz < *n) {
	*info = -9;
    } else if (*lwork < lwmin) {
	*info = -11;
    } else if (*liwork < liwmin) {
	*info = -13;
    }

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DSBEVD ", &i__1);
	goto L10;
    }

/*     Quick return if possible */

    if (*n == 0) {
	goto L10;
    }

    if (*n == 1) {
	W(1) = AB(1,1);
	if (wantz) {
	    Z(1,1) = 1.;
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

    anrm = dlansb_("M", uplo, n, kd, &AB(1,1), ldab, &WORK(1));
    iscale = 0;
    if (anrm > 0. && anrm < rmin) {
	iscale = 1;
	sigma = rmin / anrm;
    } else if (anrm > rmax) {
	iscale = 1;
	sigma = rmax / anrm;
    }
    if (iscale == 1) {
	if (lower) {
	    dlascl_("B", kd, kd, &c_b14, &sigma, n, n, &AB(1,1), ldab, 
		    info);
	} else {
	    dlascl_("Q", kd, kd, &c_b14, &sigma, n, n, &AB(1,1), ldab, 
		    info);
	}
    }

/*     Call DSBTRD to reduce symmetric band matrix to tridiagonal form. */

    inde = 1;
    indwrk = inde + *n;
    indwk2 = indwrk + *n * *n;
    llwrk2 = *lwork - indwk2 + 1;
    dsbtrd_(jobz, uplo, n, kd, &AB(1,1), ldab, &W(1), &WORK(inde), &Z(1,1), ldz, &WORK(indwrk), &iinfo);

/*     For eigenvalues only, call DSTERF.  For eigenvectors, call SSTEDC. 
*/

    if (! wantz) {
	dsterf_(n, &W(1), &WORK(inde), info);
    } else {
	dstedc_("I", n, &W(1), &WORK(inde), &WORK(indwrk), n, &WORK(indwk2), &
		llwrk2, &IWORK(1), liwork, info);
	dgemm_("N", "N", n, n, n, &c_b14, &Z(1,1), ldz, &WORK(indwrk), n,
		 &c_b21, &WORK(indwk2), n);
	dlacpy_("A", n, n, &WORK(indwk2), n, &Z(1,1), ldz);
    }

/*     If matrix was scaled, then rescale eigenvalues appropriately. */

    if (iscale == 1) {
	d__1 = 1. / sigma;
	dscal_(n, &d__1, &W(1), &c__1);
    }

L10:
    if (*lwork > 0) {
	WORK(1) = (doublereal) lwmin;
    }
    if (*liwork > 0) {
	IWORK(1) = liwmin;
    }
    return 0;

/*     End of DSBEVD */

} /* dsbevd_ */

