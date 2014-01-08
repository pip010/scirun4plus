#include "f2c.h"

/* Subroutine */ int cstedc_(char *compz, integer *n, real *d, real *e, 
	complex *z, integer *ldz, complex *work, integer *lwork, real *rwork, 
	integer *lrwork, integer *iwork, integer *liwork, integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    CSTEDC computes all eigenvalues and, optionally, eigenvectors of a   
    symmetric tridiagonal matrix using the divide and conquer method.   
    The eigenvectors of a full or band complex Hermitian matrix can also 
  
    be found if CHETRD or CHPTRD or CHBTRD has been used to reduce this   
    matrix to tridiagonal form.   

    This code makes very mild assumptions about floating point   
    arithmetic. It will work on machines with a guard digit in   
    add/subtract, or on those binary machines without guard digits   
    which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or Cray-2.   
    It could conceivably fail on hexadecimal or decimal machines   
    without guard digits, but we know of none.  See SLAED3 for details.   

    Arguments   
    =========   

    COMPZ   (input) CHARACTER*1   
            = 'N':  Compute eigenvalues only.   
            = 'I':  Compute eigenvectors of tridiagonal matrix also.   
            = 'V':  Compute eigenvectors of original Hermitian matrix   
                    also.  On entry, Z contains the unitary matrix used   
                    to reduce the original matrix to tridiagonal form.   

    N       (input) INTEGER   
            The dimension of the symmetric tridiagonal matrix.  N >= 0.   

    D       (input/output) REAL array, dimension (N)   
            On entry, the diagonal elements of the tridiagonal matrix.   
            On exit, if INFO = 0, the eigenvalues in ascending order.   

    E       (input/output) REAL array, dimension (N-1)   
            On entry, the subdiagonal elements of the tridiagonal matrix. 
  
            On exit, E has been destroyed.   

    Z       (input/output) COMPLEX array, dimension (LDZ,N)   
            On entry, if COMPZ = 'V', then Z contains the unitary   
            matrix used in the reduction to tridiagonal form.   
            On exit, if INFO = 0, then if COMPZ = 'V', Z contains the   
            orthonormal eigenvectors of the original Hermitian matrix,   
            and if COMPZ = 'I', Z contains the orthonormal eigenvectors   
            of the symmetric tridiagonal matrix.   
            If  COMPZ = 'N', then Z is not referenced.   

    LDZ     (input) INTEGER   
            The leading dimension of the array Z.  LDZ >= 1.   
            If eigenvectors are desired, then LDZ >= max(1,N).   

    WORK    (workspace/output) COMPLEX array, dimension (LWORK)   
            On exit, if LWORK > 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INTEGER   
            The dimension of the array WORK.   
            If COMPZ = 'N' or 'I', or N <= 1, LWORK must be at least 1.   
            If COMPZ = 'V' and N > 1, LWORK must be at least N*N.   

    RWORK   (workspace/output) REAL array,   
                                           dimension (LRWORK)   
            On exit, if LRWORK > 0, RWORK(1) returns the optimal LRWORK. 
  

    LRWORK  (input) INTEGER   
            The dimension of the array RWORK.   
            If COMPZ = 'N' or N <= 1, LRWORK must be at least 1.   
            If COMPZ = 'V' and N > 1, LRWORK must be at least   
                           1 + 3*N + 2*N*lg N + 3*N**2 ,   
                           where lg( N ) = smallest integer k such   
                           that 2**k >= N.   
            If COMPZ = 'I' and N > 1, LRWORK must be at least   
                           1 + 3*N + 2*N*lg N + 3*N**2 .   

    IWORK   (workspace/output) INTEGER array, dimension (LIWORK)   
            On exit, if LIWORK > 0, IWORK(1) returns the optimal LIWORK. 
  

    LIWORK  (input) INTEGER   
            The dimension of the array IWORK.   
            If COMPZ = 'N' or N <= 1, LIWORK must be at least 1.   
            If COMPZ = 'V' or N > 1,  LIWORK must be at least   
                                      6 + 6*N + 5*N*lg N.   
            If COMPZ = 'I' or N > 1,  LIWORK must be at least   
                                      2 + 5*N .   

    INFO    (output) INTEGER   
            = 0:  successful exit.   
            < 0:  if INFO = -i, the i-th argument had an illegal value.   
            > 0:  The algorithm failed to compute an eigenvalue while   
                  working on the submatrix lying in rows and columns   
                  INFO/(N+1) through mod(INFO,N+1).   

    ===================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static integer c__2 = 2;
    static real c_b12 = 0.f;
    static real c_b13 = 1.f;
    static integer c__0 = 0;
    static integer c__1 = 1;
    
    /* System generated locals */
    integer z_dim1, z_offset, i__1, i__2, i__3, i__4;
    real r__1, r__2;
    /* Builtin functions */
    double log(doublereal);
    integer pow_ii(integer *, integer *);
    double sqrt(doublereal);
    /* Local variables */
    static real tiny;
    static integer i, j, k, m;
    static real p;
    extern logical lsame_(char *, char *);
    extern /* Subroutine */ int cswap_(integer *, complex *, integer *, 
	    complex *, integer *);
    static integer lwmin;
    extern /* Subroutine */ int claed0_(integer *, integer *, real *, real *, 
	    complex *, integer *, complex *, integer *, real *, integer *, 
	    integer *);
    static integer start, ii, ll;
    extern /* Subroutine */ int clacrm_(integer *, integer *, complex *, 
	    integer *, real *, integer *, complex *, integer *, real *);
    extern doublereal slamch_(char *);
    extern /* Subroutine */ int clacpy_(char *, integer *, integer *, complex 
	    *, integer *, complex *, integer *), xerbla_(char *, 
	    integer *), slascl_(char *, integer *, integer *, real *, 
	    real *, integer *, integer *, real *, integer *, integer *), sstedc_(char *, integer *, real *, real *, real *, 
	    integer *, real *, integer *, integer *, integer *, integer *), slaset_(char *, integer *, integer *, real *, real *, 
	    real *, integer *);
    static integer liwmin, icompz;
    extern /* Subroutine */ int csteqr_(char *, integer *, real *, real *, 
	    complex *, integer *, real *, integer *);
    static real orgnrm;
    extern doublereal slanst_(char *, integer *, real *, real *);
    extern /* Subroutine */ int ssterf_(integer *, real *, real *, integer *);
    static integer lrwmin, nm1;
    extern /* Subroutine */ int ssteqr_(char *, integer *, real *, real *, 
	    real *, integer *, real *, integer *);
    static integer end, lgn;
    static real eps;



#define D(I) d[(I)-1]
#define E(I) e[(I)-1]
#define WORK(I) work[(I)-1]
#define RWORK(I) rwork[(I)-1]
#define IWORK(I) iwork[(I)-1]

#define Z(I,J) z[(I)-1 + ((J)-1)* ( *ldz)]

    *info = 0;

    if (lsame_(compz, "N")) {
	icompz = 0;
    } else if (lsame_(compz, "V")) {
	icompz = 1;
    } else if (lsame_(compz, "I")) {
	icompz = 2;
    } else {
	icompz = -1;
    }
    if (*n <= 1 || icompz <= 0) {
	lwmin = 1;
	liwmin = 1;
	lrwmin = 1;
    } else {
	lgn = (integer) (log((real) (*n)) / log(2.f));
	if (pow_ii(&c__2, &lgn) < *n) {
	    ++lgn;
	}
	if (pow_ii(&c__2, &lgn) < *n) {
	    ++lgn;
	}
	if (icompz == 1) {
	    lwmin = *n * *n;
/* Computing 2nd power */
	    i__1 = *n;
	    lrwmin = *n * 3 + 1 + (*n << 1) * lgn + i__1 * i__1 * 3;
	    liwmin = *n * 6 + 6 + *n * 5 * lgn;
	} else if (icompz == 2) {
	    lwmin = 1;
/* Computing 2nd power */
	    i__1 = *n;
	    lrwmin = *n * 3 + 1 + (*n << 1) * lgn + i__1 * i__1 * 3;
	    liwmin = *n * 5 + 2;
	}
    }
    if (icompz < 0) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*ldz < 1 || icompz > 0 && *ldz < max(1,*n)) {
	*info = -6;
    } else if (*lwork < lwmin) {
	*info = -8;
    } else if (*lrwork < lrwmin) {
	*info = -10;
    } else if (*liwork < liwmin) {
	*info = -12;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("CSTEDC", &i__1);
	goto L70;
    }

/*     Quick return if possible */

    if (*n == 0) {
	goto L70;
    }
    if (*n == 1) {
	if (icompz != 0) {
	    i__1 = z_dim1 + 1;
	    Z(1,1).r = 1.f, Z(1,1).i = 0.f;
	}
	goto L70;
    }

/*     If the following conditional clause is removed, then the routine   
       will use the Divide and Conquer routine to compute only the   
       eigenvalues, which requires (3N + 3N**2) real workspace and   
       (2 + 5N + 2N lg(N)) integer workspace.   
       Since on many architectures SSTERF is much faster than any other   
       algorithm for finding eigenvalues only, it is used here   
       as the default.   

       If COMPZ = 'N', use SSTERF to compute the eigenvalues. */

    if (icompz == 0) {
	ssterf_(n, &D(1), &E(1), info);
	goto L70;
    }

/*     If N is smaller than the minimum divide size (SMLSIZ+1), then   
       solve the problem with another solver. */

    if (*n <= 25) {
	if (icompz == 0) {
	    ssterf_(n, &D(1), &E(1), info);
	    goto L70;
	} else if (icompz == 2) {
	    csteqr_("I", n, &D(1), &E(1), &Z(1,1), ldz, &RWORK(1), info);
	    goto L70;
	} else {
	    csteqr_("V", n, &D(1), &E(1), &Z(1,1), ldz, &RWORK(1), info);
	    goto L70;
	}
    }

/*     If COMPZ = 'I', we simply call SSTEDC instead. */

    if (icompz == 2) {
	slaset_("Full", n, n, &c_b12, &c_b13, &RWORK(1), n);
	ll = *n * *n + 1;
	i__1 = *lrwork - ll + 1;
	sstedc_("I", n, &D(1), &E(1), &RWORK(1), n, &RWORK(ll), &i__1, &IWORK(
		1), liwork, info);
	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
	    i__2 = *n;
	    for (i = 1; i <= *n; ++i) {
		i__3 = i + j * z_dim1;
		i__4 = (j - 1) * *n + i;
		Z(i,j).r = RWORK((j-1)**n+i), Z(i,j).i = 0.f;
/* L10: */
	    }
/* L20: */
	}
	goto L70;
    }

/*     From now on, only option left to be handled is COMPZ = 'V',   
       i.e. ICOMPZ = 1.   

       Scale. */

    nm1 = *n - 1;
    orgnrm = slanst_("M", n, &D(1), &E(1));
    if (orgnrm == 0.f) {
	goto L70;
    }

    eps = slamch_("Epsilon");

    start = 1;

/*     while ( START <= N ) */

L30:
    if (start <= *n) {

/*     Let END be the position of the next subdiagonal entry such that
   
       E( END ) <= TINY or END = N if no such subdiagonal exists.  The
   
       matrix identified by the elements between START and END   
       constitutes an independent sub-problem. */

	end = start;
L40:
	if (end < *n) {
	    tiny = eps * sqrt((r__1 = D(end), dabs(r__1))) * sqrt((r__2 = D(
		    end + 1), dabs(r__2)));
	    if ((r__1 = E(end), dabs(r__1)) > tiny) {
		++end;
		goto L40;
	    }
	}

/*        (Sub) Problem determined.  Compute its size and solve it. */

	m = end - start + 1;
	if (m > 25) {
	    *info = 25;

/*           Scale. */

	    orgnrm = slanst_("M", &m, &D(start), &E(start));
	    slascl_("G", &c__0, &c__0, &orgnrm, &c_b13, &m, &c__1, &D(start), 
		    &m, info);
	    i__1 = m - 1;
	    i__2 = m - 1;
	    slascl_("G", &c__0, &c__0, &orgnrm, &c_b13, &i__1, &c__1, &E(
		    start), &i__2, info);

	    claed0_(n, &m, &D(start), &E(start), &Z(1,start), ldz, 
		    &WORK(1), n, &RWORK(1), &IWORK(1), info);
	    if (*info > 0) {
		*info = (*info / (m + 1) + start - 1) * (*n + 1) + *info % (m 
			+ 1) + start - 1;
		goto L70;
	    }

/*           Scale back. */

	    slascl_("G", &c__0, &c__0, &c_b13, &orgnrm, &m, &c__1, &D(start), 
		    &m, info);

	} else {
	    ssteqr_("I", &m, &D(start), &E(start), &RWORK(1), &m, &RWORK(m * 
		    m + 1), info);
	    clacrm_(n, &m, &Z(1,start), ldz, &RWORK(1), &m, &WORK(
		    1), n, &RWORK(m * m + 1));
	    clacpy_("A", n, &m, &WORK(1), n, &Z(1,start), ldz);
	    if (*info > 0) {
		*info = start * (*n + 1) + end;
		goto L70;
	    }
	}

	start = end + 1;
	goto L30;
    }

/*     endwhile   

       If the problem split any number of times, then the eigenvalues   
       will not be properly ordered.  Here we permute the eigenvalues   
       (and the associated eigenvectors) into ascending order. */

    if (m != *n) {

/*        Use Selection Sort to minimize swaps of eigenvectors */

	i__1 = *n;
	for (ii = 2; ii <= *n; ++ii) {
	    i = ii - 1;
	    k = i;
	    p = D(i);
	    i__2 = *n;
	    for (j = ii; j <= *n; ++j) {
		if (D(j) < p) {
		    k = j;
		    p = D(j);
		}
/* L50: */
	    }
	    if (k != i) {
		D(k) = D(i);
		D(i) = p;
		cswap_(n, &Z(1,i), &c__1, &Z(1,k), &
			c__1);
	    }
/* L60: */
	}
    }

L70:
    if (*lwork > 0) {
	WORK(1).r = (real) lwmin, WORK(1).i = 0.f;
    }
    if (*lrwork > 0) {
	RWORK(1) = (real) lrwmin;
    }
    if (*liwork > 0) {
	IWORK(1) = liwmin;
    }
    return 0;

/*     End of CSTEDC */

} /* cstedc_ */

