#include "f2c.h"

/* Subroutine */ int zstedc_(char *compz, integer *n, doublereal *d, 
	doublereal *e, doublecomplex *z, integer *ldz, doublecomplex *work, 
	integer *lwork, doublereal *rwork, integer *lrwork, integer *iwork, 
	integer *liwork, integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    ZSTEDC computes all eigenvalues and, optionally, eigenvectors of a   
    symmetric tridiagonal matrix using the divide and conquer method.   
    The eigenvectors of a full or band complex Hermitian matrix can also 
  
    be found if ZHETRD or ZHPTRD or ZHBTRD has been used to reduce this   
    matrix to tridiagonal form.   

    This code makes very mild assumptions about floating point   
    arithmetic. It will work on machines with a guard digit in   
    add/subtract, or on those binary machines without guard digits   
    which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or Cray-2.   
    It could conceivably fail on hexadecimal or decimal machines   
    without guard digits, but we know of none.  See DLAED3 for details.   

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

    D       (input/output) DOUBLE PRECISION array, dimension (N)   
            On entry, the diagonal elements of the tridiagonal matrix.   
            On exit, if INFO = 0, the eigenvalues in ascending order.   

    E       (input/output) DOUBLE PRECISION array, dimension (N-1)   
            On entry, the subdiagonal elements of the tridiagonal matrix. 
  
            On exit, E has been destroyed.   

    Z       (input/output) COMPLEX*16 array, dimension (LDZ,N)   
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

    WORK    (workspace/output) COMPLEX*16 array, dimension (LWORK)   
            On exit, if LWORK > 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INTEGER   
            The dimension of the array WORK.   
            If COMPZ = 'N' or 'I', or N <= 1, LWORK must be at least 1.   
            If COMPZ = 'V' and N > 1, LWORK must be at least N*N.   

    RWORK   (workspace/output) DOUBLE PRECISION array,   
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
    static doublereal c_b12 = 0.;
    static doublereal c_b13 = 1.;
    static integer c__0 = 0;
    static integer c__1 = 1;
    
    /* System generated locals */
    integer z_dim1, z_offset, i__1, i__2, i__3, i__4;
    doublereal d__1, d__2;
    /* Builtin functions */
    double log(doublereal);
    integer pow_ii(integer *, integer *);
    double sqrt(doublereal);
    /* Local variables */
    static doublereal tiny;
    static integer i, j, k, m;
    static doublereal p;
    extern logical lsame_(char *, char *);
    static integer lwmin, start;
    extern /* Subroutine */ int zswap_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), zlaed0_(integer *, integer *, 
	    doublereal *, doublereal *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublereal *, integer *, integer *);
    static integer ii, ll;
    extern doublereal dlamch_(char *);
    extern /* Subroutine */ int dlascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *), dstedc_(char *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     integer *, integer *, integer *, integer *), dlaset_(
	    char *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *), xerbla_(char *, integer *);
    extern doublereal dlanst_(char *, integer *, doublereal *, doublereal *);
    extern /* Subroutine */ int dsterf_(integer *, doublereal *, doublereal *,
	     integer *), zlacrm_(integer *, integer *, doublecomplex *, 
	    integer *, doublereal *, integer *, doublecomplex *, integer *, 
	    doublereal *);
    static integer liwmin, icompz;
    extern /* Subroutine */ int dsteqr_(char *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *), zlacpy_(char *, integer *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, integer *);
    static doublereal orgnrm;
    static integer lrwmin, nm1;
    extern /* Subroutine */ int zsteqr_(char *, integer *, doublereal *, 
	    doublereal *, doublecomplex *, integer *, doublereal *, integer *);
    static integer end, lgn;
    static doublereal eps;



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
	lgn = (integer) (log((doublereal) (*n)) / log(2.));
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
	xerbla_("ZSTEDC", &i__1);
	goto L70;
    }

/*     Quick return if possible */

    if (*n == 0) {
	goto L70;
    }
    if (*n == 1) {
	if (icompz != 0) {
	    i__1 = z_dim1 + 1;
	    Z(1,1).r = 1., Z(1,1).i = 0.;
	}
	goto L70;
    }

/*     If the following conditional clause is removed, then the routine   
       will use the Divide and Conquer routine to compute only the   
       eigenvalues, which requires (3N + 3N**2) real workspace and   
       (2 + 5N + 2N lg(N)) integer workspace.   
       Since on many architectures DSTERF is much faster than any other   
       algorithm for finding eigenvalues only, it is used here   
       as the default.   

       If COMPZ = 'N', use DSTERF to compute the eigenvalues. */

    if (icompz == 0) {
	dsterf_(n, &D(1), &E(1), info);
	goto L70;
    }

/*     If N is smaller than the minimum divide size (SMLSIZ+1), then   
       solve the problem with another solver. */

    if (*n <= 25) {
	if (icompz == 0) {
	    dsterf_(n, &D(1), &E(1), info);
	    goto L70;
	} else if (icompz == 2) {
	    zsteqr_("I", n, &D(1), &E(1), &Z(1,1), ldz, &RWORK(1), info);
	    goto L70;
	} else {
	    zsteqr_("V", n, &D(1), &E(1), &Z(1,1), ldz, &RWORK(1), info);
	    goto L70;
	}
    }

/*     If COMPZ = 'I', we simply call DSTEDC instead. */

    if (icompz == 2) {
	dlaset_("Full", n, n, &c_b12, &c_b13, &RWORK(1), n);
	ll = *n * *n + 1;
	i__1 = *lrwork - ll + 1;
	dstedc_("I", n, &D(1), &E(1), &RWORK(1), n, &RWORK(ll), &i__1, &IWORK(
		1), liwork, info);
	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
	    i__2 = *n;
	    for (i = 1; i <= *n; ++i) {
		i__3 = i + j * z_dim1;
		i__4 = (j - 1) * *n + i;
		Z(i,j).r = RWORK((j-1)**n+i), Z(i,j).i = 0.;
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
    orgnrm = dlanst_("M", n, &D(1), &E(1));
    if (orgnrm == 0.) {
	goto L70;
    }

    eps = dlamch_("Epsilon");

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
	    tiny = eps * sqrt((d__1 = D(end), abs(d__1))) * sqrt((d__2 = D(
		    end + 1), abs(d__2)));
	    if ((d__1 = E(end), abs(d__1)) > tiny) {
		++end;
		goto L40;
	    }
	}

/*        (Sub) Problem determined.  Compute its size and solve it. */

	m = end - start + 1;
	if (m > 25) {
	    *info = 25;

/*           Scale. */

	    orgnrm = dlanst_("M", &m, &D(start), &E(start));
	    dlascl_("G", &c__0, &c__0, &orgnrm, &c_b13, &m, &c__1, &D(start), 
		    &m, info);
	    i__1 = m - 1;
	    i__2 = m - 1;
	    dlascl_("G", &c__0, &c__0, &orgnrm, &c_b13, &i__1, &c__1, &E(
		    start), &i__2, info);

	    zlaed0_(n, &m, &D(start), &E(start), &Z(1,start), ldz, 
		    &WORK(1), n, &RWORK(1), &IWORK(1), info);
	    if (*info > 0) {
		*info = (*info / (m + 1) + start - 1) * (*n + 1) + *info % (m 
			+ 1) + start - 1;
		goto L70;
	    }

/*           Scale back. */

	    dlascl_("G", &c__0, &c__0, &c_b13, &orgnrm, &m, &c__1, &D(start), 
		    &m, info);

	} else {
	    dsteqr_("I", &m, &D(start), &E(start), &RWORK(1), &m, &RWORK(m * 
		    m + 1), info);
	    zlacrm_(n, &m, &Z(1,start), ldz, &RWORK(1), &m, &WORK(
		    1), n, &RWORK(m * m + 1));
	    zlacpy_("A", n, &m, &WORK(1), n, &Z(1,start), ldz);
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
		zswap_(n, &Z(1,i), &c__1, &Z(1,k), &
			c__1);
	    }
/* L60: */
	}
    }

L70:
    if (*lwork > 0) {
	WORK(1).r = (doublereal) lwmin, WORK(1).i = 0.;
    }
    if (*lrwork > 0) {
	RWORK(1) = (doublereal) lrwmin;
    }
    if (*liwork > 0) {
	IWORK(1) = liwmin;
    }
    return 0;

/*     End of ZSTEDC */

} /* zstedc_ */

