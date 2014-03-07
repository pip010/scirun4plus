#include "f2c.h"

/* Subroutine */ int zlaed0_(integer *qsiz, integer *n, doublereal *d, 
	doublereal *e, doublecomplex *q, integer *ldq, doublecomplex *qstore, 
	integer *ldqs, doublereal *rwork, integer *iwork, integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    Using the divide and conquer method, ZLAED0 computes all eigenvalues 
  
    of a symmetric tridiagonal matrix which is one diagonal block of   
    those from reducing a dense or band Hermitian matrix and   
    corresponding eigenvectors of the dense or band matrix.   

    Arguments   
    =========   

    QSIZ   (input) INTEGER   
           The dimension of the unitary matrix used to reduce   
           the full matrix to tridiagonal form.  QSIZ >= N if ICOMPQ = 1. 
  

    N      (input) INTEGER   
           The dimension of the symmetric tridiagonal matrix.  N >= 0.   

    D      (input/output) DOUBLE PRECISION array, dimension (N)   
           On entry, the diagonal elements of the tridiagonal matrix.   
           On exit, the eigenvalues in ascending order.   

    E      (input/output) DOUBLE PRECISION array, dimension (N-1)   
           On entry, the off-diagonal elements of the tridiagonal matrix. 
  
           On exit, E has been destroyed.   

    Q      (input/output) COMPLEX*16 array, dimension (LDQ,N)   
           On entry, Q must contain an QSIZ x N matrix whose columns   
           unitarily orthonormal. It is a part of the unitary matrix   
           that reduces the full dense Hermitian matrix to a   
           (reducible) symmetric tridiagonal matrix.   

    LDQ    (input) INTEGER   
           The leading dimension of the array Q.  LDQ >= max(1,N).   

    IWORK  (workspace) INTEGER array,   
           the dimension of IWORK must be at least   
                        6 + 6*N + 5*N*lg N   
                        ( lg( N ) = smallest integer k   
                                    such that 2^k >= N )   

    RWORK  (workspace) DOUBLE PRECISION array,   
                                 dimension (1 + 3*N + 2*N*lg N + 3*N**2) 
  
                          ( lg( N ) = smallest integer k   
                                      such that 2^k >= N )   

    QSTORE (workspace) COMPLEX*16 array, dimension (LDQS, N)   
           Used to store parts of   
           the eigenvector matrix when the updating matrix multiplies   
           take place.   

    LDQS   (input) INTEGER   
           The leading dimension of the array QSTORE.   
           LDQS >= max(1,N).   

    INFO   (output) INTEGER   
            = 0:  successful exit.   
            < 0:  if INFO = -i, the i-th argument had an illegal value.   
            > 0:  The algorithm failed to compute an eigenvalue while   
                  working on the submatrix lying in rows and columns   
                  INFO/(N+1) through mod(INFO,N+1).   

    ===================================================================== 
  

    Warning:      N could be as big as QSIZ!   


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static integer c__2 = 2;
    static integer c__1 = 1;
    
    /* System generated locals */
    integer q_dim1, q_offset, qstore_dim1, qstore_offset, i__1, i__2;
    doublereal d__1;
    /* Builtin functions */
    double log(doublereal);
    integer pow_ii(integer *, integer *);
    /* Local variables */
    static doublereal temp;
    static integer curr, i, j, k, iperm;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static integer indxq, iwrem, iqptr, tlvls;
    extern /* Subroutine */ int zcopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), zlaed7_(integer *, integer *, 
	    integer *, integer *, integer *, integer *, doublereal *, 
	    doublecomplex *, integer *, doublereal *, integer *, doublereal *,
	     integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublecomplex *, doublereal *, integer *, integer *)
	    ;
    static integer ll, iq, igivcl;
    extern /* Subroutine */ int xerbla_(char *, integer *), zlacrm_(
	    integer *, integer *, doublecomplex *, integer *, doublereal *, 
	    integer *, doublecomplex *, integer *, doublereal *);
    static integer igivnm, submat, curprb, subpbs, igivpt;
    extern /* Subroutine */ int dsteqr_(char *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *);
    static integer curlvl, matsiz, iprmpt, lgn, msd2, smm1, spm1, spm2;



#define D(I) d[(I)-1]
#define E(I) e[(I)-1]
#define RWORK(I) rwork[(I)-1]
#define IWORK(I) iwork[(I)-1]

#define Q(I,J) q[(I)-1 + ((J)-1)* ( *ldq)]
#define QSTORE(I,J) qstore[(I)-1 + ((J)-1)* ( *ldqs)]

    *info = 0;

/*     IF( ICOMPQ .LT. 0 .OR. ICOMPQ .GT. 2 ) THEN   
          INFO = -1   
       ELSE IF( ( ICOMPQ .EQ. 1 ) .AND. ( QSIZ .LT. MAX( 0, N ) ) )   
      $        THEN */
    if (*qsiz < max(0,*n)) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*ldq < max(1,*n)) {
	*info = -6;
    } else if (*ldqs < max(1,*n)) {
	*info = -8;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("ZLAED0", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return 0;
    }

/*     Determine the size and placement of the submatrices, and save in   
       the leading elements of IWORK. */

    IWORK(1) = *n;
    subpbs = 1;
    tlvls = 0;
L10:
    if (IWORK(subpbs) > 25) {
	for (j = subpbs; j >= 1; --j) {
	    IWORK(j * 2) = (IWORK(j) + 1) / 2;
	    IWORK((j << 1) - 1) = IWORK(j) / 2;
/* L20: */
	}
	++tlvls;
	subpbs <<= 1;
	goto L10;
    }
    i__1 = subpbs;
    for (j = 2; j <= subpbs; ++j) {
	IWORK(j) += IWORK(j - 1);
/* L30: */
    }

/*     Divide the matrix into SUBPBS submatrices of size at most SMLSIZ+1 
  
       using rank-1 modifications (cuts). */

    spm1 = subpbs - 1;
    i__1 = spm1;
    for (i = 1; i <= spm1; ++i) {
	submat = IWORK(i) + 1;
	smm1 = submat - 1;
	D(smm1) -= (d__1 = E(smm1), abs(d__1));
	D(submat) -= (d__1 = E(smm1), abs(d__1));
/* L40: */
    }

    indxq = (*n << 2) + 3;

/*     Set up workspaces for eigenvalues only/accumulate new vectors   
       routine */

    temp = log((doublereal) (*n)) / log(2.);
    lgn = (integer) temp;
    if (pow_ii(&c__2, &lgn) < *n) {
	++lgn;
    }
    if (pow_ii(&c__2, &lgn) < *n) {
	++lgn;
    }
    iprmpt = indxq + *n + 1;
    iperm = iprmpt + *n * lgn;
    iqptr = iperm + *n * lgn;
    igivpt = iqptr + *n + 2;
    igivcl = igivpt + *n * lgn;

    igivnm = 1;
    iq = igivnm + (*n << 1) * lgn;
/* Computing 2nd power */
    i__1 = *n;
    iwrem = iq + i__1 * i__1 + 1;
/*     Initialize pointers */
    i__1 = subpbs;
    for (i = 0; i <= subpbs; ++i) {
	IWORK(iprmpt + i) = 1;
	IWORK(igivpt + i) = 1;
/* L50: */
    }
    IWORK(iqptr) = 1;

/*     Solve each submatrix eigenproblem at the bottom of the divide and 
  
       conquer tree. */

    curr = 0;
    i__1 = spm1;
    for (i = 0; i <= spm1; ++i) {
	if (i == 0) {
	    submat = 1;
	    matsiz = IWORK(1);
	} else {
	    submat = IWORK(i) + 1;
	    matsiz = IWORK(i + 1) - IWORK(i);
	}
	ll = iq - 1 + IWORK(iqptr + curr);
	dsteqr_("I", &matsiz, &D(submat), &E(submat), &RWORK(ll), &matsiz, &
		RWORK(1), info);
	zlacrm_(qsiz, &matsiz, &Q(1,submat), ldq, &RWORK(ll), &
		matsiz, &QSTORE(1,submat), ldqs, &RWORK(iwrem)
		);
/* Computing 2nd power */
	i__2 = matsiz;
	IWORK(iqptr + curr + 1) = IWORK(iqptr + curr) + i__2 * i__2;
	++curr;
	if (*info > 0) {
	    *info = submat * (*n + 1) + submat + matsiz - 1;
	    return 0;
	}
	k = 1;
	i__2 = IWORK(i + 1);
	for (j = submat; j <= IWORK(i+1); ++j) {
	    IWORK(indxq + j) = k;
	    ++k;
/* L60: */
	}
/* L70: */
    }

/*     Successively merge eigensystems of adjacent submatrices   
       into eigensystem for the corresponding larger matrix.   

       while ( SUBPBS > 1 ) */

    curlvl = 1;
L80:
    if (subpbs > 1) {
	spm2 = subpbs - 2;
	i__1 = spm2;
	for (i = 0; i <= spm2; i += 2) {
	    if (i == 0) {
		submat = 1;
		matsiz = IWORK(2);
		msd2 = IWORK(1);
		curprb = 0;
	    } else {
		submat = IWORK(i) + 1;
		matsiz = IWORK(i + 2) - IWORK(i);
		msd2 = matsiz / 2;
		++curprb;
	    }

/*     Merge lower order eigensystems (of size MSD2 and MATSIZ - M
SD2)   
       into an eigensystem of size MATSIZ.  ZLAED7 handles the cas
e   
       when the eigenvectors of a full or band Hermitian matrix (w
hich   
       was reduced to tridiagonal form) are desired.   

       I am free to use Q as a valuable working space until Loop 1
50. */

	    zlaed7_(&matsiz, &msd2, qsiz, &tlvls, &curlvl, &curprb, &D(submat)
		    , &QSTORE(1,submat), ldqs, &E(submat + 
		    msd2 - 1), &IWORK(indxq + submat), &RWORK(iq), &IWORK(
		    iqptr), &IWORK(iprmpt), &IWORK(iperm), &IWORK(igivpt), &
		    IWORK(igivcl), &RWORK(igivnm), &Q(1,submat), &
		    RWORK(iwrem), &IWORK(subpbs + 1), info);
	    if (*info > 0) {
		*info = submat * (*n + 1) + submat + matsiz - 1;
		return 0;
	    }
	    IWORK(i / 2 + 1) = IWORK(i + 2);
/* L90: */
	}
	subpbs /= 2;
	++curlvl;
	goto L80;
    }

/*     end while   

       Re-merge the eigenvalues/vectors which were deflated at the final 
  
       merge step. */

    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
	j = IWORK(indxq + i);
	RWORK(i) = D(j);
	zcopy_(qsiz, &QSTORE(1,j), &c__1, &Q(1,i), 
		&c__1);
/* L100: */
    }
    dcopy_(n, &RWORK(1), &c__1, &D(1), &c__1);

    return 0;

/*     End of ZLAED0 */

} /* zlaed0_ */

