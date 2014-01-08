#include "f2c.h"

/* Subroutine */ int sgelsx_(integer *m, integer *n, integer *nrhs, real *a, 
	integer *lda, real *b, integer *ldb, integer *jpvt, real *rcond, 
	integer *rank, real *work, integer *info)
{
/*  -- LAPACK driver routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       March 31, 1993   


    Purpose   
    =======   

    SGELSX computes the minimum-norm solution to a real linear least   
    squares problem:   
        minimize || A * X - B ||   
    using a complete orthogonal factorization of A.  A is an M-by-N   
    matrix which may be rank-deficient.   

    Several right hand side vectors b and solution vectors x can be   
    handled in a single call; they are stored as the columns of the   
    M-by-NRHS right hand side matrix B and the N-by-NRHS solution   
    matrix X.   

    The routine first computes a QR factorization with column pivoting:   
        A * P = Q * [ R11 R12 ]   
                    [  0  R22 ]   
    with R11 defined as the largest leading submatrix whose estimated   
    condition number is less than 1/RCOND.  The order of R11, RANK,   
    is the effective rank of A.   

    Then, R22 is considered to be negligible, and R12 is annihilated   
    by orthogonal transformations from the right, arriving at the   
    complete orthogonal factorization:   
       A * P = Q * [ T11 0 ] * Z   
                   [  0  0 ]   
    The minimum-norm solution is then   
       X = P * Z' [ inv(T11)*Q1'*B ]   
                  [        0       ]   
    where Q1 consists of the first RANK columns of Q.   

    Arguments   
    =========   

    M       (input) INTEGER   
            The number of rows of the matrix A.  M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrix A.  N >= 0.   

    NRHS    (input) INTEGER   
            The number of right hand sides, i.e., the number of   
            columns of matrices B and X. NRHS >= 0.   

    A       (input/output) REAL array, dimension (LDA,N)   
            On entry, the M-by-N matrix A.   
            On exit, A has been overwritten by details of its   
            complete orthogonal factorization.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(1,M).   

    B       (input/output) REAL array, dimension (LDB,NRHS)   
            On entry, the M-by-NRHS right hand side matrix B.   
            On exit, the N-by-NRHS solution matrix X.   
            If m >= n and RANK = n, the residual sum-of-squares for   
            the solution in the i-th column is given by the sum of   
            squares of elements N+1:M in that column.   

    LDB     (input) INTEGER   
            The leading dimension of the array B. LDB >= max(1,M,N).   

    JPVT    (input/output) INTEGER array, dimension (N)   
            On entry, if JPVT(i) .ne. 0, the i-th column of A is an   
            initial column, otherwise it is a free column.  Before   
            the QR factorization of A, all initial columns are   
            permuted to the leading positions; only the remaining   
            free columns are moved as a result of column pivoting   
            during the factorization.   
            On exit, if JPVT(i) = k, then the i-th column of A*P   
            was the k-th column of A.   

    RCOND   (input) REAL   
            RCOND is used to determine the effective rank of A, which   
            is defined as the order of the largest leading triangular   
            submatrix R11 in the QR factorization with pivoting of A,   
            whose estimated condition number < 1/RCOND.   

    RANK    (output) INTEGER   
            The effective rank of A, i.e., the order of the submatrix   
            R11.  This is the same as the order of the submatrix T11   
            in the complete orthogonal factorization of A.   

    WORK    (workspace) REAL array, dimension   
                        (max( min(M,N)+3*N, 2*min(M,N)+NRHS )),   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   

    ===================================================================== 
  


    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static integer c__0 = 0;
    static real c_b13 = 0.f;
    static integer c__2 = 2;
    static integer c__1 = 1;
    static real c_b36 = 1.f;
    
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2;
    real r__1;
    /* Local variables */
    static real anrm, bnrm, smin, smax;
    static integer i, j, k, iascl, ibscl, ismin, ismax;
    static real c1, c2, s1, s2, t1, t2;
    extern /* Subroutine */ int strsm_(char *, char *, char *, char *, 
	    integer *, integer *, real *, real *, integer *, real *, integer *
	    ), slaic1_(integer *, integer *, 
	    real *, real *, real *, real *, real *, real *, real *), sorm2r_(
	    char *, char *, integer *, integer *, integer *, real *, integer *
	    , real *, real *, integer *, real *, integer *), 
	    slabad_(real *, real *);
    static integer mn;
    extern doublereal slamch_(char *), slange_(char *, integer *, 
	    integer *, real *, integer *, real *);
    extern /* Subroutine */ int xerbla_(char *, integer *);
    static real bignum;
    extern /* Subroutine */ int slascl_(char *, integer *, integer *, real *, 
	    real *, integer *, integer *, real *, integer *, integer *), sgeqpf_(integer *, integer *, real *, integer *, integer 
	    *, real *, real *, integer *), slaset_(char *, integer *, integer 
	    *, real *, real *, real *, integer *);
    static real sminpr, smaxpr, smlnum;
    extern /* Subroutine */ int slatzm_(char *, integer *, integer *, real *, 
	    integer *, real *, real *, real *, integer *, real *), 
	    stzrqf_(integer *, integer *, real *, integer *, real *, integer *
	    );



#define JPVT(I) jpvt[(I)-1]
#define WORK(I) work[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]
#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]

    mn = min(*m,*n);
    ismin = mn + 1;
    ismax = (mn << 1) + 1;

/*     Test the input arguments. */

    *info = 0;
    if (*m < 0) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*nrhs < 0) {
	*info = -3;
    } else if (*lda < max(1,*m)) {
	*info = -5;
    } else /* if(complicated condition) */ {
/* Computing MAX */
	i__1 = max(1,*m);
	if (*ldb < max(i__1,*n)) {
	    *info = -7;
	}
    }

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("SGELSX", &i__1);
	return 0;
    }

/*     Quick return if possible   

   Computing MIN */
    i__1 = min(*m,*n);
    if (min(i__1,*nrhs) == 0) {
	*rank = 0;
	return 0;
    }

/*     Get machine parameters */

    smlnum = slamch_("S") / slamch_("P");
    bignum = 1.f / smlnum;
    slabad_(&smlnum, &bignum);

/*     Scale A, B if max elements outside range [SMLNUM,BIGNUM] */

    anrm = slange_("M", m, n, &A(1,1), lda, &WORK(1));
    iascl = 0;
    if (anrm > 0.f && anrm < smlnum) {

/*        Scale matrix norm up to SMLNUM */

	slascl_("G", &c__0, &c__0, &anrm, &smlnum, m, n, &A(1,1), lda, 
		info);
	iascl = 1;
    } else if (anrm > bignum) {

/*        Scale matrix norm down to BIGNUM */

	slascl_("G", &c__0, &c__0, &anrm, &bignum, m, n, &A(1,1), lda, 
		info);
	iascl = 2;
    } else if (anrm == 0.f) {

/*        VISMatrix all zero. Return zero solution. */

	i__1 = max(*m,*n);
	slaset_("F", &i__1, nrhs, &c_b13, &c_b13, &B(1,1), ldb);
	*rank = 0;
	goto L100;
    }

    bnrm = slange_("M", m, nrhs, &B(1,1), ldb, &WORK(1));
    ibscl = 0;
    if (bnrm > 0.f && bnrm < smlnum) {

/*        Scale matrix norm up to SMLNUM */

	slascl_("G", &c__0, &c__0, &bnrm, &smlnum, m, nrhs, &B(1,1), ldb,
		 info);
	ibscl = 1;
    } else if (bnrm > bignum) {

/*        Scale matrix norm down to BIGNUM */

	slascl_("G", &c__0, &c__0, &bnrm, &bignum, m, nrhs, &B(1,1), ldb,
		 info);
	ibscl = 2;
    }

/*     Compute QR factorization with column pivoting of A:   
          A * P = Q * R */

    sgeqpf_(m, n, &A(1,1), lda, &JPVT(1), &WORK(1), &WORK(mn + 1), info);

/*     workspace 3*N. Details of Householder rotations stored   
       in WORK(1:MN).   

       Determine RANK using incremental condition estimation */

    WORK(ismin) = 1.f;
    WORK(ismax) = 1.f;
    smax = (r__1 = A(1,1), dabs(r__1));
    smin = smax;
    if ((r__1 = A(1,1), dabs(r__1)) == 0.f) {
	*rank = 0;
	i__1 = max(*m,*n);
	slaset_("F", &i__1, nrhs, &c_b13, &c_b13, &B(1,1), ldb);
	goto L100;
    } else {
	*rank = 1;
    }

L10:
    if (*rank < mn) {
	i = *rank + 1;
	slaic1_(&c__2, rank, &WORK(ismin), &smin, &A(1,i), &A(i,i), &sminpr, &s1, &c1);
	slaic1_(&c__1, rank, &WORK(ismax), &smax, &A(1,i), &A(i,i), &smaxpr, &s2, &c2);

	if (smaxpr * *rcond <= sminpr) {
	    i__1 = *rank;
	    for (i = 1; i <= *rank; ++i) {
		WORK(ismin + i - 1) = s1 * WORK(ismin + i - 1);
		WORK(ismax + i - 1) = s2 * WORK(ismax + i - 1);
/* L20: */
	    }
	    WORK(ismin + *rank) = c1;
	    WORK(ismax + *rank) = c2;
	    smin = sminpr;
	    smax = smaxpr;
	    ++(*rank);
	    goto L10;
	}
    }

/*     Logically partition R = [ R11 R12 ]   
                               [  0  R22 ]   
       where R11 = R(1:RANK,1:RANK)   

       [R11,R12] = [ T11, 0 ] * Y */

    if (*rank < *n) {
	stzrqf_(rank, n, &A(1,1), lda, &WORK(mn + 1), info);
    }

/*     Details of Householder rotations stored in WORK(MN+1:2*MN)   

       B(1:M,1:NRHS) := Q' * B(1:M,1:NRHS) */

    sorm2r_("Left", "Transpose", m, nrhs, &mn, &A(1,1), lda, &WORK(1), &
	    B(1,1), ldb, &WORK((mn << 1) + 1), info);

/*     workspace NRHS   

       B(1:RANK,1:NRHS) := inv(T11) * B(1:RANK,1:NRHS) */

    strsm_("Left", "Upper", "No transpose", "Non-unit", rank, nrhs, &c_b36, &
	    A(1,1), lda, &B(1,1), ldb);

    i__1 = *n;
    for (i = *rank + 1; i <= *n; ++i) {
	i__2 = *nrhs;
	for (j = 1; j <= *nrhs; ++j) {
	    B(i,j) = 0.f;
/* L30: */
	}
/* L40: */
    }

/*     B(1:N,1:NRHS) := Y' * B(1:N,1:NRHS) */

    if (*rank < *n) {
	i__1 = *rank;
	for (i = 1; i <= *rank; ++i) {
	    i__2 = *n - *rank + 1;
	    slatzm_("Left", &i__2, nrhs, &A(i,*rank+1), lda, &
		    WORK(mn + i), &B(i,1), &B(*rank+1,1), ldb,
		     &WORK((mn << 1) + 1));
/* L50: */
	}
    }

/*     workspace NRHS   

       B(1:N,1:NRHS) := P * B(1:N,1:NRHS) */

    i__1 = *nrhs;
    for (j = 1; j <= *nrhs; ++j) {
	i__2 = *n;
	for (i = 1; i <= *n; ++i) {
	    WORK((mn << 1) + i) = 1.f;
/* L60: */
	}
	i__2 = *n;
	for (i = 1; i <= *n; ++i) {
	    if (WORK((mn << 1) + i) == 1.f) {
		if (JPVT(i) != i) {
		    k = i;
		    t1 = B(k,j);
		    t2 = B(JPVT(k),j);
L70:
		    B(JPVT(k),j) = t1;
		    WORK((mn << 1) + k) = 0.f;
		    t1 = t2;
		    k = JPVT(k);
		    t2 = B(JPVT(k),j);
		    if (JPVT(k) != i) {
			goto L70;
		    }
		    B(i,j) = t1;
		    WORK((mn << 1) + k) = 0.f;
		}
	    }
/* L80: */
	}
/* L90: */
    }

/*     Undo scaling */

    if (iascl == 1) {
	slascl_("G", &c__0, &c__0, &anrm, &smlnum, n, nrhs, &B(1,1), ldb,
		 info);
	slascl_("U", &c__0, &c__0, &smlnum, &anrm, rank, rank, &A(1,1), 
		lda, info);
    } else if (iascl == 2) {
	slascl_("G", &c__0, &c__0, &anrm, &bignum, n, nrhs, &B(1,1), ldb,
		 info);
	slascl_("U", &c__0, &c__0, &bignum, &anrm, rank, rank, &A(1,1), 
		lda, info);
    }
    if (ibscl == 1) {
	slascl_("G", &c__0, &c__0, &smlnum, &bnrm, n, nrhs, &B(1,1), ldb,
		 info);
    } else if (ibscl == 2) {
	slascl_("G", &c__0, &c__0, &bignum, &bnrm, n, nrhs, &B(1,1), ldb,
		 info);
    }

L100:

    return 0;

/*     End of SGELSX */

} /* sgelsx_ */

