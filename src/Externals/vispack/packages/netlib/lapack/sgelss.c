#include "f2c.h"

/* Subroutine */ int sgelss_(integer *m, integer *n, integer *nrhs, real *a, 
	integer *lda, real *b, integer *ldb, real *s, real *rcond, integer *
	rank, real *work, integer *lwork, integer *info)
{
/*  -- LAPACK driver routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    SGELSS computes the minimum norm solution to a real linear least   
    squares problem:   

    Minimize 2-norm(| b - A*x |).   

    using the singular value decomposition (SVD) of A. A is an M-by-N   
    matrix which may be rank-deficient.   

    Several right hand side vectors b and solution vectors x can be   
    handled in a single call; they are stored as the columns of the   
    M-by-NRHS right hand side matrix B and the N-by-NRHS solution matrix 
  
    X.   

    The effective rank of A is determined by treating as zero those   
    singular values which are less than RCOND times the largest singular 
  
    value.   

    Arguments   
    =========   

    M       (input) INTEGER   
            The number of rows of the matrix A. M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrix A. N >= 0.   

    NRHS    (input) INTEGER   
            The number of right hand sides, i.e., the number of columns   
            of the matrices B and X. NRHS >= 0.   

    A       (input/output) REAL array, dimension (LDA,N)   
            On entry, the M-by-N matrix A.   
            On exit, the first min(m,n) rows of A are overwritten with   
            its right singular vectors, stored rowwise.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(1,M).   

    B       (input/output) REAL array, dimension (LDB,NRHS)   
            On entry, the M-by-NRHS right hand side matrix B.   
            On exit, B is overwritten by the N-by-NRHS solution   
            matrix X.  If m >= n and RANK = n, the residual   
            sum-of-squares for the solution in the i-th column is given   
            by the sum of squares of elements n+1:m in that column.   

    LDB     (input) INTEGER   
            The leading dimension of the array B. LDB >= max(1,max(M,N)). 
  

    S       (output) REAL array, dimension (min(M,N))   
            The singular values of A in decreasing order.   
            The condition number of A in the 2-norm = S(1)/S(min(m,n)).   

    RCOND   (input) REAL   
            RCOND is used to determine the effective rank of A.   
            Singular values S(i) <= RCOND*S(1) are treated as zero.   
            If RCOND < 0, machine precision is used instead.   

    RANK    (output) INTEGER   
            The effective rank of A, i.e., the number of singular values 
  
            which are greater than RCOND*S(1).   

    WORK    (workspace/output) REAL array, dimension (LWORK)   
            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INTEGER   
            The dimension of the array WORK. LWORK >= 1, and also:   
            LWORK >= 3*min(M,N) + max( 2*min(M,N), max(M,N), NRHS )   
            For good performance, LWORK should generally be larger.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value.   
            > 0:  the algorithm for computing the SVD failed to converge; 
  
                  if INFO = i, i off-diagonal elements of an intermediate 
  
                  bidiagonal form did not converge to zero.   

    ===================================================================== 
  


       Test the input arguments   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static integer c__6 = 6;
    static integer c_n1 = -1;
    static integer c__1 = 1;
    static integer c__0 = 0;
    static real c_b74 = 0.f;
    static real c_b108 = 1.f;
    
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3, i__4;
    real r__1;
    /* Local variables */
    static real anrm, bnrm;
    static integer itau;
    static real vdum[1];
    static integer i, iascl, ibscl, chunk;
    extern /* Subroutine */ int sgemm_(char *, char *, integer *, integer *, 
	    integer *, real *, real *, integer *, real *, integer *, real *, 
	    real *, integer *);
    static real sfmin;
    static integer minmn, maxmn;
    extern /* Subroutine */ int sgemv_(char *, integer *, integer *, real *, 
	    real *, integer *, real *, integer *, real *, real *, integer *);
    static integer itaup, itauq;
    extern /* Subroutine */ int srscl_(integer *, real *, real *, integer *);
    static integer mnthr, iwork;
    extern /* Subroutine */ int scopy_(integer *, real *, integer *, real *, 
	    integer *);
    static integer bl, ie, il;
    extern /* Subroutine */ int slabad_(real *, real *);
    static integer mm, bdspac;
    extern /* Subroutine */ int sgebrd_(integer *, integer *, real *, integer 
	    *, real *, real *, real *, real *, real *, integer *, integer *);
    extern doublereal slamch_(char *), slange_(char *, integer *, 
	    integer *, real *, integer *, real *);
    extern /* Subroutine */ int xerbla_(char *, integer *);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static real bignum;
    extern /* Subroutine */ int sgelqf_(integer *, integer *, real *, integer 
	    *, real *, real *, integer *, integer *), slascl_(char *, integer 
	    *, integer *, real *, real *, integer *, integer *, real *, 
	    integer *, integer *), sgeqrf_(integer *, integer *, real 
	    *, integer *, real *, real *, integer *, integer *), slacpy_(char 
	    *, integer *, integer *, real *, integer *, real *, integer *), slaset_(char *, integer *, integer *, real *, real *, 
	    real *, integer *), sbdsqr_(char *, integer *, integer *, 
	    integer *, integer *, real *, real *, real *, integer *, real *, 
	    integer *, real *, integer *, real *, integer *), sorgbr_(
	    char *, integer *, integer *, integer *, real *, integer *, real *
	    , real *, integer *, integer *);
    static integer ldwork;
    extern /* Subroutine */ int sormbr_(char *, char *, char *, integer *, 
	    integer *, integer *, real *, integer *, real *, real *, integer *
	    , real *, integer *, integer *);
    static integer minwrk, maxwrk;
    static real smlnum;
    extern /* Subroutine */ int sormlq_(char *, char *, integer *, integer *, 
	    integer *, real *, integer *, real *, real *, integer *, real *, 
	    integer *, integer *), sormqr_(char *, char *, 
	    integer *, integer *, integer *, real *, integer *, real *, real *
	    , integer *, real *, integer *, integer *);
    static real eps, thr;



#define S(I) s[(I)-1]
#define WORK(I) work[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]
#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]

    *info = 0;
    minmn = min(*m,*n);
    maxmn = max(*m,*n);
    mnthr = ilaenv_(&c__6, "SGELSS", " ", m, n, nrhs, &c_n1, 6L, 1L);
    if (*m < 0) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*nrhs < 0) {
	*info = -3;
    } else if (*lda < max(1,*m)) {
	*info = -5;
    } else if (*ldb < max(1,maxmn)) {
	*info = -7;
    }

/*     Compute workspace   
        (Note: Comments in the code beginning "Workspace:" describe the   
         minimal amount of workspace needed at that point in the code,   
         as well as the preferred amount for good performance.   
         NB refers to the optimal block size for the immediately   
         following subroutine, as returned by ILAENV.) */

    minwrk = 1;
    if (*info == 0 && *lwork >= 1) {
	maxwrk = 0;
	mm = *m;
	if (*m >= *n && *m >= mnthr) {

/*           Path 1a - overdetermined, with many more rows than co
lumns */

	    mm = *n;
/* Computing MAX */
	    i__1 = maxwrk, i__2 = *n + *n * ilaenv_(&c__1, "SGEQRF", " ", m, 
		    n, &c_n1, &c_n1, 6L, 1L);
	    maxwrk = max(i__1,i__2);
/* Computing MAX */
	    i__1 = maxwrk, i__2 = *n + *nrhs * ilaenv_(&c__1, "SORMQR", "LT", 
		    m, nrhs, n, &c_n1, 6L, 2L);
	    maxwrk = max(i__1,i__2);
	}
	if (*m >= *n) {

/*           Path 1 - overdetermined or exactly determined   

             Compute workspace neede for SBDSQR   

   Computing MAX */
	    i__1 = 1, i__2 = *n * 5 - 4;
	    bdspac = max(i__1,i__2);
/* Computing MAX */
	    i__1 = maxwrk, i__2 = *n * 3 + (mm + *n) * ilaenv_(&c__1, "SGEBRD"
		    , " ", &mm, n, &c_n1, &c_n1, 6L, 1L);
	    maxwrk = max(i__1,i__2);
/* Computing MAX */
	    i__1 = maxwrk, i__2 = *n * 3 + *nrhs * ilaenv_(&c__1, "SORMBR", 
		    "QLT", &mm, nrhs, n, &c_n1, 6L, 3L);
	    maxwrk = max(i__1,i__2);
/* Computing MAX */
	    i__1 = maxwrk, i__2 = *n * 3 + (*n - 1) * ilaenv_(&c__1, "SORGBR",
		     "P", n, n, n, &c_n1, 6L, 1L);
	    maxwrk = max(i__1,i__2);
	    maxwrk = max(maxwrk,bdspac);
/* Computing MAX */
	    i__1 = maxwrk, i__2 = *n * *nrhs;
	    maxwrk = max(i__1,i__2);
/* Computing MAX */
	    i__1 = *n * 3 + mm, i__2 = *n * 3 + *nrhs, i__1 = max(i__1,i__2);
	    minwrk = max(i__1,bdspac);
	    maxwrk = max(minwrk,maxwrk);
	}
	if (*n > *m) {

/*           Compute workspace neede for SBDSQR   

   Computing MAX */
	    i__1 = 1, i__2 = *m * 5 - 4;
	    bdspac = max(i__1,i__2);
/* Computing MAX */
	    i__1 = *m * 3 + *nrhs, i__2 = *m * 3 + *n, i__1 = max(i__1,i__2);
	    minwrk = max(i__1,bdspac);
	    if (*n >= mnthr) {

/*              Path 2a - underdetermined, with many more colu
mns   
                than rows */

		maxwrk = *m + *m * ilaenv_(&c__1, "SGELQF", " ", m, n, &c_n1, 
			&c_n1, 6L, 1L);
/* Computing MAX */
		i__1 = maxwrk, i__2 = *m * *m + (*m << 2) + (*m << 1) * 
			ilaenv_(&c__1, "SGEBRD", " ", m, m, &c_n1, &c_n1, 6L, 
			1L);
		maxwrk = max(i__1,i__2);
/* Computing MAX */
		i__1 = maxwrk, i__2 = *m * *m + (*m << 2) + *nrhs * ilaenv_(&
			c__1, "SORMBR", "QLT", m, nrhs, m, &c_n1, 6L, 3L);
		maxwrk = max(i__1,i__2);
/* Computing MAX */
		i__1 = maxwrk, i__2 = *m * *m + (*m << 2) + (*m - 1) * 
			ilaenv_(&c__1, "SORGBR", "P", m, m, m, &c_n1, 6L, 1L);
		maxwrk = max(i__1,i__2);
/* Computing MAX */
		i__1 = maxwrk, i__2 = *m * *m + *m + bdspac;
		maxwrk = max(i__1,i__2);
		if (*nrhs > 1) {
/* Computing MAX */
		    i__1 = maxwrk, i__2 = *m * *m + *m + *m * *nrhs;
		    maxwrk = max(i__1,i__2);
		} else {
/* Computing MAX */
		    i__1 = maxwrk, i__2 = *m * *m + (*m << 1);
		    maxwrk = max(i__1,i__2);
		}
/* Computing MAX */
		i__1 = maxwrk, i__2 = *m + *nrhs * ilaenv_(&c__1, "SORMLQ", 
			"LT", n, nrhs, m, &c_n1, 6L, 2L);
		maxwrk = max(i__1,i__2);
	    } else {

/*              Path 2 - underdetermined */

		maxwrk = *m * 3 + (*n + *m) * ilaenv_(&c__1, "SGEBRD", " ", m,
			 n, &c_n1, &c_n1, 6L, 1L);
/* Computing MAX */
		i__1 = maxwrk, i__2 = *m * 3 + *nrhs * ilaenv_(&c__1, "SORMBR"
			, "QLT", m, nrhs, m, &c_n1, 6L, 3L);
		maxwrk = max(i__1,i__2);
/* Computing MAX */
		i__1 = maxwrk, i__2 = *m * 3 + *m * ilaenv_(&c__1, "SORGBR", 
			"P", m, n, m, &c_n1, 6L, 1L);
		maxwrk = max(i__1,i__2);
		maxwrk = max(maxwrk,bdspac);
/* Computing MAX */
		i__1 = maxwrk, i__2 = *n * *nrhs;
		maxwrk = max(i__1,i__2);
	    }
	}
	maxwrk = max(minwrk,maxwrk);
	WORK(1) = (real) maxwrk;
    }

    minwrk = max(minwrk,1);
    if (*lwork < minwrk) {
	*info = -12;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("SGELSS", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*m == 0 || *n == 0) {
	*rank = 0;
	return 0;
    }

/*     Get machine parameters */

    eps = slamch_("P");
    sfmin = slamch_("S");
    smlnum = sfmin / eps;
    bignum = 1.f / smlnum;
    slabad_(&smlnum, &bignum);

/*     Scale A if max element outside range [SMLNUM,BIGNUM] */

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
	slaset_("F", &i__1, nrhs, &c_b74, &c_b74, &B(1,1), ldb);
	slaset_("F", &minmn, &c__1, &c_b74, &c_b74, &S(1), &c__1);
	*rank = 0;
	goto L70;
    }

/*     Scale B if max element outside range [SMLNUM,BIGNUM] */

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

/*     Overdetermined case */

    if (*m >= *n) {

/*        Path 1 - overdetermined or exactly determined */

	mm = *m;
	if (*m >= mnthr) {

/*           Path 1a - overdetermined, with many more rows than co
lumns */

	    mm = *n;
	    itau = 1;
	    iwork = itau + *n;

/*           Compute A=Q*R   
             (Workspace: need 2*N, prefer N+N*NB) */

	    i__1 = *lwork - iwork + 1;
	    sgeqrf_(m, n, &A(1,1), lda, &WORK(itau), &WORK(iwork), &i__1,
		     info);

/*           Multiply B by transpose(Q)   
             (Workspace: need N+NRHS, prefer N+NRHS*NB) */

	    i__1 = *lwork - iwork + 1;
	    sormqr_("L", "T", m, nrhs, n, &A(1,1), lda, &WORK(itau), &B(1,1), ldb, &WORK(iwork), &i__1, info);

/*           Zero out below R */

	    if (*n > 1) {
		i__1 = *n - 1;
		i__2 = *n - 1;
		slaset_("L", &i__1, &i__2, &c_b74, &c_b74, &A(2,1), 
			lda);
	    }
	}

	ie = 1;
	itauq = ie + *n;
	itaup = itauq + *n;
	iwork = itaup + *n;

/*        Bidiagonalize R in A   
          (Workspace: need 3*N+MM, prefer 3*N+(MM+N)*NB) */

	i__1 = *lwork - iwork + 1;
	sgebrd_(&mm, n, &A(1,1), lda, &S(1), &WORK(ie), &WORK(itauq), &
		WORK(itaup), &WORK(iwork), &i__1, info);

/*        Multiply B by transpose of left bidiagonalizing vectors of R
   
          (Workspace: need 3*N+NRHS, prefer 3*N+NRHS*NB) */

	i__1 = *lwork - iwork + 1;
	sormbr_("Q", "L", "T", &mm, nrhs, n, &A(1,1), lda, &WORK(itauq), 
		&B(1,1), ldb, &WORK(iwork), &i__1, info);

/*        Generate right bidiagonalizing vectors of R in A   
          (Workspace: need 4*N-1, prefer 3*N+(N-1)*NB) */

	i__1 = *lwork - iwork + 1;
	sorgbr_("P", n, n, n, &A(1,1), lda, &WORK(itaup), &WORK(iwork), &
		i__1, info);
	iwork = ie + *n;

/*        Perform bidiagonal QR iteration   
            multiply B by transpose of left singular vectors   
            compute right singular vectors in A   
          (Workspace: need BDSPAC) */

	sbdsqr_("U", n, n, &c__0, nrhs, &S(1), &WORK(ie), &A(1,1), lda, 
		vdum, &c__1, &B(1,1), ldb, &WORK(iwork), info);
	if (*info != 0) {
	    goto L70;
	}

/*        Multiply B by reciprocals of singular values   

   Computing MAX */
	r__1 = *rcond * S(1);
	thr = dmax(r__1,sfmin);
	if (*rcond < 0.f) {
/* Computing MAX */
	    r__1 = eps * S(1);
	    thr = dmax(r__1,sfmin);
	}
	*rank = 0;
	i__1 = *n;
	for (i = 1; i <= *n; ++i) {
	    if (S(i) > thr) {
		srscl_(nrhs, &S(i), &B(i,1), ldb);
		++(*rank);
	    } else {
		slaset_("F", &c__1, nrhs, &c_b74, &c_b74, &B(i,1), ldb);
	    }
/* L10: */
	}

/*        Multiply B by right singular vectors   
          (Workspace: need N, prefer N*NRHS) */

	if (*lwork >= *ldb * *nrhs && *nrhs > 1) {
	    sgemm_("T", "N", n, nrhs, n, &c_b108, &A(1,1), lda, &B(1,1), ldb, &c_b74, &WORK(1), ldb);
	    slacpy_("G", n, nrhs, &WORK(1), ldb, &B(1,1), ldb);
	} else if (*nrhs > 1) {
	    chunk = *lwork / *n;
	    i__1 = *nrhs;
	    i__2 = chunk;
	    for (i = 1; chunk < 0 ? i >= *nrhs : i <= *nrhs; i += chunk) {
/* Computing MIN */
		i__3 = *nrhs - i + 1;
		bl = min(i__3,chunk);
		sgemm_("T", "N", n, &bl, n, &c_b108, &A(1,1), lda, &B(1,1), ldb, &c_b74, &WORK(1), n);
		slacpy_("G", n, &bl, &WORK(1), n, &B(1,1), ldb);
/* L20: */
	    }
	} else {
	    sgemv_("T", n, n, &c_b108, &A(1,1), lda, &B(1,1), &c__1,
		     &c_b74, &WORK(1), &c__1);
	    scopy_(n, &WORK(1), &c__1, &B(1,1), &c__1);
	}

    } else /* if(complicated condition) */ {
/* Computing MAX */
	i__2 = *m, i__1 = (*m << 1) - 4, i__2 = max(i__2,i__1), i__2 = max(
		i__2,*nrhs), i__1 = *n - *m * 3;
	if (*n >= mnthr && *lwork >= (*m << 2) + *m * *m + max(i__2,i__1)) {

/*        Path 2a - underdetermined, with many more columns than r
ows   
          and sufficient workspace for an efficient algorithm */

	    ldwork = *m;
/* Computing MAX   
   Computing MAX */
	    i__3 = *m, i__4 = (*m << 1) - 4, i__3 = max(i__3,i__4), i__3 = 
		    max(i__3,*nrhs), i__4 = *n - *m * 3;
	    i__2 = (*m << 2) + *m * *lda + max(i__3,i__4), i__1 = *m * *lda + 
		    *m + *m * *nrhs;
	    if (*lwork >= max(i__2,i__1)) {
		ldwork = *lda;
	    }
	    itau = 1;
	    iwork = *m + 1;

/*        Compute A=L*Q   
          (Workspace: need 2*M, prefer M+M*NB) */

	    i__2 = *lwork - iwork + 1;
	    sgelqf_(m, n, &A(1,1), lda, &WORK(itau), &WORK(iwork), &i__2,
		     info);
	    il = iwork;

/*        Copy L to WORK(IL), zeroing out above it */

	    slacpy_("L", m, m, &A(1,1), lda, &WORK(il), &ldwork);
	    i__2 = *m - 1;
	    i__1 = *m - 1;
	    slaset_("U", &i__2, &i__1, &c_b74, &c_b74, &WORK(il + ldwork), &
		    ldwork);
	    ie = il + ldwork * *m;
	    itauq = ie + *m;
	    itaup = itauq + *m;
	    iwork = itaup + *m;

/*        Bidiagonalize L in WORK(IL)   
          (Workspace: need M*M+5*M, prefer M*M+4*M+2*M*NB) */

	    i__2 = *lwork - iwork + 1;
	    sgebrd_(m, m, &WORK(il), &ldwork, &S(1), &WORK(ie), &WORK(itauq), 
		    &WORK(itaup), &WORK(iwork), &i__2, info);

/*        Multiply B by transpose of left bidiagonalizing vectors 
of L   
          (Workspace: need M*M+4*M+NRHS, prefer M*M+4*M+NRHS*NB) 
*/

	    i__2 = *lwork - iwork + 1;
	    sormbr_("Q", "L", "T", m, nrhs, m, &WORK(il), &ldwork, &WORK(
		    itauq), &B(1,1), ldb, &WORK(iwork), &i__2, info);

/*        Generate right bidiagonalizing vectors of R in WORK(IL) 
  
          (Workspace: need M*M+5*M-1, prefer M*M+4*M+(M-1)*NB) */

	    i__2 = *lwork - iwork + 1;
	    sorgbr_("P", m, m, m, &WORK(il), &ldwork, &WORK(itaup), &WORK(
		    iwork), &i__2, info);
	    iwork = ie + *m;

/*        Perform bidiagonal QR iteration,   
             computing right singular vectors of L in WORK(IL) and
   
             multiplying B by transpose of left singular vectors 
  
          (Workspace: need M*M+M+BDSPAC) */

	    sbdsqr_("U", m, m, &c__0, nrhs, &S(1), &WORK(ie), &WORK(il), &
		    ldwork, &A(1,1), lda, &B(1,1), ldb, &WORK(iwork)
		    , info);
	    if (*info != 0) {
		goto L70;
	    }

/*        Multiply B by reciprocals of singular values   

   Computing MAX */
	    r__1 = *rcond * S(1);
	    thr = dmax(r__1,sfmin);
	    if (*rcond < 0.f) {
/* Computing MAX */
		r__1 = eps * S(1);
		thr = dmax(r__1,sfmin);
	    }
	    *rank = 0;
	    i__2 = *m;
	    for (i = 1; i <= *m; ++i) {
		if (S(i) > thr) {
		    srscl_(nrhs, &S(i), &B(i,1), ldb);
		    ++(*rank);
		} else {
		    slaset_("F", &c__1, nrhs, &c_b74, &c_b74, &B(i,1), 
			    ldb);
		}
/* L30: */
	    }
	    iwork = ie;

/*        Multiply B by right singular vectors of L in WORK(IL)   
          (Workspace: need M*M+2*M, prefer M*M+M+M*NRHS) */

	    if (*lwork >= *ldb * *nrhs + iwork - 1 && *nrhs > 1) {
		sgemm_("T", "N", m, nrhs, m, &c_b108, &WORK(il), &ldwork, &B(1,1), ldb, &c_b74, &WORK(iwork), ldb);
		slacpy_("G", m, nrhs, &WORK(iwork), ldb, &B(1,1), ldb);
	    } else if (*nrhs > 1) {
		chunk = (*lwork - iwork + 1) / *m;
		i__2 = *nrhs;
		i__1 = chunk;
		for (i = 1; chunk < 0 ? i >= *nrhs : i <= *nrhs; i += chunk) {
/* Computing MIN */
		    i__3 = *nrhs - i + 1;
		    bl = min(i__3,chunk);
		    sgemm_("T", "N", m, &bl, m, &c_b108, &WORK(il), &ldwork, &
			    B(1,i), ldb, &c_b74, &WORK(iwork), n);
		    slacpy_("G", m, &bl, &WORK(iwork), n, &B(1,1), ldb);
/* L40: */
		}
	    } else {
		sgemv_("T", m, m, &c_b108, &WORK(il), &ldwork, &B(1,1),
			 &c__1, &c_b74, &WORK(iwork), &c__1);
		scopy_(m, &WORK(iwork), &c__1, &B(1,1), &c__1);
	    }

/*        Zero out below first M rows of B */

	    i__1 = *n - *m;
	    slaset_("F", &i__1, nrhs, &c_b74, &c_b74, &B(*m+1,1), 
		    ldb);
	    iwork = itau + *m;

/*        Multiply transpose(Q) by B   
          (Workspace: need M+NRHS, prefer M+NRHS*NB) */

	    i__1 = *lwork - iwork + 1;
	    sormlq_("L", "T", n, nrhs, m, &A(1,1), lda, &WORK(itau), &B(1,1), ldb, &WORK(iwork), &i__1, info);

	} else {

/*        Path 2 - remaining underdetermined cases */

	    ie = 1;
	    itauq = ie + *m;
	    itaup = itauq + *m;
	    iwork = itaup + *m;

/*        Bidiagonalize A   
          (Workspace: need 3*M+N, prefer 3*M+(M+N)*NB) */

	    i__1 = *lwork - iwork + 1;
	    sgebrd_(m, n, &A(1,1), lda, &S(1), &WORK(ie), &WORK(itauq), &
		    WORK(itaup), &WORK(iwork), &i__1, info);

/*        Multiply B by transpose of left bidiagonalizing vectors 
  
          (Workspace: need 3*M+NRHS, prefer 3*M+NRHS*NB) */

	    i__1 = *lwork - iwork + 1;
	    sormbr_("Q", "L", "T", m, nrhs, n, &A(1,1), lda, &WORK(itauq)
		    , &B(1,1), ldb, &WORK(iwork), &i__1, info);

/*        Generate right bidiagonalizing vectors in A   
          (Workspace: need 4*M, prefer 3*M+M*NB) */

	    i__1 = *lwork - iwork + 1;
	    sorgbr_("P", m, n, m, &A(1,1), lda, &WORK(itaup), &WORK(
		    iwork), &i__1, info);
	    iwork = ie + *m;

/*        Perform bidiagonal QR iteration,   
             computing right singular vectors of A in A and   
             multiplying B by transpose of left singular vectors 
  
          (Workspace: need BDSPAC) */

	    sbdsqr_("L", m, n, &c__0, nrhs, &S(1), &WORK(ie), &A(1,1), 
		    lda, vdum, &c__1, &B(1,1), ldb, &WORK(iwork), info);
	    if (*info != 0) {
		goto L70;
	    }

/*        Multiply B by reciprocals of singular values   

   Computing MAX */
	    r__1 = *rcond * S(1);
	    thr = dmax(r__1,sfmin);
	    if (*rcond < 0.f) {
/* Computing MAX */
		r__1 = eps * S(1);
		thr = dmax(r__1,sfmin);
	    }
	    *rank = 0;
	    i__1 = *m;
	    for (i = 1; i <= *m; ++i) {
		if (S(i) > thr) {
		    srscl_(nrhs, &S(i), &B(i,1), ldb);
		    ++(*rank);
		} else {
		    slaset_("F", &c__1, nrhs, &c_b74, &c_b74, &B(i,1), 
			    ldb);
		}
/* L50: */
	    }

/*        Multiply B by right singular vectors of A   
          (Workspace: need N, prefer N*NRHS) */

	    if (*lwork >= *ldb * *nrhs && *nrhs > 1) {
		sgemm_("T", "N", n, nrhs, m, &c_b108, &A(1,1), lda, &B(1,1), ldb, &c_b74, &WORK(1), ldb);
		slacpy_("F", n, nrhs, &WORK(1), ldb, &B(1,1), ldb);
	    } else if (*nrhs > 1) {
		chunk = *lwork / *n;
		i__1 = *nrhs;
		i__2 = chunk;
		for (i = 1; chunk < 0 ? i >= *nrhs : i <= *nrhs; i += chunk) {
/* Computing MIN */
		    i__3 = *nrhs - i + 1;
		    bl = min(i__3,chunk);
		    sgemm_("T", "N", n, &bl, m, &c_b108, &A(1,1), lda, &
			    B(1,i), ldb, &c_b74, &WORK(1), n);
		    slacpy_("F", n, &bl, &WORK(1), n, &B(1,i), ldb);
/* L60: */
		}
	    } else {
		sgemv_("T", m, n, &c_b108, &A(1,1), lda, &B(1,1), &
			c__1, &c_b74, &WORK(1), &c__1);
		scopy_(n, &WORK(1), &c__1, &B(1,1), &c__1);
	    }
	}
    }

/*     Undo scaling */

    if (iascl == 1) {
	slascl_("G", &c__0, &c__0, &anrm, &smlnum, n, nrhs, &B(1,1), ldb,
		 info);
	slascl_("G", &c__0, &c__0, &smlnum, &anrm, &minmn, &c__1, &S(1), &
		minmn, info);
    } else if (iascl == 2) {
	slascl_("G", &c__0, &c__0, &anrm, &bignum, n, nrhs, &B(1,1), ldb,
		 info);
	slascl_("G", &c__0, &c__0, &bignum, &anrm, &minmn, &c__1, &S(1), &
		minmn, info);
    }
    if (ibscl == 1) {
	slascl_("G", &c__0, &c__0, &smlnum, &bnrm, n, nrhs, &B(1,1), ldb,
		 info);
    } else if (ibscl == 2) {
	slascl_("G", &c__0, &c__0, &bignum, &bnrm, n, nrhs, &B(1,1), ldb,
		 info);
    }

L70:
    WORK(1) = (real) maxwrk;
    return 0;

/*     End of SGELSS */

} /* sgelss_ */

