#include "f2c.h"

/* Subroutine */ int sgegs_(char *jobvsl, char *jobvsr, integer *n, real *a, 
	integer *lda, real *b, integer *ldb, real *alphar, real *alphai, real 
	*beta, real *vsl, integer *ldvsl, real *vsr, integer *ldvsr, real *
	work, integer *lwork, integer *info)
{
/*  -- LAPACK driver routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    SGEGS computes for a pair of N-by-N real nonsymmetric matrices A, B: 
  
    the generalized eigenvalues (alphar +/- alphai*i, beta), the real   
    Schur form (A, B), and optionally left and/or right Schur vectors   
    (VSL and VSR).   

    (If only the generalized eigenvalues are needed, use the driver SGEGV 
  
    instead.)   

    A generalized eigenvalue for a pair of matrices (A,B) is, roughly   
    speaking, a scalar w or a ratio  alpha/beta = w, such that  A - w*B   
    is singular.  It is usually represented as the pair (alpha,beta),   
    as there is a reasonable interpretation for beta=0, and even for   
    both being zero.  A good beginning reference is the book, "VISMatrix   
    Computations", by G. Golub & C. van Loan (Johns Hopkins U. Press)   

    The (generalized) Schur form of a pair of matrices is the result of   
    multiplying both matrices on the left by one orthogonal matrix and   
    both on the right by another orthogonal matrix, these two orthogonal 
  
    matrices being chosen so as to bring the pair of matrices into   
    (real) Schur form.   

    A pair of matrices A, B is in generalized real Schur form if B is   
    upper triangular with non-negative diagonal and A is block upper   
    triangular with 1-by-1 and 2-by-2 blocks.  1-by-1 blocks correspond   
    to real generalized eigenvalues, while 2-by-2 blocks of A will be   
    "standardized" by making the corresponding elements of B have the   
    form:   
            [  a  0  ]   
            [  0  b  ]   

    and the pair of corresponding 2-by-2 blocks in A and B will   
    have a complex conjugate pair of generalized eigenvalues.   

    The left and right Schur vectors are the columns of VSL and VSR,   
    respectively, where VSL and VSR are the orthogonal matrices   
    which reduce A and B to Schur form:   

    Schur form of (A,B) = ( (VSL)**T A (VSR), (VSL)**T B (VSR) )   

    Arguments   
    =========   

    JOBVSL  (input) CHARACTER*1   
            = 'N':  do not compute the left Schur vectors;   
            = 'V':  compute the left Schur vectors.   

    JOBVSR  (input) CHARACTER*1   
            = 'N':  do not compute the right Schur vectors;   
            = 'V':  compute the right Schur vectors.   

    N       (input) INTEGER   
            The order of the matrices A, B, VSL, and VSR.  N >= 0.   

    A       (input/output) REAL array, dimension (LDA, N)   
            On entry, the first of the pair of matrices whose generalized 
  
            eigenvalues and (optionally) Schur vectors are to be   
            computed.   
            On exit, the generalized Schur form of A.   
            Note: to avoid overflow, the Frobenius norm of the matrix   
            A should be less than the overflow threshold.   

    LDA     (input) INTEGER   
            The leading dimension of A.  LDA >= max(1,N).   

    B       (input/output) REAL array, dimension (LDB, N)   
            On entry, the second of the pair of matrices whose   
            generalized eigenvalues and (optionally) Schur vectors are   
            to be computed.   
            On exit, the generalized Schur form of B.   
            Note: to avoid overflow, the Frobenius norm of the matrix   
            B should be less than the overflow threshold.   

    LDB     (input) INTEGER   
            The leading dimension of B.  LDB >= max(1,N).   

    ALPHAR  (output) REAL array, dimension (N)   
    ALPHAI  (output) REAL array, dimension (N)   
    BETA    (output) REAL array, dimension (N)   
            On exit, (ALPHAR(j) + ALPHAI(j)*i)/BETA(j), j=1,...,N, will   
            be the generalized eigenvalues.  ALPHAR(j) + ALPHAI(j)*i,   
            j=1,...,N  and  BETA(j),j=1,...,N  are the diagonals of the   
            complex Schur form (A,B) that would result if the 2-by-2   
            diagonal blocks of the real Schur form of (A,B) were further 
  
            reduced to triangular form using 2-by-2 complex unitary   
            transformations.  If ALPHAI(j) is zero, then the j-th   
            eigenvalue is real; if positive, then the j-th and (j+1)-st   
            eigenvalues are a complex conjugate pair, with ALPHAI(j+1)   
            negative.   

            Note: the quotients ALPHAR(j)/BETA(j) and ALPHAI(j)/BETA(j)   
            may easily over- or underflow, and BETA(j) may even be zero. 
  
            Thus, the user should avoid naively computing the ratio   
            alpha/beta.  However, ALPHAR and ALPHAI will be always less   
            than and usually comparable with norm(A) in magnitude, and   
            BETA always less than and usually comparable with norm(B).   

    VSL     (output) REAL array, dimension (LDVSL,N)   
            If JOBVSL = 'V', VSL will contain the left Schur vectors.   
            (See "Purpose", above.)   
            Not referenced if JOBVSL = 'N'.   

    LDVSL   (input) INTEGER   
            The leading dimension of the matrix VSL. LDVSL >=1, and   
            if JOBVSL = 'V', LDVSL >= N.   

    VSR     (output) REAL array, dimension (LDVSR,N)   
            If JOBVSR = 'V', VSR will contain the right Schur vectors.   
            (See "Purpose", above.)   
            Not referenced if JOBVSR = 'N'.   

    LDVSR   (input) INTEGER   
            The leading dimension of the matrix VSR. LDVSR >= 1, and   
            if JOBVSR = 'V', LDVSR >= N.   

    WORK    (workspace/output) REAL array, dimension (LWORK)   
            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INTEGER   
            The dimension of the array WORK.  LWORK >= max(1,4*N).   
            For good performance, LWORK must generally be larger.   
            To compute the optimal value of LWORK, call ILAENV to get   
            blocksizes (for SGEQRF, SORMQR, and SORGQR.)  Then compute:   
            NB  -- MAX of the blocksizes for SGEQRF, SORMQR, and SORGQR   
            The optimal LWORK is  2*N + N*(NB+1).   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value.   
            = 1,...,N:   
                  The QZ iteration failed.  (A,B) are not in Schur   
                  form, but ALPHAR(j), ALPHAI(j), and BETA(j) should   
                  be correct for j=INFO+1,...,N.   
            > N:  errors that usually indicate LAPACK problems:   
                  =N+1: error return from SGGBAL   
                  =N+2: error return from SGEQRF   
                  =N+3: error return from SORMQR   
                  =N+4: error return from SORGQR   
                  =N+5: error return from SGGHRD   
                  =N+6: error return from SHGEQZ (other than failed   
                                                  iteration)   
                  =N+7: error return from SGGBAK (computing VSL)   
                  =N+8: error return from SGGBAK (computing VSR)   
                  =N+9: error return from SLASCL (various places)   

    ===================================================================== 
  


       Decode the input arguments   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static integer c_n1 = -1;
    static real c_b23 = 0.f;
    static real c_b24 = 1.f;
    static integer c__1 = 1;
    
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, vsl_dim1, vsl_offset, 
	    vsr_dim1, vsr_offset, i__1, i__2;
    /* Local variables */
    static real anrm, bnrm;
    static integer itau;
    extern logical lsame_(char *, char *);
    static integer ileft, iinfo, icols;
    static logical ilvsl;
    static integer iwork;
    static logical ilvsr;
    static integer irows;
    extern /* Subroutine */ int sggbak_(char *, char *, integer *, integer *, 
	    integer *, real *, real *, integer *, real *, integer *, integer *
	    ), sggbal_(char *, integer *, real *, integer *, 
	    real *, integer *, integer *, integer *, real *, real *, real *, 
	    integer *);
    static logical ilascl, ilbscl;
    extern doublereal slamch_(char *), slange_(char *, integer *, 
	    integer *, real *, integer *, real *);
    static real safmin;
    extern /* Subroutine */ int sgghrd_(char *, char *, integer *, integer *, 
	    integer *, real *, integer *, real *, integer *, real *, integer *
	    , real *, integer *, integer *), xerbla_(char *, 
	    integer *);
    static real bignum;
    extern /* Subroutine */ int slascl_(char *, integer *, integer *, real *, 
	    real *, integer *, integer *, real *, integer *, integer *);
    static integer ijobvl, iright;
    extern /* Subroutine */ int sgeqrf_(integer *, integer *, real *, integer 
	    *, real *, real *, integer *, integer *);
    static integer ijobvr;
    extern /* Subroutine */ int slacpy_(char *, integer *, integer *, real *, 
	    integer *, real *, integer *), slaset_(char *, integer *, 
	    integer *, real *, real *, real *, integer *);
    static real anrmto;
    static integer lwkmin;
    static real bnrmto;
    extern /* Subroutine */ int shgeqz_(char *, char *, char *, integer *, 
	    integer *, integer *, real *, integer *, real *, integer *, real *
	    , real *, real *, real *, integer *, real *, integer *, real *, 
	    integer *, integer *);
    static real smlnum;
    extern /* Subroutine */ int sorgqr_(integer *, integer *, integer *, real 
	    *, integer *, real *, real *, integer *, integer *);
    static integer lwkopt;
    extern /* Subroutine */ int sormqr_(char *, char *, integer *, integer *, 
	    integer *, real *, integer *, real *, real *, integer *, real *, 
	    integer *, integer *);
    static integer ihi, ilo;
    static real eps;



#define ALPHAR(I) alphar[(I)-1]
#define ALPHAI(I) alphai[(I)-1]
#define BETA(I) beta[(I)-1]
#define WORK(I) work[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]
#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]
#define VSL(I,J) vsl[(I)-1 + ((J)-1)* ( *ldvsl)]
#define VSR(I,J) vsr[(I)-1 + ((J)-1)* ( *ldvsr)]

    if (lsame_(jobvsl, "N")) {
	ijobvl = 1;
	ilvsl = FALSE_;
    } else if (lsame_(jobvsl, "V")) {
	ijobvl = 2;
	ilvsl = TRUE_;
    } else {
	ijobvl = -1;
	ilvsl = FALSE_;
    }

    if (lsame_(jobvsr, "N")) {
	ijobvr = 1;
	ilvsr = FALSE_;
    } else if (lsame_(jobvsr, "V")) {
	ijobvr = 2;
	ilvsr = TRUE_;
    } else {
	ijobvr = -1;
	ilvsr = FALSE_;
    }

/*     Test the input arguments   

   Computing MAX */
    i__1 = *n << 2;
    lwkmin = max(i__1,1);
    lwkopt = lwkmin;
    *info = 0;
    if (ijobvl <= 0) {
	*info = -1;
    } else if (ijobvr <= 0) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*lda < max(1,*n)) {
	*info = -5;
    } else if (*ldb < max(1,*n)) {
	*info = -7;
    } else if (*ldvsl < 1 || ilvsl && *ldvsl < *n) {
	*info = -12;
    } else if (*ldvsr < 1 || ilvsr && *ldvsr < *n) {
	*info = -14;
    } else if (*lwork < lwkmin) {
	*info = -16;
    }

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("SGEGS ", &i__1);
	return 0;
    }

/*     Quick return if possible */

    WORK(1) = (real) lwkopt;
    if (*n == 0) {
	return 0;
    }

/*     Get machine constants */

    eps = slamch_("E") * slamch_("B");
    safmin = slamch_("S");
    smlnum = *n * safmin / eps;
    bignum = 1.f / smlnum;

/*     Scale A if max element outside range [SMLNUM,BIGNUM] */

    anrm = slange_("M", n, n, &A(1,1), lda, &WORK(1));
    ilascl = FALSE_;
    if (anrm > 0.f && anrm < smlnum) {
	anrmto = smlnum;
	ilascl = TRUE_;
    } else if (anrm > bignum) {
	anrmto = bignum;
	ilascl = TRUE_;
    }

    if (ilascl) {
	slascl_("G", &c_n1, &c_n1, &anrm, &anrmto, n, n, &A(1,1), lda, &
		iinfo);
	if (iinfo != 0) {
	    *info = *n + 9;
	    return 0;
	}
    }

/*     Scale B if max element outside range [SMLNUM,BIGNUM] */

    bnrm = slange_("M", n, n, &B(1,1), ldb, &WORK(1));
    ilbscl = FALSE_;
    if (bnrm > 0.f && bnrm < smlnum) {
	bnrmto = smlnum;
	ilbscl = TRUE_;
    } else if (bnrm > bignum) {
	bnrmto = bignum;
	ilbscl = TRUE_;
    }

    if (ilbscl) {
	slascl_("G", &c_n1, &c_n1, &bnrm, &bnrmto, n, n, &B(1,1), ldb, &
		iinfo);
	if (iinfo != 0) {
	    *info = *n + 9;
	    return 0;
	}
    }

/*     Permute the matrix to make it more nearly triangular   
       Workspace layout:  (2*N words -- "work..." not actually used)   
          left_permutation, right_permutation, work... */

    ileft = 1;
    iright = *n + 1;
    iwork = iright + *n;
    sggbal_("P", n, &A(1,1), lda, &B(1,1), ldb, &ilo, &ihi, &WORK(
	    ileft), &WORK(iright), &WORK(iwork), &iinfo);
    if (iinfo != 0) {
	*info = *n + 1;
	goto L10;
    }

/*     Reduce B to triangular form, and initialize VSL and/or VSR   
       Workspace layout:  ("work..." must have at least N words)   
          left_permutation, right_permutation, tau, work... */

    irows = ihi + 1 - ilo;
    icols = *n + 1 - ilo;
    itau = iwork;
    iwork = itau + irows;
    i__1 = *lwork + 1 - iwork;
    sgeqrf_(&irows, &icols, &B(ilo,ilo), ldb, &WORK(itau), &WORK(
	    iwork), &i__1, &iinfo);
    if (iinfo >= 0) {
/* Computing MAX */
	i__1 = lwkopt, i__2 = (integer) WORK(iwork) + iwork - 1;
	lwkopt = max(i__1,i__2);
    }
    if (iinfo != 0) {
	*info = *n + 2;
	goto L10;
    }

    i__1 = *lwork + 1 - iwork;
    sormqr_("L", "T", &irows, &icols, &irows, &B(ilo,ilo), ldb, &
	    WORK(itau), &A(ilo,ilo), lda, &WORK(iwork), &i__1, &
	    iinfo);
    if (iinfo >= 0) {
/* Computing MAX */
	i__1 = lwkopt, i__2 = (integer) WORK(iwork) + iwork - 1;
	lwkopt = max(i__1,i__2);
    }
    if (iinfo != 0) {
	*info = *n + 3;
	goto L10;
    }

    if (ilvsl) {
	slaset_("Full", n, n, &c_b23, &c_b24, &VSL(1,1), ldvsl);
	i__1 = irows - 1;
	i__2 = irows - 1;
	slacpy_("L", &i__1, &i__2, &B(ilo+1,ilo), ldb, &VSL(ilo+1,ilo), ldvsl);
	i__1 = *lwork + 1 - iwork;
	sorgqr_(&irows, &irows, &irows, &VSL(ilo,ilo), ldvsl, &
		WORK(itau), &WORK(iwork), &i__1, &iinfo);
	if (iinfo >= 0) {
/* Computing MAX */
	    i__1 = lwkopt, i__2 = (integer) WORK(iwork) + iwork - 1;
	    lwkopt = max(i__1,i__2);
	}
	if (iinfo != 0) {
	    *info = *n + 4;
	    goto L10;
	}
    }

    if (ilvsr) {
	slaset_("Full", n, n, &c_b23, &c_b24, &VSR(1,1), ldvsr);
    }

/*     Reduce to generalized Hessenberg form */

    sgghrd_(jobvsl, jobvsr, n, &ilo, &ihi, &A(1,1), lda, &B(1,1), 
	    ldb, &VSL(1,1), ldvsl, &VSR(1,1), ldvsr, &iinfo);
    if (iinfo != 0) {
	*info = *n + 5;
	goto L10;
    }

/*     Perform QZ algorithm, computing Schur vectors if desired   
       Workspace layout:  ("work..." must have at least 1 word)   
          left_permutation, right_permutation, work... */

    iwork = itau;
    i__1 = *lwork + 1 - iwork;
    shgeqz_("S", jobvsl, jobvsr, n, &ilo, &ihi, &A(1,1), lda, &B(1,1), ldb, &ALPHAR(1), &ALPHAI(1), &BETA(1), &VSL(1,1)
	    , ldvsl, &VSR(1,1), ldvsr, &WORK(iwork), &i__1, &iinfo);
    if (iinfo >= 0) {
/* Computing MAX */
	i__1 = lwkopt, i__2 = (integer) WORK(iwork) + iwork - 1;
	lwkopt = max(i__1,i__2);
    }
    if (iinfo != 0) {
	if (iinfo > 0 && iinfo <= *n) {
	    *info = iinfo;
	} else if (iinfo > *n && iinfo <= *n << 1) {
	    *info = iinfo - *n;
	} else {
	    *info = *n + 6;
	}
	goto L10;
    }

/*     Apply permutation to VSL and VSR */

    if (ilvsl) {
	sggbak_("P", "L", n, &ilo, &ihi, &WORK(ileft), &WORK(iright), n, &VSL(1,1), ldvsl, &iinfo);
	if (iinfo != 0) {
	    *info = *n + 7;
	    goto L10;
	}
    }
    if (ilvsr) {
	sggbak_("P", "R", n, &ilo, &ihi, &WORK(ileft), &WORK(iright), n, &VSR(1,1), ldvsr, &iinfo);
	if (iinfo != 0) {
	    *info = *n + 8;
	    goto L10;
	}
    }

/*     Undo scaling */

    if (ilascl) {
	slascl_("H", &c_n1, &c_n1, &anrmto, &anrm, n, n, &A(1,1), lda, &
		iinfo);
	if (iinfo != 0) {
	    *info = *n + 9;
	    return 0;
	}
	slascl_("G", &c_n1, &c_n1, &anrmto, &anrm, n, &c__1, &ALPHAR(1), n, &
		iinfo);
	if (iinfo != 0) {
	    *info = *n + 9;
	    return 0;
	}
	slascl_("G", &c_n1, &c_n1, &anrmto, &anrm, n, &c__1, &ALPHAI(1), n, &
		iinfo);
	if (iinfo != 0) {
	    *info = *n + 9;
	    return 0;
	}
    }

    if (ilbscl) {
	slascl_("U", &c_n1, &c_n1, &bnrmto, &bnrm, n, n, &B(1,1), ldb, &
		iinfo);
	if (iinfo != 0) {
	    *info = *n + 9;
	    return 0;
	}
	slascl_("G", &c_n1, &c_n1, &bnrmto, &bnrm, n, &c__1, &BETA(1), n, &
		iinfo);
	if (iinfo != 0) {
	    *info = *n + 9;
	    return 0;
	}
    }

L10:
    WORK(1) = (real) lwkopt;

    return 0;

/*     End of SGEGS */

} /* sgegs_ */

