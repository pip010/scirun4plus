#include "f2c.h"

/* Subroutine */ int cgegs_(char *jobvsl, char *jobvsr, integer *n, complex *
	a, integer *lda, complex *b, integer *ldb, complex *alpha, complex *
	beta, complex *vsl, integer *ldvsl, complex *vsr, integer *ldvsr, 
	complex *work, integer *lwork, real *rwork, integer *info)
{
/*  -- LAPACK driver routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    SGEGS computes for a pair of N-by-N complex nonsymmetric matrices A, 
  
    B:  the generalized eigenvalues (alpha, beta), the complex Schur   
    form (A, B), and optionally left and/or right Schur vectors   
    (VSL and VSR).   

    (If only the generalized eigenvalues are needed, use the driver CGEGV 
  
    instead.)   

    A generalized eigenvalue for a pair of matrices (A,B) is, roughly   
    speaking, a scalar w or a ratio  alpha/beta = w, such that  A - w*B   
    is singular.  It is usually represented as the pair (alpha,beta),   
    as there is a reasonable interpretation for beta=0, and even for   
    both being zero.  A good beginning reference is the book, "VISMatrix   
    Computations", by G. Golub & C. van Loan (Johns Hopkins U. Press)   

    The (generalized) Schur form of a pair of matrices is the result of   
    multiplying both matrices on the left by one unitary matrix and   
    both on the right by another unitary matrix, these two unitary   
    matrices being chosen so as to bring the pair of matrices into   
    upper triangular form with the diagonal elements of B being   
    non-negative real numbers (this is also called complex Schur form.)   

    The left and right Schur vectors are the columns of VSL and VSR,   
    respectively, where VSL and VSR are the unitary matrices   
    which reduce A and B to Schur form:   

    Schur form of (A,B) = ( (VSL)**H A (VSR), (VSL)**H B (VSR) )   

    Arguments   
    =========   

    JOBVSL   (input) CHARACTER*1   
            = 'N':  do not compute the left Schur vectors;   
            = 'V':  compute the left Schur vectors.   

    JOBVSR   (input) CHARACTER*1   
            = 'N':  do not compute the right Schur vectors;   
            = 'V':  compute the right Schur vectors.   

    N       (input) INTEGER   
            The order of the matrices A, B, VSL, and VSR.  N >= 0.   

    A       (input/output) COMPLEX array, dimension (LDA, N)   
            On entry, the first of the pair of matrices whose generalized 
  
            eigenvalues and (optionally) Schur vectors are to be   
            computed.   
            On exit, the generalized Schur form of A.   

    LDA     (input) INTEGER   
            The leading dimension of A.  LDA >= max(1,N).   

    B       (input/output) COMPLEX array, dimension (LDB, N)   
            On entry, the second of the pair of matrices whose   
            generalized eigenvalues and (optionally) Schur vectors are   
            to be computed.   
            On exit, the generalized Schur form of B.   

    LDB     (input) INTEGER   
            The leading dimension of B.  LDB >= max(1,N).   

    ALPHA   (output) COMPLEX array, dimension (N)   
    BETA    (output) COMPLEX array, dimension (N)   
            On exit,  ALPHA(j)/BETA(j), j=1,...,N, will be the   
            generalized eigenvalues.  ALPHA(j), j=1,...,N  and  BETA(j), 
  
            j=1,...,N  are the diagonals of the complex Schur form (A,B) 
  
            output by CGEGS.  The  BETA(j) will be non-negative real.   

            Note: the quotients ALPHA(j)/BETA(j) may easily over- or   
            underflow, and BETA(j) may even be zero.  Thus, the user   
            should avoid naively computing the ratio alpha/beta.   
            However, ALPHA will be always less than and usually   
            comparable with norm(A) in magnitude, and BETA always less   
            than and usually comparable with norm(B).   

    VSL     (output) COMPLEX array, dimension (LDVSL,N)   
            If JOBVSL = 'V', VSL will contain the left Schur vectors.   
            (See "Purpose", above.)   
            Not referenced if JOBVSL = 'N'.   

    LDVSL   (input) INTEGER   
            The leading dimension of the matrix VSL. LDVSL >= 1, and   
            if JOBVSL = 'V', LDVSL >= N.   

    VSR     (output) COMPLEX array, dimension (LDVSR,N)   
            If JOBVSR = 'V', VSR will contain the right Schur vectors.   
            (See "Purpose", above.)   
            Not referenced if JOBVSR = 'N'.   

    LDVSR   (input) INTEGER   
            The leading dimension of the matrix VSR. LDVSR >= 1, and   
            if JOBVSR = 'V', LDVSR >= N.   

    WORK    (workspace/output) COMPLEX array, dimension (LWORK)   
            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INTEGER   
            The dimension of the array WORK.  LWORK >= max(1,2*N).   
            For good performance, LWORK must generally be larger.   
            To compute the optimal value of LWORK, call ILAENV to get   
            blocksizes (for CGEQRF, CUNMQR, and CUNGQR.)  Then compute:   
            NB  -- MAX of the blocksizes for CGEQRF, CUNMQR, and CUNGQR; 
  
            the optimal LWORK is N*(NB+1).   

    RWORK   (workspace) REAL array, dimension (3*N)   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value.   
            =1,...,N:   
                  The QZ iteration failed.  (A,B) are not in Schur   
                  form, but ALPHA(j) and BETA(j) should be correct for   
                  j=INFO+1,...,N.   
            > N:  errors that usually indicate LAPACK problems:   
                  =N+1: error return from CGGBAL   
                  =N+2: error return from CGEQRF   
                  =N+3: error return from CUNMQR   
                  =N+4: error return from CUNGQR   
                  =N+5: error return from CGGHRD   
                  =N+6: error return from CHGEQZ (other than failed   
                                                 iteration)   
                  =N+7: error return from CGGBAK (computing VSL)   
                  =N+8: error return from CGGBAK (computing VSR)   
                  =N+9: error return from CLASCL (various places)   

    ===================================================================== 
  


       Decode the input arguments   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static complex c_b1 = {0.f,0.f};
    static complex c_b2 = {1.f,0.f};
    static integer c_n1 = -1;
    static integer c__1 = 1;
    
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, vsl_dim1, vsl_offset, 
	    vsr_dim1, vsr_offset, i__1, i__2, i__3;
    /* Local variables */
    static real anrm, bnrm;
    static integer itau;
    extern logical lsame_(char *, char *);
    static integer ileft, iinfo, icols;
    static logical ilvsl;
    static integer iwork;
    static logical ilvsr;
    static integer irows;
    extern /* Subroutine */ int cggbak_(char *, char *, integer *, integer *, 
	    integer *, real *, real *, integer *, complex *, integer *, 
	    integer *), cggbal_(char *, integer *, complex *, 
	    integer *, complex *, integer *, integer *, integer *, real *, 
	    real *, real *, integer *);
    extern doublereal clange_(char *, integer *, integer *, complex *, 
	    integer *, real *);
    extern /* Subroutine */ int cgghrd_(char *, char *, integer *, integer *, 
	    integer *, complex *, integer *, complex *, integer *, complex *, 
	    integer *, complex *, integer *, integer *), 
	    clascl_(char *, integer *, integer *, real *, real *, integer *, 
	    integer *, complex *, integer *, integer *);
    static logical ilascl, ilbscl;
    extern /* Subroutine */ int cgeqrf_(integer *, integer *, complex *, 
	    integer *, complex *, complex *, integer *, integer *);
    extern doublereal slamch_(char *);
    extern /* Subroutine */ int clacpy_(char *, integer *, integer *, complex 
	    *, integer *, complex *, integer *), claset_(char *, 
	    integer *, integer *, complex *, complex *, complex *, integer *);
    static real safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *);
    static real bignum;
    extern /* Subroutine */ int chgeqz_(char *, char *, char *, integer *, 
	    integer *, integer *, complex *, integer *, complex *, integer *, 
	    complex *, complex *, complex *, integer *, complex *, integer *, 
	    complex *, integer *, real *, integer *);
    static integer ijobvl, iright, ijobvr;
    static real anrmto;
    static integer lwkmin;
    static real bnrmto;
    extern /* Subroutine */ int cungqr_(integer *, integer *, integer *, 
	    complex *, integer *, complex *, complex *, integer *, integer *),
	     cunmqr_(char *, char *, integer *, integer *, integer *, complex 
	    *, integer *, complex *, complex *, integer *, complex *, integer 
	    *, integer *);
    static real smlnum;
    static integer irwork, lwkopt, ihi, ilo;
    static real eps;



#define ALPHA(I) alpha[(I)-1]
#define BETA(I) beta[(I)-1]
#define WORK(I) work[(I)-1]
#define RWORK(I) rwork[(I)-1]

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
    i__1 = *n << 1;
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
	*info = -11;
    } else if (*ldvsr < 1 || ilvsr && *ldvsr < *n) {
	*info = -13;
    } else if (*lwork < lwkmin) {
	*info = -15;
    }

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("CGEGS ", &i__1);
	return 0;
    }

/*     Quick return if possible */

    WORK(1).r = (real) lwkopt, WORK(1).i = 0.f;
    if (*n == 0) {
	return 0;
    }

/*     Get machine constants */

    eps = slamch_("E") * slamch_("B");
    safmin = slamch_("S");
    smlnum = *n * safmin / eps;
    bignum = 1.f / smlnum;

/*     Scale A if max element outside range [SMLNUM,BIGNUM] */

    anrm = clange_("M", n, n, &A(1,1), lda, &RWORK(1));
    ilascl = FALSE_;
    if (anrm > 0.f && anrm < smlnum) {
	anrmto = smlnum;
	ilascl = TRUE_;
    } else if (anrm > bignum) {
	anrmto = bignum;
	ilascl = TRUE_;
    }

    if (ilascl) {
	clascl_("G", &c_n1, &c_n1, &anrm, &anrmto, n, n, &A(1,1), lda, &
		iinfo);
	if (iinfo != 0) {
	    *info = *n + 9;
	    return 0;
	}
    }

/*     Scale B if max element outside range [SMLNUM,BIGNUM] */

    bnrm = clange_("M", n, n, &B(1,1), ldb, &RWORK(1));
    ilbscl = FALSE_;
    if (bnrm > 0.f && bnrm < smlnum) {
	bnrmto = smlnum;
	ilbscl = TRUE_;
    } else if (bnrm > bignum) {
	bnrmto = bignum;
	ilbscl = TRUE_;
    }

    if (ilbscl) {
	clascl_("G", &c_n1, &c_n1, &bnrm, &bnrmto, n, n, &B(1,1), ldb, &
		iinfo);
	if (iinfo != 0) {
	    *info = *n + 9;
	    return 0;
	}
    }

/*     Permute the matrix to make it more nearly triangular */

    ileft = 1;
    iright = *n + 1;
    irwork = iright + *n;
    iwork = 1;
    cggbal_("P", n, &A(1,1), lda, &B(1,1), ldb, &ilo, &ihi, &RWORK(
	    ileft), &RWORK(iright), &RWORK(irwork), &iinfo);
    if (iinfo != 0) {
	*info = *n + 1;
	goto L10;
    }

/*     Reduce B to triangular form, and initialize VSL and/or VSR */

    irows = ihi + 1 - ilo;
    icols = *n + 1 - ilo;
    itau = iwork;
    iwork = itau + irows;
    i__1 = *lwork + 1 - iwork;
    cgeqrf_(&irows, &icols, &B(ilo,ilo), ldb, &WORK(itau), &WORK(
	    iwork), &i__1, &iinfo);
    if (iinfo >= 0) {
/* Computing MAX */
	i__3 = iwork;
	i__1 = lwkopt, i__2 = (integer) WORK(iwork).r + iwork - 1;
	lwkopt = max(i__1,i__2);
    }
    if (iinfo != 0) {
	*info = *n + 2;
	goto L10;
    }

    i__1 = *lwork + 1 - iwork;
    cunmqr_("L", "C", &irows, &icols, &irows, &B(ilo,ilo), ldb, &
	    WORK(itau), &A(ilo,ilo), lda, &WORK(iwork), &i__1, &
	    iinfo);
    if (iinfo >= 0) {
/* Computing MAX */
	i__3 = iwork;
	i__1 = lwkopt, i__2 = (integer) WORK(iwork).r + iwork - 1;
	lwkopt = max(i__1,i__2);
    }
    if (iinfo != 0) {
	*info = *n + 3;
	goto L10;
    }

    if (ilvsl) {
	claset_("Full", n, n, &c_b1, &c_b2, &VSL(1,1), ldvsl);
	i__1 = irows - 1;
	i__2 = irows - 1;
	clacpy_("L", &i__1, &i__2, &B(ilo+1,ilo), ldb, &VSL(ilo+1,ilo), ldvsl);
	i__1 = *lwork + 1 - iwork;
	cungqr_(&irows, &irows, &irows, &VSL(ilo,ilo), ldvsl, &
		WORK(itau), &WORK(iwork), &i__1, &iinfo);
	if (iinfo >= 0) {
/* Computing MAX */
	    i__3 = iwork;
	    i__1 = lwkopt, i__2 = (integer) WORK(iwork).r + iwork - 1;
	    lwkopt = max(i__1,i__2);
	}
	if (iinfo != 0) {
	    *info = *n + 4;
	    goto L10;
	}
    }

    if (ilvsr) {
	claset_("Full", n, n, &c_b1, &c_b2, &VSR(1,1), ldvsr);
    }

/*     Reduce to generalized Hessenberg form */

    cgghrd_(jobvsl, jobvsr, n, &ilo, &ihi, &A(1,1), lda, &B(1,1), 
	    ldb, &VSL(1,1), ldvsl, &VSR(1,1), ldvsr, &iinfo);
    if (iinfo != 0) {
	*info = *n + 5;
	goto L10;
    }

/*     Perform QZ algorithm, computing Schur vectors if desired */

    iwork = itau;
    i__1 = *lwork + 1 - iwork;
    chgeqz_("S", jobvsl, jobvsr, n, &ilo, &ihi, &A(1,1), lda, &B(1,1), ldb, &ALPHA(1), &BETA(1), &VSL(1,1), ldvsl, &
	    VSR(1,1), ldvsr, &WORK(iwork), &i__1, &RWORK(irwork), &
	    iinfo);
    if (iinfo >= 0) {
/* Computing MAX */
	i__3 = iwork;
	i__1 = lwkopt, i__2 = (integer) WORK(iwork).r + iwork - 1;
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
	cggbak_("P", "L", n, &ilo, &ihi, &RWORK(ileft), &RWORK(iright), n, &
		VSL(1,1), ldvsl, &iinfo);
	if (iinfo != 0) {
	    *info = *n + 7;
	    goto L10;
	}
    }
    if (ilvsr) {
	cggbak_("P", "R", n, &ilo, &ihi, &RWORK(ileft), &RWORK(iright), n, &
		VSR(1,1), ldvsr, &iinfo);
	if (iinfo != 0) {
	    *info = *n + 8;
	    goto L10;
	}
    }

/*     Undo scaling */

    if (ilascl) {
	clascl_("U", &c_n1, &c_n1, &anrmto, &anrm, n, n, &A(1,1), lda, &
		iinfo);
	if (iinfo != 0) {
	    *info = *n + 9;
	    return 0;
	}
	clascl_("G", &c_n1, &c_n1, &anrmto, &anrm, n, &c__1, &ALPHA(1), n, &
		iinfo);
	if (iinfo != 0) {
	    *info = *n + 9;
	    return 0;
	}
    }

    if (ilbscl) {
	clascl_("U", &c_n1, &c_n1, &bnrmto, &bnrm, n, n, &B(1,1), ldb, &
		iinfo);
	if (iinfo != 0) {
	    *info = *n + 9;
	    return 0;
	}
	clascl_("G", &c_n1, &c_n1, &bnrmto, &bnrm, n, &c__1, &BETA(1), n, &
		iinfo);
	if (iinfo != 0) {
	    *info = *n + 9;
	    return 0;
	}
    }

L10:
    WORK(1).r = (real) lwkopt, WORK(1).i = 0.f;

    return 0;

/*     End of CGEGS */

} /* cgegs_ */

