#include "f2c.h"

/* Subroutine */ int sgegv_(char *jobvl, char *jobvr, integer *n, real *a, 
	integer *lda, real *b, integer *ldb, real *alphar, real *alphai, real 
	*beta, real *vl, integer *ldvl, real *vr, integer *ldvr, real *work, 
	integer *lwork, integer *info)
{
/*  -- LAPACK driver routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    SGEGV computes for a pair of n-by-n real nonsymmetric matrices A and 
  
    B, the generalized eigenvalues (alphar +/- alphai*i, beta), and   
    optionally, the left and/or right generalized eigenvectors (VL and   
    VR).   

    A generalized eigenvalue for a pair of matrices (A,B) is, roughly   
    speaking, a scalar w or a ratio  alpha/beta = w, such that  A - w*B   
    is singular.  It is usually represented as the pair (alpha,beta),   
    as there is a reasonable interpretation for beta=0, and even for   
    both being zero.  A good beginning reference is the book, "VISMatrix   
    Computations", by G. Golub & C. van Loan (Johns Hopkins U. Press)   

    A right generalized eigenvector corresponding to a generalized   
    eigenvalue  w  for a pair of matrices (A,B) is a vector  r  such   
    that  (A - w B) r = 0 .  A left generalized eigenvector is a vector   
    l such that l**H * (A - w B) = 0, where l**H is the   
    conjugate-transpose of l.   

    Note: this routine performs "full balancing" on A and B -- see   
    "Further Details", below.   

    Arguments   
    =========   

    JOBVL   (input) CHARACTER*1   
            = 'N':  do not compute the left generalized eigenvectors;   
            = 'V':  compute the left generalized eigenvectors.   

    JOBVR   (input) CHARACTER*1   
            = 'N':  do not compute the right generalized eigenvectors;   
            = 'V':  compute the right generalized eigenvectors.   

    N       (input) INTEGER   
            The order of the matrices A, B, VL, and VR.  N >= 0.   

    A       (input/output) REAL array, dimension (LDA, N)   
            On entry, the first of the pair of matrices whose   
            generalized eigenvalues and (optionally) generalized   
            eigenvectors are to be computed.   
            On exit, the contents will have been destroyed.  (For a   
            description of the contents of A on exit, see "Further   
            Details", below.)   

    LDA     (input) INTEGER   
            The leading dimension of A.  LDA >= max(1,N).   

    B       (input/output) REAL array, dimension (LDB, N)   
            On entry, the second of the pair of matrices whose   
            generalized eigenvalues and (optionally) generalized   
            eigenvectors are to be computed.   
            On exit, the contents will have been destroyed.  (For a   
            description of the contents of B on exit, see "Further   
            Details", below.)   

    LDB     (input) INTEGER   
            The leading dimension of B.  LDB >= max(1,N).   

    ALPHAR  (output) REAL array, dimension (N)   
    ALPHAI  (output) REAL array, dimension (N)   
    BETA    (output) REAL array, dimension (N)   
            On exit, (ALPHAR(j) + ALPHAI(j)*i)/BETA(j), j=1,...,N, will   
            be the generalized eigenvalues.  If ALPHAI(j) is zero, then   
            the j-th eigenvalue is real; if positive, then the j-th and   
            (j+1)-st eigenvalues are a complex conjugate pair, with   
            ALPHAI(j+1) negative.   

            Note: the quotients ALPHAR(j)/BETA(j) and ALPHAI(j)/BETA(j)   
            may easily over- or underflow, and BETA(j) may even be zero. 
  
            Thus, the user should avoid naively computing the ratio   
            alpha/beta.  However, ALPHAR and ALPHAI will be always less   
            than and usually comparable with norm(A) in magnitude, and   
            BETA always less than and usually comparable with norm(B).   

    VL      (output) REAL array, dimension (LDVL,N)   
            If JOBVL = 'V', the left generalized eigenvectors.  (See   
            "Purpose", above.)  Real eigenvectors take one column,   
            complex take two columns, the first for the real part and   
            the second for the imaginary part.  Complex eigenvectors   
            correspond to an eigenvalue with positive imaginary part.   
            Each eigenvector will be scaled so the largest component   
            will have abs(real part) + abs(imag. part) = 1, *except*   
            that for eigenvalues with alpha=beta=0, a zero vector will   
            be returned as the corresponding eigenvector.   
            Not referenced if JOBVL = 'N'.   

    LDVL    (input) INTEGER   
            The leading dimension of the matrix VL. LDVL >= 1, and   
            if JOBVL = 'V', LDVL >= N.   

    VR      (output) REAL array, dimension (LDVR,N)   
            If JOBVL = 'V', the right generalized eigenvectors.  (See   
            "Purpose", above.)  Real eigenvectors take one column,   
            complex take two columns, the first for the real part and   
            the second for the imaginary part.  Complex eigenvectors   
            correspond to an eigenvalue with positive imaginary part.   
            Each eigenvector will be scaled so the largest component   
            will have abs(real part) + abs(imag. part) = 1, *except*   
            that for eigenvalues with alpha=beta=0, a zero vector will   
            be returned as the corresponding eigenvector.   
            Not referenced if JOBVR = 'N'.   

    LDVR    (input) INTEGER   
            The leading dimension of the matrix VR. LDVR >= 1, and   
            if JOBVR = 'V', LDVR >= N.   

    WORK    (workspace/output) REAL array, dimension (LWORK)   
            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INTEGER   
            The dimension of the array WORK.  LWORK >= max(1,8*N).   
            For good performance, LWORK must generally be larger.   
            To compute the optimal value of LWORK, call ILAENV to get   
            blocksizes (for SGEQRF, SORMQR, and SORGQR.)  Then compute:   
            NB  -- MAX of the blocksizes for SGEQRF, SORMQR, and SORGQR; 
  
            The optimal LWORK is:   
                2*N + MAX( 6*N, N*(NB+1) ).   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value.   
            = 1,...,N:   
                  The QZ iteration failed.  No eigenvectors have been   
                  calculated, but ALPHAR(j), ALPHAI(j), and BETA(j)   
                  should be correct for j=INFO+1,...,N.   
            > N:  errors that usually indicate LAPACK problems:   
                  =N+1: error return from SGGBAL   
                  =N+2: error return from SGEQRF   
                  =N+3: error return from SORMQR   
                  =N+4: error return from SORGQR   
                  =N+5: error return from SGGHRD   
                  =N+6: error return from SHGEQZ (other than failed   
                                                  iteration)   
                  =N+7: error return from STGEVC   
                  =N+8: error return from SGGBAK (computing VL)   
                  =N+9: error return from SGGBAK (computing VR)   
                  =N+10: error return from SLASCL (various calls)   

    Further Details   
    ===============   

    Balancing   
    ---------   

    This driver calls SGGBAL to both permute and scale rows and columns   
    of A and B.  The permutations PL and PR are chosen so that PL*A*PR   
    and PL*B*R will be upper triangular except for the diagonal blocks   
    A(i:j,i:j) and B(i:j,i:j), with i and j as close together as   
    possible.  The diagonal scaling matrices DL and DR are chosen so   
    that the pair  DL*PL*A*PR*DR, DL*PL*B*PR*DR have elements close to   
    one (except for the elements that start out zero.)   

    After the eigenvalues and eigenvectors of the balanced matrices   
    have been computed, SGGBAK transforms the eigenvectors back to what   
    they would have been (in perfect arithmetic) if they had not been   
    balanced.   

    Contents of A and B on Exit   
    -------- -- - --- - -- ----   

    If any eigenvectors are computed (either JOBVL='V' or JOBVR='V' or   
    both), then on exit the arrays A and B will contain the real Schur   
    form[*] of the "balanced" versions of A and B.  If no eigenvectors   
    are computed, then only the diagonal blocks will be correct.   

    [*] See SHGEQZ, SGEGS, or read the book "VISMatrix Computations",   
        by Golub & van Loan, pub. by Johns Hopkins U. Press.   

    ===================================================================== 
  


       Decode the input arguments   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static integer c_n1 = -1;
    static real c_b14 = 1.f;
    static real c_b25 = 0.f;
    static integer c__1 = 1;
    
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, vl_dim1, vl_offset, vr_dim1, 
	    vr_offset, i__1, i__2;
    real r__1, r__2, r__3, r__4;
    /* Local variables */
    static real absb, anrm, bnrm;
    static integer itau;
    static real temp;
    static logical ilvl, ilvr;
    static real anrm1, anrm2, bnrm1, bnrm2, absai, scale, absar, sbeta;
    extern logical lsame_(char *, char *);
    static integer ileft, iinfo, icols, iwork, irows, jc, in, jr;
    static real salfai;
    extern /* Subroutine */ int sggbak_(char *, char *, integer *, integer *, 
	    integer *, real *, real *, integer *, real *, integer *, integer *
	    ), sggbal_(char *, integer *, real *, integer *, 
	    real *, integer *, integer *, integer *, real *, real *, real *, 
	    integer *);
    static real salfar;
    extern doublereal slamch_(char *), slange_(char *, integer *, 
	    integer *, real *, integer *, real *);
    static real safmin;
    extern /* Subroutine */ int sgghrd_(char *, char *, integer *, integer *, 
	    integer *, real *, integer *, real *, integer *, real *, integer *
	    , real *, integer *, integer *);
    static real safmax;
    static char chtemp[1];
    static logical ldumma[1];
    extern /* Subroutine */ int slascl_(char *, integer *, integer *, real *, 
	    real *, integer *, integer *, real *, integer *, integer *), xerbla_(char *, integer *);
    static integer ijobvl, iright;
    static logical ilimit;
    extern /* Subroutine */ int sgeqrf_(integer *, integer *, real *, integer 
	    *, real *, real *, integer *, integer *);
    static integer ijobvr;
    extern /* Subroutine */ int slacpy_(char *, integer *, integer *, real *, 
	    integer *, real *, integer *), slaset_(char *, integer *, 
	    integer *, real *, real *, real *, integer *), stgevc_(
	    char *, char *, logical *, integer *, real *, integer *, real *, 
	    integer *, real *, integer *, real *, integer *, integer *, 
	    integer *, real *, integer *);
    static real onepls;
    static integer lwkmin;
    extern /* Subroutine */ int shgeqz_(char *, char *, char *, integer *, 
	    integer *, integer *, real *, integer *, real *, integer *, real *
	    , real *, real *, real *, integer *, real *, integer *, real *, 
	    integer *, integer *), sorgqr_(integer *, 
	    integer *, integer *, real *, integer *, real *, real *, integer *
	    , integer *);
    static integer lwkopt;
    extern /* Subroutine */ int sormqr_(char *, char *, integer *, integer *, 
	    integer *, real *, integer *, real *, real *, integer *, real *, 
	    integer *, integer *);
    static integer ihi, ilo;
    static real eps;
    static logical ilv;



#define ALPHAR(I) alphar[(I)-1]
#define ALPHAI(I) alphai[(I)-1]
#define BETA(I) beta[(I)-1]
#define WORK(I) work[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]
#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]
#define VL(I,J) vl[(I)-1 + ((J)-1)* ( *ldvl)]
#define VR(I,J) vr[(I)-1 + ((J)-1)* ( *ldvr)]

    if (lsame_(jobvl, "N")) {
	ijobvl = 1;
	ilvl = FALSE_;
    } else if (lsame_(jobvl, "V")) {
	ijobvl = 2;
	ilvl = TRUE_;
    } else {
	ijobvl = -1;
	ilvl = FALSE_;
    }

    if (lsame_(jobvr, "N")) {
	ijobvr = 1;
	ilvr = FALSE_;
    } else if (lsame_(jobvr, "V")) {
	ijobvr = 2;
	ilvr = TRUE_;
    } else {
	ijobvr = -1;
	ilvr = FALSE_;
    }
    ilv = ilvl || ilvr;

/*     Test the input arguments   

   Computing MAX */
    i__1 = *n << 3;
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
    } else if (*ldvl < 1 || ilvl && *ldvl < *n) {
	*info = -12;
    } else if (*ldvr < 1 || ilvr && *ldvr < *n) {
	*info = -14;
    } else if (*lwork < lwkmin) {
	*info = -16;
    }

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("SGEGV ", &i__1);
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
    safmin += safmin;
    safmax = 1.f / safmin;
    onepls = eps * 4 + 1.f;

/*     Scale A */

    anrm = slange_("M", n, n, &A(1,1), lda, &WORK(1));
    anrm1 = anrm;
    anrm2 = 1.f;
    if (anrm < 1.f) {
	if (safmax * anrm < 1.f) {
	    anrm1 = safmin;
	    anrm2 = safmax * anrm;
	}
    }

    if (anrm > 0.f) {
	slascl_("G", &c_n1, &c_n1, &anrm, &c_b14, n, n, &A(1,1), lda, &
		iinfo);
	if (iinfo != 0) {
	    *info = *n + 10;
	    return 0;
	}
    }

/*     Scale B */

    bnrm = slange_("M", n, n, &B(1,1), ldb, &WORK(1));
    bnrm1 = bnrm;
    bnrm2 = 1.f;
    if (bnrm < 1.f) {
	if (safmax * bnrm < 1.f) {
	    bnrm1 = safmin;
	    bnrm2 = safmax * bnrm;
	}
    }

    if (bnrm > 0.f) {
	slascl_("G", &c_n1, &c_n1, &bnrm, &c_b14, n, n, &B(1,1), ldb, &
		iinfo);
	if (iinfo != 0) {
	    *info = *n + 10;
	    return 0;
	}
    }

/*     Permute the matrix to make it more nearly triangular   
       Workspace layout:  (8*N words -- "work" requires 6*N words)   
          left_permutation, right_permutation, work... */

    ileft = 1;
    iright = *n + 1;
    iwork = iright + *n;
    sggbal_("B", n, &A(1,1), lda, &B(1,1), ldb, &ilo, &ihi, &WORK(
	    ileft), &WORK(iright), &WORK(iwork), &iinfo);
    if (iinfo != 0) {
	*info = *n + 1;
	goto L120;
    }

/*     Reduce B to triangular form, and initialize VL and/or VR   
       Workspace layout:  ("work..." must have at least N words)   
          left_permutation, right_permutation, tau, work... */

    irows = ihi + 1 - ilo;
    if (ilv) {
	icols = *n + 1 - ilo;
    } else {
	icols = irows;
    }
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
	goto L120;
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
	goto L120;
    }

    if (ilvl) {
	slaset_("Full", n, n, &c_b25, &c_b14, &VL(1,1), ldvl);
	i__1 = irows - 1;
	i__2 = irows - 1;
	slacpy_("L", &i__1, &i__2, &B(ilo+1,ilo), ldb, &VL(ilo+1,ilo), ldvl);
	i__1 = *lwork + 1 - iwork;
	sorgqr_(&irows, &irows, &irows, &VL(ilo,ilo), ldvl, &WORK(
		itau), &WORK(iwork), &i__1, &iinfo);
	if (iinfo >= 0) {
/* Computing MAX */
	    i__1 = lwkopt, i__2 = (integer) WORK(iwork) + iwork - 1;
	    lwkopt = max(i__1,i__2);
	}
	if (iinfo != 0) {
	    *info = *n + 4;
	    goto L120;
	}
    }

    if (ilvr) {
	slaset_("Full", n, n, &c_b25, &c_b14, &VR(1,1), ldvr);
    }

/*     Reduce to generalized Hessenberg form */

    if (ilv) {

/*        Eigenvectors requested -- work on whole matrix. */

	sgghrd_(jobvl, jobvr, n, &ilo, &ihi, &A(1,1), lda, &B(1,1), 
		ldb, &VL(1,1), ldvl, &VR(1,1), ldvr, &iinfo);
    } else {
	sgghrd_("N", "N", &irows, &c__1, &irows, &A(ilo,ilo), lda, 
		&B(ilo,ilo), ldb, &VL(1,1), ldvl, &VR(1,1), ldvr, &iinfo);
    }
    if (iinfo != 0) {
	*info = *n + 5;
	goto L120;
    }

/*     Perform QZ algorithm   
       Workspace layout:  ("work..." must have at least 1 word)   
          left_permutation, right_permutation, work... */

    iwork = itau;
    if (ilv) {
	*(unsigned char *)chtemp = 'S';
    } else {
	*(unsigned char *)chtemp = 'E';
    }
    i__1 = *lwork + 1 - iwork;
    shgeqz_(chtemp, jobvl, jobvr, n, &ilo, &ihi, &A(1,1), lda, &B(1,1), ldb, &ALPHAR(1), &ALPHAI(1), &BETA(1), &VL(1,1), 
	    ldvl, &VR(1,1), ldvr, &WORK(iwork), &i__1, &iinfo);
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
	goto L120;
    }

    if (ilv) {

/*        Compute Eigenvectors  (STGEVC requires 6*N words of workspac
e) */

	if (ilvl) {
	    if (ilvr) {
		*(unsigned char *)chtemp = 'B';
	    } else {
		*(unsigned char *)chtemp = 'L';
	    }
	} else {
	    *(unsigned char *)chtemp = 'R';
	}

	stgevc_(chtemp, "B", ldumma, n, &A(1,1), lda, &B(1,1), ldb, 
		&VL(1,1), ldvl, &VR(1,1), ldvr, n, &in, &WORK(
		iwork), &iinfo);
	if (iinfo != 0) {
	    *info = *n + 7;
	    goto L120;
	}

/*        Undo balancing on VL and VR, rescale */

	if (ilvl) {
	    sggbak_("B", "L", n, &ilo, &ihi, &WORK(ileft), &WORK(iright), n, &
		    VL(1,1), ldvl, &iinfo);
	    if (iinfo != 0) {
		*info = *n + 8;
		goto L120;
	    }
	    i__1 = *n;
	    for (jc = 1; jc <= *n; ++jc) {
		if (ALPHAI(jc) < 0.f) {
		    goto L50;
		}
		temp = 0.f;
		if (ALPHAI(jc) == 0.f) {
		    i__2 = *n;
		    for (jr = 1; jr <= *n; ++jr) {
/* Computing MAX */
			r__2 = temp, r__3 = (r__1 = VL(jr,jc), 
				dabs(r__1));
			temp = dmax(r__2,r__3);
/* L10: */
		    }
		} else {
		    i__2 = *n;
		    for (jr = 1; jr <= *n; ++jr) {
/* Computing MAX */
			r__3 = temp, r__4 = (r__1 = VL(jr,jc), 
				dabs(r__1)) + (r__2 = VL(jr,jc+1), dabs(r__2));
			temp = dmax(r__3,r__4);
/* L20: */
		    }
		}
		if (temp < safmin) {
		    goto L50;
		}
		temp = 1.f / temp;
		if (ALPHAI(jc) == 0.f) {
		    i__2 = *n;
		    for (jr = 1; jr <= *n; ++jr) {
			VL(jr,jc) *= temp;
/* L30: */
		    }
		} else {
		    i__2 = *n;
		    for (jr = 1; jr <= *n; ++jr) {
			VL(jr,jc) *= temp;
			VL(jr,jc+1) *= temp;
/* L40: */
		    }
		}
L50:
		;
	    }
	}
	if (ilvr) {
	    sggbak_("B", "R", n, &ilo, &ihi, &WORK(ileft), &WORK(iright), n, &
		    VR(1,1), ldvr, &iinfo);
	    if (iinfo != 0) {
		*info = *n + 9;
		goto L120;
	    }
	    i__1 = *n;
	    for (jc = 1; jc <= *n; ++jc) {
		if (ALPHAI(jc) < 0.f) {
		    goto L100;
		}
		temp = 0.f;
		if (ALPHAI(jc) == 0.f) {
		    i__2 = *n;
		    for (jr = 1; jr <= *n; ++jr) {
/* Computing MAX */
			r__2 = temp, r__3 = (r__1 = VR(jr,jc), 
				dabs(r__1));
			temp = dmax(r__2,r__3);
/* L60: */
		    }
		} else {
		    i__2 = *n;
		    for (jr = 1; jr <= *n; ++jr) {
/* Computing MAX */
			r__3 = temp, r__4 = (r__1 = VR(jr,jc), 
				dabs(r__1)) + (r__2 = VR(jr,jc+1), dabs(r__2));
			temp = dmax(r__3,r__4);
/* L70: */
		    }
		}
		if (temp < safmin) {
		    goto L100;
		}
		temp = 1.f / temp;
		if (ALPHAI(jc) == 0.f) {
		    i__2 = *n;
		    for (jr = 1; jr <= *n; ++jr) {
			VR(jr,jc) *= temp;
/* L80: */
		    }
		} else {
		    i__2 = *n;
		    for (jr = 1; jr <= *n; ++jr) {
			VR(jr,jc) *= temp;
			VR(jr,jc+1) *= temp;
/* L90: */
		    }
		}
L100:
		;
	    }
	}

/*        End of eigenvector calculation */

    }

/*     Undo scaling in alpha, beta   

       Note: this does not give the alpha and beta for the unscaled   
       problem.   

       Un-scaling is limited to avoid underflow in alpha and beta   
       if they are significant. */

    i__1 = *n;
    for (jc = 1; jc <= *n; ++jc) {
	absar = (r__1 = ALPHAR(jc), dabs(r__1));
	absai = (r__1 = ALPHAI(jc), dabs(r__1));
	absb = (r__1 = BETA(jc), dabs(r__1));
	salfar = anrm * ALPHAR(jc);
	salfai = anrm * ALPHAI(jc);
	sbeta = bnrm * BETA(jc);
	ilimit = FALSE_;
	scale = 1.f;

/*        Check for significant underflow in ALPHAI   

   Computing MAX */
	r__1 = safmin, r__2 = eps * absar, r__1 = max(r__1,r__2), r__2 = eps *
		 absb;
	if (dabs(salfai) < safmin && absai >= dmax(r__1,r__2)) {
	    ilimit = TRUE_;
/* Computing MAX */
	    r__1 = onepls * safmin, r__2 = anrm2 * absai;
	    scale = onepls * safmin / anrm1 / dmax(r__1,r__2);

	} else if (salfai == 0.f) {

/*           If insignificant underflow in ALPHAI, then make the 
  
             conjugate eigenvalue real. */

	    if (ALPHAI(jc) < 0.f && jc > 1) {
		ALPHAI(jc - 1) = 0.f;
	    } else if (ALPHAI(jc) > 0.f && jc < *n) {
		ALPHAI(jc + 1) = 0.f;
	    }
	}

/*        Check for significant underflow in ALPHAR   

   Computing MAX */
	r__1 = safmin, r__2 = eps * absai, r__1 = max(r__1,r__2), r__2 = eps *
		 absb;
	if (dabs(salfar) < safmin && absar >= dmax(r__1,r__2)) {
	    ilimit = TRUE_;
/* Computing MAX   
   Computing MAX */
	    r__3 = onepls * safmin, r__4 = anrm2 * absar;
	    r__1 = scale, r__2 = onepls * safmin / anrm1 / dmax(r__3,r__4);
	    scale = dmax(r__1,r__2);
	}

/*        Check for significant underflow in BETA   

   Computing MAX */
	r__1 = safmin, r__2 = eps * absar, r__1 = max(r__1,r__2), r__2 = eps *
		 absai;
	if (dabs(sbeta) < safmin && absb >= dmax(r__1,r__2)) {
	    ilimit = TRUE_;
/* Computing MAX   
   Computing MAX */
	    r__3 = onepls * safmin, r__4 = bnrm2 * absb;
	    r__1 = scale, r__2 = onepls * safmin / bnrm1 / dmax(r__3,r__4);
	    scale = dmax(r__1,r__2);
	}

/*        Check for possible overflow when limiting scaling */

	if (ilimit) {
/* Computing MAX */
	    r__1 = dabs(salfar), r__2 = dabs(salfai), r__1 = max(r__1,r__2), 
		    r__2 = dabs(sbeta);
	    temp = scale * safmin * dmax(r__1,r__2);
	    if (temp > 1.f) {
		scale /= temp;
	    }
	    if (scale < 1.f) {
		ilimit = FALSE_;
	    }
	}

/*        Recompute un-scaled ALPHAR, ALPHAI, BETA if necessary. */

	if (ilimit) {
	    salfar = scale * ALPHAR(jc) * anrm;
	    salfai = scale * ALPHAI(jc) * anrm;
	    sbeta = scale * BETA(jc) * bnrm;
	}
	ALPHAR(jc) = salfar;
	ALPHAI(jc) = salfai;
	BETA(jc) = sbeta;
/* L110: */
    }

L120:
    WORK(1) = (real) lwkopt;

    return 0;

/*     End of SGEGV */

} /* sgegv_ */

