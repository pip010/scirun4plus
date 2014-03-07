#include "f2c.h"

/* Subroutine */ int zgegv_(char *jobvl, char *jobvr, integer *n, 
	doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, 
	doublecomplex *alpha, doublecomplex *beta, doublecomplex *vl, integer 
	*ldvl, doublecomplex *vr, integer *ldvr, doublecomplex *work, integer 
	*lwork, doublereal *rwork, integer *info)
{
/*  -- LAPACK driver routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    ZGEGV computes for a pair of N-by-N complex nonsymmetric matrices A   
    and B, the generalized eigenvalues (alpha, beta), and optionally,   
    the left and/or right generalized eigenvectors (VL and VR).   

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

    A       (input/output) COMPLEX*16 array, dimension (LDA, N)   
            On entry, the first of the pair of matrices whose   
            generalized eigenvalues and (optionally) generalized   
            eigenvectors are to be computed.   
            On exit, the contents will have been destroyed.  (For a   
            description of the contents of A on exit, see "Further   
            Details", below.)   

    LDA     (input) INTEGER   
            The leading dimension of A.  LDA >= max(1,N).   

    B       (input/output) COMPLEX*16 array, dimension (LDB, N)   
            On entry, the second of the pair of matrices whose   
            generalized eigenvalues and (optionally) generalized   
            eigenvectors are to be computed.   
            On exit, the contents will have been destroyed.  (For a   
            description of the contents of B on exit, see "Further   
            Details", below.)   

    LDB     (input) INTEGER   
            The leading dimension of B.  LDB >= max(1,N).   

    ALPHA   (output) COMPLEX*16 array, dimension (N)   
    BETA    (output) COMPLEX*16 array, dimension (N)   
            On exit, ALPHA(j)/BETA(j), j=1,...,N, will be the   
            generalized eigenvalues.   

            Note: the quotients ALPHA(j)/BETA(j) may easily over- or   
            underflow, and BETA(j) may even be zero.  Thus, the user   
            should avoid naively computing the ratio alpha/beta.   
            However, ALPHA will be always less than and usually   
            comparable with norm(A) in magnitude, and BETA always less   
            than and usually comparable with norm(B).   

    VL      (output) COMPLEX*16 array, dimension (LDVL,N)   
            If JOBVL = 'V', the left generalized eigenvectors.  (See   
            "Purpose", above.)   
            Each eigenvector will be scaled so the largest component   
            will have abs(real part) + abs(imag. part) = 1, *except*   
            that for eigenvalues with alpha=beta=0, a zero vector will   
            be returned as the corresponding eigenvector.   
            Not referenced if JOBVL = 'N'.   

    LDVL    (input) INTEGER   
            The leading dimension of the matrix VL. LDVL >= 1, and   
            if JOBVL = 'V', LDVL >= N.   

    VR      (output) COMPLEX*16 array, dimension (LDVR,N)   
            If JOBVL = 'V', the right generalized eigenvectors.  (See   
            "Purpose", above.)   
            Each eigenvector will be scaled so the largest component   
            will have abs(real part) + abs(imag. part) = 1, *except*   
            that for eigenvalues with alpha=beta=0, a zero vector will   
            be returned as the corresponding eigenvector.   
            Not referenced if JOBVR = 'N'.   

    LDVR    (input) INTEGER   
            The leading dimension of the matrix VR. LDVR >= 1, and   
            if JOBVR = 'V', LDVR >= N.   

    WORK    (workspace/output) COMPLEX*16 array, dimension (LWORK)   
            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INTEGER   
            The dimension of the array WORK.  LWORK >= max(1,2*N).   
            For good performance, LWORK must generally be larger.   
            To compute the optimal value of LWORK, call ILAENV to get   
            blocksizes (for ZGEQRF, ZUNMQR, and CUNGQR.)  Then compute:   
            NB  -- MAX of the blocksizes for ZGEQRF, ZUNMQR, and CUNGQR; 
  
            The optimal LWORK is  MAX( 2*N, N*(NB+1) ).   

    RWORK   (workspace/output) DOUBLE PRECISION array, dimension (8*N)   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value.   
            =1,...,N:   
                  The QZ iteration failed.  No eigenvectors have been   
                  calculated, but ALPHA(j) and BETA(j) should be   
                  correct for j=INFO+1,...,N.   
            > N:  errors that usually indicate LAPACK problems:   
                  =N+1: error return from ZGGBAL   
                  =N+2: error return from ZGEQRF   
                  =N+3: error return from ZUNMQR   
                  =N+4: error return from ZUNGQR   
                  =N+5: error return from ZGGHRD   
                  =N+6: error return from ZHGEQZ (other than failed   
                                                 iteration)   
                  =N+7: error return from ZTGEVC   
                  =N+8: error return from ZGGBAK (computing VL)   
                  =N+9: error return from ZGGBAK (computing VR)   
                  =N+10: error return from ZLASCL (various calls)   

    Further Details   
    ===============   

    Balancing   
    ---------   

    This driver calls ZGGBAL to both permute and scale rows and columns   
    of A and B.  The permutations PL and PR are chosen so that PL*A*PR   
    and PL*B*R will be upper triangular except for the diagonal blocks   
    A(i:j,i:j) and B(i:j,i:j), with i and j as close together as   
    possible.  The diagonal scaling matrices DL and DR are chosen so   
    that the pair  DL*PL*A*PR*DR, DL*PL*B*PR*DR have elements close to   
    one (except for the elements that start out zero.)   

    After the eigenvalues and eigenvectors of the balanced matrices   
    have been computed, ZGGBAK transforms the eigenvectors back to what   
    they would have been (in perfect arithmetic) if they had not been   
    balanced.   

    Contents of A and B on Exit   
    -------- -- - --- - -- ----   

    If any eigenvectors are computed (either JOBVL='V' or JOBVR='V' or   
    both), then on exit the arrays A and B will contain the complex Schur 
  
    form[*] of the "balanced" versions of A and B.  If no eigenvectors   
    are computed, then only the diagonal blocks will be correct.   

    [*] In other words, upper triangular form.   

    ===================================================================== 
  


       Decode the input arguments   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static doublecomplex c_b1 = {0.,0.};
    static doublecomplex c_b2 = {1.,0.};
    static integer c_n1 = -1;
    static doublereal c_b16 = 1.;
    static integer c__1 = 1;
    
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, vl_dim1, vl_offset, vr_dim1, 
	    vr_offset, i__1, i__2, i__3, i__4;
    doublereal d__1, d__2, d__3, d__4;
    doublecomplex z__1, z__2;
    /* Builtin functions */
    double d_imag(doublecomplex *);
    /* Local variables */
    static doublereal absb, anrm, bnrm;
    static integer itau;
    static doublereal temp;
    static logical ilvl, ilvr;
    static doublereal anrm1, anrm2, bnrm1, bnrm2, absai, scale, absar, sbeta;
    extern logical lsame_(char *, char *);
    static integer ileft, iinfo, icols, iwork, irows, jc, in;
    extern doublereal dlamch_(char *);
    static integer jr;
    static doublereal salfai;
    extern /* Subroutine */ int zggbak_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublecomplex *,
	     integer *, integer *), zggbal_(char *, integer *,
	     doublecomplex *, integer *, doublecomplex *, integer *, integer *
	    , integer *, doublereal *, doublereal *, doublereal *, integer *);
    static doublereal salfar, safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *);
    static doublereal safmax;
    static char chtemp[1];
    static logical ldumma[1];
    extern doublereal zlange_(char *, integer *, integer *, doublecomplex *, 
	    integer *, doublereal *);
    static integer ijobvl, iright;
    static logical ilimit;
    extern /* Subroutine */ int zgghrd_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, integer *,
	     doublecomplex *, integer *, doublecomplex *, integer *, integer *
	    ), zlascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublecomplex *,
	     integer *, integer *);
    static integer ijobvr;
    extern /* Subroutine */ int zgeqrf_(integer *, integer *, doublecomplex *,
	     integer *, doublecomplex *, doublecomplex *, integer *, integer *
	    );
    static integer lwkmin;
    extern /* Subroutine */ int zlacpy_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *), 
	    zlaset_(char *, integer *, integer *, doublecomplex *, 
	    doublecomplex *, doublecomplex *, integer *), ztgevc_(
	    char *, char *, logical *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, integer *, integer *, doublecomplex *,
	     doublereal *, integer *), zhgeqz_(char *, char *,
	     char *, integer *, integer *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, integer *, doublereal *, integer *);
    static integer irwork, lwkopt;
    extern /* Subroutine */ int zungqr_(integer *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, integer *), zunmqr_(char *, char *, integer *, integer 
	    *, integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, integer *);
    static integer ihi, ilo;
    static doublereal eps;
    static logical ilv;



#define ALPHA(I) alpha[(I)-1]
#define BETA(I) beta[(I)-1]
#define WORK(I) work[(I)-1]
#define RWORK(I) rwork[(I)-1]

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
    } else if (*ldvl < 1 || ilvl && *ldvl < *n) {
	*info = -11;
    } else if (*ldvr < 1 || ilvr && *ldvr < *n) {
	*info = -13;
    } else if (*lwork < lwkmin) {
	*info = -15;
    }

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("ZGEGV ", &i__1);
	return 0;
    }

/*     Quick return if possible */

    WORK(1).r = (doublereal) lwkopt, WORK(1).i = 0.;
    if (*n == 0) {
	return 0;
    }

/*     Get machine constants */

    eps = dlamch_("E") * dlamch_("B");
    safmin = dlamch_("S");
    safmin += safmin;
    safmax = 1. / safmin;

/*     Scale A */

    anrm = zlange_("M", n, n, &A(1,1), lda, &RWORK(1));
    anrm1 = anrm;
    anrm2 = 1.;
    if (anrm < 1.) {
	if (safmax * anrm < 1.) {
	    anrm1 = safmin;
	    anrm2 = safmax * anrm;
	}
    }

    if (anrm > 0.) {
	zlascl_("G", &c_n1, &c_n1, &anrm, &c_b16, n, n, &A(1,1), lda, &
		iinfo);
	if (iinfo != 0) {
	    *info = *n + 10;
	    return 0;
	}
    }

/*     Scale B */

    bnrm = zlange_("M", n, n, &B(1,1), ldb, &RWORK(1));
    bnrm1 = bnrm;
    bnrm2 = 1.;
    if (bnrm < 1.) {
	if (safmax * bnrm < 1.) {
	    bnrm1 = safmin;
	    bnrm2 = safmax * bnrm;
	}
    }

    if (bnrm > 0.) {
	zlascl_("G", &c_n1, &c_n1, &bnrm, &c_b16, n, n, &B(1,1), ldb, &
		iinfo);
	if (iinfo != 0) {
	    *info = *n + 10;
	    return 0;
	}
    }

/*     Permute the matrix to make it more nearly triangular   
       Also "balance" the matrix. */

    ileft = 1;
    iright = *n + 1;
    irwork = iright + *n;
    zggbal_("B", n, &A(1,1), lda, &B(1,1), ldb, &ilo, &ihi, &RWORK(
	    ileft), &RWORK(iright), &RWORK(irwork), &iinfo);
    if (iinfo != 0) {
	*info = *n + 1;
	goto L80;
    }

/*     Reduce B to triangular form, and initialize VL and/or VR */

    irows = ihi + 1 - ilo;
    if (ilv) {
	icols = *n + 1 - ilo;
    } else {
	icols = irows;
    }
    itau = 1;
    iwork = itau + irows;
    i__1 = *lwork + 1 - iwork;
    zgeqrf_(&irows, &icols, &B(ilo,ilo), ldb, &WORK(itau), &WORK(
	    iwork), &i__1, &iinfo);
    if (iinfo >= 0) {
/* Computing MAX */
	i__3 = iwork;
	i__1 = lwkopt, i__2 = (integer) WORK(iwork).r + iwork - 1;
	lwkopt = max(i__1,i__2);
    }
    if (iinfo != 0) {
	*info = *n + 2;
	goto L80;
    }

    i__1 = *lwork + 1 - iwork;
    zunmqr_("L", "C", &irows, &icols, &irows, &B(ilo,ilo), ldb, &
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
	goto L80;
    }

    if (ilvl) {
	zlaset_("Full", n, n, &c_b1, &c_b2, &VL(1,1), ldvl);
	i__1 = irows - 1;
	i__2 = irows - 1;
	zlacpy_("L", &i__1, &i__2, &B(ilo+1,ilo), ldb, &VL(ilo+1,ilo), ldvl);
	i__1 = *lwork + 1 - iwork;
	zungqr_(&irows, &irows, &irows, &VL(ilo,ilo), ldvl, &WORK(
		itau), &WORK(iwork), &i__1, &iinfo);
	if (iinfo >= 0) {
/* Computing MAX */
	    i__3 = iwork;
	    i__1 = lwkopt, i__2 = (integer) WORK(iwork).r + iwork - 1;
	    lwkopt = max(i__1,i__2);
	}
	if (iinfo != 0) {
	    *info = *n + 4;
	    goto L80;
	}
    }

    if (ilvr) {
	zlaset_("Full", n, n, &c_b1, &c_b2, &VR(1,1), ldvr);
    }

/*     Reduce to generalized Hessenberg form */

    if (ilv) {

/*        Eigenvectors requested -- work on whole matrix. */

	zgghrd_(jobvl, jobvr, n, &ilo, &ihi, &A(1,1), lda, &B(1,1), 
		ldb, &VL(1,1), ldvl, &VR(1,1), ldvr, &iinfo);
    } else {
	zgghrd_("N", "N", &irows, &c__1, &irows, &A(ilo,ilo), lda, 
		&B(ilo,ilo), ldb, &VL(1,1), ldvl, &VR(1,1), ldvr, &iinfo);
    }
    if (iinfo != 0) {
	*info = *n + 5;
	goto L80;
    }

/*     Perform QZ algorithm */

    iwork = itau;
    if (ilv) {
	*(unsigned char *)chtemp = 'S';
    } else {
	*(unsigned char *)chtemp = 'E';
    }
    i__1 = *lwork + 1 - iwork;
    zhgeqz_(chtemp, jobvl, jobvr, n, &ilo, &ihi, &A(1,1), lda, &B(1,1), ldb, &ALPHA(1), &BETA(1), &VL(1,1), ldvl, &VR(1,1), ldvr, &WORK(iwork), &i__1, &RWORK(irwork), &iinfo);
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
	goto L80;
    }

    if (ilv) {

/*        Compute Eigenvectors */

	if (ilvl) {
	    if (ilvr) {
		*(unsigned char *)chtemp = 'B';
	    } else {
		*(unsigned char *)chtemp = 'L';
	    }
	} else {
	    *(unsigned char *)chtemp = 'R';
	}

	ztgevc_(chtemp, "B", ldumma, n, &A(1,1), lda, &B(1,1), ldb, 
		&VL(1,1), ldvl, &VR(1,1), ldvr, n, &in, &WORK(
		iwork), &RWORK(irwork), &iinfo);
	if (iinfo != 0) {
	    *info = *n + 7;
	    goto L80;
	}

/*        Undo balancing on VL and VR, rescale */

	if (ilvl) {
	    zggbak_("B", "L", n, &ilo, &ihi, &RWORK(ileft), &RWORK(iright), n,
		     &VL(1,1), ldvl, &iinfo);
	    if (iinfo != 0) {
		*info = *n + 8;
		goto L80;
	    }
	    i__1 = *n;
	    for (jc = 1; jc <= *n; ++jc) {
		temp = 0.;
		i__2 = *n;
		for (jr = 1; jr <= *n; ++jr) {
/* Computing MAX */
		    i__3 = jr + jc * vl_dim1;
		    d__3 = temp, d__4 = (d__1 = VL(jr,jc).r, abs(d__1)) + (
			    d__2 = d_imag(&VL(jr,jc)), abs(d__2));
		    temp = max(d__3,d__4);
/* L10: */
		}
		if (temp < safmin) {
		    goto L30;
		}
		temp = 1. / temp;
		i__2 = *n;
		for (jr = 1; jr <= *n; ++jr) {
		    i__3 = jr + jc * vl_dim1;
		    i__4 = jr + jc * vl_dim1;
		    z__1.r = temp * VL(jr,jc).r, z__1.i = temp * VL(jr,jc).i;
		    VL(jr,jc).r = z__1.r, VL(jr,jc).i = z__1.i;
/* L20: */
		}
L30:
		;
	    }
	}
	if (ilvr) {
	    zggbak_("B", "R", n, &ilo, &ihi, &RWORK(ileft), &RWORK(iright), n,
		     &VR(1,1), ldvr, &iinfo);
	    if (iinfo != 0) {
		*info = *n + 9;
		goto L80;
	    }
	    i__1 = *n;
	    for (jc = 1; jc <= *n; ++jc) {
		temp = 0.;
		i__2 = *n;
		for (jr = 1; jr <= *n; ++jr) {
/* Computing MAX */
		    i__3 = jr + jc * vr_dim1;
		    d__3 = temp, d__4 = (d__1 = VR(jr,jc).r, abs(d__1)) + (
			    d__2 = d_imag(&VR(jr,jc)), abs(d__2));
		    temp = max(d__3,d__4);
/* L40: */
		}
		if (temp < safmin) {
		    goto L60;
		}
		temp = 1. / temp;
		i__2 = *n;
		for (jr = 1; jr <= *n; ++jr) {
		    i__3 = jr + jc * vr_dim1;
		    i__4 = jr + jc * vr_dim1;
		    z__1.r = temp * VR(jr,jc).r, z__1.i = temp * VR(jr,jc).i;
		    VR(jr,jc).r = z__1.r, VR(jr,jc).i = z__1.i;
/* L50: */
		}
L60:
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
	i__2 = jc;
	absar = (d__1 = ALPHA(jc).r, abs(d__1));
	absai = (d__1 = d_imag(&ALPHA(jc)), abs(d__1));
	i__2 = jc;
	absb = (d__1 = BETA(jc).r, abs(d__1));
	i__2 = jc;
	salfar = anrm * ALPHA(jc).r;
	salfai = anrm * d_imag(&ALPHA(jc));
	i__2 = jc;
	sbeta = bnrm * BETA(jc).r;
	ilimit = FALSE_;
	scale = 1.;

/*        Check for significant underflow in imaginary part of ALPHA 
  

   Computing MAX */
	d__1 = safmin, d__2 = eps * absar, d__1 = max(d__1,d__2), d__2 = eps *
		 absb;
	if (abs(salfai) < safmin && absai >= max(d__1,d__2)) {
	    ilimit = TRUE_;
/* Computing MAX */
	    d__1 = safmin, d__2 = anrm2 * absai;
	    scale = safmin / anrm1 / max(d__1,d__2);
	}

/*        Check for significant underflow in real part of ALPHA   

   Computing MAX */
	d__1 = safmin, d__2 = eps * absai, d__1 = max(d__1,d__2), d__2 = eps *
		 absb;
	if (abs(salfar) < safmin && absar >= max(d__1,d__2)) {
	    ilimit = TRUE_;
/* Computing MAX   
   Computing MAX */
	    d__3 = safmin, d__4 = anrm2 * absar;
	    d__1 = scale, d__2 = safmin / anrm1 / max(d__3,d__4);
	    scale = max(d__1,d__2);
	}

/*        Check for significant underflow in BETA   

   Computing MAX */
	d__1 = safmin, d__2 = eps * absar, d__1 = max(d__1,d__2), d__2 = eps *
		 absai;
	if (abs(sbeta) < safmin && absb >= max(d__1,d__2)) {
	    ilimit = TRUE_;
/* Computing MAX   
   Computing MAX */
	    d__3 = safmin, d__4 = bnrm2 * absb;
	    d__1 = scale, d__2 = safmin / bnrm1 / max(d__3,d__4);
	    scale = max(d__1,d__2);
	}

/*        Check for possible overflow when limiting scaling */

	if (ilimit) {
/* Computing MAX */
	    d__1 = abs(salfar), d__2 = abs(salfai), d__1 = max(d__1,d__2), 
		    d__2 = abs(sbeta);
	    temp = scale * safmin * max(d__1,d__2);
	    if (temp > 1.) {
		scale /= temp;
	    }
	    if (scale < 1.) {
		ilimit = FALSE_;
	    }
	}

/*        Recompute un-scaled ALPHA, BETA if necessary. */

	if (ilimit) {
	    i__2 = jc;
	    salfar = scale * ALPHA(jc).r * anrm;
	    salfai = scale * d_imag(&ALPHA(jc)) * anrm;
	    i__2 = jc;
	    z__2.r = scale * BETA(jc).r, z__2.i = scale * BETA(jc).i;
	    z__1.r = bnrm * z__2.r, z__1.i = bnrm * z__2.i;
	    sbeta = z__1.r;
	}
	i__2 = jc;
	z__1.r = salfar, z__1.i = salfai;
	ALPHA(jc).r = z__1.r, ALPHA(jc).i = z__1.i;
	i__2 = jc;
	BETA(jc).r = sbeta, BETA(jc).i = 0.;
/* L70: */
    }

L80:
    WORK(1).r = (doublereal) lwkopt, WORK(1).i = 0.;

    return 0;

/*     End of ZGEGV */

} /* zgegv_ */

