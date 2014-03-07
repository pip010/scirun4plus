#include "f2c.h"

/* Subroutine */ int cgees_(char *jobvs, char *sort, L_fp select, integer *n, 
	complex *a, integer *lda, integer *sdim, complex *w, complex *vs, 
	integer *ldvs, complex *work, integer *lwork, real *rwork, logical *
	bwork, integer *info)
{
/*  -- LAPACK driver routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       March 31, 1993   


    Purpose   
    =======   

    CGEES computes for an N-by-N complex nonsymmetric matrix A, the   
    eigenvalues, the Schur form T, and, optionally, the matrix of Schur   
    vectors Z.  This gives the Schur factorization A = Z*T*(Z**H).   

    Optionally, it also orders the eigenvalues on the diagonal of the   
    Schur form so that selected eigenvalues are at the top left.   
    The leading columns of Z then form an orthonormal basis for the   
    invariant subspace corresponding to the selected eigenvalues.   
    A complex matrix is in Schur form if it is upper triangular.   

    Arguments   
    =========   

    JOBVS   (input) CHARACTER*1   
            = 'N': Schur vectors are not computed;   
            = 'V': Schur vectors are computed.   

    SORT    (input) CHARACTER*1   
            Specifies whether or not to order the eigenvalues on the   
            diagonal of the Schur form.   
            = 'N': Eigenvalues are not ordered:   
            = 'S': Eigenvalues are ordered (see SELECT).   

    SELECT  (input) LOGICAL FUNCTION of one COMPLEX argument   
            SELECT must be declared EXTERNAL in the calling subroutine.   
            If SORT = 'S', SELECT is used to select eigenvalues to order 
  
            to the top left of the Schur form.   
            IF SORT = 'N', SELECT is not referenced.   
            The eigenvalue W(j) is selected if SELECT(W(j)) is true.   

    N       (input) INTEGER   
            The order of the matrix A. N >= 0.   

    A       (input/output) COMPLEX array, dimension (LDA,N)   
            On entry, the N-by-N matrix A.   
            On exit, A has been overwritten by its Schur form T.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(1,N).   

    SDIM    (output) INTEGER   
            If SORT = 'N', SDIM = 0.   
            If SORT = 'S', SDIM = number of eigenvalues for which   
                           SELECT is true.   

    W       (output) COMPLEX array, dimension (N)   
            W contains the computed eigenvalues, in the same order that   
            they appear on the diagonal of the output Schur form T.   

    VS      (output) COMPLEX array, dimension (LDVS,N)   
            If JOBVS = 'V', VS contains the unitary matrix Z of Schur   
            vectors.   
            If JOBVS = 'N', VS is not referenced.   

    LDVS    (input) INTEGER   
            The leading dimension of the array VS.  LDVS >= 1; if   
            JOBVS = 'V', LDVS >= N.   

    WORK    (workspace/output) COMPLEX array, dimension (LWORK)   
            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INTEGER   
            The dimension of the array WORK.  LWORK >= max(1,2*N).   
            For good performance, LWORK must generally be larger.   

    RWORK   (workspace) REAL array, dimension (N)   

    BWORK   (workspace) LOGICAL array, dimension (N)   
            Not referenced if SORT = 'N'.   

    INFO    (output) INTEGER   
            = 0: successful exit   
            < 0: if INFO = -i, the i-th argument had an illegal value.   
            > 0: if INFO = i, and i is   
                 <= N:  the QR algorithm failed to compute all the   
                        eigenvalues; elements 1:ILO-1 and i+1:N of W   
                        contain those eigenvalues which have converged;   
                        if JOBVS = 'V', VS contains the matrix which   
                        reduces A to its partially converged Schur form. 
  
                 = N+1: the eigenvalues could not be reordered because   
                        some eigenvalues were too close to separate (the 
  
                        problem is very ill-conditioned);   
                 = N+2: after reordering, roundoff changed values of   
                        some complex eigenvalues so that leading   
                        eigenvalues in the Schur form no longer satisfy   
                        SELECT = .TRUE..  This could also be caused by   
                        underflow due to scaling.   

    ===================================================================== 
  


       Test the input arguments   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static integer c__1 = 1;
    static integer c__0 = 0;
    static integer c__8 = 8;
    static integer c_n1 = -1;
    static integer c__4 = 4;
    
    /* System generated locals */
    integer a_dim1, a_offset, vs_dim1, vs_offset, i__1, i__2, i__3, i__4;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    static integer ibal, maxb;
    static real anrm;
    static integer ierr, itau, iwrk, i, k;
    static real s;
    static integer icond, ieval;
    extern logical lsame_(char *, char *);
    extern /* Subroutine */ int ccopy_(integer *, complex *, integer *, 
	    complex *, integer *), cgebak_(char *, char *, integer *, integer 
	    *, integer *, real *, integer *, complex *, integer *, integer *), cgebal_(char *, integer *, complex *, integer *, 
	    integer *, integer *, real *, integer *), slabad_(real *, 
	    real *);
    static logical scalea;
    extern doublereal clange_(char *, integer *, integer *, complex *, 
	    integer *, real *);
    static real cscale;
    extern /* Subroutine */ int cgehrd_(integer *, integer *, integer *, 
	    complex *, integer *, complex *, complex *, integer *, integer *),
	     clascl_(char *, integer *, integer *, real *, real *, integer *, 
	    integer *, complex *, integer *, integer *);
    extern doublereal slamch_(char *);
    extern /* Subroutine */ int clacpy_(char *, integer *, integer *, complex 
	    *, integer *, complex *, integer *), xerbla_(char *, 
	    integer *);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static real bignum;
    extern /* Subroutine */ int chseqr_(char *, char *, integer *, integer *, 
	    integer *, complex *, integer *, complex *, complex *, integer *, 
	    complex *, integer *, integer *), cunghr_(integer 
	    *, integer *, integer *, complex *, integer *, complex *, complex 
	    *, integer *, integer *), ctrsen_(char *, char *, logical *, 
	    integer *, complex *, integer *, complex *, integer *, complex *, 
	    integer *, real *, real *, complex *, integer *, integer *);
    static integer minwrk, maxwrk;
    static real smlnum;
    static integer hswork;
    static logical wantst, wantvs;
    static integer ihi, ilo;
    static real dum[1], eps, sep;



#define W(I) w[(I)-1]
#define WORK(I) work[(I)-1]
#define RWORK(I) rwork[(I)-1]
#define BWORK(I) bwork[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]
#define VS(I,J) vs[(I)-1 + ((J)-1)* ( *ldvs)]

    *info = 0;
    wantvs = lsame_(jobvs, "V");
    wantst = lsame_(sort, "S");
    if (! wantvs && ! lsame_(jobvs, "N")) {
	*info = -1;
    } else if (! wantst && ! lsame_(sort, "N")) {
	*info = -2;
    } else if (*n < 0) {
	*info = -4;
    } else if (*lda < max(1,*n)) {
	*info = -6;
    } else if (*ldvs < 1 || wantvs && *ldvs < *n) {
	*info = -10;
    }

/*     Compute workspace   
        (Note: Comments in the code beginning "Workspace:" describe the   
         minimal amount of workspace needed at that point in the code,   
         as well as the preferred amount for good performance.   
         CWorkspace refers to complex workspace, and RWorkspace to real   
         workspace. NB refers to the optimal block size for the   
         immediately following subroutine, as returned by ILAENV.   
         HSWORK refers to the workspace preferred by CHSEQR, as   
         calculated below. HSWORK is computed assuming ILO=1 and IHI=N,   
         the worst case.) */

    minwrk = 1;
    if (*info == 0 && *lwork >= 1) {
	maxwrk = *n + *n * ilaenv_(&c__1, "CGEHRD", " ", n, &c__1, n, &c__0, 
		6L, 1L);
/* Computing MAX */
	i__1 = 1, i__2 = *n * 3;
	minwrk = max(i__1,i__2);
	if (! wantvs) {
/* Computing MAX */
	    i__1 = ilaenv_(&c__8, "CHSEQR", "SN", n, &c__1, n, &c_n1, 6L, 2L);
	    maxb = max(i__1,2);
/* Computing MIN   
   Computing MAX */
	    i__3 = 2, i__4 = ilaenv_(&c__4, "CHSEQR", "SN", n, &c__1, n, &
		    c_n1, 6L, 2L);
	    i__1 = min(maxb,*n), i__2 = max(i__3,i__4);
	    k = min(i__1,i__2);
/* Computing MAX */
	    i__1 = k * (k + 2), i__2 = *n << 1;
	    hswork = max(i__1,i__2);
/* Computing MAX */
	    i__1 = max(maxwrk,hswork);
	    maxwrk = max(i__1,1);
	} else {
/* Computing MAX */
	    i__1 = maxwrk, i__2 = *n + (*n - 1) * ilaenv_(&c__1, "CUNGHR", 
		    " ", n, &c__1, n, &c_n1, 6L, 1L);
	    maxwrk = max(i__1,i__2);
/* Computing MAX */
	    i__1 = ilaenv_(&c__8, "CHSEQR", "EN", n, &c__1, n, &c_n1, 6L, 2L);
	    maxb = max(i__1,2);
/* Computing MIN   
   Computing MAX */
	    i__3 = 2, i__4 = ilaenv_(&c__4, "CHSEQR", "EN", n, &c__1, n, &
		    c_n1, 6L, 2L);
	    i__1 = min(maxb,*n), i__2 = max(i__3,i__4);
	    k = min(i__1,i__2);
/* Computing MAX */
	    i__1 = k * (k + 2), i__2 = *n << 1;
	    hswork = max(i__1,i__2);
/* Computing MAX */
	    i__1 = max(maxwrk,hswork);
	    maxwrk = max(i__1,1);
	}
	WORK(1).r = (real) maxwrk, WORK(1).i = 0.f;
    }
    if (*lwork < minwrk) {
	*info = -12;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("CGEES ", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0) {
	*sdim = 0;
	return 0;
    }

/*     Get machine constants */

    eps = slamch_("P");
    smlnum = slamch_("S");
    bignum = 1.f / smlnum;
    slabad_(&smlnum, &bignum);
    smlnum = sqrt(smlnum) / eps;
    bignum = 1.f / smlnum;

/*     Scale A if max element outside range [SMLNUM,BIGNUM] */

    anrm = clange_("M", n, n, &A(1,1), lda, dum);
    scalea = FALSE_;
    if (anrm > 0.f && anrm < smlnum) {
	scalea = TRUE_;
	cscale = smlnum;
    } else if (anrm > bignum) {
	scalea = TRUE_;
	cscale = bignum;
    }
    if (scalea) {
	clascl_("G", &c__0, &c__0, &anrm, &cscale, n, n, &A(1,1), lda, &
		ierr);
    }

/*     Permute the matrix to make it more nearly triangular   
       (CWorkspace: none)   
       (RWorkspace: need N) */

    ibal = 1;
    cgebal_("P", n, &A(1,1), lda, &ilo, &ihi, &RWORK(ibal), &ierr);

/*     Reduce to upper Hessenberg form   
       (CWorkspace: need 2*N, prefer N+N*NB)   
       (RWorkspace: none) */

    itau = 1;
    iwrk = *n + itau;
    i__1 = *lwork - iwrk + 1;
    cgehrd_(n, &ilo, &ihi, &A(1,1), lda, &WORK(itau), &WORK(iwrk), &i__1,
	     &ierr);

    if (wantvs) {

/*        Copy Householder vectors to VS */

	clacpy_("L", n, n, &A(1,1), lda, &VS(1,1), ldvs);

/*        Generate unitary matrix in VS   
          (CWorkspace: need 2*N-1, prefer N+(N-1)*NB)   
          (RWorkspace: none) */

	i__1 = *lwork - iwrk + 1;
	cunghr_(n, &ilo, &ihi, &VS(1,1), ldvs, &WORK(itau), &WORK(iwrk),
		 &i__1, &ierr);
    }

    *sdim = 0;

/*     Perform QR iteration, accumulating Schur vectors in VS if desired 
  
       (CWorkspace: need 1, prefer HSWORK (see comments) )   
       (RWorkspace: none) */

    iwrk = itau;
    i__1 = *lwork - iwrk + 1;
    chseqr_("S", jobvs, n, &ilo, &ihi, &A(1,1), lda, &W(1), &VS(1,1), ldvs, &WORK(iwrk), &i__1, &ieval);
    if (ieval > 0) {
	*info = ieval;
    }

/*     Sort eigenvalues if desired */

    if (wantst && *info == 0) {
	if (scalea) {
	    clascl_("G", &c__0, &c__0, &cscale, &anrm, n, &c__1, &W(1), n, &
		    ierr);
	}
	i__1 = *n;
	for (i = 1; i <= *n; ++i) {
	    BWORK(i) = (*select)(&W(i));
/* L10: */
	}

/*        Reorder eigenvalues and transform Schur vectors   
          (CWorkspace: none)   
          (RWorkspace: none) */

	i__1 = *lwork - iwrk + 1;
	ctrsen_("N", jobvs, &BWORK(1), n, &A(1,1), lda, &VS(1,1), 
		ldvs, &W(1), sdim, &s, &sep, &WORK(iwrk), &i__1, &icond);
    }

    if (wantvs) {

/*        Undo balancing   
          (CWorkspace: none)   
          (RWorkspace: need N) */

	cgebak_("P", "R", n, &ilo, &ihi, &RWORK(ibal), n, &VS(1,1), 
		ldvs, &ierr);
    }

    if (scalea) {

/*        Undo scaling for the Schur form of A */

	clascl_("U", &c__0, &c__0, &cscale, &anrm, n, n, &A(1,1), lda, &
		ierr);
	i__1 = *lda + 1;
	ccopy_(n, &A(1,1), &i__1, &W(1), &c__1);
    }

    WORK(1).r = (real) maxwrk, WORK(1).i = 0.f;
    return 0;

/*     End of CGEES */

} /* cgees_ */

