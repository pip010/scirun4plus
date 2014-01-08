#include "f2c.h"

/* Subroutine */ int sgeev_(char *jobvl, char *jobvr, integer *n, real *a, 
	integer *lda, real *wr, real *wi, real *vl, integer *ldvl, real *vr, 
	integer *ldvr, real *work, integer *lwork, integer *info)
{
/*  -- LAPACK driver routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    SGEEV computes for an N-by-N real nonsymmetric matrix A, the   
    eigenvalues and, optionally, the left and/or right eigenvectors.   

    The right eigenvector v(j) of A satisfies   
                     A * v(j) = lambda(j) * v(j)   
    where lambda(j) is its eigenvalue.   
    The left eigenvector u(j) of A satisfies   
                  u(j)**H * A = lambda(j) * u(j)**H   
    where u(j)**H denotes the conjugate transpose of u(j).   

    The computed eigenvectors are normalized to have Euclidean norm   
    equal to 1 and largest component real.   

    Arguments   
    =========   

    JOBVL   (input) CHARACTER*1   
            = 'N': left eigenvectors of A are not computed;   
            = 'V': left eigenvectors of A are computed.   

    JOBVR   (input) CHARACTER*1   
            = 'N': right eigenvectors of A are not computed;   
            = 'V': right eigenvectors of A are computed.   

    N       (input) INTEGER   
            The order of the matrix A. N >= 0.   

    A       (input/output) REAL array, dimension (LDA,N)   
            On entry, the N-by-N matrix A.   
            On exit, A has been overwritten.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(1,N).   

    WR      (output) REAL array, dimension (N)   
    WI      (output) REAL array, dimension (N)   
            WR and WI contain the real and imaginary parts,   
            respectively, of the computed eigenvalues.  Complex   
            conjugate pairs of eigenvalues appear consecutively   
            with the eigenvalue having the positive imaginary part   
            first.   

    VL      (output) REAL array, dimension (LDVL,N)   
            If JOBVL = 'V', the left eigenvectors u(j) are stored one   
            after another in the columns of VL, in the same order   
            as their eigenvalues.   
            If JOBVL = 'N', VL is not referenced.   
            If the j-th eigenvalue is real, then u(j) = VL(:,j),   
            the j-th column of VL.   
            If the j-th and (j+1)-st eigenvalues form a complex   
            conjugate pair, then u(j) = VL(:,j) + i*VL(:,j+1) and   
            u(j+1) = VL(:,j) - i*VL(:,j+1).   

    LDVL    (input) INTEGER   
            The leading dimension of the array VL.  LDVL >= 1; if   
            JOBVL = 'V', LDVL >= N.   

    VR      (output) REAL array, dimension (LDVR,N)   
            If JOBVR = 'V', the right eigenvectors v(j) are stored one   
            after another in the columns of VR, in the same order   
            as their eigenvalues.   
            If JOBVR = 'N', VR is not referenced.   
            If the j-th eigenvalue is real, then v(j) = VR(:,j),   
            the j-th column of VR.   
            If the j-th and (j+1)-st eigenvalues form a complex   
            conjugate pair, then v(j) = VR(:,j) + i*VR(:,j+1) and   
            v(j+1) = VR(:,j) - i*VR(:,j+1).   

    LDVR    (input) INTEGER   
            The leading dimension of the array VR.  LDVR >= 1; if   
            JOBVR = 'V', LDVR >= N.   

    WORK    (workspace/output) REAL array, dimension (LWORK)   
            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INTEGER   
            The dimension of the array WORK.  LWORK >= max(1,3*N), and   
            if JOBVL = 'V' or JOBVR = 'V', LWORK >= 4*N.  For good   
            performance, LWORK must generally be larger.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value.   
            > 0:  if INFO = i, the QR algorithm failed to compute all the 
  
                  eigenvalues, and no eigenvectors have been computed;   
                  elements i+1:N of WR and WI contain eigenvalues which   
                  have converged.   

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
    integer a_dim1, a_offset, vl_dim1, vl_offset, vr_dim1, vr_offset, i__1, 
	    i__2, i__3, i__4;
    real r__1, r__2;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    static integer ibal;
    static char side[1];
    static integer maxb;
    static real anrm;
    static integer ierr, itau, iwrk, nout;
    extern /* Subroutine */ int srot_(integer *, real *, integer *, real *, 
	    integer *, real *, real *);
    extern doublereal snrm2_(integer *, real *, integer *);
    static integer i, k;
    static real r;
    extern logical lsame_(char *, char *);
    extern /* Subroutine */ int sscal_(integer *, real *, real *, integer *);
    extern doublereal slapy2_(real *, real *);
    static real cs;
    extern /* Subroutine */ int slabad_(real *, real *);
    static logical scalea;
    static real cscale;
    extern /* Subroutine */ int sgebak_(char *, char *, integer *, integer *, 
	    integer *, real *, integer *, real *, integer *, integer *), sgebal_(char *, integer *, real *, integer *, 
	    integer *, integer *, real *, integer *);
    static real sn;
    extern doublereal slamch_(char *), slange_(char *, integer *, 
	    integer *, real *, integer *, real *);
    extern /* Subroutine */ int sgehrd_(integer *, integer *, integer *, real 
	    *, integer *, real *, real *, integer *, integer *), xerbla_(char 
	    *, integer *);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static logical select[1];
    static real bignum;
    extern /* Subroutine */ int slascl_(char *, integer *, integer *, real *, 
	    real *, integer *, integer *, real *, integer *, integer *);
    extern integer isamax_(integer *, real *, integer *);
    extern /* Subroutine */ int slacpy_(char *, integer *, integer *, real *, 
	    integer *, real *, integer *), slartg_(real *, real *, 
	    real *, real *, real *), sorghr_(integer *, integer *, integer *, 
	    real *, integer *, real *, real *, integer *, integer *), shseqr_(
	    char *, char *, integer *, integer *, integer *, real *, integer *
	    , real *, real *, real *, integer *, real *, integer *, integer *), strevc_(char *, char *, logical *, integer *, 
	    real *, integer *, real *, integer *, real *, integer *, integer *
	    , integer *, real *, integer *);
    static integer minwrk, maxwrk;
    static logical wantvl;
    static real smlnum;
    static integer hswork;
    static logical wantvr;
    static integer ihi;
    static real scl;
    static integer ilo;
    static real dum[1], eps;



#define WR(I) wr[(I)-1]
#define WI(I) wi[(I)-1]
#define WORK(I) work[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]
#define VL(I,J) vl[(I)-1 + ((J)-1)* ( *ldvl)]
#define VR(I,J) vr[(I)-1 + ((J)-1)* ( *ldvr)]

    *info = 0;
    wantvl = lsame_(jobvl, "V");
    wantvr = lsame_(jobvr, "V");
    if (! wantvl && ! lsame_(jobvl, "N")) {
	*info = -1;
    } else if (! wantvr && ! lsame_(jobvr, "N")) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*lda < max(1,*n)) {
	*info = -5;
    } else if (*ldvl < 1 || wantvl && *ldvl < *n) {
	*info = -9;
    } else if (*ldvr < 1 || wantvr && *ldvr < *n) {
	*info = -11;
    }

/*     Compute workspace   
        (Note: Comments in the code beginning "Workspace:" describe the   
         minimal amount of workspace needed at that point in the code,   
         as well as the preferred amount for good performance.   
         NB refers to the optimal block size for the immediately   
         following subroutine, as returned by ILAENV.   
         HSWORK refers to the workspace preferred by SHSEQR, as   
         calculated below. HSWORK is computed assuming ILO=1 and IHI=N,   
         the worst case.) */

    minwrk = 1;
    if (*info == 0 && *lwork >= 1) {
	maxwrk = (*n << 1) + *n * ilaenv_(&c__1, "SGEHRD", " ", n, &c__1, n, &
		c__0, 6L, 1L);
	if (! wantvl && ! wantvr) {
/* Computing MAX */
	    i__1 = 1, i__2 = *n * 3;
	    minwrk = max(i__1,i__2);
/* Computing MAX */
	    i__1 = ilaenv_(&c__8, "SHSEQR", "EN", n, &c__1, n, &c_n1, 6L, 2L);
	    maxb = max(i__1,2);
/* Computing MIN   
   Computing MAX */
	    i__3 = 2, i__4 = ilaenv_(&c__4, "SHSEQR", "EN", n, &c__1, n, &
		    c_n1, 6L, 2L);
	    i__1 = min(maxb,*n), i__2 = max(i__3,i__4);
	    k = min(i__1,i__2);
/* Computing MAX */
	    i__1 = k * (k + 2), i__2 = *n << 1;
	    hswork = max(i__1,i__2);
/* Computing MAX */
	    i__1 = maxwrk, i__2 = *n + 1, i__1 = max(i__1,i__2), i__2 = *n + 
		    hswork;
	    maxwrk = max(i__1,i__2);
	} else {
/* Computing MAX */
	    i__1 = 1, i__2 = *n << 2;
	    minwrk = max(i__1,i__2);
/* Computing MAX */
	    i__1 = maxwrk, i__2 = (*n << 1) + (*n - 1) * ilaenv_(&c__1, "SOR"
		    "GHR", " ", n, &c__1, n, &c_n1, 6L, 1L);
	    maxwrk = max(i__1,i__2);
/* Computing MAX */
	    i__1 = ilaenv_(&c__8, "SHSEQR", "SV", n, &c__1, n, &c_n1, 6L, 2L);
	    maxb = max(i__1,2);
/* Computing MIN   
   Computing MAX */
	    i__3 = 2, i__4 = ilaenv_(&c__4, "SHSEQR", "SV", n, &c__1, n, &
		    c_n1, 6L, 2L);
	    i__1 = min(maxb,*n), i__2 = max(i__3,i__4);
	    k = min(i__1,i__2);
/* Computing MAX */
	    i__1 = k * (k + 2), i__2 = *n << 1;
	    hswork = max(i__1,i__2);
/* Computing MAX */
	    i__1 = maxwrk, i__2 = *n + 1, i__1 = max(i__1,i__2), i__2 = *n + 
		    hswork;
	    maxwrk = max(i__1,i__2);
/* Computing MAX */
	    i__1 = maxwrk, i__2 = *n << 2;
	    maxwrk = max(i__1,i__2);
	}
	WORK(1) = (real) maxwrk;
    }
    if (*lwork < minwrk) {
	*info = -13;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("SGEEV ", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0) {
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

    anrm = slange_("M", n, n, &A(1,1), lda, dum);
    scalea = FALSE_;
    if (anrm > 0.f && anrm < smlnum) {
	scalea = TRUE_;
	cscale = smlnum;
    } else if (anrm > bignum) {
	scalea = TRUE_;
	cscale = bignum;
    }
    if (scalea) {
	slascl_("G", &c__0, &c__0, &anrm, &cscale, n, n, &A(1,1), lda, &
		ierr);
    }

/*     Balance the matrix   
       (Workspace: need N) */

    ibal = 1;
    sgebal_("B", n, &A(1,1), lda, &ilo, &ihi, &WORK(ibal), &ierr);

/*     Reduce to upper Hessenberg form   
       (Workspace: need 3*N, prefer 2*N+N*NB) */

    itau = ibal + *n;
    iwrk = itau + *n;
    i__1 = *lwork - iwrk + 1;
    sgehrd_(n, &ilo, &ihi, &A(1,1), lda, &WORK(itau), &WORK(iwrk), &i__1,
	     &ierr);

    if (wantvl) {

/*        Want left eigenvectors   
          Copy Householder vectors to VL */

	*(unsigned char *)side = 'L';
	slacpy_("L", n, n, &A(1,1), lda, &VL(1,1), ldvl);

/*        Generate orthogonal matrix in VL   
          (Workspace: need 3*N-1, prefer 2*N+(N-1)*NB) */

	i__1 = *lwork - iwrk + 1;
	sorghr_(n, &ilo, &ihi, &VL(1,1), ldvl, &WORK(itau), &WORK(iwrk),
		 &i__1, &ierr);

/*        Perform QR iteration, accumulating Schur vectors in VL   
          (Workspace: need N+1, prefer N+HSWORK (see comments) ) */

	iwrk = itau;
	i__1 = *lwork - iwrk + 1;
	shseqr_("S", "V", n, &ilo, &ihi, &A(1,1), lda, &WR(1), &WI(1), &
		VL(1,1), ldvl, &WORK(iwrk), &i__1, info);

	if (wantvr) {

/*           Want left and right eigenvectors   
             Copy Schur vectors to VR */

	    *(unsigned char *)side = 'B';
	    slacpy_("F", n, n, &VL(1,1), ldvl, &VR(1,1), ldvr)
		    ;
	}

    } else if (wantvr) {

/*        Want right eigenvectors   
          Copy Householder vectors to VR */

	*(unsigned char *)side = 'R';
	slacpy_("L", n, n, &A(1,1), lda, &VR(1,1), ldvr);

/*        Generate orthogonal matrix in VR   
          (Workspace: need 3*N-1, prefer 2*N+(N-1)*NB) */

	i__1 = *lwork - iwrk + 1;
	sorghr_(n, &ilo, &ihi, &VR(1,1), ldvr, &WORK(itau), &WORK(iwrk),
		 &i__1, &ierr);

/*        Perform QR iteration, accumulating Schur vectors in VR   
          (Workspace: need N+1, prefer N+HSWORK (see comments) ) */

	iwrk = itau;
	i__1 = *lwork - iwrk + 1;
	shseqr_("S", "V", n, &ilo, &ihi, &A(1,1), lda, &WR(1), &WI(1), &
		VR(1,1), ldvr, &WORK(iwrk), &i__1, info);

    } else {

/*        Compute eigenvalues only   
          (Workspace: need N+1, prefer N+HSWORK (see comments) ) */

	iwrk = itau;
	i__1 = *lwork - iwrk + 1;
	shseqr_("E", "N", n, &ilo, &ihi, &A(1,1), lda, &WR(1), &WI(1), &
		VR(1,1), ldvr, &WORK(iwrk), &i__1, info);
    }

/*     If INFO > 0 from SHSEQR, then quit */

    if (*info > 0) {
	goto L50;
    }

    if (wantvl || wantvr) {

/*        Compute left and/or right eigenvectors   
          (Workspace: need 4*N) */

	strevc_(side, "B", select, n, &A(1,1), lda, &VL(1,1), ldvl,
		 &VR(1,1), ldvr, n, &nout, &WORK(iwrk), &ierr);
    }

    if (wantvl) {

/*        Undo balancing of left eigenvectors   
          (Workspace: need N) */

	sgebak_("B", "L", n, &ilo, &ihi, &WORK(ibal), n, &VL(1,1), ldvl,
		 &ierr);

/*        Normalize left eigenvectors and make largest component real 
*/

	i__1 = *n;
	for (i = 1; i <= *n; ++i) {
	    if (WI(i) == 0.f) {
		scl = 1.f / snrm2_(n, &VL(1,i), &c__1);
		sscal_(n, &scl, &VL(1,i), &c__1);
	    } else if (WI(i) > 0.f) {
		r__1 = snrm2_(n, &VL(1,i), &c__1);
		r__2 = snrm2_(n, &VL(1,i+1), &c__1);
		scl = 1.f / slapy2_(&r__1, &r__2);
		sscal_(n, &scl, &VL(1,i), &c__1);
		sscal_(n, &scl, &VL(1,i+1), &c__1);
		i__2 = *n;
		for (k = 1; k <= *n; ++k) {
/* Computing 2nd power */
		    r__1 = VL(k,i);
/* Computing 2nd power */
		    r__2 = VL(k,i+1);
		    WORK(iwrk + k - 1) = r__1 * r__1 + r__2 * r__2;
/* L10: */
		}
		k = isamax_(n, &WORK(iwrk), &c__1);
		slartg_(&VL(k,i), &VL(k,i+1), &cs,
			 &sn, &r);
		srot_(n, &VL(1,i), &c__1, &VL(1,i+1), &c__1, &cs, &sn);
		VL(k,i+1) = 0.f;
	    }
/* L20: */
	}
    }

    if (wantvr) {

/*        Undo balancing of right eigenvectors   
          (Workspace: need N) */

	sgebak_("B", "R", n, &ilo, &ihi, &WORK(ibal), n, &VR(1,1), ldvr,
		 &ierr);

/*        Normalize right eigenvectors and make largest component real
 */

	i__1 = *n;
	for (i = 1; i <= *n; ++i) {
	    if (WI(i) == 0.f) {
		scl = 1.f / snrm2_(n, &VR(1,i), &c__1);
		sscal_(n, &scl, &VR(1,i), &c__1);
	    } else if (WI(i) > 0.f) {
		r__1 = snrm2_(n, &VR(1,i), &c__1);
		r__2 = snrm2_(n, &VR(1,i+1), &c__1);
		scl = 1.f / slapy2_(&r__1, &r__2);
		sscal_(n, &scl, &VR(1,i), &c__1);
		sscal_(n, &scl, &VR(1,i+1), &c__1);
		i__2 = *n;
		for (k = 1; k <= *n; ++k) {
/* Computing 2nd power */
		    r__1 = VR(k,i);
/* Computing 2nd power */
		    r__2 = VR(k,i+1);
		    WORK(iwrk + k - 1) = r__1 * r__1 + r__2 * r__2;
/* L30: */
		}
		k = isamax_(n, &WORK(iwrk), &c__1);
		slartg_(&VR(k,i), &VR(k,i+1), &cs,
			 &sn, &r);
		srot_(n, &VR(1,i), &c__1, &VR(1,i+1), &c__1, &cs, &sn);
		VR(k,i+1) = 0.f;
	    }
/* L40: */
	}
    }

/*     Undo scaling if necessary */

L50:
    if (scalea) {
	i__1 = *n - *info;
/* Computing MAX */
	i__3 = *n - *info;
	i__2 = max(i__3,1);
	slascl_("G", &c__0, &c__0, &cscale, &anrm, &i__1, &c__1, &WR(*info + 
		1), &i__2, &ierr);
	i__1 = *n - *info;
/* Computing MAX */
	i__3 = *n - *info;
	i__2 = max(i__3,1);
	slascl_("G", &c__0, &c__0, &cscale, &anrm, &i__1, &c__1, &WI(*info + 
		1), &i__2, &ierr);
	if (*info > 0) {
	    i__1 = ilo - 1;
	    slascl_("G", &c__0, &c__0, &cscale, &anrm, &i__1, &c__1, &WR(1), 
		    n, &ierr);
	    i__1 = ilo - 1;
	    slascl_("G", &c__0, &c__0, &cscale, &anrm, &i__1, &c__1, &WI(1), 
		    n, &ierr);
	}
    }

    WORK(1) = (real) maxwrk;
    return 0;

/*     End of SGEEV */

} /* sgeev_ */

