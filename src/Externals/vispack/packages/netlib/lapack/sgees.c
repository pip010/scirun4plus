#include "f2c.h"

/* Subroutine */ int sgees_(char *jobvs, char *sort, L_fp select, integer *n, 
	real *a, integer *lda, integer *sdim, real *wr, real *wi, real *vs, 
	integer *ldvs, real *work, integer *lwork, logical *bwork, integer *
	info)
{
/*  -- LAPACK driver routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       March 31, 1993   


    Purpose   
    =======   

    SGEES computes for an N-by-N real nonsymmetric matrix A, the   
    eigenvalues, the real Schur form T, and, optionally, the matrix of   
    Schur vectors Z.  This gives the Schur factorization A = Z*T*(Z**T). 
  

    Optionally, it also orders the eigenvalues on the diagonal of the   
    real Schur form so that selected eigenvalues are at the top left.   
    The leading columns of Z then form an orthonormal basis for the   
    invariant subspace corresponding to the selected eigenvalues.   

    A matrix is in real Schur form if it is upper quasi-triangular with   
    1-by-1 and 2-by-2 blocks. 2-by-2 blocks will be standardized in the   
    form   
            [  a  b  ]   
            [  c  a  ]   

    where b*c < 0. The eigenvalues of such a block are a +- sqrt(bc).   

    Arguments   
    =========   

    JOBVS   (input) CHARACTER*1   
            = 'N': Schur vectors are not computed;   
            = 'V': Schur vectors are computed.   

    SORT    (input) CHARACTER*1   
            Specifies whether or not to order the eigenvalues on the   
            diagonal of the Schur form.   
            = 'N': Eigenvalues are not ordered;   
            = 'S': Eigenvalues are ordered (see SELECT).   

    SELECT  (input) LOGICAL FUNCTION of two REAL arguments   
            SELECT must be declared EXTERNAL in the calling subroutine.   
            If SORT = 'S', SELECT is used to select eigenvalues to sort   
            to the top left of the Schur form.   
            If SORT = 'N', SELECT is not referenced.   
            An eigenvalue WR(j)+sqrt(-1)*WI(j) is selected if   
            SELECT(WR(j),WI(j)) is true; i.e., if either one of a complex 
  
            conjugate pair of eigenvalues is selected, then both complex 
  
            eigenvalues are selected.   
            Note that a selected complex eigenvalue may no longer   
            satisfy SELECT(WR(j),WI(j)) = .TRUE. after ordering, since   
            ordering may change the value of complex eigenvalues   
            (especially if the eigenvalue is ill-conditioned); in this   
            case INFO is set to N+2 (see INFO below).   

    N       (input) INTEGER   
            The order of the matrix A. N >= 0.   

    A       (input/output) REAL array, dimension (LDA,N)   
            On entry, the N-by-N matrix A.   
            On exit, A has been overwritten by its real Schur form T.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(1,N).   

    SDIM    (output) INTEGER   
            If SORT = 'N', SDIM = 0.   
            If SORT = 'S', SDIM = number of eigenvalues (after sorting)   
                           for which SELECT is true. (Complex conjugate   
                           pairs for which SELECT is true for either   
                           eigenvalue count as 2.)   

    WR      (output) REAL array, dimension (N)   
    WI      (output) REAL array, dimension (N)   
            WR and WI contain the real and imaginary parts,   
            respectively, of the computed eigenvalues in the same order   
            that they appear on the diagonal of the output Schur form T. 
  
            Complex conjugate pairs of eigenvalues will appear   
            consecutively with the eigenvalue having the positive   
            imaginary part first.   

    VS      (output) REAL array, dimension (LDVS,N)   
            If JOBVS = 'V', VS contains the orthogonal matrix Z of Schur 
  
            vectors.   
            If JOBVS = 'N', VS is not referenced.   

    LDVS    (input) INTEGER   
            The leading dimension of the array VS.  LDVS >= 1; if   
            JOBVS = 'V', LDVS >= N.   

    WORK    (workspace/output) REAL array, dimension (LWORK)   
            On exit, if INFO = 0, WORK(1) contains the optimal LWORK.   

    LWORK   (input) INTEGER   
            The dimension of the array WORK.  LWORK >= max(1,3*N).   
            For good performance, LWORK must generally be larger.   

    BWORK   (workspace) LOGICAL array, dimension (N)   
            Not referenced if SORT = 'N'.   

    INFO    (output) INTEGER   
            = 0: successful exit   
            < 0: if INFO = -i, the i-th argument had an illegal value.   
            > 0: if INFO = i, and i is   
               <= N: the QR algorithm failed to compute all the   
                     eigenvalues; elements 1:ILO-1 and i+1:N of WR and WI 
  
                     contain those eigenvalues which have converged; if   
                     JOBVS = 'V', VS contains the matrix which reduces A 
  
                     to its partially converged Schur form.   
               = N+1: the eigenvalues could not be reordered because some 
  
                     eigenvalues were too close to separate (the problem 
  
                     is very ill-conditioned);   
               = N+2: after reordering, roundoff changed values of some   
                     complex eigenvalues so that leading eigenvalues in   
                     the Schur form no longer satisfy SELECT=.TRUE.  This 
  
                     could also be caused by underflow due to scaling.   

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
    static integer idum[1], ierr, itau, iwrk, inxt, i, k;
    static real s;
    static integer icond, ieval;
    extern logical lsame_(char *, char *);
    static logical cursl;
    static integer i1, i2;
    extern /* Subroutine */ int scopy_(integer *, real *, integer *, real *, 
	    integer *), sswap_(integer *, real *, integer *, real *, integer *
	    );
    static logical lst2sl;
    extern /* Subroutine */ int slabad_(real *, real *);
    static logical scalea;
    static integer ip;
    static real cscale;
    extern /* Subroutine */ int sgebak_(char *, char *, integer *, integer *, 
	    integer *, real *, integer *, real *, integer *, integer *), sgebal_(char *, integer *, real *, integer *, 
	    integer *, integer *, real *, integer *);
    extern doublereal slamch_(char *), slange_(char *, integer *, 
	    integer *, real *, integer *, real *);
    extern /* Subroutine */ int sgehrd_(integer *, integer *, integer *, real 
	    *, integer *, real *, real *, integer *, integer *), xerbla_(char 
	    *, integer *);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static real bignum;
    extern /* Subroutine */ int slascl_(char *, integer *, integer *, real *, 
	    real *, integer *, integer *, real *, integer *, integer *), slacpy_(char *, integer *, integer *, real *, integer *, 
	    real *, integer *);
    static logical lastsl;
    extern /* Subroutine */ int sorghr_(integer *, integer *, integer *, real 
	    *, integer *, real *, real *, integer *, integer *), shseqr_(char 
	    *, char *, integer *, integer *, integer *, real *, integer *, 
	    real *, real *, real *, integer *, real *, integer *, integer *);
    static integer minwrk, maxwrk;
    static real smlnum;
    static integer hswork;
    extern /* Subroutine */ int strsen_(char *, char *, logical *, integer *, 
	    real *, integer *, real *, integer *, real *, real *, integer *, 
	    real *, real *, real *, integer *, integer *, integer *, integer *
	    );
    static logical wantst, wantvs;
    static integer ihi, ilo;
    static real dum[1], eps, sep;



#define IDUM(I) idum[(I)]
#define WR(I) wr[(I)-1]
#define WI(I) wi[(I)-1]
#define WORK(I) work[(I)-1]
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
/* Computing MAX */
	i__1 = 1, i__2 = *n * 3;
	minwrk = max(i__1,i__2);
	if (! wantvs) {
/* Computing MAX */
	    i__1 = ilaenv_(&c__8, "SHSEQR", "SN", n, &c__1, n, &c_n1, 6L, 2L);
	    maxb = max(i__1,2);
/* Computing MIN   
   Computing MAX */
	    i__3 = 2, i__4 = ilaenv_(&c__4, "SHSEQR", "SN", n, &c__1, n, &
		    c_n1, 6L, 2L);
	    i__1 = min(maxb,*n), i__2 = max(i__3,i__4);
	    k = min(i__1,i__2);
/* Computing MAX */
	    i__1 = k * (k + 2), i__2 = *n << 1;
	    hswork = max(i__1,i__2);
/* Computing MAX */
	    i__1 = maxwrk, i__2 = *n + hswork, i__1 = max(i__1,i__2);
	    maxwrk = max(i__1,1);
	} else {
/* Computing MAX */
	    i__1 = maxwrk, i__2 = (*n << 1) + (*n - 1) * ilaenv_(&c__1, "SOR"
		    "GHR", " ", n, &c__1, n, &c_n1, 6L, 1L);
	    maxwrk = max(i__1,i__2);
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
	    i__1 = maxwrk, i__2 = *n + hswork, i__1 = max(i__1,i__2);
	    maxwrk = max(i__1,1);
	}
	WORK(1) = (real) maxwrk;
    }
    if (*lwork < minwrk) {
	*info = -13;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("SGEES ", &i__1);
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

/*     Permute the matrix to make it more nearly triangular   
       (Workspace: need N) */

    ibal = 1;
    sgebal_("P", n, &A(1,1), lda, &ilo, &ihi, &WORK(ibal), &ierr);

/*     Reduce to upper Hessenberg form   
       (Workspace: need 3*N, prefer 2*N+N*NB) */

    itau = *n + ibal;
    iwrk = *n + itau;
    i__1 = *lwork - iwrk + 1;
    sgehrd_(n, &ilo, &ihi, &A(1,1), lda, &WORK(itau), &WORK(iwrk), &i__1,
	     &ierr);

    if (wantvs) {

/*        Copy Householder vectors to VS */

	slacpy_("L", n, n, &A(1,1), lda, &VS(1,1), ldvs);

/*        Generate orthogonal matrix in VS   
          (Workspace: need 3*N-1, prefer 2*N+(N-1)*NB) */

	i__1 = *lwork - iwrk + 1;
	sorghr_(n, &ilo, &ihi, &VS(1,1), ldvs, &WORK(itau), &WORK(iwrk),
		 &i__1, &ierr);
    }

    *sdim = 0;

/*     Perform QR iteration, accumulating Schur vectors in VS if desired 
  
       (Workspace: need N+1, prefer N+HSWORK (see comments) ) */

    iwrk = itau;
    i__1 = *lwork - iwrk + 1;
    shseqr_("S", jobvs, n, &ilo, &ihi, &A(1,1), lda, &WR(1), &WI(1), &VS(1,1), ldvs, &WORK(iwrk), &i__1, &ieval);
    if (ieval > 0) {
	*info = ieval;
    }

/*     Sort eigenvalues if desired */

    if (wantst && *info == 0) {
	if (scalea) {
	    slascl_("G", &c__0, &c__0, &cscale, &anrm, n, &c__1, &WR(1), n, &
		    ierr);
	    slascl_("G", &c__0, &c__0, &cscale, &anrm, n, &c__1, &WI(1), n, &
		    ierr);
	}
	i__1 = *n;
	for (i = 1; i <= *n; ++i) {
	    BWORK(i) = (*select)(&WR(i), &WI(i));
/* L10: */
	}

/*        Reorder eigenvalues and transform Schur vectors   
          (Workspace: none needed) */

	i__1 = *lwork - iwrk + 1;
	strsen_("N", jobvs, &BWORK(1), n, &A(1,1), lda, &VS(1,1), 
		ldvs, &WR(1), &WI(1), sdim, &s, &sep, &WORK(iwrk), &i__1, 
		idum, &c__1, &icond);
	if (icond > 0) {
	    *info = *n + icond;
	}
    }

    if (wantvs) {

/*        Undo balancing   
          (Workspace: need N) */

	sgebak_("P", "R", n, &ilo, &ihi, &WORK(ibal), n, &VS(1,1), ldvs,
		 &ierr);
    }

    if (scalea) {

/*        Undo scaling for the Schur form of A */

	slascl_("H", &c__0, &c__0, &cscale, &anrm, n, n, &A(1,1), lda, &
		ierr);
	i__1 = *lda + 1;
	scopy_(n, &A(1,1), &i__1, &WR(1), &c__1);
	if (cscale == smlnum) {

/*           If scaling back towards underflow, adjust WI if an   
             offdiagonal element of a 2-by-2 block in the Schur fo
rm   
             underflows. */

	    if (ieval > 0) {
		i1 = ieval + 1;
		i2 = ihi - 1;
		i__1 = ilo - 1;
/* Computing MAX */
		i__3 = ilo - 1;
		i__2 = max(i__3,1);
		slascl_("G", &c__0, &c__0, &cscale, &anrm, &i__1, &c__1, &WI(
			1), &i__2, &ierr);
	    } else if (wantst) {
		i1 = 1;
		i2 = *n - 1;
	    } else {
		i1 = ilo;
		i2 = ihi - 1;
	    }
	    inxt = i1 - 1;
	    i__1 = i2;
	    for (i = i1; i <= i2; ++i) {
		if (i < inxt) {
		    goto L20;
		}
		if (WI(i) == 0.f) {
		    inxt = i + 1;
		} else {
		    if (A(i+1,i) == 0.f) {
			WI(i) = 0.f;
			WI(i + 1) = 0.f;
		    } else if (A(i+1,i) != 0.f && A(i,i+1) == 0.f) {
			WI(i) = 0.f;
			WI(i + 1) = 0.f;
			if (i > 1) {
			    i__2 = i - 1;
			    sswap_(&i__2, &A(1,i), &c__1, &A(1,i+1), &c__1);
			}
			if (*n > i + 1) {
			    i__2 = *n - i - 1;
			    sswap_(&i__2, &A(i,i+2), lda, &A(i+1,i+2), lda);
			}
			sswap_(n, &VS(1,i), &c__1, &VS(1,i+1), &c__1);
			A(i,i+1) = A(i+1,i);
			A(i+1,i) = 0.f;
		    }
		    inxt = i + 2;
		}
L20:
		;
	    }
	}

/*        Undo scaling for the imaginary part of the eigenvalues */

	i__1 = *n - ieval;
/* Computing MAX */
	i__3 = *n - ieval;
	i__2 = max(i__3,1);
	slascl_("G", &c__0, &c__0, &cscale, &anrm, &i__1, &c__1, &WI(ieval + 
		1), &i__2, &ierr);
    }

    if (wantst && *info == 0) {

/*        Check if reordering successful */

	lastsl = TRUE_;
	lst2sl = TRUE_;
	*sdim = 0;
	ip = 0;
	i__1 = *n;
	for (i = 1; i <= *n; ++i) {
	    cursl = (*select)(&WR(i), &WI(i));
	    if (WI(i) == 0.f) {
		if (cursl) {
		    ++(*sdim);
		}
		ip = 0;
		if (cursl && ! lastsl) {
		    *info = *n + 2;
		}
	    } else {
		if (ip == 1) {

/*                 Last eigenvalue of conjugate pair */

		    cursl = cursl || lastsl;
		    lastsl = cursl;
		    if (cursl) {
			*sdim += 2;
		    }
		    ip = -1;
		    if (cursl && ! lst2sl) {
			*info = *n + 2;
		    }
		} else {

/*                 First eigenvalue of conjugate pair */

		    ip = 1;
		}
	    }
	    lst2sl = lastsl;
	    lastsl = cursl;
/* L30: */
	}
    }

    WORK(1) = (real) maxwrk;
    return 0;

/*     End of SGEES */

} /* sgees_ */

