#include "f2c.h"

/* Subroutine */ int zhgeqz_(char *job, char *compq, char *compz, integer *n, 
	integer *ilo, integer *ihi, doublecomplex *a, integer *lda, 
	doublecomplex *b, integer *ldb, doublecomplex *alpha, doublecomplex *
	beta, doublecomplex *q, integer *ldq, doublecomplex *z, integer *ldz, 
	doublecomplex *work, integer *lwork, doublereal *rwork, integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    ZHGEQZ implements a single-shift version of the QZ   
    method for finding the generalized eigenvalues w(i)=ALPHA(i)/BETA(i) 
  
    of the equation   

         det( A - w(i) B ) = 0   

    If JOB='S', then the pair (A,B) is simultaneously   
    reduced to Schur form (i.e., A and B are both upper triangular) by   
    applying one unitary tranformation (usually called Q) on the left and 
  
    another (usually called Z) on the right.  The diagonal elements of   
    A are then ALPHA(1),...,ALPHA(N), and of B are BETA(1),...,BETA(N).   

    If JOB='S' and COMPQ and COMPZ are 'V' or 'I', then the unitary   
    transformations used to reduce (A,B) are accumulated into the arrays 
  
    Q and Z s.t.:   

         Q(in) A(in) Z(in)* = Q(out) A(out) Z(out)*   
         Q(in) B(in) Z(in)* = Q(out) B(out) Z(out)*   

    Ref: C.B. Moler & G.W. Stewart, "An Algorithm for Generalized VISMatrix 
  
         Eigenvalue Problems", SIAM J. Numer. Anal., 10(1973),   
         pp. 241--256.   

    Arguments   
    =========   

    JOB     (input) CHARACTER*1   
            = 'E': compute only ALPHA and BETA.  A and B will not   
                   necessarily be put into generalized Schur form.   
            = 'S': put A and B into generalized Schur form, as well   
                   as computing ALPHA and BETA.   

    COMPQ   (input) CHARACTER*1   
            = 'N': do not modify Q.   
            = 'V': multiply the array Q on the right by the conjugate   
                   transpose of the unitary tranformation that is   
                   applied to the left side of A and B to reduce them   
                   to Schur form.   
            = 'I': like COMPQ='V', except that Q will be initialized to   
                   the identity first.   

    COMPZ   (input) CHARACTER*1   
            = 'N': do not modify Z.   
            = 'V': multiply the array Z on the right by the unitary   
                   tranformation that is applied to the right side of   
                   A and B to reduce them to Schur form.   
            = 'I': like COMPZ='V', except that Z will be initialized to   
                   the identity first.   

    N       (input) INTEGER   
            The order of the matrices A, B, Q, and Z.  N >= 0.   

    ILO     (input) INTEGER   
    IHI     (input) INTEGER   
            It is assumed that A is already upper triangular in rows and 
  
            columns 1:ILO-1 and IHI+1:N.   
            1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.   

    A       (input/output) COMPLEX*16 array, dimension (LDA, N)   
            On entry, the N-by-N upper Hessenberg matrix A.  Elements   
            below the subdiagonal must be zero.   
            If JOB='S', then on exit A and B will have been   
               simultaneously reduced to upper triangular form.   
            If JOB='E', then on exit A will have been destroyed.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max( 1, N ).   

    B       (input/output) COMPLEX*16 array, dimension (LDB, N)   
            On entry, the N-by-N upper triangular matrix B.  Elements   
            below the diagonal must be zero.   
            If JOB='S', then on exit A and B will have been   
               simultaneously reduced to upper triangular form.   
            If JOB='E', then on exit B will have been destroyed.   

    LDB     (input) INTEGER   
            The leading dimension of the array B.  LDB >= max( 1, N ).   

    ALPHA   (output) COMPLEX*16 array, dimension (N)   
            The diagonal elements of A when the pair (A,B) has been   
            reduced to Schur form.  ALPHA(i)/BETA(i) i=1,...,N   
            are the generalized eigenvalues.   

    BETA    (output) COMPLEX*16 array, dimension (N)   
            The diagonal elements of B when the pair (A,B) has been   
            reduced to Schur form.  ALPHA(i)/BETA(i) i=1,...,N   
            are the generalized eigenvalues.  A and B are normalized   
            so that BETA(1),...,BETA(N) are non-negative real numbers.   

    Q       (input/output) COMPLEX*16 array, dimension (LDQ, N)   
            If COMPQ='N', then Q will not be referenced.   
            If COMPQ='V' or 'I', then the conjugate transpose of the   
               unitary transformations which are applied to A and B on   
               the left will be applied to the array Q on the right.   

    LDQ     (input) INTEGER   
            The leading dimension of the array Q.  LDQ >= 1.   
            If COMPQ='V' or 'I', then LDQ >= N.   

    Z       (input/output) COMPLEX*16 array, dimension (LDZ, N)   
            If COMPZ='N', then Z will not be referenced.   
            If COMPZ='V' or 'I', then the unitary transformations which   
               are applied to A and B on the right will be applied to the 
  
               array Z on the right.   

    LDZ     (input) INTEGER   
            The leading dimension of the array Z.  LDZ >= 1.   
            If COMPZ='V' or 'I', then LDZ >= N.   

    WORK    (workspace/output) COMPLEX*16 array, dimension (LWORK)   
            On exit, if INFO >= 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INTEGER   
            The dimension of the array WORK.  LWORK >= max(1,N).   

    RWORK   (workspace) DOUBLE PRECISION array, dimension (N)   

    INFO    (output) INTEGER   
            = 0: successful exit   
            < 0: if INFO = -i, the i-th argument had an illegal value   
            = 1,...,N: the QZ iteration did not converge.  (A,B) is not   
                       in Schur form, but ALPHA(i) and BETA(i),   
                       i=INFO+1,...,N should be correct.   
            = N+1,...,2*N: the shift calculation failed.  (A,B) is not   
                       in Schur form, but ALPHA(i) and BETA(i),   
                       i=INFO-N+1,...,N should be correct.   
            > 2*N:     various "impossible" errors.   

    Further Details   
    ===============   

    We assume that complex ABS works as long as its value is less than   
    overflow.   

    ===================================================================== 
  


       Decode JOB, COMPQ, COMPZ   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static doublecomplex c_b1 = {0.,0.};
    static doublecomplex c_b2 = {1.,0.};
    static integer c__1 = 1;
    static integer c__2 = 2;
    
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, q_dim1, q_offset, z_dim1, 
	    z_offset, i__1, i__2, i__3, i__4, i__5, i__6;
    doublereal d__1, d__2, d__3, d__4, d__5, d__6;
    doublecomplex z__1, z__2, z__3, z__4, z__5, z__6;
    /* Builtin functions */
    double z_abs(doublecomplex *);
    void d_cnjg(doublecomplex *, doublecomplex *);
    double d_imag(doublecomplex *);
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *), pow_zi(
	    doublecomplex *, doublecomplex *, integer *), z_sqrt(
	    doublecomplex *, doublecomplex *);
    /* Local variables */
    static doublereal absb, atol, btol, temp;
    extern /* Subroutine */ int zrot_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublereal *, doublecomplex *);
    static doublereal temp2, c;
    static integer j;
    static doublecomplex s, t;
    extern logical lsame_(char *, char *);
    static doublecomplex ctemp;
    static integer iiter, ilast, jiter;
    static doublereal anorm, bnorm;
    static integer maxit;
    static doublecomplex shift;
    extern /* Subroutine */ int zscal_(integer *, doublecomplex *, 
	    doublecomplex *, integer *);
    static doublereal tempr;
    static doublecomplex ctemp2, ctemp3;
    static logical ilazr2;
    static integer jc, in;
    static doublereal ascale, bscale;
    static doublecomplex u12;
    extern doublereal dlamch_(char *);
    static integer jr;
    static doublecomplex signbc;
    static doublereal safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *);
    static doublecomplex eshift;
    static logical ilschr;
    static integer icompq, ilastm;
    static doublecomplex rtdisc;
    static integer ischur;
    extern doublereal zlanhs_(char *, integer *, doublecomplex *, integer *, 
	    doublereal *);
    static logical ilazro;
    static integer icompz, ifirst;
    extern /* Subroutine */ int zlartg_(doublecomplex *, doublecomplex *, 
	    doublereal *, doublecomplex *, doublecomplex *);
    static integer ifrstm;
    extern /* Subroutine */ int zlaset_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, doublecomplex *, integer *);
    static integer istart;
    static doublecomplex ad11, ad12, ad21, ad22;
    static integer jch;
    static logical ilq, ilz;
    static doublereal ulp;
    static doublecomplex abi22;



#define ALPHA(I) alpha[(I)-1]
#define BETA(I) beta[(I)-1]
#define WORK(I) work[(I)-1]
#define RWORK(I) rwork[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]
#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]
#define Q(I,J) q[(I)-1 + ((J)-1)* ( *ldq)]
#define Z(I,J) z[(I)-1 + ((J)-1)* ( *ldz)]

    if (lsame_(job, "E")) {
	ilschr = FALSE_;
	ischur = 1;
    } else if (lsame_(job, "S")) {
	ilschr = TRUE_;
	ischur = 2;
    } else {
	ischur = 0;
    }

    if (lsame_(compq, "N")) {
	ilq = FALSE_;
	icompq = 1;
    } else if (lsame_(compq, "V")) {
	ilq = TRUE_;
	icompq = 2;
    } else if (lsame_(compq, "I")) {
	ilq = TRUE_;
	icompq = 3;
    } else {
	icompq = 0;
    }

    if (lsame_(compz, "N")) {
	ilz = FALSE_;
	icompz = 1;
    } else if (lsame_(compz, "V")) {
	ilz = TRUE_;
	icompz = 2;
    } else if (lsame_(compz, "I")) {
	ilz = TRUE_;
	icompz = 3;
    } else {
	icompz = 0;
    }

/*     Check Argument Values */

    *info = 0;
    if (ischur == 0) {
	*info = -1;
    } else if (icompq == 0) {
	*info = -2;
    } else if (icompz == 0) {
	*info = -3;
    } else if (*n < 0) {
	*info = -4;
    } else if (*ilo < 1) {
	*info = -5;
    } else if (*ihi > *n || *ihi < *ilo - 1) {
	*info = -6;
    } else if (*lda < *n) {
	*info = -8;
    } else if (*ldb < *n) {
	*info = -10;
    } else if (*ldq < 1 || ilq && *ldq < *n) {
	*info = -14;
    } else if (*ldz < 1 || ilz && *ldz < *n) {
	*info = -16;
    } else if (*lwork < max(1,*n)) {
	*info = -18;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("ZHGEQZ", &i__1);
	return 0;
    }

/*     Quick return if possible */

    WORK(1).r = 1., WORK(1).i = 0.;
    if (*n <= 0) {
	return 0;
    }

/*     Initialize Q and Z */

    if (icompq == 3) {
	zlaset_("Full", n, n, &c_b1, &c_b2, &Q(1,1), ldq);
    }
    if (icompz == 3) {
	zlaset_("Full", n, n, &c_b1, &c_b2, &Z(1,1), ldz);
    }

/*     Machine Constants */

    in = *ihi + 1 - *ilo;
    safmin = dlamch_("S");
    ulp = dlamch_("E") * dlamch_("B");
    anorm = zlanhs_("F", &in, &A(*ilo,*ilo), lda, &RWORK(1));
    bnorm = zlanhs_("F", &in, &B(*ilo,*ilo), ldb, &RWORK(1));
/* Computing MAX */
    d__1 = safmin, d__2 = ulp * anorm;
    atol = max(d__1,d__2);
/* Computing MAX */
    d__1 = safmin, d__2 = ulp * bnorm;
    btol = max(d__1,d__2);
    ascale = 1. / max(safmin,anorm);
    bscale = 1. / max(safmin,bnorm);


/*     Set Eigenvalues IHI+1:N */

    i__1 = *n;
    for (j = *ihi + 1; j <= *n; ++j) {
	absb = z_abs(&B(j,j));
	if (absb > safmin) {
	    i__2 = j + j * b_dim1;
	    z__2.r = B(j,j).r / absb, z__2.i = B(j,j).i / absb;
	    d_cnjg(&z__1, &z__2);
	    signbc.r = z__1.r, signbc.i = z__1.i;
	    i__2 = j + j * b_dim1;
	    B(j,j).r = absb, B(j,j).i = 0.;
	    if (ilschr) {
		i__2 = j - 1;
		zscal_(&i__2, &signbc, &B(1,j), &c__1);
		zscal_(&j, &signbc, &A(1,j), &c__1);
	    } else {
		i__2 = j + j * a_dim1;
		i__3 = j + j * a_dim1;
		z__1.r = A(j,j).r * signbc.r - A(j,j).i * signbc.i, z__1.i =
			 A(j,j).r * signbc.i + A(j,j).i * signbc.r;
		A(j,j).r = z__1.r, A(j,j).i = z__1.i;
	    }
	    if (ilz) {
		zscal_(n, &signbc, &Z(1,j), &c__1);
	    }
	} else {
	    i__2 = j + j * b_dim1;
	    B(j,j).r = 0., B(j,j).i = 0.;
	}
	i__2 = j;
	i__3 = j + j * a_dim1;
	ALPHA(j).r = A(j,j).r, ALPHA(j).i = A(j,j).i;
	i__2 = j;
	i__3 = j + j * b_dim1;
	BETA(j).r = B(j,j).r, BETA(j).i = B(j,j).i;
/* L10: */
    }

/*     If IHI < ILO, skip QZ steps */

    if (*ihi < *ilo) {
	goto L190;
    }

/*     MAIN QZ ITERATION LOOP   

       Initialize dynamic indices   

       Eigenvalues ILAST+1:N have been found.   
          Column operations modify rows IFRSTM:whatever   
          Row operations modify columns whatever:ILASTM   

       If only eigenvalues are being computed, then   
          IFRSTM is the row of the last splitting row above row ILAST;   
          this is always at least ILO.   
       IITER counts iterations since the last eigenvalue was found,   
          to tell when to use an extraordinary shift.   
       MAXIT is the maximum number of QZ sweeps allowed. */

    ilast = *ihi;
    if (ilschr) {
	ifrstm = 1;
	ilastm = *n;
    } else {
	ifrstm = *ilo;
	ilastm = *ihi;
    }
    iiter = 0;
    eshift.r = 0., eshift.i = 0.;
    maxit = (*ihi - *ilo + 1) * 30;

    i__1 = maxit;
    for (jiter = 1; jiter <= maxit; ++jiter) {

/*        Check for too many iterations. */

	if (jiter > maxit) {
	    goto L180;
	}

/*        Split the matrix if possible.   

          Two tests:   
             1: A(j,j-1)=0  or  j=ILO   
             2: B(j,j)=0   

          Special case: j=ILAST */

	if (ilast == *ilo) {
	    goto L60;
	} else {
	    i__2 = ilast + (ilast - 1) * a_dim1;
	    if ((d__1 = A(ilast,ilast-1).r, abs(d__1)) + (d__2 = d_imag(&A(ilast,ilast-1)), abs(d__2)) <= atol) {
		i__2 = ilast + (ilast - 1) * a_dim1;
		A(ilast,ilast-1).r = 0., A(ilast,ilast-1).i = 0.;
		goto L60;
	    }
	}

	if (z_abs(&B(ilast,ilast)) <= btol) {
	    i__2 = ilast + ilast * b_dim1;
	    B(ilast,ilast).r = 0., B(ilast,ilast).i = 0.;
	    goto L50;
	}

/*        General case: j<ILAST */

	i__2 = *ilo;
	for (j = ilast - 1; j >= *ilo; --j) {

/*           Test 1: for A(j,j-1)=0 or j=ILO */

	    if (j == *ilo) {
		ilazro = TRUE_;
	    } else {
		i__3 = j + (j - 1) * a_dim1;
		if ((d__1 = A(j,j-1).r, abs(d__1)) + (d__2 = d_imag(&A(j,j-1)), abs(d__2)) <= atol) {
		    i__3 = j + (j - 1) * a_dim1;
		    A(j,j-1).r = 0., A(j,j-1).i = 0.;
		    ilazro = TRUE_;
		} else {
		    ilazro = FALSE_;
		}
	    }

/*           Test 2: for B(j,j)=0 */

	    if (z_abs(&B(j,j)) < btol) {
		i__3 = j + j * b_dim1;
		B(j,j).r = 0., B(j,j).i = 0.;

/*              Test 1a: Check for 2 consecutive small subdiag
onals in A */

		ilazr2 = FALSE_;
		if (! ilazro) {
		    i__3 = j + (j - 1) * a_dim1;
		    i__4 = j + 1 + j * a_dim1;
		    i__5 = j + j * a_dim1;
		    if (((d__1 = A(j,j-1).r, abs(d__1)) + (d__2 = d_imag(&A(j,j-1)), abs(d__2))) * (ascale * ((
			    d__3 = A(j+1,j).r, abs(d__3)) + (d__4 = d_imag(&A(j+1,j)), abs(d__4)))) <= ((d__5 = A(j,j).r, abs(d__5)) + (d__6 = d_imag(&A(j,j)), abs(d__6))) * (ascale * atol)) {
			ilazr2 = TRUE_;
		    }
		}

/*              If both tests pass (1 & 2), i.e., the leading 
diagonal   
                element of B in the block is zero, split a 1x1
 block off   
                at the top. (I.e., at the J-th row/column) The
 leading   
                diagonal element of the remainder can also be 
zero, so   
                this may have to be done repeatedly. */

		if (ilazro || ilazr2) {
		    i__3 = ilast - 1;
		    for (jch = j; jch <= ilast-1; ++jch) {
			i__4 = jch + jch * a_dim1;
			ctemp.r = A(jch,jch).r, ctemp.i = A(jch,jch).i;
			zlartg_(&ctemp, &A(jch+1,jch), &c, &s, &
				A(jch,jch));
			i__4 = jch + 1 + jch * a_dim1;
			A(jch+1,jch).r = 0., A(jch+1,jch).i = 0.;
			i__4 = ilastm - jch;
			zrot_(&i__4, &A(jch,jch+1), lda, &A(jch+1,jch+1), lda, &c, &s);
			i__4 = ilastm - jch;
			zrot_(&i__4, &B(jch,jch+1), ldb, &B(jch+1,jch+1), ldb, &c, &s);
			if (ilq) {
			    d_cnjg(&z__1, &s);
			    zrot_(n, &Q(1,jch), &c__1, &Q(1,jch+1), &c__1, &c, &z__1);
			}
			if (ilazr2) {
			    i__4 = jch + (jch - 1) * a_dim1;
			    i__5 = jch + (jch - 1) * a_dim1;
			    z__1.r = c * A(jch,jch-1).r, z__1.i = c * A(jch,jch-1).i;
			    A(jch,jch-1).r = z__1.r, A(jch,jch-1).i = z__1.i;
			}
			ilazr2 = FALSE_;
			i__4 = jch + 1 + (jch + 1) * b_dim1;
			if ((d__1 = B(jch+1,jch+1).r, abs(d__1)) + (d__2 = d_imag(&B(jch+1,jch+1)), abs(d__2)) >= 
				btol) {
			    if (jch + 1 >= ilast) {
				goto L60;
			    } else {
				ifirst = jch + 1;
				goto L70;
			    }
			}
			i__4 = jch + 1 + (jch + 1) * b_dim1;
			B(jch+1,jch+1).r = 0., B(jch+1,jch+1).i = 0.;
/* L20: */
		    }
		    goto L50;
		} else {

/*                 Only test 2 passed -- chase the zero to
 B(ILAST,ILAST)   
                   Then process as in the case B(ILAST,ILA
ST)=0 */

		    i__3 = ilast - 1;
		    for (jch = j; jch <= ilast-1; ++jch) {
			i__4 = jch + (jch + 1) * b_dim1;
			ctemp.r = B(jch,jch+1).r, ctemp.i = B(jch,jch+1).i;
			zlartg_(&ctemp, &B(jch+1,jch+1), &c, 
				&s, &B(jch,jch+1));
			i__4 = jch + 1 + (jch + 1) * b_dim1;
			B(jch+1,jch+1).r = 0., B(jch+1,jch+1).i = 0.;
			if (jch < ilastm - 1) {
			    i__4 = ilastm - jch - 1;
			    zrot_(&i__4, &B(jch,jch+2), ldb, &
				    B(jch+1,jch+2), ldb, &c, 
				    &s);
			}
			i__4 = ilastm - jch + 2;
			zrot_(&i__4, &A(jch,jch-1), lda, &A(jch+1,jch-1), lda, &c, &s);
			if (ilq) {
			    d_cnjg(&z__1, &s);
			    zrot_(n, &Q(1,jch), &c__1, &Q(1,jch+1), &c__1, &c, &z__1);
			}
			i__4 = jch + 1 + jch * a_dim1;
			ctemp.r = A(jch+1,jch).r, ctemp.i = A(jch+1,jch).i;
			zlartg_(&ctemp, &A(jch+1,jch-1), &c, 
				&s, &A(jch+1,jch));
			i__4 = jch + 1 + (jch - 1) * a_dim1;
			A(jch+1,jch-1).r = 0., A(jch+1,jch-1).i = 0.;
			i__4 = jch + 1 - ifrstm;
			zrot_(&i__4, &A(ifrstm,jch), &c__1, &A(ifrstm,jch-1), &c__1, &c, &s);
			i__4 = jch - ifrstm;
			zrot_(&i__4, &B(ifrstm,jch), &c__1, &B(ifrstm,jch-1), &c__1, &c, &s);
			if (ilz) {
			    zrot_(n, &Z(1,jch), &c__1, &Z(1,jch-1), &c__1, &c, &s);
			}
/* L30: */
		    }
		    goto L50;
		}
	    } else if (ilazro) {

/*              Only test 1 passed -- work on J:ILAST */

		ifirst = j;
		goto L70;
	    }

/*           Neither test passed -- try next J   

   L40: */
	}

/*        (Drop-through is "impossible") */

	*info = (*n << 1) + 1;
	goto L210;

/*        B(ILAST,ILAST)=0 -- clear A(ILAST,ILAST-1) to split off a   
          1x1 block. */

L50:
	i__2 = ilast + ilast * a_dim1;
	ctemp.r = A(ilast,ilast).r, ctemp.i = A(ilast,ilast).i;
	zlartg_(&ctemp, &A(ilast,ilast-1), &c, &s, &A(ilast,ilast));
	i__2 = ilast + (ilast - 1) * a_dim1;
	A(ilast,ilast-1).r = 0., A(ilast,ilast-1).i = 0.;
	i__2 = ilast - ifrstm;
	zrot_(&i__2, &A(ifrstm,ilast), &c__1, &A(ifrstm,ilast-1), &c__1, &c, &s);
	i__2 = ilast - ifrstm;
	zrot_(&i__2, &B(ifrstm,ilast), &c__1, &B(ifrstm,ilast-1), &c__1, &c, &s);
	if (ilz) {
	    zrot_(n, &Z(1,ilast), &c__1, &Z(1,ilast-1), &c__1, &c, &s);
	}

/*        A(ILAST,ILAST-1)=0 -- Standardize B, set ALPHA and BETA */

L60:
	absb = z_abs(&B(ilast,ilast));
	if (absb > safmin) {
	    i__2 = ilast + ilast * b_dim1;
	    z__2.r = B(ilast,ilast).r / absb, z__2.i = B(ilast,ilast).i / absb;
	    d_cnjg(&z__1, &z__2);
	    signbc.r = z__1.r, signbc.i = z__1.i;
	    i__2 = ilast + ilast * b_dim1;
	    B(ilast,ilast).r = absb, B(ilast,ilast).i = 0.;
	    if (ilschr) {
		i__2 = ilast - ifrstm;
		zscal_(&i__2, &signbc, &B(ifrstm,ilast), &c__1);
		i__2 = ilast + 1 - ifrstm;
		zscal_(&i__2, &signbc, &A(ifrstm,ilast), &c__1);
	    } else {
		i__2 = ilast + ilast * a_dim1;
		i__3 = ilast + ilast * a_dim1;
		z__1.r = A(ilast,ilast).r * signbc.r - A(ilast,ilast).i * signbc.i, z__1.i =
			 A(ilast,ilast).r * signbc.i + A(ilast,ilast).i * signbc.r;
		A(ilast,ilast).r = z__1.r, A(ilast,ilast).i = z__1.i;
	    }
	    if (ilz) {
		zscal_(n, &signbc, &Z(1,ilast), &c__1);
	    }
	} else {
	    i__2 = ilast + ilast * b_dim1;
	    B(ilast,ilast).r = 0., B(ilast,ilast).i = 0.;
	}
	i__2 = ilast;
	i__3 = ilast + ilast * a_dim1;
	ALPHA(ilast).r = A(ilast,ilast).r, ALPHA(ilast).i = A(ilast,ilast).i;
	i__2 = ilast;
	i__3 = ilast + ilast * b_dim1;
	BETA(ilast).r = B(ilast,ilast).r, BETA(ilast).i = B(ilast,ilast).i;

/*        Go to next block -- exit if finished. */

	--ilast;
	if (ilast < *ilo) {
	    goto L190;
	}

/*        Reset counters */

	iiter = 0;
	eshift.r = 0., eshift.i = 0.;
	if (! ilschr) {
	    ilastm = ilast;
	    if (ifrstm > ilast) {
		ifrstm = *ilo;
	    }
	}
	goto L160;

/*        QZ step   

          This iteration only involves rows/columns IFIRST:ILAST.  We 
  
          assume IFIRST < ILAST, and that the diagonal of B is non-zer
o. */

L70:
	++iiter;
	if (! ilschr) {
	    ifrstm = ifirst;
	}

/*        Compute the Shift.   

          At this point, IFIRST < ILAST, and the diagonal elements of 
  
          B(IFIRST:ILAST,IFIRST,ILAST) are larger than BTOL (in   
          magnitude) */

	if (iiter / 10 * 10 != iiter) {

/*           The Wilkinson shift (AEP p.512), i.e., the eigenvalue
 of   
             the bottom-right 2x2 block of A inv(B) which is neare
st to   
             the bottom-right element.   

             We factor B as U*D, where U has unit diagonals, and 
  
             compute (A*inv(D))*inv(U). */

	    i__2 = ilast - 1 + ilast * b_dim1;
	    z__2.r = bscale * B(ilast-1,ilast).r, z__2.i = bscale * B(ilast-1,ilast).i;
	    i__3 = ilast + ilast * b_dim1;
	    z__3.r = bscale * B(ilast,ilast).r, z__3.i = bscale * B(ilast,ilast).i;
	    z_div(&z__1, &z__2, &z__3);
	    u12.r = z__1.r, u12.i = z__1.i;
	    i__2 = ilast - 1 + (ilast - 1) * a_dim1;
	    z__2.r = ascale * A(ilast-1,ilast-1).r, z__2.i = ascale * A(ilast-1,ilast-1).i;
	    i__3 = ilast - 1 + (ilast - 1) * b_dim1;
	    z__3.r = bscale * B(ilast-1,ilast-1).r, z__3.i = bscale * B(ilast-1,ilast-1).i;
	    z_div(&z__1, &z__2, &z__3);
	    ad11.r = z__1.r, ad11.i = z__1.i;
	    i__2 = ilast + (ilast - 1) * a_dim1;
	    z__2.r = ascale * A(ilast,ilast-1).r, z__2.i = ascale * A(ilast,ilast-1).i;
	    i__3 = ilast - 1 + (ilast - 1) * b_dim1;
	    z__3.r = bscale * B(ilast-1,ilast-1).r, z__3.i = bscale * B(ilast-1,ilast-1).i;
	    z_div(&z__1, &z__2, &z__3);
	    ad21.r = z__1.r, ad21.i = z__1.i;
	    i__2 = ilast - 1 + ilast * a_dim1;
	    z__2.r = ascale * A(ilast-1,ilast).r, z__2.i = ascale * A(ilast-1,ilast).i;
	    i__3 = ilast + ilast * b_dim1;
	    z__3.r = bscale * B(ilast,ilast).r, z__3.i = bscale * B(ilast,ilast).i;
	    z_div(&z__1, &z__2, &z__3);
	    ad12.r = z__1.r, ad12.i = z__1.i;
	    i__2 = ilast + ilast * a_dim1;
	    z__2.r = ascale * A(ilast,ilast).r, z__2.i = ascale * A(ilast,ilast).i;
	    i__3 = ilast + ilast * b_dim1;
	    z__3.r = bscale * B(ilast,ilast).r, z__3.i = bscale * B(ilast,ilast).i;
	    z_div(&z__1, &z__2, &z__3);
	    ad22.r = z__1.r, ad22.i = z__1.i;
	    z__2.r = u12.r * ad21.r - u12.i * ad21.i, z__2.i = u12.r * ad21.i 
		    + u12.i * ad21.r;
	    z__1.r = ad22.r - z__2.r, z__1.i = ad22.i - z__2.i;
	    abi22.r = z__1.r, abi22.i = z__1.i;

	    z__2.r = ad11.r + abi22.r, z__2.i = ad11.i + abi22.i;
	    z__1.r = z__2.r * .5, z__1.i = z__2.i * .5;
	    t.r = z__1.r, t.i = z__1.i;
	    pow_zi(&z__4, &t, &c__2);
	    z__5.r = ad12.r * ad21.r - ad12.i * ad21.i, z__5.i = ad12.r * 
		    ad21.i + ad12.i * ad21.r;
	    z__3.r = z__4.r + z__5.r, z__3.i = z__4.i + z__5.i;
	    z__6.r = ad11.r * ad22.r - ad11.i * ad22.i, z__6.i = ad11.r * 
		    ad22.i + ad11.i * ad22.r;
	    z__2.r = z__3.r - z__6.r, z__2.i = z__3.i - z__6.i;
	    z_sqrt(&z__1, &z__2);
	    rtdisc.r = z__1.r, rtdisc.i = z__1.i;
	    z__1.r = t.r - abi22.r, z__1.i = t.i - abi22.i;
	    z__2.r = t.r - abi22.r, z__2.i = t.i - abi22.i;
	    temp = z__1.r * rtdisc.r + d_imag(&z__2) * d_imag(&rtdisc);
	    if (temp <= 0.) {
		z__1.r = t.r + rtdisc.r, z__1.i = t.i + rtdisc.i;
		shift.r = z__1.r, shift.i = z__1.i;
	    } else {
		z__1.r = t.r - rtdisc.r, z__1.i = t.i - rtdisc.i;
		shift.r = z__1.r, shift.i = z__1.i;
	    }
	} else {

/*           Exceptional shift.  Chosen for no particularly good r
eason. */

	    i__2 = ilast - 1 + ilast * a_dim1;
	    z__4.r = ascale * A(ilast-1,ilast).r, z__4.i = ascale * A(ilast-1,ilast).i;
	    i__3 = ilast - 1 + (ilast - 1) * b_dim1;
	    z__5.r = bscale * B(ilast-1,ilast-1).r, z__5.i = bscale * B(ilast-1,ilast-1).i;
	    z_div(&z__3, &z__4, &z__5);
	    d_cnjg(&z__2, &z__3);
	    z__1.r = eshift.r + z__2.r, z__1.i = eshift.i + z__2.i;
	    eshift.r = z__1.r, eshift.i = z__1.i;
	    shift.r = eshift.r, shift.i = eshift.i;
	}

/*        Now check for two consecutive small subdiagonals. */

	i__2 = ifirst + 1;
	for (j = ilast - 1; j >= ifirst+1; --j) {
	    istart = j;
	    i__3 = j + j * a_dim1;
	    z__2.r = ascale * A(j,j).r, z__2.i = ascale * A(j,j).i;
	    i__4 = j + j * b_dim1;
	    z__4.r = bscale * B(j,j).r, z__4.i = bscale * B(j,j).i;
	    z__3.r = shift.r * z__4.r - shift.i * z__4.i, z__3.i = shift.r * 
		    z__4.i + shift.i * z__4.r;
	    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
	    ctemp.r = z__1.r, ctemp.i = z__1.i;
	    temp = (d__1 = ctemp.r, abs(d__1)) + (d__2 = d_imag(&ctemp), abs(
		    d__2));
	    i__3 = j + 1 + j * a_dim1;
	    temp2 = ascale * ((d__1 = A(j+1,j).r, abs(d__1)) + (d__2 = d_imag(&
		    A(j+1,j)), abs(d__2)));
	    tempr = max(temp,temp2);
	    if (tempr < 1. && tempr != 0.) {
		temp /= tempr;
		temp2 /= tempr;
	    }
	    i__3 = j + (j - 1) * a_dim1;
	    if (((d__1 = A(j,j-1).r, abs(d__1)) + (d__2 = d_imag(&A(j,j-1)), abs(d__2))) * temp2 <= temp * atol) {
		goto L90;
	    }
/* L80: */
	}

	istart = ifirst;
	i__2 = ifirst + ifirst * a_dim1;
	z__2.r = ascale * A(ifirst,ifirst).r, z__2.i = ascale * A(ifirst,ifirst).i;
	i__3 = ifirst + ifirst * b_dim1;
	z__4.r = bscale * B(ifirst,ifirst).r, z__4.i = bscale * B(ifirst,ifirst).i;
	z__3.r = shift.r * z__4.r - shift.i * z__4.i, z__3.i = shift.r * 
		z__4.i + shift.i * z__4.r;
	z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
	ctemp.r = z__1.r, ctemp.i = z__1.i;
L90:

/*        Do an implicit-shift QZ sweep.   

          Initial Q */

	i__2 = istart + 1 + istart * a_dim1;
	z__1.r = ascale * A(istart+1,istart).r, z__1.i = ascale * A(istart+1,istart).i;
	ctemp2.r = z__1.r, ctemp2.i = z__1.i;
	zlartg_(&ctemp, &ctemp2, &c, &s, &ctemp3);

/*        Sweep */

	i__2 = ilast - 1;
	for (j = istart; j <= ilast-1; ++j) {
	    if (j > istart) {
		i__3 = j + (j - 1) * a_dim1;
		ctemp.r = A(j,j-1).r, ctemp.i = A(j,j-1).i;
		zlartg_(&ctemp, &A(j+1,j-1), &c, &s, &A(j,j-1));
		i__3 = j + 1 + (j - 1) * a_dim1;
		A(j+1,j-1).r = 0., A(j+1,j-1).i = 0.;
	    }

	    i__3 = ilastm;
	    for (jc = j; jc <= ilastm; ++jc) {
		i__4 = j + jc * a_dim1;
		z__2.r = c * A(j,jc).r, z__2.i = c * A(j,jc).i;
		i__5 = j + 1 + jc * a_dim1;
		z__3.r = s.r * A(j+1,jc).r - s.i * A(j+1,jc).i, z__3.i = s.r * A(j+1,jc).i + s.i * A(j+1,jc).r;
		z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
		ctemp.r = z__1.r, ctemp.i = z__1.i;
		i__4 = j + 1 + jc * a_dim1;
		d_cnjg(&z__4, &s);
		z__3.r = -z__4.r, z__3.i = -z__4.i;
		i__5 = j + jc * a_dim1;
		z__2.r = z__3.r * A(j,jc).r - z__3.i * A(j,jc).i, z__2.i = 
			z__3.r * A(j,jc).i + z__3.i * A(j,jc).r;
		i__6 = j + 1 + jc * a_dim1;
		z__5.r = c * A(j+1,jc).r, z__5.i = c * A(j+1,jc).i;
		z__1.r = z__2.r + z__5.r, z__1.i = z__2.i + z__5.i;
		A(j+1,jc).r = z__1.r, A(j+1,jc).i = z__1.i;
		i__4 = j + jc * a_dim1;
		A(j,jc).r = ctemp.r, A(j,jc).i = ctemp.i;
		i__4 = j + jc * b_dim1;
		z__2.r = c * B(j,jc).r, z__2.i = c * B(j,jc).i;
		i__5 = j + 1 + jc * b_dim1;
		z__3.r = s.r * B(j+1,jc).r - s.i * B(j+1,jc).i, z__3.i = s.r * B(j+1,jc).i + s.i * B(j+1,jc).r;
		z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
		ctemp2.r = z__1.r, ctemp2.i = z__1.i;
		i__4 = j + 1 + jc * b_dim1;
		d_cnjg(&z__4, &s);
		z__3.r = -z__4.r, z__3.i = -z__4.i;
		i__5 = j + jc * b_dim1;
		z__2.r = z__3.r * B(j,jc).r - z__3.i * B(j,jc).i, z__2.i = 
			z__3.r * B(j,jc).i + z__3.i * B(j,jc).r;
		i__6 = j + 1 + jc * b_dim1;
		z__5.r = c * B(j+1,jc).r, z__5.i = c * B(j+1,jc).i;
		z__1.r = z__2.r + z__5.r, z__1.i = z__2.i + z__5.i;
		B(j+1,jc).r = z__1.r, B(j+1,jc).i = z__1.i;
		i__4 = j + jc * b_dim1;
		B(j,jc).r = ctemp2.r, B(j,jc).i = ctemp2.i;
/* L100: */
	    }
	    if (ilq) {
		i__3 = *n;
		for (jr = 1; jr <= *n; ++jr) {
		    i__4 = jr + j * q_dim1;
		    z__2.r = c * Q(jr,j).r, z__2.i = c * Q(jr,j).i;
		    d_cnjg(&z__4, &s);
		    i__5 = jr + (j + 1) * q_dim1;
		    z__3.r = z__4.r * Q(jr,j+1).r - z__4.i * Q(jr,j+1).i, z__3.i =
			     z__4.r * Q(jr,j+1).i + z__4.i * Q(jr,j+1).r;
		    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
		    ctemp.r = z__1.r, ctemp.i = z__1.i;
		    i__4 = jr + (j + 1) * q_dim1;
		    z__3.r = -s.r, z__3.i = -s.i;
		    i__5 = jr + j * q_dim1;
		    z__2.r = z__3.r * Q(jr,j).r - z__3.i * Q(jr,j).i, z__2.i =
			     z__3.r * Q(jr,j).i + z__3.i * Q(jr,j).r;
		    i__6 = jr + (j + 1) * q_dim1;
		    z__4.r = c * Q(jr,j+1).r, z__4.i = c * Q(jr,j+1).i;
		    z__1.r = z__2.r + z__4.r, z__1.i = z__2.i + z__4.i;
		    Q(jr,j+1).r = z__1.r, Q(jr,j+1).i = z__1.i;
		    i__4 = jr + j * q_dim1;
		    Q(jr,j).r = ctemp.r, Q(jr,j).i = ctemp.i;
/* L110: */
		}
	    }

	    i__3 = j + 1 + (j + 1) * b_dim1;
	    ctemp.r = B(j+1,j+1).r, ctemp.i = B(j+1,j+1).i;
	    zlartg_(&ctemp, &B(j+1,j), &c, &s, &B(j+1,j+1));
	    i__3 = j + 1 + j * b_dim1;
	    B(j+1,j).r = 0., B(j+1,j).i = 0.;

/* Computing MIN */
	    i__4 = j + 2;
	    i__3 = min(i__4,ilast);
	    for (jr = ifrstm; jr <= min(j+2,ilast); ++jr) {
		i__4 = jr + (j + 1) * a_dim1;
		z__2.r = c * A(jr,j+1).r, z__2.i = c * A(jr,j+1).i;
		i__5 = jr + j * a_dim1;
		z__3.r = s.r * A(jr,j).r - s.i * A(jr,j).i, z__3.i = s.r * A(jr,j).i + s.i * A(jr,j).r;
		z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
		ctemp.r = z__1.r, ctemp.i = z__1.i;
		i__4 = jr + j * a_dim1;
		d_cnjg(&z__4, &s);
		z__3.r = -z__4.r, z__3.i = -z__4.i;
		i__5 = jr + (j + 1) * a_dim1;
		z__2.r = z__3.r * A(jr,j+1).r - z__3.i * A(jr,j+1).i, z__2.i = 
			z__3.r * A(jr,j+1).i + z__3.i * A(jr,j+1).r;
		i__6 = jr + j * a_dim1;
		z__5.r = c * A(jr,j).r, z__5.i = c * A(jr,j).i;
		z__1.r = z__2.r + z__5.r, z__1.i = z__2.i + z__5.i;
		A(jr,j).r = z__1.r, A(jr,j).i = z__1.i;
		i__4 = jr + (j + 1) * a_dim1;
		A(jr,j+1).r = ctemp.r, A(jr,j+1).i = ctemp.i;
/* L120: */
	    }
	    i__3 = j;
	    for (jr = ifrstm; jr <= j; ++jr) {
		i__4 = jr + (j + 1) * b_dim1;
		z__2.r = c * B(jr,j+1).r, z__2.i = c * B(jr,j+1).i;
		i__5 = jr + j * b_dim1;
		z__3.r = s.r * B(jr,j).r - s.i * B(jr,j).i, z__3.i = s.r * B(jr,j).i + s.i * B(jr,j).r;
		z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
		ctemp.r = z__1.r, ctemp.i = z__1.i;
		i__4 = jr + j * b_dim1;
		d_cnjg(&z__4, &s);
		z__3.r = -z__4.r, z__3.i = -z__4.i;
		i__5 = jr + (j + 1) * b_dim1;
		z__2.r = z__3.r * B(jr,j+1).r - z__3.i * B(jr,j+1).i, z__2.i = 
			z__3.r * B(jr,j+1).i + z__3.i * B(jr,j+1).r;
		i__6 = jr + j * b_dim1;
		z__5.r = c * B(jr,j).r, z__5.i = c * B(jr,j).i;
		z__1.r = z__2.r + z__5.r, z__1.i = z__2.i + z__5.i;
		B(jr,j).r = z__1.r, B(jr,j).i = z__1.i;
		i__4 = jr + (j + 1) * b_dim1;
		B(jr,j+1).r = ctemp.r, B(jr,j+1).i = ctemp.i;
/* L130: */
	    }
	    if (ilz) {
		i__3 = *n;
		for (jr = 1; jr <= *n; ++jr) {
		    i__4 = jr + (j + 1) * z_dim1;
		    z__2.r = c * Z(jr,j+1).r, z__2.i = c * Z(jr,j+1).i;
		    i__5 = jr + j * z_dim1;
		    z__3.r = s.r * Z(jr,j).r - s.i * Z(jr,j).i, z__3.i = s.r *
			     Z(jr,j).i + s.i * Z(jr,j).r;
		    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
		    ctemp.r = z__1.r, ctemp.i = z__1.i;
		    i__4 = jr + j * z_dim1;
		    d_cnjg(&z__4, &s);
		    z__3.r = -z__4.r, z__3.i = -z__4.i;
		    i__5 = jr + (j + 1) * z_dim1;
		    z__2.r = z__3.r * Z(jr,j+1).r - z__3.i * Z(jr,j+1).i, z__2.i =
			     z__3.r * Z(jr,j+1).i + z__3.i * Z(jr,j+1).r;
		    i__6 = jr + j * z_dim1;
		    z__5.r = c * Z(jr,j).r, z__5.i = c * Z(jr,j).i;
		    z__1.r = z__2.r + z__5.r, z__1.i = z__2.i + z__5.i;
		    Z(jr,j).r = z__1.r, Z(jr,j).i = z__1.i;
		    i__4 = jr + (j + 1) * z_dim1;
		    Z(jr,j+1).r = ctemp.r, Z(jr,j+1).i = ctemp.i;
/* L140: */
		}
	    }
/* L150: */
	}

L160:

/* L170: */
	;
    }

/*     Drop-through = non-convergence */

L180:
    *info = ilast;
    goto L210;

/*     Successful completion of all QZ steps */

L190:

/*     Set Eigenvalues 1:ILO-1 */

    i__1 = *ilo - 1;
    for (j = 1; j <= *ilo-1; ++j) {
	absb = z_abs(&B(j,j));
	if (absb > safmin) {
	    i__2 = j + j * b_dim1;
	    z__2.r = B(j,j).r / absb, z__2.i = B(j,j).i / absb;
	    d_cnjg(&z__1, &z__2);
	    signbc.r = z__1.r, signbc.i = z__1.i;
	    i__2 = j + j * b_dim1;
	    B(j,j).r = absb, B(j,j).i = 0.;
	    if (ilschr) {
		i__2 = j - 1;
		zscal_(&i__2, &signbc, &B(1,j), &c__1);
		zscal_(&j, &signbc, &A(1,j), &c__1);
	    } else {
		i__2 = j + j * a_dim1;
		i__3 = j + j * a_dim1;
		z__1.r = A(j,j).r * signbc.r - A(j,j).i * signbc.i, z__1.i =
			 A(j,j).r * signbc.i + A(j,j).i * signbc.r;
		A(j,j).r = z__1.r, A(j,j).i = z__1.i;
	    }
	    if (ilz) {
		zscal_(n, &signbc, &Z(1,j), &c__1);
	    }
	} else {
	    i__2 = j + j * b_dim1;
	    B(j,j).r = 0., B(j,j).i = 0.;
	}
	i__2 = j;
	i__3 = j + j * a_dim1;
	ALPHA(j).r = A(j,j).r, ALPHA(j).i = A(j,j).i;
	i__2 = j;
	i__3 = j + j * b_dim1;
	BETA(j).r = B(j,j).r, BETA(j).i = B(j,j).i;
/* L200: */
    }

/*     Normal Termination */

    *info = 0;

/*     Exit (other than argument error) -- return optimal workspace size 
*/

L210:
    WORK(1).r = 1., WORK(1).i = 0.;
    return 0;

/*     End of ZHGEQZ */

} /* zhgeqz_ */

