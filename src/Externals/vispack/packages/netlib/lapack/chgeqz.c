#include "f2c.h"

/* Subroutine */ int chgeqz_(char *job, char *compq, char *compz, integer *n, 
	integer *ilo, integer *ihi, complex *a, integer *lda, complex *b, 
	integer *ldb, complex *alpha, complex *beta, complex *q, integer *ldq,
	 complex *z, integer *ldz, complex *work, integer *lwork, real *rwork,
	 integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    CHGEQZ implements a single-shift version of the QZ   
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

    A       (input/output) COMPLEX array, dimension (LDA, N)   
            On entry, the N-by-N upper Hessenberg matrix A.  Elements   
            below the subdiagonal must be zero.   
            If JOB='S', then on exit A and B will have been   
               simultaneously reduced to upper triangular form.   
            If JOB='E', then on exit A will have been destroyed.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max( 1, N ).   

    B       (input/output) COMPLEX array, dimension (LDB, N)   
            On entry, the N-by-N upper triangular matrix B.  Elements   
            below the diagonal must be zero.   
            If JOB='S', then on exit A and B will have been   
               simultaneously reduced to upper triangular form.   
            If JOB='E', then on exit B will have been destroyed.   

    LDB     (input) INTEGER   
            The leading dimension of the array B.  LDB >= max( 1, N ).   

    ALPHA   (output) COMPLEX array, dimension (N)   
            The diagonal elements of A when the pair (A,B) has been   
            reduced to Schur form.  ALPHA(i)/BETA(i) i=1,...,N   
            are the generalized eigenvalues.   

    BETA    (output) COMPLEX array, dimension (N)   
            The diagonal elements of B when the pair (A,B) has been   
            reduced to Schur form.  ALPHA(i)/BETA(i) i=1,...,N   
            are the generalized eigenvalues.  A and B are normalized   
            so that BETA(1),...,BETA(N) are non-negative real numbers.   

    Q       (input/output) COMPLEX array, dimension (LDQ, N)   
            If COMPQ='N', then Q will not be referenced.   
            If COMPQ='V' or 'I', then the conjugate transpose of the   
               unitary transformations which are applied to A and B on   
               the left will be applied to the array Q on the right.   

    LDQ     (input) INTEGER   
            The leading dimension of the array Q.  LDQ >= 1.   
            If COMPQ='V' or 'I', then LDQ >= N.   

    Z       (input/output) COMPLEX array, dimension (LDZ, N)   
            If COMPZ='N', then Z will not be referenced.   
            If COMPZ='V' or 'I', then the unitary transformations which   
               are applied to A and B on the right will be applied to the 
  
               array Z on the right.   

    LDZ     (input) INTEGER   
            The leading dimension of the array Z.  LDZ >= 1.   
            If COMPZ='V' or 'I', then LDZ >= N.   

    WORK    (workspace/output) COMPLEX array, dimension (LWORK)   
            On exit, if INFO >= 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INTEGER   
            The dimension of the array WORK.  LWORK >= max(1,N).   

    RWORK   (workspace) REAL array, dimension (N)   

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
    static complex c_b1 = {0.f,0.f};
    static complex c_b2 = {1.f,0.f};
    static integer c__1 = 1;
    static integer c__2 = 2;
    
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, q_dim1, q_offset, z_dim1, 
	    z_offset, i__1, i__2, i__3, i__4, i__5, i__6;
    real r__1, r__2, r__3, r__4, r__5, r__6;
    complex q__1, q__2, q__3, q__4, q__5, q__6;
    /* Builtin functions */
    double c_abs(complex *);
    void r_cnjg(complex *, complex *);
    double r_imag(complex *);
    void c_div(complex *, complex *, complex *), pow_ci(complex *, complex *, 
	    integer *), c_sqrt(complex *, complex *);
    /* Local variables */
    static real absb, atol, btol, temp;
    extern /* Subroutine */ int crot_(integer *, complex *, integer *, 
	    complex *, integer *, real *, complex *);
    static real temp2, c;
    static integer j;
    static complex s, t;
    extern /* Subroutine */ int cscal_(integer *, complex *, complex *, 
	    integer *);
    extern logical lsame_(char *, char *);
    static complex ctemp;
    static integer iiter, ilast, jiter;
    static real anorm, bnorm;
    static integer maxit;
    static complex shift;
    static real tempr;
    static complex ctemp2, ctemp3;
    static logical ilazr2;
    static integer jc, in;
    static real ascale, bscale;
    static complex u12;
    static integer jr;
    static complex signbc;
    extern doublereal slamch_(char *), clanhs_(char *, integer *, 
	    complex *, integer *, real *);
    extern /* Subroutine */ int claset_(char *, integer *, integer *, complex 
	    *, complex *, complex *, integer *), clartg_(complex *, 
	    complex *, real *, complex *, complex *);
    static real safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *);
    static complex eshift;
    static logical ilschr;
    static integer icompq, ilastm;
    static complex rtdisc;
    static integer ischur;
    static logical ilazro;
    static integer icompz, ifirst, ifrstm, istart;
    static complex ad11, ad12, ad21, ad22;
    static integer jch;
    static logical ilq, ilz;
    static real ulp;
    static complex abi22;



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
	xerbla_("CHGEQZ", &i__1);
	return 0;
    }

/*     Quick return if possible */

    WORK(1).r = 1.f, WORK(1).i = 0.f;
    if (*n <= 0) {
	return 0;
    }

/*     Initialize Q and Z */

    if (icompq == 3) {
	claset_("Full", n, n, &c_b1, &c_b2, &Q(1,1), ldq);
    }
    if (icompz == 3) {
	claset_("Full", n, n, &c_b1, &c_b2, &Z(1,1), ldz);
    }

/*     Machine Constants */

    in = *ihi + 1 - *ilo;
    safmin = slamch_("S");
    ulp = slamch_("E") * slamch_("B");
    anorm = clanhs_("F", &in, &A(*ilo,*ilo), lda, &RWORK(1));
    bnorm = clanhs_("F", &in, &B(*ilo,*ilo), ldb, &RWORK(1));
/* Computing MAX */
    r__1 = safmin, r__2 = ulp * anorm;
    atol = dmax(r__1,r__2);
/* Computing MAX */
    r__1 = safmin, r__2 = ulp * bnorm;
    btol = dmax(r__1,r__2);
    ascale = 1.f / dmax(safmin,anorm);
    bscale = 1.f / dmax(safmin,bnorm);


/*     Set Eigenvalues IHI+1:N */

    i__1 = *n;
    for (j = *ihi + 1; j <= *n; ++j) {
	absb = c_abs(&B(j,j));
	if (absb > safmin) {
	    i__2 = j + j * b_dim1;
	    q__2.r = B(j,j).r / absb, q__2.i = B(j,j).i / absb;
	    r_cnjg(&q__1, &q__2);
	    signbc.r = q__1.r, signbc.i = q__1.i;
	    i__2 = j + j * b_dim1;
	    B(j,j).r = absb, B(j,j).i = 0.f;
	    if (ilschr) {
		i__2 = j - 1;
		cscal_(&i__2, &signbc, &B(1,j), &c__1);
		cscal_(&j, &signbc, &A(1,j), &c__1);
	    } else {
		i__2 = j + j * a_dim1;
		i__3 = j + j * a_dim1;
		q__1.r = A(j,j).r * signbc.r - A(j,j).i * signbc.i, q__1.i =
			 A(j,j).r * signbc.i + A(j,j).i * signbc.r;
		A(j,j).r = q__1.r, A(j,j).i = q__1.i;
	    }
	    if (ilz) {
		cscal_(n, &signbc, &Z(1,j), &c__1);
	    }
	} else {
	    i__2 = j + j * b_dim1;
	    B(j,j).r = 0.f, B(j,j).i = 0.f;
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
    eshift.r = 0.f, eshift.i = 0.f;
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
	    if ((r__1 = A(ilast,ilast-1).r, dabs(r__1)) + (r__2 = r_imag(&A(ilast,ilast-1)), dabs(r__2)) <= atol) {
		i__2 = ilast + (ilast - 1) * a_dim1;
		A(ilast,ilast-1).r = 0.f, A(ilast,ilast-1).i = 0.f;
		goto L60;
	    }
	}

	if (c_abs(&B(ilast,ilast)) <= btol) {
	    i__2 = ilast + ilast * b_dim1;
	    B(ilast,ilast).r = 0.f, B(ilast,ilast).i = 0.f;
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
		if ((r__1 = A(j,j-1).r, dabs(r__1)) + (r__2 = r_imag(&A(j,j-1)), dabs(r__2)) <= atol) {
		    i__3 = j + (j - 1) * a_dim1;
		    A(j,j-1).r = 0.f, A(j,j-1).i = 0.f;
		    ilazro = TRUE_;
		} else {
		    ilazro = FALSE_;
		}
	    }

/*           Test 2: for B(j,j)=0 */

	    if (c_abs(&B(j,j)) < btol) {
		i__3 = j + j * b_dim1;
		B(j,j).r = 0.f, B(j,j).i = 0.f;

/*              Test 1a: Check for 2 consecutive small subdiag
onals in A */

		ilazr2 = FALSE_;
		if (! ilazro) {
		    i__3 = j + (j - 1) * a_dim1;
		    i__4 = j + 1 + j * a_dim1;
		    i__5 = j + j * a_dim1;
		    if (((r__1 = A(j,j-1).r, dabs(r__1)) + (r__2 = r_imag(&A(j,j-1)), dabs(r__2))) * (ascale * ((
			    r__3 = A(j+1,j).r, dabs(r__3)) + (r__4 = r_imag(&A(j+1,j)), dabs(r__4)))) <= ((r__5 = A(j,j).r, dabs(r__5)) + (r__6 = r_imag(&A(j,j)), dabs(r__6))) * (ascale * atol)) {
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
			clartg_(&ctemp, &A(jch+1,jch), &c, &s, &
				A(jch,jch));
			i__4 = jch + 1 + jch * a_dim1;
			A(jch+1,jch).r = 0.f, A(jch+1,jch).i = 0.f;
			i__4 = ilastm - jch;
			crot_(&i__4, &A(jch,jch+1), lda, &A(jch+1,jch+1), lda, &c, &s);
			i__4 = ilastm - jch;
			crot_(&i__4, &B(jch,jch+1), ldb, &B(jch+1,jch+1), ldb, &c, &s);
			if (ilq) {
			    r_cnjg(&q__1, &s);
			    crot_(n, &Q(1,jch), &c__1, &Q(1,jch+1), &c__1, &c, &q__1);
			}
			if (ilazr2) {
			    i__4 = jch + (jch - 1) * a_dim1;
			    i__5 = jch + (jch - 1) * a_dim1;
			    q__1.r = c * A(jch,jch-1).r, q__1.i = c * A(jch,jch-1).i;
			    A(jch,jch-1).r = q__1.r, A(jch,jch-1).i = q__1.i;
			}
			ilazr2 = FALSE_;
			i__4 = jch + 1 + (jch + 1) * b_dim1;
			if ((r__1 = B(jch+1,jch+1).r, dabs(r__1)) + (r__2 = r_imag(&
				B(jch+1,jch+1)), dabs(r__2)) 
				>= btol) {
			    if (jch + 1 >= ilast) {
				goto L60;
			    } else {
				ifirst = jch + 1;
				goto L70;
			    }
			}
			i__4 = jch + 1 + (jch + 1) * b_dim1;
			B(jch+1,jch+1).r = 0.f, B(jch+1,jch+1).i = 0.f;
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
			clartg_(&ctemp, &B(jch+1,jch+1), &c, 
				&s, &B(jch,jch+1));
			i__4 = jch + 1 + (jch + 1) * b_dim1;
			B(jch+1,jch+1).r = 0.f, B(jch+1,jch+1).i = 0.f;
			if (jch < ilastm - 1) {
			    i__4 = ilastm - jch - 1;
			    crot_(&i__4, &B(jch,jch+2), ldb, &
				    B(jch+1,jch+2), ldb, &c, 
				    &s);
			}
			i__4 = ilastm - jch + 2;
			crot_(&i__4, &A(jch,jch-1), lda, &A(jch+1,jch-1), lda, &c, &s);
			if (ilq) {
			    r_cnjg(&q__1, &s);
			    crot_(n, &Q(1,jch), &c__1, &Q(1,jch+1), &c__1, &c, &q__1);
			}
			i__4 = jch + 1 + jch * a_dim1;
			ctemp.r = A(jch+1,jch).r, ctemp.i = A(jch+1,jch).i;
			clartg_(&ctemp, &A(jch+1,jch-1), &c, 
				&s, &A(jch+1,jch));
			i__4 = jch + 1 + (jch - 1) * a_dim1;
			A(jch+1,jch-1).r = 0.f, A(jch+1,jch-1).i = 0.f;
			i__4 = jch + 1 - ifrstm;
			crot_(&i__4, &A(ifrstm,jch), &c__1, &A(ifrstm,jch-1), &c__1, &c, &s);
			i__4 = jch - ifrstm;
			crot_(&i__4, &B(ifrstm,jch), &c__1, &B(ifrstm,jch-1), &c__1, &c, &s);
			if (ilz) {
			    crot_(n, &Z(1,jch), &c__1, &Z(1,jch-1), &c__1, &c, &s);
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
	clartg_(&ctemp, &A(ilast,ilast-1), &c, &s, &A(ilast,ilast));
	i__2 = ilast + (ilast - 1) * a_dim1;
	A(ilast,ilast-1).r = 0.f, A(ilast,ilast-1).i = 0.f;
	i__2 = ilast - ifrstm;
	crot_(&i__2, &A(ifrstm,ilast), &c__1, &A(ifrstm,ilast-1), &c__1, &c, &s);
	i__2 = ilast - ifrstm;
	crot_(&i__2, &B(ifrstm,ilast), &c__1, &B(ifrstm,ilast-1), &c__1, &c, &s);
	if (ilz) {
	    crot_(n, &Z(1,ilast), &c__1, &Z(1,ilast-1), &c__1, &c, &s);
	}

/*        A(ILAST,ILAST-1)=0 -- Standardize B, set ALPHA and BETA */

L60:
	absb = c_abs(&B(ilast,ilast));
	if (absb > safmin) {
	    i__2 = ilast + ilast * b_dim1;
	    q__2.r = B(ilast,ilast).r / absb, q__2.i = B(ilast,ilast).i / absb;
	    r_cnjg(&q__1, &q__2);
	    signbc.r = q__1.r, signbc.i = q__1.i;
	    i__2 = ilast + ilast * b_dim1;
	    B(ilast,ilast).r = absb, B(ilast,ilast).i = 0.f;
	    if (ilschr) {
		i__2 = ilast - ifrstm;
		cscal_(&i__2, &signbc, &B(ifrstm,ilast), &c__1);
		i__2 = ilast + 1 - ifrstm;
		cscal_(&i__2, &signbc, &A(ifrstm,ilast), &c__1);
	    } else {
		i__2 = ilast + ilast * a_dim1;
		i__3 = ilast + ilast * a_dim1;
		q__1.r = A(ilast,ilast).r * signbc.r - A(ilast,ilast).i * signbc.i, q__1.i =
			 A(ilast,ilast).r * signbc.i + A(ilast,ilast).i * signbc.r;
		A(ilast,ilast).r = q__1.r, A(ilast,ilast).i = q__1.i;
	    }
	    if (ilz) {
		cscal_(n, &signbc, &Z(1,ilast), &c__1);
	    }
	} else {
	    i__2 = ilast + ilast * b_dim1;
	    B(ilast,ilast).r = 0.f, B(ilast,ilast).i = 0.f;
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
	eshift.r = 0.f, eshift.i = 0.f;
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
	    q__2.r = bscale * B(ilast-1,ilast).r, q__2.i = bscale * B(ilast-1,ilast).i;
	    i__3 = ilast + ilast * b_dim1;
	    q__3.r = bscale * B(ilast,ilast).r, q__3.i = bscale * B(ilast,ilast).i;
	    c_div(&q__1, &q__2, &q__3);
	    u12.r = q__1.r, u12.i = q__1.i;
	    i__2 = ilast - 1 + (ilast - 1) * a_dim1;
	    q__2.r = ascale * A(ilast-1,ilast-1).r, q__2.i = ascale * A(ilast-1,ilast-1).i;
	    i__3 = ilast - 1 + (ilast - 1) * b_dim1;
	    q__3.r = bscale * B(ilast-1,ilast-1).r, q__3.i = bscale * B(ilast-1,ilast-1).i;
	    c_div(&q__1, &q__2, &q__3);
	    ad11.r = q__1.r, ad11.i = q__1.i;
	    i__2 = ilast + (ilast - 1) * a_dim1;
	    q__2.r = ascale * A(ilast,ilast-1).r, q__2.i = ascale * A(ilast,ilast-1).i;
	    i__3 = ilast - 1 + (ilast - 1) * b_dim1;
	    q__3.r = bscale * B(ilast-1,ilast-1).r, q__3.i = bscale * B(ilast-1,ilast-1).i;
	    c_div(&q__1, &q__2, &q__3);
	    ad21.r = q__1.r, ad21.i = q__1.i;
	    i__2 = ilast - 1 + ilast * a_dim1;
	    q__2.r = ascale * A(ilast-1,ilast).r, q__2.i = ascale * A(ilast-1,ilast).i;
	    i__3 = ilast + ilast * b_dim1;
	    q__3.r = bscale * B(ilast,ilast).r, q__3.i = bscale * B(ilast,ilast).i;
	    c_div(&q__1, &q__2, &q__3);
	    ad12.r = q__1.r, ad12.i = q__1.i;
	    i__2 = ilast + ilast * a_dim1;
	    q__2.r = ascale * A(ilast,ilast).r, q__2.i = ascale * A(ilast,ilast).i;
	    i__3 = ilast + ilast * b_dim1;
	    q__3.r = bscale * B(ilast,ilast).r, q__3.i = bscale * B(ilast,ilast).i;
	    c_div(&q__1, &q__2, &q__3);
	    ad22.r = q__1.r, ad22.i = q__1.i;
	    q__2.r = u12.r * ad21.r - u12.i * ad21.i, q__2.i = u12.r * ad21.i 
		    + u12.i * ad21.r;
	    q__1.r = ad22.r - q__2.r, q__1.i = ad22.i - q__2.i;
	    abi22.r = q__1.r, abi22.i = q__1.i;

	    q__2.r = ad11.r + abi22.r, q__2.i = ad11.i + abi22.i;
	    q__1.r = q__2.r * .5f, q__1.i = q__2.i * .5f;
	    t.r = q__1.r, t.i = q__1.i;
	    pow_ci(&q__4, &t, &c__2);
	    q__5.r = ad12.r * ad21.r - ad12.i * ad21.i, q__5.i = ad12.r * 
		    ad21.i + ad12.i * ad21.r;
	    q__3.r = q__4.r + q__5.r, q__3.i = q__4.i + q__5.i;
	    q__6.r = ad11.r * ad22.r - ad11.i * ad22.i, q__6.i = ad11.r * 
		    ad22.i + ad11.i * ad22.r;
	    q__2.r = q__3.r - q__6.r, q__2.i = q__3.i - q__6.i;
	    c_sqrt(&q__1, &q__2);
	    rtdisc.r = q__1.r, rtdisc.i = q__1.i;
	    q__1.r = t.r - abi22.r, q__1.i = t.i - abi22.i;
	    q__2.r = t.r - abi22.r, q__2.i = t.i - abi22.i;
	    temp = q__1.r * rtdisc.r + r_imag(&q__2) * r_imag(&rtdisc);
	    if (temp <= 0.f) {
		q__1.r = t.r + rtdisc.r, q__1.i = t.i + rtdisc.i;
		shift.r = q__1.r, shift.i = q__1.i;
	    } else {
		q__1.r = t.r - rtdisc.r, q__1.i = t.i - rtdisc.i;
		shift.r = q__1.r, shift.i = q__1.i;
	    }
	} else {

/*           Exceptional shift.  Chosen for no particularly good r
eason. */

	    i__2 = ilast - 1 + ilast * a_dim1;
	    q__4.r = ascale * A(ilast-1,ilast).r, q__4.i = ascale * A(ilast-1,ilast).i;
	    i__3 = ilast - 1 + (ilast - 1) * b_dim1;
	    q__5.r = bscale * B(ilast-1,ilast-1).r, q__5.i = bscale * B(ilast-1,ilast-1).i;
	    c_div(&q__3, &q__4, &q__5);
	    r_cnjg(&q__2, &q__3);
	    q__1.r = eshift.r + q__2.r, q__1.i = eshift.i + q__2.i;
	    eshift.r = q__1.r, eshift.i = q__1.i;
	    shift.r = eshift.r, shift.i = eshift.i;
	}

/*        Now check for two consecutive small subdiagonals. */

	i__2 = ifirst + 1;
	for (j = ilast - 1; j >= ifirst+1; --j) {
	    istart = j;
	    i__3 = j + j * a_dim1;
	    q__2.r = ascale * A(j,j).r, q__2.i = ascale * A(j,j).i;
	    i__4 = j + j * b_dim1;
	    q__4.r = bscale * B(j,j).r, q__4.i = bscale * B(j,j).i;
	    q__3.r = shift.r * q__4.r - shift.i * q__4.i, q__3.i = shift.r * 
		    q__4.i + shift.i * q__4.r;
	    q__1.r = q__2.r - q__3.r, q__1.i = q__2.i - q__3.i;
	    ctemp.r = q__1.r, ctemp.i = q__1.i;
	    temp = (r__1 = ctemp.r, dabs(r__1)) + (r__2 = r_imag(&ctemp), 
		    dabs(r__2));
	    i__3 = j + 1 + j * a_dim1;
	    temp2 = ascale * ((r__1 = A(j+1,j).r, dabs(r__1)) + (r__2 = r_imag(
		    &A(j+1,j)), dabs(r__2)));
	    tempr = dmax(temp,temp2);
	    if (tempr < 1.f && tempr != 0.f) {
		temp /= tempr;
		temp2 /= tempr;
	    }
	    i__3 = j + (j - 1) * a_dim1;
	    if (((r__1 = A(j,j-1).r, dabs(r__1)) + (r__2 = r_imag(&A(j,j-1)), dabs(r__2))) * temp2 <= temp * atol) {
		goto L90;
	    }
/* L80: */
	}

	istart = ifirst;
	i__2 = ifirst + ifirst * a_dim1;
	q__2.r = ascale * A(ifirst,ifirst).r, q__2.i = ascale * A(ifirst,ifirst).i;
	i__3 = ifirst + ifirst * b_dim1;
	q__4.r = bscale * B(ifirst,ifirst).r, q__4.i = bscale * B(ifirst,ifirst).i;
	q__3.r = shift.r * q__4.r - shift.i * q__4.i, q__3.i = shift.r * 
		q__4.i + shift.i * q__4.r;
	q__1.r = q__2.r - q__3.r, q__1.i = q__2.i - q__3.i;
	ctemp.r = q__1.r, ctemp.i = q__1.i;
L90:

/*        Do an implicit-shift QZ sweep.   

          Initial Q */

	i__2 = istart + 1 + istart * a_dim1;
	q__1.r = ascale * A(istart+1,istart).r, q__1.i = ascale * A(istart+1,istart).i;
	ctemp2.r = q__1.r, ctemp2.i = q__1.i;
	clartg_(&ctemp, &ctemp2, &c, &s, &ctemp3);

/*        Sweep */

	i__2 = ilast - 1;
	for (j = istart; j <= ilast-1; ++j) {
	    if (j > istart) {
		i__3 = j + (j - 1) * a_dim1;
		ctemp.r = A(j,j-1).r, ctemp.i = A(j,j-1).i;
		clartg_(&ctemp, &A(j+1,j-1), &c, &s, &A(j,j-1));
		i__3 = j + 1 + (j - 1) * a_dim1;
		A(j+1,j-1).r = 0.f, A(j+1,j-1).i = 0.f;
	    }

	    i__3 = ilastm;
	    for (jc = j; jc <= ilastm; ++jc) {
		i__4 = j + jc * a_dim1;
		q__2.r = c * A(j,jc).r, q__2.i = c * A(j,jc).i;
		i__5 = j + 1 + jc * a_dim1;
		q__3.r = s.r * A(j+1,jc).r - s.i * A(j+1,jc).i, q__3.i = s.r * A(j+1,jc).i + s.i * A(j+1,jc).r;
		q__1.r = q__2.r + q__3.r, q__1.i = q__2.i + q__3.i;
		ctemp.r = q__1.r, ctemp.i = q__1.i;
		i__4 = j + 1 + jc * a_dim1;
		r_cnjg(&q__4, &s);
		q__3.r = -(doublereal)q__4.r, q__3.i = -(doublereal)q__4.i;
		i__5 = j + jc * a_dim1;
		q__2.r = q__3.r * A(j,jc).r - q__3.i * A(j,jc).i, q__2.i = 
			q__3.r * A(j,jc).i + q__3.i * A(j,jc).r;
		i__6 = j + 1 + jc * a_dim1;
		q__5.r = c * A(j+1,jc).r, q__5.i = c * A(j+1,jc).i;
		q__1.r = q__2.r + q__5.r, q__1.i = q__2.i + q__5.i;
		A(j+1,jc).r = q__1.r, A(j+1,jc).i = q__1.i;
		i__4 = j + jc * a_dim1;
		A(j,jc).r = ctemp.r, A(j,jc).i = ctemp.i;
		i__4 = j + jc * b_dim1;
		q__2.r = c * B(j,jc).r, q__2.i = c * B(j,jc).i;
		i__5 = j + 1 + jc * b_dim1;
		q__3.r = s.r * B(j+1,jc).r - s.i * B(j+1,jc).i, q__3.i = s.r * B(j+1,jc).i + s.i * B(j+1,jc).r;
		q__1.r = q__2.r + q__3.r, q__1.i = q__2.i + q__3.i;
		ctemp2.r = q__1.r, ctemp2.i = q__1.i;
		i__4 = j + 1 + jc * b_dim1;
		r_cnjg(&q__4, &s);
		q__3.r = -(doublereal)q__4.r, q__3.i = -(doublereal)q__4.i;
		i__5 = j + jc * b_dim1;
		q__2.r = q__3.r * B(j,jc).r - q__3.i * B(j,jc).i, q__2.i = 
			q__3.r * B(j,jc).i + q__3.i * B(j,jc).r;
		i__6 = j + 1 + jc * b_dim1;
		q__5.r = c * B(j+1,jc).r, q__5.i = c * B(j+1,jc).i;
		q__1.r = q__2.r + q__5.r, q__1.i = q__2.i + q__5.i;
		B(j+1,jc).r = q__1.r, B(j+1,jc).i = q__1.i;
		i__4 = j + jc * b_dim1;
		B(j,jc).r = ctemp2.r, B(j,jc).i = ctemp2.i;
/* L100: */
	    }
	    if (ilq) {
		i__3 = *n;
		for (jr = 1; jr <= *n; ++jr) {
		    i__4 = jr + j * q_dim1;
		    q__2.r = c * Q(jr,j).r, q__2.i = c * Q(jr,j).i;
		    r_cnjg(&q__4, &s);
		    i__5 = jr + (j + 1) * q_dim1;
		    q__3.r = q__4.r * Q(jr,j+1).r - q__4.i * Q(jr,j+1).i, q__3.i =
			     q__4.r * Q(jr,j+1).i + q__4.i * Q(jr,j+1).r;
		    q__1.r = q__2.r + q__3.r, q__1.i = q__2.i + q__3.i;
		    ctemp.r = q__1.r, ctemp.i = q__1.i;
		    i__4 = jr + (j + 1) * q_dim1;
		    q__3.r = -(doublereal)s.r, q__3.i = -(doublereal)s.i;
		    i__5 = jr + j * q_dim1;
		    q__2.r = q__3.r * Q(jr,j).r - q__3.i * Q(jr,j).i, q__2.i =
			     q__3.r * Q(jr,j).i + q__3.i * Q(jr,j).r;
		    i__6 = jr + (j + 1) * q_dim1;
		    q__4.r = c * Q(jr,j+1).r, q__4.i = c * Q(jr,j+1).i;
		    q__1.r = q__2.r + q__4.r, q__1.i = q__2.i + q__4.i;
		    Q(jr,j+1).r = q__1.r, Q(jr,j+1).i = q__1.i;
		    i__4 = jr + j * q_dim1;
		    Q(jr,j).r = ctemp.r, Q(jr,j).i = ctemp.i;
/* L110: */
		}
	    }

	    i__3 = j + 1 + (j + 1) * b_dim1;
	    ctemp.r = B(j+1,j+1).r, ctemp.i = B(j+1,j+1).i;
	    clartg_(&ctemp, &B(j+1,j), &c, &s, &B(j+1,j+1));
	    i__3 = j + 1 + j * b_dim1;
	    B(j+1,j).r = 0.f, B(j+1,j).i = 0.f;

/* Computing MIN */
	    i__4 = j + 2;
	    i__3 = min(i__4,ilast);
	    for (jr = ifrstm; jr <= min(j+2,ilast); ++jr) {
		i__4 = jr + (j + 1) * a_dim1;
		q__2.r = c * A(jr,j+1).r, q__2.i = c * A(jr,j+1).i;
		i__5 = jr + j * a_dim1;
		q__3.r = s.r * A(jr,j).r - s.i * A(jr,j).i, q__3.i = s.r * A(jr,j).i + s.i * A(jr,j).r;
		q__1.r = q__2.r + q__3.r, q__1.i = q__2.i + q__3.i;
		ctemp.r = q__1.r, ctemp.i = q__1.i;
		i__4 = jr + j * a_dim1;
		r_cnjg(&q__4, &s);
		q__3.r = -(doublereal)q__4.r, q__3.i = -(doublereal)q__4.i;
		i__5 = jr + (j + 1) * a_dim1;
		q__2.r = q__3.r * A(jr,j+1).r - q__3.i * A(jr,j+1).i, q__2.i = 
			q__3.r * A(jr,j+1).i + q__3.i * A(jr,j+1).r;
		i__6 = jr + j * a_dim1;
		q__5.r = c * A(jr,j).r, q__5.i = c * A(jr,j).i;
		q__1.r = q__2.r + q__5.r, q__1.i = q__2.i + q__5.i;
		A(jr,j).r = q__1.r, A(jr,j).i = q__1.i;
		i__4 = jr + (j + 1) * a_dim1;
		A(jr,j+1).r = ctemp.r, A(jr,j+1).i = ctemp.i;
/* L120: */
	    }
	    i__3 = j;
	    for (jr = ifrstm; jr <= j; ++jr) {
		i__4 = jr + (j + 1) * b_dim1;
		q__2.r = c * B(jr,j+1).r, q__2.i = c * B(jr,j+1).i;
		i__5 = jr + j * b_dim1;
		q__3.r = s.r * B(jr,j).r - s.i * B(jr,j).i, q__3.i = s.r * B(jr,j).i + s.i * B(jr,j).r;
		q__1.r = q__2.r + q__3.r, q__1.i = q__2.i + q__3.i;
		ctemp.r = q__1.r, ctemp.i = q__1.i;
		i__4 = jr + j * b_dim1;
		r_cnjg(&q__4, &s);
		q__3.r = -(doublereal)q__4.r, q__3.i = -(doublereal)q__4.i;
		i__5 = jr + (j + 1) * b_dim1;
		q__2.r = q__3.r * B(jr,j+1).r - q__3.i * B(jr,j+1).i, q__2.i = 
			q__3.r * B(jr,j+1).i + q__3.i * B(jr,j+1).r;
		i__6 = jr + j * b_dim1;
		q__5.r = c * B(jr,j).r, q__5.i = c * B(jr,j).i;
		q__1.r = q__2.r + q__5.r, q__1.i = q__2.i + q__5.i;
		B(jr,j).r = q__1.r, B(jr,j).i = q__1.i;
		i__4 = jr + (j + 1) * b_dim1;
		B(jr,j+1).r = ctemp.r, B(jr,j+1).i = ctemp.i;
/* L130: */
	    }
	    if (ilz) {
		i__3 = *n;
		for (jr = 1; jr <= *n; ++jr) {
		    i__4 = jr + (j + 1) * z_dim1;
		    q__2.r = c * Z(jr,j+1).r, q__2.i = c * Z(jr,j+1).i;
		    i__5 = jr + j * z_dim1;
		    q__3.r = s.r * Z(jr,j).r - s.i * Z(jr,j).i, q__3.i = s.r *
			     Z(jr,j).i + s.i * Z(jr,j).r;
		    q__1.r = q__2.r + q__3.r, q__1.i = q__2.i + q__3.i;
		    ctemp.r = q__1.r, ctemp.i = q__1.i;
		    i__4 = jr + j * z_dim1;
		    r_cnjg(&q__4, &s);
		    q__3.r = -(doublereal)q__4.r, q__3.i = -(doublereal)
			    q__4.i;
		    i__5 = jr + (j + 1) * z_dim1;
		    q__2.r = q__3.r * Z(jr,j+1).r - q__3.i * Z(jr,j+1).i, q__2.i =
			     q__3.r * Z(jr,j+1).i + q__3.i * Z(jr,j+1).r;
		    i__6 = jr + j * z_dim1;
		    q__5.r = c * Z(jr,j).r, q__5.i = c * Z(jr,j).i;
		    q__1.r = q__2.r + q__5.r, q__1.i = q__2.i + q__5.i;
		    Z(jr,j).r = q__1.r, Z(jr,j).i = q__1.i;
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
	absb = c_abs(&B(j,j));
	if (absb > safmin) {
	    i__2 = j + j * b_dim1;
	    q__2.r = B(j,j).r / absb, q__2.i = B(j,j).i / absb;
	    r_cnjg(&q__1, &q__2);
	    signbc.r = q__1.r, signbc.i = q__1.i;
	    i__2 = j + j * b_dim1;
	    B(j,j).r = absb, B(j,j).i = 0.f;
	    if (ilschr) {
		i__2 = j - 1;
		cscal_(&i__2, &signbc, &B(1,j), &c__1);
		cscal_(&j, &signbc, &A(1,j), &c__1);
	    } else {
		i__2 = j + j * a_dim1;
		i__3 = j + j * a_dim1;
		q__1.r = A(j,j).r * signbc.r - A(j,j).i * signbc.i, q__1.i =
			 A(j,j).r * signbc.i + A(j,j).i * signbc.r;
		A(j,j).r = q__1.r, A(j,j).i = q__1.i;
	    }
	    if (ilz) {
		cscal_(n, &signbc, &Z(1,j), &c__1);
	    }
	} else {
	    i__2 = j + j * b_dim1;
	    B(j,j).r = 0.f, B(j,j).i = 0.f;
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
    WORK(1).r = 1.f, WORK(1).i = 0.f;
    return 0;

/*     End of CHGEQZ */

} /* chgeqz_ */

