#include "f2c.h"

/* Subroutine */ int dtrsyl_(char *trana, char *tranb, integer *isgn, integer 
	*m, integer *n, doublereal *a, integer *lda, doublereal *b, integer *
	ldb, doublereal *c, integer *ldc, doublereal *scale, integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       March 31, 1993   


    Purpose   
    =======   

    DTRSYL solves the real Sylvester matrix equation:   

       op(A)*X + X*op(B) = scale*C or   
       op(A)*X - X*op(B) = scale*C,   

    where op(A) = A or A**T, and  A and B are both upper quasi-   
    triangular. A is M-by-M and B is N-by-N; the right hand side C and   
    the solution X are M-by-N; and scale is an output scale factor, set   
    <= 1 to avoid overflow in X.   

    A and B must be in Schur canonical form (as returned by DHSEQR), that 
  
    is, block upper triangular with 1-by-1 and 2-by-2 diagonal blocks;   
    each 2-by-2 diagonal block has its diagonal elements equal and its   
    off-diagonal elements of opposite sign.   

    Arguments   
    =========   

    TRANA   (input) CHARACTER*1   
            Specifies the option op(A):   
            = 'N': op(A) = A    (No transpose)   
            = 'T': op(A) = A**T (Transpose)   
            = 'C': op(A) = A**H (Conjugate transpose = Transpose)   

    TRANB   (input) CHARACTER*1   
            Specifies the option op(B):   
            = 'N': op(B) = B    (No transpose)   
            = 'T': op(B) = B**T (Transpose)   
            = 'C': op(B) = B**H (Conjugate transpose = Transpose)   

    ISGN    (input) INTEGER   
            Specifies the sign in the equation:   
            = +1: solve op(A)*X + X*op(B) = scale*C   
            = -1: solve op(A)*X - X*op(B) = scale*C   

    M       (input) INTEGER   
            The order of the matrix A, and the number of rows in the   
            matrices X and C. M >= 0.   

    N       (input) INTEGER   
            The order of the matrix B, and the number of columns in the   
            matrices X and C. N >= 0.   

    A       (input) DOUBLE PRECISION array, dimension (LDA,M)   
            The upper quasi-triangular matrix A, in Schur canonical form. 
  

    LDA     (input) INTEGER   
            The leading dimension of the array A. LDA >= max(1,M).   

    B       (input) DOUBLE PRECISION array, dimension (LDB,N)   
            The upper quasi-triangular matrix B, in Schur canonical form. 
  

    LDB     (input) INTEGER   
            The leading dimension of the array B. LDB >= max(1,N).   

    C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)   
            On entry, the M-by-N right hand side matrix C.   
            On exit, C is overwritten by the solution matrix X.   

    LDC     (input) INTEGER   
            The leading dimension of the array C. LDC >= max(1,M)   

    SCALE   (output) DOUBLE PRECISION   
            The scale factor, scale, set <= 1 to avoid overflow in X.   

    INFO    (output) INTEGER   
            = 0: successful exit   
            < 0: if INFO = -i, the i-th argument had an illegal value   
            = 1: A and B have common or very close eigenvalues; perturbed 
  
                 values were used to solve the equation (but the matrices 
  
                 A and B are unchanged).   

    ===================================================================== 
  


       Decode and Test input parameters   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static integer c__1 = 1;
    static logical c_false = FALSE_;
    static integer c__2 = 2;
    static doublereal c_b26 = 1.;
    static doublereal c_b30 = 0.;
    static logical c_true = TRUE_;
    
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, i__1, i__2, 
	    i__3, i__4;
    doublereal d__1, d__2;
    /* Local variables */
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static integer ierr;
    static doublereal smin, suml, sumr;
    static integer j, k, l;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static doublereal x[4]	/* was [2][2] */;
    extern logical lsame_(char *, char *);
    static integer knext, lnext, k1, k2, l1, l2;
    static doublereal xnorm;
    extern /* Subroutine */ int dlaln2_(logical *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     doublereal *, doublereal *, integer *, doublereal *, doublereal *
	    , doublereal *, integer *, doublereal *, doublereal *, integer *),
	     dlasy2_(logical *, logical *, integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *);
    static doublereal a11, db;
    extern /* Subroutine */ int dlabad_(doublereal *, doublereal *);
    extern doublereal dlamch_(char *), dlange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *);
    static doublereal scaloc;
    extern /* Subroutine */ int xerbla_(char *, integer *);
    static doublereal bignum;
    static logical notrna, notrnb;
    static doublereal smlnum, da11, vec[4]	/* was [2][2] */, dum[1], eps,
	     sgn;



#define X(I) x[(I)]
#define WAS(I) was[(I)]
#define DUM(I) dum[(I)]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]
#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]
#define C(I,J) c[(I)-1 + ((J)-1)* ( *ldc)]

    notrna = lsame_(trana, "N");
    notrnb = lsame_(tranb, "N");

    *info = 0;
    if (! notrna && ! lsame_(trana, "T") && ! lsame_(trana, "C")) {
	*info = -1;
    } else if (! notrnb && ! lsame_(tranb, "T") && ! lsame_(tranb, 
	    "C")) {
	*info = -2;
    } else if (*isgn != 1 && *isgn != -1) {
	*info = -3;
    } else if (*m < 0) {
	*info = -4;
    } else if (*n < 0) {
	*info = -5;
    } else if (*lda < max(1,*m)) {
	*info = -7;
    } else if (*ldb < max(1,*n)) {
	*info = -9;
    } else if (*ldc < max(1,*m)) {
	*info = -11;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DTRSYL", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*m == 0 || *n == 0) {
	return 0;
    }

/*     Set constants to control overflow */

    eps = dlamch_("P");
    smlnum = dlamch_("S");
    bignum = 1. / smlnum;
    dlabad_(&smlnum, &bignum);
    smlnum = smlnum * (doublereal) (*m * *n) / eps;
    bignum = 1. / smlnum;

/* Computing MAX */
    d__1 = smlnum, d__2 = eps * dlange_("M", m, m, &A(1,1), lda, dum)
	    , d__1 = max(d__1,d__2), d__2 = eps * dlange_("M", n, n, &B(1,1), ldb, dum);
    smin = max(d__1,d__2);

    *scale = 1.;
    sgn = (doublereal) (*isgn);

    if (notrna && notrnb) {

/*        Solve    A*X + ISGN*X*B = scale*C.   

          The (K,L)th block of X is determined starting from   
          bottom-left corner column by column by   

           A(K,K)*X(K,L) + ISGN*X(K,L)*B(L,L) = C(K,L) - R(K,L)   

          Where   
                    M                         L-1   
          R(K,L) = SUM [A(K,I)*X(I,L)] + ISGN*SUM [X(K,J)*B(J,L)].   
                  I=K+1                       J=1   

          Start column loop (index = L)   
          L1 (L2) : column index of the first (first) row of X(K,L). 
*/

	lnext = 1;
	i__1 = *n;
	for (l = 1; l <= *n; ++l) {
	    if (l < lnext) {
		goto L60;
	    }
	    if (l == *n) {
		l1 = l;
		l2 = l;
	    } else {
		if (B(l+1,l) != 0.) {
		    l1 = l;
		    l2 = l + 1;
		    lnext = l + 2;
		} else {
		    l1 = l;
		    l2 = l;
		    lnext = l + 1;
		}
	    }

/*           Start row loop (index = K)   
             K1 (K2): row index of the first (last) row of X(K,L).
 */

	    knext = *m;
	    for (k = *m; k >= 1; --k) {
		if (k > knext) {
		    goto L50;
		}
		if (k == 1) {
		    k1 = k;
		    k2 = k;
		} else {
		    if (A(k,k-1) != 0.) {
			k1 = k - 1;
			k2 = k;
			knext = k - 2;
		    } else {
			k1 = k;
			k2 = k;
			knext = k - 1;
		    }
		}

		if (l1 == l2 && k1 == k2) {
		    i__2 = *m - k1;
/* Computing MIN */
		    i__3 = k1 + 1;
/* Computing MIN */
		    i__4 = k1 + 1;
		    suml = ddot_(&i__2, &A(k1,min(k1+1,*m)), lda, &
			    C(min(k1+1,*m),l1), &c__1);
		    i__2 = l1 - 1;
		    sumr = ddot_(&i__2, &C(k1,1), ldc, &B(1,l1), &c__1);
		    vec[0] = C(k1,l1) - (suml + sgn * sumr);
		    scaloc = 1.;

		    a11 = A(k1,k1) + sgn * B(l1,l1);
		    da11 = abs(a11);
		    if (da11 <= smin) {
			a11 = smin;
			da11 = smin;
			*info = 1;
		    }
		    db = abs(vec[0]);
		    if (da11 < 1. && db > 1.) {
			if (db > bignum * da11) {
			    scaloc = 1. / db;
			}
		    }
		    X(0) = vec[0] * scaloc / a11;

		    if (scaloc != 1.) {
			i__2 = *n;
			for (j = 1; j <= *n; ++j) {
			    dscal_(m, &scaloc, &C(1,j), &c__1);
/* L10: */
			}
			*scale *= scaloc;
		    }
		    C(k1,l1) = X(0);

		} else if (l1 == l2 && k1 != k2) {

		    i__2 = *m - k2;
/* Computing MIN */
		    i__3 = k2 + 1;
/* Computing MIN */
		    i__4 = k2 + 1;
		    suml = ddot_(&i__2, &A(k1,min(k2+1,*m)), lda, &
			    C(min(k2+1,*m),l1), &c__1);
		    i__2 = l1 - 1;
		    sumr = ddot_(&i__2, &C(k1,1), ldc, &B(1,l1), &c__1);
		    vec[0] = C(k1,l1) - (suml + sgn * sumr);

		    i__2 = *m - k2;
/* Computing MIN */
		    i__3 = k2 + 1;
/* Computing MIN */
		    i__4 = k2 + 1;
		    suml = ddot_(&i__2, &A(k2,min(k2+1,*m)), lda, &
			    C(min(k2+1,*m),l1), &c__1);
		    i__2 = l1 - 1;
		    sumr = ddot_(&i__2, &C(k2,1), ldc, &B(1,l1), &c__1);
		    vec[1] = C(k2,l1) - (suml + sgn * sumr);

		    d__1 = -sgn * B(l1,l1);
		    dlaln2_(&c_false, &c__2, &c__1, &smin, &c_b26, &A(k1,k1), lda, &c_b26, &c_b26, vec, &c__2, &d__1,
			     &c_b30, x, &c__2, &scaloc, &xnorm, &ierr);
		    if (ierr != 0) {
			*info = 1;
		    }

		    if (scaloc != 1.) {
			i__2 = *n;
			for (j = 1; j <= *n; ++j) {
			    dscal_(m, &scaloc, &C(1,j), &c__1);
/* L20: */
			}
			*scale *= scaloc;
		    }
		    C(k1,l1) = X(0);
		    C(k2,l1) = X(1);

		} else if (l1 != l2 && k1 == k2) {

		    i__2 = *m - k1;
/* Computing MIN */
		    i__3 = k1 + 1;
/* Computing MIN */
		    i__4 = k1 + 1;
		    suml = ddot_(&i__2, &A(k1,min(k1+1,*m)), lda, &
			    C(min(k1+1,*m),l1), &c__1);
		    i__2 = l1 - 1;
		    sumr = ddot_(&i__2, &C(k1,1), ldc, &B(1,l1), &c__1);
		    vec[0] = sgn * (C(k1,l1) - (suml + sgn * sumr))
			    ;

		    i__2 = *m - k1;
/* Computing MIN */
		    i__3 = k1 + 1;
/* Computing MIN */
		    i__4 = k1 + 1;
		    suml = ddot_(&i__2, &A(k1,min(k1+1,*m)), lda, &
			    C(min(k1+1,*m),l2), &c__1);
		    i__2 = l1 - 1;
		    sumr = ddot_(&i__2, &C(k1,1), ldc, &B(1,l2), &c__1);
		    vec[1] = sgn * (C(k1,l2) - (suml + sgn * sumr))
			    ;

		    d__1 = -sgn * A(k1,k1);
		    dlaln2_(&c_true, &c__2, &c__1, &smin, &c_b26, &B(l1,l1), ldb, &c_b26, &c_b26, vec, &c__2, &d__1, 
			    &c_b30, x, &c__2, &scaloc, &xnorm, &ierr);
		    if (ierr != 0) {
			*info = 1;
		    }

		    if (scaloc != 1.) {
			i__2 = *n;
			for (j = 1; j <= *n; ++j) {
			    dscal_(m, &scaloc, &C(1,j), &c__1);
/* L30: */
			}
			*scale *= scaloc;
		    }
		    C(k1,l1) = X(0);
		    C(k1,l2) = X(1);

		} else if (l1 != l2 && k1 != k2) {

		    i__2 = *m - k2;
/* Computing MIN */
		    i__3 = k2 + 1;
/* Computing MIN */
		    i__4 = k2 + 1;
		    suml = ddot_(&i__2, &A(k1,min(k2+1,*m)), lda, &
			    C(min(k2+1,*m),l1), &c__1);
		    i__2 = l1 - 1;
		    sumr = ddot_(&i__2, &C(k1,1), ldc, &B(1,l1), &c__1);
		    vec[0] = C(k1,l1) - (suml + sgn * sumr);

		    i__2 = *m - k2;
/* Computing MIN */
		    i__3 = k2 + 1;
/* Computing MIN */
		    i__4 = k2 + 1;
		    suml = ddot_(&i__2, &A(k1,min(k2+1,*m)), lda, &
			    C(min(k2+1,*m),l2), &c__1);
		    i__2 = l1 - 1;
		    sumr = ddot_(&i__2, &C(k1,1), ldc, &B(1,l2), &c__1);
		    vec[2] = C(k1,l2) - (suml + sgn * sumr);

		    i__2 = *m - k2;
/* Computing MIN */
		    i__3 = k2 + 1;
/* Computing MIN */
		    i__4 = k2 + 1;
		    suml = ddot_(&i__2, &A(k2,min(k2+1,*m)), lda, &
			    C(min(k2+1,*m),l1), &c__1);
		    i__2 = l1 - 1;
		    sumr = ddot_(&i__2, &C(k2,1), ldc, &B(1,l1), &c__1);
		    vec[1] = C(k2,l1) - (suml + sgn * sumr);

		    i__2 = *m - k2;
/* Computing MIN */
		    i__3 = k2 + 1;
/* Computing MIN */
		    i__4 = k2 + 1;
		    suml = ddot_(&i__2, &A(k2,min(k2+1,*m)), lda, &
			    C(min(k2+1,*m),l2), &c__1);
		    i__2 = l1 - 1;
		    sumr = ddot_(&i__2, &C(k2,1), ldc, &B(1,l2), &c__1);
		    vec[3] = C(k2,l2) - (suml + sgn * sumr);

		    dlasy2_(&c_false, &c_false, isgn, &c__2, &c__2, &A(k1,k1), lda, &B(l1,l1), ldb, vec,
			     &c__2, &scaloc, x, &c__2, &xnorm, &ierr);
		    if (ierr != 0) {
			*info = 1;
		    }

		    if (scaloc != 1.) {
			i__2 = *n;
			for (j = 1; j <= *n; ++j) {
			    dscal_(m, &scaloc, &C(1,j), &c__1);
/* L40: */
			}
			*scale *= scaloc;
		    }
		    C(k1,l1) = X(0);
		    C(k1,l2) = X(2);
		    C(k2,l1) = X(1);
		    C(k2,l2) = X(3);
		}

L50:
		;
	    }

L60:
	    ;
	}

    } else if (! notrna && notrnb) {

/*        Solve    A' *X + ISGN*X*B = scale*C.   

          The (K,L)th block of X is determined starting from   
          upper-left corner column by column by   

            A(K,K)'*X(K,L) + ISGN*X(K,L)*B(L,L) = C(K,L) - R(K,L)   

          Where   
                     K-1                        L-1   
            R(K,L) = SUM [A(I,K)'*X(I,L)] +ISGN*SUM [X(K,J)*B(J,L)]   
                     I=1                        J=1   

          Start column loop (index = L)   
          L1 (L2): column index of the first (last) row of X(K,L) */

	lnext = 1;
	i__1 = *n;
	for (l = 1; l <= *n; ++l) {
	    if (l < lnext) {
		goto L120;
	    }
	    if (l == *n) {
		l1 = l;
		l2 = l;
	    } else {
		if (B(l+1,l) != 0.) {
		    l1 = l;
		    l2 = l + 1;
		    lnext = l + 2;
		} else {
		    l1 = l;
		    l2 = l;
		    lnext = l + 1;
		}
	    }

/*           Start row loop (index = K)   
             K1 (K2): row index of the first (last) row of X(K,L) 
*/

	    knext = 1;
	    i__2 = *m;
	    for (k = 1; k <= *m; ++k) {
		if (k < knext) {
		    goto L110;
		}
		if (k == *m) {
		    k1 = k;
		    k2 = k;
		} else {
		    if (A(k+1,k) != 0.) {
			k1 = k;
			k2 = k + 1;
			knext = k + 2;
		    } else {
			k1 = k;
			k2 = k;
			knext = k + 1;
		    }
		}

		if (l1 == l2 && k1 == k2) {
		    i__3 = k1 - 1;
		    suml = ddot_(&i__3, &A(1,k1), &c__1, &C(1,l1), &c__1);
		    i__3 = l1 - 1;
		    sumr = ddot_(&i__3, &C(k1,1), ldc, &B(1,l1), &c__1);
		    vec[0] = C(k1,l1) - (suml + sgn * sumr);
		    scaloc = 1.;

		    a11 = A(k1,k1) + sgn * B(l1,l1);
		    da11 = abs(a11);
		    if (da11 <= smin) {
			a11 = smin;
			da11 = smin;
			*info = 1;
		    }
		    db = abs(vec[0]);
		    if (da11 < 1. && db > 1.) {
			if (db > bignum * da11) {
			    scaloc = 1. / db;
			}
		    }
		    X(0) = vec[0] * scaloc / a11;

		    if (scaloc != 1.) {
			i__3 = *n;
			for (j = 1; j <= *n; ++j) {
			    dscal_(m, &scaloc, &C(1,j), &c__1);
/* L70: */
			}
			*scale *= scaloc;
		    }
		    C(k1,l1) = X(0);

		} else if (l1 == l2 && k1 != k2) {

		    i__3 = k1 - 1;
		    suml = ddot_(&i__3, &A(1,k1), &c__1, &C(1,l1), &c__1);
		    i__3 = l1 - 1;
		    sumr = ddot_(&i__3, &C(k1,1), ldc, &B(1,l1), &c__1);
		    vec[0] = C(k1,l1) - (suml + sgn * sumr);

		    i__3 = k1 - 1;
		    suml = ddot_(&i__3, &A(1,k2), &c__1, &C(1,l1), &c__1);
		    i__3 = l1 - 1;
		    sumr = ddot_(&i__3, &C(k2,1), ldc, &B(1,l1), &c__1);
		    vec[1] = C(k2,l1) - (suml + sgn * sumr);

		    d__1 = -sgn * B(l1,l1);
		    dlaln2_(&c_true, &c__2, &c__1, &smin, &c_b26, &A(k1,k1), lda, &c_b26, &c_b26, vec, &c__2, &d__1, 
			    &c_b30, x, &c__2, &scaloc, &xnorm, &ierr);
		    if (ierr != 0) {
			*info = 1;
		    }

		    if (scaloc != 1.) {
			i__3 = *n;
			for (j = 1; j <= *n; ++j) {
			    dscal_(m, &scaloc, &C(1,j), &c__1);
/* L80: */
			}
			*scale *= scaloc;
		    }
		    C(k1,l1) = X(0);
		    C(k2,l1) = X(1);

		} else if (l1 != l2 && k1 == k2) {

		    i__3 = k1 - 1;
		    suml = ddot_(&i__3, &A(1,k1), &c__1, &C(1,l1), &c__1);
		    i__3 = l1 - 1;
		    sumr = ddot_(&i__3, &C(k1,1), ldc, &B(1,l1), &c__1);
		    vec[0] = sgn * (C(k1,l1) - (suml + sgn * sumr))
			    ;

		    i__3 = k1 - 1;
		    suml = ddot_(&i__3, &A(1,k1), &c__1, &C(1,l2), &c__1);
		    i__3 = l1 - 1;
		    sumr = ddot_(&i__3, &C(k1,1), ldc, &B(1,l2), &c__1);
		    vec[1] = sgn * (C(k1,l2) - (suml + sgn * sumr))
			    ;

		    d__1 = -sgn * A(k1,k1);
		    dlaln2_(&c_true, &c__2, &c__1, &smin, &c_b26, &B(l1,l1), ldb, &c_b26, &c_b26, vec, &c__2, &d__1, 
			    &c_b30, x, &c__2, &scaloc, &xnorm, &ierr);
		    if (ierr != 0) {
			*info = 1;
		    }

		    if (scaloc != 1.) {
			i__3 = *n;
			for (j = 1; j <= *n; ++j) {
			    dscal_(m, &scaloc, &C(1,j), &c__1);
/* L90: */
			}
			*scale *= scaloc;
		    }
		    C(k1,l1) = X(0);
		    C(k1,l2) = X(1);

		} else if (l1 != l2 && k1 != k2) {

		    i__3 = k1 - 1;
		    suml = ddot_(&i__3, &A(1,k1), &c__1, &C(1,l1), &c__1);
		    i__3 = l1 - 1;
		    sumr = ddot_(&i__3, &C(k1,1), ldc, &B(1,l1), &c__1);
		    vec[0] = C(k1,l1) - (suml + sgn * sumr);

		    i__3 = k1 - 1;
		    suml = ddot_(&i__3, &A(1,k1), &c__1, &C(1,l2), &c__1);
		    i__3 = l1 - 1;
		    sumr = ddot_(&i__3, &C(k1,1), ldc, &B(1,l2), &c__1);
		    vec[2] = C(k1,l2) - (suml + sgn * sumr);

		    i__3 = k1 - 1;
		    suml = ddot_(&i__3, &A(1,k2), &c__1, &C(1,l1), &c__1);
		    i__3 = l1 - 1;
		    sumr = ddot_(&i__3, &C(k2,1), ldc, &B(1,l1), &c__1);
		    vec[1] = C(k2,l1) - (suml + sgn * sumr);

		    i__3 = k1 - 1;
		    suml = ddot_(&i__3, &A(1,k2), &c__1, &C(1,l2), &c__1);
		    i__3 = l1 - 1;
		    sumr = ddot_(&i__3, &C(k2,1), ldc, &B(1,l2), &c__1);
		    vec[3] = C(k2,l2) - (suml + sgn * sumr);

		    dlasy2_(&c_true, &c_false, isgn, &c__2, &c__2, &A(k1,k1), lda, &B(l1,l1), ldb, vec, &
			    c__2, &scaloc, x, &c__2, &xnorm, &ierr);
		    if (ierr != 0) {
			*info = 1;
		    }

		    if (scaloc != 1.) {
			i__3 = *n;
			for (j = 1; j <= *n; ++j) {
			    dscal_(m, &scaloc, &C(1,j), &c__1);
/* L100: */
			}
			*scale *= scaloc;
		    }
		    C(k1,l1) = X(0);
		    C(k1,l2) = X(2);
		    C(k2,l1) = X(1);
		    C(k2,l2) = X(3);
		}

L110:
		;
	    }
L120:
	    ;
	}

    } else if (! notrna && ! notrnb) {

/*        Solve    A'*X + ISGN*X*B' = scale*C.   

          The (K,L)th block of X is determined starting from   
          top-right corner column by column by   

             A(K,K)'*X(K,L) + ISGN*X(K,L)*B(L,L)' = C(K,L) - R(K,L)   

          Where   
                       K-1                          N   
              R(K,L) = SUM [A(I,K)'*X(I,L)] + ISGN*SUM [X(K,J)*B(L,J)'
].   
                       I=1                        J=L+1   

          Start column loop (index = L)   
          L1 (L2): column index of the first (last) row of X(K,L) */

	lnext = *n;
	for (l = *n; l >= 1; --l) {
	    if (l > lnext) {
		goto L180;
	    }
	    if (l == 1) {
		l1 = l;
		l2 = l;
	    } else {
		if (B(l,l-1) != 0.) {
		    l1 = l - 1;
		    l2 = l;
		    lnext = l - 2;
		} else {
		    l1 = l;
		    l2 = l;
		    lnext = l - 1;
		}
	    }

/*           Start row loop (index = K)   
             K1 (K2): row index of the first (last) row of X(K,L) 
*/

	    knext = 1;
	    i__1 = *m;
	    for (k = 1; k <= *m; ++k) {
		if (k < knext) {
		    goto L170;
		}
		if (k == *m) {
		    k1 = k;
		    k2 = k;
		} else {
		    if (A(k+1,k) != 0.) {
			k1 = k;
			k2 = k + 1;
			knext = k + 2;
		    } else {
			k1 = k;
			k2 = k;
			knext = k + 1;
		    }
		}

		if (l1 == l2 && k1 == k2) {
		    i__2 = k1 - 1;
		    suml = ddot_(&i__2, &A(1,k1), &c__1, &C(1,l1), &c__1);
		    i__2 = *n - l1;
/* Computing MIN */
		    i__3 = l1 + 1;
/* Computing MIN */
		    i__4 = l1 + 1;
		    sumr = ddot_(&i__2, &C(k1,min(l1+1,*n)), ldc, &
			    B(l1,min(l1+1,*n)), ldb);
		    vec[0] = C(k1,l1) - (suml + sgn * sumr);
		    scaloc = 1.;

		    a11 = A(k1,k1) + sgn * B(l1,l1);
		    da11 = abs(a11);
		    if (da11 <= smin) {
			a11 = smin;
			da11 = smin;
			*info = 1;
		    }
		    db = abs(vec[0]);
		    if (da11 < 1. && db > 1.) {
			if (db > bignum * da11) {
			    scaloc = 1. / db;
			}
		    }
		    X(0) = vec[0] * scaloc / a11;

		    if (scaloc != 1.) {
			i__2 = *n;
			for (j = 1; j <= *n; ++j) {
			    dscal_(m, &scaloc, &C(1,j), &c__1);
/* L130: */
			}
			*scale *= scaloc;
		    }
		    C(k1,l1) = X(0);

		} else if (l1 == l2 && k1 != k2) {

		    i__2 = k1 - 1;
		    suml = ddot_(&i__2, &A(1,k1), &c__1, &C(1,l1), &c__1);
		    i__2 = *n - l2;
/* Computing MIN */
		    i__3 = l2 + 1;
/* Computing MIN */
		    i__4 = l2 + 1;
		    sumr = ddot_(&i__2, &C(k1,min(l2+1,*n)), ldc, &
			    B(l1,min(l2+1,*n)), ldb);
		    vec[0] = C(k1,l1) - (suml + sgn * sumr);

		    i__2 = k1 - 1;
		    suml = ddot_(&i__2, &A(1,k2), &c__1, &C(1,l1), &c__1);
		    i__2 = *n - l2;
/* Computing MIN */
		    i__3 = l2 + 1;
/* Computing MIN */
		    i__4 = l2 + 1;
		    sumr = ddot_(&i__2, &C(k2,min(l2+1,*n)), ldc, &
			    B(l1,min(l2+1,*n)), ldb);
		    vec[1] = C(k2,l1) - (suml + sgn * sumr);

		    d__1 = -sgn * B(l1,l1);
		    dlaln2_(&c_true, &c__2, &c__1, &smin, &c_b26, &A(k1,k1), lda, &c_b26, &c_b26, vec, &c__2, &d__1, 
			    &c_b30, x, &c__2, &scaloc, &xnorm, &ierr);
		    if (ierr != 0) {
			*info = 1;
		    }

		    if (scaloc != 1.) {
			i__2 = *n;
			for (j = 1; j <= *n; ++j) {
			    dscal_(m, &scaloc, &C(1,j), &c__1);
/* L140: */
			}
			*scale *= scaloc;
		    }
		    C(k1,l1) = X(0);
		    C(k2,l1) = X(1);

		} else if (l1 != l2 && k1 == k2) {

		    i__2 = k1 - 1;
		    suml = ddot_(&i__2, &A(1,k1), &c__1, &C(1,l1), &c__1);
		    i__2 = *n - l2;
/* Computing MIN */
		    i__3 = l2 + 1;
/* Computing MIN */
		    i__4 = l2 + 1;
		    sumr = ddot_(&i__2, &C(k1,min(l2+1,*n)), ldc, &
			    B(l1,min(l2+1,*n)), ldb);
		    vec[0] = sgn * (C(k1,l1) - (suml + sgn * sumr))
			    ;

		    i__2 = k1 - 1;
		    suml = ddot_(&i__2, &A(1,k1), &c__1, &C(1,l2), &c__1);
		    i__2 = *n - l2;
/* Computing MIN */
		    i__3 = l2 + 1;
/* Computing MIN */
		    i__4 = l2 + 1;
		    sumr = ddot_(&i__2, &C(k1,min(l2+1,*n)), ldc, &
			    B(l2,min(l2+1,*n)), ldb);
		    vec[1] = sgn * (C(k1,l2) - (suml + sgn * sumr))
			    ;

		    d__1 = -sgn * A(k1,k1);
		    dlaln2_(&c_false, &c__2, &c__1, &smin, &c_b26, &B(l1,l1), ldb, &c_b26, &c_b26, vec, &c__2, &d__1,
			     &c_b30, x, &c__2, &scaloc, &xnorm, &ierr);
		    if (ierr != 0) {
			*info = 1;
		    }

		    if (scaloc != 1.) {
			i__2 = *n;
			for (j = 1; j <= *n; ++j) {
			    dscal_(m, &scaloc, &C(1,j), &c__1);
/* L150: */
			}
			*scale *= scaloc;
		    }
		    C(k1,l1) = X(0);
		    C(k1,l2) = X(1);

		} else if (l1 != l2 && k1 != k2) {

		    i__2 = k1 - 1;
		    suml = ddot_(&i__2, &A(1,k1), &c__1, &C(1,l1), &c__1);
		    i__2 = *n - l2;
/* Computing MIN */
		    i__3 = l2 + 1;
/* Computing MIN */
		    i__4 = l2 + 1;
		    sumr = ddot_(&i__2, &C(k1,min(l2+1,*n)), ldc, &
			    B(l1,min(l2+1,*n)), ldb);
		    vec[0] = C(k1,l1) - (suml + sgn * sumr);

		    i__2 = k1 - 1;
		    suml = ddot_(&i__2, &A(1,k1), &c__1, &C(1,l2), &c__1);
		    i__2 = *n - l2;
/* Computing MIN */
		    i__3 = l2 + 1;
/* Computing MIN */
		    i__4 = l2 + 1;
		    sumr = ddot_(&i__2, &C(k1,min(l2+1,*n)), ldc, &
			    B(l2,min(l2+1,*n)), ldb);
		    vec[2] = C(k1,l2) - (suml + sgn * sumr);

		    i__2 = k1 - 1;
		    suml = ddot_(&i__2, &A(1,k2), &c__1, &C(1,l1), &c__1);
		    i__2 = *n - l2;
/* Computing MIN */
		    i__3 = l2 + 1;
/* Computing MIN */
		    i__4 = l2 + 1;
		    sumr = ddot_(&i__2, &C(k2,min(l2+1,*n)), ldc, &
			    B(l1,min(l2+1,*n)), ldb);
		    vec[1] = C(k2,l1) - (suml + sgn * sumr);

		    i__2 = k1 - 1;
		    suml = ddot_(&i__2, &A(1,k2), &c__1, &C(1,l2), &c__1);
		    i__2 = *n - l2;
/* Computing MIN */
		    i__3 = l2 + 1;
/* Computing MIN */
		    i__4 = l2 + 1;
		    sumr = ddot_(&i__2, &C(k2,min(l2+1,*n)), ldc, &
			    B(l2,min(l2+1,*n)), ldb);
		    vec[3] = C(k2,l2) - (suml + sgn * sumr);

		    dlasy2_(&c_true, &c_true, isgn, &c__2, &c__2, &A(k1,k1), lda, &B(l1,l1), ldb, vec, &
			    c__2, &scaloc, x, &c__2, &xnorm, &ierr);
		    if (ierr != 0) {
			*info = 1;
		    }

		    if (scaloc != 1.) {
			i__2 = *n;
			for (j = 1; j <= *n; ++j) {
			    dscal_(m, &scaloc, &C(1,j), &c__1);
/* L160: */
			}
			*scale *= scaloc;
		    }
		    C(k1,l1) = X(0);
		    C(k1,l2) = X(2);
		    C(k2,l1) = X(1);
		    C(k2,l2) = X(3);
		}

L170:
		;
	    }
L180:
	    ;
	}

    } else if (notrna && ! notrnb) {

/*        Solve    A*X + ISGN*X*B' = scale*C.   

          The (K,L)th block of X is determined starting from   
          bottom-right corner column by column by   

              A(K,K)*X(K,L) + ISGN*X(K,L)*B(L,L)' = C(K,L) - R(K,L)   

          Where   
                        M                          N   
              R(K,L) = SUM [A(K,I)*X(I,L)] + ISGN*SUM [X(K,J)*B(L,J)']
.   
                      I=K+1                      J=L+1   

          Start column loop (index = L)   
          L1 (L2): column index of the first (last) row of X(K,L) */

	lnext = *n;
	for (l = *n; l >= 1; --l) {
	    if (l > lnext) {
		goto L240;
	    }
	    if (l == 1) {
		l1 = l;
		l2 = l;
	    } else {
		if (B(l,l-1) != 0.) {
		    l1 = l - 1;
		    l2 = l;
		    lnext = l - 2;
		} else {
		    l1 = l;
		    l2 = l;
		    lnext = l - 1;
		}
	    }

/*           Start row loop (index = K)   
             K1 (K2): row index of the first (last) row of X(K,L) 
*/

	    knext = *m;
	    for (k = *m; k >= 1; --k) {
		if (k > knext) {
		    goto L230;
		}
		if (k == 1) {
		    k1 = k;
		    k2 = k;
		} else {
		    if (A(k,k-1) != 0.) {
			k1 = k - 1;
			k2 = k;
			knext = k - 2;
		    } else {
			k1 = k;
			k2 = k;
			knext = k - 1;
		    }
		}

		if (l1 == l2 && k1 == k2) {
		    i__1 = *m - k1;
/* Computing MIN */
		    i__2 = k1 + 1;
/* Computing MIN */
		    i__3 = k1 + 1;
		    suml = ddot_(&i__1, &A(k1,min(k1+1,*m)), lda, &
			    C(min(k1+1,*m),l1), &c__1);
		    i__1 = *n - l1;
/* Computing MIN */
		    i__2 = l1 + 1;
/* Computing MIN */
		    i__3 = l1 + 1;
		    sumr = ddot_(&i__1, &C(k1,min(l1+1,*n)), ldc, &
			    B(l1,min(l1+1,*n)), ldb);
		    vec[0] = C(k1,l1) - (suml + sgn * sumr);
		    scaloc = 1.;

		    a11 = A(k1,k1) + sgn * B(l1,l1);
		    da11 = abs(a11);
		    if (da11 <= smin) {
			a11 = smin;
			da11 = smin;
			*info = 1;
		    }
		    db = abs(vec[0]);
		    if (da11 < 1. && db > 1.) {
			if (db > bignum * da11) {
			    scaloc = 1. / db;
			}
		    }
		    X(0) = vec[0] * scaloc / a11;

		    if (scaloc != 1.) {
			i__1 = *n;
			for (j = 1; j <= *n; ++j) {
			    dscal_(m, &scaloc, &C(1,j), &c__1);
/* L190: */
			}
			*scale *= scaloc;
		    }
		    C(k1,l1) = X(0);

		} else if (l1 == l2 && k1 != k2) {

		    i__1 = *m - k2;
/* Computing MIN */
		    i__2 = k2 + 1;
/* Computing MIN */
		    i__3 = k2 + 1;
		    suml = ddot_(&i__1, &A(k1,min(k2+1,*m)), lda, &
			    C(min(k2+1,*m),l1), &c__1);
		    i__1 = *n - l2;
/* Computing MIN */
		    i__2 = l2 + 1;
/* Computing MIN */
		    i__3 = l2 + 1;
		    sumr = ddot_(&i__1, &C(k1,min(l2+1,*n)), ldc, &
			    B(l1,min(l2+1,*n)), ldb);
		    vec[0] = C(k1,l1) - (suml + sgn * sumr);

		    i__1 = *m - k2;
/* Computing MIN */
		    i__2 = k2 + 1;
/* Computing MIN */
		    i__3 = k2 + 1;
		    suml = ddot_(&i__1, &A(k2,min(k2+1,*m)), lda, &
			    C(min(k2+1,*m),l1), &c__1);
		    i__1 = *n - l2;
/* Computing MIN */
		    i__2 = l2 + 1;
/* Computing MIN */
		    i__3 = l2 + 1;
		    sumr = ddot_(&i__1, &C(k2,min(l2+1,*n)), ldc, &
			    B(l1,min(l2+1,*n)), ldb);
		    vec[1] = C(k2,l1) - (suml + sgn * sumr);

		    d__1 = -sgn * B(l1,l1);
		    dlaln2_(&c_false, &c__2, &c__1, &smin, &c_b26, &A(k1,k1), lda, &c_b26, &c_b26, vec, &c__2, &d__1,
			     &c_b30, x, &c__2, &scaloc, &xnorm, &ierr);
		    if (ierr != 0) {
			*info = 1;
		    }

		    if (scaloc != 1.) {
			i__1 = *n;
			for (j = 1; j <= *n; ++j) {
			    dscal_(m, &scaloc, &C(1,j), &c__1);
/* L200: */
			}
			*scale *= scaloc;
		    }
		    C(k1,l1) = X(0);
		    C(k2,l1) = X(1);

		} else if (l1 != l2 && k1 == k2) {

		    i__1 = *m - k1;
/* Computing MIN */
		    i__2 = k1 + 1;
/* Computing MIN */
		    i__3 = k1 + 1;
		    suml = ddot_(&i__1, &A(k1,min(k1+1,*m)), lda, &
			    C(min(k1+1,*m),l1), &c__1);
		    i__1 = *n - l2;
/* Computing MIN */
		    i__2 = l2 + 1;
/* Computing MIN */
		    i__3 = l2 + 1;
		    sumr = ddot_(&i__1, &C(k1,min(l2+1,*n)), ldc, &
			    B(l1,min(l2+1,*n)), ldb);
		    vec[0] = sgn * (C(k1,l1) - (suml + sgn * sumr))
			    ;

		    i__1 = *m - k1;
/* Computing MIN */
		    i__2 = k1 + 1;
/* Computing MIN */
		    i__3 = k1 + 1;
		    suml = ddot_(&i__1, &A(k1,min(k1+1,*m)), lda, &
			    C(min(k1+1,*m),l2), &c__1);
		    i__1 = *n - l2;
/* Computing MIN */
		    i__2 = l2 + 1;
/* Computing MIN */
		    i__3 = l2 + 1;
		    sumr = ddot_(&i__1, &C(k1,min(l2+1,*n)), ldc, &
			    B(l2,min(l2+1,*n)), ldb);
		    vec[1] = sgn * (C(k1,l2) - (suml + sgn * sumr))
			    ;

		    d__1 = -sgn * A(k1,k1);
		    dlaln2_(&c_false, &c__2, &c__1, &smin, &c_b26, &B(l1,l1), ldb, &c_b26, &c_b26, vec, &c__2, &d__1,
			     &c_b30, x, &c__2, &scaloc, &xnorm, &ierr);
		    if (ierr != 0) {
			*info = 1;
		    }

		    if (scaloc != 1.) {
			i__1 = *n;
			for (j = 1; j <= *n; ++j) {
			    dscal_(m, &scaloc, &C(1,j), &c__1);
/* L210: */
			}
			*scale *= scaloc;
		    }
		    C(k1,l1) = X(0);
		    C(k1,l2) = X(1);

		} else if (l1 != l2 && k1 != k2) {

		    i__1 = *m - k2;
/* Computing MIN */
		    i__2 = k2 + 1;
/* Computing MIN */
		    i__3 = k2 + 1;
		    suml = ddot_(&i__1, &A(k1,min(k2+1,*m)), lda, &
			    C(min(k2+1,*m),l1), &c__1);
		    i__1 = *n - l2;
/* Computing MIN */
		    i__2 = l2 + 1;
/* Computing MIN */
		    i__3 = l2 + 1;
		    sumr = ddot_(&i__1, &C(k1,min(l2+1,*n)), ldc, &
			    B(l1,min(l2+1,*n)), ldb);
		    vec[0] = C(k1,l1) - (suml + sgn * sumr);

		    i__1 = *m - k2;
/* Computing MIN */
		    i__2 = k2 + 1;
/* Computing MIN */
		    i__3 = k2 + 1;
		    suml = ddot_(&i__1, &A(k1,min(k2+1,*m)), lda, &
			    C(min(k2+1,*m),l2), &c__1);
		    i__1 = *n - l2;
/* Computing MIN */
		    i__2 = l2 + 1;
/* Computing MIN */
		    i__3 = l2 + 1;
		    sumr = ddot_(&i__1, &C(k1,min(l2+1,*n)), ldc, &
			    B(l2,min(l2+1,*n)), ldb);
		    vec[2] = C(k1,l2) - (suml + sgn * sumr);

		    i__1 = *m - k2;
/* Computing MIN */
		    i__2 = k2 + 1;
/* Computing MIN */
		    i__3 = k2 + 1;
		    suml = ddot_(&i__1, &A(k2,min(k2+1,*m)), lda, &
			    C(min(k2+1,*m),l1), &c__1);
		    i__1 = *n - l2;
/* Computing MIN */
		    i__2 = l2 + 1;
/* Computing MIN */
		    i__3 = l2 + 1;
		    sumr = ddot_(&i__1, &C(k2,min(l2+1,*n)), ldc, &
			    B(l1,min(l2+1,*n)), ldb);
		    vec[1] = C(k2,l1) - (suml + sgn * sumr);

		    i__1 = *m - k2;
/* Computing MIN */
		    i__2 = k2 + 1;
/* Computing MIN */
		    i__3 = k2 + 1;
		    suml = ddot_(&i__1, &A(k2,min(k2+1,*m)), lda, &
			    C(min(k2+1,*m),l2), &c__1);
		    i__1 = *n - l2;
/* Computing MIN */
		    i__2 = l2 + 1;
/* Computing MIN */
		    i__3 = l2 + 1;
		    sumr = ddot_(&i__1, &C(k2,min(l2+1,*n)), ldc, &
			    B(l2,min(l2+1,*n)), ldb);
		    vec[3] = C(k2,l2) - (suml + sgn * sumr);

		    dlasy2_(&c_false, &c_true, isgn, &c__2, &c__2, &A(k1,k1), lda, &B(l1,l1), ldb, vec, &
			    c__2, &scaloc, x, &c__2, &xnorm, &ierr);
		    if (ierr != 0) {
			*info = 1;
		    }

		    if (scaloc != 1.) {
			i__1 = *n;
			for (j = 1; j <= *n; ++j) {
			    dscal_(m, &scaloc, &C(1,j), &c__1);
/* L220: */
			}
			*scale *= scaloc;
		    }
		    C(k1,l1) = X(0);
		    C(k1,l2) = X(2);
		    C(k2,l1) = X(1);
		    C(k2,l2) = X(3);
		}

L230:
		;
	    }
L240:
	    ;
	}

    }

    return 0;

/*     End of DTRSYL */

} /* dtrsyl_ */

