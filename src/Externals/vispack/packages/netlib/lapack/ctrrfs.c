#include "f2c.h"

/* Subroutine */ int ctrrfs_(char *uplo, char *trans, char *diag, integer *n, 
	integer *nrhs, complex *a, integer *lda, complex *b, integer *ldb, 
	complex *x, integer *ldx, real *ferr, real *berr, complex *work, real 
	*rwork, integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    CTRRFS provides error bounds and backward error estimates for the   
    solution to a system of linear equations with a triangular   
    coefficient matrix.   

    The solution matrix X must be computed by CTRTRS or some other   
    means before entering this routine.  CTRRFS does not do iterative   
    refinement because doing so cannot improve the backward error.   

    Arguments   
    =========   

    UPLO    (input) CHARACTER*1   
            = 'U':  A is upper triangular;   
            = 'L':  A is lower triangular.   

    TRANS   (input) CHARACTER*1   
            Specifies the form of the system of equations:   
            = 'N':  A * X = B     (No transpose)   
            = 'T':  A**T * X = B  (Transpose)   
            = 'C':  A**H * X = B  (Conjugate transpose)   

    DIAG    (input) CHARACTER*1   
            = 'N':  A is non-unit triangular;   
            = 'U':  A is unit triangular.   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    NRHS    (input) INTEGER   
            The number of right hand sides, i.e., the number of columns   
            of the matrices B and X.  NRHS >= 0.   

    A       (input) COMPLEX array, dimension (LDA,N)   
            The triangular matrix A.  If UPLO = 'U', the leading N-by-N   
            upper triangular part of the array A contains the upper   
            triangular matrix, and the strictly lower triangular part of 
  
            A is not referenced.  If UPLO = 'L', the leading N-by-N lower 
  
            triangular part of the array A contains the lower triangular 
  
            matrix, and the strictly upper triangular part of A is not   
            referenced.  If DIAG = 'U', the diagonal elements of A are   
            also not referenced and are assumed to be 1.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(1,N).   

    B       (input) COMPLEX array, dimension (LDB,NRHS)   
            The right hand side matrix B.   

    LDB     (input) INTEGER   
            The leading dimension of the array B.  LDB >= max(1,N).   

    X       (input) COMPLEX array, dimension (LDX,NRHS)   
            The solution matrix X.   

    LDX     (input) INTEGER   
            The leading dimension of the array X.  LDX >= max(1,N).   

    FERR    (output) REAL array, dimension (NRHS)   
            The estimated forward error bound for each solution vector   
            X(j) (the j-th column of the solution matrix X).   
            If XTRUE is the true solution corresponding to X(j), FERR(j) 
  
            is an estimated upper bound for the magnitude of the largest 
  
            element in (X(j) - XTRUE) divided by the magnitude of the   
            largest element in X(j).  The estimate is as reliable as   
            the estimate for RCOND, and is almost always a slight   
            overestimate of the true error.   

    BERR    (output) REAL array, dimension (NRHS)   
            The componentwise relative backward error of each solution   
            vector X(j) (i.e., the smallest relative change in   
            any element of A or B that makes X(j) an exact solution).   

    WORK    (workspace) COMPLEX array, dimension (2*N)   

    RWORK   (workspace) REAL array, dimension (N)   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   

    ===================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static integer c__1 = 1;
    
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, x_dim1, x_offset, i__1, i__2, 
	    i__3, i__4, i__5;
    real r__1, r__2, r__3, r__4;
    complex q__1;
    /* Builtin functions */
    double r_imag(complex *);
    /* Local variables */
    static integer kase;
    static real safe1, safe2;
    static integer i, j, k;
    static real s;
    extern logical lsame_(char *, char *);
    extern /* Subroutine */ int ccopy_(integer *, complex *, integer *, 
	    complex *, integer *), caxpy_(integer *, complex *, complex *, 
	    integer *, complex *, integer *);
    static logical upper;
    extern /* Subroutine */ int ctrmv_(char *, char *, char *, integer *, 
	    complex *, integer *, complex *, integer *), ctrsv_(char *, char *, char *, integer *, complex *, 
	    integer *, complex *, integer *), clacon_(
	    integer *, complex *, complex *, real *, integer *);
    static real xk;
    extern doublereal slamch_(char *);
    static integer nz;
    static real safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *);
    static logical notran;
    static char transn[1], transt[1];
    static logical nounit;
    static real lstres, eps;



#define FERR(I) ferr[(I)-1]
#define BERR(I) berr[(I)-1]
#define WORK(I) work[(I)-1]
#define RWORK(I) rwork[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]
#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]
#define X(I,J) x[(I)-1 + ((J)-1)* ( *ldx)]

    *info = 0;
    upper = lsame_(uplo, "U");
    notran = lsame_(trans, "N");
    nounit = lsame_(diag, "N");

    if (! upper && ! lsame_(uplo, "L")) {
	*info = -1;
    } else if (! notran && ! lsame_(trans, "T") && ! lsame_(trans, 
	    "C")) {
	*info = -2;
    } else if (! nounit && ! lsame_(diag, "U")) {
	*info = -3;
    } else if (*n < 0) {
	*info = -4;
    } else if (*nrhs < 0) {
	*info = -5;
    } else if (*lda < max(1,*n)) {
	*info = -7;
    } else if (*ldb < max(1,*n)) {
	*info = -9;
    } else if (*ldx < max(1,*n)) {
	*info = -11;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("CTRRFS", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0 || *nrhs == 0) {
	i__1 = *nrhs;
	for (j = 1; j <= *nrhs; ++j) {
	    FERR(j) = 0.f;
	    BERR(j) = 0.f;
/* L10: */
	}
	return 0;
    }

    if (notran) {
	*(unsigned char *)transn = 'N';
	*(unsigned char *)transt = 'C';
    } else {
	*(unsigned char *)transn = 'C';
	*(unsigned char *)transt = 'N';
    }

/*     NZ = maximum number of nonzero elements in each row of A, plus 1 */

    nz = *n + 1;
    eps = slamch_("Epsilon");
    safmin = slamch_("Safe minimum");
    safe1 = nz * safmin;
    safe2 = safe1 / eps;

/*     Do for each right hand side */

    i__1 = *nrhs;
    for (j = 1; j <= *nrhs; ++j) {

/*        Compute residual R = B - op(A) * X,   
          where op(A) = A, A**T, or A**H, depending on TRANS. */

	ccopy_(n, &X(1,j), &c__1, &WORK(1), &c__1);
	ctrmv_(uplo, trans, diag, n, &A(1,1), lda, &WORK(1), &c__1);
	q__1.r = -1.f, q__1.i = 0.f;
	caxpy_(n, &q__1, &B(1,j), &c__1, &WORK(1), &c__1);

/*        Compute componentwise relative backward error from formula 
  

          max(i) ( abs(R(i)) / ( abs(op(A))*abs(X) + abs(B) )(i) )   

          where abs(Z) is the componentwise absolute value of the matr
ix   
          or vector Z.  If the i-th component of the denominator is le
ss   
          than SAFE2, then SAFE1 is added to the i-th components of th
e   
          numerator and denominator before dividing. */

	i__2 = *n;
	for (i = 1; i <= *n; ++i) {
	    i__3 = i + j * b_dim1;
	    RWORK(i) = (r__1 = B(i,j).r, dabs(r__1)) + (r__2 = r_imag(&B(i,j)), dabs(r__2));
/* L20: */
	}

	if (notran) {

/*           Compute abs(A)*abs(X) + abs(B). */

	    if (upper) {
		if (nounit) {
		    i__2 = *n;
		    for (k = 1; k <= *n; ++k) {
			i__3 = k + j * x_dim1;
			xk = (r__1 = X(k,j).r, dabs(r__1)) + (r__2 = r_imag(&
				X(k,j)), dabs(r__2));
			i__3 = k;
			for (i = 1; i <= k; ++i) {
			    i__4 = i + k * a_dim1;
			    RWORK(i) += ((r__1 = A(i,k).r, dabs(r__1)) + (
				    r__2 = r_imag(&A(i,k)), dabs(
				    r__2))) * xk;
/* L30: */
			}
/* L40: */
		    }
		} else {
		    i__2 = *n;
		    for (k = 1; k <= *n; ++k) {
			i__3 = k + j * x_dim1;
			xk = (r__1 = X(k,j).r, dabs(r__1)) + (r__2 = r_imag(&
				X(k,j)), dabs(r__2));
			i__3 = k - 1;
			for (i = 1; i <= k-1; ++i) {
			    i__4 = i + k * a_dim1;
			    RWORK(i) += ((r__1 = A(i,k).r, dabs(r__1)) + (
				    r__2 = r_imag(&A(i,k)), dabs(
				    r__2))) * xk;
/* L50: */
			}
			RWORK(k) += xk;
/* L60: */
		    }
		}
	    } else {
		if (nounit) {
		    i__2 = *n;
		    for (k = 1; k <= *n; ++k) {
			i__3 = k + j * x_dim1;
			xk = (r__1 = X(k,j).r, dabs(r__1)) + (r__2 = r_imag(&
				X(k,j)), dabs(r__2));
			i__3 = *n;
			for (i = k; i <= *n; ++i) {
			    i__4 = i + k * a_dim1;
			    RWORK(i) += ((r__1 = A(i,k).r, dabs(r__1)) + (
				    r__2 = r_imag(&A(i,k)), dabs(
				    r__2))) * xk;
/* L70: */
			}
/* L80: */
		    }
		} else {
		    i__2 = *n;
		    for (k = 1; k <= *n; ++k) {
			i__3 = k + j * x_dim1;
			xk = (r__1 = X(k,j).r, dabs(r__1)) + (r__2 = r_imag(&
				X(k,j)), dabs(r__2));
			i__3 = *n;
			for (i = k + 1; i <= *n; ++i) {
			    i__4 = i + k * a_dim1;
			    RWORK(i) += ((r__1 = A(i,k).r, dabs(r__1)) + (
				    r__2 = r_imag(&A(i,k)), dabs(
				    r__2))) * xk;
/* L90: */
			}
			RWORK(k) += xk;
/* L100: */
		    }
		}
	    }
	} else {

/*           Compute abs(A**H)*abs(X) + abs(B). */

	    if (upper) {
		if (nounit) {
		    i__2 = *n;
		    for (k = 1; k <= *n; ++k) {
			s = 0.f;
			i__3 = k;
			for (i = 1; i <= k; ++i) {
			    i__4 = i + k * a_dim1;
			    i__5 = i + j * x_dim1;
			    s += ((r__1 = A(i,k).r, dabs(r__1)) + (r__2 = 
				    r_imag(&A(i,k)), dabs(r__2))) *
				     ((r__3 = X(i,j).r, dabs(r__3)) + (r__4 =
				     r_imag(&X(i,j)), dabs(r__4)));
/* L110: */
			}
			RWORK(k) += s;
/* L120: */
		    }
		} else {
		    i__2 = *n;
		    for (k = 1; k <= *n; ++k) {
			i__3 = k + j * x_dim1;
			s = (r__1 = X(k,j).r, dabs(r__1)) + (r__2 = r_imag(&
				X(k,j)), dabs(r__2));
			i__3 = k - 1;
			for (i = 1; i <= k-1; ++i) {
			    i__4 = i + k * a_dim1;
			    i__5 = i + j * x_dim1;
			    s += ((r__1 = A(i,k).r, dabs(r__1)) + (r__2 = 
				    r_imag(&A(i,k)), dabs(r__2))) *
				     ((r__3 = X(i,j).r, dabs(r__3)) + (r__4 =
				     r_imag(&X(i,j)), dabs(r__4)));
/* L130: */
			}
			RWORK(k) += s;
/* L140: */
		    }
		}
	    } else {
		if (nounit) {
		    i__2 = *n;
		    for (k = 1; k <= *n; ++k) {
			s = 0.f;
			i__3 = *n;
			for (i = k; i <= *n; ++i) {
			    i__4 = i + k * a_dim1;
			    i__5 = i + j * x_dim1;
			    s += ((r__1 = A(i,k).r, dabs(r__1)) + (r__2 = 
				    r_imag(&A(i,k)), dabs(r__2))) *
				     ((r__3 = X(i,j).r, dabs(r__3)) + (r__4 =
				     r_imag(&X(i,j)), dabs(r__4)));
/* L150: */
			}
			RWORK(k) += s;
/* L160: */
		    }
		} else {
		    i__2 = *n;
		    for (k = 1; k <= *n; ++k) {
			i__3 = k + j * x_dim1;
			s = (r__1 = X(k,j).r, dabs(r__1)) + (r__2 = r_imag(&
				X(k,j)), dabs(r__2));
			i__3 = *n;
			for (i = k + 1; i <= *n; ++i) {
			    i__4 = i + k * a_dim1;
			    i__5 = i + j * x_dim1;
			    s += ((r__1 = A(i,k).r, dabs(r__1)) + (r__2 = 
				    r_imag(&A(i,k)), dabs(r__2))) *
				     ((r__3 = X(i,j).r, dabs(r__3)) + (r__4 =
				     r_imag(&X(i,j)), dabs(r__4)));
/* L170: */
			}
			RWORK(k) += s;
/* L180: */
		    }
		}
	    }
	}
	s = 0.f;
	i__2 = *n;
	for (i = 1; i <= *n; ++i) {
	    if (RWORK(i) > safe2) {
/* Computing MAX */
		i__3 = i;
		r__3 = s, r__4 = ((r__1 = WORK(i).r, dabs(r__1)) + (r__2 = 
			r_imag(&WORK(i)), dabs(r__2))) / RWORK(i);
		s = dmax(r__3,r__4);
	    } else {
/* Computing MAX */
		i__3 = i;
		r__3 = s, r__4 = ((r__1 = WORK(i).r, dabs(r__1)) + (r__2 = 
			r_imag(&WORK(i)), dabs(r__2)) + safe1) / (RWORK(i) + 
			safe1);
		s = dmax(r__3,r__4);
	    }
/* L190: */
	}
	BERR(j) = s;

/*        Bound error from formula   

          norm(X - XTRUE) / norm(X) .le. FERR =   
          norm( abs(inv(op(A)))*   
             ( abs(R) + NZ*EPS*( abs(op(A))*abs(X)+abs(B) ))) / norm(X
)   

          where   
            norm(Z) is the magnitude of the largest component of Z   
            inv(op(A)) is the inverse of op(A)   
            abs(Z) is the componentwise absolute value of the matrix o
r   
               vector Z   
            NZ is the maximum number of nonzeros in any row of A, plus
 1   
            EPS is machine epsilon   

          The i-th component of abs(R)+NZ*EPS*(abs(op(A))*abs(X)+abs(B
))   
          is incremented by SAFE1 if the i-th component of   
          abs(op(A))*abs(X) + abs(B) is less than SAFE2.   

          Use CLACON to estimate the infinity-norm of the matrix   
             inv(op(A)) * diag(W),   
          where W = abs(R) + NZ*EPS*( abs(op(A))*abs(X)+abs(B) ))) */

	i__2 = *n;
	for (i = 1; i <= *n; ++i) {
	    if (RWORK(i) > safe2) {
		i__3 = i;
		RWORK(i) = (r__1 = WORK(i).r, dabs(r__1)) + (r__2 = r_imag(
			&WORK(i)), dabs(r__2)) + nz * eps * RWORK(i);
	    } else {
		i__3 = i;
		RWORK(i) = (r__1 = WORK(i).r, dabs(r__1)) + (r__2 = r_imag(
			&WORK(i)), dabs(r__2)) + nz * eps * RWORK(i) + safe1;
	    }
/* L200: */
	}

	kase = 0;
L210:
	clacon_(n, &WORK(*n + 1), &WORK(1), &FERR(j), &kase);
	if (kase != 0) {
	    if (kase == 1) {

/*              Multiply by diag(W)*inv(op(A)**H). */

		ctrsv_(uplo, transt, diag, n, &A(1,1), lda, &WORK(1), &
			c__1);
		i__2 = *n;
		for (i = 1; i <= *n; ++i) {
		    i__3 = i;
		    i__4 = i;
		    i__5 = i;
		    q__1.r = RWORK(i) * WORK(i).r, q__1.i = RWORK(i) 
			    * WORK(i).i;
		    WORK(i).r = q__1.r, WORK(i).i = q__1.i;
/* L220: */
		}
	    } else {

/*              Multiply by inv(op(A))*diag(W). */

		i__2 = *n;
		for (i = 1; i <= *n; ++i) {
		    i__3 = i;
		    i__4 = i;
		    i__5 = i;
		    q__1.r = RWORK(i) * WORK(i).r, q__1.i = RWORK(i) 
			    * WORK(i).i;
		    WORK(i).r = q__1.r, WORK(i).i = q__1.i;
/* L230: */
		}
		ctrsv_(uplo, transn, diag, n, &A(1,1), lda, &WORK(1), &
			c__1);
	    }
	    goto L210;
	}

/*        Normalize error. */

	lstres = 0.f;
	i__2 = *n;
	for (i = 1; i <= *n; ++i) {
/* Computing MAX */
	    i__3 = i + j * x_dim1;
	    r__3 = lstres, r__4 = (r__1 = X(i,j).r, dabs(r__1)) + (r__2 = 
		    r_imag(&X(i,j)), dabs(r__2));
	    lstres = dmax(r__3,r__4);
/* L240: */
	}
	if (lstres != 0.f) {
	    FERR(j) /= lstres;
	}

/* L250: */
    }

    return 0;

/*     End of CTRRFS */

} /* ctrrfs_ */

