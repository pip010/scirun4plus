#include "f2c.h"

/* Subroutine */ int ztprfs_(char *uplo, char *trans, char *diag, integer *n, 
	integer *nrhs, doublecomplex *ap, doublecomplex *b, integer *ldb, 
	doublecomplex *x, integer *ldx, doublereal *ferr, doublereal *berr, 
	doublecomplex *work, doublereal *rwork, integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    ZTPRFS provides error bounds and backward error estimates for the   
    solution to a system of linear equations with a triangular packed   
    coefficient matrix.   

    The solution matrix X must be computed by ZTPTRS or some other   
    means before entering this routine.  ZTPRFS does not do iterative   
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

    AP      (input) COMPLEX*16 array, dimension (N*(N+1)/2)   
            The upper or lower triangular matrix A, packed columnwise in 
  
            a linear array.  The j-th column of A is stored in the array 
  
            AP as follows:   
            if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;   
            if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.   
            If DIAG = 'U', the diagonal elements of A are not referenced 
  
            and are assumed to be 1.   

    B       (input) COMPLEX*16 array, dimension (LDB,NRHS)   
            The right hand side matrix B.   

    LDB     (input) INTEGER   
            The leading dimension of the array B.  LDB >= max(1,N).   

    X       (input) COMPLEX*16 array, dimension (LDX,NRHS)   
            The solution matrix X.   

    LDX     (input) INTEGER   
            The leading dimension of the array X.  LDX >= max(1,N).   

    FERR    (output) DOUBLE PRECISION array, dimension (NRHS)   
            The estimated forward error bound for each solution vector   
            X(j) (the j-th column of the solution matrix X).   
            If XTRUE is the true solution corresponding to X(j), FERR(j) 
  
            is an estimated upper bound for the magnitude of the largest 
  
            element in (X(j) - XTRUE) divided by the magnitude of the   
            largest element in X(j).  The estimate is as reliable as   
            the estimate for RCOND, and is almost always a slight   
            overestimate of the true error.   

    BERR    (output) DOUBLE PRECISION array, dimension (NRHS)   
            The componentwise relative backward error of each solution   
            vector X(j) (i.e., the smallest relative change in   
            any element of A or B that makes X(j) an exact solution).   

    WORK    (workspace) COMPLEX*16 array, dimension (2*N)   

    RWORK   (workspace) DOUBLE PRECISION array, dimension (N)   

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
    integer b_dim1, b_offset, x_dim1, x_offset, i__1, i__2, i__3, i__4, i__5;
    doublereal d__1, d__2, d__3, d__4;
    doublecomplex z__1;
    /* Builtin functions */
    double d_imag(doublecomplex *);
    /* Local variables */
    static integer kase;
    static doublereal safe1, safe2;
    static integer i, j, k;
    static doublereal s;
    extern logical lsame_(char *, char *);
    static logical upper;
    extern /* Subroutine */ int zcopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), zaxpy_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *), ztpmv_(
	    char *, char *, char *, integer *, doublecomplex *, doublecomplex 
	    *, integer *), ztpsv_(char *, char *, 
	    char *, integer *, doublecomplex *, doublecomplex *, integer *);
    static integer kc;
    extern doublereal dlamch_(char *);
    static doublereal xk;
    static integer nz;
    static doublereal safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *), zlacon_(
	    integer *, doublecomplex *, doublecomplex *, doublereal *, 
	    integer *);
    static logical notran;
    static char transn[1], transt[1];
    static logical nounit;
    static doublereal lstres, eps;



#define AP(I) ap[(I)-1]
#define FERR(I) ferr[(I)-1]
#define BERR(I) berr[(I)-1]
#define WORK(I) work[(I)-1]
#define RWORK(I) rwork[(I)-1]

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
    } else if (*ldb < max(1,*n)) {
	*info = -8;
    } else if (*ldx < max(1,*n)) {
	*info = -10;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("ZTPRFS", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0 || *nrhs == 0) {
	i__1 = *nrhs;
	for (j = 1; j <= *nrhs; ++j) {
	    FERR(j) = 0.;
	    BERR(j) = 0.;
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
    eps = dlamch_("Epsilon");
    safmin = dlamch_("Safe minimum");
    safe1 = nz * safmin;
    safe2 = safe1 / eps;

/*     Do for each right hand side */

    i__1 = *nrhs;
    for (j = 1; j <= *nrhs; ++j) {

/*        Compute residual R = B - op(A) * X,   
          where op(A) = A, A**T, or A**H, depending on TRANS. */

	zcopy_(n, &X(1,j), &c__1, &WORK(1), &c__1);
	ztpmv_(uplo, trans, diag, n, &AP(1), &WORK(1), &c__1);
	z__1.r = -1., z__1.i = 0.;
	zaxpy_(n, &z__1, &B(1,j), &c__1, &WORK(1), &c__1);

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
	    RWORK(i) = (d__1 = B(i,j).r, abs(d__1)) + (d__2 = d_imag(&B(i,j)), abs(d__2));
/* L20: */
	}

	if (notran) {

/*           Compute abs(A)*abs(X) + abs(B). */

	    if (upper) {
		kc = 1;
		if (nounit) {
		    i__2 = *n;
		    for (k = 1; k <= *n; ++k) {
			i__3 = k + j * x_dim1;
			xk = (d__1 = X(k,j).r, abs(d__1)) + (d__2 = d_imag(&
				X(k,j)), abs(d__2));
			i__3 = k;
			for (i = 1; i <= k; ++i) {
			    i__4 = kc + i - 1;
			    RWORK(i) += ((d__1 = AP(kc+i-1).r, abs(d__1)) + (
				    d__2 = d_imag(&AP(kc + i - 1)), abs(d__2))
				    ) * xk;
/* L30: */
			}
			kc += k;
/* L40: */
		    }
		} else {
		    i__2 = *n;
		    for (k = 1; k <= *n; ++k) {
			i__3 = k + j * x_dim1;
			xk = (d__1 = X(k,j).r, abs(d__1)) + (d__2 = d_imag(&
				X(k,j)), abs(d__2));
			i__3 = k - 1;
			for (i = 1; i <= k-1; ++i) {
			    i__4 = kc + i - 1;
			    RWORK(i) += ((d__1 = AP(kc+i-1).r, abs(d__1)) + (
				    d__2 = d_imag(&AP(kc + i - 1)), abs(d__2))
				    ) * xk;
/* L50: */
			}
			RWORK(k) += xk;
			kc += k;
/* L60: */
		    }
		}
	    } else {
		kc = 1;
		if (nounit) {
		    i__2 = *n;
		    for (k = 1; k <= *n; ++k) {
			i__3 = k + j * x_dim1;
			xk = (d__1 = X(k,j).r, abs(d__1)) + (d__2 = d_imag(&
				X(k,j)), abs(d__2));
			i__3 = *n;
			for (i = k; i <= *n; ++i) {
			    i__4 = kc + i - k;
			    RWORK(i) += ((d__1 = AP(kc+i-k).r, abs(d__1)) + (
				    d__2 = d_imag(&AP(kc + i - k)), abs(d__2))
				    ) * xk;
/* L70: */
			}
			kc = kc + *n - k + 1;
/* L80: */
		    }
		} else {
		    i__2 = *n;
		    for (k = 1; k <= *n; ++k) {
			i__3 = k + j * x_dim1;
			xk = (d__1 = X(k,j).r, abs(d__1)) + (d__2 = d_imag(&
				X(k,j)), abs(d__2));
			i__3 = *n;
			for (i = k + 1; i <= *n; ++i) {
			    i__4 = kc + i - k;
			    RWORK(i) += ((d__1 = AP(kc+i-k).r, abs(d__1)) + (
				    d__2 = d_imag(&AP(kc + i - k)), abs(d__2))
				    ) * xk;
/* L90: */
			}
			RWORK(k) += xk;
			kc = kc + *n - k + 1;
/* L100: */
		    }
		}
	    }
	} else {

/*           Compute abs(A**H)*abs(X) + abs(B). */

	    if (upper) {
		kc = 1;
		if (nounit) {
		    i__2 = *n;
		    for (k = 1; k <= *n; ++k) {
			s = 0.;
			i__3 = k;
			for (i = 1; i <= k; ++i) {
			    i__4 = kc + i - 1;
			    i__5 = i + j * x_dim1;
			    s += ((d__1 = AP(kc+i-1).r, abs(d__1)) + (d__2 = 
				    d_imag(&AP(kc + i - 1)), abs(d__2))) * ((
				    d__3 = X(i,j).r, abs(d__3)) + (d__4 = 
				    d_imag(&X(i,j)), abs(d__4)));
/* L110: */
			}
			RWORK(k) += s;
			kc += k;
/* L120: */
		    }
		} else {
		    i__2 = *n;
		    for (k = 1; k <= *n; ++k) {
			i__3 = k + j * x_dim1;
			s = (d__1 = X(k,j).r, abs(d__1)) + (d__2 = d_imag(&X(k,j)), abs(d__2));
			i__3 = k - 1;
			for (i = 1; i <= k-1; ++i) {
			    i__4 = kc + i - 1;
			    i__5 = i + j * x_dim1;
			    s += ((d__1 = AP(kc+i-1).r, abs(d__1)) + (d__2 = 
				    d_imag(&AP(kc + i - 1)), abs(d__2))) * ((
				    d__3 = X(i,j).r, abs(d__3)) + (d__4 = 
				    d_imag(&X(i,j)), abs(d__4)));
/* L130: */
			}
			RWORK(k) += s;
			kc += k;
/* L140: */
		    }
		}
	    } else {
		kc = 1;
		if (nounit) {
		    i__2 = *n;
		    for (k = 1; k <= *n; ++k) {
			s = 0.;
			i__3 = *n;
			for (i = k; i <= *n; ++i) {
			    i__4 = kc + i - k;
			    i__5 = i + j * x_dim1;
			    s += ((d__1 = AP(kc+i-k).r, abs(d__1)) + (d__2 = 
				    d_imag(&AP(kc + i - k)), abs(d__2))) * ((
				    d__3 = X(i,j).r, abs(d__3)) + (d__4 = 
				    d_imag(&X(i,j)), abs(d__4)));
/* L150: */
			}
			RWORK(k) += s;
			kc = kc + *n - k + 1;
/* L160: */
		    }
		} else {
		    i__2 = *n;
		    for (k = 1; k <= *n; ++k) {
			i__3 = k + j * x_dim1;
			s = (d__1 = X(k,j).r, abs(d__1)) + (d__2 = d_imag(&X(k,j)), abs(d__2));
			i__3 = *n;
			for (i = k + 1; i <= *n; ++i) {
			    i__4 = kc + i - k;
			    i__5 = i + j * x_dim1;
			    s += ((d__1 = AP(kc+i-k).r, abs(d__1)) + (d__2 = 
				    d_imag(&AP(kc + i - k)), abs(d__2))) * ((
				    d__3 = X(i,j).r, abs(d__3)) + (d__4 = 
				    d_imag(&X(i,j)), abs(d__4)));
/* L170: */
			}
			RWORK(k) += s;
			kc = kc + *n - k + 1;
/* L180: */
		    }
		}
	    }
	}
	s = 0.;
	i__2 = *n;
	for (i = 1; i <= *n; ++i) {
	    if (RWORK(i) > safe2) {
/* Computing MAX */
		i__3 = i;
		d__3 = s, d__4 = ((d__1 = WORK(i).r, abs(d__1)) + (d__2 = 
			d_imag(&WORK(i)), abs(d__2))) / RWORK(i);
		s = max(d__3,d__4);
	    } else {
/* Computing MAX */
		i__3 = i;
		d__3 = s, d__4 = ((d__1 = WORK(i).r, abs(d__1)) + (d__2 = 
			d_imag(&WORK(i)), abs(d__2)) + safe1) / (RWORK(i) + 
			safe1);
		s = max(d__3,d__4);
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

          Use ZLACON to estimate the infinity-norm of the matrix   
             inv(op(A)) * diag(W),   
          where W = abs(R) + NZ*EPS*( abs(op(A))*abs(X)+abs(B) ))) */

	i__2 = *n;
	for (i = 1; i <= *n; ++i) {
	    if (RWORK(i) > safe2) {
		i__3 = i;
		RWORK(i) = (d__1 = WORK(i).r, abs(d__1)) + (d__2 = d_imag(&
			WORK(i)), abs(d__2)) + nz * eps * RWORK(i);
	    } else {
		i__3 = i;
		RWORK(i) = (d__1 = WORK(i).r, abs(d__1)) + (d__2 = d_imag(&
			WORK(i)), abs(d__2)) + nz * eps * RWORK(i) + safe1;
	    }
/* L200: */
	}

	kase = 0;
L210:
	zlacon_(n, &WORK(*n + 1), &WORK(1), &FERR(j), &kase);
	if (kase != 0) {
	    if (kase == 1) {

/*              Multiply by diag(W)*inv(op(A)**H). */

		ztpsv_(uplo, transt, diag, n, &AP(1), &WORK(1), &c__1);
		i__2 = *n;
		for (i = 1; i <= *n; ++i) {
		    i__3 = i;
		    i__4 = i;
		    i__5 = i;
		    z__1.r = RWORK(i) * WORK(i).r, z__1.i = RWORK(i) 
			    * WORK(i).i;
		    WORK(i).r = z__1.r, WORK(i).i = z__1.i;
/* L220: */
		}
	    } else {

/*              Multiply by inv(op(A))*diag(W). */

		i__2 = *n;
		for (i = 1; i <= *n; ++i) {
		    i__3 = i;
		    i__4 = i;
		    i__5 = i;
		    z__1.r = RWORK(i) * WORK(i).r, z__1.i = RWORK(i) 
			    * WORK(i).i;
		    WORK(i).r = z__1.r, WORK(i).i = z__1.i;
/* L230: */
		}
		ztpsv_(uplo, transn, diag, n, &AP(1), &WORK(1), &c__1);
	    }
	    goto L210;
	}

/*        Normalize error. */

	lstres = 0.;
	i__2 = *n;
	for (i = 1; i <= *n; ++i) {
/* Computing MAX */
	    i__3 = i + j * x_dim1;
	    d__3 = lstres, d__4 = (d__1 = X(i,j).r, abs(d__1)) + (d__2 = 
		    d_imag(&X(i,j)), abs(d__2));
	    lstres = max(d__3,d__4);
/* L240: */
	}
	if (lstres != 0.) {
	    FERR(j) /= lstres;
	}

/* L250: */
    }

    return 0;

/*     End of ZTPRFS */

} /* ztprfs_ */

