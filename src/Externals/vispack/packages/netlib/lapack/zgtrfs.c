#include "f2c.h"

/* Subroutine */ int zgtrfs_(char *trans, integer *n, integer *nrhs, 
	doublecomplex *dl, doublecomplex *d, doublecomplex *du, doublecomplex 
	*dlf, doublecomplex *df, doublecomplex *duf, doublecomplex *du2, 
	integer *ipiv, doublecomplex *b, integer *ldb, doublecomplex *x, 
	integer *ldx, doublereal *ferr, doublereal *berr, doublecomplex *work,
	 doublereal *rwork, integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    ZGTRFS improves the computed solution to a system of linear   
    equations when the coefficient matrix is tridiagonal, and provides   
    error bounds and backward error estimates for the solution.   

    Arguments   
    =========   

    TRANS   (input) CHARACTER*1   
            Specifies the form of the system of equations:   
            = 'N':  A * X = B     (No transpose)   
            = 'T':  A**T * X = B  (Transpose)   
            = 'C':  A**H * X = B  (Conjugate transpose)   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    NRHS    (input) INTEGER   
            The number of right hand sides, i.e., the number of columns   
            of the matrix B.  NRHS >= 0.   

    DL      (input) COMPLEX*16 array, dimension (N-1)   
            The (n-1) subdiagonal elements of A.   

    D       (input) COMPLEX*16 array, dimension (N)   
            The diagonal elements of A.   

    DU      (input) COMPLEX*16 array, dimension (N-1)   
            The (n-1) superdiagonal elements of A.   

    DLF     (input) COMPLEX*16 array, dimension (N-1)   
            The (n-1) multipliers that define the matrix L from the   
            LU factorization of A as computed by ZGTTRF.   

    DF      (input) COMPLEX*16 array, dimension (N)   
            The n diagonal elements of the upper triangular matrix U from 
  
            the LU factorization of A.   

    DUF     (input) COMPLEX*16 array, dimension (N-1)   
            The (n-1) elements of the first superdiagonal of U.   

    DU2     (input) COMPLEX*16 array, dimension (N-2)   
            The (n-2) elements of the second superdiagonal of U.   

    IPIV    (input) INTEGER array, dimension (N)   
            The pivot indices; for 1 <= i <= n, row i of the matrix was   
            interchanged with row IPIV(i).  IPIV(i) will always be either 
  
            i or i+1; IPIV(i) = i indicates a row interchange was not   
            required.   

    B       (input) COMPLEX*16 array, dimension (LDB,NRHS)   
            The right hand side matrix B.   

    LDB     (input) INTEGER   
            The leading dimension of the array B.  LDB >= max(1,N).   

    X       (input/output) COMPLEX*16 array, dimension (LDX,NRHS)   
            On entry, the solution matrix X, as computed by ZGTTRS.   
            On exit, the improved solution matrix X.   

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

    Internal Parameters   
    ===================   

    ITMAX is the maximum number of steps of iterative refinement.   

    ===================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static integer c__1 = 1;
    static doublereal c_b18 = -1.;
    static doublereal c_b19 = 1.;
    static doublecomplex c_b26 = {1.,0.};
    
    /* System generated locals */
    integer b_dim1, b_offset, x_dim1, x_offset, i__1, i__2, i__3, i__4, i__5, 
	    i__6, i__7, i__8, i__9;
    doublereal d__1, d__2, d__3, d__4, d__5, d__6, d__7, d__8, d__9, d__10, 
	    d__11, d__12, d__13, d__14;
    doublecomplex z__1;
    /* Builtin functions */
    double d_imag(doublecomplex *);
    /* Local variables */
    static integer kase;
    static doublereal safe1, safe2;
    static integer i, j;
    static doublereal s;
    extern logical lsame_(char *, char *);
    static integer count;
    extern /* Subroutine */ int zcopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), zaxpy_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    extern doublereal dlamch_(char *);
    static integer nz;
    static doublereal safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *), zlacon_(
	    integer *, doublecomplex *, doublecomplex *, doublereal *, 
	    integer *), zlagtm_(char *, integer *, integer *, doublereal *, 
	    doublecomplex *, doublecomplex *, doublecomplex *, doublecomplex *
	    , integer *, doublereal *, doublecomplex *, integer *);
    static logical notran;
    static char transn[1], transt[1];
    static doublereal lstres;
    extern /* Subroutine */ int zgttrs_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, doublecomplex *, doublecomplex *
	    , integer *, doublecomplex *, integer *, integer *);
    static doublereal eps;



#define DL(I) dl[(I)-1]
#define D(I) d[(I)-1]
#define DU(I) du[(I)-1]
#define DLF(I) dlf[(I)-1]
#define DF(I) df[(I)-1]
#define DUF(I) duf[(I)-1]
#define DU2(I) du2[(I)-1]
#define IPIV(I) ipiv[(I)-1]
#define FERR(I) ferr[(I)-1]
#define BERR(I) berr[(I)-1]
#define WORK(I) work[(I)-1]
#define RWORK(I) rwork[(I)-1]

#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]
#define X(I,J) x[(I)-1 + ((J)-1)* ( *ldx)]

    *info = 0;
    notran = lsame_(trans, "N");
    if (! notran && ! lsame_(trans, "T") && ! lsame_(trans, "C")) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*nrhs < 0) {
	*info = -3;
    } else if (*ldb < max(1,*n)) {
	*info = -13;
    } else if (*ldx < max(1,*n)) {
	*info = -15;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("ZGTRFS", &i__1);
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

    nz = 4;
    eps = dlamch_("Epsilon");
    safmin = dlamch_("Safe minimum");
    safe1 = nz * safmin;
    safe2 = safe1 / eps;

/*     Do for each right hand side */

    i__1 = *nrhs;
    for (j = 1; j <= *nrhs; ++j) {

	count = 1;
	lstres = 3.;
L20:

/*        Loop until stopping criterion is satisfied.   

          Compute residual R = B - op(A) * X,   
          where op(A) = A, A**T, or A**H, depending on TRANS. */

	zcopy_(n, &B(1,j), &c__1, &WORK(1), &c__1);
	zlagtm_(trans, n, &c__1, &c_b18, &DL(1), &D(1), &DU(1), &X(1,j), ldx, &c_b19, &WORK(1), n);

/*        Compute abs(op(A))*abs(x) + abs(b) for use in the backward 
  
          error bound. */

	if (notran) {
	    if (*n == 1) {
		i__2 = j * b_dim1 + 1;
		i__3 = j * x_dim1 + 1;
		RWORK(1) = (d__1 = B(1,j).r, abs(d__1)) + (d__2 = d_imag(&B(1,j)), abs(d__2)) + ((d__3 = D(1).r, abs(
			d__3)) + (d__4 = d_imag(&D(1)), abs(d__4))) * ((d__5 =
			 X(1,j).r, abs(d__5)) + (d__6 = d_imag(&X(1,j)), abs(d__6)));
	    } else {
		i__2 = j * b_dim1 + 1;
		i__3 = j * x_dim1 + 1;
		i__4 = j * x_dim1 + 2;
		RWORK(1) = (d__1 = B(1,j).r, abs(d__1)) + (d__2 = d_imag(&B(1,j)), abs(d__2)) + ((d__3 = D(1).r, abs(
			d__3)) + (d__4 = d_imag(&D(1)), abs(d__4))) * ((d__5 =
			 X(1,j).r, abs(d__5)) + (d__6 = d_imag(&X(1,j)), abs(d__6))) + ((d__7 = DU(1).r, abs(d__7)) + (
			d__8 = d_imag(&DU(1)), abs(d__8))) * ((d__9 = X(2,j)
			.r, abs(d__9)) + (d__10 = d_imag(&X(2,j)), 
			abs(d__10)));
		i__2 = *n - 1;
		for (i = 2; i <= *n-1; ++i) {
		    i__3 = i + j * b_dim1;
		    i__4 = i - 1;
		    i__5 = i - 1 + j * x_dim1;
		    i__6 = i;
		    i__7 = i + j * x_dim1;
		    i__8 = i;
		    i__9 = i + 1 + j * x_dim1;
		    RWORK(i) = (d__1 = B(i,j).r, abs(d__1)) + (d__2 = d_imag(
			    &B(i,j)), abs(d__2)) + ((d__3 = DL(
			    i-1).r, abs(d__3)) + (d__4 = d_imag(&DL(i - 1)), 
			    abs(d__4))) * ((d__5 = X(i-1,j).r, abs(d__5)) + (
			    d__6 = d_imag(&X(i-1,j)), abs(d__6)))
			     + ((d__7 = D(i).r, abs(d__7)) + (d__8 = 
			    d_imag(&D(i)), abs(d__8))) * ((d__9 = X(i,j).r, 
			    abs(d__9)) + (d__10 = d_imag(&X(i,j)), 
			    abs(d__10))) + ((d__11 = DU(i).r, abs(d__11)) 
			    + (d__12 = d_imag(&DU(i)), abs(d__12))) * ((d__13 
			    = X(i+1,j).r, abs(d__13)) + (d__14 = d_imag(&X(i+1,j)), abs(d__14)));
/* L30: */
		}
		i__2 = *n + j * b_dim1;
		i__3 = *n - 1;
		i__4 = *n - 1 + j * x_dim1;
		i__5 = *n;
		i__6 = *n + j * x_dim1;
		RWORK(*n) = (d__1 = B(*n,j).r, abs(d__1)) + (d__2 = d_imag(&B(*n,j)), abs(d__2)) + ((d__3 = DL(*n-1).r, 
			abs(d__3)) + (d__4 = d_imag(&DL(*n - 1)), abs(d__4))) 
			* ((d__5 = X(*n-1,j).r, abs(d__5)) + (d__6 = d_imag(&X(*n-1,j)), abs(d__6))) + ((d__7 = D(*n)
			.r, abs(d__7)) + (d__8 = d_imag(&D(*n)), abs(d__8))) *
			 ((d__9 = X(*n,j).r, abs(d__9)) + (d__10 = d_imag(&X(*n,j)), abs(d__10)));
	    }
	} else {
	    if (*n == 1) {
		i__2 = j * b_dim1 + 1;
		i__3 = j * x_dim1 + 1;
		RWORK(1) = (d__1 = B(1,j).r, abs(d__1)) + (d__2 = d_imag(&B(1,j)), abs(d__2)) + ((d__3 = D(1).r, abs(
			d__3)) + (d__4 = d_imag(&D(1)), abs(d__4))) * ((d__5 =
			 X(1,j).r, abs(d__5)) + (d__6 = d_imag(&X(1,j)), abs(d__6)));
	    } else {
		i__2 = j * b_dim1 + 1;
		i__3 = j * x_dim1 + 1;
		i__4 = j * x_dim1 + 2;
		RWORK(1) = (d__1 = B(1,j).r, abs(d__1)) + (d__2 = d_imag(&B(1,j)), abs(d__2)) + ((d__3 = D(1).r, abs(
			d__3)) + (d__4 = d_imag(&D(1)), abs(d__4))) * ((d__5 =
			 X(1,j).r, abs(d__5)) + (d__6 = d_imag(&X(1,j)), abs(d__6))) + ((d__7 = DL(1).r, abs(d__7)) + (
			d__8 = d_imag(&DL(1)), abs(d__8))) * ((d__9 = X(2,j)
			.r, abs(d__9)) + (d__10 = d_imag(&X(2,j)), 
			abs(d__10)));
		i__2 = *n - 1;
		for (i = 2; i <= *n-1; ++i) {
		    i__3 = i + j * b_dim1;
		    i__4 = i - 1;
		    i__5 = i - 1 + j * x_dim1;
		    i__6 = i;
		    i__7 = i + j * x_dim1;
		    i__8 = i;
		    i__9 = i + 1 + j * x_dim1;
		    RWORK(i) = (d__1 = B(i,j).r, abs(d__1)) + (d__2 = d_imag(
			    &B(i,j)), abs(d__2)) + ((d__3 = DU(
			    i-1).r, abs(d__3)) + (d__4 = d_imag(&DU(i - 1)), 
			    abs(d__4))) * ((d__5 = X(i-1,j).r, abs(d__5)) + (
			    d__6 = d_imag(&X(i-1,j)), abs(d__6)))
			     + ((d__7 = D(i).r, abs(d__7)) + (d__8 = 
			    d_imag(&D(i)), abs(d__8))) * ((d__9 = X(i,j).r, 
			    abs(d__9)) + (d__10 = d_imag(&X(i,j)), 
			    abs(d__10))) + ((d__11 = DL(i).r, abs(d__11)) 
			    + (d__12 = d_imag(&DL(i)), abs(d__12))) * ((d__13 
			    = X(i+1,j).r, abs(d__13)) + (d__14 = d_imag(&X(i+1,j)), abs(d__14)));
/* L40: */
		}
		i__2 = *n + j * b_dim1;
		i__3 = *n - 1;
		i__4 = *n - 1 + j * x_dim1;
		i__5 = *n;
		i__6 = *n + j * x_dim1;
		RWORK(*n) = (d__1 = B(*n,j).r, abs(d__1)) + (d__2 = d_imag(&B(*n,j)), abs(d__2)) + ((d__3 = DU(*n-1).r, 
			abs(d__3)) + (d__4 = d_imag(&DU(*n - 1)), abs(d__4))) 
			* ((d__5 = X(*n-1,j).r, abs(d__5)) + (d__6 = d_imag(&X(*n-1,j)), abs(d__6))) + ((d__7 = D(*n)
			.r, abs(d__7)) + (d__8 = d_imag(&D(*n)), abs(d__8))) *
			 ((d__9 = X(*n,j).r, abs(d__9)) + (d__10 = d_imag(&X(*n,j)), abs(d__10)));
	    }
	}

/*        Compute componentwise relative backward error from formula 
  

          max(i) ( abs(R(i)) / ( abs(op(A))*abs(X) + abs(B) )(i) )   

          where abs(Z) is the componentwise absolute value of the matr
ix   
          or vector Z.  If the i-th component of the denominator is le
ss   
          than SAFE2, then SAFE1 is added to the i-th components of th
e   
          numerator and denominator before dividing. */

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
/* L50: */
	}
	BERR(j) = s;

/*        Test stopping criterion. Continue iterating if   
             1) The residual BERR(J) is larger than machine epsilon, a
nd   
             2) BERR(J) decreased by at least a factor of 2 during the
   
                last iteration, and   
             3) At most ITMAX iterations tried. */

	if (BERR(j) > eps && BERR(j) * 2. <= lstres && count <= 5) {

/*           Update solution and try again. */

	    zgttrs_(trans, n, &c__1, &DLF(1), &DF(1), &DUF(1), &DU2(1), &IPIV(
		    1), &WORK(1), n, info);
	    zaxpy_(n, &c_b26, &WORK(1), &c__1, &X(1,j), &c__1);
	    lstres = BERR(j);
	    ++count;
	    goto L20;
	}

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
/* L60: */
	}

	kase = 0;
L70:
	zlacon_(n, &WORK(*n + 1), &WORK(1), &FERR(j), &kase);
	if (kase != 0) {
	    if (kase == 1) {

/*              Multiply by diag(W)*inv(op(A)**H). */

		zgttrs_(transt, n, &c__1, &DLF(1), &DF(1), &DUF(1), &DU2(1), &
			IPIV(1), &WORK(1), n, info);
		i__2 = *n;
		for (i = 1; i <= *n; ++i) {
		    i__3 = i;
		    i__4 = i;
		    i__5 = i;
		    z__1.r = RWORK(i) * WORK(i).r, z__1.i = RWORK(i) 
			    * WORK(i).i;
		    WORK(i).r = z__1.r, WORK(i).i = z__1.i;
/* L80: */
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
/* L90: */
		}
		zgttrs_(transn, n, &c__1, &DLF(1), &DF(1), &DUF(1), &DU2(1), &
			IPIV(1), &WORK(1), n, info);
	    }
	    goto L70;
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
/* L100: */
	}
	if (lstres != 0.) {
	    FERR(j) /= lstres;
	}

/* L110: */
    }

    return 0;

/*     End of ZGTRFS */

} /* zgtrfs_ */

