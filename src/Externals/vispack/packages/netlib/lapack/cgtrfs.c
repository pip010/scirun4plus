#include "f2c.h"

/* Subroutine */ int cgtrfs_(char *trans, integer *n, integer *nrhs, complex *
	dl, complex *d, complex *du, complex *dlf, complex *df, complex *duf, 
	complex *du2, integer *ipiv, complex *b, integer *ldb, complex *x, 
	integer *ldx, real *ferr, real *berr, complex *work, real *rwork, 
	integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    CGTRFS improves the computed solution to a system of linear   
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

    DL      (input) COMPLEX array, dimension (N-1)   
            The (n-1) subdiagonal elements of A.   

    D       (input) COMPLEX array, dimension (N)   
            The diagonal elements of A.   

    DU      (input) COMPLEX array, dimension (N-1)   
            The (n-1) superdiagonal elements of A.   

    DLF     (input) COMPLEX array, dimension (N-1)   
            The (n-1) multipliers that define the matrix L from the   
            LU factorization of A as computed by CGTTRF.   

    DF      (input) COMPLEX array, dimension (N)   
            The n diagonal elements of the upper triangular matrix U from 
  
            the LU factorization of A.   

    DUF     (input) COMPLEX array, dimension (N-1)   
            The (n-1) elements of the first superdiagonal of U.   

    DU2     (input) COMPLEX array, dimension (N-2)   
            The (n-2) elements of the second superdiagonal of U.   

    IPIV    (input) INTEGER array, dimension (N)   
            The pivot indices; for 1 <= i <= n, row i of the matrix was   
            interchanged with row IPIV(i).  IPIV(i) will always be either 
  
            i or i+1; IPIV(i) = i indicates a row interchange was not   
            required.   

    B       (input) COMPLEX array, dimension (LDB,NRHS)   
            The right hand side matrix B.   

    LDB     (input) INTEGER   
            The leading dimension of the array B.  LDB >= max(1,N).   

    X       (input/output) COMPLEX array, dimension (LDX,NRHS)   
            On entry, the solution matrix X, as computed by CGTTRS.   
            On exit, the improved solution matrix X.   

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

    Internal Parameters   
    ===================   

    ITMAX is the maximum number of steps of iterative refinement.   

    ===================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static integer c__1 = 1;
    static real c_b18 = -1.f;
    static real c_b19 = 1.f;
    static complex c_b26 = {1.f,0.f};
    
    /* System generated locals */
    integer b_dim1, b_offset, x_dim1, x_offset, i__1, i__2, i__3, i__4, i__5, 
	    i__6, i__7, i__8, i__9;
    real r__1, r__2, r__3, r__4, r__5, r__6, r__7, r__8, r__9, r__10, r__11, 
	    r__12, r__13, r__14;
    complex q__1;
    /* Builtin functions */
    double r_imag(complex *);
    /* Local variables */
    static integer kase;
    static real safe1, safe2;
    static integer i, j;
    static real s;
    extern logical lsame_(char *, char *);
    extern /* Subroutine */ int ccopy_(integer *, complex *, integer *, 
	    complex *, integer *), caxpy_(integer *, complex *, complex *, 
	    integer *, complex *, integer *);
    static integer count;
    extern /* Subroutine */ int clacon_(integer *, complex *, complex *, real 
	    *, integer *), clagtm_(char *, integer *, integer *, real *, 
	    complex *, complex *, complex *, complex *, integer *, real *, 
	    complex *, integer *);
    static integer nz;
    extern doublereal slamch_(char *);
    static real safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *);
    static logical notran;
    static char transn[1];
    extern /* Subroutine */ int cgttrs_(char *, integer *, integer *, complex 
	    *, complex *, complex *, complex *, integer *, complex *, integer 
	    *, integer *);
    static char transt[1];
    static real lstres, eps;



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
	xerbla_("CGTRFS", &i__1);
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

    nz = 4;
    eps = slamch_("Epsilon");
    safmin = slamch_("Safe minimum");
    safe1 = nz * safmin;
    safe2 = safe1 / eps;

/*     Do for each right hand side */

    i__1 = *nrhs;
    for (j = 1; j <= *nrhs; ++j) {

	count = 1;
	lstres = 3.f;
L20:

/*        Loop until stopping criterion is satisfied.   

          Compute residual R = B - op(A) * X,   
          where op(A) = A, A**T, or A**H, depending on TRANS. */

	ccopy_(n, &B(1,j), &c__1, &WORK(1), &c__1);
	clagtm_(trans, n, &c__1, &c_b18, &DL(1), &D(1), &DU(1), &X(1,j), ldx, &c_b19, &WORK(1), n);

/*        Compute abs(op(A))*abs(x) + abs(b) for use in the backward 
  
          error bound. */

	if (notran) {
	    if (*n == 1) {
		i__2 = j * b_dim1 + 1;
		i__3 = j * x_dim1 + 1;
		RWORK(1) = (r__1 = B(1,j).r, dabs(r__1)) + (r__2 = r_imag(&B(1,j)), dabs(r__2)) + ((r__3 = D(1).r, dabs(
			r__3)) + (r__4 = r_imag(&D(1)), dabs(r__4))) * ((r__5 
			= X(1,j).r, dabs(r__5)) + (r__6 = r_imag(&X(1,j)), dabs(r__6)));
	    } else {
		i__2 = j * b_dim1 + 1;
		i__3 = j * x_dim1 + 1;
		i__4 = j * x_dim1 + 2;
		RWORK(1) = (r__1 = B(1,j).r, dabs(r__1)) + (r__2 = r_imag(&B(1,j)), dabs(r__2)) + ((r__3 = D(1).r, dabs(
			r__3)) + (r__4 = r_imag(&D(1)), dabs(r__4))) * ((r__5 
			= X(1,j).r, dabs(r__5)) + (r__6 = r_imag(&X(1,j)), dabs(r__6))) + ((r__7 = DU(1).r, dabs(
			r__7)) + (r__8 = r_imag(&DU(1)), dabs(r__8))) * ((
			r__9 = X(2,j).r, dabs(r__9)) + (r__10 = r_imag(&X(2,j)), dabs(r__10)));
		i__2 = *n - 1;
		for (i = 2; i <= *n-1; ++i) {
		    i__3 = i + j * b_dim1;
		    i__4 = i - 1;
		    i__5 = i - 1 + j * x_dim1;
		    i__6 = i;
		    i__7 = i + j * x_dim1;
		    i__8 = i;
		    i__9 = i + 1 + j * x_dim1;
		    RWORK(i) = (r__1 = B(i,j).r, dabs(r__1)) + (r__2 = 
			    r_imag(&B(i,j)), dabs(r__2)) + ((r__3 =
			     DL(i-1).r, dabs(r__3)) + (r__4 = r_imag(&DL(i - 
			    1)), dabs(r__4))) * ((r__5 = X(i-1,j).r, dabs(r__5)
			    ) + (r__6 = r_imag(&X(i-1,j)), dabs(
			    r__6))) + ((r__7 = D(i).r, dabs(r__7)) + (r__8 
			    = r_imag(&D(i)), dabs(r__8))) * ((r__9 = X(i,j)
			    .r, dabs(r__9)) + (r__10 = r_imag(&X(i,j)), dabs(r__10))) + ((r__11 = DU(i).r, 
			    dabs(r__11)) + (r__12 = r_imag(&DU(i)), dabs(
			    r__12))) * ((r__13 = X(i+1,j).r, dabs(r__13)) + (
			    r__14 = r_imag(&X(i+1,j)), dabs(
			    r__14)));
/* L30: */
		}
		i__2 = *n + j * b_dim1;
		i__3 = *n - 1;
		i__4 = *n - 1 + j * x_dim1;
		i__5 = *n;
		i__6 = *n + j * x_dim1;
		RWORK(*n) = (r__1 = B(*n,j).r, dabs(r__1)) + (r__2 = r_imag(&
			B(*n,j)), dabs(r__2)) + ((r__3 = DL(*n-1)
			.r, dabs(r__3)) + (r__4 = r_imag(&DL(*n - 1)), dabs(
			r__4))) * ((r__5 = X(*n-1,j).r, dabs(r__5)) + (r__6 = 
			r_imag(&X(*n-1,j)), dabs(r__6))) + ((
			r__7 = D(*n).r, dabs(r__7)) + (r__8 = r_imag(&D(*n))
			, dabs(r__8))) * ((r__9 = X(*n,j).r, dabs(r__9)) + (
			r__10 = r_imag(&X(*n,j)), dabs(r__10)));
	    }
	} else {
	    if (*n == 1) {
		i__2 = j * b_dim1 + 1;
		i__3 = j * x_dim1 + 1;
		RWORK(1) = (r__1 = B(1,j).r, dabs(r__1)) + (r__2 = r_imag(&B(1,j)), dabs(r__2)) + ((r__3 = D(1).r, dabs(
			r__3)) + (r__4 = r_imag(&D(1)), dabs(r__4))) * ((r__5 
			= X(1,j).r, dabs(r__5)) + (r__6 = r_imag(&X(1,j)), dabs(r__6)));
	    } else {
		i__2 = j * b_dim1 + 1;
		i__3 = j * x_dim1 + 1;
		i__4 = j * x_dim1 + 2;
		RWORK(1) = (r__1 = B(1,j).r, dabs(r__1)) + (r__2 = r_imag(&B(1,j)), dabs(r__2)) + ((r__3 = D(1).r, dabs(
			r__3)) + (r__4 = r_imag(&D(1)), dabs(r__4))) * ((r__5 
			= X(1,j).r, dabs(r__5)) + (r__6 = r_imag(&X(1,j)), dabs(r__6))) + ((r__7 = DL(1).r, dabs(
			r__7)) + (r__8 = r_imag(&DL(1)), dabs(r__8))) * ((
			r__9 = X(2,j).r, dabs(r__9)) + (r__10 = r_imag(&X(2,j)), dabs(r__10)));
		i__2 = *n - 1;
		for (i = 2; i <= *n-1; ++i) {
		    i__3 = i + j * b_dim1;
		    i__4 = i - 1;
		    i__5 = i - 1 + j * x_dim1;
		    i__6 = i;
		    i__7 = i + j * x_dim1;
		    i__8 = i;
		    i__9 = i + 1 + j * x_dim1;
		    RWORK(i) = (r__1 = B(i,j).r, dabs(r__1)) + (r__2 = 
			    r_imag(&B(i,j)), dabs(r__2)) + ((r__3 =
			     DU(i-1).r, dabs(r__3)) + (r__4 = r_imag(&DU(i - 
			    1)), dabs(r__4))) * ((r__5 = X(i-1,j).r, dabs(r__5)
			    ) + (r__6 = r_imag(&X(i-1,j)), dabs(
			    r__6))) + ((r__7 = D(i).r, dabs(r__7)) + (r__8 
			    = r_imag(&D(i)), dabs(r__8))) * ((r__9 = X(i,j)
			    .r, dabs(r__9)) + (r__10 = r_imag(&X(i,j)), dabs(r__10))) + ((r__11 = DL(i).r, 
			    dabs(r__11)) + (r__12 = r_imag(&DL(i)), dabs(
			    r__12))) * ((r__13 = X(i+1,j).r, dabs(r__13)) + (
			    r__14 = r_imag(&X(i+1,j)), dabs(
			    r__14)));
/* L40: */
		}
		i__2 = *n + j * b_dim1;
		i__3 = *n - 1;
		i__4 = *n - 1 + j * x_dim1;
		i__5 = *n;
		i__6 = *n + j * x_dim1;
		RWORK(*n) = (r__1 = B(*n,j).r, dabs(r__1)) + (r__2 = r_imag(&
			B(*n,j)), dabs(r__2)) + ((r__3 = DU(*n-1)
			.r, dabs(r__3)) + (r__4 = r_imag(&DU(*n - 1)), dabs(
			r__4))) * ((r__5 = X(*n-1,j).r, dabs(r__5)) + (r__6 = 
			r_imag(&X(*n-1,j)), dabs(r__6))) + ((
			r__7 = D(*n).r, dabs(r__7)) + (r__8 = r_imag(&D(*n))
			, dabs(r__8))) * ((r__9 = X(*n,j).r, dabs(r__9)) + (
			r__10 = r_imag(&X(*n,j)), dabs(r__10)));
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
/* L50: */
	}
	BERR(j) = s;

/*        Test stopping criterion. Continue iterating if   
             1) The residual BERR(J) is larger than machine epsilon, a
nd   
             2) BERR(J) decreased by at least a factor of 2 during the
   
                last iteration, and   
             3) At most ITMAX iterations tried. */

	if (BERR(j) > eps && BERR(j) * 2.f <= lstres && count <= 5) {

/*           Update solution and try again. */

	    cgttrs_(trans, n, &c__1, &DLF(1), &DF(1), &DUF(1), &DU2(1), &IPIV(
		    1), &WORK(1), n, info);
	    caxpy_(n, &c_b26, &WORK(1), &c__1, &X(1,j), &c__1);
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
/* L60: */
	}

	kase = 0;
L70:
	clacon_(n, &WORK(*n + 1), &WORK(1), &FERR(j), &kase);
	if (kase != 0) {
	    if (kase == 1) {

/*              Multiply by diag(W)*inv(op(A)**H). */

		cgttrs_(transt, n, &c__1, &DLF(1), &DF(1), &DUF(1), &DU2(1), &
			IPIV(1), &WORK(1), n, info);
		i__2 = *n;
		for (i = 1; i <= *n; ++i) {
		    i__3 = i;
		    i__4 = i;
		    i__5 = i;
		    q__1.r = RWORK(i) * WORK(i).r, q__1.i = RWORK(i) 
			    * WORK(i).i;
		    WORK(i).r = q__1.r, WORK(i).i = q__1.i;
/* L80: */
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
/* L90: */
		}
		cgttrs_(transn, n, &c__1, &DLF(1), &DF(1), &DUF(1), &DU2(1), &
			IPIV(1), &WORK(1), n, info);
	    }
	    goto L70;
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
/* L100: */
	}
	if (lstres != 0.f) {
	    FERR(j) /= lstres;
	}

/* L110: */
    }

    return 0;

/*     End of CGTRFS */

} /* cgtrfs_ */

