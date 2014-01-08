#include "f2c.h"

/* Subroutine */ int sgtrfs_(char *trans, integer *n, integer *nrhs, real *dl,
	 real *d, real *du, real *dlf, real *df, real *duf, real *du2, 
	integer *ipiv, real *b, integer *ldb, real *x, integer *ldx, real *
	ferr, real *berr, real *work, integer *iwork, integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    SGTRFS improves the computed solution to a system of linear   
    equations when the coefficient matrix is tridiagonal, and provides   
    error bounds and backward error estimates for the solution.   

    Arguments   
    =========   

    TRANS   (input) CHARACTER*1   
            Specifies the form of the system of equations:   
            = 'N':  A * X = B     (No transpose)   
            = 'T':  A**T * X = B  (Transpose)   
            = 'C':  A**H * X = B  (Conjugate transpose = Transpose)   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    NRHS    (input) INTEGER   
            The number of right hand sides, i.e., the number of columns   
            of the matrix B.  NRHS >= 0.   

    DL      (input) REAL array, dimension (N-1)   
            The (n-1) subdiagonal elements of A.   

    D       (input) REAL array, dimension (N)   
            The diagonal elements of A.   

    DU      (input) REAL array, dimension (N-1)   
            The (n-1) superdiagonal elements of A.   

    DLF     (input) REAL array, dimension (N-1)   
            The (n-1) multipliers that define the matrix L from the   
            LU factorization of A as computed by SGTTRF.   

    DF      (input) REAL array, dimension (N)   
            The n diagonal elements of the upper triangular matrix U from 
  
            the LU factorization of A.   

    DUF     (input) REAL array, dimension (N-1)   
            The (n-1) elements of the first superdiagonal of U.   

    DU2     (input) REAL array, dimension (N-2)   
            The (n-2) elements of the second superdiagonal of U.   

    IPIV    (input) INTEGER array, dimension (N)   
            The pivot indices; for 1 <= i <= n, row i of the matrix was   
            interchanged with row IPIV(i).  IPIV(i) will always be either 
  
            i or i+1; IPIV(i) = i indicates a row interchange was not   
            required.   

    B       (input) REAL array, dimension (LDB,NRHS)   
            The right hand side matrix B.   

    LDB     (input) INTEGER   
            The leading dimension of the array B.  LDB >= max(1,N).   

    X       (input/output) REAL array, dimension (LDX,NRHS)   
            On entry, the solution matrix X, as computed by SGTTRS.   
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

    WORK    (workspace) REAL array, dimension (3*N)   

    IWORK   (workspace) INTEGER array, dimension (N)   

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
    
    /* System generated locals */
    integer b_dim1, b_offset, x_dim1, x_offset, i__1, i__2;
    real r__1, r__2, r__3, r__4;
    /* Local variables */
    static integer kase;
    static real safe1, safe2;
    static integer i, j;
    static real s;
    extern logical lsame_(char *, char *);
    static integer count;
    extern /* Subroutine */ int scopy_(integer *, real *, integer *, real *, 
	    integer *), saxpy_(integer *, real *, real *, integer *, real *, 
	    integer *);
    extern doublereal slamch_(char *);
    static integer nz;
    static real safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *), slacon_(
	    integer *, real *, real *, integer *, real *, integer *), slagtm_(
	    char *, integer *, integer *, real *, real *, real *, real *, 
	    real *, integer *, real *, real *, integer *);
    static logical notran;
    static char transn[1], transt[1];
    static real lstres;
    extern /* Subroutine */ int sgttrs_(char *, integer *, integer *, real *, 
	    real *, real *, real *, integer *, real *, integer *, integer *);
    static real eps;



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
#define IWORK(I) iwork[(I)-1]

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
	xerbla_("SGTRFS", &i__1);
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
	*(unsigned char *)transt = 'T';
    } else {
	*(unsigned char *)transn = 'T';
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

	scopy_(n, &B(1,j), &c__1, &WORK(*n + 1), &c__1);
	slagtm_(trans, n, &c__1, &c_b18, &DL(1), &D(1), &DU(1), &X(1,j), ldx, &c_b19, &WORK(*n + 1), n);

/*        Compute abs(op(A))*abs(x) + abs(b) for use in the backward 
  
          error bound. */

	if (notran) {
	    if (*n == 1) {
		WORK(1) = (r__1 = B(1,j), dabs(r__1)) + (r__2 = D(
			1) * X(1,j), dabs(r__2));
	    } else {
		WORK(1) = (r__1 = B(1,j), dabs(r__1)) + (r__2 = D(
			1) * X(1,j), dabs(r__2)) + (r__3 = DU(1) * 
			X(2,j), dabs(r__3));
		i__2 = *n - 1;
		for (i = 2; i <= *n-1; ++i) {
		    WORK(i) = (r__1 = B(i,j), dabs(r__1)) + (r__2 =
			     DL(i - 1) * X(i-1,j), dabs(r__2)) + 
			    (r__3 = D(i) * X(i,j), dabs(r__3)) + (
			    r__4 = DU(i) * X(i+1,j), dabs(r__4));
/* L30: */
		}
		WORK(*n) = (r__1 = B(*n,j), dabs(r__1)) + (r__2 = 
			DL(*n - 1) * X(*n-1,j), dabs(r__2)) + (
			r__3 = D(*n) * X(*n,j), dabs(r__3));
	    }
	} else {
	    if (*n == 1) {
		WORK(1) = (r__1 = B(1,j), dabs(r__1)) + (r__2 = D(
			1) * X(1,j), dabs(r__2));
	    } else {
		WORK(1) = (r__1 = B(1,j), dabs(r__1)) + (r__2 = D(
			1) * X(1,j), dabs(r__2)) + (r__3 = DL(1) * 
			X(2,j), dabs(r__3));
		i__2 = *n - 1;
		for (i = 2; i <= *n-1; ++i) {
		    WORK(i) = (r__1 = B(i,j), dabs(r__1)) + (r__2 =
			     DU(i - 1) * X(i-1,j), dabs(r__2)) + 
			    (r__3 = D(i) * X(i,j), dabs(r__3)) + (
			    r__4 = DL(i) * X(i+1,j), dabs(r__4));
/* L40: */
		}
		WORK(*n) = (r__1 = B(*n,j), dabs(r__1)) + (r__2 = 
			DU(*n - 1) * X(*n-1,j), dabs(r__2)) + (
			r__3 = D(*n) * X(*n,j), dabs(r__3));
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
	    if (WORK(i) > safe2) {
/* Computing MAX */
		r__2 = s, r__3 = (r__1 = WORK(*n + i), dabs(r__1)) / WORK(i);
		s = dmax(r__2,r__3);
	    } else {
/* Computing MAX */
		r__2 = s, r__3 = ((r__1 = WORK(*n + i), dabs(r__1)) + safe1) /
			 (WORK(i) + safe1);
		s = dmax(r__2,r__3);
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

	    sgttrs_(trans, n, &c__1, &DLF(1), &DF(1), &DUF(1), &DU2(1), &IPIV(
		    1), &WORK(*n + 1), n, info);
	    saxpy_(n, &c_b19, &WORK(*n + 1), &c__1, &X(1,j), &c__1)
		    ;
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

          Use SLACON to estimate the infinity-norm of the matrix   
             inv(op(A)) * diag(W),   
          where W = abs(R) + NZ*EPS*( abs(op(A))*abs(X)+abs(B) ))) */

	i__2 = *n;
	for (i = 1; i <= *n; ++i) {
	    if (WORK(i) > safe2) {
		WORK(i) = (r__1 = WORK(*n + i), dabs(r__1)) + nz * eps * WORK(
			i);
	    } else {
		WORK(i) = (r__1 = WORK(*n + i), dabs(r__1)) + nz * eps * WORK(
			i) + safe1;
	    }
/* L60: */
	}

	kase = 0;
L70:
	slacon_(n, &WORK((*n << 1) + 1), &WORK(*n + 1), &IWORK(1), &FERR(j), &
		kase);
	if (kase != 0) {
	    if (kase == 1) {

/*              Multiply by diag(W)*inv(op(A)**T). */

		sgttrs_(transt, n, &c__1, &DLF(1), &DF(1), &DUF(1), &DU2(1), &
			IPIV(1), &WORK(*n + 1), n, info);
		i__2 = *n;
		for (i = 1; i <= *n; ++i) {
		    WORK(*n + i) = WORK(i) * WORK(*n + i);
/* L80: */
		}
	    } else {

/*              Multiply by inv(op(A))*diag(W). */

		i__2 = *n;
		for (i = 1; i <= *n; ++i) {
		    WORK(*n + i) = WORK(i) * WORK(*n + i);
/* L90: */
		}
		sgttrs_(transn, n, &c__1, &DLF(1), &DF(1), &DUF(1), &DU2(1), &
			IPIV(1), &WORK(*n + 1), n, info);
	    }
	    goto L70;
	}

/*        Normalize error. */

	lstres = 0.f;
	i__2 = *n;
	for (i = 1; i <= *n; ++i) {
/* Computing MAX */
	    r__2 = lstres, r__3 = (r__1 = X(i,j), dabs(r__1));
	    lstres = dmax(r__2,r__3);
/* L100: */
	}
	if (lstres != 0.f) {
	    FERR(j) /= lstres;
	}

/* L110: */
    }

    return 0;

/*     End of SGTRFS */

} /* sgtrfs_ */

