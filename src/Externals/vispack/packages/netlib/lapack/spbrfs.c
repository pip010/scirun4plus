#include "f2c.h"

/* Subroutine */ int spbrfs_(char *uplo, integer *n, integer *kd, integer *
	nrhs, real *ab, integer *ldab, real *afb, integer *ldafb, real *b, 
	integer *ldb, real *x, integer *ldx, real *ferr, real *berr, real *
	work, integer *iwork, integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    SPBRFS improves the computed solution to a system of linear   
    equations when the coefficient matrix is symmetric positive definite 
  
    and banded, and provides error bounds and backward error estimates   
    for the solution.   

    Arguments   
    =========   

    UPLO    (input) CHARACTER*1   
            = 'U':  Upper triangle of A is stored;   
            = 'L':  Lower triangle of A is stored.   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    KD      (input) INTEGER   
            The number of superdiagonals of the matrix A if UPLO = 'U',   
            or the number of subdiagonals if UPLO = 'L'.  KD >= 0.   

    NRHS    (input) INTEGER   
            The number of right hand sides, i.e., the number of columns   
            of the matrices B and X.  NRHS >= 0.   

    AB      (input) REAL array, dimension (LDAB,N)   
            The upper or lower triangle of the symmetric band matrix A,   
            stored in the first KD+1 rows of the array.  The j-th column 
  
            of A is stored in the j-th column of the array AB as follows: 
  
            if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j; 
  
            if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd). 
  

    LDAB    (input) INTEGER   
            The leading dimension of the array AB.  LDAB >= KD+1.   

    AFB     (input) REAL array, dimension (LDAFB,N)   
            The triangular factor U or L from the Cholesky factorization 
  
            A = U**T*U or A = L*L**T of the band matrix A as computed by 
  
            SPBTRF, in the same storage format as A (see AB).   

    LDAFB   (input) INTEGER   
            The leading dimension of the array AFB.  LDAFB >= KD+1.   

    B       (input) REAL array, dimension (LDB,NRHS)   
            The right hand side matrix B.   

    LDB     (input) INTEGER   
            The leading dimension of the array B.  LDB >= max(1,N).   

    X       (input/output) REAL array, dimension (LDX,NRHS)   
            On entry, the solution matrix X, as computed by SPBTRS.   
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
    static real c_b12 = -1.f;
    static real c_b14 = 1.f;
    
    /* System generated locals */
    integer ab_dim1, ab_offset, afb_dim1, afb_offset, b_dim1, b_offset, 
	    x_dim1, x_offset, i__1, i__2, i__3, i__4, i__5;
    real r__1, r__2, r__3;
    /* Local variables */
    static integer kase;
    static real safe1, safe2;
    static integer i, j, k, l;
    static real s;
    extern logical lsame_(char *, char *);
    static integer count;
    extern /* Subroutine */ int ssbmv_(char *, integer *, integer *, real *, 
	    real *, integer *, real *, integer *, real *, real *, integer *);
    static logical upper;
    extern /* Subroutine */ int scopy_(integer *, real *, integer *, real *, 
	    integer *), saxpy_(integer *, real *, real *, integer *, real *, 
	    integer *);
    static real xk;
    extern doublereal slamch_(char *);
    static integer nz;
    static real safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *), slacon_(
	    integer *, real *, real *, integer *, real *, integer *);
    static real lstres;
    extern /* Subroutine */ int spbtrs_(char *, integer *, integer *, integer 
	    *, real *, integer *, real *, integer *, integer *);
    static real eps;



#define FERR(I) ferr[(I)-1]
#define BERR(I) berr[(I)-1]
#define WORK(I) work[(I)-1]
#define IWORK(I) iwork[(I)-1]

#define AB(I,J) ab[(I)-1 + ((J)-1)* ( *ldab)]
#define AFB(I,J) afb[(I)-1 + ((J)-1)* ( *ldafb)]
#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]
#define X(I,J) x[(I)-1 + ((J)-1)* ( *ldx)]

    *info = 0;
    upper = lsame_(uplo, "U");
    if (! upper && ! lsame_(uplo, "L")) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*kd < 0) {
	*info = -3;
    } else if (*nrhs < 0) {
	*info = -4;
    } else if (*ldab < *kd + 1) {
	*info = -6;
    } else if (*ldafb < *kd + 1) {
	*info = -8;
    } else if (*ldb < max(1,*n)) {
	*info = -10;
    } else if (*ldx < max(1,*n)) {
	*info = -12;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("SPBRFS", &i__1);
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

/*     NZ = maximum number of nonzero elements in each row of A, plus 1   

   Computing MIN */
    i__1 = *n + 1, i__2 = (*kd << 1) + 2;
    nz = min(i__1,i__2);
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

          Compute residual R = B - A * X */

	scopy_(n, &B(1,j), &c__1, &WORK(*n + 1), &c__1);
	ssbmv_(uplo, n, kd, &c_b12, &AB(1,1), ldab, &X(1,j), 
		&c__1, &c_b14, &WORK(*n + 1), &c__1);

/*        Compute componentwise relative backward error from formula 
  

          max(i) ( abs(R(i)) / ( abs(A)*abs(X) + abs(B) )(i) )   

          where abs(Z) is the componentwise absolute value of the matr
ix   
          or vector Z.  If the i-th component of the denominator is le
ss   
          than SAFE2, then SAFE1 is added to the i-th components of th
e   
          numerator and denominator before dividing. */

	i__2 = *n;
	for (i = 1; i <= *n; ++i) {
	    WORK(i) = (r__1 = B(i,j), dabs(r__1));
/* L30: */
	}

/*        Compute abs(A)*abs(X) + abs(B). */

	if (upper) {
	    i__2 = *n;
	    for (k = 1; k <= *n; ++k) {
		s = 0.f;
		xk = (r__1 = X(k,j), dabs(r__1));
		l = *kd + 1 - k;
/* Computing MAX */
		i__3 = 1, i__4 = k - *kd;
		i__5 = k - 1;
		for (i = max(1,k-*kd); i <= k-1; ++i) {
		    WORK(i) += (r__1 = AB(l+i,k), dabs(r__1)) * 
			    xk;
		    s += (r__1 = AB(l+i,k), dabs(r__1)) * (r__2 
			    = X(i,j), dabs(r__2));
/* L40: */
		}
		WORK(k) = WORK(k) + (r__1 = AB(*kd+1,k), dabs(
			r__1)) * xk + s;
/* L50: */
	    }
	} else {
	    i__2 = *n;
	    for (k = 1; k <= *n; ++k) {
		s = 0.f;
		xk = (r__1 = X(k,j), dabs(r__1));
		WORK(k) += (r__1 = AB(1,k), dabs(r__1)) * xk;
		l = 1 - k;
/* Computing MIN */
		i__3 = *n, i__4 = k + *kd;
		i__5 = min(i__3,i__4);
		for (i = k + 1; i <= min(*n,k+*kd); ++i) {
		    WORK(i) += (r__1 = AB(l+i,k), dabs(r__1)) * 
			    xk;
		    s += (r__1 = AB(l+i,k), dabs(r__1)) * (r__2 
			    = X(i,j), dabs(r__2));
/* L60: */
		}
		WORK(k) += s;
/* L70: */
	    }
	}
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
/* L80: */
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

	    spbtrs_(uplo, n, kd, &c__1, &AFB(1,1), ldafb, &WORK(*n + 1)
		    , n, info);
	    saxpy_(n, &c_b14, &WORK(*n + 1), &c__1, &X(1,j), &c__1)
		    ;
	    lstres = BERR(j);
	    ++count;
	    goto L20;
	}

/*        Bound error from formula   

          norm(X - XTRUE) / norm(X) .le. FERR =   
          norm( abs(inv(A))*   
             ( abs(R) + NZ*EPS*( abs(A)*abs(X)+abs(B) ))) / norm(X)   

          where   
            norm(Z) is the magnitude of the largest component of Z   
            inv(A) is the inverse of A   
            abs(Z) is the componentwise absolute value of the matrix o
r   
               vector Z   
            NZ is the maximum number of nonzeros in any row of A, plus
 1   
            EPS is machine epsilon   

          The i-th component of abs(R)+NZ*EPS*(abs(A)*abs(X)+abs(B)) 
  
          is incremented by SAFE1 if the i-th component of   
          abs(A)*abs(X) + abs(B) is less than SAFE2.   

          Use SLACON to estimate the infinity-norm of the matrix   
             inv(A) * diag(W),   
          where W = abs(R) + NZ*EPS*( abs(A)*abs(X)+abs(B) ))) */

	i__2 = *n;
	for (i = 1; i <= *n; ++i) {
	    if (WORK(i) > safe2) {
		WORK(i) = (r__1 = WORK(*n + i), dabs(r__1)) + nz * eps * WORK(
			i);
	    } else {
		WORK(i) = (r__1 = WORK(*n + i), dabs(r__1)) + nz * eps * WORK(
			i) + safe1;
	    }
/* L90: */
	}

	kase = 0;
L100:
	slacon_(n, &WORK((*n << 1) + 1), &WORK(*n + 1), &IWORK(1), &FERR(j), &
		kase);
	if (kase != 0) {
	    if (kase == 1) {

/*              Multiply by diag(W)*inv(A'). */

		spbtrs_(uplo, n, kd, &c__1, &AFB(1,1), ldafb, &WORK(*n 
			+ 1), n, info);
		i__2 = *n;
		for (i = 1; i <= *n; ++i) {
		    WORK(*n + i) *= WORK(i);
/* L110: */
		}
	    } else if (kase == 2) {

/*              Multiply by inv(A)*diag(W). */

		i__2 = *n;
		for (i = 1; i <= *n; ++i) {
		    WORK(*n + i) *= WORK(i);
/* L120: */
		}
		spbtrs_(uplo, n, kd, &c__1, &AFB(1,1), ldafb, &WORK(*n 
			+ 1), n, info);
	    }
	    goto L100;
	}

/*        Normalize error. */

	lstres = 0.f;
	i__2 = *n;
	for (i = 1; i <= *n; ++i) {
/* Computing MAX */
	    r__2 = lstres, r__3 = (r__1 = X(i,j), dabs(r__1));
	    lstres = dmax(r__2,r__3);
/* L130: */
	}
	if (lstres != 0.f) {
	    FERR(j) /= lstres;
	}

/* L140: */
    }

    return 0;

/*     End of SPBRFS */

} /* spbrfs_ */

