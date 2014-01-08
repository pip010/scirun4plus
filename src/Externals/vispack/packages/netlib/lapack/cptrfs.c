#include "f2c.h"

/* Subroutine */ int cptrfs_(char *uplo, integer *n, integer *nrhs, real *d, 
	complex *e, real *df, complex *ef, complex *b, integer *ldb, complex *
	x, integer *ldx, real *ferr, real *berr, complex *work, real *rwork, 
	integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    CPTRFS improves the computed solution to a system of linear   
    equations when the coefficient matrix is Hermitian positive definite 
  
    and tridiagonal, and provides error bounds and backward error   
    estimates for the solution.   

    Arguments   
    =========   

    UPLO    (input) CHARACTER*1   
            Specifies whether the superdiagonal or the subdiagonal of the 
  
            tridiagonal matrix A is stored and the form of the   
            factorization:   
            = 'U':  E is the superdiagonal of A, and A = U**H*D*U;   
            = 'L':  E is the subdiagonal of A, and A = L*D*L**H.   
            (The two forms are equivalent if A is real.)   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    NRHS    (input) INTEGER   
            The number of right hand sides, i.e., the number of columns   
            of the matrix B.  NRHS >= 0.   

    D       (input) REAL array, dimension (N)   
            The n real diagonal elements of the tridiagonal matrix A.   

    E       (input) COMPLEX array, dimension (N-1)   
            The (n-1) off-diagonal elements of the tridiagonal matrix A   
            (see UPLO).   

    DF      (input) REAL array, dimension (N)   
            The n diagonal elements of the diagonal matrix D from   
            the factorization computed by CPTTRF.   

    EF      (input) COMPLEX array, dimension (N-1)   
            The (n-1) off-diagonal elements of the unit bidiagonal   
            factor U or L from the factorization computed by CPTTRF   
            (see UPLO).   

    B       (input) COMPLEX array, dimension (LDB,NRHS)   
            The right hand side matrix B.   

    LDB     (input) INTEGER   
            The leading dimension of the array B.  LDB >= max(1,N).   

    X       (input/output) COMPLEX array, dimension (LDX,NRHS)   
            On entry, the solution matrix X, as computed by CPTTRS.   
            On exit, the improved solution matrix X.   

    LDX     (input) INTEGER   
            The leading dimension of the array X.  LDX >= max(1,N).   

    FERR    (output) REAL array, dimension (NRHS)   
            The forward error bound for each solution vector   
            X(j) (the j-th column of the solution matrix X).   
            If XTRUE is the true solution corresponding to X(j), FERR(j) 
  
            is an estimated upper bound for the magnitude of the largest 
  
            element in (X(j) - XTRUE) divided by the magnitude of the   
            largest element in X(j).   

    BERR    (output) REAL array, dimension (NRHS)   
            The componentwise relative backward error of each solution   
            vector X(j) (i.e., the smallest relative change in   
            any element of A or B that makes X(j) an exact solution).   

    WORK    (workspace) COMPLEX array, dimension (N)   

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
    static complex c_b16 = {1.f,0.f};
    
    /* System generated locals */
    integer b_dim1, b_offset, x_dim1, x_offset, i__1, i__2, i__3, i__4, i__5, 
	    i__6;
    real r__1, r__2, r__3, r__4, r__5, r__6, r__7, r__8, r__9, r__10, r__11, 
	    r__12;
    complex q__1, q__2, q__3;
    /* Builtin functions */
    double r_imag(complex *);
    void r_cnjg(complex *, complex *);
    double c_abs(complex *);
    /* Local variables */
    static real safe1, safe2;
    static integer i, j;
    static real s;
    extern logical lsame_(char *, char *);
    extern /* Subroutine */ int caxpy_(integer *, complex *, complex *, 
	    integer *, complex *, integer *);
    static integer count;
    static logical upper;
    static complex bi, cx, dx, ex;
    static integer ix;
    extern doublereal slamch_(char *);
    static integer nz;
    static real safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *);
    extern integer isamax_(integer *, real *, integer *);
    static real lstres;
    extern /* Subroutine */ int cpttrs_(char *, integer *, integer *, real *, 
	    complex *, complex *, integer *, integer *);
    static real eps;



#define D(I) d[(I)-1]
#define E(I) e[(I)-1]
#define DF(I) df[(I)-1]
#define EF(I) ef[(I)-1]
#define FERR(I) ferr[(I)-1]
#define BERR(I) berr[(I)-1]
#define WORK(I) work[(I)-1]
#define RWORK(I) rwork[(I)-1]

#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]
#define X(I,J) x[(I)-1 + ((J)-1)* ( *ldx)]

    *info = 0;
    upper = lsame_(uplo, "U");
    if (! upper && ! lsame_(uplo, "L")) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*nrhs < 0) {
	*info = -3;
    } else if (*ldb < max(1,*n)) {
	*info = -9;
    } else if (*ldx < max(1,*n)) {
	*info = -11;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("CPTRFS", &i__1);
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

          Compute residual R = B - A * X.  Also compute   
          abs(A)*abs(x) + abs(b) for use in the backward error bound. 
*/

	if (upper) {
	    if (*n == 1) {
		i__2 = j * b_dim1 + 1;
		bi.r = B(1,j).r, bi.i = B(1,j).i;
		i__2 = j * x_dim1 + 1;
		q__1.r = D(1) * X(1,j).r, q__1.i = D(1) * X(1,j).i;
		dx.r = q__1.r, dx.i = q__1.i;
		q__1.r = bi.r - dx.r, q__1.i = bi.i - dx.i;
		WORK(1).r = q__1.r, WORK(1).i = q__1.i;
		RWORK(1) = (r__1 = bi.r, dabs(r__1)) + (r__2 = r_imag(&bi), 
			dabs(r__2)) + ((r__3 = dx.r, dabs(r__3)) + (r__4 = 
			r_imag(&dx), dabs(r__4)));
	    } else {
		i__2 = j * b_dim1 + 1;
		bi.r = B(1,j).r, bi.i = B(1,j).i;
		i__2 = j * x_dim1 + 1;
		q__1.r = D(1) * X(1,j).r, q__1.i = D(1) * X(1,j).i;
		dx.r = q__1.r, dx.i = q__1.i;
		i__2 = j * x_dim1 + 2;
		q__1.r = E(1).r * X(2,j).r - E(1).i * X(2,j).i, q__1.i = E(
			1).r * X(2,j).i + E(1).i * X(2,j).r;
		ex.r = q__1.r, ex.i = q__1.i;
		q__2.r = bi.r - dx.r, q__2.i = bi.i - dx.i;
		q__1.r = q__2.r - ex.r, q__1.i = q__2.i - ex.i;
		WORK(1).r = q__1.r, WORK(1).i = q__1.i;
		i__2 = j * x_dim1 + 2;
		RWORK(1) = (r__1 = bi.r, dabs(r__1)) + (r__2 = r_imag(&bi), 
			dabs(r__2)) + ((r__3 = dx.r, dabs(r__3)) + (r__4 = 
			r_imag(&dx), dabs(r__4))) + ((r__5 = E(1).r, dabs(
			r__5)) + (r__6 = r_imag(&E(1)), dabs(r__6))) * ((r__7 
			= X(2,j).r, dabs(r__7)) + (r__8 = r_imag(&X(2,j)), dabs(r__8)));
		i__2 = *n - 1;
		for (i = 2; i <= *n-1; ++i) {
		    i__3 = i + j * b_dim1;
		    bi.r = B(i,j).r, bi.i = B(i,j).i;
		    r_cnjg(&q__2, &E(i - 1));
		    i__3 = i - 1 + j * x_dim1;
		    q__1.r = q__2.r * X(i-1,j).r - q__2.i * X(i-1,j).i, q__1.i =
			     q__2.r * X(i-1,j).i + q__2.i * X(i-1,j).r;
		    cx.r = q__1.r, cx.i = q__1.i;
		    i__3 = i;
		    i__4 = i + j * x_dim1;
		    q__1.r = D(i) * X(i,j).r, q__1.i = D(i) * X(i,j)
			    .i;
		    dx.r = q__1.r, dx.i = q__1.i;
		    i__3 = i;
		    i__4 = i + 1 + j * x_dim1;
		    q__1.r = E(i).r * X(i+1,j).r - E(i).i * X(i+1,j).i, 
			    q__1.i = E(i).r * X(i+1,j).i + E(i).i * X(i+1,j).r;
		    ex.r = q__1.r, ex.i = q__1.i;
		    i__3 = i;
		    q__3.r = bi.r - cx.r, q__3.i = bi.i - cx.i;
		    q__2.r = q__3.r - dx.r, q__2.i = q__3.i - dx.i;
		    q__1.r = q__2.r - ex.r, q__1.i = q__2.i - ex.i;
		    WORK(i).r = q__1.r, WORK(i).i = q__1.i;
		    i__3 = i - 1;
		    i__4 = i - 1 + j * x_dim1;
		    i__5 = i;
		    i__6 = i + 1 + j * x_dim1;
		    RWORK(i) = (r__1 = bi.r, dabs(r__1)) + (r__2 = r_imag(&bi)
			    , dabs(r__2)) + ((r__3 = E(i-1).r, dabs(r__3)) + 
			    (r__4 = r_imag(&E(i - 1)), dabs(r__4))) * ((r__5 =
			     X(i-1,j).r, dabs(r__5)) + (r__6 = r_imag(&X(i-1,j)), dabs(r__6))) + ((r__7 = dx.r, 
			    dabs(r__7)) + (r__8 = r_imag(&dx), dabs(r__8))) + 
			    ((r__9 = E(i).r, dabs(r__9)) + (r__10 = r_imag(
			    &E(i)), dabs(r__10))) * ((r__11 = X(i+1,j).r, dabs(
			    r__11)) + (r__12 = r_imag(&X(i+1,j)),
			     dabs(r__12)));
/* L30: */
		}
		i__2 = *n + j * b_dim1;
		bi.r = B(*n,j).r, bi.i = B(*n,j).i;
		r_cnjg(&q__2, &E(*n - 1));
		i__2 = *n - 1 + j * x_dim1;
		q__1.r = q__2.r * X(*n-1,j).r - q__2.i * X(*n-1,j).i, q__1.i = 
			q__2.r * X(*n-1,j).i + q__2.i * X(*n-1,j).r;
		cx.r = q__1.r, cx.i = q__1.i;
		i__2 = *n;
		i__3 = *n + j * x_dim1;
		q__1.r = D(*n) * X(*n,j).r, q__1.i = D(*n) * X(*n,j).i;
		dx.r = q__1.r, dx.i = q__1.i;
		i__2 = *n;
		q__2.r = bi.r - cx.r, q__2.i = bi.i - cx.i;
		q__1.r = q__2.r - dx.r, q__1.i = q__2.i - dx.i;
		WORK(*n).r = q__1.r, WORK(*n).i = q__1.i;
		i__2 = *n - 1;
		i__3 = *n - 1 + j * x_dim1;
		RWORK(*n) = (r__1 = bi.r, dabs(r__1)) + (r__2 = r_imag(&bi), 
			dabs(r__2)) + ((r__3 = E(*n-1).r, dabs(r__3)) + (r__4 
			= r_imag(&E(*n - 1)), dabs(r__4))) * ((r__5 = X(*n-1,j)
			.r, dabs(r__5)) + (r__6 = r_imag(&X(*n-1,j)), dabs(r__6))) + ((r__7 = dx.r, dabs(r__7)) + 
			(r__8 = r_imag(&dx), dabs(r__8)));
	    }
	} else {
	    if (*n == 1) {
		i__2 = j * b_dim1 + 1;
		bi.r = B(1,j).r, bi.i = B(1,j).i;
		i__2 = j * x_dim1 + 1;
		q__1.r = D(1) * X(1,j).r, q__1.i = D(1) * X(1,j).i;
		dx.r = q__1.r, dx.i = q__1.i;
		q__1.r = bi.r - dx.r, q__1.i = bi.i - dx.i;
		WORK(1).r = q__1.r, WORK(1).i = q__1.i;
		RWORK(1) = (r__1 = bi.r, dabs(r__1)) + (r__2 = r_imag(&bi), 
			dabs(r__2)) + ((r__3 = dx.r, dabs(r__3)) + (r__4 = 
			r_imag(&dx), dabs(r__4)));
	    } else {
		i__2 = j * b_dim1 + 1;
		bi.r = B(1,j).r, bi.i = B(1,j).i;
		i__2 = j * x_dim1 + 1;
		q__1.r = D(1) * X(1,j).r, q__1.i = D(1) * X(1,j).i;
		dx.r = q__1.r, dx.i = q__1.i;
		r_cnjg(&q__2, &E(1));
		i__2 = j * x_dim1 + 2;
		q__1.r = q__2.r * X(2,j).r - q__2.i * X(2,j).i, q__1.i = 
			q__2.r * X(2,j).i + q__2.i * X(2,j).r;
		ex.r = q__1.r, ex.i = q__1.i;
		q__2.r = bi.r - dx.r, q__2.i = bi.i - dx.i;
		q__1.r = q__2.r - ex.r, q__1.i = q__2.i - ex.i;
		WORK(1).r = q__1.r, WORK(1).i = q__1.i;
		i__2 = j * x_dim1 + 2;
		RWORK(1) = (r__1 = bi.r, dabs(r__1)) + (r__2 = r_imag(&bi), 
			dabs(r__2)) + ((r__3 = dx.r, dabs(r__3)) + (r__4 = 
			r_imag(&dx), dabs(r__4))) + ((r__5 = E(1).r, dabs(
			r__5)) + (r__6 = r_imag(&E(1)), dabs(r__6))) * ((r__7 
			= X(2,j).r, dabs(r__7)) + (r__8 = r_imag(&X(2,j)), dabs(r__8)));
		i__2 = *n - 1;
		for (i = 2; i <= *n-1; ++i) {
		    i__3 = i + j * b_dim1;
		    bi.r = B(i,j).r, bi.i = B(i,j).i;
		    i__3 = i - 1;
		    i__4 = i - 1 + j * x_dim1;
		    q__1.r = E(i-1).r * X(i-1,j).r - E(i-1).i * X(i-1,j).i, 
			    q__1.i = E(i-1).r * X(i-1,j).i + E(i-1).i * X(i-1,j).r;
		    cx.r = q__1.r, cx.i = q__1.i;
		    i__3 = i;
		    i__4 = i + j * x_dim1;
		    q__1.r = D(i) * X(i,j).r, q__1.i = D(i) * X(i,j)
			    .i;
		    dx.r = q__1.r, dx.i = q__1.i;
		    r_cnjg(&q__2, &E(i));
		    i__3 = i + 1 + j * x_dim1;
		    q__1.r = q__2.r * X(i+1,j).r - q__2.i * X(i+1,j).i, q__1.i =
			     q__2.r * X(i+1,j).i + q__2.i * X(i+1,j).r;
		    ex.r = q__1.r, ex.i = q__1.i;
		    i__3 = i;
		    q__3.r = bi.r - cx.r, q__3.i = bi.i - cx.i;
		    q__2.r = q__3.r - dx.r, q__2.i = q__3.i - dx.i;
		    q__1.r = q__2.r - ex.r, q__1.i = q__2.i - ex.i;
		    WORK(i).r = q__1.r, WORK(i).i = q__1.i;
		    i__3 = i - 1;
		    i__4 = i - 1 + j * x_dim1;
		    i__5 = i;
		    i__6 = i + 1 + j * x_dim1;
		    RWORK(i) = (r__1 = bi.r, dabs(r__1)) + (r__2 = r_imag(&bi)
			    , dabs(r__2)) + ((r__3 = E(i-1).r, dabs(r__3)) + 
			    (r__4 = r_imag(&E(i - 1)), dabs(r__4))) * ((r__5 =
			     X(i-1,j).r, dabs(r__5)) + (r__6 = r_imag(&X(i-1,j)), dabs(r__6))) + ((r__7 = dx.r, 
			    dabs(r__7)) + (r__8 = r_imag(&dx), dabs(r__8))) + 
			    ((r__9 = E(i).r, dabs(r__9)) + (r__10 = r_imag(
			    &E(i)), dabs(r__10))) * ((r__11 = X(i+1,j).r, dabs(
			    r__11)) + (r__12 = r_imag(&X(i+1,j)),
			     dabs(r__12)));
/* L40: */
		}
		i__2 = *n + j * b_dim1;
		bi.r = B(*n,j).r, bi.i = B(*n,j).i;
		i__2 = *n - 1;
		i__3 = *n - 1 + j * x_dim1;
		q__1.r = E(*n-1).r * X(*n-1,j).r - E(*n-1).i * X(*n-1,j).i, 
			q__1.i = E(*n-1).r * X(*n-1,j).i + E(*n-1).i * X(*n-1,j)
			.r;
		cx.r = q__1.r, cx.i = q__1.i;
		i__2 = *n;
		i__3 = *n + j * x_dim1;
		q__1.r = D(*n) * X(*n,j).r, q__1.i = D(*n) * X(*n,j).i;
		dx.r = q__1.r, dx.i = q__1.i;
		i__2 = *n;
		q__2.r = bi.r - cx.r, q__2.i = bi.i - cx.i;
		q__1.r = q__2.r - dx.r, q__1.i = q__2.i - dx.i;
		WORK(*n).r = q__1.r, WORK(*n).i = q__1.i;
		i__2 = *n - 1;
		i__3 = *n - 1 + j * x_dim1;
		RWORK(*n) = (r__1 = bi.r, dabs(r__1)) + (r__2 = r_imag(&bi), 
			dabs(r__2)) + ((r__3 = E(*n-1).r, dabs(r__3)) + (r__4 
			= r_imag(&E(*n - 1)), dabs(r__4))) * ((r__5 = X(*n-1,j)
			.r, dabs(r__5)) + (r__6 = r_imag(&X(*n-1,j)), dabs(r__6))) + ((r__7 = dx.r, dabs(r__7)) + 
			(r__8 = r_imag(&dx), dabs(r__8)));
	    }
	}

/*        Compute componentwise relative backward error from formula 
  

          max(i) ( abs(R(i)) / ( abs(A)*abs(X) + abs(B) )(i) )   

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

	    cpttrs_(uplo, n, &c__1, &DF(1), &EF(1), &WORK(1), n, info);
	    caxpy_(n, &c_b16, &WORK(1), &c__1, &X(1,j), &c__1);
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
          abs(A)*abs(X) + abs(B) is less than SAFE2. */

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
	ix = isamax_(n, &RWORK(1), &c__1);
	FERR(j) = RWORK(ix);

/*        Estimate the norm of inv(A).   

          Solve M(A) * x = e, where M(A) = (m(i,j)) is given by   

             m(i,j) =  abs(A(i,j)), i = j,   
             m(i,j) = -abs(A(i,j)), i .ne. j,   

          and e = [ 1, 1, ..., 1 ]'.  Note M(A) = M(L)*D*M(L)'.   

          Solve M(L) * x = e. */

	RWORK(1) = 1.f;
	i__2 = *n;
	for (i = 2; i <= *n; ++i) {
	    RWORK(i) = RWORK(i - 1) * c_abs(&EF(i - 1)) + 1.f;
/* L70: */
	}

/*        Solve D * M(L)' * x = b. */

	RWORK(*n) /= DF(*n);
	for (i = *n - 1; i >= 1; --i) {
	    RWORK(i) = RWORK(i) / DF(i) + RWORK(i + 1) * c_abs(&EF(i));
/* L80: */
	}

/*        Compute norm(inv(A)) = max(x(i)), 1<=i<=n. */

	ix = isamax_(n, &RWORK(1), &c__1);
	FERR(j) *= (r__1 = RWORK(ix), dabs(r__1));

/*        Normalize error. */

	lstres = 0.f;
	i__2 = *n;
	for (i = 1; i <= *n; ++i) {
/* Computing MAX */
	    r__1 = lstres, r__2 = c_abs(&X(i,j));
	    lstres = dmax(r__1,r__2);
/* L90: */
	}
	if (lstres != 0.f) {
	    FERR(j) /= lstres;
	}

/* L100: */
    }

    return 0;

/*     End of CPTRFS */

} /* cptrfs_ */

