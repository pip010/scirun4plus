#include "f2c.h"

/* Subroutine */ int clagtm_(char *trans, integer *n, integer *nrhs, real *
	alpha, complex *dl, complex *d, complex *du, complex *x, integer *ldx,
	 real *beta, complex *b, integer *ldb)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    CLAGTM performs a matrix-vector product of the form   

       B := alpha * A * X + beta * B   

    where A is a tridiagonal matrix of order N, B and X are N by NRHS   
    matrices, and alpha and beta are real scalars, each of which may be   
    0., 1., or -1.   

    Arguments   
    =========   

    TRANS   (input) CHARACTER   
            Specifies the operation applied to A.   
            = 'N':  No transpose, B := alpha * A * X + beta * B   
            = 'T':  Transpose,    B := alpha * A**T * X + beta * B   
            = 'C':  Conjugate transpose, B := alpha * A**H * X + beta * B 
  

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    NRHS    (input) INTEGER   
            The number of right hand sides, i.e., the number of columns   
            of the matrices X and B.   

    ALPHA   (input) REAL   
            The scalar alpha.  ALPHA must be 0., 1., or -1.; otherwise,   
            it is assumed to be 0.   

    DL      (input) COMPLEX array, dimension (N-1)   
            The (n-1) sub-diagonal elements of T.   

    D       (input) COMPLEX array, dimension (N)   
            The diagonal elements of T.   

    DU      (input) COMPLEX array, dimension (N-1)   
            The (n-1) super-diagonal elements of T.   

    X       (input) COMPLEX array, dimension (LDX,NRHS)   
            The N by NRHS matrix X.   
    LDX     (input) INTEGER   
            The leading dimension of the array X.  LDX >= max(N,1).   

    BETA    (input) REAL   
            The scalar beta.  BETA must be 0., 1., or -1.; otherwise,   
            it is assumed to be 1.   

    B       (input/output) COMPLEX array, dimension (LDB,NRHS)   
            On entry, the N by NRHS matrix B.   
            On exit, B is overwritten by the matrix expression   
            B := alpha * A * X + beta * B.   

    LDB     (input) INTEGER   
            The leading dimension of the array B.  LDB >= max(N,1).   

    ===================================================================== 
  


    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    integer b_dim1, b_offset, x_dim1, x_offset, i__1, i__2, i__3, i__4, i__5, 
	    i__6, i__7, i__8, i__9, i__10;
    complex q__1, q__2, q__3, q__4, q__5, q__6, q__7, q__8, q__9;
    /* Builtin functions */
    void r_cnjg(complex *, complex *);
    /* Local variables */
    static integer i, j;
    extern logical lsame_(char *, char *);


#define DL(I) dl[(I)-1]
#define D(I) d[(I)-1]
#define DU(I) du[(I)-1]

#define X(I,J) x[(I)-1 + ((J)-1)* ( *ldx)]
#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]

    if (*n == 0) {
	return 0;
    }

/*     Multiply B by BETA if BETA.NE.1. */

    if (*beta == 0.f) {
	i__1 = *nrhs;
	for (j = 1; j <= *nrhs; ++j) {
	    i__2 = *n;
	    for (i = 1; i <= *n; ++i) {
		i__3 = i + j * b_dim1;
		B(i,j).r = 0.f, B(i,j).i = 0.f;
/* L10: */
	    }
/* L20: */
	}
    } else if (*beta == -1.f) {
	i__1 = *nrhs;
	for (j = 1; j <= *nrhs; ++j) {
	    i__2 = *n;
	    for (i = 1; i <= *n; ++i) {
		i__3 = i + j * b_dim1;
		i__4 = i + j * b_dim1;
		q__1.r = -(doublereal)B(i,j).r, q__1.i = -(doublereal)B(i,j)
			.i;
		B(i,j).r = q__1.r, B(i,j).i = q__1.i;
/* L30: */
	    }
/* L40: */
	}
    }

    if (*alpha == 1.f) {
	if (lsame_(trans, "N")) {

/*           Compute B := B + A*X */

	    i__1 = *nrhs;
	    for (j = 1; j <= *nrhs; ++j) {
		if (*n == 1) {
		    i__2 = j * b_dim1 + 1;
		    i__3 = j * b_dim1 + 1;
		    i__4 = j * x_dim1 + 1;
		    q__2.r = D(1).r * X(1,j).r - D(1).i * X(1,j).i, q__2.i =
			     D(1).r * X(1,j).i + D(1).i * X(1,j).r;
		    q__1.r = B(1,j).r + q__2.r, q__1.i = B(1,j).i + q__2.i;
		    B(1,j).r = q__1.r, B(1,j).i = q__1.i;
		} else {
		    i__2 = j * b_dim1 + 1;
		    i__3 = j * b_dim1 + 1;
		    i__4 = j * x_dim1 + 1;
		    q__3.r = D(1).r * X(1,j).r - D(1).i * X(1,j).i, q__3.i =
			     D(1).r * X(1,j).i + D(1).i * X(1,j).r;
		    q__2.r = B(1,j).r + q__3.r, q__2.i = B(1,j).i + q__3.i;
		    i__5 = j * x_dim1 + 2;
		    q__4.r = DU(1).r * X(2,j).r - DU(1).i * X(2,j).i, 
			    q__4.i = DU(1).r * X(2,j).i + DU(1).i * X(2,j)
			    .r;
		    q__1.r = q__2.r + q__4.r, q__1.i = q__2.i + q__4.i;
		    B(1,j).r = q__1.r, B(1,j).i = q__1.i;
		    i__2 = *n + j * b_dim1;
		    i__3 = *n + j * b_dim1;
		    i__4 = *n - 1;
		    i__5 = *n - 1 + j * x_dim1;
		    q__3.r = DL(*n-1).r * X(*n-1,j).r - DL(*n-1).i * X(*n-1,j).i, 
			    q__3.i = DL(*n-1).r * X(*n-1,j).i + DL(*n-1).i * X(*n-1,j).r;
		    q__2.r = B(*n,j).r + q__3.r, q__2.i = B(*n,j).i + q__3.i;
		    i__6 = *n;
		    i__7 = *n + j * x_dim1;
		    q__4.r = D(*n).r * X(*n,j).r - D(*n).i * X(*n,j).i, 
			    q__4.i = D(*n).r * X(*n,j).i + D(*n).i * X(*n,j).r;
		    q__1.r = q__2.r + q__4.r, q__1.i = q__2.i + q__4.i;
		    B(*n,j).r = q__1.r, B(*n,j).i = q__1.i;
		    i__2 = *n - 1;
		    for (i = 2; i <= *n-1; ++i) {
			i__3 = i + j * b_dim1;
			i__4 = i + j * b_dim1;
			i__5 = i - 1;
			i__6 = i - 1 + j * x_dim1;
			q__4.r = DL(i-1).r * X(i-1,j).r - DL(i-1).i * X(i-1,j)
				.i, q__4.i = DL(i-1).r * X(i-1,j).i + DL(i-1)
				.i * X(i-1,j).r;
			q__3.r = B(i,j).r + q__4.r, q__3.i = B(i,j).i + 
				q__4.i;
			i__7 = i;
			i__8 = i + j * x_dim1;
			q__5.r = D(i).r * X(i,j).r - D(i).i * X(i,j)
				.i, q__5.i = D(i).r * X(i,j).i + D(i)
				.i * X(i,j).r;
			q__2.r = q__3.r + q__5.r, q__2.i = q__3.i + q__5.i;
			i__9 = i;
			i__10 = i + 1 + j * x_dim1;
			q__6.r = DU(i).r * X(i+1,j).r - DU(i).i * X(i+1,j).i, q__6.i = DU(i).r * X(i+1,j).i + 
				DU(i).i * X(i+1,j).r;
			q__1.r = q__2.r + q__6.r, q__1.i = q__2.i + q__6.i;
			B(i,j).r = q__1.r, B(i,j).i = q__1.i;
/* L50: */
		    }
		}
/* L60: */
	    }
	} else if (lsame_(trans, "T")) {

/*           Compute B := B + A**T * X */

	    i__1 = *nrhs;
	    for (j = 1; j <= *nrhs; ++j) {
		if (*n == 1) {
		    i__2 = j * b_dim1 + 1;
		    i__3 = j * b_dim1 + 1;
		    i__4 = j * x_dim1 + 1;
		    q__2.r = D(1).r * X(1,j).r - D(1).i * X(1,j).i, q__2.i =
			     D(1).r * X(1,j).i + D(1).i * X(1,j).r;
		    q__1.r = B(1,j).r + q__2.r, q__1.i = B(1,j).i + q__2.i;
		    B(1,j).r = q__1.r, B(1,j).i = q__1.i;
		} else {
		    i__2 = j * b_dim1 + 1;
		    i__3 = j * b_dim1 + 1;
		    i__4 = j * x_dim1 + 1;
		    q__3.r = D(1).r * X(1,j).r - D(1).i * X(1,j).i, q__3.i =
			     D(1).r * X(1,j).i + D(1).i * X(1,j).r;
		    q__2.r = B(1,j).r + q__3.r, q__2.i = B(1,j).i + q__3.i;
		    i__5 = j * x_dim1 + 2;
		    q__4.r = DL(1).r * X(2,j).r - DL(1).i * X(2,j).i, 
			    q__4.i = DL(1).r * X(2,j).i + DL(1).i * X(2,j)
			    .r;
		    q__1.r = q__2.r + q__4.r, q__1.i = q__2.i + q__4.i;
		    B(1,j).r = q__1.r, B(1,j).i = q__1.i;
		    i__2 = *n + j * b_dim1;
		    i__3 = *n + j * b_dim1;
		    i__4 = *n - 1;
		    i__5 = *n - 1 + j * x_dim1;
		    q__3.r = DU(*n-1).r * X(*n-1,j).r - DU(*n-1).i * X(*n-1,j).i, 
			    q__3.i = DU(*n-1).r * X(*n-1,j).i + DU(*n-1).i * X(*n-1,j).r;
		    q__2.r = B(*n,j).r + q__3.r, q__2.i = B(*n,j).i + q__3.i;
		    i__6 = *n;
		    i__7 = *n + j * x_dim1;
		    q__4.r = D(*n).r * X(*n,j).r - D(*n).i * X(*n,j).i, 
			    q__4.i = D(*n).r * X(*n,j).i + D(*n).i * X(*n,j).r;
		    q__1.r = q__2.r + q__4.r, q__1.i = q__2.i + q__4.i;
		    B(*n,j).r = q__1.r, B(*n,j).i = q__1.i;
		    i__2 = *n - 1;
		    for (i = 2; i <= *n-1; ++i) {
			i__3 = i + j * b_dim1;
			i__4 = i + j * b_dim1;
			i__5 = i - 1;
			i__6 = i - 1 + j * x_dim1;
			q__4.r = DU(i-1).r * X(i-1,j).r - DU(i-1).i * X(i-1,j)
				.i, q__4.i = DU(i-1).r * X(i-1,j).i + DU(i-1)
				.i * X(i-1,j).r;
			q__3.r = B(i,j).r + q__4.r, q__3.i = B(i,j).i + 
				q__4.i;
			i__7 = i;
			i__8 = i + j * x_dim1;
			q__5.r = D(i).r * X(i,j).r - D(i).i * X(i,j)
				.i, q__5.i = D(i).r * X(i,j).i + D(i)
				.i * X(i,j).r;
			q__2.r = q__3.r + q__5.r, q__2.i = q__3.i + q__5.i;
			i__9 = i;
			i__10 = i + 1 + j * x_dim1;
			q__6.r = DL(i).r * X(i+1,j).r - DL(i).i * X(i+1,j).i, q__6.i = DL(i).r * X(i+1,j).i + 
				DL(i).i * X(i+1,j).r;
			q__1.r = q__2.r + q__6.r, q__1.i = q__2.i + q__6.i;
			B(i,j).r = q__1.r, B(i,j).i = q__1.i;
/* L70: */
		    }
		}
/* L80: */
	    }
	} else if (lsame_(trans, "C")) {

/*           Compute B := B + A**H * X */

	    i__1 = *nrhs;
	    for (j = 1; j <= *nrhs; ++j) {
		if (*n == 1) {
		    i__2 = j * b_dim1 + 1;
		    i__3 = j * b_dim1 + 1;
		    r_cnjg(&q__3, &D(1));
		    i__4 = j * x_dim1 + 1;
		    q__2.r = q__3.r * X(1,j).r - q__3.i * X(1,j).i, q__2.i =
			     q__3.r * X(1,j).i + q__3.i * X(1,j).r;
		    q__1.r = B(1,j).r + q__2.r, q__1.i = B(1,j).i + q__2.i;
		    B(1,j).r = q__1.r, B(1,j).i = q__1.i;
		} else {
		    i__2 = j * b_dim1 + 1;
		    i__3 = j * b_dim1 + 1;
		    r_cnjg(&q__4, &D(1));
		    i__4 = j * x_dim1 + 1;
		    q__3.r = q__4.r * X(1,j).r - q__4.i * X(1,j).i, q__3.i =
			     q__4.r * X(1,j).i + q__4.i * X(1,j).r;
		    q__2.r = B(1,j).r + q__3.r, q__2.i = B(1,j).i + q__3.i;
		    r_cnjg(&q__6, &DL(1));
		    i__5 = j * x_dim1 + 2;
		    q__5.r = q__6.r * X(2,j).r - q__6.i * X(2,j).i, q__5.i =
			     q__6.r * X(2,j).i + q__6.i * X(2,j).r;
		    q__1.r = q__2.r + q__5.r, q__1.i = q__2.i + q__5.i;
		    B(1,j).r = q__1.r, B(1,j).i = q__1.i;
		    i__2 = *n + j * b_dim1;
		    i__3 = *n + j * b_dim1;
		    r_cnjg(&q__4, &DU(*n - 1));
		    i__4 = *n - 1 + j * x_dim1;
		    q__3.r = q__4.r * X(*n-1,j).r - q__4.i * X(*n-1,j).i, q__3.i =
			     q__4.r * X(*n-1,j).i + q__4.i * X(*n-1,j).r;
		    q__2.r = B(*n,j).r + q__3.r, q__2.i = B(*n,j).i + q__3.i;
		    r_cnjg(&q__6, &D(*n));
		    i__5 = *n + j * x_dim1;
		    q__5.r = q__6.r * X(*n,j).r - q__6.i * X(*n,j).i, q__5.i =
			     q__6.r * X(*n,j).i + q__6.i * X(*n,j).r;
		    q__1.r = q__2.r + q__5.r, q__1.i = q__2.i + q__5.i;
		    B(*n,j).r = q__1.r, B(*n,j).i = q__1.i;
		    i__2 = *n - 1;
		    for (i = 2; i <= *n-1; ++i) {
			i__3 = i + j * b_dim1;
			i__4 = i + j * b_dim1;
			r_cnjg(&q__5, &DU(i - 1));
			i__5 = i - 1 + j * x_dim1;
			q__4.r = q__5.r * X(i-1,j).r - q__5.i * X(i-1,j).i, 
				q__4.i = q__5.r * X(i-1,j).i + q__5.i * X(i-1,j)
				.r;
			q__3.r = B(i,j).r + q__4.r, q__3.i = B(i,j).i + 
				q__4.i;
			r_cnjg(&q__7, &D(i));
			i__6 = i + j * x_dim1;
			q__6.r = q__7.r * X(i,j).r - q__7.i * X(i,j).i, 
				q__6.i = q__7.r * X(i,j).i + q__7.i * X(i,j)
				.r;
			q__2.r = q__3.r + q__6.r, q__2.i = q__3.i + q__6.i;
			r_cnjg(&q__9, &DL(i));
			i__7 = i + 1 + j * x_dim1;
			q__8.r = q__9.r * X(i+1,j).r - q__9.i * X(i+1,j).i, 
				q__8.i = q__9.r * X(i+1,j).i + q__9.i * X(i+1,j)
				.r;
			q__1.r = q__2.r + q__8.r, q__1.i = q__2.i + q__8.i;
			B(i,j).r = q__1.r, B(i,j).i = q__1.i;
/* L90: */
		    }
		}
/* L100: */
	    }
	}
    } else if (*alpha == -1.f) {
	if (lsame_(trans, "N")) {

/*           Compute B := B - A*X */

	    i__1 = *nrhs;
	    for (j = 1; j <= *nrhs; ++j) {
		if (*n == 1) {
		    i__2 = j * b_dim1 + 1;
		    i__3 = j * b_dim1 + 1;
		    i__4 = j * x_dim1 + 1;
		    q__2.r = D(1).r * X(1,j).r - D(1).i * X(1,j).i, q__2.i =
			     D(1).r * X(1,j).i + D(1).i * X(1,j).r;
		    q__1.r = B(1,j).r - q__2.r, q__1.i = B(1,j).i - q__2.i;
		    B(1,j).r = q__1.r, B(1,j).i = q__1.i;
		} else {
		    i__2 = j * b_dim1 + 1;
		    i__3 = j * b_dim1 + 1;
		    i__4 = j * x_dim1 + 1;
		    q__3.r = D(1).r * X(1,j).r - D(1).i * X(1,j).i, q__3.i =
			     D(1).r * X(1,j).i + D(1).i * X(1,j).r;
		    q__2.r = B(1,j).r - q__3.r, q__2.i = B(1,j).i - q__3.i;
		    i__5 = j * x_dim1 + 2;
		    q__4.r = DU(1).r * X(2,j).r - DU(1).i * X(2,j).i, 
			    q__4.i = DU(1).r * X(2,j).i + DU(1).i * X(2,j)
			    .r;
		    q__1.r = q__2.r - q__4.r, q__1.i = q__2.i - q__4.i;
		    B(1,j).r = q__1.r, B(1,j).i = q__1.i;
		    i__2 = *n + j * b_dim1;
		    i__3 = *n + j * b_dim1;
		    i__4 = *n - 1;
		    i__5 = *n - 1 + j * x_dim1;
		    q__3.r = DL(*n-1).r * X(*n-1,j).r - DL(*n-1).i * X(*n-1,j).i, 
			    q__3.i = DL(*n-1).r * X(*n-1,j).i + DL(*n-1).i * X(*n-1,j).r;
		    q__2.r = B(*n,j).r - q__3.r, q__2.i = B(*n,j).i - q__3.i;
		    i__6 = *n;
		    i__7 = *n + j * x_dim1;
		    q__4.r = D(*n).r * X(*n,j).r - D(*n).i * X(*n,j).i, 
			    q__4.i = D(*n).r * X(*n,j).i + D(*n).i * X(*n,j).r;
		    q__1.r = q__2.r - q__4.r, q__1.i = q__2.i - q__4.i;
		    B(*n,j).r = q__1.r, B(*n,j).i = q__1.i;
		    i__2 = *n - 1;
		    for (i = 2; i <= *n-1; ++i) {
			i__3 = i + j * b_dim1;
			i__4 = i + j * b_dim1;
			i__5 = i - 1;
			i__6 = i - 1 + j * x_dim1;
			q__4.r = DL(i-1).r * X(i-1,j).r - DL(i-1).i * X(i-1,j)
				.i, q__4.i = DL(i-1).r * X(i-1,j).i + DL(i-1)
				.i * X(i-1,j).r;
			q__3.r = B(i,j).r - q__4.r, q__3.i = B(i,j).i - 
				q__4.i;
			i__7 = i;
			i__8 = i + j * x_dim1;
			q__5.r = D(i).r * X(i,j).r - D(i).i * X(i,j)
				.i, q__5.i = D(i).r * X(i,j).i + D(i)
				.i * X(i,j).r;
			q__2.r = q__3.r - q__5.r, q__2.i = q__3.i - q__5.i;
			i__9 = i;
			i__10 = i + 1 + j * x_dim1;
			q__6.r = DU(i).r * X(i+1,j).r - DU(i).i * X(i+1,j).i, q__6.i = DU(i).r * X(i+1,j).i + 
				DU(i).i * X(i+1,j).r;
			q__1.r = q__2.r - q__6.r, q__1.i = q__2.i - q__6.i;
			B(i,j).r = q__1.r, B(i,j).i = q__1.i;
/* L110: */
		    }
		}
/* L120: */
	    }
	} else if (lsame_(trans, "T")) {

/*           Compute B := B - A'*X */

	    i__1 = *nrhs;
	    for (j = 1; j <= *nrhs; ++j) {
		if (*n == 1) {
		    i__2 = j * b_dim1 + 1;
		    i__3 = j * b_dim1 + 1;
		    i__4 = j * x_dim1 + 1;
		    q__2.r = D(1).r * X(1,j).r - D(1).i * X(1,j).i, q__2.i =
			     D(1).r * X(1,j).i + D(1).i * X(1,j).r;
		    q__1.r = B(1,j).r - q__2.r, q__1.i = B(1,j).i - q__2.i;
		    B(1,j).r = q__1.r, B(1,j).i = q__1.i;
		} else {
		    i__2 = j * b_dim1 + 1;
		    i__3 = j * b_dim1 + 1;
		    i__4 = j * x_dim1 + 1;
		    q__3.r = D(1).r * X(1,j).r - D(1).i * X(1,j).i, q__3.i =
			     D(1).r * X(1,j).i + D(1).i * X(1,j).r;
		    q__2.r = B(1,j).r - q__3.r, q__2.i = B(1,j).i - q__3.i;
		    i__5 = j * x_dim1 + 2;
		    q__4.r = DL(1).r * X(2,j).r - DL(1).i * X(2,j).i, 
			    q__4.i = DL(1).r * X(2,j).i + DL(1).i * X(2,j)
			    .r;
		    q__1.r = q__2.r - q__4.r, q__1.i = q__2.i - q__4.i;
		    B(1,j).r = q__1.r, B(1,j).i = q__1.i;
		    i__2 = *n + j * b_dim1;
		    i__3 = *n + j * b_dim1;
		    i__4 = *n - 1;
		    i__5 = *n - 1 + j * x_dim1;
		    q__3.r = DU(*n-1).r * X(*n-1,j).r - DU(*n-1).i * X(*n-1,j).i, 
			    q__3.i = DU(*n-1).r * X(*n-1,j).i + DU(*n-1).i * X(*n-1,j).r;
		    q__2.r = B(*n,j).r - q__3.r, q__2.i = B(*n,j).i - q__3.i;
		    i__6 = *n;
		    i__7 = *n + j * x_dim1;
		    q__4.r = D(*n).r * X(*n,j).r - D(*n).i * X(*n,j).i, 
			    q__4.i = D(*n).r * X(*n,j).i + D(*n).i * X(*n,j).r;
		    q__1.r = q__2.r - q__4.r, q__1.i = q__2.i - q__4.i;
		    B(*n,j).r = q__1.r, B(*n,j).i = q__1.i;
		    i__2 = *n - 1;
		    for (i = 2; i <= *n-1; ++i) {
			i__3 = i + j * b_dim1;
			i__4 = i + j * b_dim1;
			i__5 = i - 1;
			i__6 = i - 1 + j * x_dim1;
			q__4.r = DU(i-1).r * X(i-1,j).r - DU(i-1).i * X(i-1,j)
				.i, q__4.i = DU(i-1).r * X(i-1,j).i + DU(i-1)
				.i * X(i-1,j).r;
			q__3.r = B(i,j).r - q__4.r, q__3.i = B(i,j).i - 
				q__4.i;
			i__7 = i;
			i__8 = i + j * x_dim1;
			q__5.r = D(i).r * X(i,j).r - D(i).i * X(i,j)
				.i, q__5.i = D(i).r * X(i,j).i + D(i)
				.i * X(i,j).r;
			q__2.r = q__3.r - q__5.r, q__2.i = q__3.i - q__5.i;
			i__9 = i;
			i__10 = i + 1 + j * x_dim1;
			q__6.r = DL(i).r * X(i+1,j).r - DL(i).i * X(i+1,j).i, q__6.i = DL(i).r * X(i+1,j).i + 
				DL(i).i * X(i+1,j).r;
			q__1.r = q__2.r - q__6.r, q__1.i = q__2.i - q__6.i;
			B(i,j).r = q__1.r, B(i,j).i = q__1.i;
/* L130: */
		    }
		}
/* L140: */
	    }
	} else if (lsame_(trans, "C")) {

/*           Compute B := B - A'*X */

	    i__1 = *nrhs;
	    for (j = 1; j <= *nrhs; ++j) {
		if (*n == 1) {
		    i__2 = j * b_dim1 + 1;
		    i__3 = j * b_dim1 + 1;
		    r_cnjg(&q__3, &D(1));
		    i__4 = j * x_dim1 + 1;
		    q__2.r = q__3.r * X(1,j).r - q__3.i * X(1,j).i, q__2.i =
			     q__3.r * X(1,j).i + q__3.i * X(1,j).r;
		    q__1.r = B(1,j).r - q__2.r, q__1.i = B(1,j).i - q__2.i;
		    B(1,j).r = q__1.r, B(1,j).i = q__1.i;
		} else {
		    i__2 = j * b_dim1 + 1;
		    i__3 = j * b_dim1 + 1;
		    r_cnjg(&q__4, &D(1));
		    i__4 = j * x_dim1 + 1;
		    q__3.r = q__4.r * X(1,j).r - q__4.i * X(1,j).i, q__3.i =
			     q__4.r * X(1,j).i + q__4.i * X(1,j).r;
		    q__2.r = B(1,j).r - q__3.r, q__2.i = B(1,j).i - q__3.i;
		    r_cnjg(&q__6, &DL(1));
		    i__5 = j * x_dim1 + 2;
		    q__5.r = q__6.r * X(2,j).r - q__6.i * X(2,j).i, q__5.i =
			     q__6.r * X(2,j).i + q__6.i * X(2,j).r;
		    q__1.r = q__2.r - q__5.r, q__1.i = q__2.i - q__5.i;
		    B(1,j).r = q__1.r, B(1,j).i = q__1.i;
		    i__2 = *n + j * b_dim1;
		    i__3 = *n + j * b_dim1;
		    r_cnjg(&q__4, &DU(*n - 1));
		    i__4 = *n - 1 + j * x_dim1;
		    q__3.r = q__4.r * X(*n-1,j).r - q__4.i * X(*n-1,j).i, q__3.i =
			     q__4.r * X(*n-1,j).i + q__4.i * X(*n-1,j).r;
		    q__2.r = B(*n,j).r - q__3.r, q__2.i = B(*n,j).i - q__3.i;
		    r_cnjg(&q__6, &D(*n));
		    i__5 = *n + j * x_dim1;
		    q__5.r = q__6.r * X(*n,j).r - q__6.i * X(*n,j).i, q__5.i =
			     q__6.r * X(*n,j).i + q__6.i * X(*n,j).r;
		    q__1.r = q__2.r - q__5.r, q__1.i = q__2.i - q__5.i;
		    B(*n,j).r = q__1.r, B(*n,j).i = q__1.i;
		    i__2 = *n - 1;
		    for (i = 2; i <= *n-1; ++i) {
			i__3 = i + j * b_dim1;
			i__4 = i + j * b_dim1;
			r_cnjg(&q__5, &DU(i - 1));
			i__5 = i - 1 + j * x_dim1;
			q__4.r = q__5.r * X(i-1,j).r - q__5.i * X(i-1,j).i, 
				q__4.i = q__5.r * X(i-1,j).i + q__5.i * X(i-1,j)
				.r;
			q__3.r = B(i,j).r - q__4.r, q__3.i = B(i,j).i - 
				q__4.i;
			r_cnjg(&q__7, &D(i));
			i__6 = i + j * x_dim1;
			q__6.r = q__7.r * X(i,j).r - q__7.i * X(i,j).i, 
				q__6.i = q__7.r * X(i,j).i + q__7.i * X(i,j)
				.r;
			q__2.r = q__3.r - q__6.r, q__2.i = q__3.i - q__6.i;
			r_cnjg(&q__9, &DL(i));
			i__7 = i + 1 + j * x_dim1;
			q__8.r = q__9.r * X(i+1,j).r - q__9.i * X(i+1,j).i, 
				q__8.i = q__9.r * X(i+1,j).i + q__9.i * X(i+1,j)
				.r;
			q__1.r = q__2.r - q__8.r, q__1.i = q__2.i - q__8.i;
			B(i,j).r = q__1.r, B(i,j).i = q__1.i;
/* L150: */
		    }
		}
/* L160: */
	    }
	}
    }
    return 0;

/*     End of CLAGTM */

} /* clagtm_ */

