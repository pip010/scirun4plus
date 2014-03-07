#include "f2c.h"

/* Subroutine */ int zpttrs_(char *uplo, integer *n, integer *nrhs, 
	doublereal *d, doublecomplex *e, doublecomplex *b, integer *ldb, 
	integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       March 31, 1993   


    Purpose   
    =======   

    ZPTTRS solves a system of linear equations A * X = B with a   
    Hermitian positive definite tridiagonal matrix A using the   
    factorization A = U**H*D*U or A = L*D*L**H computed by ZPTTRF.   

    Arguments   
    =========   

    UPLO    (input) CHARACTER*1   
            Specifies whether the superdiagonal or the subdiagonal   
            of the tridiagonal matrix A is stored and the form of the   
            factorization:   
            = 'U':  E is the superdiagonal of U, and A = U'*D*U;   
            = 'L':  E is the subdiagonal of L, and A = L*D*L'.   
            (The two forms are equivalent if A is real.)   

    N       (input) INTEGER   
            The order of the tridiagonal matrix A.  N >= 0.   

    NRHS    (input) INTEGER   
            The number of right hand sides, i.e., the number of columns   
            of the matrix B.  NRHS >= 0.   

    D       (input) DOUBLE PRECISION array, dimension (N)   
            The n diagonal elements of the diagonal matrix D from the   
            factorization computed by ZPTTRF.   

    E       (input) COMPLEX*16 array, dimension (N-1)   
            The (n-1) off-diagonal elements of the unit bidiagonal   
            factor U or L from the factorization computed by ZPTTRF   
            (see UPLO).   

    B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)   
            On entry, the right hand side matrix B.   
            On exit, the solution matrix X.   

    LDB     (input) INTEGER   
            The leading dimension of the array B.  LDB >= max(1,N).   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   

    ===================================================================== 
  


       Test the input arguments.   

    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    integer b_dim1, b_offset, i__1, i__2, i__3, i__4, i__5, i__6;
    doublecomplex z__1, z__2, z__3, z__4;
    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);
    /* Local variables */
    static integer i, j;
    extern logical lsame_(char *, char *);
    static logical upper;
    extern /* Subroutine */ int xerbla_(char *, integer *);


#define D(I) d[(I)-1]
#define E(I) e[(I)-1]

#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]

    *info = 0;
    upper = lsame_(uplo, "U");
    if (! upper && ! lsame_(uplo, "L")) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*nrhs < 0) {
	*info = -3;
    } else if (*ldb < max(1,*n)) {
	*info = -7;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("ZPTTRS", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return 0;
    }

    if (upper) {

/*        Solve A * X = B using the factorization A = U'*D*U,   
          overwriting each right hand side vector with its solution. 
*/

	i__1 = *nrhs;
	for (j = 1; j <= *nrhs; ++j) {

/*           Solve U' * x = b. */

	    i__2 = *n;
	    for (i = 2; i <= *n; ++i) {
		i__3 = i + j * b_dim1;
		i__4 = i + j * b_dim1;
		i__5 = i - 1 + j * b_dim1;
		d_cnjg(&z__3, &E(i - 1));
		z__2.r = B(i-1,j).r * z__3.r - B(i-1,j).i * z__3.i, z__2.i = B(i-1,j).r * z__3.i + B(i-1,j).i * z__3.r;
		z__1.r = B(i,j).r - z__2.r, z__1.i = B(i,j).i - z__2.i;
		B(i,j).r = z__1.r, B(i,j).i = z__1.i;
/* L10: */
	    }

/*           Solve D * U * x = b. */

	    i__2 = *n + j * b_dim1;
	    i__3 = *n + j * b_dim1;
	    i__4 = *n;
	    z__1.r = B(*n,j).r / D(*n), z__1.i = B(*n,j).i / D(*n);
	    B(*n,j).r = z__1.r, B(*n,j).i = z__1.i;
	    for (i = *n - 1; i >= 1; --i) {
		i__2 = i + j * b_dim1;
		i__3 = i + j * b_dim1;
		i__4 = i;
		z__2.r = B(i,j).r / D(i), z__2.i = B(i,j).i / D(i);
		i__5 = i + 1 + j * b_dim1;
		i__6 = i;
		z__3.r = B(i+1,j).r * E(i).r - B(i+1,j).i * E(i).i, 
			z__3.i = B(i+1,j).r * E(i).i + B(i+1,j).i * E(i)
			.r;
		z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
		B(i,j).r = z__1.r, B(i,j).i = z__1.i;
/* L20: */
	    }
/* L30: */
	}
    } else {

/*        Solve A * X = B using the factorization A = L*D*L',   
          overwriting each right hand side vector with its solution. 
*/

	i__1 = *nrhs;
	for (j = 1; j <= *nrhs; ++j) {

/*           Solve L * x = b. */

	    i__2 = *n;
	    for (i = 2; i <= *n; ++i) {
		i__3 = i + j * b_dim1;
		i__4 = i + j * b_dim1;
		i__5 = i - 1 + j * b_dim1;
		i__6 = i - 1;
		z__2.r = B(i-1,j).r * E(i-1).r - B(i-1,j).i * E(i-1).i, 
			z__2.i = B(i-1,j).r * E(i-1).i + B(i-1,j).i * E(i-1)
			.r;
		z__1.r = B(i,j).r - z__2.r, z__1.i = B(i,j).i - z__2.i;
		B(i,j).r = z__1.r, B(i,j).i = z__1.i;
/* L40: */
	    }

/*           Solve D * L' * x = b. */

	    i__2 = *n + j * b_dim1;
	    i__3 = *n + j * b_dim1;
	    i__4 = *n;
	    z__1.r = B(*n,j).r / D(*n), z__1.i = B(*n,j).i / D(*n);
	    B(*n,j).r = z__1.r, B(*n,j).i = z__1.i;
	    for (i = *n - 1; i >= 1; --i) {
		i__2 = i + j * b_dim1;
		i__3 = i + j * b_dim1;
		i__4 = i;
		z__2.r = B(i,j).r / D(i), z__2.i = B(i,j).i / D(i);
		i__5 = i + 1 + j * b_dim1;
		d_cnjg(&z__4, &E(i));
		z__3.r = B(i+1,j).r * z__4.r - B(i+1,j).i * z__4.i, z__3.i = B(i+1,j).r * z__4.i + B(i+1,j).i * z__4.r;
		z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
		B(i,j).r = z__1.r, B(i,j).i = z__1.i;
/* L50: */
	    }
/* L60: */
	}
    }

    return 0;

/*     End of ZPTTRS */

} /* zpttrs_ */

