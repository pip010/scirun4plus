#include "f2c.h"

/* Subroutine */ int zgttrs_(char *trans, integer *n, integer *nrhs, 
	doublecomplex *dl, doublecomplex *d, doublecomplex *du, doublecomplex 
	*du2, integer *ipiv, doublecomplex *b, integer *ldb, integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    ZGTTRS solves one of the systems of equations   
       A * X = B,  A**T * X = B,  or  A**H * X = B,   
    with a tridiagonal matrix A using the LU factorization computed   
    by ZGTTRF.   

    Arguments   
    =========   

    TRANS   (input) CHARACTER   
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
            The (n-1) multipliers that define the matrix L from the   
            LU factorization of A.   

    D       (input) COMPLEX*16 array, dimension (N)   
            The n diagonal elements of the upper triangular matrix U from 
  
            the LU factorization of A.   

    DU      (input) COMPLEX*16 array, dimension (N-1)   
            The (n-1) elements of the first superdiagonal of U.   

    DU2     (input) COMPLEX*16 array, dimension (N-2)   
            The (n-2) elements of the second superdiagonal of U.   

    IPIV    (input) INTEGER array, dimension (N)   
            The pivot indices; for 1 <= i <= n, row i of the matrix was   
            interchanged with row IPIV(i).  IPIV(i) will always be either 
  
            i or i+1; IPIV(i) = i indicates a row interchange was not   
            required.   

    B       (input/output) COMPLEX*16 array, dimension (LDB,NRHS)   
            On entry, the right hand side matrix B.   
            On exit, B is overwritten by the solution matrix X.   

    LDB     (input) INTEGER   
            The leading dimension of the array B.  LDB >= max(1,N).   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   

    ===================================================================== 
  


    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    integer b_dim1, b_offset, i__1, i__2, i__3, i__4, i__5, i__6, i__7, i__8;
    doublecomplex z__1, z__2, z__3, z__4, z__5, z__6, z__7, z__8;
    /* Builtin functions */
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *), d_cnjg(
	    doublecomplex *, doublecomplex *);
    /* Local variables */
    static doublecomplex temp;
    static integer i, j;
    extern logical lsame_(char *, char *);
    extern /* Subroutine */ int xerbla_(char *, integer *);
    static logical notran;


#define DL(I) dl[(I)-1]
#define D(I) d[(I)-1]
#define DU(I) du[(I)-1]
#define DU2(I) du2[(I)-1]
#define IPIV(I) ipiv[(I)-1]

#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]

    *info = 0;
    notran = lsame_(trans, "N");
    if (! notran && ! lsame_(trans, "T") && ! lsame_(trans, "C")) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*nrhs < 0) {
	*info = -3;
    } else if (*ldb < max(*n,1)) {
	*info = -10;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("ZGTTRS", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0 || *nrhs == 0) {
	return 0;
    }

    if (notran) {

/*        Solve A*X = B using the LU factorization of A,   
          overwriting each right hand side vector with its solution. 
*/

	i__1 = *nrhs;
	for (j = 1; j <= *nrhs; ++j) {

/*           Solve L*x = b. */

	    i__2 = *n - 1;
	    for (i = 1; i <= *n-1; ++i) {
		if (IPIV(i) == i) {
		    i__3 = i + 1 + j * b_dim1;
		    i__4 = i + 1 + j * b_dim1;
		    i__5 = i;
		    i__6 = i + j * b_dim1;
		    z__2.r = DL(i).r * B(i,j).r - DL(i).i * B(i,j).i, 
			    z__2.i = DL(i).r * B(i,j).i + DL(i).i * B(i,j).r;
		    z__1.r = B(i+1,j).r - z__2.r, z__1.i = B(i+1,j).i - z__2.i;
		    B(i+1,j).r = z__1.r, B(i+1,j).i = z__1.i;
		} else {
		    i__3 = i + j * b_dim1;
		    temp.r = B(i,j).r, temp.i = B(i,j).i;
		    i__3 = i + j * b_dim1;
		    i__4 = i + 1 + j * b_dim1;
		    B(i,j).r = B(i+1,j).r, B(i,j).i = B(i+1,j).i;
		    i__3 = i + 1 + j * b_dim1;
		    i__4 = i;
		    i__5 = i + j * b_dim1;
		    z__2.r = DL(i).r * B(i,j).r - DL(i).i * B(i,j).i, 
			    z__2.i = DL(i).r * B(i,j).i + DL(i).i * B(i,j).r;
		    z__1.r = temp.r - z__2.r, z__1.i = temp.i - z__2.i;
		    B(i+1,j).r = z__1.r, B(i+1,j).i = z__1.i;
		}
/* L10: */
	    }

/*           Solve U*x = b. */

	    i__2 = *n + j * b_dim1;
	    z_div(&z__1, &B(*n,j), &D(*n));
	    B(*n,j).r = z__1.r, B(*n,j).i = z__1.i;
	    if (*n > 1) {
		i__2 = *n - 1 + j * b_dim1;
		i__3 = *n - 1 + j * b_dim1;
		i__4 = *n - 1;
		i__5 = *n + j * b_dim1;
		z__3.r = DU(*n-1).r * B(*n,j).r - DU(*n-1).i * B(*n,j).i, 
			z__3.i = DU(*n-1).r * B(*n,j).i + DU(*n-1).i * B(*n,j)
			.r;
		z__2.r = B(*n-1,j).r - z__3.r, z__2.i = B(*n-1,j).i - z__3.i;
		z_div(&z__1, &z__2, &D(*n - 1));
		B(*n-1,j).r = z__1.r, B(*n-1,j).i = z__1.i;
	    }
	    for (i = *n - 2; i >= 1; --i) {
		i__2 = i + j * b_dim1;
		i__3 = i + j * b_dim1;
		i__4 = i;
		i__5 = i + 1 + j * b_dim1;
		z__4.r = DU(i).r * B(i+1,j).r - DU(i).i * B(i+1,j).i, 
			z__4.i = DU(i).r * B(i+1,j).i + DU(i).i * B(i+1,j)
			.r;
		z__3.r = B(i,j).r - z__4.r, z__3.i = B(i,j).i - z__4.i;
		i__6 = i;
		i__7 = i + 2 + j * b_dim1;
		z__5.r = DU2(i).r * B(i+2,j).r - DU2(i).i * B(i+2,j).i, 
			z__5.i = DU2(i).r * B(i+2,j).i + DU2(i).i * B(i+2,j).r;
		z__2.r = z__3.r - z__5.r, z__2.i = z__3.i - z__5.i;
		z_div(&z__1, &z__2, &D(i));
		B(i,j).r = z__1.r, B(i,j).i = z__1.i;
/* L20: */
	    }
/* L30: */
	}
    } else if (lsame_(trans, "T")) {

/*        Solve A**T * X = B. */

	i__1 = *nrhs;
	for (j = 1; j <= *nrhs; ++j) {

/*           Solve U**T * x = b. */

	    i__2 = j * b_dim1 + 1;
	    z_div(&z__1, &B(1,j), &D(1));
	    B(1,j).r = z__1.r, B(1,j).i = z__1.i;
	    if (*n > 1) {
		i__2 = j * b_dim1 + 2;
		i__3 = j * b_dim1 + 2;
		i__4 = j * b_dim1 + 1;
		z__3.r = DU(1).r * B(1,j).r - DU(1).i * B(1,j).i, z__3.i = 
			DU(1).r * B(1,j).i + DU(1).i * B(1,j).r;
		z__2.r = B(2,j).r - z__3.r, z__2.i = B(2,j).i - z__3.i;
		z_div(&z__1, &z__2, &D(2));
		B(2,j).r = z__1.r, B(2,j).i = z__1.i;
	    }
	    i__2 = *n;
	    for (i = 3; i <= *n; ++i) {
		i__3 = i + j * b_dim1;
		i__4 = i + j * b_dim1;
		i__5 = i - 1;
		i__6 = i - 1 + j * b_dim1;
		z__4.r = DU(i-1).r * B(i-1,j).r - DU(i-1).i * B(i-1,j).i, 
			z__4.i = DU(i-1).r * B(i-1,j).i + DU(i-1).i * B(i-1,j)
			.r;
		z__3.r = B(i,j).r - z__4.r, z__3.i = B(i,j).i - z__4.i;
		i__7 = i - 2;
		i__8 = i - 2 + j * b_dim1;
		z__5.r = DU2(i-2).r * B(i-2,j).r - DU2(i-2).i * B(i-2,j).i, 
			z__5.i = DU2(i-2).r * B(i-2,j).i + DU2(i-2).i * B(i-2,j).r;
		z__2.r = z__3.r - z__5.r, z__2.i = z__3.i - z__5.i;
		z_div(&z__1, &z__2, &D(i));
		B(i,j).r = z__1.r, B(i,j).i = z__1.i;
/* L40: */
	    }

/*           Solve L**T * x = b. */

	    for (i = *n - 1; i >= 1; --i) {
		if (IPIV(i) == i) {
		    i__2 = i + j * b_dim1;
		    i__3 = i + j * b_dim1;
		    i__4 = i;
		    i__5 = i + 1 + j * b_dim1;
		    z__2.r = DL(i).r * B(i+1,j).r - DL(i).i * B(i+1,j).i, 
			    z__2.i = DL(i).r * B(i+1,j).i + DL(i).i * B(i+1,j).r;
		    z__1.r = B(i,j).r - z__2.r, z__1.i = B(i,j).i - z__2.i;
		    B(i,j).r = z__1.r, B(i,j).i = z__1.i;
		} else {
		    i__2 = i + 1 + j * b_dim1;
		    temp.r = B(i+1,j).r, temp.i = B(i+1,j).i;
		    i__2 = i + 1 + j * b_dim1;
		    i__3 = i + j * b_dim1;
		    i__4 = i;
		    z__2.r = DL(i).r * temp.r - DL(i).i * temp.i, 
			    z__2.i = DL(i).r * temp.i + DL(i).i * 
			    temp.r;
		    z__1.r = B(i,j).r - z__2.r, z__1.i = B(i,j).i - z__2.i;
		    B(i+1,j).r = z__1.r, B(i+1,j).i = z__1.i;
		    i__2 = i + j * b_dim1;
		    B(i,j).r = temp.r, B(i,j).i = temp.i;
		}
/* L50: */
	    }
/* L60: */
	}
    } else {

/*        Solve A**H * X = B. */

	i__1 = *nrhs;
	for (j = 1; j <= *nrhs; ++j) {

/*           Solve U**H * x = b. */

	    i__2 = j * b_dim1 + 1;
	    d_cnjg(&z__2, &D(1));
	    z_div(&z__1, &B(1,j), &z__2);
	    B(1,j).r = z__1.r, B(1,j).i = z__1.i;
	    if (*n > 1) {
		i__2 = j * b_dim1 + 2;
		i__3 = j * b_dim1 + 2;
		d_cnjg(&z__4, &DU(1));
		i__4 = j * b_dim1 + 1;
		z__3.r = z__4.r * B(1,j).r - z__4.i * B(1,j).i, z__3.i = 
			z__4.r * B(1,j).i + z__4.i * B(1,j).r;
		z__2.r = B(2,j).r - z__3.r, z__2.i = B(2,j).i - z__3.i;
		d_cnjg(&z__5, &D(2));
		z_div(&z__1, &z__2, &z__5);
		B(2,j).r = z__1.r, B(2,j).i = z__1.i;
	    }
	    i__2 = *n;
	    for (i = 3; i <= *n; ++i) {
		i__3 = i + j * b_dim1;
		i__4 = i + j * b_dim1;
		d_cnjg(&z__5, &DU(i - 1));
		i__5 = i - 1 + j * b_dim1;
		z__4.r = z__5.r * B(i-1,j).r - z__5.i * B(i-1,j).i, z__4.i = 
			z__5.r * B(i-1,j).i + z__5.i * B(i-1,j).r;
		z__3.r = B(i,j).r - z__4.r, z__3.i = B(i,j).i - z__4.i;
		d_cnjg(&z__7, &DU2(i - 2));
		i__6 = i - 2 + j * b_dim1;
		z__6.r = z__7.r * B(i-2,j).r - z__7.i * B(i-2,j).i, z__6.i = 
			z__7.r * B(i-2,j).i + z__7.i * B(i-2,j).r;
		z__2.r = z__3.r - z__6.r, z__2.i = z__3.i - z__6.i;
		d_cnjg(&z__8, &D(i));
		z_div(&z__1, &z__2, &z__8);
		B(i,j).r = z__1.r, B(i,j).i = z__1.i;
/* L70: */
	    }

/*           Solve L**H * x = b. */

	    for (i = *n - 1; i >= 1; --i) {
		if (IPIV(i) == i) {
		    i__2 = i + j * b_dim1;
		    i__3 = i + j * b_dim1;
		    d_cnjg(&z__3, &DL(i));
		    i__4 = i + 1 + j * b_dim1;
		    z__2.r = z__3.r * B(i+1,j).r - z__3.i * B(i+1,j).i, z__2.i =
			     z__3.r * B(i+1,j).i + z__3.i * B(i+1,j).r;
		    z__1.r = B(i,j).r - z__2.r, z__1.i = B(i,j).i - z__2.i;
		    B(i,j).r = z__1.r, B(i,j).i = z__1.i;
		} else {
		    i__2 = i + 1 + j * b_dim1;
		    temp.r = B(i+1,j).r, temp.i = B(i+1,j).i;
		    i__2 = i + 1 + j * b_dim1;
		    i__3 = i + j * b_dim1;
		    d_cnjg(&z__3, &DL(i));
		    z__2.r = z__3.r * temp.r - z__3.i * temp.i, z__2.i = 
			    z__3.r * temp.i + z__3.i * temp.r;
		    z__1.r = B(i,j).r - z__2.r, z__1.i = B(i,j).i - z__2.i;
		    B(i+1,j).r = z__1.r, B(i+1,j).i = z__1.i;
		    i__2 = i + j * b_dim1;
		    B(i,j).r = temp.r, B(i,j).i = temp.i;
		}
/* L80: */
	    }
/* L90: */
	}
    }

/*     End of ZGTTRS */

    return 0;
} /* zgttrs_ */

