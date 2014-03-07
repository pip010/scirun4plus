#include "f2c.h"

/* Subroutine */ int dgttrs_(char *trans, integer *n, integer *nrhs, 
	doublereal *dl, doublereal *d, doublereal *du, doublereal *du2, 
	integer *ipiv, doublereal *b, integer *ldb, integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DGTTRS solves one of the systems of equations   
       A*X = B  or  A'*X = B,   
    with a tridiagonal matrix A using the LU factorization computed   
    by DGTTRF.   

    Arguments   
    =========   

    TRANS   (input) CHARACTER   
            Specifies the form of the system of equations:   
            = 'N':  A * X = B  (No transpose)   
            = 'T':  A'* X = B  (Transpose)   
            = 'C':  A'* X = B  (Conjugate transpose = Transpose)   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    NRHS    (input) INTEGER   
            The number of right hand sides, i.e., the number of columns   
            of the matrix B.  NRHS >= 0.   

    DL      (input) DOUBLE PRECISION array, dimension (N-1)   
            The (n-1) multipliers that define the matrix L from the   
            LU factorization of A.   

    D       (input) DOUBLE PRECISION array, dimension (N)   
            The n diagonal elements of the upper triangular matrix U from 
  
            the LU factorization of A.   

    DU      (input) DOUBLE PRECISION array, dimension (N-1)   
            The (n-1) elements of the first superdiagonal of U.   

    DU2     (input) DOUBLE PRECISION array, dimension (N-2)   
            The (n-2) elements of the second superdiagonal of U.   

    IPIV    (input) INTEGER array, dimension (N)   
            The pivot indices; for 1 <= i <= n, row i of the matrix was   
            interchanged with row IPIV(i).  IPIV(i) will always be either 
  
            i or i+1; IPIV(i) = i indicates a row interchange was not   
            required.   

    B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)   
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
    integer b_dim1, b_offset, i__1, i__2;
    /* Local variables */
    static doublereal temp;
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
	xerbla_("DGTTRS", &i__1);
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
		    B(i+1,j) -= DL(i) * B(i,j);
		} else {
		    temp = B(i,j);
		    B(i,j) = B(i+1,j);
		    B(i+1,j) = temp - DL(i) * B(i,j);
		}
/* L10: */
	    }

/*           Solve U*x = b. */

	    B(*n,j) /= D(*n);
	    if (*n > 1) {
		B(*n-1,j) = (B(*n-1,j) - DU(*n - 1) 
			* B(*n,j)) / D(*n - 1);
	    }
	    for (i = *n - 2; i >= 1; --i) {
		B(i,j) = (B(i,j) - DU(i) * B(i+1,j) - DU2(i) * B(i+2,j)) / D(i);
/* L20: */
	    }
/* L30: */
	}
    } else {

/*        Solve A' * X = B. */

	i__1 = *nrhs;
	for (j = 1; j <= *nrhs; ++j) {

/*           Solve U'*x = b. */

	    B(1,j) /= D(1);
	    if (*n > 1) {
		B(2,j) = (B(2,j) - DU(1) * B(1,j)) / D(2);
	    }
	    i__2 = *n;
	    for (i = 3; i <= *n; ++i) {
		B(i,j) = (B(i,j) - DU(i - 1) * B(i-1,j) - DU2(i - 2) * B(i-2,j)) / 
			D(i);
/* L40: */
	    }

/*           Solve L'*x = b. */

	    for (i = *n - 1; i >= 1; --i) {
		if (IPIV(i) == i) {
		    B(i,j) -= DL(i) * B(i+1,j);
		} else {
		    temp = B(i+1,j);
		    B(i+1,j) = B(i,j) - DL(i) * temp;
		    B(i,j) = temp;
		}
/* L50: */
	    }
/* L60: */
	}
    }

/*     End of DGTTRS */

    return 0;
} /* dgttrs_ */

