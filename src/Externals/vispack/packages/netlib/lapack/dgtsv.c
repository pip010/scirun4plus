#include "f2c.h"

/* Subroutine */ int dgtsv_(integer *n, integer *nrhs, doublereal *dl, 
	doublereal *d, doublereal *du, doublereal *b, integer *ldb, integer *
	info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DGTSV  solves the equation   

       A*X = B,   

    where A is an N-by-N tridiagonal matrix, by Gaussian elimination with 
  
    partial pivoting.   

    Note that the equation  A'*X = B  may be solved by interchanging the 
  
    order of the arguments DU and DL.   

    Arguments   
    =========   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    NRHS    (input) INTEGER   
            The number of right hand sides, i.e., the number of columns   
            of the matrix B.  NRHS >= 0.   

    DL      (input/output) DOUBLE PRECISION array, dimension (N-1)   
            On entry, DL must contain the (n-1) subdiagonal elements of   
            A.   
            On exit, DL is overwritten by the (n-2) elements of the   
            second superdiagonal of the upper triangular matrix U from   
            the LU factorization of A, in DL(1), ..., DL(n-2).   

    D       (input/output) DOUBLE PRECISION array, dimension (N)   
            On entry, D must contain the diagonal elements of A.   
            On exit, D is overwritten by the n diagonal elements of U.   

    DU      (input/output) DOUBLE PRECISION array, dimension (N-1)   
            On entry, DU must contain the (n-1) superdiagonal elements   
            of A.   
            On exit, DU is overwritten by the (n-1) elements of the first 
  
            superdiagonal of U.   

    B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)   
            On entry, the N-by-NRHS right hand side matrix B.   
            On exit, if INFO = 0, the N-by-NRHS solution matrix X.   

    LDB     (input) INTEGER   
            The leading dimension of the array B.  LDB >= max(1,N).   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   
            > 0:  if INFO = i, U(i,i) is exactly zero, and the solution   
                  has not been computed.  The factorization has not been 
  
                  completed unless i = N.   

    ===================================================================== 
  


    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    integer b_dim1, b_offset, i__1, i__2;
    doublereal d__1, d__2;
    /* Local variables */
    static doublereal temp, mult;
    static integer j, k;
    extern /* Subroutine */ int xerbla_(char *, integer *);


#define DL(I) dl[(I)-1]
#define D(I) d[(I)-1]
#define DU(I) du[(I)-1]

#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]

    *info = 0;
    if (*n < 0) {
	*info = -1;
    } else if (*nrhs < 0) {
	*info = -2;
    } else if (*ldb < max(1,*n)) {
	*info = -7;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DGTSV ", &i__1);
	return 0;
    }

    if (*n == 0) {
	return 0;
    }

    i__1 = *n - 1;
    for (k = 1; k <= *n-1; ++k) {
	if (DL(k) == 0.) {

/*           Subdiagonal is zero, no elimination is required. */

	    if (D(k) == 0.) {

/*              Diagonal is zero: set INFO = K and return; a u
nique   
                solution can not be found. */

		*info = k;
		return 0;
	    }
	} else if ((d__1 = D(k), abs(d__1)) >= (d__2 = DL(k), abs(d__2))) {

/*           No row interchange required */

	    mult = DL(k) / D(k);
	    D(k + 1) -= mult * DU(k);
	    i__2 = *nrhs;
	    for (j = 1; j <= *nrhs; ++j) {
		B(k+1,j) -= mult * B(k,j);
/* L10: */
	    }
	    if (k < *n - 1) {
		DL(k) = 0.;
	    }
	} else {

/*           Interchange rows K and K+1 */

	    mult = D(k) / DL(k);
	    D(k) = DL(k);
	    temp = D(k + 1);
	    D(k + 1) = DU(k) - mult * temp;
	    if (k < *n - 1) {
		DL(k) = DU(k + 1);
		DU(k + 1) = -mult * DL(k);
	    }
	    DU(k) = temp;
	    i__2 = *nrhs;
	    for (j = 1; j <= *nrhs; ++j) {
		temp = B(k,j);
		B(k,j) = B(k+1,j);
		B(k+1,j) = temp - mult * B(k+1,j);
/* L20: */
	    }
	}
/* L30: */
    }
    if (D(*n) == 0.) {
	*info = *n;
	return 0;
    }

/*     Back solve with the matrix U from the factorization. */

    i__1 = *nrhs;
    for (j = 1; j <= *nrhs; ++j) {
	B(*n,j) /= D(*n);
	if (*n > 1) {
	    B(*n-1,j) = (B(*n-1,j) - DU(*n - 1) * B(*n,j)) / D(*n - 1);
	}
	for (k = *n - 2; k >= 1; --k) {
	    B(k,j) = (B(k,j) - DU(k) * B(k+1,j) - DL(k) * B(k+2,j)) / D(k);
/* L40: */
	}
/* L50: */
    }

    return 0;

/*     End of DGTSV */

} /* dgtsv_ */

