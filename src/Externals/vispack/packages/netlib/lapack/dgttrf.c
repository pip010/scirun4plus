#include "f2c.h"

/* Subroutine */ int dgttrf_(integer *n, doublereal *dl, doublereal *d, 
	doublereal *du, doublereal *du2, integer *ipiv, integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DGTTRF computes an LU factorization of a real tridiagonal matrix A   
    using elimination with partial pivoting and row interchanges.   

    The factorization has the form   
       A = L * U   
    where L is a product of permutation and unit lower bidiagonal   
    matrices and U is upper triangular with nonzeros in only the main   
    diagonal and first two superdiagonals.   

    Arguments   
    =========   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    DL      (input/output) DOUBLE PRECISION array, dimension (N-1)   
            On entry, DL must contain the (n-1) subdiagonal elements of   
            A.   
            On exit, DL is overwritten by the (n-1) multipliers that   
            define the matrix L from the LU factorization of A.   

    D       (input/output) DOUBLE PRECISION array, dimension (N)   
            On entry, D must contain the diagonal elements of A.   
            On exit, D is overwritten by the n diagonal elements of the   
            upper triangular matrix U from the LU factorization of A.   

    DU      (input/output) DOUBLE PRECISION array, dimension (N-1)   
            On entry, DU must contain the (n-1) superdiagonal elements   
            of A.   
            On exit, DU is overwritten by the (n-1) elements of the first 
  
            superdiagonal of U.   

    DU2     (output) DOUBLE PRECISION array, dimension (N-2)   
            On exit, DU2 is overwritten by the (n-2) elements of the   
            second superdiagonal of U.   

    IPIV    (output) INTEGER array, dimension (N)   
            The pivot indices; for 1 <= i <= n, row i of the matrix was   
            interchanged with row IPIV(i).  IPIV(i) will always be either 
  
            i or i+1; IPIV(i) = i indicates a row interchange was not   
            required.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   
            > 0:  if INFO = i, U(i,i) is exactly zero. The factorization 
  
                  has been completed, but the factor U is exactly   
                  singular, and division by zero will occur if it is used 
  
                  to solve a system of equations.   

    ===================================================================== 
  


    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;
    /* Local variables */
    static doublereal fact, temp;
    static integer i;
    extern /* Subroutine */ int xerbla_(char *, integer *);


#define IPIV(I) ipiv[(I)-1]
#define DU2(I) du2[(I)-1]
#define DU(I) du[(I)-1]
#define D(I) d[(I)-1]
#define DL(I) dl[(I)-1]


    *info = 0;
    if (*n < 0) {
	*info = -1;
	i__1 = -(*info);
	xerbla_("DGTTRF", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return 0;
    }

/*     Initialize IPIV(i) = i */

    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
	IPIV(i) = i;
/* L10: */
    }

    i__1 = *n - 1;
    for (i = 1; i <= *n-1; ++i) {
	if (DL(i) == 0.) {

/*           Subdiagonal is zero, no elimination is required. */

	    if (D(i) == 0. && *info == 0) {
		*info = i;
	    }
	    if (i < *n - 1) {
		DU2(i) = 0.;
	    }
	} else if ((d__1 = D(i), abs(d__1)) >= (d__2 = DL(i), abs(d__2))) {

/*           No row interchange required, eliminate DL(I) */

	    fact = DL(i) / D(i);
	    DL(i) = fact;
	    D(i + 1) -= fact * DU(i);
	    if (i < *n - 1) {
		DU2(i) = 0.;
	    }
	} else {

/*           Interchange rows I and I+1, eliminate DL(I) */

	    fact = D(i) / DL(i);
	    D(i) = DL(i);
	    DL(i) = fact;
	    temp = DU(i);
	    DU(i) = D(i + 1);
	    D(i + 1) = temp - fact * D(i + 1);
	    if (i < *n - 1) {
		DU2(i) = DU(i + 1);
		DU(i + 1) = -fact * DU(i + 1);
	    }
	    ++IPIV(i);
	}
/* L20: */
    }
    if (D(*n) == 0. && *info == 0) {
	*info = *n;
	return 0;
    }

    return 0;

/*     End of DGTTRF */

} /* dgttrf_ */

