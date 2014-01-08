#include "f2c.h"

/* Subroutine */ int zlaswp_(integer *n, doublecomplex *a, integer *lda, 
	integer *k1, integer *k2, integer *ipiv, integer *incx)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    ZLASWP performs a series of row interchanges on the matrix A.   
    One row interchange is initiated for each of rows K1 through K2 of A. 
  

    Arguments   
    =========   

    N       (input) INTEGER   
            The number of columns of the matrix A.   

    A       (input/output) COMPLEX*16 array, dimension (LDA,N)   
            On entry, the matrix of column dimension N to which the row   
            interchanges will be applied.   
            On exit, the permuted matrix.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.   

    K1      (input) INTEGER   
            The first element of IPIV for which a row interchange will   
            be done.   

    K2      (input) INTEGER   
            The last element of IPIV for which a row interchange will   
            be done.   

    IPIV    (input) INTEGER array, dimension (M*abs(INCX))   
            The vector of pivot indices.  Only the elements in positions 
  
            K1 through K2 of IPIV are accessed.   
            IPIV(K) = L implies rows K and L are to be interchanged.   

    INCX    (input) INTEGER   
            The increment between successive values of IPIV.  If IPIV   
            is negative, the pivots are applied in reverse order.   

   ===================================================================== 
  


       Interchange row I with row IPIV(I) for each of rows K1 through K2. 
  

    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    integer a_dim1, a_offset, i__1;
    /* Local variables */
    static integer i;
    extern /* Subroutine */ int zswap_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    static integer ip, ix;


#define IPIV(I) ipiv[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    if (*incx == 0) {
	return 0;
    }
    if (*incx > 0) {
	ix = *k1;
    } else {
	ix = (1 - *k2) * *incx + 1;
    }
    if (*incx == 1) {
	i__1 = *k2;
	for (i = *k1; i <= *k2; ++i) {
	    ip = IPIV(i);
	    if (ip != i) {
		zswap_(n, &A(i,1), lda, &A(ip,1), lda);
	    }
/* L10: */
	}
    } else if (*incx > 1) {
	i__1 = *k2;
	for (i = *k1; i <= *k2; ++i) {
	    ip = IPIV(ix);
	    if (ip != i) {
		zswap_(n, &A(i,1), lda, &A(ip,1), lda);
	    }
	    ix += *incx;
/* L20: */
	}
    } else if (*incx < 0) {
	i__1 = *k1;
	for (i = *k2; i >= *k1; --i) {
	    ip = IPIV(ix);
	    if (ip != i) {
		zswap_(n, &A(i,1), lda, &A(ip,1), lda);
	    }
	    ix += *incx;
/* L30: */
	}
    }

    return 0;

/*     End of ZLASWP */

} /* zlaswp_ */

