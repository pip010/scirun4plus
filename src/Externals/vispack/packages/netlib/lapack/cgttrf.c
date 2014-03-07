#include "f2c.h"

/* Subroutine */ int cgttrf_(integer *n, complex *dl, complex *d, complex *du,
	 complex *du2, integer *ipiv, integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    CGTTRF computes an LU factorization of a complex tridiagonal matrix A 
  
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

    DL      (input/output) COMPLEX array, dimension (N-1)   
            On entry, DL must contain the (n-1) subdiagonal elements of   
            A.   
            On exit, DL is overwritten by the (n-1) multipliers that   
            define the matrix L from the LU factorization of A.   

    D       (input/output) COMPLEX array, dimension (N)   
            On entry, D must contain the diagonal elements of A.   
            On exit, D is overwritten by the n diagonal elements of the   
            upper triangular matrix U from the LU factorization of A.   

    DU      (input/output) COMPLEX array, dimension (N-1)   
            On entry, DU must contain the (n-1) superdiagonal elements   
            of A.   
            On exit, DU is overwritten by the (n-1) elements of the first 
  
            superdiagonal of U.   

    DU2     (output) COMPLEX array, dimension (N-2)   
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
    integer i__1, i__2, i__3, i__4;
    real r__1, r__2, r__3, r__4;
    complex q__1, q__2;
    /* Builtin functions */
    double r_imag(complex *);
    void c_div(complex *, complex *, complex *);
    /* Local variables */
    static complex fact, temp;
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
	xerbla_("CGTTRF", &i__1);
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
	i__2 = i;
	if (DL(i).r == 0.f && DL(i).i == 0.f) {

/*           Subdiagonal is zero, no elimination is required. */

	    i__2 = i;
	    if (D(i).r == 0.f && D(i).i == 0.f && *info == 0) {
		*info = i;
	    }
	    if (i < *n - 1) {
		i__2 = i;
		DU2(i).r = 0.f, DU2(i).i = 0.f;
	    }
	} else /* if(complicated condition) */ {
	    i__2 = i;
	    i__3 = i;
	    if ((r__1 = D(i).r, dabs(r__1)) + (r__2 = r_imag(&D(i)), dabs(
		    r__2)) >= (r__3 = DL(i).r, dabs(r__3)) + (r__4 = 
		    r_imag(&DL(i)), dabs(r__4))) {

/*           No row interchange required, eliminate DL(I) */

		c_div(&q__1, &DL(i), &D(i));
		fact.r = q__1.r, fact.i = q__1.i;
		i__2 = i;
		DL(i).r = fact.r, DL(i).i = fact.i;
		i__2 = i + 1;
		i__3 = i + 1;
		i__4 = i;
		q__2.r = fact.r * DU(i).r - fact.i * DU(i).i, q__2.i = 
			fact.r * DU(i).i + fact.i * DU(i).r;
		q__1.r = D(i+1).r - q__2.r, q__1.i = D(i+1).i - q__2.i;
		D(i+1).r = q__1.r, D(i+1).i = q__1.i;
		if (i < *n - 1) {
		    i__2 = i;
		    DU2(i).r = 0.f, DU2(i).i = 0.f;
		}
	    } else {

/*           Interchange rows I and I+1, eliminate DL(I) */

		c_div(&q__1, &D(i), &DL(i));
		fact.r = q__1.r, fact.i = q__1.i;
		i__2 = i;
		i__3 = i;
		D(i).r = DL(i).r, D(i).i = DL(i).i;
		i__2 = i;
		DL(i).r = fact.r, DL(i).i = fact.i;
		i__2 = i;
		temp.r = DU(i).r, temp.i = DU(i).i;
		i__2 = i;
		i__3 = i + 1;
		DU(i).r = D(i+1).r, DU(i).i = D(i+1).i;
		i__2 = i + 1;
		i__3 = i + 1;
		q__2.r = fact.r * D(i+1).r - fact.i * D(i+1).i, q__2.i = 
			fact.r * D(i+1).i + fact.i * D(i+1).r;
		q__1.r = temp.r - q__2.r, q__1.i = temp.i - q__2.i;
		D(i+1).r = q__1.r, D(i+1).i = q__1.i;
		if (i < *n - 1) {
		    i__2 = i;
		    i__3 = i + 1;
		    DU2(i).r = DU(i+1).r, DU2(i).i = DU(i+1).i;
		    i__2 = i + 1;
		    q__2.r = -(doublereal)fact.r, q__2.i = -(doublereal)
			    fact.i;
		    i__3 = i + 1;
		    q__1.r = q__2.r * DU(i+1).r - q__2.i * DU(i+1).i, 
			    q__1.i = q__2.r * DU(i+1).i + q__2.i * DU(i+1)
			    .r;
		    DU(i+1).r = q__1.r, DU(i+1).i = q__1.i;
		}
		++IPIV(i);
	    }
	}
/* L20: */
    }
    i__1 = *n;
    if (D(*n).r == 0.f && D(*n).i == 0.f && *info == 0) {
	*info = *n;
	return 0;
    }

    return 0;

/*     End of CGTTRF */

} /* cgttrf_ */

