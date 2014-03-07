#include "f2c.h"

/* Subroutine */ int zgelq2_(integer *m, integer *n, doublecomplex *a, 
	integer *lda, doublecomplex *tau, doublecomplex *work, integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    ZGELQ2 computes an LQ factorization of a complex m by n matrix A:   
    A = L * Q.   

    Arguments   
    =========   

    M       (input) INTEGER   
            The number of rows of the matrix A.  M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrix A.  N >= 0.   

    A       (input/output) COMPLEX*16 array, dimension (LDA,N)   
            On entry, the m by n matrix A.   
            On exit, the elements on and below the diagonal of the array 
  
            contain the m by min(m,n) lower trapezoidal matrix L (L is   
            lower triangular if m <= n); the elements above the diagonal, 
  
            with the array TAU, represent the unitary matrix Q as a   
            product of elementary reflectors (see Further Details).   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(1,M).   

    TAU     (output) COMPLEX*16 array, dimension (min(M,N))   
            The scalar factors of the elementary reflectors (see Further 
  
            Details).   

    WORK    (workspace) COMPLEX*16 array, dimension (M)   

    INFO    (output) INTEGER   
            = 0: successful exit   
            < 0: if INFO = -i, the i-th argument had an illegal value   

    Further Details   
    ===============   

    The matrix Q is represented as a product of elementary reflectors   

       Q = H(k)' . . . H(2)' H(1)', where k = min(m,n).   

    Each H(i) has the form   

       H(i) = I - tau * v * v'   

    where tau is a complex scalar, and v is a complex vector with   
    v(1:i-1) = 0 and v(i) = 1; conjg(v(i+1:n)) is stored on exit in   
    A(i,i+1:n), and tau in TAU(i).   

    ===================================================================== 
  


       Test the input arguments   

    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    /* Local variables */
    static integer i, k;
    static doublecomplex alpha;
    extern /* Subroutine */ int zlarf_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, doublecomplex *), xerbla_(char *, integer *), zlarfg_(integer *, doublecomplex *, doublecomplex *, 
	    integer *, doublecomplex *), zlacgv_(integer *, doublecomplex *, 
	    integer *);


#define TAU(I) tau[(I)-1]
#define WORK(I) work[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    *info = 0;
    if (*m < 0) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*lda < max(1,*m)) {
	*info = -4;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("ZGELQ2", &i__1);
	return 0;
    }

    k = min(*m,*n);

    i__1 = k;
    for (i = 1; i <= k; ++i) {

/*        Generate elementary reflector H(i) to annihilate A(i,i+1:n) 
*/

	i__2 = *n - i + 1;
	zlacgv_(&i__2, &A(i,i), lda);
	i__2 = i + i * a_dim1;
	alpha.r = A(i,i).r, alpha.i = A(i,i).i;
	i__2 = *n - i + 1;
/* Computing MIN */
	i__3 = i + 1;
	zlarfg_(&i__2, &alpha, &A(i,min(i+1,*n)), lda, &TAU(i));
	if (i < *m) {

/*           Apply H(i) to A(i+1:m,i:n) from the right */

	    i__2 = i + i * a_dim1;
	    A(i,i).r = 1., A(i,i).i = 0.;
	    i__2 = *m - i;
	    i__3 = *n - i + 1;
	    zlarf_("Right", &i__2, &i__3, &A(i,i), lda, &TAU(i), &
		    A(i+1,i), lda, &WORK(1));
	}
	i__2 = i + i * a_dim1;
	A(i,i).r = alpha.r, A(i,i).i = alpha.i;
	i__2 = *n - i + 1;
	zlacgv_(&i__2, &A(i,i), lda);
/* L10: */
    }
    return 0;

/*     End of ZGELQ2 */

} /* zgelq2_ */

