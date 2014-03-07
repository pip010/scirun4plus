#include "f2c.h"

/* Subroutine */ int zgeql2_(integer *m, integer *n, doublecomplex *a, 
	integer *lda, doublecomplex *tau, doublecomplex *work, integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    ZGEQL2 computes a QL factorization of a complex m by n matrix A:   
    A = Q * L.   

    Arguments   
    =========   

    M       (input) INTEGER   
            The number of rows of the matrix A.  M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrix A.  N >= 0.   

    A       (input/output) COMPLEX*16 array, dimension (LDA,N)   
            On entry, the m by n matrix A.   
            On exit, if m >= n, the lower triangle of the subarray   
            A(m-n+1:m,1:n) contains the n by n lower triangular matrix L; 
  
            if m <= n, the elements on and below the (n-m)-th   
            superdiagonal contain the m by n lower trapezoidal matrix L; 
  
            the remaining elements, with the array TAU, represent the   
            unitary matrix Q as a product of elementary reflectors   
            (see Further Details).   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(1,M).   

    TAU     (output) COMPLEX*16 array, dimension (min(M,N))   
            The scalar factors of the elementary reflectors (see Further 
  
            Details).   

    WORK    (workspace) COMPLEX*16 array, dimension (N)   

    INFO    (output) INTEGER   
            = 0: successful exit   
            < 0: if INFO = -i, the i-th argument had an illegal value   

    Further Details   
    ===============   

    The matrix Q is represented as a product of elementary reflectors   

       Q = H(k) . . . H(2) H(1), where k = min(m,n).   

    Each H(i) has the form   

       H(i) = I - tau * v * v'   

    where tau is a complex scalar, and v is a complex vector with   
    v(m-k+i+1:m) = 0 and v(m-k+i) = 1; v(1:m-k+i-1) is stored on exit in 
  
    A(1:m-k+i-1,n-k+i), and tau in TAU(i).   

    ===================================================================== 
  


       Test the input arguments   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static integer c__1 = 1;
    
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    doublecomplex z__1;
    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);
    /* Local variables */
    static integer i, k;
    static doublecomplex alpha;
    extern /* Subroutine */ int zlarf_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, doublecomplex *), xerbla_(char *, integer *), zlarfg_(integer *, doublecomplex *, doublecomplex *, 
	    integer *, doublecomplex *);



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
	xerbla_("ZGEQL2", &i__1);
	return 0;
    }

    k = min(*m,*n);

    for (i = k; i >= 1; --i) {

/*        Generate elementary reflector H(i) to annihilate   
          A(1:m-k+i-1,n-k+i) */

	i__1 = *m - k + i + (*n - k + i) * a_dim1;
	alpha.r = A(*m-k+i,*n-k+i).r, alpha.i = A(*m-k+i,*n-k+i).i;
	i__1 = *m - k + i;
	zlarfg_(&i__1, &alpha, &A(1,*n-k+i), &c__1, &TAU(i));

/*        Apply H(i)' to A(1:m-k+i,1:n-k+i-1) from the left */

	i__1 = *m - k + i + (*n - k + i) * a_dim1;
	A(*m-k+i,*n-k+i).r = 1., A(*m-k+i,*n-k+i).i = 0.;
	i__1 = *m - k + i;
	i__2 = *n - k + i - 1;
	d_cnjg(&z__1, &TAU(i));
	zlarf_("Left", &i__1, &i__2, &A(1,*n-k+i), &c__1, &
		z__1, &A(1,1), lda, &WORK(1));
	i__1 = *m - k + i + (*n - k + i) * a_dim1;
	A(*m-k+i,*n-k+i).r = alpha.r, A(*m-k+i,*n-k+i).i = alpha.i;
/* L10: */
    }
    return 0;

/*     End of ZGEQL2 */

} /* zgeql2_ */

