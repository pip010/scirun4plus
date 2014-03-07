#include "f2c.h"

/* Subroutine */ int zungl2_(integer *m, integer *n, integer *k, 
	doublecomplex *a, integer *lda, doublecomplex *tau, doublecomplex *
	work, integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    ZUNGL2 generates an m-by-n complex matrix Q with orthonormal rows,   
    which is defined as the first m rows of a product of k elementary   
    reflectors of order n   

          Q  =  H(k)' . . . H(2)' H(1)'   

    as returned by ZGELQF.   

    Arguments   
    =========   

    M       (input) INTEGER   
            The number of rows of the matrix Q. M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrix Q. N >= M.   

    K       (input) INTEGER   
            The number of elementary reflectors whose product defines the 
  
            matrix Q. M >= K >= 0.   

    A       (input/output) COMPLEX*16 array, dimension (LDA,N)   
            On entry, the i-th row must contain the vector which defines 
  
            the elementary reflector H(i), for i = 1,2,...,k, as returned 
  
            by ZGELQF in the first k rows of its array argument A.   
            On exit, the m by n matrix Q.   

    LDA     (input) INTEGER   
            The first dimension of the array A. LDA >= max(1,M).   

    TAU     (input) COMPLEX*16 array, dimension (K)   
            TAU(i) must contain the scalar factor of the elementary   
            reflector H(i), as returned by ZGELQF.   

    WORK    (workspace) COMPLEX*16 array, dimension (M)   

    INFO    (output) INTEGER   
            = 0: successful exit   
            < 0: if INFO = -i, the i-th argument has an illegal value   

    ===================================================================== 
  


       Test the input arguments   

    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublecomplex z__1, z__2;
    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);
    /* Local variables */
    static integer i, j, l;
    extern /* Subroutine */ int zscal_(integer *, doublecomplex *, 
	    doublecomplex *, integer *), zlarf_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, doublecomplex *), xerbla_(char *, integer *), zlacgv_(integer *, doublecomplex *, integer *);


#define TAU(I) tau[(I)-1]
#define WORK(I) work[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    *info = 0;
    if (*m < 0) {
	*info = -1;
    } else if (*n < *m) {
	*info = -2;
    } else if (*k < 0 || *k > *m) {
	*info = -3;
    } else if (*lda < max(1,*m)) {
	*info = -5;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("ZUNGL2", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*m <= 0) {
	return 0;
    }

    if (*k < *m) {

/*        Initialise rows k+1:m to rows of the unit matrix */

	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
	    i__2 = *m;
	    for (l = *k + 1; l <= *m; ++l) {
		i__3 = l + j * a_dim1;
		A(l,j).r = 0., A(l,j).i = 0.;
/* L10: */
	    }
	    if (j > *k && j <= *m) {
		i__2 = j + j * a_dim1;
		A(j,j).r = 1., A(j,j).i = 0.;
	    }
/* L20: */
	}
    }

    for (i = *k; i >= 1; --i) {

/*        Apply H(i)' to A(i:m,i:n) from the right */

	if (i < *n) {
	    i__1 = *n - i;
	    zlacgv_(&i__1, &A(i,i+1), lda);
	    if (i < *m) {
		i__1 = i + i * a_dim1;
		A(i,i).r = 1., A(i,i).i = 0.;
		i__1 = *m - i;
		i__2 = *n - i + 1;
		d_cnjg(&z__1, &TAU(i));
		zlarf_("Right", &i__1, &i__2, &A(i,i), lda, &z__1, 
			&A(i+1,i), lda, &WORK(1));
	    }
	    i__1 = *n - i;
	    i__2 = i;
	    z__1.r = -TAU(i).r, z__1.i = -TAU(i).i;
	    zscal_(&i__1, &z__1, &A(i,i+1), lda);
	    i__1 = *n - i;
	    zlacgv_(&i__1, &A(i,i+1), lda);
	}
	i__1 = i + i * a_dim1;
	d_cnjg(&z__2, &TAU(i));
	z__1.r = 1. - z__2.r, z__1.i = 0. - z__2.i;
	A(i,i).r = z__1.r, A(i,i).i = z__1.i;

/*        Set A(1:i-1,i) to zero */

	i__1 = i - 1;
	for (l = 1; l <= i-1; ++l) {
	    i__2 = i + l * a_dim1;
	    A(i,l).r = 0., A(i,l).i = 0.;
/* L30: */
	}
/* L40: */
    }
    return 0;

/*     End of ZUNGL2 */

} /* zungl2_ */

