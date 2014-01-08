#include "f2c.h"

/* Subroutine */ int zungtr_(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, doublecomplex *tau, doublecomplex *work, integer *lwork,
	 integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    ZUNGTR generates a complex unitary matrix Q which is defined as the   
    product of n-1 elementary reflectors of order N, as returned by   
    ZHETRD:   

    if UPLO = 'U', Q = H(n-1) . . . H(2) H(1),   

    if UPLO = 'L', Q = H(1) H(2) . . . H(n-1).   

    Arguments   
    =========   

    UPLO    (input) CHARACTER*1   
            = 'U': Upper triangle of A contains elementary reflectors   
                   from ZHETRD;   
            = 'L': Lower triangle of A contains elementary reflectors   
                   from ZHETRD.   

    N       (input) INTEGER   
            The order of the matrix Q. N >= 0.   

    A       (input/output) COMPLEX*16 array, dimension (LDA,N)   
            On entry, the vectors which define the elementary reflectors, 
  
            as returned by ZHETRD.   
            On exit, the N-by-N unitary matrix Q.   

    LDA     (input) INTEGER   
            The leading dimension of the array A. LDA >= N.   

    TAU     (input) COMPLEX*16 array, dimension (N-1)   
            TAU(i) must contain the scalar factor of the elementary   
            reflector H(i), as returned by ZHETRD.   

    WORK    (workspace/output) COMPLEX*16 array, dimension (LWORK)   
            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INTEGER   
            The dimension of the array WORK. LWORK >= N-1.   
            For optimum performance LWORK >= (N-1)*NB, where NB is   
            the optimal blocksize.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   

    ===================================================================== 
  


       Test the input arguments   

    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;
    /* Local variables */
    static integer i, j;
    extern logical lsame_(char *, char *);
    static integer iinfo;
    static logical upper;
    extern /* Subroutine */ int xerbla_(char *, integer *), zungql_(
	    integer *, integer *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, integer *), zungqr_(
	    integer *, integer *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, integer *);


#define TAU(I) tau[(I)-1]
#define WORK(I) work[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    *info = 0;
    upper = lsame_(uplo, "U");
    if (! upper && ! lsame_(uplo, "L")) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*lda < max(1,*n)) {
	*info = -4;
    } else /* if(complicated condition) */ {
/* Computing MAX */
	i__1 = 1, i__2 = *n - 1;
	if (*lwork < max(i__1,i__2)) {
	    *info = -7;
	}
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("ZUNGTR", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0) {
	WORK(1).r = 1., WORK(1).i = 0.;
	return 0;
    }

    if (upper) {

/*        Q was determined by a call to ZHETRD with UPLO = 'U'   

          Shift the vectors which define the elementary reflectors one
   
          column to the left, and set the last row and column of Q to 
  
          those of the unit matrix */

	i__1 = *n - 1;
	for (j = 1; j <= *n-1; ++j) {
	    i__2 = j - 1;
	    for (i = 1; i <= j-1; ++i) {
		i__3 = i + j * a_dim1;
		i__4 = i + (j + 1) * a_dim1;
		A(i,j).r = A(i,j+1).r, A(i,j).i = A(i,j+1).i;
/* L10: */
	    }
	    i__2 = *n + j * a_dim1;
	    A(*n,j).r = 0., A(*n,j).i = 0.;
/* L20: */
	}
	i__1 = *n - 1;
	for (i = 1; i <= *n-1; ++i) {
	    i__2 = i + *n * a_dim1;
	    A(i,*n).r = 0., A(i,*n).i = 0.;
/* L30: */
	}
	i__1 = *n + *n * a_dim1;
	A(*n,*n).r = 1., A(*n,*n).i = 0.;

/*        Generate Q(1:n-1,1:n-1) */

	i__1 = *n - 1;
	i__2 = *n - 1;
	i__3 = *n - 1;
	zungql_(&i__1, &i__2, &i__3, &A(1,1), lda, &TAU(1), &WORK(1), 
		lwork, &iinfo);

    } else {

/*        Q was determined by a call to ZHETRD with UPLO = 'L'.   

          Shift the vectors which define the elementary reflectors one
   
          column to the right, and set the first row and column of Q t
o   
          those of the unit matrix */

	for (j = *n; j >= 2; --j) {
	    i__1 = j * a_dim1 + 1;
	    A(1,j).r = 0., A(1,j).i = 0.;
	    i__1 = *n;
	    for (i = j + 1; i <= *n; ++i) {
		i__2 = i + j * a_dim1;
		i__3 = i + (j - 1) * a_dim1;
		A(i,j).r = A(i,j-1).r, A(i,j).i = A(i,j-1).i;
/* L40: */
	    }
/* L50: */
	}
	i__1 = a_dim1 + 1;
	A(1,1).r = 1., A(1,1).i = 0.;
	i__1 = *n;
	for (i = 2; i <= *n; ++i) {
	    i__2 = i + a_dim1;
	    A(i,1).r = 0., A(i,1).i = 0.;
/* L60: */
	}
	if (*n > 1) {

/*           Generate Q(2:n,2:n) */

	    i__1 = *n - 1;
	    i__2 = *n - 1;
	    i__3 = *n - 1;
	    zungqr_(&i__1, &i__2, &i__3, &A(2,2), lda, &TAU(1), 
		    &WORK(1), lwork, &iinfo);
	}
    }
    return 0;

/*     End of ZUNGTR */

} /* zungtr_ */

