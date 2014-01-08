#include "f2c.h"

/* Subroutine */ int zupgtr_(char *uplo, integer *n, doublecomplex *ap, 
	doublecomplex *tau, doublecomplex *q, integer *ldq, doublecomplex *
	work, integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    ZUPGTR generates a complex unitary matrix Q which is defined as the   
    product of n-1 elementary reflectors H(i) of order n, as returned by 
  
    ZHPTRD using packed storage:   

    if UPLO = 'U', Q = H(n-1) . . . H(2) H(1),   

    if UPLO = 'L', Q = H(1) H(2) . . . H(n-1).   

    Arguments   
    =========   

    UPLO    (input) CHARACTER*1   
            = 'U': Upper triangular packed storage used in previous   
                   call to ZHPTRD;   
            = 'L': Lower triangular packed storage used in previous   
                   call to ZHPTRD.   

    N       (input) INTEGER   
            The order of the matrix Q. N >= 0.   

    AP      (input) COMPLEX*16 array, dimension (N*(N+1)/2)   
            The vectors which define the elementary reflectors, as   
            returned by ZHPTRD.   

    TAU     (input) COMPLEX*16 array, dimension (N-1)   
            TAU(i) must contain the scalar factor of the elementary   
            reflector H(i), as returned by ZHPTRD.   

    Q       (output) COMPLEX*16 array, dimension (LDQ,N)   
            The N-by-N unitary matrix Q.   

    LDQ     (input) INTEGER   
            The leading dimension of the array Q. LDQ >= max(1,N).   

    WORK    (workspace) COMPLEX*16 array, dimension (N-1)   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   

    ===================================================================== 
  


       Test the input arguments   

    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    integer q_dim1, q_offset, i__1, i__2, i__3, i__4;
    /* Local variables */
    static integer i, j;
    extern logical lsame_(char *, char *);
    static integer iinfo;
    static logical upper;
    extern /* Subroutine */ int zung2l_(integer *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *);
    static integer ij;
    extern /* Subroutine */ int zung2r_(integer *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *), xerbla_(char *, integer *);


#define AP(I) ap[(I)-1]
#define TAU(I) tau[(I)-1]
#define WORK(I) work[(I)-1]

#define Q(I,J) q[(I)-1 + ((J)-1)* ( *ldq)]

    *info = 0;
    upper = lsame_(uplo, "U");
    if (! upper && ! lsame_(uplo, "L")) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*ldq < max(1,*n)) {
	*info = -6;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("ZUPGTR", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return 0;
    }

    if (upper) {

/*        Q was determined by a call to ZHPTRD with UPLO = 'U'   

          Unpack the vectors which define the elementary reflectors an
d   
          set the last row and column of Q equal to those of the unit 
  
          matrix */

	ij = 2;
	i__1 = *n - 1;
	for (j = 1; j <= *n-1; ++j) {
	    i__2 = j - 1;
	    for (i = 1; i <= j-1; ++i) {
		i__3 = i + j * q_dim1;
		i__4 = ij;
		Q(i,j).r = AP(ij).r, Q(i,j).i = AP(ij).i;
		++ij;
/* L10: */
	    }
	    ij += 2;
	    i__2 = *n + j * q_dim1;
	    Q(*n,j).r = 0., Q(*n,j).i = 0.;
/* L20: */
	}
	i__1 = *n - 1;
	for (i = 1; i <= *n-1; ++i) {
	    i__2 = i + *n * q_dim1;
	    Q(i,*n).r = 0., Q(i,*n).i = 0.;
/* L30: */
	}
	i__1 = *n + *n * q_dim1;
	Q(*n,*n).r = 1., Q(*n,*n).i = 0.;

/*        Generate Q(1:n-1,1:n-1) */

	i__1 = *n - 1;
	i__2 = *n - 1;
	i__3 = *n - 1;
	zung2l_(&i__1, &i__2, &i__3, &Q(1,1), ldq, &TAU(1), &WORK(1), &
		iinfo);

    } else {

/*        Q was determined by a call to ZHPTRD with UPLO = 'L'.   

          Unpack the vectors which define the elementary reflectors an
d   
          set the first row and column of Q equal to those of the unit
   
          matrix */

	i__1 = q_dim1 + 1;
	Q(1,1).r = 1., Q(1,1).i = 0.;
	i__1 = *n;
	for (i = 2; i <= *n; ++i) {
	    i__2 = i + q_dim1;
	    Q(i,1).r = 0., Q(i,1).i = 0.;
/* L40: */
	}
	ij = 3;
	i__1 = *n;
	for (j = 2; j <= *n; ++j) {
	    i__2 = j * q_dim1 + 1;
	    Q(1,j).r = 0., Q(1,j).i = 0.;
	    i__2 = *n;
	    for (i = j + 1; i <= *n; ++i) {
		i__3 = i + j * q_dim1;
		i__4 = ij;
		Q(i,j).r = AP(ij).r, Q(i,j).i = AP(ij).i;
		++ij;
/* L50: */
	    }
	    ij += 2;
/* L60: */
	}
	if (*n > 1) {

/*           Generate Q(2:n,2:n) */

	    i__1 = *n - 1;
	    i__2 = *n - 1;
	    i__3 = *n - 1;
	    zung2r_(&i__1, &i__2, &i__3, &Q(2,2), ldq, &TAU(1), 
		    &WORK(1), &iinfo);
	}
    }
    return 0;

/*     End of ZUPGTR */

} /* zupgtr_ */

