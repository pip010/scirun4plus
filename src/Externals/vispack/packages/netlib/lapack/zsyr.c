#include "f2c.h"

/* Subroutine */ int zsyr_(char *uplo, integer *n, doublecomplex *alpha, 
	doublecomplex *x, integer *incx, doublecomplex *a, integer *lda)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    ZSYR   performs the symmetric rank 1 operation   

       A := alpha*x*( x' ) + A,   

    where alpha is a complex scalar, x is an n element vector and A is an 
  
    n by n symmetric matrix.   

    Arguments   
    ==========   

    UPLO   - CHARACTER*1   
             On entry, UPLO specifies whether the upper or lower   
             triangular part of the array A is to be referenced as   
             follows:   

                UPLO = 'U' or 'u'   Only the upper triangular part of A   
                                    is to be referenced.   

                UPLO = 'L' or 'l'   Only the lower triangular part of A   
                                    is to be referenced.   

             Unchanged on exit.   

    N      - INTEGER   
             On entry, N specifies the order of the matrix A.   
             N must be at least zero.   
             Unchanged on exit.   

    ALPHA  - COMPLEX*16   
             On entry, ALPHA specifies the scalar alpha.   
             Unchanged on exit.   

    X      - COMPLEX*16 array, dimension at least   
             ( 1 + ( N - 1 )*abs( INCX ) ).   
             Before entry, the incremented array X must contain the N-   
             element vector x.   
             Unchanged on exit.   

    INCX   - INTEGER   
             On entry, INCX specifies the increment for the elements of   
             X. INCX must not be zero.   
             Unchanged on exit.   

    A      - COMPLEX*16 array, dimension ( LDA, N )   
             Before entry, with  UPLO = 'U' or 'u', the leading n by n   
             upper triangular part of the array A must contain the upper 
  
             triangular part of the symmetric matrix and the strictly   
             lower triangular part of A is not referenced. On exit, the   
             upper triangular part of the array A is overwritten by the   
             upper triangular part of the updated matrix.   
             Before entry, with UPLO = 'L' or 'l', the leading n by n   
             lower triangular part of the array A must contain the lower 
  
             triangular part of the symmetric matrix and the strictly   
             upper triangular part of A is not referenced. On exit, the   
             lower triangular part of the array A is overwritten by the   
             lower triangular part of the updated matrix.   

    LDA    - INTEGER   
             On entry, LDA specifies the first dimension of A as declared 
  
             in the calling (sub) program. LDA must be at least   
             max( 1, N ).   
             Unchanged on exit.   

   ===================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5;
    doublecomplex z__1, z__2;
    /* Local variables */
    static integer info;
    static doublecomplex temp;
    static integer i, j;
    extern logical lsame_(char *, char *);
    static integer ix, jx, kx;
    extern /* Subroutine */ int xerbla_(char *, integer *);


#define X(I) x[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    info = 0;
    if (! lsame_(uplo, "U") && ! lsame_(uplo, "L")) {
	info = 1;
    } else if (*n < 0) {
	info = 2;
    } else if (*incx == 0) {
	info = 5;
    } else if (*lda < max(1,*n)) {
	info = 7;
    }
    if (info != 0) {
	xerbla_("ZSYR  ", &info);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0 || alpha->r == 0. && alpha->i == 0.) {
	return 0;
    }

/*     Set the start point in X if the increment is not unity. */

    if (*incx <= 0) {
	kx = 1 - (*n - 1) * *incx;
    } else if (*incx != 1) {
	kx = 1;
    }

/*     Start the operations. In this version the elements of A are   
       accessed sequentially with one pass through the triangular part   
       of A. */

    if (lsame_(uplo, "U")) {

/*        Form  A  when A is stored in upper triangle. */

	if (*incx == 1) {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		i__2 = j;
		if (X(j).r != 0. || X(j).i != 0.) {
		    i__2 = j;
		    z__1.r = alpha->r * X(j).r - alpha->i * X(j).i, 
			    z__1.i = alpha->r * X(j).i + alpha->i * X(j)
			    .r;
		    temp.r = z__1.r, temp.i = z__1.i;
		    i__2 = j;
		    for (i = 1; i <= j; ++i) {
			i__3 = i + j * a_dim1;
			i__4 = i + j * a_dim1;
			i__5 = i;
			z__2.r = X(i).r * temp.r - X(i).i * temp.i, 
				z__2.i = X(i).r * temp.i + X(i).i * 
				temp.r;
			z__1.r = A(i,j).r + z__2.r, z__1.i = A(i,j).i + 
				z__2.i;
			A(i,j).r = z__1.r, A(i,j).i = z__1.i;
/* L10: */
		    }
		}
/* L20: */
	    }
	} else {
	    jx = kx;
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		i__2 = jx;
		if (X(jx).r != 0. || X(jx).i != 0.) {
		    i__2 = jx;
		    z__1.r = alpha->r * X(jx).r - alpha->i * X(jx).i, 
			    z__1.i = alpha->r * X(jx).i + alpha->i * X(jx)
			    .r;
		    temp.r = z__1.r, temp.i = z__1.i;
		    ix = kx;
		    i__2 = j;
		    for (i = 1; i <= j; ++i) {
			i__3 = i + j * a_dim1;
			i__4 = i + j * a_dim1;
			i__5 = ix;
			z__2.r = X(ix).r * temp.r - X(ix).i * temp.i, 
				z__2.i = X(ix).r * temp.i + X(ix).i * 
				temp.r;
			z__1.r = A(i,j).r + z__2.r, z__1.i = A(i,j).i + 
				z__2.i;
			A(i,j).r = z__1.r, A(i,j).i = z__1.i;
			ix += *incx;
/* L30: */
		    }
		}
		jx += *incx;
/* L40: */
	    }
	}
    } else {

/*        Form  A  when A is stored in lower triangle. */

	if (*incx == 1) {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		i__2 = j;
		if (X(j).r != 0. || X(j).i != 0.) {
		    i__2 = j;
		    z__1.r = alpha->r * X(j).r - alpha->i * X(j).i, 
			    z__1.i = alpha->r * X(j).i + alpha->i * X(j)
			    .r;
		    temp.r = z__1.r, temp.i = z__1.i;
		    i__2 = *n;
		    for (i = j; i <= *n; ++i) {
			i__3 = i + j * a_dim1;
			i__4 = i + j * a_dim1;
			i__5 = i;
			z__2.r = X(i).r * temp.r - X(i).i * temp.i, 
				z__2.i = X(i).r * temp.i + X(i).i * 
				temp.r;
			z__1.r = A(i,j).r + z__2.r, z__1.i = A(i,j).i + 
				z__2.i;
			A(i,j).r = z__1.r, A(i,j).i = z__1.i;
/* L50: */
		    }
		}
/* L60: */
	    }
	} else {
	    jx = kx;
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		i__2 = jx;
		if (X(jx).r != 0. || X(jx).i != 0.) {
		    i__2 = jx;
		    z__1.r = alpha->r * X(jx).r - alpha->i * X(jx).i, 
			    z__1.i = alpha->r * X(jx).i + alpha->i * X(jx)
			    .r;
		    temp.r = z__1.r, temp.i = z__1.i;
		    ix = jx;
		    i__2 = *n;
		    for (i = j; i <= *n; ++i) {
			i__3 = i + j * a_dim1;
			i__4 = i + j * a_dim1;
			i__5 = ix;
			z__2.r = X(ix).r * temp.r - X(ix).i * temp.i, 
				z__2.i = X(ix).r * temp.i + X(ix).i * 
				temp.r;
			z__1.r = A(i,j).r + z__2.r, z__1.i = A(i,j).i + 
				z__2.i;
			A(i,j).r = z__1.r, A(i,j).i = z__1.i;
			ix += *incx;
/* L70: */
		    }
		}
		jx += *incx;
/* L80: */
	    }
	}
    }

    return 0;

/*     End of ZSYR */

} /* zsyr_ */

