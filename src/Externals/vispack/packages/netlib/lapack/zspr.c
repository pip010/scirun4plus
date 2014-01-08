#include "f2c.h"

/* Subroutine */ int zspr_(char *uplo, integer *n, doublecomplex *alpha, 
	doublecomplex *x, integer *incx, doublecomplex *ap)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    ZSPR    performs the symmetric rank 1 operation   

       A := alpha*x*conjg( x' ) + A,   

    where alpha is a complex scalar, x is an n element vector and A is an 
  
    n by n symmetric matrix, supplied in packed form.   

    Arguments   
    ==========   

    UPLO   - CHARACTER*1   
             On entry, UPLO specifies whether the upper or lower   
             triangular part of the matrix A is supplied in the packed   
             array AP as follows:   

                UPLO = 'U' or 'u'   The upper triangular part of A is   
                                    supplied in AP.   

                UPLO = 'L' or 'l'   The lower triangular part of A is   
                                    supplied in AP.   

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

    AP     - COMPLEX*16 array, dimension at least   
             ( ( N*( N + 1 ) )/2 ).   
             Before entry, with  UPLO = 'U' or 'u', the array AP must   
             contain the upper triangular part of the symmetric matrix   
             packed sequentially, column by column, so that AP( 1 )   
             contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 1, 2 )   
             and a( 2, 2 ) respectively, and so on. On exit, the array   
             AP is overwritten by the upper triangular part of the   
             updated matrix.   
             Before entry, with UPLO = 'L' or 'l', the array AP must   
             contain the lower triangular part of the symmetric matrix   
             packed sequentially, column by column, so that AP( 1 )   
             contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 2, 1 )   
             and a( 3, 1 ) respectively, and so on. On exit, the array   
             AP is overwritten by the lower triangular part of the   
             updated matrix.   
             Note that the imaginary parts of the diagonal elements need 
  
             not be set, they are assumed to be zero, and on exit they   
             are set to zero.   

   ===================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5;
    doublecomplex z__1, z__2;
    /* Local variables */
    static integer info;
    static doublecomplex temp;
    static integer i, j, k;
    extern logical lsame_(char *, char *);
    static integer kk, ix, jx, kx;
    extern /* Subroutine */ int xerbla_(char *, integer *);


#define AP(I) ap[(I)-1]
#define X(I) x[(I)-1]


    info = 0;
    if (! lsame_(uplo, "U") && ! lsame_(uplo, "L")) {
	info = 1;
    } else if (*n < 0) {
	info = 2;
    } else if (*incx == 0) {
	info = 5;
    }
    if (info != 0) {
	xerbla_("ZSPR  ", &info);
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

/*     Start the operations. In this version the elements of the array AP 
  
       are accessed sequentially with one pass through AP. */

    kk = 1;
    if (lsame_(uplo, "U")) {

/*        Form  A  when upper triangle is stored in AP. */

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
		    k = kk;
		    i__2 = j - 1;
		    for (i = 1; i <= j-1; ++i) {
			i__3 = k;
			i__4 = k;
			i__5 = i;
			z__2.r = X(i).r * temp.r - X(i).i * temp.i, 
				z__2.i = X(i).r * temp.i + X(i).i * 
				temp.r;
			z__1.r = AP(k).r + z__2.r, z__1.i = AP(k).i + 
				z__2.i;
			AP(k).r = z__1.r, AP(k).i = z__1.i;
			++k;
/* L10: */
		    }
		    i__2 = kk + j - 1;
		    i__3 = kk + j - 1;
		    i__4 = j;
		    z__2.r = X(j).r * temp.r - X(j).i * temp.i, z__2.i =
			     X(j).r * temp.i + X(j).i * temp.r;
		    z__1.r = AP(kk+j-1).r + z__2.r, z__1.i = AP(kk+j-1).i + 
			    z__2.i;
		    AP(kk+j-1).r = z__1.r, AP(kk+j-1).i = z__1.i;
		} else {
		    i__2 = kk + j - 1;
		    i__3 = kk + j - 1;
		    AP(kk+j-1).r = AP(kk+j-1).r, AP(kk+j-1).i = AP(kk+j-1).i;
		}
		kk += j;
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
		    i__2 = kk + j - 2;
		    for (k = kk; k <= kk+j-2; ++k) {
			i__3 = k;
			i__4 = k;
			i__5 = ix;
			z__2.r = X(ix).r * temp.r - X(ix).i * temp.i, 
				z__2.i = X(ix).r * temp.i + X(ix).i * 
				temp.r;
			z__1.r = AP(k).r + z__2.r, z__1.i = AP(k).i + 
				z__2.i;
			AP(k).r = z__1.r, AP(k).i = z__1.i;
			ix += *incx;
/* L30: */
		    }
		    i__2 = kk + j - 1;
		    i__3 = kk + j - 1;
		    i__4 = jx;
		    z__2.r = X(jx).r * temp.r - X(jx).i * temp.i, z__2.i =
			     X(jx).r * temp.i + X(jx).i * temp.r;
		    z__1.r = AP(kk+j-1).r + z__2.r, z__1.i = AP(kk+j-1).i + 
			    z__2.i;
		    AP(kk+j-1).r = z__1.r, AP(kk+j-1).i = z__1.i;
		} else {
		    i__2 = kk + j - 1;
		    i__3 = kk + j - 1;
		    AP(kk+j-1).r = AP(kk+j-1).r, AP(kk+j-1).i = AP(kk+j-1).i;
		}
		jx += *incx;
		kk += j;
/* L40: */
	    }
	}
    } else {

/*        Form  A  when lower triangle is stored in AP. */

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
		    i__2 = kk;
		    i__3 = kk;
		    i__4 = j;
		    z__2.r = temp.r * X(j).r - temp.i * X(j).i, z__2.i =
			     temp.r * X(j).i + temp.i * X(j).r;
		    z__1.r = AP(kk).r + z__2.r, z__1.i = AP(kk).i + 
			    z__2.i;
		    AP(kk).r = z__1.r, AP(kk).i = z__1.i;
		    k = kk + 1;
		    i__2 = *n;
		    for (i = j + 1; i <= *n; ++i) {
			i__3 = k;
			i__4 = k;
			i__5 = i;
			z__2.r = X(i).r * temp.r - X(i).i * temp.i, 
				z__2.i = X(i).r * temp.i + X(i).i * 
				temp.r;
			z__1.r = AP(k).r + z__2.r, z__1.i = AP(k).i + 
				z__2.i;
			AP(k).r = z__1.r, AP(k).i = z__1.i;
			++k;
/* L50: */
		    }
		} else {
		    i__2 = kk;
		    i__3 = kk;
		    AP(kk).r = AP(kk).r, AP(kk).i = AP(kk).i;
		}
		kk = kk + *n - j + 1;
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
		    i__2 = kk;
		    i__3 = kk;
		    i__4 = jx;
		    z__2.r = temp.r * X(jx).r - temp.i * X(jx).i, z__2.i =
			     temp.r * X(jx).i + temp.i * X(jx).r;
		    z__1.r = AP(kk).r + z__2.r, z__1.i = AP(kk).i + 
			    z__2.i;
		    AP(kk).r = z__1.r, AP(kk).i = z__1.i;
		    ix = jx;
		    i__2 = kk + *n - j;
		    for (k = kk + 1; k <= kk+*n-j; ++k) {
			ix += *incx;
			i__3 = k;
			i__4 = k;
			i__5 = ix;
			z__2.r = X(ix).r * temp.r - X(ix).i * temp.i, 
				z__2.i = X(ix).r * temp.i + X(ix).i * 
				temp.r;
			z__1.r = AP(k).r + z__2.r, z__1.i = AP(k).i + 
				z__2.i;
			AP(k).r = z__1.r, AP(k).i = z__1.i;
/* L70: */
		    }
		} else {
		    i__2 = kk;
		    i__3 = kk;
		    AP(kk).r = AP(kk).r, AP(kk).i = AP(kk).i;
		}
		jx += *incx;
		kk = kk + *n - j + 1;
/* L80: */
	    }
	}
    }

    return 0;

/*     End of ZSPR */

} /* zspr_ */

