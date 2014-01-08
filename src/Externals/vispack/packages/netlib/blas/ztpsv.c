
/*  -- translated by f2c (version 19940927).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Subroutine */ int ztpsv_(char *uplo, char *trans, char *diag, integer *n, 
	doublecomplex *ap, doublecomplex *x, integer *incx)
{


    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *), d_cnjg(
	    doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer info;
    static doublecomplex temp;
    static integer i, j, k;
    extern logical lsame_(char *, char *);
    static integer kk, ix, jx, kx;
    extern /* Subroutine */ int xerbla_(char *, integer *);
    static logical noconj, nounit;


/*  Purpose   
    =======   

    ZTPSV  solves one of the systems of equations   

       A*x = b,   or   A'*x = b,   or   conjg( A' )*x = b,   

    where b and x are n element vectors and A is an n by n unit, or   
    non-unit, upper or lower triangular matrix, supplied in packed form. 
  

    No test for singularity or near-singularity is included in this   
    routine. Such tests must be performed before calling this routine.   

    Parameters   
    ==========   

    UPLO   - CHARACTER*1.   
             On entry, UPLO specifies whether the matrix is an upper or   
             lower triangular matrix as follows:   

                UPLO = 'U' or 'u'   A is an upper triangular matrix.   

                UPLO = 'L' or 'l'   A is a lower triangular matrix.   

             Unchanged on exit.   

    TRANS  - CHARACTER*1.   
             On entry, TRANS specifies the equations to be solved as   
             follows:   

                TRANS = 'N' or 'n'   A*x = b.   

                TRANS = 'T' or 't'   A'*x = b.   

                TRANS = 'C' or 'c'   conjg( A' )*x = b.   

             Unchanged on exit.   

    DIAG   - CHARACTER*1.   
             On entry, DIAG specifies whether or not A is unit   
             triangular as follows:   

                DIAG = 'U' or 'u'   A is assumed to be unit triangular.   

                DIAG = 'N' or 'n'   A is not assumed to be unit   
                                    triangular.   

             Unchanged on exit.   

    N      - INTEGER.   
             On entry, N specifies the order of the matrix A.   
             N must be at least zero.   
             Unchanged on exit.   

    AP     - COMPLEX*16       array of DIMENSION at least   
             ( ( n*( n + 1 ) )/2 ).   
             Before entry with  UPLO = 'U' or 'u', the array AP must   
             contain the upper triangular matrix packed sequentially,   
             column by column, so that AP( 1 ) contains a( 1, 1 ),   
             AP( 2 ) and AP( 3 ) contain a( 1, 2 ) and a( 2, 2 )   
             respectively, and so on.   
             Before entry with UPLO = 'L' or 'l', the array AP must   
             contain the lower triangular matrix packed sequentially,   
             column by column, so that AP( 1 ) contains a( 1, 1 ),   
             AP( 2 ) and AP( 3 ) contain a( 2, 1 ) and a( 3, 1 )   
             respectively, and so on.   
             Note that when  DIAG = 'U' or 'u', the diagonal elements of 
  
             A are not referenced, but are assumed to be unity.   
             Unchanged on exit.   

    X      - COMPLEX*16       array of dimension at least   
             ( 1 + ( n - 1 )*abs( INCX ) ).   
             Before entry, the incremented array X must contain the n   
             element right-hand side vector b. On exit, X is overwritten 
  
             with the solution vector x.   

    INCX   - INTEGER.   
             On entry, INCX specifies the increment for the elements of   
             X. INCX must not be zero.   
             Unchanged on exit.   


    Level 2 Blas routine.   

    -- Written on 22-October-1986.   
       Jack Dongarra, Argonne National Lab.   
       Jeremy Du Croz, Nag Central Office.   
       Sven Hammarling, Nag Central Office.   
       Richard Hanson, Sandia National Labs.   



       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
#define X(I) x[(I)-1]
#define AP(I) ap[(I)-1]


    info = 0;
    if (! lsame_(uplo, "U") && ! lsame_(uplo, "L")) {
	info = 1;
    } else if (! lsame_(trans, "N") && ! lsame_(trans, "T") &&
	     ! lsame_(trans, "C")) {
	info = 2;
    } else if (! lsame_(diag, "U") && ! lsame_(diag, "N")) {
	info = 3;
    } else if (*n < 0) {
	info = 4;
    } else if (*incx == 0) {
	info = 7;
    }
    if (info != 0) {
	xerbla_("ZTPSV ", &info);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0) {
	return 0;
    }

    noconj = lsame_(trans, "T");
    nounit = lsame_(diag, "N");

/*     Set up the start point in X if the increment is not unity. This   
       will be  ( N - 1 )*INCX  too small for descending loops. */

    if (*incx <= 0) {
	kx = 1 - (*n - 1) * *incx;
    } else if (*incx != 1) {
	kx = 1;
    }

/*     Start the operations. In this version the elements of AP are   
       accessed sequentially with one pass through AP. */

    if (lsame_(trans, "N")) {

/*        Form  x := inv( A )*x. */

	if (lsame_(uplo, "U")) {
	    kk = *n * (*n + 1) / 2;
	    if (*incx == 1) {
		for (j = *n; j >= 1; --j) {
		    i__1 = j;
		    if (X(j).r != 0. || X(j).i != 0.) {
			if (nounit) {
			    i__1 = j;
			    z_div(&z__1, &X(j), &AP(kk));
			    X(j).r = z__1.r, X(j).i = z__1.i;
			}
			i__1 = j;
			temp.r = X(j).r, temp.i = X(j).i;
			k = kk - 1;
			for (i = j - 1; i >= 1; --i) {
			    i__1 = i;
			    i__2 = i;
			    i__3 = k;
			    z__2.r = temp.r * AP(k).r - temp.i * AP(k)
				    .i, z__2.i = temp.r * AP(k).i + temp.i 
				    * AP(k).r;
			    z__1.r = X(i).r - z__2.r, z__1.i = X(i).i - 
				    z__2.i;
			    X(i).r = z__1.r, X(i).i = z__1.i;
			    --k;
/* L10: */
			}
		    }
		    kk -= j;
/* L20: */
		}
	    } else {
		jx = kx + (*n - 1) * *incx;
		for (j = *n; j >= 1; --j) {
		    i__1 = jx;
		    if (X(jx).r != 0. || X(jx).i != 0.) {
			if (nounit) {
			    i__1 = jx;
			    z_div(&z__1, &X(jx), &AP(kk));
			    X(jx).r = z__1.r, X(jx).i = z__1.i;
			}
			i__1 = jx;
			temp.r = X(jx).r, temp.i = X(jx).i;
			ix = jx;
			i__1 = kk - j + 1;
			for (k = kk - 1; k >= kk-j+1; --k) {
			    ix -= *incx;
			    i__2 = ix;
			    i__3 = ix;
			    i__4 = k;
			    z__2.r = temp.r * AP(k).r - temp.i * AP(k)
				    .i, z__2.i = temp.r * AP(k).i + temp.i 
				    * AP(k).r;
			    z__1.r = X(ix).r - z__2.r, z__1.i = X(ix).i - 
				    z__2.i;
			    X(ix).r = z__1.r, X(ix).i = z__1.i;
/* L30: */
			}
		    }
		    jx -= *incx;
		    kk -= j;
/* L40: */
		}
	    }
	} else {
	    kk = 1;
	    if (*incx == 1) {
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
		    i__2 = j;
		    if (X(j).r != 0. || X(j).i != 0.) {
			if (nounit) {
			    i__2 = j;
			    z_div(&z__1, &X(j), &AP(kk));
			    X(j).r = z__1.r, X(j).i = z__1.i;
			}
			i__2 = j;
			temp.r = X(j).r, temp.i = X(j).i;
			k = kk + 1;
			i__2 = *n;
			for (i = j + 1; i <= *n; ++i) {
			    i__3 = i;
			    i__4 = i;
			    i__5 = k;
			    z__2.r = temp.r * AP(k).r - temp.i * AP(k)
				    .i, z__2.i = temp.r * AP(k).i + temp.i 
				    * AP(k).r;
			    z__1.r = X(i).r - z__2.r, z__1.i = X(i).i - 
				    z__2.i;
			    X(i).r = z__1.r, X(i).i = z__1.i;
			    ++k;
/* L50: */
			}
		    }
		    kk += *n - j + 1;
/* L60: */
		}
	    } else {
		jx = kx;
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
		    i__2 = jx;
		    if (X(jx).r != 0. || X(jx).i != 0.) {
			if (nounit) {
			    i__2 = jx;
			    z_div(&z__1, &X(jx), &AP(kk));
			    X(jx).r = z__1.r, X(jx).i = z__1.i;
			}
			i__2 = jx;
			temp.r = X(jx).r, temp.i = X(jx).i;
			ix = jx;
			i__2 = kk + *n - j;
			for (k = kk + 1; k <= kk+*n-j; ++k) {
			    ix += *incx;
			    i__3 = ix;
			    i__4 = ix;
			    i__5 = k;
			    z__2.r = temp.r * AP(k).r - temp.i * AP(k)
				    .i, z__2.i = temp.r * AP(k).i + temp.i 
				    * AP(k).r;
			    z__1.r = X(ix).r - z__2.r, z__1.i = X(ix).i - 
				    z__2.i;
			    X(ix).r = z__1.r, X(ix).i = z__1.i;
/* L70: */
			}
		    }
		    jx += *incx;
		    kk += *n - j + 1;
/* L80: */
		}
	    }
	}
    } else {

/*        Form  x := inv( A' )*x  or  x := inv( conjg( A' ) )*x. */

	if (lsame_(uplo, "U")) {
	    kk = 1;
	    if (*incx == 1) {
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
		    i__2 = j;
		    temp.r = X(j).r, temp.i = X(j).i;
		    k = kk;
		    if (noconj) {
			i__2 = j - 1;
			for (i = 1; i <= j-1; ++i) {
			    i__3 = k;
			    i__4 = i;
			    z__2.r = AP(k).r * X(i).r - AP(k).i * X(
				    i).i, z__2.i = AP(k).r * X(i).i 
				    + AP(k).i * X(i).r;
			    z__1.r = temp.r - z__2.r, z__1.i = temp.i - 
				    z__2.i;
			    temp.r = z__1.r, temp.i = z__1.i;
			    ++k;
/* L90: */
			}
			if (nounit) {
			    z_div(&z__1, &temp, &AP(kk + j - 1));
			    temp.r = z__1.r, temp.i = z__1.i;
			}
		    } else {
			i__2 = j - 1;
			for (i = 1; i <= j-1; ++i) {
			    d_cnjg(&z__3, &AP(k));
			    i__3 = i;
			    z__2.r = z__3.r * X(i).r - z__3.i * X(i).i, 
				    z__2.i = z__3.r * X(i).i + z__3.i * X(
				    i).r;
			    z__1.r = temp.r - z__2.r, z__1.i = temp.i - 
				    z__2.i;
			    temp.r = z__1.r, temp.i = z__1.i;
			    ++k;
/* L100: */
			}
			if (nounit) {
			    d_cnjg(&z__2, &AP(kk + j - 1));
			    z_div(&z__1, &temp, &z__2);
			    temp.r = z__1.r, temp.i = z__1.i;
			}
		    }
		    i__2 = j;
		    X(j).r = temp.r, X(j).i = temp.i;
		    kk += j;
/* L110: */
		}
	    } else {
		jx = kx;
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
		    i__2 = jx;
		    temp.r = X(jx).r, temp.i = X(jx).i;
		    ix = kx;
		    if (noconj) {
			i__2 = kk + j - 2;
			for (k = kk; k <= kk+j-2; ++k) {
			    i__3 = k;
			    i__4 = ix;
			    z__2.r = AP(k).r * X(ix).r - AP(k).i * X(
				    ix).i, z__2.i = AP(k).r * X(ix).i 
				    + AP(k).i * X(ix).r;
			    z__1.r = temp.r - z__2.r, z__1.i = temp.i - 
				    z__2.i;
			    temp.r = z__1.r, temp.i = z__1.i;
			    ix += *incx;
/* L120: */
			}
			if (nounit) {
			    z_div(&z__1, &temp, &AP(kk + j - 1));
			    temp.r = z__1.r, temp.i = z__1.i;
			}
		    } else {
			i__2 = kk + j - 2;
			for (k = kk; k <= kk+j-2; ++k) {
			    d_cnjg(&z__3, &AP(k));
			    i__3 = ix;
			    z__2.r = z__3.r * X(ix).r - z__3.i * X(ix).i, 
				    z__2.i = z__3.r * X(ix).i + z__3.i * X(
				    ix).r;
			    z__1.r = temp.r - z__2.r, z__1.i = temp.i - 
				    z__2.i;
			    temp.r = z__1.r, temp.i = z__1.i;
			    ix += *incx;
/* L130: */
			}
			if (nounit) {
			    d_cnjg(&z__2, &AP(kk + j - 1));
			    z_div(&z__1, &temp, &z__2);
			    temp.r = z__1.r, temp.i = z__1.i;
			}
		    }
		    i__2 = jx;
		    X(jx).r = temp.r, X(jx).i = temp.i;
		    jx += *incx;
		    kk += j;
/* L140: */
		}
	    }
	} else {
	    kk = *n * (*n + 1) / 2;
	    if (*incx == 1) {
		for (j = *n; j >= 1; --j) {
		    i__1 = j;
		    temp.r = X(j).r, temp.i = X(j).i;
		    k = kk;
		    if (noconj) {
			i__1 = j + 1;
			for (i = *n; i >= j+1; --i) {
			    i__2 = k;
			    i__3 = i;
			    z__2.r = AP(k).r * X(i).r - AP(k).i * X(
				    i).i, z__2.i = AP(k).r * X(i).i 
				    + AP(k).i * X(i).r;
			    z__1.r = temp.r - z__2.r, z__1.i = temp.i - 
				    z__2.i;
			    temp.r = z__1.r, temp.i = z__1.i;
			    --k;
/* L150: */
			}
			if (nounit) {
			    z_div(&z__1, &temp, &AP(kk - *n + j));
			    temp.r = z__1.r, temp.i = z__1.i;
			}
		    } else {
			i__1 = j + 1;
			for (i = *n; i >= j+1; --i) {
			    d_cnjg(&z__3, &AP(k));
			    i__2 = i;
			    z__2.r = z__3.r * X(i).r - z__3.i * X(i).i, 
				    z__2.i = z__3.r * X(i).i + z__3.i * X(
				    i).r;
			    z__1.r = temp.r - z__2.r, z__1.i = temp.i - 
				    z__2.i;
			    temp.r = z__1.r, temp.i = z__1.i;
			    --k;
/* L160: */
			}
			if (nounit) {
			    d_cnjg(&z__2, &AP(kk - *n + j));
			    z_div(&z__1, &temp, &z__2);
			    temp.r = z__1.r, temp.i = z__1.i;
			}
		    }
		    i__1 = j;
		    X(j).r = temp.r, X(j).i = temp.i;
		    kk -= *n - j + 1;
/* L170: */
		}
	    } else {
		kx += (*n - 1) * *incx;
		jx = kx;
		for (j = *n; j >= 1; --j) {
		    i__1 = jx;
		    temp.r = X(jx).r, temp.i = X(jx).i;
		    ix = kx;
		    if (noconj) {
			i__1 = kk - (*n - (j + 1));
			for (k = kk; k >= kk-(*n-(j+1)); --k) {
			    i__2 = k;
			    i__3 = ix;
			    z__2.r = AP(k).r * X(ix).r - AP(k).i * X(
				    ix).i, z__2.i = AP(k).r * X(ix).i 
				    + AP(k).i * X(ix).r;
			    z__1.r = temp.r - z__2.r, z__1.i = temp.i - 
				    z__2.i;
			    temp.r = z__1.r, temp.i = z__1.i;
			    ix -= *incx;
/* L180: */
			}
			if (nounit) {
			    z_div(&z__1, &temp, &AP(kk - *n + j));
			    temp.r = z__1.r, temp.i = z__1.i;
			}
		    } else {
			i__1 = kk - (*n - (j + 1));
			for (k = kk; k >= kk-(*n-(j+1)); --k) {
			    d_cnjg(&z__3, &AP(k));
			    i__2 = ix;
			    z__2.r = z__3.r * X(ix).r - z__3.i * X(ix).i, 
				    z__2.i = z__3.r * X(ix).i + z__3.i * X(
				    ix).r;
			    z__1.r = temp.r - z__2.r, z__1.i = temp.i - 
				    z__2.i;
			    temp.r = z__1.r, temp.i = z__1.i;
			    ix -= *incx;
/* L190: */
			}
			if (nounit) {
			    d_cnjg(&z__2, &AP(kk - *n + j));
			    z_div(&z__1, &temp, &z__2);
			    temp.r = z__1.r, temp.i = z__1.i;
			}
		    }
		    i__1 = jx;
		    X(jx).r = temp.r, X(jx).i = temp.i;
		    jx -= *incx;
		    kk -= *n - j + 1;
/* L200: */
		}
	    }
	}
    }

    return 0;

/*     End of ZTPSV . */

} /* ztpsv_ */

