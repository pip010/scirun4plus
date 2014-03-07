#include "f2c.h"

/* Subroutine */ int zlasr_(char *side, char *pivot, char *direct, integer *m,
	 integer *n, doublereal *c, doublereal *s, doublecomplex *a, integer *
	lda)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    ZLASR   performs the transformation   

       A := P*A,   when SIDE = 'L' or 'l'  (  Left-hand side )   

       A := A*P',  when SIDE = 'R' or 'r'  ( Right-hand side )   

    where A is an m by n complex matrix and P is an orthogonal matrix,   
    consisting of a sequence of plane rotations determined by the   
    parameters PIVOT and DIRECT as follows ( z = m when SIDE = 'L' or 'l' 
  
    and z = n when SIDE = 'R' or 'r' ):   

    When  DIRECT = 'F' or 'f'  ( Forward sequence ) then   

       P = P( z - 1 )*...*P( 2 )*P( 1 ),   

    and when DIRECT = 'B' or 'b'  ( Backward sequence ) then   

       P = P( 1 )*P( 2 )*...*P( z - 1 ),   

    where  P( k ) is a plane rotation matrix for the following planes:   

       when  PIVOT = 'V' or 'v'  ( Variable pivot ),   
          the plane ( k, k + 1 )   

       when  PIVOT = 'T' or 't'  ( Top pivot ),   
          the plane ( 1, k + 1 )   

       when  PIVOT = 'B' or 'b'  ( Bottom pivot ),   
          the plane ( k, z )   

    c( k ) and s( k )  must contain the  cosine and sine that define the 
  
    matrix  P( k ).  The two by two plane rotation part of the matrix   
    P( k ), R( k ), is assumed to be of the form   

       R( k ) = (  c( k )  s( k ) ).   
                ( -s( k )  c( k ) )   

    Arguments   
    =========   

    SIDE    (input) CHARACTER*1   
            Specifies whether the plane rotation matrix P is applied to   
            A on the left or the right.   
            = 'L':  Left, compute A := P*A   
            = 'R':  Right, compute A:= A*P'   

    DIRECT  (input) CHARACTER*1   
            Specifies whether P is a forward or backward sequence of   
            plane rotations.   
            = 'F':  Forward, P = P( z - 1 )*...*P( 2 )*P( 1 )   
            = 'B':  Backward, P = P( 1 )*P( 2 )*...*P( z - 1 )   

    PIVOT   (input) CHARACTER*1   
            Specifies the plane for which P(k) is a plane rotation   
            matrix.   
            = 'V':  Variable pivot, the plane (k,k+1)   
            = 'T':  Top pivot, the plane (1,k+1)   
            = 'B':  Bottom pivot, the plane (k,z)   

    M       (input) INTEGER   
            The number of rows of the matrix A.  If m <= 1, an immediate 
  
            return is effected.   

    N       (input) INTEGER   
            The number of columns of the matrix A.  If n <= 1, an   
            immediate return is effected.   

    C, S    (input) DOUBLE PRECISION arrays, dimension   
                    (M-1) if SIDE = 'L'   
                    (N-1) if SIDE = 'R'   
            c(k) and s(k) contain the cosine and sine that define the   
            matrix P(k).  The two by two plane rotation part of the   
            matrix P(k), R(k), is assumed to be of the form   
            R( k ) = (  c( k )  s( k ) ).   
                     ( -s( k )  c( k ) )   

    A       (input/output) COMPLEX*16 array, dimension (LDA,N)   
            The m by n matrix A.  On exit, A is overwritten by P*A if   
            SIDE = 'R' or by A*P' if SIDE = 'L'.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(1,M).   

    ===================================================================== 
  


       Test the input parameters   

    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;
    doublecomplex z__1, z__2, z__3;
    /* Local variables */
    static integer info;
    static doublecomplex temp;
    static integer i, j;
    extern logical lsame_(char *, char *);
    static doublereal ctemp, stemp;
    extern /* Subroutine */ int xerbla_(char *, integer *);


#define C(I) c[(I)-1]
#define S(I) s[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    info = 0;
    if (! (lsame_(side, "L") || lsame_(side, "R"))) {
	info = 1;
    } else if (! (lsame_(pivot, "V") || lsame_(pivot, "T") || 
	    lsame_(pivot, "B"))) {
	info = 2;
    } else if (! (lsame_(direct, "F") || lsame_(direct, "B")))
	     {
	info = 3;
    } else if (*m < 0) {
	info = 4;
    } else if (*n < 0) {
	info = 5;
    } else if (*lda < max(1,*m)) {
	info = 9;
    }
    if (info != 0) {
	xerbla_("ZLASR ", &info);
	return 0;
    }

/*     Quick return if possible */

    if (*m == 0 || *n == 0) {
	return 0;
    }
    if (lsame_(side, "L")) {

/*        Form  P * A */

	if (lsame_(pivot, "V")) {
	    if (lsame_(direct, "F")) {
		i__1 = *m - 1;
		for (j = 1; j <= *m-1; ++j) {
		    ctemp = C(j);
		    stemp = S(j);
		    if (ctemp != 1. || stemp != 0.) {
			i__2 = *n;
			for (i = 1; i <= *n; ++i) {
			    i__3 = j + 1 + i * a_dim1;
			    temp.r = A(j+1,i).r, temp.i = A(j+1,i).i;
			    i__3 = j + 1 + i * a_dim1;
			    z__2.r = ctemp * temp.r, z__2.i = ctemp * temp.i;
			    i__4 = j + i * a_dim1;
			    z__3.r = stemp * A(j,i).r, z__3.i = stemp * A(j,i).i;
			    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - 
				    z__3.i;
			    A(j+1,i).r = z__1.r, A(j+1,i).i = z__1.i;
			    i__3 = j + i * a_dim1;
			    z__2.r = stemp * temp.r, z__2.i = stemp * temp.i;
			    i__4 = j + i * a_dim1;
			    z__3.r = ctemp * A(j,i).r, z__3.i = ctemp * A(j,i).i;
			    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + 
				    z__3.i;
			    A(j,i).r = z__1.r, A(j,i).i = z__1.i;
/* L10: */
			}
		    }
/* L20: */
		}
	    } else if (lsame_(direct, "B")) {
		for (j = *m - 1; j >= 1; --j) {
		    ctemp = C(j);
		    stemp = S(j);
		    if (ctemp != 1. || stemp != 0.) {
			i__1 = *n;
			for (i = 1; i <= *n; ++i) {
			    i__2 = j + 1 + i * a_dim1;
			    temp.r = A(j+1,i).r, temp.i = A(j+1,i).i;
			    i__2 = j + 1 + i * a_dim1;
			    z__2.r = ctemp * temp.r, z__2.i = ctemp * temp.i;
			    i__3 = j + i * a_dim1;
			    z__3.r = stemp * A(j,i).r, z__3.i = stemp * A(j,i).i;
			    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - 
				    z__3.i;
			    A(j+1,i).r = z__1.r, A(j+1,i).i = z__1.i;
			    i__2 = j + i * a_dim1;
			    z__2.r = stemp * temp.r, z__2.i = stemp * temp.i;
			    i__3 = j + i * a_dim1;
			    z__3.r = ctemp * A(j,i).r, z__3.i = ctemp * A(j,i).i;
			    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + 
				    z__3.i;
			    A(j,i).r = z__1.r, A(j,i).i = z__1.i;
/* L30: */
			}
		    }
/* L40: */
		}
	    }
	} else if (lsame_(pivot, "T")) {
	    if (lsame_(direct, "F")) {
		i__1 = *m;
		for (j = 2; j <= *m; ++j) {
		    ctemp = C(j - 1);
		    stemp = S(j - 1);
		    if (ctemp != 1. || stemp != 0.) {
			i__2 = *n;
			for (i = 1; i <= *n; ++i) {
			    i__3 = j + i * a_dim1;
			    temp.r = A(j,i).r, temp.i = A(j,i).i;
			    i__3 = j + i * a_dim1;
			    z__2.r = ctemp * temp.r, z__2.i = ctemp * temp.i;
			    i__4 = i * a_dim1 + 1;
			    z__3.r = stemp * A(1,i).r, z__3.i = stemp * A(1,i).i;
			    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - 
				    z__3.i;
			    A(j,i).r = z__1.r, A(j,i).i = z__1.i;
			    i__3 = i * a_dim1 + 1;
			    z__2.r = stemp * temp.r, z__2.i = stemp * temp.i;
			    i__4 = i * a_dim1 + 1;
			    z__3.r = ctemp * A(1,i).r, z__3.i = ctemp * A(1,i).i;
			    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + 
				    z__3.i;
			    A(1,i).r = z__1.r, A(1,i).i = z__1.i;
/* L50: */
			}
		    }
/* L60: */
		}
	    } else if (lsame_(direct, "B")) {
		for (j = *m; j >= 2; --j) {
		    ctemp = C(j - 1);
		    stemp = S(j - 1);
		    if (ctemp != 1. || stemp != 0.) {
			i__1 = *n;
			for (i = 1; i <= *n; ++i) {
			    i__2 = j + i * a_dim1;
			    temp.r = A(j,i).r, temp.i = A(j,i).i;
			    i__2 = j + i * a_dim1;
			    z__2.r = ctemp * temp.r, z__2.i = ctemp * temp.i;
			    i__3 = i * a_dim1 + 1;
			    z__3.r = stemp * A(1,i).r, z__3.i = stemp * A(1,i).i;
			    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - 
				    z__3.i;
			    A(j,i).r = z__1.r, A(j,i).i = z__1.i;
			    i__2 = i * a_dim1 + 1;
			    z__2.r = stemp * temp.r, z__2.i = stemp * temp.i;
			    i__3 = i * a_dim1 + 1;
			    z__3.r = ctemp * A(1,i).r, z__3.i = ctemp * A(1,i).i;
			    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + 
				    z__3.i;
			    A(1,i).r = z__1.r, A(1,i).i = z__1.i;
/* L70: */
			}
		    }
/* L80: */
		}
	    }
	} else if (lsame_(pivot, "B")) {
	    if (lsame_(direct, "F")) {
		i__1 = *m - 1;
		for (j = 1; j <= *m-1; ++j) {
		    ctemp = C(j);
		    stemp = S(j);
		    if (ctemp != 1. || stemp != 0.) {
			i__2 = *n;
			for (i = 1; i <= *n; ++i) {
			    i__3 = j + i * a_dim1;
			    temp.r = A(j,i).r, temp.i = A(j,i).i;
			    i__3 = j + i * a_dim1;
			    i__4 = *m + i * a_dim1;
			    z__2.r = stemp * A(*m,i).r, z__2.i = stemp * A(*m,i).i;
			    z__3.r = ctemp * temp.r, z__3.i = ctemp * temp.i;
			    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + 
				    z__3.i;
			    A(j,i).r = z__1.r, A(j,i).i = z__1.i;
			    i__3 = *m + i * a_dim1;
			    i__4 = *m + i * a_dim1;
			    z__2.r = ctemp * A(*m,i).r, z__2.i = ctemp * A(*m,i).i;
			    z__3.r = stemp * temp.r, z__3.i = stemp * temp.i;
			    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - 
				    z__3.i;
			    A(*m,i).r = z__1.r, A(*m,i).i = z__1.i;
/* L90: */
			}
		    }
/* L100: */
		}
	    } else if (lsame_(direct, "B")) {
		for (j = *m - 1; j >= 1; --j) {
		    ctemp = C(j);
		    stemp = S(j);
		    if (ctemp != 1. || stemp != 0.) {
			i__1 = *n;
			for (i = 1; i <= *n; ++i) {
			    i__2 = j + i * a_dim1;
			    temp.r = A(j,i).r, temp.i = A(j,i).i;
			    i__2 = j + i * a_dim1;
			    i__3 = *m + i * a_dim1;
			    z__2.r = stemp * A(*m,i).r, z__2.i = stemp * A(*m,i).i;
			    z__3.r = ctemp * temp.r, z__3.i = ctemp * temp.i;
			    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + 
				    z__3.i;
			    A(j,i).r = z__1.r, A(j,i).i = z__1.i;
			    i__2 = *m + i * a_dim1;
			    i__3 = *m + i * a_dim1;
			    z__2.r = ctemp * A(*m,i).r, z__2.i = ctemp * A(*m,i).i;
			    z__3.r = stemp * temp.r, z__3.i = stemp * temp.i;
			    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - 
				    z__3.i;
			    A(*m,i).r = z__1.r, A(*m,i).i = z__1.i;
/* L110: */
			}
		    }
/* L120: */
		}
	    }
	}
    } else if (lsame_(side, "R")) {

/*        Form A * P' */

	if (lsame_(pivot, "V")) {
	    if (lsame_(direct, "F")) {
		i__1 = *n - 1;
		for (j = 1; j <= *n-1; ++j) {
		    ctemp = C(j);
		    stemp = S(j);
		    if (ctemp != 1. || stemp != 0.) {
			i__2 = *m;
			for (i = 1; i <= *m; ++i) {
			    i__3 = i + (j + 1) * a_dim1;
			    temp.r = A(i,j+1).r, temp.i = A(i,j+1).i;
			    i__3 = i + (j + 1) * a_dim1;
			    z__2.r = ctemp * temp.r, z__2.i = ctemp * temp.i;
			    i__4 = i + j * a_dim1;
			    z__3.r = stemp * A(i,j).r, z__3.i = stemp * A(i,j).i;
			    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - 
				    z__3.i;
			    A(i,j+1).r = z__1.r, A(i,j+1).i = z__1.i;
			    i__3 = i + j * a_dim1;
			    z__2.r = stemp * temp.r, z__2.i = stemp * temp.i;
			    i__4 = i + j * a_dim1;
			    z__3.r = ctemp * A(i,j).r, z__3.i = ctemp * A(i,j).i;
			    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + 
				    z__3.i;
			    A(i,j).r = z__1.r, A(i,j).i = z__1.i;
/* L130: */
			}
		    }
/* L140: */
		}
	    } else if (lsame_(direct, "B")) {
		for (j = *n - 1; j >= 1; --j) {
		    ctemp = C(j);
		    stemp = S(j);
		    if (ctemp != 1. || stemp != 0.) {
			i__1 = *m;
			for (i = 1; i <= *m; ++i) {
			    i__2 = i + (j + 1) * a_dim1;
			    temp.r = A(i,j+1).r, temp.i = A(i,j+1).i;
			    i__2 = i + (j + 1) * a_dim1;
			    z__2.r = ctemp * temp.r, z__2.i = ctemp * temp.i;
			    i__3 = i + j * a_dim1;
			    z__3.r = stemp * A(i,j).r, z__3.i = stemp * A(i,j).i;
			    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - 
				    z__3.i;
			    A(i,j+1).r = z__1.r, A(i,j+1).i = z__1.i;
			    i__2 = i + j * a_dim1;
			    z__2.r = stemp * temp.r, z__2.i = stemp * temp.i;
			    i__3 = i + j * a_dim1;
			    z__3.r = ctemp * A(i,j).r, z__3.i = ctemp * A(i,j).i;
			    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + 
				    z__3.i;
			    A(i,j).r = z__1.r, A(i,j).i = z__1.i;
/* L150: */
			}
		    }
/* L160: */
		}
	    }
	} else if (lsame_(pivot, "T")) {
	    if (lsame_(direct, "F")) {
		i__1 = *n;
		for (j = 2; j <= *n; ++j) {
		    ctemp = C(j - 1);
		    stemp = S(j - 1);
		    if (ctemp != 1. || stemp != 0.) {
			i__2 = *m;
			for (i = 1; i <= *m; ++i) {
			    i__3 = i + j * a_dim1;
			    temp.r = A(i,j).r, temp.i = A(i,j).i;
			    i__3 = i + j * a_dim1;
			    z__2.r = ctemp * temp.r, z__2.i = ctemp * temp.i;
			    i__4 = i + a_dim1;
			    z__3.r = stemp * A(i,1).r, z__3.i = stemp * A(i,1).i;
			    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - 
				    z__3.i;
			    A(i,j).r = z__1.r, A(i,j).i = z__1.i;
			    i__3 = i + a_dim1;
			    z__2.r = stemp * temp.r, z__2.i = stemp * temp.i;
			    i__4 = i + a_dim1;
			    z__3.r = ctemp * A(i,1).r, z__3.i = ctemp * A(i,1).i;
			    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + 
				    z__3.i;
			    A(i,1).r = z__1.r, A(i,1).i = z__1.i;
/* L170: */
			}
		    }
/* L180: */
		}
	    } else if (lsame_(direct, "B")) {
		for (j = *n; j >= 2; --j) {
		    ctemp = C(j - 1);
		    stemp = S(j - 1);
		    if (ctemp != 1. || stemp != 0.) {
			i__1 = *m;
			for (i = 1; i <= *m; ++i) {
			    i__2 = i + j * a_dim1;
			    temp.r = A(i,j).r, temp.i = A(i,j).i;
			    i__2 = i + j * a_dim1;
			    z__2.r = ctemp * temp.r, z__2.i = ctemp * temp.i;
			    i__3 = i + a_dim1;
			    z__3.r = stemp * A(i,1).r, z__3.i = stemp * A(i,1).i;
			    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - 
				    z__3.i;
			    A(i,j).r = z__1.r, A(i,j).i = z__1.i;
			    i__2 = i + a_dim1;
			    z__2.r = stemp * temp.r, z__2.i = stemp * temp.i;
			    i__3 = i + a_dim1;
			    z__3.r = ctemp * A(i,1).r, z__3.i = ctemp * A(i,1).i;
			    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + 
				    z__3.i;
			    A(i,1).r = z__1.r, A(i,1).i = z__1.i;
/* L190: */
			}
		    }
/* L200: */
		}
	    }
	} else if (lsame_(pivot, "B")) {
	    if (lsame_(direct, "F")) {
		i__1 = *n - 1;
		for (j = 1; j <= *n-1; ++j) {
		    ctemp = C(j);
		    stemp = S(j);
		    if (ctemp != 1. || stemp != 0.) {
			i__2 = *m;
			for (i = 1; i <= *m; ++i) {
			    i__3 = i + j * a_dim1;
			    temp.r = A(i,j).r, temp.i = A(i,j).i;
			    i__3 = i + j * a_dim1;
			    i__4 = i + *n * a_dim1;
			    z__2.r = stemp * A(i,*n).r, z__2.i = stemp * A(i,*n).i;
			    z__3.r = ctemp * temp.r, z__3.i = ctemp * temp.i;
			    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + 
				    z__3.i;
			    A(i,j).r = z__1.r, A(i,j).i = z__1.i;
			    i__3 = i + *n * a_dim1;
			    i__4 = i + *n * a_dim1;
			    z__2.r = ctemp * A(i,*n).r, z__2.i = ctemp * A(i,*n).i;
			    z__3.r = stemp * temp.r, z__3.i = stemp * temp.i;
			    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - 
				    z__3.i;
			    A(i,*n).r = z__1.r, A(i,*n).i = z__1.i;
/* L210: */
			}
		    }
/* L220: */
		}
	    } else if (lsame_(direct, "B")) {
		for (j = *n - 1; j >= 1; --j) {
		    ctemp = C(j);
		    stemp = S(j);
		    if (ctemp != 1. || stemp != 0.) {
			i__1 = *m;
			for (i = 1; i <= *m; ++i) {
			    i__2 = i + j * a_dim1;
			    temp.r = A(i,j).r, temp.i = A(i,j).i;
			    i__2 = i + j * a_dim1;
			    i__3 = i + *n * a_dim1;
			    z__2.r = stemp * A(i,*n).r, z__2.i = stemp * A(i,*n).i;
			    z__3.r = ctemp * temp.r, z__3.i = ctemp * temp.i;
			    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + 
				    z__3.i;
			    A(i,j).r = z__1.r, A(i,j).i = z__1.i;
			    i__2 = i + *n * a_dim1;
			    i__3 = i + *n * a_dim1;
			    z__2.r = ctemp * A(i,*n).r, z__2.i = ctemp * A(i,*n).i;
			    z__3.r = stemp * temp.r, z__3.i = stemp * temp.i;
			    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - 
				    z__3.i;
			    A(i,*n).r = z__1.r, A(i,*n).i = z__1.i;
/* L230: */
			}
		    }
/* L240: */
		}
	    }
	}
    }

    return 0;

/*     End of ZLASR */

} /* zlasr_ */

