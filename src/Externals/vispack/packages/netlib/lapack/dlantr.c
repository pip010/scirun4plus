#include "f2c.h"

doublereal dlantr_(char *norm, char *uplo, char *diag, integer *m, integer *n,
	 doublereal *a, integer *lda, doublereal *work)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    DLANTR  returns the value of the one norm,  or the Frobenius norm, or 
  
    the  infinity norm,  or the  element of  largest absolute value  of a 
  
    trapezoidal or triangular matrix A.   

    Description   
    ===========   

    DLANTR returns the value   

       DLANTR = ( max(abs(A(i,j))), NORM = 'M' or 'm'   
                (   
                ( norm1(A),         NORM = '1', 'O' or 'o'   
                (   
                ( normI(A),         NORM = 'I' or 'i'   
                (   
                ( normF(A),         NORM = 'F', 'f', 'E' or 'e'   

    where  norm1  denotes the  one norm of a matrix (maximum column sum), 
  
    normI  denotes the  infinity norm  of a matrix  (maximum row sum) and 
  
    normF  denotes the  Frobenius norm of a matrix (square root of sum of 
  
    squares).  Note that  max(abs(A(i,j)))  is not a  matrix norm.   

    Arguments   
    =========   

    NORM    (input) CHARACTER*1   
            Specifies the value to be returned in DLANTR as described   
            above.   

    UPLO    (input) CHARACTER*1   
            Specifies whether the matrix A is upper or lower trapezoidal. 
  
            = 'U':  Upper trapezoidal   
            = 'L':  Lower trapezoidal   
            Note that A is triangular instead of trapezoidal if M = N.   

    DIAG    (input) CHARACTER*1   
            Specifies whether or not the matrix A has unit diagonal.   
            = 'N':  Non-unit diagonal   
            = 'U':  Unit diagonal   

    M       (input) INTEGER   
            The number of rows of the matrix A.  M >= 0, and if   
            UPLO = 'U', M <= N.  When M = 0, DLANTR is set to zero.   

    N       (input) INTEGER   
            The number of columns of the matrix A.  N >= 0, and if   
            UPLO = 'L', N <= M.  When N = 0, DLANTR is set to zero.   

    A       (input) DOUBLE PRECISION array, dimension (LDA,N)   
            The trapezoidal matrix A (A is triangular if M = N).   
            If UPLO = 'U', the leading m by n upper trapezoidal part of   
            the array A contains the upper trapezoidal matrix, and the   
            strictly lower triangular part of A is not referenced.   
            If UPLO = 'L', the leading m by n lower trapezoidal part of   
            the array A contains the lower trapezoidal matrix, and the   
            strictly upper triangular part of A is not referenced.  Note 
  
            that when DIAG = 'U', the diagonal elements of A are not   
            referenced and are assumed to be one.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(M,1).   

    WORK    (workspace) DOUBLE PRECISION array, dimension (LWORK),   
            where LWORK >= M when NORM = 'I'; otherwise, WORK is not   
            referenced.   

   ===================================================================== 
  


    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static integer c__1 = 1;
    
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;
    doublereal ret_val, d__1, d__2, d__3;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    static integer i, j;
    static doublereal scale;
    static logical udiag;
    extern logical lsame_(char *, char *);
    static doublereal value;
    extern /* Subroutine */ int dlassq_(integer *, doublereal *, integer *, 
	    doublereal *, doublereal *);
    static doublereal sum;



#define WORK(I) work[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    if (min(*m,*n) == 0) {
	value = 0.;
    } else if (lsame_(norm, "M")) {

/*        Find max(abs(A(i,j))). */

	if (lsame_(diag, "U")) {
	    value = 1.;
	    if (lsame_(uplo, "U")) {
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
/* Computing MIN */
		    i__3 = *m, i__4 = j - 1;
		    i__2 = min(i__3,i__4);
		    for (i = 1; i <= min(*m,j-1); ++i) {
/* Computing MAX */
			d__2 = value, d__3 = (d__1 = A(i,j), abs(
				d__1));
			value = max(d__2,d__3);
/* L10: */
		    }
/* L20: */
		}
	    } else {
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
		    i__2 = *m;
		    for (i = j + 1; i <= *m; ++i) {
/* Computing MAX */
			d__2 = value, d__3 = (d__1 = A(i,j), abs(
				d__1));
			value = max(d__2,d__3);
/* L30: */
		    }
/* L40: */
		}
	    }
	} else {
	    value = 0.;
	    if (lsame_(uplo, "U")) {
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
		    i__2 = min(*m,j);
		    for (i = 1; i <= min(*m,j); ++i) {
/* Computing MAX */
			d__2 = value, d__3 = (d__1 = A(i,j), abs(
				d__1));
			value = max(d__2,d__3);
/* L50: */
		    }
/* L60: */
		}
	    } else {
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
		    i__2 = *m;
		    for (i = j; i <= *m; ++i) {
/* Computing MAX */
			d__2 = value, d__3 = (d__1 = A(i,j), abs(
				d__1));
			value = max(d__2,d__3);
/* L70: */
		    }
/* L80: */
		}
	    }
	}
    } else if (lsame_(norm, "O") || *(unsigned char *)norm == '1') {

/*        Find norm1(A). */

	value = 0.;
	udiag = lsame_(diag, "U");
	if (lsame_(uplo, "U")) {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		if (udiag && j <= *m) {
		    sum = 1.;
		    i__2 = j - 1;
		    for (i = 1; i <= j-1; ++i) {
			sum += (d__1 = A(i,j), abs(d__1));
/* L90: */
		    }
		} else {
		    sum = 0.;
		    i__2 = min(*m,j);
		    for (i = 1; i <= min(*m,j); ++i) {
			sum += (d__1 = A(i,j), abs(d__1));
/* L100: */
		    }
		}
		value = max(value,sum);
/* L110: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		if (udiag) {
		    sum = 1.;
		    i__2 = *m;
		    for (i = j + 1; i <= *m; ++i) {
			sum += (d__1 = A(i,j), abs(d__1));
/* L120: */
		    }
		} else {
		    sum = 0.;
		    i__2 = *m;
		    for (i = j; i <= *m; ++i) {
			sum += (d__1 = A(i,j), abs(d__1));
/* L130: */
		    }
		}
		value = max(value,sum);
/* L140: */
	    }
	}
    } else if (lsame_(norm, "I")) {

/*        Find normI(A). */

	if (lsame_(uplo, "U")) {
	    if (lsame_(diag, "U")) {
		i__1 = *m;
		for (i = 1; i <= *m; ++i) {
		    WORK(i) = 1.;
/* L150: */
		}
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
/* Computing MIN */
		    i__3 = *m, i__4 = j - 1;
		    i__2 = min(i__3,i__4);
		    for (i = 1; i <= min(*m,j-1); ++i) {
			WORK(i) += (d__1 = A(i,j), abs(d__1));
/* L160: */
		    }
/* L170: */
		}
	    } else {
		i__1 = *m;
		for (i = 1; i <= *m; ++i) {
		    WORK(i) = 0.;
/* L180: */
		}
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
		    i__2 = min(*m,j);
		    for (i = 1; i <= min(*m,j); ++i) {
			WORK(i) += (d__1 = A(i,j), abs(d__1));
/* L190: */
		    }
/* L200: */
		}
	    }
	} else {
	    if (lsame_(diag, "U")) {
		i__1 = *n;
		for (i = 1; i <= *n; ++i) {
		    WORK(i) = 1.;
/* L210: */
		}
		i__1 = *m;
		for (i = *n + 1; i <= *m; ++i) {
		    WORK(i) = 0.;
/* L220: */
		}
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
		    i__2 = *m;
		    for (i = j + 1; i <= *m; ++i) {
			WORK(i) += (d__1 = A(i,j), abs(d__1));
/* L230: */
		    }
/* L240: */
		}
	    } else {
		i__1 = *m;
		for (i = 1; i <= *m; ++i) {
		    WORK(i) = 0.;
/* L250: */
		}
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
		    i__2 = *m;
		    for (i = j; i <= *m; ++i) {
			WORK(i) += (d__1 = A(i,j), abs(d__1));
/* L260: */
		    }
/* L270: */
		}
	    }
	}
	value = 0.;
	i__1 = *m;
	for (i = 1; i <= *m; ++i) {
/* Computing MAX */
	    d__1 = value, d__2 = WORK(i);
	    value = max(d__1,d__2);
/* L280: */
	}
    } else if (lsame_(norm, "F") || lsame_(norm, "E")) {

/*        Find normF(A). */

	if (lsame_(uplo, "U")) {
	    if (lsame_(diag, "U")) {
		scale = 1.;
		sum = (doublereal) min(*m,*n);
		i__1 = *n;
		for (j = 2; j <= *n; ++j) {
/* Computing MIN */
		    i__3 = *m, i__4 = j - 1;
		    i__2 = min(i__3,i__4);
		    dlassq_(&i__2, &A(1,j), &c__1, &scale, &sum);
/* L290: */
		}
	    } else {
		scale = 0.;
		sum = 1.;
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
		    i__2 = min(*m,j);
		    dlassq_(&i__2, &A(1,j), &c__1, &scale, &sum);
/* L300: */
		}
	    }
	} else {
	    if (lsame_(diag, "U")) {
		scale = 1.;
		sum = (doublereal) min(*m,*n);
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
		    i__2 = *m - j;
/* Computing MIN */
		    i__3 = *m, i__4 = j + 1;
		    dlassq_(&i__2, &A(min(*m,j+1),j), &c__1, &
			    scale, &sum);
/* L310: */
		}
	    } else {
		scale = 0.;
		sum = 1.;
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
		    i__2 = *m - j + 1;
		    dlassq_(&i__2, &A(j,j), &c__1, &scale, &sum);
/* L320: */
		}
	    }
	}
	value = scale * sqrt(sum);
    }

    ret_val = value;
    return ret_val;

/*     End of DLANTR */

} /* dlantr_ */

