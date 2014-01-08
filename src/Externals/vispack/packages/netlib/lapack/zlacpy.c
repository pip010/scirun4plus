#include "f2c.h"

/* Subroutine */ int zlacpy_(char *uplo, integer *m, integer *n, 
	doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   


    Purpose   
    =======   

    ZLACPY copies all or part of a two-dimensional matrix A to another   
    matrix B.   

    Arguments   
    =========   

    UPLO    (input) CHARACTER*1   
            Specifies the part of the matrix A to be copied to B.   
            = 'U':      Upper triangular part   
            = 'L':      Lower triangular part   
            Otherwise:  All of the matrix A   

    M       (input) INTEGER   
            The number of rows of the matrix A.  M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrix A.  N >= 0.   

    A       (input) COMPLEX*16 array, dimension (LDA,N)   
            The m by n matrix A.  If UPLO = 'U', only the upper trapezium 
  
            is accessed; if UPLO = 'L', only the lower trapezium is   
            accessed.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(1,M).   

    B       (output) COMPLEX*16 array, dimension (LDB,N)   
            On exit, B = A in the locations specified by UPLO.   

    LDB     (input) INTEGER   
            The leading dimension of the array B.  LDB >= max(1,M).   

    ===================================================================== 
  


    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3, i__4;
    /* Local variables */
    static integer i, j;
    extern logical lsame_(char *, char *);



#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]
#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]

    if (lsame_(uplo, "U")) {
	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
	    i__2 = min(j,*m);
	    for (i = 1; i <= min(j,*m); ++i) {
		i__3 = i + j * b_dim1;
		i__4 = i + j * a_dim1;
		B(i,j).r = A(i,j).r, B(i,j).i = A(i,j).i;
/* L10: */
	    }
/* L20: */
	}

    } else if (lsame_(uplo, "L")) {
	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
	    i__2 = *m;
	    for (i = j; i <= *m; ++i) {
		i__3 = i + j * b_dim1;
		i__4 = i + j * a_dim1;
		B(i,j).r = A(i,j).r, B(i,j).i = A(i,j).i;
/* L30: */
	    }
/* L40: */
	}

    } else {
	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
	    i__2 = *m;
	    for (i = 1; i <= *m; ++i) {
		i__3 = i + j * b_dim1;
		i__4 = i + j * a_dim1;
		B(i,j).r = A(i,j).r, B(i,j).i = A(i,j).i;
/* L50: */
	    }
/* L60: */
	}
    }

    return 0;

/*     End of ZLACPY */

} /* zlacpy_ */

