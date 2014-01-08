#include "f2c.h"

/* Subroutine */ int zlacrm_(integer *m, integer *n, doublecomplex *a, 
	integer *lda, doublereal *b, integer *ldb, doublecomplex *c, integer *
	ldc, doublereal *rwork)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    ZLACRM performs a very simple matrix-matrix multiplication:   
             C := A * B,   
    where A is M by N and complex; B is N by N and real;   
    C is M by N and complex.   

    Arguments   
    =========   

    M       (input) INTEGER   
            The number of rows of the matrix A and of the matrix C.   
            M >= 0.   

    N       (input) INTEGER   
            The number of columns and rows of the matrix B and   
            the number of columns of the matrix C.   
            N >= 0.   

    A       (input) COMPLEX*16 array, dimension (LDA, N)   
            A contains the M by N matrix A.   

    LDA     (input) INTEGER   
            The leading dimension of the array A. LDA >=max(1,M).   

    B       (input) DOUBLE PRECISION array, dimension (LDB, N)   
            B contains the N by N matrix B.   

    LDB     (input) INTEGER   
            The leading dimension of the array B. LDB >=max(1,N).   

    C       (input) COMPLEX*16 array, dimension (LDC, N)   
            C contains the M by N matrix C.   

    LDC     (input) INTEGER   
            The leading dimension of the array C. LDC >=max(1,N).   

    RWORK   (workspace) DOUBLE PRECISION array, dimension (2*M*N)   

    ===================================================================== 
  


       Quick return if possible.   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static doublereal c_b6 = 1.;
    static doublereal c_b7 = 0.;
    
    /* System generated locals */
    integer b_dim1, b_offset, a_dim1, a_offset, c_dim1, c_offset, i__1, i__2, 
	    i__3, i__4, i__5;
    doublereal d__1;
    doublecomplex z__1;
    /* Builtin functions */
    double d_imag(doublecomplex *);
    /* Local variables */
    static integer i, j, l;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *);



#define RWORK(I) rwork[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]
#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]
#define C(I,J) c[(I)-1 + ((J)-1)* ( *ldc)]

    if (*m == 0 || *n == 0) {
	return 0;
    }

    i__1 = *n;
    for (j = 1; j <= *n; ++j) {
	i__2 = *m;
	for (i = 1; i <= *m; ++i) {
	    i__3 = i + j * a_dim1;
	    RWORK((j - 1) * *m + i) = A(i,j).r;
/* L10: */
	}
/* L20: */
    }

    l = *m * *n + 1;
    dgemm_("N", "N", m, n, n, &c_b6, &RWORK(1), m, &B(1,1), ldb, &c_b7, &
	    RWORK(l), m);
    i__1 = *n;
    for (j = 1; j <= *n; ++j) {
	i__2 = *m;
	for (i = 1; i <= *m; ++i) {
	    i__3 = i + j * c_dim1;
	    i__4 = l + (j - 1) * *m + i - 1;
	    C(i,j).r = RWORK(l+(j-1)**m+i-1), C(i,j).i = 0.;
/* L30: */
	}
/* L40: */
    }

    i__1 = *n;
    for (j = 1; j <= *n; ++j) {
	i__2 = *m;
	for (i = 1; i <= *m; ++i) {
	    RWORK((j - 1) * *m + i) = d_imag(&A(i,j));
/* L50: */
	}
/* L60: */
    }
    dgemm_("N", "N", m, n, n, &c_b6, &RWORK(1), m, &B(1,1), ldb, &c_b7, &
	    RWORK(l), m);
    i__1 = *n;
    for (j = 1; j <= *n; ++j) {
	i__2 = *m;
	for (i = 1; i <= *m; ++i) {
	    i__3 = i + j * c_dim1;
	    i__4 = i + j * c_dim1;
	    d__1 = C(i,j).r;
	    i__5 = l + (j - 1) * *m + i - 1;
	    z__1.r = d__1, z__1.i = RWORK(l+(j-1)**m+i-1);
	    C(i,j).r = z__1.r, C(i,j).i = z__1.i;
/* L70: */
	}
/* L80: */
    }

    return 0;

/*     End of ZLACRM */

} /* zlacrm_ */

