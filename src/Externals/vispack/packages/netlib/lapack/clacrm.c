#include "f2c.h"

/* Subroutine */ int clacrm_(integer *m, integer *n, complex *a, integer *lda,
	 real *b, integer *ldb, complex *c, integer *ldc, real *rwork)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    CLACRM performs a very simple matrix-matrix multiplication:   
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

    A       (input) COMPLEX array, dimension (LDA, N)   
            A contains the M by N matrix A.   

    LDA     (input) INTEGER   
            The leading dimension of the array A. LDA >=max(1,M).   

    B       (input) REAL array, dimension (LDB, N)   
            B contains the N by N matrix B.   

    LDB     (input) INTEGER   
            The leading dimension of the array B. LDB >=max(1,N).   

    C       (input) COMPLEX array, dimension (LDC, N)   
            C contains the M by N matrix C.   

    LDC     (input) INTEGER   
            The leading dimension of the array C. LDC >=max(1,N).   

    RWORK   (workspace) REAL array, dimension (2*M*N)   

    ===================================================================== 
  


       Quick return if possible.   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static real c_b6 = 1.f;
    static real c_b7 = 0.f;
    
    /* System generated locals */
    integer b_dim1, b_offset, a_dim1, a_offset, c_dim1, c_offset, i__1, i__2, 
	    i__3, i__4, i__5;
    doublereal d__1;
    complex q__1;
    /* Builtin functions */
    double r_imag(complex *);
    /* Local variables */
    static integer i, j, l;
    extern /* Subroutine */ int sgemm_(char *, char *, integer *, integer *, 
	    integer *, real *, real *, integer *, real *, integer *, real *, 
	    real *, integer *);



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
    sgemm_("N", "N", m, n, n, &c_b6, &RWORK(1), m, &B(1,1), ldb, &c_b7, &
	    RWORK(l), m);
    i__1 = *n;
    for (j = 1; j <= *n; ++j) {
	i__2 = *m;
	for (i = 1; i <= *m; ++i) {
	    i__3 = i + j * c_dim1;
	    i__4 = l + (j - 1) * *m + i - 1;
	    C(i,j).r = RWORK(l+(j-1)**m+i-1), C(i,j).i = 0.f;
/* L30: */
	}
/* L40: */
    }

    i__1 = *n;
    for (j = 1; j <= *n; ++j) {
	i__2 = *m;
	for (i = 1; i <= *m; ++i) {
	    RWORK((j - 1) * *m + i) = r_imag(&A(i,j));
/* L50: */
	}
/* L60: */
    }
    sgemm_("N", "N", m, n, n, &c_b6, &RWORK(1), m, &B(1,1), ldb, &c_b7, &
	    RWORK(l), m);
    i__1 = *n;
    for (j = 1; j <= *n; ++j) {
	i__2 = *m;
	for (i = 1; i <= *m; ++i) {
	    i__3 = i + j * c_dim1;
	    i__4 = i + j * c_dim1;
	    d__1 = C(i,j).r;
	    i__5 = l + (j - 1) * *m + i - 1;
	    q__1.r = d__1, q__1.i = RWORK(l+(j-1)**m+i-1);
	    C(i,j).r = q__1.r, C(i,j).i = q__1.i;
/* L70: */
	}
/* L80: */
    }

    return 0;

/*     End of CLACRM */

} /* clacrm_ */

