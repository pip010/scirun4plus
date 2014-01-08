#include "f2c.h"

/* Subroutine */ int zhegs2_(integer *itype, char *uplo, integer *n, 
	doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, 
	integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    ZHEGS2 reduces a complex Hermitian-definite generalized   
    eigenproblem to standard form.   

    If ITYPE = 1, the problem is A*x = lambda*B*x,   
    and A is overwritten by inv(U')*A*inv(U) or inv(L)*A*inv(L')   

    If ITYPE = 2 or 3, the problem is A*B*x = lambda*x or   
    B*A*x = lambda*x, and A is overwritten by U*A*U` or L'*A*L.   

    B must have been previously factorized as U'*U or L*L' by ZPOTRF.   

    Arguments   
    =========   

    ITYPE   (input) INTEGER   
            = 1: compute inv(U')*A*inv(U) or inv(L)*A*inv(L');   
            = 2 or 3: compute U*A*U' or L'*A*L.   

    UPLO    (input) CHARACTER   
            Specifies whether the upper or lower triangular part of the   
            Hermitian matrix A is stored, and how B has been factorized. 
  
            = 'U':  Upper triangular   
            = 'L':  Lower triangular   

    N       (input) INTEGER   
            The order of the matrices A and B.  N >= 0.   

    A       (input/output) COMPLEX*16 array, dimension (LDA,N)   
            On entry, the Hermitian matrix A.  If UPLO = 'U', the leading 
  
            n by n upper triangular part of A contains the upper   
            triangular part of the matrix A, and the strictly lower   
            triangular part of A is not referenced.  If UPLO = 'L', the   
            leading n by n lower triangular part of A contains the lower 
  
            triangular part of the matrix A, and the strictly upper   
            triangular part of A is not referenced.   

            On exit, if INFO = 0, the transformed matrix, stored in the   
            same format as A.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(1,N).   

    B       (input) COMPLEX*16 array, dimension (LDB,N)   
            The triangular factor from the Cholesky factorization of B,   
            as returned by ZPOTRF.   

    LDB     (input) INTEGER   
            The leading dimension of the array B.  LDB >= max(1,N).   

    INFO    (output) INTEGER   
            = 0:  successful exit.   
            < 0:  if INFO = -i, the i-th argument had an illegal value.   

    ===================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static doublecomplex c_b1 = {1.,0.};
    static integer c__1 = 1;
    
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2;
    doublereal d__1, d__2;
    doublecomplex z__1;
    /* Local variables */
    extern /* Subroutine */ int zher2_(char *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    static integer k;
    extern logical lsame_(char *, char *);
    static logical upper;
    extern /* Subroutine */ int zaxpy_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *), ztrmv_(
	    char *, char *, char *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), ztrsv_(char *
	    , char *, char *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    static doublecomplex ct;
    extern /* Subroutine */ int xerbla_(char *, integer *), zdscal_(
	    integer *, doublereal *, doublecomplex *, integer *), zlacgv_(
	    integer *, doublecomplex *, integer *);
    static doublereal akk, bkk;




#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]
#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]

    *info = 0;
    upper = lsame_(uplo, "U");
    if (*itype < 1 || *itype > 3) {
	*info = -1;
    } else if (! upper && ! lsame_(uplo, "L")) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*lda < max(1,*n)) {
	*info = -5;
    } else if (*ldb < max(1,*n)) {
	*info = -7;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("ZHEGS2", &i__1);
	return 0;
    }

    if (*itype == 1) {
	if (upper) {

/*           Compute inv(U')*A*inv(U) */

	    i__1 = *n;
	    for (k = 1; k <= *n; ++k) {

/*              Update the upper triangle of A(k:n,k:n) */

		i__2 = k + k * a_dim1;
		akk = A(k,k).r;
		i__2 = k + k * b_dim1;
		bkk = B(k,k).r;
/* Computing 2nd power */
		d__1 = bkk;
		akk /= d__1 * d__1;
		i__2 = k + k * a_dim1;
		A(k,k).r = akk, A(k,k).i = 0.;
		if (k < *n) {
		    i__2 = *n - k;
		    d__1 = 1. / bkk;
		    zdscal_(&i__2, &d__1, &A(k,k+1), lda);
		    d__1 = akk * -.5;
		    ct.r = d__1, ct.i = 0.;
		    i__2 = *n - k;
		    zlacgv_(&i__2, &A(k,k+1), lda);
		    i__2 = *n - k;
		    zlacgv_(&i__2, &B(k,k+1), ldb);
		    i__2 = *n - k;
		    zaxpy_(&i__2, &ct, &B(k,k+1), ldb, &A(k,k+1), lda);
		    i__2 = *n - k;
		    z__1.r = -1., z__1.i = 0.;
		    zher2_(uplo, &i__2, &z__1, &A(k,k+1), lda, 
			    &B(k,k+1), ldb, &A(k+1,k+1), lda);
		    i__2 = *n - k;
		    zaxpy_(&i__2, &ct, &B(k,k+1), ldb, &A(k,k+1), lda);
		    i__2 = *n - k;
		    zlacgv_(&i__2, &B(k,k+1), ldb);
		    i__2 = *n - k;
		    ztrsv_(uplo, "Conjugate transpose", "Non-unit", &i__2, &B(k+1,k+1), ldb, &A(k,k+1), lda);
		    i__2 = *n - k;
		    zlacgv_(&i__2, &A(k,k+1), lda);
		}
/* L10: */
	    }
	} else {

/*           Compute inv(L)*A*inv(L') */

	    i__1 = *n;
	    for (k = 1; k <= *n; ++k) {

/*              Update the lower triangle of A(k:n,k:n) */

		i__2 = k + k * a_dim1;
		akk = A(k,k).r;
		i__2 = k + k * b_dim1;
		bkk = B(k,k).r;
/* Computing 2nd power */
		d__1 = bkk;
		akk /= d__1 * d__1;
		i__2 = k + k * a_dim1;
		A(k,k).r = akk, A(k,k).i = 0.;
		if (k < *n) {
		    i__2 = *n - k;
		    d__1 = 1. / bkk;
		    zdscal_(&i__2, &d__1, &A(k+1,k), &c__1);
		    d__1 = akk * -.5;
		    ct.r = d__1, ct.i = 0.;
		    i__2 = *n - k;
		    zaxpy_(&i__2, &ct, &B(k+1,k), &c__1, &A(k+1,k), &c__1);
		    i__2 = *n - k;
		    z__1.r = -1., z__1.i = 0.;
		    zher2_(uplo, &i__2, &z__1, &A(k+1,k), &c__1, 
			    &B(k+1,k), &c__1, &A(k+1,k+1), lda);
		    i__2 = *n - k;
		    zaxpy_(&i__2, &ct, &B(k+1,k), &c__1, &A(k+1,k), &c__1);
		    i__2 = *n - k;
		    ztrsv_(uplo, "No transpose", "Non-unit", &i__2, &B(k+1,k+1), ldb, &A(k+1,k), 
			    &c__1);
		}
/* L20: */
	    }
	}
    } else {
	if (upper) {

/*           Compute U*A*U' */

	    i__1 = *n;
	    for (k = 1; k <= *n; ++k) {

/*              Update the upper triangle of A(1:k,1:k) */

		i__2 = k + k * a_dim1;
		akk = A(k,k).r;
		i__2 = k + k * b_dim1;
		bkk = B(k,k).r;
		i__2 = k - 1;
		ztrmv_(uplo, "No transpose", "Non-unit", &i__2, &B(1,1), 
			ldb, &A(1,k), &c__1);
		d__1 = akk * .5;
		ct.r = d__1, ct.i = 0.;
		i__2 = k - 1;
		zaxpy_(&i__2, &ct, &B(1,k), &c__1, &A(1,k), &c__1);
		i__2 = k - 1;
		zher2_(uplo, &i__2, &c_b1, &A(1,k), &c__1, &B(1,k), &c__1, &A(1,1), lda);
		i__2 = k - 1;
		zaxpy_(&i__2, &ct, &B(1,k), &c__1, &A(1,k), &c__1);
		i__2 = k - 1;
		zdscal_(&i__2, &bkk, &A(1,k), &c__1);
		i__2 = k + k * a_dim1;
/* Computing 2nd power */
		d__2 = bkk;
		d__1 = akk * (d__2 * d__2);
		A(k,k).r = d__1, A(k,k).i = 0.;
/* L30: */
	    }
	} else {

/*           Compute L'*A*L */

	    i__1 = *n;
	    for (k = 1; k <= *n; ++k) {

/*              Update the lower triangle of A(1:k,1:k) */

		i__2 = k + k * a_dim1;
		akk = A(k,k).r;
		i__2 = k + k * b_dim1;
		bkk = B(k,k).r;
		i__2 = k - 1;
		zlacgv_(&i__2, &A(k,1), lda);
		i__2 = k - 1;
		ztrmv_(uplo, "Conjugate transpose", "Non-unit", &i__2, &B(1,1), ldb, &A(k,1), lda);
		d__1 = akk * .5;
		ct.r = d__1, ct.i = 0.;
		i__2 = k - 1;
		zlacgv_(&i__2, &B(k,1), ldb);
		i__2 = k - 1;
		zaxpy_(&i__2, &ct, &B(k,1), ldb, &A(k,1), lda);
		i__2 = k - 1;
		zher2_(uplo, &i__2, &c_b1, &A(k,1), lda, &B(k,1)
			, ldb, &A(1,1), lda);
		i__2 = k - 1;
		zaxpy_(&i__2, &ct, &B(k,1), ldb, &A(k,1), lda);
		i__2 = k - 1;
		zlacgv_(&i__2, &B(k,1), ldb);
		i__2 = k - 1;
		zdscal_(&i__2, &bkk, &A(k,1), lda);
		i__2 = k - 1;
		zlacgv_(&i__2, &A(k,1), lda);
		i__2 = k + k * a_dim1;
/* Computing 2nd power */
		d__2 = bkk;
		d__1 = akk * (d__2 * d__2);
		A(k,k).r = d__1, A(k,k).i = 0.;
/* L40: */
	    }
	}
    }
    return 0;

/*     End of ZHEGS2 */

} /* zhegs2_ */

