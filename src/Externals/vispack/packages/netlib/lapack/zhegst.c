#include "f2c.h"

/* Subroutine */ int zhegst_(integer *itype, char *uplo, integer *n, 
	doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, 
	integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    ZHEGST reduces a complex Hermitian-definite generalized   
    eigenproblem to standard form.   

    If ITYPE = 1, the problem is A*x = lambda*B*x,   
    and A is overwritten by inv(U**H)*A*inv(U) or inv(L)*A*inv(L**H)   

    If ITYPE = 2 or 3, the problem is A*B*x = lambda*x or   
    B*A*x = lambda*x, and A is overwritten by U*A*U**H or L**H*A*L.   

    B must have been previously factorized as U**H*U or L*L**H by ZPOTRF. 
  

    Arguments   
    =========   

    ITYPE   (input) INTEGER   
            = 1: compute inv(U**H)*A*inv(U) or inv(L)*A*inv(L**H);   
            = 2 or 3: compute U*A*U**H or L**H*A*L.   

    UPLO    (input) CHARACTER   
            = 'U':  Upper triangle of A is stored and B is factored as   
                    U**H*U;   
            = 'L':  Lower triangle of A is stored and B is factored as   
                    L*L**H.   

    N       (input) INTEGER   
            The order of the matrices A and B.  N >= 0.   

    A       (input/output) COMPLEX*16 array, dimension (LDA,N)   
            On entry, the Hermitian matrix A.  If UPLO = 'U', the leading 
  
            N-by-N upper triangular part of A contains the upper   
            triangular part of the matrix A, and the strictly lower   
            triangular part of A is not referenced.  If UPLO = 'L', the   
            leading N-by-N lower triangular part of A contains the lower 
  
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
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   

    ===================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static doublecomplex c_b1 = {1.,0.};
    static doublecomplex c_b2 = {.5,0.};
    static integer c__1 = 1;
    static integer c_n1 = -1;
    static doublereal c_b18 = 1.;
    
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3;
    doublecomplex z__1;
    /* Local variables */
    static integer k;
    extern logical lsame_(char *, char *);
    extern /* Subroutine */ int zhemm_(char *, char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, doublecomplex *, integer *);
    static logical upper;
    extern /* Subroutine */ int ztrmm_(char *, char *, char *, char *, 
	    integer *, integer *, doublecomplex *, doublecomplex *, integer *,
	     doublecomplex *, integer *), 
	    ztrsm_(char *, char *, char *, char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *), zhegs2_(integer *, 
	    char *, integer *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, integer *), zher2k_(char *, char *, integer *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublereal *, doublecomplex *, 
	    integer *);
    static integer kb, nb;
    extern /* Subroutine */ int xerbla_(char *, integer *);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);




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
	xerbla_("ZHEGST", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return 0;
    }

/*     Determine the block size for this environment. */

    nb = ilaenv_(&c__1, "ZHEGST", uplo, n, &c_n1, &c_n1, &c_n1, 6L, 1L);

    if (nb <= 1 || nb >= *n) {

/*        Use unblocked code */

	zhegs2_(itype, uplo, n, &A(1,1), lda, &B(1,1), ldb, info);
    } else {

/*        Use blocked code */

	if (*itype == 1) {
	    if (upper) {

/*              Compute inv(U')*A*inv(U) */

		i__1 = *n;
		i__2 = nb;
		for (k = 1; nb < 0 ? k >= *n : k <= *n; k += nb) {
/* Computing MIN */
		    i__3 = *n - k + 1;
		    kb = min(i__3,nb);

/*                 Update the upper triangle of A(k:n,k:n)
 */

		    zhegs2_(itype, uplo, &kb, &A(k,k), lda, &B(k,k), ldb, info);
		    if (k + kb <= *n) {
			i__3 = *n - k - kb + 1;
			ztrsm_("Left", uplo, "Conjugate transpose", "Non-unit"
				, &kb, &i__3, &c_b1, &B(k,k), ldb, 
				&A(k,k+kb), lda);
			i__3 = *n - k - kb + 1;
			z__1.r = -.5, z__1.i = 0.;
			zhemm_("Left", uplo, &kb, &i__3, &z__1, &A(k,k), lda, &B(k,k+kb), ldb, 
				&c_b1, &A(k,k+kb), lda)
				;
			i__3 = *n - k - kb + 1;
			z__1.r = -1., z__1.i = 0.;
			zher2k_(uplo, "Conjugate transpose", &i__3, &kb, &
				z__1, &A(k,k+kb), lda, &B(k,k+kb), ldb, &c_b18, &A(k+kb,k+kb), lda);
			i__3 = *n - k - kb + 1;
			z__1.r = -.5, z__1.i = 0.;
			zhemm_("Left", uplo, &kb, &i__3, &z__1, &A(k,k), lda, &B(k,k+kb), ldb, 
				&c_b1, &A(k,k+kb), lda)
				;
			i__3 = *n - k - kb + 1;
			ztrsm_("Right", uplo, "No transpose", "Non-unit", &kb,
				 &i__3, &c_b1, &B(k+kb,k+kb),
				 ldb, &A(k,k+kb), lda);
		    }
/* L10: */
		}
	    } else {

/*              Compute inv(L)*A*inv(L') */

		i__2 = *n;
		i__1 = nb;
		for (k = 1; nb < 0 ? k >= *n : k <= *n; k += nb) {
/* Computing MIN */
		    i__3 = *n - k + 1;
		    kb = min(i__3,nb);

/*                 Update the lower triangle of A(k:n,k:n)
 */

		    zhegs2_(itype, uplo, &kb, &A(k,k), lda, &B(k,k), ldb, info);
		    if (k + kb <= *n) {
			i__3 = *n - k - kb + 1;
			ztrsm_("Right", uplo, "Conjugate transpose", "Non-un"
				"it", &i__3, &kb, &c_b1, &B(k,k), 
				ldb, &A(k+kb,k), lda);
			i__3 = *n - k - kb + 1;
			z__1.r = -.5, z__1.i = 0.;
			zhemm_("Right", uplo, &i__3, &kb, &z__1, &A(k,k), lda, &B(k+kb,k), ldb, &
				c_b1, &A(k+kb,k), lda);
			i__3 = *n - k - kb + 1;
			z__1.r = -1., z__1.i = 0.;
			zher2k_(uplo, "No transpose", &i__3, &kb, &z__1, &A(k+kb,k), lda, &B(k+kb,k), ldb, &c_b18, &A(k+kb,k+kb), lda);
			i__3 = *n - k - kb + 1;
			z__1.r = -.5, z__1.i = 0.;
			zhemm_("Right", uplo, &i__3, &kb, &z__1, &A(k,k), lda, &B(k+kb,k), ldb, &
				c_b1, &A(k+kb,k), lda);
			i__3 = *n - k - kb + 1;
			ztrsm_("Left", uplo, "No transpose", "Non-unit", &
				i__3, &kb, &c_b1, &B(k+kb,k+kb), ldb, &A(k+kb,k), lda);
		    }
/* L20: */
		}
	    }
	} else {
	    if (upper) {

/*              Compute U*A*U' */

		i__1 = *n;
		i__2 = nb;
		for (k = 1; nb < 0 ? k >= *n : k <= *n; k += nb) {
/* Computing MIN */
		    i__3 = *n - k + 1;
		    kb = min(i__3,nb);

/*                 Update the upper triangle of A(1:k+kb-1
,1:k+kb-1) */

		    i__3 = k - 1;
		    ztrmm_("Left", uplo, "No transpose", "Non-unit", &i__3, &
			    kb, &c_b1, &B(1,1), ldb, &A(1,k), 
			    lda);
		    i__3 = k - 1;
		    zhemm_("Right", uplo, &i__3, &kb, &c_b2, &A(k,k), lda, &B(1,k), ldb, &c_b1, &A(1,k), lda);
		    i__3 = k - 1;
		    zher2k_(uplo, "No transpose", &i__3, &kb, &c_b1, &A(1,k), lda, &B(1,k), ldb, &c_b18,
			     &A(1,1), lda);
		    i__3 = k - 1;
		    zhemm_("Right", uplo, &i__3, &kb, &c_b2, &A(k,k), lda, &B(1,k), ldb, &c_b1, &A(1,k), lda);
		    i__3 = k - 1;
		    ztrmm_("Right", uplo, "Conjugate transpose", "Non-unit", &
			    i__3, &kb, &c_b1, &B(k,k), ldb, &A(1,k), lda);
		    zhegs2_(itype, uplo, &kb, &A(k,k), lda, &B(k,k), ldb, info);
/* L30: */
		}
	    } else {

/*              Compute L'*A*L */

		i__2 = *n;
		i__1 = nb;
		for (k = 1; nb < 0 ? k >= *n : k <= *n; k += nb) {
/* Computing MIN */
		    i__3 = *n - k + 1;
		    kb = min(i__3,nb);

/*                 Update the lower triangle of A(1:k+kb-1
,1:k+kb-1) */

		    i__3 = k - 1;
		    ztrmm_("Right", uplo, "No transpose", "Non-unit", &kb, &
			    i__3, &c_b1, &B(1,1), ldb, &A(k,1), 
			    lda);
		    i__3 = k - 1;
		    zhemm_("Left", uplo, &kb, &i__3, &c_b2, &A(k,k)
			    , lda, &B(k,1), ldb, &c_b1, &A(k,1),
			     lda);
		    i__3 = k - 1;
		    zher2k_(uplo, "Conjugate transpose", &i__3, &kb, &c_b1, &
			    A(k,1), lda, &B(k,1), ldb, &c_b18, &
			    A(1,1), lda);
		    i__3 = k - 1;
		    zhemm_("Left", uplo, &kb, &i__3, &c_b2, &A(k,k)
			    , lda, &B(k,1), ldb, &c_b1, &A(k,1),
			     lda);
		    i__3 = k - 1;
		    ztrmm_("Left", uplo, "Conjugate transpose", "Non-unit", &
			    kb, &i__3, &c_b1, &B(k,k), ldb, &A(k,1), lda);
		    zhegs2_(itype, uplo, &kb, &A(k,k), lda, &B(k,k), ldb, info);
/* L40: */
		}
	    }
	}
    }
    return 0;

/*     End of ZHEGST */

} /* zhegst_ */

