#include "f2c.h"

/* Subroutine */ int ssysvx_(char *fact, char *uplo, integer *n, integer *
	nrhs, real *a, integer *lda, real *af, integer *ldaf, integer *ipiv, 
	real *b, integer *ldb, real *x, integer *ldx, real *rcond, real *ferr,
	 real *berr, real *work, integer *lwork, integer *iwork, integer *
	info)
{
/*  -- LAPACK driver routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    SSYSVX uses the diagonal pivoting factorization to compute the   
    solution to a real system of linear equations A * X = B,   
    where A is an N-by-N symmetric matrix and X and B are N-by-NRHS   
    matrices.   

    Error bounds on the solution and a condition estimate are also   
    provided.   

    Description   
    ===========   

    The following steps are performed:   

    1. If FACT = 'N', the diagonal pivoting method is used to factor A.   
       The form of the factorization is   
          A = U * D * U**T,  if UPLO = 'U', or   
          A = L * D * L**T,  if UPLO = 'L',   
       where U (or L) is a product of permutation and unit upper (lower) 
  
       triangular matrices, and D is symmetric and block diagonal with   
       1-by-1 and 2-by-2 diagonal blocks.   

    2. The factored form of A is used to estimate the condition number   
       of the matrix A.  If the reciprocal of the condition number is   
       less than machine precision, steps 3 and 4 are skipped.   

    3. The system of equations is solved for X using the factored form   
       of A.   

    4. Iterative refinement is applied to improve the computed solution   
       matrix and calculate error bounds and backward error estimates   
       for it.   

    Arguments   
    =========   

    FACT    (input) CHARACTER*1   
            Specifies whether or not the factored form of A has been   
            supplied on entry.   
            = 'F':  On entry, AF and IPIV contain the factored form of   
                    A.  AF and IPIV will not be modified.   
            = 'N':  The matrix A will be copied to AF and factored.   

    UPLO    (input) CHARACTER*1   
            = 'U':  Upper triangle of A is stored;   
            = 'L':  Lower triangle of A is stored.   

    N       (input) INTEGER   
            The number of linear equations, i.e., the order of the   
            matrix A.  N >= 0.   

    NRHS    (input) INTEGER   
            The number of right hand sides, i.e., the number of columns   
            of the matrices B and X.  NRHS >= 0.   

    A       (input) REAL array, dimension (LDA,N)   
            The symmetric matrix A.  If UPLO = 'U', the leading N-by-N   
            upper triangular part of A contains the upper triangular part 
  
            of the matrix A, and the strictly lower triangular part of A 
  
            is not referenced.  If UPLO = 'L', the leading N-by-N lower   
            triangular part of A contains the lower triangular part of   
            the matrix A, and the strictly upper triangular part of A is 
  
            not referenced.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(1,N).   

    AF      (input or output) REAL array, dimension (LDAF,N)   
            If FACT = 'F', then AF is an input argument and on entry   
            contains the block diagonal matrix D and the multipliers used 
  
            to obtain the factor U or L from the factorization   
            A = U*D*U**T or A = L*D*L**T as computed by SSYTRF.   

            If FACT = 'N', then AF is an output argument and on exit   
            returns the block diagonal matrix D and the multipliers used 
  
            to obtain the factor U or L from the factorization   
            A = U*D*U**T or A = L*D*L**T.   

    LDAF    (input) INTEGER   
            The leading dimension of the array AF.  LDAF >= max(1,N).   

    IPIV    (input or output) INTEGER array, dimension (N)   
            If FACT = 'F', then IPIV is an input argument and on entry   
            contains details of the interchanges and the block structure 
  
            of D, as determined by SSYTRF.   
            If IPIV(k) > 0, then rows and columns k and IPIV(k) were   
            interchanged and D(k,k) is a 1-by-1 diagonal block.   
            If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0, then rows and   
            columns k-1 and -IPIV(k) were interchanged and D(k-1:k,k-1:k) 
  
            is a 2-by-2 diagonal block.  If UPLO = 'L' and IPIV(k) =   
            IPIV(k+1) < 0, then rows and columns k+1 and -IPIV(k) were   
            interchanged and D(k:k+1,k:k+1) is a 2-by-2 diagonal block.   

            If FACT = 'N', then IPIV is an output argument and on exit   
            contains details of the interchanges and the block structure 
  
            of D, as determined by SSYTRF.   

    B       (input) REAL array, dimension (LDB,NRHS)   
            The N-by-NRHS right hand side matrix B.   

    LDB     (input) INTEGER   
            The leading dimension of the array B.  LDB >= max(1,N).   

    X       (output) REAL array, dimension (LDX,NRHS)   
            If INFO = 0, the N-by-NRHS solution matrix X.   

    LDX     (input) INTEGER   
            The leading dimension of the array X.  LDX >= max(1,N).   

    RCOND   (output) REAL   
            The estimate of the reciprocal condition number of the matrix 
  
            A.  If RCOND is less than the machine precision (in   
            particular, if RCOND = 0), the matrix is singular to working 
  
            precision.  This condition is indicated by a return code of   
            INFO > 0, and the solution and error bounds are not computed. 
  

    FERR    (output) REAL array, dimension (NRHS)   
            The estimated forward error bound for each solution vector   
            X(j) (the j-th column of the solution matrix X).   
            If XTRUE is the true solution corresponding to X(j), FERR(j) 
  
            is an estimated upper bound for the magnitude of the largest 
  
            element in (X(j) - XTRUE) divided by the magnitude of the   
            largest element in X(j).  The estimate is as reliable as   
            the estimate for RCOND, and is almost always a slight   
            overestimate of the true error.   

    BERR    (output) REAL array, dimension (NRHS)   
            The componentwise relative backward error of each solution   
            vector X(j) (i.e., the smallest relative change in   
            any element of A or B that makes X(j) an exact solution).   

    WORK    (workspace/output) REAL array, dimension (LWORK)   
            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INTEGER   
            The length of WORK.  LWORK >= 3*N, and for best performance   
            LWORK >= N*NB, where NB is the optimal blocksize for   
            SSYTRF.   

    IWORK   (workspace) INTEGER array, dimension (N)   

    INFO    (output) INTEGER   
            = 0: successful exit   
            < 0: if INFO = -i, the i-th argument had an illegal value   
            > 0: if INFO = i, and i is   
                 <= N: D(i,i) is exactly zero.  The factorization has   
                       has been completed, but the block diagonal   
                       matrix D is exactly singular, so the solution and 
  
                       error bounds could not be computed.   
                 = N+1: the block diagonal matrix D is nonsingular, but   
                       RCOND is less than machine precision.  The   
                       factorization has been completed, but the matrix   
                       is singular to working precision, so the solution 
  
                       and error bounds have not been computed.   

    ===================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    integer a_dim1, a_offset, af_dim1, af_offset, b_dim1, b_offset, x_dim1, 
	    x_offset, i__1;
    /* Local variables */
    extern logical lsame_(char *, char *);
    static real anorm;
    extern doublereal slamch_(char *);
    static logical nofact;
    extern /* Subroutine */ int xerbla_(char *, integer *), slacpy_(
	    char *, integer *, integer *, real *, integer *, real *, integer *
	    );
    extern doublereal slansy_(char *, char *, integer *, real *, integer *, 
	    real *);
    extern /* Subroutine */ int ssycon_(char *, integer *, real *, integer *, 
	    integer *, real *, real *, real *, integer *, integer *), 
	    ssyrfs_(char *, integer *, integer *, real *, integer *, real *, 
	    integer *, integer *, real *, integer *, real *, integer *, real *
	    , real *, real *, integer *, integer *), ssytrf_(char *, 
	    integer *, real *, integer *, integer *, real *, integer *, 
	    integer *), ssytrs_(char *, integer *, integer *, real *, 
	    integer *, integer *, real *, integer *, integer *);


#define IPIV(I) ipiv[(I)-1]
#define FERR(I) ferr[(I)-1]
#define BERR(I) berr[(I)-1]
#define WORK(I) work[(I)-1]
#define IWORK(I) iwork[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]
#define AF(I,J) af[(I)-1 + ((J)-1)* ( *ldaf)]
#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]
#define X(I,J) x[(I)-1 + ((J)-1)* ( *ldx)]

    *info = 0;
    nofact = lsame_(fact, "N");
    if (! nofact && ! lsame_(fact, "F")) {
	*info = -1;
    } else if (! lsame_(uplo, "U") && ! lsame_(uplo, "L")) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*nrhs < 0) {
	*info = -4;
    } else if (*lda < max(1,*n)) {
	*info = -6;
    } else if (*ldaf < max(1,*n)) {
	*info = -8;
    } else if (*ldb < max(1,*n)) {
	*info = -11;
    } else if (*ldx < max(1,*n)) {
	*info = -13;
    } else if (*lwork < *n << 1) {
	*info = -18;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("SSYSVX", &i__1);
	return 0;
    }

    if (nofact) {

/*        Compute the factorization A = U*D*U' or A = L*D*L'. */

	slacpy_(uplo, n, n, &A(1,1), lda, &AF(1,1), ldaf);
	ssytrf_(uplo, n, &AF(1,1), ldaf, &IPIV(1), &WORK(1), lwork, 
		info);

/*        Return if INFO is non-zero. */

	if (*info != 0) {
	    if (*info > 0) {
		*rcond = 0.f;
	    }
	    return 0;
	}
    }

/*     Compute the norm of the matrix A. */

    anorm = slansy_("I", uplo, n, &A(1,1), lda, &WORK(1));

/*     Compute the reciprocal of the condition number of A. */

    ssycon_(uplo, n, &AF(1,1), ldaf, &IPIV(1), &anorm, rcond, &WORK(1), 
	    &IWORK(1), info);

/*     Return if the matrix is singular to working precision. */

    if (*rcond < slamch_("Epsilon")) {
	*info = *n + 1;
	return 0;
    }

/*     Compute the solution vectors X. */

    slacpy_("Full", n, nrhs, &B(1,1), ldb, &X(1,1), ldx);
    ssytrs_(uplo, n, nrhs, &AF(1,1), ldaf, &IPIV(1), &X(1,1), ldx, 
	    info);

/*     Use iterative refinement to improve the computed solutions and   
       compute error bounds and backward error estimates for them. */

    ssyrfs_(uplo, n, nrhs, &A(1,1), lda, &AF(1,1), ldaf, &IPIV(1), 
	    &B(1,1), ldb, &X(1,1), ldx, &FERR(1), &BERR(1), &WORK(1)
	    , &IWORK(1), info);

    return 0;

/*     End of SSYSVX */

} /* ssysvx_ */

