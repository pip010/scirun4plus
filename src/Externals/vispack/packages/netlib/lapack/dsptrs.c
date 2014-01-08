#include "f2c.h"

/* Subroutine */ int dsptrs_(char *uplo, integer *n, integer *nrhs, 
	doublereal *ap, integer *ipiv, doublereal *b, integer *ldb, integer *
	info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       March 31, 1993   


    Purpose   
    =======   

    DSPTRS solves a system of linear equations A*X = B with a real   
    symmetric matrix A stored in packed format using the factorization   
    A = U*D*U**T or A = L*D*L**T computed by DSPTRF.   

    Arguments   
    =========   

    UPLO    (input) CHARACTER*1   
            Specifies whether the details of the factorization are stored 
  
            as an upper or lower triangular matrix.   
            = 'U':  Upper triangular, form is A = U*D*U**T;   
            = 'L':  Lower triangular, form is A = L*D*L**T.   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    NRHS    (input) INTEGER   
            The number of right hand sides, i.e., the number of columns   
            of the matrix B.  NRHS >= 0.   

    AP      (input) DOUBLE PRECISION array, dimension (N*(N+1)/2)   
            The block diagonal matrix D and the multipliers used to   
            obtain the factor U or L as computed by DSPTRF, stored as a   
            packed triangular matrix.   

    IPIV    (input) INTEGER array, dimension (N)   
            Details of the interchanges and the block structure of D   
            as determined by DSPTRF.   

    B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)   
            On entry, the right hand side matrix B.   
            On exit, the solution matrix X.   

    LDB     (input) INTEGER   
            The leading dimension of the array B.  LDB >= max(1,N).   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0: if INFO = -i, the i-th argument had an illegal value   

    ===================================================================== 
  


    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static doublereal c_b7 = -1.;
    static integer c__1 = 1;
    static doublereal c_b19 = 1.;
    
    /* System generated locals */
    integer b_dim1, b_offset, i__1;
    doublereal d__1;
    /* Local variables */
    extern /* Subroutine */ int dger_(integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static doublereal akm1k;
    static integer j, k;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    extern logical lsame_(char *, char *);
    static doublereal denom;
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *), dswap_(integer *, 
	    doublereal *, integer *, doublereal *, integer *);
    static logical upper;
    static doublereal ak, bk;
    static integer kc, kp;
    extern /* Subroutine */ int xerbla_(char *, integer *);
    static doublereal akm1, bkm1;



#define AP(I) ap[(I)-1]
#define IPIV(I) ipiv[(I)-1]

#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]

    *info = 0;
    upper = lsame_(uplo, "U");
    if (! upper && ! lsame_(uplo, "L")) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*nrhs < 0) {
	*info = -3;
    } else if (*ldb < max(1,*n)) {
	*info = -7;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DSPTRS", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0 || *nrhs == 0) {
	return 0;
    }

    if (upper) {

/*        Solve A*X = B, where A = U*D*U'.   

          First solve U*D*X = B, overwriting B with X.   

          K is the main loop index, decreasing from N to 1 in steps of
   
          1 or 2, depending on the size of the diagonal blocks. */

	k = *n;
	kc = *n * (*n + 1) / 2 + 1;
L10:

/*        If K < 1, exit from loop. */

	if (k < 1) {
	    goto L30;
	}

	kc -= k;
	if (IPIV(k) > 0) {

/*           1 x 1 diagonal block   

             Interchange rows K and IPIV(K). */

	    kp = IPIV(k);
	    if (kp != k) {
		dswap_(nrhs, &B(k,1), ldb, &B(kp,1), ldb);
	    }

/*           Multiply by inv(U(K)), where U(K) is the transformati
on   
             stored in column K of A. */

	    i__1 = k - 1;
	    dger_(&i__1, nrhs, &c_b7, &AP(kc), &c__1, &B(k,1), ldb, &B(1,1), ldb);

/*           Multiply by the inverse of the diagonal block. */

	    d__1 = 1. / AP(kc + k - 1);
	    dscal_(nrhs, &d__1, &B(k,1), ldb);
	    --k;
	} else {

/*           2 x 2 diagonal block   

             Interchange rows K-1 and -IPIV(K). */

	    kp = -IPIV(k);
	    if (kp != k - 1) {
		dswap_(nrhs, &B(k-1,1), ldb, &B(kp,1), ldb);
	    }

/*           Multiply by inv(U(K)), where U(K) is the transformati
on   
             stored in columns K-1 and K of A. */

	    i__1 = k - 2;
	    dger_(&i__1, nrhs, &c_b7, &AP(kc), &c__1, &B(k,1), ldb, &B(1,1), ldb);
	    i__1 = k - 2;
	    dger_(&i__1, nrhs, &c_b7, &AP(kc - (k - 1)), &c__1, &B(k-1,1), ldb, &B(1,1), ldb);

/*           Multiply by the inverse of the diagonal block. */

	    akm1k = AP(kc + k - 2);
	    akm1 = AP(kc - 1) / akm1k;
	    ak = AP(kc + k - 1) / akm1k;
	    denom = akm1 * ak - 1.;
	    i__1 = *nrhs;
	    for (j = 1; j <= *nrhs; ++j) {
		bkm1 = B(k-1,j) / akm1k;
		bk = B(k,j) / akm1k;
		B(k-1,j) = (ak * bkm1 - bk) / denom;
		B(k,j) = (akm1 * bk - bkm1) / denom;
/* L20: */
	    }
	    kc = kc - k + 1;
	    k += -2;
	}

	goto L10;
L30:

/*        Next solve U'*X = B, overwriting B with X.   

          K is the main loop index, increasing from 1 to N in steps of
   
          1 or 2, depending on the size of the diagonal blocks. */

	k = 1;
	kc = 1;
L40:

/*        If K > N, exit from loop. */

	if (k > *n) {
	    goto L50;
	}

	if (IPIV(k) > 0) {

/*           1 x 1 diagonal block   

             Multiply by inv(U'(K)), where U(K) is the transformat
ion   
             stored in column K of A. */

	    i__1 = k - 1;
	    dgemv_("Transpose", &i__1, nrhs, &c_b7, &B(1,1), ldb, &AP(kc)
		    , &c__1, &c_b19, &B(k,1), ldb);

/*           Interchange rows K and IPIV(K). */

	    kp = IPIV(k);
	    if (kp != k) {
		dswap_(nrhs, &B(k,1), ldb, &B(kp,1), ldb);
	    }
	    kc += k;
	    ++k;
	} else {

/*           2 x 2 diagonal block   

             Multiply by inv(U'(K+1)), where U(K+1) is the transfo
rmation   
             stored in columns K and K+1 of A. */

	    i__1 = k - 1;
	    dgemv_("Transpose", &i__1, nrhs, &c_b7, &B(1,1), ldb, &AP(kc)
		    , &c__1, &c_b19, &B(k,1), ldb);
	    i__1 = k - 1;
	    dgemv_("Transpose", &i__1, nrhs, &c_b7, &B(1,1), ldb, &AP(kc 
		    + k), &c__1, &c_b19, &B(k+1,1), ldb);

/*           Interchange rows K and -IPIV(K). */

	    kp = -IPIV(k);
	    if (kp != k) {
		dswap_(nrhs, &B(k,1), ldb, &B(kp,1), ldb);
	    }
	    kc = kc + (k << 1) + 1;
	    k += 2;
	}

	goto L40;
L50:

	;
    } else {

/*        Solve A*X = B, where A = L*D*L'.   

          First solve L*D*X = B, overwriting B with X.   

          K is the main loop index, increasing from 1 to N in steps of
   
          1 or 2, depending on the size of the diagonal blocks. */

	k = 1;
	kc = 1;
L60:

/*        If K > N, exit from loop. */

	if (k > *n) {
	    goto L80;
	}

	if (IPIV(k) > 0) {

/*           1 x 1 diagonal block   

             Interchange rows K and IPIV(K). */

	    kp = IPIV(k);
	    if (kp != k) {
		dswap_(nrhs, &B(k,1), ldb, &B(kp,1), ldb);
	    }

/*           Multiply by inv(L(K)), where L(K) is the transformati
on   
             stored in column K of A. */

	    if (k < *n) {
		i__1 = *n - k;
		dger_(&i__1, nrhs, &c_b7, &AP(kc + 1), &c__1, &B(k,1), 
			ldb, &B(k+1,1), ldb);
	    }

/*           Multiply by the inverse of the diagonal block. */

	    d__1 = 1. / AP(kc);
	    dscal_(nrhs, &d__1, &B(k,1), ldb);
	    kc = kc + *n - k + 1;
	    ++k;
	} else {

/*           2 x 2 diagonal block   

             Interchange rows K+1 and -IPIV(K). */

	    kp = -IPIV(k);
	    if (kp != k + 1) {
		dswap_(nrhs, &B(k+1,1), ldb, &B(kp,1), ldb);
	    }

/*           Multiply by inv(L(K)), where L(K) is the transformati
on   
             stored in columns K and K+1 of A. */

	    if (k < *n - 1) {
		i__1 = *n - k - 1;
		dger_(&i__1, nrhs, &c_b7, &AP(kc + 2), &c__1, &B(k,1), 
			ldb, &B(k+2,1), ldb);
		i__1 = *n - k - 1;
		dger_(&i__1, nrhs, &c_b7, &AP(kc + *n - k + 2), &c__1, &B(k+1,1), ldb, &B(k+2,1), ldb);
	    }

/*           Multiply by the inverse of the diagonal block. */

	    akm1k = AP(kc + 1);
	    akm1 = AP(kc) / akm1k;
	    ak = AP(kc + *n - k + 1) / akm1k;
	    denom = akm1 * ak - 1.;
	    i__1 = *nrhs;
	    for (j = 1; j <= *nrhs; ++j) {
		bkm1 = B(k,j) / akm1k;
		bk = B(k+1,j) / akm1k;
		B(k,j) = (ak * bkm1 - bk) / denom;
		B(k+1,j) = (akm1 * bk - bkm1) / denom;
/* L70: */
	    }
	    kc = kc + (*n - k << 1) + 1;
	    k += 2;
	}

	goto L60;
L80:

/*        Next solve L'*X = B, overwriting B with X.   

          K is the main loop index, decreasing from N to 1 in steps of
   
          1 or 2, depending on the size of the diagonal blocks. */

	k = *n;
	kc = *n * (*n + 1) / 2 + 1;
L90:

/*        If K < 1, exit from loop. */

	if (k < 1) {
	    goto L100;
	}

	kc -= *n - k + 1;
	if (IPIV(k) > 0) {

/*           1 x 1 diagonal block   

             Multiply by inv(L'(K)), where L(K) is the transformat
ion   
             stored in column K of A. */

	    if (k < *n) {
		i__1 = *n - k;
		dgemv_("Transpose", &i__1, nrhs, &c_b7, &B(k+1,1), 
			ldb, &AP(kc + 1), &c__1, &c_b19, &B(k,1), ldb);
	    }

/*           Interchange rows K and IPIV(K). */

	    kp = IPIV(k);
	    if (kp != k) {
		dswap_(nrhs, &B(k,1), ldb, &B(kp,1), ldb);
	    }
	    --k;
	} else {

/*           2 x 2 diagonal block   

             Multiply by inv(L'(K-1)), where L(K-1) is the transfo
rmation   
             stored in columns K-1 and K of A. */

	    if (k < *n) {
		i__1 = *n - k;
		dgemv_("Transpose", &i__1, nrhs, &c_b7, &B(k+1,1), 
			ldb, &AP(kc + 1), &c__1, &c_b19, &B(k,1), ldb);
		i__1 = *n - k;
		dgemv_("Transpose", &i__1, nrhs, &c_b7, &B(k+1,1), 
			ldb, &AP(kc - (*n - k)), &c__1, &c_b19, &B(k-1,1), ldb);
	    }

/*           Interchange rows K and -IPIV(K). */

	    kp = -IPIV(k);
	    if (kp != k) {
		dswap_(nrhs, &B(k,1), ldb, &B(kp,1), ldb);
	    }
	    kc -= *n - k + 2;
	    k += -2;
	}

	goto L90;
L100:
	;
    }

    return 0;

/*     End of DSPTRS */

} /* dsptrs_ */

