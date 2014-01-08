#include "f2c.h"

/* Subroutine */ int dsptri_(char *uplo, integer *n, doublereal *ap, integer *
	ipiv, doublereal *work, integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       March 31, 1993   


    Purpose   
    =======   

    DSPTRI computes the inverse of a real symmetric indefinite matrix   
    A in packed storage using the factorization A = U*D*U**T or   
    A = L*D*L**T computed by DSPTRF.   

    Arguments   
    =========   

    UPLO    (input) CHARACTER*1   
            Specifies whether the details of the factorization are stored 
  
            as an upper or lower triangular matrix.   
            = 'U':  Upper triangular, form is A = U*D*U**T;   
            = 'L':  Lower triangular, form is A = L*D*L**T.   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    AP      (input/output) DOUBLE PRECISION array, dimension (N*(N+1)/2) 
  
            On entry, the block diagonal matrix D and the multipliers   
            used to obtain the factor U or L as computed by DSPTRF,   
            stored as a packed triangular matrix.   

            On exit, if INFO = 0, the (symmetric) inverse of the original 
  
            matrix, stored as a packed triangular matrix. The j-th column 
  
            of inv(A) is stored in the array AP as follows:   
            if UPLO = 'U', AP(i + (j-1)*j/2) = inv(A)(i,j) for 1<=i<=j;   
            if UPLO = 'L',   
               AP(i + (j-1)*(2n-j)/2) = inv(A)(i,j) for j<=i<=n.   

    IPIV    (input) INTEGER array, dimension (N)   
            Details of the interchanges and the block structure of D   
            as determined by DSPTRF.   

    WORK    (workspace) DOUBLE PRECISION array, dimension (N)   

    INFO    (output) INTEGER   
            = 0: successful exit   
            < 0: if INFO = -i, the i-th argument had an illegal value   
            > 0: if INFO = i, D(i,i) = 0; the matrix is singular and its 
  
                 inverse could not be computed.   

    ===================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static integer c__1 = 1;
    static doublereal c_b11 = -1.;
    static doublereal c_b13 = 0.;
    
    /* System generated locals */
    integer i__1;
    doublereal d__1;
    /* Local variables */
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static doublereal temp, akkp1, d;
    static integer j, k;
    static doublereal t;
    extern logical lsame_(char *, char *);
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dswap_(integer *, doublereal *, integer 
	    *, doublereal *, integer *);
    static integer kstep;
    extern /* Subroutine */ int dspmv_(char *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     integer *);
    static logical upper;
    static doublereal ak;
    static integer kc, kp, kx;
    extern /* Subroutine */ int xerbla_(char *, integer *);
    static integer kcnext, kpc, npp;
    static doublereal akp1;



#define WORK(I) work[(I)-1]
#define IPIV(I) ipiv[(I)-1]
#define AP(I) ap[(I)-1]


    *info = 0;
    upper = lsame_(uplo, "U");
    if (! upper && ! lsame_(uplo, "L")) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DSPTRI", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return 0;
    }

/*     Check that the diagonal matrix D is nonsingular. */

    if (upper) {

/*        Upper triangular storage: examine D from bottom to top */

	kp = *n * (*n + 1) / 2;
	for (*info = *n; *info >= 1; --(*info)) {
	    if (IPIV(*info) > 0 && AP(kp) == 0.) {
		return 0;
	    }
	    kp -= *info;
/* L10: */
	}
    } else {

/*        Lower triangular storage: examine D from top to bottom. */

	kp = 1;
	i__1 = *n;
	for (*info = 1; *info <= i__1; ++(*info)) {
	    if (IPIV(*info) > 0 && AP(kp) == 0.) {
		return 0;
	    }
	    kp = kp + *n - *info + 1;
/* L20: */
	}
    }
    *info = 0;

    if (upper) {

/*        Compute inv(A) from the factorization A = U*D*U'.   

          K is the main loop index, increasing from 1 to N in steps of
   
          1 or 2, depending on the size of the diagonal blocks. */

	k = 1;
	kc = 1;
L30:

/*        If K > N, exit from loop. */

	if (k > *n) {
	    goto L50;
	}

	kcnext = kc + k;
	if (IPIV(k) > 0) {

/*           1 x 1 diagonal block   

             Invert the diagonal block. */

	    AP(kc + k - 1) = 1. / AP(kc + k - 1);

/*           Compute column K of the inverse. */

	    if (k > 1) {
		i__1 = k - 1;
		dcopy_(&i__1, &AP(kc), &c__1, &WORK(1), &c__1);
		i__1 = k - 1;
		dspmv_(uplo, &i__1, &c_b11, &AP(1), &WORK(1), &c__1, &c_b13, &
			AP(kc), &c__1);
		i__1 = k - 1;
		AP(kc + k - 1) -= ddot_(&i__1, &WORK(1), &c__1, &AP(kc), &
			c__1);
	    }
	    kstep = 1;
	} else {

/*           2 x 2 diagonal block   

             Invert the diagonal block. */

	    t = (d__1 = AP(kcnext + k - 1), abs(d__1));
	    ak = AP(kc + k - 1) / t;
	    akp1 = AP(kcnext + k) / t;
	    akkp1 = AP(kcnext + k - 1) / t;
	    d = t * (ak * akp1 - 1.);
	    AP(kc + k - 1) = akp1 / d;
	    AP(kcnext + k) = ak / d;
	    AP(kcnext + k - 1) = -akkp1 / d;

/*           Compute columns K and K+1 of the inverse. */

	    if (k > 1) {
		i__1 = k - 1;
		dcopy_(&i__1, &AP(kc), &c__1, &WORK(1), &c__1);
		i__1 = k - 1;
		dspmv_(uplo, &i__1, &c_b11, &AP(1), &WORK(1), &c__1, &c_b13, &
			AP(kc), &c__1);
		i__1 = k - 1;
		AP(kc + k - 1) -= ddot_(&i__1, &WORK(1), &c__1, &AP(kc), &
			c__1);
		i__1 = k - 1;
		AP(kcnext + k - 1) -= ddot_(&i__1, &AP(kc), &c__1, &AP(kcnext)
			, &c__1);
		i__1 = k - 1;
		dcopy_(&i__1, &AP(kcnext), &c__1, &WORK(1), &c__1);
		i__1 = k - 1;
		dspmv_(uplo, &i__1, &c_b11, &AP(1), &WORK(1), &c__1, &c_b13, &
			AP(kcnext), &c__1);
		i__1 = k - 1;
		AP(kcnext + k) -= ddot_(&i__1, &WORK(1), &c__1, &AP(kcnext), &
			c__1);
	    }
	    kstep = 2;
	    kcnext = kcnext + k + 1;
	}

	kp = (i__1 = IPIV(k), abs(i__1));
	if (kp != k) {

/*           Interchange rows and columns K and KP in the leading 
  
             submatrix A(1:k+1,1:k+1) */

	    kpc = (kp - 1) * kp / 2 + 1;
	    i__1 = kp - 1;
	    dswap_(&i__1, &AP(kc), &c__1, &AP(kpc), &c__1);
	    kx = kpc + kp - 1;
	    i__1 = k - 1;
	    for (j = kp + 1; j <= k-1; ++j) {
		kx = kx + j - 1;
		temp = AP(kc + j - 1);
		AP(kc + j - 1) = AP(kx);
		AP(kx) = temp;
/* L40: */
	    }
	    temp = AP(kc + k - 1);
	    AP(kc + k - 1) = AP(kpc + kp - 1);
	    AP(kpc + kp - 1) = temp;
	    if (kstep == 2) {
		temp = AP(kc + k + k - 1);
		AP(kc + k + k - 1) = AP(kc + k + kp - 1);
		AP(kc + k + kp - 1) = temp;
	    }
	}

	k += kstep;
	kc = kcnext;
	goto L30;
L50:

	;
    } else {

/*        Compute inv(A) from the factorization A = L*D*L'.   

          K is the main loop index, increasing from 1 to N in steps of
   
          1 or 2, depending on the size of the diagonal blocks. */

	npp = *n * (*n + 1) / 2;
	k = *n;
	kc = npp;
L60:

/*        If K < 1, exit from loop. */

	if (k < 1) {
	    goto L80;
	}

	kcnext = kc - (*n - k + 2);
	if (IPIV(k) > 0) {

/*           1 x 1 diagonal block   

             Invert the diagonal block. */

	    AP(kc) = 1. / AP(kc);

/*           Compute column K of the inverse. */

	    if (k < *n) {
		i__1 = *n - k;
		dcopy_(&i__1, &AP(kc + 1), &c__1, &WORK(1), &c__1);
		i__1 = *n - k;
		dspmv_(uplo, &i__1, &c_b11, &AP(kc + *n - k + 1), &WORK(1), &
			c__1, &c_b13, &AP(kc + 1), &c__1);
		i__1 = *n - k;
		AP(kc) -= ddot_(&i__1, &WORK(1), &c__1, &AP(kc + 1), &c__1);
	    }
	    kstep = 1;
	} else {

/*           2 x 2 diagonal block   

             Invert the diagonal block. */

	    t = (d__1 = AP(kcnext + 1), abs(d__1));
	    ak = AP(kcnext) / t;
	    akp1 = AP(kc) / t;
	    akkp1 = AP(kcnext + 1) / t;
	    d = t * (ak * akp1 - 1.);
	    AP(kcnext) = akp1 / d;
	    AP(kc) = ak / d;
	    AP(kcnext + 1) = -akkp1 / d;

/*           Compute columns K-1 and K of the inverse. */

	    if (k < *n) {
		i__1 = *n - k;
		dcopy_(&i__1, &AP(kc + 1), &c__1, &WORK(1), &c__1);
		i__1 = *n - k;
		dspmv_(uplo, &i__1, &c_b11, &AP(kc + (*n - k + 1)), &WORK(1), 
			&c__1, &c_b13, &AP(kc + 1), &c__1);
		i__1 = *n - k;
		AP(kc) -= ddot_(&i__1, &WORK(1), &c__1, &AP(kc + 1), &c__1);
		i__1 = *n - k;
		AP(kcnext + 1) -= ddot_(&i__1, &AP(kc + 1), &c__1, &AP(kcnext 
			+ 2), &c__1);
		i__1 = *n - k;
		dcopy_(&i__1, &AP(kcnext + 2), &c__1, &WORK(1), &c__1);
		i__1 = *n - k;
		dspmv_(uplo, &i__1, &c_b11, &AP(kc + (*n - k + 1)), &WORK(1), 
			&c__1, &c_b13, &AP(kcnext + 2), &c__1);
		i__1 = *n - k;
		AP(kcnext) -= ddot_(&i__1, &WORK(1), &c__1, &AP(kcnext + 2), &
			c__1);
	    }
	    kstep = 2;
	    kcnext -= *n - k + 3;
	}

	kp = (i__1 = IPIV(k), abs(i__1));
	if (kp != k) {

/*           Interchange rows and columns K and KP in the trailing
   
             submatrix A(k-1:n,k-1:n) */

	    kpc = npp - (*n - kp + 1) * (*n - kp + 2) / 2 + 1;
	    if (kp < *n) {
		i__1 = *n - kp;
		dswap_(&i__1, &AP(kc + kp - k + 1), &c__1, &AP(kpc + 1), &
			c__1);
	    }
	    kx = kc + kp - k;
	    i__1 = kp - 1;
	    for (j = k + 1; j <= kp-1; ++j) {
		kx = kx + *n - j + 1;
		temp = AP(kc + j - k);
		AP(kc + j - k) = AP(kx);
		AP(kx) = temp;
/* L70: */
	    }
	    temp = AP(kc);
	    AP(kc) = AP(kpc);
	    AP(kpc) = temp;
	    if (kstep == 2) {
		temp = AP(kc - *n + k - 1);
		AP(kc - *n + k - 1) = AP(kc - *n + kp - 1);
		AP(kc - *n + kp - 1) = temp;
	    }
	}

	k -= kstep;
	kc = kcnext;
	goto L60;
L80:
	;
    }

    return 0;

/*     End of DSPTRI */

} /* dsptri_ */

