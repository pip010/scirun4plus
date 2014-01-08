#include "f2c.h"

/* Subroutine */ int zsytri_(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, integer *ipiv, doublecomplex *work, integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    ZSYTRI computes the inverse of a complex symmetric indefinite matrix 
  
    A using the factorization A = U*D*U**T or A = L*D*L**T computed by   
    ZSYTRF.   

    Arguments   
    =========   

    UPLO    (input) CHARACTER*1   
            Specifies whether the details of the factorization are stored 
  
            as an upper or lower triangular matrix.   
            = 'U':  Upper triangular, form is A = U*D*U**T;   
            = 'L':  Lower triangular, form is A = L*D*L**T.   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    A       (input/output) COMPLEX*16 array, dimension (LDA,N)   
            On entry, the block diagonal matrix D and the multipliers   
            used to obtain the factor U or L as computed by ZSYTRF.   

            On exit, if INFO = 0, the (symmetric) inverse of the original 
  
            matrix.  If UPLO = 'U', the upper triangular part of the   
            inverse is formed and the part of A below the diagonal is not 
  
            referenced; if UPLO = 'L' the lower triangular part of the   
            inverse is formed and the part of A above the diagonal is   
            not referenced.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(1,N).   

    IPIV    (input) INTEGER array, dimension (N)   
            Details of the interchanges and the block structure of D   
            as determined by ZSYTRF.   

    WORK    (workspace) COMPLEX*16 array, dimension (2*N)   

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
    static doublecomplex c_b1 = {1.,0.};
    static doublecomplex c_b2 = {0.,0.};
    static integer c__1 = 1;
    
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublecomplex z__1, z__2, z__3;
    /* Builtin functions */
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *);
    /* Local variables */
    static doublecomplex temp, akkp1, d;
    static integer k;
    static doublecomplex t;
    extern logical lsame_(char *, char *);
    static integer kstep;
    static logical upper;
    extern /* Subroutine */ int zcopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    extern /* Double Complex */ VOID zdotu_(doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    extern /* Subroutine */ int zswap_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), zsymv_(char *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, doublecomplex *, integer *);
    static doublecomplex ak;
    static integer kp;
    extern /* Subroutine */ int xerbla_(char *, integer *);
    static doublecomplex akp1;



#define IPIV(I) ipiv[(I)-1]
#define WORK(I) work[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    *info = 0;
    upper = lsame_(uplo, "U");
    if (! upper && ! lsame_(uplo, "L")) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*lda < max(1,*n)) {
	*info = -4;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("ZSYTRI", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return 0;
    }

/*     Check that the diagonal matrix D is nonsingular. */

    if (upper) {

/*        Upper triangular storage: examine D from bottom to top */

	for (*info = *n; *info >= 1; --(*info)) {
	    i__1 = *info + *info * a_dim1;
	    if (IPIV(*info) > 0 && (A(*info,*info).r == 0. && A(*info,*info).i == 0.)) {
		return 0;
	    }
/* L10: */
	}
    } else {

/*        Lower triangular storage: examine D from top to bottom. */

	i__1 = *n;
	for (*info = 1; *info <= i__1; ++(*info)) {
	    i__2 = *info + *info * a_dim1;
	    if (IPIV(*info) > 0 && (A(*info,*info).r == 0. && A(*info,*info).i == 0.)) {
		return 0;
	    }
/* L20: */
	}
    }
    *info = 0;

    if (upper) {

/*        Compute inv(A) from the factorization A = U*D*U'.   

          K is the main loop index, increasing from 1 to N in steps of
   
          1 or 2, depending on the size of the diagonal blocks. */

	k = 1;
L30:

/*        If K > N, exit from loop. */

	if (k > *n) {
	    goto L40;
	}

	if (IPIV(k) > 0) {

/*           1 x 1 diagonal block   

             Invert the diagonal block. */

	    i__1 = k + k * a_dim1;
	    z_div(&z__1, &c_b1, &A(k,k));
	    A(k,k).r = z__1.r, A(k,k).i = z__1.i;

/*           Compute column K of the inverse. */

	    if (k > 1) {
		i__1 = k - 1;
		zcopy_(&i__1, &A(1,k), &c__1, &WORK(1), &c__1);
		i__1 = k - 1;
		z__1.r = -1., z__1.i = 0.;
		zsymv_(uplo, &i__1, &z__1, &A(1,1), lda, &WORK(1), &c__1,
			 &c_b2, &A(1,k), &c__1);
		i__1 = k + k * a_dim1;
		i__2 = k + k * a_dim1;
		i__3 = k - 1;
		zdotu_(&z__2, &i__3, &WORK(1), &c__1, &A(1,k), &
			c__1);
		z__1.r = A(k,k).r - z__2.r, z__1.i = A(k,k).i - z__2.i;
		A(k,k).r = z__1.r, A(k,k).i = z__1.i;
	    }
	    kstep = 1;
	} else {

/*           2 x 2 diagonal block   

             Invert the diagonal block. */

	    i__1 = k + (k + 1) * a_dim1;
	    t.r = A(k,k+1).r, t.i = A(k,k+1).i;
	    z_div(&z__1, &A(k,k), &t);
	    ak.r = z__1.r, ak.i = z__1.i;
	    z_div(&z__1, &A(k+1,k+1), &t);
	    akp1.r = z__1.r, akp1.i = z__1.i;
	    z_div(&z__1, &A(k,k+1), &t);
	    akkp1.r = z__1.r, akkp1.i = z__1.i;
	    z__3.r = ak.r * akp1.r - ak.i * akp1.i, z__3.i = ak.r * akp1.i + 
		    ak.i * akp1.r;
	    z__2.r = z__3.r - 1., z__2.i = z__3.i + 0.;
	    z__1.r = t.r * z__2.r - t.i * z__2.i, z__1.i = t.r * z__2.i + t.i 
		    * z__2.r;
	    d.r = z__1.r, d.i = z__1.i;
	    i__1 = k + k * a_dim1;
	    z_div(&z__1, &akp1, &d);
	    A(k,k).r = z__1.r, A(k,k).i = z__1.i;
	    i__1 = k + 1 + (k + 1) * a_dim1;
	    z_div(&z__1, &ak, &d);
	    A(k+1,k+1).r = z__1.r, A(k+1,k+1).i = z__1.i;
	    i__1 = k + (k + 1) * a_dim1;
	    z__2.r = -akkp1.r, z__2.i = -akkp1.i;
	    z_div(&z__1, &z__2, &d);
	    A(k,k+1).r = z__1.r, A(k,k+1).i = z__1.i;

/*           Compute columns K and K+1 of the inverse. */

	    if (k > 1) {
		i__1 = k - 1;
		zcopy_(&i__1, &A(1,k), &c__1, &WORK(1), &c__1);
		i__1 = k - 1;
		z__1.r = -1., z__1.i = 0.;
		zsymv_(uplo, &i__1, &z__1, &A(1,1), lda, &WORK(1), &c__1,
			 &c_b2, &A(1,k), &c__1);
		i__1 = k + k * a_dim1;
		i__2 = k + k * a_dim1;
		i__3 = k - 1;
		zdotu_(&z__2, &i__3, &WORK(1), &c__1, &A(1,k), &
			c__1);
		z__1.r = A(k,k).r - z__2.r, z__1.i = A(k,k).i - z__2.i;
		A(k,k).r = z__1.r, A(k,k).i = z__1.i;
		i__1 = k + (k + 1) * a_dim1;
		i__2 = k + (k + 1) * a_dim1;
		i__3 = k - 1;
		zdotu_(&z__2, &i__3, &A(1,k), &c__1, &A(1,k+1), &c__1);
		z__1.r = A(k,k+1).r - z__2.r, z__1.i = A(k,k+1).i - z__2.i;
		A(k,k+1).r = z__1.r, A(k,k+1).i = z__1.i;
		i__1 = k - 1;
		zcopy_(&i__1, &A(1,k+1), &c__1, &WORK(1), &
			c__1);
		i__1 = k - 1;
		z__1.r = -1., z__1.i = 0.;
		zsymv_(uplo, &i__1, &z__1, &A(1,1), lda, &WORK(1), &c__1,
			 &c_b2, &A(1,k+1), &c__1);
		i__1 = k + 1 + (k + 1) * a_dim1;
		i__2 = k + 1 + (k + 1) * a_dim1;
		i__3 = k - 1;
		zdotu_(&z__2, &i__3, &WORK(1), &c__1, &A(1,k+1)
			, &c__1);
		z__1.r = A(k+1,k+1).r - z__2.r, z__1.i = A(k+1,k+1).i - z__2.i;
		A(k+1,k+1).r = z__1.r, A(k+1,k+1).i = z__1.i;
	    }
	    kstep = 2;
	}

	kp = (i__1 = IPIV(k), abs(i__1));
	if (kp != k) {

/*           Interchange rows and columns K and KP in the leading 
  
             submatrix A(1:k+1,1:k+1) */

	    i__1 = kp - 1;
	    zswap_(&i__1, &A(1,k), &c__1, &A(1,kp), &
		    c__1);
	    i__1 = k - kp - 1;
	    zswap_(&i__1, &A(kp+1,k), &c__1, &A(kp,kp+1), lda);
	    i__1 = k + k * a_dim1;
	    temp.r = A(k,k).r, temp.i = A(k,k).i;
	    i__1 = k + k * a_dim1;
	    i__2 = kp + kp * a_dim1;
	    A(k,k).r = A(kp,kp).r, A(k,k).i = A(kp,kp).i;
	    i__1 = kp + kp * a_dim1;
	    A(kp,kp).r = temp.r, A(kp,kp).i = temp.i;
	    if (kstep == 2) {
		i__1 = k + (k + 1) * a_dim1;
		temp.r = A(k,k+1).r, temp.i = A(k,k+1).i;
		i__1 = k + (k + 1) * a_dim1;
		i__2 = kp + (k + 1) * a_dim1;
		A(k,k+1).r = A(kp,k+1).r, A(k,k+1).i = A(kp,k+1).i;
		i__1 = kp + (k + 1) * a_dim1;
		A(kp,k+1).r = temp.r, A(kp,k+1).i = temp.i;
	    }
	}

	k += kstep;
	goto L30;
L40:

	;
    } else {

/*        Compute inv(A) from the factorization A = L*D*L'.   

          K is the main loop index, increasing from 1 to N in steps of
   
          1 or 2, depending on the size of the diagonal blocks. */

	k = *n;
L50:

/*        If K < 1, exit from loop. */

	if (k < 1) {
	    goto L60;
	}

	if (IPIV(k) > 0) {

/*           1 x 1 diagonal block   

             Invert the diagonal block. */

	    i__1 = k + k * a_dim1;
	    z_div(&z__1, &c_b1, &A(k,k));
	    A(k,k).r = z__1.r, A(k,k).i = z__1.i;

/*           Compute column K of the inverse. */

	    if (k < *n) {
		i__1 = *n - k;
		zcopy_(&i__1, &A(k+1,k), &c__1, &WORK(1), &c__1);
		i__1 = *n - k;
		z__1.r = -1., z__1.i = 0.;
		zsymv_(uplo, &i__1, &z__1, &A(k+1,k+1), lda, 
			&WORK(1), &c__1, &c_b2, &A(k+1,k), &c__1);
		i__1 = k + k * a_dim1;
		i__2 = k + k * a_dim1;
		i__3 = *n - k;
		zdotu_(&z__2, &i__3, &WORK(1), &c__1, &A(k+1,k), 
			&c__1);
		z__1.r = A(k,k).r - z__2.r, z__1.i = A(k,k).i - z__2.i;
		A(k,k).r = z__1.r, A(k,k).i = z__1.i;
	    }
	    kstep = 1;
	} else {

/*           2 x 2 diagonal block   

             Invert the diagonal block. */

	    i__1 = k + (k - 1) * a_dim1;
	    t.r = A(k,k-1).r, t.i = A(k,k-1).i;
	    z_div(&z__1, &A(k-1,k-1), &t);
	    ak.r = z__1.r, ak.i = z__1.i;
	    z_div(&z__1, &A(k,k), &t);
	    akp1.r = z__1.r, akp1.i = z__1.i;
	    z_div(&z__1, &A(k,k-1), &t);
	    akkp1.r = z__1.r, akkp1.i = z__1.i;
	    z__3.r = ak.r * akp1.r - ak.i * akp1.i, z__3.i = ak.r * akp1.i + 
		    ak.i * akp1.r;
	    z__2.r = z__3.r - 1., z__2.i = z__3.i + 0.;
	    z__1.r = t.r * z__2.r - t.i * z__2.i, z__1.i = t.r * z__2.i + t.i 
		    * z__2.r;
	    d.r = z__1.r, d.i = z__1.i;
	    i__1 = k - 1 + (k - 1) * a_dim1;
	    z_div(&z__1, &akp1, &d);
	    A(k-1,k-1).r = z__1.r, A(k-1,k-1).i = z__1.i;
	    i__1 = k + k * a_dim1;
	    z_div(&z__1, &ak, &d);
	    A(k,k).r = z__1.r, A(k,k).i = z__1.i;
	    i__1 = k + (k - 1) * a_dim1;
	    z__2.r = -akkp1.r, z__2.i = -akkp1.i;
	    z_div(&z__1, &z__2, &d);
	    A(k,k-1).r = z__1.r, A(k,k-1).i = z__1.i;

/*           Compute columns K-1 and K of the inverse. */

	    if (k < *n) {
		i__1 = *n - k;
		zcopy_(&i__1, &A(k+1,k), &c__1, &WORK(1), &c__1);
		i__1 = *n - k;
		z__1.r = -1., z__1.i = 0.;
		zsymv_(uplo, &i__1, &z__1, &A(k+1,k+1), lda, 
			&WORK(1), &c__1, &c_b2, &A(k+1,k), &c__1);
		i__1 = k + k * a_dim1;
		i__2 = k + k * a_dim1;
		i__3 = *n - k;
		zdotu_(&z__2, &i__3, &WORK(1), &c__1, &A(k+1,k), 
			&c__1);
		z__1.r = A(k,k).r - z__2.r, z__1.i = A(k,k).i - z__2.i;
		A(k,k).r = z__1.r, A(k,k).i = z__1.i;
		i__1 = k + (k - 1) * a_dim1;
		i__2 = k + (k - 1) * a_dim1;
		i__3 = *n - k;
		zdotu_(&z__2, &i__3, &A(k+1,k), &c__1, &A(k+1,k-1), &c__1);
		z__1.r = A(k,k-1).r - z__2.r, z__1.i = A(k,k-1).i - z__2.i;
		A(k,k-1).r = z__1.r, A(k,k-1).i = z__1.i;
		i__1 = *n - k;
		zcopy_(&i__1, &A(k+1,k-1), &c__1, &WORK(1), &
			c__1);
		i__1 = *n - k;
		z__1.r = -1., z__1.i = 0.;
		zsymv_(uplo, &i__1, &z__1, &A(k+1,k+1), lda, 
			&WORK(1), &c__1, &c_b2, &A(k+1,k-1), 
			&c__1);
		i__1 = k - 1 + (k - 1) * a_dim1;
		i__2 = k - 1 + (k - 1) * a_dim1;
		i__3 = *n - k;
		zdotu_(&z__2, &i__3, &WORK(1), &c__1, &A(k+1,k-1), &c__1);
		z__1.r = A(k-1,k-1).r - z__2.r, z__1.i = A(k-1,k-1).i - z__2.i;
		A(k-1,k-1).r = z__1.r, A(k-1,k-1).i = z__1.i;
	    }
	    kstep = 2;
	}

	kp = (i__1 = IPIV(k), abs(i__1));
	if (kp != k) {

/*           Interchange rows and columns K and KP in the trailing
   
             submatrix A(k-1:n,k-1:n) */

	    if (kp < *n) {
		i__1 = *n - kp;
		zswap_(&i__1, &A(kp+1,k), &c__1, &A(kp+1,kp), &c__1);
	    }
	    i__1 = kp - k - 1;
	    zswap_(&i__1, &A(k+1,k), &c__1, &A(kp,k+1), lda);
	    i__1 = k + k * a_dim1;
	    temp.r = A(k,k).r, temp.i = A(k,k).i;
	    i__1 = k + k * a_dim1;
	    i__2 = kp + kp * a_dim1;
	    A(k,k).r = A(kp,kp).r, A(k,k).i = A(kp,kp).i;
	    i__1 = kp + kp * a_dim1;
	    A(kp,kp).r = temp.r, A(kp,kp).i = temp.i;
	    if (kstep == 2) {
		i__1 = k + (k - 1) * a_dim1;
		temp.r = A(k,k-1).r, temp.i = A(k,k-1).i;
		i__1 = k + (k - 1) * a_dim1;
		i__2 = kp + (k - 1) * a_dim1;
		A(k,k-1).r = A(kp,k-1).r, A(k,k-1).i = A(kp,k-1).i;
		i__1 = kp + (k - 1) * a_dim1;
		A(kp,k-1).r = temp.r, A(kp,k-1).i = temp.i;
	    }
	}

	k -= kstep;
	goto L50;
L60:
	;
    }

    return 0;

/*     End of ZSYTRI */

} /* zsytri_ */

