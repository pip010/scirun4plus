#include "f2c.h"

/* Subroutine */ int zpteqr_(char *compz, integer *n, doublereal *d, 
	doublereal *e, doublecomplex *z, integer *ldz, doublereal *work, 
	integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    ZPTEQR computes all eigenvalues and, optionally, eigenvectors of a   
    symmetric positive definite tridiagonal matrix by first factoring the 
  
    matrix using DPTTRF and then calling ZBDSQR to compute the singular   
    values of the bidiagonal factor.   

    This routine computes the eigenvalues of the positive definite   
    tridiagonal matrix to high relative accuracy.  This means that if the 
  
    eigenvalues range over many orders of magnitude in size, then the   
    small eigenvalues and corresponding eigenvectors will be computed   
    more accurately than, for example, with the standard QR method.   

    The eigenvectors of a full or band positive definite Hermitian matrix 
  
    can also be found if ZHETRD, ZHPTRD, or ZHBTRD has been used to   
    reduce this matrix to tridiagonal form.  (The reduction to   
    tridiagonal form, however, may preclude the possibility of obtaining 
  
    high relative accuracy in the small eigenvalues of the original   
    matrix, if these eigenvalues range over many orders of magnitude.)   

    Arguments   
    =========   

    COMPZ   (input) CHARACTER*1   
            = 'N':  Compute eigenvalues only.   
            = 'V':  Compute eigenvectors of original Hermitian   
                    matrix also.  VISArray Z contains the unitary matrix   
                    used to reduce the original matrix to tridiagonal   
                    form.   
            = 'I':  Compute eigenvectors of tridiagonal matrix also.   

    N       (input) INTEGER   
            The order of the matrix.  N >= 0.   

    D       (input/output) DOUBLE PRECISION array, dimension (N)   
            On entry, the n diagonal elements of the tridiagonal matrix. 
  
            On normal exit, D contains the eigenvalues, in descending   
            order.   

    E       (input/output) DOUBLE PRECISION array, dimension (N-1)   
            On entry, the (n-1) subdiagonal elements of the tridiagonal   
            matrix.   
            On exit, E has been destroyed.   

    Z       (input/output) COMPLEX*16 array, dimension (LDZ, N)   
            On entry, if COMPZ = 'V', the unitary matrix used in the   
            reduction to tridiagonal form.   
            On exit, if COMPZ = 'V', the orthonormal eigenvectors of the 
  
            original Hermitian matrix;   
            if COMPZ = 'I', the orthonormal eigenvectors of the   
            tridiagonal matrix.   
            If INFO > 0 on exit, Z contains the eigenvectors associated   
            with only the stored eigenvalues.   
            If  COMPZ = 'N', then Z is not referenced.   

    LDZ     (input) INTEGER   
            The leading dimension of the array Z.  LDZ >= 1, and if   
            COMPZ = 'V' or 'I', LDZ >= max(1,N).   

    WORK    (workspace) DOUBLE PRECISION array, dimension (LWORK)   
            If  COMPZ = 'N', then LWORK = 2*N   
            If  COMPZ = 'V' or 'I', then LWORK = MAX(1,4*N-4)   

    INFO    (output) INTEGER   
            = 0:  successful exit.   
            < 0:  if INFO = -i, the i-th argument had an illegal value.   
            > 0:  if INFO = i, and i is:   
                  <= N  the Cholesky factorization of the matrix could   
                        not be performed because the i-th principal minor 
  
                        was not positive definite.   
                  > N   the SVD algorithm failed to converge;   
                        if INFO = N+i, i off-diagonal elements of the   
                        bidiagonal factor did not converge to zero.   

    ==================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static doublecomplex c_b1 = {0.,0.};
    static doublecomplex c_b2 = {1.,0.};
    static integer c__0 = 0;
    static integer c__1 = 1;
    
    /* System generated locals */
    integer z_dim1, z_offset, i__1;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    static doublecomplex c[1]	/* was [1][1] */;
    static integer i;
    extern logical lsame_(char *, char *);
    static doublecomplex vt[1]	/* was [1][1] */;
    extern /* Subroutine */ int xerbla_(char *, integer *);
    static integer icompz;
    extern /* Subroutine */ int zlaset_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, doublecomplex *, integer *), dpttrf_(integer *, doublereal *, doublereal *, integer *)
	    , zbdsqr_(char *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublereal *, integer *);
    static integer nru;



#define C(I) c[(I)]
#define WAS(I) was[(I)]
#define VT(I) vt[(I)]
#define WAS(I) was[(I)]
#define D(I) d[(I)-1]
#define E(I) e[(I)-1]
#define WORK(I) work[(I)-1]

#define Z(I,J) z[(I)-1 + ((J)-1)* ( *ldz)]

    *info = 0;

    if (lsame_(compz, "N")) {
	icompz = 0;
    } else if (lsame_(compz, "V")) {
	icompz = 1;
    } else if (lsame_(compz, "I")) {
	icompz = 2;
    } else {
	icompz = -1;
    }
    if (icompz < 0) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*ldz < 1 || icompz > 0 && *ldz < max(1,*n)) {
	*info = -6;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("ZPTEQR", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return 0;
    }

    if (*n == 1) {
	if (icompz > 0) {
	    i__1 = z_dim1 + 1;
	    Z(1,1).r = 1., Z(1,1).i = 0.;
	}
	return 0;
    }
    if (icompz == 2) {
	zlaset_("Full", n, n, &c_b1, &c_b2, &Z(1,1), ldz);
    }

/*     Call DPTTRF to factor the matrix. */

    dpttrf_(n, &D(1), &E(1), info);
    if (*info != 0) {
	return 0;
    }
    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
	D(i) = sqrt(D(i));
/* L10: */
    }
    i__1 = *n - 1;
    for (i = 1; i <= *n-1; ++i) {
	E(i) *= D(i);
/* L20: */
    }

/*     Call ZBDSQR to compute the singular values/vectors of the   
       bidiagonal factor. */

    if (icompz > 0) {
	nru = *n;
    } else {
	nru = 0;
    }
    zbdsqr_("Lower", n, &c__0, &nru, &c__0, &D(1), &E(1), vt, &c__1, &Z(1,1), ldz, c, &c__1, &WORK(1), info);

/*     Square the singular values. */

    if (*info == 0) {
	i__1 = *n;
	for (i = 1; i <= *n; ++i) {
	    D(i) *= D(i);
/* L30: */
	}
    } else {
	*info = *n + *info;
    }

    return 0;

/*     End of ZPTEQR */

} /* zpteqr_ */

