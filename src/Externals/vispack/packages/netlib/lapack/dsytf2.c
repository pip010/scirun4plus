#include "f2c.h"

/* Subroutine */ int dsytf2_(char *uplo, integer *n, doublereal *a, integer *
	lda, integer *ipiv, integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    DSYTF2 computes the factorization of a real symmetric matrix A using 
  
    the Bunch-Kaufman diagonal pivoting method:   

       A = U*D*U'  or  A = L*D*L'   

    where U (or L) is a product of permutation and unit upper (lower)   
    triangular matrices, U' is the transpose of U, and D is symmetric and 
  
    block diagonal with 1-by-1 and 2-by-2 diagonal blocks.   

    This is the unblocked version of the algorithm, calling Level 2 BLAS. 
  

    Arguments   
    =========   

    UPLO    (input) CHARACTER*1   
            Specifies whether the upper or lower triangular part of the   
            symmetric matrix A is stored:   
            = 'U':  Upper triangular   
            = 'L':  Lower triangular   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the symmetric matrix A.  If UPLO = 'U', the leading 
  
            n-by-n upper triangular part of A contains the upper   
            triangular part of the matrix A, and the strictly lower   
            triangular part of A is not referenced.  If UPLO = 'L', the   
            leading n-by-n lower triangular part of A contains the lower 
  
            triangular part of the matrix A, and the strictly upper   
            triangular part of A is not referenced.   

            On exit, the block diagonal matrix D and the multipliers used 
  
            to obtain the factor U or L (see below for further details). 
  

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(1,N).   

    IPIV    (output) INTEGER array, dimension (N)   
            Details of the interchanges and the block structure of D.   
            If IPIV(k) > 0, then rows and columns k and IPIV(k) were   
            interchanged and D(k,k) is a 1-by-1 diagonal block.   
            If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0, then rows and   
            columns k-1 and -IPIV(k) were interchanged and D(k-1:k,k-1:k) 
  
            is a 2-by-2 diagonal block.  If UPLO = 'L' and IPIV(k) =   
            IPIV(k+1) < 0, then rows and columns k+1 and -IPIV(k) were   
            interchanged and D(k:k+1,k:k+1) is a 2-by-2 diagonal block.   

    INFO    (output) INTEGER   
            = 0: successful exit   
            < 0: if INFO = -k, the k-th argument had an illegal value   
            > 0: if INFO = k, D(k,k) is exactly zero.  The factorization 
  
                 has been completed, but the block diagonal matrix D is   
                 exactly singular, and division by zero will occur if it 
  
                 is used to solve a system of equations.   

    Further Details   
    ===============   

    If UPLO = 'U', then A = U*D*U', where   
       U = P(n)*U(n)* ... *P(k)U(k)* ...,   
    i.e., U is a product of terms P(k)*U(k), where k decreases from n to 
  
    1 in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1   
    and 2-by-2 diagonal blocks D(k).  P(k) is a permutation matrix as   
    defined by IPIV(k), and U(k) is a unit upper triangular matrix, such 
  
    that if the diagonal block D(k) is of order s (s = 1 or 2), then   

               (   I    v    0   )   k-s   
       U(k) =  (   0    I    0   )   s   
               (   0    0    I   )   n-k   
                  k-s   s   n-k   

    If s = 1, D(k) overwrites A(k,k), and v overwrites A(1:k-1,k).   
    If s = 2, the upper triangle of D(k) overwrites A(k-1,k-1), A(k-1,k), 
  
    and A(k,k), and v overwrites A(1:k-2,k-1:k).   

    If UPLO = 'L', then A = L*D*L', where   
       L = P(1)*L(1)* ... *P(k)*L(k)* ...,   
    i.e., L is a product of terms P(k)*L(k), where k increases from 1 to 
  
    n in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1   
    and 2-by-2 diagonal blocks D(k).  P(k) is a permutation matrix as   
    defined by IPIV(k), and L(k) is a unit lower triangular matrix, such 
  
    that if the diagonal block D(k) is of order s (s = 1 or 2), then   

               (   I    0     0   )  k-1   
       L(k) =  (   0    I     0   )  s   
               (   0    v     I   )  n-k-s+1   
                  k-1   s  n-k-s+1   

    If s = 1, D(k) overwrites A(k,k), and v overwrites A(k+1:n,k).   
    If s = 2, the lower triangle of D(k) overwrites A(k,k), A(k+1,k),   
    and A(k+1,k+1), and v overwrites A(k+2:n,k:k+1).   

    ===================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static integer c__1 = 1;
    
    /* System generated locals */
    integer a_dim1, a_offset, i__1;
    doublereal d__1, d__2, d__3;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    static integer imax, jmax;
    extern /* Subroutine */ int drot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *), dsyr_(char *
	    , integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *);
    static doublereal c;
    static integer k;
    static doublereal s, t, alpha;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    extern logical lsame_(char *, char *);
    extern /* Subroutine */ int dswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static integer kstep;
    static logical upper;
    static doublereal r1, r2;
    extern /* Subroutine */ int dlaev2_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *);
    static integer kk, kp;
    static doublereal absakk;
    extern integer idamax_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int xerbla_(char *, integer *);
    static doublereal colmax, rowmax;



#define IPIV(I) ipiv[(I)-1]

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
	xerbla_("DSYTF2", &i__1);
	return 0;
    }

/*     Initialize ALPHA for use in choosing pivot block size. */

    alpha = (sqrt(17.) + 1.) / 8.;

    if (upper) {

/*        Factorize A as U*D*U' using the upper triangle of A   

          K is the main loop index, decreasing from N to 1 in steps of
   
          1 or 2 */

	k = *n;
L10:

/*        If K < 1, exit from loop */

	if (k < 1) {
	    goto L30;
	}
	kstep = 1;

/*        Determine rows and columns to be interchanged and whether   
          a 1-by-1 or 2-by-2 pivot block will be used */

	absakk = (d__1 = A(k,k), abs(d__1));

/*        IMAX is the row-index of the largest off-diagonal element in
   
          column K, and COLMAX is its absolute value */

	if (k > 1) {
	    i__1 = k - 1;
	    imax = idamax_(&i__1, &A(1,k), &c__1);
	    colmax = (d__1 = A(imax,k), abs(d__1));
	} else {
	    colmax = 0.;
	}

	if (max(absakk,colmax) == 0.) {

/*           Column K is zero: set INFO and continue */

	    if (*info == 0) {
		*info = k;
	    }
	    kp = k;
	} else {
	    if (absakk >= alpha * colmax) {

/*              no interchange, use 1-by-1 pivot block */

		kp = k;
	    } else {

/*              JMAX is the column-index of the largest off-di
agonal   
                element in row IMAX, and ROWMAX is its absolut
e value */

		i__1 = k - imax;
		jmax = imax + idamax_(&i__1, &A(imax,imax+1), 
			lda);
		rowmax = (d__1 = A(imax,jmax), abs(d__1));
		if (imax > 1) {
		    i__1 = imax - 1;
		    jmax = idamax_(&i__1, &A(1,imax), &c__1);
/* Computing MAX */
		    d__2 = rowmax, d__3 = (d__1 = A(jmax,imax), 
			    abs(d__1));
		    rowmax = max(d__2,d__3);
		}

		if (absakk >= alpha * colmax * (colmax / rowmax)) {

/*                 no interchange, use 1-by-1 pivot block 
*/

		    kp = k;
		} else if ((d__1 = A(imax,imax), abs(d__1)) >= 
			alpha * rowmax) {

/*                 interchange rows and columns K and IMAX
, use 1-by-1   
                   pivot block */

		    kp = imax;
		} else {

/*                 interchange rows and columns K-1 and IM
AX, use 2-by-2   
                   pivot block */

		    kp = imax;
		    kstep = 2;
		}
	    }

	    kk = k - kstep + 1;
	    if (kp != kk) {

/*              Interchange rows and columns KK and KP in the 
leading   
                submatrix A(1:k,1:k) */

		i__1 = kp - 1;
		dswap_(&i__1, &A(1,kk), &c__1, &A(1,kp),
			 &c__1);
		i__1 = kk - kp - 1;
		dswap_(&i__1, &A(kp+1,kk), &c__1, &A(kp,kp+1), lda);
		t = A(kk,kk);
		A(kk,kk) = A(kp,kp);
		A(kp,kp) = t;
		if (kstep == 2) {
		    t = A(k-1,k);
		    A(k-1,k) = A(kp,k);
		    A(kp,k) = t;
		}
	    }

/*           Update the leading submatrix */

	    if (kstep == 1) {

/*              1-by-1 pivot block D(k): column k now holds   

                W(k) = U(k)*D(k)   

                where U(k) is the k-th column of U   

                Perform a rank-1 update of A(1:k-1,1:k-1) as 
  

                A := A - U(k)*D(k)*U(k)' = A - W(k)*1/D(k)*W(k
)' */

		r1 = 1. / A(k,k);
		i__1 = k - 1;
		d__1 = -r1;
		dsyr_(uplo, &i__1, &d__1, &A(1,k), &c__1, &A(1,1), lda);

/*              Store U(k) in column k */

		i__1 = k - 1;
		dscal_(&i__1, &r1, &A(1,k), &c__1);
	    } else {

/*              2-by-2 pivot block D(k): columns k and k-1 now
 hold   

                ( W(k-1) W(k) ) = ( U(k-1) U(k) )*D(k)   

                where U(k) and U(k-1) are the k-th and (k-1)-t
h columns   
                of U   

                Perform a rank-2 update of A(1:k-2,1:k-2) as 
  

                A := A - ( U(k-1) U(k) )*D(k)*( U(k-1) U(k) )'
   
                   = A - ( W(k-1) W(k) )*inv(D(k))*( W(k-1) W(
k) )'   

                Convert this to two rank-1 updates by using th
e eigen-   
                decomposition of D(k) */

		dlaev2_(&A(k-1,k-1), &A(k-1,k), 
			&A(k,k), &r1, &r2, &c, &s);
		r1 = 1. / r1;
		r2 = 1. / r2;
		i__1 = k - 2;
		drot_(&i__1, &A(1,k-1), &c__1, &A(1,k), &c__1, &c, &s);
		i__1 = k - 2;
		d__1 = -r1;
		dsyr_(uplo, &i__1, &d__1, &A(1,k-1), &c__1, &A(1,1), lda);
		i__1 = k - 2;
		d__1 = -r2;
		dsyr_(uplo, &i__1, &d__1, &A(1,k), &c__1, &A(1,1), lda);

/*              Store U(k) and U(k-1) in columns k and k-1 */

		i__1 = k - 2;
		dscal_(&i__1, &r1, &A(1,k-1), &c__1);
		i__1 = k - 2;
		dscal_(&i__1, &r2, &A(1,k), &c__1);
		i__1 = k - 2;
		d__1 = -s;
		drot_(&i__1, &A(1,k-1), &c__1, &A(1,k), &c__1, &c, &d__1);
	    }
	}

/*        Store details of the interchanges in IPIV */

	if (kstep == 1) {
	    IPIV(k) = kp;
	} else {
	    IPIV(k) = -kp;
	    IPIV(k - 1) = -kp;
	}

/*        Decrease K and return to the start of the main loop */

	k -= kstep;
	goto L10;

    } else {

/*        Factorize A as L*D*L' using the lower triangle of A   

          K is the main loop index, increasing from 1 to N in steps of
   
          1 or 2 */

	k = 1;
L20:

/*        If K > N, exit from loop */

	if (k > *n) {
	    goto L30;
	}
	kstep = 1;

/*        Determine rows and columns to be interchanged and whether   
          a 1-by-1 or 2-by-2 pivot block will be used */

	absakk = (d__1 = A(k,k), abs(d__1));

/*        IMAX is the row-index of the largest off-diagonal element in
   
          column K, and COLMAX is its absolute value */

	if (k < *n) {
	    i__1 = *n - k;
	    imax = k + idamax_(&i__1, &A(k+1,k), &c__1);
	    colmax = (d__1 = A(imax,k), abs(d__1));
	} else {
	    colmax = 0.;
	}

	if (max(absakk,colmax) == 0.) {

/*           Column K is zero: set INFO and continue */

	    if (*info == 0) {
		*info = k;
	    }
	    kp = k;
	} else {
	    if (absakk >= alpha * colmax) {

/*              no interchange, use 1-by-1 pivot block */

		kp = k;
	    } else {

/*              JMAX is the column-index of the largest off-di
agonal   
                element in row IMAX, and ROWMAX is its absolut
e value */

		i__1 = imax - k;
		jmax = k - 1 + idamax_(&i__1, &A(imax,k), lda);
		rowmax = (d__1 = A(imax,jmax), abs(d__1));
		if (imax < *n) {
		    i__1 = *n - imax;
		    jmax = imax + idamax_(&i__1, &A(imax+1,imax),
			     &c__1);
/* Computing MAX */
		    d__2 = rowmax, d__3 = (d__1 = A(jmax,imax), 
			    abs(d__1));
		    rowmax = max(d__2,d__3);
		}

		if (absakk >= alpha * colmax * (colmax / rowmax)) {

/*                 no interchange, use 1-by-1 pivot block 
*/

		    kp = k;
		} else if ((d__1 = A(imax,imax), abs(d__1)) >= 
			alpha * rowmax) {

/*                 interchange rows and columns K and IMAX
, use 1-by-1   
                   pivot block */

		    kp = imax;
		} else {

/*                 interchange rows and columns K+1 and IM
AX, use 2-by-2   
                   pivot block */

		    kp = imax;
		    kstep = 2;
		}
	    }

	    kk = k + kstep - 1;
	    if (kp != kk) {

/*              Interchange rows and columns KK and KP in the 
trailing   
                submatrix A(k:n,k:n) */

		if (kp < *n) {
		    i__1 = *n - kp;
		    dswap_(&i__1, &A(kp+1,kk), &c__1, &A(kp+1,kp), &c__1);
		}
		i__1 = kp - kk - 1;
		dswap_(&i__1, &A(kk+1,kk), &c__1, &A(kp,kk+1), lda);
		t = A(kk,kk);
		A(kk,kk) = A(kp,kp);
		A(kp,kp) = t;
		if (kstep == 2) {
		    t = A(k+1,k);
		    A(k+1,k) = A(kp,k);
		    A(kp,k) = t;
		}
	    }

/*           Update the trailing submatrix */

	    if (kstep == 1) {

/*              1-by-1 pivot block D(k): column k now holds   

                W(k) = L(k)*D(k)   

                where L(k) is the k-th column of L */

		if (k < *n) {

/*                 Perform a rank-1 update of A(k+1:n,k+1:
n) as   

                   A := A - L(k)*D(k)*L(k)' = A - W(k)*(1/
D(k))*W(k)' */

		    r1 = 1. / A(k,k);
		    i__1 = *n - k;
		    d__1 = -r1;
		    dsyr_(uplo, &i__1, &d__1, &A(k+1,k), &c__1, &
			    A(k+1,k+1), lda);

/*                 Store L(k) in column K */

		    i__1 = *n - k;
		    dscal_(&i__1, &r1, &A(k+1,k), &c__1);
		}
	    } else {

/*              2-by-2 pivot block D(k): columns K and K+1 now
 hold   

                ( W(k) W(k+1) ) = ( L(k) L(k+1) )*D(k)   

                where L(k) and L(k+1) are the k-th and (k+1)-t
h columns   
                of L */

		if (k < *n - 1) {

/*                 Perform a rank-2 update of A(k+2:n,k+2:
n) as   

                   A := A - ( L(k) L(k+1) )*D(k)*( L(k) L(
k+1) )'   
                      = A - ( W(k) W(k+1) )*inv(D(k))*( W(
k) W(k+1) )'   

                   Convert this to two rank-1 updates by u
sing the eigen-   
                   decomposition of D(k) */

		    dlaev2_(&A(k,k), &A(k+1,k), &A(k+1,k+1), &r1, &r2, &c, &s);
		    r1 = 1. / r1;
		    r2 = 1. / r2;
		    i__1 = *n - k - 1;
		    drot_(&i__1, &A(k+2,k), &c__1, &A(k+2,k+1), &c__1, &c, &s);
		    i__1 = *n - k - 1;
		    d__1 = -r1;
		    dsyr_(uplo, &i__1, &d__1, &A(k+2,k), &c__1, &
			    A(k+2,k+2), lda);
		    i__1 = *n - k - 1;
		    d__1 = -r2;
		    dsyr_(uplo, &i__1, &d__1, &A(k+2,k+1), &
			    c__1, &A(k+2,k+2), lda);

/*                 Store L(k) and L(k+1) in columns k and 
k+1 */

		    i__1 = *n - k - 1;
		    dscal_(&i__1, &r1, &A(k+2,k), &c__1);
		    i__1 = *n - k - 1;
		    dscal_(&i__1, &r2, &A(k+2,k+1), &c__1);
		    i__1 = *n - k - 1;
		    d__1 = -s;
		    drot_(&i__1, &A(k+2,k), &c__1, &A(k+2,k+1), &c__1, &c, &d__1);
		}
	    }
	}

/*        Store details of the interchanges in IPIV */

	if (kstep == 1) {
	    IPIV(k) = kp;
	} else {
	    IPIV(k) = -kp;
	    IPIV(k + 1) = -kp;
	}

/*        Increase K and return to the start of the main loop */

	k += kstep;
	goto L20;

    }

L30:
    return 0;

/*     End of DSYTF2 */

} /* dsytf2_ */

