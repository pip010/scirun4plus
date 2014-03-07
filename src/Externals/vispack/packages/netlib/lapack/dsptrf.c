#include "f2c.h"

/* Subroutine */ int dsptrf_(char *uplo, integer *n, doublereal *ap, integer *
	ipiv, integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       March 31, 1993   


    Purpose   
    =======   

    DSPTRF computes the factorization of a real symmetric matrix A stored 
  
    in packed format using the Bunch-Kaufman diagonal pivoting method:   

       A = U*D*U**T  or  A = L*D*L**T   

    where U (or L) is a product of permutation and unit upper (lower)   
    triangular matrices, and D is symmetric and block diagonal with   
    1-by-1 and 2-by-2 diagonal blocks.   

    Arguments   
    =========   

    UPLO    (input) CHARACTER*1   
            = 'U':  Upper triangle of A is stored;   
            = 'L':  Lower triangle of A is stored.   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    AP      (input/output) DOUBLE PRECISION array, dimension (N*(N+1)/2) 
  
            On entry, the upper or lower triangle of the symmetric matrix 
  
            A, packed columnwise in a linear array.  The j-th column of A 
  
            is stored in the array AP as follows:   
            if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;   
            if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.   

            On exit, the block diagonal matrix D and the multipliers used 
  
            to obtain the factor U or L, stored as a packed triangular   
            matrix overwriting A (see below for further details).   

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
            < 0: if INFO = -i, the i-th argument had an illegal value   
            > 0: if INFO = i, D(i,i) is exactly zero.  The factorization 
  
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
    integer i__1;
    doublereal d__1, d__2, d__3;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    static integer imax, jmax;
    extern /* Subroutine */ int drot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *), dspr_(char *
	    , integer *, doublereal *, doublereal *, integer *, doublereal *);
    static doublereal c;
    static integer j, k;
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
    static integer kc, kk, kp;
    static doublereal absakk;
    static integer kx;
    extern integer idamax_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int xerbla_(char *, integer *);
    static doublereal colmax, rowmax;
    static integer knc, kpc, npp;



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
	xerbla_("DSPTRF", &i__1);
	return 0;
    }

/*     Initialize ALPHA for use in choosing pivot block size. */

    alpha = (sqrt(17.) + 1.) / 8.;

    if (upper) {

/*        Factorize A as U*D*U' using the upper triangle of A   

          K is the main loop index, decreasing from N to 1 in steps of
   
          1 or 2 */

	k = *n;
	kc = (*n - 1) * *n / 2 + 1;
L10:
	knc = kc;

/*        If K < 1, exit from loop */

	if (k < 1) {
	    goto L70;
	}
	kstep = 1;

/*        Determine rows and columns to be interchanged and whether   
          a 1-by-1 or 2-by-2 pivot block will be used */

	absakk = (d__1 = AP(kc + k - 1), abs(d__1));

/*        IMAX is the row-index of the largest off-diagonal element in
   
          column K, and COLMAX is its absolute value */

	if (k > 1) {
	    i__1 = k - 1;
	    imax = idamax_(&i__1, &AP(kc), &c__1);
	    colmax = (d__1 = AP(kc + imax - 1), abs(d__1));
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

		rowmax = 0.;
		jmax = imax;
		kx = imax * (imax + 1) / 2 + imax;
		i__1 = k;
		for (j = imax + 1; j <= k; ++j) {
		    if ((d__1 = AP(kx), abs(d__1)) > rowmax) {
			rowmax = (d__1 = AP(kx), abs(d__1));
			jmax = j;
		    }
		    kx += j;
/* L20: */
		}
		kpc = (imax - 1) * imax / 2 + 1;
		if (imax > 1) {
		    i__1 = imax - 1;
		    jmax = idamax_(&i__1, &AP(kpc), &c__1);
/* Computing MAX */
		    d__2 = rowmax, d__3 = (d__1 = AP(kpc + jmax - 1), abs(
			    d__1));
		    rowmax = max(d__2,d__3);
		}

		if (absakk >= alpha * colmax * (colmax / rowmax)) {

/*                 no interchange, use 1-by-1 pivot block 
*/

		    kp = k;
		} else if ((d__1 = AP(kpc + imax - 1), abs(d__1)) >= alpha * 
			rowmax) {

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
	    if (kstep == 2) {
		knc = knc - k + 1;
	    }
	    if (kp != kk) {

/*              Interchange rows and columns KK and KP in the 
leading   
                submatrix A(1:k,1:k) */

		i__1 = kp - 1;
		dswap_(&i__1, &AP(knc), &c__1, &AP(kpc), &c__1);
		kx = kpc + kp - 1;
		i__1 = kk - 1;
		for (j = kp + 1; j <= kk-1; ++j) {
		    kx = kx + j - 1;
		    t = AP(knc + j - 1);
		    AP(knc + j - 1) = AP(kx);
		    AP(kx) = t;
/* L30: */
		}
		t = AP(knc + kk - 1);
		AP(knc + kk - 1) = AP(kpc + kp - 1);
		AP(kpc + kp - 1) = t;
		if (kstep == 2) {
		    t = AP(kc + k - 2);
		    AP(kc + k - 2) = AP(kc + kp - 1);
		    AP(kc + kp - 1) = t;
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

		r1 = 1. / AP(kc + k - 1);
		i__1 = k - 1;
		d__1 = -r1;
		dspr_(uplo, &i__1, &d__1, &AP(kc), &c__1, &AP(1));

/*              Store U(k) in column k */

		i__1 = k - 1;
		dscal_(&i__1, &r1, &AP(kc), &c__1);
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

		dlaev2_(&AP(kc - 1), &AP(kc + k - 2), &AP(kc + k - 1), &r1, &
			r2, &c, &s);
		r1 = 1. / r1;
		r2 = 1. / r2;
		i__1 = k - 2;
		drot_(&i__1, &AP(knc), &c__1, &AP(kc), &c__1, &c, &s);
		i__1 = k - 2;
		d__1 = -r1;
		dspr_(uplo, &i__1, &d__1, &AP(knc), &c__1, &AP(1));
		i__1 = k - 2;
		d__1 = -r2;
		dspr_(uplo, &i__1, &d__1, &AP(kc), &c__1, &AP(1));

/*              Store U(k) and U(k-1) in columns k and k-1 */

		i__1 = k - 2;
		dscal_(&i__1, &r1, &AP(knc), &c__1);
		i__1 = k - 2;
		dscal_(&i__1, &r2, &AP(kc), &c__1);
		i__1 = k - 2;
		d__1 = -s;
		drot_(&i__1, &AP(knc), &c__1, &AP(kc), &c__1, &c, &d__1);
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
	kc = knc - k;
	goto L10;

    } else {

/*        Factorize A as L*D*L' using the lower triangle of A   

          K is the main loop index, increasing from 1 to N in steps of
   
          1 or 2 */

	k = 1;
	kc = 1;
	npp = *n * (*n + 1) / 2;
L40:
	knc = kc;

/*        If K > N, exit from loop */

	if (k > *n) {
	    goto L70;
	}
	kstep = 1;

/*        Determine rows and columns to be interchanged and whether   
          a 1-by-1 or 2-by-2 pivot block will be used */

	absakk = (d__1 = AP(kc), abs(d__1));

/*        IMAX is the row-index of the largest off-diagonal element in
   
          column K, and COLMAX is its absolute value */

	if (k < *n) {
	    i__1 = *n - k;
	    imax = k + idamax_(&i__1, &AP(kc + 1), &c__1);
	    colmax = (d__1 = AP(kc + imax - k), abs(d__1));
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

		rowmax = 0.;
		kx = kc + imax - k;
		i__1 = imax - 1;
		for (j = k; j <= imax-1; ++j) {
		    if ((d__1 = AP(kx), abs(d__1)) > rowmax) {
			rowmax = (d__1 = AP(kx), abs(d__1));
			jmax = j;
		    }
		    kx = kx + *n - j;
/* L50: */
		}
		kpc = npp - (*n - imax + 1) * (*n - imax + 2) / 2 + 1;
		if (imax < *n) {
		    i__1 = *n - imax;
		    jmax = imax + idamax_(&i__1, &AP(kpc + 1), &c__1);
/* Computing MAX */
		    d__2 = rowmax, d__3 = (d__1 = AP(kpc + jmax - imax), abs(
			    d__1));
		    rowmax = max(d__2,d__3);
		}

		if (absakk >= alpha * colmax * (colmax / rowmax)) {

/*                 no interchange, use 1-by-1 pivot block 
*/

		    kp = k;
		} else if ((d__1 = AP(kpc), abs(d__1)) >= alpha * rowmax) {

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
	    if (kstep == 2) {
		knc = knc + *n - k + 1;
	    }
	    if (kp != kk) {

/*              Interchange rows and columns KK and KP in the 
trailing   
                submatrix A(k:n,k:n) */

		if (kp < *n) {
		    i__1 = *n - kp;
		    dswap_(&i__1, &AP(knc + kp - kk + 1), &c__1, &AP(kpc + 1),
			     &c__1);
		}
		kx = knc + kp - kk;
		i__1 = kp - 1;
		for (j = kk + 1; j <= kp-1; ++j) {
		    kx = kx + *n - j + 1;
		    t = AP(knc + j - kk);
		    AP(knc + j - kk) = AP(kx);
		    AP(kx) = t;
/* L60: */
		}
		t = AP(knc);
		AP(knc) = AP(kpc);
		AP(kpc) = t;
		if (kstep == 2) {
		    t = AP(kc + 1);
		    AP(kc + 1) = AP(kc + kp - k);
		    AP(kc + kp - k) = t;
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

		    r1 = 1. / AP(kc);
		    i__1 = *n - k;
		    d__1 = -r1;
		    dspr_(uplo, &i__1, &d__1, &AP(kc + 1), &c__1, &AP(kc + *n 
			    - k + 1));

/*                 Store L(k) in column K */

		    i__1 = *n - k;
		    dscal_(&i__1, &r1, &AP(kc + 1), &c__1);
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

		    dlaev2_(&AP(kc), &AP(kc + 1), &AP(knc), &r1, &r2, &c, &s);
		    r1 = 1. / r1;
		    r2 = 1. / r2;
		    i__1 = *n - k - 1;
		    drot_(&i__1, &AP(kc + 2), &c__1, &AP(knc + 1), &c__1, &c, 
			    &s);
		    i__1 = *n - k - 1;
		    d__1 = -r1;
		    dspr_(uplo, &i__1, &d__1, &AP(kc + 2), &c__1, &AP(knc + *
			    n - k));
		    i__1 = *n - k - 1;
		    d__1 = -r2;
		    dspr_(uplo, &i__1, &d__1, &AP(knc + 1), &c__1, &AP(knc + *
			    n - k));

/*                 Store L(k) and L(k+1) in columns k and 
k+1 */

		    i__1 = *n - k - 1;
		    dscal_(&i__1, &r1, &AP(kc + 2), &c__1);
		    i__1 = *n - k - 1;
		    dscal_(&i__1, &r2, &AP(knc + 1), &c__1);
		    i__1 = *n - k - 1;
		    d__1 = -s;
		    drot_(&i__1, &AP(kc + 2), &c__1, &AP(knc + 1), &c__1, &c, 
			    &d__1);
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
	kc = knc + *n - k + 2;
	goto L40;

    }

L70:
    return 0;

/*     End of DSPTRF */

} /* dsptrf_ */

