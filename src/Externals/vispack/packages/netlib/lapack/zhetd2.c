#include "f2c.h"

/* Subroutine */ int zhetd2_(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, doublereal *d, doublereal *e, doublecomplex *tau, 
	integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    ZHETD2 reduces a complex Hermitian matrix A to real symmetric   
    tridiagonal form T by a unitary similarity transformation:   
    Q' * A * Q = T.   

    Arguments   
    =========   

    UPLO    (input) CHARACTER*1   
            Specifies whether the upper or lower triangular part of the   
            Hermitian matrix A is stored:   
            = 'U':  Upper triangular   
            = 'L':  Lower triangular   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    A       (input/output) COMPLEX*16 array, dimension (LDA,N)   
            On entry, the Hermitian matrix A.  If UPLO = 'U', the leading 
  
            n-by-n upper triangular part of A contains the upper   
            triangular part of the matrix A, and the strictly lower   
            triangular part of A is not referenced.  If UPLO = 'L', the   
            leading n-by-n lower triangular part of A contains the lower 
  
            triangular part of the matrix A, and the strictly upper   
            triangular part of A is not referenced.   
            On exit, if UPLO = 'U', the diagonal and first superdiagonal 
  
            of A are overwritten by the corresponding elements of the   
            tridiagonal matrix T, and the elements above the first   
            superdiagonal, with the array TAU, represent the unitary   
            matrix Q as a product of elementary reflectors; if UPLO   
            = 'L', the diagonal and first subdiagonal of A are over-   
            written by the corresponding elements of the tridiagonal   
            matrix T, and the elements below the first subdiagonal, with 
  
            the array TAU, represent the unitary matrix Q as a product   
            of elementary reflectors. See Further Details.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(1,N).   

    D       (output) DOUBLE PRECISION array, dimension (N)   
            The diagonal elements of the tridiagonal matrix T:   
            D(i) = A(i,i).   

    E       (output) DOUBLE PRECISION array, dimension (N-1)   
            The off-diagonal elements of the tridiagonal matrix T:   
            E(i) = A(i,i+1) if UPLO = 'U', E(i) = A(i+1,i) if UPLO = 'L'. 
  

    TAU     (output) COMPLEX*16 array, dimension (N-1)   
            The scalar factors of the elementary reflectors (see Further 
  
            Details).   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value.   

    Further Details   
    ===============   

    If UPLO = 'U', the matrix Q is represented as a product of elementary 
  
    reflectors   

       Q = H(n-1) . . . H(2) H(1).   

    Each H(i) has the form   

       H(i) = I - tau * v * v'   

    where tau is a complex scalar, and v is a complex vector with   
    v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in   
    A(1:i-1,i+1), and tau in TAU(i).   

    If UPLO = 'L', the matrix Q is represented as a product of elementary 
  
    reflectors   

       Q = H(1) H(2) . . . H(n-1).   

    Each H(i) has the form   

       H(i) = I - tau * v * v'   

    where tau is a complex scalar, and v is a complex vector with   
    v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in A(i+2:n,i), 
  
    and tau in TAU(i).   

    The contents of A on exit are illustrated by the following examples   
    with n = 5:   

    if UPLO = 'U':                       if UPLO = 'L':   

      (  d   e   v2  v3  v4 )              (  d                  )   
      (      d   e   v3  v4 )              (  e   d              )   
      (          d   e   v4 )              (  v1  e   d          )   
      (              d   e  )              (  v1  v2  e   d      )   
      (                  d  )              (  v1  v2  v3  e   d  )   

    where d and e denote diagonal and off-diagonal elements of T, and vi 
  
    denotes an element of the vector defining H(i).   

    ===================================================================== 
  


       Test the input parameters   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static doublecomplex c_b2 = {0.,0.};
    static integer c__1 = 1;
    
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublereal d__1;
    doublecomplex z__1, z__2, z__3, z__4;
    /* Local variables */
    static doublecomplex taui;
    extern /* Subroutine */ int zher2_(char *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    static integer i;
    static doublecomplex alpha;
    extern logical lsame_(char *, char *);
    extern /* Double Complex */ VOID zdotc_(doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    extern /* Subroutine */ int zhemv_(char *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, doublecomplex *, integer *);
    static logical upper;
    extern /* Subroutine */ int zaxpy_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *), xerbla_(
	    char *, integer *), zlarfg_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *);



#define D(I) d[(I)-1]
#define E(I) e[(I)-1]
#define TAU(I) tau[(I)-1]

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
	xerbla_("ZHETD2", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*n <= 0) {
	return 0;
    }

    if (upper) {

/*        Reduce the upper triangle of A */

	i__1 = *n + *n * a_dim1;
	i__2 = *n + *n * a_dim1;
	d__1 = A(*n,*n).r;
	A(*n,*n).r = d__1, A(*n,*n).i = 0.;
	for (i = *n - 1; i >= 1; --i) {

/*           Generate elementary reflector H(i) = I - tau * v * v'
   
             to annihilate A(1:i-1,i+1) */

	    i__1 = i + (i + 1) * a_dim1;
	    alpha.r = A(i,i+1).r, alpha.i = A(i,i+1).i;
	    zlarfg_(&i, &alpha, &A(1,i+1), &c__1, &taui);
	    i__1 = i;
	    E(i) = alpha.r;

	    if (taui.r != 0. || taui.i != 0.) {

/*              Apply H(i) from both sides to A(1:i,1:i) */

		i__1 = i + (i + 1) * a_dim1;
		A(i,i+1).r = 1., A(i,i+1).i = 0.;

/*              Compute  x := tau * A * v  storing x in TAU(1:
i) */

		zhemv_(uplo, &i, &taui, &A(1,1), lda, &A(1,i+1), &c__1, &c_b2, &TAU(1), &c__1);

/*              Compute  w := x - 1/2 * tau * (x'*v) * v */

		z__3.r = -.5, z__3.i = 0.;
		z__2.r = z__3.r * taui.r - z__3.i * taui.i, z__2.i = z__3.r * 
			taui.i + z__3.i * taui.r;
		zdotc_(&z__4, &i, &TAU(1), &c__1, &A(1,i+1), &
			c__1);
		z__1.r = z__2.r * z__4.r - z__2.i * z__4.i, z__1.i = z__2.r * 
			z__4.i + z__2.i * z__4.r;
		alpha.r = z__1.r, alpha.i = z__1.i;
		zaxpy_(&i, &alpha, &A(1,i+1), &c__1, &TAU(1), &
			c__1);

/*              Apply the transformation as a rank-2 update: 
  
                   A := A - v * w' - w * v' */

		z__1.r = -1., z__1.i = 0.;
		zher2_(uplo, &i, &z__1, &A(1,i+1), &c__1, &TAU(
			1), &c__1, &A(1,1), lda);

	    }
	    i__1 = i + (i + 1) * a_dim1;
	    i__2 = i;
	    A(i,i+1).r = E(i), A(i,i+1).i = 0.;
	    i__1 = i + 1;
	    i__2 = i + 1 + (i + 1) * a_dim1;
	    D(i+1) = A(i+1,i+1).r;
	    i__1 = i;
	    TAU(i).r = taui.r, TAU(i).i = taui.i;
/* L10: */
	}
	i__1 = a_dim1 + 1;
	D(1) = A(1,1).r;
    } else {

/*        Reduce the lower triangle of A */

	i__1 = a_dim1 + 1;
	i__2 = a_dim1 + 1;
	d__1 = A(1,1).r;
	A(1,1).r = d__1, A(1,1).i = 0.;
	i__1 = *n - 1;
	for (i = 1; i <= *n-1; ++i) {

/*           Generate elementary reflector H(i) = I - tau * v * v'
   
             to annihilate A(i+2:n,i) */

	    i__2 = i + 1 + i * a_dim1;
	    alpha.r = A(i+1,i).r, alpha.i = A(i+1,i).i;
	    i__2 = *n - i;
/* Computing MIN */
	    i__3 = i + 2;
	    zlarfg_(&i__2, &alpha, &A(min(i+2,*n),i), &c__1, &
		    taui);
	    i__2 = i;
	    E(i) = alpha.r;

	    if (taui.r != 0. || taui.i != 0.) {

/*              Apply H(i) from both sides to A(i+1:n,i+1:n) 
*/

		i__2 = i + 1 + i * a_dim1;
		A(i+1,i).r = 1., A(i+1,i).i = 0.;

/*              Compute  x := tau * A * v  storing y in TAU(i:
n-1) */

		i__2 = *n - i;
		zhemv_(uplo, &i__2, &taui, &A(i+1,i+1), lda, 
			&A(i+1,i), &c__1, &c_b2, &TAU(i), &c__1);

/*              Compute  w := x - 1/2 * tau * (x'*v) * v */

		z__3.r = -.5, z__3.i = 0.;
		z__2.r = z__3.r * taui.r - z__3.i * taui.i, z__2.i = z__3.r * 
			taui.i + z__3.i * taui.r;
		i__2 = *n - i;
		zdotc_(&z__4, &i__2, &TAU(i), &c__1, &A(i+1,i), &
			c__1);
		z__1.r = z__2.r * z__4.r - z__2.i * z__4.i, z__1.i = z__2.r * 
			z__4.i + z__2.i * z__4.r;
		alpha.r = z__1.r, alpha.i = z__1.i;
		i__2 = *n - i;
		zaxpy_(&i__2, &alpha, &A(i+1,i), &c__1, &TAU(i), 
			&c__1);

/*              Apply the transformation as a rank-2 update: 
  
                   A := A - v * w' - w * v' */

		i__2 = *n - i;
		z__1.r = -1., z__1.i = 0.;
		zher2_(uplo, &i__2, &z__1, &A(i+1,i), &c__1, &
			TAU(i), &c__1, &A(i+1,i+1), lda);

	    }
	    i__2 = i + 1 + i * a_dim1;
	    i__3 = i;
	    A(i+1,i).r = E(i), A(i+1,i).i = 0.;
	    i__2 = i;
	    i__3 = i + i * a_dim1;
	    D(i) = A(i,i).r;
	    i__2 = i;
	    TAU(i).r = taui.r, TAU(i).i = taui.i;
/* L20: */
	}
	i__1 = *n;
	i__2 = *n + *n * a_dim1;
	D(*n) = A(*n,*n).r;
    }

    return 0;

/*     End of ZHETD2 */

} /* zhetd2_ */

