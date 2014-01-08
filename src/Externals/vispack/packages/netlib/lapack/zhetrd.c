#include "f2c.h"

/* Subroutine */ int zhetrd_(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, doublereal *d, doublereal *e, doublecomplex *tau, 
	doublecomplex *work, integer *lwork, integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    ZHETRD reduces a complex Hermitian matrix A to real symmetric   
    tridiagonal form T by a unitary similarity transformation:   
    Q**H * A * Q = T.   

    Arguments   
    =========   

    UPLO    (input) CHARACTER*1   
            = 'U':  Upper triangle of A is stored;   
            = 'L':  Lower triangle of A is stored.   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    A       (input/output) COMPLEX*16 array, dimension (LDA,N)   
            On entry, the Hermitian matrix A.  If UPLO = 'U', the leading 
  
            N-by-N upper triangular part of A contains the upper   
            triangular part of the matrix A, and the strictly lower   
            triangular part of A is not referenced.  If UPLO = 'L', the   
            leading N-by-N lower triangular part of A contains the lower 
  
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

    WORK    (workspace/output) COMPLEX*16 array, dimension (LWORK)   
            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INTEGER   
            The dimension of the array WORK.  LWORK >= 1.   
            For optimum performance LWORK >= N*NB, where NB is the   
            optimal blocksize.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   

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
    static integer c__1 = 1;
    static integer c_n1 = -1;
    static integer c__3 = 3;
    static integer c__2 = 2;
    static doublereal c_b23 = 1.;
    
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5;
    doublecomplex z__1;
    /* Local variables */
    static integer i, j;
    extern logical lsame_(char *, char *);
    static integer nbmin, iinfo;
    static logical upper;
    extern /* Subroutine */ int zhetd2_(char *, integer *, doublecomplex *, 
	    integer *, doublereal *, doublereal *, doublecomplex *, integer *), zher2k_(char *, char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublereal *, doublecomplex *, integer *);
    static integer nb, kk, nx;
    extern /* Subroutine */ int xerbla_(char *, integer *);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int zlatrd_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublereal *, doublecomplex *, 
	    doublecomplex *, integer *);
    static integer ldwork, iws;



#define D(I) d[(I)-1]
#define E(I) e[(I)-1]
#define TAU(I) tau[(I)-1]
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
    } else if (*lwork < 1) {
	*info = -9;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("ZHETRD", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0) {
	WORK(1).r = 1., WORK(1).i = 0.;
	return 0;
    }

/*     Determine the block size. */

    nb = ilaenv_(&c__1, "ZHETRD", uplo, n, &c_n1, &c_n1, &c_n1, 6L, 1L);
    nx = *n;
    iws = 1;
    if (nb > 1 && nb < *n) {

/*        Determine when to cross over from blocked to unblocked code 
  
          (last block is always handled by unblocked code).   

   Computing MAX */
	i__1 = nb, i__2 = ilaenv_(&c__3, "ZHETRD", uplo, n, &c_n1, &c_n1, &
		c_n1, 6L, 1L);
	nx = max(i__1,i__2);
	if (nx < *n) {

/*           Determine if workspace is large enough for blocked co
de. */

	    ldwork = *n;
	    iws = ldwork * nb;
	    if (*lwork < iws) {

/*              Not enough workspace to use optimal NB:  deter
mine the   
                minimum value of NB, and reduce NB or force us
e of   
                unblocked code by setting NX = N.   

   Computing MAX */
		i__1 = *lwork / ldwork;
		nb = max(i__1,1);
		nbmin = ilaenv_(&c__2, "ZHETRD", uplo, n, &c_n1, &c_n1, &c_n1,
			 6L, 1L);
		if (nb < nbmin) {
		    nx = *n;
		}
	    }
	} else {
	    nx = *n;
	}
    } else {
	nb = 1;
    }

    if (upper) {

/*        Reduce the upper triangle of A.   
          Columns 1:kk are handled by the unblocked method. */

	kk = *n - (*n - nx + nb - 1) / nb * nb;
	i__1 = kk + 1;
	i__2 = -nb;
	for (i = *n - nb + 1; -nb < 0 ? i >= kk+1 : i <= kk+1; i += -nb) {

/*           Reduce columns i:i+nb-1 to tridiagonal form and form 
the   
             matrix W which is needed to update the unreduced part
 of   
             the matrix */

	    i__3 = i + nb - 1;
	    zlatrd_(uplo, &i__3, &nb, &A(1,1), lda, &E(1), &TAU(1), &
		    WORK(1), &ldwork);

/*           Update the unreduced submatrix A(1:i-1,1:i-1), using 
an   
             update of the form:  A := A - V*W' - W*V' */

	    i__3 = i - 1;
	    z__1.r = -1., z__1.i = 0.;
	    zher2k_(uplo, "No transpose", &i__3, &nb, &z__1, &A(1,i), lda, &WORK(1), &ldwork, &c_b23, &A(1,1), lda);

/*           Copy superdiagonal elements back into A, and diagonal
   
             elements into D */

	    i__3 = i + nb - 1;
	    for (j = i; j <= i+nb-1; ++j) {
		i__4 = j - 1 + j * a_dim1;
		i__5 = j - 1;
		A(j-1,j).r = E(j-1), A(j-1,j).i = 0.;
		i__4 = j;
		i__5 = j + j * a_dim1;
		D(j) = A(j,j).r;
/* L10: */
	    }
/* L20: */
	}

/*        Use unblocked code to reduce the last or only block */

	zhetd2_(uplo, &kk, &A(1,1), lda, &D(1), &E(1), &TAU(1), &iinfo);
    } else {

/*        Reduce the lower triangle of A */

	i__2 = *n - nx;
	i__1 = nb;
	for (i = 1; nb < 0 ? i >= *n-nx : i <= *n-nx; i += nb) {

/*           Reduce columns i:i+nb-1 to tridiagonal form and form 
the   
             matrix W which is needed to update the unreduced part
 of   
             the matrix */

	    i__3 = *n - i + 1;
	    zlatrd_(uplo, &i__3, &nb, &A(i,i), lda, &E(i), &TAU(i),
		     &WORK(1), &ldwork);

/*           Update the unreduced submatrix A(i+nb:n,i+nb:n), usin
g   
             an update of the form:  A := A - V*W' - W*V' */

	    i__3 = *n - i - nb + 1;
	    z__1.r = -1., z__1.i = 0.;
	    zher2k_(uplo, "No transpose", &i__3, &nb, &z__1, &A(i+nb,i), lda, &WORK(nb + 1), &ldwork, &c_b23, &A(i+nb,i+nb), lda);

/*           Copy subdiagonal elements back into A, and diagonal 
  
             elements into D */

	    i__3 = i + nb - 1;
	    for (j = i; j <= i+nb-1; ++j) {
		i__4 = j + 1 + j * a_dim1;
		i__5 = j;
		A(j+1,j).r = E(j), A(j+1,j).i = 0.;
		i__4 = j;
		i__5 = j + j * a_dim1;
		D(j) = A(j,j).r;
/* L30: */
	    }
/* L40: */
	}

/*        Use unblocked code to reduce the last or only block */

	i__1 = *n - i + 1;
	zhetd2_(uplo, &i__1, &A(i,i), lda, &D(i), &E(i), &TAU(i), &
		iinfo);
    }

    WORK(1).r = (doublereal) iws, WORK(1).i = 0.;
    return 0;

/*     End of ZHETRD */

} /* zhetrd_ */

