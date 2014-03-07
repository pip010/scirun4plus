#include "f2c.h"

/* Subroutine */ int zlatrd_(char *uplo, integer *n, integer *nb, 
	doublecomplex *a, integer *lda, doublereal *e, doublecomplex *tau, 
	doublecomplex *w, integer *ldw)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    ZLATRD reduces NB rows and columns of a complex Hermitian matrix A to 
  
    Hermitian tridiagonal form by a unitary similarity   
    transformation Q' * A * Q, and returns the matrices V and W which are 
  
    needed to apply the transformation to the unreduced part of A.   

    If UPLO = 'U', ZLATRD reduces the last NB rows and columns of a   
    matrix, of which the upper triangle is supplied;   
    if UPLO = 'L', ZLATRD reduces the first NB rows and columns of a   
    matrix, of which the lower triangle is supplied.   

    This is an auxiliary routine called by ZHETRD.   

    Arguments   
    =========   

    UPLO    (input) CHARACTER   
            Specifies whether the upper or lower triangular part of the   
            Hermitian matrix A is stored:   
            = 'U': Upper triangular   
            = 'L': Lower triangular   

    N       (input) INTEGER   
            The order of the matrix A.   

    NB      (input) INTEGER   
            The number of rows and columns to be reduced.   

    A       (input/output) COMPLEX*16 array, dimension (LDA,N)   
            On entry, the Hermitian matrix A.  If UPLO = 'U', the leading 
  
            n-by-n upper triangular part of A contains the upper   
            triangular part of the matrix A, and the strictly lower   
            triangular part of A is not referenced.  If UPLO = 'L', the   
            leading n-by-n lower triangular part of A contains the lower 
  
            triangular part of the matrix A, and the strictly upper   
            triangular part of A is not referenced.   
            On exit:   
            if UPLO = 'U', the last NB columns have been reduced to   
              tridiagonal form, with the diagonal elements overwriting   
              the diagonal elements of A; the elements above the diagonal 
  
              with the array TAU, represent the unitary matrix Q as a   
              product of elementary reflectors;   
            if UPLO = 'L', the first NB columns have been reduced to   
              tridiagonal form, with the diagonal elements overwriting   
              the diagonal elements of A; the elements below the diagonal 
  
              with the array TAU, represent the  unitary matrix Q as a   
              product of elementary reflectors.   
            See Further Details.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(1,N).   

    E       (output) DOUBLE PRECISION array, dimension (N-1)   
            If UPLO = 'U', E(n-nb:n-1) contains the superdiagonal   
            elements of the last NB columns of the reduced matrix;   
            if UPLO = 'L', E(1:nb) contains the subdiagonal elements of   
            the first NB columns of the reduced matrix.   

    TAU     (output) COMPLEX*16 array, dimension (N-1)   
            The scalar factors of the elementary reflectors, stored in   
            TAU(n-nb:n-1) if UPLO = 'U', and in TAU(1:nb) if UPLO = 'L'. 
  
            See Further Details.   

    W       (output) COMPLEX*16 array, dimension (LDW,NB)   
            The n-by-nb matrix W required to update the unreduced part   
            of A.   

    LDW     (input) INTEGER   
            The leading dimension of the array W. LDW >= max(1,N).   

    Further Details   
    ===============   

    If UPLO = 'U', the matrix Q is represented as a product of elementary 
  
    reflectors   

       Q = H(n) H(n-1) . . . H(n-nb+1).   

    Each H(i) has the form   

       H(i) = I - tau * v * v'   

    where tau is a complex scalar, and v is a complex vector with   
    v(i:n) = 0 and v(i-1) = 1; v(1:i-1) is stored on exit in A(1:i-1,i), 
  
    and tau in TAU(i-1).   

    If UPLO = 'L', the matrix Q is represented as a product of elementary 
  
    reflectors   

       Q = H(1) H(2) . . . H(nb).   

    Each H(i) has the form   

       H(i) = I - tau * v * v'   

    where tau is a complex scalar, and v is a complex vector with   
    v(1:i) = 0 and v(i+1) = 1; v(i+1:n) is stored on exit in A(i+1:n,i), 
  
    and tau in TAU(i).   

    The elements of the vectors v together form the n-by-nb matrix V   
    which is needed, with W, to apply the transformation to the unreduced 
  
    part of the matrix, using a Hermitian rank-2k update of the form:   
    A := A - V*W' - W*V'.   

    The contents of A on exit are illustrated by the following examples   
    with n = 5 and nb = 2:   

    if UPLO = 'U':                       if UPLO = 'L':   

      (  a   a   a   v4  v5 )              (  d                  )   
      (      a   a   v4  v5 )              (  1   d              )   
      (          a   1   v5 )              (  v1  1   a          )   
      (              d   1  )              (  v1  v2  a   a      )   
      (                  d  )              (  v1  v2  a   a   a  )   

    where d denotes a diagonal element of the reduced matrix, a denotes   
    an element of the original matrix that is unchanged, and vi denotes   
    an element of the vector defining H(i).   

    ===================================================================== 
  


       Quick return if possible   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static doublecomplex c_b1 = {0.,0.};
    static doublecomplex c_b2 = {1.,0.};
    static integer c__1 = 1;
    
    /* System generated locals */
    integer a_dim1, a_offset, w_dim1, w_offset, i__1, i__2, i__3;
    doublereal d__1;
    doublecomplex z__1, z__2, z__3, z__4;
    /* Local variables */
    static integer i;
    static doublecomplex alpha;
    extern logical lsame_(char *, char *);
    extern /* Subroutine */ int zscal_(integer *, doublecomplex *, 
	    doublecomplex *, integer *);
    extern /* Double Complex */ VOID zdotc_(doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    extern /* Subroutine */ int zgemv_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, doublecomplex *, integer *), 
	    zhemv_(char *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, integer *), zaxpy_(integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *);
    static integer iw;
    extern /* Subroutine */ int zlarfg_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *), zlacgv_(integer *, 
	    doublecomplex *, integer *);



#define E(I) e[(I)-1]
#define TAU(I) tau[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]
#define W(I,J) w[(I)-1 + ((J)-1)* ( *ldw)]

    if (*n <= 0) {
	return 0;
    }

    if (lsame_(uplo, "U")) {

/*        Reduce last NB columns of upper triangle */

	i__1 = *n - *nb + 1;
	for (i = *n; i >= *n-*nb+1; --i) {
	    iw = i - *n + *nb;
	    if (i < *n) {

/*              Update A(1:i,i) */

		i__2 = i + i * a_dim1;
		i__3 = i + i * a_dim1;
		d__1 = A(i,i).r;
		A(i,i).r = d__1, A(i,i).i = 0.;
		i__2 = *n - i;
		zlacgv_(&i__2, &W(i,iw+1), ldw);
		i__2 = *n - i;
		z__1.r = -1., z__1.i = 0.;
		zgemv_("No transpose", &i, &i__2, &z__1, &A(1,i+1), lda, &W(i,iw+1), ldw, &c_b2, &A(1,i), &c__1);
		i__2 = *n - i;
		zlacgv_(&i__2, &W(i,iw+1), ldw);
		i__2 = *n - i;
		zlacgv_(&i__2, &A(i,i+1), lda);
		i__2 = *n - i;
		z__1.r = -1., z__1.i = 0.;
		zgemv_("No transpose", &i, &i__2, &z__1, &W(1,iw+1), ldw, &A(i,i+1), lda, &c_b2, &A(1,i), &c__1);
		i__2 = *n - i;
		zlacgv_(&i__2, &A(i,i+1), lda);
		i__2 = i + i * a_dim1;
		i__3 = i + i * a_dim1;
		d__1 = A(i,i).r;
		A(i,i).r = d__1, A(i,i).i = 0.;
	    }
	    if (i > 1) {

/*              Generate elementary reflector H(i) to annihila
te   
                A(1:i-2,i) */

		i__2 = i - 1 + i * a_dim1;
		alpha.r = A(i-1,i).r, alpha.i = A(i-1,i).i;
		i__2 = i - 1;
		zlarfg_(&i__2, &alpha, &A(1,i), &c__1, &TAU(i - 1))
			;
		i__2 = i - 1;
		E(i-1) = alpha.r;
		i__2 = i - 1 + i * a_dim1;
		A(i-1,i).r = 1., A(i-1,i).i = 0.;

/*              Compute W(1:i-1,i) */

		i__2 = i - 1;
		zhemv_("Upper", &i__2, &c_b2, &A(1,1), lda, &A(1,i), &c__1, &c_b1, &W(1,iw), &c__1);
		if (i < *n) {
		    i__2 = i - 1;
		    i__3 = *n - i;
		    zgemv_("Conjugate transpose", &i__2, &i__3, &c_b2, &W(1,iw+1), ldw, &A(1,i), &
			    c__1, &c_b1, &W(i+1,iw), &c__1);
		    i__2 = i - 1;
		    i__3 = *n - i;
		    z__1.r = -1., z__1.i = 0.;
		    zgemv_("No transpose", &i__2, &i__3, &z__1, &A(1,i+1), lda, &W(i+1,iw), &c__1, 
			    &c_b2, &W(1,iw), &c__1);
		    i__2 = i - 1;
		    i__3 = *n - i;
		    zgemv_("Conjugate transpose", &i__2, &i__3, &c_b2, &A(1,i+1), lda, &A(1,i), &
			    c__1, &c_b1, &W(i+1,iw), &c__1);
		    i__2 = i - 1;
		    i__3 = *n - i;
		    z__1.r = -1., z__1.i = 0.;
		    zgemv_("No transpose", &i__2, &i__3, &z__1, &W(1,iw+1), ldw, &W(i+1,iw), &c__1, 
			    &c_b2, &W(1,iw), &c__1);
		}
		i__2 = i - 1;
		zscal_(&i__2, &TAU(i - 1), &W(1,iw), &c__1);
		z__3.r = -.5, z__3.i = 0.;
		i__2 = i - 1;
		z__2.r = z__3.r * TAU(i-1).r - z__3.i * TAU(i-1).i, z__2.i =
			 z__3.r * TAU(i-1).i + z__3.i * TAU(i-1).r;
		i__3 = i - 1;
		zdotc_(&z__4, &i__3, &W(1,iw), &c__1, &A(1,i), &c__1);
		z__1.r = z__2.r * z__4.r - z__2.i * z__4.i, z__1.i = z__2.r * 
			z__4.i + z__2.i * z__4.r;
		alpha.r = z__1.r, alpha.i = z__1.i;
		i__2 = i - 1;
		zaxpy_(&i__2, &alpha, &A(1,i), &c__1, &W(1,iw), &c__1);
	    }

/* L10: */
	}
    } else {

/*        Reduce first NB columns of lower triangle */

	i__1 = *nb;
	for (i = 1; i <= *nb; ++i) {

/*           Update A(i:n,i) */

	    i__2 = i + i * a_dim1;
	    i__3 = i + i * a_dim1;
	    d__1 = A(i,i).r;
	    A(i,i).r = d__1, A(i,i).i = 0.;
	    i__2 = i - 1;
	    zlacgv_(&i__2, &W(i,1), ldw);
	    i__2 = *n - i + 1;
	    i__3 = i - 1;
	    z__1.r = -1., z__1.i = 0.;
	    zgemv_("No transpose", &i__2, &i__3, &z__1, &A(i,1), lda, &
		    W(i,1), ldw, &c_b2, &A(i,i), &c__1)
		    ;
	    i__2 = i - 1;
	    zlacgv_(&i__2, &W(i,1), ldw);
	    i__2 = i - 1;
	    zlacgv_(&i__2, &A(i,1), lda);
	    i__2 = *n - i + 1;
	    i__3 = i - 1;
	    z__1.r = -1., z__1.i = 0.;
	    zgemv_("No transpose", &i__2, &i__3, &z__1, &W(i,1), ldw, &
		    A(i,1), lda, &c_b2, &A(i,i), &c__1)
		    ;
	    i__2 = i - 1;
	    zlacgv_(&i__2, &A(i,1), lda);
	    i__2 = i + i * a_dim1;
	    i__3 = i + i * a_dim1;
	    d__1 = A(i,i).r;
	    A(i,i).r = d__1, A(i,i).i = 0.;
	    if (i < *n) {

/*              Generate elementary reflector H(i) to annihila
te   
                A(i+2:n,i) */

		i__2 = i + 1 + i * a_dim1;
		alpha.r = A(i+1,i).r, alpha.i = A(i+1,i).i;
		i__2 = *n - i;
/* Computing MIN */
		i__3 = i + 2;
		zlarfg_(&i__2, &alpha, &A(min(i+2,*n),i), &c__1, &
			TAU(i));
		i__2 = i;
		E(i) = alpha.r;
		i__2 = i + 1 + i * a_dim1;
		A(i+1,i).r = 1., A(i+1,i).i = 0.;

/*              Compute W(i+1:n,i) */

		i__2 = *n - i;
		zhemv_("Lower", &i__2, &c_b2, &A(i+1,i+1), 
			lda, &A(i+1,i), &c__1, &c_b1, &W(i+1,i), &c__1);
		i__2 = *n - i;
		i__3 = i - 1;
		zgemv_("Conjugate transpose", &i__2, &i__3, &c_b2, &W(i+1,1), ldw, &A(i+1,i), &c__1, &c_b1, &
			W(1,i), &c__1);
		i__2 = *n - i;
		i__3 = i - 1;
		z__1.r = -1., z__1.i = 0.;
		zgemv_("No transpose", &i__2, &i__3, &z__1, &A(i+1,1)
			, lda, &W(1,i), &c__1, &c_b2, &W(i+1,i), &c__1);
		i__2 = *n - i;
		i__3 = i - 1;
		zgemv_("Conjugate transpose", &i__2, &i__3, &c_b2, &A(i+1,1), lda, &A(i+1,i), &c__1, &c_b1, &
			W(1,i), &c__1);
		i__2 = *n - i;
		i__3 = i - 1;
		z__1.r = -1., z__1.i = 0.;
		zgemv_("No transpose", &i__2, &i__3, &z__1, &W(i+1,1)
			, ldw, &W(1,i), &c__1, &c_b2, &W(i+1,i), &c__1);
		i__2 = *n - i;
		zscal_(&i__2, &TAU(i), &W(i+1,i), &c__1);
		z__3.r = -.5, z__3.i = 0.;
		i__2 = i;
		z__2.r = z__3.r * TAU(i).r - z__3.i * TAU(i).i, z__2.i =
			 z__3.r * TAU(i).i + z__3.i * TAU(i).r;
		i__3 = *n - i;
		zdotc_(&z__4, &i__3, &W(i+1,i), &c__1, &A(i+1,i), &c__1);
		z__1.r = z__2.r * z__4.r - z__2.i * z__4.i, z__1.i = z__2.r * 
			z__4.i + z__2.i * z__4.r;
		alpha.r = z__1.r, alpha.i = z__1.i;
		i__2 = *n - i;
		zaxpy_(&i__2, &alpha, &A(i+1,i), &c__1, &W(i+1,i), &c__1);
	    }

/* L20: */
	}
    }

    return 0;

/*     End of ZLATRD */

} /* zlatrd_ */

