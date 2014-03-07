#include "f2c.h"

/* Subroutine */ int clahrd_(integer *n, integer *k, integer *nb, complex *a, 
	integer *lda, complex *tau, complex *t, integer *ldt, complex *y, 
	integer *ldy)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    CLAHRD reduces the first NB columns of a complex general n-by-(n-k+1) 
  
    matrix A so that elements below the k-th subdiagonal are zero. The   
    reduction is performed by a unitary similarity transformation   
    Q' * A * Q. The routine returns the matrices V and T which determine 
  
    Q as a block reflector I - V*T*V', and also the matrix Y = A * V * T. 
  

    This is an auxiliary routine called by CGEHRD.   

    Arguments   
    =========   

    N       (input) INTEGER   
            The order of the matrix A.   

    K       (input) INTEGER   
            The offset for the reduction. Elements below the k-th   
            subdiagonal in the first NB columns are reduced to zero.   

    NB      (input) INTEGER   
            The number of columns to be reduced.   

    A       (input/output) COMPLEX array, dimension (LDA,N-K+1)   
            On entry, the n-by-(n-k+1) general matrix A.   
            On exit, the elements on and above the k-th subdiagonal in   
            the first NB columns are overwritten with the corresponding   
            elements of the reduced matrix; the elements below the k-th   
            subdiagonal, with the array TAU, represent the matrix Q as a 
  
            product of elementary reflectors. The other columns of A are 
  
            unchanged. See Further Details.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(1,N).   

    TAU     (output) COMPLEX array, dimension (NB)   
            The scalar factors of the elementary reflectors. See Further 
  
            Details.   

    T       (output) COMPLEX array, dimension (NB,NB)   
            The upper triangular matrix T.   

    LDT     (input) INTEGER   
            The leading dimension of the array T.  LDT >= NB.   

    Y       (output) COMPLEX array, dimension (LDY,NB)   
            The n-by-nb matrix Y.   

    LDY     (input) INTEGER   
            The leading dimension of the array Y. LDY >= max(1,N).   

    Further Details   
    ===============   

    The matrix Q is represented as a product of nb elementary reflectors 
  

       Q = H(1) H(2) . . . H(nb).   

    Each H(i) has the form   

       H(i) = I - tau * v * v'   

    where tau is a complex scalar, and v is a complex vector with   
    v(1:i+k-1) = 0, v(i+k) = 1; v(i+k+1:n) is stored on exit in   
    A(i+k+1:n,i), and tau in TAU(i).   

    The elements of the vectors v together form the (n-k+1)-by-nb matrix 
  
    V which is needed, with T and Y, to apply the transformation to the   
    unreduced part of the matrix, using an update of the form:   
    A := (I - V*T*V') * (A - Y*V').   

    The contents of A on exit are illustrated by the following example   
    with n = 7, k = 3 and nb = 2:   

       ( a   h   a   a   a )   
       ( a   h   a   a   a )   
       ( a   h   a   a   a )   
       ( h   h   a   a   a )   
       ( v1  h   a   a   a )   
       ( v1  v2  a   a   a )   
       ( v1  v2  a   a   a )   

    where a denotes an element of the original matrix A, h denotes a   
    modified element of the upper Hessenberg matrix H, and vi denotes an 
  
    element of the vector defining H(i).   

    ===================================================================== 
  


       Quick return if possible   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static complex c_b1 = {0.f,0.f};
    static complex c_b2 = {1.f,0.f};
    static integer c__1 = 1;
    
    /* System generated locals */
    integer a_dim1, a_offset, t_dim1, t_offset, y_dim1, y_offset, i__1, i__2, 
	    i__3;
    complex q__1;
    /* Local variables */
    static integer i;
    extern /* Subroutine */ int cscal_(integer *, complex *, complex *, 
	    integer *), cgemv_(char *, integer *, integer *, complex *, 
	    complex *, integer *, complex *, integer *, complex *, complex *, 
	    integer *), ccopy_(integer *, complex *, integer *, 
	    complex *, integer *), caxpy_(integer *, complex *, complex *, 
	    integer *, complex *, integer *), ctrmv_(char *, char *, char *, 
	    integer *, complex *, integer *, complex *, integer *);
    static complex ei;
    extern /* Subroutine */ int clarfg_(integer *, complex *, complex *, 
	    integer *, complex *), clacgv_(integer *, complex *, integer *);



#define TAU(I) tau[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]
#define T(I,J) t[(I)-1 + ((J)-1)* ( *ldt)]
#define Y(I,J) y[(I)-1 + ((J)-1)* ( *ldy)]

    if (*n <= 1) {
	return 0;
    }

    i__1 = *nb;
    for (i = 1; i <= *nb; ++i) {
	if (i > 1) {

/*           Update A(1:n,i)   

             Compute i-th column of A - Y * V' */

	    i__2 = i - 1;
	    clacgv_(&i__2, &A(*k+i-1,1), lda);
	    i__2 = i - 1;
	    q__1.r = -1.f, q__1.i = 0.f;
	    cgemv_("No transpose", n, &i__2, &q__1, &Y(1,1), ldy, &A(*k+i-1,1), lda, &c_b2, &A(1,i), &c__1);
	    i__2 = i - 1;
	    clacgv_(&i__2, &A(*k+i-1,1), lda);

/*           Apply I - V * T' * V' to this column (call it b) from
 the   
             left, using the last column of T as workspace   

             Let  V = ( V1 )   and   b = ( b1 )   (first I-1 rows)
   
                      ( V2 )             ( b2 )   

             where V1 is unit lower triangular   

             w := V1' * b1 */

	    i__2 = i - 1;
	    ccopy_(&i__2, &A(*k+1,i), &c__1, &T(1,*nb)
		    , &c__1);
	    i__2 = i - 1;
	    ctrmv_("Lower", "Conjugate transpose", "Unit", &i__2, &A(*k+1,1), lda, &T(1,*nb), &c__1);

/*           w := w + V2'*b2 */

	    i__2 = *n - *k - i + 1;
	    i__3 = i - 1;
	    cgemv_("Conjugate transpose", &i__2, &i__3, &c_b2, &A(*k+i,1), lda, &A(*k+i,i), &c__1, &c_b2, &T(1,*nb), &c__1);

/*           w := T'*w */

	    i__2 = i - 1;
	    ctrmv_("Upper", "Conjugate transpose", "Non-unit", &i__2, &T(1,1), ldt, &T(1,*nb), &c__1);

/*           b2 := b2 - V2*w */

	    i__2 = *n - *k - i + 1;
	    i__3 = i - 1;
	    q__1.r = -1.f, q__1.i = 0.f;
	    cgemv_("No transpose", &i__2, &i__3, &q__1, &A(*k+i,1), 
		    lda, &T(1,*nb), &c__1, &c_b2, &A(*k+i,i), &c__1);

/*           b1 := b1 - V1*w */

	    i__2 = i - 1;
	    ctrmv_("Lower", "No transpose", "Unit", &i__2, &A(*k+1,1)
		    , lda, &T(1,*nb), &c__1);
	    i__2 = i - 1;
	    q__1.r = -1.f, q__1.i = 0.f;
	    caxpy_(&i__2, &q__1, &T(1,*nb), &c__1, &A(*k+1,i), &c__1);

	    i__2 = *k + i - 1 + (i - 1) * a_dim1;
	    A(*k+i-1,i-1).r = ei.r, A(*k+i-1,i-1).i = ei.i;
	}

/*        Generate the elementary reflector H(i) to annihilate   
          A(k+i+1:n,i) */

	i__2 = *k + i + i * a_dim1;
	ei.r = A(*k+i,i).r, ei.i = A(*k+i,i).i;
	i__2 = *n - *k - i + 1;
/* Computing MIN */
	i__3 = *k + i + 1;
	clarfg_(&i__2, &ei, &A(min(*k+i+1,*n),i), &c__1, &TAU(i));
	i__2 = *k + i + i * a_dim1;
	A(*k+i,i).r = 1.f, A(*k+i,i).i = 0.f;

/*        Compute  Y(1:n,i) */

	i__2 = *n - *k - i + 1;
	cgemv_("No transpose", n, &i__2, &c_b2, &A(1,i+1), lda,
		 &A(*k+i,i), &c__1, &c_b1, &Y(1,i), &
		c__1);
	i__2 = *n - *k - i + 1;
	i__3 = i - 1;
	cgemv_("Conjugate transpose", &i__2, &i__3, &c_b2, &A(*k+i,1)
		, lda, &A(*k+i,i), &c__1, &c_b1, &T(1,i), &c__1);
	i__2 = i - 1;
	q__1.r = -1.f, q__1.i = 0.f;
	cgemv_("No transpose", n, &i__2, &q__1, &Y(1,1), ldy, &T(1,i), &c__1, &c_b2, &Y(1,i), &c__1);
	cscal_(n, &TAU(i), &Y(1,i), &c__1);

/*        Compute T(1:i,i) */

	i__2 = i - 1;
	i__3 = i;
	q__1.r = -(doublereal)TAU(i).r, q__1.i = -(doublereal)TAU(i).i;
	cscal_(&i__2, &q__1, &T(1,i), &c__1);
	i__2 = i - 1;
	ctrmv_("Upper", "No transpose", "Non-unit", &i__2, &T(1,1), ldt, 
		&T(1,i), &c__1);
	i__2 = i + i * t_dim1;
	i__3 = i;
	T(i,i).r = TAU(i).r, T(i,i).i = TAU(i).i;

/* L10: */
    }
    i__1 = *k + *nb + *nb * a_dim1;
    A(*k+*nb,*nb).r = ei.r, A(*k+*nb,*nb).i = ei.i;

    return 0;

/*     End of CLAHRD */

} /* clahrd_ */

