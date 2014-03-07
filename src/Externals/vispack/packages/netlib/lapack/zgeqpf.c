#include "f2c.h"

/* Subroutine */ int zgeqpf_(integer *m, integer *n, doublecomplex *a, 
	integer *lda, integer *jpvt, doublecomplex *tau, doublecomplex *work, 
	doublereal *rwork, integer *info)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       March 31, 1993   


    Purpose   
    =======   

    ZGEQPF computes a QR factorization with column pivoting of a   
    complex M-by-N matrix A: A*P = Q*R.   

    Arguments   
    =========   

    M       (input) INTEGER   
            The number of rows of the matrix A. M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrix A. N >= 0   

    A       (input/output) COMPLEX*16 array, dimension (LDA,N)   
            On entry, the M-by-N matrix A.   
            On exit, the upper triangle of the array contains the   
            min(M,N)-by-N upper triangular matrix R; the elements   
            below the diagonal, together with the array TAU,   
            represent the orthogonal matrix Q as a product of   
            min(m,n) elementary reflectors.   

    LDA     (input) INTEGER   
            The leading dimension of the array A. LDA >= max(1,M).   

    JPVT    (input/output) INTEGER array, dimension (N)   
            On entry, if JPVT(i) .ne. 0, the i-th column of A is permuted 
  
            to the front of A*P (a leading column); if JPVT(i) = 0,   
            the i-th column of A is a free column.   
            On exit, if JPVT(i) = k, then the i-th column of A*P   
            was the k-th column of A.   

    TAU     (output) COMPLEX*16 array, dimension (min(M,N))   
            The scalar factors of the elementary reflectors.   

    WORK    (workspace) COMPLEX*16 array, dimension (N)   

    RWORK   (workspace) DOUBLE PRECISION array, dimension (2*N)   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   

    Further Details   
    ===============   

    The matrix Q is represented as a product of elementary reflectors   

       Q = H(1) H(2) . . . H(n)   

    Each H(i) has the form   

       H = I - tau * v * v'   

    where tau is a complex scalar, and v is a complex vector with   
    v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i). 
  

    The matrix P is represented in jpvt as follows: If   
       jpvt(j) = i   
    then the jth column of P is the ith canonical unit vector.   

    ===================================================================== 
  


       Test the input arguments   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static integer c__1 = 1;
    
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublereal d__1;
    doublecomplex z__1;
    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);
    double z_abs(doublecomplex *), sqrt(doublereal);
    /* Local variables */
    static doublereal temp, temp2;
    static integer i, j, itemp;
    extern /* Subroutine */ int zlarf_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, doublecomplex *), zswap_(integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *), zgeqr2_(
	    integer *, integer *, doublecomplex *, integer *, doublecomplex *,
	     doublecomplex *, integer *);
    extern doublereal dznrm2_(integer *, doublecomplex *, integer *);
    static integer ma, mn;
    extern /* Subroutine */ int zunm2r_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    extern integer idamax_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int xerbla_(char *, integer *), zlarfg_(
	    integer *, doublecomplex *, doublecomplex *, integer *, 
	    doublecomplex *);
    static doublecomplex aii;
    static integer pvt;



#define JPVT(I) jpvt[(I)-1]
#define TAU(I) tau[(I)-1]
#define WORK(I) work[(I)-1]
#define RWORK(I) rwork[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    *info = 0;
    if (*m < 0) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*lda < max(1,*m)) {
	*info = -4;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("ZGEQPF", &i__1);
	return 0;
    }

    mn = min(*m,*n);

/*     Move initial columns up front */

    itemp = 1;
    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
	if (JPVT(i) != 0) {
	    if (i != itemp) {
		zswap_(m, &A(1,i), &c__1, &A(1,itemp), &
			c__1);
		JPVT(i) = JPVT(itemp);
		JPVT(itemp) = i;
	    } else {
		JPVT(i) = i;
	    }
	    ++itemp;
	} else {
	    JPVT(i) = i;
	}
/* L10: */
    }
    --itemp;

/*     Compute the QR factorization and update remaining columns */

    if (itemp > 0) {
	ma = min(itemp,*m);
	zgeqr2_(m, &ma, &A(1,1), lda, &TAU(1), &WORK(1), info);
	if (ma < *n) {
	    i__1 = *n - ma;
	    zunm2r_("Left", "Conjugate transpose", m, &i__1, &ma, &A(1,1)
		    , lda, &TAU(1), &A(1,ma+1), lda, &WORK(1), 
		    info);
	}
    }

    if (itemp < mn) {

/*        Initialize partial column norms. The first n elements of   
          work store the exact column norms. */

	i__1 = *n;
	for (i = itemp + 1; i <= *n; ++i) {
	    i__2 = *m - itemp;
	    RWORK(i) = dznrm2_(&i__2, &A(itemp+1,i), &c__1);
	    RWORK(*n + i) = RWORK(i);
/* L20: */
	}

/*        Compute factorization */

	i__1 = mn;
	for (i = itemp + 1; i <= mn; ++i) {

/*           Determine ith pivot column and swap if necessary */

	    i__2 = *n - i + 1;
	    pvt = i - 1 + idamax_(&i__2, &RWORK(i), &c__1);

	    if (pvt != i) {
		zswap_(m, &A(1,pvt), &c__1, &A(1,i), &
			c__1);
		itemp = JPVT(pvt);
		JPVT(pvt) = JPVT(i);
		JPVT(i) = itemp;
		RWORK(pvt) = RWORK(i);
		RWORK(*n + pvt) = RWORK(*n + i);
	    }

/*           Generate elementary reflector H(i) */

	    i__2 = i + i * a_dim1;
	    aii.r = A(i,i).r, aii.i = A(i,i).i;
	    i__2 = *m - i + 1;
/* Computing MIN */
	    i__3 = i + 1;
	    zlarfg_(&i__2, &aii, &A(min(i+1,*m),i), &c__1, &TAU(i)
		    );
	    i__2 = i + i * a_dim1;
	    A(i,i).r = aii.r, A(i,i).i = aii.i;

	    if (i < *n) {

/*              Apply H(i) to A(i:m,i+1:n) from the left */

		i__2 = i + i * a_dim1;
		aii.r = A(i,i).r, aii.i = A(i,i).i;
		i__2 = i + i * a_dim1;
		A(i,i).r = 1., A(i,i).i = 0.;
		i__2 = *m - i + 1;
		i__3 = *n - i;
		d_cnjg(&z__1, &TAU(i));
		zlarf_("Left", &i__2, &i__3, &A(i,i), &c__1, &z__1,
			 &A(i,i+1), lda, &WORK(1));
		i__2 = i + i * a_dim1;
		A(i,i).r = aii.r, A(i,i).i = aii.i;
	    }

/*           Update partial column norms */

	    i__2 = *n;
	    for (j = i + 1; j <= *n; ++j) {
		if (RWORK(j) != 0.) {
/* Computing 2nd power */
		    d__1 = z_abs(&A(i,j)) / RWORK(j);
		    temp = 1. - d__1 * d__1;
		    temp = max(temp,0.);
/* Computing 2nd power */
		    d__1 = RWORK(j) / RWORK(*n + j);
		    temp2 = temp * .05 * (d__1 * d__1) + 1.;
		    if (temp2 == 1.) {
			if (*m - i > 0) {
			    i__3 = *m - i;
			    RWORK(j) = dznrm2_(&i__3, &A(i+1,j), 
				    &c__1);
			    RWORK(*n + j) = RWORK(j);
			} else {
			    RWORK(j) = 0.;
			    RWORK(*n + j) = 0.;
			}
		    } else {
			RWORK(j) *= sqrt(temp);
		    }
		}
/* L30: */
	    }

/* L40: */
	}
    }
    return 0;

/*     End of ZGEQPF */

} /* zgeqpf_ */

