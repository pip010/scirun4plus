#include "f2c.h"

/* Subroutine */ int cunglq_(integer *m, integer *n, integer *k, complex *a, 
	integer *lda, complex *tau, complex *work, integer *lwork, integer *
	info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    CUNGLQ generates an M-by-N complex matrix Q with orthonormal rows,   
    which is defined as the first M rows of a product of K elementary   
    reflectors of order N   

          Q  =  H(k)' . . . H(2)' H(1)'   

    as returned by CGELQF.   

    Arguments   
    =========   

    M       (input) INTEGER   
            The number of rows of the matrix Q. M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrix Q. N >= M.   

    K       (input) INTEGER   
            The number of elementary reflectors whose product defines the 
  
            matrix Q. M >= K >= 0.   

    A       (input/output) COMPLEX array, dimension (LDA,N)   
            On entry, the i-th row must contain the vector which defines 
  
            the elementary reflector H(i), for i = 1,2,...,k, as returned 
  
            by CGELQF in the first k rows of its array argument A.   
            On exit, the M-by-N matrix Q.   

    LDA     (input) INTEGER   
            The first dimension of the array A. LDA >= max(1,M).   

    TAU     (input) COMPLEX array, dimension (K)   
            TAU(i) must contain the scalar factor of the elementary   
            reflector H(i), as returned by CGELQF.   

    WORK    (workspace/output) COMPLEX array, dimension (LWORK)   
            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INTEGER   
            The dimension of the array WORK. LWORK >= max(1,M).   
            For optimum performance LWORK >= M*NB, where NB is   
            the optimal blocksize.   

    INFO    (output) INTEGER   
            = 0:  successful exit;   
            < 0:  if INFO = -i, the i-th argument has an illegal value   

    ===================================================================== 
  


       Test the input arguments   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static integer c__1 = 1;
    static integer c_n1 = -1;
    static integer c__3 = 3;
    static integer c__2 = 2;
    
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;
    /* Local variables */
    static integer i, j, l, nbmin, iinfo;
    extern /* Subroutine */ int cungl2_(integer *, integer *, integer *, 
	    complex *, integer *, complex *, complex *, integer *);
    static integer ib, nb, ki, kk;
    extern /* Subroutine */ int clarfb_(char *, char *, char *, char *, 
	    integer *, integer *, integer *, complex *, integer *, complex *, 
	    integer *, complex *, integer *, complex *, integer *);
    static integer nx;
    extern /* Subroutine */ int clarft_(char *, char *, integer *, integer *, 
	    complex *, integer *, complex *, complex *, integer *), xerbla_(char *, integer *);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static integer ldwork, iws;



#define TAU(I) tau[(I)-1]
#define WORK(I) work[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    *info = 0;
    if (*m < 0) {
	*info = -1;
    } else if (*n < *m) {
	*info = -2;
    } else if (*k < 0 || *k > *m) {
	*info = -3;
    } else if (*lda < max(1,*m)) {
	*info = -5;
    } else if (*lwork < max(1,*m)) {
	*info = -8;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("CUNGLQ", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*m <= 0) {
	WORK(1).r = 1.f, WORK(1).i = 0.f;
	return 0;
    }

/*     Determine the block size. */

    nb = ilaenv_(&c__1, "CUNGLQ", " ", m, n, k, &c_n1, 6L, 1L);
    nbmin = 2;
    nx = 0;
    iws = *m;
    if (nb > 1 && nb < *k) {

/*        Determine when to cross over from blocked to unblocked code.
   

   Computing MAX */
	i__1 = 0, i__2 = ilaenv_(&c__3, "CUNGLQ", " ", m, n, k, &c_n1, 6L, 1L)
		;
	nx = max(i__1,i__2);
	if (nx < *k) {

/*           Determine if workspace is large enough for blocked co
de. */

	    ldwork = *m;
	    iws = ldwork * nb;
	    if (*lwork < iws) {

/*              Not enough workspace to use optimal NB:  reduc
e NB and   
                determine the minimum value of NB. */

		nb = *lwork / ldwork;
/* Computing MAX */
		i__1 = 2, i__2 = ilaenv_(&c__2, "CUNGLQ", " ", m, n, k, &c_n1,
			 6L, 1L);
		nbmin = max(i__1,i__2);
	    }
	}
    }

    if (nb >= nbmin && nb < *k && nx < *k) {

/*        Use blocked code after the last block.   
          The first kk rows are handled by the block method. */

	ki = (*k - nx - 1) / nb * nb;
/* Computing MIN */
	i__1 = *k, i__2 = ki + nb;
	kk = min(i__1,i__2);

/*        Set A(kk+1:m,1:kk) to zero. */

	i__1 = kk;
	for (j = 1; j <= kk; ++j) {
	    i__2 = *m;
	    for (i = kk + 1; i <= *m; ++i) {
		i__3 = i + j * a_dim1;
		A(i,j).r = 0.f, A(i,j).i = 0.f;
/* L10: */
	    }
/* L20: */
	}
    } else {
	kk = 0;
    }

/*     Use unblocked code for the last or only block. */

    if (kk < *m) {
	i__1 = *m - kk;
	i__2 = *n - kk;
	i__3 = *k - kk;
	cungl2_(&i__1, &i__2, &i__3, &A(kk+1,kk+1), lda, &
		TAU(kk + 1), &WORK(1), &iinfo);
    }

    if (kk > 0) {

/*        Use blocked code */

	i__1 = -nb;
	for (i = ki + 1; -nb < 0 ? i >= 1 : i <= 1; i += -nb) {
/* Computing MIN */
	    i__2 = nb, i__3 = *k - i + 1;
	    ib = min(i__2,i__3);
	    if (i + ib <= *m) {

/*              Form the triangular factor of the block reflec
tor   
                H = H(i) H(i+1) . . . H(i+ib-1) */

		i__2 = *n - i + 1;
		clarft_("Forward", "Rowwise", &i__2, &ib, &A(i,i), 
			lda, &TAU(i), &WORK(1), &ldwork);

/*              Apply H' to A(i+ib:m,i:n) from the right */

		i__2 = *m - i - ib + 1;
		i__3 = *n - i + 1;
		clarfb_("Right", "Conjugate transpose", "Forward", "Rowwise", 
			&i__2, &i__3, &ib, &A(i,i), lda, &WORK(1), 
			&ldwork, &A(i+ib,i), lda, &WORK(ib + 1), 
			&ldwork);
	    }

/*           Apply H' to columns i:n of current block */

	    i__2 = *n - i + 1;
	    cungl2_(&ib, &i__2, &ib, &A(i,i), lda, &TAU(i), &WORK(
		    1), &iinfo);

/*           Set columns 1:i-1 of current block to zero */

	    i__2 = i - 1;
	    for (j = 1; j <= i-1; ++j) {
		i__3 = i + ib - 1;
		for (l = i; l <= i+ib-1; ++l) {
		    i__4 = l + j * a_dim1;
		    A(l,j).r = 0.f, A(l,j).i = 0.f;
/* L30: */
		}
/* L40: */
	    }
/* L50: */
	}
    }

    WORK(1).r = (real) iws, WORK(1).i = 0.f;
    return 0;

/*     End of CUNGLQ */

} /* cunglq_ */

