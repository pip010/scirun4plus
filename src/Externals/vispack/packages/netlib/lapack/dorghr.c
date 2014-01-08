#include "f2c.h"

/* Subroutine */ int dorghr_(integer *n, integer *ilo, integer *ihi, 
	doublereal *a, integer *lda, doublereal *tau, doublereal *work, 
	integer *lwork, integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DORGHR generates a real orthogonal matrix Q which is defined as the   
    product of IHI-ILO elementary reflectors of order N, as returned by   
    DGEHRD:   

    Q = H(ilo) H(ilo+1) . . . H(ihi-1).   

    Arguments   
    =========   

    N       (input) INTEGER   
            The order of the matrix Q. N >= 0.   

    ILO     (input) INTEGER   
    IHI     (input) INTEGER   
            ILO and IHI must have the same values as in the previous call 
  
            of DGEHRD. Q is equal to the unit matrix except in the   
            submatrix Q(ilo+1:ihi,ilo+1:ihi).   
            1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.   

    A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the vectors which define the elementary reflectors, 
  
            as returned by DGEHRD.   
            On exit, the N-by-N orthogonal matrix Q.   

    LDA     (input) INTEGER   
            The leading dimension of the array A. LDA >= max(1,N).   

    TAU     (input) DOUBLE PRECISION array, dimension (N-1)   
            TAU(i) must contain the scalar factor of the elementary   
            reflector H(i), as returned by DGEHRD.   

    WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK) 
  
            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INTEGER   
            The dimension of the array WORK. LWORK >= IHI-ILO.   
            For optimum performance LWORK >= (IHI-ILO)*NB, where NB is   
            the optimal blocksize.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   

    ===================================================================== 
  


       Test the input arguments   

    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    /* Local variables */
    static integer i, j, iinfo, nh;
    extern /* Subroutine */ int xerbla_(char *, integer *), dorgqr_(
	    integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, integer *);


#define TAU(I) tau[(I)-1]
#define WORK(I) work[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    *info = 0;
    if (*n < 0) {
	*info = -1;
    } else if (*ilo < 1 || *ilo > max(1,*n)) {
	*info = -2;
    } else if (*ihi < min(*ilo,*n) || *ihi > *n) {
	*info = -3;
    } else if (*lda < max(1,*n)) {
	*info = -5;
    } else /* if(complicated condition) */ {
/* Computing MAX */
	i__1 = 1, i__2 = *ihi - *ilo;
	if (*lwork < max(i__1,i__2)) {
	    *info = -8;
	}
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DORGHR", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0) {
	WORK(1) = 1.;
	return 0;
    }

/*     Shift the vectors which define the elementary reflectors one   
       column to the right, and set the first ilo and the last n-ihi   
       rows and columns to those of the unit matrix */

    i__1 = *ilo + 1;
    for (j = *ihi; j >= *ilo+1; --j) {
	i__2 = j - 1;
	for (i = 1; i <= j-1; ++i) {
	    A(i,j) = 0.;
/* L10: */
	}
	i__2 = *ihi;
	for (i = j + 1; i <= *ihi; ++i) {
	    A(i,j) = A(i,j-1);
/* L20: */
	}
	i__2 = *n;
	for (i = *ihi + 1; i <= *n; ++i) {
	    A(i,j) = 0.;
/* L30: */
	}
/* L40: */
    }
    i__1 = *ilo;
    for (j = 1; j <= *ilo; ++j) {
	i__2 = *n;
	for (i = 1; i <= *n; ++i) {
	    A(i,j) = 0.;
/* L50: */
	}
	A(j,j) = 1.;
/* L60: */
    }
    i__1 = *n;
    for (j = *ihi + 1; j <= *n; ++j) {
	i__2 = *n;
	for (i = 1; i <= *n; ++i) {
	    A(i,j) = 0.;
/* L70: */
	}
	A(j,j) = 1.;
/* L80: */
    }

    nh = *ihi - *ilo;
    if (nh > 0) {

/*        Generate Q(ilo+1:ihi,ilo+1:ihi) */

	dorgqr_(&nh, &nh, &nh, &A(*ilo+1,*ilo+1), lda, &TAU(*
		ilo), &WORK(1), lwork, &iinfo);
    }
    return 0;

/*     End of DORGHR */

} /* dorghr_ */

