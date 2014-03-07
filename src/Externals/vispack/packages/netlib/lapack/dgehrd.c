#include "f2c.h"

/* Subroutine */ int dgehrd_(integer *n, integer *ilo, integer *ihi, 
	doublereal *a, integer *lda, doublereal *tau, doublereal *work, 
	integer *lwork, integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DGEHRD reduces a real general matrix A to upper Hessenberg form H by 
  
    an orthogonal similarity transformation:  Q' * A * Q = H .   

    Arguments   
    =========   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    ILO     (input) INTEGER   
    IHI     (input) INTEGER   
            It is assumed that A is already upper triangular in rows   
            and columns 1:ILO-1 and IHI+1:N. ILO and IHI are normally   
            set by a previous call to DGEBAL; otherwise they should be   
            set to 1 and N respectively. See Further Details.   
            1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.   

    A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the N-by-N general matrix to be reduced.   
            On exit, the upper triangle and the first subdiagonal of A   
            are overwritten with the upper Hessenberg matrix H, and the   
            elements below the first subdiagonal, with the array TAU,   
            represent the orthogonal matrix Q as a product of elementary 
  
            reflectors. See Further Details.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(1,N).   

    TAU     (output) DOUBLE PRECISION array, dimension (N-1)   
            The scalar factors of the elementary reflectors (see Further 
  
            Details). Elements 1:ILO-1 and IHI:N-1 of TAU are set to   
            zero.   

    WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK) 
  
            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INTEGER   
            The length of the array WORK.  LWORK >= max(1,N).   
            For optimum performance LWORK >= N*NB, where NB is the   
            optimal blocksize.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value.   

    Further Details   
    ===============   

    The matrix Q is represented as a product of (ihi-ilo) elementary   
    reflectors   

       Q = H(ilo) H(ilo+1) . . . H(ihi-1).   

    Each H(i) has the form   

       H(i) = I - tau * v * v'   

    where tau is a real scalar, and v is a real vector with   
    v(1:i) = 0, v(i+1) = 1 and v(ihi+1:n) = 0; v(i+2:ihi) is stored on   
    exit in A(i+2:ihi,i), and tau in TAU(i).   

    The contents of A are illustrated by the following example, with   
    n = 7, ilo = 2 and ihi = 6:   

    on entry,                        on exit,   

    ( a   a   a   a   a   a   a )    (  a   a   h   h   h   h   a )   
    (     a   a   a   a   a   a )    (      a   h   h   h   h   a )   
    (     a   a   a   a   a   a )    (      h   h   h   h   h   h )   
    (     a   a   a   a   a   a )    (      v2  h   h   h   h   h )   
    (     a   a   a   a   a   a )    (      v2  v3  h   h   h   h )   
    (     a   a   a   a   a   a )    (      v2  v3  v4  h   h   h )   
    (                         a )    (                          a )   

    where a denotes an element of the original matrix A, h denotes a   
    modified element of the upper Hessenberg matrix H, and vi denotes an 
  
    element of the vector defining H(i).   

    ===================================================================== 
  


       Test the input parameters   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static integer c__1 = 1;
    static integer c_n1 = -1;
    static integer c__3 = 3;
    static integer c__2 = 2;
    static integer c__65 = 65;
    static doublereal c_b21 = -1.;
    static doublereal c_b22 = 1.;
    
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;
    /* Local variables */
    static integer i;
    static doublereal t[4160]	/* was [65][64] */;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *);
    static integer nbmin, iinfo;
    extern /* Subroutine */ int dgehd2_(integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *);
    static integer ib;
    static doublereal ei;
    static integer nb, nh;
    extern /* Subroutine */ int dlarfb_(char *, char *, char *, char *, 
	    integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *), dlahrd_(integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *);
    static integer nx;
    extern /* Subroutine */ int xerbla_(char *, integer *);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static integer ldwork, iws;



#define T(I) t[(I)]
#define WAS(I) was[(I)]
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
    } else if (*lwork < max(1,*n)) {
	*info = -8;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DGEHRD", &i__1);
	return 0;
    }

/*     Set elements 1:ILO-1 and IHI:N-1 of TAU to zero */

    i__1 = *ilo - 1;
    for (i = 1; i <= *ilo-1; ++i) {
	TAU(i) = 0.;
/* L10: */
    }
    i__1 = *n - 1;
    for (i = max(1,*ihi); i <= *n-1; ++i) {
	TAU(i) = 0.;
/* L20: */
    }

/*     Quick return if possible */

    nh = *ihi - *ilo + 1;
    if (nh <= 1) {
	WORK(1) = 1.;
	return 0;
    }

/*     Determine the block size.   

   Computing MIN */
    i__1 = 64, i__2 = ilaenv_(&c__1, "DGEHRD", " ", n, ilo, ihi, &c_n1, 6L, 
	    1L);
    nb = min(i__1,i__2);
    nbmin = 2;
    iws = 1;
    if (nb > 1 && nb < nh) {

/*        Determine when to cross over from blocked to unblocked code 
  
          (last block is always handled by unblocked code).   

   Computing MAX */
	i__1 = nb, i__2 = ilaenv_(&c__3, "DGEHRD", " ", n, ilo, ihi, &c_n1, 
		6L, 1L);
	nx = max(i__1,i__2);
	if (nx < nh) {

/*           Determine if workspace is large enough for blocked co
de. */

	    iws = *n * nb;
	    if (*lwork < iws) {

/*              Not enough workspace to use optimal NB:  deter
mine the   
                minimum value of NB, and reduce NB or force us
e of   
                unblocked code.   

   Computing MAX */
		i__1 = 2, i__2 = ilaenv_(&c__2, "DGEHRD", " ", n, ilo, ihi, &
			c_n1, 6L, 1L);
		nbmin = max(i__1,i__2);
		if (*lwork >= *n * nbmin) {
		    nb = *lwork / *n;
		} else {
		    nb = 1;
		}
	    }
	}
    }
    ldwork = *n;

    if (nb < nbmin || nb >= nh) {

/*        Use unblocked code below */

	i = *ilo;

    } else {

/*        Use blocked code */

	i__1 = *ihi - 1 - nx;
	i__2 = nb;
	for (i = *ilo; nb < 0 ? i >= *ihi-1-nx : i <= *ihi-1-nx; i += nb) {
/* Computing MIN */
	    i__3 = nb, i__4 = *ihi - i;
	    ib = min(i__3,i__4);

/*           Reduce columns i:i+ib-1 to Hessenberg form, returning
 the   
             matrices V and T of the block reflector H = I - V*T*V
'   
             which performs the reduction, and also the matrix Y =
 A*V*T */

	    dlahrd_(ihi, &i, &ib, &A(1,i), lda, &TAU(i), t, &c__65,
		     &WORK(1), &ldwork);

/*           Apply the block reflector H to A(1:ihi,i+ib:ihi) from
 the   
             right, computing  A := A - Y * V'. V(i+ib,ib-1) must 
be set   
             to 1. */

	    ei = A(i+ib,i+ib-1);
	    A(i+ib,i+ib-1) = 1.;
	    i__3 = *ihi - i - ib + 1;
	    dgemm_("No transpose", "Transpose", ihi, &i__3, &ib, &c_b21, &
		    WORK(1), &ldwork, &A(i+ib,i), lda, &c_b22, &
		    A(1,i+ib), lda);
	    A(i+ib,i+ib-1) = ei;

/*           Apply the block reflector H to A(i+1:ihi,i+ib:n) from
 the   
             left */

	    i__3 = *ihi - i;
	    i__4 = *n - i - ib + 1;
	    dlarfb_("Left", "Transpose", "Forward", "Columnwise", &i__3, &
		    i__4, &ib, &A(i+1,i), lda, t, &c__65, &A(i+1,i+ib), lda, &WORK(1), &ldwork);
/* L30: */
	}
    }

/*     Use unblocked code to reduce the rest of the matrix */

    dgehd2_(n, &i, ihi, &A(1,1), lda, &TAU(1), &WORK(1), &iinfo);
    WORK(1) = (doublereal) iws;

    return 0;

/*     End of DGEHRD */

} /* dgehrd_ */

