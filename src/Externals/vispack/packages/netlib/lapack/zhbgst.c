#include "f2c.h"

/* Subroutine */ int zhbgst_(char *vect, char *uplo, integer *n, integer *ka, 
	integer *kb, doublecomplex *ab, integer *ldab, doublecomplex *bb, 
	integer *ldbb, doublecomplex *x, integer *ldx, doublecomplex *work, 
	doublereal *rwork, integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    ZHBGST reduces a complex Hermitian-definite banded generalized   
    eigenproblem  A*x = lambda*B*x  to standard form  C*y = lambda*y,   
    such that C has the same bandwidth as A.   

    B must have been previously factorized as S**H*S by ZPBSTF, using a   
    split Cholesky factorization. A is overwritten by C = X**H*A*X, where 
  
    X = S**(-1)*Q and Q is a unitary matrix chosen to preserve the   
    bandwidth of A.   

    Arguments   
    =========   

    VECT    (input) CHARACTER*1   
            = 'N':  do not form the transformation matrix X;   
            = 'V':  form X.   

    UPLO    (input) CHARACTER*1   
            = 'U':  Upper triangle of A is stored;   
            = 'L':  Lower triangle of A is stored.   

    N       (input) INTEGER   
            The order of the matrices A and B.  N >= 0.   

    KA      (input) INTEGER   
            The number of superdiagonals of the matrix A if UPLO = 'U',   
            or the number of subdiagonals if UPLO = 'L'.  KA >= 0.   

    KB      (input) INTEGER   
            The number of superdiagonals of the matrix B if UPLO = 'U',   
            or the number of subdiagonals if UPLO = 'L'.  KA >= KB >= 0. 
  

    AB      (input/output) COMPLEX*16 array, dimension (LDAB,N)   
            On entry, the upper or lower triangle of the Hermitian band   
            matrix A, stored in the first ka+1 rows of the array.  The   
            j-th column of A is stored in the j-th column of the array AB 
  
            as follows:   
            if UPLO = 'U', AB(ka+1+i-j,j) = A(i,j) for max(1,j-ka)<=i<=j; 
  
            if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+ka). 
  

            On exit, the transformed matrix X**H*A*X, stored in the same 
  
            format as A.   

    LDAB    (input) INTEGER   
            The leading dimension of the array AB.  LDAB >= KA+1.   

    BB      (input) COMPLEX*16 array, dimension (LDBB,N)   
            The banded factor S from the split Cholesky factorization of 
  
            B, as returned by ZPBSTF, stored in the first kb+1 rows of   
            the array.   

    LDBB    (input) INTEGER   
            The leading dimension of the array BB.  LDBB >= KB+1.   

    X       (output) COMPLEX*16 array, dimension (LDX,N)   
            If VECT = 'V', the n-by-n matrix X.   
            If VECT = 'N', the array X is not referenced.   

    LDX     (input) INTEGER   
            The leading dimension of the array X.   
            LDX >= max(1,N) if VECT = 'V'; LDX >= 1 otherwise.   

    WORK    (workspace) COMPLEX*16 array, dimension (N)   

    RWORK   (workspace) DOUBLE PRECISION array, dimension (N)   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value.   

    ===================================================================== 
  


       Test the input parameters   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static doublecomplex c_b1 = {0.,0.};
    static doublecomplex c_b2 = {1.,0.};
    static integer c__1 = 1;
    
    /* System generated locals */
    integer ab_dim1, ab_offset, bb_dim1, bb_offset, x_dim1, x_offset, i__1, 
	    i__2, i__3, i__4, i__5, i__6, i__7, i__8;
    doublereal d__1;
    doublecomplex z__1, z__2, z__3, z__4, z__5, z__6, z__7, z__8, z__9, z__10;
    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);
    /* Local variables */
    static integer inca;
    extern /* Subroutine */ int zrot_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublereal *, doublecomplex *);
    static integer i, j, k, l, m;
    static doublecomplex t;
    extern logical lsame_(char *, char *);
    extern /* Subroutine */ int zgerc_(integer *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    static integer i0, i1;
    static logical upper;
    static integer i2, j1, j2;
    extern /* Subroutine */ int zgeru_(integer *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    static logical wantx;
    extern /* Subroutine */ int zlar2v_(integer *, doublecomplex *, 
	    doublecomplex *, doublecomplex *, integer *, doublereal *, 
	    doublecomplex *, integer *);
    static doublecomplex ra;
    static integer nr, nx;
    extern /* Subroutine */ int xerbla_(char *, integer *), zdscal_(
	    integer *, doublereal *, doublecomplex *, integer *);
    static logical update;
    extern /* Subroutine */ int zlacgv_(integer *, doublecomplex *, integer *)
	    ;
    static integer ka1, kb1;
    extern /* Subroutine */ int zlaset_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, doublecomplex *, integer *), zlartg_(doublecomplex *, doublecomplex *, doublereal *, 
	    doublecomplex *, doublecomplex *);
    static doublecomplex ra1;
    extern /* Subroutine */ int zlargv_(integer *, doublecomplex *, integer *,
	     doublecomplex *, integer *, doublereal *, integer *);
    static integer j1t, j2t;
    extern /* Subroutine */ int zlartv_(integer *, doublecomplex *, integer *,
	     doublecomplex *, integer *, doublereal *, doublecomplex *, 
	    integer *);
    static doublereal bii;
    static integer kbt, nrt;



#define WORK(I) work[(I)-1]
#define RWORK(I) rwork[(I)-1]

#define AB(I,J) ab[(I)-1 + ((J)-1)* ( *ldab)]
#define BB(I,J) bb[(I)-1 + ((J)-1)* ( *ldbb)]
#define X(I,J) x[(I)-1 + ((J)-1)* ( *ldx)]

    wantx = lsame_(vect, "V");
    upper = lsame_(uplo, "U");
    ka1 = *ka + 1;
    kb1 = *kb + 1;
    *info = 0;
    if (! wantx && ! lsame_(vect, "N")) {
	*info = -1;
    } else if (! upper && ! lsame_(uplo, "L")) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*ka < 0) {
	*info = -4;
    } else if (*kb < 0) {
	*info = -5;
    } else if (*ldab < *ka + 1) {
	*info = -7;
    } else if (*ldbb < *kb + 1) {
	*info = -9;
    } else if (*ldx < 1 || wantx && *ldx < max(1,*n)) {
	*info = -11;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("ZHBGST", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return 0;
    }

    inca = *ldab * ka1;

/*     Initialize X to the unit matrix, if needed */

    if (wantx) {
	zlaset_("Full", n, n, &c_b1, &c_b2, &X(1,1), ldx);
    }

/*     Set M to the splitting point m. It must be the same value as is   
       used in ZPBSTF. The chosen value allows the arrays WORK and RWORK 
  
       to be of dimension (N). */

    m = (*n + *kb) / 2;

/*     The routine works in two phases, corresponding to the two halves   
       of the split Cholesky factorization of B as S**H*S where   

       S = ( U    )   
           ( M  L )   

       with U upper triangular of order m, and L lower triangular of   
       order n-m. S has the same bandwidth as B.   

       S is treated as a product of elementary matrices:   

       S = S(m)*S(m-1)*...*S(2)*S(1)*S(m+1)*S(m+2)*...*S(n-1)*S(n)   

       where S(i) is determined by the i-th row of S.   

       In phase 1, the index i takes the values n, n-1, ... , m+1;   
       in phase 2, it takes the values 1, 2, ... , m.   

       For each value of i, the current matrix A is updated by forming   
       inv(S(i))**H*A*inv(S(i)). This creates a triangular bulge outside 
  
       the band of A. The bulge is then pushed down toward the bottom of 
  
       A in phase 1, and up toward the top of A in phase 2, by applying   
       plane rotations.   

       There are kb*(kb+1)/2 elements in the bulge, but at most 2*kb-1   
       of them are linearly independent, so annihilating a bulge requires 
  
       only 2*kb-1 plane rotations. The rotations are divided into a 1st 
  
       set of kb-1 rotations, and a 2nd set of kb rotations.   

       Wherever possible, rotations are generated and applied in vector   
       operations of length NR between the indices J1 and J2 (sometimes   
       replaced by modified values NRT, J1T or J2T).   

       The real cosines and complex sines of the rotations are stored in 
  
       the arrays RWORK and WORK, those of the 1st set in elements   
       2:m-kb-1, and those of the 2nd set in elements m-kb+1:n.   

       The bulges are not formed explicitly; nonzero elements outside the 
  
       band are created only when they are required for generating new   
       rotations; they are stored in the array WORK, in positions where   
       they are later overwritten by the sines of the rotations which   
       annihilate them.   

       **************************** Phase 1 ***************************** 
  

       The logical structure of this phase is:   

       UPDATE = .TRUE.   
       DO I = N, M + 1, -1   
          use S(i) to update A and create a new bulge   
          apply rotations to push all bulges KA positions downward   
       END DO   
       UPDATE = .FALSE.   
       DO I = M + KA + 1, N - 1   
          apply rotations to push all bulges KA positions downward   
       END DO   

       To avoid duplicating code, the two loops are merged. */

    update = TRUE_;
    i = *n + 1;
L10:
    if (update) {
	--i;
/* Computing MIN */
	i__1 = *kb, i__2 = i - 1;
	kbt = min(i__1,i__2);
	i0 = i - 1;
/* Computing MIN */
	i__1 = *n, i__2 = i + *ka;
	i1 = min(i__1,i__2);
	i2 = i - kbt + ka1;
	if (i < m + 1) {
	    update = FALSE_;
	    ++i;
	    i0 = m;
	    if (*ka == 0) {
		goto L480;
	    }
	    goto L10;
	}
    } else {
	i += *ka;
	if (i > *n - 1) {
	    goto L480;
	}
    }

    if (upper) {

/*        Transform A, working with the upper triangle */

	if (update) {

/*           Form  inv(S(i))**H * A * inv(S(i)) */

	    i__1 = kb1 + i * bb_dim1;
	    bii = BB(kb1,i).r;
	    i__1 = ka1 + i * ab_dim1;
	    i__2 = ka1 + i * ab_dim1;
	    d__1 = AB(ka1,i).r / bii / bii;
	    AB(ka1,i).r = d__1, AB(ka1,i).i = 0.;
	    i__1 = i1;
	    for (j = i + 1; j <= i1; ++j) {
		i__2 = i - j + ka1 + j * ab_dim1;
		i__3 = i - j + ka1 + j * ab_dim1;
		z__1.r = AB(i-j+ka1,j).r / bii, z__1.i = AB(i-j+ka1,j).i / bii;
		AB(i-j+ka1,j).r = z__1.r, AB(i-j+ka1,j).i = z__1.i;
/* L20: */
	    }
/* Computing MAX */
	    i__1 = 1, i__2 = i - *ka;
	    i__3 = i - 1;
	    for (j = max(1,i-*ka); j <= i-1; ++j) {
		i__1 = j - i + ka1 + i * ab_dim1;
		i__2 = j - i + ka1 + i * ab_dim1;
		z__1.r = AB(j-i+ka1,i).r / bii, z__1.i = AB(j-i+ka1,i).i / bii;
		AB(j-i+ka1,i).r = z__1.r, AB(j-i+ka1,i).i = z__1.i;
/* L30: */
	    }
	    i__3 = i - 1;
	    for (k = i - kbt; k <= i-1; ++k) {
		i__1 = k;
		for (j = i - kbt; j <= k; ++j) {
		    i__2 = j - k + ka1 + k * ab_dim1;
		    i__4 = j - k + ka1 + k * ab_dim1;
		    i__5 = j - i + kb1 + i * bb_dim1;
		    d_cnjg(&z__5, &AB(k-i+ka1,i));
		    z__4.r = BB(j-i+kb1,i).r * z__5.r - BB(j-i+kb1,i).i * z__5.i, 
			    z__4.i = BB(j-i+kb1,i).r * z__5.i + BB(j-i+kb1,i).i * 
			    z__5.r;
		    z__3.r = AB(j-k+ka1,k).r - z__4.r, z__3.i = AB(j-k+ka1,k).i - 
			    z__4.i;
		    d_cnjg(&z__7, &BB(k-i+kb1,i));
		    i__6 = j - i + ka1 + i * ab_dim1;
		    z__6.r = z__7.r * AB(j-i+ka1,i).r - z__7.i * AB(j-i+ka1,i).i, 
			    z__6.i = z__7.r * AB(j-i+ka1,i).i + z__7.i * AB(j-i+ka1,i)
			    .r;
		    z__2.r = z__3.r - z__6.r, z__2.i = z__3.i - z__6.i;
		    i__7 = ka1 + i * ab_dim1;
		    d__1 = AB(ka1,i).r;
		    i__8 = j - i + kb1 + i * bb_dim1;
		    z__9.r = d__1 * BB(j-i+kb1,i).r, z__9.i = d__1 * BB(j-i+kb1,i).i;
		    d_cnjg(&z__10, &BB(k-i+kb1,i));
		    z__8.r = z__9.r * z__10.r - z__9.i * z__10.i, z__8.i = 
			    z__9.r * z__10.i + z__9.i * z__10.r;
		    z__1.r = z__2.r + z__8.r, z__1.i = z__2.i + z__8.i;
		    AB(j-k+ka1,k).r = z__1.r, AB(j-k+ka1,k).i = z__1.i;
/* L40: */
		}
/* Computing MAX */
		i__1 = 1, i__2 = i - *ka;
		i__4 = i - kbt - 1;
		for (j = max(1,i-*ka); j <= i-kbt-1; ++j) {
		    i__1 = j - k + ka1 + k * ab_dim1;
		    i__2 = j - k + ka1 + k * ab_dim1;
		    d_cnjg(&z__3, &BB(k-i+kb1,i));
		    i__5 = j - i + ka1 + i * ab_dim1;
		    z__2.r = z__3.r * AB(j-i+ka1,i).r - z__3.i * AB(j-i+ka1,i).i, 
			    z__2.i = z__3.r * AB(j-i+ka1,i).i + z__3.i * AB(j-i+ka1,i)
			    .r;
		    z__1.r = AB(j-k+ka1,k).r - z__2.r, z__1.i = AB(j-k+ka1,k).i - 
			    z__2.i;
		    AB(j-k+ka1,k).r = z__1.r, AB(j-k+ka1,k).i = z__1.i;
/* L50: */
		}
/* L60: */
	    }
	    i__3 = i1;
	    for (j = i; j <= i1; ++j) {
/* Computing MAX */
		i__4 = j - *ka, i__1 = i - kbt;
		i__2 = i - 1;
		for (k = max(j-*ka,i-kbt); k <= i-1; ++k) {
		    i__4 = k - j + ka1 + j * ab_dim1;
		    i__1 = k - j + ka1 + j * ab_dim1;
		    i__5 = k - i + kb1 + i * bb_dim1;
		    i__6 = i - j + ka1 + j * ab_dim1;
		    z__2.r = BB(k-i+kb1,i).r * AB(i-j+ka1,j).r - BB(k-i+kb1,i).i * AB(i-j+ka1,j)
			    .i, z__2.i = BB(k-i+kb1,i).r * AB(i-j+ka1,j).i + BB(k-i+kb1,i).i 
			    * AB(i-j+ka1,j).r;
		    z__1.r = AB(k-j+ka1,j).r - z__2.r, z__1.i = AB(k-j+ka1,j).i - 
			    z__2.i;
		    AB(k-j+ka1,j).r = z__1.r, AB(k-j+ka1,j).i = z__1.i;
/* L70: */
		}
/* L80: */
	    }

	    if (wantx) {

/*              post-multiply X by inv(S(i)) */

		i__3 = *n - m;
		d__1 = 1. / bii;
		zdscal_(&i__3, &d__1, &X(m+1,i), &c__1);
		if (kbt > 0) {
		    i__3 = *n - m;
		    z__1.r = -1., z__1.i = 0.;
		    zgerc_(&i__3, &kbt, &z__1, &X(m+1,i), &c__1, 
			    &BB(kb1-kbt,i), &c__1, &X(m+1,i-kbt), ldx);
		}
	    }

/*           store a(i,i1) in RA1 for use in next loop over K */

	    i__3 = i - i1 + ka1 + i1 * ab_dim1;
	    ra1.r = AB(i-i1+ka1,i1).r, ra1.i = AB(i-i1+ka1,i1).i;
	}

/*        Generate and apply vectors of rotations to chase all the   
          existing bulges KA positions down toward the bottom of the 
  
          band */

	i__3 = *kb - 1;
	for (k = 1; k <= *kb-1; ++k) {
	    if (update) {

/*              Determine the rotations which would annihilate
 the bulge   
                which has in theory just been created */

		if (i - k + *ka < *n && i - k > 1) {

/*                 generate rotation to annihilate a(i,i-k
+ka+1) */

		    zlartg_(&AB(k+1,i-k+*ka), &ra1, &
			    RWORK(i - k + *ka - m), &WORK(i - k + *ka - m), &
			    ra);

/*                 create nonzero element a(i-k,i-k+ka+1) 
outside the   
                   band and store it in WORK(i-k) */

		    i__2 = kb1 - k + i * bb_dim1;
		    z__2.r = -BB(kb1-k,i).r, z__2.i = -BB(kb1-k,i).i;
		    z__1.r = z__2.r * ra1.r - z__2.i * ra1.i, z__1.i = z__2.r 
			    * ra1.i + z__2.i * ra1.r;
		    t.r = z__1.r, t.i = z__1.i;
		    i__2 = i - k;
		    i__4 = i - k + *ka - m;
		    z__2.r = RWORK(i-k+*ka-m) * t.r, z__2.i = RWORK(i-k+*ka-m) * t.i;
		    d_cnjg(&z__4, &WORK(i - k + *ka - m));
		    i__1 = (i - k + *ka) * ab_dim1 + 1;
		    z__3.r = z__4.r * AB(1,i-k+*ka).r - z__4.i * AB(1,i-k+*ka).i, 
			    z__3.i = z__4.r * AB(1,i-k+*ka).i + z__4.i * AB(1,i-k+*ka)
			    .r;
		    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
		    WORK(i-k).r = z__1.r, WORK(i-k).i = z__1.i;
		    i__2 = (i - k + *ka) * ab_dim1 + 1;
		    i__4 = i - k + *ka - m;
		    z__2.r = WORK(i-k+*ka-m).r * t.r - WORK(i-k+*ka-m).i * t.i, z__2.i =
			     WORK(i-k+*ka-m).r * t.i + WORK(i-k+*ka-m).i * t.r;
		    i__1 = i - k + *ka - m;
		    i__5 = (i - k + *ka) * ab_dim1 + 1;
		    z__3.r = RWORK(i-k+*ka-m) * AB(1,i-k+*ka).r, z__3.i = RWORK(i-k+*ka-m) * 
			    AB(1,i-k+*ka).i;
		    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
		    AB(1,i-k+*ka).r = z__1.r, AB(1,i-k+*ka).i = z__1.i;
		    ra1.r = ra.r, ra1.i = ra.i;
		}
	    }
/* Computing MAX */
	    i__2 = 1, i__4 = k - i0 + 2;
	    j2 = i - k - 1 + max(i__2,i__4) * ka1;
	    nr = (*n - j2 + *ka) / ka1;
	    j1 = j2 + (nr - 1) * ka1;
	    if (update) {
/* Computing MAX */
		i__2 = j2, i__4 = i + (*ka << 1) - k + 1;
		j2t = max(i__2,i__4);
	    } else {
		j2t = j2;
	    }
	    nrt = (*n - j2t + *ka) / ka1;
	    i__2 = j1;
	    i__4 = ka1;
	    for (j = j2t; ka1 < 0 ? j >= j1 : j <= j1; j += ka1) {

/*              create nonzero element a(j-ka,j+1) outside the
 band   
                and store it in WORK(j-m) */

		i__1 = j - m;
		i__5 = j - m;
		i__6 = (j + 1) * ab_dim1 + 1;
		z__1.r = WORK(j-m).r * AB(1,j+1).r - WORK(j-m).i * AB(1,j+1)
			.i, z__1.i = WORK(j-m).r * AB(1,j+1).i + WORK(j-m).i 
			* AB(1,j+1).r;
		WORK(j-m).r = z__1.r, WORK(j-m).i = z__1.i;
		i__1 = (j + 1) * ab_dim1 + 1;
		i__5 = j - m;
		i__6 = (j + 1) * ab_dim1 + 1;
		z__1.r = RWORK(j-m) * AB(1,j+1).r, z__1.i = RWORK(j-m) * AB(1,j+1).i;
		AB(1,j+1).r = z__1.r, AB(1,j+1).i = z__1.i;
/* L90: */
	    }

/*           generate rotations in 1st set to annihilate elements 
which   
             have been created outside the band */

	    if (nrt > 0) {
		zlargv_(&nrt, &AB(1,j2t), &inca, &WORK(j2t - m), &
			ka1, &RWORK(j2t - m), &ka1);
	    }
	    if (nr > 0) {

/*              apply rotations in 1st set from the right */

		i__4 = *ka - 1;
		for (l = 1; l <= *ka-1; ++l) {
		    zlartv_(&nr, &AB(ka1-l,j2), &inca, &AB(*ka-l,j2+1), &inca, &RWORK(j2 - m), 
			    &WORK(j2 - m), &ka1);
/* L100: */
		}

/*              apply rotations in 1st set from both sides to 
diagonal   
                blocks */

		zlar2v_(&nr, &AB(ka1,j2), &AB(ka1,j2+1), &AB(*ka,j2+1), &inca, &
			RWORK(j2 - m), &WORK(j2 - m), &ka1);

		zlacgv_(&nr, &WORK(j2 - m), &ka1);
	    }

/*           start applying rotations in 1st set from the left */

	    i__4 = *kb - k + 1;
	    for (l = *ka - 1; l >= *kb-k+1; --l) {
		nrt = (*n - j2 + l) / ka1;
		if (nrt > 0) {
		    zlartv_(&nrt, &AB(l,j2+ka1-l), &inca, &
			    AB(l+1,j2+ka1-l), &inca, &
			    RWORK(j2 - m), &WORK(j2 - m), &ka1);
		}
/* L110: */
	    }

	    if (wantx) {

/*              post-multiply X by product of rotations in 1st
 set */

		i__4 = j1;
		i__2 = ka1;
		for (j = j2; ka1 < 0 ? j >= j1 : j <= j1; j += ka1) {
		    i__1 = *n - m;
		    d_cnjg(&z__1, &WORK(j - m));
		    zrot_(&i__1, &X(m+1,j), &c__1, &X(m+1,j+1), &c__1, &RWORK(j - m), &z__1);
/* L120: */
		}
	    }
/* L130: */
	}

	if (update) {
	    if (i2 <= *n && kbt > 0) {

/*              create nonzero element a(i-kbt,i-kbt+ka+1) out
side the   
                band and store it in WORK(i-kbt) */

		i__3 = i - kbt;
		i__2 = kb1 - kbt + i * bb_dim1;
		z__2.r = -BB(kb1-kbt,i).r, z__2.i = -BB(kb1-kbt,i).i;
		z__1.r = z__2.r * ra1.r - z__2.i * ra1.i, z__1.i = z__2.r * 
			ra1.i + z__2.i * ra1.r;
		WORK(i-kbt).r = z__1.r, WORK(i-kbt).i = z__1.i;
	    }
	}

	for (k = *kb; k >= 1; --k) {
	    if (update) {
/* Computing MAX */
		i__3 = 2, i__2 = k - i0 + 1;
		j2 = i - k - 1 + max(i__3,i__2) * ka1;
	    } else {
/* Computing MAX */
		i__3 = 1, i__2 = k - i0 + 1;
		j2 = i - k - 1 + max(i__3,i__2) * ka1;
	    }

/*           finish applying rotations in 2nd set from the left */

	    for (l = *kb - k; l >= 1; --l) {
		nrt = (*n - j2 + *ka + l) / ka1;
		if (nrt > 0) {
		    zlartv_(&nrt, &AB(l,j2-l+1), &inca, &AB(l+1,j2-l+1), &inca, &RWORK(j2 
			    - *ka), &WORK(j2 - *ka), &ka1);
		}
/* L140: */
	    }
	    nr = (*n - j2 + *ka) / ka1;
	    j1 = j2 + (nr - 1) * ka1;
	    i__3 = j2;
	    i__2 = -ka1;
	    for (j = j1; -ka1 < 0 ? j >= j2 : j <= j2; j += -ka1) {
		i__4 = j;
		i__1 = j - *ka;
		WORK(j).r = WORK(j-*ka).r, WORK(j).i = WORK(j-*ka).i;
		RWORK(j) = RWORK(j - *ka);
/* L150: */
	    }
	    i__2 = j1;
	    i__3 = ka1;
	    for (j = j2; ka1 < 0 ? j >= j1 : j <= j1; j += ka1) {

/*              create nonzero element a(j-ka,j+1) outside the
 band   
                and store it in WORK(j) */

		i__4 = j;
		i__1 = j;
		i__5 = (j + 1) * ab_dim1 + 1;
		z__1.r = WORK(j).r * AB(1,j+1).r - WORK(j).i * AB(1,j+1)
			.i, z__1.i = WORK(j).r * AB(1,j+1).i + WORK(j).i 
			* AB(1,j+1).r;
		WORK(j).r = z__1.r, WORK(j).i = z__1.i;
		i__4 = (j + 1) * ab_dim1 + 1;
		i__1 = j;
		i__5 = (j + 1) * ab_dim1 + 1;
		z__1.r = RWORK(j) * AB(1,j+1).r, z__1.i = RWORK(j) * AB(1,j+1).i;
		AB(1,j+1).r = z__1.r, AB(1,j+1).i = z__1.i;
/* L160: */
	    }
	    if (update) {
		if (i - k < *n - *ka && k <= kbt) {
		    i__3 = i - k + *ka;
		    i__2 = i - k;
		    WORK(i-k+*ka).r = WORK(i-k).r, WORK(i-k+*ka).i = WORK(i-k).i;
		}
	    }
/* L170: */
	}

	for (k = *kb; k >= 1; --k) {
/* Computing MAX */
	    i__3 = 1, i__2 = k - i0 + 1;
	    j2 = i - k - 1 + max(i__3,i__2) * ka1;
	    nr = (*n - j2 + *ka) / ka1;
	    j1 = j2 + (nr - 1) * ka1;
	    if (nr > 0) {

/*              generate rotations in 2nd set to annihilate el
ements   
                which have been created outside the band */

		zlargv_(&nr, &AB(1,j2), &inca, &WORK(j2), &ka1, &
			RWORK(j2), &ka1);

/*              apply rotations in 2nd set from the right */

		i__3 = *ka - 1;
		for (l = 1; l <= *ka-1; ++l) {
		    zlartv_(&nr, &AB(ka1-l,j2), &inca, &AB(*ka-l,j2+1), &inca, &RWORK(j2), &
			    WORK(j2), &ka1);
/* L180: */
		}

/*              apply rotations in 2nd set from both sides to 
diagonal   
                blocks */

		zlar2v_(&nr, &AB(ka1,j2), &AB(ka1,j2+1), &AB(*ka,j2+1), &inca, &
			RWORK(j2), &WORK(j2), &ka1);

		zlacgv_(&nr, &WORK(j2), &ka1);
	    }

/*           start applying rotations in 2nd set from the left */

	    i__3 = *kb - k + 1;
	    for (l = *ka - 1; l >= *kb-k+1; --l) {
		nrt = (*n - j2 + l) / ka1;
		if (nrt > 0) {
		    zlartv_(&nrt, &AB(l,j2+ka1-l), &inca, &
			    AB(l+1,j2+ka1-l), &inca, &
			    RWORK(j2), &WORK(j2), &ka1);
		}
/* L190: */
	    }

	    if (wantx) {

/*              post-multiply X by product of rotations in 2nd
 set */

		i__3 = j1;
		i__2 = ka1;
		for (j = j2; ka1 < 0 ? j >= j1 : j <= j1; j += ka1) {
		    i__4 = *n - m;
		    d_cnjg(&z__1, &WORK(j));
		    zrot_(&i__4, &X(m+1,j), &c__1, &X(m+1,j+1), &c__1, &RWORK(j), &z__1);
/* L200: */
		}
	    }
/* L210: */
	}

	i__2 = *kb - 1;
	for (k = 1; k <= *kb-1; ++k) {
/* Computing MAX */
	    i__3 = 1, i__4 = k - i0 + 2;
	    j2 = i - k - 1 + max(i__3,i__4) * ka1;

/*           finish applying rotations in 1st set from the left */

	    for (l = *kb - k; l >= 1; --l) {
		nrt = (*n - j2 + l) / ka1;
		if (nrt > 0) {
		    zlartv_(&nrt, &AB(l,j2+ka1-l), &inca, &
			    AB(l+1,j2+ka1-l), &inca, &
			    RWORK(j2 - m), &WORK(j2 - m), &ka1);
		}
/* L220: */
	    }
/* L230: */
	}

	i__2 = i2 + *ka;
	for (j = *n - 1; j >= i2+*ka; --j) {
	    RWORK(j - m) = RWORK(j - *ka - m);
	    i__3 = j - m;
	    i__4 = j - *ka - m;
	    WORK(j-m).r = WORK(j-*ka-m).r, WORK(j-m).i = WORK(j-*ka-m).i;
/* L240: */
	}

    } else {

/*        Transform A, working with the lower triangle */

	if (update) {

/*           Form  inv(S(i))**H * A * inv(S(i)) */

	    i__2 = i * bb_dim1 + 1;
	    bii = BB(1,i).r;
	    i__2 = i * ab_dim1 + 1;
	    i__3 = i * ab_dim1 + 1;
	    d__1 = AB(1,i).r / bii / bii;
	    AB(1,i).r = d__1, AB(1,i).i = 0.;
	    i__2 = i1;
	    for (j = i + 1; j <= i1; ++j) {
		i__3 = j - i + 1 + i * ab_dim1;
		i__4 = j - i + 1 + i * ab_dim1;
		z__1.r = AB(j-i+1,i).r / bii, z__1.i = AB(j-i+1,i).i / bii;
		AB(j-i+1,i).r = z__1.r, AB(j-i+1,i).i = z__1.i;
/* L250: */
	    }
/* Computing MAX */
	    i__2 = 1, i__3 = i - *ka;
	    i__4 = i - 1;
	    for (j = max(1,i-*ka); j <= i-1; ++j) {
		i__2 = i - j + 1 + j * ab_dim1;
		i__3 = i - j + 1 + j * ab_dim1;
		z__1.r = AB(i-j+1,j).r / bii, z__1.i = AB(i-j+1,j).i / bii;
		AB(i-j+1,j).r = z__1.r, AB(i-j+1,j).i = z__1.i;
/* L260: */
	    }
	    i__4 = i - 1;
	    for (k = i - kbt; k <= i-1; ++k) {
		i__2 = k;
		for (j = i - kbt; j <= k; ++j) {
		    i__3 = k - j + 1 + j * ab_dim1;
		    i__1 = k - j + 1 + j * ab_dim1;
		    i__5 = i - j + 1 + j * bb_dim1;
		    d_cnjg(&z__5, &AB(i-k+1,k));
		    z__4.r = BB(i-j+1,j).r * z__5.r - BB(i-j+1,j).i * z__5.i, 
			    z__4.i = BB(i-j+1,j).r * z__5.i + BB(i-j+1,j).i * 
			    z__5.r;
		    z__3.r = AB(k-j+1,j).r - z__4.r, z__3.i = AB(k-j+1,j).i - 
			    z__4.i;
		    d_cnjg(&z__7, &BB(i-k+1,k));
		    i__6 = i - j + 1 + j * ab_dim1;
		    z__6.r = z__7.r * AB(i-j+1,j).r - z__7.i * AB(i-j+1,j).i, 
			    z__6.i = z__7.r * AB(i-j+1,j).i + z__7.i * AB(i-j+1,j)
			    .r;
		    z__2.r = z__3.r - z__6.r, z__2.i = z__3.i - z__6.i;
		    i__7 = i * ab_dim1 + 1;
		    d__1 = AB(1,i).r;
		    i__8 = i - j + 1 + j * bb_dim1;
		    z__9.r = d__1 * BB(i-j+1,j).r, z__9.i = d__1 * BB(i-j+1,j).i;
		    d_cnjg(&z__10, &BB(i-k+1,k));
		    z__8.r = z__9.r * z__10.r - z__9.i * z__10.i, z__8.i = 
			    z__9.r * z__10.i + z__9.i * z__10.r;
		    z__1.r = z__2.r + z__8.r, z__1.i = z__2.i + z__8.i;
		    AB(k-j+1,j).r = z__1.r, AB(k-j+1,j).i = z__1.i;
/* L270: */
		}
/* Computing MAX */
		i__2 = 1, i__3 = i - *ka;
		i__1 = i - kbt - 1;
		for (j = max(1,i-*ka); j <= i-kbt-1; ++j) {
		    i__2 = k - j + 1 + j * ab_dim1;
		    i__3 = k - j + 1 + j * ab_dim1;
		    d_cnjg(&z__3, &BB(i-k+1,k));
		    i__5 = i - j + 1 + j * ab_dim1;
		    z__2.r = z__3.r * AB(i-j+1,j).r - z__3.i * AB(i-j+1,j).i, 
			    z__2.i = z__3.r * AB(i-j+1,j).i + z__3.i * AB(i-j+1,j)
			    .r;
		    z__1.r = AB(k-j+1,j).r - z__2.r, z__1.i = AB(k-j+1,j).i - 
			    z__2.i;
		    AB(k-j+1,j).r = z__1.r, AB(k-j+1,j).i = z__1.i;
/* L280: */
		}
/* L290: */
	    }
	    i__4 = i1;
	    for (j = i; j <= i1; ++j) {
/* Computing MAX */
		i__1 = j - *ka, i__2 = i - kbt;
		i__3 = i - 1;
		for (k = max(j-*ka,i-kbt); k <= i-1; ++k) {
		    i__1 = j - k + 1 + k * ab_dim1;
		    i__2 = j - k + 1 + k * ab_dim1;
		    i__5 = i - k + 1 + k * bb_dim1;
		    i__6 = j - i + 1 + i * ab_dim1;
		    z__2.r = BB(i-k+1,k).r * AB(j-i+1,i).r - BB(i-k+1,k).i * AB(j-i+1,i)
			    .i, z__2.i = BB(i-k+1,k).r * AB(j-i+1,i).i + BB(i-k+1,k).i 
			    * AB(j-i+1,i).r;
		    z__1.r = AB(j-k+1,k).r - z__2.r, z__1.i = AB(j-k+1,k).i - 
			    z__2.i;
		    AB(j-k+1,k).r = z__1.r, AB(j-k+1,k).i = z__1.i;
/* L300: */
		}
/* L310: */
	    }

	    if (wantx) {

/*              post-multiply X by inv(S(i)) */

		i__4 = *n - m;
		d__1 = 1. / bii;
		zdscal_(&i__4, &d__1, &X(m+1,i), &c__1);
		if (kbt > 0) {
		    i__4 = *n - m;
		    z__1.r = -1., z__1.i = 0.;
		    i__3 = *ldbb - 1;
		    zgeru_(&i__4, &kbt, &z__1, &X(m+1,i), &c__1, 
			    &BB(kbt+1,i-kbt), &i__3, &X(m+1,i-kbt), ldx);
		}
	    }

/*           store a(i1,i) in RA1 for use in next loop over K */

	    i__4 = i1 - i + 1 + i * ab_dim1;
	    ra1.r = AB(i1-i+1,i).r, ra1.i = AB(i1-i+1,i).i;
	}

/*        Generate and apply vectors of rotations to chase all the   
          existing bulges KA positions down toward the bottom of the 
  
          band */

	i__4 = *kb - 1;
	for (k = 1; k <= *kb-1; ++k) {
	    if (update) {

/*              Determine the rotations which would annihilate
 the bulge   
                which has in theory just been created */

		if (i - k + *ka < *n && i - k > 1) {

/*                 generate rotation to annihilate a(i-k+k
a+1,i) */

		    zlartg_(&AB(ka1-k,i), &ra1, &RWORK(i - k + *
			    ka - m), &WORK(i - k + *ka - m), &ra);

/*                 create nonzero element a(i-k+ka+1,i-k) 
outside the   
                   band and store it in WORK(i-k) */

		    i__3 = k + 1 + (i - k) * bb_dim1;
		    z__2.r = -BB(k+1,i-k).r, z__2.i = -BB(k+1,i-k).i;
		    z__1.r = z__2.r * ra1.r - z__2.i * ra1.i, z__1.i = z__2.r 
			    * ra1.i + z__2.i * ra1.r;
		    t.r = z__1.r, t.i = z__1.i;
		    i__3 = i - k;
		    i__1 = i - k + *ka - m;
		    z__2.r = RWORK(i-k+*ka-m) * t.r, z__2.i = RWORK(i-k+*ka-m) * t.i;
		    d_cnjg(&z__4, &WORK(i - k + *ka - m));
		    i__2 = ka1 + (i - k) * ab_dim1;
		    z__3.r = z__4.r * AB(ka1,i-k).r - z__4.i * AB(ka1,i-k).i, 
			    z__3.i = z__4.r * AB(ka1,i-k).i + z__4.i * AB(ka1,i-k)
			    .r;
		    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
		    WORK(i-k).r = z__1.r, WORK(i-k).i = z__1.i;
		    i__3 = ka1 + (i - k) * ab_dim1;
		    i__1 = i - k + *ka - m;
		    z__2.r = WORK(i-k+*ka-m).r * t.r - WORK(i-k+*ka-m).i * t.i, z__2.i =
			     WORK(i-k+*ka-m).r * t.i + WORK(i-k+*ka-m).i * t.r;
		    i__2 = i - k + *ka - m;
		    i__5 = ka1 + (i - k) * ab_dim1;
		    z__3.r = RWORK(i-k+*ka-m) * AB(ka1,i-k).r, z__3.i = RWORK(i-k+*ka-m) * 
			    AB(ka1,i-k).i;
		    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
		    AB(ka1,i-k).r = z__1.r, AB(ka1,i-k).i = z__1.i;
		    ra1.r = ra.r, ra1.i = ra.i;
		}
	    }
/* Computing MAX */
	    i__3 = 1, i__1 = k - i0 + 2;
	    j2 = i - k - 1 + max(i__3,i__1) * ka1;
	    nr = (*n - j2 + *ka) / ka1;
	    j1 = j2 + (nr - 1) * ka1;
	    if (update) {
/* Computing MAX */
		i__3 = j2, i__1 = i + (*ka << 1) - k + 1;
		j2t = max(i__3,i__1);
	    } else {
		j2t = j2;
	    }
	    nrt = (*n - j2t + *ka) / ka1;
	    i__3 = j1;
	    i__1 = ka1;
	    for (j = j2t; ka1 < 0 ? j >= j1 : j <= j1; j += ka1) {

/*              create nonzero element a(j+1,j-ka) outside the
 band   
                and store it in WORK(j-m) */

		i__2 = j - m;
		i__5 = j - m;
		i__6 = ka1 + (j - *ka + 1) * ab_dim1;
		z__1.r = WORK(j-m).r * AB(ka1,j-*ka+1).r - WORK(j-m).i * AB(ka1,j-*ka+1)
			.i, z__1.i = WORK(j-m).r * AB(ka1,j-*ka+1).i + WORK(j-m).i 
			* AB(ka1,j-*ka+1).r;
		WORK(j-m).r = z__1.r, WORK(j-m).i = z__1.i;
		i__2 = ka1 + (j - *ka + 1) * ab_dim1;
		i__5 = j - m;
		i__6 = ka1 + (j - *ka + 1) * ab_dim1;
		z__1.r = RWORK(j-m) * AB(ka1,j-*ka+1).r, z__1.i = RWORK(j-m) * AB(ka1,j-*ka+1).i;
		AB(ka1,j-*ka+1).r = z__1.r, AB(ka1,j-*ka+1).i = z__1.i;
/* L320: */
	    }

/*           generate rotations in 1st set to annihilate elements 
which   
             have been created outside the band */

	    if (nrt > 0) {
		zlargv_(&nrt, &AB(ka1,j2t-*ka), &inca, &WORK(
			j2t - m), &ka1, &RWORK(j2t - m), &ka1);
	    }
	    if (nr > 0) {

/*              apply rotations in 1st set from the left */

		i__1 = *ka - 1;
		for (l = 1; l <= *ka-1; ++l) {
		    zlartv_(&nr, &AB(l+1,j2-l), &inca, &AB(l+2,j2-l), &inca, &RWORK(j2 - m)
			    , &WORK(j2 - m), &ka1);
/* L330: */
		}

/*              apply rotations in 1st set from both sides to 
diagonal   
                blocks */

		zlar2v_(&nr, &AB(1,j2), &AB(1,j2+1), &AB(2,j2), &inca, &RWORK(j2 - m), &
			WORK(j2 - m), &ka1);

		zlacgv_(&nr, &WORK(j2 - m), &ka1);
	    }

/*           start applying rotations in 1st set from the right */

	    i__1 = *kb - k + 1;
	    for (l = *ka - 1; l >= *kb-k+1; --l) {
		nrt = (*n - j2 + l) / ka1;
		if (nrt > 0) {
		    zlartv_(&nrt, &AB(ka1-l+1,j2), &inca, &AB(ka1-l,j2+1), &inca, &RWORK(j2 - 
			    m), &WORK(j2 - m), &ka1);
		}
/* L340: */
	    }

	    if (wantx) {

/*              post-multiply X by product of rotations in 1st
 set */

		i__1 = j1;
		i__3 = ka1;
		for (j = j2; ka1 < 0 ? j >= j1 : j <= j1; j += ka1) {
		    i__2 = *n - m;
		    zrot_(&i__2, &X(m+1,j), &c__1, &X(m+1,j+1), &c__1, &RWORK(j - m), &WORK(j - m)
			    );
/* L350: */
		}
	    }
/* L360: */
	}

	if (update) {
	    if (i2 <= *n && kbt > 0) {

/*              create nonzero element a(i-kbt+ka+1,i-kbt) out
side the   
                band and store it in WORK(i-kbt) */

		i__4 = i - kbt;
		i__3 = kbt + 1 + (i - kbt) * bb_dim1;
		z__2.r = -BB(kbt+1,i-kbt).r, z__2.i = -BB(kbt+1,i-kbt).i;
		z__1.r = z__2.r * ra1.r - z__2.i * ra1.i, z__1.i = z__2.r * 
			ra1.i + z__2.i * ra1.r;
		WORK(i-kbt).r = z__1.r, WORK(i-kbt).i = z__1.i;
	    }
	}

	for (k = *kb; k >= 1; --k) {
	    if (update) {
/* Computing MAX */
		i__4 = 2, i__3 = k - i0 + 1;
		j2 = i - k - 1 + max(i__4,i__3) * ka1;
	    } else {
/* Computing MAX */
		i__4 = 1, i__3 = k - i0 + 1;
		j2 = i - k - 1 + max(i__4,i__3) * ka1;
	    }

/*           finish applying rotations in 2nd set from the right 
*/

	    for (l = *kb - k; l >= 1; --l) {
		nrt = (*n - j2 + *ka + l) / ka1;
		if (nrt > 0) {
		    zlartv_(&nrt, &AB(ka1-l+1,j2-*ka), &
			    inca, &AB(ka1-l,j2-*ka+1), &
			    inca, &RWORK(j2 - *ka), &WORK(j2 - *ka), &ka1);
		}
/* L370: */
	    }
	    nr = (*n - j2 + *ka) / ka1;
	    j1 = j2 + (nr - 1) * ka1;
	    i__4 = j2;
	    i__3 = -ka1;
	    for (j = j1; -ka1 < 0 ? j >= j2 : j <= j2; j += -ka1) {
		i__1 = j;
		i__2 = j - *ka;
		WORK(j).r = WORK(j-*ka).r, WORK(j).i = WORK(j-*ka).i;
		RWORK(j) = RWORK(j - *ka);
/* L380: */
	    }
	    i__3 = j1;
	    i__4 = ka1;
	    for (j = j2; ka1 < 0 ? j >= j1 : j <= j1; j += ka1) {

/*              create nonzero element a(j+1,j-ka) outside the
 band   
                and store it in WORK(j) */

		i__1 = j;
		i__2 = j;
		i__5 = ka1 + (j - *ka + 1) * ab_dim1;
		z__1.r = WORK(j).r * AB(ka1,j-*ka+1).r - WORK(j).i * AB(ka1,j-*ka+1)
			.i, z__1.i = WORK(j).r * AB(ka1,j-*ka+1).i + WORK(j).i 
			* AB(ka1,j-*ka+1).r;
		WORK(j).r = z__1.r, WORK(j).i = z__1.i;
		i__1 = ka1 + (j - *ka + 1) * ab_dim1;
		i__2 = j;
		i__5 = ka1 + (j - *ka + 1) * ab_dim1;
		z__1.r = RWORK(j) * AB(ka1,j-*ka+1).r, z__1.i = RWORK(j) * AB(ka1,j-*ka+1).i;
		AB(ka1,j-*ka+1).r = z__1.r, AB(ka1,j-*ka+1).i = z__1.i;
/* L390: */
	    }
	    if (update) {
		if (i - k < *n - *ka && k <= kbt) {
		    i__4 = i - k + *ka;
		    i__3 = i - k;
		    WORK(i-k+*ka).r = WORK(i-k).r, WORK(i-k+*ka).i = WORK(i-k).i;
		}
	    }
/* L400: */
	}

	for (k = *kb; k >= 1; --k) {
/* Computing MAX */
	    i__4 = 1, i__3 = k - i0 + 1;
	    j2 = i - k - 1 + max(i__4,i__3) * ka1;
	    nr = (*n - j2 + *ka) / ka1;
	    j1 = j2 + (nr - 1) * ka1;
	    if (nr > 0) {

/*              generate rotations in 2nd set to annihilate el
ements   
                which have been created outside the band */

		zlargv_(&nr, &AB(ka1,j2-*ka), &inca, &WORK(j2)
			, &ka1, &RWORK(j2), &ka1);

/*              apply rotations in 2nd set from the left */

		i__4 = *ka - 1;
		for (l = 1; l <= *ka-1; ++l) {
		    zlartv_(&nr, &AB(l+1,j2-l), &inca, &AB(l+2,j2-l), &inca, &RWORK(j2), &
			    WORK(j2), &ka1);
/* L410: */
		}

/*              apply rotations in 2nd set from both sides to 
diagonal   
                blocks */

		zlar2v_(&nr, &AB(1,j2), &AB(1,j2+1), &AB(2,j2), &inca, &RWORK(j2), &WORK(
			j2), &ka1);

		zlacgv_(&nr, &WORK(j2), &ka1);
	    }

/*           start applying rotations in 2nd set from the right */

	    i__4 = *kb - k + 1;
	    for (l = *ka - 1; l >= *kb-k+1; --l) {
		nrt = (*n - j2 + l) / ka1;
		if (nrt > 0) {
		    zlartv_(&nrt, &AB(ka1-l+1,j2), &inca, &AB(ka1-l,j2+1), &inca, &RWORK(j2), 
			    &WORK(j2), &ka1);
		}
/* L420: */
	    }

	    if (wantx) {

/*              post-multiply X by product of rotations in 2nd
 set */

		i__4 = j1;
		i__3 = ka1;
		for (j = j2; ka1 < 0 ? j >= j1 : j <= j1; j += ka1) {
		    i__1 = *n - m;
		    zrot_(&i__1, &X(m+1,j), &c__1, &X(m+1,j+1), &c__1, &RWORK(j), &WORK(j));
/* L430: */
		}
	    }
/* L440: */
	}

	i__3 = *kb - 1;
	for (k = 1; k <= *kb-1; ++k) {
/* Computing MAX */
	    i__4 = 1, i__1 = k - i0 + 2;
	    j2 = i - k - 1 + max(i__4,i__1) * ka1;

/*           finish applying rotations in 1st set from the right 
*/

	    for (l = *kb - k; l >= 1; --l) {
		nrt = (*n - j2 + l) / ka1;
		if (nrt > 0) {
		    zlartv_(&nrt, &AB(ka1-l+1,j2), &inca, &AB(ka1-l,j2+1), &inca, &RWORK(j2 - 
			    m), &WORK(j2 - m), &ka1);
		}
/* L450: */
	    }
/* L460: */
	}

	i__3 = i2 + *ka;
	for (j = *n - 1; j >= i2+*ka; --j) {
	    RWORK(j - m) = RWORK(j - *ka - m);
	    i__4 = j - m;
	    i__1 = j - *ka - m;
	    WORK(j-m).r = WORK(j-*ka-m).r, WORK(j-m).i = WORK(j-*ka-m).i;
/* L470: */
	}

    }

    goto L10;

L480:

/*     **************************** Phase 2 ***************************** 
  

       The logical structure of this phase is:   

       UPDATE = .TRUE.   
       DO I = 1, M   
          use S(i) to update A and create a new bulge   
          apply rotations to push all bulges KA positions upward   
       END DO   
       UPDATE = .FALSE.   
       DO I = M - KA - 1, 2, -1   
          apply rotations to push all bulges KA positions upward   
       END DO   

       To avoid duplicating code, the two loops are merged. */

    update = TRUE_;
    i = 0;
L490:
    if (update) {
	++i;
/* Computing MIN */
	i__3 = *kb, i__4 = m - i;
	kbt = min(i__3,i__4);
	i0 = i + 1;
/* Computing MAX */
	i__3 = 1, i__4 = i - *ka;
	i1 = max(i__3,i__4);
	i2 = i + kbt - ka1;
	if (i > m) {
	    update = FALSE_;
	    --i;
	    i0 = m + 1;
	    if (*ka == 0) {
		return 0;
	    }
	    goto L490;
	}
    } else {
	i -= *ka;
	if (i < 2) {
	    return 0;
	}
    }

    if (i < m - kbt) {
	nx = m;
    } else {
	nx = *n;
    }

    if (upper) {

/*        Transform A, working with the upper triangle */

	if (update) {

/*           Form  inv(S(i))**H * A * inv(S(i)) */

	    i__3 = kb1 + i * bb_dim1;
	    bii = BB(kb1,i).r;
	    i__3 = ka1 + i * ab_dim1;
	    i__4 = ka1 + i * ab_dim1;
	    d__1 = AB(ka1,i).r / bii / bii;
	    AB(ka1,i).r = d__1, AB(ka1,i).i = 0.;
	    i__3 = i - 1;
	    for (j = i1; j <= i-1; ++j) {
		i__4 = j - i + ka1 + i * ab_dim1;
		i__1 = j - i + ka1 + i * ab_dim1;
		z__1.r = AB(j-i+ka1,i).r / bii, z__1.i = AB(j-i+ka1,i).i / bii;
		AB(j-i+ka1,i).r = z__1.r, AB(j-i+ka1,i).i = z__1.i;
/* L500: */
	    }
/* Computing MIN */
	    i__4 = *n, i__1 = i + *ka;
	    i__3 = min(i__4,i__1);
	    for (j = i + 1; j <= min(*n,i+*ka); ++j) {
		i__4 = i - j + ka1 + j * ab_dim1;
		i__1 = i - j + ka1 + j * ab_dim1;
		z__1.r = AB(i-j+ka1,j).r / bii, z__1.i = AB(i-j+ka1,j).i / bii;
		AB(i-j+ka1,j).r = z__1.r, AB(i-j+ka1,j).i = z__1.i;
/* L510: */
	    }
	    i__3 = i + kbt;
	    for (k = i + 1; k <= i+kbt; ++k) {
		i__4 = i + kbt;
		for (j = k; j <= i+kbt; ++j) {
		    i__1 = k - j + ka1 + j * ab_dim1;
		    i__2 = k - j + ka1 + j * ab_dim1;
		    i__5 = i - j + kb1 + j * bb_dim1;
		    d_cnjg(&z__5, &AB(i-k+ka1,k));
		    z__4.r = BB(i-j+kb1,j).r * z__5.r - BB(i-j+kb1,j).i * z__5.i, 
			    z__4.i = BB(i-j+kb1,j).r * z__5.i + BB(i-j+kb1,j).i * 
			    z__5.r;
		    z__3.r = AB(k-j+ka1,j).r - z__4.r, z__3.i = AB(k-j+ka1,j).i - 
			    z__4.i;
		    d_cnjg(&z__7, &BB(i-k+kb1,k));
		    i__6 = i - j + ka1 + j * ab_dim1;
		    z__6.r = z__7.r * AB(i-j+ka1,j).r - z__7.i * AB(i-j+ka1,j).i, 
			    z__6.i = z__7.r * AB(i-j+ka1,j).i + z__7.i * AB(i-j+ka1,j)
			    .r;
		    z__2.r = z__3.r - z__6.r, z__2.i = z__3.i - z__6.i;
		    i__7 = ka1 + i * ab_dim1;
		    d__1 = AB(ka1,i).r;
		    i__8 = i - j + kb1 + j * bb_dim1;
		    z__9.r = d__1 * BB(i-j+kb1,j).r, z__9.i = d__1 * BB(i-j+kb1,j).i;
		    d_cnjg(&z__10, &BB(i-k+kb1,k));
		    z__8.r = z__9.r * z__10.r - z__9.i * z__10.i, z__8.i = 
			    z__9.r * z__10.i + z__9.i * z__10.r;
		    z__1.r = z__2.r + z__8.r, z__1.i = z__2.i + z__8.i;
		    AB(k-j+ka1,j).r = z__1.r, AB(k-j+ka1,j).i = z__1.i;
/* L520: */
		}
/* Computing MIN */
		i__1 = *n, i__2 = i + *ka;
		i__4 = min(i__1,i__2);
		for (j = i + kbt + 1; j <= min(*n,i+*ka); ++j) {
		    i__1 = k - j + ka1 + j * ab_dim1;
		    i__2 = k - j + ka1 + j * ab_dim1;
		    d_cnjg(&z__3, &BB(i-k+kb1,k));
		    i__5 = i - j + ka1 + j * ab_dim1;
		    z__2.r = z__3.r * AB(i-j+ka1,j).r - z__3.i * AB(i-j+ka1,j).i, 
			    z__2.i = z__3.r * AB(i-j+ka1,j).i + z__3.i * AB(i-j+ka1,j)
			    .r;
		    z__1.r = AB(k-j+ka1,j).r - z__2.r, z__1.i = AB(k-j+ka1,j).i - 
			    z__2.i;
		    AB(k-j+ka1,j).r = z__1.r, AB(k-j+ka1,j).i = z__1.i;
/* L530: */
		}
/* L540: */
	    }
	    i__3 = i;
	    for (j = i1; j <= i; ++j) {
/* Computing MIN */
		i__1 = j + *ka, i__2 = i + kbt;
		i__4 = min(i__1,i__2);
		for (k = i + 1; k <= min(j+*ka,i+kbt); ++k) {
		    i__1 = j - k + ka1 + k * ab_dim1;
		    i__2 = j - k + ka1 + k * ab_dim1;
		    i__5 = i - k + kb1 + k * bb_dim1;
		    i__6 = j - i + ka1 + i * ab_dim1;
		    z__2.r = BB(i-k+kb1,k).r * AB(j-i+ka1,i).r - BB(i-k+kb1,k).i * AB(j-i+ka1,i)
			    .i, z__2.i = BB(i-k+kb1,k).r * AB(j-i+ka1,i).i + BB(i-k+kb1,k).i 
			    * AB(j-i+ka1,i).r;
		    z__1.r = AB(j-k+ka1,k).r - z__2.r, z__1.i = AB(j-k+ka1,k).i - 
			    z__2.i;
		    AB(j-k+ka1,k).r = z__1.r, AB(j-k+ka1,k).i = z__1.i;
/* L550: */
		}
/* L560: */
	    }

	    if (wantx) {

/*              post-multiply X by inv(S(i)) */

		d__1 = 1. / bii;
		zdscal_(&nx, &d__1, &X(1,i), &c__1);
		if (kbt > 0) {
		    z__1.r = -1., z__1.i = 0.;
		    i__3 = *ldbb - 1;
		    zgeru_(&nx, &kbt, &z__1, &X(1,i), &c__1, &BB(*kb,i+1), &i__3, &X(1,i+1), ldx);
		}
	    }

/*           store a(i1,i) in RA1 for use in next loop over K */

	    i__3 = i1 - i + ka1 + i * ab_dim1;
	    ra1.r = AB(i1-i+ka1,i).r, ra1.i = AB(i1-i+ka1,i).i;
	}

/*        Generate and apply vectors of rotations to chase all the   
          existing bulges KA positions up toward the top of the band 
*/

	i__3 = *kb - 1;
	for (k = 1; k <= *kb-1; ++k) {
	    if (update) {

/*              Determine the rotations which would annihilate
 the bulge   
                which has in theory just been created */

		if (i + k - ka1 > 0 && i + k < m) {

/*                 generate rotation to annihilate a(i+k-k
a-1,i) */

		    zlartg_(&AB(k+1,i), &ra1, &RWORK(i + k - *
			    ka), &WORK(i + k - *ka), &ra);

/*                 create nonzero element a(i+k-ka-1,i+k) 
outside the   
                   band and store it in WORK(m-kb+i+k) */

		    i__4 = kb1 - k + (i + k) * bb_dim1;
		    z__2.r = -BB(kb1-k,i+k).r, z__2.i = -BB(kb1-k,i+k).i;
		    z__1.r = z__2.r * ra1.r - z__2.i * ra1.i, z__1.i = z__2.r 
			    * ra1.i + z__2.i * ra1.r;
		    t.r = z__1.r, t.i = z__1.i;
		    i__4 = m - *kb + i + k;
		    i__1 = i + k - *ka;
		    z__2.r = RWORK(i+k-*ka) * t.r, z__2.i = RWORK(i+k-*ka) * t.i;
		    d_cnjg(&z__4, &WORK(i + k - *ka));
		    i__2 = (i + k) * ab_dim1 + 1;
		    z__3.r = z__4.r * AB(1,i+k).r - z__4.i * AB(1,i+k).i, 
			    z__3.i = z__4.r * AB(1,i+k).i + z__4.i * AB(1,i+k)
			    .r;
		    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
		    WORK(m-*kb+i+k).r = z__1.r, WORK(m-*kb+i+k).i = z__1.i;
		    i__4 = (i + k) * ab_dim1 + 1;
		    i__1 = i + k - *ka;
		    z__2.r = WORK(i+k-*ka).r * t.r - WORK(i+k-*ka).i * t.i, z__2.i =
			     WORK(i+k-*ka).r * t.i + WORK(i+k-*ka).i * t.r;
		    i__2 = i + k - *ka;
		    i__5 = (i + k) * ab_dim1 + 1;
		    z__3.r = RWORK(i+k-*ka) * AB(1,i+k).r, z__3.i = RWORK(i+k-*ka) * 
			    AB(1,i+k).i;
		    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
		    AB(1,i+k).r = z__1.r, AB(1,i+k).i = z__1.i;
		    ra1.r = ra.r, ra1.i = ra.i;
		}
	    }
/* Computing MAX */
	    i__4 = 1, i__1 = k + i0 - m + 1;
	    j2 = i + k + 1 - max(i__4,i__1) * ka1;
	    nr = (j2 + *ka - 1) / ka1;
	    j1 = j2 - (nr - 1) * ka1;
	    if (update) {
/* Computing MIN */
		i__4 = j2, i__1 = i - (*ka << 1) + k - 1;
		j2t = min(i__4,i__1);
	    } else {
		j2t = j2;
	    }
	    nrt = (j2t + *ka - 1) / ka1;
	    i__4 = j2t;
	    i__1 = ka1;
	    for (j = j1; ka1 < 0 ? j >= j2t : j <= j2t; j += ka1) {

/*              create nonzero element a(j-1,j+ka) outside the
 band   
                and store it in WORK(j) */

		i__2 = j;
		i__5 = j;
		i__6 = (j + *ka - 1) * ab_dim1 + 1;
		z__1.r = WORK(j).r * AB(1,j+*ka-1).r - WORK(j).i * AB(1,j+*ka-1)
			.i, z__1.i = WORK(j).r * AB(1,j+*ka-1).i + WORK(j).i 
			* AB(1,j+*ka-1).r;
		WORK(j).r = z__1.r, WORK(j).i = z__1.i;
		i__2 = (j + *ka - 1) * ab_dim1 + 1;
		i__5 = j;
		i__6 = (j + *ka - 1) * ab_dim1 + 1;
		z__1.r = RWORK(j) * AB(1,j+*ka-1).r, z__1.i = RWORK(j) * AB(1,j+*ka-1).i;
		AB(1,j+*ka-1).r = z__1.r, AB(1,j+*ka-1).i = z__1.i;
/* L570: */
	    }

/*           generate rotations in 1st set to annihilate elements 
which   
             have been created outside the band */

	    if (nrt > 0) {
		zlargv_(&nrt, &AB(1,j1+*ka), &inca, &WORK(j1),
			 &ka1, &RWORK(j1), &ka1);
	    }
	    if (nr > 0) {

/*              apply rotations in 1st set from the left */

		i__1 = *ka - 1;
		for (l = 1; l <= *ka-1; ++l) {
		    zlartv_(&nr, &AB(ka1-l,j1+l), &inca, &
			    AB(*ka-l,j1+l), &inca, &RWORK(
			    j1), &WORK(j1), &ka1);
/* L580: */
		}

/*              apply rotations in 1st set from both sides to 
diagonal   
                blocks */

		zlar2v_(&nr, &AB(ka1,j1), &AB(ka1,j1-1), &AB(*ka,j1), &inca, &RWORK(j1), 
			&WORK(j1), &ka1);

		zlacgv_(&nr, &WORK(j1), &ka1);
	    }

/*           start applying rotations in 1st set from the right */

	    i__1 = *kb - k + 1;
	    for (l = *ka - 1; l >= *kb-k+1; --l) {
		nrt = (j2 + l - 1) / ka1;
		j1t = j2 - (nrt - 1) * ka1;
		if (nrt > 0) {
		    zlartv_(&nrt, &AB(l,j1t), &inca, &AB(l+1,j1t-1), &inca, &RWORK(j1t), &WORK(
			    j1t), &ka1);
		}
/* L590: */
	    }

	    if (wantx) {

/*              post-multiply X by product of rotations in 1st
 set */

		i__1 = j2;
		i__4 = ka1;
		for (j = j1; ka1 < 0 ? j >= j2 : j <= j2; j += ka1) {
		    zrot_(&nx, &X(1,j), &c__1, &X(1,j-1), &c__1, &RWORK(j), &WORK(j));
/* L600: */
		}
	    }
/* L610: */
	}

	if (update) {
	    if (i2 > 0 && kbt > 0) {

/*              create nonzero element a(i+kbt-ka-1,i+kbt) out
side the   
                band and store it in WORK(m-kb+i+kbt) */

		i__3 = m - *kb + i + kbt;
		i__4 = kb1 - kbt + (i + kbt) * bb_dim1;
		z__2.r = -BB(kb1-kbt,i+kbt).r, z__2.i = -BB(kb1-kbt,i+kbt).i;
		z__1.r = z__2.r * ra1.r - z__2.i * ra1.i, z__1.i = z__2.r * 
			ra1.i + z__2.i * ra1.r;
		WORK(m-*kb+i+kbt).r = z__1.r, WORK(m-*kb+i+kbt).i = z__1.i;
	    }
	}

	for (k = *kb; k >= 1; --k) {
	    if (update) {
/* Computing MAX */
		i__3 = 2, i__4 = k + i0 - m;
		j2 = i + k + 1 - max(i__3,i__4) * ka1;
	    } else {
/* Computing MAX */
		i__3 = 1, i__4 = k + i0 - m;
		j2 = i + k + 1 - max(i__3,i__4) * ka1;
	    }

/*           finish applying rotations in 2nd set from the right 
*/

	    for (l = *kb - k; l >= 1; --l) {
		nrt = (j2 + *ka + l - 1) / ka1;
		j1t = j2 - (nrt - 1) * ka1;
		if (nrt > 0) {
		    zlartv_(&nrt, &AB(l,j1t+*ka), &inca, &AB(l+1,j1t+*ka-1), &inca, &RWORK(
			    m - *kb + j1t + *ka), &WORK(m - *kb + j1t + *ka), 
			    &ka1);
		}
/* L620: */
	    }
	    nr = (j2 + *ka - 1) / ka1;
	    j1 = j2 - (nr - 1) * ka1;
	    i__3 = j2;
	    i__4 = ka1;
	    for (j = j1; ka1 < 0 ? j >= j2 : j <= j2; j += ka1) {
		i__1 = m - *kb + j;
		i__2 = m - *kb + j + *ka;
		WORK(m-*kb+j).r = WORK(m-*kb+j+*ka).r, WORK(m-*kb+j).i = WORK(m-*kb+j+*ka).i;
		RWORK(m - *kb + j) = RWORK(m - *kb + j + *ka);
/* L630: */
	    }
	    i__4 = j2;
	    i__3 = ka1;
	    for (j = j1; ka1 < 0 ? j >= j2 : j <= j2; j += ka1) {

/*              create nonzero element a(j-1,j+ka) outside the
 band   
                and store it in WORK(m-kb+j) */

		i__1 = m - *kb + j;
		i__2 = m - *kb + j;
		i__5 = (j + *ka - 1) * ab_dim1 + 1;
		z__1.r = WORK(m-*kb+j).r * AB(1,j+*ka-1).r - WORK(m-*kb+j).i * AB(1,j+*ka-1)
			.i, z__1.i = WORK(m-*kb+j).r * AB(1,j+*ka-1).i + WORK(m-*kb+j).i 
			* AB(1,j+*ka-1).r;
		WORK(m-*kb+j).r = z__1.r, WORK(m-*kb+j).i = z__1.i;
		i__1 = (j + *ka - 1) * ab_dim1 + 1;
		i__2 = m - *kb + j;
		i__5 = (j + *ka - 1) * ab_dim1 + 1;
		z__1.r = RWORK(m-*kb+j) * AB(1,j+*ka-1).r, z__1.i = RWORK(m-*kb+j) * AB(1,j+*ka-1).i;
		AB(1,j+*ka-1).r = z__1.r, AB(1,j+*ka-1).i = z__1.i;
/* L640: */
	    }
	    if (update) {
		if (i + k > ka1 && k <= kbt) {
		    i__3 = m - *kb + i + k - *ka;
		    i__4 = m - *kb + i + k;
		    WORK(m-*kb+i+k-*ka).r = WORK(m-*kb+i+k).r, WORK(m-*kb+i+k-*ka).i = WORK(m-*kb+i+k).i;
		}
	    }
/* L650: */
	}

	for (k = *kb; k >= 1; --k) {
/* Computing MAX */
	    i__3 = 1, i__4 = k + i0 - m;
	    j2 = i + k + 1 - max(i__3,i__4) * ka1;
	    nr = (j2 + *ka - 1) / ka1;
	    j1 = j2 - (nr - 1) * ka1;
	    if (nr > 0) {

/*              generate rotations in 2nd set to annihilate el
ements   
                which have been created outside the band */

		zlargv_(&nr, &AB(1,j1+*ka), &inca, &WORK(m - *
			kb + j1), &ka1, &RWORK(m - *kb + j1), &ka1);

/*              apply rotations in 2nd set from the left */

		i__3 = *ka - 1;
		for (l = 1; l <= *ka-1; ++l) {
		    zlartv_(&nr, &AB(ka1-l,j1+l), &inca, &
			    AB(*ka-l,j1+l), &inca, &RWORK(m 
			    - *kb + j1), &WORK(m - *kb + j1), &ka1);
/* L660: */
		}

/*              apply rotations in 2nd set from both sides to 
diagonal   
                blocks */

		zlar2v_(&nr, &AB(ka1,j1), &AB(ka1,j1-1), &AB(*ka,j1), &inca, &RWORK(m - *
			kb + j1), &WORK(m - *kb + j1), &ka1);

		zlacgv_(&nr, &WORK(m - *kb + j1), &ka1);
	    }

/*           start applying rotations in 2nd set from the right */

	    i__3 = *kb - k + 1;
	    for (l = *ka - 1; l >= *kb-k+1; --l) {
		nrt = (j2 + l - 1) / ka1;
		j1t = j2 - (nrt - 1) * ka1;
		if (nrt > 0) {
		    zlartv_(&nrt, &AB(l,j1t), &inca, &AB(l+1,j1t-1), &inca, &RWORK(m - *kb + j1t),
			     &WORK(m - *kb + j1t), &ka1);
		}
/* L670: */
	    }

	    if (wantx) {

/*              post-multiply X by product of rotations in 2nd
 set */

		i__3 = j2;
		i__4 = ka1;
		for (j = j1; ka1 < 0 ? j >= j2 : j <= j2; j += ka1) {
		    zrot_(&nx, &X(1,j), &c__1, &X(1,j-1), &c__1, &RWORK(m - *kb + j), &WORK(m - *kb + 
			    j));
/* L680: */
		}
	    }
/* L690: */
	}

	i__4 = *kb - 1;
	for (k = 1; k <= *kb-1; ++k) {
/* Computing MAX */
	    i__3 = 1, i__1 = k + i0 - m + 1;
	    j2 = i + k + 1 - max(i__3,i__1) * ka1;

/*           finish applying rotations in 1st set from the right 
*/

	    for (l = *kb - k; l >= 1; --l) {
		nrt = (j2 + l - 1) / ka1;
		j1t = j2 - (nrt - 1) * ka1;
		if (nrt > 0) {
		    zlartv_(&nrt, &AB(l,j1t), &inca, &AB(l+1,j1t-1), &inca, &RWORK(j1t), &WORK(
			    j1t), &ka1);
		}
/* L700: */
	    }
/* L710: */
	}

	i__4 = i2 - *ka;
	for (j = 2; j <= i2-*ka; ++j) {
	    RWORK(j) = RWORK(j + *ka);
	    i__3 = j;
	    i__1 = j + *ka;
	    WORK(j).r = WORK(j+*ka).r, WORK(j).i = WORK(j+*ka).i;
/* L720: */
	}

    } else {

/*        Transform A, working with the lower triangle */

	if (update) {

/*           Form  inv(S(i))**H * A * inv(S(i)) */

	    i__4 = i * bb_dim1 + 1;
	    bii = BB(1,i).r;
	    i__4 = i * ab_dim1 + 1;
	    i__3 = i * ab_dim1 + 1;
	    d__1 = AB(1,i).r / bii / bii;
	    AB(1,i).r = d__1, AB(1,i).i = 0.;
	    i__4 = i - 1;
	    for (j = i1; j <= i-1; ++j) {
		i__3 = i - j + 1 + j * ab_dim1;
		i__1 = i - j + 1 + j * ab_dim1;
		z__1.r = AB(i-j+1,j).r / bii, z__1.i = AB(i-j+1,j).i / bii;
		AB(i-j+1,j).r = z__1.r, AB(i-j+1,j).i = z__1.i;
/* L730: */
	    }
/* Computing MIN */
	    i__3 = *n, i__1 = i + *ka;
	    i__4 = min(i__3,i__1);
	    for (j = i + 1; j <= min(*n,i+*ka); ++j) {
		i__3 = j - i + 1 + i * ab_dim1;
		i__1 = j - i + 1 + i * ab_dim1;
		z__1.r = AB(j-i+1,i).r / bii, z__1.i = AB(j-i+1,i).i / bii;
		AB(j-i+1,i).r = z__1.r, AB(j-i+1,i).i = z__1.i;
/* L740: */
	    }
	    i__4 = i + kbt;
	    for (k = i + 1; k <= i+kbt; ++k) {
		i__3 = i + kbt;
		for (j = k; j <= i+kbt; ++j) {
		    i__1 = j - k + 1 + k * ab_dim1;
		    i__2 = j - k + 1 + k * ab_dim1;
		    i__5 = j - i + 1 + i * bb_dim1;
		    d_cnjg(&z__5, &AB(k-i+1,i));
		    z__4.r = BB(j-i+1,i).r * z__5.r - BB(j-i+1,i).i * z__5.i, 
			    z__4.i = BB(j-i+1,i).r * z__5.i + BB(j-i+1,i).i * 
			    z__5.r;
		    z__3.r = AB(j-k+1,k).r - z__4.r, z__3.i = AB(j-k+1,k).i - 
			    z__4.i;
		    d_cnjg(&z__7, &BB(k-i+1,i));
		    i__6 = j - i + 1 + i * ab_dim1;
		    z__6.r = z__7.r * AB(j-i+1,i).r - z__7.i * AB(j-i+1,i).i, 
			    z__6.i = z__7.r * AB(j-i+1,i).i + z__7.i * AB(j-i+1,i)
			    .r;
		    z__2.r = z__3.r - z__6.r, z__2.i = z__3.i - z__6.i;
		    i__7 = i * ab_dim1 + 1;
		    d__1 = AB(1,i).r;
		    i__8 = j - i + 1 + i * bb_dim1;
		    z__9.r = d__1 * BB(j-i+1,i).r, z__9.i = d__1 * BB(j-i+1,i).i;
		    d_cnjg(&z__10, &BB(k-i+1,i));
		    z__8.r = z__9.r * z__10.r - z__9.i * z__10.i, z__8.i = 
			    z__9.r * z__10.i + z__9.i * z__10.r;
		    z__1.r = z__2.r + z__8.r, z__1.i = z__2.i + z__8.i;
		    AB(j-k+1,k).r = z__1.r, AB(j-k+1,k).i = z__1.i;
/* L750: */
		}
/* Computing MIN */
		i__1 = *n, i__2 = i + *ka;
		i__3 = min(i__1,i__2);
		for (j = i + kbt + 1; j <= min(*n,i+*ka); ++j) {
		    i__1 = j - k + 1 + k * ab_dim1;
		    i__2 = j - k + 1 + k * ab_dim1;
		    d_cnjg(&z__3, &BB(k-i+1,i));
		    i__5 = j - i + 1 + i * ab_dim1;
		    z__2.r = z__3.r * AB(j-i+1,i).r - z__3.i * AB(j-i+1,i).i, 
			    z__2.i = z__3.r * AB(j-i+1,i).i + z__3.i * AB(j-i+1,i)
			    .r;
		    z__1.r = AB(j-k+1,k).r - z__2.r, z__1.i = AB(j-k+1,k).i - 
			    z__2.i;
		    AB(j-k+1,k).r = z__1.r, AB(j-k+1,k).i = z__1.i;
/* L760: */
		}
/* L770: */
	    }
	    i__4 = i;
	    for (j = i1; j <= i; ++j) {
/* Computing MIN */
		i__1 = j + *ka, i__2 = i + kbt;
		i__3 = min(i__1,i__2);
		for (k = i + 1; k <= min(j+*ka,i+kbt); ++k) {
		    i__1 = k - j + 1 + j * ab_dim1;
		    i__2 = k - j + 1 + j * ab_dim1;
		    i__5 = k - i + 1 + i * bb_dim1;
		    i__6 = i - j + 1 + j * ab_dim1;
		    z__2.r = BB(k-i+1,i).r * AB(i-j+1,j).r - BB(k-i+1,i).i * AB(i-j+1,j)
			    .i, z__2.i = BB(k-i+1,i).r * AB(i-j+1,j).i + BB(k-i+1,i).i 
			    * AB(i-j+1,j).r;
		    z__1.r = AB(k-j+1,j).r - z__2.r, z__1.i = AB(k-j+1,j).i - 
			    z__2.i;
		    AB(k-j+1,j).r = z__1.r, AB(k-j+1,j).i = z__1.i;
/* L780: */
		}
/* L790: */
	    }

	    if (wantx) {

/*              post-multiply X by inv(S(i)) */

		d__1 = 1. / bii;
		zdscal_(&nx, &d__1, &X(1,i), &c__1);
		if (kbt > 0) {
		    z__1.r = -1., z__1.i = 0.;
		    zgerc_(&nx, &kbt, &z__1, &X(1,i), &c__1, &BB(2,i), &c__1, &X(1,i+1), 
			    ldx);
		}
	    }

/*           store a(i,i1) in RA1 for use in next loop over K */

	    i__4 = i - i1 + 1 + i1 * ab_dim1;
	    ra1.r = AB(i-i1+1,i1).r, ra1.i = AB(i-i1+1,i1).i;
	}

/*        Generate and apply vectors of rotations to chase all the   
          existing bulges KA positions up toward the top of the band 
*/

	i__4 = *kb - 1;
	for (k = 1; k <= *kb-1; ++k) {
	    if (update) {

/*              Determine the rotations which would annihilate
 the bulge   
                which has in theory just been created */

		if (i + k - ka1 > 0 && i + k < m) {

/*                 generate rotation to annihilate a(i,i+k
-ka-1) */

		    zlartg_(&AB(ka1-k,i+k-*ka), &ra1, &
			    RWORK(i + k - *ka), &WORK(i + k - *ka), &ra);

/*                 create nonzero element a(i+k,i+k-ka-1) 
outside the   
                   band and store it in WORK(m-kb+i+k) */

		    i__3 = k + 1 + i * bb_dim1;
		    z__2.r = -BB(k+1,i).r, z__2.i = -BB(k+1,i).i;
		    z__1.r = z__2.r * ra1.r - z__2.i * ra1.i, z__1.i = z__2.r 
			    * ra1.i + z__2.i * ra1.r;
		    t.r = z__1.r, t.i = z__1.i;
		    i__3 = m - *kb + i + k;
		    i__1 = i + k - *ka;
		    z__2.r = RWORK(i+k-*ka) * t.r, z__2.i = RWORK(i+k-*ka) * t.i;
		    d_cnjg(&z__4, &WORK(i + k - *ka));
		    i__2 = ka1 + (i + k - *ka) * ab_dim1;
		    z__3.r = z__4.r * AB(ka1,i+k-*ka).r - z__4.i * AB(ka1,i+k-*ka).i, 
			    z__3.i = z__4.r * AB(ka1,i+k-*ka).i + z__4.i * AB(ka1,i+k-*ka)
			    .r;
		    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
		    WORK(m-*kb+i+k).r = z__1.r, WORK(m-*kb+i+k).i = z__1.i;
		    i__3 = ka1 + (i + k - *ka) * ab_dim1;
		    i__1 = i + k - *ka;
		    z__2.r = WORK(i+k-*ka).r * t.r - WORK(i+k-*ka).i * t.i, z__2.i =
			     WORK(i+k-*ka).r * t.i + WORK(i+k-*ka).i * t.r;
		    i__2 = i + k - *ka;
		    i__5 = ka1 + (i + k - *ka) * ab_dim1;
		    z__3.r = RWORK(i+k-*ka) * AB(ka1,i+k-*ka).r, z__3.i = RWORK(i+k-*ka) * 
			    AB(ka1,i+k-*ka).i;
		    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
		    AB(ka1,i+k-*ka).r = z__1.r, AB(ka1,i+k-*ka).i = z__1.i;
		    ra1.r = ra.r, ra1.i = ra.i;
		}
	    }
/* Computing MAX */
	    i__3 = 1, i__1 = k + i0 - m + 1;
	    j2 = i + k + 1 - max(i__3,i__1) * ka1;
	    nr = (j2 + *ka - 1) / ka1;
	    j1 = j2 - (nr - 1) * ka1;
	    if (update) {
/* Computing MIN */
		i__3 = j2, i__1 = i - (*ka << 1) + k - 1;
		j2t = min(i__3,i__1);
	    } else {
		j2t = j2;
	    }
	    nrt = (j2t + *ka - 1) / ka1;
	    i__3 = j2t;
	    i__1 = ka1;
	    for (j = j1; ka1 < 0 ? j >= j2t : j <= j2t; j += ka1) {

/*              create nonzero element a(j+ka,j-1) outside the
 band   
                and store it in WORK(j) */

		i__2 = j;
		i__5 = j;
		i__6 = ka1 + (j - 1) * ab_dim1;
		z__1.r = WORK(j).r * AB(ka1,j-1).r - WORK(j).i * AB(ka1,j-1)
			.i, z__1.i = WORK(j).r * AB(ka1,j-1).i + WORK(j).i 
			* AB(ka1,j-1).r;
		WORK(j).r = z__1.r, WORK(j).i = z__1.i;
		i__2 = ka1 + (j - 1) * ab_dim1;
		i__5 = j;
		i__6 = ka1 + (j - 1) * ab_dim1;
		z__1.r = RWORK(j) * AB(ka1,j-1).r, z__1.i = RWORK(j) * AB(ka1,j-1).i;
		AB(ka1,j-1).r = z__1.r, AB(ka1,j-1).i = z__1.i;
/* L800: */
	    }

/*           generate rotations in 1st set to annihilate elements 
which   
             have been created outside the band */

	    if (nrt > 0) {
		zlargv_(&nrt, &AB(ka1,j1), &inca, &WORK(j1), &ka1,
			 &RWORK(j1), &ka1);
	    }
	    if (nr > 0) {

/*              apply rotations in 1st set from the right */

		i__1 = *ka - 1;
		for (l = 1; l <= *ka-1; ++l) {
		    zlartv_(&nr, &AB(l+1,j1), &inca, &AB(l+2,j1-1), &inca, &RWORK(j1), &WORK(
			    j1), &ka1);
/* L810: */
		}

/*              apply rotations in 1st set from both sides to 
diagonal   
                blocks */

		zlar2v_(&nr, &AB(1,j1), &AB(1,j1-1), &AB(2,j1-1), &inca, &RWORK(j1), &
			WORK(j1), &ka1);

		zlacgv_(&nr, &WORK(j1), &ka1);
	    }

/*           start applying rotations in 1st set from the left */

	    i__1 = *kb - k + 1;
	    for (l = *ka - 1; l >= *kb-k+1; --l) {
		nrt = (j2 + l - 1) / ka1;
		j1t = j2 - (nrt - 1) * ka1;
		if (nrt > 0) {
		    zlartv_(&nrt, &AB(ka1-l+1,j1t-ka1+l)
			    , &inca, &AB(ka1-l,j1t-ka1+l),
			     &inca, &RWORK(j1t), &WORK(j1t), &ka1);
		}
/* L820: */
	    }

	    if (wantx) {

/*              post-multiply X by product of rotations in 1st
 set */

		i__1 = j2;
		i__3 = ka1;
		for (j = j1; ka1 < 0 ? j >= j2 : j <= j2; j += ka1) {
		    d_cnjg(&z__1, &WORK(j));
		    zrot_(&nx, &X(1,j), &c__1, &X(1,j-1), &c__1, &RWORK(j), &z__1);
/* L830: */
		}
	    }
/* L840: */
	}

	if (update) {
	    if (i2 > 0 && kbt > 0) {

/*              create nonzero element a(i+kbt,i+kbt-ka-1) out
side the   
                band and store it in WORK(m-kb+i+kbt) */

		i__4 = m - *kb + i + kbt;
		i__3 = kbt + 1 + i * bb_dim1;
		z__2.r = -BB(kbt+1,i).r, z__2.i = -BB(kbt+1,i).i;
		z__1.r = z__2.r * ra1.r - z__2.i * ra1.i, z__1.i = z__2.r * 
			ra1.i + z__2.i * ra1.r;
		WORK(m-*kb+i+kbt).r = z__1.r, WORK(m-*kb+i+kbt).i = z__1.i;
	    }
	}

	for (k = *kb; k >= 1; --k) {
	    if (update) {
/* Computing MAX */
		i__4 = 2, i__3 = k + i0 - m;
		j2 = i + k + 1 - max(i__4,i__3) * ka1;
	    } else {
/* Computing MAX */
		i__4 = 1, i__3 = k + i0 - m;
		j2 = i + k + 1 - max(i__4,i__3) * ka1;
	    }

/*           finish applying rotations in 2nd set from the left */

	    for (l = *kb - k; l >= 1; --l) {
		nrt = (j2 + *ka + l - 1) / ka1;
		j1t = j2 - (nrt - 1) * ka1;
		if (nrt > 0) {
		    zlartv_(&nrt, &AB(ka1-l+1,j1t+l-1), 
			    &inca, &AB(ka1-l,j1t+l-1), &
			    inca, &RWORK(m - *kb + j1t + *ka), &WORK(m - *kb 
			    + j1t + *ka), &ka1);
		}
/* L850: */
	    }
	    nr = (j2 + *ka - 1) / ka1;
	    j1 = j2 - (nr - 1) * ka1;
	    i__4 = j2;
	    i__3 = ka1;
	    for (j = j1; ka1 < 0 ? j >= j2 : j <= j2; j += ka1) {
		i__1 = m - *kb + j;
		i__2 = m - *kb + j + *ka;
		WORK(m-*kb+j).r = WORK(m-*kb+j+*ka).r, WORK(m-*kb+j).i = WORK(m-*kb+j+*ka).i;
		RWORK(m - *kb + j) = RWORK(m - *kb + j + *ka);
/* L860: */
	    }
	    i__3 = j2;
	    i__4 = ka1;
	    for (j = j1; ka1 < 0 ? j >= j2 : j <= j2; j += ka1) {

/*              create nonzero element a(j+ka,j-1) outside the
 band   
                and store it in WORK(m-kb+j) */

		i__1 = m - *kb + j;
		i__2 = m - *kb + j;
		i__5 = ka1 + (j - 1) * ab_dim1;
		z__1.r = WORK(m-*kb+j).r * AB(ka1,j-1).r - WORK(m-*kb+j).i * AB(ka1,j-1)
			.i, z__1.i = WORK(m-*kb+j).r * AB(ka1,j-1).i + WORK(m-*kb+j).i 
			* AB(ka1,j-1).r;
		WORK(m-*kb+j).r = z__1.r, WORK(m-*kb+j).i = z__1.i;
		i__1 = ka1 + (j - 1) * ab_dim1;
		i__2 = m - *kb + j;
		i__5 = ka1 + (j - 1) * ab_dim1;
		z__1.r = RWORK(m-*kb+j) * AB(ka1,j-1).r, z__1.i = RWORK(m-*kb+j) * AB(ka1,j-1).i;
		AB(ka1,j-1).r = z__1.r, AB(ka1,j-1).i = z__1.i;
/* L870: */
	    }
	    if (update) {
		if (i + k > ka1 && k <= kbt) {
		    i__4 = m - *kb + i + k - *ka;
		    i__3 = m - *kb + i + k;
		    WORK(m-*kb+i+k-*ka).r = WORK(m-*kb+i+k).r, WORK(m-*kb+i+k-*ka).i = WORK(m-*kb+i+k).i;
		}
	    }
/* L880: */
	}

	for (k = *kb; k >= 1; --k) {
/* Computing MAX */
	    i__4 = 1, i__3 = k + i0 - m;
	    j2 = i + k + 1 - max(i__4,i__3) * ka1;
	    nr = (j2 + *ka - 1) / ka1;
	    j1 = j2 - (nr - 1) * ka1;
	    if (nr > 0) {

/*              generate rotations in 2nd set to annihilate el
ements   
                which have been created outside the band */

		zlargv_(&nr, &AB(ka1,j1), &inca, &WORK(m - *kb + 
			j1), &ka1, &RWORK(m - *kb + j1), &ka1);

/*              apply rotations in 2nd set from the right */

		i__4 = *ka - 1;
		for (l = 1; l <= *ka-1; ++l) {
		    zlartv_(&nr, &AB(l+1,j1), &inca, &AB(l+2,j1-1), &inca, &RWORK(m - *kb + j1)
			    , &WORK(m - *kb + j1), &ka1);
/* L890: */
		}

/*              apply rotations in 2nd set from both sides to 
diagonal   
                blocks */

		zlar2v_(&nr, &AB(1,j1), &AB(1,j1-1), &AB(2,j1-1), &inca, &RWORK(m - *
			kb + j1), &WORK(m - *kb + j1), &ka1);

		zlacgv_(&nr, &WORK(m - *kb + j1), &ka1);
	    }

/*           start applying rotations in 2nd set from the left */

	    i__4 = *kb - k + 1;
	    for (l = *ka - 1; l >= *kb-k+1; --l) {
		nrt = (j2 + l - 1) / ka1;
		j1t = j2 - (nrt - 1) * ka1;
		if (nrt > 0) {
		    zlartv_(&nrt, &AB(ka1-l+1,j1t-ka1+l)
			    , &inca, &AB(ka1-l,j1t-ka1+l),
			     &inca, &RWORK(m - *kb + j1t), &WORK(m - *kb + 
			    j1t), &ka1);
		}
/* L900: */
	    }

	    if (wantx) {

/*              post-multiply X by product of rotations in 2nd
 set */

		i__4 = j2;
		i__3 = ka1;
		for (j = j1; ka1 < 0 ? j >= j2 : j <= j2; j += ka1) {
		    d_cnjg(&z__1, &WORK(m - *kb + j));
		    zrot_(&nx, &X(1,j), &c__1, &X(1,j-1), &c__1, &RWORK(m - *kb + j), &z__1);
/* L910: */
		}
	    }
/* L920: */
	}

	i__3 = *kb - 1;
	for (k = 1; k <= *kb-1; ++k) {
/* Computing MAX */
	    i__4 = 1, i__1 = k + i0 - m + 1;
	    j2 = i + k + 1 - max(i__4,i__1) * ka1;

/*           finish applying rotations in 1st set from the left */

	    for (l = *kb - k; l >= 1; --l) {
		nrt = (j2 + l - 1) / ka1;
		j1t = j2 - (nrt - 1) * ka1;
		if (nrt > 0) {
		    zlartv_(&nrt, &AB(ka1-l+1,j1t-ka1+l)
			    , &inca, &AB(ka1-l,j1t-ka1+l),
			     &inca, &RWORK(j1t), &WORK(j1t), &ka1);
		}
/* L930: */
	    }
/* L940: */
	}

	i__3 = i2 - *ka;
	for (j = 2; j <= i2-*ka; ++j) {
	    RWORK(j) = RWORK(j + *ka);
	    i__4 = j;
	    i__1 = j + *ka;
	    WORK(j).r = WORK(j+*ka).r, WORK(j).i = WORK(j+*ka).i;
/* L950: */
	}

    }

    goto L490;

/*     End of ZHBGST */

} /* zhbgst_ */

