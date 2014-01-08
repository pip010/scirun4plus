#include "f2c.h"

/* Subroutine */ int zhbtrd_(char *vect, char *uplo, integer *n, integer *kd, 
	doublecomplex *ab, integer *ldab, doublereal *d, doublereal *e, 
	doublecomplex *q, integer *ldq, doublecomplex *work, integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    ZHBTRD reduces a complex Hermitian band matrix A to real symmetric   
    tridiagonal form T by a unitary similarity transformation:   
    Q**H * A * Q = T.   

    Arguments   
    =========   

    VECT    (input) CHARACTER*1   
            = 'N':  do not form Q;   
            = 'V':  form Q;   
            = 'U':  update a matrix X, by forming X*Q.   

    UPLO    (input) CHARACTER*1   
            = 'U':  Upper triangle of A is stored;   
            = 'L':  Lower triangle of A is stored.   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    KD      (input) INTEGER   
            The number of superdiagonals of the matrix A if UPLO = 'U',   
            or the number of subdiagonals if UPLO = 'L'.  KD >= 0.   

    AB      (input/output) COMPLEX*16 array, dimension (LDAB,N)   
            On entry, the upper or lower triangle of the Hermitian band   
            matrix A, stored in the first KD+1 rows of the array.  The   
            j-th column of A is stored in the j-th column of the array AB 
  
            as follows:   
            if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j; 
  
            if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd). 
  
            On exit, the diagonal elements of AB are overwritten by the   
            diagonal elements of the tridiagonal matrix T; if KD > 0, the 
  
            elements on the first superdiagonal (if UPLO = 'U') or the   
            first subdiagonal (if UPLO = 'L') are overwritten by the   
            offdiagonal elements of T; the rest of AB is overwritten by   
            values generated during the reduction.   

    LDAB    (input) INTEGER   
            The leading dimension of the array AB.  LDAB >= KD+1.   

    D       (output) DOUBLE PRECISION array, dimension (N)   
            The diagonal elements of the tridiagonal matrix T.   

    E       (output) DOUBLE PRECISION array, dimension (N-1)   
            The off-diagonal elements of the tridiagonal matrix T:   
            E(i) = T(i,i+1) if UPLO = 'U'; E(i) = T(i+1,i) if UPLO = 'L'. 
  

    Q       (input/output) COMPLEX*16 array, dimension (LDQ,N)   
            On entry, if VECT = 'U', then Q must contain an N-by-N   
            matrix X; if VECT = 'N' or 'V', then Q need not be set.   

            On exit:   
            if VECT = 'V', Q contains the N-by-N unitary matrix Q;   
            if VECT = 'U', Q contains the product X*Q;   
            if VECT = 'N', the array Q is not referenced.   

    LDQ     (input) INTEGER   
            The leading dimension of the array Q.   
            LDQ >= 1, and LDQ >= N if VECT = 'V' or 'U'.   

    WORK    (workspace) COMPLEX*16 array, dimension (N)   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   

    ===================================================================== 
  


       Test the input parameters   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static doublecomplex c_b1 = {0.,0.};
    static doublecomplex c_b2 = {1.,0.};
    static integer c__1 = 1;
    
    /* System generated locals */
    integer ab_dim1, ab_offset, q_dim1, q_offset, i__1, i__2, i__3, i__4, 
	    i__5, i__6;
    doublereal d__1;
    doublecomplex z__1;
    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);
    double z_abs(doublecomplex *);
    /* Local variables */
    static integer inca;
    static doublereal abst;
    static doublecomplex temp;
    extern /* Subroutine */ int zrot_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublereal *, doublecomplex *);
    static integer i, j, k, l;
    static doublecomplex t;
    extern logical lsame_(char *, char *);
    extern /* Subroutine */ int zscal_(integer *, doublecomplex *, 
	    doublecomplex *, integer *);
    static logical initq, wantq, upper;
    static integer j1, j2;
    extern /* Subroutine */ int zlar2v_(integer *, doublecomplex *, 
	    doublecomplex *, doublecomplex *, integer *, doublereal *, 
	    doublecomplex *, integer *);
    static integer nr;
    extern /* Subroutine */ int xerbla_(char *, integer *), zlacgv_(
	    integer *, doublecomplex *, integer *);
    static integer kd1;
    extern /* Subroutine */ int zlaset_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, doublecomplex *, integer *), zlartg_(doublecomplex *, doublecomplex *, doublereal *, 
	    doublecomplex *, doublecomplex *), zlargv_(integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublereal *, integer *), zlartv_(integer *, doublecomplex *, 
	    integer *, doublecomplex *, integer *, doublereal *, 
	    doublecomplex *, integer *);
    static integer kdn, nrt;



#define D(I) d[(I)-1]
#define E(I) e[(I)-1]
#define WORK(I) work[(I)-1]

#define AB(I,J) ab[(I)-1 + ((J)-1)* ( *ldab)]
#define Q(I,J) q[(I)-1 + ((J)-1)* ( *ldq)]

    initq = lsame_(vect, "V");
    wantq = initq || lsame_(vect, "U");
    upper = lsame_(uplo, "U");
    kd1 = *kd + 1;
    *info = 0;
    if (! wantq && ! lsame_(vect, "N")) {
	*info = -1;
    } else if (! upper && ! lsame_(uplo, "L")) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*kd < 0) {
	*info = -4;
    } else if (*ldab < kd1) {
	*info = -6;
    } else if (*ldq < max(1,*n) && wantq) {
	*info = -10;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("ZHBTRD", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return 0;
    }

/*     Initialize Q to the unit matrix, if needed */

    if (initq) {
	zlaset_("Full", n, n, &c_b1, &c_b2, &Q(1,1), ldq);
    }

/*     Wherever possible, plane rotations are generated and applied in   
       vector operations of length NR over the index set J1:J2:KD1.   

       The real cosines and complex sines of the plane rotations are   
       stored in the arrays D and WORK. */

    inca = kd1 * *ldab;
/* Computing MIN */
    i__1 = *n - 1;
    kdn = min(i__1,*kd);
    if (upper) {

	if (*kd > 1) {

/*           Reduce to complex Hermitian tridiagonal form, working
 with   
             the upper triangle */

	    nr = 0;
	    j1 = kdn + 2;
	    j2 = 1;

	    i__1 = kd1 + ab_dim1;
	    i__2 = kd1 + ab_dim1;
	    d__1 = AB(kd1,1).r;
	    AB(kd1,1).r = d__1, AB(kd1,1).i = 0.;
	    i__1 = *n - 2;
	    for (i = 1; i <= *n-2; ++i) {

/*              Reduce i-th row of matrix to tridiagonal form 
*/

		for (k = kdn + 1; k >= 2; --k) {
		    j1 += kdn;
		    j2 += kdn;

		    if (nr > 0) {

/*                    generate plane rotations to anni
hilate nonzero   
                      elements which have been created
 outside the band */

			zlargv_(&nr, &AB(1,j1-1), &inca, &
				WORK(j1), &kd1, &D(j1), &kd1);

/*                    apply rotations from the right 
*/

			i__2 = *kd - 1;
			for (l = 1; l <= *kd-1; ++l) {
			    zlartv_(&nr, &AB(l+1,j1-1), &
				    inca, &AB(l,j1), &inca, &D(j1)
				    , &WORK(j1), &kd1);
/* L10: */
			}
		    }

		    if (k > 2) {
			if (k <= *n - i + 1) {

/*                       generate plane rotation t
o annihilate a(i,i+k-1)   
                         within the band */

			    zlartg_(&AB(*kd-k+3,i+k-2), 
				    &AB(*kd-k+2,i+k-1), 
				    &D(i + k - 1), &WORK(i + k - 1), &temp);
			    i__2 = *kd - k + 3 + (i + k - 2) * ab_dim1;
			    AB(*kd-k+3,i+k-2).r = temp.r, AB(*kd-k+3,i+k-2).i = temp.i;

/*                       apply rotation from the r
ight */

			    i__2 = k - 3;
			    zrot_(&i__2, &AB(*kd-k+4,i+k-2), &c__1, &AB(*kd-k+3,i+k-1), &c__1, &D(i + k - 1), &
				    WORK(i + k - 1));
			}
			++nr;
			j1 = j1 - kdn - 1;
		    }

/*                 apply plane rotations from both sides t
o diagonal   
                   blocks */

		    if (nr > 0) {
			zlar2v_(&nr, &AB(kd1,j1-1), &AB(kd1,j1), &AB(*kd,j1), &inca,
				 &D(j1), &WORK(j1), &kd1);
		    }

/*                 apply plane rotations from the left */

		    zlacgv_(&nr, &WORK(j1), &kd1);
		    i__2 = *kd - 1;
		    for (l = 1; l <= *kd-1; ++l) {
			if (j2 + l > *n) {
			    nrt = nr - 1;
			} else {
			    nrt = nr;
			}
			if (nrt > 0) {
			    zlartv_(&nrt, &AB(*kd-l,j1+l), &
				    inca, &AB(*kd-l+1,j1+l), &inca, &D(j1), &WORK(j1), &kd1);
			}
/* L20: */
		    }

		    if (wantq) {

/*                    accumulate product of plane rota
tions in Q */

			i__2 = j2;
			i__3 = kd1;
			for (j = j1; i__3 < 0 ? j >= i__2 : j <= i__2; j += 
				i__3) {
			    d_cnjg(&z__1, &WORK(j));
			    zrot_(n, &Q(1,j-1), &c__1, &Q(1,j), &c__1, &D(j), &z__1);
/* L30: */
			}
		    }

		    if (j2 + kdn > *n) {

/*                    adjust J2 to keep within the bou
nds of the matrix */

			--nr;
			j2 = j2 - kdn - 1;
		    }

		    i__3 = j2;
		    i__2 = kd1;
		    for (j = j1; kd1 < 0 ? j >= j2 : j <= j2; j += kd1) 
			    {

/*                    create nonzero element a(j-1,j+k
d) outside the band   
                      and store it in WORK */

			i__4 = j + *kd;
			i__5 = j;
			i__6 = (j + *kd) * ab_dim1 + 1;
			z__1.r = WORK(j).r * AB(1,j+*kd).r - WORK(j).i * 
				AB(1,j+*kd).i, z__1.i = WORK(j).r * AB(1,j+*kd)
				.i + WORK(j).i * AB(1,j+*kd).r;
			WORK(j+*kd).r = z__1.r, WORK(j+*kd).i = z__1.i;
			i__4 = (j + *kd) * ab_dim1 + 1;
			i__5 = j;
			i__6 = (j + *kd) * ab_dim1 + 1;
			z__1.r = D(j) * AB(1,j+*kd).r, z__1.i = D(j) * AB(1,j+*kd).i;
			AB(1,j+*kd).r = z__1.r, AB(1,j+*kd).i = z__1.i;
/* L40: */
		    }
/* L50: */
		}
/* L60: */
	    }
	}

	if (*kd > 0) {

/*           make off-diagonal elements real and copy them to E */

	    i__1 = *n - 1;
	    for (i = 1; i <= *n-1; ++i) {
		i__2 = *kd + (i + 1) * ab_dim1;
		t.r = AB(*kd,i+1).r, t.i = AB(*kd,i+1).i;
		abst = z_abs(&t);
		i__2 = *kd + (i + 1) * ab_dim1;
		AB(*kd,i+1).r = abst, AB(*kd,i+1).i = 0.;
		E(i) = abst;
		if (abst != 0.) {
		    z__1.r = t.r / abst, z__1.i = t.i / abst;
		    t.r = z__1.r, t.i = z__1.i;
		} else {
		    t.r = 1., t.i = 0.;
		}
		if (i < *n - 1) {
		    i__2 = *kd + (i + 2) * ab_dim1;
		    i__3 = *kd + (i + 2) * ab_dim1;
		    z__1.r = AB(*kd,i+2).r * t.r - AB(*kd,i+2).i * t.i, z__1.i = AB(*kd,i+2).r * t.i + AB(*kd,i+2).i * t.r;
		    AB(*kd,i+2).r = z__1.r, AB(*kd,i+2).i = z__1.i;
		}
		if (wantq) {
		    d_cnjg(&z__1, &t);
		    zscal_(n, &z__1, &Q(1,i+1), &c__1);
		}
/* L70: */
	    }
	} else {

/*           set E to zero if original matrix was diagonal */

	    i__1 = *n - 1;
	    for (i = 1; i <= *n-1; ++i) {
		E(i) = 0.;
/* L80: */
	    }
	}

/*        copy diagonal elements to D */

	i__1 = *n;
	for (i = 1; i <= *n; ++i) {
	    i__2 = i;
	    i__3 = kd1 + i * ab_dim1;
	    D(i) = AB(kd1,i).r;
/* L90: */
	}

    } else {

	if (*kd > 1) {

/*           Reduce to complex Hermitian tridiagonal form, working
 with   
             the lower triangle */

	    nr = 0;
	    j1 = kdn + 2;
	    j2 = 1;

	    i__1 = ab_dim1 + 1;
	    i__2 = ab_dim1 + 1;
	    d__1 = AB(1,1).r;
	    AB(1,1).r = d__1, AB(1,1).i = 0.;
	    i__1 = *n - 2;
	    for (i = 1; i <= *n-2; ++i) {

/*              Reduce i-th column of matrix to tridiagonal fo
rm */

		for (k = kdn + 1; k >= 2; --k) {
		    j1 += kdn;
		    j2 += kdn;

		    if (nr > 0) {

/*                    generate plane rotations to anni
hilate nonzero   
                      elements which have been created
 outside the band */

			zlargv_(&nr, &AB(kd1,j1-kd1), &inca, &
				WORK(j1), &kd1, &D(j1), &kd1);

/*                    apply plane rotations from one s
ide */

			i__2 = *kd - 1;
			for (l = 1; l <= *kd-1; ++l) {
			    zlartv_(&nr, &AB(kd1-l,j1-kd1+l), &inca, &AB(kd1-l+1,j1-kd1+l), &inca, &D(j1), &WORK(
				    j1), &kd1);
/* L100: */
			}
		    }

		    if (k > 2) {
			if (k <= *n - i + 1) {

/*                       generate plane rotation t
o annihilate a(i+k-1,i)   
                         within the band */

			    zlartg_(&AB(k-1,i), &AB(k,i), &D(i + k - 1), &WORK(i + k - 1),
				     &temp);
			    i__2 = k - 1 + i * ab_dim1;
			    AB(k-1,i).r = temp.r, AB(k-1,i).i = temp.i;

/*                       apply rotation from the l
eft */

			    i__2 = k - 3;
			    i__3 = *ldab - 1;
			    i__4 = *ldab - 1;
			    zrot_(&i__2, &AB(k-2,i+1), &
				    i__3, &AB(k-1,i+1), &
				    i__4, &D(i + k - 1), &WORK(i + k - 1));
			}
			++nr;
			j1 = j1 - kdn - 1;
		    }

/*                 apply plane rotations from both sides t
o diagonal   
                   blocks */

		    if (nr > 0) {
			zlar2v_(&nr, &AB(1,j1-1), &AB(1,j1), &AB(2,j1-1), &
				inca, &D(j1), &WORK(j1), &kd1);
		    }

/*                 apply plane rotations from the right */

		    zlacgv_(&nr, &WORK(j1), &kd1);
		    i__2 = *kd - 1;
		    for (l = 1; l <= *kd-1; ++l) {
			if (j2 + l > *n) {
			    nrt = nr - 1;
			} else {
			    nrt = nr;
			}
			if (nrt > 0) {
			    zlartv_(&nrt, &AB(l+2,j1-1), &
				    inca, &AB(l+1,j1), &inca, &
				    D(j1), &WORK(j1), &kd1);
			}
/* L110: */
		    }

		    if (wantq) {

/*                    accumulate product of plane rota
tions in Q */

			i__2 = j2;
			i__3 = kd1;
			for (j = j1; i__3 < 0 ? j >= i__2 : j <= i__2; j += 
				i__3) {
			    zrot_(n, &Q(1,j-1), &c__1, &Q(1,j), &c__1, &D(j), &WORK(j));
/* L120: */
			}
		    }

		    if (j2 + kdn > *n) {

/*                    adjust J2 to keep within the bou
nds of the matrix */

			--nr;
			j2 = j2 - kdn - 1;
		    }

		    i__3 = j2;
		    i__2 = kd1;
		    for (j = j1; kd1 < 0 ? j >= j2 : j <= j2; j += kd1) 
			    {

/*                    create nonzero element a(j+kd,j-
1) outside the   
                      band and store it in WORK */

			i__4 = j + *kd;
			i__5 = j;
			i__6 = kd1 + j * ab_dim1;
			z__1.r = WORK(j).r * AB(kd1,j).r - WORK(j).i * 
				AB(kd1,j).i, z__1.i = WORK(j).r * AB(kd1,j)
				.i + WORK(j).i * AB(kd1,j).r;
			WORK(j+*kd).r = z__1.r, WORK(j+*kd).i = z__1.i;
			i__4 = kd1 + j * ab_dim1;
			i__5 = j;
			i__6 = kd1 + j * ab_dim1;
			z__1.r = D(j) * AB(kd1,j).r, z__1.i = D(j) * AB(kd1,j).i;
			AB(kd1,j).r = z__1.r, AB(kd1,j).i = z__1.i;
/* L130: */
		    }
/* L140: */
		}
/* L150: */
	    }
	}

	if (*kd > 0) {

/*           make off-diagonal elements real and copy them to E */

	    i__1 = *n - 1;
	    for (i = 1; i <= *n-1; ++i) {
		i__2 = i * ab_dim1 + 2;
		t.r = AB(2,i).r, t.i = AB(2,i).i;
		abst = z_abs(&t);
		i__2 = i * ab_dim1 + 2;
		AB(2,i).r = abst, AB(2,i).i = 0.;
		E(i) = abst;
		if (abst != 0.) {
		    z__1.r = t.r / abst, z__1.i = t.i / abst;
		    t.r = z__1.r, t.i = z__1.i;
		} else {
		    t.r = 1., t.i = 0.;
		}
		if (i < *n - 1) {
		    i__2 = (i + 1) * ab_dim1 + 2;
		    i__3 = (i + 1) * ab_dim1 + 2;
		    z__1.r = AB(2,i+1).r * t.r - AB(2,i+1).i * t.i, z__1.i = AB(2,i+1).r * t.i + AB(2,i+1).i * t.r;
		    AB(2,i+1).r = z__1.r, AB(2,i+1).i = z__1.i;
		}
		if (wantq) {
		    zscal_(n, &t, &Q(1,i+1), &c__1);
		}
/* L160: */
	    }
	} else {

/*           set E to zero if original matrix was diagonal */

	    i__1 = *n - 1;
	    for (i = 1; i <= *n-1; ++i) {
		E(i) = 0.;
/* L170: */
	    }
	}

/*        copy diagonal elements to D */

	i__1 = *n;
	for (i = 1; i <= *n; ++i) {
	    i__2 = i;
	    i__3 = i * ab_dim1 + 1;
	    D(i) = AB(1,i).r;
/* L180: */
	}
    }

    return 0;

/*     End of ZHBTRD */

} /* zhbtrd_ */

