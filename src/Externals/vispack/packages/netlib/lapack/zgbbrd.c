#include "f2c.h"

/* Subroutine */ int zgbbrd_(char *vect, integer *m, integer *n, integer *ncc,
	 integer *kl, integer *ku, doublecomplex *ab, integer *ldab, 
	doublereal *d, doublereal *e, doublecomplex *q, integer *ldq, 
	doublecomplex *pt, integer *ldpt, doublecomplex *c, integer *ldc, 
	doublecomplex *work, doublereal *rwork, integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    ZGBBRD reduces a complex general m-by-n band matrix A to real upper   
    bidiagonal form B by a unitary transformation: Q' * A * P = B.   

    The routine computes B, and optionally forms Q or P', or computes   
    Q'*C for a given matrix C.   

    Arguments   
    =========   

    VECT    (input) CHARACTER*1   
            Specifies whether or not the matrices Q and P' are to be   
            formed.   
            = 'N': do not form Q or P';   
            = 'Q': form Q only;   
            = 'P': form P' only;   
            = 'B': form both.   

    M       (input) INTEGER   
            The number of rows of the matrix A.  M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrix A.  N >= 0.   

    NCC     (input) INTEGER   
            The number of columns of the matrix C.  NCC >= 0.   

    KL      (input) INTEGER   
            The number of subdiagonals of the matrix A. KL >= 0.   

    KU      (input) INTEGER   
            The number of superdiagonals of the matrix A. KU >= 0.   

    AB      (input/output) COMPLEX*16 array, dimension (LDAB,N)   
            On entry, the m-by-n band matrix A, stored in rows 1 to   
            KL+KU+1. The j-th column of A is stored in the j-th column of 
  
            the array AB as follows:   
            AB(ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl).   
            On exit, A is overwritten by values generated during the   
            reduction.   

    LDAB    (input) INTEGER   
            The leading dimension of the array A. LDAB >= KL+KU+1.   

    D       (output) DOUBLE PRECISION array, dimension (min(M,N))   
            The diagonal elements of the bidiagonal matrix B.   

    E       (output) DOUBLE PRECISION array, dimension (min(M,N)-1)   
            The superdiagonal elements of the bidiagonal matrix B.   

    Q       (output) COMPLEX*16 array, dimension (LDQ,M)   
            If VECT = 'Q' or 'B', the m-by-m unitary matrix Q.   
            If VECT = 'N' or 'P', the array Q is not referenced.   

    LDQ     (input) INTEGER   
            The leading dimension of the array Q.   
            LDQ >= max(1,M) if VECT = 'Q' or 'B'; LDQ >= 1 otherwise.   

    PT      (output) COMPLEX*16 array, dimension (LDPT,N)   
            If VECT = 'P' or 'B', the n-by-n unitary matrix P'.   
            If VECT = 'N' or 'Q', the array PT is not referenced.   

    LDPT    (input) INTEGER   
            The leading dimension of the array PT.   
            LDPT >= max(1,N) if VECT = 'P' or 'B'; LDPT >= 1 otherwise.   

    C       (input/output) COMPLEX*16 array, dimension (LDC,NCC)   
            On entry, an m-by-ncc matrix C.   
            On exit, C is overwritten by Q'*C.   
            C is not referenced if NCC = 0.   

    LDC     (input) INTEGER   
            The leading dimension of the array C.   
            LDC >= max(1,M) if NCC > 0; LDC >= 1 if NCC = 0.   

    WORK    (workspace) COMPLEX*16 array, dimension (max(M,N))   

    RWORK   (workspace) DOUBLE PRECISION array, dimension (max(M,N))   

    INFO    (output) INTEGER   
            = 0:  successful exit.   
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
    integer ab_dim1, ab_offset, c_dim1, c_offset, pt_dim1, pt_offset, q_dim1, 
	    q_offset, i__1, i__2, i__3, i__4, i__5, i__6, i__7;
    doublecomplex z__1, z__2, z__3;
    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);
    double z_abs(doublecomplex *);
    /* Local variables */
    static integer inca;
    static doublereal abst;
    extern /* Subroutine */ int zrot_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublereal *, doublecomplex *);
    static integer i, j, l;
    static doublecomplex t;
    extern logical lsame_(char *, char *);
    static logical wantb, wantc;
    extern /* Subroutine */ int zscal_(integer *, doublecomplex *, 
	    doublecomplex *, integer *);
    static integer minmn;
    static logical wantq;
    static integer j1, j2, kb;
    static doublecomplex ra, rb;
    static doublereal rc;
    static integer kk, ml, nr, mu;
    static doublecomplex rs;
    extern /* Subroutine */ int xerbla_(char *, integer *);
    static integer kb1;
    extern /* Subroutine */ int zlaset_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, doublecomplex *, integer *), zlartg_(doublecomplex *, doublecomplex *, doublereal *, 
	    doublecomplex *, doublecomplex *), zlargv_(integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublereal *, integer *);
    static integer ml0;
    static logical wantpt;
    static integer mu0;
    extern /* Subroutine */ int zlartv_(integer *, doublecomplex *, integer *,
	     doublecomplex *, integer *, doublereal *, doublecomplex *, 
	    integer *);
    static integer klm, kun, nrt, klu1;



#define D(I) d[(I)-1]
#define E(I) e[(I)-1]
#define WORK(I) work[(I)-1]
#define RWORK(I) rwork[(I)-1]

#define AB(I,J) ab[(I)-1 + ((J)-1)* ( *ldab)]
#define Q(I,J) q[(I)-1 + ((J)-1)* ( *ldq)]
#define PT(I,J) pt[(I)-1 + ((J)-1)* ( *ldpt)]
#define C(I,J) c[(I)-1 + ((J)-1)* ( *ldc)]

    wantb = lsame_(vect, "B");
    wantq = lsame_(vect, "Q") || wantb;
    wantpt = lsame_(vect, "P") || wantb;
    wantc = *ncc > 0;
    klu1 = *kl + *ku + 1;
    *info = 0;
    if (! wantq && ! wantpt && ! lsame_(vect, "N")) {
	*info = -1;
    } else if (*m < 0) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*ncc < 0) {
	*info = -4;
    } else if (*kl < 0) {
	*info = -5;
    } else if (*ku < 0) {
	*info = -6;
    } else if (*ldab < klu1) {
	*info = -8;
    } else if (*ldq < 1 || wantq && *ldq < max(1,*m)) {
	*info = -12;
    } else if (*ldpt < 1 || wantpt && *ldpt < max(1,*n)) {
	*info = -14;
    } else if (*ldc < 1 || wantc && *ldc < max(1,*m)) {
	*info = -16;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("ZGBBRD", &i__1);
	return 0;
    }

/*     Initialize Q and P' to the unit matrix, if needed */

    if (wantq) {
	zlaset_("Full", m, m, &c_b1, &c_b2, &Q(1,1), ldq);
    }
    if (wantpt) {
	zlaset_("Full", n, n, &c_b1, &c_b2, &PT(1,1), ldpt);
    }

/*     Quick return if possible. */

    if (*m == 0 || *n == 0) {
	return 0;
    }

    minmn = min(*m,*n);

    if (*kl + *ku > 1) {

/*        Reduce to upper bidiagonal form if KU > 0; if KU = 0, reduce
   
          first to lower bidiagonal form and then transform to upper 
  
          bidiagonal */

	if (*ku > 0) {
	    ml0 = 1;
	    mu0 = 2;
	} else {
	    ml0 = 2;
	    mu0 = 1;
	}

/*        Wherever possible, plane rotations are generated and applied
 in   
          vector operations of length NR over the index set J1:J2:KLU1
.   

          The complex sines of the plane rotations are stored in WORK,
   
          and the real cosines in RWORK.   

   Computing MIN */
	i__1 = *m - 1;
	klm = min(i__1,*kl);
/* Computing MIN */
	i__1 = *n - 1;
	kun = min(i__1,*ku);
	kb = klm + kun;
	kb1 = kb + 1;
	inca = kb1 * *ldab;
	nr = 0;
	j1 = klm + 2;
	j2 = 1 - kun;

	i__1 = minmn;
	for (i = 1; i <= minmn; ++i) {

/*           Reduce i-th column and i-th row of matrix to bidiagon
al form */

	    ml = klm + 1;
	    mu = kun + 1;
	    i__2 = kb;
	    for (kk = 1; kk <= kb; ++kk) {
		j1 += kb;
		j2 += kb;

/*              generate plane rotations to annihilate nonzero
 elements   
                which have been created below the band */

		if (nr > 0) {
		    zlargv_(&nr, &AB(klu1,j1-klm-1), &inca, 
			    &WORK(j1), &kb1, &RWORK(j1), &kb1);
		}

/*              apply plane rotations from the left */

		i__3 = kb;
		for (l = 1; l <= kb; ++l) {
		    if (j2 - klm + l - 1 > *n) {
			nrt = nr - 1;
		    } else {
			nrt = nr;
		    }
		    if (nrt > 0) {
			zlartv_(&nrt, &AB(klu1-l,j1-klm+l-1), &inca, &AB(klu1-l+1,j1-klm+l-1), &inca, &RWORK(j1), &WORK(
				j1), &kb1);
		    }
/* L10: */
		}

		if (ml > ml0) {
		    if (ml <= *m - i + 1) {

/*                    generate plane rotation to annih
ilate a(i+ml-1,i)   
                      within the band, and apply rotat
ion from the left */

			zlartg_(&AB(*ku+ml-1,i), &AB(*ku+ml,i), &RWORK(i + ml - 1), &WORK(i + 
				ml - 1), &ra);
			i__3 = *ku + ml - 1 + i * ab_dim1;
			AB(*ku+ml-1,i).r = ra.r, AB(*ku+ml-1,i).i = ra.i;
			if (i < *n) {
/* Computing MIN */
			    i__4 = *ku + ml - 2, i__5 = *n - i;
			    i__3 = min(i__4,i__5);
			    i__6 = *ldab - 1;
			    i__7 = *ldab - 1;
			    zrot_(&i__3, &AB(*ku+ml-2,i+1)
				    , &i__6, &AB(*ku+ml-1,i+1), &i__7, &RWORK(i + ml - 1), &
				    WORK(i + ml - 1));
			}
		    }
		    ++nr;
		    j1 -= kb1;
		}

		if (wantq) {

/*                 accumulate product of plane rotations i
n Q */

		    i__3 = j2;
		    i__4 = kb1;
		    for (j = j1; kb1 < 0 ? j >= j2 : j <= j2; j += kb1) 
			    {
			d_cnjg(&z__1, &WORK(j));
			zrot_(m, &Q(1,j-1), &c__1, &Q(1,j), &c__1, &RWORK(j), &z__1);
/* L20: */
		    }
		}

		if (wantc) {

/*                 apply plane rotations to C */

		    i__4 = j2;
		    i__3 = kb1;
		    for (j = j1; kb1 < 0 ? j >= j2 : j <= j2; j += kb1) 
			    {
			zrot_(ncc, &C(j-1,1), ldc, &C(j,1), 
				ldc, &RWORK(j), &WORK(j));
/* L30: */
		    }
		}

		if (j2 + kun > *n) {

/*                 adjust J2 to keep within the bounds of 
the matrix */

		    --nr;
		    j2 -= kb1;
		}

		i__3 = j2;
		i__4 = kb1;
		for (j = j1; kb1 < 0 ? j >= j2 : j <= j2; j += kb1) {

/*                 create nonzero element a(j-1,j+ku) abov
e the band   
                   and store it in WORK(n+1:2*n) */

		    i__5 = j + kun;
		    i__6 = j;
		    i__7 = (j + kun) * ab_dim1 + 1;
		    z__1.r = WORK(j).r * AB(1,j+kun).r - WORK(j).i * AB(1,j+kun).i, z__1.i = WORK(j).r * AB(1,j+kun).i + 
			    WORK(j).i * AB(1,j+kun).r;
		    WORK(j+kun).r = z__1.r, WORK(j+kun).i = z__1.i;
		    i__5 = (j + kun) * ab_dim1 + 1;
		    i__6 = j;
		    i__7 = (j + kun) * ab_dim1 + 1;
		    z__1.r = RWORK(j) * AB(1,j+kun).r, z__1.i = RWORK(j) * 
			    AB(1,j+kun).i;
		    AB(1,j+kun).r = z__1.r, AB(1,j+kun).i = z__1.i;
/* L40: */
		}

/*              generate plane rotations to annihilate nonzero
 elements   
                which have been generated above the band */

		if (nr > 0) {
		    zlargv_(&nr, &AB(1,j1+kun-1), &inca, &
			    WORK(j1 + kun), &kb1, &RWORK(j1 + kun), &kb1);
		}

/*              apply plane rotations from the right */

		i__4 = kb;
		for (l = 1; l <= kb; ++l) {
		    if (j2 + l - 1 > *m) {
			nrt = nr - 1;
		    } else {
			nrt = nr;
		    }
		    if (nrt > 0) {
			zlartv_(&nrt, &AB(l+1,j1+kun-1), &
				inca, &AB(l,j1+kun), &inca, &
				RWORK(j1 + kun), &WORK(j1 + kun), &kb1);
		    }
/* L50: */
		}

		if (ml == ml0 && mu > mu0) {
		    if (mu <= *n - i + 1) {

/*                    generate plane rotation to annih
ilate a(i,i+mu-1)   
                      within the band, and apply rotat
ion from the right */

			zlartg_(&AB(*ku-mu+3,i+mu-2), &
				AB(*ku-mu+2,i+mu-1), &
				RWORK(i + mu - 1), &WORK(i + mu - 1), &ra);
			i__4 = *ku - mu + 3 + (i + mu - 2) * ab_dim1;
			AB(*ku-mu+3,i+mu-2).r = ra.r, AB(*ku-mu+3,i+mu-2).i = ra.i;
/* Computing MIN */
			i__3 = *kl + mu - 2, i__5 = *m - i;
			i__4 = min(i__3,i__5);
			zrot_(&i__4, &AB(*ku-mu+4,i+mu-2), &c__1, &AB(*ku-mu+3,i+mu-1), &c__1, &RWORK(i + mu - 1), &
				WORK(i + mu - 1));
		    }
		    ++nr;
		    j1 -= kb1;
		}

		if (wantpt) {

/*                 accumulate product of plane rotations i
n P' */

		    i__4 = j2;
		    i__3 = kb1;
		    for (j = j1; kb1 < 0 ? j >= j2 : j <= j2; j += kb1) 
			    {
			d_cnjg(&z__1, &WORK(j + kun));
			zrot_(n, &PT(j+kun-1,1), ldpt, &PT(j+kun,1), ldpt, &RWORK(j + kun), &z__1);
/* L60: */
		    }
		}

		if (j2 + kb > *m) {

/*                 adjust J2 to keep within the bounds of 
the matrix */

		    --nr;
		    j2 -= kb1;
		}

		i__3 = j2;
		i__4 = kb1;
		for (j = j1; kb1 < 0 ? j >= j2 : j <= j2; j += kb1) {

/*                 create nonzero element a(j+kl+ku,j+ku-1
) below the   
                   band and store it in WORK(1:n) */

		    i__5 = j + kb;
		    i__6 = j + kun;
		    i__7 = klu1 + (j + kun) * ab_dim1;
		    z__1.r = WORK(j+kun).r * AB(klu1,j+kun).r - WORK(j+kun).i * AB(klu1,j+kun).i, z__1.i = WORK(j+kun).r * AB(klu1,j+kun).i + 
			    WORK(j+kun).i * AB(klu1,j+kun).r;
		    WORK(j+kb).r = z__1.r, WORK(j+kb).i = z__1.i;
		    i__5 = klu1 + (j + kun) * ab_dim1;
		    i__6 = j + kun;
		    i__7 = klu1 + (j + kun) * ab_dim1;
		    z__1.r = RWORK(j+kun) * AB(klu1,j+kun).r, z__1.i = RWORK(j+kun) * 
			    AB(klu1,j+kun).i;
		    AB(klu1,j+kun).r = z__1.r, AB(klu1,j+kun).i = z__1.i;
/* L70: */
		}

		if (ml > ml0) {
		    --ml;
		} else {
		    --mu;
		}
/* L80: */
	    }
/* L90: */
	}
    }

    if (*ku == 0 && *kl > 0) {

/*        A has been reduced to complex lower bidiagonal form   

          Transform lower bidiagonal form to upper bidiagonal by apply
ing   
          plane rotations from the left, overwriting superdiagonal   
          elements on subdiagonal elements   

   Computing MIN */
	i__2 = *m - 1;
	i__1 = min(i__2,*n);
	for (i = 1; i <= min(*m-1,*n); ++i) {
	    zlartg_(&AB(1,i), &AB(2,i), &rc, &rs, &ra)
		    ;
	    i__2 = i * ab_dim1 + 1;
	    AB(1,i).r = ra.r, AB(1,i).i = ra.i;
	    if (i < *n) {
		i__2 = i * ab_dim1 + 2;
		i__4 = (i + 1) * ab_dim1 + 1;
		z__1.r = rs.r * AB(1,i+1).r - rs.i * AB(1,i+1).i, z__1.i = rs.r 
			* AB(1,i+1).i + rs.i * AB(1,i+1).r;
		AB(2,i).r = z__1.r, AB(2,i).i = z__1.i;
		i__2 = (i + 1) * ab_dim1 + 1;
		i__4 = (i + 1) * ab_dim1 + 1;
		z__1.r = rc * AB(1,i+1).r, z__1.i = rc * AB(1,i+1).i;
		AB(1,i+1).r = z__1.r, AB(1,i+1).i = z__1.i;
	    }
	    if (wantq) {
		d_cnjg(&z__1, &rs);
		zrot_(m, &Q(1,i), &c__1, &Q(1,i+1), 
			&c__1, &rc, &z__1);
	    }
	    if (wantc) {
		zrot_(ncc, &C(i,1), ldc, &C(i+1,1), ldc, &rc, 
			&rs);
	    }
/* L100: */
	}
    } else {

/*        A has been reduced to complex upper bidiagonal form or is   
          diagonal */

	if (*ku > 0 && *m < *n) {

/*           Annihilate a(m,m+1) by applying plane rotations from 
the   
             right */

	    i__1 = *ku + (*m + 1) * ab_dim1;
	    rb.r = AB(*ku,*m+1).r, rb.i = AB(*ku,*m+1).i;
	    for (i = *m; i >= 1; --i) {
		zlartg_(&AB(*ku+1,i), &rb, &rc, &rs, &ra);
		i__1 = *ku + 1 + i * ab_dim1;
		AB(*ku+1,i).r = ra.r, AB(*ku+1,i).i = ra.i;
		if (i > 1) {
		    d_cnjg(&z__3, &rs);
		    z__2.r = -z__3.r, z__2.i = -z__3.i;
		    i__1 = *ku + i * ab_dim1;
		    z__1.r = z__2.r * AB(*ku,i).r - z__2.i * AB(*ku,i).i, 
			    z__1.i = z__2.r * AB(*ku,i).i + z__2.i * AB(*ku,i)
			    .r;
		    rb.r = z__1.r, rb.i = z__1.i;
		    i__1 = *ku + i * ab_dim1;
		    i__2 = *ku + i * ab_dim1;
		    z__1.r = rc * AB(*ku,i).r, z__1.i = rc * AB(*ku,i).i;
		    AB(*ku,i).r = z__1.r, AB(*ku,i).i = z__1.i;
		}
		if (wantpt) {
		    d_cnjg(&z__1, &rs);
		    zrot_(n, &PT(i,1), ldpt, &PT(*m+1,1), 
			    ldpt, &rc, &z__1);
		}
/* L110: */
	    }
	}
    }

/*     Make diagonal and superdiagonal elements real, storing them in D   
       and E */

    i__1 = *ku + 1 + ab_dim1;
    t.r = AB(*ku+1,1).r, t.i = AB(*ku+1,1).i;
    i__1 = minmn;
    for (i = 1; i <= minmn; ++i) {
	abst = z_abs(&t);
	D(i) = abst;
	if (abst != 0.) {
	    z__1.r = t.r / abst, z__1.i = t.i / abst;
	    t.r = z__1.r, t.i = z__1.i;
	} else {
	    t.r = 1., t.i = 0.;
	}
	if (wantq) {
	    zscal_(m, &t, &Q(1,i), &c__1);
	}
	if (wantc) {
	    d_cnjg(&z__1, &t);
	    zscal_(ncc, &z__1, &C(i,1), ldc);
	}
	if (i < minmn) {
	    if (*ku == 0 && *kl == 0) {
		E(i) = 0.;
		i__2 = (i + 1) * ab_dim1 + 1;
		t.r = AB(1,i+1).r, t.i = AB(1,i+1).i;
	    } else {
		if (*ku == 0) {
		    i__2 = i * ab_dim1 + 2;
		    d_cnjg(&z__2, &t);
		    z__1.r = AB(2,i).r * z__2.r - AB(2,i).i * z__2.i, 
			    z__1.i = AB(2,i).r * z__2.i + AB(2,i).i * 
			    z__2.r;
		    t.r = z__1.r, t.i = z__1.i;
		} else {
		    i__2 = *ku + (i + 1) * ab_dim1;
		    d_cnjg(&z__2, &t);
		    z__1.r = AB(*ku,i+1).r * z__2.r - AB(*ku,i+1).i * z__2.i, 
			    z__1.i = AB(*ku,i+1).r * z__2.i + AB(*ku,i+1).i * 
			    z__2.r;
		    t.r = z__1.r, t.i = z__1.i;
		}
		abst = z_abs(&t);
		E(i) = abst;
		if (abst != 0.) {
		    z__1.r = t.r / abst, z__1.i = t.i / abst;
		    t.r = z__1.r, t.i = z__1.i;
		} else {
		    t.r = 1., t.i = 0.;
		}
		if (wantpt) {
		    zscal_(n, &t, &PT(i+1,1), ldpt);
		}
		i__2 = *ku + 1 + (i + 1) * ab_dim1;
		d_cnjg(&z__2, &t);
		z__1.r = AB(*ku+1,i+1).r * z__2.r - AB(*ku+1,i+1).i * z__2.i, z__1.i = 
			AB(*ku+1,i+1).r * z__2.i + AB(*ku+1,i+1).i * z__2.r;
		t.r = z__1.r, t.i = z__1.i;
	    }
	}
/* L120: */
    }
    return 0;

/*     End of ZGBBRD */

} /* zgbbrd_ */

