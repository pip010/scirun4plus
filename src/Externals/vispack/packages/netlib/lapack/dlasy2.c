#include "f2c.h"

/* Subroutine */ int dlasy2_(logical *ltranl, logical *ltranr, integer *isgn, 
	integer *n1, integer *n2, doublereal *tl, integer *ldtl, doublereal *
	tr, integer *ldtr, doublereal *b, integer *ldb, doublereal *scale, 
	doublereal *x, integer *ldx, doublereal *xnorm, integer *info)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    DLASY2 solves for the N1 by N2 matrix X, 1 <= N1,N2 <= 2, in   

           op(TL)*X + ISGN*X*op(TR) = SCALE*B,   

    where TL is N1 by N1, TR is N2 by N2, B is N1 by N2, and ISGN = 1 or 
  
    -1.  op(T) = T or T', where T' denotes the transpose of T.   

    Arguments   
    =========   

    LTRANL  (input) LOGICAL   
            On entry, LTRANL specifies the op(TL):   
               = .FALSE., op(TL) = TL,   
               = .TRUE., op(TL) = TL'.   

    LTRANR  (input) LOGICAL   
            On entry, LTRANR specifies the op(TR):   
              = .FALSE., op(TR) = TR,   
              = .TRUE., op(TR) = TR'.   

    ISGN    (input) INTEGER   
            On entry, ISGN specifies the sign of the equation   
            as described before. ISGN may only be 1 or -1.   

    N1      (input) INTEGER   
            On entry, N1 specifies the order of matrix TL.   
            N1 may only be 0, 1 or 2.   

    N2      (input) INTEGER   
            On entry, N2 specifies the order of matrix TR.   
            N2 may only be 0, 1 or 2.   

    TL      (input) DOUBLE PRECISION array, dimension (LDTL,2)   
            On entry, TL contains an N1 by N1 matrix.   

    LDTL    (input) INTEGER   
            The leading dimension of the matrix TL. LDTL >= max(1,N1).   

    TR      (input) DOUBLE PRECISION array, dimension (LDTR,2)   
            On entry, TR contains an N2 by N2 matrix.   

    LDTR    (input) INTEGER   
            The leading dimension of the matrix TR. LDTR >= max(1,N2).   

    B       (input) DOUBLE PRECISION array, dimension (LDB,2)   
            On entry, the N1 by N2 matrix B contains the right-hand   
            side of the equation.   

    LDB     (input) INTEGER   
            The leading dimension of the matrix B. LDB >= max(1,N1).   

    SCALE   (output) DOUBLE PRECISION   
            On exit, SCALE contains the scale factor. SCALE is chosen   
            less than or equal to 1 to prevent the solution overflowing. 
  

    X       (output) DOUBLE PRECISION array, dimension (LDX,2)   
            On exit, X contains the N1 by N2 solution.   

    LDX     (input) INTEGER   
            The leading dimension of the matrix X. LDX >= max(1,N1).   

    XNORM   (output) DOUBLE PRECISION   
            On exit, XNORM is the infinity-norm of the solution.   

    INFO    (output) INTEGER   
            On exit, INFO is set to   
               0: successful exit.   
               1: TL and TR have too close eigenvalues, so TL or   
                  TR is perturbed to get a nonsingular equation.   
            NOTE: In the interests of speed, this routine does not   
                  check the inputs for errors.   

   ===================================================================== 
  

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static integer c__4 = 4;
    static integer c__1 = 1;
    static integer c__16 = 16;
    static integer c__0 = 0;
    
    /* Initialized data */
    static integer locu12[4] = { 3,4,1,2 };
    static integer locl21[4] = { 2,1,4,3 };
    static integer locu22[4] = { 4,3,2,1 };
    static logical xswpiv[4] = { FALSE_,FALSE_,TRUE_,TRUE_ };
    static logical bswpiv[4] = { FALSE_,TRUE_,FALSE_,TRUE_ };
    /* System generated locals */
    integer b_dim1, b_offset, tl_dim1, tl_offset, tr_dim1, tr_offset, x_dim1, 
	    x_offset;
    doublereal d__1, d__2, d__3, d__4, d__5, d__6, d__7, d__8;
    /* Local variables */
    static doublereal btmp[4], smin;
    static integer ipiv;
    static doublereal temp;
    static integer jpiv[4];
    static doublereal xmax;
    static integer ipsv, jpsv, i, j, k;
    static logical bswap;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dswap_(integer *, doublereal *, integer 
	    *, doublereal *, integer *);
    static logical xswap;
    static doublereal x2[2], l21, u11, u12;
    static integer ip, jp;
    static doublereal u22, t16[16]	/* was [4][4] */;
    extern doublereal dlamch_(char *);
    extern integer idamax_(integer *, doublereal *, integer *);
    static doublereal smlnum, gam, bet, eps, sgn, tmp[4], tau1;



#define LOCU12(I) locu12[(I)]
#define LOCL21(I) locl21[(I)]
#define LOCU22(I) locu22[(I)]
#define BTMP(I) btmp[(I)]
#define JPIV(I) jpiv[(I)]
#define X2(I) x2[(I)]
#define TMP(I) tmp[(I)]

#define TL(I,J) tl[(I)-1 + ((J)-1)* ( *ldtl)]
#define TR(I,J) tr[(I)-1 + ((J)-1)* ( *ldtr)]
#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]
#define X(I,J) x[(I)-1 + ((J)-1)* ( *ldx)]


/*     Do not check the input parameters for errors */

    *info = 0;

/*     Quick return if possible */

    if (*n1 == 0 || *n2 == 0) {
	return 0;
    }

/*     Set constants to control overflow */

    eps = dlamch_("P");
    smlnum = dlamch_("S") / eps;
    sgn = (doublereal) (*isgn);

    k = *n1 + *n1 + *n2 - 2;
    switch (k) {
	case 1:  goto L10;
	case 2:  goto L20;
	case 3:  goto L30;
	case 4:  goto L50;
    }

/*     1 by 1: TL11*X + SGN*X*TR11 = B11 */

L10:
    tau1 = TL(1,1) + sgn * TR(1,1);
    bet = abs(tau1);
    if (bet <= smlnum) {
	tau1 = smlnum;
	bet = smlnum;
	*info = 1;
    }

    *scale = 1.;
    gam = (d__1 = B(1,1), abs(d__1));
    if (smlnum * gam > bet) {
	*scale = 1. / gam;
    }

    X(1,1) = B(1,1) * *scale / tau1;
    *xnorm = (d__1 = X(1,1), abs(d__1));
    return 0;

/*     1 by 2:   
       TL11*[X11 X12] + ISGN*[X11 X12]*op[TR11 TR12]  = [B11 B12]   
                                         [TR21 TR22] */

L20:

/* Computing MAX   
   Computing MAX */
    d__7 = (d__1 = TL(1,1), abs(d__1)), d__8 = (d__2 = TR(1,1)
	    , abs(d__2)), d__7 = max(d__7,d__8), d__8 = (d__3 = TR(1,2), abs(d__3)), d__7 = max(d__7,d__8), d__8 = (d__4 = TR(2,1), abs(d__4)), d__7 = max(d__7,d__8), d__8 = (d__5 = 
	    TR(2,2), abs(d__5));
    d__6 = eps * max(d__7,d__8);
    smin = max(d__6,smlnum);
    TMP(0) = TL(1,1) + sgn * TR(1,1);
    TMP(3) = TL(1,1) + sgn * TR(2,2);
    if (*ltranr) {
	TMP(1) = sgn * TR(2,1);
	TMP(2) = sgn * TR(1,2);
    } else {
	TMP(1) = sgn * TR(1,2);
	TMP(2) = sgn * TR(2,1);
    }
    BTMP(0) = B(1,1);
    BTMP(1) = B(1,2);
    goto L40;

/*     2 by 1:   
            op[TL11 TL12]*[X11] + ISGN* [X11]*TR11  = [B11]   
              [TL21 TL22] [X21]         [X21]         [B21] */

L30:
/* Computing MAX   
   Computing MAX */
    d__7 = (d__1 = TR(1,1), abs(d__1)), d__8 = (d__2 = TL(1,1)
	    , abs(d__2)), d__7 = max(d__7,d__8), d__8 = (d__3 = TL(1,2), abs(d__3)), d__7 = max(d__7,d__8), d__8 = (d__4 = TL(2,1), abs(d__4)), d__7 = max(d__7,d__8), d__8 = (d__5 = 
	    TL(2,2), abs(d__5));
    d__6 = eps * max(d__7,d__8);
    smin = max(d__6,smlnum);
    TMP(0) = TL(1,1) + sgn * TR(1,1);
    TMP(3) = TL(2,2) + sgn * TR(1,1);
    if (*ltranl) {
	TMP(1) = TL(1,2);
	TMP(2) = TL(2,1);
    } else {
	TMP(1) = TL(2,1);
	TMP(2) = TL(1,2);
    }
    BTMP(0) = B(1,1);
    BTMP(1) = B(2,1);
L40:

/*     Solve 2 by 2 system using complete pivoting.   
       Set pivots less than SMIN to SMIN. */

    ipiv = idamax_(&c__4, tmp, &c__1);
    u11 = TMP(ipiv - 1);
    if (abs(u11) <= smin) {
	*info = 1;
	u11 = smin;
    }
    u12 = TMP(LOCU12(ipiv - 1) - 1);
    l21 = TMP(LOCL21(ipiv - 1) - 1) / u11;
    u22 = TMP(LOCU22(ipiv - 1) - 1) - u12 * l21;
    xswap = xswpiv[ipiv - 1];
    bswap = bswpiv[ipiv - 1];
    if (abs(u22) <= smin) {
	*info = 1;
	u22 = smin;
    }
    if (bswap) {
	temp = BTMP(1);
	BTMP(1) = BTMP(0) - l21 * temp;
	BTMP(0) = temp;
    } else {
	BTMP(1) -= l21 * BTMP(0);
    }
    *scale = 1.;
    if (smlnum * 2. * abs(BTMP(1)) > abs(u22) || smlnum * 2. * abs(BTMP(0)) > 
	    abs(u11)) {
/* Computing MAX */
	d__1 = abs(BTMP(0)), d__2 = abs(BTMP(1));
	*scale = .5 / max(d__1,d__2);
	BTMP(0) *= *scale;
	BTMP(1) *= *scale;
    }
    X2(1) = BTMP(1) / u22;
    X2(0) = BTMP(0) / u11 - u12 / u11 * X2(1);
    if (xswap) {
	temp = X2(1);
	X2(1) = X2(0);
	X2(0) = temp;
    }
    X(1,1) = X2(0);
    if (*n1 == 1) {
	X(1,2) = X2(1);
	*xnorm = (d__1 = X(1,1), abs(d__1)) + (d__2 = X(1,2), abs(d__2));
    } else {
	X(2,1) = X2(1);
/* Computing MAX */
	d__3 = (d__1 = X(1,1), abs(d__1)), d__4 = (d__2 = X(2,1)
		, abs(d__2));
	*xnorm = max(d__3,d__4);
    }
    return 0;

/*     2 by 2:   
       op[TL11 TL12]*[X11 X12] +ISGN* [X11 X12]*op[TR11 TR12] = [B11 B12] 
  
         [TL21 TL22] [X21 X22]        [X21 X22]   [TR21 TR22]   [B21 B22] 
  

       Solve equivalent 4 by 4 system using complete pivoting.   
       Set pivots less than SMIN to SMIN. */

L50:
/* Computing MAX */
    d__5 = (d__1 = TR(1,1), abs(d__1)), d__6 = (d__2 = TR(1,2), abs(d__2)), d__5 = max(d__5,d__6), d__6 = (d__3 = TR(2,1), abs(d__3)), d__5 = max(d__5,d__6), d__6 = (d__4 = 
	    TR(2,2), abs(d__4));
    smin = max(d__5,d__6);
/* Computing MAX */
    d__5 = smin, d__6 = (d__1 = TL(1,1), abs(d__1)), d__5 = max(d__5,
	    d__6), d__6 = (d__2 = TL(1,2), abs(d__2)), d__5 = 
	    max(d__5,d__6), d__6 = (d__3 = TL(2,1), abs(d__3)), d__5 =
	     max(d__5,d__6), d__6 = (d__4 = TL(2,2), abs(d__4))
	    ;
    smin = max(d__5,d__6);
/* Computing MAX */
    d__1 = eps * smin;
    smin = max(d__1,smlnum);
    BTMP(0) = 0.;
    dcopy_(&c__16, btmp, &c__0, t16, &c__1);
    t16[0] = TL(1,1) + sgn * TR(1,1);
    t16[5] = TL(2,2) + sgn * TR(1,1);
    t16[10] = TL(1,1) + sgn * TR(2,2);
    t16[15] = TL(2,2) + sgn * TR(2,2);
    if (*ltranl) {
	t16[4] = TL(2,1);
	t16[1] = TL(1,2);
	t16[14] = TL(2,1);
	t16[11] = TL(1,2);
    } else {
	t16[4] = TL(1,2);
	t16[1] = TL(2,1);
	t16[14] = TL(1,2);
	t16[11] = TL(2,1);
    }
    if (*ltranr) {
	t16[8] = sgn * TR(1,2);
	t16[13] = sgn * TR(1,2);
	t16[2] = sgn * TR(2,1);
	t16[7] = sgn * TR(2,1);
    } else {
	t16[8] = sgn * TR(2,1);
	t16[13] = sgn * TR(2,1);
	t16[2] = sgn * TR(1,2);
	t16[7] = sgn * TR(1,2);
    }
    BTMP(0) = B(1,1);
    BTMP(1) = B(2,1);
    BTMP(2) = B(1,2);
    BTMP(3) = B(2,2);

/*     Perform elimination */

    for (i = 1; i <= 3; ++i) {
	xmax = 0.;
	for (ip = i; ip <= 4; ++ip) {
	    for (jp = i; jp <= 4; ++jp) {
		if ((d__1 = t16[ip + (jp << 2) - 5], abs(d__1)) >= xmax) {
		    xmax = (d__1 = t16[ip + (jp << 2) - 5], abs(d__1));
		    ipsv = ip;
		    jpsv = jp;
		}
/* L60: */
	    }
/* L70: */
	}
	if (ipsv != i) {
	    dswap_(&c__4, &t16[ipsv - 1], &c__4, &t16[i - 1], &c__4);
	    temp = BTMP(i - 1);
	    BTMP(i - 1) = BTMP(ipsv - 1);
	    BTMP(ipsv - 1) = temp;
	}
	if (jpsv != i) {
	    dswap_(&c__4, &t16[(jpsv << 2) - 4], &c__1, &t16[(i << 2) - 4], &
		    c__1);
	}
	JPIV(i - 1) = jpsv;
	if ((d__1 = t16[i + (i << 2) - 5], abs(d__1)) < smin) {
	    *info = 1;
	    t16[i + (i << 2) - 5] = smin;
	}
	for (j = i + 1; j <= 4; ++j) {
	    t16[j + (i << 2) - 5] /= t16[i + (i << 2) - 5];
	    BTMP(j - 1) -= t16[j + (i << 2) - 5] * BTMP(i - 1);
	    for (k = i + 1; k <= 4; ++k) {
		t16[j + (k << 2) - 5] -= t16[j + (i << 2) - 5] * t16[i + (k <<
			 2) - 5];
/* L80: */
	    }
/* L90: */
	}
/* L100: */
    }
    if (abs(t16[15]) < smin) {
	t16[15] = smin;
    }
    *scale = 1.;
    if (smlnum * 8. * abs(BTMP(0)) > abs(t16[0]) || smlnum * 8. * abs(BTMP(1))
	     > abs(t16[5]) || smlnum * 8. * abs(BTMP(2)) > abs(t16[10]) || 
	    smlnum * 8. * abs(BTMP(3)) > abs(t16[15])) {
/* Computing MAX */
	d__1 = abs(BTMP(0)), d__2 = abs(BTMP(1)), d__1 = max(d__1,d__2), d__2 
		= abs(BTMP(2)), d__1 = max(d__1,d__2), d__2 = abs(BTMP(3));
	*scale = .125 / max(d__1,d__2);
	BTMP(0) *= *scale;
	BTMP(1) *= *scale;
	BTMP(2) *= *scale;
	BTMP(3) *= *scale;
    }
    for (i = 1; i <= 4; ++i) {
	k = 5 - i;
	temp = 1. / t16[k + (k << 2) - 5];
	TMP(k - 1) = BTMP(k - 1) * temp;
	for (j = k + 1; j <= 4; ++j) {
	    TMP(k - 1) -= temp * t16[k + (j << 2) - 5] * TMP(j - 1);
/* L110: */
	}
/* L120: */
    }
    for (i = 1; i <= 3; ++i) {
	if (JPIV(4 - i - 1) != 4 - i) {
	    temp = TMP(4 - i - 1);
	    TMP(4 - i - 1) = TMP(JPIV(4 - i - 1) - 1);
	    TMP(JPIV(4 - i - 1) - 1) = temp;
	}
/* L130: */
    }
    X(1,1) = TMP(0);
    X(2,1) = TMP(1);
    X(1,2) = TMP(2);
    X(2,2) = TMP(3);
/* Computing MAX */
    d__1 = abs(TMP(0)) + abs(TMP(2)), d__2 = abs(TMP(1)) + abs(TMP(3));
    *xnorm = max(d__1,d__2);
    return 0;

/*     End of DLASY2 */

} /* dlasy2_ */

