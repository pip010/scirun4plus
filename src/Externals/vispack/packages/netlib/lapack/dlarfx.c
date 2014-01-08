#include "f2c.h"

/* Subroutine */ int dlarfx_(char *side, integer *m, integer *n, doublereal *
	v, doublereal *tau, doublereal *c, integer *ldc, doublereal *work)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   


    Purpose   
    =======   

    DLARFX applies a real elementary reflector H to a real m by n   
    matrix C, from either the left or the right. H is represented in the 
  
    form   

          H = I - tau * v * v'   

    where tau is a real scalar and v is a real vector.   

    If tau = 0, then H is taken to be the unit matrix   

    This version uses inline code if H has order < 11.   

    Arguments   
    =========   

    SIDE    (input) CHARACTER*1   
            = 'L': form  H * C   
            = 'R': form  C * H   

    M       (input) INTEGER   
            The number of rows of the matrix C.   

    N       (input) INTEGER   
            The number of columns of the matrix C.   

    V       (input) DOUBLE PRECISION array, dimension (M) if SIDE = 'L'   
                                       or (N) if SIDE = 'R'   
            The vector v in the representation of H.   

    TAU     (input) DOUBLE PRECISION   
            The value tau in the representation of H.   

    C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)   
            On entry, the m by n matrix C.   
            On exit, C is overwritten by the matrix H * C if SIDE = 'L', 
  
            or C * H if SIDE = 'R'.   

    LDC     (input) INTEGER   
            The leading dimension of the array C. LDA >= (1,M).   

    WORK    (workspace) DOUBLE PRECISION array, dimension   
                        (N) if SIDE = 'L'   
                        or (M) if SIDE = 'R'   
            WORK is not referenced if H has order < 11.   

    ===================================================================== 
  


    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static doublereal c_b14 = 1.;
    static integer c__1 = 1;
    static doublereal c_b16 = 0.;
    
    /* System generated locals */
    integer c_dim1, c_offset, i__1;
    doublereal d__1;
    /* Local variables */
    extern /* Subroutine */ int dger_(integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static integer j;
    extern logical lsame_(char *, char *);
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *);
    static doublereal t1, t2, t3, t4, t5, t6, t7, t8, t9, v1, v2, v3, v4, v5, 
	    v6, v7, v8, v9, t10, v10, sum;



#define V(I) v[(I)-1]
#define WORK(I) work[(I)-1]

#define C(I,J) c[(I)-1 + ((J)-1)* ( *ldc)]

    if (*tau == 0.) {
	return 0;
    }
    if (lsame_(side, "L")) {

/*        Form  H * C, where H has order m. */

	switch (*m) {
	    case 1:  goto L10;
	    case 2:  goto L30;
	    case 3:  goto L50;
	    case 4:  goto L70;
	    case 5:  goto L90;
	    case 6:  goto L110;
	    case 7:  goto L130;
	    case 8:  goto L150;
	    case 9:  goto L170;
	    case 10:  goto L190;
	}

/*        Code for general M   

          w := C'*v */

	dgemv_("Transpose", m, n, &c_b14, &C(1,1), ldc, &V(1), &c__1, &
		c_b16, &WORK(1), &c__1);

/*        C := C - tau * v * w' */

	d__1 = -(*tau);
	dger_(m, n, &d__1, &V(1), &c__1, &WORK(1), &c__1, &C(1,1), ldc);
	goto L410;
L10:

/*        Special code for 1 x 1 Householder */

	t1 = 1. - *tau * V(1) * V(1);
	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
	    C(1,j) = t1 * C(1,j);
/* L20: */
	}
	goto L410;
L30:

/*        Special code for 2 x 2 Householder */

	v1 = V(1);
	t1 = *tau * v1;
	v2 = V(2);
	t2 = *tau * v2;
	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
	    sum = v1 * C(1,j) + v2 * C(2,j);
	    C(1,j) -= sum * t1;
	    C(2,j) -= sum * t2;
/* L40: */
	}
	goto L410;
L50:

/*        Special code for 3 x 3 Householder */

	v1 = V(1);
	t1 = *tau * v1;
	v2 = V(2);
	t2 = *tau * v2;
	v3 = V(3);
	t3 = *tau * v3;
	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
	    sum = v1 * C(1,j) + v2 * C(2,j) + v3 * C(3,j);
	    C(1,j) -= sum * t1;
	    C(2,j) -= sum * t2;
	    C(3,j) -= sum * t3;
/* L60: */
	}
	goto L410;
L70:

/*        Special code for 4 x 4 Householder */

	v1 = V(1);
	t1 = *tau * v1;
	v2 = V(2);
	t2 = *tau * v2;
	v3 = V(3);
	t3 = *tau * v3;
	v4 = V(4);
	t4 = *tau * v4;
	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
	    sum = v1 * C(1,j) + v2 * C(2,j) + v3 * C(3,j) + v4 * C(4,j);
	    C(1,j) -= sum * t1;
	    C(2,j) -= sum * t2;
	    C(3,j) -= sum * t3;
	    C(4,j) -= sum * t4;
/* L80: */
	}
	goto L410;
L90:

/*        Special code for 5 x 5 Householder */

	v1 = V(1);
	t1 = *tau * v1;
	v2 = V(2);
	t2 = *tau * v2;
	v3 = V(3);
	t3 = *tau * v3;
	v4 = V(4);
	t4 = *tau * v4;
	v5 = V(5);
	t5 = *tau * v5;
	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
	    sum = v1 * C(1,j) + v2 * C(2,j) + v3 * C(3,j) + v4 * C(4,j) + v5 * C(5,j);
	    C(1,j) -= sum * t1;
	    C(2,j) -= sum * t2;
	    C(3,j) -= sum * t3;
	    C(4,j) -= sum * t4;
	    C(5,j) -= sum * t5;
/* L100: */
	}
	goto L410;
L110:

/*        Special code for 6 x 6 Householder */

	v1 = V(1);
	t1 = *tau * v1;
	v2 = V(2);
	t2 = *tau * v2;
	v3 = V(3);
	t3 = *tau * v3;
	v4 = V(4);
	t4 = *tau * v4;
	v5 = V(5);
	t5 = *tau * v5;
	v6 = V(6);
	t6 = *tau * v6;
	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
	    sum = v1 * C(1,j) + v2 * C(2,j) + v3 * C(3,j) + v4 * C(4,j) + v5 * C(5,j) + v6 * C(6,j);
	    C(1,j) -= sum * t1;
	    C(2,j) -= sum * t2;
	    C(3,j) -= sum * t3;
	    C(4,j) -= sum * t4;
	    C(5,j) -= sum * t5;
	    C(6,j) -= sum * t6;
/* L120: */
	}
	goto L410;
L130:

/*        Special code for 7 x 7 Householder */

	v1 = V(1);
	t1 = *tau * v1;
	v2 = V(2);
	t2 = *tau * v2;
	v3 = V(3);
	t3 = *tau * v3;
	v4 = V(4);
	t4 = *tau * v4;
	v5 = V(5);
	t5 = *tau * v5;
	v6 = V(6);
	t6 = *tau * v6;
	v7 = V(7);
	t7 = *tau * v7;
	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
	    sum = v1 * C(1,j) + v2 * C(2,j) + v3 * C(3,j) + v4 * C(4,j) + v5 * C(5,j) + v6 * C(6,j) + v7 * C(7,j);
	    C(1,j) -= sum * t1;
	    C(2,j) -= sum * t2;
	    C(3,j) -= sum * t3;
	    C(4,j) -= sum * t4;
	    C(5,j) -= sum * t5;
	    C(6,j) -= sum * t6;
	    C(7,j) -= sum * t7;
/* L140: */
	}
	goto L410;
L150:

/*        Special code for 8 x 8 Householder */

	v1 = V(1);
	t1 = *tau * v1;
	v2 = V(2);
	t2 = *tau * v2;
	v3 = V(3);
	t3 = *tau * v3;
	v4 = V(4);
	t4 = *tau * v4;
	v5 = V(5);
	t5 = *tau * v5;
	v6 = V(6);
	t6 = *tau * v6;
	v7 = V(7);
	t7 = *tau * v7;
	v8 = V(8);
	t8 = *tau * v8;
	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
	    sum = v1 * C(1,j) + v2 * C(2,j) + v3 * C(3,j) + v4 * C(4,j) + v5 * C(5,j) + v6 * C(6,j) + v7 * C(7,j) + 
		    v8 * C(8,j);
	    C(1,j) -= sum * t1;
	    C(2,j) -= sum * t2;
	    C(3,j) -= sum * t3;
	    C(4,j) -= sum * t4;
	    C(5,j) -= sum * t5;
	    C(6,j) -= sum * t6;
	    C(7,j) -= sum * t7;
	    C(8,j) -= sum * t8;
/* L160: */
	}
	goto L410;
L170:

/*        Special code for 9 x 9 Householder */

	v1 = V(1);
	t1 = *tau * v1;
	v2 = V(2);
	t2 = *tau * v2;
	v3 = V(3);
	t3 = *tau * v3;
	v4 = V(4);
	t4 = *tau * v4;
	v5 = V(5);
	t5 = *tau * v5;
	v6 = V(6);
	t6 = *tau * v6;
	v7 = V(7);
	t7 = *tau * v7;
	v8 = V(8);
	t8 = *tau * v8;
	v9 = V(9);
	t9 = *tau * v9;
	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
	    sum = v1 * C(1,j) + v2 * C(2,j) + v3 * C(3,j) + v4 * C(4,j) + v5 * C(5,j) + v6 * C(6,j) + v7 * C(7,j) + 
		    v8 * C(8,j) + v9 * C(9,j);
	    C(1,j) -= sum * t1;
	    C(2,j) -= sum * t2;
	    C(3,j) -= sum * t3;
	    C(4,j) -= sum * t4;
	    C(5,j) -= sum * t5;
	    C(6,j) -= sum * t6;
	    C(7,j) -= sum * t7;
	    C(8,j) -= sum * t8;
	    C(9,j) -= sum * t9;
/* L180: */
	}
	goto L410;
L190:

/*        Special code for 10 x 10 Householder */

	v1 = V(1);
	t1 = *tau * v1;
	v2 = V(2);
	t2 = *tau * v2;
	v3 = V(3);
	t3 = *tau * v3;
	v4 = V(4);
	t4 = *tau * v4;
	v5 = V(5);
	t5 = *tau * v5;
	v6 = V(6);
	t6 = *tau * v6;
	v7 = V(7);
	t7 = *tau * v7;
	v8 = V(8);
	t8 = *tau * v8;
	v9 = V(9);
	t9 = *tau * v9;
	v10 = V(10);
	t10 = *tau * v10;
	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
	    sum = v1 * C(1,j) + v2 * C(2,j) + v3 * C(3,j) + v4 * C(4,j) + v5 * C(5,j) + v6 * C(6,j) + v7 * C(7,j) + 
		    v8 * C(8,j) + v9 * C(9,j) + v10 * C(10,j);
	    C(1,j) -= sum * t1;
	    C(2,j) -= sum * t2;
	    C(3,j) -= sum * t3;
	    C(4,j) -= sum * t4;
	    C(5,j) -= sum * t5;
	    C(6,j) -= sum * t6;
	    C(7,j) -= sum * t7;
	    C(8,j) -= sum * t8;
	    C(9,j) -= sum * t9;
	    C(10,j) -= sum * t10;
/* L200: */
	}
	goto L410;
    } else {

/*        Form  C * H, where H has order n. */

	switch (*n) {
	    case 1:  goto L210;
	    case 2:  goto L230;
	    case 3:  goto L250;
	    case 4:  goto L270;
	    case 5:  goto L290;
	    case 6:  goto L310;
	    case 7:  goto L330;
	    case 8:  goto L350;
	    case 9:  goto L370;
	    case 10:  goto L390;
	}

/*        Code for general N   

          w := C * v */

	dgemv_("No transpose", m, n, &c_b14, &C(1,1), ldc, &V(1), &c__1, 
		&c_b16, &WORK(1), &c__1);

/*        C := C - tau * w * v' */

	d__1 = -(*tau);
	dger_(m, n, &d__1, &WORK(1), &c__1, &V(1), &c__1, &C(1,1), ldc);
	goto L410;
L210:

/*        Special code for 1 x 1 Householder */

	t1 = 1. - *tau * V(1) * V(1);
	i__1 = *m;
	for (j = 1; j <= *m; ++j) {
	    C(j,1) = t1 * C(j,1);
/* L220: */
	}
	goto L410;
L230:

/*        Special code for 2 x 2 Householder */

	v1 = V(1);
	t1 = *tau * v1;
	v2 = V(2);
	t2 = *tau * v2;
	i__1 = *m;
	for (j = 1; j <= *m; ++j) {
	    sum = v1 * C(j,1) + v2 * C(j,2);
	    C(j,1) -= sum * t1;
	    C(j,2) -= sum * t2;
/* L240: */
	}
	goto L410;
L250:

/*        Special code for 3 x 3 Householder */

	v1 = V(1);
	t1 = *tau * v1;
	v2 = V(2);
	t2 = *tau * v2;
	v3 = V(3);
	t3 = *tau * v3;
	i__1 = *m;
	for (j = 1; j <= *m; ++j) {
	    sum = v1 * C(j,1) + v2 * C(j,2) + v3 * C(j,3);
	    C(j,1) -= sum * t1;
	    C(j,2) -= sum * t2;
	    C(j,3) -= sum * t3;
/* L260: */
	}
	goto L410;
L270:

/*        Special code for 4 x 4 Householder */

	v1 = V(1);
	t1 = *tau * v1;
	v2 = V(2);
	t2 = *tau * v2;
	v3 = V(3);
	t3 = *tau * v3;
	v4 = V(4);
	t4 = *tau * v4;
	i__1 = *m;
	for (j = 1; j <= *m; ++j) {
	    sum = v1 * C(j,1) + v2 * C(j,2) + v3 * C(j,3) + v4 * C(j,4);
	    C(j,1) -= sum * t1;
	    C(j,2) -= sum * t2;
	    C(j,3) -= sum * t3;
	    C(j,4) -= sum * t4;
/* L280: */
	}
	goto L410;
L290:

/*        Special code for 5 x 5 Householder */

	v1 = V(1);
	t1 = *tau * v1;
	v2 = V(2);
	t2 = *tau * v2;
	v3 = V(3);
	t3 = *tau * v3;
	v4 = V(4);
	t4 = *tau * v4;
	v5 = V(5);
	t5 = *tau * v5;
	i__1 = *m;
	for (j = 1; j <= *m; ++j) {
	    sum = v1 * C(j,1) + v2 * C(j,2) + v3 * C(j,3) + v4 * C(j,4) + v5 * C(j,5);
	    C(j,1) -= sum * t1;
	    C(j,2) -= sum * t2;
	    C(j,3) -= sum * t3;
	    C(j,4) -= sum * t4;
	    C(j,5) -= sum * t5;
/* L300: */
	}
	goto L410;
L310:

/*        Special code for 6 x 6 Householder */

	v1 = V(1);
	t1 = *tau * v1;
	v2 = V(2);
	t2 = *tau * v2;
	v3 = V(3);
	t3 = *tau * v3;
	v4 = V(4);
	t4 = *tau * v4;
	v5 = V(5);
	t5 = *tau * v5;
	v6 = V(6);
	t6 = *tau * v6;
	i__1 = *m;
	for (j = 1; j <= *m; ++j) {
	    sum = v1 * C(j,1) + v2 * C(j,2) + v3 * C(j,3) + v4 * C(j,4) + v5 * C(j,5) + v6 * C(j,6);
	    C(j,1) -= sum * t1;
	    C(j,2) -= sum * t2;
	    C(j,3) -= sum * t3;
	    C(j,4) -= sum * t4;
	    C(j,5) -= sum * t5;
	    C(j,6) -= sum * t6;
/* L320: */
	}
	goto L410;
L330:

/*        Special code for 7 x 7 Householder */

	v1 = V(1);
	t1 = *tau * v1;
	v2 = V(2);
	t2 = *tau * v2;
	v3 = V(3);
	t3 = *tau * v3;
	v4 = V(4);
	t4 = *tau * v4;
	v5 = V(5);
	t5 = *tau * v5;
	v6 = V(6);
	t6 = *tau * v6;
	v7 = V(7);
	t7 = *tau * v7;
	i__1 = *m;
	for (j = 1; j <= *m; ++j) {
	    sum = v1 * C(j,1) + v2 * C(j,2) + v3 * C(j,3) + v4 * C(j,4) + v5 * C(j,5) + v6 * C(j,6) + v7 * C(j,7);
	    C(j,1) -= sum * t1;
	    C(j,2) -= sum * t2;
	    C(j,3) -= sum * t3;
	    C(j,4) -= sum * t4;
	    C(j,5) -= sum * t5;
	    C(j,6) -= sum * t6;
	    C(j,7) -= sum * t7;
/* L340: */
	}
	goto L410;
L350:

/*        Special code for 8 x 8 Householder */

	v1 = V(1);
	t1 = *tau * v1;
	v2 = V(2);
	t2 = *tau * v2;
	v3 = V(3);
	t3 = *tau * v3;
	v4 = V(4);
	t4 = *tau * v4;
	v5 = V(5);
	t5 = *tau * v5;
	v6 = V(6);
	t6 = *tau * v6;
	v7 = V(7);
	t7 = *tau * v7;
	v8 = V(8);
	t8 = *tau * v8;
	i__1 = *m;
	for (j = 1; j <= *m; ++j) {
	    sum = v1 * C(j,1) + v2 * C(j,2) + v3 * C(j,3) + v4 * C(j,4) + v5 * C(j,5) + v6 * C(j,6) + v7 * C(j,7) + v8 * C(j,8);
	    C(j,1) -= sum * t1;
	    C(j,2) -= sum * t2;
	    C(j,3) -= sum * t3;
	    C(j,4) -= sum * t4;
	    C(j,5) -= sum * t5;
	    C(j,6) -= sum * t6;
	    C(j,7) -= sum * t7;
	    C(j,8) -= sum * t8;
/* L360: */
	}
	goto L410;
L370:

/*        Special code for 9 x 9 Householder */

	v1 = V(1);
	t1 = *tau * v1;
	v2 = V(2);
	t2 = *tau * v2;
	v3 = V(3);
	t3 = *tau * v3;
	v4 = V(4);
	t4 = *tau * v4;
	v5 = V(5);
	t5 = *tau * v5;
	v6 = V(6);
	t6 = *tau * v6;
	v7 = V(7);
	t7 = *tau * v7;
	v8 = V(8);
	t8 = *tau * v8;
	v9 = V(9);
	t9 = *tau * v9;
	i__1 = *m;
	for (j = 1; j <= *m; ++j) {
	    sum = v1 * C(j,1) + v2 * C(j,2) + v3 * C(j,3) + v4 * C(j,4) + v5 * C(j,5) + v6 * C(j,6) + v7 * C(j,7) + v8 * C(j,8) + v9 * C(j,9);
	    C(j,1) -= sum * t1;
	    C(j,2) -= sum * t2;
	    C(j,3) -= sum * t3;
	    C(j,4) -= sum * t4;
	    C(j,5) -= sum * t5;
	    C(j,6) -= sum * t6;
	    C(j,7) -= sum * t7;
	    C(j,8) -= sum * t8;
	    C(j,9) -= sum * t9;
/* L380: */
	}
	goto L410;
L390:

/*        Special code for 10 x 10 Householder */

	v1 = V(1);
	t1 = *tau * v1;
	v2 = V(2);
	t2 = *tau * v2;
	v3 = V(3);
	t3 = *tau * v3;
	v4 = V(4);
	t4 = *tau * v4;
	v5 = V(5);
	t5 = *tau * v5;
	v6 = V(6);
	t6 = *tau * v6;
	v7 = V(7);
	t7 = *tau * v7;
	v8 = V(8);
	t8 = *tau * v8;
	v9 = V(9);
	t9 = *tau * v9;
	v10 = V(10);
	t10 = *tau * v10;
	i__1 = *m;
	for (j = 1; j <= *m; ++j) {
	    sum = v1 * C(j,1) + v2 * C(j,2) + v3 * C(j,3) + v4 * C(j,4) + v5 * C(j,5) + v6 * C(j,6) + v7 * C(j,7) + v8 * C(j,8) + v9 * C(j,9) 
		    + v10 * C(j,10);
	    C(j,1) -= sum * t1;
	    C(j,2) -= sum * t2;
	    C(j,3) -= sum * t3;
	    C(j,4) -= sum * t4;
	    C(j,5) -= sum * t5;
	    C(j,6) -= sum * t6;
	    C(j,7) -= sum * t7;
	    C(j,8) -= sum * t8;
	    C(j,9) -= sum * t9;
	    C(j,10) -= sum * t10;
/* L400: */
	}
	goto L410;
    }
L410:
    return 0;

/*     End of DLARFX */

} /* dlarfx_ */

