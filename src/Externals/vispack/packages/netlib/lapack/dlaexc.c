#include "f2c.h"

/* Subroutine */ int dlaexc_(logical *wantq, integer *n, doublereal *t, 
	integer *ldt, doublereal *q, integer *ldq, integer *j1, integer *n1, 
	integer *n2, doublereal *work, integer *info)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   


    Purpose   
    =======   

    DLAEXC swaps adjacent diagonal blocks T11 and T22 of order 1 or 2 in 
  
    an upper quasi-triangular matrix T by an orthogonal similarity   
    transformation.   

    T must be in Schur canonical form, that is, block upper triangular   
    with 1-by-1 and 2-by-2 diagonal blocks; each 2-by-2 diagonal block   
    has its diagonal elemnts equal and its off-diagonal elements of   
    opposite sign.   

    Arguments   
    =========   

    WANTQ   (input) LOGICAL   
            = .TRUE. : accumulate the transformation in the matrix Q;   
            = .FALSE.: do not accumulate the transformation.   

    N       (input) INTEGER   
            The order of the matrix T. N >= 0.   

    T       (input/output) DOUBLE PRECISION array, dimension (LDT,N)   
            On entry, the upper quasi-triangular matrix T, in Schur   
            canonical form.   
            On exit, the updated matrix T, again in Schur canonical form. 
  

    LDT     (input)  INTEGER   
            The leading dimension of the array T. LDT >= max(1,N).   

    Q       (input/output) DOUBLE PRECISION array, dimension (LDQ,N)   
            On entry, if WANTQ is .TRUE., the orthogonal matrix Q.   
            On exit, if WANTQ is .TRUE., the updated matrix Q.   
            If WANTQ is .FALSE., Q is not referenced.   

    LDQ     (input) INTEGER   
            The leading dimension of the array Q.   
            LDQ >= 1; and if WANTQ is .TRUE., LDQ >= N.   

    J1      (input) INTEGER   
            The index of the first row of the first block T11.   

    N1      (input) INTEGER   
            The order of the first block T11. N1 = 0, 1 or 2.   

    N2      (input) INTEGER   
            The order of the second block T22. N2 = 0, 1 or 2.   

    WORK    (workspace) DOUBLE PRECISION array, dimension (N)   

    INFO    (output) INTEGER   
            = 0: successful exit   
            = 1: the transformed matrix T would be too far from Schur   
                 form; the blocks are not swapped and T and Q are   
                 unchanged.   

    ===================================================================== 
  


    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static integer c__1 = 1;
    static integer c__4 = 4;
    static logical c_false = FALSE_;
    static integer c_n1 = -1;
    static integer c__2 = 2;
    static integer c__3 = 3;
    
    /* System generated locals */
    integer q_dim1, q_offset, t_dim1, t_offset, i__1;
    doublereal d__1, d__2, d__3;
    /* Local variables */
    static integer ierr;
    static doublereal temp;
    extern /* Subroutine */ int drot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *);
    static doublereal d[16]	/* was [4][4] */;
    static integer k;
    static doublereal u[3], scale, x[4]	/* was [2][2] */, dnorm;
    static integer j2, j3, j4;
    static doublereal xnorm, u1[3], u2[3];
    extern /* Subroutine */ int dlanv2_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *), dlasy2_(
	    logical *, logical *, integer *, integer *, integer *, doublereal 
	    *, integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *);
    static integer nd;
    static doublereal cs, t11, t22;
    extern doublereal dlamch_(char *);
    static doublereal t33;
    extern doublereal dlange_(char *, integer *, integer *, doublereal *, 
	    integer *, doublereal *);
    extern /* Subroutine */ int dlarfg_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *);
    static doublereal sn;
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *), 
	    dlartg_(doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *), dlarfx_(char *, integer *, integer *, doublereal *,
	     doublereal *, doublereal *, integer *, doublereal *);
    static doublereal thresh, smlnum, wi1, wi2, wr1, wr2, eps, tau, tau1, 
	    tau2;



#define D(I) d[(I)]
#define WAS(I) was[(I)]
#define U(I) u[(I)]
#define U1(I) u1[(I)]
#define U2(I) u2[(I)]
#define WORK(I) work[(I)-1]

#define T(I,J) t[(I)-1 + ((J)-1)* ( *ldt)]
#define Q(I,J) q[(I)-1 + ((J)-1)* ( *ldq)]

    *info = 0;

/*     Quick return if possible */

    if (*n == 0 || *n1 == 0 || *n2 == 0) {
	return 0;
    }
    if (*j1 + *n1 > *n) {
	return 0;
    }

    j2 = *j1 + 1;
    j3 = *j1 + 2;
    j4 = *j1 + 3;

    if (*n1 == 1 && *n2 == 1) {

/*        Swap two 1-by-1 blocks. */

	t11 = T(*j1,*j1);
	t22 = T(j2,j2);

/*        Determine the transformation to perform the interchange. */

	d__1 = t22 - t11;
	dlartg_(&T(*j1,j2), &d__1, &cs, &sn, &temp);

/*        Apply transformation to the matrix T. */

	if (j3 <= *n) {
	    i__1 = *n - *j1 - 1;
	    drot_(&i__1, &T(*j1,j3), ldt, &T(j2,j3), 
		    ldt, &cs, &sn);
	}
	i__1 = *j1 - 1;
	drot_(&i__1, &T(1,*j1), &c__1, &T(1,j2), &c__1, 
		&cs, &sn);

	T(*j1,*j1) = t22;
	T(j2,j2) = t11;

	if (*wantq) {

/*           Accumulate transformation in the matrix Q. */

	    drot_(n, &Q(1,*j1), &c__1, &Q(1,j2), &c__1, 
		    &cs, &sn);
	}

    } else {

/*        Swapping involves at least one 2-by-2 block.   

          Copy the diagonal block of order N1+N2 to the local array D 
  
          and compute its norm. */

	nd = *n1 + *n2;
	dlacpy_("Full", &nd, &nd, &T(*j1,*j1), ldt, d, &c__4);
	dnorm = dlange_("Max", &nd, &nd, d, &c__4, &WORK(1));

/*        Compute machine-dependent threshold for test for accepting 
  
          swap. */

	eps = dlamch_("P");
	smlnum = dlamch_("S") / eps;
/* Computing MAX */
	d__1 = eps * 10. * dnorm;
	thresh = max(d__1,smlnum);

/*        Solve T11*X - X*T22 = scale*T12 for X. */

	dlasy2_(&c_false, &c_false, &c_n1, n1, n2, d, &c__4, &D(*n1 + 1 + (*
		n1 + 1 << 2) - 5), &c__4, &D((*n1 + 1 << 2) - 4), &c__4, &
		scale, x, &c__2, &xnorm, &ierr);

/*        Swap the adjacent diagonal blocks. */

	k = *n1 + *n1 + *n2 - 3;
	switch (k) {
	    case 1:  goto L10;
	    case 2:  goto L20;
	    case 3:  goto L30;
	}

L10:

/*        N1 = 1, N2 = 2: generate elementary reflector H so that:   

          ( scale, X11, X12 ) H = ( 0, 0, * ) */

	U(0) = scale;
	U(1) = x[0];
	U(2) = x[2];
	dlarfg_(&c__3, &U(2), u, &c__1, &tau);
	U(2) = 1.;
	t11 = T(*j1,*j1);

/*        Perform swap provisionally on diagonal block in D. */

	dlarfx_("L", &c__3, &c__3, u, &tau, d, &c__4, &WORK(1));
	dlarfx_("R", &c__3, &c__3, u, &tau, d, &c__4, &WORK(1));

/*        Test whether to reject swap.   

   Computing MAX */
	d__2 = abs(D(2)), d__3 = abs(D(6)), d__2 = max(d__2,d__3), d__3 = (
		d__1 = D(10) - t11, abs(d__1));
	if (max(d__2,d__3) > thresh) {
	    goto L50;
	}

/*        Accept swap: apply transformation to the entire matrix T. */

	i__1 = *n - *j1 + 1;
	dlarfx_("L", &c__3, &i__1, u, &tau, &T(*j1,*j1), ldt, &
		WORK(1));
	dlarfx_("R", &j2, &c__3, u, &tau, &T(1,*j1), ldt, &WORK(1));

	T(j3,*j1) = 0.;
	T(j3,j2) = 0.;
	T(j3,j3) = t11;

	if (*wantq) {

/*           Accumulate transformation in the matrix Q. */

	    dlarfx_("R", n, &c__3, u, &tau, &Q(1,*j1), ldq, &WORK(
		    1));
	}
	goto L40;

L20:

/*        N1 = 2, N2 = 1: generate elementary reflector H so that:   

          H (  -X11 ) = ( * )   
            (  -X21 ) = ( 0 )   
            ( scale ) = ( 0 ) */

	U(0) = -x[0];
	U(1) = -x[1];
	U(2) = scale;
	dlarfg_(&c__3, u, &U(1), &c__1, &tau);
	U(0) = 1.;
	t33 = T(j3,j3);

/*        Perform swap provisionally on diagonal block in D. */

	dlarfx_("L", &c__3, &c__3, u, &tau, d, &c__4, &WORK(1));
	dlarfx_("R", &c__3, &c__3, u, &tau, d, &c__4, &WORK(1));

/*        Test whether to reject swap.   

   Computing MAX */
	d__2 = abs(D(1)), d__3 = abs(D(2)), d__2 = max(d__2,d__3), d__3 = (
		d__1 = D(0) - t33, abs(d__1));
	if (max(d__2,d__3) > thresh) {
	    goto L50;
	}

/*        Accept swap: apply transformation to the entire matrix T. */

	dlarfx_("R", &j3, &c__3, u, &tau, &T(1,*j1), ldt, &WORK(1));
	i__1 = *n - *j1;
	dlarfx_("L", &c__3, &i__1, u, &tau, &T(*j1,j2), ldt, &WORK(
		1));

	T(*j1,*j1) = t33;
	T(j2,*j1) = 0.;
	T(j3,*j1) = 0.;

	if (*wantq) {

/*           Accumulate transformation in the matrix Q. */

	    dlarfx_("R", n, &c__3, u, &tau, &Q(1,*j1), ldq, &WORK(
		    1));
	}
	goto L40;

L30:

/*        N1 = 2, N2 = 2: generate elementary reflectors H(1) and H(2)
 so   
          that:   

          H(2) H(1) (  -X11  -X12 ) = (  *  * )   
                    (  -X21  -X22 )   (  0  * )   
                    ( scale    0  )   (  0  0 )   
                    (    0  scale )   (  0  0 ) */

	U1(0) = -x[0];
	U1(1) = -x[1];
	U1(2) = scale;
	dlarfg_(&c__3, u1, &U1(1), &c__1, &tau1);
	U1(0) = 1.;

	temp = -tau1 * (x[2] + U1(1) * x[3]);
	U2(0) = -temp * U1(1) - x[3];
	U2(1) = -temp * U1(2);
	U2(2) = scale;
	dlarfg_(&c__3, u2, &U2(1), &c__1, &tau2);
	U2(0) = 1.;

/*        Perform swap provisionally on diagonal block in D. */

	dlarfx_("L", &c__3, &c__4, u1, &tau1, d, &c__4, &WORK(1));
	dlarfx_("R", &c__4, &c__3, u1, &tau1, d, &c__4, &WORK(1));
	dlarfx_("L", &c__3, &c__4, u2, &tau2, &D(1), &c__4, &WORK(1));
	dlarfx_("R", &c__4, &c__3, u2, &tau2, &D(4), &c__4, &WORK(1));

/*        Test whether to reject swap.   

   Computing MAX */
	d__1 = abs(D(2)), d__2 = abs(D(6)), d__1 = max(d__1,d__2), d__2 = abs(
		D(3)), d__1 = max(d__1,d__2), d__2 = abs(D(7));
	if (max(d__1,d__2) > thresh) {
	    goto L50;
	}

/*        Accept swap: apply transformation to the entire matrix T. */

	i__1 = *n - *j1 + 1;
	dlarfx_("L", &c__3, &i__1, u1, &tau1, &T(*j1,*j1), ldt, &
		WORK(1));
	dlarfx_("R", &j4, &c__3, u1, &tau1, &T(1,*j1), ldt, &WORK(
		1));
	i__1 = *n - *j1 + 1;
	dlarfx_("L", &c__3, &i__1, u2, &tau2, &T(j2,*j1), ldt, &
		WORK(1));
	dlarfx_("R", &j4, &c__3, u2, &tau2, &T(1,j2), ldt, &WORK(1)
		);

	T(j3,*j1) = 0.;
	T(j3,j2) = 0.;
	T(j4,*j1) = 0.;
	T(j4,j2) = 0.;

	if (*wantq) {

/*           Accumulate transformation in the matrix Q. */

	    dlarfx_("R", n, &c__3, u1, &tau1, &Q(1,*j1), ldq, &
		    WORK(1));
	    dlarfx_("R", n, &c__3, u2, &tau2, &Q(1,j2), ldq, &WORK(
		    1));
	}

L40:

	if (*n2 == 2) {

/*           Standardize new 2-by-2 block T11 */

	    dlanv2_(&T(*j1,*j1), &T(*j1,j2), &T(j2,*j1), &T(j2,j2), &wr1, &wi1, &wr2, &
		    wi2, &cs, &sn);
	    i__1 = *n - *j1 - 1;
	    drot_(&i__1, &T(*j1,*j1+2), ldt, &T(j2,*j1+2), ldt, &cs, &sn);
	    i__1 = *j1 - 1;
	    drot_(&i__1, &T(1,*j1), &c__1, &T(1,j2), &
		    c__1, &cs, &sn);
	    if (*wantq) {
		drot_(n, &Q(1,*j1), &c__1, &Q(1,j2), &
			c__1, &cs, &sn);
	    }
	}

	if (*n1 == 2) {

/*           Standardize new 2-by-2 block T22 */

	    j3 = *j1 + *n2;
	    j4 = j3 + 1;
	    dlanv2_(&T(j3,j3), &T(j3,j4), &T(j4,j3), &T(j4,j4), &wr1, &wi1, &wr2, &wi2, &
		    cs, &sn);
	    if (j3 + 2 <= *n) {
		i__1 = *n - j3 - 1;
		drot_(&i__1, &T(j3,j3+2), ldt, &T(j4,j3+2), ldt, &cs, &sn);
	    }
	    i__1 = j3 - 1;
	    drot_(&i__1, &T(1,j3), &c__1, &T(1,j4), &
		    c__1, &cs, &sn);
	    if (*wantq) {
		drot_(n, &Q(1,j3), &c__1, &Q(1,j4), &
			c__1, &cs, &sn);
	    }
	}

    }
    return 0;

/*     Exit with INFO = 1 if swap was rejected. */

L50:
    *info = 1;
    return 0;

/*     End of DLAEXC */

} /* dlaexc_ */

