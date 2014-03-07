#include "f2c.h"

/* Subroutine */ int ctgevc_(char *side, char *howmny, logical *select, 
	integer *n, complex *a, integer *lda, complex *b, integer *ldb, 
	complex *vl, integer *ldvl, complex *vr, integer *ldvr, integer *mm, 
	integer *m, complex *work, real *rwork, integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   



    Purpose   
    =======   

    CTGEVC computes some or all of the right and/or left generalized   
    eigenvectors of a pair of complex upper triangular matrices (A,B).   

    The right generalized eigenvector x and the left generalized   
    eigenvector y of (A,B) corresponding to a generalized eigenvalue   
    w are defined by:   

            (A - wB) * x = 0  and  y**H * (A - wB) = 0   

    where y**H denotes the conjugate tranpose of y.   

    If an eigenvalue w is determined by zero diagonal elements of both A 
  
    and B, a unit vector is returned as the corresponding eigenvector.   

    If all eigenvectors are requested, the routine may either return   
    the matrices X and/or Y of right or left eigenvectors of (A,B), or   
    the products Z*X and/or Q*Y, where Z and Q are input unitary   
    matrices.  If (A,B) was obtained from the generalized Schur   
    factorization of an original pair of matrices   
       (A0,B0) = (Q*A*Z**H,Q*B*Z**H),   
    then Z*X and Q*Y are the matrices of right or left eigenvectors of   
    A.   

    Arguments   
    =========   

    SIDE    (input) CHARACTER*1   
            = 'R': compute right eigenvectors only;   
            = 'L': compute left eigenvectors only;   
            = 'B': compute both right and left eigenvectors.   

    HOWMNY  (input) CHARACTER*1   
            = 'A': compute all right and/or left eigenvectors;   
            = 'B': compute all right and/or left eigenvectors, and   
                   backtransform them using the input matrices supplied   
                   in VR and/or VL;   
            = 'S': compute selected right and/or left eigenvectors,   
                   specified by the logical array SELECT.   

    SELECT  (input) LOGICAL array, dimension (N)   
            If HOWMNY='S', SELECT specifies the eigenvectors to be   
            computed.   
            If HOWMNY='A' or 'B', SELECT is not referenced.   
            To select the eigenvector corresponding to the j-th   
            eigenvalue, SELECT(j) must be set to .TRUE..   

    N       (input) INTEGER   
            The order of the matrices A and B.  N >= 0.   

    A       (input) COMPLEX array, dimension (LDA,N)   
            The upper triangular matrix A.   

    LDA     (input) INTEGER   
            The leading dimension of array A.  LDA >= max(1,N).   

    B       (input) COMPLEX array, dimension (LDB,N)   
            The upper triangular matrix B.  B must have real diagonal   
            elements.   

    LDB     (input) INTEGER   
            The leading dimension of array B.  LDB >= max(1,N).   

    VL      (input/output) COMPLEX array, dimension (LDVL,MM)   
            On entry, if SIDE = 'L' or 'B' and HOWMNY = 'B', VL must   
            contain an N-by-N matrix Q (usually the unitary matrix Q   
            of left Schur vectors returned by CHGEQZ).   
            On exit, if SIDE = 'L' or 'B', VL contains:   
            if HOWMNY = 'A', the matrix Y of left eigenvectors of (A,B); 
  
            if HOWMNY = 'B', the matrix Q*Y;   
            if HOWMNY = 'S', the left eigenvectors of (A,B) specified by 
  
                        SELECT, stored consecutively in the columns of   
                        VL, in the same order as their eigenvalues.   
            If SIDE = 'R', VL is not referenced.   

    LDVL    (input) INTEGER   
            The leading dimension of array VL.   
            LDVL >= max(1,N) if SIDE = 'L' or 'B'; LDVL >= 1 otherwise.   

    VR      (input/output) COMPLEX array, dimension (LDVR,MM)   
            On entry, if SIDE = 'R' or 'B' and HOWMNY = 'B', VR must   
            contain an N-by-N matrix Q (usually the unitary matrix Z   
            of right Schur vectors returned by CHGEQZ).   
            On exit, if SIDE = 'R' or 'B', VR contains:   
            if HOWMNY = 'A', the matrix X of right eigenvectors of (A,B); 
  
            if HOWMNY = 'B', the matrix Z*X;   
            if HOWMNY = 'S', the right eigenvectors of (A,B) specified by 
  
                        SELECT, stored consecutively in the columns of   
                        VR, in the same order as their eigenvalues.   
            If SIDE = 'L', VR is not referenced.   

    LDVR    (input) INTEGER   
            The leading dimension of the array VR.   
            LDVR >= max(1,N) if SIDE = 'R' or 'B'; LDVR >= 1 otherwise.   

    MM      (input) INTEGER   
            The leading dimension of the array VR.   
            LDVR >= max(1,N) if SIDE = 'R' or 'B'; LDVR >= 1 otherwise.   

    MM      (input) INTEGER   
            The number of columns in the arrays VL and/or VR. MM >= M.   

    M       (output) INTEGER   
            The number of columns in the arrays VL and/or VR actually   
            used to store the eigenvectors.  If HOWMNY = 'A' or 'B', M   
            is set to N.  Each selected eigenvector occupies one column. 
  

    WORK    (workspace) COMPLEX array, dimension (2*N)   

    RWORK   (workspace) REAL array, dimension (2*N)   

    INFO    (output) INTEGER   
            = 0:  successful exit.   
            < 0:  if INFO = -i, the i-th argument had an illegal value.   

    ===================================================================== 
  


       Decode and Test the input parameters   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static complex c_b1 = {0.f,0.f};
    static complex c_b2 = {1.f,0.f};
    static integer c__1 = 1;
    
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, vl_dim1, vl_offset, vr_dim1, 
	    vr_offset, i__1, i__2, i__3, i__4, i__5;
    real r__1, r__2, r__3, r__4, r__5, r__6;
    complex q__1, q__2, q__3, q__4;
    /* Builtin functions */
    double r_imag(complex *);
    void r_cnjg(complex *, complex *), c_div(complex *, complex *, complex *);
    /* Local variables */
    static integer ibeg, ieig, iend;
    static real dmin__;
    static integer isrc;
    static real temp;
    static complex suma, sumb;
    static real xmax;
    static complex d;
    static integer i, j;
    static real scale;
    static logical ilall;
    static integer iside;
    static real sbeta;
    extern logical lsame_(char *, char *);
    extern /* Subroutine */ int cgemv_(char *, integer *, integer *, complex *
	    , complex *, integer *, complex *, integer *, complex *, complex *
	    , integer *);
    static real small;
    static logical compl;
    static real anorm, bnorm;
    static logical compr;
    static complex ca, cb;
    static logical ilbbad;
    static real acoefa;
    static integer je;
    static real bcoefa, acoeff;
    static complex bcoeff;
    static logical ilback;
    static integer im;
    extern /* Subroutine */ int slabad_(real *, real *);
    static real ascale, bscale;
    static integer jr;
    extern doublereal slamch_(char *);
    static complex salpha;
    static real safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *);
    static real bignum;
    static logical ilcomp;
    static integer ihwmny;
    static real big;
    static logical lsa, lsb;
    static real ulp;
    static complex sum;



#define SELECT(I) select[(I)-1]
#define WORK(I) work[(I)-1]
#define RWORK(I) rwork[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]
#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]
#define VL(I,J) vl[(I)-1 + ((J)-1)* ( *ldvl)]
#define VR(I,J) vr[(I)-1 + ((J)-1)* ( *ldvr)]

    if (lsame_(howmny, "A")) {
	ihwmny = 1;
	ilall = TRUE_;
	ilback = FALSE_;
    } else if (lsame_(howmny, "S")) {
	ihwmny = 2;
	ilall = FALSE_;
	ilback = FALSE_;
    } else if (lsame_(howmny, "B") || lsame_(howmny, "T")) {
	ihwmny = 3;
	ilall = TRUE_;
	ilback = TRUE_;
    } else {
	ihwmny = -1;
    }

    if (lsame_(side, "R")) {
	iside = 1;
	compl = FALSE_;
	compr = TRUE_;
    } else if (lsame_(side, "L")) {
	iside = 2;
	compl = TRUE_;
	compr = FALSE_;
    } else if (lsame_(side, "B")) {
	iside = 3;
	compl = TRUE_;
	compr = TRUE_;
    } else {
	iside = -1;
    }

/*     Count the number of eigenvectors */

    if (! ilall) {
	im = 0;
	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
	    if (SELECT(j)) {
		++im;
	    }
/* L10: */
	}
    } else {
	im = *n;
    }

/*     Check diagonal of B */

    ilbbad = FALSE_;
    i__1 = *n;
    for (j = 1; j <= *n; ++j) {
	if (r_imag(&B(j,j)) != 0.f) {
	    ilbbad = TRUE_;
	}
/* L20: */
    }

    *info = 0;
    if (iside < 0) {
	*info = -1;
    } else if (ihwmny < 0) {
	*info = -2;
    } else if (*n < 0) {
	*info = -4;
    } else if (*lda < max(1,*n)) {
	*info = -6;
    } else if (ilbbad) {
	*info = -7;
    } else if (*ldb < max(1,*n)) {
	*info = -8;
    } else if (compl && *ldvl < *n || *ldvl < 1) {
	*info = -10;
    } else if (compr && *ldvr < *n || *ldvr < 1) {
	*info = -12;
    } else if (*mm < im) {
	*info = -13;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("CTGEVC", &i__1);
	return 0;
    }

/*     Quick return if possible */

    *m = im;
    if (*n == 0) {
	return 0;
    }

/*     Machine Constants */

    safmin = slamch_("Safe minimum");
    big = 1.f / safmin;
    slabad_(&safmin, &big);
    ulp = slamch_("Epsilon") * slamch_("Base");
    small = safmin * *n / ulp;
    big = 1.f / small;
    bignum = 1.f / (safmin * *n);

/*     Compute the 1-norm of each column of the strictly upper triangular 
  
       part of A and B to check for possible overflow in the triangular   
       solver. */

    i__1 = a_dim1 + 1;
    anorm = (r__1 = A(1,1).r, dabs(r__1)) + (r__2 = r_imag(&A(1,1)), 
	    dabs(r__2));
    i__1 = b_dim1 + 1;
    bnorm = (r__1 = B(1,1).r, dabs(r__1)) + (r__2 = r_imag(&B(1,1)), 
	    dabs(r__2));
    RWORK(1) = 0.f;
    RWORK(*n + 1) = 0.f;
    i__1 = *n;
    for (j = 2; j <= *n; ++j) {
	RWORK(j) = 0.f;
	RWORK(*n + j) = 0.f;
	i__2 = j - 1;
	for (i = 1; i <= j-1; ++i) {
	    i__3 = i + j * a_dim1;
	    RWORK(j) += (r__1 = A(i,j).r, dabs(r__1)) + (r__2 = r_imag(&A(i,j)), dabs(r__2));
	    i__3 = i + j * b_dim1;
	    RWORK(*n + j) += (r__1 = B(i,j).r, dabs(r__1)) + (r__2 = r_imag(&
		    B(i,j)), dabs(r__2));
/* L30: */
	}
/* Computing MAX */
	i__2 = j + j * a_dim1;
	r__3 = anorm, r__4 = RWORK(j) + ((r__1 = A(j,j).r, dabs(r__1)) + (
		r__2 = r_imag(&A(j,j)), dabs(r__2)));
	anorm = dmax(r__3,r__4);
/* Computing MAX */
	i__2 = j + j * b_dim1;
	r__3 = bnorm, r__4 = RWORK(*n + j) + ((r__1 = B(j,j).r, dabs(r__1)) 
		+ (r__2 = r_imag(&B(j,j)), dabs(r__2)));
	bnorm = dmax(r__3,r__4);
/* L40: */
    }

    ascale = 1.f / dmax(anorm,safmin);
    bscale = 1.f / dmax(bnorm,safmin);

/*     Left eigenvectors */

    if (compl) {
	ieig = 0;

/*        Main loop over eigenvalues */

	i__1 = *n;
	for (je = 1; je <= *n; ++je) {
	    if (ilall) {
		ilcomp = TRUE_;
	    } else {
		ilcomp = SELECT(je);
	    }
	    if (ilcomp) {
		++ieig;

		i__2 = je + je * a_dim1;
		i__3 = je + je * b_dim1;
		if ((r__1 = A(je,je).r, dabs(r__1)) + (r__2 = r_imag(&A(je,je)), dabs(r__2)) <= safmin && (r__3 = B(je,je).r, dabs(r__3)) <= safmin) {

/*                 Singular matrix pencil -- return unit e
igenvector */

		    i__2 = *n;
		    for (jr = 1; jr <= *n; ++jr) {
			i__3 = jr + ieig * vl_dim1;
			VL(jr,ieig).r = 0.f, VL(jr,ieig).i = 0.f;
/* L50: */
		    }
		    i__2 = ieig + ieig * vl_dim1;
		    VL(ieig,ieig).r = 1.f, VL(ieig,ieig).i = 0.f;
		    goto L140;
		}

/*              Non-singular eigenvalue:   
                Compute coefficients  a  and  b  in   
                     H   
                   y  ( a A - b B ) = 0   

   Computing MAX */
		i__2 = je + je * a_dim1;
		i__3 = je + je * b_dim1;
		r__4 = ((r__1 = A(je,je).r, dabs(r__1)) + (r__2 = r_imag(&A(je,je)), dabs(r__2))) * ascale, r__5 = (r__3 =
			 B(je,je).r, dabs(r__3)) * bscale, r__4 = max(r__4,
			r__5);
		temp = 1.f / dmax(r__4,safmin);
		i__2 = je + je * a_dim1;
		q__2.r = temp * A(je,je).r, q__2.i = temp * A(je,je).i;
		q__1.r = ascale * q__2.r, q__1.i = ascale * q__2.i;
		salpha.r = q__1.r, salpha.i = q__1.i;
		i__2 = je + je * b_dim1;
		sbeta = temp * B(je,je).r * bscale;
		acoeff = sbeta * ascale;
		q__1.r = bscale * salpha.r, q__1.i = bscale * salpha.i;
		bcoeff.r = q__1.r, bcoeff.i = q__1.i;

/*              Scale to avoid underflow */

		lsa = dabs(sbeta) >= safmin && dabs(acoeff) < small;
		lsb = (r__1 = salpha.r, dabs(r__1)) + (r__2 = r_imag(&salpha),
			 dabs(r__2)) >= safmin && (r__3 = bcoeff.r, dabs(r__3)
			) + (r__4 = r_imag(&bcoeff), dabs(r__4)) < small;

		scale = 1.f;
		if (lsa) {
		    scale = small / dabs(sbeta) * dmin(anorm,big);
		}
		if (lsb) {
/* Computing MAX */
		    r__3 = scale, r__4 = small / ((r__1 = salpha.r, dabs(r__1)
			    ) + (r__2 = r_imag(&salpha), dabs(r__2))) * dmin(
			    bnorm,big);
		    scale = dmax(r__3,r__4);
		}
		if (lsa || lsb) {
/* Computing MIN   
   Computing MAX */
		    r__5 = 1.f, r__6 = dabs(acoeff), r__5 = max(r__5,r__6), 
			    r__6 = (r__1 = bcoeff.r, dabs(r__1)) + (r__2 = 
			    r_imag(&bcoeff), dabs(r__2));
		    r__3 = scale, r__4 = 1.f / (safmin * dmax(r__5,r__6));
		    scale = dmin(r__3,r__4);
		    if (lsa) {
			acoeff = ascale * (scale * sbeta);
		    } else {
			acoeff = scale * acoeff;
		    }
		    if (lsb) {
			q__2.r = scale * salpha.r, q__2.i = scale * salpha.i;
			q__1.r = bscale * q__2.r, q__1.i = bscale * q__2.i;
			bcoeff.r = q__1.r, bcoeff.i = q__1.i;
		    } else {
			q__1.r = scale * bcoeff.r, q__1.i = scale * bcoeff.i;
			bcoeff.r = q__1.r, bcoeff.i = q__1.i;
		    }
		}

		acoefa = dabs(acoeff);
		bcoefa = (r__1 = bcoeff.r, dabs(r__1)) + (r__2 = r_imag(&
			bcoeff), dabs(r__2));
		xmax = 1.f;
		i__2 = *n;
		for (jr = 1; jr <= *n; ++jr) {
		    i__3 = jr;
		    WORK(jr).r = 0.f, WORK(jr).i = 0.f;
/* L60: */
		}
		i__2 = je;
		WORK(je).r = 1.f, WORK(je).i = 0.f;
/* Computing MAX */
		r__1 = ulp * acoefa * anorm, r__2 = ulp * bcoefa * bnorm, 
			r__1 = max(r__1,r__2);
		dmin__ = dmax(r__1,safmin);

/*                                              H   
                Triangular solve of  (a A - b B)  y = 0   

                                        H   
                (rowwise in  (a A - b B) , or columnwise in a 
A - b B) */

		i__2 = *n;
		for (j = je + 1; j <= *n; ++j) {

/*                 Compute   
                         j-1   
                   SUM = sum  conjg( a*A(k,j) - b*B(k,j) )
*x(k)   
                         k=je   
                   (Scale if necessary) */

		    temp = 1.f / xmax;
		    if (acoefa * RWORK(j) + bcoefa * RWORK(*n + j) > bignum * 
			    temp) {
			i__3 = j - 1;
			for (jr = je; jr <= j-1; ++jr) {
			    i__4 = jr;
			    i__5 = jr;
			    q__1.r = temp * WORK(jr).r, q__1.i = temp * 
				    WORK(jr).i;
			    WORK(jr).r = q__1.r, WORK(jr).i = q__1.i;
/* L70: */
			}
			xmax = 1.f;
		    }
		    suma.r = 0.f, suma.i = 0.f;
		    sumb.r = 0.f, sumb.i = 0.f;

		    i__3 = j - 1;
		    for (jr = je; jr <= j-1; ++jr) {
			r_cnjg(&q__3, &A(jr,j));
			i__4 = jr;
			q__2.r = q__3.r * WORK(jr).r - q__3.i * WORK(jr)
				.i, q__2.i = q__3.r * WORK(jr).i + q__3.i * 
				WORK(jr).r;
			q__1.r = suma.r + q__2.r, q__1.i = suma.i + q__2.i;
			suma.r = q__1.r, suma.i = q__1.i;
			r_cnjg(&q__3, &B(jr,j));
			i__4 = jr;
			q__2.r = q__3.r * WORK(jr).r - q__3.i * WORK(jr)
				.i, q__2.i = q__3.r * WORK(jr).i + q__3.i * 
				WORK(jr).r;
			q__1.r = sumb.r + q__2.r, q__1.i = sumb.i + q__2.i;
			sumb.r = q__1.r, sumb.i = q__1.i;
/* L80: */
		    }
		    q__2.r = acoeff * suma.r, q__2.i = acoeff * suma.i;
		    r_cnjg(&q__4, &bcoeff);
		    q__3.r = q__4.r * sumb.r - q__4.i * sumb.i, q__3.i = 
			    q__4.r * sumb.i + q__4.i * sumb.r;
		    q__1.r = q__2.r - q__3.r, q__1.i = q__2.i - q__3.i;
		    sum.r = q__1.r, sum.i = q__1.i;

/*                 Form x(j) = - SUM / conjg( a*A(j,j) - b
*B(j,j) )   

                   with scaling and perturbation of the de
nominator */

		    i__3 = j + j * a_dim1;
		    q__3.r = acoeff * A(j,j).r, q__3.i = acoeff * A(j,j).i;
		    i__4 = j + j * b_dim1;
		    q__4.r = bcoeff.r * B(j,j).r - bcoeff.i * B(j,j).i, 
			    q__4.i = bcoeff.r * B(j,j).i + bcoeff.i * B(j,j)
			    .r;
		    q__2.r = q__3.r - q__4.r, q__2.i = q__3.i - q__4.i;
		    r_cnjg(&q__1, &q__2);
		    d.r = q__1.r, d.i = q__1.i;
		    if ((r__1 = d.r, dabs(r__1)) + (r__2 = r_imag(&d), dabs(
			    r__2)) <= dmin__) {
			q__1.r = dmin__, q__1.i = 0.f;
			d.r = q__1.r, d.i = q__1.i;
		    }

		    if ((r__1 = d.r, dabs(r__1)) + (r__2 = r_imag(&d), dabs(
			    r__2)) < 1.f) {
			if ((r__1 = sum.r, dabs(r__1)) + (r__2 = r_imag(&sum),
				 dabs(r__2)) >= bignum * ((r__3 = d.r, dabs(
				r__3)) + (r__4 = r_imag(&d), dabs(r__4)))) {
			    temp = 1.f / ((r__1 = sum.r, dabs(r__1)) + (r__2 =
				     r_imag(&sum), dabs(r__2)));
			    i__3 = j - 1;
			    for (jr = je; jr <= j-1; ++jr) {
				i__4 = jr;
				i__5 = jr;
				q__1.r = temp * WORK(jr).r, q__1.i = temp * 
					WORK(jr).i;
				WORK(jr).r = q__1.r, WORK(jr).i = q__1.i;
/* L90: */
			    }
			    xmax = temp * xmax;
			    q__1.r = temp * sum.r, q__1.i = temp * sum.i;
			    sum.r = q__1.r, sum.i = q__1.i;
			}
		    }
		    i__3 = j;
		    q__2.r = -(doublereal)sum.r, q__2.i = -(doublereal)sum.i;
		    c_div(&q__1, &q__2, &d);
		    WORK(j).r = q__1.r, WORK(j).i = q__1.i;
/* Computing MAX */
		    i__3 = j;
		    r__3 = xmax, r__4 = (r__1 = WORK(j).r, dabs(r__1)) + (
			    r__2 = r_imag(&WORK(j)), dabs(r__2));
		    xmax = dmax(r__3,r__4);
/* L100: */
		}

/*              Back transform eigenvector if HOWMNY='B'. */

		if (ilback) {
		    i__2 = *n + 1 - je;
		    cgemv_("N", n, &i__2, &c_b2, &VL(1,je), ldvl, 
			    &WORK(je), &c__1, &c_b1, &WORK(*n + 1), &c__1)
			    ;
		    isrc = 2;
		    ibeg = 1;
		} else {
		    isrc = 1;
		    ibeg = je;
		}

/*              Copy and scale eigenvector into column of VL 
*/

		xmax = 0.f;
		i__2 = *n;
		for (jr = ibeg; jr <= *n; ++jr) {
/* Computing MAX */
		    i__3 = (isrc - 1) * *n + jr;
		    r__3 = xmax, r__4 = (r__1 = WORK((isrc-1)**n+jr).r, dabs(r__1)) + (
			    r__2 = r_imag(&WORK((isrc - 1) * *n + jr)), dabs(
			    r__2));
		    xmax = dmax(r__3,r__4);
/* L110: */
		}

		if (xmax > safmin) {
		    temp = 1.f / xmax;
		    i__2 = *n;
		    for (jr = ibeg; jr <= *n; ++jr) {
			i__3 = jr + ieig * vl_dim1;
			i__4 = (isrc - 1) * *n + jr;
			q__1.r = temp * WORK((isrc-1)**n+jr).r, q__1.i = temp * WORK(
				(isrc-1)**n+jr).i;
			VL(jr,ieig).r = q__1.r, VL(jr,ieig).i = q__1.i;
/* L120: */
		    }
		} else {
		    ibeg = *n + 1;
		}

		i__2 = ibeg - 1;
		for (jr = 1; jr <= ibeg-1; ++jr) {
		    i__3 = jr + ieig * vl_dim1;
		    VL(jr,ieig).r = 0.f, VL(jr,ieig).i = 0.f;
/* L130: */
		}

	    }
L140:
	    ;
	}
    }

/*     Right eigenvectors */

    if (compr) {
	ieig = im + 1;

/*        Main loop over eigenvalues */

	for (je = *n; je >= 1; --je) {
	    if (ilall) {
		ilcomp = TRUE_;
	    } else {
		ilcomp = SELECT(je);
	    }
	    if (ilcomp) {
		--ieig;

		i__1 = je + je * a_dim1;
		i__2 = je + je * b_dim1;
		if ((r__1 = A(je,je).r, dabs(r__1)) + (r__2 = r_imag(&A(je,je)), dabs(r__2)) <= safmin && (r__3 = B(je,je).r, dabs(r__3)) <= safmin) {

/*                 Singular matrix pencil -- return unit e
igenvector */

		    i__1 = *n;
		    for (jr = 1; jr <= *n; ++jr) {
			i__2 = jr + ieig * vr_dim1;
			VR(jr,ieig).r = 0.f, VR(jr,ieig).i = 0.f;
/* L150: */
		    }
		    i__1 = ieig + ieig * vr_dim1;
		    VR(ieig,ieig).r = 1.f, VR(ieig,ieig).i = 0.f;
		    goto L250;
		}

/*              Non-singular eigenvalue:   
                Compute coefficients  a  and  b  in   

                ( a A - b B ) x  = 0   

   Computing MAX */
		i__1 = je + je * a_dim1;
		i__2 = je + je * b_dim1;
		r__4 = ((r__1 = A(je,je).r, dabs(r__1)) + (r__2 = r_imag(&A(je,je)), dabs(r__2))) * ascale, r__5 = (r__3 =
			 B(je,je).r, dabs(r__3)) * bscale, r__4 = max(r__4,
			r__5);
		temp = 1.f / dmax(r__4,safmin);
		i__1 = je + je * a_dim1;
		q__2.r = temp * A(je,je).r, q__2.i = temp * A(je,je).i;
		q__1.r = ascale * q__2.r, q__1.i = ascale * q__2.i;
		salpha.r = q__1.r, salpha.i = q__1.i;
		i__1 = je + je * b_dim1;
		sbeta = temp * B(je,je).r * bscale;
		acoeff = sbeta * ascale;
		q__1.r = bscale * salpha.r, q__1.i = bscale * salpha.i;
		bcoeff.r = q__1.r, bcoeff.i = q__1.i;

/*              Scale to avoid underflow */

		lsa = dabs(sbeta) >= safmin && dabs(acoeff) < small;
		lsb = (r__1 = salpha.r, dabs(r__1)) + (r__2 = r_imag(&salpha),
			 dabs(r__2)) >= safmin && (r__3 = bcoeff.r, dabs(r__3)
			) + (r__4 = r_imag(&bcoeff), dabs(r__4)) < small;

		scale = 1.f;
		if (lsa) {
		    scale = small / dabs(sbeta) * dmin(anorm,big);
		}
		if (lsb) {
/* Computing MAX */
		    r__3 = scale, r__4 = small / ((r__1 = salpha.r, dabs(r__1)
			    ) + (r__2 = r_imag(&salpha), dabs(r__2))) * dmin(
			    bnorm,big);
		    scale = dmax(r__3,r__4);
		}
		if (lsa || lsb) {
/* Computing MIN   
   Computing MAX */
		    r__5 = 1.f, r__6 = dabs(acoeff), r__5 = max(r__5,r__6), 
			    r__6 = (r__1 = bcoeff.r, dabs(r__1)) + (r__2 = 
			    r_imag(&bcoeff), dabs(r__2));
		    r__3 = scale, r__4 = 1.f / (safmin * dmax(r__5,r__6));
		    scale = dmin(r__3,r__4);
		    if (lsa) {
			acoeff = ascale * (scale * sbeta);
		    } else {
			acoeff = scale * acoeff;
		    }
		    if (lsb) {
			q__2.r = scale * salpha.r, q__2.i = scale * salpha.i;
			q__1.r = bscale * q__2.r, q__1.i = bscale * q__2.i;
			bcoeff.r = q__1.r, bcoeff.i = q__1.i;
		    } else {
			q__1.r = scale * bcoeff.r, q__1.i = scale * bcoeff.i;
			bcoeff.r = q__1.r, bcoeff.i = q__1.i;
		    }
		}

		acoefa = dabs(acoeff);
		bcoefa = (r__1 = bcoeff.r, dabs(r__1)) + (r__2 = r_imag(&
			bcoeff), dabs(r__2));
		xmax = 1.f;
		i__1 = *n;
		for (jr = 1; jr <= *n; ++jr) {
		    i__2 = jr;
		    WORK(jr).r = 0.f, WORK(jr).i = 0.f;
/* L160: */
		}
		i__1 = je;
		WORK(je).r = 1.f, WORK(je).i = 0.f;
/* Computing MAX */
		r__1 = ulp * acoefa * anorm, r__2 = ulp * bcoefa * bnorm, 
			r__1 = max(r__1,r__2);
		dmin__ = dmax(r__1,safmin);

/*              Triangular solve of  (a A - b B) x = 0  (colum
nwise)   

                WORK(1:j-1) contains sums w,   
                WORK(j+1:JE) contains x */

		i__1 = je - 1;
		for (jr = 1; jr <= je-1; ++jr) {
		    i__2 = jr;
		    i__3 = jr + je * a_dim1;
		    q__2.r = acoeff * A(jr,je).r, q__2.i = acoeff * A(jr,je).i;
		    i__4 = jr + je * b_dim1;
		    q__3.r = bcoeff.r * B(jr,je).r - bcoeff.i * B(jr,je).i, 
			    q__3.i = bcoeff.r * B(jr,je).i + bcoeff.i * B(jr,je)
			    .r;
		    q__1.r = q__2.r - q__3.r, q__1.i = q__2.i - q__3.i;
		    WORK(jr).r = q__1.r, WORK(jr).i = q__1.i;
/* L170: */
		}
		i__1 = je;
		WORK(je).r = 1.f, WORK(je).i = 0.f;

		for (j = je - 1; j >= 1; --j) {

/*                 Form x(j) := - w(j) / d   
                   with scaling and perturbation of the de
nominator */

		    i__1 = j + j * a_dim1;
		    q__2.r = acoeff * A(j,j).r, q__2.i = acoeff * A(j,j).i;
		    i__2 = j + j * b_dim1;
		    q__3.r = bcoeff.r * B(j,j).r - bcoeff.i * B(j,j).i, 
			    q__3.i = bcoeff.r * B(j,j).i + bcoeff.i * B(j,j)
			    .r;
		    q__1.r = q__2.r - q__3.r, q__1.i = q__2.i - q__3.i;
		    d.r = q__1.r, d.i = q__1.i;
		    if ((r__1 = d.r, dabs(r__1)) + (r__2 = r_imag(&d), dabs(
			    r__2)) <= dmin__) {
			q__1.r = dmin__, q__1.i = 0.f;
			d.r = q__1.r, d.i = q__1.i;
		    }

		    if ((r__1 = d.r, dabs(r__1)) + (r__2 = r_imag(&d), dabs(
			    r__2)) < 1.f) {
			i__1 = j;
			if ((r__1 = WORK(j).r, dabs(r__1)) + (r__2 = 
				r_imag(&WORK(j)), dabs(r__2)) >= bignum * ((
				r__3 = d.r, dabs(r__3)) + (r__4 = r_imag(&d), 
				dabs(r__4)))) {
			    i__1 = j;
			    temp = 1.f / ((r__1 = WORK(j).r, dabs(r__1)) + 
				    (r__2 = r_imag(&WORK(j)), dabs(r__2)));
			    i__1 = je;
			    for (jr = 1; jr <= je; ++jr) {
				i__2 = jr;
				i__3 = jr;
				q__1.r = temp * WORK(jr).r, q__1.i = temp * 
					WORK(jr).i;
				WORK(jr).r = q__1.r, WORK(jr).i = q__1.i;
/* L180: */
			    }
			}
		    }

		    i__1 = j;
		    i__2 = j;
		    q__2.r = -(doublereal)WORK(j).r, q__2.i = -(doublereal)
			    WORK(j).i;
		    c_div(&q__1, &q__2, &d);
		    WORK(j).r = q__1.r, WORK(j).i = q__1.i;

		    if (j > 1) {

/*                    w = w + x(j)*(a A(*,j) - b B(*,j
) ) with scaling */

			i__1 = j;
			if ((r__1 = WORK(j).r, dabs(r__1)) + (r__2 = 
				r_imag(&WORK(j)), dabs(r__2)) > 1.f) {
			    i__1 = j;
			    temp = 1.f / ((r__1 = WORK(j).r, dabs(r__1)) + 
				    (r__2 = r_imag(&WORK(j)), dabs(r__2)));
			    if (acoefa * RWORK(j) + bcoefa * RWORK(*n + j) >= 
				    bignum * temp) {
				i__1 = je;
				for (jr = 1; jr <= je; ++jr) {
				    i__2 = jr;
				    i__3 = jr;
				    q__1.r = temp * WORK(jr).r, q__1.i = 
					    temp * WORK(jr).i;
				    WORK(jr).r = q__1.r, WORK(jr).i = 
					    q__1.i;
/* L190: */
				}
			    }
			}

			i__1 = j;
			q__1.r = acoeff * WORK(j).r, q__1.i = acoeff * 
				WORK(j).i;
			ca.r = q__1.r, ca.i = q__1.i;
			i__1 = j;
			q__1.r = bcoeff.r * WORK(j).r - bcoeff.i * WORK(
				j).i, q__1.i = bcoeff.r * WORK(j).i + 
				bcoeff.i * WORK(j).r;
			cb.r = q__1.r, cb.i = q__1.i;
			i__1 = j - 1;
			for (jr = 1; jr <= j-1; ++jr) {
			    i__2 = jr;
			    i__3 = jr;
			    i__4 = jr + j * a_dim1;
			    q__3.r = ca.r * A(jr,j).r - ca.i * A(jr,j).i, 
				    q__3.i = ca.r * A(jr,j).i + ca.i * A(jr,j)
				    .r;
			    q__2.r = WORK(jr).r + q__3.r, q__2.i = WORK(
				    jr).i + q__3.i;
			    i__5 = jr + j * b_dim1;
			    q__4.r = cb.r * B(jr,j).r - cb.i * B(jr,j).i, 
				    q__4.i = cb.r * B(jr,j).i + cb.i * B(jr,j)
				    .r;
			    q__1.r = q__2.r - q__4.r, q__1.i = q__2.i - 
				    q__4.i;
			    WORK(jr).r = q__1.r, WORK(jr).i = q__1.i;
/* L200: */
			}
		    }
/* L210: */
		}

/*              Back transform eigenvector if HOWMNY='B'. */

		if (ilback) {
		    cgemv_("N", n, &je, &c_b2, &VR(1,1), ldvr, &WORK(1),
			     &c__1, &c_b1, &WORK(*n + 1), &c__1);
		    isrc = 2;
		    iend = *n;
		} else {
		    isrc = 1;
		    iend = je;
		}

/*              Copy and scale eigenvector into column of VR 
*/

		xmax = 0.f;
		i__1 = iend;
		for (jr = 1; jr <= iend; ++jr) {
/* Computing MAX */
		    i__2 = (isrc - 1) * *n + jr;
		    r__3 = xmax, r__4 = (r__1 = WORK((isrc-1)**n+jr).r, dabs(r__1)) + (
			    r__2 = r_imag(&WORK((isrc - 1) * *n + jr)), dabs(
			    r__2));
		    xmax = dmax(r__3,r__4);
/* L220: */
		}

		if (xmax > safmin) {
		    temp = 1.f / xmax;
		    i__1 = iend;
		    for (jr = 1; jr <= iend; ++jr) {
			i__2 = jr + ieig * vr_dim1;
			i__3 = (isrc - 1) * *n + jr;
			q__1.r = temp * WORK((isrc-1)**n+jr).r, q__1.i = temp * WORK(
				(isrc-1)**n+jr).i;
			VR(jr,ieig).r = q__1.r, VR(jr,ieig).i = q__1.i;
/* L230: */
		    }
		} else {
		    iend = 0;
		}

		i__1 = *n;
		for (jr = iend + 1; jr <= *n; ++jr) {
		    i__2 = jr + ieig * vr_dim1;
		    VR(jr,ieig).r = 0.f, VR(jr,ieig).i = 0.f;
/* L240: */
		}

	    }
L250:
	    ;
	}
    }

    return 0;

/*     End of CTGEVC */

} /* ctgevc_ */

