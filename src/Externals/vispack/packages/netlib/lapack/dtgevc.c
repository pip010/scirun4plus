#include "f2c.h"

/* Subroutine */ int dtgevc_(char *side, char *howmny, logical *select, 
	integer *n, doublereal *a, integer *lda, doublereal *b, integer *ldb, 
	doublereal *vl, integer *ldvl, doublereal *vr, integer *ldvr, integer 
	*mm, integer *m, doublereal *work, integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   



    Purpose   
    =======   

    DTGEVC computes some or all of the right and/or left generalized   
    eigenvectors of a pair of real upper triangular matrices (A,B).   

    The right generalized eigenvector x and the left generalized   
    eigenvector y of (A,B) corresponding to a generalized eigenvalue   
    w are defined by:   

            (A - wB) * x = 0  and  y**H * (A - wB) = 0   

    where y**H denotes the conjugate tranpose of y.   

    If an eigenvalue w is determined by zero diagonal elements of both A 
  
    and B, a unit vector is returned as the corresponding eigenvector.   

    If all eigenvectors are requested, the routine may either return   
    the matrices X and/or Y of right or left eigenvectors of (A,B), or   
    the products Z*X and/or Q*Y, where Z and Q are input orthogonal   
    matrices.  If (A,B) was obtained from the generalized real-Schur   
    factorization of an original pair of matrices   
       (A0,B0) = (Q*A*Z**H,Q*B*Z**H),   
    then Z*X and Q*Y are the matrices of right or left eigenvectors of   
    A.   

    A must be block upper triangular, with 1-by-1 and 2-by-2 diagonal   
    blocks.  Corresponding to each 2-by-2 diagonal block is a complex   
    conjugate pair of eigenvalues and eigenvectors; only one   
    eigenvector of the pair is computed, namely the one corresponding   
    to the eigenvalue with positive imaginary part.   

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
            To select the real eigenvector corresponding to the real   
            eigenvalue w(j), SELECT(j) must be set to .TRUE.  To select   
            the complex eigenvector corresponding to a complex conjugate 
  
            pair w(j) and w(j+1), either SELECT(j) or SELECT(j+1) must   
            be set to .TRUE..   

    N       (input) INTEGER   
            The order of the matrices A and B.  N >= 0.   

    A       (input) DOUBLE PRECISION array, dimension (LDA,N)   
            The upper quasi-triangular matrix A.   

    LDA     (input) INTEGER   
            The leading dimension of array A.  LDA >= max(1, N).   

    B       (input) DOUBLE PRECISION array, dimension (LDB,N)   
            The upper triangular matrix B.  If A has a 2-by-2 diagonal   
            block, then the corresponding 2-by-2 block of B must be   
            diagonal with positive elements.   

    LDB     (input) INTEGER   
            The leading dimension of array B.  LDB >= max(1,N).   

    VL      (input/output) DOUBLE PRECISION array, dimension (LDVL,MM)   
            On entry, if SIDE = 'L' or 'B' and HOWMNY = 'B', VL must   
            contain an N-by-N matrix Q (usually the orthogonal matrix Q   
            of left Schur vectors returned by DHGEQZ).   
            On exit, if SIDE = 'L' or 'B', VL contains:   
            if HOWMNY = 'A', the matrix Y of left eigenvectors of (A,B); 
  
            if HOWMNY = 'B', the matrix Q*Y;   
            if HOWMNY = 'S', the left eigenvectors of (A,B) specified by 
  
                        SELECT, stored consecutively in the columns of   
                        VL, in the same order as their eigenvalues.   
            If SIDE = 'R', VL is not referenced.   

            A complex eigenvector corresponding to a complex eigenvalue   
            is stored in two consecutive columns, the first holding the   
            real part, and the second the imaginary part.   

    LDVL    (input) INTEGER   
            The leading dimension of array VL.   
            LDVL >= max(1,N) if SIDE = 'L' or 'B'; LDVL >= 1 otherwise.   

    VR      (input/output) COMPLEX*16 array, dimension (LDVR,MM)   
            On entry, if SIDE = 'R' or 'B' and HOWMNY = 'B', VR must   
            contain an N-by-N matrix Q (usually the orthogonal matrix Z   
            of right Schur vectors returned by DHGEQZ).   
            On exit, if SIDE = 'R' or 'B', VR contains:   
            if HOWMNY = 'A', the matrix X of right eigenvectors of (A,B); 
  
            if HOWMNY = 'B', the matrix Z*X;   
            if HOWMNY = 'S', the right eigenvectors of (A,B) specified by 
  
                        SELECT, stored consecutively in the columns of   
                        VR, in the same order as their eigenvalues.   
            If SIDE = 'L', VR is not referenced.   

            A complex eigenvector corresponding to a complex eigenvalue   
            is stored in two consecutive columns, the first holding the   
            real part and the second the imaginary part.   

    LDVR    (input) INTEGER   
            The leading dimension of the array VR.   
            LDVR >= max(1,N) if SIDE = 'R' or 'B'; LDVR >= 1 otherwise.   

    MM      (input) INTEGER   
            The number of columns in the arrays VL and/or VR. MM >= M.   

    M       (output) INTEGER   
            The number of columns in the arrays VL and/or VR actually   
            used to store the eigenvectors.  If HOWMNY = 'A' or 'B', M   
            is set to N.  Each selected real eigenvector occupies one   
            column and each selected complex eigenvector occupies two   
            columns.   

    WORK    (workspace) DOUBLE PRECISION array, dimension (6*N)   

    INFO    (output) INTEGER   
            = 0:  successful exit.   
            < 0:  if INFO = -i, the i-th argument had an illegal value.   
            > 0:  the 2-by-2 block (INFO:INFO+1) does not have a complex 
  
                  eigenvalue.   

    Further Details   
    ===============   

    Allocation of workspace:   
    ---------- -- ---------   

       WORK( j ) = 1-norm of j-th column of A, above the diagonal   
       WORK( N+j ) = 1-norm of j-th column of B, above the diagonal   
       WORK( 2*N+1:3*N ) = real part of eigenvector   
       WORK( 3*N+1:4*N ) = imaginary part of eigenvector   
       WORK( 4*N+1:5*N ) = real part of back-transformed eigenvector   
       WORK( 5*N+1:6*N ) = imaginary part of back-transformed eigenvector 
  

    Rowwise vs. columnwise solution methods:   
    ------- --  ---------- -------- -------   

    Finding a generalized eigenvector consists basically of solving the   
    singular triangular system   

     (A - w B) x = 0     (for right) or:   (A - w B)**H y = 0  (for left) 
  

    Consider finding the i-th right eigenvector (assume all eigenvalues   
    are real). The equation to be solved is:   
         n                   i   
    0 = sum  C(j,k) v(k)  = sum  C(j,k) v(k)     for j = i,. . .,1   
        k=j                 k=j   

    where  C = (A - w B)  (The components v(i+1:n) are 0.)   

    The "rowwise" method is:   

    (1)  v(i) := 1   
    for j = i-1,. . .,1:   
                            i   
        (2) compute  s = - sum C(j,k) v(k)   and   
                          k=j+1   

        (3) v(j) := s / C(j,j)   

    Step 2 is sometimes called the "dot product" step, since it is an   
    inner product between the j-th row and the portion of the eigenvector 
  
    that has been computed so far.   

    The "columnwise" method consists basically in doing the sums   
    for all the rows in parallel.  As each v(j) is computed, the   
    contribution of v(j) times the j-th column of C is added to the   
    partial sums.  Since FORTRAN arrays are stored columnwise, this has   
    the advantage that at each step, the elements of C that are accessed 
  
    are adjacent to one another, whereas with the rowwise method, the   
    elements accessed at a step are spaced LDA (and LDB) words apart.   

    When finding left eigenvectors, the matrix in question is the   
    transpose of the one in storage, so the rowwise method then   
    actually accesses columns of A and B at each step, and so is the   
    preferred method.   

    ===================================================================== 
  


       Decode and Test the input parameters   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static logical c_true = TRUE_;
    static integer c__2 = 2;
    static doublereal c_b34 = 1.;
    static integer c__1 = 1;
    static doublereal c_b36 = 0.;
    static logical c_false = FALSE_;
    
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, vl_dim1, vl_offset, vr_dim1, 
	    vr_offset, i__1, i__2, i__3, i__4, i__5;
    doublereal d__1, d__2, d__3, d__4, d__5, d__6;
    /* Local variables */
    static integer ibeg, ieig, iend;
    static doublereal dmin__, temp, suma[4]	/* was [2][2] */, sumb[4]	
	    /* was [2][2] */, xmax;
    extern /* Subroutine */ int dlag2_(doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *);
    static doublereal cim2a, cim2b, cre2a, cre2b, temp2, bdiag[2];
    static integer i, j;
    static doublereal acoef, scale;
    static logical ilall;
    static integer iside;
    static doublereal sbeta;
    extern logical lsame_(char *, char *);
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *);
    static logical il2by2;
    static integer iinfo;
    static doublereal small;
    static logical compl;
    static doublereal anorm, bnorm;
    static logical compr;
    extern /* Subroutine */ int dlaln2_(logical *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     doublereal *, doublereal *, integer *, doublereal *, doublereal *
	    , doublereal *, integer *, doublereal *, doublereal *, integer *);
    static doublereal temp2i;
    extern /* Subroutine */ int dlabad_(doublereal *, doublereal *);
    static doublereal temp2r;
    static integer ja;
    static logical ilabad, ilbbad;
    static integer jc, je, na;
    static doublereal acoefa, bcoefa, cimaga, cimagb;
    static logical ilback;
    static integer im;
    static doublereal bcoefi, ascale, bscale, creala;
    static integer jr;
    static doublereal crealb;
    extern doublereal dlamch_(char *);
    static doublereal bcoefr;
    static integer jw, nw;
    static doublereal salfar, safmin;
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *);
    static doublereal xscale, bignum;
    extern /* Subroutine */ int xerbla_(char *, integer *);
    static logical ilcomp, ilcplx;
    static integer ihwmny;
    static doublereal big;
    static logical lsa, lsb;
    static doublereal ulp, sum[4]	/* was [2][2] */;



#define BDIAG(I) bdiag[(I)]
#define SELECT(I) select[(I)-1]
#define WORK(I) work[(I)-1]

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
	ilall = TRUE_;
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

/*     Count the number of eigenvectors to be computed */

    if (! ilall) {
	im = 0;
	ilcplx = FALSE_;
	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
	    if (ilcplx) {
		ilcplx = FALSE_;
		goto L10;
	    }
	    if (j < *n) {
		if (A(j+1,j) != 0.) {
		    ilcplx = TRUE_;
		}
	    }
	    if (ilcplx) {
		if (SELECT(j) || SELECT(j + 1)) {
		    im += 2;
		}
	    } else {
		if (SELECT(j)) {
		    ++im;
		}
	    }
L10:
	    ;
	}
    } else {
	im = *n;
    }

/*     Check 2-by-2 diagonal blocks of A, B */

    ilabad = FALSE_;
    ilbbad = FALSE_;
    i__1 = *n - 1;
    for (j = 1; j <= *n-1; ++j) {
	if (A(j+1,j) != 0.) {
	    if (B(j,j) == 0. || B(j+1,j+1) == 0. 
		    || B(j,j+1) != 0.) {
		ilbbad = TRUE_;
	    }
	    if (j < *n - 1) {
		if (A(j+2,j+1) != 0.) {
		    ilabad = TRUE_;
		}
	    }
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
    } else if (ilabad) {
	*info = -5;
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
	xerbla_("DTGEVC", &i__1);
	return 0;
    }

/*     Quick return if possible */

    *m = im;
    if (*n == 0) {
	return 0;
    }

/*     Machine Constants */

    safmin = dlamch_("Safe minimum");
    big = 1. / safmin;
    dlabad_(&safmin, &big);
    ulp = dlamch_("Epsilon") * dlamch_("Base");
    small = safmin * *n / ulp;
    big = 1. / small;
    bignum = 1. / (safmin * *n);

/*     Compute the 1-norm of each column of the strictly upper triangular 
  
       part (i.e., excluding all elements belonging to the diagonal   
       blocks) of A and B to check for possible overflow in the   
       triangular solver. */

    anorm = (d__1 = A(1,1), abs(d__1));
    if (*n > 1) {
	anorm += (d__1 = A(2,1), abs(d__1));
    }
    bnorm = (d__1 = B(1,1), abs(d__1));
    WORK(1) = 0.;
    WORK(*n + 1) = 0.;

    i__1 = *n;
    for (j = 2; j <= *n; ++j) {
	temp = 0.;
	temp2 = 0.;
	if (A(j,j-1) == 0.) {
	    iend = j - 1;
	} else {
	    iend = j - 2;
	}
	i__2 = iend;
	for (i = 1; i <= iend; ++i) {
	    temp += (d__1 = A(i,j), abs(d__1));
	    temp2 += (d__1 = B(i,j), abs(d__1));
/* L30: */
	}
	WORK(j) = temp;
	WORK(*n + j) = temp2;
/* Computing MIN */
	i__3 = j + 1;
	i__2 = min(i__3,*n);
	for (i = iend + 1; i <= min(j+1,*n); ++i) {
	    temp += (d__1 = A(i,j), abs(d__1));
	    temp2 += (d__1 = B(i,j), abs(d__1));
/* L40: */
	}
	anorm = max(anorm,temp);
	bnorm = max(bnorm,temp2);
/* L50: */
    }

    ascale = 1. / max(anorm,safmin);
    bscale = 1. / max(bnorm,safmin);

/*     Left eigenvectors */

    if (compl) {
	ieig = 0;

/*        Main loop over eigenvalues */

	ilcplx = FALSE_;
	i__1 = *n;
	for (je = 1; je <= *n; ++je) {

/*           Skip this iteration if (a) HOWMNY='S' and SELECT=.FAL
SE., or   
             (b) this would be the second of a complex pair.   
             Check for complex eigenvalue, so as to be sure of whi
ch   
             entry(-ies) of SELECT to look at. */

	    if (ilcplx) {
		ilcplx = FALSE_;
		goto L220;
	    }
	    nw = 1;
	    if (je < *n) {
		if (A(je+1,je) != 0.) {
		    ilcplx = TRUE_;
		    nw = 2;
		}
	    }
	    if (ilall) {
		ilcomp = TRUE_;
	    } else if (ilcplx) {
		ilcomp = SELECT(je) || SELECT(je + 1);
	    } else {
		ilcomp = SELECT(je);
	    }
	    if (! ilcomp) {
		goto L220;
	    }

/*           Decide if (a) singular pencil, (b) real eigenvalue, o
r   
             (c) complex eigenvalue. */

	    if (! ilcplx) {
		if ((d__1 = A(je,je), abs(d__1)) <= safmin && (
			d__2 = B(je,je), abs(d__2)) <= safmin) {

/*                 Singular matrix pencil -- return unit e
igenvector */

		    ++ieig;
		    i__2 = *n;
		    for (jr = 1; jr <= *n; ++jr) {
			VL(jr,ieig) = 0.;
/* L60: */
		    }
		    VL(ieig,ieig) = 1.;
		    goto L220;
		}
	    }

/*           Clear vector */

	    i__2 = nw * *n;
	    for (jr = 1; jr <= nw**n; ++jr) {
		WORK((*n << 1) + jr) = 0.;
/* L70: */
	    }
/*                                                 T   
             Compute coefficients in  ( a A - b B )  y = 0   
                a  is  ACOEF   
                b  is  BCOEFR + i*BCOEFI */

	    if (! ilcplx) {

/*              Real eigenvalue   

   Computing MAX */
		d__3 = (d__1 = A(je,je), abs(d__1)) * ascale, d__4 
			= (d__2 = B(je,je), abs(d__2)) * bscale, 
			d__3 = max(d__3,d__4);
		temp = 1. / max(d__3,safmin);
		salfar = temp * A(je,je) * ascale;
		sbeta = temp * B(je,je) * bscale;
		acoef = sbeta * ascale;
		bcoefr = salfar * bscale;
		bcoefi = 0.;

/*              Scale to avoid underflow */

		scale = 1.;
		lsa = abs(sbeta) >= safmin && abs(acoef) < small;
		lsb = abs(salfar) >= safmin && abs(bcoefr) < small;
		if (lsa) {
		    scale = small / abs(sbeta) * min(anorm,big);
		}
		if (lsb) {
/* Computing MAX */
		    d__1 = scale, d__2 = small / abs(salfar) * min(bnorm,big);
		    scale = max(d__1,d__2);
		}
		if (lsa || lsb) {
/* Computing MIN   
   Computing MAX */
		    d__3 = 1., d__4 = abs(acoef), d__3 = max(d__3,d__4), d__4 
			    = abs(bcoefr);
		    d__1 = scale, d__2 = 1. / (safmin * max(d__3,d__4));
		    scale = min(d__1,d__2);
		    if (lsa) {
			acoef = ascale * (scale * sbeta);
		    } else {
			acoef = scale * acoef;
		    }
		    if (lsb) {
			bcoefr = bscale * (scale * salfar);
		    } else {
			bcoefr = scale * bcoefr;
		    }
		}
		acoefa = abs(acoef);
		bcoefa = abs(bcoefr);

/*              First component is 1 */

		WORK((*n << 1) + je) = 1.;
		xmax = 1.;
	    } else {

/*              Complex eigenvalue */

		d__1 = safmin * 100.;
		dlag2_(&A(je,je), lda, &B(je,je), ldb, &
			d__1, &acoef, &temp, &bcoefr, &temp2, &bcoefi);
		bcoefi = -bcoefi;
		if (bcoefi == 0.) {
		    *info = je;
		    return 0;
		}

/*              Scale to avoid over/underflow */

		acoefa = abs(acoef);
		bcoefa = abs(bcoefr) + abs(bcoefi);
		scale = 1.;
		if (acoefa * ulp < safmin && acoefa >= safmin) {
		    scale = safmin / ulp / acoefa;
		}
		if (bcoefa * ulp < safmin && bcoefa >= safmin) {
/* Computing MAX */
		    d__1 = scale, d__2 = safmin / ulp / bcoefa;
		    scale = max(d__1,d__2);
		}
		if (safmin * acoefa > ascale) {
		    scale = ascale / (safmin * acoefa);
		}
		if (safmin * bcoefa > bscale) {
/* Computing MIN */
		    d__1 = scale, d__2 = bscale / (safmin * bcoefa);
		    scale = min(d__1,d__2);
		}
		if (scale != 1.) {
		    acoef = scale * acoef;
		    acoefa = abs(acoef);
		    bcoefr = scale * bcoefr;
		    bcoefi = scale * bcoefi;
		    bcoefa = abs(bcoefr) + abs(bcoefi);
		}

/*              Compute first two components of eigenvector */

		temp = acoef * A(je+1,je);
		temp2r = acoef * A(je,je) - bcoefr * B(je,je);
		temp2i = -bcoefi * B(je,je);
		if (abs(temp) > abs(temp2r) + abs(temp2i)) {
		    WORK((*n << 1) + je) = 1.;
		    WORK(*n * 3 + je) = 0.;
		    WORK((*n << 1) + je + 1) = -temp2r / temp;
		    WORK(*n * 3 + je + 1) = -temp2i / temp;
		} else {
		    WORK((*n << 1) + je + 1) = 1.;
		    WORK(*n * 3 + je + 1) = 0.;
		    temp = acoef * A(je,je+1);
		    WORK((*n << 1) + je) = (bcoefr * B(je+1,je+1) - acoef * A(je+1,je+1)) /
			     temp;
		    WORK(*n * 3 + je) = bcoefi * B(je+1,je+1)
			     / temp;
		}
/* Computing MAX */
		d__5 = (d__1 = WORK((*n << 1) + je), abs(d__1)) + (d__2 = 
			WORK(*n * 3 + je), abs(d__2)), d__6 = (d__3 = WORK((*
			n << 1) + je + 1), abs(d__3)) + (d__4 = WORK(*n * 3 + 
			je + 1), abs(d__4));
		xmax = max(d__5,d__6);
	    }

/* Computing MAX */
	    d__1 = ulp * acoefa * anorm, d__2 = ulp * bcoefa * bnorm, d__1 = 
		    max(d__1,d__2);
	    dmin__ = max(d__1,safmin);

/*                                           T   
             Triangular solve of  (a A - b B)  y = 0   

                                     T   
             (rowwise in  (a A - b B) , or columnwise in (a A - b 
B) ) */

	    il2by2 = FALSE_;

	    i__2 = *n;
	    for (j = je + nw; j <= *n; ++j) {
		if (il2by2) {
		    il2by2 = FALSE_;
		    goto L160;
		}

		na = 1;
		BDIAG(0) = B(j,j);
		if (j < *n) {
		    if (A(j+1,j) != 0.) {
			il2by2 = TRUE_;
			BDIAG(1) = B(j+1,j+1);
			na = 2;
		    }
		}

/*              Check whether scaling is necessary for dot pro
ducts */

		xscale = 1. / max(1.,xmax);
/* Computing MAX */
		d__1 = WORK(j), d__2 = WORK(*n + j), d__1 = max(d__1,d__2), 
			d__2 = acoefa * WORK(j) + bcoefa * WORK(*n + j);
		temp = max(d__1,d__2);
		if (il2by2) {
/* Computing MAX */
		    d__1 = temp, d__2 = WORK(j + 1), d__1 = max(d__1,d__2), 
			    d__2 = WORK(*n + j + 1), d__1 = max(d__1,d__2), 
			    d__2 = acoefa * WORK(j + 1) + bcoefa * WORK(*n + 
			    j + 1);
		    temp = max(d__1,d__2);
		}
		if (temp > bignum * xscale) {
		    i__3 = nw - 1;
		    for (jw = 0; jw <= nw-1; ++jw) {
			i__4 = j - 1;
			for (jr = je; jr <= j-1; ++jr) {
			    WORK((jw + 2) * *n + jr) = xscale * WORK((jw + 2) 
				    * *n + jr);
/* L80: */
			}
/* L90: */
		    }
		    xmax *= xscale;
		}

/*              Compute dot products   

                      j-1   
                SUM = sum  conjg( a*A(k,j) - b*B(k,j) )*x(k) 
  
                      k=je   

                To reduce the op count, this is done as   

                _        j-1                  _        j-1   
                a*conjg( sum  A(k,j)*x(k) ) - b*conjg( sum  B(
k,j)*x(k) )   
                         k=je                          k=je   

                which may cause underflow problems if A or B a
re close   
                to underflow.  (E.g., less than SMALL.)   


                A series of compiler directives to defeat vect
orization   
                for the next loop   

   $PL$ CMCHAR=' '   
   DIR$          NEXTSCALAR   
   $DIR          SCALAR   
   DIR$          NEXT SCALAR   
   VD$L          NOVECTOR   
   DEC$          NOVECTOR   
   VD$           NOVECTOR   
   VDIR          NOVECTOR   
   VOCL          LOOP,SCALAR   
   IBM           PREFER SCALAR   
   $PL$ CMCHAR='*' */

		i__3 = nw;
		for (jw = 1; jw <= nw; ++jw) {

/* $PL$ CMCHAR=' '   
   DIR$             NEXTSCALAR   
   $DIR             SCALAR   
   DIR$             NEXT SCALAR   
   VD$L             NOVECTOR   
   DEC$             NOVECTOR   
   VD$              NOVECTOR   
   VDIR             NOVECTOR   
   VOCL             LOOP,SCALAR   
   IBM              PREFER SCALAR   
   $PL$ CMCHAR='*' */

		    i__4 = na;
		    for (ja = 1; ja <= na; ++ja) {
			suma[ja + (jw << 1) - 3] = 0.;
			sumb[ja + (jw << 1) - 3] = 0.;

			i__5 = j - 1;
			for (jr = je; jr <= j-1; ++jr) {
			    suma[ja + (jw << 1) - 3] += A(jr,j+ja-1) * WORK((jw + 1) * *n + jr);
			    sumb[ja + (jw << 1) - 3] += B(jr,j+ja-1) * WORK((jw + 1) * *n + jr);
/* L100: */
			}
/* L110: */
		    }
/* L120: */
		}

/* $PL$ CMCHAR=' '   
   DIR$          NEXTSCALAR   
   $DIR          SCALAR   
   DIR$          NEXT SCALAR   
   VD$L          NOVECTOR   
   DEC$          NOVECTOR   
   VD$           NOVECTOR   
   VDIR          NOVECTOR   
   VOCL          LOOP,SCALAR   
   IBM           PREFER SCALAR   
   $PL$ CMCHAR='*' */

		i__3 = na;
		for (ja = 1; ja <= na; ++ja) {
		    if (ilcplx) {
			sum[ja - 1] = -acoef * suma[ja - 1] + bcoefr * sumb[
				ja - 1] - bcoefi * sumb[ja + 1];
			sum[ja + 1] = -acoef * suma[ja + 1] + bcoefr * sumb[
				ja + 1] + bcoefi * sumb[ja - 1];
		    } else {
			sum[ja - 1] = -acoef * suma[ja - 1] + bcoefr * sumb[
				ja - 1];
		    }
/* L130: */
		}

/*                                  T   
                Solve  ( a A - b B )  y = SUM(,)   
                with scaling and perturbation of the denominat
or */

		dlaln2_(&c_true, &na, &nw, &dmin__, &acoef, &A(j,j)
			, lda, bdiag, &BDIAG(1), sum, &c__2, &bcoefr, &bcoefi,
			 &WORK((*n << 1) + j), n, &scale, &temp, &iinfo);
		if (scale < 1.) {
		    i__3 = nw - 1;
		    for (jw = 0; jw <= nw-1; ++jw) {
			i__4 = j - 1;
			for (jr = je; jr <= j-1; ++jr) {
			    WORK((jw + 2) * *n + jr) = scale * WORK((jw + 2) *
				     *n + jr);
/* L140: */
			}
/* L150: */
		    }
		    xmax = scale * xmax;
		}
		xmax = max(xmax,temp);
L160:
		;
	    }

/*           Copy eigenvector to VL, back transforming if   
             HOWMNY='B'. */

	    ++ieig;
	    if (ilback) {
		i__2 = nw - 1;
		for (jw = 0; jw <= nw-1; ++jw) {
		    i__3 = *n + 1 - je;
		    dgemv_("N", n, &i__3, &c_b34, &VL(1,je), ldvl,
			     &WORK((jw + 2) * *n + je), &c__1, &c_b36, &WORK((
			    jw + 4) * *n + 1), &c__1);
/* L170: */
		}
		dlacpy_(" ", n, &nw, &WORK((*n << 2) + 1), n, &VL(1,je), ldvl);
		ibeg = 1;
	    } else {
		dlacpy_(" ", n, &nw, &WORK((*n << 1) + 1), n, &VL(1,ieig), ldvl);
		ibeg = je;
	    }

/*           Scale eigenvector */

	    xmax = 0.;
	    if (ilcplx) {
		i__2 = *n;
		for (j = ibeg; j <= *n; ++j) {
/* Computing MAX */
		    d__3 = xmax, d__4 = (d__1 = VL(j,ieig), abs(
			    d__1)) + (d__2 = VL(j,ieig+1), 
			    abs(d__2));
		    xmax = max(d__3,d__4);
/* L180: */
		}
	    } else {
		i__2 = *n;
		for (j = ibeg; j <= *n; ++j) {
/* Computing MAX */
		    d__2 = xmax, d__3 = (d__1 = VL(j,ieig), abs(
			    d__1));
		    xmax = max(d__2,d__3);
/* L190: */
		}
	    }

	    if (xmax > safmin) {
		xscale = 1. / xmax;

		i__2 = nw - 1;
		for (jw = 0; jw <= nw-1; ++jw) {
		    i__3 = *n;
		    for (jr = ibeg; jr <= *n; ++jr) {
			VL(jr,ieig+jw) = xscale * VL(jr,ieig+jw);
/* L200: */
		    }
/* L210: */
		}
	    }
	    ieig = ieig + nw - 1;

L220:
	    ;
	}
    }

/*     Right eigenvectors */

    if (compr) {
	ieig = im + 1;

/*        Main loop over eigenvalues */

	ilcplx = FALSE_;
	for (je = *n; je >= 1; --je) {

/*           Skip this iteration if (a) HOWMNY='S' and SELECT=.FAL
SE., or   
             (b) this would be the second of a complex pair.   
             Check for complex eigenvalue, so as to be sure of whi
ch   
             entry(-ies) of SELECT to look at -- if complex, SELEC
T(JE)   
             or SELECT(JE-1).   
             If this is a complex pair, the 2-by-2 diagonal block 
  
             corresponding to the eigenvalue is in rows/columns JE
-1:JE */

	    if (ilcplx) {
		ilcplx = FALSE_;
		goto L500;
	    }
	    nw = 1;
	    if (je > 1) {
		if (A(je,je-1) != 0.) {
		    ilcplx = TRUE_;
		    nw = 2;
		}
	    }
	    if (ilall) {
		ilcomp = TRUE_;
	    } else if (ilcplx) {
		ilcomp = SELECT(je) || SELECT(je - 1);
	    } else {
		ilcomp = SELECT(je);
	    }
	    if (! ilcomp) {
		goto L500;
	    }

/*           Decide if (a) singular pencil, (b) real eigenvalue, o
r   
             (c) complex eigenvalue. */

	    if (! ilcplx) {
		if ((d__1 = A(je,je), abs(d__1)) <= safmin && (
			d__2 = B(je,je), abs(d__2)) <= safmin) {

/*                 Singular matrix pencil -- unit eigenvec
tor */

		    --ieig;
		    i__1 = *n;
		    for (jr = 1; jr <= *n; ++jr) {
			VR(jr,ieig) = 0.;
/* L230: */
		    }
		    VR(ieig,ieig) = 1.;
		    goto L500;
		}
	    }

/*           Clear vector */

	    i__1 = nw - 1;
	    for (jw = 0; jw <= nw-1; ++jw) {
		i__2 = *n;
		for (jr = 1; jr <= *n; ++jr) {
		    WORK((jw + 2) * *n + jr) = 0.;
/* L240: */
		}
/* L250: */
	    }

/*           Compute coefficients in  ( a A - b B ) x = 0   
                a  is  ACOEF   
                b  is  BCOEFR + i*BCOEFI */

	    if (! ilcplx) {

/*              Real eigenvalue   

   Computing MAX */
		d__3 = (d__1 = A(je,je), abs(d__1)) * ascale, d__4 
			= (d__2 = B(je,je), abs(d__2)) * bscale, 
			d__3 = max(d__3,d__4);
		temp = 1. / max(d__3,safmin);
		salfar = temp * A(je,je) * ascale;
		sbeta = temp * B(je,je) * bscale;
		acoef = sbeta * ascale;
		bcoefr = salfar * bscale;
		bcoefi = 0.;

/*              Scale to avoid underflow */

		scale = 1.;
		lsa = abs(sbeta) >= safmin && abs(acoef) < small;
		lsb = abs(salfar) >= safmin && abs(bcoefr) < small;
		if (lsa) {
		    scale = small / abs(sbeta) * min(anorm,big);
		}
		if (lsb) {
/* Computing MAX */
		    d__1 = scale, d__2 = small / abs(salfar) * min(bnorm,big);
		    scale = max(d__1,d__2);
		}
		if (lsa || lsb) {
/* Computing MIN   
   Computing MAX */
		    d__3 = 1., d__4 = abs(acoef), d__3 = max(d__3,d__4), d__4 
			    = abs(bcoefr);
		    d__1 = scale, d__2 = 1. / (safmin * max(d__3,d__4));
		    scale = min(d__1,d__2);
		    if (lsa) {
			acoef = ascale * (scale * sbeta);
		    } else {
			acoef = scale * acoef;
		    }
		    if (lsb) {
			bcoefr = bscale * (scale * salfar);
		    } else {
			bcoefr = scale * bcoefr;
		    }
		}
		acoefa = abs(acoef);
		bcoefa = abs(bcoefr);

/*              First component is 1 */

		WORK((*n << 1) + je) = 1.;
		xmax = 1.;

/*              Compute contribution from column JE of A and B
 to sum   
                (See "Further Details", above.) */

		i__1 = je - 1;
		for (jr = 1; jr <= je-1; ++jr) {
		    WORK((*n << 1) + jr) = bcoefr * B(jr,je) - 
			    acoef * A(jr,je);
/* L260: */
		}
	    } else {

/*              Complex eigenvalue */

		d__1 = safmin * 100.;
		dlag2_(&A(je-1,je-1), lda, &B(je-1,je-1), ldb, &d__1, &acoef, &temp, &bcoefr, &
			temp2, &bcoefi);
		if (bcoefi == 0.) {
		    *info = je - 1;
		    return 0;
		}

/*              Scale to avoid over/underflow */

		acoefa = abs(acoef);
		bcoefa = abs(bcoefr) + abs(bcoefi);
		scale = 1.;
		if (acoefa * ulp < safmin && acoefa >= safmin) {
		    scale = safmin / ulp / acoefa;
		}
		if (bcoefa * ulp < safmin && bcoefa >= safmin) {
/* Computing MAX */
		    d__1 = scale, d__2 = safmin / ulp / bcoefa;
		    scale = max(d__1,d__2);
		}
		if (safmin * acoefa > ascale) {
		    scale = ascale / (safmin * acoefa);
		}
		if (safmin * bcoefa > bscale) {
/* Computing MIN */
		    d__1 = scale, d__2 = bscale / (safmin * bcoefa);
		    scale = min(d__1,d__2);
		}
		if (scale != 1.) {
		    acoef = scale * acoef;
		    acoefa = abs(acoef);
		    bcoefr = scale * bcoefr;
		    bcoefi = scale * bcoefi;
		    bcoefa = abs(bcoefr) + abs(bcoefi);
		}

/*              Compute first two components of eigenvector   
                and contribution to sums */

		temp = acoef * A(je,je-1);
		temp2r = acoef * A(je,je) - bcoefr * B(je,je);
		temp2i = -bcoefi * B(je,je);
		if (abs(temp) >= abs(temp2r) + abs(temp2i)) {
		    WORK((*n << 1) + je) = 1.;
		    WORK(*n * 3 + je) = 0.;
		    WORK((*n << 1) + je - 1) = -temp2r / temp;
		    WORK(*n * 3 + je - 1) = -temp2i / temp;
		} else {
		    WORK((*n << 1) + je - 1) = 1.;
		    WORK(*n * 3 + je - 1) = 0.;
		    temp = acoef * A(je-1,je);
		    WORK((*n << 1) + je) = (bcoefr * B(je-1,je-1) - acoef * A(je-1,je-1)) /
			     temp;
		    WORK(*n * 3 + je) = bcoefi * B(je-1,je-1)
			     / temp;
		}

/* Computing MAX */
		d__5 = (d__1 = WORK((*n << 1) + je), abs(d__1)) + (d__2 = 
			WORK(*n * 3 + je), abs(d__2)), d__6 = (d__3 = WORK((*
			n << 1) + je - 1), abs(d__3)) + (d__4 = WORK(*n * 3 + 
			je - 1), abs(d__4));
		xmax = max(d__5,d__6);

/*              Compute contribution from columns JE and JE-1 
  
                of A and B to the sums. */

		creala = acoef * WORK((*n << 1) + je - 1);
		cimaga = acoef * WORK(*n * 3 + je - 1);
		crealb = bcoefr * WORK((*n << 1) + je - 1) - bcoefi * WORK(*n 
			* 3 + je - 1);
		cimagb = bcoefi * WORK((*n << 1) + je - 1) + bcoefr * WORK(*n 
			* 3 + je - 1);
		cre2a = acoef * WORK((*n << 1) + je);
		cim2a = acoef * WORK(*n * 3 + je);
		cre2b = bcoefr * WORK((*n << 1) + je) - bcoefi * WORK(*n * 3 
			+ je);
		cim2b = bcoefi * WORK((*n << 1) + je) + bcoefr * WORK(*n * 3 
			+ je);
		i__1 = je - 2;
		for (jr = 1; jr <= je-2; ++jr) {
		    WORK((*n << 1) + jr) = -creala * A(jr,je-1)
			     + crealb * B(jr,je-1) - cre2a * A(jr,je) + cre2b * B(jr,je);
		    WORK(*n * 3 + jr) = -cimaga * A(jr,je-1) + 
			    cimagb * B(jr,je-1) - cim2a * A(jr,je) + cim2b * B(jr,je);
/* L270: */
		}
	    }

/* Computing MAX */
	    d__1 = ulp * acoefa * anorm, d__2 = ulp * bcoefa * bnorm, d__1 = 
		    max(d__1,d__2);
	    dmin__ = max(d__1,safmin);

/*           Columnwise triangular solve of  (a A - b B)  x = 0 */

	    il2by2 = FALSE_;
	    for (j = je - nw; j >= 1; --j) {

/*              If a 2-by-2 block, is in position j-1:j, wait 
until   
                next iteration to process it (when it will be 
j:j+1) */

		if (! il2by2 && j > 1) {
		    if (A(j,j-1) != 0.) {
			il2by2 = TRUE_;
			goto L370;
		    }
		}
		BDIAG(0) = B(j,j);
		if (il2by2) {
		    na = 2;
		    BDIAG(1) = B(j+1,j+1);
		} else {
		    na = 1;
		}

/*              Compute x(j) (and x(j+1), if 2-by-2 block) */

		dlaln2_(&c_false, &na, &nw, &dmin__, &acoef, &A(j,j), lda, bdiag, &BDIAG(1), &WORK((*n << 1) + j), 
			n, &bcoefr, &bcoefi, sum, &c__2, &scale, &temp, &
			iinfo);
		if (scale < 1.) {

		    i__1 = nw - 1;
		    for (jw = 0; jw <= nw-1; ++jw) {
			i__2 = je;
			for (jr = 1; jr <= je; ++jr) {
			    WORK((jw + 2) * *n + jr) = scale * WORK((jw + 2) *
				     *n + jr);
/* L280: */
			}
/* L290: */
		    }
		}
/* Computing MAX */
		d__1 = scale * xmax;
		xmax = max(d__1,temp);

		i__1 = nw;
		for (jw = 1; jw <= nw; ++jw) {
		    i__2 = na;
		    for (ja = 1; ja <= na; ++ja) {
			WORK((jw + 1) * *n + j + ja - 1) = sum[ja + (jw << 1) 
				- 3];
/* L300: */
		    }
/* L310: */
		}

/*              w = w + x(j)*(a A(*,j) - b B(*,j) ) with scali
ng */

		if (j > 1) {

/*                 Check whether scaling is necessary for 
sum. */

		    xscale = 1. / max(1.,xmax);
		    temp = acoefa * WORK(j) + bcoefa * WORK(*n + j);
		    if (il2by2) {
/* Computing MAX */
			d__1 = temp, d__2 = acoefa * WORK(j + 1) + bcoefa * 
				WORK(*n + j + 1);
			temp = max(d__1,d__2);
		    }
/* Computing MAX */
		    d__1 = max(temp,acoefa);
		    temp = max(d__1,bcoefa);
		    if (temp > bignum * xscale) {

			i__1 = nw - 1;
			for (jw = 0; jw <= nw-1; ++jw) {
			    i__2 = je;
			    for (jr = 1; jr <= je; ++jr) {
				WORK((jw + 2) * *n + jr) = xscale * WORK((jw 
					+ 2) * *n + jr);
/* L320: */
			    }
/* L330: */
			}
			xmax *= xscale;
		    }

/*                 Compute the contributions of the off-di
agonals of   
                   column j (and j+1, if 2-by-2 block) of 
A and B to the   
                   sums. */


		    i__1 = na;
		    for (ja = 1; ja <= na; ++ja) {
			if (ilcplx) {
			    creala = acoef * WORK((*n << 1) + j + ja - 1);
			    cimaga = acoef * WORK(*n * 3 + j + ja - 1);
			    crealb = bcoefr * WORK((*n << 1) + j + ja - 1) - 
				    bcoefi * WORK(*n * 3 + j + ja - 1);
			    cimagb = bcoefi * WORK((*n << 1) + j + ja - 1) + 
				    bcoefr * WORK(*n * 3 + j + ja - 1);
			    i__2 = j - 1;
			    for (jr = 1; jr <= j-1; ++jr) {
				WORK((*n << 1) + jr) = WORK((*n << 1) + jr) - 
					creala * A(jr,j+ja-1)
					 + crealb * B(jr,j+ja-1);
				WORK(*n * 3 + jr) = WORK(*n * 3 + jr) - 
					cimaga * A(jr,j+ja-1)
					 + cimagb * B(jr,j+ja-1);
/* L340: */
			    }
			} else {
			    creala = acoef * WORK((*n << 1) + j + ja - 1);
			    crealb = bcoefr * WORK((*n << 1) + j + ja - 1);
			    i__2 = j - 1;
			    for (jr = 1; jr <= j-1; ++jr) {
				WORK((*n << 1) + jr) = WORK((*n << 1) + jr) - 
					creala * A(jr,j+ja-1)
					 + crealb * B(jr,j+ja-1);
/* L350: */
			    }
			}
/* L360: */
		    }
		}

		il2by2 = FALSE_;
L370:
		;
	    }

/*           Copy eigenvector to VR, back transforming if   
             HOWMNY='B'. */

	    ieig -= nw;
	    if (ilback) {

		i__1 = nw - 1;
		for (jw = 0; jw <= nw-1; ++jw) {
		    i__2 = *n;
		    for (jr = 1; jr <= *n; ++jr) {
			WORK((jw + 4) * *n + jr) = WORK((jw + 2) * *n + 1) * 
				VR(jr,1);
/* L380: */
		    }

/*                 A series of compiler directives to defe
at   
                   vectorization for the next loop */


		    i__2 = je;
		    for (jc = 2; jc <= je; ++jc) {
			i__3 = *n;
			for (jr = 1; jr <= *n; ++jr) {
			    WORK((jw + 4) * *n + jr) += WORK((jw + 2) * *n + 
				    jc) * VR(jr,jc);
/* L390: */
			}
/* L400: */
		    }
/* L410: */
		}

		i__1 = nw - 1;
		for (jw = 0; jw <= nw-1; ++jw) {
		    i__2 = *n;
		    for (jr = 1; jr <= *n; ++jr) {
			VR(jr,ieig+jw) = WORK((jw + 4) * *n + 
				jr);
/* L420: */
		    }
/* L430: */
		}

		iend = *n;
	    } else {
		i__1 = nw - 1;
		for (jw = 0; jw <= nw-1; ++jw) {
		    i__2 = *n;
		    for (jr = 1; jr <= *n; ++jr) {
			VR(jr,ieig+jw) = WORK((jw + 2) * *n + 
				jr);
/* L440: */
		    }
/* L450: */
		}

		iend = je;
	    }

/*           Scale eigenvector */

	    xmax = 0.;
	    if (ilcplx) {
		i__1 = iend;
		for (j = 1; j <= iend; ++j) {
/* Computing MAX */
		    d__3 = xmax, d__4 = (d__1 = VR(j,ieig), abs(
			    d__1)) + (d__2 = VR(j,ieig+1), 
			    abs(d__2));
		    xmax = max(d__3,d__4);
/* L460: */
		}
	    } else {
		i__1 = iend;
		for (j = 1; j <= iend; ++j) {
/* Computing MAX */
		    d__2 = xmax, d__3 = (d__1 = VR(j,ieig), abs(
			    d__1));
		    xmax = max(d__2,d__3);
/* L470: */
		}
	    }

	    if (xmax > safmin) {
		xscale = 1. / xmax;
		i__1 = nw - 1;
		for (jw = 0; jw <= nw-1; ++jw) {
		    i__2 = iend;
		    for (jr = 1; jr <= iend; ++jr) {
			VR(jr,ieig+jw) = xscale * VR(jr,ieig+jw);
/* L480: */
		    }
/* L490: */
		}
	    }
L500:
	    ;
	}
    }

    return 0;

/*     End of DTGEVC */

} /* dtgevc_ */

