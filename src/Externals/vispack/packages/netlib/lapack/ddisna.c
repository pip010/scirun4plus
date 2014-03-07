#include "f2c.h"

/* Subroutine */ int ddisna_(char *job, integer *m, integer *n, doublereal *d,
	 doublereal *sep, integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DDISNA computes the reciprocal condition numbers for the eigenvectors 
  
    of a real symmetric or complex Hermitian matrix or for the left or   
    right singular vectors of a general m-by-n matrix. The reciprocal   
    condition number is the 'gap' between the corresponding eigenvalue or 
  
    singular value and the nearest other one.   

    The bound on the error, measured by angle in radians, in the I-th   
    computed vector is given by   

           DLAMCH( 'E' ) * ( ANORM / SEP( I ) )   

    where ANORM = 2-norm(A) = max( abs( D(j) ) ).  SEP(I) is not allowed 
  
    to be smaller than DLAMCH( 'E' )*ANORM in order to limit the size of 
  
    the error bound.   

    DDISNA may also be used to compute error bounds for eigenvectors of   
    the generalized symmetric definite eigenproblem.   

    Arguments   
    =========   

    JOB     (input) CHARACTER*1   
            Specifies for which problem the reciprocal condition numbers 
  
            should be computed:   
            = 'E':  the eigenvectors of a symmetric/Hermitian matrix;   
            = 'L':  the left singular vectors of a general matrix;   
            = 'R':  the right singular vectors of a general matrix.   

    M       (input) INTEGER   
            The number of rows of the matrix. M >= 0.   

    N       (input) INTEGER   
            If JOB = 'L' or 'R', the number of columns of the matrix,   
            in which case N >= 0. Ignored if JOB = 'E'.   

    D       (input) DOUBLE PRECISION array, dimension (M) if JOB = 'E'   
                                dimension (min(M,N)) if JOB = 'L' or 'R' 
  
            The eigenvalues (if JOB = 'E') or singular values (if JOB =   
            'L' or 'R') of the matrix, in either increasing or decreasing 
  
            order. If singular values, they must be non-negative.   

    SEP     (output) DOUBLE PRECISION array, dimension (M) if JOB = 'E'   
                                 dimension (min(M,N)) if JOB = 'L' or 'R' 
  
            The reciprocal condition numbers of the vectors.   

    INFO    (output) INTEGER   
            = 0:  successful exit.   
            < 0:  if INFO = -i, the i-th argument had an illegal value.   

    ===================================================================== 
  


       Test the input arguments   

    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2, d__3;
    /* Local variables */
    static logical decr, left, incr, sing;
    static integer i, k;
    static logical eigen;
    extern logical lsame_(char *, char *);
    static doublereal anorm;
    static logical right;
    extern doublereal dlamch_(char *);
    static doublereal oldgap, safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *);
    static doublereal newgap, thresh, eps;


#define SEP(I) sep[(I)-1]
#define D(I) d[(I)-1]


    *info = 0;
    eigen = lsame_(job, "E");
    left = lsame_(job, "L");
    right = lsame_(job, "R");
    sing = left || right;
    if (eigen) {
	k = *m;
    } else if (sing) {
	k = min(*m,*n);
    }
    if (! eigen && ! sing) {
	*info = -1;
    } else if (*m < 0) {
	*info = -2;
    } else if (k < 0) {
	*info = -3;
    } else {
	incr = TRUE_;
	decr = TRUE_;
	i__1 = k - 1;
	for (i = 1; i <= k-1; ++i) {
	    if (incr) {
		incr = incr && D(i) <= D(i + 1);
	    }
	    if (decr) {
		decr = decr && D(i) >= D(i + 1);
	    }
/* L10: */
	}
	if (sing && k > 0) {
	    if (incr) {
		incr = incr && 0. <= D(1);
	    }
	    if (decr) {
		decr = decr && D(k) >= 0.;
	    }
	}
	if (! (incr || decr)) {
	    *info = -4;
	}
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DDISNA", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (k == 0) {
	return 0;
    }

/*     Compute reciprocal condition numbers */

    if (k == 1) {
	SEP(1) = dlamch_("O");
    } else {
	oldgap = (d__1 = D(2) - D(1), abs(d__1));
	SEP(1) = oldgap;
	i__1 = k - 1;
	for (i = 2; i <= k-1; ++i) {
	    newgap = (d__1 = D(i + 1) - D(i), abs(d__1));
	    SEP(i) = min(oldgap,newgap);
	    oldgap = newgap;
/* L20: */
	}
	SEP(k) = oldgap;
    }
    if (sing) {
	if (left && *m > *n || right && *m < *n) {
	    if (incr) {
		SEP(1) = min(SEP(1),D(1));
	    }
	    if (decr) {
/* Computing MIN */
		d__1 = SEP(k), d__2 = D(k);
		SEP(k) = min(d__1,d__2);
	    }
	}
    }

/*     Ensure that reciprocal condition numbers are not less than   
       threshold, in order to limit the size of the error bound */

    eps = dlamch_("E");
    safmin = dlamch_("S");
/* Computing MAX */
    d__2 = abs(D(1)), d__3 = (d__1 = D(k), abs(d__1));
    anorm = max(d__2,d__3);
    if (anorm == 0.) {
	thresh = eps;
    } else {
/* Computing MAX */
	d__1 = eps * anorm;
	thresh = max(d__1,safmin);
    }
    i__1 = k;
    for (i = 1; i <= k; ++i) {
/* Computing MAX */
	d__1 = SEP(i);
	SEP(i) = max(d__1,thresh);
/* L30: */
    }

    return 0;

/*     End of DDISNA */

} /* ddisna_ */

