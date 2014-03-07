#include "f2c.h"

/* Subroutine */ int dlaed9_(integer *k, integer *kstart, integer *kstop, 
	integer *n, doublereal *d, doublereal *q, integer *ldq, doublereal *
	rho, doublereal *dlamda, doublereal *w, doublereal *s, integer *lds, 
	integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Oak Ridge National Lab, Argonne National Lab, 
  
       Courant Institute, NAG Ltd., and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DLAED9 finds the roots of the secular equation, as defined by the   
    values in D, Z, and RHO, between KSTART and KSTOP.  It makes the   
    appropriate calls to DLAED4 and then stores the new matrix of   
    eigenvectors for use in calculating the next level of Z vectors.   

    Arguments   
    =========   

    K       (input) INTEGER   
            The number of terms in the rational function to be solved by 
  
            DLAED4.  K >= 0.   

    KSTART  (input) INTEGER   
    KSTOP   (input) INTEGER   
            The updated eigenvalues Lambda(I), KSTART <= I <= KSTOP   
            are to be computed.  1 <= KSTART <= KSTOP <= K.   

    N       (input) INTEGER   
            The number of rows and columns in the Q matrix.   
            N >= K (delation may result in N > K).   

    D       (output) DOUBLE PRECISION array, dimension (N)   
            D(I) contains the updated eigenvalues   
            for KSTART <= I <= KSTOP.   

    Q       (workspace) DOUBLE PRECISION array, dimension (LDQ,N)   

    LDQ     (input) INTEGER   
            The leading dimension of the array Q.  LDQ >= max( 1, N ).   

    RHO     (input) DOUBLE PRECISION   
            The value of the parameter in the rank one update equation.   
            RHO >= 0 required.   

    DLAMDA  (input) DOUBLE PRECISION array, dimension (K)   
            The first K elements of this array contain the old roots   
            of the deflated updating problem.  These are the poles   
            of the secular equation.   

    W       (input) DOUBLE PRECISION array, dimension (K)   
            The first K elements of this array contain the components   
            of the deflation-adjusted updating vector.   

    S       (output) DOUBLE PRECISION array, dimension (LDS, K)   
            Will contain the eigenvectors of the repaired matrix which   
            will be stored for subsequent Z vector calculation and   
            multiplied by the previously accumulated eigenvectors   
            to update the system.   

    LDS     (input) INTEGER   
            The leading dimension of S.  LDS >= max( 1, K ).   

    INFO    (output) INTEGER   
            = 0:  successful exit.   
            < 0:  if INFO = -i, the i-th argument had an illegal value.   
            > 0:  if INFO = 1, an eigenvalue did not converge   

    ===================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static integer c__1 = 1;
    
    /* System generated locals */
    integer q_dim1, q_offset, s_dim1, s_offset, i__1, i__2;
    doublereal d__1;
    /* Builtin functions */
    double sqrt(doublereal), d_sign(doublereal *, doublereal *);
    /* Local variables */
    static doublereal temp;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    static integer i, j;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dlaed4_(integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *);
    extern doublereal dlamc3_(doublereal *, doublereal *);
    extern /* Subroutine */ int xerbla_(char *, integer *);



#define D(I) d[(I)-1]
#define DLAMDA(I) dlamda[(I)-1]
#define W(I) w[(I)-1]

#define Q(I,J) q[(I)-1 + ((J)-1)* ( *ldq)]
#define S(I,J) s[(I)-1 + ((J)-1)* ( *lds)]

    *info = 0;

    if (*k < 0) {
	*info = -1;
    } else if (*kstart < 1 || *kstart > max(1,*k)) {
	*info = -2;
    } else if (max(1,*kstop) < *kstart || *kstop > max(1,*k)) {
	*info = -3;
    } else if (*n < *k) {
	*info = -4;
    } else if (*ldq < max(1,*k)) {
	*info = -7;
    } else if (*lds < max(1,*k)) {
	*info = -12;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DLAED9", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*k == 0) {
	return 0;
    }

/*     Modify values DLAMDA(i) to make sure all DLAMDA(i)-DLAMDA(j) can   
       be computed with high relative accuracy (barring over/underflow). 
  
       This is a problem on machines without a guard digit in   
       add/subtract (Cray XMP, Cray YMP, Cray C 90 and Cray 2).   
       The following code replaces DLAMDA(I) by 2*DLAMDA(I)-DLAMDA(I),   
       which on any of these machines zeros out the bottommost   
       bit of DLAMDA(I) if it is 1; this makes the subsequent   
       subtractions DLAMDA(I)-DLAMDA(J) unproblematic when cancellation   
       occurs. On binary machines with a guard digit (almost all   
       machines) it does not change DLAMDA(I) at all. On hexadecimal   
       and decimal machines with a guard digit, it slightly   
       changes the bottommost bits of DLAMDA(I). It does not account   
       for hexadecimal or decimal machines without guard digits   
       (we know of none). We use a subroutine call to compute   
       2*DLAMBDA(I) to prevent optimizing compilers from eliminating   
       this code. */

    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
	DLAMDA(i) = dlamc3_(&DLAMDA(i), &DLAMDA(i)) - DLAMDA(i);
/* L10: */
    }

    i__1 = *kstop;
    for (j = *kstart; j <= *kstop; ++j) {
	dlaed4_(k, &j, &DLAMDA(1), &W(1), &Q(1,j), rho, &D(j), 
		info);

/*        If the zero finder fails, the computation is terminated. */

	if (*info != 0) {
	    goto L120;
	}
/* L20: */
    }

    if (*k == 1 || *k == 2) {
	i__1 = *k;
	for (i = 1; i <= *k; ++i) {
	    i__2 = *k;
	    for (j = 1; j <= *k; ++j) {
		S(j,i) = Q(j,i);
/* L30: */
	    }
/* L40: */
	}
	goto L120;
    }

/*     Compute updated W. */

    dcopy_(k, &W(1), &c__1, &S(1,1), &c__1);

/*     Initialize W(I) = Q(I,I) */

    i__1 = *ldq + 1;
    dcopy_(k, &Q(1,1), &i__1, &W(1), &c__1);
    i__1 = *k;
    for (j = 1; j <= *k; ++j) {
	i__2 = j - 1;
	for (i = 1; i <= j-1; ++i) {
	    W(i) *= Q(i,j) / (DLAMDA(i) - DLAMDA(j));
/* L50: */
	}
	i__2 = *k;
	for (i = j + 1; i <= *k; ++i) {
	    W(i) *= Q(i,j) / (DLAMDA(i) - DLAMDA(j));
/* L60: */
	}
/* L70: */
    }
    i__1 = *k;
    for (i = 1; i <= *k; ++i) {
	d__1 = sqrt(-W(i));
	W(i) = d_sign(&d__1, &S(i,1));
/* L80: */
    }

/*     Compute eigenvectors of the modified rank-1 modification. */

    i__1 = *k;
    for (j = 1; j <= *k; ++j) {
	i__2 = *k;
	for (i = 1; i <= *k; ++i) {
	    Q(i,j) = W(i) / Q(i,j);
/* L90: */
	}
	temp = dnrm2_(k, &Q(1,j), &c__1);
	i__2 = *k;
	for (i = 1; i <= *k; ++i) {
	    S(i,j) = Q(i,j) / temp;
/* L100: */
	}
/* L110: */
    }

L120:
    return 0;

/*     End of DLAED9 */

} /* dlaed9_ */

