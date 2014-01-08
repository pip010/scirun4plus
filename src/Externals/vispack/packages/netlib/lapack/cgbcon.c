#include "f2c.h"

/* Subroutine */ int cgbcon_(char *norm, integer *n, integer *kl, integer *ku,
	 complex *ab, integer *ldab, integer *ipiv, real *anorm, real *rcond, 
	complex *work, real *rwork, integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    CGBCON estimates the reciprocal of the condition number of a complex 
  
    general band matrix A, in either the 1-norm or the infinity-norm,   
    using the LU factorization computed by CGBTRF.   

    An estimate is obtained for norm(inv(A)), and the reciprocal of the   
    condition number is computed as   
       RCOND = 1 / ( norm(A) * norm(inv(A)) ).   

    Arguments   
    =========   

    NORM    (input) CHARACTER*1   
            Specifies whether the 1-norm condition number or the   
            infinity-norm condition number is required:   
            = '1' or 'O':  1-norm;   
            = 'I':         Infinity-norm.   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    KL      (input) INTEGER   
            The number of subdiagonals within the band of A.  KL >= 0.   

    KU      (input) INTEGER   
            The number of superdiagonals within the band of A.  KU >= 0. 
  

    AB      (input) COMPLEX array, dimension (LDAB,N)   
            Details of the LU factorization of the band matrix A, as   
            computed by CGBTRF.  U is stored as an upper triangular band 
  
            matrix with KL+KU superdiagonals in rows 1 to KL+KU+1, and   
            the multipliers used during the factorization are stored in   
            rows KL+KU+2 to 2*KL+KU+1.   

    LDAB    (input) INTEGER   
            The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.   

    IPIV    (input) INTEGER array, dimension (N)   
            The pivot indices; for 1 <= i <= N, row i of the matrix was   
            interchanged with row IPIV(i).   

    ANORM   (input) REAL   
            If NORM = '1' or 'O', the 1-norm of the original matrix A.   
            If NORM = 'I', the infinity-norm of the original matrix A.   

    RCOND   (output) REAL   
            The reciprocal of the condition number of the matrix A,   
            computed as RCOND = 1/(norm(A) * norm(inv(A))).   

    WORK    (workspace) COMPLEX array, dimension (2*N)   

    RWORK   (workspace) REAL array, dimension (N)   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0: if INFO = -i, the i-th argument had an illegal value   

    ===================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static integer c__1 = 1;
    
    /* System generated locals */
    integer ab_dim1, ab_offset, i__1, i__2, i__3;
    real r__1, r__2;
    complex q__1, q__2;
    /* Builtin functions */
    double r_imag(complex *);
    /* Local variables */
    static integer kase, kase1, j;
    static complex t;
    static real scale;
    extern /* Complex */ VOID cdotc_(complex *, integer *, complex *, integer 
	    *, complex *, integer *);
    extern logical lsame_(char *, char *);
    extern /* Subroutine */ int caxpy_(integer *, complex *, complex *, 
	    integer *, complex *, integer *);
    static logical lnoti;
    static integer kd, lm, jp;
    extern /* Subroutine */ int clacon_(integer *, complex *, complex *, real 
	    *, integer *);
    static integer ix;
    extern integer icamax_(integer *, complex *, integer *);
    extern doublereal slamch_(char *);
    extern /* Subroutine */ int clatbs_(char *, char *, char *, char *, 
	    integer *, integer *, complex *, integer *, complex *, real *, 
	    real *, integer *), xerbla_(char *
	    , integer *);
    static real ainvnm;
    extern /* Subroutine */ int csrscl_(integer *, real *, complex *, integer 
	    *);
    static logical onenrm;
    static char normin[1];
    static real smlnum;



#define IPIV(I) ipiv[(I)-1]
#define WORK(I) work[(I)-1]
#define RWORK(I) rwork[(I)-1]

#define AB(I,J) ab[(I)-1 + ((J)-1)* ( *ldab)]

    *info = 0;
    onenrm = *(unsigned char *)norm == '1' || lsame_(norm, "O");
    if (! onenrm && ! lsame_(norm, "I")) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*kl < 0) {
	*info = -3;
    } else if (*ku < 0) {
	*info = -4;
    } else if (*ldab < (*kl << 1) + *ku + 1) {
	*info = -6;
    } else if (*anorm < 0.f) {
	*info = -8;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("CGBCON", &i__1);
	return 0;
    }

/*     Quick return if possible */

    *rcond = 0.f;
    if (*n == 0) {
	*rcond = 1.f;
	return 0;
    } else if (*anorm == 0.f) {
	return 0;
    }

    smlnum = slamch_("Safe minimum");

/*     Estimate the norm of inv(A). */

    ainvnm = 0.f;
    *(unsigned char *)normin = 'N';
    if (onenrm) {
	kase1 = 1;
    } else {
	kase1 = 2;
    }
    kd = *kl + *ku + 1;
    lnoti = *kl > 0;
    kase = 0;
L10:
    clacon_(n, &WORK(*n + 1), &WORK(1), &ainvnm, &kase);
    if (kase != 0) {
	if (kase == kase1) {

/*           Multiply by inv(L). */

	    if (lnoti) {
		i__1 = *n - 1;
		for (j = 1; j <= *n-1; ++j) {
/* Computing MIN */
		    i__2 = *kl, i__3 = *n - j;
		    lm = min(i__2,i__3);
		    jp = IPIV(j);
		    i__2 = jp;
		    t.r = WORK(jp).r, t.i = WORK(jp).i;
		    if (jp != j) {
			i__2 = jp;
			i__3 = j;
			WORK(jp).r = WORK(j).r, WORK(jp).i = WORK(j)
				.i;
			i__2 = j;
			WORK(j).r = t.r, WORK(j).i = t.i;
		    }
		    q__1.r = -(doublereal)t.r, q__1.i = -(doublereal)t.i;
		    caxpy_(&lm, &q__1, &AB(kd+1,j), &c__1, &
			    WORK(j + 1), &c__1);
/* L20: */
		}
	    }

/*           Multiply by inv(U). */

	    i__1 = *kl + *ku;
	    clatbs_("Upper", "No transpose", "Non-unit", normin, n, &i__1, &
		    AB(1,1), ldab, &WORK(1), &scale, &RWORK(1), info);
	} else {

/*           Multiply by inv(U'). */

	    i__1 = *kl + *ku;
	    clatbs_("Upper", "Conjugate transpose", "Non-unit", normin, n, &
		    i__1, &AB(1,1), ldab, &WORK(1), &scale, &RWORK(1), 
		    info);

/*           Multiply by inv(L'). */

	    if (lnoti) {
		for (j = *n - 1; j >= 1; --j) {
/* Computing MIN */
		    i__1 = *kl, i__2 = *n - j;
		    lm = min(i__1,i__2);
		    i__1 = j;
		    i__2 = j;
		    cdotc_(&q__2, &lm, &AB(kd+1,j), &c__1, &
			    WORK(j + 1), &c__1);
		    q__1.r = WORK(j).r - q__2.r, q__1.i = WORK(j).i - 
			    q__2.i;
		    WORK(j).r = q__1.r, WORK(j).i = q__1.i;
		    jp = IPIV(j);
		    if (jp != j) {
			i__1 = jp;
			t.r = WORK(jp).r, t.i = WORK(jp).i;
			i__1 = jp;
			i__2 = j;
			WORK(jp).r = WORK(j).r, WORK(jp).i = WORK(j)
				.i;
			i__1 = j;
			WORK(j).r = t.r, WORK(j).i = t.i;
		    }
/* L30: */
		}
	    }
	}

/*        Divide X by 1/SCALE if doing so will not cause overflow. */

	*(unsigned char *)normin = 'Y';
	if (scale != 1.f) {
	    ix = icamax_(n, &WORK(1), &c__1);
	    i__1 = ix;
	    if (scale < ((r__1 = WORK(ix).r, dabs(r__1)) + (r__2 = r_imag(&
		    WORK(ix)), dabs(r__2))) * smlnum || scale == 0.f) {
		goto L40;
	    }
	    csrscl_(n, &scale, &WORK(1), &c__1);
	}
	goto L10;
    }

/*     Compute the estimate of the reciprocal condition number. */

    if (ainvnm != 0.f) {
	*rcond = 1.f / ainvnm / *anorm;
    }

L40:
    return 0;

/*     End of CGBCON */

} /* cgbcon_ */

