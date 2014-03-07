#include "f2c.h"

/* Subroutine */ int sggbal_(char *job, integer *n, real *a, integer *lda, 
	real *b, integer *ldb, integer *ilo, integer *ihi, real *lscale, real 
	*rscale, real *work, integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    SGGBAL balances a pair of general real matrices (A,B).  This   
    involves, first, permuting A and B by similarity transformations to   
    isolate eigenvalues in the first 1 to ILO$-$1 and last IHI+1 to N   
    elements on the diagonal; and second, applying a diagonal similarity 
  
    transformation to rows and columns ILO to IHI to make the rows   
    and columns as close in norm as possible. Both steps are optional.   

    Balancing may reduce the 1-norm of the matrices, and improve the   
    accuracy of the computed eigenvalues and/or eigenvectors in the   
    generalized eigenvalue problem A*x = lambda*B*x.   

    Arguments   
    =========   

    JOB     (input) CHARACTER*1   
            Specifies the operations to be performed on A and B:   
            = 'N':  none:  simply set ILO = 1, IHI = N, LSCALE(I) = 1.0   
                    and RSCALE(I) = 1.0 for i = 1,...,N.   
            = 'P':  permute only;   
            = 'S':  scale only;   
            = 'B':  both permute and scale.   

    N       (input) INTEGER   
            The order of the matrices A and B.  N >= 0.   

    A       (input/output) REAL array, dimension (LDA,N)   
            On entry, the input matrix A.   
            On exit,  A is overwritten by the balanced matrix.   
            If JOB = 'N', A is not referenced.   

    LDA     (input) INTEGER   
            The leading dimension of the array A. LDA >= max(1,N).   

    B       (input/output) REAL array, dimension (LDB,N)   
            On entry, the input matrix B.   
            On exit,  B is overwritten by the balanced matrix.   
            If JOB = 'N', B is not referenced.   

    LDB     (input) INTEGER   
            The leading dimension of the array B. LDB >= max(1,N).   

    ILO     (output) INTEGER   
    IHI     (output) INTEGER   
            ILO and IHI are set to integers such that on exit   
            A(i,j) = 0 and B(i,j) = 0 if i > j and   
            j = 1,...,ILO-1 or i = IHI+1,...,N.   
            If JOB = 'N' or 'S', ILO = 1 and IHI = N.   

    LSCALE  (output) REAL array, dimension (N)   
            Details of the permutations and scaling factors applied   
            to the left side of A and B.  If P(j) is the index of the   
            row interchanged with row j, and D(j)   
            is the scaling factor applied to row j, then   
              LSCALE(j) = P(j)    for J = 1,...,ILO-1   
                        = D(j)    for J = ILO,...,IHI   
                        = P(j)    for J = IHI+1,...,N.   
            The order in which the interchanges are made is N to IHI+1,   
            then 1 to ILO-1.   

    RSCALE  (output) REAL array, dimension (N)   
            Details of the permutations and scaling factors applied   
            to the right side of A and B.  If P(j) is the index of the   
            column interchanged with column j, and D(j)   
            is the scaling factor applied to column j, then   
              LSCALE(j) = P(j)    for J = 1,...,ILO-1   
                        = D(j)    for J = ILO,...,IHI   
                        = P(j)    for J = IHI+1,...,N.   
            The order in which the interchanges are made is N to IHI+1,   
            then 1 to ILO-1.   

    WORK    (workspace) REAL array, dimension (6*N)   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value.   

    Further Details   
    ===============   

    See R.C. WARD, Balancing the generalized eigenvalue problem,   
                   SIAM J. Sci. Stat. Comp. 2 (1981), 141-152.   

    ===================================================================== 
  


       Test the input parameters   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static integer c__1 = 1;
    static real c_b34 = 10.f;
    static real c_b70 = .5f;
    
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3;
    real r__1, r__2, r__3;
    /* Builtin functions */
    double r_lg10(real *), r_sign(real *, real *), pow_ri(real *, integer *);
    /* Local variables */
    static integer lcab;
    static real beta, coef;
    static integer irab, lrab;
    static real basl, cmax;
    extern doublereal sdot_(integer *, real *, integer *, real *, integer *);
    static real coef2, coef5;
    static integer i, j, k, l, m;
    static real gamma, t, alpha;
    extern logical lsame_(char *, char *);
    extern /* Subroutine */ int sscal_(integer *, real *, real *, integer *);
    static real sfmin, sfmax;
    static integer iflow;
    extern /* Subroutine */ int sswap_(integer *, real *, integer *, real *, 
	    integer *);
    static integer kount;
    extern /* Subroutine */ int saxpy_(integer *, real *, real *, integer *, 
	    real *, integer *);
    static integer jc;
    static real ta, tb, tc;
    static integer ir, it;
    static real ew;
    static integer nr;
    static real pgamma;
    extern doublereal slamch_(char *);
    extern /* Subroutine */ int xerbla_(char *, integer *);
    extern integer isamax_(integer *, real *, integer *);
    static integer lsfmin, lsfmax, ip1, jp1, lm1;
    static real cab, rab, ewc, cor, sum;
    static integer nrp2, icab;



#define LSCALE(I) lscale[(I)-1]
#define RSCALE(I) rscale[(I)-1]
#define WORK(I) work[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]
#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]

    *info = 0;
    if (! lsame_(job, "N") && ! lsame_(job, "P") && ! lsame_(
	    job, "S") && ! lsame_(job, "B")) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*lda < max(1,*n)) {
	*info = -4;
    } else if (*ldb < max(1,*n)) {
	*info = -5;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("SGGBAL", &i__1);
	return 0;
    }

    k = 1;
    l = *n;

/*     Quick return if possible */

    if (*n == 0) {
	return 0;
    }

    if (lsame_(job, "N")) {
	*ilo = 1;
	*ihi = *n;
	i__1 = *n;
	for (i = 1; i <= *n; ++i) {
	    LSCALE(i) = 1.f;
	    RSCALE(i) = 1.f;
/* L10: */
	}
	return 0;
    }

    if (k == l) {
	*ilo = 1;
	*ihi = 1;
	LSCALE(1) = 1.f;
	RSCALE(1) = 1.f;
	return 0;
    }

    if (lsame_(job, "S")) {
	goto L190;
    }

    goto L30;

/*     Permute the matrices A and B to isolate the eigenvalues.   

       Find row with one nonzero in columns 1 through L */

L20:
    l = lm1;
    if (l != 1) {
	goto L30;
    }

    RSCALE(1) = 1.f;
    LSCALE(1) = 1.f;
    goto L190;

L30:
    lm1 = l - 1;
    for (i = l; i >= 1; --i) {
	i__1 = lm1;
	for (j = 1; j <= lm1; ++j) {
	    jp1 = j + 1;
	    if (A(i,j) != 0.f || B(i,j) != 0.f) {
		goto L50;
	    }
/* L40: */
	}
	j = l;
	goto L70;

L50:
	i__1 = l;
	for (j = jp1; j <= l; ++j) {
	    if (A(i,j) != 0.f || B(i,j) != 0.f) {
		goto L80;
	    }
/* L60: */
	}
	j = jp1 - 1;

L70:
	m = l;
	iflow = 1;
	goto L160;
L80:
	;
    }
    goto L100;

/*     Find column with one nonzero in rows K through N */

L90:
    ++k;

L100:
    i__1 = l;
    for (j = k; j <= l; ++j) {
	i__2 = lm1;
	for (i = k; i <= lm1; ++i) {
	    ip1 = i + 1;
	    if (A(i,j) != 0.f || B(i,j) != 0.f) {
		goto L120;
	    }
/* L110: */
	}
	i = l;
	goto L140;
L120:
	i__2 = l;
	for (i = ip1; i <= l; ++i) {
	    if (A(i,j) != 0.f || B(i,j) != 0.f) {
		goto L150;
	    }
/* L130: */
	}
	i = ip1 - 1;
L140:
	m = k;
	iflow = 2;
	goto L160;
L150:
	;
    }
    goto L190;

/*     Permute rows M and I */

L160:
    LSCALE(m) = (real) i;
    if (i == m) {
	goto L170;
    }
    i__1 = *n - k + 1;
    sswap_(&i__1, &A(i,k), lda, &A(m,k), lda);
    i__1 = *n - k + 1;
    sswap_(&i__1, &B(i,k), ldb, &B(m,k), ldb);

/*     Permute columns M and J */

L170:
    RSCALE(m) = (real) j;
    if (j == m) {
	goto L180;
    }
    sswap_(&l, &A(1,j), &c__1, &A(1,m), &c__1);
    sswap_(&l, &B(1,j), &c__1, &B(1,m), &c__1);

L180:
    switch (iflow) {
	case 1:  goto L20;
	case 2:  goto L90;
    }

L190:
    *ilo = k;
    *ihi = l;

    if (*ilo == *ihi) {
	return 0;
    }

    if (lsame_(job, "P")) {
	return 0;
    }

/*     Balance the submatrix in rows ILO to IHI. */

    nr = *ihi - *ilo + 1;
    i__1 = *ihi;
    for (i = *ilo; i <= *ihi; ++i) {
	RSCALE(i) = 0.f;
	LSCALE(i) = 0.f;

	WORK(i) = 0.f;
	WORK(i + *n) = 0.f;
	WORK(i + (*n << 1)) = 0.f;
	WORK(i + *n * 3) = 0.f;
	WORK(i + (*n << 2)) = 0.f;
	WORK(i + *n * 5) = 0.f;
/* L200: */
    }

/*     Compute right side vector in resulting linear equations */

    basl = r_lg10(&c_b34);
    i__1 = *ihi;
    for (i = *ilo; i <= *ihi; ++i) {
	i__2 = *ihi;
	for (j = *ilo; j <= *ihi; ++j) {
	    tb = B(i,j);
	    ta = A(i,j);
	    if (ta == 0.f) {
		goto L210;
	    }
	    r__1 = dabs(ta);
	    ta = r_lg10(&r__1) / basl;
L210:
	    if (tb == 0.f) {
		goto L220;
	    }
	    r__1 = dabs(tb);
	    tb = r_lg10(&r__1) / basl;
L220:
	    WORK(i + (*n << 2)) = WORK(i + (*n << 2)) - ta - tb;
	    WORK(j + *n * 5) = WORK(j + *n * 5) - ta - tb;
/* L230: */
	}
/* L240: */
    }

    coef = 1.f / (real) (nr << 1);
    coef2 = coef * coef;
    coef5 = coef2 * .5f;
    nrp2 = nr + 2;
    beta = 0.f;
    it = 1;

/*     Start generalized conjugate gradient iteration */

L250:

    gamma = sdot_(&nr, &WORK(*ilo + (*n << 2)), &c__1, &WORK(*ilo + (*n << 2))
	    , &c__1) + sdot_(&nr, &WORK(*ilo + *n * 5), &c__1, &WORK(*ilo + *
	    n * 5), &c__1);

    ew = 0.f;
    ewc = 0.f;
    i__1 = *ihi;
    for (i = *ilo; i <= *ihi; ++i) {
	ew += WORK(i + (*n << 2));
	ewc += WORK(i + *n * 5);
/* L260: */
    }

/* Computing 2nd power */
    r__1 = ew;
/* Computing 2nd power */
    r__2 = ewc;
/* Computing 2nd power */
    r__3 = ew - ewc;
    gamma = coef * gamma - coef2 * (r__1 * r__1 + r__2 * r__2) - coef5 * (
	    r__3 * r__3);
    if (gamma == 0.f) {
	goto L350;
    }
    if (it != 1) {
	beta = gamma / pgamma;
    }
    t = coef5 * (ewc - ew * 3.f);
    tc = coef5 * (ew - ewc * 3.f);

    sscal_(&nr, &beta, &WORK(*ilo), &c__1);
    sscal_(&nr, &beta, &WORK(*ilo + *n), &c__1);

    saxpy_(&nr, &coef, &WORK(*ilo + (*n << 2)), &c__1, &WORK(*ilo + *n), &
	    c__1);
    saxpy_(&nr, &coef, &WORK(*ilo + *n * 5), &c__1, &WORK(*ilo), &c__1);

    i__1 = *ihi;
    for (i = *ilo; i <= *ihi; ++i) {
	WORK(i) += tc;
	WORK(i + *n) += t;
/* L270: */
    }

/*     Apply matrix to vector */

    i__1 = *ihi;
    for (i = *ilo; i <= *ihi; ++i) {
	kount = 0;
	sum = 0.f;
	i__2 = *ihi;
	for (j = *ilo; j <= *ihi; ++j) {
	    if (A(i,j) == 0.f) {
		goto L280;
	    }
	    ++kount;
	    sum += WORK(j);
L280:
	    if (B(i,j) == 0.f) {
		goto L290;
	    }
	    ++kount;
	    sum += WORK(j);
L290:
	    ;
	}
	WORK(i + (*n << 1)) = (real) kount * WORK(i + *n) + sum;
/* L300: */
    }

    i__1 = *ihi;
    for (j = *ilo; j <= *ihi; ++j) {
	kount = 0;
	sum = 0.f;
	i__2 = *ihi;
	for (i = *ilo; i <= *ihi; ++i) {
	    if (A(i,j) == 0.f) {
		goto L310;
	    }
	    ++kount;
	    sum += WORK(i + *n);
L310:
	    if (B(i,j) == 0.f) {
		goto L320;
	    }
	    ++kount;
	    sum += WORK(i + *n);
L320:
	    ;
	}
	WORK(j + *n * 3) = (real) kount * WORK(j) + sum;
/* L330: */
    }

    sum = sdot_(&nr, &WORK(*ilo + *n), &c__1, &WORK(*ilo + (*n << 1)), &c__1) 
	    + sdot_(&nr, &WORK(*ilo), &c__1, &WORK(*ilo + *n * 3), &c__1);
    alpha = gamma / sum;

/*     Determine correction to current iteration */

    cmax = 0.f;
    i__1 = *ihi;
    for (i = *ilo; i <= *ihi; ++i) {
	cor = alpha * WORK(i + *n);
	if (dabs(cor) > cmax) {
	    cmax = dabs(cor);
	}
	LSCALE(i) += cor;
	cor = alpha * WORK(i);
	if (dabs(cor) > cmax) {
	    cmax = dabs(cor);
	}
	RSCALE(i) += cor;
/* L340: */
    }
    if (cmax < .5f) {
	goto L350;
    }

    r__1 = -(doublereal)alpha;
    saxpy_(&nr, &r__1, &WORK(*ilo + (*n << 1)), &c__1, &WORK(*ilo + (*n << 2))
	    , &c__1);
    r__1 = -(doublereal)alpha;
    saxpy_(&nr, &r__1, &WORK(*ilo + *n * 3), &c__1, &WORK(*ilo + *n * 5), &
	    c__1);

    pgamma = gamma;
    ++it;
    if (it <= nrp2) {
	goto L250;
    }

/*     End generalized conjugate gradient iteration */

L350:
    sfmin = slamch_("S");
    sfmax = 1.f / sfmin;
    lsfmin = (integer) (r_lg10(&sfmin) / basl + 1.f);
    lsfmax = (integer) (r_lg10(&sfmax) / basl);
    i__1 = *ihi;
    for (i = *ilo; i <= *ihi; ++i) {
	i__2 = *n - *ilo + 1;
	irab = isamax_(&i__2, &A(i,*ilo), lda);
	rab = (r__1 = A(i,irab+*ilo-1), dabs(r__1));
	i__2 = *n - *ilo + 1;
	irab = isamax_(&i__2, &B(i,*ilo), lda);
/* Computing MAX */
	r__2 = rab, r__3 = (r__1 = B(i,irab+*ilo-1), dabs(
		r__1));
	rab = dmax(r__2,r__3);
	r__1 = rab + sfmin;
	lrab = (integer) (r_lg10(&r__1) / basl + 1.f);
	ir = LSCALE(i) + r_sign(&c_b70, &LSCALE(i));
/* Computing MIN */
	i__2 = max(ir,lsfmin), i__2 = min(i__2,lsfmax), i__3 = lsfmax - lrab;
	ir = min(i__2,i__3);
	LSCALE(i) = pow_ri(&c_b34, &ir);
	icab = isamax_(ihi, &A(1,i), &c__1);
	cab = (r__1 = A(icab,i), dabs(r__1));
	icab = isamax_(ihi, &B(1,i), &c__1);
/* Computing MAX */
	r__2 = cab, r__3 = (r__1 = B(icab,i), dabs(r__1));
	cab = dmax(r__2,r__3);
	r__1 = cab + sfmin;
	lcab = (integer) (r_lg10(&r__1) / basl + 1.f);
	jc = RSCALE(i) + r_sign(&c_b70, &RSCALE(i));
/* Computing MIN */
	i__2 = max(jc,lsfmin), i__2 = min(i__2,lsfmax), i__3 = lsfmax - lcab;
	jc = min(i__2,i__3);
	RSCALE(i) = pow_ri(&c_b34, &jc);
/* L360: */
    }

/*     Row scaling of matrices A and B */

    i__1 = *ihi;
    for (i = *ilo; i <= *ihi; ++i) {
	i__2 = *n - *ilo + 1;
	sscal_(&i__2, &LSCALE(i), &A(i,*ilo), lda);
	i__2 = *n - *ilo + 1;
	sscal_(&i__2, &LSCALE(i), &B(i,*ilo), ldb);
/* L370: */
    }

/*     Column scaling of matrices A and B */

    i__1 = *ihi;
    for (j = *ilo; j <= *ihi; ++j) {
	sscal_(ihi, &RSCALE(j), &A(1,j), &c__1);
	sscal_(ihi, &RSCALE(j), &B(1,j), &c__1);
/* L380: */
    }

    return 0;

/*     End of SGGBAL */

} /* sggbal_ */

