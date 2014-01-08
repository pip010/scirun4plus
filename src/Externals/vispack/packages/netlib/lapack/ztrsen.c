#include "f2c.h"

/* Subroutine */ int ztrsen_(char *job, char *compq, logical *select, integer 
	*n, doublecomplex *t, integer *ldt, doublecomplex *q, integer *ldq, 
	doublecomplex *w, integer *m, doublereal *s, doublereal *sep, 
	doublecomplex *work, integer *lwork, integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       March 31, 1993   


    Purpose   
    =======   

    ZTRSEN reorders the Schur factorization of a complex matrix   
    A = Q*T*Q**H, so that a selected cluster of eigenvalues appears in   
    the leading positions on the diagonal of the upper triangular matrix 
  
    T, and the leading columns of Q form an orthonormal basis of the   
    corresponding right invariant subspace.   

    Optionally the routine computes the reciprocal condition numbers of   
    the cluster of eigenvalues and/or the invariant subspace.   

    Arguments   
    =========   

    JOB     (input) CHARACTER*1   
            Specifies whether condition numbers are required for the   
            cluster of eigenvalues (S) or the invariant subspace (SEP):   
            = 'N': none;   
            = 'E': for eigenvalues only (S);   
            = 'V': for invariant subspace only (SEP);   
            = 'B': for both eigenvalues and invariant subspace (S and   
                   SEP).   

    COMPQ   (input) CHARACTER*1   
            = 'V': update the matrix Q of Schur vectors;   
            = 'N': do not update Q.   

    SELECT  (input) LOGICAL array, dimension (N)   
            SELECT specifies the eigenvalues in the selected cluster. To 
  
            select the j-th eigenvalue, SELECT(j) must be set to .TRUE.. 
  

    N       (input) INTEGER   
            The order of the matrix T. N >= 0.   

    T       (input/output) COMPLEX*16 array, dimension (LDT,N)   
            On entry, the upper triangular matrix T.   
            On exit, T is overwritten by the reordered matrix T, with the 
  
            selected eigenvalues as the leading diagonal elements.   

    LDT     (input) INTEGER   
            The leading dimension of the array T. LDT >= max(1,N).   

    Q       (input/output) COMPLEX*16 array, dimension (LDQ,N)   
            On entry, if COMPQ = 'V', the matrix Q of Schur vectors.   
            On exit, if COMPQ = 'V', Q has been postmultiplied by the   
            unitary transformation matrix which reorders T; the leading M 
  
            columns of Q form an orthonormal basis for the specified   
            invariant subspace.   
            If COMPQ = 'N', Q is not referenced.   

    LDQ     (input) INTEGER   
            The leading dimension of the array Q.   
            LDQ >= 1; and if COMPQ = 'V', LDQ >= N.   

    W       (output) COMPLEX*16   
            The reordered eigenvalues of T, in the same order as they   
            appear on the diagonal of T.   

    M       (output) INTEGER   
            The dimension of the specified invariant subspace.   
            0 <= M <= N.   

    S       (output) DOUBLE PRECISION   
            If JOB = 'E' or 'B', S is a lower bound on the reciprocal   
            condition number for the selected cluster of eigenvalues.   
            S cannot underestimate the true reciprocal condition number   
            by more than a factor of sqrt(N). If M = 0 or N, S = 1.   
            If JOB = 'N' or 'V', S is not referenced.   

    SEP     (output) DOUBLE PRECISION   
            If JOB = 'V' or 'B', SEP is the estimated reciprocal   
            condition number of the specified invariant subspace. If   
            M = 0 or N, SEP = norm(T).   
            If JOB = 'N' or 'E', SEP is not referenced.   

    WORK    (workspace) COMPLEX*16 array, dimension (LWORK)   
            If JOB = 'N', WORK is not referenced.   

    LWORK   (input) INTEGER   
            The dimension of the array WORK.   
            If JOB = 'N', LWORK >= 1;   
            if JOB = 'E', LWORK = M*(N-M);   
            if JOB = 'V' or 'B', LWORK >= 2*M*(N-M).   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   

    Further Details   
    ===============   

    ZTRSEN first collects the selected eigenvalues by computing a unitary 
  
    transformation Z to move them to the top left corner of T. In other   
    words, the selected eigenvalues are the eigenvalues of T11 in:   

                  Z'*T*Z = ( T11 T12 ) n1   
                           (  0  T22 ) n2   
                              n1  n2   

    where N = n1+n2 and Z' means the conjugate transpose of Z. The first 
  
    n1 columns of Z span the specified invariant subspace of T.   

    If T has been obtained from the Schur factorization of a matrix   
    A = Q*T*Q', then the reordered Schur factorization of A is given by   
    A = (Q*Z)*(Z'*T*Z)*(Q*Z)', and the first n1 columns of Q*Z span the   
    corresponding invariant subspace of A.   

    The reciprocal condition number of the average of the eigenvalues of 
  
    T11 may be returned in S. S lies between 0 (very badly conditioned)   
    and 1 (very well conditioned). It is computed as follows. First we   
    compute R so that   

                           P = ( I  R ) n1   
                               ( 0  0 ) n2   
                                 n1 n2   

    is the projector on the invariant subspace associated with T11.   
    R is the solution of the Sylvester equation:   

                          T11*R - R*T22 = T12.   

    Let F-norm(M) denote the Frobenius-norm of M and 2-norm(M) denote   
    the two-norm of M. Then S is computed as the lower bound   

                        (1 + F-norm(R)**2)**(-1/2)   

    on the reciprocal of 2-norm(P), the true reciprocal condition number. 
  
    S cannot underestimate 1 / 2-norm(P) by more than a factor of   
    sqrt(N).   

    An approximate error bound for the computed average of the   
    eigenvalues of T11 is   

                           EPS * norm(T) / S   

    where EPS is the machine precision.   

    The reciprocal condition number of the right invariant subspace   
    spanned by the first n1 columns of Z (or of Q*Z) is returned in SEP. 
  
    SEP is defined as the separation of T11 and T22:   

                       sep( T11, T22 ) = sigma-min( C )   

    where sigma-min(C) is the smallest singular value of the   
    n1*n2-by-n1*n2 matrix   

       C  = kprod( I(n2), T11 ) - kprod( transpose(T22), I(n1) )   

    I(m) is an m by m identity matrix, and kprod denotes the Kronecker   
    product. We estimate sigma-min(C) by the reciprocal of an estimate of 
  
    the 1-norm of inverse(C). The true reciprocal 1-norm of inverse(C)   
    cannot differ from sigma-min(C) by more than a factor of sqrt(n1*n2). 
  

    When SEP is small, small changes in T can cause large changes in   
    the invariant subspace. An approximate bound on the maximum angular   
    error in the computed right invariant subspace is   

                        EPS * norm(T) / SEP   

    ===================================================================== 
  


       Decode and test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static integer c_n1 = -1;
    
    /* System generated locals */
    integer q_dim1, q_offset, t_dim1, t_offset, i__1, i__2, i__3;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    static integer kase, ierr, k;
    static doublereal scale;
    extern logical lsame_(char *, char *);
    static logical wantq, wants;
    static doublereal rnorm;
    static integer n1, n2;
    static doublereal rwork[1];
    static integer nn, ks;
    extern /* Subroutine */ int xerbla_(char *, integer *);
    extern doublereal zlange_(char *, integer *, integer *, doublecomplex *, 
	    integer *, doublereal *);
    static logical wantbh;
    extern /* Subroutine */ int zlacon_(integer *, doublecomplex *, 
	    doublecomplex *, doublereal *, integer *), zlacpy_(char *, 
	    integer *, integer *, doublecomplex *, integer *, doublecomplex *,
	     integer *);
    static logical wantsp;
    extern /* Subroutine */ int ztrexc_(char *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, integer *, integer *, integer *, 
	    integer *), ztrsyl_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, integer *,
	     doublecomplex *, integer *, doublereal *, integer *);
    static doublereal est;



#define RWORK(I) rwork[(I)]
#define SELECT(I) select[(I)-1]
#define W(I) w[(I)-1]
#define WORK(I) work[(I)-1]

#define T(I,J) t[(I)-1 + ((J)-1)* ( *ldt)]
#define Q(I,J) q[(I)-1 + ((J)-1)* ( *ldq)]

    wantbh = lsame_(job, "B");
    wants = lsame_(job, "E") || wantbh;
    wantsp = lsame_(job, "V") || wantbh;
    wantq = lsame_(compq, "V");

/*     Set M to the number of selected eigenvalues. */

    *m = 0;
    i__1 = *n;
    for (k = 1; k <= *n; ++k) {
	if (SELECT(k)) {
	    ++(*m);
	}
/* L10: */
    }

    n1 = *m;
    n2 = *n - *m;
    nn = n1 * n2;

    *info = 0;
    if (! lsame_(job, "N") && ! wants && ! wantsp) {
	*info = -1;
    } else if (! lsame_(compq, "N") && ! wantq) {
	*info = -2;
    } else if (*n < 0) {
	*info = -4;
    } else if (*ldt < max(1,*n)) {
	*info = -6;
    } else if (*ldq < 1 || wantq && *ldq < *n) {
	*info = -8;
    } else if (*lwork < 1 || wants && ! wantsp && *lwork < nn || wantsp && *
	    lwork < nn << 1) {
	*info = -14;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("ZTRSEN", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*m == *n || *m == 0) {
	if (wants) {
	    *s = 1.;
	}
	if (wantsp) {
	    *sep = zlange_("1", n, n, &T(1,1), ldt, rwork);
	}
	goto L40;
    }

/*     Collect the selected eigenvalues at the top left corner of T. */

    ks = 0;
    i__1 = *n;
    for (k = 1; k <= *n; ++k) {
	if (SELECT(k)) {
	    ++ks;

/*           Swap the K-th eigenvalue to position KS. */

	    if (k != ks) {
		ztrexc_(compq, n, &T(1,1), ldt, &Q(1,1), ldq, &k, &
			ks, &ierr);
	    }
	}
/* L20: */
    }

    if (wants) {

/*        Solve the Sylvester equation for R:   

             T11*R - R*T22 = scale*T12 */

	zlacpy_("F", &n1, &n2, &T(1,n1+1), ldt, &WORK(1), &n1);
	ztrsyl_("N", "N", &c_n1, &n1, &n2, &T(1,1), ldt, &T(n1+1,n1+1), ldt, &WORK(1), &n1, &scale, &ierr);

/*        Estimate the reciprocal of the condition number of the clust
er   
          of eigenvalues. */

	rnorm = zlange_("F", &n1, &n2, &WORK(1), &n1, rwork);
	if (rnorm == 0.) {
	    *s = 1.;
	} else {
	    *s = scale / (sqrt(scale * scale / rnorm + rnorm) * sqrt(rnorm));
	}
    }

    if (wantsp) {

/*        Estimate sep(T11,T22). */

	est = 0.;
	kase = 0;
L30:
	zlacon_(&nn, &WORK(nn + 1), &WORK(1), &est, &kase);
	if (kase != 0) {
	    if (kase == 1) {

/*              Solve T11*R - R*T22 = scale*X. */

		ztrsyl_("N", "N", &c_n1, &n1, &n2, &T(1,1), ldt, &T(n1+1,n1+1), ldt, &WORK(1), &n1, &scale, &
			ierr);
	    } else {

/*              Solve T11'*R - R*T22' = scale*X. */

		ztrsyl_("C", "C", &c_n1, &n1, &n2, &T(1,1), ldt, &T(n1+1,n1+1), ldt, &WORK(1), &n1, &scale, &
			ierr);
	    }
	    goto L30;
	}

	*sep = scale / est;
    }

L40:

/*     Copy reordered eigenvalues to W. */

    i__1 = *n;
    for (k = 1; k <= *n; ++k) {
	i__2 = k;
	i__3 = k + k * t_dim1;
	W(k).r = T(k,k).r, W(k).i = T(k,k).i;
/* L50: */
    }
    return 0;

/*     End of ZTRSEN */

} /* ztrsen_ */

