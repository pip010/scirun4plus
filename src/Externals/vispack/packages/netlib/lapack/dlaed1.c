#include "f2c.h"

/* Subroutine */ int dlaed1_(integer *n, doublereal *d, doublereal *q, 
	integer *ldq, integer *indxq, doublereal *rho, integer *cutpnt, 
	doublereal *work, integer *iwork, integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DLAED1 computes the updated eigensystem of a diagonal   
    matrix after modification by a rank-one symmetric matrix.  This   
    routine is used only for the eigenproblem which requires all   
    eigenvalues and eigenvectors of a tridiagonal matrix.  DLAED7 handles 
  
    the case in which eigenvalues only or eigenvalues and eigenvectors   
    of a full symmetric matrix (which was reduced to tridiagonal form)   
    are desired.   

      T = Q(in) ( D(in) + RHO * Z*Z' ) Q'(in) = Q(out) * D(out) * Q'(out) 
  

       where Z = Q'u, u is a vector of length N with ones in the   
       CUTPNT and CUTPNT + 1 th elements and zeros elsewhere.   

       The eigenvectors of the original matrix are stored in Q, and the   
       eigenvalues are in D.  The algorithm consists of three stages:   

          The first stage consists of deflating the size of the problem   
          when there are multiple eigenvalues or if there is a zero in   
          the Z vector.  For each such occurence the dimension of the   
          secular equation problem is reduced by one.  This stage is   
          performed by the routine DLAED2.   

          The second stage consists of calculating the updated   
          eigenvalues. This is done by finding the roots of the secular   
          equation via the routine DLAED4 (as called by SLAED3).   
          This routine also calculates the eigenvectors of the current   
          problem.   

          The final stage consists of computing the updated eigenvectors 
  
          directly using the updated eigenvalues.  The eigenvectors for   
          the current problem are multiplied with the eigenvectors from   
          the overall problem.   

    Arguments   
    =========   

    N      (input) INTEGER   
           The dimension of the symmetric tridiagonal matrix.  N >= 0.   

    D      (input/output) DOUBLE PRECISION array, dimension (N)   
           On entry, the eigenvalues of the rank-1-perturbed matrix.   
           On exit, the eigenvalues of the repaired matrix.   

    Q      (input/output) DOUBLE PRECISION array, dimension (LDQ,N)   
           On entry, the eigenvectors of the rank-1-perturbed matrix.   
           On exit, the eigenvectors of the repaired tridiagonal matrix. 
  

    LDQ    (input) INTEGER   
           The leading dimension of the array Q.  LDQ >= max(1,N).   

    INDXQ  (input/output) INTEGER array, dimension (N)   
           On entry, the permutation which separately sorts the two   
           subproblems in D into ascending order.   
           On exit, the permutation which will reintegrate the   
           subproblems back into sorted order,   
           i.e. D( INDXQ( I = 1, N ) ) will be in ascending order.   

    RHO    (input) DOUBLE PRECISION   
           The subdiagonal entry used to create the rank-1 modification. 
  

    CUTPNT (input) INTEGER   
           The location of the last eigenvalue in the leading sub-matrix. 
  
           min(1,N) <= CUTPNT <= N.   

    WORK   (workspace) DOUBLE PRECISION array, dimension (3*N+2*N**2)   

    IWORK  (workspace) INTEGER array, dimension (4*N)   

    INFO   (output) INTEGER   
            = 0:  successful exit.   
            < 0:  if INFO = -i, the i-th argument had an illegal value.   
            > 0:  if INFO = 1, an eigenvalue did not converge   

    ===================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static integer c__1 = 1;
    static integer c_n1 = -1;
    
    /* System generated locals */
    integer q_dim1, q_offset, i__1;
    /* Local variables */
    static integer indx, i, k, indxc;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static integer indxp;
    extern /* Subroutine */ int dlaed2_(integer *, integer *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    doublereal *, integer *, integer *, integer *, integer *), 
	    dlaed3_(integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, integer *);
    static integer n1, n2, idlmda, is, iw, iz;
    extern /* Subroutine */ int dlamrg_(integer *, integer *, doublereal *, 
	    integer *, integer *, integer *), xerbla_(char *, integer *);
    static integer coltyp, iq2, ldq2, zpp1;



#define D(I) d[(I)-1]
#define INDXQ(I) indxq[(I)-1]
#define WORK(I) work[(I)-1]
#define IWORK(I) iwork[(I)-1]

#define Q(I,J) q[(I)-1 + ((J)-1)* ( *ldq)]

    *info = 0;

    if (*n < 0) {
	*info = -1;
    } else if (*ldq < max(1,*n)) {
	*info = -4;
    } else if (min(1,*n) > *cutpnt || *n < *cutpnt) {
	*info = -7;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DLAED1", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return 0;
    }

/*     The following values are for bookkeeping purposes only.  They are 
  
       integer pointers which indicate the portion of the workspace   
       used by a particular array in DLAED2 and SLAED3. */

    ldq2 = *n;

    iz = 1;
    idlmda = iz + *n;
    iw = idlmda + *n;
    iq2 = iw + *n;
    is = iq2 + *n * ldq2;

    indx = 1;
    indxc = indx + *n;
    coltyp = indxc + *n;
    indxp = coltyp + *n;


/*     Form the z-vector which consists of the last row of Q_1 and the   
       first row of Q_2. */

    dcopy_(cutpnt, &Q(*cutpnt,1), ldq, &WORK(iz), &c__1);
    zpp1 = *cutpnt + 1;
    i__1 = *n - *cutpnt;
    dcopy_(&i__1, &Q(zpp1,zpp1), ldq, &WORK(iz + *cutpnt), &c__1);

/*     Deflate eigenvalues. */

    dlaed2_(&k, n, &D(1), &Q(1,1), ldq, &INDXQ(1), rho, cutpnt, &WORK(iz)
	    , &WORK(idlmda), &WORK(iq2), &ldq2, &IWORK(indxc), &WORK(iw), &
	    IWORK(indxp), &IWORK(indx), &IWORK(coltyp), info);
    if (*info != 0) {
	goto L20;
    }

/*     Solve Secular Equation. */

    if (k != 0) {
	dlaed3_(&k, &c__1, &k, n, &D(1), &Q(1,1), ldq, rho, cutpnt, &
		WORK(idlmda), &WORK(iq2), &ldq2, &IWORK(indxc), &IWORK(coltyp)
		, &WORK(iw), &WORK(is), &k, info);
	if (*info != 0) {
	    goto L20;
	}

/*     Prepare the INDXQ sorting permutation. */

	n1 = k;
	n2 = *n - k;
	dlamrg_(&n1, &n2, &D(1), &c__1, &c_n1, &INDXQ(1));
    } else {
	i__1 = *n;
	for (i = 1; i <= *n; ++i) {
	    INDXQ(i) = i;
/* L10: */
	}
    }

L20:
    return 0;

/*     End of DLAED1 */

} /* dlaed1_ */

