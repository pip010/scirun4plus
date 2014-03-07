#include "f2c.h"

/* Subroutine */ int zhpev_(char *jobz, char *uplo, integer *n, doublecomplex 
	*ap, doublereal *w, doublecomplex *z, integer *ldz, doublecomplex *
	work, doublereal *rwork, integer *info)
{
/*  -- LAPACK driver routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       March 31, 1993   


    Purpose   
    =======   

    ZHPEV computes all the eigenvalues and, optionally, eigenvectors of a 
  
    complex Hermitian matrix in packed storage.   

    Arguments   
    =========   

    JOBZ    (input) CHARACTER*1   
            = 'N':  Compute eigenvalues only;   
            = 'V':  Compute eigenvalues and eigenvectors.   

    UPLO    (input) CHARACTER*1   
            = 'U':  Upper triangle of A is stored;   
            = 'L':  Lower triangle of A is stored.   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    AP      (input/output) COMPLEX*16 array, dimension (N*(N+1)/2)   
            On entry, the upper or lower triangle of the Hermitian matrix 
  
            A, packed columnwise in a linear array.  The j-th column of A 
  
            is stored in the array AP as follows:   
            if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;   
            if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n. 
  

            On exit, AP is overwritten by values generated during the   
            reduction to tridiagonal form.  If UPLO = 'U', the diagonal   
            and first superdiagonal of the tridiagonal matrix T overwrite 
  
            the corresponding elements of A, and if UPLO = 'L', the   
            diagonal and first subdiagonal of T overwrite the   
            corresponding elements of A.   

    W       (output) DOUBLE PRECISION array, dimension (N)   
            If INFO = 0, the eigenvalues in ascending order.   

    Z       (output) COMPLEX*16 array, dimension (LDZ, N)   
            If JOBZ = 'V', then if INFO = 0, Z contains the orthonormal   
            eigenvectors of the matrix A, with the i-th column of Z   
            holding the eigenvector associated with W(i).   
            If JOBZ = 'N', then Z is not referenced.   

    LDZ     (input) INTEGER   
            The leading dimension of the array Z.  LDZ >= 1, and if   
            JOBZ = 'V', LDZ >= max(1,N).   

    WORK    (workspace) COMPLEX*16 array, dimension (max(1, 2*N-1))   

    RWORK   (workspace) DOUBLE PRECISION array, dimension (max(1, 3*N-2)) 
  

    INFO    (output) INTEGER   
            = 0:  successful exit.   
            < 0:  if INFO = -i, the i-th argument had an illegal value.   
            > 0:  if INFO = i, the algorithm failed to converge; i   
                  off-diagonal elements of an intermediate tridiagonal   
                  form did not converge to zero.   

    ===================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static integer c__1 = 1;
    
    /* System generated locals */
    integer z_dim1, z_offset, i__1;
    doublereal d__1;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    static integer inde;
    static doublereal anrm;
    static integer imax;
    static doublereal rmin, rmax;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static doublereal sigma;
    extern logical lsame_(char *, char *);
    static integer iinfo;
    static logical wantz;
    extern doublereal dlamch_(char *);
    static integer iscale;
    static doublereal safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *), zdscal_(
	    integer *, doublereal *, doublecomplex *, integer *);
    static doublereal bignum;
    static integer indtau;
    extern /* Subroutine */ int dsterf_(integer *, doublereal *, doublereal *,
	     integer *);
    extern doublereal zlanhp_(char *, char *, integer *, doublecomplex *, 
	    doublereal *);
    static integer indrwk, indwrk;
    static doublereal smlnum;
    extern /* Subroutine */ int zhptrd_(char *, integer *, doublecomplex *, 
	    doublereal *, doublereal *, doublecomplex *, integer *), 
	    zsteqr_(char *, integer *, doublereal *, doublereal *, 
	    doublecomplex *, integer *, doublereal *, integer *), 
	    zupgtr_(char *, integer *, doublecomplex *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    static doublereal eps;



#define AP(I) ap[(I)-1]
#define W(I) w[(I)-1]
#define WORK(I) work[(I)-1]
#define RWORK(I) rwork[(I)-1]

#define Z(I,J) z[(I)-1 + ((J)-1)* ( *ldz)]

    wantz = lsame_(jobz, "V");

    *info = 0;
    if (! (wantz || lsame_(jobz, "N"))) {
	*info = -1;
    } else if (! (lsame_(uplo, "L") || lsame_(uplo, "U"))) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*ldz < 1 || wantz && *ldz < *n) {
	*info = -7;
    }

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("ZHPEV ", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return 0;
    }

    if (*n == 1) {
	W(1) = AP(1).r;
	RWORK(1) = 1.;
	if (wantz) {
	    i__1 = z_dim1 + 1;
	    Z(1,1).r = 1., Z(1,1).i = 0.;
	}
	return 0;
    }

/*     Get machine constants. */

    safmin = dlamch_("Safe minimum");
    eps = dlamch_("Precision");
    smlnum = safmin / eps;
    bignum = 1. / smlnum;
    rmin = sqrt(smlnum);
    rmax = sqrt(bignum);

/*     Scale matrix to allowable range, if necessary. */

    anrm = zlanhp_("M", uplo, n, &AP(1), &RWORK(1));
    iscale = 0;
    if (anrm > 0. && anrm < rmin) {
	iscale = 1;
	sigma = rmin / anrm;
    } else if (anrm > rmax) {
	iscale = 1;
	sigma = rmax / anrm;
    }
    if (iscale == 1) {
	i__1 = *n * (*n + 1) / 2;
	zdscal_(&i__1, &sigma, &AP(1), &c__1);
    }

/*     Call ZHPTRD to reduce Hermitian packed matrix to tridiagonal form. 
*/

    inde = 1;
    indtau = 1;
    zhptrd_(uplo, n, &AP(1), &W(1), &RWORK(inde), &WORK(indtau), &iinfo);

/*     For eigenvalues only, call DSTERF.  For eigenvectors, first call   
       ZUPGTR to generate the orthogonal matrix, then call ZSTEQR. */

    if (! wantz) {
	dsterf_(n, &W(1), &RWORK(inde), info);
    } else {
	indwrk = indtau + *n;
	zupgtr_(uplo, n, &AP(1), &WORK(indtau), &Z(1,1), ldz, &WORK(
		indwrk), &iinfo);
	indrwk = inde + *n;
	zsteqr_(jobz, n, &W(1), &RWORK(inde), &Z(1,1), ldz, &RWORK(
		indrwk), info);
    }

/*     If matrix was scaled, then rescale eigenvalues appropriately. */

    if (iscale == 1) {
	if (*info == 0) {
	    imax = *n;
	} else {
	    imax = *info - 1;
	}
	d__1 = 1. / sigma;
	dscal_(&imax, &d__1, &W(1), &c__1);
    }

    return 0;

/*     End of ZHPEV */

} /* zhpev_ */

