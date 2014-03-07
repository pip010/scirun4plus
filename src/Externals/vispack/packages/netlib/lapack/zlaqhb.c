#include "f2c.h"

/* Subroutine */ int zlaqhb_(char *uplo, integer *n, integer *kd, 
	doublecomplex *ab, integer *ldab, doublereal *s, doublereal *scond, 
	doublereal *amax, char *equed)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    ZLAQHB equilibrates a symmetric band matrix A using the scaling   
    factors in the vector S.   

    Arguments   
    =========   

    UPLO    (input) CHARACTER*1   
            Specifies whether the upper or lower triangular part of the   
            symmetric matrix A is stored.   
            = 'U':  Upper triangular   
            = 'L':  Lower triangular   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    KD      (input) INTEGER   
            The number of super-diagonals of the matrix A if UPLO = 'U', 
  
            or the number of sub-diagonals if UPLO = 'L'.  KD >= 0.   

    AB      (input/output) COMPLEX*16 array, dimension (LDAB,N)   
            On entry, the upper or lower triangle of the symmetric band   
            matrix A, stored in the first KD+1 rows of the array.  The   
            j-th column of A is stored in the j-th column of the array AB 
  
            as follows:   
            if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j; 
  
            if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd). 
  

            On exit, if INFO = 0, the triangular factor U or L from the   
            Cholesky factorization A = U'*U or A = L*L' of the band   
            matrix A, in the same storage format as A.   

    LDAB    (input) INTEGER   
            The leading dimension of the array AB.  LDAB >= KD+1.   

    S       (output) DOUBLE PRECISION array, dimension (N)   
            The scale factors for A.   

    SCOND   (input) DOUBLE PRECISION   
            Ratio of the smallest S(i) to the largest S(i).   

    AMAX    (input) DOUBLE PRECISION   
            Absolute value of largest matrix entry.   

    EQUED   (output) CHARACTER*1   
            Specifies whether or not equilibration was done.   
            = 'N':  No equilibration.   
            = 'Y':  Equilibration was done, i.e., A has been replaced by 
  
                    diag(S) * A * diag(S).   

    Internal Parameters   
    ===================   

    THRESH is a threshold value used to decide if scaling should be done 
  
    based on the ratio of the scaling factors.  If SCOND < THRESH,   
    scaling is done.   

    LARGE and SMALL are threshold values used to decide if scaling should 
  
    be done based on the absolute size of the largest matrix element.   
    If AMAX > LARGE or AMAX < SMALL, scaling is done.   

    ===================================================================== 
  


       Quick return if possible   

    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    integer ab_dim1, ab_offset, i__1, i__2, i__3, i__4;
    doublereal d__1;
    doublecomplex z__1;
    /* Local variables */
    static integer i, j;
    static doublereal large;
    extern logical lsame_(char *, char *);
    static doublereal small, cj;
    extern doublereal dlamch_(char *);


#define S(I) s[(I)-1]

#define AB(I,J) ab[(I)-1 + ((J)-1)* ( *ldab)]

    if (*n <= 0) {
	*(unsigned char *)equed = 'N';
	return 0;
    }

/*     Initialize LARGE and SMALL. */

    small = dlamch_("Safe minimum") / dlamch_("Precision");
    large = 1. / small;

    if (*scond >= .1 && *amax >= small && *amax <= large) {

/*        No equilibration */

	*(unsigned char *)equed = 'N';
    } else {

/*        Replace A by diag(S) * A * diag(S). */

	if (lsame_(uplo, "U")) {

/*           Upper triangle of A is stored in band format. */

	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		cj = S(j);
/* Computing MAX */
		i__2 = 1, i__3 = j - *kd;
		i__4 = j - 1;
		for (i = max(1,j-*kd); i <= j-1; ++i) {
		    i__2 = *kd + 1 + i - j + j * ab_dim1;
		    d__1 = cj * S(i);
		    i__3 = *kd + 1 + i - j + j * ab_dim1;
		    z__1.r = d__1 * AB(*kd+1+i-j,j).r, z__1.i = d__1 * AB(*kd+1+i-j,j).i;
		    AB(*kd+1+i-j,j).r = z__1.r, AB(*kd+1+i-j,j).i = z__1.i;
/* L10: */
		}
		i__4 = *kd + 1 + j * ab_dim1;
		i__2 = *kd + 1 + j * ab_dim1;
		d__1 = cj * cj * AB(*kd+1,j).r;
		AB(*kd+1,j).r = d__1, AB(*kd+1,j).i = 0.;
/* L20: */
	    }
	} else {

/*           Lower triangle of A is stored. */

	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		cj = S(j);
		i__4 = j * ab_dim1 + 1;
		i__2 = j * ab_dim1 + 1;
		d__1 = cj * cj * AB(1,j).r;
		AB(1,j).r = d__1, AB(1,j).i = 0.;
/* Computing MIN */
		i__2 = *n, i__3 = j + *kd;
		i__4 = min(i__2,i__3);
		for (i = j + 1; i <= min(*n,j+*kd); ++i) {
		    i__2 = i + 1 - j + j * ab_dim1;
		    d__1 = cj * S(i);
		    i__3 = i + 1 - j + j * ab_dim1;
		    z__1.r = d__1 * AB(i+1-j,j).r, z__1.i = d__1 * AB(i+1-j,j).i;
		    AB(i+1-j,j).r = z__1.r, AB(i+1-j,j).i = z__1.i;
/* L30: */
		}
/* L40: */
	    }
	}
	*(unsigned char *)equed = 'Y';
    }

    return 0;

/*     End of ZLAQHB */

} /* zlaqhb_ */

