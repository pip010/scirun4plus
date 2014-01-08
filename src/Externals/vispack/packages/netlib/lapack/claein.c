#include "f2c.h"

/* Subroutine */ int claein_(logical *rightv, logical *noinit, integer *n, 
	complex *h, integer *ldh, complex *w, complex *v, complex *b, integer 
	*ldb, real *rwork, real *eps3, real *smlnum, integer *info)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    CLAEIN uses inverse iteration to find a right or left eigenvector   
    corresponding to the eigenvalue W of a complex upper Hessenberg   
    matrix H.   

    Arguments   
    =========   

    RIGHTV   (input) LOGICAL   
            = .TRUE. : compute right eigenvector;   
            = .FALSE.: compute left eigenvector.   

    NOINIT   (input) LOGICAL   
            = .TRUE. : no initial vector supplied in V   
            = .FALSE.: initial vector supplied in V.   

    N       (input) INTEGER   
            The order of the matrix H.  N >= 0.   

    H       (input) COMPLEX array, dimension (LDH,N)   
            The upper Hessenberg matrix H.   

    LDH     (input) INTEGER   
            The leading dimension of the array H.  LDH >= max(1,N).   

    W       (input) COMPLEX   
            The eigenvalue of H whose corresponding right or left   
            eigenvector is to be computed.   

    V       (input/output) COMPLEX array, dimension (N)   
            On entry, if NOINIT = .FALSE., V must contain a starting   
            vector for inverse iteration; otherwise V need not be set.   
            On exit, V contains the computed eigenvector, normalized so   
            that the component of largest magnitude has magnitude 1; here 
  
            the magnitude of a complex number (x,y) is taken to be   
            |x| + |y|.   

    B       (workspace) COMPLEX array, dimension (LDB,N)   

    LDB     (input) INTEGER   
            The leading dimension of the array B.  LDB >= max(1,N).   

    RWORK   (workspace) REAL array, dimension (N)   

    EPS3    (input) REAL   
            A small machine-dependent value which is used to perturb   
            close eigenvalues, and to replace zero pivots.   

    SMLNUM  (input) REAL   
            A machine-dependent value close to the underflow threshold.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            = 1:  inverse iteration did not converge; V is set to the   
                  last iterate.   

    ===================================================================== 
  


    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static integer c__1 = 1;
    
    /* System generated locals */
    integer b_dim1, b_offset, h_dim1, h_offset, i__1, i__2, i__3, i__4, i__5;
    real r__1, r__2, r__3, r__4;
    doublereal d__1;
    complex q__1, q__2;
    /* Builtin functions */
    double sqrt(doublereal), r_imag(complex *);
    /* Local variables */
    static integer ierr;
    static complex temp;
    static integer i, j;
    static real scale;
    static complex x;
    static char trans[1];
    static real rtemp, rootn, vnorm;
    extern doublereal scnrm2_(integer *, complex *, integer *);
    static complex ei, ej;
    extern integer icamax_(integer *, complex *, integer *);
    extern /* Complex */ VOID cladiv_(complex *, complex *, complex *);
    extern /* Subroutine */ int csscal_(integer *, real *, complex *, integer 
	    *), clatrs_(char *, char *, char *, char *, integer *, complex *, 
	    integer *, complex *, real *, real *, integer *);
    extern doublereal scasum_(integer *, complex *, integer *);
    static char normin[1];
    static real nrmsml, growto;
    static integer its;



#define V(I) v[(I)-1]
#define RWORK(I) rwork[(I)-1]

#define H(I,J) h[(I)-1 + ((J)-1)* ( *ldh)]
#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]

    *info = 0;

/*     GROWTO is the threshold used in the acceptance test for an   
       eigenvector. */

    rootn = sqrt((real) (*n));
    growto = .1f / rootn;
/* Computing MAX */
    r__1 = 1.f, r__2 = *eps3 * rootn;
    nrmsml = dmax(r__1,r__2) * *smlnum;

/*     Form B = H - W*I (except that the subdiagonal elements are not   
       stored). */

    i__1 = *n;
    for (j = 1; j <= *n; ++j) {
	i__2 = j - 1;
	for (i = 1; i <= j-1; ++i) {
	    i__3 = i + j * b_dim1;
	    i__4 = i + j * h_dim1;
	    B(i,j).r = H(i,j).r, B(i,j).i = H(i,j).i;
/* L10: */
	}
	i__2 = j + j * b_dim1;
	i__3 = j + j * h_dim1;
	q__1.r = H(j,j).r - w->r, q__1.i = H(j,j).i - w->i;
	B(j,j).r = q__1.r, B(j,j).i = q__1.i;
/* L20: */
    }

    if (*noinit) {

/*        Initialize V. */

	i__1 = *n;
	for (i = 1; i <= *n; ++i) {
	    i__2 = i;
	    V(i).r = *eps3, V(i).i = 0.f;
/* L30: */
	}
    } else {

/*        Scale supplied initial vector. */

	vnorm = scnrm2_(n, &V(1), &c__1);
	r__1 = *eps3 * rootn / dmax(vnorm,nrmsml);
	csscal_(n, &r__1, &V(1), &c__1);
    }

    if (*rightv) {

/*        LU decomposition with partial pivoting of B, replacing zero 
  
          pivots by EPS3. */

	i__1 = *n - 1;
	for (i = 1; i <= *n-1; ++i) {
	    i__2 = i + 1 + i * h_dim1;
	    ei.r = H(i+1,i).r, ei.i = H(i+1,i).i;
	    i__2 = i + i * b_dim1;
	    if ((r__1 = B(i,i).r, dabs(r__1)) + (r__2 = r_imag(&B(i,i)), dabs(r__2)) < (r__3 = ei.r, dabs(r__3)) + (r__4 
		    = r_imag(&ei), dabs(r__4))) {

/*              Interchange rows and eliminate. */

		cladiv_(&q__1, &B(i,i), &ei);
		x.r = q__1.r, x.i = q__1.i;
		i__2 = i + i * b_dim1;
		B(i,i).r = ei.r, B(i,i).i = ei.i;
		i__2 = *n;
		for (j = i + 1; j <= *n; ++j) {
		    i__3 = i + 1 + j * b_dim1;
		    temp.r = B(i+1,j).r, temp.i = B(i+1,j).i;
		    i__3 = i + 1 + j * b_dim1;
		    i__4 = i + j * b_dim1;
		    q__2.r = x.r * temp.r - x.i * temp.i, q__2.i = x.r * 
			    temp.i + x.i * temp.r;
		    q__1.r = B(i,j).r - q__2.r, q__1.i = B(i,j).i - q__2.i;
		    B(i+1,j).r = q__1.r, B(i+1,j).i = q__1.i;
		    i__3 = i + j * b_dim1;
		    B(i,j).r = temp.r, B(i,j).i = temp.i;
/* L40: */
		}
	    } else {

/*              Eliminate without interchange. */

		i__2 = i + i * b_dim1;
		if (B(i,i).r == 0.f && B(i,i).i == 0.f) {
		    i__3 = i + i * b_dim1;
		    B(i,i).r = *eps3, B(i,i).i = 0.f;
		}
		cladiv_(&q__1, &ei, &B(i,i));
		x.r = q__1.r, x.i = q__1.i;
		if (x.r != 0.f || x.i != 0.f) {
		    i__2 = *n;
		    for (j = i + 1; j <= *n; ++j) {
			i__3 = i + 1 + j * b_dim1;
			i__4 = i + 1 + j * b_dim1;
			i__5 = i + j * b_dim1;
			q__2.r = x.r * B(i,j).r - x.i * B(i,j).i, q__2.i = 
				x.r * B(i,j).i + x.i * B(i,j).r;
			q__1.r = B(i+1,j).r - q__2.r, q__1.i = B(i+1,j).i - 
				q__2.i;
			B(i+1,j).r = q__1.r, B(i+1,j).i = q__1.i;
/* L50: */
		    }
		}
	    }
/* L60: */
	}
	i__1 = *n + *n * b_dim1;
	if (B(*n,*n).r == 0.f && B(*n,*n).i == 0.f) {
	    i__2 = *n + *n * b_dim1;
	    B(*n,*n).r = *eps3, B(*n,*n).i = 0.f;
	}

	*(unsigned char *)trans = 'N';

    } else {

/*        UL decomposition with partial pivoting of B, replacing zero 
  
          pivots by EPS3. */

	for (j = *n; j >= 2; --j) {
	    i__1 = j + (j - 1) * h_dim1;
	    ej.r = H(j,j-1).r, ej.i = H(j,j-1).i;
	    i__1 = j + j * b_dim1;
	    if ((r__1 = B(j,j).r, dabs(r__1)) + (r__2 = r_imag(&B(j,j)), dabs(r__2)) < (r__3 = ej.r, dabs(r__3)) + (r__4 
		    = r_imag(&ej), dabs(r__4))) {

/*              Interchange columns and eliminate. */

		cladiv_(&q__1, &B(j,j), &ej);
		x.r = q__1.r, x.i = q__1.i;
		i__1 = j + j * b_dim1;
		B(j,j).r = ej.r, B(j,j).i = ej.i;
		i__1 = j - 1;
		for (i = 1; i <= j-1; ++i) {
		    i__2 = i + (j - 1) * b_dim1;
		    temp.r = B(i,j-1).r, temp.i = B(i,j-1).i;
		    i__2 = i + (j - 1) * b_dim1;
		    i__3 = i + j * b_dim1;
		    q__2.r = x.r * temp.r - x.i * temp.i, q__2.i = x.r * 
			    temp.i + x.i * temp.r;
		    q__1.r = B(i,j).r - q__2.r, q__1.i = B(i,j).i - q__2.i;
		    B(i,j-1).r = q__1.r, B(i,j-1).i = q__1.i;
		    i__2 = i + j * b_dim1;
		    B(i,j).r = temp.r, B(i,j).i = temp.i;
/* L70: */
		}
	    } else {

/*              Eliminate without interchange. */

		i__1 = j + j * b_dim1;
		if (B(j,j).r == 0.f && B(j,j).i == 0.f) {
		    i__2 = j + j * b_dim1;
		    B(j,j).r = *eps3, B(j,j).i = 0.f;
		}
		cladiv_(&q__1, &ej, &B(j,j));
		x.r = q__1.r, x.i = q__1.i;
		if (x.r != 0.f || x.i != 0.f) {
		    i__1 = j - 1;
		    for (i = 1; i <= j-1; ++i) {
			i__2 = i + (j - 1) * b_dim1;
			i__3 = i + (j - 1) * b_dim1;
			i__4 = i + j * b_dim1;
			q__2.r = x.r * B(i,j).r - x.i * B(i,j).i, q__2.i = 
				x.r * B(i,j).i + x.i * B(i,j).r;
			q__1.r = B(i,j-1).r - q__2.r, q__1.i = B(i,j-1).i - 
				q__2.i;
			B(i,j-1).r = q__1.r, B(i,j-1).i = q__1.i;
/* L80: */
		    }
		}
	    }
/* L90: */
	}
	i__1 = b_dim1 + 1;
	if (B(1,1).r == 0.f && B(1,1).i == 0.f) {
	    i__2 = b_dim1 + 1;
	    B(1,1).r = *eps3, B(1,1).i = 0.f;
	}

	*(unsigned char *)trans = 'C';

    }

    *(unsigned char *)normin = 'N';
    i__1 = *n;
    for (its = 1; its <= *n; ++its) {

/*        Solve U*x = scale*v for a right eigenvector   
            or U'*x = scale*v for a left eigenvector,   
          overwriting x on v. */

	clatrs_("Upper", trans, "Nonunit", normin, n, &B(1,1), ldb, &V(1)
		, &scale, &RWORK(1), &ierr);
	*(unsigned char *)normin = 'Y';

/*        Test for sufficient growth in the norm of v. */

	vnorm = scasum_(n, &V(1), &c__1);
	if (vnorm >= growto * scale) {
	    goto L120;
	}

/*        Choose new orthogonal starting vector and try again. */

	rtemp = *eps3 / (rootn + 1.f);
	V(1).r = *eps3, V(1).i = 0.f;
	i__2 = *n;
	for (i = 2; i <= *n; ++i) {
	    i__3 = i;
	    V(i).r = rtemp, V(i).i = 0.f;
/* L100: */
	}
	i__2 = *n - its + 1;
	i__3 = *n - its + 1;
	d__1 = *eps3 * rootn;
	q__1.r = V(*n-its+1).r - d__1, q__1.i = V(*n-its+1).i;
	V(*n-its+1).r = q__1.r, V(*n-its+1).i = q__1.i;
/* L110: */
    }

/*     Failure to find eigenvector in N iterations. */

    *info = 1;

L120:

/*     Normalize eigenvector. */

    i = icamax_(n, &V(1), &c__1);
    i__1 = i;
    r__3 = 1.f / ((r__1 = V(i).r, dabs(r__1)) + (r__2 = r_imag(&V(i)), 
	    dabs(r__2)));
    csscal_(n, &r__3, &V(1), &c__1);

    return 0;

/*     End of CLAEIN */

} /* claein_ */

