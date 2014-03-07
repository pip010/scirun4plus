
/*  -- translated by f2c (version 19940927).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Subroutine */ int zherk_(char *uplo, char *trans, integer *n, integer *k, 
	doublereal *alpha, doublecomplex *a, integer *lda, doublereal *beta, 
	doublecomplex *c, integer *ldc)
{


    /* System generated locals */
    integer a_dim1, a_offset, c_dim1, c_offset, i__1, i__2, i__3, i__4, i__5, 
	    i__6;
    doublereal d__1;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer info;
    static doublecomplex temp;
    static integer i, j, l;
    extern logical lsame_(char *, char *);
    static integer nrowa;
    static doublereal rtemp;
    static logical upper;
    extern /* Subroutine */ int xerbla_(char *, integer *);


/*  Purpose   
    =======   

    ZHERK  performs one of the hermitian rank k operations   

       C := alpha*A*conjg( A' ) + beta*C,   

    or   

       C := alpha*conjg( A' )*A + beta*C,   

    where  alpha and beta  are  real scalars,  C is an  n by n  hermitian 
  
    matrix and  A  is an  n by k  matrix in the  first case and a  k by n 
  
    matrix in the second case.   

    Parameters   
    ==========   

    UPLO   - CHARACTER*1.   
             On  entry,   UPLO  specifies  whether  the  upper  or  lower 
  
             triangular  part  of the  array  C  is to be  referenced  as 
  
             follows:   

                UPLO = 'U' or 'u'   Only the  upper triangular part of  C 
  
                                    is to be referenced.   

                UPLO = 'L' or 'l'   Only the  lower triangular part of  C 
  
                                    is to be referenced.   

             Unchanged on exit.   

    TRANS  - CHARACTER*1.   
             On entry,  TRANS  specifies the operation to be performed as 
  
             follows:   

                TRANS = 'N' or 'n'   C := alpha*A*conjg( A' ) + beta*C.   

                TRANS = 'C' or 'c'   C := alpha*conjg( A' )*A + beta*C.   

             Unchanged on exit.   

    N      - INTEGER.   
             On entry,  N specifies the order of the matrix C.  N must be 
  
             at least zero.   
             Unchanged on exit.   

    K      - INTEGER.   
             On entry with  TRANS = 'N' or 'n',  K  specifies  the number 
  
             of  columns   of  the   matrix   A,   and  on   entry   with 
  
             TRANS = 'C' or 'c',  K  specifies  the number of rows of the 
  
             matrix A.  K must be at least zero.   
             Unchanged on exit.   

    ALPHA  - DOUBLE PRECISION            .   
             On entry, ALPHA specifies the scalar alpha.   
             Unchanged on exit.   

    A      - COMPLEX*16       array of DIMENSION ( LDA, ka ), where ka is 
  
             k  when  TRANS = 'N' or 'n',  and is  n  otherwise.   
             Before entry with  TRANS = 'N' or 'n',  the  leading  n by k 
  
             part of the array  A  must contain the matrix  A,  otherwise 
  
             the leading  k by n  part of the array  A  must contain  the 
  
             matrix A.   
             Unchanged on exit.   

    LDA    - INTEGER.   
             On entry, LDA specifies the first dimension of A as declared 
  
             in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n' 
  
             then  LDA must be at least  max( 1, n ), otherwise  LDA must 
  
             be at least  max( 1, k ).   
             Unchanged on exit.   

    BETA   - DOUBLE PRECISION.   
             On entry, BETA specifies the scalar beta.   
             Unchanged on exit.   

    C      - COMPLEX*16          array of DIMENSION ( LDC, n ).   
             Before entry  with  UPLO = 'U' or 'u',  the leading  n by n 
  
             upper triangular part of the array C must contain the upper 
  
             triangular part  of the  hermitian matrix  and the strictly 
  
             lower triangular part of C is not referenced.  On exit, the 
  
             upper triangular part of the array  C is overwritten by the 
  
             upper triangular part of the updated matrix.   
             Before entry  with  UPLO = 'L' or 'l',  the leading  n by n 
  
             lower triangular part of the array C must contain the lower 
  
             triangular part  of the  hermitian matrix  and the strictly 
  
             upper triangular part of C is not referenced.  On exit, the 
  
             lower triangular part of the array  C is overwritten by the 
  
             lower triangular part of the updated matrix.   
             Note that the imaginary parts of the diagonal elements need 
  
             not be set,  they are assumed to be zero,  and on exit they 
  
             are set to zero.   

    LDC    - INTEGER.   
             On entry, LDC specifies the first dimension of C as declared 
  
             in  the  calling  (sub)  program.   LDC  must  be  at  least 
  
             max( 1, n ).   
             Unchanged on exit.   


    Level 3 Blas routine.   

    -- Written on 8-February-1989.   
       Jack Dongarra, Argonne National Laboratory.   
       Iain Duff, AERE Harwell.   
       Jeremy Du Croz, Numerical Algorithms Group Ltd.   
       Sven Hammarling, Numerical Algorithms Group Ltd.   

    -- Modified 8-Nov-93 to set C(J,J) to DBLE( C(J,J) ) when BETA = 1.   
       Ed Anderson, Cray Research Inc.   



       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]
#define C(I,J) c[(I)-1 + ((J)-1)* ( *ldc)]

    if (lsame_(trans, "N")) {
	nrowa = *n;
    } else {
	nrowa = *k;
    }
    upper = lsame_(uplo, "U");

    info = 0;
    if (! upper && ! lsame_(uplo, "L")) {
	info = 1;
    } else if (! lsame_(trans, "N") && ! lsame_(trans, "C")) {
	info = 2;
    } else if (*n < 0) {
	info = 3;
    } else if (*k < 0) {
	info = 4;
    } else if (*lda < max(1,nrowa)) {
	info = 7;
    } else if (*ldc < max(1,*n)) {
	info = 10;
    }
    if (info != 0) {
	xerbla_("ZHERK ", &info);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0 || (*alpha == 0. || *k == 0) && *beta == 1.) {
	return 0;
    }

/*     And when  alpha.eq.zero. */

    if (*alpha == 0.) {
	if (upper) {
	    if (*beta == 0.) {
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
		    i__2 = j;
		    for (i = 1; i <= j; ++i) {
			i__3 = i + j * c_dim1;
			C(i,j).r = 0., C(i,j).i = 0.;
/* L10: */
		    }
/* L20: */
		}
	    } else {
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
		    i__2 = j - 1;
		    for (i = 1; i <= j-1; ++i) {
			i__3 = i + j * c_dim1;
			i__4 = i + j * c_dim1;
			z__1.r = *beta * C(i,j).r, z__1.i = *beta * C(i,j)
				.i;
			C(i,j).r = z__1.r, C(i,j).i = z__1.i;
/* L30: */
		    }
		    i__2 = j + j * c_dim1;
		    i__3 = j + j * c_dim1;
		    d__1 = *beta * C(j,j).r;
		    C(j,j).r = d__1, C(j,j).i = 0.;
/* L40: */
		}
	    }
	} else {
	    if (*beta == 0.) {
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
		    i__2 = *n;
		    for (i = j; i <= *n; ++i) {
			i__3 = i + j * c_dim1;
			C(i,j).r = 0., C(i,j).i = 0.;
/* L50: */
		    }
/* L60: */
		}
	    } else {
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
		    i__2 = j + j * c_dim1;
		    i__3 = j + j * c_dim1;
		    d__1 = *beta * C(j,j).r;
		    C(j,j).r = d__1, C(j,j).i = 0.;
		    i__2 = *n;
		    for (i = j + 1; i <= *n; ++i) {
			i__3 = i + j * c_dim1;
			i__4 = i + j * c_dim1;
			z__1.r = *beta * C(i,j).r, z__1.i = *beta * C(i,j)
				.i;
			C(i,j).r = z__1.r, C(i,j).i = z__1.i;
/* L70: */
		    }
/* L80: */
		}
	    }
	}
	return 0;
    }

/*     Start the operations. */

    if (lsame_(trans, "N")) {

/*        Form  C := alpha*A*conjg( A' ) + beta*C. */

	if (upper) {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		if (*beta == 0.) {
		    i__2 = j;
		    for (i = 1; i <= j; ++i) {
			i__3 = i + j * c_dim1;
			C(i,j).r = 0., C(i,j).i = 0.;
/* L90: */
		    }
		} else if (*beta != 1.) {
		    i__2 = j - 1;
		    for (i = 1; i <= j-1; ++i) {
			i__3 = i + j * c_dim1;
			i__4 = i + j * c_dim1;
			z__1.r = *beta * C(i,j).r, z__1.i = *beta * C(i,j)
				.i;
			C(i,j).r = z__1.r, C(i,j).i = z__1.i;
/* L100: */
		    }
		    i__2 = j + j * c_dim1;
		    i__3 = j + j * c_dim1;
		    d__1 = *beta * C(j,j).r;
		    C(j,j).r = d__1, C(j,j).i = 0.;
		} else {
		    i__2 = j + j * c_dim1;
		    i__3 = j + j * c_dim1;
		    d__1 = C(j,j).r;
		    C(j,j).r = d__1, C(j,j).i = 0.;
		}
		i__2 = *k;
		for (l = 1; l <= *k; ++l) {
		    i__3 = j + l * a_dim1;
		    if (A(j,l).r != 0. || A(j,l).i != 0.) {
			d_cnjg(&z__2, &A(j,l));
			z__1.r = *alpha * z__2.r, z__1.i = *alpha * z__2.i;
			temp.r = z__1.r, temp.i = z__1.i;
			i__3 = j - 1;
			for (i = 1; i <= j-1; ++i) {
			    i__4 = i + j * c_dim1;
			    i__5 = i + j * c_dim1;
			    i__6 = i + l * a_dim1;
			    z__2.r = temp.r * A(i,l).r - temp.i * A(i,l).i, 
				    z__2.i = temp.r * A(i,l).i + temp.i * A(i,l).r;
			    z__1.r = C(i,j).r + z__2.r, z__1.i = C(i,j).i + 
				    z__2.i;
			    C(i,j).r = z__1.r, C(i,j).i = z__1.i;
/* L110: */
			}
			i__3 = j + j * c_dim1;
			i__4 = j + j * c_dim1;
			i__5 = i + l * a_dim1;
			z__1.r = temp.r * A(i,l).r - temp.i * A(i,l).i, 
				z__1.i = temp.r * A(i,l).i + temp.i * A(i,l)
				.r;
			d__1 = C(j,j).r + z__1.r;
			C(j,j).r = d__1, C(j,j).i = 0.;
		    }
/* L120: */
		}
/* L130: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		if (*beta == 0.) {
		    i__2 = *n;
		    for (i = j; i <= *n; ++i) {
			i__3 = i + j * c_dim1;
			C(i,j).r = 0., C(i,j).i = 0.;
/* L140: */
		    }
		} else if (*beta != 1.) {
		    i__2 = j + j * c_dim1;
		    i__3 = j + j * c_dim1;
		    d__1 = *beta * C(j,j).r;
		    C(j,j).r = d__1, C(j,j).i = 0.;
		    i__2 = *n;
		    for (i = j + 1; i <= *n; ++i) {
			i__3 = i + j * c_dim1;
			i__4 = i + j * c_dim1;
			z__1.r = *beta * C(i,j).r, z__1.i = *beta * C(i,j)
				.i;
			C(i,j).r = z__1.r, C(i,j).i = z__1.i;
/* L150: */
		    }
		} else {
		    i__2 = j + j * c_dim1;
		    i__3 = j + j * c_dim1;
		    d__1 = C(j,j).r;
		    C(j,j).r = d__1, C(j,j).i = 0.;
		}
		i__2 = *k;
		for (l = 1; l <= *k; ++l) {
		    i__3 = j + l * a_dim1;
		    if (A(j,l).r != 0. || A(j,l).i != 0.) {
			d_cnjg(&z__2, &A(j,l));
			z__1.r = *alpha * z__2.r, z__1.i = *alpha * z__2.i;
			temp.r = z__1.r, temp.i = z__1.i;
			i__3 = j + j * c_dim1;
			i__4 = j + j * c_dim1;
			i__5 = j + l * a_dim1;
			z__1.r = temp.r * A(j,l).r - temp.i * A(j,l).i, 
				z__1.i = temp.r * A(j,l).i + temp.i * A(j,l)
				.r;
			d__1 = C(j,j).r + z__1.r;
			C(j,j).r = d__1, C(j,j).i = 0.;
			i__3 = *n;
			for (i = j + 1; i <= *n; ++i) {
			    i__4 = i + j * c_dim1;
			    i__5 = i + j * c_dim1;
			    i__6 = i + l * a_dim1;
			    z__2.r = temp.r * A(i,l).r - temp.i * A(i,l).i, 
				    z__2.i = temp.r * A(i,l).i + temp.i * A(i,l).r;
			    z__1.r = C(i,j).r + z__2.r, z__1.i = C(i,j).i + 
				    z__2.i;
			    C(i,j).r = z__1.r, C(i,j).i = z__1.i;
/* L160: */
			}
		    }
/* L170: */
		}
/* L180: */
	    }
	}
    } else {

/*        Form  C := alpha*conjg( A' )*A + beta*C. */

	if (upper) {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		i__2 = j - 1;
		for (i = 1; i <= j-1; ++i) {
		    temp.r = 0., temp.i = 0.;
		    i__3 = *k;
		    for (l = 1; l <= *k; ++l) {
			d_cnjg(&z__3, &A(l,i));
			i__4 = l + j * a_dim1;
			z__2.r = z__3.r * A(l,j).r - z__3.i * A(l,j).i, 
				z__2.i = z__3.r * A(l,j).i + z__3.i * A(l,j)
				.r;
			z__1.r = temp.r + z__2.r, z__1.i = temp.i + z__2.i;
			temp.r = z__1.r, temp.i = z__1.i;
/* L190: */
		    }
		    if (*beta == 0.) {
			i__3 = i + j * c_dim1;
			z__1.r = *alpha * temp.r, z__1.i = *alpha * temp.i;
			C(i,j).r = z__1.r, C(i,j).i = z__1.i;
		    } else {
			i__3 = i + j * c_dim1;
			z__2.r = *alpha * temp.r, z__2.i = *alpha * temp.i;
			i__4 = i + j * c_dim1;
			z__3.r = *beta * C(i,j).r, z__3.i = *beta * C(i,j)
				.i;
			z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
			C(i,j).r = z__1.r, C(i,j).i = z__1.i;
		    }
/* L200: */
		}
		rtemp = 0.;
		i__2 = *k;
		for (l = 1; l <= *k; ++l) {
		    d_cnjg(&z__3, &A(l,j));
		    i__3 = l + j * a_dim1;
		    z__2.r = z__3.r * A(l,j).r - z__3.i * A(l,j).i, z__2.i =
			     z__3.r * A(l,j).i + z__3.i * A(l,j).r;
		    z__1.r = rtemp + z__2.r, z__1.i = z__2.i;
		    rtemp = z__1.r;
/* L210: */
		}
		if (*beta == 0.) {
		    i__2 = j + j * c_dim1;
		    d__1 = *alpha * rtemp;
		    C(j,j).r = d__1, C(j,j).i = 0.;
		} else {
		    i__2 = j + j * c_dim1;
		    i__3 = j + j * c_dim1;
		    d__1 = *alpha * rtemp + *beta * C(j,j).r;
		    C(j,j).r = d__1, C(j,j).i = 0.;
		}
/* L220: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		rtemp = 0.;
		i__2 = *k;
		for (l = 1; l <= *k; ++l) {
		    d_cnjg(&z__3, &A(l,j));
		    i__3 = l + j * a_dim1;
		    z__2.r = z__3.r * A(l,j).r - z__3.i * A(l,j).i, z__2.i =
			     z__3.r * A(l,j).i + z__3.i * A(l,j).r;
		    z__1.r = rtemp + z__2.r, z__1.i = z__2.i;
		    rtemp = z__1.r;
/* L230: */
		}
		if (*beta == 0.) {
		    i__2 = j + j * c_dim1;
		    d__1 = *alpha * rtemp;
		    C(j,j).r = d__1, C(j,j).i = 0.;
		} else {
		    i__2 = j + j * c_dim1;
		    i__3 = j + j * c_dim1;
		    d__1 = *alpha * rtemp + *beta * C(j,j).r;
		    C(j,j).r = d__1, C(j,j).i = 0.;
		}
		i__2 = *n;
		for (i = j + 1; i <= *n; ++i) {
		    temp.r = 0., temp.i = 0.;
		    i__3 = *k;
		    for (l = 1; l <= *k; ++l) {
			d_cnjg(&z__3, &A(l,i));
			i__4 = l + j * a_dim1;
			z__2.r = z__3.r * A(l,j).r - z__3.i * A(l,j).i, 
				z__2.i = z__3.r * A(l,j).i + z__3.i * A(l,j)
				.r;
			z__1.r = temp.r + z__2.r, z__1.i = temp.i + z__2.i;
			temp.r = z__1.r, temp.i = z__1.i;
/* L240: */
		    }
		    if (*beta == 0.) {
			i__3 = i + j * c_dim1;
			z__1.r = *alpha * temp.r, z__1.i = *alpha * temp.i;
			C(i,j).r = z__1.r, C(i,j).i = z__1.i;
		    } else {
			i__3 = i + j * c_dim1;
			z__2.r = *alpha * temp.r, z__2.i = *alpha * temp.i;
			i__4 = i + j * c_dim1;
			z__3.r = *beta * C(i,j).r, z__3.i = *beta * C(i,j)
				.i;
			z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
			C(i,j).r = z__1.r, C(i,j).i = z__1.i;
		    }
/* L250: */
		}
/* L260: */
	    }
	}
    }

    return 0;

/*     End of ZHERK . */

} /* zherk_ */

