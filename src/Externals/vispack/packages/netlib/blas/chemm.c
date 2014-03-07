
/*  -- translated by f2c (version 19940927).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Subroutine */ int chemm_(char *side, char *uplo, integer *m, integer *n, 
	complex *alpha, complex *a, integer *lda, complex *b, integer *ldb, 
	complex *beta, complex *c, integer *ldc)
{


    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, i__1, i__2, 
	    i__3, i__4, i__5, i__6;
    doublereal d__1;
    complex q__1, q__2, q__3, q__4, q__5;

    /* Builtin functions */
    void r_cnjg(complex *, complex *);

    /* Local variables */
    static integer info;
    static complex temp1, temp2;
    static integer i, j, k;
    extern logical lsame_(char *, char *);
    static integer nrowa;
    static logical upper;
    extern /* Subroutine */ int xerbla_(char *, integer *);


/*  Purpose   
    =======   

    CHEMM  performs one of the matrix-matrix operations   

       C := alpha*A*B + beta*C,   

    or   

       C := alpha*B*A + beta*C,   

    where alpha and beta are scalars, A is an hermitian matrix and  B and 
  
    C are m by n matrices.   

    Parameters   
    ==========   

    SIDE   - CHARACTER*1.   
             On entry,  SIDE  specifies whether  the  hermitian matrix  A 
  
             appears on the  left or right  in the  operation as follows: 
  

                SIDE = 'L' or 'l'   C := alpha*A*B + beta*C,   

                SIDE = 'R' or 'r'   C := alpha*B*A + beta*C,   

             Unchanged on exit.   

    UPLO   - CHARACTER*1.   
             On  entry,   UPLO  specifies  whether  the  upper  or  lower 
  
             triangular  part  of  the  hermitian  matrix   A  is  to  be 
  
             referenced as follows:   

                UPLO = 'U' or 'u'   Only the upper triangular part of the 
  
                                    hermitian matrix is to be referenced. 
  

                UPLO = 'L' or 'l'   Only the lower triangular part of the 
  
                                    hermitian matrix is to be referenced. 
  

             Unchanged on exit.   

    M      - INTEGER.   
             On entry,  M  specifies the number of rows of the matrix  C. 
  
             M  must be at least zero.   
             Unchanged on exit.   

    N      - INTEGER.   
             On entry, N specifies the number of columns of the matrix C. 
  
             N  must be at least zero.   
             Unchanged on exit.   

    ALPHA  - COMPLEX         .   
             On entry, ALPHA specifies the scalar alpha.   
             Unchanged on exit.   

    A      - COMPLEX          array of DIMENSION ( LDA, ka ), where ka is 
  
             m  when  SIDE = 'L' or 'l'  and is n  otherwise.   
             Before entry  with  SIDE = 'L' or 'l',  the  m by m  part of 
  
             the array  A  must contain the  hermitian matrix,  such that 
  
             when  UPLO = 'U' or 'u', the leading m by m upper triangular 
  
             part of the array  A  must contain the upper triangular part 
  
             of the  hermitian matrix and the  strictly  lower triangular 
  
             part of  A  is not referenced,  and when  UPLO = 'L' or 'l', 
  
             the leading  m by m  lower triangular part  of the  array  A 
  
             must  contain  the  lower triangular part  of the  hermitian 
  
             matrix and the  strictly upper triangular part of  A  is not 
  
             referenced.   
             Before entry  with  SIDE = 'R' or 'r',  the  n by n  part of 
  
             the array  A  must contain the  hermitian matrix,  such that 
  
             when  UPLO = 'U' or 'u', the leading n by n upper triangular 
  
             part of the array  A  must contain the upper triangular part 
  
             of the  hermitian matrix and the  strictly  lower triangular 
  
             part of  A  is not referenced,  and when  UPLO = 'L' or 'l', 
  
             the leading  n by n  lower triangular part  of the  array  A 
  
             must  contain  the  lower triangular part  of the  hermitian 
  
             matrix and the  strictly upper triangular part of  A  is not 
  
             referenced.   
             Note that the imaginary parts  of the diagonal elements need 
  
             not be set, they are assumed to be zero.   
             Unchanged on exit.   

    LDA    - INTEGER.   
             On entry, LDA specifies the first dimension of A as declared 
  
             in the  calling (sub) program. When  SIDE = 'L' or 'l'  then 
  
             LDA must be at least  max( 1, m ), otherwise  LDA must be at 
  
             least max( 1, n ).   
             Unchanged on exit.   

    B      - COMPLEX          array of DIMENSION ( LDB, n ).   
             Before entry, the leading  m by n part of the array  B  must 
  
             contain the matrix B.   
             Unchanged on exit.   

    LDB    - INTEGER.   
             On entry, LDB specifies the first dimension of B as declared 
  
             in  the  calling  (sub)  program.   LDB  must  be  at  least 
  
             max( 1, m ).   
             Unchanged on exit.   

    BETA   - COMPLEX         .   
             On entry,  BETA  specifies the scalar  beta.  When  BETA  is 
  
             supplied as zero then C need not be set on input.   
             Unchanged on exit.   

    C      - COMPLEX          array of DIMENSION ( LDC, n ).   
             Before entry, the leading  m by n  part of the array  C must 
  
             contain the matrix  C,  except when  beta  is zero, in which 
  
             case C need not be set on entry.   
             On exit, the array  C  is overwritten by the  m by n updated 
  
             matrix.   

    LDC    - INTEGER.   
             On entry, LDC specifies the first dimension of C as declared 
  
             in  the  calling  (sub)  program.   LDC  must  be  at  least 
  
             max( 1, m ).   
             Unchanged on exit.   


    Level 3 Blas routine.   

    -- Written on 8-February-1989.   
       Jack Dongarra, Argonne National Laboratory.   
       Iain Duff, AERE Harwell.   
       Jeremy Du Croz, Numerical Algorithms Group Ltd.   
       Sven Hammarling, Numerical Algorithms Group Ltd.   



       Set NROWA as the number of rows of A.   

    
   Parameter adjustments   
       Function Body */

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]
#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]
#define C(I,J) c[(I)-1 + ((J)-1)* ( *ldc)]

    if (lsame_(side, "L")) {
	nrowa = *m;
    } else {
	nrowa = *n;
    }
    upper = lsame_(uplo, "U");

/*     Test the input parameters. */

    info = 0;
    if (! lsame_(side, "L") && ! lsame_(side, "R")) {
	info = 1;
    } else if (! upper && ! lsame_(uplo, "L")) {
	info = 2;
    } else if (*m < 0) {
	info = 3;
    } else if (*n < 0) {
	info = 4;
    } else if (*lda < max(1,nrowa)) {
	info = 7;
    } else if (*ldb < max(1,*m)) {
	info = 9;
    } else if (*ldc < max(1,*m)) {
	info = 12;
    }
    if (info != 0) {
	xerbla_("CHEMM ", &info);
	return 0;
    }

/*     Quick return if possible. */

    if (*m == 0 || *n == 0 || alpha->r == 0.f && alpha->i == 0.f && (beta->r 
	    == 1.f && beta->i == 0.f)) {
	return 0;
    }

/*     And when  alpha.eq.zero. */

    if (alpha->r == 0.f && alpha->i == 0.f) {
	if (beta->r == 0.f && beta->i == 0.f) {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		i__2 = *m;
		for (i = 1; i <= *m; ++i) {
		    i__3 = i + j * c_dim1;
		    C(i,j).r = 0.f, C(i,j).i = 0.f;
/* L10: */
		}
/* L20: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		i__2 = *m;
		for (i = 1; i <= *m; ++i) {
		    i__3 = i + j * c_dim1;
		    i__4 = i + j * c_dim1;
		    q__1.r = beta->r * C(i,j).r - beta->i * C(i,j).i, 
			    q__1.i = beta->r * C(i,j).i + beta->i * C(i,j)
			    .r;
		    C(i,j).r = q__1.r, C(i,j).i = q__1.i;
/* L30: */
		}
/* L40: */
	    }
	}
	return 0;
    }

/*     Start the operations. */

    if (lsame_(side, "L")) {

/*        Form  C := alpha*A*B + beta*C. */

	if (upper) {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		i__2 = *m;
		for (i = 1; i <= *m; ++i) {
		    i__3 = i + j * b_dim1;
		    q__1.r = alpha->r * B(i,j).r - alpha->i * B(i,j).i, 
			    q__1.i = alpha->r * B(i,j).i + alpha->i * B(i,j)
			    .r;
		    temp1.r = q__1.r, temp1.i = q__1.i;
		    temp2.r = 0.f, temp2.i = 0.f;
		    i__3 = i - 1;
		    for (k = 1; k <= i-1; ++k) {
			i__4 = k + j * c_dim1;
			i__5 = k + j * c_dim1;
			i__6 = k + i * a_dim1;
			q__2.r = temp1.r * A(k,i).r - temp1.i * A(k,i).i, 
				q__2.i = temp1.r * A(k,i).i + temp1.i * A(k,i).r;
			q__1.r = C(k,j).r + q__2.r, q__1.i = C(k,j).i + 
				q__2.i;
			C(k,j).r = q__1.r, C(k,j).i = q__1.i;
			i__4 = k + j * b_dim1;
			r_cnjg(&q__3, &A(k,i));
			q__2.r = B(k,j).r * q__3.r - B(k,j).i * q__3.i, 
				q__2.i = B(k,j).r * q__3.i + B(k,j).i * 
				q__3.r;
			q__1.r = temp2.r + q__2.r, q__1.i = temp2.i + q__2.i;
			temp2.r = q__1.r, temp2.i = q__1.i;
/* L50: */
		    }
		    if (beta->r == 0.f && beta->i == 0.f) {
			i__3 = i + j * c_dim1;
			i__4 = i + i * a_dim1;
			d__1 = A(i,i).r;
			q__2.r = d__1 * temp1.r, q__2.i = d__1 * temp1.i;
			q__3.r = alpha->r * temp2.r - alpha->i * temp2.i, 
				q__3.i = alpha->r * temp2.i + alpha->i * 
				temp2.r;
			q__1.r = q__2.r + q__3.r, q__1.i = q__2.i + q__3.i;
			C(i,j).r = q__1.r, C(i,j).i = q__1.i;
		    } else {
			i__3 = i + j * c_dim1;
			i__4 = i + j * c_dim1;
			q__3.r = beta->r * C(i,j).r - beta->i * C(i,j).i, 
				q__3.i = beta->r * C(i,j).i + beta->i * C(i,j).r;
			i__5 = i + i * a_dim1;
			d__1 = A(i,i).r;
			q__4.r = d__1 * temp1.r, q__4.i = d__1 * temp1.i;
			q__2.r = q__3.r + q__4.r, q__2.i = q__3.i + q__4.i;
			q__5.r = alpha->r * temp2.r - alpha->i * temp2.i, 
				q__5.i = alpha->r * temp2.i + alpha->i * 
				temp2.r;
			q__1.r = q__2.r + q__5.r, q__1.i = q__2.i + q__5.i;
			C(i,j).r = q__1.r, C(i,j).i = q__1.i;
		    }
/* L60: */
		}
/* L70: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		for (i = *m; i >= 1; --i) {
		    i__2 = i + j * b_dim1;
		    q__1.r = alpha->r * B(i,j).r - alpha->i * B(i,j).i, 
			    q__1.i = alpha->r * B(i,j).i + alpha->i * B(i,j)
			    .r;
		    temp1.r = q__1.r, temp1.i = q__1.i;
		    temp2.r = 0.f, temp2.i = 0.f;
		    i__2 = *m;
		    for (k = i + 1; k <= *m; ++k) {
			i__3 = k + j * c_dim1;
			i__4 = k + j * c_dim1;
			i__5 = k + i * a_dim1;
			q__2.r = temp1.r * A(k,i).r - temp1.i * A(k,i).i, 
				q__2.i = temp1.r * A(k,i).i + temp1.i * A(k,i).r;
			q__1.r = C(k,j).r + q__2.r, q__1.i = C(k,j).i + 
				q__2.i;
			C(k,j).r = q__1.r, C(k,j).i = q__1.i;
			i__3 = k + j * b_dim1;
			r_cnjg(&q__3, &A(k,i));
			q__2.r = B(k,j).r * q__3.r - B(k,j).i * q__3.i, 
				q__2.i = B(k,j).r * q__3.i + B(k,j).i * 
				q__3.r;
			q__1.r = temp2.r + q__2.r, q__1.i = temp2.i + q__2.i;
			temp2.r = q__1.r, temp2.i = q__1.i;
/* L80: */
		    }
		    if (beta->r == 0.f && beta->i == 0.f) {
			i__2 = i + j * c_dim1;
			i__3 = i + i * a_dim1;
			d__1 = A(i,i).r;
			q__2.r = d__1 * temp1.r, q__2.i = d__1 * temp1.i;
			q__3.r = alpha->r * temp2.r - alpha->i * temp2.i, 
				q__3.i = alpha->r * temp2.i + alpha->i * 
				temp2.r;
			q__1.r = q__2.r + q__3.r, q__1.i = q__2.i + q__3.i;
			C(i,j).r = q__1.r, C(i,j).i = q__1.i;
		    } else {
			i__2 = i + j * c_dim1;
			i__3 = i + j * c_dim1;
			q__3.r = beta->r * C(i,j).r - beta->i * C(i,j).i, 
				q__3.i = beta->r * C(i,j).i + beta->i * C(i,j).r;
			i__4 = i + i * a_dim1;
			d__1 = A(i,i).r;
			q__4.r = d__1 * temp1.r, q__4.i = d__1 * temp1.i;
			q__2.r = q__3.r + q__4.r, q__2.i = q__3.i + q__4.i;
			q__5.r = alpha->r * temp2.r - alpha->i * temp2.i, 
				q__5.i = alpha->r * temp2.i + alpha->i * 
				temp2.r;
			q__1.r = q__2.r + q__5.r, q__1.i = q__2.i + q__5.i;
			C(i,j).r = q__1.r, C(i,j).i = q__1.i;
		    }
/* L90: */
		}
/* L100: */
	    }
	}
    } else {

/*        Form  C := alpha*B*A + beta*C. */

	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
	    i__2 = j + j * a_dim1;
	    d__1 = A(j,j).r;
	    q__1.r = d__1 * alpha->r, q__1.i = d__1 * alpha->i;
	    temp1.r = q__1.r, temp1.i = q__1.i;
	    if (beta->r == 0.f && beta->i == 0.f) {
		i__2 = *m;
		for (i = 1; i <= *m; ++i) {
		    i__3 = i + j * c_dim1;
		    i__4 = i + j * b_dim1;
		    q__1.r = temp1.r * B(i,j).r - temp1.i * B(i,j).i, 
			    q__1.i = temp1.r * B(i,j).i + temp1.i * B(i,j)
			    .r;
		    C(i,j).r = q__1.r, C(i,j).i = q__1.i;
/* L110: */
		}
	    } else {
		i__2 = *m;
		for (i = 1; i <= *m; ++i) {
		    i__3 = i + j * c_dim1;
		    i__4 = i + j * c_dim1;
		    q__2.r = beta->r * C(i,j).r - beta->i * C(i,j).i, 
			    q__2.i = beta->r * C(i,j).i + beta->i * C(i,j)
			    .r;
		    i__5 = i + j * b_dim1;
		    q__3.r = temp1.r * B(i,j).r - temp1.i * B(i,j).i, 
			    q__3.i = temp1.r * B(i,j).i + temp1.i * B(i,j)
			    .r;
		    q__1.r = q__2.r + q__3.r, q__1.i = q__2.i + q__3.i;
		    C(i,j).r = q__1.r, C(i,j).i = q__1.i;
/* L120: */
		}
	    }
	    i__2 = j - 1;
	    for (k = 1; k <= j-1; ++k) {
		if (upper) {
		    i__3 = k + j * a_dim1;
		    q__1.r = alpha->r * A(k,j).r - alpha->i * A(k,j).i, 
			    q__1.i = alpha->r * A(k,j).i + alpha->i * A(k,j)
			    .r;
		    temp1.r = q__1.r, temp1.i = q__1.i;
		} else {
		    r_cnjg(&q__2, &A(j,k));
		    q__1.r = alpha->r * q__2.r - alpha->i * q__2.i, q__1.i = 
			    alpha->r * q__2.i + alpha->i * q__2.r;
		    temp1.r = q__1.r, temp1.i = q__1.i;
		}
		i__3 = *m;
		for (i = 1; i <= *m; ++i) {
		    i__4 = i + j * c_dim1;
		    i__5 = i + j * c_dim1;
		    i__6 = i + k * b_dim1;
		    q__2.r = temp1.r * B(i,k).r - temp1.i * B(i,k).i, 
			    q__2.i = temp1.r * B(i,k).i + temp1.i * B(i,k)
			    .r;
		    q__1.r = C(i,j).r + q__2.r, q__1.i = C(i,j).i + q__2.i;
		    C(i,j).r = q__1.r, C(i,j).i = q__1.i;
/* L130: */
		}
/* L140: */
	    }
	    i__2 = *n;
	    for (k = j + 1; k <= *n; ++k) {
		if (upper) {
		    r_cnjg(&q__2, &A(j,k));
		    q__1.r = alpha->r * q__2.r - alpha->i * q__2.i, q__1.i = 
			    alpha->r * q__2.i + alpha->i * q__2.r;
		    temp1.r = q__1.r, temp1.i = q__1.i;
		} else {
		    i__3 = k + j * a_dim1;
		    q__1.r = alpha->r * A(k,j).r - alpha->i * A(k,j).i, 
			    q__1.i = alpha->r * A(k,j).i + alpha->i * A(k,j)
			    .r;
		    temp1.r = q__1.r, temp1.i = q__1.i;
		}
		i__3 = *m;
		for (i = 1; i <= *m; ++i) {
		    i__4 = i + j * c_dim1;
		    i__5 = i + j * c_dim1;
		    i__6 = i + k * b_dim1;
		    q__2.r = temp1.r * B(i,k).r - temp1.i * B(i,k).i, 
			    q__2.i = temp1.r * B(i,k).i + temp1.i * B(i,k)
			    .r;
		    q__1.r = C(i,j).r + q__2.r, q__1.i = C(i,j).i + q__2.i;
		    C(i,j).r = q__1.r, C(i,j).i = q__1.i;
/* L150: */
		}
/* L160: */
	    }
/* L170: */
	}
    }

    return 0;

/*     End of CHEMM . */

} /* chemm_ */

