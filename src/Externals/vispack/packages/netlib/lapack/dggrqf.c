#include "f2c.h"

/* Subroutine */ int dggrqf_(integer *m, integer *p, integer *n, doublereal *
	a, integer *lda, doublereal *taua, doublereal *b, integer *ldb, 
	doublereal *taub, doublereal *work, integer *lwork, integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DGGRQF computes a generalized RQ factorization of an M-by-N matrix A 
  
    and a P-by-N matrix B:   

                A = R*Q,        B = Z*T*Q,   

    where Q is an N-by-N orthogonal matrix, Z is a P-by-P orthogonal   
    matrix, and R and T assume one of the forms:   

    if M <= N,  R = ( 0  R12 ) M,   or if M > N,  R = ( R11 ) M-N,   
                     N-M  M                           ( R21 ) N   
                                                         N   

    where R12 or R21 is upper triangular, and   

    if P >= N,  T = ( T11 ) N  ,   or if P < N,  T = ( T11  T12 ) P,   
                    (  0  ) P-N                         P   N-P   
                       N   

    where T11 is upper triangular.   

    In particular, if B is square and nonsingular, the GRQ factorization 
  
    of A and B implicitly gives the RQ factorization of A*inv(B):   

                 A*inv(B) = (R*inv(T))*Z'   

    where inv(B) denotes the inverse of the matrix B, and Z' denotes the 
  
    transpose of the matrix Z.   

    Arguments   
    =========   

    M       (input) INTEGER   
            The number of rows of the matrix A.  M >= 0.   

    P       (input) INTEGER   
            The number of rows of the matrix B.  P >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrices A and B. N >= 0.   

    A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the M-by-N matrix A.   
            On exit, if M <= N, the upper triangle of the subarray   
            A(1:M,N-M+1:N) contains the M-by-M upper triangular matrix R; 
  
            if M > N, the elements on and above the (M-N)-th subdiagonal 
  
            contain the M-by-N upper trapezoidal matrix R; the remaining 
  
            elements, with the array TAUA, represent the orthogonal   
            matrix Q as a product of elementary reflectors (see Further   
            Details).   

    LDA     (input) INTEGER   
            The leading dimension of the array A. LDA >= max(1,M).   

    TAUA    (output) DOUBLE PRECISION array, dimension (min(M,N))   
            The scalar factors of the elementary reflectors which   
            represent the orthogonal matrix Q (see Further Details).   

    B       (input/output) DOUBLE PRECISION array, dimension (LDB,N)   
            On entry, the P-by-N matrix B.   
            On exit, the elements on and above the diagonal of the array 
  
            contain the min(P,N)-by-N upper trapezoidal matrix T (T is   
            upper triangular if P >= N); the elements below the diagonal, 
  
            with the array TAUB, represent the orthogonal matrix Z as a   
            product of elementary reflectors (see Further Details).   

    LDB     (input) INTEGER   
            The leading dimension of the array B. LDB >= max(1,P).   

    TAUB    (output) DOUBLE PRECISION array, dimension (min(P,N))   
            The scalar factors of the elementary reflectors which   
            represent the orthogonal matrix Z (see Further Details).   

    WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK) 
  
            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INTEGER   
            The dimension of the array WORK. LWORK >= max(1,N,M,P).   
            For optimum performance LWORK >= max(N,M,P)*max(NB1,NB2,NB3), 
  
            where NB1 is the optimal blocksize for the RQ factorization   
            of an M-by-N matrix, NB2 is the optimal blocksize for the   
            QR factorization of a P-by-N matrix, and NB3 is the optimal   
            blocksize for a call of DORMRQ.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INF0= -i, the i-th argument had an illegal value.   

    Further Details   
    ===============   

    The matrix Q is represented as a product of elementary reflectors   

       Q = H(1) H(2) . . . H(k), where k = min(m,n).   

    Each H(i) has the form   

       H(i) = I - taua * v * v'   

    where taua is a real scalar, and v is a real vector with   
    v(n-k+i+1:n) = 0 and v(n-k+i) = 1; v(1:n-k+i-1) is stored on exit in 
  
    A(m-k+i,1:n-k+i-1), and taua in TAUA(i).   
    To form Q explicitly, use LAPACK subroutine DORGRQ.   
    To use Q to update another matrix, use LAPACK subroutine DORMRQ.   

    The matrix Z is represented as a product of elementary reflectors   

       Z = H(1) H(2) . . . H(k), where k = min(p,n).   

    Each H(i) has the form   

       H(i) = I - taub * v * v'   

    where taub is a real scalar, and v is a real vector with   
    v(1:i-1) = 0 and v(i) = 1; v(i+1:p) is stored on exit in B(i+1:p,i), 
  
    and taub in TAUB(i).   
    To form Z explicitly, use LAPACK subroutine DORGQR.   
    To use Z to update another matrix, use LAPACK subroutine DORMQR.   

    ===================================================================== 
  


       Test the input parameters   

    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3;
    /* Local variables */
    static integer lopt;
    extern /* Subroutine */ int dgeqrf_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, integer *), 
	    dgerqf_(integer *, integer *, doublereal *, integer *, doublereal 
	    *, doublereal *, integer *, integer *), xerbla_(char *, integer *), dormrq_(char *, char *, integer *, integer *, integer *,
	     doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, integer *);


#define TAUA(I) taua[(I)-1]
#define TAUB(I) taub[(I)-1]
#define WORK(I) work[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]
#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]

    *info = 0;
    if (*m < 0) {
	*info = -1;
    } else if (*p < 0) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*lda < max(1,*m)) {
	*info = -5;
    } else if (*ldb < max(1,*p)) {
	*info = -8;
    } else /* if(complicated condition) */ {
/* Computing MAX */
	i__1 = max(1,*m), i__1 = max(i__1,*p);
	if (*lwork < max(i__1,*n)) {
	    *info = -11;
	}
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DGGRQF", &i__1);
	return 0;
    }

/*     RQ factorization of M-by-N matrix A: A = R*Q */

    dgerqf_(m, n, &A(1,1), lda, &TAUA(1), &WORK(1), lwork, info);
    lopt = (integer) WORK(1);

/*     Update B := B*Q' */

    i__1 = min(*m,*n);
/* Computing MAX */
    i__2 = 1, i__3 = *m - *n + 1;
    dormrq_("Right", "Transpose", p, n, &i__1, &A(max(1,*m-*n+1),1), 
	    lda, &TAUA(1), &B(1,1), ldb, &WORK(1), lwork, info);
/* Computing MAX */
    i__1 = lopt, i__2 = (integer) WORK(1);
    lopt = max(i__1,i__2);

/*     QR factorization of P-by-N matrix B: B = Z*T */

    dgeqrf_(p, n, &B(1,1), ldb, &TAUB(1), &WORK(1), lwork, info);
/* Computing MAX */
    i__1 = lopt, i__2 = (integer) WORK(1);
    WORK(1) = (doublereal) max(i__1,i__2);

    return 0;

/*     End of DGGRQF */

} /* dggrqf_ */

