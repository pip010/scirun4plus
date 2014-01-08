#include "f2c.h"

/* Subroutine */ int zggsvd_(char *jobu, char *jobv, char *jobq, integer *m, 
	integer *n, integer *p, integer *k, integer *l, doublecomplex *a, 
	integer *lda, doublecomplex *b, integer *ldb, doublereal *alpha, 
	doublereal *beta, doublecomplex *u, integer *ldu, doublecomplex *v, 
	integer *ldv, doublecomplex *q, integer *ldq, doublecomplex *work, 
	doublereal *rwork, integer *iwork, integer *info)
{
/*  -- LAPACK driver routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    ZGGSVD computes the generalized singular value decomposition (GSVD)   
    of an M-by-N complex matrix A and P-by-N complex matrix B:   

          U'*A*Q = D1*( 0 R ),    V'*B*Q = D2*( 0 R )   

    where U, V and Q are unitary matrices, and Z' means the conjugate   
    transpose of Z.  Let K+L = the effective numerical rank of the   
    matrix (A',B')', then R is a (K+L)-by-(K+L) nonsingular upper   
    triangular matrix, D1 and D2 are M-by-(K+L) and P-by-(K+L) "diagonal" 
  
    matrices and of the following structures, respectively:   

    If M-K-L >= 0,   

                        K  L   
           D1 =     K ( I  0 )   
                    L ( 0  C )   
                M-K-L ( 0  0 )   

                      K  L   
           D2 =   L ( 0  S )   
                P-L ( 0  0 )   

                    N-K-L  K    L   
      ( 0 R ) = K (  0   R11  R12 )   
                L (  0    0   R22 )   
    where   

      C = diag( ALPHA(K+1), ... , ALPHA(K+L) ),   
      S = diag( BETA(K+1),  ... , BETA(K+L) ),   
      C**2 + S**2 = I.   

      R is stored in A(1:K+L,N-K-L+1:N) on exit.   

    If M-K-L < 0,   

                      K M-K K+L-M   
           D1 =   K ( I  0    0   )   
                M-K ( 0  C    0   )   

                        K M-K K+L-M   
           D2 =   M-K ( 0  S    0  )   
                K+L-M ( 0  0    I  )   
                  P-L ( 0  0    0  )   

                       N-K-L  K   M-K  K+L-M   
      ( 0 R ) =     K ( 0    R11  R12  R13  )   
                  M-K ( 0     0   R22  R23  )   
                K+L-M ( 0     0    0   R33  )   

    where   

      C = diag( ALPHA(K+1), ... , ALPHA(M) ),   
      S = diag( BETA(K+1),  ... , BETA(M) ),   
      C**2 + S**2 = I.   

      (R11 R12 R13 ) is stored in A(1:M, N-K-L+1:N), and R33 is stored   
      ( 0  R22 R23 )   
      in B(M-K+1:L,N+M-K-L+1:N) on exit.   

    The routine computes C, S, R, and optionally the unitary   
    transformation matrices U, V and Q.   

    In particular, if B is an N-by-N nonsingular matrix, then the GSVD of 
  
    A and B implicitly gives the SVD of A*inv(B):   
                         A*inv(B) = U*(D1*inv(D2))*V'.   
    If ( A',B')' has orthnormal columns, then the GSVD of A and B is also 
  
    equal to the CS decomposition of A and B. Furthermore, the GSVD can   
    be used to derive the solution of the eigenvalue problem:   
                         A'*A x = lambda* B'*B x.   
    In some literature, the GSVD of A and B is presented in the form   
                     U'*A*X = ( 0 D1 ),   V'*B*X = ( 0 D2 )   
    where U and V are orthogonal and X is nonsingular, and D1 and D2 are 
  
    ``diagonal''.  The former GSVD form can be converted to the latter   
    form by taking the nonsingular matrix X as   

                          X = Q*(  I   0    )   
                                (  0 inv(R) )   

    Arguments   
    =========   

    JOBU    (input) CHARACTER*1   
            = 'U':  Unitary matrix U is computed;   
            = 'N':  U is not computed.   

    JOBV    (input) CHARACTER*1   
            = 'V':  Unitary matrix V is computed;   
            = 'N':  V is not computed.   

    JOBQ    (input) CHARACTER*1   
            = 'Q':  Unitary matrix Q is computed;   
            = 'N':  Q is not computed.   

    M       (input) INTEGER   
            The number of rows of the matrix A.  M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrices A and B.  N >= 0.   

    P       (input) INTEGER   
            The number of rows of the matrix B.  P >= 0.   

    K       (output) INTEGER   
    L       (output) INTEGER   
            On exit, K and L specify the dimension of the subblocks   
            described in Purpose.   
            K + L = effective numerical rank of (A',B')'.   

    A       (input/output) COMPLEX*16 array, dimension (LDA,N)   
            On entry, the M-by-N matrix A.   
            On exit, A contains the triangular matrix R, or part of R.   
            See Purpose for details.   

    LDA     (input) INTEGER   
            The leading dimension of the array A. LDA >= max(1,M).   

    B       (input/output) COMPLEX*16 array, dimension (LDB,N)   
            On entry, the P-by-N matrix B.   
            On exit, B contains part of the triangular matrix R if   
            M-K-L < 0.  See Purpose for details.   

    LDB     (input) INTEGER   
            The leading dimension of the array B. LDB >= max(1,P).   

    ALPHA   (output) DOUBLE PRECISION array, dimension (N)   
    BETA    (output) DOUBLE PRECISION array, dimension (N)   
            On exit, ALPHA and BETA contain the generalized singular   
            value pairs of A and B;   
              ALPHA(1:K) = 1,   
              BETA(1:K)  = 0,   
            and if M-K-L >= 0,   
              ALPHA(K+1:K+L) = C,   
              BETA(K+1:K+L)  = S,   
            or if M-K-L < 0,   
              ALPHA(K+1:M)= C, ALPHA(M+1:K+L)= 0   
              BETA(K+1:M) = S, BETA(M+1:K+L) = 1   
            and   
              ALPHA(K+L+1:N) = 0   
              BETA(K+L+1:N)  = 0   

    U       (output) COMPLEX*16 array, dimension (LDU,M)   
            If JOBU = 'U', U contains the M-by-M unitary matrix U.   
            If JOBU = 'N', U is not referenced.   

    LDU     (input) INTEGER   
            The leading dimension of the array U. LDU >= max(1,M) if   
            JOBU = 'U'; LDU >= 1 otherwise.   

    V       (output) COMPLEX*16 array, dimension (LDV,P)   
            If JOBV = 'V', V contains the P-by-P unitary matrix V.   
            If JOBV = 'N', V is not referenced.   

    LDV     (input) INTEGER   
            The leading dimension of the array V. LDV >= max(1,P) if   
            JOBV = 'V'; LDV >= 1 otherwise.   

    Q       (output) COMPLEX*16 array, dimension (LDQ,N)   
            If JOBQ = 'Q', Q contains the N-by-N unitary matrix Q.   
            If JOBQ = 'N', Q is not referenced.   

    LDQ     (input) INTEGER   
            The leading dimension of the array Q. LDQ >= max(1,N) if   
            JOBQ = 'Q'; LDQ >= 1 otherwise.   

    WORK    (workspace) COMPLEX*16 array, dimension (max(3*N,M,P)+N)   

    RWORK   (workspace) DOUBLE PRECISION array, dimension (2*N)   

    IWORK   (workspace) INTEGER array, dimension (N)   

    INFO    (output)INTEGER   
            = 0:  successful exit.   
            < 0:  if INFO = -i, the i-th argument had an illegal value.   
            > 0:  if INFO = 1, the Jacobi-type procedure failed to   
                  converge.  For further details, see subroutine ZTGSJA. 
  

    Internal Parameters   
    ===================   

    TOLA    DOUBLE PRECISION   
    TOLB    DOUBLE PRECISION   
            TOLA and TOLB are the thresholds to determine the effective   
            rank of (A',B')'. Generally, they are set to   
                     TOLA = MAX(M,N)*norm(A)*MAZHEPS,   
                     TOLB = MAX(P,N)*norm(B)*MAZHEPS.   
            The size of TOLA and TOLB may affect the size of backward   
            errors of the decomposition.   

    ===================================================================== 
  


       Decode and test the input parameters   

    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, q_dim1, q_offset, u_dim1, 
	    u_offset, v_dim1, v_offset, i__1;
    /* Local variables */
    static doublereal tola, tolb, unfl;
    extern logical lsame_(char *, char *);
    static doublereal anorm, bnorm;
    static logical wantq, wantu, wantv;
    extern doublereal dlamch_(char *);
    static integer ncycle;
    extern /* Subroutine */ int xerbla_(char *, integer *);
    extern doublereal zlange_(char *, integer *, integer *, doublecomplex *, 
	    integer *, doublereal *);
    extern /* Subroutine */ int ztgsja_(char *, char *, char *, integer *, 
	    integer *, integer *, integer *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, integer *), 
	    zggsvp_(char *, char *, char *, integer *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublecomplex *,
	     integer *, doublecomplex *, integer *, doublecomplex *, integer *
	    , integer *, doublereal *, doublecomplex *, doublecomplex *, 
	    integer *);
    static doublereal ulp;


#define ALPHA(I) alpha[(I)-1]
#define BETA(I) beta[(I)-1]
#define WORK(I) work[(I)-1]
#define RWORK(I) rwork[(I)-1]
#define IWORK(I) iwork[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]
#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]
#define U(I,J) u[(I)-1 + ((J)-1)* ( *ldu)]
#define V(I,J) v[(I)-1 + ((J)-1)* ( *ldv)]
#define Q(I,J) q[(I)-1 + ((J)-1)* ( *ldq)]

    wantu = lsame_(jobu, "U");
    wantv = lsame_(jobv, "V");
    wantq = lsame_(jobq, "Q");

    *info = 0;
    if (! (wantu || lsame_(jobu, "N"))) {
	*info = -1;
    } else if (! (wantv || lsame_(jobv, "N"))) {
	*info = -2;
    } else if (! (wantq || lsame_(jobq, "N"))) {
	*info = -3;
    } else if (*m < 0) {
	*info = -4;
    } else if (*n < 0) {
	*info = -5;
    } else if (*p < 0) {
	*info = -6;
    } else if (*lda < max(1,*m)) {
	*info = -10;
    } else if (*ldb < max(1,*p)) {
	*info = -12;
    } else if (*ldu < 1 || wantu && *ldu < *m) {
	*info = -16;
    } else if (*ldv < 1 || wantv && *ldv < *p) {
	*info = -18;
    } else if (*ldq < 1 || wantq && *ldq < *n) {
	*info = -20;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("ZGGSVD", &i__1);
	return 0;
    }

/*     Compute the Frobenius norm of matrices A and B */

    anorm = zlange_("1", m, n, &A(1,1), lda, &RWORK(1));
    bnorm = zlange_("1", p, n, &B(1,1), ldb, &RWORK(1));

/*     Get machine precision and set up threshold for determining   
       the effective numerical rank of the matrices A and B. */

    ulp = dlamch_("Precision");
    unfl = dlamch_("Safe Minimum");
    tola = max(*m,*n) * max(anorm,unfl) * ulp;
    tolb = max(*p,*n) * max(bnorm,unfl) * ulp;

    zggsvp_(jobu, jobv, jobq, m, p, n, &A(1,1), lda, &B(1,1), ldb, &
	    tola, &tolb, k, l, &U(1,1), ldu, &V(1,1), ldv, &Q(1,1), ldq, &IWORK(1), &RWORK(1), &WORK(1), &WORK(*n + 1), 
	    info);

/*     Compute the GSVD of two upper "triangular" matrices */

    ztgsja_(jobu, jobv, jobq, m, p, n, k, l, &A(1,1), lda, &B(1,1), 
	    ldb, &tola, &tolb, &ALPHA(1), &BETA(1), &U(1,1), ldu, &V(1,1), ldv, &Q(1,1), ldq, &WORK(1), &ncycle, info);

    return 0;

/*     End of ZGGSVD */

} /* zggsvd_ */

