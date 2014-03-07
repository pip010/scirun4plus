#include "f2c.h"

/* Subroutine */ int zlarft_(char *direct, char *storev, integer *n, integer *
	k, doublecomplex *v, integer *ldv, doublecomplex *tau, doublecomplex *
	t, integer *ldt)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    ZLARFT forms the triangular factor T of a complex block reflector H   
    of order n, which is defined as a product of k elementary reflectors. 
  

    If DIRECT = 'F', H = H(1) H(2) . . . H(k) and T is upper triangular; 
  

    If DIRECT = 'B', H = H(k) . . . H(2) H(1) and T is lower triangular. 
  

    If STOREV = 'C', the vector which defines the elementary reflector   
    H(i) is stored in the i-th column of the array V, and   

       H  =  I - V * T * V'   

    If STOREV = 'R', the vector which defines the elementary reflector   
    H(i) is stored in the i-th row of the array V, and   

       H  =  I - V' * T * V   

    Arguments   
    =========   

    DIRECT  (input) CHARACTER*1   
            Specifies the order in which the elementary reflectors are   
            multiplied to form the block reflector:   
            = 'F': H = H(1) H(2) . . . H(k) (Forward)   
            = 'B': H = H(k) . . . H(2) H(1) (Backward)   

    STOREV  (input) CHARACTER*1   
            Specifies how the vectors which define the elementary   
            reflectors are stored (see also Further Details):   
            = 'C': columnwise   
            = 'R': rowwise   

    N       (input) INTEGER   
            The order of the block reflector H. N >= 0.   

    K       (input) INTEGER   
            The order of the triangular factor T (= the number of   
            elementary reflectors). K >= 1.   

    V       (input/output) COMPLEX*16 array, dimension   
                                 (LDV,K) if STOREV = 'C'   
                                 (LDV,N) if STOREV = 'R'   
            The matrix V. See further details.   

    LDV     (input) INTEGER   
            The leading dimension of the array V.   
            If STOREV = 'C', LDV >= max(1,N); if STOREV = 'R', LDV >= K. 
  

    TAU     (input) COMPLEX*16 array, dimension (K)   
            TAU(i) must contain the scalar factor of the elementary   
            reflector H(i).   

    T       (output) COMPLEX*16 array, dimension (LDT,K)   
            The k by k triangular factor T of the block reflector.   
            If DIRECT = 'F', T is upper triangular; if DIRECT = 'B', T is 
  
            lower triangular. The rest of the array is not used.   

    LDT     (input) INTEGER   
            The leading dimension of the array T. LDT >= K.   

    Further Details   
    ===============   

    The shape of the matrix V and the storage of the vectors which define 
  
    the H(i) is best illustrated by the following example with n = 5 and 
  
    k = 3. The elements equal to 1 are not stored; the corresponding   
    array elements are modified but restored on exit. The rest of the   
    array is not used.   

    DIRECT = 'F' and STOREV = 'C':         DIRECT = 'F' and STOREV = 'R': 
  

                 V = (  1       )                 V = (  1 v1 v1 v1 v1 ) 
  
                     ( v1  1    )                     (     1 v2 v2 v2 ) 
  
                     ( v1 v2  1 )                     (        1 v3 v3 ) 
  
                     ( v1 v2 v3 )   
                     ( v1 v2 v3 )   

    DIRECT = 'B' and STOREV = 'C':         DIRECT = 'B' and STOREV = 'R': 
  

                 V = ( v1 v2 v3 )                 V = ( v1 v1  1       ) 
  
                     ( v1 v2 v3 )                     ( v2 v2 v2  1    ) 
  
                     (  1 v2 v3 )                     ( v3 v3 v3 v3  1 ) 
  
                     (     1 v3 )   
                     (        1 )   

    ===================================================================== 
  


       Quick return if possible   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static doublecomplex c_b2 = {0.,0.};
    static integer c__1 = 1;
    
    /* System generated locals */
    integer t_dim1, t_offset, v_dim1, v_offset, i__1, i__2, i__3, i__4;
    doublecomplex z__1;
    /* Local variables */
    static integer i, j;
    extern logical lsame_(char *, char *);
    extern /* Subroutine */ int zgemv_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, doublecomplex *, integer *), 
	    ztrmv_(char *, char *, char *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, integer *), 
	    zlacgv_(integer *, doublecomplex *, integer *);
    static doublecomplex vii;



#define TAU(I) tau[(I)-1]

#define V(I,J) v[(I)-1 + ((J)-1)* ( *ldv)]
#define T(I,J) t[(I)-1 + ((J)-1)* ( *ldt)]

    if (*n == 0) {
	return 0;
    }

    if (lsame_(direct, "F")) {
	i__1 = *k;
	for (i = 1; i <= *k; ++i) {
	    i__2 = i;
	    if (TAU(i).r == 0. && TAU(i).i == 0.) {

/*              H(i)  =  I */

		i__2 = i;
		for (j = 1; j <= i; ++j) {
		    i__3 = j + i * t_dim1;
		    T(j,i).r = 0., T(j,i).i = 0.;
/* L10: */
		}
	    } else {

/*              general case */

		i__2 = i + i * v_dim1;
		vii.r = V(i,i).r, vii.i = V(i,i).i;
		i__2 = i + i * v_dim1;
		V(i,i).r = 1., V(i,i).i = 0.;
		if (lsame_(storev, "C")) {

/*                 T(1:i-1,i) := - tau(i) * V(i:n,1:i-1)' 
* V(i:n,i) */

		    i__2 = *n - i + 1;
		    i__3 = i - 1;
		    i__4 = i;
		    z__1.r = -TAU(i).r, z__1.i = -TAU(i).i;
		    zgemv_("Conjugate transpose", &i__2, &i__3, &z__1, &V(i,1), ldv, &V(i,i), &c__1, &c_b2, &
			    T(1,i), &c__1);
		} else {

/*                 T(1:i-1,i) := - tau(i) * V(1:i-1,i:n) *
 V(i,i:n)' */

		    if (i < *n) {
			i__2 = *n - i;
			zlacgv_(&i__2, &V(i,i+1), ldv);
		    }
		    i__2 = i - 1;
		    i__3 = *n - i + 1;
		    i__4 = i;
		    z__1.r = -TAU(i).r, z__1.i = -TAU(i).i;
		    zgemv_("No transpose", &i__2, &i__3, &z__1, &V(1,i), ldv, &V(i,i), ldv, &c_b2, &T(1,i), &c__1);
		    if (i < *n) {
			i__2 = *n - i;
			zlacgv_(&i__2, &V(i,i+1), ldv);
		    }
		}
		i__2 = i + i * v_dim1;
		V(i,i).r = vii.r, V(i,i).i = vii.i;

/*              T(1:i-1,i) := T(1:i-1,1:i-1) * T(1:i-1,i) */

		i__2 = i - 1;
		ztrmv_("Upper", "No transpose", "Non-unit", &i__2, &T(1,1), ldt, &T(1,i), &c__1);
		i__2 = i + i * t_dim1;
		i__3 = i;
		T(i,i).r = TAU(i).r, T(i,i).i = TAU(i).i;
	    }
/* L20: */
	}
    } else {
	for (i = *k; i >= 1; --i) {
	    i__1 = i;
	    if (TAU(i).r == 0. && TAU(i).i == 0.) {

/*              H(i)  =  I */

		i__1 = *k;
		for (j = i; j <= *k; ++j) {
		    i__2 = j + i * t_dim1;
		    T(j,i).r = 0., T(j,i).i = 0.;
/* L30: */
		}
	    } else {

/*              general case */

		if (i < *k) {
		    if (lsame_(storev, "C")) {
			i__1 = *n - *k + i + i * v_dim1;
			vii.r = V(*n-*k+i,i).r, vii.i = V(*n-*k+i,i).i;
			i__1 = *n - *k + i + i * v_dim1;
			V(*n-*k+i,i).r = 1., V(*n-*k+i,i).i = 0.;

/*                    T(i+1:k,i) :=   
                              - tau(i) * V(1:n-k+i,i+1
:k)' * V(1:n-k+i,i) */

			i__1 = *n - *k + i;
			i__2 = *k - i;
			i__3 = i;
			z__1.r = -TAU(i).r, z__1.i = -TAU(i).i;
			zgemv_("Conjugate transpose", &i__1, &i__2, &z__1, &V(1,i+1), ldv, &V(1,i)
				, &c__1, &c_b2, &T(i+1,i), &c__1);
			i__1 = *n - *k + i + i * v_dim1;
			V(*n-*k+i,i).r = vii.r, V(*n-*k+i,i).i = vii.i;
		    } else {
			i__1 = i + (*n - *k + i) * v_dim1;
			vii.r = V(i,*n-*k+i).r, vii.i = V(i,*n-*k+i).i;
			i__1 = i + (*n - *k + i) * v_dim1;
			V(i,*n-*k+i).r = 1., V(i,*n-*k+i).i = 0.;

/*                    T(i+1:k,i) :=   
                              - tau(i) * V(i+1:k,1:n-k
+i) * V(i,1:n-k+i)' */

			i__1 = *n - *k + i - 1;
			zlacgv_(&i__1, &V(i,1), ldv);
			i__1 = *k - i;
			i__2 = *n - *k + i;
			i__3 = i;
			z__1.r = -TAU(i).r, z__1.i = -TAU(i).i;
			zgemv_("No transpose", &i__1, &i__2, &z__1, &V(i+1,1), ldv, &V(i,1), ldv, &c_b2, &
				T(i+1,i), &c__1);
			i__1 = *n - *k + i - 1;
			zlacgv_(&i__1, &V(i,1), ldv);
			i__1 = i + (*n - *k + i) * v_dim1;
			V(i,*n-*k+i).r = vii.r, V(i,*n-*k+i).i = vii.i;
		    }

/*                 T(i+1:k,i) := T(i+1:k,i+1:k) * T(i+1:k,
i) */

		    i__1 = *k - i;
		    ztrmv_("Lower", "No transpose", "Non-unit", &i__1, &T(i+1,i+1), ldt, &T(i+1,i)
			    , &c__1);
		}
		i__1 = i + i * t_dim1;
		i__2 = i;
		T(i,i).r = TAU(i).r, T(i,i).i = TAU(i).i;
	    }
/* L40: */
	}
    }
    return 0;

/*     End of ZLARFT */

} /* zlarft_ */

