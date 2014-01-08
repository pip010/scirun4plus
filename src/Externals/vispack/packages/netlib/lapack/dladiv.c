#include "f2c.h"

/* Subroutine */ int dladiv_(doublereal *a, doublereal *b, doublereal *c, 
	doublereal *d, doublereal *p, doublereal *q)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    DLADIV performs complex division in  real arithmetic   

                          a + i*b   
               p + i*q = ---------   
                          c + i*d   

    The algorithm is due to Robert L. Smith and can be found   
    in D. Knuth, The art of Computer Programming, Vol.2, p.195   

    Arguments   
    =========   

    A       (input) DOUBLE PRECISION   
    B       (input) DOUBLE PRECISION   
    C       (input) DOUBLE PRECISION   
    D       (input) DOUBLE PRECISION   
            The scalars a, b, c, and d in the above expression.   

    P       (output) DOUBLE PRECISION   
    Q       (output) DOUBLE PRECISION   
            The scalars p and q in the above expression.   

    ===================================================================== 
*/
    static doublereal e, f;



    if (abs(*d) < abs(*c)) {
	e = *d / *c;
	f = *c + *d * e;
	*p = (*a + *b * e) / f;
	*q = (*b - *a * e) / f;
    } else {
	e = *c / *d;
	f = *d + *c * e;
	*p = (*b + *a * e) / f;
	*q = (-(*a) + *b * e) / f;
    }

    return 0;

/*     End of DLADIV */

} /* dladiv_ */

