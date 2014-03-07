#include "nrutil.h"
// you will need to put a declaration for each c routine you plan to use
matrix_SHARE void gaussj(float **a, int n, float **b, int m);
matrix_SHARE void svdcmp(float **a, int m, int n, float w[], float **v);
matrix_SHARE float pythag(float a, float b);
matrix_SHARE void ludcmp(float **a, int n, int *indx, float *d);
matrix_SHARE void qrdcmp(float **a, int n, float *c, float *d, int *sing);
