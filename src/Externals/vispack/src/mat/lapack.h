#ifndef __LAPACK_H
#define __LAPACK_H

#include "matrix.h"
#include "f2c.h"


//Declarations for all lapack functions to be used
#ifndef VISMATRIXDOUBLE
extern "C" int sgesvd_(char *, char *, integer *, integer *, real *,\
		       integer *, real *, real *, integer *, real *,\
		       integer *, real *, integer *, integer *); //SVD
extern "C" int slacpy_(char *, integer *, integer *, real *, integer *,\
		       real *, integer *); //LAPACK copy for QR
extern "C" int sorgqr_(integer *, integer *, integer *, real *, integer *,\
		       real *, real *, integer *, integer *); //make Q for QR
extern "C" int sgeqrf_(integer *, integer *, real *, integer *, real *,\
		       real *, integer *, integer *); //QR
extern "C" int sgetrf_(integer *, integer *, real *, integer *, integer *,\
		       integer *); //LU
extern "C" int sgesv_(integer *, integer *, real *, integer *, integer *,\
		      real *, integer *, integer *); //solves linear equations
#else
extern "C" int dgesvd_(char *, char *, integer *, integer *, doublereal *,\
		       integer *, doublereal *, doublereal *, integer *, doublereal *,\
		       integer *, doublereal *, integer *, integer *); //SVD
extern "C" int dlacpy_(char *, integer *, integer *, doublereal *, integer *,\
		       doublereal *, integer *); //LAPACK copy for QR
extern "C" int dorgqr_(integer *, integer *, integer *, doublereal *, integer *,\
		       doublereal *, doublereal *, integer *, integer *); //make Q for QR
extern "C" int dgeqrf_(integer *, integer *, doublereal *, integer *, doublereal *,\
		       doublereal *, integer *, integer *); //QR
extern "C" int dgetrf_(integer *, integer *, doublereal *, integer *, integer *,\
		       integer *); //LU
extern "C" int dgesv_(integer *, integer *, doublereal *, integer *, integer *,\
		      doublereal *, integer *, integer *); //solves linear equations
#endif //VISMATRIXDOUBLE

#endif





