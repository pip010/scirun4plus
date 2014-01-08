#ifndef __DEFINES_H__
#define __DEFINES_H__

#ifdef _WIN32
#pragma warning( disable: 4244 )
#pragma warning( disable: 4251 )
#endif

#define EPSILON 1.0e-6
#define PI 3.141592
#define RADIANS2DEGREES 180.0/PI
#define DEGREES2RADIANS PI/180.0

#define SUPER_EPSILON 1.0e-10

#define INITIAL_LAMBDA 0.1
#define MIN_LAMBDA     1.0e-15
#define MAX_LAMBDA     1.0e15

#define MAX_VALUE      3.402823466e+38F

#define sign(a) ((a>=0) ? 1.0 : -1.0)

//-----------------------------------------
// for defining 2D or 3D
//-----------------------------------------

#define NORMALIZED_IDEAL_DISTANCE_6 (0.5/cos(30.0*DEGREES2RADIANS))
#define NORMALIZED_IDEAL_DISTANCE_2 0.5

//#define TWO_D
#define THREE_D

#ifdef TWO_D
#define DIMENSION 2
#define NUM_NEIGHBORS 2

#endif

#ifdef THREE_D
#define DIMENSION 3
#define NUM_NEIGHBORS 6// back to six

#endif

#include <features/mtxlib.h>

typedef vec<DIMENSION> vector_type;
typedef matrix<DIMENSION,DIMENSION> matrix_type;

enum { MPU, MLS };
enum { QUADRATIC, QUADRATIC_LINEAR, QUADRATIC_CONSTANT,
CUBIC_QUADRATIC_LINEAR };


#endif // __DEFINES_H__
