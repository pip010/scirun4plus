#ifndef __MTXLIB_H__
#define __MTXLIB_H__

#include <iostream>

static inline float DegToRad(float a) { return a*0.01745329252f;};
static inline float RadToDeg(float a) { return a*57.29577951f;};


////////////////////////////////////////////////////////////////////////
//                           VECTOR CLASS                             //
////////////////////////////////////////////////////////////////////////


template <int dim> class vec 
{
public:
  //--------------------------------------------------------------------
  // constructors and assignment
  vec() {};
  vec(float rhs);
  vec(float inX, float inY);
  vec(float inX, float inY, float inZ);
  vec(float inX, float inY, float inZ, float inW);
  vec(const vec &rhs) { *this = rhs; };

  vec& operator = (const vec &rhs);
  vec& operator = (float rhs);

  void set(float inX, float inY); 
  void set(float inX, float inY, float inZ);
  void set(float inX, float inY, float inZ, float inW);

  //--------------------------------------------------------------------
  // array indexing
  float& operator [] (unsigned int i) { return x[i]; };
  const float& operator [] (unsigned int i) const { return x[i]; };
  float& operator () (unsigned int i) { return x[i]; };
  const float& operator () (unsigned int i) const { return x[i]; };

  //--------------------------------------------------------------------
  // basic arithmetic
  vec& operator += (const vec &rhs);
  vec& operator -= (const vec &rhs);
  vec& operator *= (float rhs);
  vec& operator *= (const vec &rhs);
  vec& operator /= (float rhs);
  vec& operator /= (const vec &rhs);

  //--------------------------------------------------------------------
  // vector methods
  float length(float epsilon=0.0) const;
  float lengthSqr() const;
  vec& normalize();
  vec& random();


  float  x[dim];
};

//----------------------------------------------------------------------
// basic arithmetic

template <int dim>
vec<dim> operator + (const vec<dim> &a, const vec<dim> &b); 

template <int dim>
vec<dim> operator - (const vec<dim> &a, const vec<dim> &b);

template <int dim>
vec<dim> operator - (const vec<dim> &a);

template <int dim>
vec<dim> operator * (const vec<dim> &v, float f); 

template <int dim>
vec<dim> operator * (float f, const vec<dim> &v);

template <int dim>
vec<dim> operator * (const vec<dim> &a, const vec<dim> &b); 

template <int dim>
vec<dim> operator / (const vec<dim> &v, float f); 

template <int dim>
vec<dim> operator / (const vec<dim> &a, const vec<dim> &b); 

//----------------------------------------------------------------------
// equal ops

template <int dim>
bool operator == (const vec<dim> &a, const vec<dim> &b);

template <int dim>
bool operator != (const vec<dim> &a, const vec<dim> &b);


//----------------------------------------------------------------------
// vector methods
template <int dim>
float DotProduct(const vec<dim> &a, const vec<dim> &b);

vec<3> CrossProduct(const vec<3> &a, const vec<3> &b);


template <int dim>
std::ostream& operator << (std::ostream& outs, 
                           const vec<dim>& source);



////////////////////////////////////////////////////////////////////////
//                           MATRIX CLASS                             //
////////////////////////////////////////////////////////////////////////


template <int r, int c> class matrix 
{
public:
  //--------------------------------------------------------------------
  // constructors and assignment
  matrix() {};
  matrix(float v);
  matrix(const vec<2> &r0, const vec<2> &r1);
  matrix(const vec<3> &r0, const vec<3> &r1, const vec<3> &r2);
  matrix(const vec<3> &r0, const vec<3> &r1, 
         const vec<3> &r2, const vec<3> &r3);
  matrix(const vec<4> &r0, const vec<4> &r1, 
         const vec<4> &r2, const vec<4> &r3);
  matrix(const matrix &rhs) { *this = rhs; };

  matrix& operator = (const matrix &rhs);
  matrix& operator = (float rhs);

  void set(const vec<2> &r0, const vec<2> &r1); 
  void set(const vec<3> &r0, const vec<3> &r1, const vec<3> &r2);
  void set(const vec<3> &r0, const vec<3> &r1, 
           const vec<3> &r2, const vec<3> &r3);
  void set(const vec<2> &r0, const vec<2> &r1, 
           const vec<2> &r2, const vec<2> &r3);
  void set(const vec<4> &r0, const vec<4> &r1, 
           const vec<4> &r2, const vec<4> &r3);

  //--------------------------------------------------------------------
  // array indexing
  float& operator () (unsigned int i, unsigned int j) 
  { return x[i*c+j]; };
  const float& operator () (unsigned int i, unsigned int j) const 
  { return x[i*c+j]; };

  //--------------------------------------------------------------------
  // basic arithmetic
  matrix& operator += (const matrix &rhs);
  matrix& operator -= (const matrix &rhs);
  matrix& operator *= (float rhs);
  matrix& operator /= (float rhs);

  //--------------------------------------------------------------------
  // vector methods
  matrix& identity();
  matrix& transpose();
  matrix& invert();
  float trace() const;
  float norm() const;
 
  float x[r*c];
};

//----------------------------------------------------------------------
// basic arithmetic

template <int r, int c>
matrix<r,c> operator + (const matrix<r,c> &a, const matrix<r,c> &b); 

template <int r, int c>
matrix<r,c> operator - (const matrix<r,c> &a, const matrix<r,c> &b);

template <int r, int c>
matrix<r,c> operator - (const matrix<r,c> &a);

//template <int r, int c, int c2>
//matrix<r,c2> operator * (const matrix<r,c> &a, const matrix<c,c2> &b );

template <int r, int c>
matrix<r,r> operator * (const matrix<r,c> &a, const matrix<c,r> &b );

template <int r, int c>
matrix<r,c> operator * (const matrix<r,c> &a, float f); 

template <int r, int c>
matrix<r,c> operator * (float f, const matrix<r,c> &a);

template <int r, int c>
vec<r> operator * (const matrix<r,c> &m, const vec<c> &v);

template <int r, int c>
vec<c> operator * (const vec<r> &v, const matrix<r,c> &m);

template <int r, int c>
matrix<r,c> operator / (const matrix<r,c> &a, float f);


//----------------------------------------------------------------------
// equal ops

template <int r, int c>
bool operator == (const matrix<r,c> &a, const matrix<r,c> &b);

template <int r, int c>
bool operator != (const matrix<r,c> &a, const matrix<r,c> &b);


template <int r, int c>
std::ostream& operator << (std::ostream& outs, 
                           const matrix<r,c>& source);


//----------------------------------------------------------------------
// matrix ops

template <int r, int c>
matrix<r,c> DirectProduct(const vec<r> &v, const vec<c> &vT);


////////////////////////////////////////////////////////////////////////
//                         QUATERNION CLASS                           //
////////////////////////////////////////////////////////////////////////

class quaternion
{
 public:
	float x,y,z,w;
	quaternion(){ init(0.0,0.0,0.0,1.0); };

  inline void init(float xv, float yv, float zv, float wv) 
  { x=xv; y=yv; z=zv; w=wv; };

public:
  friend quaternion operator*(const quaternion &a,
                              const quaternion &b) 
  {
	  quaternion c;
	  c.w = a.w*b.w - a.x*b.x - a.y*b.y - a.z*b.z; 
	  c.x = a.w*b.x + a.x*b.w + a.y*b.z - a.z*b.y; 
	  c.y = a.w*b.y - a.x*b.z + a.y*b.w + a.z*b.x; 
	  c.z = a.w*b.z + a.x*b.y - a.y*b.x + a.z*b.w; 
	  return c;
  };
};

quaternion RotationArc(const vec<3> &v0, const vec<3> &v1);
void       MakeMatrix44(float *m, const vec<3> &v,
                        const quaternion &q);
quaternion VirtualTrackBall(const vec<3> &cop,  const vec<3> &cor,
                            const vec<3> &dir1, const vec<3> &dir2);

vec<3> PlaneLineIntersection(const vec<3> &n, float d,
                             const vec<3> &p1, const vec<3> &p2);

std::ostream& operator << (std::ostream& outs, 
                           const quaternion& source);


#include "mtxlib.T"

#endif // __MTXLIB_H__
