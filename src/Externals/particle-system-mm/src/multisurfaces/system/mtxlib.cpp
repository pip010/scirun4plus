#include <cmath>
#include <features/mtxlib.h>

#include <iostream>

#ifdef _WIN32
#pragma warning( disable: 4244 )
#endif


//-------------------------------------------------------------------------
// function   : Constructors
// description: 
//-------------------------------------------------------------------------
template<>
vec<2>::vec( float inX, float inY )
{
  x[0] = inX; x[1] = inY;    
}

template<>
vec<3>::vec( float inX, float inY, float inZ )
{
  x[0] = inX; x[1] = inY; x[2] = inZ; 
}

template<>
vec<4>::vec( float inX, float inY, float inZ, float inW )
{
  x[0] = inX; x[1] = inY; x[2] = inZ; x[3] = inW; 
}

template<>
void vec<2>::set( float inX, float inY )
{
  x[0] = inX; x[1] = inY;    
}

template<>
void vec<3>::set( float inX, float inY, float inZ )
{
  x[0] = inX; x[1] = inY; x[2] = inZ;  
}

template<>
void vec<4>::set( float inX, float inY, float inZ, float inW )
{
  x[0] = inX; x[1] = inY; x[2] = inZ; x[3] = inW;  
}

//-------------------------------------------------------------------------
// function   : vector methods
// description: 
//-------------------------------------------------------------------------

vec<3> CrossProduct( const vec<3> &a, const vec<3> &b )
{
  return vec<3>( a[1]*b[2] - a[2]*b[1],
                 a[2]*b[0] - a[0]*b[2],
                 a[0]*b[1] - a[1]*b[0] );
}

template<>
float DotProduct(const vec<4> &a, const vec<4> &b)
{
  return a(0)*b(0) + a(1)*b(1) + a(2)*b(2) + a(3)*b(3);
}

template<>
float DotProduct(const vec<3> &a, const vec<3> &b)
{
return a(0)*b(0) + a(1)*b(1) + a(2)*b(2);
}

/************************************************************************/


//-------------------------------------------------------------------------
// function   : Constructors
// description: 
//-------------------------------------------------------------------------

template<>
matrix<2,2>::matrix( const vec<2> &r0, const vec<2> &r1 )
{ 
  x[0] = r0[0]; x[1] = r0[1];
  x[2] = r1[0]; x[3] = r1[1];
} 

template<>
matrix<3,3>::matrix( const vec<3> &r0, const vec<3> &r1, const vec<3> &r2 )
{
  x[0] = r0[0]; x[1] = r0[1]; x[2] = r0[2];
  x[3] = r1[0]; x[4] = r1[1]; x[5] = r1[2];
  x[6] = r2[0]; x[7] = r2[1]; x[8] = r2[2];
} 

template<>
matrix<4,3>::matrix( const vec<3> &r0, const vec<3> &r1, 
                     const vec<3> &r2, const vec<3> &r3 )
{
  x[0] = r0[0]; x[1]  = r0[1]; x[2]  = r0[2];
  x[3] = r1[0]; x[4]  = r1[1]; x[5]  = r1[2];
  x[6] = r2[0]; x[7]  = r2[1]; x[8]  = r2[2];
  x[9] = r3[0]; x[10] = r3[1]; x[11] = r3[2];
} 

template<>
matrix<4,4>::matrix( const vec<4> &r0, const vec<4> &r1, 
                     const vec<4> &r2, const vec<4> &r3 )
{
  x[0] = r0[0]; x[1] = r0[1]; x[2] = r0[2]; x[3] = r0[3];
  x[4] = r1[0]; x[5] = r1[1]; x[6] = r1[2]; x[7] = r1[3];
  x[8] = r2[0]; x[9] = r2[1]; x[10] = r2[2]; x[11] = r2[3];
  x[12] = r3[0]; x[13] = r3[1]; x[14] = r3[2]; x[15] = r3[3];
} 

template<>
void matrix<2,2>::set( const vec<2> &r0, const vec<2> &r1 )
{
  x[0] = r0[0]; x[1] = r0[1];
  x[2] = r1[0]; x[3] = r1[1];
}

template<>
void matrix<3,3>::set(const vec<3> &r0, const vec<3> &r1, const vec<3> &r2)
{
  x[0] = r0[0]; x[1] = r0[1]; x[2] = r0[2];
  x[3] = r1[0]; x[4] = r1[1]; x[5] = r1[2];
  x[6] = r2[0]; x[7] = r2[1]; x[8] = r2[2];
}



template<>
void matrix<4,2>::set( const vec<2> &r0, const vec<2> &r1, 
                       const vec<2> &r2, const vec<2> &r3 )
{ 
  x[0] = r0[0]; x[1] = r0[1];
  x[2] = r1[0]; x[3] = r1[1];
  x[4] = r2[0]; x[5] = r2[1];
  x[6] = r3[0]; x[7] = r3[1];
} 

template<>
void matrix<4,3>::set( const vec<3> &r0, const vec<3> &r1,
			      const vec<3> &r2, const vec<3> &r3 )
{
  x[0] = r0[0]; x[1]  = r0[1]; x[2]  = r0[2];
  x[3] = r1[0]; x[4]  = r1[1]; x[5]  = r1[2];
  x[6] = r2[0]; x[7]  = r2[1]; x[8]  = r2[2];
  x[9] = r3[0]; x[10] = r3[1]; x[11] = r3[2];
} 

template<>
void matrix<4,4>::set( const vec<4> &r0, const vec<4> &r1, 
                       const vec<4> &r2, const vec<4> &r3 )
{
  x[0] = r0[0]; x[1] = r0[1]; x[2] = r0[2]; x[3] = r0[3];
  x[4] = r1[0]; x[5] = r1[1]; x[6] = r1[2]; x[7] = r1[3];
  x[8] = r2[0]; x[9] = r2[1]; x[10] = r2[2]; x[11] = r2[3];
  x[12] = r3[0]; x[13] = r3[1]; x[14] = r3[2]; x[15] = r3[3];
}

//-------------------------------------------------------------------------
// function   : basic arthmetic
// description: 
//-------------------------------------------------------------------------

template<>
matrix<2,2>& 
matrix<2,2>::operator*=(const matrix<2,2> &rhs)
{
  matrix<2,2> lhs(*this);
  
  x[0] = lhs(0,0)*rhs(0,0) + lhs(0,1)*rhs(1,0);
  x[1] = lhs(0,0)*rhs(0,1) + lhs(0,1)*rhs(1,1);
  
  x[2] = lhs(1,0)*rhs(0,0) + lhs(1,1)*rhs(1,0);
  x[3] = lhs(1,0)*rhs(0,1) + lhs(1,1)*rhs(1,1);

  return *this;
}


template<>
matrix<2,2> operator*( const matrix<2,2> &a, const matrix<2,2> &b )
{
  matrix<2,2> ret(a);
  ret *= b;
  return ret;
}


template<>
matrix<3,3>& 
matrix<3,3>::operator*=(const matrix<3,3> &rhs)
{
  matrix<3,3> lhs(*this);
  
  x[0] = lhs(0,0)*rhs(0,0) + lhs(0,1)*rhs(1,0) + lhs(0,2)*rhs(2,0);
  x[1] = lhs(0,0)*rhs(0,1) + lhs(0,1)*rhs(1,1) + lhs(0,2)*rhs(2,1);
  x[2] = lhs(0,0)*rhs(0,2) + lhs(0,1)*rhs(1,2) + lhs(0,2)*rhs(2,2);
  
  x[3] = lhs(1,0)*rhs(0,0) + lhs(1,1)*rhs(1,0) + lhs(1,2)*rhs(2,0);
  x[4] = lhs(1,0)*rhs(0,1) + lhs(1,1)*rhs(1,1) + lhs(1,2)*rhs(2,1);
  x[5] = lhs(1,0)*rhs(0,2) + lhs(1,1)*rhs(1,2) + lhs(1,2)*rhs(2,2);
  
  x[6] = lhs(2,0)*rhs(0,0) + lhs(2,1)*rhs(1,0) + lhs(2,2)*rhs(2,0);
  x[7] = lhs(2,0)*rhs(0,1) + lhs(2,1)*rhs(1,1) + lhs(2,2)*rhs(2,1);
  x[8] = lhs(2,0)*rhs(0,2) + lhs(2,1)*rhs(1,2) + lhs(2,2)*rhs(2,2);

  return *this;
}


template<>
matrix<3,3> 
operator*( const matrix<3,3> &a, const matrix<3,3> &b )
{
  matrix<3,3> ret(a);
  ret *= b;
  return ret;
}

template <>
vec<4> operator * ( const matrix<4,4> &m, const vec<4> &v )
{
	float v0 = v(0);
	float v1 = v(1);
	float v2 = v(2);
	float v3 = v(3);
	
  vec<4> ret(
  m(0,0)*v0 + m(0,1)*v1 + m(0,2)*v2 + m(0,3)*v3,
  m(1,0)*v0 + m(1,1)*v1 + m(1,2)*v2 + m(1,3)*v3,
  m(2,0)*v0 + m(2,1)*v1 + m(2,2)*v2 + m(2,3)*v3,
  m(3,0)*v0 + m(3,1)*v1 + m(3,2)*v2 + m(3,3)*v3
  );
  return ret;
}

template <>
vec<4> operator * ( const vec<4> &v, const matrix<4,4> &m )
{
	float v0 = v(0);
	float v1 = v(1);
	float v2 = v(2);
	float v3 = v(3);
	
  vec<4> ret(
  m(0,0)*v0 + m(0,1)*v1 + m(0,2)*v2 + m(0,3)*v3,
  m(1,0)*v0 + m(1,1)*v1 + m(1,2)*v2 + m(1,3)*v3,
  m(2,0)*v0 + m(2,1)*v1 + m(2,2)*v2 + m(2,3)*v3,
  m(3,0)*v0 + m(3,1)*v1 + m(3,2)*v2 + m(3,3)*v3
  );
  return ret;
}

template <>
vec<3> operator * ( const vec<4> &v, const matrix<4,3> &m )
{
	float v0 = v(0);
	float v1 = v(1);
	float v2 = v(2);
	float v3 = v(3);
	
  vec<3> ret(
  m(0,0)*v0 + m(0,1)*v1 + m(0,2)*v2 + m(0,3)*v3,
  m(1,0)*v0 + m(1,1)*v1 + m(1,2)*v2 + m(1,3)*v3,
  m(2,0)*v0 + m(2,1)*v1 + m(2,2)*v2 + m(2,3)*v3
  );
  return ret;
}


//-------------------------------------------------------------------------
// function   : matrix methods
// description: 
//-------------------------------------------------------------------------

template<>
matrix<2,2>& matrix<2,2>::invert()
{
  matrix<2,2> m(*this);

  float det = m(0,0)*m(1,1) - m(0,1)*m(1,0);
  
  if ( !det )
  {
    std::cout << "matrix<2,2>::invert() : determinant equal to zero!\n";
    return identity();
  }
  
  (*this)(1,1) =  m(0,0);
  (*this)(0,1) = -m(0,1);
  (*this)(1,0) = -m(1,0);
  (*this)(0,0) =  m(1,1);

  (*this) /= det;
  return *this;
}

template<>
matrix<3,3>& matrix<3,3>::invert()
{
  matrix<3,3> m(*this);

  // get the inverse determinate
  float idet = m(0,0)*m(1,1)*m(2,2) - m(0,0)*m(1,2)*m(2,1) -
	             m(0,1)*m(1,0)*m(2,2) + m(0,1)*m(1,2)*m(2,0) +
	             m(0,2)*m(1,0)*m(2,1) - m(0,2)*m(1,1)*m(2,0);

  if ( !idet )
  {
    std::cout << "matrix<3,3>::invert() : determinant equal to zero!\n";
	  return identity();
  }

  idet = 1.0f / idet;

  // build the inverse matrix
  (*this)(0,0) = idet * (m(1,1)*m(2,2) - m(1,2)*m(2,1));
  (*this)(0,1) = idet * (m(0,2)*m(2,1) - m(0,1)*m(2,2));
  (*this)(0,2) = idet * (m(0,1)*m(1,2) - m(0,2)*m(1,1));
  (*this)(1,0) = idet * (m(1,2)*m(2,0) - m(1,0)*m(2,2));
  (*this)(1,1) = idet * (m(0,0)*m(2,2) - m(0,2)*m(2,0));
  (*this)(1,2) = idet * (m(0,2)*m(1,0) - m(0,0)*m(1,2));
  (*this)(2,0) = idet * (m(1,0)*m(2,1) - m(1,1)*m(2,0));
  (*this)(2,1) = idet * (m(0,1)*m(2,0) - m(0,0)*m(2,1));
  (*this)(2,2) = idet * (m(0,0)*m(1,1) - m(0,1)*m(1,0));
        
  return *this;
}

template<>
matrix<4,4>& matrix<4,4>::invert() 
{
  matrix<4,4> a(*this);
  matrix<4,4> b; b.identity();

  unsigned int R, C;
  unsigned int cc;
  unsigned int rowMax; // Points to max abs value row in this column
  unsigned int row;
  float tmp;

  // Go through columns
  for (C=0; C<4; C++)
  {

  // Find the row with max value in this column
  rowMax = C;
  for (R=C+1; R<4; R++)
    if (fabs(a(R,C)) > fabs(a(rowMax,C)))
      rowMax = R;

  // If the max value here is 0, we can't invert.  Return identity.
  if (a(rowMax,C) == 0.0F)
    return(identity());

  // Swap row "rowMax" with row "c"
  for (cc=0; cc<4; cc++)
  {
    tmp = a(C,cc);
    a(C,cc) = a(rowMax,cc);
    a(rowMax,cc) = tmp;
    tmp = b(C,cc);
    b(C,cc) = b(rowMax,cc);
    b(rowMax,cc) = tmp;
  }

  // Now everything we do is on row "c".
  // Set the max cell to 1 by dividing the entire row by that value
  tmp = a(C,C);
  for (cc=0; cc<4; cc++)
  {
    a(C,cc) /= tmp;
    b(C,cc) /= tmp;
  }

  // Now do the other rows, so that this column only has a 1 and 0's
  for (row = 0; row < 4; row++)
  {
    if (row != C)
    {
    tmp = a(row,C);
    for (cc=0; cc<4; cc++)
    {
      a(row,cc) -= a(C,cc) * tmp;
      b(row,cc) -= b(C,cc) * tmp;
    }
    }
  }

  }

  *this = b;

  return *this;
}

/************************************************************************/

quaternion RotationArc(const vec<3> &v0, const vec<3> &v1)
{
  quaternion q;
  // v0.normalize(); 
  // v1.normalize();  // If vector is already unit length then why do it again?
  vec<3> c = CrossProduct(v0,v1);
  float  d = DotProduct(v0,v1);
  float  s = (float)sqrt((1+d)*2);
  q.x = c(0) / s;
  q.y = c(1) / s;
  q.z = c(2) / s;
  q.w = s /2.0f;
  return q;
}

void MakeMatrix44(float *m, const vec<3> &v, const quaternion &q) 
{
	// a lean function for filling m[16] with
	// a 4x4 transformation matrix based on 
	// translation v and rotation q
	// This routine would likely be used by an opengl
	// programmer calling glmultmatrix()
	m[0] = 1-2*(q.y*q.y+q.z*q.z);
	m[1] = 2*(q.x*q.y+q.w*q.z)  ;
	m[2] = 2*(q.x*q.z-q.w*q.y)  ;
	m[3] = 0            ;
	m[4] = 2*(q.x*q.y-q.w*q.z)  ;	
	m[5] = 1-2*(q.x*q.x+q.z*q.z);
	m[6] = 2*(q.y*q.z+q.w*q.x)  ;
	m[7] = 0    ;   
	m[8] = 2*(q.x*q.z+q.w*q.y)  ; 
	m[9] = 2*(q.y*q.z-q.w*q.x)  ; 
	m[10]= 1-2*(q.x*q.x+q.y*q.y);
	m[11]= 0    ; 
	m[12] = v(0) ;
	m[13] = v(1) ;
	m[14] = v(2) ;
	m[15] = 1.0f;
}

quaternion VirtualTrackBall( const vec<3> &cop,
                             const vec<3> &cor,
                             const vec<3> &dir1,
                             const vec<3> &dir2) 
{
	// Implement track ball functionality to spin stuf on the screen
	//  cop   center of projection
	//  cor   center of rotation
	//  dir1  old mouse direction 
	//  dir2  new mouse direction
	// pretend there is a sphere around cor.  Then find the points
	// where dir1 and dir2 intersect that sphere.  Find the
	// rotation that takes the first point to the second.
	// If the line "dir" doesn't hit the sphere then take the closest
	// point on the sphere to that line.
	float m;
	// compute plane 
	vec<3> nrml = cor - cop;
	float fudgefactor = 1.0f/(nrml.length() * 0.25f); // since trackball proportional to distance from cop
	nrml.normalize();
	float dist = -DotProduct(nrml,cor);
	vec<3> u= PlaneLineIntersection(nrml,dist,cop,cop+dir1);
	u=u-cor;
	u=u*fudgefactor;
	m= u.length();
	if(m>1) 
    {
        u=u*1.0f/m;
    }
	else 
    {
		u=u - (nrml * (float)sqrt(1-m*m));
	}
	vec<3> v= PlaneLineIntersection(nrml,dist,cop,cop+dir2);
	v=v-cor;
	v=v*fudgefactor;
	m= v.length();
	if(m>1) 
    {
        v=v*1.0f/m;
    }
	else 
    {
		v=v - (nrml * (float)sqrt(1-m*m));
	}
	return RotationArc(u,v);
}

vec<3> PlaneLineIntersection( const vec<3> &n, float d,
                              const vec<3> &p1, const vec<3> &p2)
{
  // returns the point where the line p1-p2 intersects the plane n&d
  vec<3> dif  = p2-p1;
  float dn= DotProduct(n,dif);
  float t = -(d+DotProduct(n,p1) )/dn;
  return p1 + (dif*t);
}

//------------------------------------------------------------------------
// Function    : operator << 
// Description : output function for the quaternion class
//------------------------------------------------------------------------
std::ostream& operator << ( std::ostream& outs, const quaternion& source ) 
{
  outs << "(" << source.x << ", " << source.y << ", " << source.z <<
    ", " << source.w << ")" << std::endl;
  
  return outs; 
}




