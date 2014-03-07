#include <cmath>
#include <iostream>

#ifdef _WIN32
#pragma warning( disable: 4244 )
#endif


//-------------------------------------------------------------------------
// function   : Constructors
// description: 
//-------------------------------------------------------------------------
template <int dim>
vec<dim>::vec( float rhs )
{
  for ( int i = 0; i < dim; i++ )
    x[i] = rhs;
}

template <int dim>
vec<dim>& vec<dim>::operator = ( const vec &rhs )
{
  for ( int i = 0; i < dim; i++ )
    x[i] = rhs[i];
    
  return *this;
}

template <int dim>
vec<dim>& vec<dim>::operator = ( float rhs )
{
  for ( int i = 0; i < dim; i++ )
    x[i] = rhs;
    
  return *this;
}

//-------------------------------------------------------------------------
// function   : basic arthmetic
// description: 
//-------------------------------------------------------------------------
template <int dim> 
vec<dim>& vec<dim>::operator += ( const vec<dim> &rhs )
{
  for ( int i = 0; i < dim; i++ )
    x[i] += rhs[i];
    
  return *this;
}

template <int dim>
vec<dim> operator + ( const vec<dim> &a, const vec<dim> &b )
{
  vec<dim> ret(a);
  ret += b;
  return ret; 
}

template <int dim> 
vec<dim>& vec<dim>::operator -= ( const vec<dim> &rhs )
{
  for ( int i = 0; i < dim; i++ )
    x[i] -= rhs[i];
    
  return *this;
}

template <int dim>
vec<dim> operator - ( const vec<dim> &a, const vec<dim> &b )
{
  vec<dim> ret(a);
  ret -= b;
  return ret; 
}

template <int dim> 
vec<dim>& vec<dim>::operator *= ( const vec<dim> &rhs )
{
  for ( int i = 0; i < dim; i++ )
    x[i] *= rhs[i];
    
  return *this;
}

template <int dim> 
vec<dim>& vec<dim>::operator *= ( float rhs )
{
  for ( int i = 0; i < dim; i++ )
    x[i] *= rhs;
    
  return *this;
}

template <int dim>
vec<dim> operator * ( const vec<dim> &v, float f )
{
  vec<dim> ret(v);
  ret *= f;
  return ret; 
}

template <int dim>
vec<dim> operator * ( float f, const vec<dim> &v )
{
  vec<dim> ret(v);
  ret *= f;
  return ret; 
}

template <int dim>
vec<dim> operator * ( const vec<dim> &a, const vec<dim> &b )
{
  vec<dim> ret(a);
  ret *= b;
  return ret; 
}

template <int dim> 
vec<dim>& vec<dim>::operator /= ( const vec<dim> &rhs )
{
  for ( int i = 0; i < dim; i++ )
    x[i] /= rhs[i];
    
  return *this;
}

template <int dim> 
vec<dim>& vec<dim>::operator /= ( float rhs )
{
  for ( int i = 0; i < dim; i++ )
    x[i] /= rhs;
    
  return *this;
}

template <int dim>
vec<dim> operator / ( const vec<dim> &v, float f )
{
  vec<dim> ret(v);
  ret /= f;
  return ret; 
}

template <int dim>
vec<dim> operator / ( const vec<dim> &a, const vec<dim> &b )
{
  vec<dim> ret(a);
  ret /= b;
  return ret; 
}

template <int dim>
vec<dim> operator - (const vec<dim> &a)
{
  vec<dim> ret(a);
  for ( int i = 0; i < dim; i++ )
    ret[i] = -ret[i];
  return ret;
}

//-------------------------------------------------------------------------
// function   : equal ops
// description: 
//-------------------------------------------------------------------------
template <int dim>
bool operator == ( const vec<dim> &a, const vec<dim> &b )
{
  for ( int i = 0; i < dim; i++ )
    if ( a[i] != b[i] )
      return false;
  return true; 
}

template <int dim>
bool operator != ( const vec<dim> &a, const vec<dim> &b )
{
  for ( int i = 0; i < dim; i++ )
    if ( a[i] != b[i] )
      return true;
  return false; 
}

//-------------------------------------------------------------------------
// function   : vector methods
// description: 
//-------------------------------------------------------------------------
template <int dim> 
float vec<dim>::length( float epsilon ) const
{
  float len = 0.0;
  for ( int i = 0; i < dim; i++ )
    len += x[i]*x[i];
    
  return sqrt( len + epsilon );
}

template <int dim> 
float vec<dim>::lengthSqr() const
{
  float len = 0.0;
  for ( int i = 0; i < dim; i++ )
    len += x[i]*x[i];
    
  return len;
}

template <int dim>
vec<dim>& vec<dim>::normalize()
{
  float len = length();
  
  if ( len > 0.0 )
    len = 1.0 / len;
  else
    len = 0.0;
    
  for ( int i = 0; i < dim; i++ )
    x[i] *= len;
    
  return *this;
}

template <int dim>
vec<dim>& vec<dim>::random()
{
  for ( int i = 0; i < dim; i++ )
    x[i] = rand();
    
  return *this;
}

//------------------------------------------------------------------------
// Function    : operator << 
// Description : output function for the vec class
//------------------------------------------------------------------------
template <int dim>
std::ostream& operator << ( std::ostream& outs, const vec<dim>& source ) 
{
  outs << "(";
  for ( int i = 0; i < dim-1; i++ )
    outs << source[i] << ", ";
  outs << source[dim-1] << ")";
  
  return outs; 
}


/************************************************************************/


//-------------------------------------------------------------------------
// function   : Constructors
// description: 
//-------------------------------------------------------------------------
template <int r, int c> 
matrix<r,c>::matrix( float f )
{
  for ( int i = 0; i < r*c; i++ )
    x[i] = f;
}

template <int r, int c>
matrix<r,c>& matrix<r,c>::operator = ( const matrix<r,c> &rhs )
{
  for ( int i = 0; i < r*c; i++ )
    x[i] = rhs.x[i];
    
  return *this;
}

template <int r, int c>
matrix<r,c>& matrix<r,c>::operator = ( float rhs )
{
  for ( int i = 0; i < r*c; i++ )
    x[i] = rhs;
    
  return *this;
}

//-------------------------------------------------------------------------
// function   : basic arthmetic
// description: 
//-------------------------------------------------------------------------
template <int r, int c>
matrix<r,c>& matrix<r,c>::operator += ( const matrix<r,c> &rhs )
{
  for ( int i = 0; i < r*c; i++ )
    x[i] += rhs.x[i];
    
  return *this;
}

template <int r, int c>
matrix<r,c> operator + ( const matrix<r,c> &a, const matrix<r,c> &b )
{
  matrix<r,c> ret(a);
  ret += b;
  return ret; 
}

template <int r, int c>
matrix<r,c>& matrix<r,c>::operator -= ( const matrix<r,c> &rhs )
{
  for ( int i = 0; i < r*c; i++ )
    x[i] -= rhs.x[i];
    
  return *this;
}

template <int r, int c>
matrix<r,c> operator - ( const matrix<r,c> &a, const matrix<r,c> &b )
{
  matrix<r,c> ret(a);
  ret -= b;
  return ret; 
}

template <int r, int c>
matrix<r,c>& matrix<r,c>::operator *= ( float rhs )
{
  for ( int i = 0; i < r*c; i++ )
    x[i] *= rhs;
    
  return *this;
}



// Adding done buy Pierre D. Georgel on 08/01/2005
/*
template <int r, int c, int c2>
matrix<r,c2> operator * (const matrix<r,c> &a, const matrix<c,c2> &b )
{
  matrix<r,c2> ret(0);  // create the matrix full of zero values
  for ( int i=0; i < r; i++ )
    for ( int j=0; j < c2; j++ )
      for ( int j2=0; j2 < c; j2++)
        ret( i, j) += a ( i, j2) * b ( j2, j);
  return ret;
}
*/


template <int r, int c>
matrix<r,r> operator * (const matrix<r,c> &a, const matrix<c,r> &b )
{
  matrix<r,r> ret(0);  // create the matrix full of zero values
  for ( int i=0; i < r; i++ )
    for ( int j=0; j < c; j++ )
      for ( int j2=0; j2 < r; j2++)
        ret(i,j) += a(i,j2) * b(j2,j);
  return ret;
}


template <int r, int c>
matrix<r,c> operator * ( float f, const matrix<r,c> &a )
{
  matrix<r,c> ret(a);
  ret *= f;
    
  return ret;
}

template <int r, int c>
matrix<r,c> operator * ( const matrix<r,c> &a, float f )
{
  matrix<r,c> ret(a);
  ret *= f;
    
  return ret;
}

template <int r, int c>
vec<r> operator * ( const matrix<r,c> &m, const vec<c> &v )
{
  vec<r> ret(0.0);
  for ( int i = 0; i < r; i++ )
    for ( int j = 0; j < c; j++ )
      ret(i) += m(i,j)*v(j);
  return ret;
}

template <int r, int c>
vec<c> operator * ( const vec<r> &v, const matrix<r,c> &m )
{
  vec<c> ret(0.0);
  for ( int i = 0; i < c; i++ )
    for ( int j = 0; j < r; j++ )
      ret(i) += v(j)*m(j,i);
  return ret;
}

template <int dim>
float DotProduct(const vec<dim> &a, const vec<dim> &b)
{
  float ret = 0.0;
  for ( int i = 0; i < dim; i++ )
    ret += a[i]*b[i];
  return ret;
}

template <int r, int c>
matrix<r,c> operator / ( const matrix<r,c> &a, float f )
{
  matrix<r,c> ret(a);
  ret /= f;
    
  return ret;
}

template <int r, int c>
matrix<r,c>& matrix<r,c>::operator /= ( float rhs )
{
  for ( int i = 0; i < r*c; i++ )
    x[i] /= rhs;
    
  return *this;
}

template <int r, int c>
matrix<r,c> operator - (const matrix<r,c> &a)
{
  matrix<r,c> ret(a);
  for ( int i = 0; i < r*c; i++ )
    ret.x[i] = -ret.x[i];
  return ret;
}

//-------------------------------------------------------------------------
// function   : matrix methods
// description: 
//-------------------------------------------------------------------------
template <int r, int c>
matrix<r,c>& matrix<r,c>::identity()
{
#ifdef _DEBUG
  if ( r != c )
  {
    std::cout << "matrix<>::identity() : ERROR --> matrix must be square!\n";
    return *this;
  }
#endif

  (*this) = 0.0;
  for ( int i = 0; i < r; i++ )
    (*this)(i,i) = 1.0;
  return *this;
}

template <int r, int c>
matrix<r,c>& matrix<r,c>::transpose()
{
#ifdef _DEBUG
  if ( r != c )
  {
    std::cout << "matrix<>::transpose() : ERROR --> matrix must be square!\n";
    return *this;
  }
#endif

  float tmp;
  for ( int i = 0; i < r; i++ )
    for ( int j = i+1; j < c; j++ )
    {
      tmp = (*this)(i,j);
      (*this)(i,j) = (*this)(j,i);
      (*this)(j,i) = tmp;
    }
  return *this;
}

template <int r, int c>
float matrix<r,c>::trace() const
{
#ifdef _DEBUG
  if ( r != c )
  {
    std::cout << "matrix<>::trace() : ERROR --> matrix must be square!\n";
    return 0.0;
  }
#endif

  float t = 0.0;
  for ( int i = 0; i < r; i++ )
    t += (*this)(i,i);
  return t;
}

template <int r, int c>
float matrix<r,c>::norm() const
{
#ifdef _DEBUG
  if ( r != c )
  {
    std::cout << "matrix<>::norm() : ERROR --> matrix must be square!\n";
    return 0.0;
  }
#endif

  float n = 0.0;
  for ( int i = 0; i < r; i++ )
    for ( int j = 0; j < r; j++ )
      n += (*this)(i,j)*(*this)(i,j);
  return sqrt( n );
}

template <int r, int c>
matrix<r,c> DirectProduct( const vec<r> &v, const vec<c> &vT )
{
  matrix<r,c> ret;
  for ( int i = 0; i < r; i++ )
    for ( int j = 0; j < c; j++ )
      ret(i,j) = v(i)*vT(j);
  return ret;
}


//-------------------------------------------------------------------------
// function   : equal ops
// description: 
//-------------------------------------------------------------------------
template <int r, int c>
bool operator == ( const matrix<r,c> &a, const matrix<r,c> &b )
{
  for ( int i = 0; i < r*c; i++ )
    if ( a.x[i] != b.x[i] )
      return false;
  return true; 
}

template <int r, int c>
bool operator != ( const matrix<r,c> &a, const matrix<r,c> &b )
{
  for ( int i = 0; i < r*c; i++ )
    if ( a.x[i] != b.x[i] )
      return true;
  return false; 
}

//------------------------------------------------------------------------
// Function    : operator << 
// Description : output function for the matrix class
//------------------------------------------------------------------------
template <int r, int c>
std::ostream& operator << ( std::ostream& outs, const matrix<r,c>& source ) 
{
  for ( int i = 0; i < r; i++ )
  {
    outs << "|";
    for ( int j = 0; j < c-1; j++ )
      outs << source.x[i*c+j] << ", ";
    outs << source.x[i*c+c-1] << "|\n";
  }
    
  return outs; 
}


