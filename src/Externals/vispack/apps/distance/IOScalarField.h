//----------------------------------------------------------------------//
// FILE        : IOScalarField.h                                         
// DESCRIPTION : This function takes a set of indicator functions
//               (ReconstructionField) and computes the
//               inside/outside scalar field
//
//               Assumes that the indicator functions are [0.5,1.0] for
//               material and [0.0, 0.5) for not material. 
//----------------------------------------------------------------------//

#ifndef __IO_SCALAR_FIELD_H__
#define __IO_SCALAR_FIELD_H__

#include <cstdlib>
#include <iostream>
#include "ScalarField.h"
#include "ReconstructionField.h"
#include "mtxlib.h"

#define EPSILON_SQRD_DEFAULT 1.0e-12

/************************************************************************/
// NOTE: this is written for 3D only right now!
/************************************************************************/

template<int dim> class IOScalarField : public ScalarField<dim>
{
public:
  IOScalarField(const char *filename,
                int kernel_type=APPROXIMATING);
  ~IOScalarField();

  void setEpsilon(double epsilon)
  {
    _epsilon_sqrd = epsilon*epsilon;
  }

  enum {INTERPOLATING, APPROXIMATING};

  inline void setMainIndicatorFunction(int i)
  { if ( i >=  _num_indicators )
    { std::cout << "PROBLEM: main indicator number too big!" << std::endl;
      exit( 1 ); }
    _main_indicator = i; };

  void computeScalarFieldParams(const vec<dim> &pos, 
                                ScalarFieldParams<dim>&params,
                                bool computeHessian=false) const;

  virtual bool checkBounds(const vec<dim> &pos) const
  { return(_indicators[_main_indicator]->checkBounds(pos));}
  virtual vec<dim> boundingBoxLow() const
  { return(_indicators[_main_indicator]->boundingBoxLow()); };
  virtual vec<dim> boundingBoxHigh() const
  { return(_indicators[_main_indicator]->boundingBoxHigh()); };

  int xdim () const {return(_indicators[_main_indicator]->xdim());}
  int ydim () const {return(_indicators[_main_indicator]->ydim());}
  int zdim () const {return(_indicators[_main_indicator]->zdim());}

protected:
  ScalarField<dim> **_indicators;
  int _num_indicators, _main_indicator;

  double _epsilon_sqrd;

  float interp_max(const float* const f, int num_values) const;
  void  interp_max_D(const float* const f, const vec<dim>* const fx, 
                     int num_values,  vec<dim> &D) const;
  void interp_max_DD(const float* const f, const vec<dim>* const fx,
                     const matrix<dim,dim>* const fxx, int num_values,
                     matrix<dim,dim> &DD) const;
  void interp_max_D_DD(const float* const f, const vec<dim>* const fx,
                       const matrix<dim,dim>* const fxx, int num_values,
                       vec<dim> &D, matrix<dim,dim> &DD) const;

  float sign_function(float a) const;
  void sign_function_D(float a, const vec<dim> &aD, vec<dim> &D) const;
  void sign_function_DD(float a, const vec<dim> &aD,
                        const matrix<dim,dim> &aDD, matrix<dim,dim> &DD) const;
  void sign_function_D_DD(float a, const vec<dim> &aD,
                          const matrix<dim,dim> &aDD, 
                          vec<dim> &D, matrix<dim,dim> &DD) const;

  float power_function(float a, int b) const;

};

//------------------------------------------------------------------------
// Function    : constructor and destructor
// Description : 
//------------------------------------------------------------------------
template <int dim>
IOScalarField<dim>::IOScalarField( const char *filename,
                                   int kernel_type ) :
  ScalarField<dim>()
{
  _epsilon_sqrd = EPSILON_SQRD_DEFAULT;
  // open the file containing the indicator function file names
  FILE* in = fopen( filename, "r" );
  if ( !in )
  { 
    std::cout << "Error opening indicator filename file " << std::endl;
    exit ( 1 );
  }

  // read in the number of indicator functions
  fscanf( in, "%i\n", &_num_indicators );
  if ( _num_indicators == 0 )
  {
    std::cout << "PROBLEM: Zero indicator functions specified!" <<
      std::endl;
    exit( 1 );
  }

  // construct ReconstructionFields for each indiator file
  _indicators = new ScalarField<dim>*[_num_indicators];
  char buffer[300];
  float scale_x, scale_y, scale_z;
  for ( int i = 0; i < _num_indicators; i++ )
  {
    fscanf( in, "%f %f %f %s\n", &scale_x, &scale_y, &scale_z, buffer );
    std::cout << buffer << std::endl;

    switch ( kernel_type )
    {
      case INTERPOLATING:
      _indicators[i] =
        new ReconstructionField<dim>( buffer,
                                      ReconstructionField<dim>::INTERPOLATING,
                                      scale_x, scale_y, scale_z );
      break;

      case APPROXIMATING:
      _indicators[i] =
        new ReconstructionField<dim>( buffer,
                                      ReconstructionField<dim>::APPROXIMATING,
                                      scale_x, scale_y, scale_z );
      break;

      default:
      std::cout << "No such kernel type" << std::endl;
      exit( 1 );
      break;
    }
  }

  // initialize the main indicator function to be the first one
  _main_indicator = 0;
}

template <int dim>
IOScalarField<dim>::~IOScalarField()
{
  for ( int i = 0; i < _num_indicators; i++ )
    delete _indicators[i];
  delete [] _indicators;
}

//------------------------------------------------------------------------
// Function    : interp_max()
// Description : 
//------------------------------------------------------------------------
template <int dim>
float IOScalarField<dim>::interp_max( const float* const f,
                                      int num_values ) const
{
  float max_value = 0.0;
  float val;
  for ( int i = 0; i < num_values; i++ )
  {
    val = f[i];
    for ( int j = 1; j < num_values; j++ )
      val *= 1.0 + sign_function(f[i] - f[(i+j)%num_values]);

    max_value += val;
  }

  max_value /= power_function( 2.0, (num_values-1) );
  return max_value;
}

template <int dim>
void IOScalarField<dim>::interp_max_D( const float* const f, 
                                       const vec<dim>* const fx, 
                                       int num_values, vec<dim> &D ) const
{
  D = 0.0;
  vec<dim> tmpv1, tmpv2, tmpv3;
  for ( int i = 0; i < num_values; i++ )
  {
    tmpv1 = fx[i];
    for ( int j = 1; j < num_values; j++ )
      tmpv1 *= 1.0 + sign_function(f[i] - f[(i+j)%num_values]);

    tmpv2 = 0.0;
    for ( int j = 1; j < num_values; j++ )
    {
      sign_function_D( (f[i] - f[(i+j)%num_values]),
                       (fx[i] - fx[(i+j)%num_values]), tmpv3 );
      for ( int k = 1; k < num_values; k++ )
      {
        if ( k == j ) continue;
        tmpv3 *= 1.0 + sign_function(f[i] - f[(i+k)%num_values]); 
      }
      tmpv2 += tmpv3;
    }
    tmpv2 *= f[i];

    D += tmpv1 + tmpv2;
  }

  D /= power_function( 2.0, (num_values-1) );
}

template <int dim>
void IOScalarField<dim>::interp_max_DD( const float* const f, 
                                        const vec<dim>* const fx, 
                                        const matrix<dim,dim>* const fxx,
                                        int num_values, 
                                        matrix<dim,dim> &DD ) const
{
  DD = 0.0;
  matrix<dim,dim> tmpm1, tmpm2, tmpm3, tmpm4;
  vec<dim> tmpv1, tmpv2, tmpv3;
  for ( int i = 0; i < num_values; i++ )
  {
    tmpm1 = fxx[i];
    for ( int j = 1; j < num_values; j++ )
      tmpm1 *= 1.0 + sign_function(f[i] - f[(i+j)%num_values]);

    tmpv1 = 0.0;
    for ( int j = 1; j < num_values; j++ )
    {
      sign_function_D( (f[i] - f[(i+j)%num_values]),
                       (fx[i] - fx[(i+j)%num_values]), tmpv2 );
      for ( int k = 1; k < num_values; k++ )
      {
        if ( k == j ) continue;
        tmpv2 *= 1.0 + sign_function(f[i] - f[(i+k)%num_values]); 
      }
      tmpv1 += tmpv2;
    }
    tmpm2 = 2.0*DirectProduct( fx[i], tmpv1 );

    tmpm3 = 0.0;
    for ( int j = 1; j < num_values; j++ )
    {
      sign_function_DD( (f[i] - f[(i+j)%num_values]),
                        (fx[i] - fx[(i+j)%num_values]), 
                        (fxx[i] - fxx[(i+j)%num_values]), tmpm4 );
      for ( int k = 1; k < num_values; k++ )
      {
        if ( k == j ) continue;
        tmpm4 *= 1.0 + sign_function(f[i] - f[(i+k)%num_values]); 
      }
      tmpm3 += tmpm4;
    }
    tmpm3 *= f[i];

    tmpm4 = 0.0;
    for ( int j = 1; j < num_values; j++ )
    { 
      sign_function_D( (f[i] - f[(i+j)%num_values]),
                       (fx[i] - fx[(i+j)%num_values]), tmpv1 );
      tmpv3 = 0.0;
      for ( int k = 1; k < num_values; k++ )
      {
        if ( k == j ) continue;
        sign_function_D( (f[i] - f[(i+k)%num_values]),
                         (fx[i] - fx[(i+k)%num_values]), tmpv2 );
        for ( int l = 1; l < num_values; l++ )
        {
          if ( (l == k) || (l == j) ) continue;
          tmpv2 *= 1.0 + sign_function(f[i] - f[(i+l)%num_values]); 
        }
        tmpv3 += tmpv2;
      }
      tmpm4 += DirectProduct( tmpv1, tmpv3 );
    }
    tmpm4 *= f[i];

    DD += tmpm1 + tmpm2 + tmpm3 + tmpm4;
  }

  DD /= power_function( 2.0, (num_values-1) );
}

template <int dim>
void IOScalarField<dim>::interp_max_D_DD( const float* const f, 
                                          const vec<dim>* const fx, 
                                          const matrix<dim,dim>* const fxx,
                                          int num_values, 
                                          vec<dim> &D,
                                          matrix<dim,dim> &DD ) const
{
  D = 0.0;
  DD = 0.0;
  vec<dim> tmpv0, tmpv1, tmpv2, tmpv3, tmpv4, tmpv5, tmpv6;
  matrix<dim,dim> tmpm1, tmpm2, tmpm3, tmpm4, tmpm5;
  float tmpf1, tmpf2;
  for ( int i = 0; i < num_values; i++ )
  {
    tmpv0 = fx[i];
    tmpm1 = fxx[i];
    tmpm3 = 0.0;
    tmpm4 = 0.0;
    tmpv1 = 0.0;
    for ( int j = 1; j < num_values; j++ )
    {
      tmpf1 = sign_function(f[i] - f[(i+j)%num_values]);
      tmpm1 *= 1.0 + tmpf1;
      tmpv0 *= 1.0 + tmpf1;

      sign_function_D_DD( (f[i] - f[(i+j)%num_values]),
                          (fx[i] - fx[(i+j)%num_values]), 
                          (fxx[i] - fxx[(i+j)%num_values]), 
                          tmpv2, tmpm5 );

      tmpv5 = 0.0;
      for ( int k = 1; k < num_values; k++ )
      {
        if ( k == j ) continue;

        tmpf2 = sign_function(f[i] - f[(i+k)%num_values]);
        tmpv2 *= 1.0 + tmpf2; 
        tmpm5 *= 1.0 + tmpf2; 

        sign_function_D( (f[i] - f[(i+k)%num_values]),
                         (fx[i] - fx[(i+k)%num_values]), tmpv4 );
        for ( int l = 1; l < num_values; l++ )
        {
          if ( (l == k) || (l == j) ) continue;
          tmpv4 *= 1.0 + sign_function(f[i] - f[(i+l)%num_values]); 
        }
        tmpv5 += tmpv4;
      }
      tmpv1 += tmpv2;
      tmpm3 += tmpm5;  
      tmpm4 += DirectProduct( tmpv2, tmpv5 );
    }
    tmpm2 = 2.0*DirectProduct( fx[i], tmpv1 );
    tmpm3 *= f[i];
    tmpm4 *= f[i];

    D += tmpv0 + f[i]*tmpv1;
    DD += tmpm1 + tmpm2 + tmpm3 + tmpm4;
  }

  D /= power_function( 2.0, (num_values-1) );
  DD /= power_function( 2.0, (num_values-1) );
}



template <int dim>
float IOScalarField<dim>::sign_function( float a ) const
{
  return ( a / sqrt( a*a + _epsilon_sqrd ));
}

template <int dim>
void IOScalarField<dim>::sign_function_D( float a, const vec<dim> &aD,
                                          vec<dim> &D ) const
{
  float isqrt_a_sqrd = 1.0/sqrt( a*a + _epsilon_sqrd );
  D = aD * 
    ( isqrt_a_sqrd - a*a*(isqrt_a_sqrd*isqrt_a_sqrd*isqrt_a_sqrd) );
}

template <int dim>
void IOScalarField<dim>::sign_function_DD( float a, const vec<dim> &aD,
                                           const matrix<dim,dim> &aDD, 
                                           matrix<dim,dim> &DD ) const
{
  float isqrt_a_sqrd = 1.0/sqrt( a*a + _epsilon_sqrd );
  float i3 = isqrt_a_sqrd*isqrt_a_sqrd*isqrt_a_sqrd;
  float i5 = i3*isqrt_a_sqrd*isqrt_a_sqrd;

  DD = aDD*(isqrt_a_sqrd - a*a*i3) - 
    3.0*DirectProduct(aD,aD)*(a*i3 - a*a*a*i5);
}

template <int dim>
void IOScalarField<dim>::sign_function_D_DD( float a, const vec<dim> &aD,
                                             const matrix<dim,dim> &aDD, 
                                             vec<dim> &D,
                                             matrix<dim,dim> &DD ) const
{
  float a_sqrd = a*a;
  float isqrt_a_sqrd = 1.0/sqrt( a_sqrd + _epsilon_sqrd );
  float i3 = isqrt_a_sqrd*isqrt_a_sqrd*isqrt_a_sqrd;
  float i5 = i3*isqrt_a_sqrd*isqrt_a_sqrd;

  D = aD * (isqrt_a_sqrd - a_sqrd*i3);

  DD = (isqrt_a_sqrd - a_sqrd*i3) * aDD - 
    3.0*(a*i3 - a_sqrd*a*i5)*DirectProduct(aD,aD);
}

template <int dim>
float IOScalarField<dim>::power_function( float a, int b ) const
{
  float val = 1.0;
  for ( int i = 0; i < b; i++ )
    val *= a;
  return val;
}

#endif // __IO_SCALAR_FIELD_H__
