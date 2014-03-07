#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <IOScalarField.h>

//------------------------------------------------------------------------
// Function    : constructor and destructor
// Description : 
//------------------------------------------------------------------------
IOScalarField::IOScalarField( const char *filename, int kernel_type,
                              int main_indicator ) :
  Surface()
{
  std::ifstream in;
  in.open(filename);

  if (! in)
  {
    // Bail entirely out if file processing fails.
    std::cout << "Error opening indicator file: " << filename << std::endl;
    exit ( 1 );
  }

  std::string header;
  if (getline(in, header, '\n'))
  {
    std::istringstream iss(header);
    iss.exceptions(std::ifstream::failbit | std::ifstream::badbit);
    try
    {
      iss >> _num_indicators;
    }
    catch (...)
    {
      std::cout << "File formatting error: first line of file "
                << filename
                << " must contain number of materials"
                << std::endl;
      exit( 1 );
    }
  }

  if ( _num_indicators == 0 )
  {
    std::cout << "PROBLEM: Zero indicator functions specified!" <<
      std::endl;
    exit( 1 );
  }

  // check if we need to inquire about the isovalue
  if ( _num_indicators == 1 )
  {
    std::cout << "Enter isovalue: ";
    std::cin >> _isovalue;
  }
  else {
    _isovalue = 0.0;
  }

  // construct ReconstructionFields for each indicator file
  _indicators = new Surface*[_num_indicators];
  float scale_x = 0.0f, scale_y = 0.0f, scale_z = 0.0f;
  bool newIndicatorFile = false;

  for ( int i = 0; i < _num_indicators; ++i )
  {
    std::string line;
    std::string buffer;
    if (getline(in, line, '\n'))
    {
      if (newIndicatorFile)
      {
        buffer = line;
      }
      else
      {
        std::istringstream lss(line);
        lss.exceptions(std::ifstream::failbit | std::ifstream::badbit);

        try
        {
          // old-style indicator file line:
          // scale_x scale_y scale_x filename
          lss >> scale_x >> scale_y >> scale_z >> buffer;
        }
        catch(...) {
          newIndicatorFile = true;
        }
        
        if (newIndicatorFile)
        {
          buffer = lss.str();
        }
      }
    }
    else {
      // Bail entirely out if file processing fails.
      std::cout << "Missing line in " << filename
                << ". Exiting" << endl;
          exit ( 1 );
    }

    std::cout << "ScalarField read file: " << buffer << std::endl;

    switch ( kernel_type )
    {
      case INTERPOLATING:
      _indicators[i] =
        new ScalarField( buffer.c_str(), ScalarField::CATMULLROM,
                         scale_x, scale_y, scale_z );
      break;

      case APPROXIMATING:
      _indicators[i] =
        new ScalarField( buffer.c_str(), ScalarField::BSPLINE,
                         scale_x, scale_y, scale_z );
      break;

      default:
      std::cout << "No such kernel type" << std::endl;
      exit( 1 );
      break;
    }
  }

  in.close();

  // initialize the main indicator function to be the first one
  _main_indicator = main_indicator;

  _indicators[_main_indicator]->domain( _start, _end );

//   // check if the user wants to specify the epsilon value
//   char yesorno;
//   do
//   {
//     cout << "Do you want to enter an epsilon value (default = 1.0e-4)" <<
//       " (y|n)?  ";
//     cin >> yesorno;
//   } while ( (yesorno != 'y') && (yesorno != 'n') );

  
  char yesorno='n';
  if ( yesorno == 'y' )
  {
    float epsilon;
    
    std::cout << "   epsilon : ";
    std::cin >> epsilon;
    _epsilon_sqrd = epsilon*epsilon;
  }
  else
    _epsilon_sqrd = static_cast<float>(EPSILON_SQRD);  
}

IOScalarField::~IOScalarField()
{
  for ( int i = 0; i < _num_indicators; i++ )
    delete _indicators[i];
  delete [] _indicators;
}

//------------------------------------------------------------------------
// Function    : interp_max()
// Description : 
//------------------------------------------------------------------------
float IOScalarField::interp_max( const float* const f, int num_values ) const
{
  float max_value = 0.0;
  float val;
  for ( int i = 0; i < num_values; i++ )
  {
    val = f[i];
    for ( int j = 1; j < num_values; j++ )
      val *= 1.0f + sign_function(f[i] - f[(i+j)%num_values]);

    max_value += val;
  }

  max_value /= power_function( 2.0, (num_values-1) );
  return max_value;
}

void IOScalarField::interp_max_D( const float* const f, 
                                  const vec<3>* const fx, 
                                  int num_values, vec<3> &D ) const
{
  D = 0.0;
  vec<3> tmpv1, tmpv2, tmpv3;
  for ( int i = 0; i < num_values; i++ )
  {
    tmpv1 = fx[i];
    for ( int j = 1; j < num_values; j++ )
      tmpv1 *= 1.0f + sign_function(f[i] - f[(i+j)%num_values]);

    tmpv2 = 0.0;
    for ( int j = 1; j < num_values; j++ )
    {
      sign_function_D( (f[i] - f[(i+j)%num_values]),
                       (fx[i] - fx[(i+j)%num_values]), tmpv3 );
      for ( int k = 1; k < num_values; k++ )
      {
        if ( k == j ) continue;
        tmpv3 *= 1.0f + sign_function(f[i] - f[(i+k)%num_values]); 
      }
      tmpv2 += tmpv3;
    }
    tmpv2 *= f[i];

    D += tmpv1 + tmpv2;
  }

  D /= power_function( 2.0, (num_values-1) );
}

void IOScalarField::interp_max_DD( const float* const f, 
                                   const vec<3>* const fx, 
                                   const matrix<3,3>* const fxx,
                                   int num_values, 
                                   matrix<3,3> &DD ) const
{
  DD = 0.0;
  matrix<3,3> tmpm1, tmpm2, tmpm3, tmpm4;
  vec<3> tmpv1, tmpv2, tmpv3;
  for ( int i = 0; i < num_values; i++ )
  {
    tmpm1 = fxx[i];
    for ( int j = 1; j < num_values; j++ )
      tmpm1 *= 1.0f + sign_function(f[i] - f[(i+j)%num_values]);

    tmpv1 = 0.0;
    for ( int j = 1; j < num_values; j++ )
    {
      sign_function_D( (f[i] - f[(i+j)%num_values]),
                       (fx[i] - fx[(i+j)%num_values]), tmpv2 );
      for ( int k = 1; k < num_values; k++ )
      {
        if ( k == j ) continue;
        tmpv2 *= 1.0f + sign_function(f[i] - f[(i+k)%num_values]); 
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
        tmpm4 *= 1.0f + sign_function(f[i] - f[(i+k)%num_values]); 
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
          tmpv2 *= 1.0f + sign_function(f[i] - f[(i+l)%num_values]); 
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

void IOScalarField::interp_max_D_DD( const float* const f, 
                                     const vec<3>* const fx, 
                                     const matrix<3,3>* const fxx,
                                     int num_values, 
                                     vec<3> &D, matrix<3,3> &DD ) const
{
  D = 0.0;
  DD = 0.0;
  vec<3> tmpv0, tmpv1, tmpv2, tmpv3, tmpv4, tmpv5, tmpv6;
  matrix<3,3> tmpm1, tmpm2, tmpm3, tmpm4, tmpm5;
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
      tmpm1 *= 1.0f + tmpf1;
      tmpv0 *= 1.0f + tmpf1;

      sign_function_D_DD( (f[i] - f[(i+j)%num_values]),
                          (fx[i] - fx[(i+j)%num_values]), 
                          (fxx[i] - fxx[(i+j)%num_values]), 
                          tmpv2, tmpm5 );

      tmpv5 = 0.0;
      for ( int k = 1; k < num_values; k++ )
      {
        if ( k == j ) continue;

        tmpf2 = sign_function(f[i] - f[(i+k)%num_values]);
        tmpv2 *= 1.0f + tmpf2; 
        tmpm5 *= 1.0f + tmpf2; 

        sign_function_D( (f[i] - f[(i+k)%num_values]),
                         (fx[i] - fx[(i+k)%num_values]), tmpv4 );
        for ( int l = 1; l < num_values; l++ )
        {
          if ( (l == k) || (l == j) ) continue;
          tmpv4 *= 1.0f + sign_function(f[i] - f[(i+l)%num_values]); 
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



float IOScalarField::sign_function( float a ) const
{
  return ( a / sqrt( a*a + _epsilon_sqrd ));
}

void IOScalarField::sign_function_D( float a, const vec<3> &aD,
                                     vec<3> &D ) const
{
  float isqrt_a_sqrd = 1.0f/sqrt( a*a + _epsilon_sqrd );
  D = aD * 
    ( isqrt_a_sqrd - a*a*(isqrt_a_sqrd*isqrt_a_sqrd*isqrt_a_sqrd) );
}

void IOScalarField::sign_function_DD( float a, const vec<3> &aD,
                                      const matrix<3,3> &aDD, 
                                      matrix<3,3> &DD ) const
{
  float isqrt_a_sqrd = 1.0f/sqrt( a*a + _epsilon_sqrd );
  float i3 = isqrt_a_sqrd*isqrt_a_sqrd*isqrt_a_sqrd;
  float i5 = i3*isqrt_a_sqrd*isqrt_a_sqrd;

  DD = aDD*(isqrt_a_sqrd - a*a*i3) - 
    3.0*DirectProduct(aD,aD)*(a*i3 - a*a*a*i5);
}

void IOScalarField::sign_function_D_DD( float a, const vec<3> &aD,
                                        const matrix<3,3> &aDD, 
                                        vec<3> &D,
                                        matrix<3,3> &DD ) const
{
  float a_sqrd = a*a;
  float isqrt_a_sqrd = 1.0f/sqrt( a_sqrd + _epsilon_sqrd );
  float i3 = isqrt_a_sqrd*isqrt_a_sqrd*isqrt_a_sqrd;
  float i5 = i3*isqrt_a_sqrd*isqrt_a_sqrd;

  D = aD * (isqrt_a_sqrd - a_sqrd*i3);

  DD = aDD*(isqrt_a_sqrd - a_sqrd*i3) - 
    3.0*DirectProduct(aD,aD)*(a*i3 - a_sqrd*a*i5);
}

float IOScalarField::power_function( float a, int b ) const
{
  float val = 1.0;
  for ( int i = 0; i < b; i++ )
    val *= a;
  return val;
}

//------------------------------------------------------------------------
// Function    : computeSurfacePointParams()
// Description : 
//------------------------------------------------------------------------
bool IOScalarField::computeSurfacePointParams(
  const vec<3> &pos, SurfacePointParams &params,
  bool computeGradient, bool computeHessian ) 
{
  // store the values for the main indicator function
  _indicators[_main_indicator]->computeSurfacePointParams( pos, params,
                                                           computeGradient,
                                                           computeHessian );
  float f_1 = params._F;
  vec<3> fx_1 = params._Fx;
  matrix<3,3> fxx_1 = params._Fxx;
  
  float *f = new float[_num_indicators-1];
  vec<3> *fx = new vec<3>[_num_indicators-1];
  matrix<3,3> *fxx = new matrix<3,3>[_num_indicators-1];
  int counter=0;
  for ( int i = 0; i < _num_indicators; i++ )
  {
    if ( i == _main_indicator )
      continue;

    _indicators[i]->computeSurfacePointParams( pos, params, 
                                               computeGradient, 
                                               computeHessian );
    f[counter] = params._F;
    fx[counter] = params._Fx;
    fxx[counter] = params._Fxx;

    ++counter;
   }
  
  // have the combined F, the smallest sf, and Fx from the first
  //   surface

  params._F = f_1 - interp_max( f, (_num_indicators-1) ) - _isovalue;

  if ( computeGradient && !computeHessian )
  {
    interp_max_D( f, fx, (_num_indicators-1), params._Fx );
    params._Fx = fx_1 - params._Fx;
  }

  if ( computeHessian )
  {
    interp_max_D_DD( f, fx, fxx, (_num_indicators-1), 
                     params._Fx, params._Fxx );
    params._Fx = fx_1 - params._Fx;
    params._Fxx = fxx_1 - params._Fxx;
  }

  delete [] f;
  delete [] fx;
  delete [] fxx;

  return true;
}