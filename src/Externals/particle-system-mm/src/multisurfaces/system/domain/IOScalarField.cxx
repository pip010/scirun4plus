#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <system/domain/ScalarField.h>
#include <system/domain/IOScalarField.h>

using namespace std;
using namespace particle_sys;

//------------------------------------------------------------------------
// Function    : constructor
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
  float scale_x = 1.0f, scale_y = 1.0f, scale_z = 1.0f;
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

  // set the domain bounds
  _indicators[_main_indicator]->bounds( _d_start, _d_end );

}

IOScalarField::~IOScalarField()
{
  for ( int i = 0; i < _num_indicators; i++ )
    delete _indicators[i];
  delete [] _indicators;
}

//------------------------------------------------------------------------
// Function    : maxSF
// Description : 
//------------------------------------------------------------------------
float IOScalarField::maxSF() const
{
  float max_sf = _indicators[0]->maxSF();
  for ( int i = 1; i < _num_indicators; i++ )
    max_sf = max( _indicators[i]->maxSF(), max_sf );
  return max_sf;
}

//------------------------------------------------------------------------
// Function    : printFunctionValues
// Description : 
//------------------------------------------------------------------------
void IOScalarField::printFunctionValues( const vector_type &pos )
{
  SurfacePointParams params;
  cout << "main indicator number : " << _main_indicator << endl;
  for ( int i = 0; i < _num_indicators; i++ )
  {
    _indicators[i]->computeSurfacePointParams( pos,params,false );
    cout << "function " << i << " : " << params._F << "  " <<
      params._Fx << endl;
  }
  computeSurfacePointParams( pos, params, false );
  cout << "IO function : " << params._F << "  " << 
    params._Fx << endl;
}

//------------------------------------------------------------------------
// Function    : computeSurfacePointParams()
// Description : 
//------------------------------------------------------------------------
bool IOScalarField::computeSurfacePointParams( const vector_type &pos, 
                                               SurfacePointParams 
                                               &params,
                                               bool computeHessian ) const
{  
  // store the values for the main indicator function
  _indicators[_main_indicator]->computeSurfacePointParams( pos, params,
                                                           computeHessian );
  float f_1 = params._F;
  float sf = params._sf;
  vector_type fx_1 = params._Fx;

  float *f = new float[_num_indicators-1];
  vector_type *fx = new vector_type[_num_indicators-1];
  int counter=0;
  for ( int i = 0; i < _num_indicators; i++ )
  {
    if ( i == _main_indicator )
      continue;

    _indicators[i]->computeSurfacePointParams( pos, params, 
                                               computeHessian );
    f[counter] = params._F;
    fx[counter] = params._Fx;
    sf = min( params._sf, sf );

    ++counter;
   }
  
  // have the combined F, the smallest sf, and Fx from the first
  //   surface

  params._F = f_1 - interp_max( f, (_num_indicators-1) );

  interp_max_D( f, fx, (_num_indicators-1), params._Fx );
  params._Fx = fx_1 - params._Fx;

  params._sf = sf;

  delete [] f;
  delete [] fx;

  return true;
}

//------------------------------------------------------------------------
// Function    : interp_max()
// Description : 
//------------------------------------------------------------------------
float IOScalarField::interp_max( const float *f, int num_values ) const
{
  float max_value = 0.0;
  float val;
  for ( int i = 0; i < num_values; i++ )
  {
    val = f[i];
    for ( int j = 1; j < num_values; j++ )
      val *= 1.0 + sign_function( f[i] - f[(i+j)%num_values] );

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

float IOScalarField::sign_function( float a ) const
{
  return ( a / sqrt( a*a + EPSILON_SQRD ));
}

void IOScalarField::sign_function_D( float a, const vec<3> &aD,
                                     vec<3> &D ) const
{
  float isqrt_a_sqrd = 1.0/sqrt( a*a + EPSILON_SQRD );
  D = aD * 
    ( isqrt_a_sqrd - a*a*(isqrt_a_sqrd*isqrt_a_sqrd*isqrt_a_sqrd) );
}


float IOScalarField::power_function( float a, int b ) const
{
  float val = 1.0;
  for ( int i = 0; i < b; i++ )
    val *= a;
  return val;
}



////------------------------------------------------------------------------
//// Function    : approx_max()
//// Description : 
////------------------------------------------------------------------------
//float IOScalarField::approx_max( float a, float b ) const
//{
//  return 
//    ( 0.5*(a + b + sqrt( (a-b)*(a-b) + SUPER_EPSILON*SUPER_EPSILON )) );
//
//}
//
//void IOScalarField::approx_max_D( float a, float b, const vector_type &aD,
//                                  const vector_type &bD, vector_type &D ) const
//{
//  D = 0.5*(aD + bD + 
//    ((a-b)/sqrt((a-b)*(a-b) + SUPER_EPSILON*SUPER_EPSILON))*(aD - bD));
//}

////------------------------------------------------------------------------
//// Function    : erf_max()
//// Description : 
////------------------------------------------------------------------------
//float IOScalarField::erf_max( float a, float b ) const
//{
//  return 
//    ( 0.5*( (1+erff((b-a)/SUPER_EPSILON))*b + 
//            (1-erff((b-a)/SUPER_EPSILON))*a ) );
//
//}
//
//void IOScalarField::erf_max_D( float a, float b, const vector_type &aD,
//                               const vector_type &bD, vector_type &D ) const
//{
//  float b_a = b - a;
//  float b_a_epsilon = b_a / SUPER_EPSILON;
//
//  D = ((2.0/sqrt(PI))*exp( -b_a_epsilon*b_a_epsilon ) * (b-a) +
//       erff( b_a_epsilon)) * (bD - aD) + bD + aD;
//  D *= 0.5;
//}
