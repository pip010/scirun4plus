#include <math.h>
#include <iostream>
#include <system/domain/Energy.h>

#ifdef _WIN32
#pragma warning( disable : 4244 4305 )
#endif

using namespace particle_sys;
using namespace std;


/************************************************************************/
//                        ENERGY BASE CLASS                             //
/************************************************************************/






/************************************************************************/
//                           COTAN ENERGY                               //
/************************************************************************/

//------------------------------------------------------------------------
// Function    : Constructor
// Description : 
//------------------------------------------------------------------------
CotanEnergy::CotanEnergy(float nid) : Energy(nid)
{
  // compute the ideal energy
  vector_type d( 
    sqrt(_normalized_ideal_distance*_normalized_ideal_distance/3.0) );

  initializeParameters( d );

  _ideal_energy = solveEnergy();
}

//------------------------------------------------------------------------
// Function    : initializeParameters
// Description : 
//------------------------------------------------------------------------
void CotanEnergy::initializeParameters( const vector_type &d_ij )
{
  _d = d_ij.length();
  _d_ij_unit = d_ij / _d;
}

//------------------------------------------------------------------------
// Function    : CotanEnergy::solveEnergy(float distance)
// Description : gives you the scalar energy between two particles
//------------------------------------------------------------------------
float CotanEnergy::solveEnergy()  
{
  if ( _d >= 1.0 )
    return 0.0;

  float theta = _d * PI * 0.5;

  _e = cos(theta)/(EPSILON + sin(theta)) + theta - PI*0.5;
  return _e;
}

//------------------------------------------------------------------------
// Function    : CotanEnergy::solveForce(float distance)
// Description : gives you the scalar force between two particles
//               --> NOTE : this force is in the negative r_ij direction!
//------------------------------------------------------------------------
vector_type CotanEnergy::solveForce()  
{
  if ( _d >= 1.0 )
  {
    vector_type ret(0.0);
    return ret;
  }

  _f = 1.0 / (EPSILON + sin( _d * PI*0.5));
  _f *= -_f;
  ++_f;
  _f *= PI*0.5;

  return ( _f * _d_ij_unit );
}

//------------------------------------------------------------------------
// Function    : CotanEnergy::solveYank(float distance)
// Description : gives you the acceleration between two particles
//               NOTE: assumes that you have already computed the force!!!
//------------------------------------------------------------------------
matrix_type CotanEnergy::solveYank()  
{
  //if ( (_r_ij_over_sigma1 == 1.0) && (_r_ij_over_sigma2 == 1.0) )
  //{
  //  matrix_type ret(0.0);
  //  return ret;
  //}

  //float theta = _r_ij_over_sigma1 * PI * 0.5;
  //float sin_theta = sin(theta);
  //float yank1 = 2.0*(cos(theta)/ (EPSILON + sin_theta*sin_theta*sin_theta));
  //yank1 *= (PI/(2.0*_sigma1))*(PI/(2.0*_sigma1));

  //theta = _r_ij_over_sigma2 * PI * 0.5;
  //sin_theta = sin(theta);
  //float yank2 = 2.0*(cos(theta)/ (EPSILON + sin_theta*sin_theta*sin_theta));
  //yank2 *= (PI/(2.0*_sigma2))*(PI/(2.0*_sigma2));

  //_y = 0.5*(yank1 + yank2);

  //matrix_type rxr = DirectProduct( _r_ij_unit, _r_ij_unit );
  //matrix_type I; I.identity();

  //return ( _y*rxr + _f*(I-rxr) );



  matrix_type tmp;
  return tmp.identity();
}


/************************************************************************/
//                         RADIAL ENERGY                                //
/************************************************************************/

//------------------------------------------------------------------------
// Function    : Constructor
// Description : 
//------------------------------------------------------------------------
RadialEnergy::RadialEnergy(float nid) : Energy(nid)
{
  // compute the ideal energy
  vector_type d( 
    sqrt(_normalized_ideal_distance*_normalized_ideal_distance/3.0) );

  initializeParameters( d );

  _ideal_energy = solveEnergy();
}

//------------------------------------------------------------------------
// Function    : initializeParameters
// Description : 
//------------------------------------------------------------------------
void RadialEnergy::initializeParameters( const vector_type &d_ij )
{
  _d = d_ij.length();
  _d_ij_unit = d_ij / _d;
}

//------------------------------------------------------------------------
// Function    : RadialEnergy::solveEnergy(float distance)
// Description : gives you the scalar energy between two particles
//------------------------------------------------------------------------
float RadialEnergy::solveEnergy()  
{
  if ( _d >= 1.0 )
    return 0.0;

  _e = (1.0-_d) / (_d+EPSILON);
  _e *= _e;

  return _e;
}

//------------------------------------------------------------------------
// Function    : RadialEnergy::solveForce(float distance)
// Description : gives you the scalar force between two particles
//               --> NOTE : this force is in the negative r_ij direction!
//------------------------------------------------------------------------
vector_type RadialEnergy::solveForce() 
{
  if ( _d >= 1.0 )
  {
    vector_type ret(0.0);
    return ret;
  }

  float _f = (-2.0 / (_d*_d*_d + EPSILON)) * (1.0-_d);


  return (_f * _d_ij_unit);
}




