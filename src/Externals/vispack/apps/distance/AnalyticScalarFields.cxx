#include <iostream>
#include <cstdlib>
#include "AnalyticScalarFields.h"

using namespace std;



void SpherePlane::computeScalarFieldParams( const vec<3> &pos, 
                                          ScalarFieldParams<3> &params,
                                          bool computeHessian )
                                        const
{
  mtx_datatype x = pos[0] - 0.5, y = pos[1] - 0.5, z = pos[2] - 0.5;
  mtx_datatype sphere_dist = x*x + 9.0*y*y + z*z - _radius*_radius, 
    plane_dist = (_z - z);
  
  if (sphere_dist > plane_dist)
    {
      params._F = plane_dist;
      params._Fx[0] = 0.0;
      params._Fx[1] = 0.0;
      params._Fx[2] = 1.0;
      if ( computeHessian )
	{
	  params._Fxx = 0.0f;
	}
    }
  else
    {
      params._F = sphere_dist;
      params._Fx[0] = 2.0*x;
      params._Fx[1] = 2.0*y;
      params._Fx[2] = 2.0*z;

      // the second derivative (Hessian matrix)
      if ( computeHessian )
	{
	  params._Fxx.identity();
	  params._Fxx *= 2.0;
	}
    }
}

  



/************************************************************************/
//                            SPHERE                                    //
/************************************************************************/


//------------------------------------------------------------------------
// Function    : computeScalarFieldParams()
// Description : Use to compute/set F Fx H at the particle pos.
//------------------------------------------------------------------------
template <>
void Sphere<2>::computeScalarFieldParams( const vec<2> &pos, 
                                          ScalarFieldParams<2> &params,
                                          bool computeHessian )
                                        const
{
  params._F = (pos[0])*(pos[0]) + (pos[1])*(pos[1]) - 
    _radius*_radius;  

  // the first derivative
  params._Fx[0] = 2.0*pos[0];
  params._Fx[1] = 2.0*pos[1];

  // the second derivative (Hessian matrix)
  if ( computeHessian )
  {
    params._Fxx.identity();
    params._Fxx *= 2.0;
  }
}

template <>
void Sphere<3>::computeScalarFieldParams( const vec<3> &pos, 
                                          ScalarFieldParams<3> &params,
                                          bool computeHessian )
                                        const
{
  params._F = (pos[0])*(pos[0]) + (pos[1])*(pos[1]) + (pos[2])*(pos[2]) - 
    _radius*_radius;  

  // the first derivative
  params._Fx[0] = 2.0*pos[0];
  params._Fx[1] = 2.0*pos[1];
  params._Fx[2] = 2.0*pos[2];

  // the second derivative (Hessian matrix)
  if ( computeHessian )
  {
    params._Fxx.identity();
    params._Fxx *= 2.0;
  }
}

/************************************************************************/
//                            TORUS                                     //
/************************************************************************/

//------------------------------------------------------------------------
// Function    : computeScalarFieldParams()
// Description : Use to compute/set F Fx H at the particle pos.
//------------------------------------------------------------------------
template <>
void Torus<2>::computeScalarFieldParams( const vec<2> &pos, 
                                         ScalarFieldParams<2> &params,
                                         bool computeHessian )
                                       const
{
  mtx_datatype x = pos[0],
        y = pos[1];
  mtx_datatype xy_sqrt = sqrt(x*x + y*y);
  
  params._F = (xy_sqrt - _R)*(xy_sqrt - _R) - _r*_r;

  // the first derivative
  params._Fx[0] = 2.0 * (xy_sqrt - _R) * (x / xy_sqrt);
  params._Fx[1] = 2.0 * (xy_sqrt - _R) * (y / xy_sqrt);

  // the second derivative (Hessian matrix)
  if ( computeHessian )
  {
    mtx_datatype xy_sqrt_cubed = xy_sqrt*xy_sqrt*xy_sqrt;
    mtx_datatype xy = xy_sqrt*xy_sqrt;

    params._Fxx = 0.0;

    params._Fxx(0,0) = 2.0*x*x/(xy) +
      2.0*(xy_sqrt - _R)/xy_sqrt -
      2.0*x*x*(xy_sqrt - _R)/(xy_sqrt_cubed);
    params._Fxx(1,1) = 2.0*y*y/(xy) +
      2.0*(xy_sqrt - _R)/xy_sqrt -
      2.0*y*y*(xy_sqrt - _R)/(xy_sqrt_cubed);
  
    params._Fxx(0,1) = 2.0*x*y/(xy) -
      2.0*x*y*(xy_sqrt - _R)/(xy_sqrt_cubed);
    params._Fxx(1,0) = params._Fxx(0,1);
  }
}

template <>
void Torus<3>::computeScalarFieldParams( const vec<3> &pos, 
                                         ScalarFieldParams<3> &params,
                                         bool computeHessian )
                                       const
{
  mtx_datatype x = pos[0],
        y = pos[1],
        z = pos[2];
  mtx_datatype xy_sqrt = sqrt(x*x + y*y);
  
  params._F = z*z + (xy_sqrt - _R)*(xy_sqrt - _R) - _r*_r;

  // the first derivative
  params._Fx[0] = 2.0 * (xy_sqrt - _R) * (x / xy_sqrt);
  params._Fx[1] = 2.0 * (xy_sqrt - _R) * (y / xy_sqrt);
  params._Fx[2] = 2.0 * z;  

  // the second derivative (Hessian matrix)
  if ( computeHessian )
  {
    mtx_datatype xy_sqrt_cubed = xy_sqrt*xy_sqrt*xy_sqrt;
    mtx_datatype xy = xy_sqrt*xy_sqrt;

    params._Fxx = 0.0;

    params._Fxx(0,0) = 2.0*x*x/(xy) +
      2.0*(xy_sqrt - _R)/xy_sqrt -
      2.0*x*x*(xy_sqrt - _R)/(xy_sqrt_cubed);
    params._Fxx(1,1) = 2.0*y*y/(xy) +
      2.0*(xy_sqrt - _R)/xy_sqrt -
      2.0*y*y*(xy_sqrt - _R)/(xy_sqrt_cubed);
    params._Fxx(2,2) = 2.0;

  
    params._Fxx(0,1) = 2.0*x*y/(xy) -
      2.0*x*y*(xy_sqrt - _R)/(xy_sqrt_cubed);
    params._Fxx(1,0) = params._Fxx(0,1);
  }
}

/************************************************************************/
//                         QUARTIC CUBE                                 //
/************************************************************************/

//------------------------------------------------------------------------
// Function    : computeScalarFieldParams
// Description : Use to compute/set F Fx H at the particle pos.
//------------------------------------------------------------------------
template <>
void QuarticCube<2>::computeScalarFieldParams( const vec<2> &pos, 
                                               ScalarFieldParams<2> &params,
                                               bool computeHessian )
                                             const
{
  mtx_datatype isovalue = 0.1;
  mtx_datatype radius_sqrd = _radius*_radius;

  params._F = pow(pos[0],4) + pow(pos[1],4) -
      10.0*radius_sqrd*pos[0]*pos[0] -
      10.0*radius_sqrd*pos[1]*pos[1] - isovalue;

  // the first derivative
  params._Fx[0] = 4.0*pow(pos[0],3) - 20.0*radius_sqrd*pos[0];
  params._Fx[1] = 4.0*pow(pos[1],3) - 20.0*radius_sqrd*pos[1];
  
  // the second derivative (Hessian matrix)
  if ( computeHessian )
  {
    params._Fxx = 0.0;
    params._Fxx(0,0) = 12.0*pos[0]*pos[0] - 20.0*radius_sqrd;
    params._Fxx(1,1) = 12.0*pos[1]*pos[1] - 20.0*radius_sqrd;
  }
}

template <>
void QuarticCube<3>::computeScalarFieldParams( const vec<3> &pos, 
                                               ScalarFieldParams<3> &params,
                                               bool computeHessian )
                                             const
{
  mtx_datatype isovalue = 0.1;
  mtx_datatype radius_sqrd = _radius*_radius;

  params._F = pow(pos[0],4) + pow(pos[1],4) + pow(pos[2],4) -
      10.0*radius_sqrd*pos[0]*pos[0] -
      10.0*radius_sqrd*pos[1]*pos[1] -
      10.0*radius_sqrd*pos[2]*pos[2] - isovalue;

  // the first derivative
  params._Fx[0] = 4.0*pow(pos[0],3) - 20.0*radius_sqrd*pos[0];
  params._Fx[1] = 4.0*pow(pos[1],3) - 20.0*radius_sqrd*pos[1];
  params._Fx[2] = 4.0*pow(pos[2],3) - 20.0*radius_sqrd*pos[2];
  
  // the second derivative (Hessian matrix)
  if ( computeHessian )
  {
    params._Fxx = 0.0;
    params._Fxx(0,0) = 12.0*pos[0]*pos[0] - 20.0*radius_sqrd;
    params._Fxx(1,1) = 12.0*pos[1]*pos[1] - 20.0*radius_sqrd;
    params._Fxx(2,2) = 12.0*pos[2]*pos[2] - 20.0*radius_sqrd;
  }
}
