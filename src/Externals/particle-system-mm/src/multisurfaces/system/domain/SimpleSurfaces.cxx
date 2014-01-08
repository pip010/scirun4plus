#include <iostream>
#include <cstdlib>
#include <system/domain/SimpleSurfaces.h>

// 4244 -- conversion of double to float warning
#ifdef _WIN32
#pragma warning( disable : 4244 )
#pragma warning( disable : 4305 )
#endif

using namespace particle_sys;
using namespace std;

/************************************************************************/
//                            SPHERE                                    //
/************************************************************************/

//------------------------------------------------------------------------
// Function    : SimpleSphere::SphereIF()
// Description : initializing constructor
//------------------------------------------------------------------------
SimpleSphere::SimpleSphere( float radius ) : Surface()
{
  _radius = radius;
}

//------------------------------------------------------------------------
// Function    : SimpleSphere::computeSurfacePointParams()
// Description : Use to compute/set F Fx H at the particle pos.
//------------------------------------------------------------------------
bool SimpleSphere::computeSurfacePointParams( const vector_type &pos, 
                                              SurfacePointParams &params,
                                              bool computeHessian )
                                        const
{
#ifdef TWO_D
  params._F = (pos[0])*(pos[0]) + (pos[1])*(pos[1]) - 
    _radius*_radius - _isovalue;  

  // the first derivative
  params._Fx[0] = 2.0*pos[0];
  params._Fx[1] = 2.0*pos[1];

  // the second derivative (Hessian matrix)
  params._Fxx.identity();
  params._Fxx *= 2.0;
#endif

#ifdef THREE_D

  params._F = (pos[0])*(pos[0]) + (pos[1])*(pos[1]) + (pos[2])*(pos[2]) - 
    _radius*_radius - _isovalue;  

  // the first derivative
  params._Fx[0] = 2.0*pos[0];
  params._Fx[1] = 2.0*pos[1];
  params._Fx[2] = 2.0*pos[2];

  // the second derivative (Hessian matrix)
  params._Fxx.identity();
  params._Fxx *= 2.0;
#endif

  return true;
}

/************************************************************************/
//                            TORUS                                     //
/************************************************************************/

//------------------------------------------------------------------------
// Function    : SimpleTorus::init()
// Description : initialize the sphere
//------------------------------------------------------------------------
SimpleTorus::SimpleTorus( float r, float R ) : Surface()
{
  _r = r; _R = R;
}

//------------------------------------------------------------------------
// Function    : SimpleTorus::computeF_Fx_H()
// Description : Use to compute/set F Fx H at the particle pos.
//------------------------------------------------------------------------
bool SimpleTorus::computeSurfacePointParams( const vector_type &pos, 
                                             SurfacePointParams &params,
                                             bool computeHessian )
                                       const
{
#ifdef TWO_D
  float x = pos[0],
        y = pos[1];
  float xy_sqrt = sqrt(x*x + y*y);
  
  params._F = (xy_sqrt - _R)*(xy_sqrt - _R) - _r*_r;

  // the first derivative
  params._Fx[0] = 2.0 * (xy_sqrt - _R) * (x / xy_sqrt);
  params._Fx[1] = 2.0 * (xy_sqrt - _R) * (y / xy_sqrt);

  // the second derivative (Hessian matrix)
  float xy_sqrt_cubed = xy_sqrt*xy_sqrt*xy_sqrt;
  float xy = xy_sqrt*xy_sqrt;

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
#endif

#ifdef THREE_D
  float x = pos[0],
        y = pos[1],
        z = pos[2];
  float xy_sqrt = sqrt(x*x + y*y);
  
  params._F = z*z + (xy_sqrt - _R)*(xy_sqrt - _R) - _r*_r;

  // the first derivative
  params._Fx[0] = 2.0 * (xy_sqrt - _R) * (x / xy_sqrt);
  params._Fx[1] = 2.0 * (xy_sqrt - _R) * (y / xy_sqrt);
  params._Fx[2] = 2.0 * z;  

  // the second derivative (Hessian matrix)
  float xy_sqrt_cubed = xy_sqrt*xy_sqrt*xy_sqrt;
  float xy = xy_sqrt*xy_sqrt;

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
#endif

  return true;
}

/************************************************************************/
//                         HOLLOW CUBE                                  //
/************************************************************************/

//------------------------------------------------------------------------
// Function    : SimpleQuarticCube::HollowCubeIF()
// Description : initializing constructor
//------------------------------------------------------------------------
SimpleQuarticCube::SimpleQuarticCube( float radius ) : Surface()
{
  _radius = radius;
  _isovalue = 0.0;
}

//------------------------------------------------------------------------
// Function    : SimpleQuarticCube::computeF_Fx_H()
// Description : Use to compute/set F Fx H at the particle pos.
//------------------------------------------------------------------------
bool SimpleQuarticCube::computeSurfacePointParams( const vector_type &pos, 
                                                   SurfacePointParams &params,
                                                   bool computeHessian )
                                             const
{
#ifdef TWO_D
  float radius_sqrd = _radius*_radius;

  params._F = pow(pos[0],4) + pow(pos[1],4) -
      10.0*radius_sqrd*pos[0]*pos[0] -
      10.0*radius_sqrd*pos[1]*pos[1] -_isovalue;

  // the first derivative
  params._Fx[0] = 4.0*pow(pos[0],3) - 20.0*radius_sqrd*pos[0];
  params._Fx[1] = 4.0*pow(pos[1],3) - 20.0*radius_sqrd*pos[1];
  
  // the second derivative (Hessian matrix)
  params._Fxx = 0.0;
  params._Fxx(0,0) = 12.0*pos[0]*pos[0] - 20.0*radius_sqrd;
  params._Fxx(1,1) = 12.0*pos[1]*pos[1] - 20.0*radius_sqrd;
#endif

#ifdef THREE_D
  float radius_sqrd = _radius*_radius;

  params._F = pow(pos[0],4) + pow(pos[1],4) + pow(pos[2],4) -
      10.0*radius_sqrd*pos[0]*pos[0] -
      10.0*radius_sqrd*pos[1]*pos[1] -
      10.0*radius_sqrd*pos[2]*pos[2] - _isovalue;

  // the first derivative
  params._Fx[0] = 4.0*pow(pos[0],3) - 20.0*radius_sqrd*pos[0];
  params._Fx[1] = 4.0*pow(pos[1],3) - 20.0*radius_sqrd*pos[1];
  params._Fx[2] = 4.0*pow(pos[2],3) - 20.0*radius_sqrd*pos[2];
  
  // the second derivative (Hessian matrix)
  params._Fxx = 0.0;
  params._Fxx(0,0) = 12.0*pos[0]*pos[0] - 20.0*radius_sqrd;
  params._Fxx(1,1) = 12.0*pos[1]*pos[1] - 20.0*radius_sqrd;
  params._Fxx(2,2) = 12.0*pos[2]*pos[2] - 20.0*radius_sqrd;
#endif

  return true;
}



/************************************************************************/
//                            ELLIPSE                                   //
/************************************************************************/

//------------------------------------------------------------------------
// Function    : SimpleEllipse::SphereIF()
// Description : initializing constructor
//------------------------------------------------------------------------
SimpleEllipse::SimpleEllipse( float a, float b, float c, float d, float e, float f ) 
: Surface()
{
  _a = a;  _b = b;  _c = c;  _d = d;  _e = e;  _f = f; 
}

//------------------------------------------------------------------------
// Function    : SimpleEllipse::computeSurfacePointParams()
// Description : Use to compute/set F Fx H at the particle pos.
//------------------------------------------------------------------------
bool SimpleEllipse::computeSurfacePointParams( const vector_type &pos, 
                                               SurfacePointParams &params,
                                               bool computeHessian )
                                        const
{
#ifdef TWO_D
  params._F = _a*(pos[0])*(pos[0]) + _b*(pos[0])*(pos[1]) + 
    _c*(pos[1])*(pos[1]) + _d*pos[0] + _e*pos[1] + _f - _isovalue;  

  // the first derivative
  params._Fx[0] = 2.0*_a*pos[0] + _b*pos[1] + _d;
  params._Fx[1] = _b*pos[0] + 2.0*_c*pos[1] + _e;

  // the second derivative (Hessian matrix)
  params._Fxx(0,0) = 2.0*_a;
  params._Fxx(0,1) = params._Fxx(1,0) = _b;
  params._Fxx(1,1) = 2.0*_c;
#endif

#ifdef THREE_D
#endif

  return true;
}



/************************************************************************/
//                            COSNTANT                                  //
/************************************************************************/

//------------------------------------------------------------------------
// Function    : Constant::computeSurfacePointParams()
// Description : Use to compute/set F Fx H at the particle pos.
//------------------------------------------------------------------------
bool Constant::computeSurfacePointParams( const vector_type &pos, 
                                          SurfacePointParams &params,
                                          bool computeHessian ) const
{
  params._F = 0.0;

  params._Fx = 0.0;
   
  params._Fxx.identity();

  return true;
}
