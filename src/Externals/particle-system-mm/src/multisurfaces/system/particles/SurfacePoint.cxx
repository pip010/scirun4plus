#include <cstdlib>
#include <iostream>
#include <system/domain/Domain.h>
#include <system/particles/SurfacePoint.h>

using namespace particle_sys;
using namespace std;

#pragma warning( disable : 4800 )

//------------------------------------------------------------------------
// Function    : updateSurfaceParameters
// Description : query the Surface to get the parameters at this point
//               location
//------------------------------------------------------------------------
bool SurfacePoint::updateSurfaceParameters()
{
  if ( !_domain->insideDomain( _position ) )
  {
    cout << "SurfacePoint::updateSurfaceParameters() --> " <<
      "Outside Domain!!!\n";
    return false;
  }

  // have the surface compute the surface parameters
  if ( !_domain->computeSurfacePointParams( _position, _surf_params,
                                            (bool)_angular_density ) )
    return false;

  // compute the curvature matrix from the Hessian
  _domain->computeCurvature( _surf_params );

  // set the Point's normal to be in the direction of the gradient
  _normal = _surf_params._Fx;
  _normal.normalize();

  return true;
}

void SurfacePoint::onlyQuerySurface()
{
  _domain->computeSurfacePointParams( _position, _surf_params,
                                            (bool)_angular_density );
}

//------------------------------------------------------------------------
// Function    : operator =
// Description : the assignment operator
//------------------------------------------------------------------------
const SurfacePoint& SurfacePoint::operator = (const SurfacePoint &that)
{
  // copy the Point variables
  (*this).pointAssignmentOp( that );

  _surf_params = that._surf_params;
  _max_sf = that._max_sf;
  _min_sf = that._min_sf;

  return *this;
}

void SurfacePoint::surfacePointAssignmentOp(const SurfacePoint &that) 
{ 
  *this = that; 
}

//------------------------------------------------------------------------
// Function    : operator << 
// Description : output function for the SurfacePoint class
//------------------------------------------------------------------------
ostream& operator << ( ostream& outs, const SurfacePoint& source ) 
{
  outs << (Point)source << "              " << 
    "<F,Fx,Fxx> : <" << source.F() <<
    ", (" << source.Fx() << "),\n" << source.Fxx();
  
  return outs; 
}











  

