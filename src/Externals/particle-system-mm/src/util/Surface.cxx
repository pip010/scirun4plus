#include <iostream>
#include <cstdlib>
#include <Surface.h>

using namespace std;

/************************************************************************/
//                       SURFACE BASE CLASS                             //
/************************************************************************/

bool Surface::inBounds( const vec<3> &pos ) const
{
  if ( (pos(0) >= _start(0)) && (pos(0) <= _end(0)) &&
       (pos(1) >= _start(1)) && (pos(1) <= _end(1)) &&
       (pos(2) >= _start(2)) && (pos(2) <= _end(2)) )
    return true;
  else
    return false;
}

//------------------------------------------------------------------------
// Function    : Surface::computeCurvature()
// Description : compute the curvature at the point pos -- uses the algo
//               from gordon's vis03 paper
//------------------------------------------------------------------------
float Surface::computeCurvature( SurfacePointParams &params ) const
{
  // get the normalized gradient at this point on the surface
  // TODO: i'm not sure why this is negative -- will have to look into 
  //       this later!!!!
  vec<3> n = -params._Fx;
  n.normalize();

  // get the tangent plane projection matrix
  matrix<3,3> P; P.identity();
  P -= DirectProduct( n, n ); 

  // store the curvature matrix that has been projected into the 
  //   tangent plane
  matrix<3,3> G = (P*params._Fxx) / 
    (params._Fx.length()+1.0e-6f); // positive bc inside is lower values!

  // TODO : not sure if we want to set the curvature matrix in this step
  //        or the next (do want the sigma values or not?)
  //params._curvature = G;

  // isolate the curvature values
  G = G*P; 
  params._curvature_mag = computeCurvatureMagnitude( G, params );

  return params._curvature_mag;
}

//------------------------------------------------------------------------
// Function    : Surface::computeCurvatureMagnitude()
// Description : compute the magnitude of the curvature --> this is an
//               overloaded function to deal with 2D or 3D
//------------------------------------------------------------------------
float Surface::computeCurvatureMagnitude( const matrix<3,3> &G,
                                          SurfacePointParams &params ) const
{
  matrix<3,3> GT(G); 
  GT.transpose();
  
  // solve for kappa1 and kappa2
  float T = G.trace();
  float F = sqrt((G*GT).trace());

  vec<2> curv;
  curv[0] = ( T + sqrt(abs(2.0f * F * F - T * T))) / 2.0f;  
  curv[1] = ( T - sqrt(abs(2.0f * F * F - T * T))) / 2.0f;

  // take the abs() of the kappa values because curvature can be negative,
  //   and we only want the magnitude! (remember, negative curvature 
  //   indicates concavities --> ++ is concave, +- is saddle, -- is sink)
  //return ( sqrt(curv[0]*curv[0] + curv[1]*curv[1]) );
  return ( max( fabs(curv[0]), fabs(curv[1]) ) );
}
