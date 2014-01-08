#include <iostream>
#include <cstdlib>
#include <float.h>
#include "defines.h"
#include "Surface.h"

#ifdef _WIN32
#pragma warning( disable : 4244 )
#endif

using namespace std;

/************************************************************************/
//                       SURFACE BASE CLASS                             //
/************************************************************************/

//------------------------------------------------------------------------
// Function    : Surface::Surface()
// Description : constructor
//------------------------------------------------------------------------
Surface::Surface()
{
  _isovalue = 0.0;

  _d_start = -1.0;
  _d_end = 1.0;
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
  vector_type n = -params._Fx;
  n.normalize();

  // get the tangent plane projection matrix
  matrix_type P; P.identity();
  P -= DirectProduct( n, n ); 

  // store the curvature matrix that has been projected into the 
  //   tangent plane
  matrix_type G = (P*params._Fxx) / 
    (params._Fx.length()+EPSILON); // positive bc inside is lower values!

  // TODO : not sure if we want to set the curvature matrix in this step
  //        or the next (do want the sigma values or not?)
  //params._curvature = G;

  // isolate the curvature values
  G = G*P; 
  params._curvature = G;
  params._curvature_mag = computeCurvatureMagnitude( G, params );

  return params._curvature_mag;
}

float Surface::computeSignedCurvature( SurfacePointParams &params ) const
{
  // get the normalized gradient at this point on the surface
  // TODO: i'm not sure why this is negative -- will have to look into 
  //       this later!!!!
  vector_type n = -params._Fx;
  n.normalize();

  // get the tangent plane projection matrix
  matrix_type P; P.identity();
  P -= DirectProduct( n, n ); 

  // store the curvature matrix that has been projected into the 
  //   tangent plane
  matrix_type G = (P*params._Fxx) / 
    (params._Fx.length()+EPSILON); // positive bc inside is lower values!

  // TODO : not sure if we want to set the curvature matrix in this step
  //        or the next (do want the sigma values or not?)
  //params._curvature = G;

  // isolate the curvature values
  G = G*P; 
  params._curvature = G;
  params._curvature_mag = computeSignedCurvatureMagnitude( G, params );

  return params._curvature_mag;
}

//------------------------------------------------------------------------
// Function    : Surface::computeCurvatureMagnitude()
// Description : compute the magnitude of the curvature --> this is an
//               overloaded function to deal with 2D or 3D
//------------------------------------------------------------------------
float Surface::computeCurvatureMagnitude( const matrix<2,2> &G,
                                          SurfacePointParams &params) const
{
  // the trace and norm SHOULD be the same for 2D --> maybe?
  return fabs( G.trace() );
}

float Surface::computeCurvatureMagnitude( const matrix<3,3> &G,
                                          SurfacePointParams &params ) const
{
  matrix<3,3> GT(G); 
  GT.transpose();
  
  // solve for kappa1 and kappa2
  float T = G.trace();
  float F = sqrt((G*GT).trace());

  vec<2> curv;
  curv[0] = ( T + sqrt(abs(2.0 * F * F - T * T))) / 2.0;  
  curv[1] = ( T - sqrt(abs(2.0 * F * F - T * T))) / 2.0;

  params._kappa1 = curv[0];
  params._kappa2 = curv[1];

  computeEigen( G, params );

  // take the abs() of the kappa values because curvature can be negative,
  //   and we only want the magnitude! (remember, negative curvature 
  //   indicates concavities --> ++ is concave, +- is saddle, -- is sink)
  return ( sqrt(curv[0]*curv[0] + curv[1]*curv[1]) );
}

float Surface::computeSignedCurvatureMagnitude( const matrix<2,2> &G,
                                                SurfacePointParams &params ) 
  const
{
  // the trace and norm SHOULD be the same for 2D --> maybe?
  return G.trace();
}

float Surface::computeSignedCurvatureMagnitude( const matrix<3,3> &G,
                                                SurfacePointParams &params ) 
  const
{
  matrix<3,3> GT(G); 
  GT.transpose();
  
  // solve for kappa1 and kappa2
  float T = G.trace();
  float F = sqrt((G*GT).trace());

  vec<2> curv;
  curv[0] = ( T + sqrt(abs(2.0 * F * F - T * T))) / 2.0;  
  curv[1] = ( T - sqrt(abs(2.0 * F * F - T * T))) / 2.0;

  params._kappa1 = curv[0];
  params._kappa2 = curv[1];

  // take the abs() of the kappa values because curvature can be negative,
  //   and we only want the magnitude! (remember, negative curvature 
  //   indicates concavities --> ++ is concave, +- is saddle, -- is sink)
  return ( curv[0] + curv[1] );
}

//------------------------------------------------------------------------
// Function    : Surface::computeEigen()
// Description : 
//------------------------------------------------------------------------
void Surface::computeEigen( const matrix<2,2> &G,
                            SurfacePointParams &params) const
{
}

void Surface::computeEigen( const matrix<3,3> &G,
                            SurfacePointParams &params) const
{
  float a = G(0,0);
  float b = G(0,1);
  float c = G(0,2);
  float d = G(1,1);
  float e = G(1,2);

  float lambda = params._kappa1;

  float x3 = 1.0;
  float x2 = (b*c + e*(lambda-a))/((lambda-d)*(lambda-a)-b*b);
  float x1 = (b*x2 + c)/(lambda-a);

#ifdef THREE_D
  params._eigenvector1.set( x1, x2, x3 );
#endif
  params._eigenvector1.normalize();
}
