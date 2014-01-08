#include <iostream>
#include <cstdlib>
#include <float.h>
#include "ScalarField.h"

using namespace std;

/************************************************************************/
//                     SCALAR FIELD BASE CLASS                          //
/************************************************************************/

//------------------------------------------------------------------------
// Function    : computeCurvatureMagnitude()
// Description : compute the magnitude of the curvature 
//------------------------------------------------------------------------
template <>
mtx_datatype ScalarField<2>::computeCurvatureMagnitude( const matrix<2,2> &G,
                                                 ScalarFieldParams<2>
                                                 &params ) const
{
  // the trace and norm SHOULD be the same for 2D 
  return fabs( G.trace() );
}

template <>
mtx_datatype ScalarField<3>::computeCurvatureMagnitude( const matrix<3,3> &G,
                                                 ScalarFieldParams<3>
                                                 &params ) const
{
  matrix<3,3> GT(G); 
  GT.transpose();
  
  // solve for kappa1 and kappa2
  mtx_datatype T = G.trace();
  mtx_datatype F = sqrt((G*GT).trace());

  vec<2> curv;
  curv[0] = ( T + sqrt(abs(2.0 * F * F - T * T))) / 2.0;  
  curv[1] = ( T - sqrt(abs(2.0 * F * F - T * T))) / 2.0;

  params._kappa1 = curv[0];
  params._kappa2 = curv[1];

  // take the abs() of the kappa values because curvature can be negative,
  //   and we only want the magnitude! (remember, negative curvature 
  //   indicates concavities --> ++ is concave, +- is saddle, -- is sink)
  return ( sqrt(curv[0]*curv[0] + curv[1]*curv[1]) );
}

