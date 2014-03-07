#include <cstdlib>
#include "IOScalarField.h"

using namespace std;

//------------------------------------------------------------------------
// Function    : computeSurfacePointParams()
// Description : 
//------------------------------------------------------------------------
template<>
void IOScalarField<3>::computeScalarFieldParams(
  const vec<3> &pos, ScalarFieldParams<3> &params,
  bool computeHessian ) const
{
  // store the values for the main indicator function
  _indicators[_main_indicator]->computeScalarFieldParams( pos,
                                                          params,
                                                          computeHessian );
  float f_1 = params._F;
  vec<3> fx_1 = params._Fx;
  matrix<3,3> fxx_1;
  if ( computeHessian )
    fxx_1 = params._Fxx;

  float *f = new float[_num_indicators-1];
  vec<3> *fx = new vec<3>[_num_indicators-1];
  matrix<3,3> *fxx;
  if ( computeHessian )
    fxx = new matrix<3,3>[_num_indicators-1];
  int counter=0;
  for ( int i = 0; i < _num_indicators; i++ )
  {
    if ( i == _main_indicator )
      continue;

    _indicators[i]->computeScalarFieldParams( pos, params, 
                                              computeHessian );
    f[counter] = params._F;
    fx[counter] = params._Fx;

    if ( computeHessian )
      fxx[counter] = params._Fxx;

    ++counter;
   }

  params._F = f_1 - interp_max( f, (_num_indicators-1) );

  if ( !computeHessian )
  {
    interp_max_D( f, fx, (_num_indicators-1), params._Fx );
    params._Fx = fx_1 - params._Fx;
  }
  else
  {
    interp_max_D_DD( f, fx, fxx, (_num_indicators-1), 
                     params._Fx, params._Fxx );
    params._Fx = fx_1 - params._Fx;
    params._Fxx = fxx_1 - params._Fxx;
  }

  delete [] f;
  delete [] fx;
  if ( computeHessian )
    delete [] fxx;
}






