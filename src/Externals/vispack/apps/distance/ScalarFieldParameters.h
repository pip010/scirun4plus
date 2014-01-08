#ifndef __SCALAR_FIELD_PARAMETERS_H__
#define __SCALAR_FIELD_PARAMETERS_H__

#include "mtxlib.h"

template <int dim> class ScalarFieldParams
{
public:
  float           _F, _curvature_mag, _kappa1, _kappa2;
  vec<dim>        _Fx;
  matrix<dim,dim> _Fxx, _curvature;

  ScalarFieldParams() {};
  const ScalarFieldParams& operator = (const ScalarFieldParams &that)
  { _F=that._F; _Fx=that._Fx; _Fxx=that._Fxx; 
    _curvature=that._curvature; _curvature_mag=that._curvature_mag;
    _kappa1=_kappa1; _kappa2=_kappa2;
    return *this; };
};


#endif // __SCALAR_FIELD_PARAMETERS_H__
