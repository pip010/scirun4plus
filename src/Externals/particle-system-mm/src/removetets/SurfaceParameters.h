#ifndef __SURFACE_PARAMETERS_H__
#define __SURFACE_PARAMETERS_H__

#include "mtxlib.h"
#include "defines.h"

struct SurfacePointParams
{
  float       _F, _curvature_mag, _sf;
  vector_type _Fx;
  matrix_type _Fxx, _curvature;
  float _kappa1, _kappa2;
  vector_type _eigenvector1;

  SurfacePointParams() {};
  const SurfacePointParams& operator = (const SurfacePointParams &that)
  { _F=that._F; _Fx=that._Fx; _Fxx=that._Fxx; 
  _curvature=that._curvature; _curvature_mag=that._curvature_mag;
  _kappa1 = that._kappa1; _kappa2 = that._kappa2;
  _eigenvector1 = that._eigenvector1;
  _sf = that._sf;
  return *this; };
};

#endif // __SURFACE_PARAMETERS_H__
