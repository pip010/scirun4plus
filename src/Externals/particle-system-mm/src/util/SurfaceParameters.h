#ifndef __SURFACE_PARAMETERS_H__
#define __SURFACE_PARAMETERS_H__

#include <mtxlib.h>
#include <defines.h>

struct SurfacePointParams
{
  float       _F, _curvature_mag;
  vector_type _Fx;
  matrix_type _Fxx;

  SurfacePointParams() {};
  const SurfacePointParams& operator = (const SurfacePointParams &that)
  { _F=that._F; _curvature_mag = that._curvature_mag;
  _Fx = that._Fx; _Fxx = that._Fxx;
  return *this; };
};

#endif // __SURFACE_PARAMETERS_H__
