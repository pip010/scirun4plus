//----------------------------------------------------------------------//
// FILE        : ScalarField.h                                                
// DESCRIPTION : Reconstructs a scalar field using spline kernels. 
//----------------------------------------------------------------------//


#ifndef __SCALAR_FIELD_H__
#define __SCALAR_FIELD_H__

#include <utilExports.h>

#include <Surface.h>
#include <mtxlib.h>
#include <multiDarrays.h>
#include <Kernel.h>
#include <defines.h>

#include <vector>

#if defined(_WIN32) && defined(_USRDLL)
template class MM_util_SHARE array3D<float>;
#endif

class MM_util_SHARE ScalarField : public Surface
{
public:
  ScalarField(const char *filename, int k, float sx=1.0,
              float sy=1.0, float sz=1.0);
  ~ScalarField();

  enum { BSPLINE, CATMULLROM, LINEAR };

  bool computeSurfacePointParams(const vector_type &pos, 
                                 SurfacePointParams &params,
                                 bool computeGradient=false,
                                 bool computeHessian=false);

  int xdim() const { return _xdim; }
  int ydim() const { return _ydim; }
  int zdim() const { return _zdim; }
  
  DataCenter xcenter() const { return _centers[0]; }
  DataCenter ycenter() const { return _centers[1]; }
  DataCenter zcenter() const { return _centers[2]; }

private:
  int _xdim, _ydim, _zdim;
  array3D<float> _field;
  Kernel *_kernel, *_kernelD, *_kernelDD;

  float _scale[3];
  DataCenter _centers[3];

  void readFile(const char *filename);

  // a simple linear interpolation
  inline float interpolate(float t, float a, float b) const
  { return ( (1.0f - t)*a + t*b ); };

};

#endif // __SCALAR_FIELD_H__
