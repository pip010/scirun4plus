//----------------------------------------------------------------------//
// FILE        : ScalarField.h                                                
// DESCRIPTION : 
//----------------------------------------------------------------------//


#ifndef __SCALAR_FIELD_H__
#define __SCALAR_FIELD_H__

#include "Surface.h"
#include "mtxlib.h"
#include "multiDarrays.h"
#include "defines.h"
#include "Kernel.h"

class ScalarField : public Surface
{
public:
  ScalarField(const char *df_file_base,  
              float isovalue, int kernel=BSPLINE,
              float scale_x=1.0, 
              float scale_y=1.0, 
              float scale_z=1.0);
  ~ScalarField();

  enum { BSPLINE, CATMULLROM };

  bool computeSurfacePointParams(const vector_type &pos, 
                                 SurfacePointParams &params,
                                 bool computeGradient=true,
                                 bool computeHessian=false) const;  

 
private:
  array3D<float> _field, _sf;
  int _xdim, _ydim, _zdim;
  float _scale[3];
  Kernel *_kernel, *_kernelD;

  void readFile(const char *df_file_base, 
                int num_header_lines=2);
                
  void readNrrdFile( const char *df_file_base);

  //--------------------------------------------------------------------
  // a simple linear interpolation
  inline float interpolate(float t, float a, float b) const
  { return ( (1.0 - t)*a + t*b ); };

};


#endif // __SCALAR_FIELD_H__
