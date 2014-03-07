//----------------------------------------------------------------------//
// FILE        : ScalarField.h                                                
// DESCRIPTION : 
//----------------------------------------------------------------------//


#ifndef __SCALAR_FIELD_H__
#define __SCALAR_FIELD_H__

#include <system/systemExports.h>

#include <system/domain/Surface.h>
#include <features/mtxlib.h>
#include <features/multiDarrays.h>
#include <system/defines.h>
#include <system/domain/Kernel.h>

#include <constants.h>

namespace particle_sys 
{
  class System_SHARE ScalarField : public Surface
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
                                   bool computeHessian) const;  

    float maxSF() const { return _max_sf; };
  
    DataCenter xcenter() const { return _centers[0]; }
    DataCenter ycenter() const { return _centers[1]; }
    DataCenter zcenter() const { return _centers[2]; }

  private:
#ifdef THREE_D
    array3D<float> _field, _sf;
    int _xdim, _ydim, _zdim;
    float _scale[3];
    DataCenter _centers[3];

#else
    array2D<float> _field;
    int _xdim, _ydim;
#endif
    Kernel *_kernel, *_kernelD;

    float _min_sf, _max_sf;

    void readFile(const char *df_file_base, 
                  int num_header_lines=2);
                  
    void readNrrdFile( const char *df_file_base);

    //--------------------------------------------------------------------
    // a simple linear interpolation
    inline float interpolate(float t, float a, float b) const
    { return ( (1.0 - t)*a + t*b ); };

  };

} // namespace particle_sys

#endif // __SCALAR_FIELD_H__
