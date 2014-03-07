//----------------------------------------------------------------------//
// FILE        : DistanceFieldwSF.h                                                
// DESCRIPTION : This surface reads in a distance field volume and the 
//               first derivatives, along with a sizing field.
//----------------------------------------------------------------------//

#ifndef __DISTANCE_FIELD_W_SF_H__
#define __DISTANCE_FIELD_W_SF_H__

#include <system/systemExports.h>

#include <system/domain/Surface.h>
#include <features/mtxlib.h>
#include <features/multiDarrays.h>
#include <system/defines.h>

namespace particle_sys 
{
  class System_SHARE DistanceFieldwSF : public Surface
  {
  public:
    DistanceFieldwSF(const char *df_file_base, float &max_sf, 
                     int &xdim, int &ydim, int &zdim);
    ~DistanceFieldwSF() {};

    //--------------------------------------------------------------------
    // evaluate the surface parameters at this position
    bool computeSurfacePointParams(const vector_type &pos, 
                                   SurfacePointParams &params,
                                   bool computeHessian) const;

  private:
    int _xdim, _ydim, _zdim;
    array3D<float> _field, _sf;
    array4D<float> _field_d; // first derivatives

    float _min_sf, _max_sf;


    //--------------------------------------------------------------------
    // read in either a distance field and create the volumes
    //    ---> the num_header_lines specifies the numer of lines in the 
    //         files before the data starts
    void readDistanceFieldFile(const char *df_file_base, 
                               int num_header_lines=2);

    //--------------------------------------------------------------------
    // a simple linear interpolation
    inline float interpolate(float t, float a, float b) const
    { return ( (1.0 - t)*a + t*b ); };
  };

 

} // namespace particle_sys




#endif // __DISTANCE_FIELD_W_SF_H__
