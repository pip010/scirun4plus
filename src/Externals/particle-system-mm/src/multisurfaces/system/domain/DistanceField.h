//----------------------------------------------------------------------//
// FILE        : DistanceField.h                                                
// DESCRIPTION : This surface reads in a distance field volume, and the 
//               first and second derivatives, and then uses the volumes
//               as the surface rep.
//----------------------------------------------------------------------//

#ifndef __DISTANCE_FIELD_H__
#define __DISTANCE_FIELD_H__

#include <system/systemExports.h>

#include <system/domain/Surface.h>
#include <features/mtxlib.h>
#include <features/multiDarrays.h>
#include <system/defines.h>

namespace particle_sys 
{
  class System_SHARE DistanceField : public Surface
  {
  public:
    DistanceField(int xdim, int ydim, const char *df_file_base);
    DistanceField(int xdim, int ydim, int zdim, const char *df_file_base);
    ~DistanceField() {};

    //--------------------------------------------------------------------
    // evaluate the surface parameters at this position
    bool computeSurfacePointParams(const vector_type &pos, 
                                   SurfacePointParams &params,
                                   bool computeHessian) const;

  private:
#ifdef TWO_D
    //int _xdim, _ydim;
    //array3D<float> _field;
    //array3D<float> _field_d; // first derivatives
    //array3D<float> _field_dd; // second derivatives
#endif

#ifdef THREE_D
    int _xdim, _ydim, _zdim;
    array3D<float> _field;
    array4D<float> _field_d; // first derivatives
    array4D<float> _field_dd; // second derivatives
#endif


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




#endif // __DISTANCE_FIELD_H__
