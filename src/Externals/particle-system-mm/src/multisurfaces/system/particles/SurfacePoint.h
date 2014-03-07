//----------------------------------------------------------------------//
// FILE        : SurfacePoint.h                                                
// DESCRIPTION : The SurfacePoint class extends the Point class to
//               define Points that are associated with a Surface. This
//               class defines methods for querying the Surface for the
//               Surface parameters at this location. Also, defines the
//               Point normal to be the normalized Surface gradient. 
//               When a SurfacePoint is constructed, a point to the
//               associated Surface must be passed into the constructor.
//----------------------------------------------------------------------//

#ifndef __SURFACE_POINT_H__ 
#define __SURFACE_POINT_H__ 

#include <system/systemExports.h>

#include <features/mtxlib.h>
#include <system/defines.h>
#include <system/particles/Point.h>
#include <system/domain/SurfaceParameters.h>

namespace particle_sys 
{ 
  class System_SHARE SurfacePoint : public Point
  {
  public:
    SurfacePoint(Domain *d, ParticleSystem *sys):Point(d,sys) 
    { _max_sf = 1.0e6; _min_sf = 0.0; };
    SurfacePoint(Domain *d, ParticleSystem *sys, const vector_type &pos)
      :Point(d, sys, pos) { _max_sf = 1.0e6; _min_sf = 0.0; };
    SurfacePoint(const SurfacePoint &that):Point() { *this = that; };
    ~SurfacePoint() {};

    const SurfacePoint& operator = (const SurfacePoint &that);

    //--------------------------------------------------------------------
    // getter functions for the basic Surface parameters of this Point
    float F() const { return _surf_params._F; };
    const vector_type& Fx() const { return _surf_params._Fx; };
    const matrix_type& Fxx() const { return _surf_params._Fxx; };
    const matrix_type& curvature() const { return _surf_params._curvature; };
    float curvature_mag() const { return _surf_params._curvature_mag; };
    float kappa1() const { return _surf_params._kappa1; };
    float kappa2() const { return _surf_params._kappa2; };
    const vector_type& eigenvector1() const { return _surf_params._eigenvector1; };
    
    inline float sf() const 
    { return max(_min_sf, min(_surf_params._sf, _max_sf)); };
    inline void max_sf(float msf) { _max_sf = msf; };
    inline void min_sf(float msf) { _min_sf = msf; };

    const SurfacePointParams& surf_params() { return _surf_params; };

    //--------------------------------------------------------------------
    // compute the surface parameters at this Point location --> this
    //   will query the Surface for the F, Fx, and Fxx values
    //   
    // addition to handle some FE stuff -- will return false if for
    //   some reason this position does not fall within the defined
    //   domain (ie, the domain is not a rectangle!)
    bool updateSurfaceParameters();

    // a timing tester
    void onlyQuerySurface();

  protected:
    SurfacePointParams _surf_params;
    float _angular_density;
    float _max_sf, _min_sf;

    SurfacePoint():Point() {};
    void surfacePointAssignmentOp(const SurfacePoint &that); 

  private:

  };  

} // namespace particle_sys 

// print function
std::ostream& operator << (std::ostream& outs, 
                           const particle_sys::SurfacePoint& source);
#endif // __SURFACE_POINT_H__ 
