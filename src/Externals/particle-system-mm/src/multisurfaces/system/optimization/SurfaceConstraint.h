//----------------------------------------------------------------------//
// FILE        : SurfaceConstraint.h                                                
// DESCRIPTION : The SurfaceConstraint projects a DynamicSurfacePoint
//               back onto the Surface.
//----------------------------------------------------------------------//

#ifndef __SURFACE_CONSTRAINT_H__
#define __SURFACE_CONSTRAINT_H__

#include <system/systemExports.h>

#include <system/optimization/Optimization.h>

namespace particle_sys 
{
  class DynamicSurfacePoint;

  class System_SHARE SurfaceConstraint : public Optimization 
  {
  public:
    SurfaceConstraint(float F_threshold=1.0) {_F_threshold=F_threshold;};
    ~SurfaceConstraint(){};

    //--------------------------------------------------------------------
    // need to initialize the Points by first putting them onto the 
    //   surface
    void init(custom_class::svector<DynamicSurfacePoint*> &points,
              int num_iterations=1);

    //--------------------------------------------------------------------
    // the function is called to perform one iteration of optimization
    void optimize(custom_class::svector<DynamicSurfacePoint*> &points)
    { init(points); _optimized=true; };

    //--------------------------------------------------------------------
    // projects a DynamicSurfacePoint onto the Surface, and updates the
    //   Point's position
    void projectOntoSurface(DynamicSurfacePoint *point, float t = 0.000001);

  private:
    float _F_threshold;
    
  };

} // namespace particle_sys

#endif // __SURFACE_CONSTRAINT_H__
