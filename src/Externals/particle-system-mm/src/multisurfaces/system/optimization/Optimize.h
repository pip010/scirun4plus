//----------------------------------------------------------------------//
// FILE        : Optimize.h                                                
// DESCRIPTION : Optimize is the temporal component of the particle
//               system. This class is a container for the different
//               attributes of the System that we are striving to
//               optimize (such as distribution, density, etc.).
//               --> Currently, the optimization process works by first
//                   trying to optimize the first Optimization, then
//                   the second, etc. It is a nested loop.
//----------------------------------------------------------------------//

#ifndef __OPTIMIZE_H__
#define __OPTIMIZE_H__

#include <system/systemExports.h>

#include <features/svector.h>
#include <system/optimization/Optimization.h>

namespace particle_sys 
{
  class DynamicSurfacePoint;

  class Optimize : public Optimization
  {
  public:
    Optimize(Optimization **ops, int num_ops);
    virtual ~Optimize();

    //--------------------------------------------------------------------
    // initialize the Points by putting them onto the surface
    void init(custom_class::svector<DynamicSurfacePoint*> &points,
              int num_iterations=1);

    //--------------------------------------------------------------------
    // perform one step of optimization
    void optimize(custom_class::svector<DynamicSurfacePoint*> &points);

  private:
    Optimization **_optimizations;
    int            _num_ops;
    
  };

} // namespace particle_sys

#endif // __OPTIMIZE_H__
