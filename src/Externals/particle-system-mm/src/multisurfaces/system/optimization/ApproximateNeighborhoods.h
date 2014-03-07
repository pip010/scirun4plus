//----------------------------------------------------------------------//
// FILE        : ApproximateNeighborhoods.h                                                
// DESCRIPTION : This is an Optimization that only updates the 
//               particles' neighborhood list every so many iterations.
//----------------------------------------------------------------------//

#ifndef __APPROXIMATE_NEIGHBORHOODS_H__
#define __APPROXIMATE_NEIGHBORHOODS_H__

#include <system/systemExports.h>

#include <features/svector.h>
#include <system/optimization/Optimization.h>

namespace particle_sys 
{
  class DynamicSurfacePoint;

  class ApproximateNeighborhoods : public Optimization
  {
  public:
    ApproximateNeighborhoods(int intervals=0);
    ~ApproximateNeighborhoods() {};

    //--------------------------------------------------------------------
    // the function is called to perform one iteration of optimization
    void optimize(custom_class::svector<DynamicSurfacePoint*> &points);

  protected:
    int _countdown, _num_iterations, _intervals;

  };
} // namespace particle_sys

#endif // __APPROXIMATE_NEIGHBORHOODS_H__
