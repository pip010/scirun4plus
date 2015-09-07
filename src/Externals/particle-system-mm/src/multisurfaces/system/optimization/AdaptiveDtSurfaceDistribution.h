//----------------------------------------------------------------------//
// FILE        : AdaptiveDtSurfaceDistribution.h                                                
// DESCRIPTION : This is an Optimization that distributes Points across
//               a Surface based on inter-Point repulsion forces, and 
//               integrates with a modified LM method adaptive dt.
//----------------------------------------------------------------------//

#ifndef __ADAPTIVE_DT_SURFACE_DISTRIBUTION_H__
#define __ADAPTIVE_DT_SURFACE_DISTRIBUTION_H__

#include <system/systemExports.h>

#include <features/svector.h>
#include <system/optimization/SurfaceConstraint.h>
#include <system/optimization/Optimization.h>


namespace particle_sys 
{
  class DynamicSurfacePoint;

  class AdaptiveDtSurfaceDistribution : public Optimization
  {
  public:
    AdaptiveDtSurfaceDistribution(float F_threshold=0.0001);
    ~AdaptiveDtSurfaceDistribution();

    //--------------------------------------------------------------------
    // need to initialize the Points by first putting them onto the 
    //   surface
    void init(custom_class::svector<DynamicSurfacePoint*> &points,
              int num_iterations=1);

    //--------------------------------------------------------------------
    // the function is called to perform one iteration of optimization
    void optimize(custom_class::svector<DynamicSurfacePoint*> &points);

// protected:
//    void parallelKernel(custom_class::svector<DynamicSurfacePoint*> &points, unsigned int start,  unsigned int  end);

  private:
    SurfaceConstraint *_constraint;
    float _prev_global_energy;
    int _num_iterations;

    float _F_threshold;
    std::vector<int> delP;

    //--------------------------------------------------------------------
    // check if the system is done moving
    bool doneMoving(custom_class::svector<DynamicSurfacePoint*> &points);
    
    void doMoving(custom_class::svector<DynamicSurfacePoint*> &points);
    void doEnergyCalc(custom_class::svector<DynamicSurfacePoint*> &points, unsigned int start, unsigned int end);
    
  };

} // namespace particle_sys

#endif // __ADAPTIVE_DT_SURFACE_DISTRIBUTION_H__
