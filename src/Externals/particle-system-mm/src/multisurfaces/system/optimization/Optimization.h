//----------------------------------------------------------------------//
// FILE        : Optimization.h                                                
// DESCRIPTION : This is the interaface for the Optimizor. Defines how
//               all optimizations will interact with the Optimizor.
//----------------------------------------------------------------------//

#ifndef __OPTIMIZATION_H__
#define __OPTIMIZATION_H__

#include <system/systemExports.h>

#include <features/svector.h>

namespace particle_sys 
{
  class DynamicSurfacePoint;

  class System_SHARE Optimization
  {
  public:
    Optimization() { _optimized=true; };
    virtual ~Optimization(){};

    //--------------------------------------------------------------------
    // do any neccessary initializations to the System
    virtual 
      void init(custom_class::svector<DynamicSurfacePoint*> &points,
                int num_iterations=1) {};

    //--------------------------------------------------------------------
    // the function is called to perform one iteration of optimization
    virtual 
      void 
      optimize(custom_class::svector<DynamicSurfacePoint*> &points) = 0;

    //--------------------------------------------------------------------
    // returns whether or not this optimization is optimal
    inline bool optimized() const { return _optimized; };

    virtual void set_bool(bool t) { };

  protected:
    bool _optimized;
  };

} // namespace particle_sys

#endif // __OPTIMIZATION_H__
