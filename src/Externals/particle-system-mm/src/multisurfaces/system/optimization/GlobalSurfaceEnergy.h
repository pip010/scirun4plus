//----------------------------------------------------------------------//
// FILE        : GlobalSurfaceEnergy.h                                                
// DESCRIPTION : This is an Optimization that strives to reach the 
//               ideal energy of the system. The ideal energy is based 
//               on the ideal energy of the Points, and is achieved
//               by splitting or killing Points across a surface.
//
//               There are two flavors of GlobalSurfaceEnergy -- the 
//               first is the Non Adaptive, which is used with a 
//               static binnning structure, or with global radius of 
//               influences for the Points. The other is Adaptive, 
//               which takes into account the problem of accumulating 
//               forces with adaptive radii and an adaptive binning
//               structure.
//               NOTE --> I haven't tackled the Adaptive class yet!!!!
//----------------------------------------------------------------------//

#ifndef __GLOBAL_SURFACE_ENERGY_H__
#define __GLOBAL_SURFACE_ENERGY_H__

#include <system/systemExports.h>

#include <features/svector.h>
#include <system/optimization/SurfaceConstraint.h>
#include <system/optimization/SurfacePopulationController.h>
#include <system/optimization/Optimization.h>

namespace particle_sys 
{
  class DynamicSurfacePoint;

  class GlobalSurfaceEnergy : public Optimization, 
    public SurfacePopulationController
  {
  public:
    GlobalSurfaceEnergy(bool biased_splitting_dying=true,
                        float threshold=100.0);
    ~GlobalSurfaceEnergy();

    //--------------------------------------------------------------------
    // the function is called to perform one iteration of optimization
    void optimize(custom_class::svector<DynamicSurfacePoint*> &points);

    void set_bool(bool t) { _first_time_through = t; };

  protected:
    SurfaceConstraint *_constraint;
    bool _biased_splitting_dying;
    bool _first_time_through;

    //--------------------------------------------------------------------
    // biased splitting and dying 
    virtual void 
      biasedOptimize(custom_class::svector<DynamicSurfacePoint*> &points)=0;

    //--------------------------------------------------------------------
    // random splitting and dying 
    virtual void 
      randomOptimize(custom_class::svector<DynamicSurfacePoint*> &points)=0;

  };

  /**********************************************************************/
  //                 Non Adaptive Updating Scheme                       //
  /**********************************************************************/
  class GlobalSurfaceEnergyNA : public GlobalSurfaceEnergy
  {
  public:
    GlobalSurfaceEnergyNA(bool biased_splitting_dying=true,
                          float threshold=100.0):
        GlobalSurfaceEnergy(biased_splitting_dying,threshold) {};
    ~GlobalSurfaceEnergyNA() {};

  private:

    //--------------------------------------------------------------------
    // biased splitting and dying 
    void biasedOptimize(custom_class::svector<DynamicSurfacePoint*> &points);

    //--------------------------------------------------------------------
    // random splitting and dying 
    void randomOptimize(custom_class::svector<DynamicSurfacePoint*> &points);

  };

} // namespace particle_sys

#endif // __GLOBAL_SURFACE_ENERGY_H__
