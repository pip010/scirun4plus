//----------------------------------------------------------------------//
// FILE        : SurfacePopulationController.h                                                
// DESCRIPTION : This is a simple class that has functions for splitting
//               and deleting particles. This class constrains points
//               that are split to the surface.
//----------------------------------------------------------------------//

#ifndef __SURFACE_POPULATION_CONTROLLER_H__
#define __SURFACE_POPULATION_CONTROLLER_H__

#include <system/systemExports.h>

#include <features/svector.h>
#include <features/mtxlib.h>
#include <system/defines.h>
#include <system/optimization/SurfaceConstraint.h>

namespace particle_sys 
{
  class DynamicSurfacePoint;

  class System_SHARE SurfacePopulationController
  {
  public:
    SurfacePopulationController(float threshold=100.0);
    ~SurfacePopulationController();

    //--------------------------------------------------------------------
    // for deleting a splitting one or many points
    void deletePoints(int num_delete,
                      custom_class::svector<DynamicSurfacePoint*> &points);
    void splitPoints( int num_split,
                      custom_class::svector<DynamicSurfacePoint*> &points);

    void deleteAPoint(int index, DynamicSurfacePoint *point);
    void splitAPoint( int index, DynamicSurfacePoint *point);

    //--------------------------------------------------------------------
    // split every point
    void splitEveryPoint(custom_class::svector<DynamicSurfacePoint*> &points,
                         int num_splits=1.0);
    void splitEveryPointIntoFour(custom_class::svector<DynamicSurfacePoint*> 
                                 &points);
    void splitEveryPointIntoFourWSF(custom_class::svector<DynamicSurfacePoint*> 
                                    &points);

    //--------------------------------------------------------------------
    // split the point into 4 points located in a square around the
    //   original position (or, into two if this is 2D)
    void splitIntoFour(int index, DynamicSurfacePoint* point);

    void splitIntoFourWSF(int index, DynamicSurfacePoint* point);

  private:
    SurfaceConstraint *_constraint;
    float _threshold;

    //--------------------------------------------------------------------
    // take the point and try to move it to the new location and project
    //    back onto the surface -- return true if the point is still
    //    inside the domain and should be added to the system
    bool moveAndProject(DynamicSurfacePoint *point,
                        vector_type new_pos);

    //--------------------------------------------------------------------
    // compute a tangent plane vector
    void computeTangentPlaneVector(const vec<3> &normal, 
                                   vec<3> &tangent);
    void computeTangentPlaneVector(const vec<2> &normal, 
                                   vec<2> &tangent);

  };

} // namespace particle_sys

#endif // __SURFACE_POPULATION_CONTROLLER_H__
