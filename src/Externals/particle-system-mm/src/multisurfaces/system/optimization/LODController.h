//----------------------------------------------------------------------//
// FILE        : LODController.h                                                
// DESCRIPTION : This is an Optimization that creates a course to fine
//               set of particles. The op starts by distributing a set
//               of large particles, and gradually decreases the size of
//               the particles.
//----------------------------------------------------------------------//

#ifndef __LOD_DISTRIBUTOR_H__
#define __LOD_DISTRIBUTOR_H__

#include <system/systemExports.h>

#include <features/svector.h>
#include <system/optimization/SurfaceConstraint.h>
#include <system/optimization/SurfacePopulationController.h>
#include <system/optimization/Optimization.h>

namespace particle_sys 
{
  class DynamicSurfacePoint;
  class Domain;
  class ParticleSystem;

  class LODController : public Optimization, 
    public SurfacePopulationController
  {
  public:
    LODController(Domain *d, ParticleSystem *ps,
                  custom_class::svector<DynamicSurfacePoint*> &p,
                  float initial_sf, float max_surface_sf,
                  int num_additional_ps=0, 
                  ParticleSystem **additional_ps=NULL,
                  int init_num_pts=100,
                  float threshold=100.0,
                  const char* basename=NULL,
                  int modulo=1);
    ~LODController();

    //--------------------------------------------------------------------
    // the function is called to perform one iteration of optimization
    void optimize(custom_class::svector<DynamicSurfacePoint*> &points);

  protected:
    SurfaceConstraint *_constraint;
    Domain *_domain;
    float _min_sf, _max_surface_sf;
    int _modulo;

    int _num_additional_ps;
    ParticleSystem **_additional_ps;

    //--------------------------------------------------------------------
    // initialize the particle system with a course neighborhood and big
    //   particles
    void buildInitialParticlesAndNeighborhood(
      ParticleSystem *ps,
      custom_class::svector<DynamicSurfacePoint*> &points,
      int init_num_pts, const char* basename);

    //--------------------------------------------------------------------
    // do a pass of level of detail reduction
    void reduceLOD(custom_class::svector<DynamicSurfacePoint*> &points);

    void initializePointsWithMesh(const char* basename, 
      custom_class::svector<DynamicSurfacePoint*> &points,
      ParticleSystem *ps);

  };
} // namespace particle_sys

#endif // __LOD_DISTRIBUTOR_H__
