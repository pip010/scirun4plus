//----------------------------------------------------------------------//
// FILE        : ParticleSystem.h                                                
// DESCRIPTION : ParticleSystem is the container for the dynamic particle 
//               system. A ParticleSystem is composed of objects 
//               (DynamicSurfacePoints), spatial (Domain), and temporal
//               (Optimize). ParticleSystem instatinates and builds the 
//               complete package.
//----------------------------------------------------------------------//

#ifndef __PARTICLE_SYSTEM_H__
#define __PARTICLE_SYSTEM_H__

#include <system/systemExports.h>

#include <features/svector.h>
#include <features/mtxlib.h>
#include <system/defines.h>
#include <system/domain/Surface.h>
#include <system/domain/Domain.h>
#include <system/optimization/Optimize.h>
#include <system/particles/DynamicSurfacePoint.h>

namespace particle_sys 
{
  class System_SHARE ParticleSystem
  {
  public:
    ParticleSystem();
    ~ParticleSystem();

    //--------------------------------------------------------------------
    // set the variables for this object
    void init(Domain *d, Optimize *o,
      custom_class::svector<DynamicSurfacePoint*> &p);
    void init(Domain *d, Optimize *o);
    void init(custom_class::svector<DynamicSurfacePoint*> &p);

    //--------------------------------------------------------------------
    // perform one step of optimization of the system
    inline void optimize() { _optimizer->optimize(_points); };
    
    //--------------------------------------------------------------------
    // return true if the system is at an optimal state, false otherwise
    inline bool optimized(){ return _optimizer->optimized(); };

    //--------------------------------------------------------------------
    // sets the input variables to the dimensions of the (assumed!) 
    //   rectangular domain
    inline void domain(vector_type &start, vector_type &end) const
    { _domain->domain( start, end ); };
    inline void domain(vector_type &start, vector_type &end,
                       float &subd_w) const
    { _domain->domain( start, end, subd_w ); };

    inline Domain* domain() const
    { return _domain; };

    //--------------------------------------------------------------------
    // returns a constant reference to the array of Points --> to be
    //   used for observing the Points, such as for rendering
    inline custom_class::svector<DynamicSurfacePoint*>& 
      points() { return _points; };
    inline DynamicSurfacePoint* point(int i) const
    { return _points[i]; };

    //--------------------------------------------------------------------
    // returns a constant reference to the a Surface --> to be
    //   used for observing the Surface, such as for rendering
    inline const Surface* surface() const { return _domain->surface(); };

    //--------------------------------------------------------------------
    // functions for sending down new components (the old ones are deleted)
    inline void surface(Surface *s) const { _domain->surface(s); };
    inline void energy(Energy *e) const { _domain->energy(e); };
    inline void optimizer(Optimize *o) 
    { delete _optimizer; _optimizer = o; };

    //--------------------------------------------------------------------
    // rendering and selection functions for interacting with the 
    //   Surface
//    inline void surfaceRender() const 
//    { _domain->surfaceRender(); };
    inline void surfaceSelect(int x, int y) const 
    { _domain->surfaceSelect(x,y); };

    void surfaceIsovalue(float iso);
    inline float  surfaceIsovalue() const
    { return _domain->surfaceIsovalue(); };

    //--------------------------------------------------------------------
    // delete and add a Point from/to the System
    void removePoint(int i);
    void addPoint(DynamicSurfacePoint *point);

    void reinitializePoints();
    void resetLambdas();
    void offsetPoints(float offset_amt);

    void splitEveryPoint();

    //--------------------------------------------------------------------
    // go through the points and get rid of ones that have a suface
    //   scalar value below some threshold
    void cleanUpSystem(float threshold);

    //--------------------------------------------------------------------
    // print some stuff
    void print() const;

    //--------------------------------------------------------------------
    // functions for interacting with the system
    inline float angular_density() const
    { if ( !_points.empty() ) return _points[0]->angular_density();
      else return 0.0; };
    inline void angular_density(float ad) 
    { for (unsigned i=0; i<_points.size(); i++) 
        _points[i]->angular_density(ad); };
    inline float planar_seperation() const
    { if ( !_points.empty() ) return _points[0]->planar_seperation();
      else return 1.0; };
    void planar_seperation(float ps);

    //--------------------------------------------------------------------
    // render the bounding box for the domain
//    inline void domainRender() const
//    { _domain->render(); };

    //--------------------------------------------------------------------
    // write out the point positions to a file
    void writePointFile(const char *filename) const;

    void writeEpsilonSampleFile(const char *filename) const;

    //--------------------------------------------------------------------
    // number of particles
    inline int numParticles() const
    { return _points.size(); };

    //--------------------------------------------------------------------
    // some timing functions
    void timeEnergyForceComputations();
    void timeSurfaceUpdates();

  private:
    Domain *_domain;
    Optimize *_optimizer;
    custom_class::svector<DynamicSurfacePoint*> _points;
  };

} // namespace particle_sys

#endif // __PARTICLE_SYSTEM_H__
