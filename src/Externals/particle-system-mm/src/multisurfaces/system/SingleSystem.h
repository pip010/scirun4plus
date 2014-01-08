//----------------------------------------------------------------------//
// FILE        : SingleSystem.h                                                
// DESCRIPTION : SingleSystem is just one system.
//----------------------------------------------------------------------//

#ifndef __SINGLE_SYSTEM_H__
#define __SINGLE_SYSTEM_H__

#include <system/systemExports.h>

#include <features/svector.h>
#include <features/mtxlib.h>
#include <defines.h>
#include <domain/Surface.h>
#include <particles/DynamicSurfacePoint.h>
#include <ParticleSystem.h>

namespace particle_sys 
{
  class System_SHARE SingleSystem 
  {
  public:
    SingleSystem(const char *param_file);
    ~SingleSystem();

    //--------------------------------------------------------------------
    // perform one step of optimization of the system
    void optimize();
    
    //--------------------------------------------------------------------
    // return true if the system is at an optimal state, false otherwise
    inline bool optimized() const
    { return _particle_sys->optimized(); };

    //--------------------------------------------------------------------
    // sets the input variables to the dimensions of the (assumed!) 
    //   rectangular domain
    inline void domain(vector_type &start, vector_type &end) const
    { _particle_sys->domain( start, end ); };

    //--------------------------------------------------------------------
    // returns a constant reference to the array of Points --> to be
    //   used for observing the Points, such as for rendering
    inline const custom_class::svector<DynamicSurfacePoint*>& 
      points() const { return _particle_sys->points(); };

    //--------------------------------------------------------------------
    // returns a constant reference to the a Surface --> to be
    //   used for observing the Surface, such as for rendering
    inline const Surface* surface() const 
    { return _particle_sys->surface(); };

    inline void reinitializePoints() const
    { _particle_sys->reinitializePoints(); };
    inline void offsetPoints(float offset_amt) const
    { _particle_sys->offsetPoints( offset_amt ); };

    inline void resetLambdas() const
    { _particle_sys->resetLambdas(); };

    inline void splitEveryPoint() const
    { _particle_sys->splitEveryPoint(); };

    //--------------------------------------------------------------------
    // functions for interacting with the system
    inline float angular_density() const
    { return _particle_sys->angular_density(); };
    inline void angular_density(float ad) const
    { _particle_sys->angular_density(ad); };
    inline float planar_seperation() const
    { return _particle_sys->planar_seperation(); };
    inline void planar_seperation(float ps) const
    { _particle_sys->planar_seperation(ps); };

    //--------------------------------------------------------------------
    // render the bounding box for the domain
    inline void domainRender() const
    { _particle_sys->domainRender(); };

    inline void surfaceRender() const
    { _particle_sys->surface()->render(); };
    inline void surfaceSelect(int x, int y) const
    { const_cast<Surface*>(_particle_sys->surface())->select(x,y); };

    inline void surfaceIsovalue(float iv) const
    { const_cast<Surface*>(_particle_sys->surface())->isovalue(iv); };


    //--------------------------------------------------------------------
    // print some stuff
    inline void print() const
    { _particle_sys->print(); };

    //--------------------------------------------------------------------
    // write out the point positions to a file
    inline void writePointFile(const char *filename) const
    { _particle_sys->writePointFile( filename ); };

    //--------------------------------------------------------------------
    // generate a Point Set Surface
    void generatePointSetSurface(int type);

  private:
    ParticleSystem *_particle_sys;
    //ofstream ***_out_file;
    //vector_type _original_pos[800], _previous_pos[800];

    char _mpu_weighting_function;

    void readParamFile(const char *param_file, Surface* &s, Energy* &e,
                       vec<3> &start, vec<3> &end, float &subd_w,
                       float &planar_sep, float &angular_den, 
                       bool &global_influence_radius,
                       float &influence_r, int &init_num_pts,
                       Optimization** &ops, int &num_ops, 
                       bool &biased_splitting_dying); 
    void readParamFile(const char *param_file, Surface* &s, Energy* &e,
                       vec<2> &start, vec<2> &end, float &subd_w,
                       float &planar_sep, float &angular_den, 
                       bool &global_influence_radius,
                       float &influence_r, int &init_num_pts,
                       Optimization** &ops, int &num_ops, 
                       bool &biased_splitting_dying); 
  };

} // namespace particle_sys

#endif // __SINGLE_SYSTEM_H__
