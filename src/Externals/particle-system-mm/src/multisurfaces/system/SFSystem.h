//----------------------------------------------------------------------//
// FILE        : SFSystem.h                                                
// DESCRIPTION : SFSystem is just one system. 
//----------------------------------------------------------------------//

#ifndef __SF_SYSTEM_H__
#define __SF_SYSTEM_H__

#include <system/systemExports.h>

#include <features/svector.h>
#include <features/mtxlib.h>
#include <system/defines.h>
#include <system/domain/Surface.h>
#include <system/domain/Energy.h>
#include <system/domain/IOScalarField.h>
#include <system/particles/DynamicSurfacePoint.h>
#include <system/ParticleSystem.h>

namespace particle_sys 
{
  typedef struct 
  {
    ParticleSystem *ps;
    Surface *surf;
    float max_sf;
    int num_materials;
    int *materials;
    char mesh_filename[300];
    vec<3> color;
    bool init_with_ptcl, optimize_ptcl;
  } IntersectionPackage;

  class System_SHARE SFSystem 
  {
  public:
    SFSystem(const char *param_file);
    ~SFSystem();

    //--------------------------------------------------------------------
    // perform one step of optimization of the system
    inline void optimize()
    { if (_intersection[_current_intersection].ps->points().size() > 0) _intersection[_current_intersection].ps->optimize();}


    
    //--------------------------------------------------------------------
    // return true if the system is at an optimal state, false otherwise
    inline bool optimized()
    { return _intersection[_current_intersection].ps->optimized(); };

    //--------------------------------------------------------------------
    // sets the input variables to the dimensions of the (assumed!) 
    //   rectangular domain
    inline void domain(vector_type &start, vector_type &end) const
    { _intersection[0].ps->domain( start, end ); };
//    inline void domainRender() const 
//    { _intersection[0].ps->domainRender(); };

    //--------------------------------------------------------------------
    // returns a constant reference to the array of Points --> to be
    //   used for observing the Points, such as for rendering
    inline const custom_class::svector<DynamicSurfacePoint*>& 
      points(int i) const
    { return _intersection[i].ps->points(); };

    inline const vec<3>& color(int i) const
    { return _intersection[i].color; };

    inline int numberOfSystems() const 
    { return /*_num_intersections*/(_current_intersection+1); };
    inline int totalNumberOfSystems() const 
    { return _num_intersections; };

    //--------------------------------------------------------------------
    // returns a constant reference to the a Surface --> to be
    //   used for observing the Surface, such as for rendering
    inline const Surface* surface(int i) const 
    { return _intersection[i].surf; };

    //--------------------------------------------------------------------
    // print some stuff
    void print() const;

    //--------------------------------------------------------------------
    // write out the point positions to a file
    void writePointFile(int i) const;

    void freezeIntersection();

    inline void resetLambdas()
    { _intersection[_current_intersection].ps->resetLambdas(); };

    void cleanDoubleJunction();

//    void selectParticle(int x, int y);
//    void renderSelectedParticle();

  private:
    IntersectionPackage *_intersection;
    IOScalarField *_io_field;
    Surface **_material, **_IOmaterial;
    int _num_materials, _num_intersections, _current_intersection;
    float _initial_sf;

    int _num_qj, _num_tj, _num_dj;
    int *_qj, *_tj, *_dj;

    vector_type _start, _end;

    char _file_name[300];
    int _init_num_points;

    vector_type _selected_pt;

//    void processHits(GLint hits, GLuint buffer[], 
//                     float &min_z, int &min_index);

    void readIntersectionFile( const char *i_file );


    void readParamFile(const char *param_file, char *s, int &e,
                       float &max_sf, float &initial_sf,
                       int &init_num_pts,
                       int &num_surfaces, int &num_intersections,
                       char *i_file); 

    void initializePointsWithMesh(const char* basename, 
      custom_class::svector<DynamicSurfacePoint*> &points,
      ParticleSystem *ps, Domain *domain, float max_surface_sf,
      bool init_with_ptcl, int modulo=1);

    
  };

} // namespace particle_sys

#endif // __SF_SYSTEM_H__
