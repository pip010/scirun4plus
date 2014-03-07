//----------------------------------------------------------------------//
// FILE        : Domain.h                                                
// DESCRIPTION : This is the spatial component of the system. Domain
//               is a container for things that are defined in space 
//               --> Surface, Neighborhood, Energy 
//----------------------------------------------------------------------//

#ifndef __DOMAIN_H__
#define __DOMAIN_H__

#include <system/systemExports.h>

#include <features/mtxlib.h>
#include <system/defines.h>
#include <features/svector.h>
#include <system/domain/SurfaceParameters.h>
#include <system/domain/Energy.h>
#include <system/domain/Neighborhood.h>
#include <system/domain/ApproximateStaticNeighborhood.h>
#include <system/domain/StaticNeighborhood.h>
#include <system/particles/DynamicSurfacePoint.h>
#include <system/domain/Surface.h>

namespace particle_sys 
{
  // here is where we will define the neighborhood type for the system
  typedef ApproximateStaticNeighborhood<DynamicSurfacePoint> neighborhood_type;
  //typedef StaticNeighborhood<DynamicSurfacePoint> neighborhood_type;


  class System_SHARE Domain
  {
  public:
    Domain(Energy *e, Surface *s, vector_type &start, vector_type &end, 
           float subd_w, int num_neighbors, float nid);
    Domain(Energy *e, Surface *s, vector_type &start, vector_type &end,  
           int num_neighbors, float nid);
    ~Domain();

    bool insideDomain(const vec<2> &pos); 
    bool insideDomain(const vec<3> &pos); 

    inline void domain(vector_type &start, vector_type &end) const
    { start = _start; end = _end; };
    inline void domain(vector_type &start, vector_type &end,
                       float &subd_w) const
    { start = _start; end = _end; subd_w = _subd_w; };

    inline float subd_width() const
    { return _subd_w; };

    inline void populateDomain(const 
                               custom_class::svector<DynamicSurfacePoint*>
                               &points)
    { _neighborhood->populateNeighborhood(points); };

//    void render() const;

    //--------------------------------------------------------------------
    // remove and add a Point into the Domain -- tell Neighborhood to 
    //   update it's data structures
    inline void removePoint(DynamicSurfacePoint *point)
    { _neighborhood->removeFromNeighborhood( point ); }

    inline void addPoint(DynamicSurfacePoint *point)
    { _neighborhood->addToNeighborhood( point ); }

    //--------------------------------------------------------------------
    // tell the Neighborhood that a Point is getting ready to move to a
    //   new position, so that Neighborhood can update it's data
    //   structures
    inline void movingToNewPosition(DynamicSurfacePoint* const point,
                                    const vector_type &new_pos)
    { _neighborhood->movingLocation(point,new_pos); };

    //--------------------------------------------------------------------
    // functions for interacting with a Surface
    inline bool computeSurfacePointParams(const vector_type &pos,
                                          SurfacePointParams &params,
                                          bool computeHessian) 
    { return _surface->computeSurfacePointParams(pos,params,computeHessian); };
    inline float computeCurvature(SurfacePointParams &params) 
    { return _surface->computeCurvature(params); };
//    inline void surfaceRender() const 
//    { _surface->render(); };
    inline void surfaceSelect(int x, int y) const 
    { _surface->select(x,y); };
    inline void surfaceIsovalue(float iv) const
    { _surface->isovalue(iv); };
    inline const Surface* surface() const
    { return _surface; };

    inline float surfaceIsovalue() const
    { return _surface->isovalue(); };

    inline void surface(Surface *surface)
    { delete _surface; _surface = surface; };
    inline void energy(Energy *energy)
    { delete _energy; _energy = energy; };

    //--------------------------------------------------------------------
    // functions for interacting with a Neighborhood
    inline void determineNeighborhood(const DynamicSurfacePoint *point) 
    { _neighborhood->determineNeighborhood(point); };
    inline DynamicSurfacePoint* nextNeighbor(const
                                             DynamicSurfacePoint *point) 
    { return _neighborhood->nextNeighbor(point); };

    inline Neighborhood<DynamicSurfacePoint>* neighborhood() const 
    { return _neighborhood; };
    inline void neighborhood(Neighborhood<DynamicSurfacePoint>* n,
                             float sub_d)
    { _neighborhood = n; _subd_w = sub_d; };

    //--------------------------------------------------------------------
    // functions for interacting with an Energy
    inline void initializeEnergy(const vector_type &d_ij)
    { _energy->initializeParameters(d_ij);};
    inline float solveEnergy() const 
    { return _energy->solveEnergy(); };
    inline vector_type solveForce() const 
    { return _energy->solveForce(); };
    inline matrix_type solveYank() const 
    { return _energy->solveYank(); };
    inline float idealEnergy() const
    { return (float)_num_neighbors*_energy->idealEnergy(); }

    inline float normalizedIdealDistance() const
    { return _normalized_ideal_distance; };


    inline int numNeighbors() const
    { return _num_neighbors; };

    

  private:
    Energy *_energy;
    Surface *_surface;
    Neighborhood<DynamicSurfacePoint> *_neighborhood;
    vector_type _start, _end;
    float _subd_w;

    int _num_neighbors;
    float _normalized_ideal_distance;

    //float _ideal_interparticle_energy;



  };

} // namespace particle_sys

// print function
System_SHARE std::ostream& operator << (std::ostream& outs,
                                        const particle_sys::Domain* source);

#endif // __DOMAIN_H__
