//----------------------------------------------------------------------//
// FILE        : ApproximateStaticNeighborhood.h                                                
// DESCRIPTION : The ApproximateStaticNeighborhood implements the 
//               Neighborhood interface, and contains a static spatial 
//               binning structure. However, this neighborhood only 
//               updates the neighbor list of each particle every so
//               many iterations.
//----------------------------------------------------------------------//

#ifndef __APPROXIMATE_STATIC_NEIGHBORHOOD_H__ 
#define __APPROXIMATE_STATIC_NEIGHBORHOOD_H__ 

#include <system/systemExports.h>

#include <features/mtxlib.h>
#include <system/defines.h>
#include <features/svector.h>
#include <system/domain/Neighborhood.h>
#include <features/Bin.h>

namespace particle_sys 
{ 
  template <class T>
  class ApproximateStaticNeighborhood : public Neighborhood<T>
  {
  public:
    ApproximateStaticNeighborhood(float bin_w, const vector_type &start, 
                                 const vector_type &end);
    ~ApproximateStaticNeighborhood() { delete _bins; };

    //--------------------------------------------------------------------
    // pass a vector of Points, and put the Points into the appropriate
    //   Bins
    void populateNeighborhood(const custom_class::svector<T*> points);

    inline void clearNeighborhood() { _bins->clearBins(); };

    //--------------------------------------------------------------------
    // takes a DynamicSurfacePoint and creates the variables that will
    //   be used to return a stream of neighboring Points
    //   --> this must be called before nextNeighbor() is ever called for
    //       a specific Point query!
    void determineNeighborhood(const T *point);

    //--------------------------------------------------------------------
    // after determineNeighborhood() has been called for a Point, this
    //   function will return a stream of neighboring Points, and then
    //   returns NULL when the stream is done (or, if 
    //   determineNeighborhood() wasn't called first!) 
    T* nextNeighbor(const T *point);

    //--------------------------------------------------------------------
    // add or remove a Point from the Neighborhood (to be called when
    //   Points are split or deleted)
    void addToNeighborhood(T *point);
    void removeFromNeighborhood(T *point);    

    //--------------------------------------------------------------------
    // a Point is changing locations, so update the Bins 
    void movingLocation(T* const point, const vector_type &new_pos);

  private:
    BinningStructure<DIMENSION,T*> *_bins;
    const T *_calling_point;
    int _active_bin, _calling_point_bin, 
        _active_neighbor, _calling_point_bin_neighbor;
    
    custom_class::svector<T*> _points;

    int _num_neighbors, _neighbor_index;

    //--------------------------------------------------------------------
    // update the neighbor list of the object
    void updateNeighborhoodList(const T *point);

    T* nextNeighborForList(const T *point);


  };  

} // namespace particle_sys 

#include "ApproximateStaticNeighborhood.T"

#endif // __APPROXIMATE_STATIC_NEIGHBORHOOD_H__ 
