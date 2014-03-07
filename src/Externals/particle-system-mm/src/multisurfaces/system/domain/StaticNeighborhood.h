//----------------------------------------------------------------------//
// FILE        : StaticNeighborhood.h                                                
// DESCRIPTION : The StaticNeighborhood implements the Neighborhood
//               interface, and contains a static spatial binning
//               structure.
//----------------------------------------------------------------------//

#ifndef __STATIC_NEIGHBORHOOD_H__ 
#define __STATIC_NEIGHBORHOOD_H__ 

#include <features/mtxlib.h>
#include <system/defines.h>
#include <features/svector.h>
#include <system/domain/Neighborhood.h>
#include <features/Bin.h>

namespace particle_sys 
{ 
  template <class T>
  class StaticNeighborhood : public Neighborhood<T>
  {
  public:
    StaticNeighborhood(float bin_w, const vector_type &start, 
                       const vector_type &end);
    ~StaticNeighborhood() { delete _bins; };

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
    //   determineNeighborhood() wasn't called first!) --> does a check
    //   that determineNeighborhood() was called by comparing the
    //   'active' Point with the one that is calling this function!
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

  };  

} // namespace particle_sys 

#include "StaticNeighborhood.T"

#endif // __STATIC_NEIGHBORHOOD_H__ 
