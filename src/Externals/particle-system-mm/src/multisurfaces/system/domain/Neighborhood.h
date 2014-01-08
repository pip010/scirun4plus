//----------------------------------------------------------------------//
// FILE        : Neighborhood.h                                                
// DESCRIPTION : The Neighborhood is the interface between the Points 
//               and the SpatialBinningStructure. It returns a 
//               successive list of Points that are neighbors.
//
//               The neighborhood is templated on Point type -- thus,
//               it stores a neighborhood to objects that have a 
//               position() function.
//----------------------------------------------------------------------//

#ifndef __NEIGHBORHOOD_H__ 
#define __NEIGHBORHOOD_H__ 

#include <system/systemExports.h>

#include <features/mtxlib.h>
#include <system/defines.h>
#include <features/svector.h>

namespace particle_sys 
{ 
  template <class T>
  class Neighborhood
  {
  public:
    Neighborhood() { _storage_int = 0; };
    virtual ~Neighborhood(){};

    //--------------------------------------------------------------------
    // initially populates the spatial structure
    virtual void populateNeighborhood( const custom_class::svector<T*> 
                                       points) = 0;

    //--------------------------------------------------------------------
    // clear all the objects stored in the Neighborhood
    virtual void clearNeighborhood() = 0;

    //--------------------------------------------------------------------
    // does any sort of initialialize for subsequent calls to 
    //   nextNeighbor()
    virtual void determineNeighborhood(const T *point) = 0;

    //--------------------------------------------------------------------
    // returns a stream of Points that are indeed neighbors, then returns
    //   NULL when all neighbors have been returned
    virtual T* nextNeighbor(const T *point) = 0;

    //--------------------------------------------------------------------
    // add or remove a Point from the Neighborhood (to be called when
    //   Points are split or deleted)
    virtual void addToNeighborhood(T *point) = 0;
    virtual void removeFromNeighborhood(T *point) = 0;    

    //--------------------------------------------------------------------
    // a Point is changing locations, so update the spatial structure 
    virtual void movingLocation(T* const point, const vector_type &new_pos) = 0;

    //--------------------------------------------------------------------
    // this a function that allows an integer to be passed into the 
    //   Neighborhood structure if need be....
    inline void storageInt(int si) { _storage_int = si; };
    inline int storageInt() { return _storage_int; };

  protected:
    //--------------------------------------------------------------------
    // check if these two Points are neighbors based on the distance
    //   between them, and their normals
    virtual bool neighbors(const T *point, const T *neighbor) const;

    int _storage_int;
  };  

} // namespace particle_sys 



// print function
template <class T>
std::ostream& operator << (std::ostream& outs, 
                           const particle_sys::Neighborhood<T>* source);

#include "Neighborhood.T"

#endif // __NEIGHBORHOOD_H__ 
