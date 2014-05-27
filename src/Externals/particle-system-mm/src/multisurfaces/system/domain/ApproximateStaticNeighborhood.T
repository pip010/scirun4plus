#include <cstdlib>
#include <iostream>




using namespace custom_class;
using namespace std;



//------------------------------------------------------------------------
// Function    : Constructor
// Description : 
//------------------------------------------------------------------------
template <class T>
ApproximateStaticNeighborhood<T>::ApproximateStaticNeighborhood( float bin_w, 
                                                  const vector_type &start, 
                                                  const vector_type &end ) : 
                                             Neighborhood<T>()
{
  _bins = new BinningStructure<DIMENSION,T*>
  ( bin_w, start, end );
  _bins->buildBinNeighbors();

  _calling_point = NULL;
  _calling_point_bin = 0;
  _active_bin = 0;
  _active_neighbor = _calling_point_bin_neighbor = -1;
  
  Neighborhood<T>::_storage_int = 0;
}

//------------------------------------------------------------------------
// Function    : populateNeighborhood
// Description : go through the list of Points and add them to the 
//               appropriate Bin
//------------------------------------------------------------------------
template <class T>
void ApproximateStaticNeighborhood<T>::populateNeighborhood(const svector<T*> 
                                                           points)
{
  for ( unsigned int i = 0; i < points.size(); i++ )
    _bins->addElementToBin( _bins->whichBin(points[i]->position()), 
                            points[i] );

}

//------------------------------------------------------------------------
// Function    : removeFrom/addToNeighborhood
// Description : remove or add the Point from or to the neighborhood
//               (to be called when Points are split or deleted)
//------------------------------------------------------------------------
template <class T>
void ApproximateStaticNeighborhood<T>::removeFromNeighborhood( T *point )
{
lock_guard<mutex> guardian(lock);
  _bins->removeElementFromBin( _bins->whichBin(point->position()),
                               point );
                               
  // set flag so that neighborhoods are updated!!
  Neighborhood<T>::_storage_int = -1;
}

template <class T>
void ApproximateStaticNeighborhood<T>::addToNeighborhood( T *point )
{
  _bins->addElementToBin( _bins->whichBin(point->position()), 
                          point );
                          
  // set flag so that neighborhoods are updated!!
  Neighborhood<T>::_storage_int = -1;
}

//------------------------------------------------------------------------
// Function    : movingLocation
// Description : a Point is moving to a new location, so we need
//               to remove the pointer to the Point from its old bin,
//               and add the pointer to the new bin
//------------------------------------------------------------------------
template <class T>
void ApproximateStaticNeighborhood<T>::movingLocation( T* const point,
                                                      const vector_type &new_pos )
{

  //lock_guard<mutex> guardian(lock);

  int old_bin = _bins->whichBin( point->position() );
  int new_bin = _bins->whichBin( new_pos );

  if ( new_bin != old_bin )
  {
	//cout << " BM ";
    _bins->removeElementFromBin( old_bin, point );
    _bins->addElementToBin( new_bin, point );
  }
}

//------------------------------------------------------------------------
// Function    : determineNeighborhood
// Description : 
//------------------------------------------------------------------------
template <class T>
void ApproximateStaticNeighborhood<T>::determineNeighborhood( const T *point )
{
  // reset the neighbor index and the active point
  _calling_point = point;
  _calling_point_bin = _active_bin = 
    _bins->whichBin( _calling_point->position() );
  _active_neighbor = _calling_point_bin_neighbor = -1;
  
  // update the neighborhood lists if we have reached the end of the 
  //   countdown
  if ( Neighborhood<T>::_storage_int <= 0 )
    updateNeighborhoodList( point );
    
  _neighbor_index = 0;
  _num_neighbors = const_cast<T*>(point)->neighbors()->size();
}

//------------------------------------------------------------------------
// Function    : updateNeighborhoodList
// Description : 
//------------------------------------------------------------------------
template <class T>
void ApproximateStaticNeighborhood<T>::updateNeighborhoodList( const T *point )
{

  // get the location of the neighbor list
  svector<T*> *list = const_cast<T*>(point)->neighbors();
  list->clear();
  
  T* neighbor = nextNeighborForList( point );
  while ( neighbor )
  {
    list->push_back( neighbor );
    neighbor = nextNeighborForList( point );
  }
}

//------------------------------------------------------------------------
// Function    : nextNeighborForList
// Description : this function will iterate through all the potential 
//               neighbors, return those that qualify as a valid
//               neighbor (ie. are within the active Point's radius of
//               influence, and their normals are less than 90 degrees
//               apart)
//               --> return NULL when the end of the neighbor list is
//                   reached, or if determineNeighborhood() was not
//                   called before calling this function (such that 
//                   the active_point doesn't match the Point that called
//                   this function)
//------------------------------------------------------------------------
template <class T>
T* ApproximateStaticNeighborhood<T>::nextNeighborForList( const T *point )
{
  // make sure that determineNeighborhood was called for this Point
  if ( point != _calling_point )
  {
    cout << "ApproximateStaticNeighborhood::nextNeighborForList() : Need to call " <<
      "determineNeighborhood() first!\n";
    return NULL;
  }

  // get the next potential neighbor
  T *next_neighbor;

  // go to the next Point in the active bin
  ++_active_neighbor;

  // first, check if the active bin still has Points to return
  if ( _active_neighbor == 
       _bins->getBinsNumElements( _active_bin ) )
  {
    // go to the next neighboring bin
    ++_calling_point_bin_neighbor;

    // but first, need to make sure we are not at the end of the
    //   neighboring bin list!
    if ( _calling_point_bin_neighbor == 
         _bins->getBinsNumNeighbors( _calling_point_bin ) )
    {
      _calling_point = NULL;
      return NULL;
    }

    // get the next bin
    _active_bin = 
      _bins->getBinsNeighborAt( _calling_point_bin, 
                                _calling_point_bin_neighbor );
    _active_neighbor = 0;

    // check if this neighboring bin is empty
    while ( _bins->isBinEmpty( _active_bin ) )
    {
      // go to the next neighboring bin
      ++_calling_point_bin_neighbor;

      // but first, need to make sure we are not at the end of the
      //   neighboring bin list!
      if ( _calling_point_bin_neighbor == 
          _bins->getBinsNumNeighbors( _calling_point_bin ) )
      {
        _calling_point = NULL;
        return NULL;
      }

      // get the next bin
      _active_bin = 
        _bins->getBinsNeighborAt( _calling_point_bin, 
                                  _calling_point_bin_neighbor );
    }
  }

  // now, _active_bin is the index of a non-empty bin, and 
  //   _active_neighbor is the index of the next neighbor in the 
  //   bin

  // get the next neighbor
  next_neighbor = _bins->getBinsElementAt( _active_bin, 
                                           _active_neighbor );
  
  // check if this neighbor is actually the active Point
  if ( next_neighbor == point )
    return nextNeighborForList( point );

  // now, check if this neighbor is really a neighbor, based on the 
  //   the distance between them and their normals
  if ( this->neighbors( point, next_neighbor ) )
    return next_neighbor;
  else
    return nextNeighborForList( point );
}

//------------------------------------------------------------------------
// Function    : nextNeighbor
// Description : 
//------------------------------------------------------------------------
template <class T>
T* ApproximateStaticNeighborhood<T>::nextNeighbor( const T *point )
{

  //cout << _neighbor_index << "  " << _num_neighbors << endl;
  if ( _neighbor_index == _num_neighbors )
    return NULL;
    
  _neighbor_index++;
  
  if( !  const_cast<T*>(point)->neighbors()->isValid(_neighbor_index) )
	return NULL;
      
  return ( const_cast<T*>(point)->neighbors()->at( _neighbor_index ) );
}
