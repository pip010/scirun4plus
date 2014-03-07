#include <cstdlib>
#include <iostream>

using namespace std;
using namespace particle_sys;

//------------------------------------------------------------------------
// Function    : neighbors
// Description : check to see if these two Points are neighbors based
//               on the distance between them (must be less than the
//               radius of influence of the active Point) and their
//               normals (the angle between them must be less than 
//               90 degrees)
//------------------------------------------------------------------------
template <class T>
bool Neighborhood<T>::neighbors( const T *point, const T *neighbor ) const
{
  // first check the distance between these Points
  vector_type d_ij; point->distance( neighbor, d_ij );
  if ( d_ij.length() >= 1.0 )
    return false;
  
  //vector_type r_ij = point->position() - neighbor->position();
  //if ( (r_ij.length() >= point->influence_radius()) &&  
  //     (r_ij.length() >= neighbor->influence_radius()) )
  //  return false;

  // now, check if their normals are more than 90 degrees apart
  //if ( DotProduct(point->normal(),neighbor->normal()) <= 0.0 )
  //  return false;

  return true;
}

//------------------------------------------------------------------------
// Function    : operator << 
// Description : output function for the Neighborhood class
//------------------------------------------------------------------------
template <class T>
ostream& operator << ( ostream& outs, const Neighborhood<T>* source ) 
{  
  return outs; 
}





  

