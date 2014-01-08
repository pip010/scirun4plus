#include <cstdlib>
#include <iostream>
#include <system/particles/Point.h>

using namespace particle_sys;
using namespace std;

//------------------------------------------------------------------------
// Function    : operator << 
// Description : output function for the Point class
//------------------------------------------------------------------------
ostream& operator << ( ostream& outs, const Point& source ) 
{
  outs << "<position,normal,radius> : <" << source.position() <<
    "), (" << source.normal() << "), " << source.radius() << ">\n"; 
  
  return outs; 
}











  

