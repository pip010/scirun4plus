#include <cstdlib>
#include <iostream>
#include <system/domain/StaticNeighborhood.h>
#include <system/domain/ApproximateStaticNeighborhood.h>
#include <system/domain/Domain.h>

using namespace particle_sys;
using namespace std;

//------------------------------------------------------------------------
// Function    : Constructor and Destructor
// Description : 
//------------------------------------------------------------------------
Domain::Domain( Energy *e, Surface *s, vector_type &start, 
                vector_type &end, float subd_w, int num_neighbors,
                float nid )
{
  _energy = e;
  _surface = s;

  _start = start; _end = end; _subd_w = subd_w;
  _neighborhood = new neighborhood_type( _subd_w, _start, _end );

  _num_neighbors = num_neighbors;
  _normalized_ideal_distance = nid;
}

Domain::Domain( Energy *e, Surface *s, vector_type &start, 
                vector_type &end, int num_neighbors,
                float nid )
{
  _energy = e;
  _surface = s;

  _start = start; _end = end; 
  _neighborhood = NULL;

  _num_neighbors = num_neighbors;
  _normalized_ideal_distance = nid;
}

Domain::~Domain()
{
  delete _neighborhood;
  delete _energy;
  delete _surface;
}

//------------------------------------------------------------------------
// Function    : insideDomain()
// Description : determine whether a position is inside or outside of the
//               domain 
//               --> overloaded this function to account for 2D or 3D
//------------------------------------------------------------------------
bool Domain::insideDomain( const vec<2> &pos )
{
  // need to do it this way to make sure that nan's don't pass!!
  if ( (pos[0] <= _end[0])   && (pos[1] <= _end[1]) &&
       (pos[0] >= _start[0]) && (pos[1] >= _start[1]) ) 
    return true;
  return false;
}

bool Domain::insideDomain( const vec<3> &pos )
{
  // need to do it this way to make sure that nan's don't pass!!
  if ( (pos[0] <= _end[0])   && (pos[1] <= _end[1])   && (pos[2] <= _end[2]) &&
       (pos[0] >= _start[0]) && (pos[1] >= _start[1]) && (pos[2] >= _start[2]) ) 
    return true;
  return false;
}

#if 0
//------------------------------------------------------------------------
// Function    : render()
// Description : render the bounding box for the domain
//------------------------------------------------------------------------
void Domain::render() const
{
  glColor3f( 0.0, 0.0, 0.0 );
  glLineWidth( 1.0 );
  glBegin( GL_LINE_STRIP );
  {
    glColor3f( 1.0, 0.0, 0.0 );
    glVertex3f( _start(0), _start(1), _start(2) );
    glVertex3f( _end(0),   _start(1), _start(2) );
    glColor3f( 0.0, 0.0, 0.0 );
    glVertex3f( _end(0),   _end(1),   _start(2) );
    glVertex3f( _start(0), _end(1),   _start(2) );
    glColor3f( 0.0, 1.0, 0.0 );
    glVertex3f( _start(0), _start(1), _start(2) );
  } glEnd();
  glBegin( GL_LINE_LOOP );
  {
    glColor3f( 0.0, 0.0, 0.0 );
    glVertex3f( _start(0), _start(1), _end(2) );
    glVertex3f( _end(0),   _start(1), _end(2) );
    glVertex3f( _end(0),   _end(1),   _end(2) );
    glVertex3f( _start(0), _end(1),   _end(2) );
  } glEnd();
  glBegin( GL_LINES );
  {
    glColor3f( 0.0, 0.0, 1.0 );
    glVertex3f( _start(0), _start(1), _start(2) );
    glVertex3f( _start(0), _start(1), _end(2)   );

    glColor3f( 0.0, 0.0, 0.0 );
    glVertex3f( _end(0),   _start(1), _start(2) );
    glVertex3f( _end(0),   _start(1), _end(2)   );

    glVertex3f( _end(0),   _end(1),   _start(2) );
    glVertex3f( _end(0),   _end(1),   _end(2)   );

    glVertex3f( _start(0), _end(1),   _start(2) );
    glVertex3f( _start(0), _end(1),   _end(2)   );
  }glEnd();
}
#endif

//------------------------------------------------------------------------
// Function    : operator << 
// Description : output function for the Domain class
//------------------------------------------------------------------------
ostream& operator << ( ostream& outs, const Domain* source ) 
{
  outs << "Neighborhood : " << endl;
  outs << "               " << source->neighborhood() << endl;
  
  return outs; 
}





  

