#include <cstdlib>
#include <iostream>
#include <LFS.h>

#include <stdio.h>

using namespace std;
using namespace custom_class;

//------------------------------------------------------------------------
// Function    : Constructor and Destructor
// Description : 
//------------------------------------------------------------------------
LocalFeatureSize::LocalFeatureSize( const char *filename,
                                    const vec<3> &s, const vec<3> &e,
                                    int xdim, int ydim, int zdim ) 
{  
  _xdim = xdim; _ydim = ydim; _zdim = zdim;
  _start = s;
  _end = e;

  // create the octree and set the domain dimensions
  _octree = new OctreeNode<Point>;
  _octree->setDimensions( _start, _end );

  // read in the medial axis points and populate the octree
  readInMAPoints( filename );
}

LocalFeatureSize::~LocalFeatureSize()
{
  delete _octree;
  
  for ( int i = 0; i < _medial_axis_pts.size(); i++ ) {
    cout << "deleting at idx: " << i << endl;
    delete _medial_axis_pts[i];
  }
}

//------------------------------------------------------------------------
// Function    : lfs()
// Description : return the local feature size at a position -- ie.
//   distance to the closest medial axis point
//------------------------------------------------------------------------
float LocalFeatureSize::lfs( const vector_type &pos ) const
{
  Point pt;
  pt.x = pos;

  return _octree->distanceToClosestPoint( &pt );
}

//------------------------------------------------------------------------
// Function    : addMedialAxisPoint()
// Description : add a point to the medial axis list, and send the
//   pointer to the point into the octree
//------------------------------------------------------------------------
void LocalFeatureSize::addMedialAxisPoint( const vector_type &pos,
                                           const vector_type &n )
{
  Point *pt = new Point;
  pt->x = pos;
  pt->n = n;
  pt->got_normal(true);
  _medial_axis_pts.push_back( pt );

  Point *opt = new Point;
  opt->x = pos;
  opt->n = n;
  opt->got_normal(true);

  _octree->addPoint( opt );
}

//------------------------------------------------------------------------
// Function    : readInMAPoints()
// Description : 
//------------------------------------------------------------------------
void LocalFeatureSize::readInMAPoints( const char *filename )
{
   FILE *in_ma_ptcl  = fopen( filename, "r" );
  if ( !in_ma_ptcl  )
  {
    cout << "LocalFeatureSize:: PROBLEM, no medial axis!!!" 
         << " Couldn't open: " << filename << endl;
    exit( 1 );
  }
  
  cout << "Reading the set of medial points..." << endl;

  // read in the particle file
  int num_ma_pts;
  fscanf( in_ma_ptcl, "%i\n", &num_ma_pts );

  cout << "     Number of medial axis points = " << num_ma_pts << endl;

  vector_type pos, n;
  float r;
  for ( int i = 0; i < num_ma_pts; i++ )
  {
    fscanf( in_ma_ptcl, "%f %f %f %f %f %f %f\n", 
      &pos[0], &pos[1], &pos[2], &n[0], &n[1], &n[2], &r );

    if ( !(i%1) )
    {
      //cout << i << endl;
      addMedialAxisPoint( pos, n );
    }
  }
  cout << "added " << _medial_axis_pts.size() << " medial axis points."
       << endl;

  cout << "    ... DONE." << endl;
  fclose( in_ma_ptcl );
}

