#include <cstdlib>
#include <iostream>
#include <list>
#include <math.h>
#include <system/particles/DynamicSurfacePoint.h>
#include <system/domain/Domain.h>
#include <system/ParticleSystem.h>
#include <system/optimization/SurfacePopulationController.h>

#ifdef _WIN32
#pragma warning( disable : 4244 4018 )
#endif

using namespace particle_sys;
using namespace custom_class;
using namespace std;

//------------------------------------------------------------------------
// Function    : Constructor and Destructor
// Description : 
//------------------------------------------------------------------------
SurfacePopulationController::SurfacePopulationController( float threshold )
{
  _constraint = new SurfaceConstraint;
  _threshold = threshold;
}

SurfacePopulationController::~SurfacePopulationController()
{
  delete _constraint;
}

//------------------------------------------------------------------------
// Function    : deletePoints
// Description : delete some Points -- for now just do this randomly
//------------------------------------------------------------------------
void SurfacePopulationController::deletePoints( int num_delete,
                                 svector<DynamicSurfacePoint*> &points )
{
  // randomly pick Points to delete
  list<int> deleteList;
  unsigned i;

  for ( int j = 0; j < num_delete; j++ )
    deleteList.push_front( (int)((float)(points.size()-1) *
                                  rand()/(float)RAND_MAX) );

  // go through the list and remove the Points
  unsigned int index;
  while ( deleteList.size() )
  {
    index = deleteList.front();
    if ( index >= points.size() )
      index = points.size() - 1;
    deleteAPoint( index, points[index] );

    deleteList.pop_front();
  }

  // reset the lambda values of these Points
  for ( i = 0; i < points.size(); i++ )
    points[i]->resetLambda();
}

//------------------------------------------------------------------------
// Function    : splitPoints
// Description : split some Points -- for now just do this randomly
//------------------------------------------------------------------------
void SurfacePopulationController::splitPoints( int num_add,
                                    svector<DynamicSurfacePoint*> &points )
{
  // randomly split some Points
  list<int> splitList;
  unsigned i;

  num_add = min( (int)num_add, (int)points.size() );
  for ( int j = 0; j < num_add; j++ )
    splitList.push_front( (int)((float)(points.size()-1) *
                                 rand()/(float)RAND_MAX) );

  // go through the list and split the Points
  int index;
  while ( splitList.size() )
  {
    // pick a Point to split
    index = splitList.front();
    splitList.pop_front();

    // split the Point and add the new Point to the System
    splitAPoint( index, points[index] );    
  }

  // reset the lambda values of these Points
  for ( i = 0; i < points.size(); i++ )
    points[i]->resetLambda();

  cout << "    ... DONE\n";
}

//------------------------------------------------------------------------
// Function    : splitEveryPoint
// Description : 
//------------------------------------------------------------------------
void SurfacePopulationController::splitEveryPoint( svector<DynamicSurfacePoint*> 
                                           &points,
                                           int num_splits )
{
  int original_num_points = points.size();
  for ( int i = 0; i < original_num_points; i++ )
    for ( int j = 0; j < num_splits; j++ )
      splitAPoint( i, points[i] );
}

void SurfacePopulationController::splitEveryPointIntoFour( 
  svector<DynamicSurfacePoint*> &points )
{
  int original_num_points = points.size();
  for ( int i = 0; i < original_num_points; i++ )
    splitIntoFour( i, points[i] );
}

void SurfacePopulationController::splitEveryPointIntoFourWSF( 
  svector<DynamicSurfacePoint*> &points )
{
  int original_num_points = points.size();
  for ( int i = 0; i < original_num_points; i++ )
    splitIntoFourWSF( i, points[i] );
}

//------------------------------------------------------------------------
// Function    : deleteAPoint
// Description : delete the Point at the specified index
//------------------------------------------------------------------------
void SurfacePopulationController::deleteAPoint( int index, 
                                        DynamicSurfacePoint *point )
{
  (point->system())->removePoint( index );
}

//------------------------------------------------------------------------
// Function    : splitAPoint
// Description : split the Point at the specified index
//------------------------------------------------------------------------
void SurfacePopulationController::splitAPoint( int index,
                                       DynamicSurfacePoint *point )
{
  // split the Point 
  point->resetLambda(); 
  DynamicSurfacePoint *split_point = new DynamicSurfacePoint( *point );
  split_point->temporary( true );

  // store the orginial position of Point
  vector_type original_position = split_point->position();

  // position the new point between the original point, 
  // and its farthest neighbor.

  //! find farthest neighbor.
  double far_d2 = -1.0;
  const DynamicSurfacePoint *farthest = 0;
  // initialize the Neighborhood and get the first neighbor
  point->domain()->determineNeighborhood(point);
  // the current neighboring Point
  const DynamicSurfacePoint *cur = point->domain()->nextNeighbor(point);

  // iterate through all of the neighboring points
  vector_type d;
  while(cur) {
    // compute the scaled vector between these two points
    point->distance(cur, d);
    double d2 = d.lengthSqr();
    if (d2 > far_d2) {
      far_d2 = d2;
      farthest = cur;
    }
    cur = point->domain()->nextNeighbor(point);
  }


  if (! farthest) {
    delete split_point;
    return;
  }

  //! place the new point halfway between the this and the farthest.
  point->distance(farthest, d);
  d = d * 0.5;
  split_point->move(original_position + d);
  _constraint->projectOntoSurface(split_point);
  split_point->temporary( false );
  split_point->system()->addPoint( split_point );
}

//------------------------------------------------------------------------
// Function    : splitIntoFour
// Description : split a point into four points, and position them in a
//               square in the local tangent plane of the original point
//------------------------------------------------------------------------
void SurfacePopulationController::splitIntoFour( int index,
                                                 DynamicSurfacePoint *point )
{
  // get the original position of the center particle
  vector_type pos = point->position();

  // make an offset vector in the tangent plane
  vector_type tangent; 
  computeTangentPlaneVector( point->normal(), tangent );
  vector_type offset = 0.25*point->planar_seperation() * tangent;

#ifdef THREE_D
  // create the orthogonal offset vector
  vector_type ortho_offset = CrossProduct( point->normal(), tangent );
  ortho_offset *= 0.35*point->planar_seperation();

  // create three points
  DynamicSurfacePoint *split_point = new DynamicSurfacePoint( *point );
  split_point->temporary( true );
  if ( moveAndProject(split_point, pos-offset) )
  {
    split_point->temporary( false );
    split_point->system()->addPoint( split_point );
  }
  else
    delete split_point;

  split_point = new DynamicSurfacePoint( *point );
  split_point->temporary( true );
  if ( moveAndProject(split_point, pos+offset) )
  {
    split_point->temporary( false );
    split_point->system()->addPoint( split_point );
  }
  else
    delete split_point;

  split_point = new DynamicSurfacePoint( *point );
  split_point->temporary( true );
  if ( moveAndProject(split_point, pos+ortho_offset) )
  {
    split_point->temporary( false );
    split_point->system()->addPoint( split_point );
  }
  else
    delete split_point;

  // and move the original point
  if ( !moveAndProject(point, pos-ortho_offset) )
    point->system()->removePoint( index );

#endif

#ifdef TWO_D

#endif

}

//------------------------------------------------------------------------
// Function    : splitIntoFourWSF
// Description : split a point into four points, and position them in a
//               square in the local tangent plane of the original point
//------------------------------------------------------------------------
void SurfacePopulationController::splitIntoFourWSF( int index,
                                                    DynamicSurfacePoint *point )
{
  // get the original position of the center particle
  vector_type pos = point->position();

  // make an offset vector in the tangent plane
  vector_type tangent; 
  computeTangentPlaneVector( point->normal(), tangent );
  vector_type offset = 0.5*point->sf() * tangent;

#ifdef THREE_D
  // create the orthogonal offset vector
  vector_type ortho_offset = CrossProduct( point->normal(), tangent );
  ortho_offset *= 0.7*point->sf();

  // create three points
  DynamicSurfacePoint *split_point = new DynamicSurfacePoint( *point );
  split_point->temporary( true );
  if ( moveAndProject(split_point, pos-offset) )
  {
    split_point->temporary( false );
    split_point->system()->addPoint( split_point );
  }
  else
    delete split_point;

  split_point = new DynamicSurfacePoint( *point );
  split_point->temporary( true );
  if ( moveAndProject(split_point, pos+offset) )
  {
    split_point->temporary( false );
    split_point->system()->addPoint( split_point );
  }
  else
    delete split_point;

  split_point = new DynamicSurfacePoint( *point );
  split_point->temporary( true );
  if ( moveAndProject(split_point, pos+ortho_offset) )
  {
    split_point->temporary( false );
    split_point->system()->addPoint( split_point );
  }
  else
    delete split_point;

  // and move the original point
  if ( !moveAndProject(point, pos-ortho_offset) )
    point->system()->removePoint( index );

#endif

#ifdef TWO_D

#endif

}

//------------------------------------------------------------------------
// Function    : moveProjectAndAdd
// Description : 
//------------------------------------------------------------------------
void SurfacePopulationController::computeTangentPlaneVector(
  const vec<3> &normal, vec<3> &tangent )
{
  // get a vector that is orthonormal to the normal and the y-axis
  vec<3> y( 0.0, 1.0, 0.0 );

  // first, make sure we can actually take the cross product!
  if ( normal == y )
    tangent.set( 0.0, 0.0, 1.0 );
  else if ( normal == -y )
    tangent.set( 0.0, 0.0, -1.0 );

  else
  {
    tangent = CrossProduct( y, normal );
    tangent.normalize();
  }
}

void SurfacePopulationController::computeTangentPlaneVector(
  const vec<2> &normal, vec<2> &tangent )
{
  // first, check for the special cases:
  vec<2> x(1.0, 0.0), y(0.0, 1.0);
  if ( normal == y )
    tangent = x;
  else if ( normal == -y )
    tangent = -x;
  else if ( normal == x )
    tangent = -y;
  else if ( normal == -x )
    tangent = y;

  else
  {
    // otherwise, compute the k vector
    tangent = y - DotProduct( y, normal ) * normal;
    tangent.normalize();
  }
}

//------------------------------------------------------------------------
// Function    : moveAndProject
// Description : 
//------------------------------------------------------------------------
bool SurfacePopulationController::moveAndProject( DynamicSurfacePoint *point,
                                                  vector_type new_pos )
{
  // try to move the point to the new location
  point->move( new_pos );

  // if this new position is outside of the domain, kill the point
  if ( point->moved_outside() )
    return false;
  
  // otherwise, project the point back onto the surface
  _constraint->projectOntoSurface( point );
  
  // if this projected position is outside of the domain, kill the point
  if ( point->moved_outside() || (fabs(point->F()) > _threshold) )
    return false;

  return true;
}






