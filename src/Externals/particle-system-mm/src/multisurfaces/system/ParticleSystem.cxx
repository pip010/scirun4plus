#include <cstdlib>
#include <iostream>
#include <fstream>
#include <time.h>
#include <system/optimization/GlobalSurfaceEnergy.h>
#include <system/ParticleSystem.h>

#ifdef _WIN32
#pragma warning ( disable: 4305 4018 )
#endif

using namespace particle_sys;
using namespace custom_class;
using namespace std;

//------------------------------------------------------------------------
// Function    : Constructor and Destructor
// Description : 
//------------------------------------------------------------------------
ParticleSystem::ParticleSystem()
{
  _domain = NULL;
  _optimizer = NULL;
  _points.clear();
}

ParticleSystem::~ParticleSystem()
{
  delete _domain;
  delete _optimizer;
  for ( unsigned i = 0; i < _points.size(); i++ )
    delete _points[i];
}

//------------------------------------------------------------------------
// Function    : init
// Description : initialize the variables
//------------------------------------------------------------------------
void ParticleSystem::init( Domain *d, Optimize *o, 
                           svector<DynamicSurfacePoint*> &p )
{
  _domain = d;
  _optimizer = o;
  _points = p;

  // initialize the optimizations
  _optimizer->init( _points, 5 );

  cout << "Done intializing" << endl;
}

void ParticleSystem::init( Domain *d, Optimize *o )
{
  _domain = d;
  _optimizer = o;
}

void ParticleSystem::init( svector<DynamicSurfacePoint*> &p)
{
  _points = p;

  // initialize the optimizations
  _optimizer->init( _points, 5 );
}

//------------------------------------------------------------------------
// Function    : removePoint
// Description : remove the Point at index i from the ParticleSystem
//------------------------------------------------------------------------
void ParticleSystem::removePoint( int i )
{
  // remove Point from the Point's Domain
  _domain->removePoint( _points[i] );

  // delete this Point
  delete _points[i];

  // if we are not at the end of the Point list, need to reset
  //   this element's pointer
  if ( unsigned(i) != (_points.size()-1) )
    _points[i] = _points.back();

  // and remove the last element in the list
  _points.pop_back();
}

//------------------------------------------------------------------------
// Function    : addPoint
// Description : add a Point to the ParticleSystem
//------------------------------------------------------------------------
void ParticleSystem::addPoint( DynamicSurfacePoint *point )
{
  // add Point to the Point's Domain
  _domain->addPoint( point );

  // and add the Point to the ParticleSystem
  _points.push_back( point );
}
 
//------------------------------------------------------------------------
// Function    : splitEveryPoint
// Description : 
//------------------------------------------------------------------------
void ParticleSystem::splitEveryPoint()
{
  GlobalSurfaceEnergyNA gse;
  gse.splitEveryPoint( _points );
}

//------------------------------------------------------------------------
// Function    : cleanUpSystem
// Description : delete all the points that do not have a surface 
//               scalar value within +-threshold
//------------------------------------------------------------------------
void ParticleSystem::cleanUpSystem( float threshold )
{
  int num_pts = _points.size();
  for ( int i = 0; i < num_pts; i++ )
  {
    if ( fabs(_points[i]->F()) > threshold )
    {
      removePoint( i-- );
      --num_pts;
    }
  }
}

//------------------------------------------------------------------------
// Function    : offsetPoints
// Description : offset the Points along their normals, and reset their
//               lambda values
//------------------------------------------------------------------------
void ParticleSystem::offsetPoints( float offset_amt )
{
  for ( int i = 0; i < _points.size(); i++ )
  {
    _points[i]->offset( offset_amt );

    if ( _points[i]->moved_outside() )
    {
      removePoint( i );
      --i;
    }
    else
      _points[i]->resetLambda();
  }
}

//------------------------------------------------------------------------
// Function    : surfaceIsovalue
// Description : change the isovalue and offset the points to the new
//               surface
//------------------------------------------------------------------------
void ParticleSystem::surfaceIsovalue( float iso )
{
  // update the isovalue
  _domain->surfaceIsovalue( iso );

  SurfaceConstraint proj;
  for ( int i = 0; i < _points.size(); i++ )
  {
    _points[i]->updateSurfaceParameters();

    // project the point back onto the surface
    proj.projectOntoSurface( _points[i] );

    // and check that the point is still inside the domain
    if ( _points[i]->moved_outside() )
    {
      removePoint( i );
      --i;
    }
    else
      _points[i]->resetLambda();
  }
}

//------------------------------------------------------------------------
// Function    : reinitializePoints
// Description : project the points back onto the surface (and also
//               update their surface values)
//------------------------------------------------------------------------
void ParticleSystem::reinitializePoints()
{
  cout << "Projecting the points onto the surface ... ";

  // and stick the points back onto the surface
  _optimizer->init( _points, 1 );

  cout << "DONE\n";
}

//------------------------------------------------------------------------
// Function    : resetLambdas
// Description : reset all the Points' lambda values
//------------------------------------------------------------------------
void ParticleSystem::resetLambdas()
{
  for ( int i = 0; i < _points.size(); i++ )
    _points[i]->resetLambda();
}

//------------------------------------------------------------------------
// Function    : planar_seperation
// Description : set the new planar seperation in the particles, and 
//               if this is a bigger or smaller than the current, check
//               if we should rebuild the neighborhood
//------------------------------------------------------------------------
void ParticleSystem::planar_seperation( float ps )
{
  // update the variable in each particle
  for ( unsigned i = 0; i < _points.size(); i++ ) 
    _points[i]->planar_seperation( ps );

  // now check if we need to create a new neighborhood
  float new_subd_w = ps / _domain->normalizedIdealDistance(); 
  float old_subd_w = _domain->subd_width();

  if ( (old_subd_w < new_subd_w) || 
       ((0.5*old_subd_w) >= new_subd_w) )
  {
    // the size of the particles has increased or is less than half
    //    the size of the old, so make a new neighborhood
    cout << "Making a new neighborhood, <old,new> = " <<
      old_subd_w << " , " << new_subd_w << endl;

    // delete the old neighborhood and create a new one
    delete _domain->neighborhood();
    vector_type d_start, d_end; 
    _domain->domain( d_start, d_end );
    Neighborhood<DynamicSurfacePoint>* n = 
      new neighborhood_type( new_subd_w, d_start, d_end );

    // pass the new neighborhood to the domain
    _domain->neighborhood( n, new_subd_w );

    // repopulate 
    n->populateNeighborhood( _points );
  }
}
 
//------------------------------------------------------------------------
// Function    : print
// Description : 
//------------------------------------------------------------------------
void ParticleSystem::print() const
{
  float min_r = 100.0;
  float max_r = 0.0;
  float r;
  for ( int i = 0; i < _points.size(); i++ )
  {
    r = 2.0*_points[i]->radius();
    if ( r < min_r )
      min_r = r;
    else if ( r > max_r )
      max_r = r;
  }

  cout << "min radius= " << min_r << "   max radius = " << max_r << endl;
  cout << "Num Points = " << _points.size() << endl << endl;

  //_points[7]->printNeighbors();
  //_points[200]->printNeighbors();
  //_points[400]->printNeighbors();
  //_points[603]->printNeighbors();

  float max_f = 0.0, min_f = 100.0;
  float ave_f = 0.0;
  float f;
  for ( int i = 0; i < _points.size(); i++ )
  {
    //cout << _points[i]->F() << "  " << _points[i]->sf() << endl;
    f = abs(_points[i]->F());
   
    if ( f > max_f )
      max_f = f;
    else if ( f < min_f )
      min_f = f;
    ave_f += log10( f );
  }
  ave_f /= (float)_points.size();
  cout << endl;
  cout << "max_f = " << max_f << "    min_f = " << min_f << endl;
  cout << "ave log10 f = " << ave_f << endl << endl;

  //for ( int i = 0; i < _points.size(); i++ )
  //  cout << _points[i]->lambda() << "  ";
  //cout << endl << endl;
}

//------------------------------------------------------------------------
// Function    : writePointFile
// Description : 
//------------------------------------------------------------------------
void ParticleSystem::writePointFile( const char *filename ) const
{
  if ( !filename )
  {
    cout << "ParticleSystem():: No output filename given." << endl;
    return;
  }

  // open the file
  ofstream out( filename );

  // check that it opened okay
  if ( !out )
  { 
    cout << "ParticleSystem()::Error opening output file" << endl;
    return;
  }

  cout << "Writing point file..." << endl;

  //out << _points.size() << endl;

  vector_type pos;
  for ( int i = 0; i < _points.size(); i++ )
  {
    pos = _points[i]->position();

#ifdef TWO_D
    out << pos(0) << " " << pos(1) << endl;
#endif
#ifdef THREE_D
    out << pos(0) << " " << pos(1) << " " << pos(2) << endl;
#endif
  }  

//  vector_type n;
//  for ( int i = 0; i < _points.size(); i++ )
//  {
//    n = _points[i]->normal();
//
//#ifdef TWO_D
//    out << n(0) << " " << n(1) << endl;
//#endif
//#ifdef THREE_D
//    out << n(0) << " " << n(1) << " " << n(2) << endl;
//#endif
//  }  

  out.close();

  cout << "          ... DONE" << endl;
}

//------------------------------------------------------------------------
// Function    : writeEpsilonSampleFile
// Description : 
//------------------------------------------------------------------------
void ParticleSystem::writeEpsilonSampleFile( const char *filename ) const
{
  if ( !filename )
  {
    cout << "ParticleSystem():: No output filename given." << endl;
    return;
  }

  // open the file
  ofstream out( filename );

  // check that it opened okay
  if ( !out )
  { 
    cout << "ParticleSystem()::Error opening output file" << endl;
    return;
  }

  cout << "Writing epsilon sample file " << filename << endl;

  vector_type pos;
  int counter = 0;
  float closest_pt, sf;
  for ( int i = 0; i < _points.size(); i++ )
  {
    _points[i]->epsilonSample( closest_pt, sf );
  
    if ( closest_pt > sf ) ++counter;

    //out << closest_pt << " " << sf << endl;
  }  

  out.close();
  cout << "          ... DONE" << endl;

  cout << "   " << counter << " points not an epsilon sample out of " <<
    _points.size() << endl;
}

//------------------------------------------------------------------------
// Function    :  
// Description : some timing functions
//------------------------------------------------------------------------
void ParticleSystem::timeEnergyForceComputations()
{
  clock_t _start_time, _end_time;
  double _total_time;
  _start_time = clock();

  float tmpf;
  vector_type tmpv;
  for ( int i = 0; i < _points.size(); i++ )
    _points[i]->computeEnergyForce( tmpf, tmpv );

  _end_time = clock();
  _total_time = 
    (double)(_end_time - _start_time)/(double)CLOCKS_PER_SEC;

  cout << "Time to compute " << _points.size() << " interparicle forces = "
    << _total_time << ",   per particle = " << 
    (_total_time/(float)_points.size()) << endl;
}

void ParticleSystem::timeSurfaceUpdates()
{
  clock_t _start_time, _end_time;
  double _total_time;
  _start_time = clock();

  for ( int i = 0; i < _points.size(); i++ )
    _points[i]->onlyQuerySurface();

  _end_time = clock();
  _total_time = 
    (double)(_end_time - _start_time)/(double)CLOCKS_PER_SEC;

  cout << "Time to compute " << _points.size() << " surface parameters = "
    << _total_time << ",   per particle = " << 
    (_total_time/(float)_points.size()) << endl;
}



  

