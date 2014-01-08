#include <cstdlib>
#include <iostream>
#include <list>
#include <math.h>
#include <fstream>
#include <string>
#include <features/mtxlib.h>
#include <system/defines.h>
#include <system/particles/DynamicSurfacePoint.h>
#include <system/domain/Domain.h>
#include <system/ParticleSystem.h>
#include <system/domain/StaticNeighborhood.h>
#include <system/domain/ApproximateStaticNeighborhood.h>
#include <system/domain/Neighborhood.h>
#include <system/optimization/GlobalSurfaceEnergy.h>
#include <system/optimization/LODController.h>

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
LODController::LODController( Domain *d, ParticleSystem *ps,
                              custom_class::svector<DynamicSurfacePoint*> &p,
                              float initial_sf, float max_surface_sf, 
                              int num_additional_ps,
                              ParticleSystem **additional_ps,
                              int init_num_pts, float threshold,
                              const char* basename,
                              int modulo ) :
  Optimization(), SurfacePopulationController(threshold)
{
  _constraint = new SurfaceConstraint;
  _domain = d;

  _modulo = modulo;

  _min_sf = initial_sf;
  _max_surface_sf = max_surface_sf;

  // store pointers to the additional particle systems
  _num_additional_ps = num_additional_ps;
  _additional_ps = additional_ps;

  buildInitialParticlesAndNeighborhood( ps, p, init_num_pts, basename );
}

LODController::~LODController()
{
  delete _constraint;
  delete [] _additional_ps;
}

//------------------------------------------------------------------------
// Function    : optimize
// Description : 
//------------------------------------------------------------------------
void LODController::optimize( svector<DynamicSurfacePoint*> &points )
{
  _optimized = false;

  if ( _min_sf == 0.0 )
    _optimized = true;
  else
  {
    cout << endl << "Reducing LOD" << endl;
    reduceLOD( points );
  }
}

//------------------------------------------------------------------------
// Function    : buildInitialParticlesAndNeighborhood
// Description : build a course neighborhood with just a few particles
//               to begin with
//------------------------------------------------------------------------
void LODController::buildInitialParticlesAndNeighborhood(
  ParticleSystem *ps, svector<DynamicSurfacePoint*> &points, 
  int init_num_pts, const char* basename )
{ 
  //******************************
  // figure out the initial sizes
  //******************************

  // get the size of the domain
  vector_type d_start, d_end;
  _domain->domain( d_start, d_end );
  vector_type diff = d_end - d_start;

  // compute the initial radius and bin width
  float radius = max(_min_sf, _max_surface_sf) / 
    _domain->normalizedIdealDistance();

  cout << "Binning structure radius = " << radius << endl;


  //******************************
  // create the neighborhood
  //******************************
  
  Neighborhood<DynamicSurfacePoint>* n = 
    new neighborhood_type( radius, d_start, d_end );

  _domain->neighborhood( n, radius );


  //******************************
  // make the list of points
  //******************************

  if ( !basename )
  {
    points.resize( init_num_pts );
    vector_type pos;

    for ( int i = 0; i < init_num_pts; i++ )
    {
      // make a random position for the point in the domain
      pos.random(); 
      pos /= (float)RAND_MAX;
      pos *= diff; 
      pos += d_start;

      points[i] = new DynamicSurfacePoint( _domain, ps, pos );
      points[i]->using_sf();
      points[i]->max_sf( _max_surface_sf );
      points[i]->min_sf( _min_sf );

      // initialize the SurfaceParams
      points[i]->updateSurfaceParameters();
    }
  }
  else
    initializePointsWithMesh( basename, points, ps );


  //******************************
  // populate the domain
  //******************************
  _domain->populateDomain( points );
}

//------------------------------------------------------------------------
// Function    : initializePointsWithMesh
// Description : 
//------------------------------------------------------------------------
void LODController::initializePointsWithMesh( const char *basename,
                                              svector<DynamicSurfacePoint*> 
                                              &points,
                                              ParticleSystem *ps )
{
  cout << "Initializing with a mesh..." << endl;
  char filename[300]; sprintf( filename, "%s.m", basename );

  // open the mesh file
  ifstream in( filename );
  if ( !in )
  {
    cout << "Error reading mesh file " << filename << endl;
    exit( 1 );
  }

  // read the mesh data into the temp arrays
  string buffer;
  DynamicSurfacePoint *point;
  points.resize( 0 );
  vector_type pos;
  int counter = 0;
  while ( in.peek() != EOF )
  {
    getline( in, buffer, ' ' ); 

    if ( buffer == "Vertex" )
    {
      // read in the index and the second space
      getline( in, buffer, ' ' );     
      getline( in, buffer, ' ' );

      // now get the vertex coordinates
      getline( in, buffer, ' ' ); pos[0] = atof( buffer.c_str() );
      getline( in, buffer, ' ' ); pos[1] = atof( buffer.c_str() );
      getline( in, buffer );      
#ifdef THREE_D
      pos[2] = atof( buffer.c_str() );
#endif

      if ( !(counter % _modulo) )
      {
        point = new DynamicSurfacePoint( _domain, ps, pos );
        point->using_sf();
        point->max_sf( _max_surface_sf );
        point->min_sf( _min_sf );

        // initialize the SurfaceParams
        point->updateSurfaceParameters();

        points.push_back( point );
      }

      ++counter;
    }
    else
      getline( in, buffer );
  }

  in.close();
}

//------------------------------------------------------------------------
// Function    : reduceLOD
// Description : 
//------------------------------------------------------------------------
void LODController::reduceLOD(custom_class::svector<DynamicSurfacePoint*> 
                              &points)
{
  // compute the new minimum particle size
  //  --> if the current min is already less than one, then just 
  //      go ahead and set the new min to be zero!
  float new_min_sf;
  if ( _min_sf < 0.5 )
    new_min_sf = 0.0;
  else
    new_min_sf = _min_sf / 2.0;

  cout << "   new LOD:" << new_min_sf << endl;

  //-----------------------------------------------------------------
  // first, check if the new size is bigger than the surface max size
  if ( new_min_sf >= _max_surface_sf )
  {
    // get the new size of the bins
    float bin_width = new_min_sf / _domain->normalizedIdealDistance();

    // delete the current Neighborhood structure and make new one
    //   with smaller bins
    delete _domain->neighborhood();

    vector_type d_start, d_end; 
    _domain->domain( d_start, d_end );
    Neighborhood<DynamicSurfacePoint>* n = 
      new neighborhood_type( bin_width, d_start, d_end );

    _domain->neighborhood( n, bin_width );

    // repopulate with the current batch of particles
    n->populateNeighborhood( points );

    // and add any other particle systems
    for ( int i = 0; i < _num_additional_ps; i++ )
      n->populateNeighborhood( _additional_ps[i]->points() );

    // set the particles' new radius and min size
    for ( int i = 0; i < points.size(); i++ )
    {
      points[i]->min_sf( new_min_sf );
      //points[i]->min_sf( max(0.1f,new_min_sf) );

      // just for rendering and debugging
      points[i]->radius( 0.5*points[i]->sf() );
    }

    // split the current particles into smaller ones
#ifdef THREE_D
    splitEveryPointIntoFourWSF( points );
#else
    splitEveryPoint( points );
#endif
  }

  //-----------------------------------------------------------
  // we know the new size is smaller than the max size we can
  //   have, so check if still need to make a new binning 
  //   structure (ie. the last size was still bigger)
  else if ( _min_sf > _max_surface_sf )
  {
    // get the new size of the bins
    float bin_width = _max_surface_sf / _domain->normalizedIdealDistance();

    // delete the current Neighborhood structure and make new one
    //   with smaller bins
    delete _domain->neighborhood();

    vector_type d_start, d_end; 
    _domain->domain( d_start, d_end );
    Neighborhood<DynamicSurfacePoint>* n = 
      new neighborhood_type( bin_width, d_start, d_end );

    _domain->neighborhood( n, bin_width );

    // repopulate with the current batch of particles
    n->populateNeighborhood( points );

    // and add any other particle systems
    for ( int i = 0; i < _num_additional_ps; i++ )
      n->populateNeighborhood( _additional_ps[i]->points() );

    // split the points 
    int original_num_pts = points.size();
    for ( int i = 0; i < original_num_pts; i++ )
    {
      points[i]->min_sf( new_min_sf );
      //points[i]->min_sf( max(0.1f,new_min_sf) );

      // just for rendering and debugging
      points[i]->radius( 0.5*points[i]->sf() );

      // split into 4 if we can
      if ( points[i]->sf() == new_min_sf )
#ifdef THREE_D
        splitIntoFourWSF( i, points[i] );
#else
      splitAPoint( i, points[i] );
#endif
        

      //// otherwise, just split
      //else
      //  splitAPoint( i, points[i] );
    }

    // do a pass to check if we should split some particles
    GlobalSurfaceEnergyNA gse( true );
    gse.optimize( points );
  }

  //-----------------------------------------------------------
  // finally, we know that the previous and current size is
  //   smaller than the max surface size (no new neighborhood!)
  else
  {
    int original_num_pts = points.size();
    for ( int i = 0; i < original_num_pts; i++ )
    {
      points[i]->min_sf( new_min_sf );
      //points[i]->min_sf( max(0.1f,new_min_sf) );

      // just for rendering and debugging
      points[i]->radius( 0.5*points[i]->sf() );

      // split into 4 if we can
      if ( points[i]->sf() == new_min_sf )
#ifdef THREE_D
        splitIntoFourWSF( i, points[i] );
#else
        splitAPoint( i, points[i] );
#endif

      //// otherwise, just split
      //else
      //  splitAPoint( i, points[i] );
    }

    // do a pass to check if we should split some particles
    GlobalSurfaceEnergyNA gse( true );
    gse.optimize( points );
  }

  _min_sf = new_min_sf;

  cout << "   number of points = " << points.size() << endl;
}





