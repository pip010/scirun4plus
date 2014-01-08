#include <cstdlib>
#include <iostream>
#include <fstream>
#include <system/domain/Domain.h>
#include <system/optimization/Optimize.h>
#include <system/domain/SimpleSurfaces.h>
//#include <system/domain/MorseFunctions.h>
#include <system/domain/DistanceField.h>
//#include <system/domain/TransformationTester.h>
#include <system/optimization/AdaptiveDtSurfaceDistribution.h>
//#include <system/optimization/LMSurfaceDistribution.h>
#include <system/optimization/GlobalSurfaceEnergy.h>
//#include <system/optimization/ConstNumParticles.h>
#include <system/optimization/ApproximateNeighborhoods.h>
#include <system/optimization/LODController.h>
#include <system/LODSystem.h>

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
LODSystem::LODSystem( const char *param_file )
{
  _particle_sys = new ParticleSystem();

  // domain variables
  Surface *surface = NULL;
  Energy *energy = NULL;
  vector_type start(0.0), end(1.0);

  // point variables
  float feature_size;
  int init_num_points=10;
  bool global_influence_radius;

  // optimizer variables
  Optimization **ops=NULL;
  ops = new Optimization*[4]; 
  int num_ops = 0;
  ops[num_ops++] = new ApproximateNeighborhoods( 6 );

  bool biased_splitting_dying;

  _threshold = 0.1;

  // now, pass all the variables in by reference to the file reader
  readParamFile( param_file, surface, energy, start, end,
                 _planar_seperation, _angular_density, 
                 global_influence_radius, feature_size,
                 init_num_points, ops, num_ops, biased_splitting_dying );

  // create the surface density optimization
  ops[num_ops++] = new GlobalSurfaceEnergyNA( biased_splitting_dying,
                                              _threshold );

  // create the domain
  Domain *domain = new Domain( energy, surface, start, end );

  // create the LOD controller
  svector<DynamicSurfacePoint*> points;
  ops[num_ops++] = new LODController( domain, _particle_sys,
                                      points, _planar_seperation,
                                      _angular_density, 
                                      feature_size, init_num_points, 0.25,
                                      _threshold );

  // create the optimizer
  Optimize *optimizer = new Optimize( ops, num_ops );
  
  // and initialize the system
  _particle_sys->init( domain, optimizer, points );
}

LODSystem::~LODSystem()
{
  delete _particle_sys;
}

//------------------------------------------------------------------------
// Function    : readParamFile
// Description : read in the parameter file, and set the Domain 
//               variables
//------------------------------------------------------------------------
void LODSystem::readParamFile( const char *param_file, Surface* &s, 
                            Energy* &e, 
                            vec<3> &start, vec<3> &end, 
                            float &planar_sep, 
                            float &angular_den, 
                            bool &global_influence_radius,
                            float &feature_size,
                            int &init_num_pts, Optimization** &ops,
                            int &num_ops, bool &biased_splitting_dying )
{
  // read in the file
  FILE *file = fopen( param_file, "r" );
  if ( !file )
  {
    cout << "Parameter File Could Not Be Opened...\n";
    return;
  }

  cout << "READING in parameter file " << param_file << endl;

  char v[100], value_s[100];

  // ENERGY
  fscanf( file, "%s%s", v, value_s );
  if ( !strncmp( value_s, "cotan", 5 ) )
  {
    e = new CotanEnergy();
    cout << "     ENERGY : cotan\n";
  }
  else if ( !strncmp( value_s, "radial", 6 ) )
  {
    e = new RadialEnergy();
    cout << "     ENERGY : radial\n";
  }

  // SURFACE
  bool fe_surface = false;
  fscanf( file, "%s%s", v, value_s );
  if ( !strncmp( value_s, "quarticcube", 11 ) )
  {
    float r = 0.26;
    s = new SimpleQuarticCube( r );
    cout << "     SURFACE : quartic cube\n";
  }
  else if ( !strncmp( value_s, "sphere", 6 ) )
  {
    float r = 0.9;
    s = new SimpleSphere( r );
    cout << "     SURFACE : sphere\n";
  }
  else if ( !strncmp( value_s, "torus", 5 ) )
  {
    float r = 0.2;
    float R = 0.7;
    s = new SimpleTorus( r, R );
    cout << "     SURFACE : torus\n";
  }
  else if ( !strncmp( value_s, "morse", 5 ) )
  {
    float f = 2.0;
    s = new MarschnerLobb( f );
    cout << "     SURFACE : morse\n";
  }

  else if ( !strncmp( value_s, "transtester", 11 ) )
  {
    s = new TransformationTester();
    cout << "     SURFACE : transformation tester\n";
  }

  else if ( !strncmp( value_s, "distancefield", 13 ) )
  {
    fscanf( file, "%s%s", v, value_s );
    int xdim, ydim, zdim;
    fscanf( file, "%s %i,%i,%i", v, &xdim, &ydim, &zdim );
    s = new DistanceField( xdim, ydim, zdim, value_s );

    start(0) = start(1) = start(2) = 0.0;
    end(0) = (float)(xdim-1); 
    end(1) = (float)(ydim-1); 
    end(2) = (float)(zdim-1);

    cout << "     SURFACE : distance field\n";
    cout << "     DOMAIN : [" << start(0) << "," << start(1) << "," <<
      start(2) << "] x [" << end(0) << "," << end(1) << "," << end(2) <<
      "]\n";

    cout << "     MAX VALUE : " << s->max_isovalue() << endl;
    cout << "     MIN VALUE : " << s->min_isovalue() << endl;

    fe_surface = true;
  }

  // DOMAIN DIMENSIONS
  if ( !fe_surface )
  {
    fscanf( file, "%s%f,%f,%f", v, &start(0), &start(1), &start(2) );
    fscanf( file, "%s%f,%f,%f", v, &end(0), &end(1), &end(2) );
    cout << "     DOMAIN : [" << start(0) << "," << start(1) << "," <<
      start(2) << "] x [" << end(0) << "," << end(1) << "," << end(2) <<
      "]\n";
  }

  // POINT VARIABLES
  fscanf( file, "%s%s", v, value_s );
  if ( !strncmp( value_s, "global", 6 ) )
  {
    global_influence_radius = true;
    cout << "     INFLUENCE RADIUS TYPE : global\n";
  }
  else if ( !strncmp( value_s, "adaptive", 8 ) )
  {
    global_influence_radius = false;
    cout << "     INFLUENCE RADIUS TYPE : adaptive\n";
  }

  fscanf( file, "%s%f", v, &planar_sep );
  cout << "     PLANAR SEPERATION : " << planar_sep << endl;

  fscanf( file, "%s%f", v, &angular_den );
  cout << "     ANGULAR DENSITY : " << angular_den << endl;

  fscanf( file, "%s%f", v, &feature_size );
  cout << "     FEATURE SIZE : " << feature_size << endl;
  
  // LODSystem INITIALIZATION
  fscanf( file, "%s%i", v, &init_num_pts );
  cout << "     INITIAL NUMBER OF POINTS : " << init_num_pts << endl;

  // OPTIMIZATIONS
  fscanf( file, "%s%s", v, value_s );
  if ( !strncmp( value_s, "adaptive_dt", 11 ) )
  {
    float dt;
    fscanf( file, "%s%f", v, &dt );

    ops[num_ops++] = new AdaptiveDtSurfaceDistribution(_threshold);
    cout << "     DISTRIBUTOR : adaptive dt\n";
  }
  else if ( !strncmp( value_s, "lm", 2 ) )
  {
    float dt;
    fscanf( file, "%s%f", v, &dt );

    ops[num_ops++] = new LMSurfaceDistribution();
    cout << "     DISTRIBUTOR : lm method" << endl;
  }

  fscanf( file, "%s%s", v, value_s );
  if ( !strncmp( value_s, "biased", 6 ) )
  {
    biased_splitting_dying = true;
    cout << "     SPLITTING and DYING : biased\n";
  }
  else if  ( !strncmp( value_s, "random", 6 ) )
  {
    biased_splitting_dying = false;
    cout << "     SPLITTING and DYING : random\n";
  }


  fclose( file );
}



