#include <cstdlib>
#include <iostream>
#include <fstream>
#include <domain/Domain.h>
#include <optimization/Optimize.h>
#include <domain/SimpleSurfaces.h>
#include <domain/DistanceField.h>
#include <optimization/AdaptiveDtSurfaceDistribution.h>
#include <optimization/LMSurfaceDistribution.h>
#include <optimization/GlobalDtSurfaceDistribution.h>
#include <optimization/JacobianGlobalDtSurfaceDistribution.h>
#include <optimization/GlobalSurfaceEnergy.h>
#include <optimization/ConstNumParticles.h>
#include <optimization/ApproximateNeighborhoods.h>
#include <domain/MPUSurface.h>
#include <domain/MLSSurface.h>
#include <SingleSystem.h>

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
SingleSystem::SingleSystem( const char *param_file )
{
  _particle_sys = new ParticleSystem();

  // domain variables
  Surface *surface = NULL;
  Energy *energy = NULL;
  vector_type start(0.0), end(1.0);
  float subdomain_w=0.0;

  // point variables
  float planar_seperation=1.0, angular_density=0.0, influence_radius=1.0;
  int init_num_points=10;
  bool global_influence_radius;

  // optimizer variables
  Optimization **ops=NULL;
  ops = new Optimization*[3]; // for now, max of 3 optimations
  int num_ops = 0;

  bool biased_splitting_dying;

  // now, pass all the variables in by reference to the file reader
  readParamFile( param_file, surface, energy, start, end, subdomain_w,
                 planar_seperation, angular_density, 
                 global_influence_radius, influence_radius,
                 init_num_points, ops, num_ops, biased_splitting_dying );

  // create the domain
  Domain *domain = new Domain( energy, surface, start, end, subdomain_w );

  // create the points
  svector<DynamicSurfacePoint*> points;
  points.resize( init_num_points );
  vector_type pos, diff=end-start;
  float r = /*25.0*/0.9;
  for ( int i = 0; i < init_num_points; i++ )
  {
    // make a random position for the point in the domain
    pos.random(); 
    pos /= (float)RAND_MAX;
#ifdef THREE_D

    pos *= diff; 
    pos += start;
    pos.normalize();
    pos *=r;

    //theta = pos(0)*2.0*PI;
    //phi = pos(1)*2.0*PI;
    //pos(0) = r*cos(theta)*sin(phi);
    //pos(1) = r*sin(theta)*sin(phi);
    //pos(2) = r*cos(phi);

    //pos += 1.0;
    //pos *= 0.5;

#endif

    
    points[i] = new DynamicSurfacePoint( domain, _particle_sys, pos );

    points[i]->storageInt( i );
    /*_original_pos[i] = pos;
    _previous_pos[i] = pos;*/

    // set the other variables
    points[i]->angular_density( angular_density );
    points[i]->planar_seperation( planar_seperation );

    // initialize the SurfaceParams
    points[i]->updateSurfaceParameters();
  }

  // populate the domain
  domain->populateDomain( points );

  // create the optimizations
  //ops[num_ops++] = new GlobalSurfaceEnergyNA( biased_splitting_dying );
  ops[num_ops++] = new ConstNumParticles( init_num_points );
  Optimize *optimizer = new Optimize( ops, num_ops );
  
  _particle_sys->init( domain, optimizer, points );

  //// open the files
  //_out_file = new ofstream**[3];
  //for ( int i = 0; i < 3; i++ )
  //  _out_file[i] = new ofstream*[4];

  //// the distance from the orgins
  //_out_file[0][0] = new ofstream( "dist_from_origin_0.txt" );
  //_out_file[0][1] = new ofstream( "dist_from_origin_1.txt" );
  //_out_file[0][2] = new ofstream( "dist_from_origin_2.txt" );
  //_out_file[0][3] = new ofstream( "dist_from_origin_3.txt" );

  //// movement each iteration
  //_out_file[1][0] = new ofstream( "dist_from_previous_0.txt" );
  //_out_file[1][1] = new ofstream( "dist_from_previous_1.txt" );
  //_out_file[1][2] = new ofstream( "dist_from_previous_2.txt" );
  //_out_file[1][3] = new ofstream( "dist_from_previous_3.txt" );

  //// distance from center point
  //_out_file[2][0] = new ofstream( "dist_from_center_0.txt" );
  //_out_file[2][1] = new ofstream( "dist_from_center_1.txt" );
  //_out_file[2][2] = new ofstream( "dist_from_center_2.txt" );
  //_out_file[2][3] = new ofstream( "dist_from_center_3.txt" );
}

SingleSystem::~SingleSystem()
{
  //_out_file[0][0]->close();
  //_out_file[0][1]->close();
  //_out_file[0][2]->close();
  //_out_file[0][3]->close();

  //_out_file[1][0]->close();
  //_out_file[1][1]->close();
  //_out_file[1][2]->close();
  //_out_file[1][3]->close();

  //_out_file[2][0]->close();
  //_out_file[2][1]->close();
  //_out_file[2][2]->close();
  //_out_file[2][3]->close();

  //delete _out_file[0][0];
  //delete _out_file[0][1];
  //delete _out_file[0][2];
  //delete _out_file[0][3];

  //delete _out_file[1][0];
  //delete _out_file[1][1];
  //delete _out_file[1][2];
  //delete _out_file[1][3];

  //delete _out_file[2][0];
  //delete _out_file[2][1];
  //delete _out_file[2][2];
  //delete _out_file[2][3];

  //delete [] _out_file[0];
  //delete [] _out_file[1];
  //delete [] _out_file[2];

  //delete [] _out_file;

  delete _particle_sys;
}

//------------------------------------------------------------------------
// Function    : Optimize
// Description : 
//------------------------------------------------------------------------
void SingleSystem::optimize()
{
  _particle_sys->optimize(); 

  //// write the neighborhood stuff out to the files
  //vector_type v;

  //// distance from original positions
  //v = _particle_sys->point(256)->position() - _original_pos[256];
  //*_out_file[0][0] << v.length() << " ";
  //v = _particle_sys->point(720)->position() - _original_pos[720];
  //*_out_file[0][0] << v.length() << " ";
  //v = _particle_sys->point(97)->position() - _original_pos[97];
  //*_out_file[0][0] << v.length() << " ";
  //v = _particle_sys->point(182)->position() - _original_pos[182];
  //*_out_file[0][0] << v.length() << " ";
  //v = _particle_sys->point(351)->position() - _original_pos[351];
  //*_out_file[0][0] << v.length() << " ";
  //v = _particle_sys->point(798)->position() - _original_pos[798];
  //*_out_file[0][0] << v.length() << " " << endl;

  //v = _particle_sys->point(757)->position() - _original_pos[757];
  //*_out_file[0][1] << v.length() << " ";
  //v = _particle_sys->point(661)->position() - _original_pos[661];
  //*_out_file[0][1] << v.length() << " ";
  //v = _particle_sys->point(561)->position() - _original_pos[561];
  //*_out_file[0][1] << v.length() << " ";
  //v = _particle_sys->point(50)->position() - _original_pos[50];
  //*_out_file[0][1] << v.length() << " ";
  //v = _particle_sys->point(336)->position() - _original_pos[336];
  //*_out_file[0][1] << v.length() << " ";
  //v = _particle_sys->point(91)->position() - _original_pos[91];
  //*_out_file[0][1] << v.length() << " " << endl;

  //v = _particle_sys->point(509)->position() - _original_pos[509];
  //*_out_file[0][2] << v.length() << " ";
  //v = _particle_sys->point(769)->position() - _original_pos[769];
  //*_out_file[0][2] << v.length() << " ";
  //v = _particle_sys->point(599)->position() - _original_pos[599];
  //*_out_file[0][2] << v.length() << " ";
  //v = _particle_sys->point(424)->position() - _original_pos[424];
  //*_out_file[0][2] << v.length() << " ";
  //v = _particle_sys->point(173)->position() - _original_pos[173];
  //*_out_file[0][2] << v.length() << " ";
  //v = _particle_sys->point(380)->position() - _original_pos[380];
  //*_out_file[0][2] << v.length() << " " << endl;

  //v = _particle_sys->point(265)->position() - _original_pos[265];
  //*_out_file[0][3] << v.length() << " ";
  //v = _particle_sys->point(546)->position() - _original_pos[546];
  //*_out_file[0][3] << v.length() << " ";
  //v = _particle_sys->point(487)->position() - _original_pos[487];
  //*_out_file[0][3] << v.length() << " ";
  //v = _particle_sys->point(465)->position() - _original_pos[465];
  //*_out_file[0][3] << v.length() << " ";
  //v = _particle_sys->point(193)->position() - _original_pos[193];
  //*_out_file[0][3] << v.length() << " ";
  //v = _particle_sys->point(701)->position() - _original_pos[701];
  //*_out_file[0][3] << v.length() << " " << endl;


  //// distance moved durning iteration
  //v = _particle_sys->point(256)->position() - _previous_pos[256];
  //*_out_file[1][0] << v.length() << " ";
  //v = _particle_sys->point(720)->position() - _previous_pos[720];
  //*_out_file[1][0] << v.length() << " ";
  //v = _particle_sys->point(97)->position() - _previous_pos[97];
  //*_out_file[1][0] << v.length() << " ";
  //v = _particle_sys->point(182)->position() - _previous_pos[182];
  //*_out_file[1][0] << v.length() << " ";
  //v = _particle_sys->point(351)->position() - _previous_pos[351];
  //*_out_file[1][0] << v.length() << " ";
  //v = _particle_sys->point(798)->position() - _previous_pos[798];
  //*_out_file[1][0] << v.length() << " " << endl;
  //_previous_pos[256] = _particle_sys->point(256)->position();
  //_previous_pos[720] = _particle_sys->point(720)->position();
  //_previous_pos[97]  = _particle_sys->point(97)->position();
  //_previous_pos[182] = _particle_sys->point(182)->position();
  //_previous_pos[351] = _particle_sys->point(351)->position();
  //_previous_pos[798] = _particle_sys->point(798)->position();

  //v = _particle_sys->point(757)->position() - _previous_pos[757];
  //*_out_file[1][1] << v.length() << " ";
  //v = _particle_sys->point(661)->position() - _previous_pos[661];
  //*_out_file[1][1] << v.length() << " ";
  //v = _particle_sys->point(561)->position() - _previous_pos[561];
  //*_out_file[1][1] << v.length() << " ";
  //v = _particle_sys->point(50)->position() - _previous_pos[50];
  //*_out_file[1][1] << v.length() << " ";
  //v = _particle_sys->point(336)->position() - _previous_pos[336];
  //*_out_file[1][1] << v.length() << " ";
  //v = _particle_sys->point(91)->position() - _previous_pos[91];
  //*_out_file[1][1] << v.length() << " " << endl;
  //_previous_pos[757] = _particle_sys->point(757)->position();
  //_previous_pos[661] = _particle_sys->point(661)->position();
  //_previous_pos[561] = _particle_sys->point(561)->position();
  //_previous_pos[50]  = _particle_sys->point(50)->position();
  //_previous_pos[336] = _particle_sys->point(336)->position();
  //_previous_pos[91]  = _particle_sys->point(91)->position();

  //v = _particle_sys->point(509)->position() - _previous_pos[509];
  //*_out_file[1][2] << v.length() << " ";
  //v = _particle_sys->point(769)->position() - _previous_pos[769];
  //*_out_file[1][2] << v.length() << " ";
  //v = _particle_sys->point(599)->position() - _previous_pos[599];
  //*_out_file[1][2] << v.length() << " ";
  //v = _particle_sys->point(424)->position() - _previous_pos[424];
  //*_out_file[1][2] << v.length() << " ";
  //v = _particle_sys->point(173)->position() - _previous_pos[173];
  //*_out_file[1][2] << v.length() << " ";
  //v = _particle_sys->point(380)->position() - _previous_pos[380];
  //*_out_file[1][2] << v.length() << " " << endl;
  //_previous_pos[509] = _particle_sys->point(509)->position();
  //_previous_pos[769] = _particle_sys->point(769)->position();
  //_previous_pos[599] = _particle_sys->point(599)->position();
  //_previous_pos[424] = _particle_sys->point(424)->position();
  //_previous_pos[173] = _particle_sys->point(173)->position();
  //_previous_pos[380] = _particle_sys->point(380)->position();


  //v = _particle_sys->point(265)->position() - _previous_pos[265];
  //*_out_file[1][3] << v.length() << " ";
  //v = _particle_sys->point(546)->position() - _previous_pos[546];
  //*_out_file[1][3] << v.length() << " ";
  //v = _particle_sys->point(487)->position() - _previous_pos[487];
  //*_out_file[1][3] << v.length() << " ";
  //v = _particle_sys->point(465)->position() - _previous_pos[465];
  //*_out_file[1][3] << v.length() << " ";
  //v = _particle_sys->point(193)->position() - _previous_pos[193];
  //*_out_file[1][3] << v.length() << " ";
  //v = _particle_sys->point(701)->position() - _previous_pos[701];
  //*_out_file[1][3] << v.length() << " " << endl;
  //_previous_pos[265] = _particle_sys->point(265)->position();
  //_previous_pos[546] = _particle_sys->point(546)->position();
  //_previous_pos[487] = _particle_sys->point(487)->position();
  //_previous_pos[465] = _particle_sys->point(465)->position();
  //_previous_pos[193] = _particle_sys->point(193)->position();
  //_previous_pos[701] = _particle_sys->point(701)->position();


  //// distance from center point
  //v = _particle_sys->point(256)->position() - 
  //  _particle_sys->point(7)->position();
  //*_out_file[2][0] << v.length() << " ";
  //v = _particle_sys->point(720)->position() - 
  //  _particle_sys->point(7)->position();
  //*_out_file[2][0] << v.length() << " ";
  //v = _particle_sys->point(97)->position() - 
  //  _particle_sys->point(7)->position();
  //*_out_file[2][0] << v.length() << " ";
  //v = _particle_sys->point(182)->position() - 
  //  _particle_sys->point(7)->position();
  //*_out_file[2][0] << v.length() << " ";
  //v = _particle_sys->point(351)->position() - 
  //  _particle_sys->point(7)->position();
  //*_out_file[2][0] << v.length() << " ";
  //v = _particle_sys->point(798)->position() - 
  //  _particle_sys->point(7)->position();
  //*_out_file[2][0] << v.length() << " " << endl;

  //v = _particle_sys->point(757)->position() - 
  //  _particle_sys->point(200)->position();
  //*_out_file[2][1] << v.length() << " ";
  //v = _particle_sys->point(661)->position() - 
  //  _particle_sys->point(200)->position();
  //*_out_file[2][1] << v.length() << " ";
  //v = _particle_sys->point(561)->position() - 
  //  _particle_sys->point(200)->position();
  //*_out_file[2][1] << v.length() << " ";
  //v = _particle_sys->point(50)->position() - 
  //  _particle_sys->point(200)->position();
  //*_out_file[2][1] << v.length() << " ";
  //v = _particle_sys->point(336)->position() - 
  //  _particle_sys->point(200)->position();
  //*_out_file[2][1] << v.length() << " ";
  //v = _particle_sys->point(91)->position() - 
  //  _particle_sys->point(200)->position();
  //*_out_file[2][1] << v.length() << " " << endl;

  //v = _particle_sys->point(509)->position() - 
  //  _particle_sys->point(400)->position();
  //*_out_file[2][2] << v.length() << " ";
  //v = _particle_sys->point(769)->position() - 
  //  _particle_sys->point(400)->position();
  //*_out_file[2][2] << v.length() << " ";
  //v = _particle_sys->point(599)->position() - 
  //  _particle_sys->point(400)->position();
  //*_out_file[2][2] << v.length() << " ";
  //v = _particle_sys->point(424)->position() - 
  //  _particle_sys->point(400)->position();
  //*_out_file[2][2] << v.length() << " ";
  //v = _particle_sys->point(173)->position() - 
  //  _particle_sys->point(400)->position();
  //*_out_file[2][2] << v.length() << " ";
  //v = _particle_sys->point(380)->position() - 
  //  _particle_sys->point(400)->position();
  //*_out_file[2][2] << v.length() << " " << endl;

  //v = _particle_sys->point(265)->position() - 
  //  _particle_sys->point(603)->position();
  //*_out_file[2][3] << v.length() << " ";
  //v = _particle_sys->point(546)->position() - 
  //  _particle_sys->point(603)->position();
  //*_out_file[2][3] << v.length() << " ";
  //v = _particle_sys->point(487)->position() - 
  //  _particle_sys->point(603)->position();
  //*_out_file[2][3] << v.length() << " ";
  //v = _particle_sys->point(465)->position() - 
  //  _particle_sys->point(603)->position();
  //*_out_file[2][3] << v.length() << " ";
  //v = _particle_sys->point(193)->position() - 
  //  _particle_sys->point(603)->position();
  //*_out_file[2][3] << v.length() << " ";
  //v = _particle_sys->point(701)->position() - 
  //  _particle_sys->point(603)->position();
  //*_out_file[2][3] << v.length() << " " << endl;

}

//------------------------------------------------------------------------
// Function    : generateMPUSurface
// Description : 
//------------------------------------------------------------------------
void SingleSystem::generatePointSetSurface( int type )
{
  cout << "Generating the point set surface ... ";

  // create the point set surface
  //
  // get the info needed by the Neighborhood in the point set surface
  vector_type start, end;
  float subd_w;
  _particle_sys->domain( start, end, subd_w );

  // instantiate a new point set surface
  Surface *tmp;
  switch ( type )
  {
  case MPU:
    tmp =  new MPUSurface( start, end, subd_w, 
                           _mpu_weighting_function );
    break;

  case MLS: 
    tmp = new MLSSurface( start, end, subd_w, 
                          _mpu_weighting_function );
    break;
  }
        
  // pass the point set surface the Points from which to generate 
  //   the surface from
  tmp->generateSurfaceFromPoints( _particle_sys->points() );

  // and update the domain with this point set surface
  _particle_sys->surface( tmp );  

  // reset the Points lambdas
  _particle_sys->resetLambdas();

  cout << "DONE\n";
}

//------------------------------------------------------------------------
// Function    : readParamFile
// Description : read in the parameter file, and set the Domain 
//               variables
//------------------------------------------------------------------------
void SingleSystem::readParamFile( const char *param_file, Surface* &s, 
                            Energy* &e, 
                            vec<3> &start, vec<3> &end, 
                            float &subd_w, float &planar_sep, 
                            float &angular_den, 
                            bool &global_influence_radius,
                            float &influence_r,
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

  else if ( !strncmp( value_s, "distancefield", 13 ) )
  {
    fscanf( file, "%s%s", v, value_s );
    int xdim, ydim, zdim;
    fscanf( file, "%s %i,%i,%i", v, &xdim, &ydim, &zdim );
    s = new DistanceField( xdim, ydim, zdim, value_s );
    /*s = new DistanceField( 100, 100, 100, 
      "C:/home/research/volumes/data/sphere.vol",
      "C:/home/research/volumes/data/sphere_dx.vol",
      "C:/home/research/volumes/data/sphere_dy.vol",
      "C:/home/research/volumes/data/sphere_dz.vol" );*/

    start(0) = start(1) = start(2) = 0.0;
    end(0) = (float)(xdim-1); 
    end(1) = (float)(ydim-1); 
    end(2) = (float)(zdim-1);
    //end(0) = 99.0; end(1) = 99.0; end(2) = 99.0;


    cout << "     SURFACE : distance field\n";
    cout << "     DOMAIN : [" << start(0) << "," << start(1) << "," <<
      start(2) << "] x [" << end(0) << "," << end(1) << "," << end(2) <<
      "]\n";

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

  // SPATIAL BIN WIDTH
  fscanf( file, "%s%f", v, &subd_w );
  cout << "     BIN WIDTH : " << subd_w << endl;

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

  fscanf( file, "%s%f", v, &influence_r );
  cout << "     RADIUS OF INFLUENCE : " << influence_r << endl;

  fscanf( file, "%s%f", v, &planar_sep );
  cout << "     PLANAR SEPERATION : " << planar_sep << endl;

  fscanf( file, "%s%f", v, &angular_den );
  cout << "     ANGULAR DENSITY : " << angular_den << endl;
  
  // SingleSystem INITIALIZATION
  fscanf( file, "%s%i", v, &init_num_pts );
  cout << "     INITIAL NUMBER OF POINTS : " << init_num_pts << endl;

  // OPTIMIZATIONS
  fscanf( file, "%s%s", v, value_s );
  if ( !strncmp( value_s, "adaptive_dt", 11 ) )
  {
    float dt;
    fscanf( file, "%s%f", v, &dt );

    ops[num_ops++] = new AdaptiveDtSurfaceDistribution();
    cout << "     DISTRIBUTOR : adaptive dt\n";
  }
  else if ( !strncmp( value_s, "lm", 2 ) )
  {
    float dt;
    fscanf( file, "%s%f", v, &dt );

    ops[num_ops++] = new LMSurfaceDistribution();
    cout << "     DISTRIBUTOR : lm method" << endl;
  }
  else if ( !strncmp( value_s, "global_dt", 9 ) )
  {
    float dt;
    fscanf( file, "%s%f", v, &dt );

    ops[num_ops++] = new GlobalDtSurfaceDistribution( dt );
    cout << "     DISTRIBUTOR : global dt = " << dt << endl;
  }
  else if ( !strncmp( value_s, "jacobian_global_dt", 18 ) )
  {
    float dt;
    fscanf( file, "%s%f", v, &dt );

    ops[num_ops++] = new JacobianGlobalDtSurfaceDistribution( dt );
    cout << "     DISTRIBUTOR : jacobian global dt = " << dt << endl;
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

void SingleSystem::readParamFile( const char *param_file, Surface* &s, 
                            Energy* &e, 
                            vec<2> &start, vec<2> &end, 
                            float &subd_w, float &planar_sep, 
                            float &angular_den, 
                            bool &global_influence_radius,
                            float &influence_r,
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

  char v[30], value_s[30];

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
  fscanf( file, "%s%s", v, value_s );
  if ( !strncmp( value_s, "quarticcube", 11 ) )
  {
    float r = 0.26;
    //r *= 10.0;
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
  else if ( !strncmp( value_s, "ellipse", 7 ) )
  {
    s = new SimpleEllipse(3.0,0.0,5.0,0.0,0.0,-0.8);
    cout << "     SURFACE : ellipse\n";
  }

  // DOMAIN DIMENSIONS
  fscanf( file, "%s%f,%f,%f", v, &start(0), &start(1) );
  fscanf( file, "%s%f,%f,%f", v, &end(0), &end(1) );
  cout << "     DOMAIN : [" << start << "] x [" << end << "]\n";

  // SPATIAL BIN WIDTH
  fscanf( file, "%s%f", v, &subd_w );
  cout << "     BIN WIDTH : " << subd_w << endl;

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

  fscanf( file, "%s%f", v, &influence_r );
  cout << "     RADIUS OF INFLUENCE : " << influence_r << endl;

  fscanf( file, "%s%f", v, &planar_sep );
  cout << "     PLANAR SEPERATION : " << planar_sep << endl;

  fscanf( file, "%s%f", v, &angular_den );
  cout << "     ANGULAR DENSITY : " << angular_den << endl;
  
  // SingleSystem INITIALIZATION
  fscanf( file, "%s%i", v, &init_num_pts );
  cout << "     INITIAL NUMBER OF POINTS : " << init_num_pts << endl;

  // OPTIMIZATIONS
  fscanf( file, "%s%s", v, value_s );
  if ( !strncmp( value_s, "adaptive_dt", 11 ) )
  {
    float dt;
    fscanf( file, "%s%f", v, &dt );

    ops[num_ops++] = new AdaptiveDtSurfaceDistribution();
    cout << "     DISTRIBUTOR : adaptive dt\n";
  }
  else if ( !strncmp( value_s, "lm", 2 ) )
  {
    float dt;
    fscanf( file, "%s%f", v, &dt );

    ops[num_ops++] = new LMSurfaceDistribution();
    cout << "     DISTRIBUTOR : lm method" << endl;
  }
  else if ( !strncmp( value_s, "global_dt", 9 ) )
  {
    float dt;
    fscanf( file, "%s%f", v, &dt );

    ops[num_ops++] = new GlobalDtSurfaceDistribution( dt );
    cout << "     DISTRIBUTOR : global dt = " << dt << endl;
  }
  else if ( !strncmp( value_s, "jacobian_global_dt", 18 ) )
  {
    float dt;
    fscanf( file, "%s%f", v, &dt );

    ops[num_ops++] = new JacobianGlobalDtSurfaceDistribution( dt );
    cout << "     DISTRIBUTOR : jacobian global dt = " << dt << endl;
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

  // MPU weighting function
  fscanf( file, "%s%s", v, value_s );
  if  ( !strncmp( value_s, "radialcubed", 11 ) )
  {
    _mpu_weighting_function = 't';
    cout << "     MPU weighting function : radial cubed\n";
  }
  else if ( !strncmp( value_s, "radial", 6 ) )
  {
    _mpu_weighting_function = 'r';
    cout << "     MPU weighting function : radial squared\n";
  }
  else if  ( !strncmp( value_s, "cotan", 5 ) )
  {
    _mpu_weighting_function = 'c';
    cout << "     MPU weighting function : cotan\n";
  }
  else if  ( !strncmp( value_s, "gaussian", 8 ) )
  {
    _mpu_weighting_function = 'g';
    cout << "     MPU weighting function : gaussian\n";
  }
  else if  ( !strncmp( value_s, "cubicspline", 11 ) )
  {
    _mpu_weighting_function = 's';
    cout << "     MPU weighting function : cubic spline\n";
  }


  fclose( file );
}

  

