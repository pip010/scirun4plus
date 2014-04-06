#if 0
  #ifdef _WIN32
  #include <windows.h>
  #endif
  #include <GL/glut.h>
#endif
#include <algorithm>
#include <locale>
#include <functional>
#include <cstdlib>
#include <ctype.h>
#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include <system/domain/Domain.h>
#include <system/optimization/Optimize.h>
#include <system/domain/DistanceFieldwSF.h>
#include <system/domain/ScalarField.h>
#include <system/domain/Intersection.h>
#include <system/domain/SurfaceParameters.h>
#include <system/optimization/AdaptiveDtSurfaceDistribution.h>
#include <system/optimization/SurfaceConstraint.h>
#include <system/optimization/GlobalSurfaceEnergy.h>
#include <system/optimization/ApproximateNeighborhoods.h>
#include <system/optimization/LODController.h>
#include <system/SFSystem.h>

#ifdef _WIN32
#pragma warning ( disable: 4305 4018 )
#endif

using namespace particle_sys;
using namespace custom_class;
using namespace std;

enum { COTAN, RADIAL };

//////////////////////////////////////////////////////////////////////////

//void SFSystem::optimize()
//{
//   cout << "About to optimize!!!" << endl;
//  if (_intersection[_current_intersection].ps->points().size() > 0) 
//     _intersection[_current_intersection].ps->optimize();
//}

//////////////////////////////////////////////////////////////////////////


//------------------------------------------------------------------------
// Function    : Constructor and Destructor
// Description : 
//------------------------------------------------------------------------
SFSystem::SFSystem( const char *param_file )
{
  // domain variables
  float max_allowed_sf;

  _selected_pt = 0.0;

  // now, pass all the variables in by reference to the file reader
  int energy;
  char i_file[350];
  readParamFile( param_file, _file_name, energy, max_allowed_sf,
                 _initial_sf,
                 _init_num_points, _num_materials, _num_intersections, 
                 i_file );

  // create the IO function
  _io_field = new IOScalarField( _file_name, IOScalarField::INTERPOLATING );
  _io_field->bounds( _start, _end );

  cout << "     DOMAIN : " << _start << " x " << _end << endl;

  // initialize the intersections
  _intersection = new IntersectionPackage[_num_intersections];

  readIntersectionFile( i_file );

  // create the intersection surfaces
  for ( int i = 0; i < _num_intersections; i++ )
  {
    // set the name of the mesh file (for the LOD controller)
    if ( !_intersection[i].init_with_ptcl )
    {
      if ( _num_materials != 1 )
        sprintf( _intersection[i].mesh_filename, "%s_%i", _file_name, 
                _intersection[i].materials[0]);
      else
        sprintf( _intersection[i].mesh_filename, "%s_%i", _file_name, i );
    }


    // instantiate the particle system
    _intersection[i].ps = new ParticleSystem();

    // create the boundary surfaces
    _intersection[i].surf = new Intersection( _io_field, 
      _intersection[i].num_materials, _intersection[i].materials );

    _intersection[i].max_sf = 
        min( max_allowed_sf, _intersection[i].surf->maxSF() ); 
    cout << "max sf: " << _intersection[i].max_sf << endl;

    // create the domain
    Domain *domain;
    if ( _intersection[i].num_materials == 1 )
    {
      Energy *e;
      switch ( energy )
      {
      case COTAN:
        e = new CotanEnergy( NORMALIZED_IDEAL_DISTANCE_6 );
        break;

      case RADIAL:
        e = new RadialEnergy( NORMALIZED_IDEAL_DISTANCE_6 );
        break;
      }
      domain = new Domain( e, _intersection[i].surf, 
                           _start, _end, 6, 
                           NORMALIZED_IDEAL_DISTANCE_6 );
    }
    else if ( _intersection[i].num_materials == 2 )
    {
      Energy *e;
      switch ( energy )
      {
      case COTAN:
        e = new CotanEnergy( NORMALIZED_IDEAL_DISTANCE_6 );
        break;

      case RADIAL:
        e = new RadialEnergy( NORMALIZED_IDEAL_DISTANCE_6 );
        break;
      }
      domain = new Domain( e, _intersection[i].surf, 
                           _start, _end, 6, 
                           NORMALIZED_IDEAL_DISTANCE_6 );
    }
    else if ( _intersection[i].num_materials == 3 )
    {
      Energy *e;
      switch ( energy )
      {
      case COTAN:
        e = new CotanEnergy( NORMALIZED_IDEAL_DISTANCE_2 );
        break;

      case RADIAL:
        e = new RadialEnergy( NORMALIZED_IDEAL_DISTANCE_2 );
        break;
      }
      domain = new Domain( e, _intersection[i].surf, 
                           _start, _end, 2,
                           NORMALIZED_IDEAL_DISTANCE_2 );
    }
    else if ( _intersection[i].num_materials == 4 )
    {
      Energy *e;
      switch ( energy )
      {
      case COTAN:
        e = new CotanEnergy( 0.0 );
        break;

      case RADIAL:
        e = new RadialEnergy( 0.0 );
        break;
      }
      domain = new Domain( e, _intersection[i].surf, 
                           _start, _end, 0,
                           NORMALIZED_IDEAL_DISTANCE_2 );
    }
//PIP
// decreased the precision by 1 to all, try relaxing the algo
    // optimizer variables
    float threshold = 10.0e-5;
    int modulo = 4;
    if ( _intersection[i].num_materials == 1 )
      threshold = 1.0e-5;
    else if ( _intersection[i].num_materials == 2 )
    {
      threshold = 1.0e-4;
      //modulo = 2;
    }	
    else if ( _intersection[i].num_materials == 3 )
      threshold = 1.0e-4;
    else if ( _intersection[i].num_materials == 4 )
    {
      threshold = 1.0e-3;
      modulo = 2;
    }
    
    // create the point data structure
    svector<DynamicSurfacePoint*> points;

    // get the size of the domain
    vector_type d_start, d_end;
    domain->domain( d_start, d_end );
    vector_type diff = d_end - d_start;

    // compute the initial radius and bin width
    float radius = _intersection[i].max_sf / 
      domain->normalizedIdealDistance();

    cout << "Binning structure radius = " << radius << endl;


    //******************************
    // create the neighborhood
    //******************************
    
    // find number of bins
    // if past memory threshold, multiply binning radius by appropriate factor to reduce # of bins
    // TODO: find a more versatile fix for this problem
    radius = radius*4;
    cout << "creating Neighborhood " << endl;
    Neighborhood<DynamicSurfacePoint>* n = 
      new neighborhood_type( radius, d_start, d_end );
    cout << "created Neighborhood " << endl;
    domain->neighborhood( n, radius );

    initializePointsWithMesh( _intersection[i].mesh_filename, points, 
                              _intersection[i].ps, domain, 
                              _intersection[i].max_sf,
                              _intersection[i].init_with_ptcl,
                              modulo );
    cout << "done initializePointsWithMesh" << endl;
    
    //******************************
    // populate the domain
    //******************************
    domain->populateDomain( points );
    cout << "done populateDomain" << endl;
    // optimizations
    Optimization **ops=NULL;
    ops = new Optimization*[4]; 
    int num_ops = 0;

    if ( _intersection[i].optimize_ptcl )
    {
      ops[num_ops++] = new ApproximateNeighborhoods( 6 );
      
      if ( _intersection[i].num_materials == 4 )
        ops[num_ops++] = new SurfaceConstraint(threshold);
      else
        ops[num_ops++] = new AdaptiveDtSurfaceDistribution(threshold);
        
      ops[num_ops++] = new GlobalSurfaceEnergyNA();
    }

cout << "DEBUG num_materials: " << _intersection[i].num_materials << endl;
    
    // create the optimizer
    Optimize *optimizer = new Optimize( ops, num_ops );
    cout << "created optimizer" << endl;
    // and initialize the system
    _intersection[i].ps->init( domain, optimizer, points );
    cout << "intersection init" << endl;
  }

  _current_intersection = -1;
  freezeIntersection();
  
  cout << "DEBUG point size: " <<  _intersection[_current_intersection].ps->points().size()<< endl;
}

SFSystem::~SFSystem()
{
  for ( int i = 0; i < _num_intersections; i++ )
  {
    delete _intersection[i].ps;
    delete [] _intersection[i].materials;
  }
  delete [] _intersection;
  delete _io_field;

  if ( _qj ) delete [] _qj;
  if ( _tj ) delete [] _tj;
  if ( _dj ) delete [] _dj;
}

//------------------------------------------------------------------------
// Function    : initializePointsWithMesh
// Description : 
//------------------------------------------------------------------------
void SFSystem::initializePointsWithMesh( const char *basename,
                                         svector<DynamicSurfacePoint*> 
                                         &points,
                                         ParticleSystem *ps, 
                                         Domain *domain,
                                         float max_surface_sf,
                                         bool init_with_ptcl,
                                         int modulo )
{
  cout << "Initializing with a mesh..." << endl;
  cout << "\t\t" << basename << endl;
  if ( !init_with_ptcl )
  {
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
        getline( in, buffer );      pos[2] = atof( buffer.c_str() );

        if ( !(counter % modulo) )
        {
          point = new DynamicSurfacePoint( domain, ps, pos );
          point->using_sf();
          point->max_sf( max_surface_sf );
          point->min_sf( 0.0 );

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
  else
  {
    // open the file
    ifstream in( basename, ifstream::in  );

    // check that it opened okay
    if ( !in )
    { 
      cout << "Error opening input file " << basename << endl;
      exit(1);
    }

    cout << "Reading particle file..." << endl;

    string buffer;
    in >> buffer;
    int num_particles = atoi( buffer.c_str() );
    points.resize( num_particles );

    cout << "  number of particles = " << num_particles << endl;

    DynamicSurfacePoint *point;
    vector_type pos, norm;
    float r;
    for ( int i = 0; i < num_particles; i++ )
    {
      in >> buffer; pos(0)  = atof( buffer.c_str() );
      in >> buffer; pos(1)  = atof( buffer.c_str() );
      in >> buffer; pos(2)  = atof( buffer.c_str() );
      in >> buffer; norm(0) = atof( buffer.c_str() );
      in >> buffer; norm(1) = atof( buffer.c_str() );
      in >> buffer; norm(2) = atof( buffer.c_str() );
      in >> buffer; r       = atof( buffer.c_str() );  

      point = new DynamicSurfacePoint( domain, ps, pos );
      point->using_sf();
      point->max_sf( max_surface_sf );
      point->min_sf( 0.0 );

      point->radius( r );
      point->normal( norm );

      // initialize the SurfaceParams
      point->updateSurfaceParameters();

      points[i] = point;
    }
    in.close();
  }
}

//------------------------------------------------------------------------
// Function    : cleanDoubleJunction
// Description : 
//------------------------------------------------------------------------
void SFSystem::cleanDoubleJunction()
{
  // make sure we are on a double junction
  if ( _intersection[_current_intersection].num_materials != 2 )
  {
    cout << "Current particle system is not sampling a double junction!"
      << endl;
    return;
  }

  int num_deleted = 0;

  // now, go through the quad and triple-junctions and find out if
  //   this double junction has any points on a quad or triple-junction
  float threshold = 1.0e-4;
  SurfacePointParams params;
  int num_ptcl;
  vec<3> pos;
  for ( int i = 0; i < _num_qj; i++ )
  {
    num_ptcl = _intersection[_current_intersection].ps->numParticles();
    for ( int p = 0; p < num_ptcl; p++ )
    {
      pos = _intersection[_current_intersection].ps->point(p)->position();
      _intersection[_qj[i]].surf->computeSurfacePointParams(
        pos, params, false );

      /*if ( params._F < threshold ) 
      {
        ++num_deleted;
        _intersection[_current_intersection].ps->removePoint( p );
        --p;
      }*/
    }
  }

  for ( int i = 0; i < _num_tj; i++ )
  {
    num_ptcl = _intersection[_current_intersection].ps->numParticles();
    for ( int p = 0; p < num_ptcl; p++ )
    {
      pos = _intersection[_current_intersection].ps->point(p)->position();
      _intersection[_tj[i]].surf->computeSurfacePointParams(
        pos, params, false );

      if ( params._F < threshold ) 
      {
        ++num_deleted;
        _intersection[_current_intersection].ps->removePoint( p );
        --p;
        --num_ptcl;
      }
    }
  }
  
  cout << "Number of particles removed = " << num_deleted << endl;
}

//------------------------------------------------------------------------
// Function    : freezeIntersection
// Description : 
//------------------------------------------------------------------------
void SFSystem::freezeIntersection()
{
  //if ( _current_intersection >= 0 )
  //  cleanDoubleJunction();

  // check if we are already at the last intersection
  if ( _current_intersection == (_num_intersections-1) )
    return;

  ++_current_intersection;

  // add the points from this system and all previous into the 
  //   neighborhood structure of this system
  for ( int i = (_current_intersection-1); i >=0; i-- )
  {
    // for quad junctions, don't add any 
    if ( _intersection[_current_intersection].num_materials == 4 )
      continue;

    // for triple junctions, only add the quad junctions
    if ( _intersection[_current_intersection].num_materials == 3 )
    {
      if ( _intersection[i].num_materials == 4 )
        _intersection[_current_intersection].ps->domain()->populateDomain( 
          _intersection[i].ps->points() );
      continue;
    }

    // for double junctions, only add the quad and triple junctions
    if ( (_intersection[i].num_materials == 3) || 
         (_intersection[i].num_materials == 4) )
      _intersection[_current_intersection].ps->domain()->populateDomain( 
        _intersection[i].ps->points() );
  }

  //cleanDoubleJunction();
}

//------------------------------------------------------------------------
// Function    : print
// Description : 
//------------------------------------------------------------------------
void SFSystem::print() const
{
  for ( int i = 0; i <= _current_intersection; i++ )
  {
    _intersection[i].ps->print();
    cout << endl;
  }
}

//------------------------------------------------------------------------
// Function    : writePointFile
// Description : 
//------------------------------------------------------------------------
void SFSystem::writePointFile(int i) const
{
  char filename[300];

  // write out the points for triangulating
  vector_type pos;
//  for ( int i = 0; i < _num_intersections; i++ )
  {
    sprintf( filename, "%s_%i.pts", _file_name, i );

    FILE *file = fopen( filename, "w" );
    if ( !file )
    {
      cout << "Points File Could Not Be Opened...\n";
      return;
    }

    cout << "Writing point file " << i << " ..." << endl;

    for ( int j = 0; j < _intersection[i].ps->numParticles(); j++ )
    {
      pos = _intersection[i].ps->point(j)->position();
#ifdef TWO_D
      fprintf( file, "%f %f\n", pos(0), pos(1) );
#else
      fprintf( file, "%f %f %f\n", pos(0), pos(1), pos(2) );
#endif
    }

    fclose( file );
    cout << "       DONE." << endl;
  }
  
  // write out the particle files
  vector_type n;
  float r;
//  for ( int i = 0; i < _num_intersections; i++ )
  {
    sprintf( filename, "%s_%i.ptcl", _file_name, i );

    FILE *file = fopen( filename, "w" );
    if ( !file )
    {
      cout << "Particle File Could Not Be Opened...\n";
      return;
    }

    cout << "Writing particle file " << i << " ..." << endl;

    fprintf( file, "%i\n", _intersection[i].ps->numParticles() );

    for ( int j = 0; j < _intersection[i].ps->numParticles(); j++ )
    {
      pos = _intersection[i].ps->point(j)->position();
      n = _intersection[i].ps->point(j)->normal();
      r = _intersection[i].ps->point(j)->radius();
#ifdef TWO_D
      fprintf( file, "%f %f %f %f %f\n", pos(0), pos(1), n(0), n(1), r );
#else
      fprintf( file, "%f %f %f %f %f %f %f\n", pos(0), pos(1), pos(2), 
        n(0), n(1), n(2), r );
#endif
    }
    fclose( file );
    cout << "       DONE." << endl;
  }
}



// trim from start
static inline std::string &ltrim(std::string &s) {
	s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(::isspace))));
	return s;
}

// trim from end
static inline std::string &rtrim(std::string &s) {
	s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(::isspace))).base(), s.end());
	return s;
}

// trim from both ends
static inline std::string &trim(std::string &s) {
	return ltrim(rtrim(s));
}

//------------------------------------------------------------------------
// Function    : readParamFile
// Description : read in the parameter file, and set the Domain 
//               variables
//------------------------------------------------------------------------
void SFSystem::readParamFile( const char *param_file, char *s, 
                              int &e, float &max_sf, 
                              float &initial_sf,
                              int &init_num_pts, int &num_surfaces,
                              int &num_intersections,
                              char *i_file )
{
  // read in the file
  FILE *file = fopen( param_file, "r" );
  if ( !file )
  {
    cout << "Parameter File Could Not Be Opened...\n";
    return;
  }

  cout << "READING in parameter file " << param_file << endl;

  char v[300], value_s[300];

  // ENERGY
  fscanf( file, "%s%s", v, value_s );
  if ( !strncmp( value_s, "cotan", 5 ) )
  {
    e = COTAN;
    cout << "     ENERGY : cotan\n";
  }
  else if ( !strncmp( value_s, "radial", 6 ) )
  {
    e = RADIAL;
    cout << "     ENERGY : radial\n";
  }

  fscanf( file, "%s%i", v, &num_surfaces );
  cout << "     NUMBER OF SURFACES : " << num_surfaces << endl;

  fscanf( file, "%s%i", v, &num_intersections );
  cout << "     NUMBER OF INTERSECTIONS : " << num_intersections << endl;

  // intersection file
  //fscanf( file, "%s%s", v, value_s );
  fscanf( file, "%s", v );
  fgets( value_s, 300, file );
	std::string value_str(value_s);
	value_str = trim(value_str);
  sprintf( i_file, "%s", value_str.c_str() );

  // filenames file
//	fscanf( file, "%s%s", v, value_s );
  fscanf( file, "%s", v );
  fgets( value_s, 300, file );
	value_str = std::string(value_s);
	value_str = trim(value_str);
  sprintf( s, "%s", value_str.c_str() );
  
  // SFSystem INITIALIZATION
  fscanf( file, "%s%i", v, &init_num_pts );
  cout << "     INITIAL NUMBER OF POINTS : " << init_num_pts << endl;

  fscanf( file, "%s%f", v, &max_sf );
  cout << "     MAX ALLOWED SF : " << max_sf << endl;

  fscanf( file, "%s%f", v, &initial_sf );
  cout << "     INITIAL SF : " << initial_sf << endl;

  fclose( file );
}

//------------------------------------------------------------------------
// Function    : readIntersectionFile
// Description : 
//------------------------------------------------------------------------
void SFSystem::readIntersectionFile( const char *i_file )
{
  // read in the file
  FILE *file = fopen( i_file, "r" );
  if ( !file )
  {
    cout << "Intersection File: " << i_file 
         << " Could Not Be Opened..." << endl;
    return;
  }

  cout << "READING in intersection file " << i_file << endl;

  // initialize the intersection counter
  int c = 0;

  // num quad junctions
  fscanf( file, "%i\n", &_num_qj );
  cout << "   num quad junctions = " << _num_qj << endl;

  char yesno;
  if ( _num_qj )
  {
    _qj = new int[_num_qj];
    int m1, m2, m3, m4;
    for ( int i = 0; i < _num_qj; i++ )
    {
      fscanf( file, "%i %i %i %i %c", &m1, &m2, &m3, &m4, &yesno );
      _qj[i] = c;
      _intersection[c].num_materials = 4;
      _intersection[c].materials = new int[_intersection[c].num_materials];
      _intersection[c].materials[0] = m1;
      _intersection[c].materials[1] = m2;
      _intersection[c].materials[2] = m3;
      _intersection[c].materials[3] = m4;
      _intersection[c].color.set( 0.1f, 0.1f, 0.1f );

      if ( yesno == 'y' )
      {
        _intersection[c].init_with_ptcl = true;
        fscanf( file, "%s %c\n", _intersection[c].mesh_filename, &yesno );
        if ( yesno == 'y' )
          _intersection[c].optimize_ptcl = true;
        else
          _intersection[c].optimize_ptcl = false;
      }
      else
      {
        _intersection[c].init_with_ptcl = false;
        _intersection[c].optimize_ptcl = true;
        fscanf( file, "%\n" );
      }
      
      c++;
    }
  }
  else
    _qj = NULL;

  // num triple junctions
  fscanf( file, "%i\n", &_num_tj );
  cout << "   num triple junctions = " << _num_tj << endl;

  if ( _num_tj )
  {
    _tj = new int[_num_tj];
    int m1, m2, m3;
    for ( int i = 0; i < _num_tj; i++ )
    {
      fscanf( file, "%i %i %i %c", &m1, &m2, &m3, &yesno );
      _tj[i] = c;
      _intersection[c].num_materials = 3;
      _intersection[c].materials = new int[_intersection[c].num_materials];
      _intersection[c].materials[0] = m1;
      _intersection[c].materials[1] = m2;
      _intersection[c].materials[2] = m3;

      _intersection[c].color.random();
      _intersection[c].color.normalize();
      _intersection[c].color *= 1.5f;

      if ( yesno == 'y' )
      {
        _intersection[c].init_with_ptcl = true;
        fscanf( file, "%s %c\n", _intersection[c].mesh_filename, &yesno );
        if ( yesno == 'y' )
          _intersection[c].optimize_ptcl = true;
        else
          _intersection[c].optimize_ptcl = false;
      }
      else
      {
        _intersection[c].init_with_ptcl = false;
        _intersection[c].optimize_ptcl = true;
        fscanf( file, "%\n" );
      }

      c++;
    }
  }
  else
    _tj = NULL;

  // num double junctions
  fscanf( file, "%i\n", &_num_dj );
  cout << "   num double junctions = " << _num_dj << endl;

  if ( _num_dj )
  {
    _dj = new int[_num_dj];
    int m1, m2;
    for ( int i = 0; i < _num_dj; i++ )
    {
      fscanf( file, "%i %i %c", &m1, &m2, &yesno );
      _dj[i] = c;
      _intersection[c].num_materials = 2;
      _intersection[c].materials = new int[_intersection[c].num_materials];
      _intersection[c].materials[0] = m1;
      _intersection[c].materials[1] = m2;
      
      _intersection[c].color.random();
      _intersection[c].color.random();
      _intersection[c].color.normalize();
      _intersection[c].color *= 1.5f;

      if ( yesno == 'y' )
      {
        _intersection[c].init_with_ptcl = true;
        fscanf( file, "%s %c\n", _intersection[c].mesh_filename, &yesno );
        if ( yesno == 'y' )
          _intersection[c].optimize_ptcl = true;
        else
          _intersection[c].optimize_ptcl = false;
      }
      else
      {
        _intersection[c].init_with_ptcl = false;
        _intersection[c].optimize_ptcl = true;
        fscanf( file, "%\n" );
      }

      c++;
    }
  }
  else
    _dj = NULL;


  fclose( file );
}

#if 0
//------------------------------------------------------------------------
// select: 
//------------------------------------------------------------------------
#define BUFSIZE 512
void SFSystem::selectParticle( int x, int y )
{
  float min_z=MAX_VALUE;
  for ( int i = 0; i <= _current_intersection; i++ )
  {
    GLint  viewport[4];
    glGetIntegerv( GL_VIEWPORT, viewport );

    GLuint selectBuf[BUFSIZE];
    glSelectBuffer( BUFSIZE, selectBuf );
    glRenderMode( GL_SELECT );

    glInitNames();
    glPushName( 0 );

    // save the current gl states
    glPushAttrib( GL_ENABLE_BIT | GL_LIGHTING_BIT | GL_TRANSFORM_BIT );
    glDisable( GL_LIGHTING );

    GLdouble projection[16];
    glGetDoublev( GL_PROJECTION_MATRIX, projection );

    glMatrixMode( GL_PROJECTION );
    glPushMatrix();
    glLoadIdentity();

    // need to select the picking region (in window coordinates!)
    gluPickMatrix( (GLdouble)x, (GLdouble)( viewport[3] - y ), 2.0, 2.0,
                    viewport );

    glMultMatrixd( projection );
    float s = 1.0;
    vector_type pos;
    for ( int j = 0; j < _intersection[i].ps->numParticles(); j++ )
    {
      pos = _intersection[i].ps->point(j)->position();
      glLoadName( j );

      glBegin( GL_TRIANGLES );
      {
        glVertex3f( pos(0)-0.25, pos(1)-0.25, pos(2) );
        glVertex3f( pos(0)+0.25, pos(1)-0.25, pos(2) );
        glVertex3f( pos(0),      pos(1)+0.25, pos(2) );
      }
      glEnd();
    }

    glPopName();

    glMatrixMode( GL_PROJECTION );
    glPopMatrix();
    glFlush();

    GLint hits = glRenderMode( GL_RENDER );

    float mz;
    int mi;
    if ( hits == 0 )
      cout << endl << "NO HITS -- please select again..." << endl;
    else if ( hits == -1 )
      cout << endl<< "OpenGL selection error -- please select again...";
    else
    {
      processHits( hits, selectBuf, mz, mi );
      if ( mz < min_z )
      {
        min_z = mz;
        _selected_pt = _intersection[i].ps->point(mi)->position();

        cout << min_z << " " << mi << endl;
      }
    }
    glPopAttrib();  
  }

  cout << "Selected position: " << _selected_pt << endl;

  // print out some info
  _io_field->printFunctionValues( _selected_pt );

}

void SFSystem::processHits( GLint hits, GLuint buffer[],
                            float &min_z, int &min_index )
{
  cout << "Processing hits..." << endl;
  float *min_z_vals = new float[hits];
  int   *names = new int[hits];
  GLuint *ptr = (GLuint*)buffer;

  for ( int i = 0; i < hits; i++ )
  {
    ++ptr;  // number of names in hit

    min_z_vals[i] = (float)(*ptr / pow(2.0,32)); ++ptr;
    ++ptr;  // max z value

    names[i] = *ptr;
    ptr++;
  }

  // now, determine which primitive got hit first (the closest)
  int smallest_z = 0;
  for (int i = 1; i < hits; i++ )
  {
    if ( min_z_vals[smallest_z] > min_z_vals[i] )
      smallest_z = i;
  }

  min_z = min_z_vals[smallest_z];
  min_index = names[smallest_z];

  // clean up memory
  delete [] names;
  delete [] min_z_vals;
}

void SFSystem::renderSelectedParticle()
{
  // save the current gl states
  glPushAttrib( GL_ENABLE_BIT | GL_LIGHTING_BIT | GL_TRANSFORM_BIT );
  glDisable( GL_LIGHTING );

  glColor3f( 1.0, 0.0, 0.0 );
 
  float offset = 0.25;
  glPushMatrix();
  {
    glTranslatef( _selected_pt(0), 
                  _selected_pt(1), 
                  _selected_pt(2) );
    glutSolidCube( 2.0*offset );

  }
  glPopMatrix();

  glPopAttrib();  
}
#endif
