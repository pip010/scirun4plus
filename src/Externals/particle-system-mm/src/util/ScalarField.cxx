#include <iostream>
#include <cstdlib>
#include <ScalarField.h>
#include <NrrdIO.h>

#include <string>

using namespace std;

//------------------------------------------------------------------------
// Function    : Constructor and Destructor
// Description : initializing constructor
//------------------------------------------------------------------------
ScalarField::ScalarField( const char *filename, int k,
                          float sx, float sy, float sz ) : Surface()
{
  std::string filenameString(filename);
  std::string::size_type pos = filenameString.rfind(".nrrd");
  _centers[0] = UNKNOWN;
  _centers[1] = UNKNOWN;
  _centers[2] = UNKNOWN;
  
  if (pos == std::string::npos)
  {
    _scale[0] = sx;
    _scale[1] = sy;
    _scale[2] = sz;

    // read a raw binary file
    readFile( filename );
  }
  else
  {
    // If reading nrrd, ignore spacing arguments
    // and use nrrd spacing instead.
    // Could there be a more elegant way to handle this?
    NrrdIO::read_nrrd(filename, _field,
                      _scale[0], _scale[1], _scale[2],
                      _xdim, _ydim, _zdim,
                      _centers[0], _centers[1], _centers[2]);
  }

  switch ( k )
  {
    case BSPLINE:
    _kernel  = new CubicBSpline();
    _kernelD = new CubicBSplineD();
    _kernelDD = new CubicBSplineDD();
    break;

    case CATMULLROM:
    _kernel  = new CatmullRom();
    _kernelD = new CatmullRomD();
    _kernelDD = new CatmullRomDD();
    break;

    case LINEAR:
    _kernel = NULL;
    _kernelD = NULL;
    _kernelDD = NULL;

    default:
    _kernel  = new CubicBSpline();
    _kernelD = new CubicBSplineD();
    _kernelDD = new CubicBSplineDD();
    break;
  }

  _start = 0.0;
  _end.set( _scale[0]*(float)(_xdim-1),
            _scale[1]*(float)(_ydim-1),
            _scale[2]*(float)(_zdim-1) );
}

ScalarField::~ScalarField()
{
  if ( _kernel ) delete _kernel;
  if ( _kernelD ) delete _kernelD;
  if ( _kernelDD ) delete _kernelDD;
}

//------------------------------------------------------------------------
// Function    : computeSurfacePointParams()
// Description : Use to compute/set F Fx H at the particle pos.
//------------------------------------------------------------------------
bool ScalarField::computeSurfacePointParams( const vector_type &pos, 
                                             SurfacePointParams &params,
                                             bool computeGradient,
                                             bool computeHessian )
                              
{
#ifdef THREE_D
  // get the lower left corner of the bounding voxel
  //   NOTE: if the position is exactly on the high end boundary, make 
  //         sure to not run outside of the domain!
  int llc_x = min( (int)floor(pos(0)/_scale[0]), _xdim-2 );
  int llc_y = min( (int)floor(pos(1)/_scale[1]), _ydim-2 );
  int llc_z = min( (int)floor(pos(2)/_scale[2]), _zdim-2 );

  // first, do the interpolation in the front xy slab
  float t_x = (pos(0)/_scale[0] - (float)llc_x);
  float t_y = (pos(1)/_scale[1] - (float)llc_y);
  float t_z = (pos(2)/_scale[2] - (float)llc_z);

  if ( !_kernel )
  {
    float b_interp = interpolate(t_x, _field(llc_x,llc_y,llc_z),
                                 _field(llc_x+1,llc_y,llc_z));
                                 
    float t_interp = interpolate(t_x, _field(llc_x,llc_y+1,llc_z),
                                 _field(llc_x+1,llc_y+1,llc_z));

    float f_interp = interpolate(t_y, b_interp, t_interp);

    // second, do the interpolation in the back xy slab
    b_interp = interpolate(t_x, _field(llc_x,llc_y,llc_z+1),
                           _field(llc_x+1,llc_y,llc_z+1));
    t_interp = interpolate(t_x, _field(llc_x,llc_y+1,llc_z+1), 
                           _field(llc_x+1,llc_y+1,llc_z+1));

    float bk_interp = interpolate(t_y, b_interp, t_interp);

    params._F = interpolate(t_z, f_interp, bk_interp);
    
    return true;
  }

  vec<4> tau_x( t_x*t_x*t_x, t_x*t_x, t_x, 1 );
  vec<4> tau_y( t_y*t_y*t_y, t_y*t_y, t_y, 1 );
  vec<4> tau_z( t_z*t_z*t_z, t_z*t_z, t_z, 1 );

  vec<3> tau_xD( tau_x(1), tau_x(2), tau_x(3) );
  vec<3> tau_yD( tau_y(1), tau_y(2), tau_y(3) );
  vec<3> tau_zD( tau_z(1), tau_z(2), tau_z(3) );

  vec<2> tau_xDD( tau_x(2), tau_x(3) );
  vec<2> tau_yDD( tau_y(2), tau_y(3) );
  vec<2> tau_zDD( tau_z(2), tau_z(3) );

  // assuming a 4^3 support
  vec<4> x, y, z, yDx, yDy, yDz, zDx, zDy, zDz,
    yDxx, yDxy, yDxz, yDyy, yDyz, yDzz,
    zDxx, zDxy, zDxz, zDyy, zDyz, zDzz;
  for ( int slice = -1; slice <= 2; slice++ )
  {
    for ( int row = -1; row <= 2; row++ )
    {
      // set the x values
      for ( int element = -1; element <= 2; element++ )
      {
        x[element+1] = _field( min(max(llc_x+element,0),(_xdim-1)), 
                               min(max(llc_y+row,0),(_ydim-1)), 
                               min(max(llc_z+slice,0),(_zdim-1)) );
      }
      y[row+1]   = _kernel->w( x, tau_x );

      if ( computeGradient || computeHessian )
      {
        yDx[row+1] = _kernelD->w( x, tau_xD );
        yDy[row+1] = y[row+1];
        yDz[row+1] = y[row+1];
      }

      if ( computeHessian )
      {
        yDxx[row+1] = _kernelDD->w( x, tau_xDD );
        yDxy[row+1] = yDx[row+1];
        yDxz[row+1] = yDx[row+1];
        yDyy[row+1] = yDy[row+1];
        yDyz[row+1] = yDy[row+1];
        yDzz[row+1] = yDy[row+1];
      }
    }
    z[slice+1]  = _kernel->w( y, tau_y );

    if ( computeGradient )
    {
      zDx[slice+1] = _kernel->w( yDx, tau_y );
      zDy[slice+1] = _kernelD->w( yDy, tau_yD );
      zDz[slice+1] = _kernel->w( yDz, tau_y );
    }

    if ( computeHessian )
    {
      zDxx[slice+1] = _kernel->w( yDxx, tau_y );
      zDxy[slice+1] = _kernelD->w( yDxy, tau_yD );
      zDxz[slice+1] = _kernel->w( yDxz, tau_y );
      zDyy[slice+1] = _kernelDD->w( yDyy, tau_yDD );
      zDyz[slice+1] = _kernelD->w( yDyz, tau_yD );
      zDzz[slice+1] = _kernel->w( yDzz, tau_y );
    }
  }

  params._F = _kernel->w( z, tau_z );

  if ( computeGradient )
  {
    params._Fx(0) = _kernel->w( zDx, tau_z );
    params._Fx(1) = _kernel->w( zDy, tau_z );
    params._Fx(2) = _kernelD->w( zDz, tau_zD );
  }

  if ( computeHessian )
  {
    params._Fxx(0,0) = _kernel->w( zDxx, tau_z );
    params._Fxx(0,1) = _kernel->w( zDxy, tau_z );
    params._Fxx(0,2) = _kernelD->w( zDxz, tau_zD );
    params._Fxx(1,1) = _kernel->w( zDyy, tau_z );
    params._Fxx(1,2) = _kernelD->w( zDyz, tau_zD );
    params._Fxx(2,2) = _kernelDD->w( zDzz, tau_zDD );

    params._Fxx(1,0) = params._Fxx(0,1);
    params._Fxx(2,0) = params._Fxx(0,2);
    params._Fxx(2,1) = params._Fxx(1,2);
  }
#else
  cout << "ScalarField: PROBLEM, haven't implemented 2D!!!!" << endl;
#endif

  return true;
}

//------------------------------------------------------------------------
// Function    : readDataFile()
// Description : 
//------------------------------------------------------------------------
void ScalarField::readFile( const char *filename )
{
#ifdef THREE_D
  cout << "Reading in the data file ..." << endl;

  // first open one file to get the dimensions
  FILE *in_file_field = fopen( filename, "r" );
  if ( !in_file_field )
  {
    cout << "ScalarField::Error reading field file " << filename << endl;
    exit( 1 );
  }

  char tmp[100];
  fgets( tmp, 100, in_file_field );
  fscanf( in_file_field, " %i %i %i", &_xdim, &_ydim, &_zdim );

  fclose( in_file_field );

  // resize the field
  _field.resize( _xdim, _ydim, _zdim );

  // open the files
  in_file_field = fopen( filename, "rb" );
  if ( !in_file_field )
  {
    cout << "ScalarField::Error reading field file." << endl;
    exit( 1 );
  }
  
  // read in and discard the headers
  for ( int i = 0; i < 2; i++ )
    fgets( tmp, 100, in_file_field );

  // read in the data values
  for ( int k = 0; k < _zdim; k++ )  
    for ( int j = 0; j < _ydim; j++ )
      for ( int i = 0; i < _xdim; i++ )
        // the distance field
        fread( &_field(i,j,k), sizeof(float), 1, 
               in_file_field );

  // close the files
  fclose( in_file_field );
  cout << "              DONE" << endl;
#else 
  cout << "ScalarField: PROBLEM, haven't implemented 2D!!!!" << endl;
#endif
}

