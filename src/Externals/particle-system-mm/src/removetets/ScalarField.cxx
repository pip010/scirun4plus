#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <NrrdIO.h>
#include "ScalarField.h"


#define SF_VALUE 0.25

using namespace std;

//------------------------------------------------------------------------
// Function    : Constructor and Destructor
// Description : initializing constructor
//------------------------------------------------------------------------
ScalarField::ScalarField( const char *df_file_base,
                          float isovalue, int kernel,
                          float scale_x, 
                          float scale_y, 
                          float scale_z) : Surface()                         
{
  std::string filenameString(df_file_base);
  std::string::size_type pos = filenameString.rfind(".nrrd");
  if (pos == std::string::npos)
  {

    // axis scaling
    _scale[0] = scale_x;
    _scale[1] = scale_y;
    _scale[2] = scale_z;

    _isovalue = isovalue;

    readFile( df_file_base );
    }
  else
  {
    // If reading nrrd, ignore spacing arguments
    // and use nrrd spacing instead.
    // Could there be a more elegant way to handle this?

    _isovalue = isovalue;
    
    readNrrdFile( df_file_base );
    
  }


  switch ( kernel )
  {
    case BSPLINE:
    _kernel  = new CubicBSpline();
    _kernelD = new CubicBSplineD();
    break;

    case CATMULLROM:
    _kernel  = new CatmullRom();
    _kernelD = new CatmullRomD();
    break;

    default:
    _kernel  = new CubicBSpline();
    _kernelD = new CubicBSplineD();
    break;
  }  
}

ScalarField::~ScalarField()
{
  delete _kernel;
  delete _kernelD;
}

//------------------------------------------------------------------------
// Function    : computeSurfacePointParams()
// Description : Use to compute/set F Fx H at the particle pos.
//------------------------------------------------------------------------
bool ScalarField::computeSurfacePointParams( const vector_type &pos, 
                                             SurfacePointParams &params,
                                             bool computeGradient,
                                             bool computeHessian ) const
                              
{
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

  vec<4> tau_x( t_x*t_x*t_x, t_x*t_x, t_x, 1 );
  vec<4> tau_y( t_y*t_y*t_y, t_y*t_y, t_y, 1 );
  vec<4> tau_z( t_z*t_z*t_z, t_z*t_z, t_z, 1 );

  vec<3> tau_xD( tau_x(1), tau_x(2), tau_x(3) );
  vec<3> tau_yD( tau_y(1), tau_y(2), tau_y(3) );
  vec<3> tau_zD( tau_z(1), tau_z(2), tau_z(3) );

  // assuming a 4^3 support
  vec<4> x, y, z, yDx, yDy, yDz, zDx, zDy, zDz;
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
    }
    z[slice+1]  = _kernel->w( y, tau_y );

    if ( computeGradient || computeHessian )
    {
      zDx[slice+1] = _kernel->w( yDx, tau_y );
      zDy[slice+1] = _kernelD->w( yDy, tau_yD );
      zDz[slice+1] = _kernel->w( yDz, tau_y );
    }
  }

  params._F = _kernel->w( z, tau_z ) - _isovalue;

  if( computeGradient || computeHessian )
  {
    params._Fx(0) = _kernel->w( zDx, tau_z );
    params._Fx(1) = _kernel->w( zDy, tau_z );
    params._Fx(2) = _kernelD->w( zDz, tau_zD );
  }

  return true;
}

//------------------------------------------------------------------------
// Function    : readDataFile()
// Description : 
//------------------------------------------------------------------------
void ScalarField::readFile( const char *df_file_base,
                            int num_header_lines )
{
  cout << "Reading in the .vol data file ..." << endl;

  // create the file names
  char df_file[150]; sprintf( df_file, "%s", df_file_base );

  // first open one file to get the dimensions
  FILE *in_file_field = fopen( df_file, "r" );
  if ( !in_file_field )
  {
    cout << "ScalarField::Error reading field file " << df_file << endl;
    exit( 1 );
  }

  char tmp[100];
  fgets( tmp, 100, in_file_field );
  _d_start = 0.0;
  fscanf( in_file_field, " %i %i %i", &_xdim, &_ydim, &_zdim );
  _d_end.set( _scale[0]*(float)(_xdim-1), 
              _scale[1]*(float)(_ydim-1), 
              _scale[2]*(float)(_zdim-1) );
  fclose( in_file_field );

  cout<<"scale: "<< _scale[0]<<" "<< _scale[1]<<" "<< _scale[2]<<endl;
  cout<<"dims: "<<_xdim<<" "<< _ydim<<" "<< _zdim<<endl;
  
  // resize the field
  _field.resize( _xdim, _ydim, _zdim );

  // open the files
  in_file_field = fopen( df_file, "rb" );
  if ( !in_file_field )
  {
    cout << "ScalarField::Error reading field file." << endl;
    exit( 1 );
  }

  // read in and discard the headers
  for ( int i = 0; i < 2; i++ )
  {
    fgets( tmp, 100, in_file_field );
  }

  // read in the data values
  for ( int k = 0; k < _zdim; k++ )  
    for ( int j = 0; j < _ydim; j++ )
      for ( int i = 0; i < _xdim; i++ )
      {
        // the distance field
        fread( &_field(i,j,k), sizeof(float), 1, 
               in_file_field );
        if ( _field(i,j,k) < _min_isovalue )
          _min_isovalue = _field(i,j,k);
        if ( _field(i,j,k) > _max_isovalue )
          _max_isovalue = _field(i,j,k);
      }


  // close the files
  fclose( in_file_field );

  cout << "              DONE" << endl;
}



//------------------------------------------------------------------------
// Function    : readNrrdDataFile()
// Description : 
//------------------------------------------------------------------------
void ScalarField::readNrrdFile( const char *df_file_base)
{
  cout << "Reading in the .nrrd data file ..." << endl;

  // create the file names
  char df_file[150]; sprintf( df_file, "%s", df_file_base );

  
  // first open one file to get the dimensions
  NrrdIO::read_nrrd(df_file, _field, _scale[0], _scale[1], _scale[2], _xdim, _ydim, _zdim);
  cout<<"scale: "<< _scale[0]<<" "<< _scale[1]<<" "<< _scale[2]<<endl;
  cout<<"dims: "<<_xdim<<" "<< _ydim<<" "<< _zdim<<endl;
  
  _d_start = 0.0;
  _d_end.set( _scale[0]*(float)(_xdim-1), 
              _scale[1]*(float)(_ydim-1), 
              _scale[2]*(float)(_zdim-1) );



  

  
  // check the data values
  for ( int k = 0; k < _zdim; k++ ) 
    for ( int j = 0; j < _ydim; j++ )
      for ( int i = 0; i < _xdim; i++ )
      {
        if ( _field(i,j,k) < _min_isovalue )
          _min_isovalue = _field(i,j,k);
        if ( _field(i,j,k) > _max_isovalue )
          _max_isovalue = _field(i,j,k);
      }

   


  cout << "              DONE" << endl;
}


