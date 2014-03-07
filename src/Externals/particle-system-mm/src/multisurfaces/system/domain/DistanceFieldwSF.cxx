#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>

typedef unsigned char byte;

#include <system/domain/DistanceFieldwSF.h>

using namespace std;
using namespace particle_sys;

//------------------------------------------------------------------------
// Function    : Constructor and Destructor
// Description : 
//------------------------------------------------------------------------
DistanceFieldwSF::DistanceFieldwSF( const char *df_file_base,
                                    float &max_sf,
                                    int &xdim, int &ydim, int &zdim) : Surface()
{
  _min_sf = 1.0e6;
  _max_sf = 0.0;

  readDistanceFieldFile( df_file_base );

  _isovalue = 0.0;

  max_sf = _max_sf;
  xdim = _xdim; ydim = _ydim; zdim = _zdim;
}

//------------------------------------------------------------------------
// Function    : computeSurfacePointParams()
// Description : 
//               NOTE --> shouldn't break this up any more! from a 
//                        software engineering perspective, we'd want
//                        to reorder the fields, so that we could just
//                        pass a pointer to a voxel and interpolate in
//                        a single function. but, then this would be
//                        a horrible memory lay out!!!!!
//------------------------------------------------------------------------
bool DistanceFieldwSF::computeSurfacePointParams( const vector_type &pos, 
                                                  SurfacePointParams 
                                                  &params,
                                                  bool computeHessian ) const
{  
  // get the lower left corner of the bounding voxel
  //   NOTE: if the position is exactly on the high end boundary, make 
  //         sure to not run outside of the domain!
  int llc_x = min( (int)floor(pos(0)), _xdim-2 );
  int llc_y = min( (int)floor(pos(1)), _ydim-2 );
#ifdef THREE_D
  int llc_z = min( (int)floor(pos(2)), _zdim-2 );
#else
  int llc_z = 0;
#endif

  // first, do the interpolation in the front xy slab
  float x_offset = pos(0) - (float)llc_x;
  float y_offset = pos(1) - (float)llc_y;
#ifdef THREE_D
  float z_offset = pos(2) - (float)llc_z; 
#endif

  // the field
  float b_interp = interpolate(x_offset, _field(llc_x,llc_y,llc_z),
                                         _field(llc_x+1,llc_y,llc_z));
                                 
  float t_interp = interpolate(x_offset, _field(llc_x,llc_y+1,llc_z),
                                         _field(llc_x+1,llc_y+1,llc_z));

#ifdef THREE_D
  float f_interp = interpolate(y_offset, b_interp, t_interp);

  // second, do the interpolation in the back xy slab
  b_interp = interpolate(x_offset, _field(llc_x,llc_y,llc_z+1),
                                   _field(llc_x+1,llc_y,llc_z+1));
  t_interp = interpolate(x_offset, _field(llc_x,llc_y+1,llc_z+1), 
                                   _field(llc_x+1,llc_y+1,llc_z+1));

  float bk_interp = interpolate(y_offset, b_interp, t_interp);

  params._F = interpolate(z_offset, f_interp, bk_interp) - _isovalue;
#else
  params._F = interpolate(y_offset, b_interp, t_interp) - _isovalue;
#endif


  // sizing field
  b_interp = interpolate(x_offset, _sf(llc_x,llc_y,llc_z),
                                   _sf(llc_x+1,llc_y,llc_z));
                                 
  t_interp = interpolate(x_offset, _sf(llc_x,llc_y+1,llc_z),
                                   _sf(llc_x+1,llc_y+1,llc_z));

#ifdef THREE_D
  f_interp = interpolate(y_offset, b_interp, t_interp);

  // second, do the interpolation in the back xy slab
  b_interp = interpolate(x_offset, _sf(llc_x,llc_y,llc_z+1),
                                   _sf(llc_x+1,llc_y,llc_z+1));
  t_interp = interpolate(x_offset, _sf(llc_x,llc_y+1,llc_z+1), 
                                   _sf(llc_x+1,llc_y+1,llc_z+1));

  bk_interp = interpolate(y_offset, b_interp, t_interp);

  params._sf = interpolate(z_offset, f_interp, bk_interp);
#else
  params._sf = interpolate(y_offset, b_interp, t_interp);
#endif



  // first derivative in the x direction
  b_interp = interpolate(x_offset, _field_d(llc_x,llc_y,llc_z,0),
                                   _field_d(llc_x+1,llc_y,llc_z,0));
  t_interp = interpolate(x_offset, _field_d(llc_x,llc_y+1,llc_z,0),
                                   _field_d(llc_x+1,llc_y+1,llc_z,0));

#ifdef THREE_D
  f_interp = interpolate(y_offset, b_interp, t_interp);

  // second, do the interpolation in the back xy slab
  b_interp = interpolate(x_offset, _field_d(llc_x,llc_y,llc_z+1,0),
                                   _field_d(llc_x+1,llc_y,llc_z+1,0));
  t_interp = interpolate(x_offset, _field_d(llc_x,llc_y+1,llc_z+1,0),
                                   _field_d(llc_x+1,llc_y+1,llc_z+1,0));

  bk_interp = interpolate(y_offset, b_interp, t_interp);

  params._Fx[0] = interpolate(z_offset, f_interp, bk_interp);
#else
  params._Fx[0] = interpolate(y_offset, b_interp, t_interp);
#endif



  // first derivative in the y direction
  b_interp = interpolate(x_offset, _field_d(llc_x,llc_y,llc_z,1),
                                   _field_d(llc_x+1,llc_y,llc_z,1));
  t_interp = interpolate(x_offset, _field_d(llc_x,llc_y+1,llc_z,1),
                                   _field_d(llc_x+1,llc_y+1,llc_z,1));

#ifdef THREE_D
  f_interp = interpolate(y_offset, b_interp, t_interp);

  // second, do the interpolation in the back xy slab
  b_interp = interpolate(x_offset, _field_d(llc_x,llc_y,llc_z+1,1),
                                   _field_d(llc_x+1,llc_y,llc_z+1,1));
  t_interp = interpolate(x_offset, _field_d(llc_x,llc_y+1,llc_z+1,1),
                                   _field_d(llc_x+1,llc_y+1,llc_z+1,1));

  bk_interp = interpolate(y_offset, b_interp, t_interp);

  params._Fx[1] = interpolate(z_offset, f_interp, bk_interp);
#else
  params._Fx[1] = interpolate(y_offset, b_interp, t_interp);
#endif


#ifdef THREE_D
  // first derivative in the z direction
  b_interp = interpolate(x_offset, _field_d(llc_x,llc_y,llc_z,2),
                                   _field_d(llc_x+1,llc_y,llc_z,2));
  t_interp = interpolate(x_offset, _field_d(llc_x,llc_y+1,llc_z,2),
                                   _field_d(llc_x+1,llc_y+1,llc_z,2));

  f_interp = interpolate(y_offset, b_interp, t_interp);

  // second, do the interpolation in the back xy slab
  b_interp = interpolate(x_offset, _field_d(llc_x,llc_y,llc_z+1,2),
                                   _field_d(llc_x+1,llc_y,llc_z+1,2));
  t_interp = interpolate(x_offset, _field_d(llc_x,llc_y+1,llc_z+1,2),
                                   _field_d(llc_x+1,llc_y+1,llc_z+1,2));

  bk_interp = interpolate(y_offset, b_interp, t_interp);

  params._Fx[2] = interpolate(z_offset, f_interp, bk_interp);
#endif

  //params._F = -params._F;
  //params._Fx = -params._Fx;

  return true;
}

//------------------------------------------------------------------------
// Function    : readDataFile()
// Description : 
//------------------------------------------------------------------------
void DistanceFieldwSF::readDistanceFieldFile( const char *df_file_base,
                                              int num_header_lines )
{
  cout << "Reading in the distance field files ..." << endl;

  _min_isovalue = 1.0e6;
  _max_isovalue = -1.0e6;

  // create the file names
  char df_file[150];     sprintf( df_file, "%s.vol", df_file_base );
  char sf_file[150];     sprintf( sf_file, "%s_sf.vol", df_file_base );
  char df_dx_file[150];  sprintf( df_dx_file, "%s_dx.vol", df_file_base );
  char df_dy_file[150];  sprintf( df_dy_file, "%s_dy.vol", df_file_base );
#ifdef THREE_D
  char df_dz_file[150];  sprintf( df_dz_file, "%s_dz.vol", df_file_base );
#endif

  // first open one file to get the dimensions
  FILE *in_file_field = fopen( df_file, "r" );
  if ( !in_file_field )
  {
    cout << "DistanceField::Error reading field file " << df_file << endl;
    exit( 1 );
  }

  char tmp[100];
  fgets( tmp, 100, in_file_field );
  fscanf( in_file_field, " %i %i %i", &_xdim, &_ydim, &_zdim );

  fclose( in_file_field );

  _field.resize( _xdim, _ydim, _zdim );
  _sf.resize( _xdim, _ydim, _zdim );
#ifdef THREE_D
  _field_d.resize( _xdim, _ydim, _zdim, 3 );
#else
  _field_d.resize( _xdim, _ydim, _zdim, 2 );
#endif


  // open the files
  in_file_field = fopen( df_file, "rb" );
  if ( !in_file_field )
  {
    cout << "DistanceFieldwSF::Error reading field file." << endl;
    exit( 1 );
  }

  FILE *in_file_sf = fopen( sf_file, "rb" );
  if ( !in_file_sf )
  {
    cout << "DistanceFieldwSF::Error reading sizing field file." << endl;
    exit( 1 );
  }

  // first derivatives
  FILE *in_file_field_dx = fopen( df_dx_file, "rb" );
  if ( !in_file_field_dx )
  {
    cout << "DistanceFieldwSF::Error reading field dx file." << endl;
    exit( 1 );
  }

  FILE *in_file_field_dy = fopen( df_dy_file, "rb" );
  if ( !in_file_field_dy )
  {
    cout << "DistanceFieldwSF::Error reading field dy file." << endl;
    exit( 1 );
  }

#ifdef THREE_D
  FILE *in_file_field_dz = fopen( df_dz_file, "rb" );
  if ( !in_file_field_dz )
  {
    cout << "DistanceFieldwSF::Error reading field dz file." << endl;
    exit( 1 );
  }
#endif
  
  // read in and discard the headers
  for ( int i = 0; i < num_header_lines; i++ )
  {
    fgets( tmp, 100, in_file_field );
    fgets( tmp, 100, in_file_sf );
    fgets( tmp, 100, in_file_field_dx );
    fgets( tmp, 100, in_file_field_dy );
#ifdef THREE_D
    fgets( tmp, 100, in_file_field_dz );
#endif
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

        // read in the sizing field
        fread( &_sf(i,j,k), sizeof(float), 1, in_file_sf );

        if ( _sf(i,j,k) == 0.0 )
          _sf(i,j,k) = 0.5;

        if ( _sf(i,j,k) != -1.0 )
        {
          if ( _sf(i,j,k) < _min_sf ) _min_sf = _sf(i,j,k);
          if ( _sf(i,j,k) > _max_sf ) _max_sf = _sf(i,j,k);
        }

        // the first derivatives
        fread( &_field_d(i,j,k,0), sizeof(float), 1, in_file_field_dx );
        fread( &_field_d(i,j,k,1), sizeof(float), 1, in_file_field_dy );
#ifdef THREE_D
        fread( &_field_d(i,j,k,2), sizeof(float), 1, in_file_field_dz );
#endif
      }

  cout << "min and max sf = " << _min_sf << " " << _max_sf << endl;

  // close the files
  fclose( in_file_field );
  fclose( in_file_field_dx );
  fclose( in_file_field_dy );
#ifdef THREE_D
  fclose( in_file_field_dz );
#endif

  cout << "              DONE" << endl;
}

