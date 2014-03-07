#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>

typedef unsigned char byte;

#include <system/domain/DistanceField.h>

using namespace std;
using namespace particle_sys;

//------------------------------------------------------------------------
// Function    : Constructor and Destructor
// Description : 
//------------------------------------------------------------------------
DistanceField::DistanceField( int xdim, int ydim, int zdim,
                              const char *df_file_base ) : Surface()
{
#ifdef THREE_D
  _xdim = xdim; _ydim = ydim; _zdim = zdim;
  _field.resize( _xdim, _ydim, _zdim );
  _field_d.resize( _xdim, _ydim, _zdim, 3 );
  _field_dd.resize( _xdim, _ydim, _zdim, 6 );

  readDistanceFieldFile( df_file_base );

  _isovalue = 0.0;
#endif
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
bool DistanceField::computeSurfacePointParams( const vector_type &pos, 
                                               SurfacePointParams 
                                               &params,
                                               bool computeHessian ) const
{
#ifdef TWO_D
#endif

#ifdef THREE_D

  // get the lower left corner of the bounding voxel
  //   NOTE: if the position is exactly on the high end boundary, make 
  //         sure to not run outside of the domain!
  int llc_x = min( (int)floor(pos(0)), _xdim-2 );
  int llc_y = min( (int)floor(pos(1)), _ydim-2 );
  int llc_z = min( (int)floor(pos(2)), _zdim-2 );

  // first, do the interpolation in the front xy slab
  float x_offset = pos(0) - (float)llc_x;
  float y_offset = pos(1) - (float)llc_y;
  float z_offset = pos(2) - (float)llc_z; 

  // the field
  float b_interp = interpolate(x_offset, _field(llc_x,llc_y,llc_z),
                                         _field(llc_x+1,llc_y,llc_z));
                                 
  float t_interp = interpolate(x_offset, _field(llc_x,llc_y+1,llc_z),
                                         _field(llc_x+1,llc_y+1,llc_z));

  float f_interp = interpolate(y_offset, b_interp, t_interp);

  // second, do the interpolation in the back xy slab
  b_interp = interpolate(x_offset, _field(llc_x,llc_y,llc_z+1),
                                   _field(llc_x+1,llc_y,llc_z+1));
  t_interp = interpolate(x_offset, _field(llc_x,llc_y+1,llc_z+1), 
                                   _field(llc_x+1,llc_y+1,llc_z+1));

  float bk_interp = interpolate(y_offset, b_interp, t_interp);

  params._F = -interpolate(z_offset, f_interp, bk_interp) - _isovalue;

  // first derivative in the x direction
  b_interp = interpolate(x_offset, _field_d(llc_x,llc_y,llc_z,0),
                                   _field_d(llc_x+1,llc_y,llc_z,0));
  t_interp = interpolate(x_offset, _field_d(llc_x,llc_y+1,llc_z,0),
                                   _field_d(llc_x+1,llc_y+1,llc_z,0));

  f_interp = interpolate(y_offset, b_interp, t_interp);

  // second, do the interpolation in the back xy slab
  b_interp = interpolate(x_offset, _field_d(llc_x,llc_y,llc_z+1,0),
                                   _field_d(llc_x+1,llc_y,llc_z+1,0));
  t_interp = interpolate(x_offset, _field_d(llc_x,llc_y+1,llc_z+1,0),
                                   _field_d(llc_x+1,llc_y+1,llc_z+1,0));

  bk_interp = interpolate(y_offset, b_interp, t_interp);

  params._Fx[0] = -interpolate(z_offset, f_interp, bk_interp);



  // first derivative in the y direction
  b_interp = interpolate(x_offset, _field_d(llc_x,llc_y,llc_z,1),
                                   _field_d(llc_x+1,llc_y,llc_z,1));
  t_interp = interpolate(x_offset, _field_d(llc_x,llc_y+1,llc_z,1),
                                   _field_d(llc_x+1,llc_y+1,llc_z,1));

  f_interp = interpolate(y_offset, b_interp, t_interp);

  // second, do the interpolation in the back xy slab
  b_interp = interpolate(x_offset, _field_d(llc_x,llc_y,llc_z+1,1),
                                   _field_d(llc_x+1,llc_y,llc_z+1,1));
  t_interp = interpolate(x_offset, _field_d(llc_x,llc_y+1,llc_z+1,1),
                                   _field_d(llc_x+1,llc_y+1,llc_z+1,1));

  bk_interp = interpolate(y_offset, b_interp, t_interp);

  params._Fx[1] = -interpolate(z_offset, f_interp, bk_interp);



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

  params._Fx[2] = -interpolate(z_offset, f_interp, bk_interp);



  // second derivative in the xx direction
  b_interp = interpolate(x_offset, _field_dd(llc_x,llc_y,llc_z,0),
                                   _field_dd(llc_x+1,llc_y,llc_z,0));
  t_interp = interpolate(x_offset, _field_dd(llc_x,llc_y+1,llc_z,0),
                                   _field_dd(llc_x+1,llc_y+1,llc_z,0));

  f_interp = interpolate(y_offset, b_interp, t_interp);

  b_interp = interpolate(x_offset, _field_dd(llc_x,llc_y,llc_z+1,0),
                                   _field_dd(llc_x+1,llc_y,llc_z+1,0));
  t_interp = interpolate(x_offset, _field_dd(llc_x,llc_y+1,llc_z+1,0),
                                   _field_dd(llc_x+1,llc_y+1,llc_z+1,0));

  bk_interp = interpolate(y_offset, b_interp, t_interp);

  params._Fxx(0,0) = interpolate(z_offset, f_interp, bk_interp);


  // second derivative in the xy direction
  b_interp = interpolate(x_offset, _field_dd(llc_x,llc_y,llc_z,1),
                                   _field_dd(llc_x+1,llc_y,llc_z,1));
  t_interp = interpolate(x_offset, _field_dd(llc_x,llc_y+1,llc_z,1),
                                   _field_dd(llc_x+1,llc_y+1,llc_z,1));

  f_interp = interpolate(y_offset, b_interp, t_interp);

  b_interp = interpolate(x_offset, _field_dd(llc_x,llc_y,llc_z+1,1),
                                   _field_dd(llc_x+1,llc_y,llc_z+1,1));
  t_interp = interpolate(x_offset, _field_dd(llc_x,llc_y+1,llc_z+1,1),
                                   _field_dd(llc_x+1,llc_y+1,llc_z+1,1));

  bk_interp = interpolate(y_offset, b_interp, t_interp);

  params._Fxx(0,1) = params._Fxx(1,0) = 
    interpolate(z_offset, f_interp, bk_interp);


  // second derivative in the xz direction
  b_interp = interpolate(x_offset, _field_dd(llc_x,llc_y,llc_z,2),
                                   _field_dd(llc_x+1,llc_y,llc_z,2));
  t_interp = interpolate(x_offset, _field_dd(llc_x,llc_y+1,llc_z,2),
                                   _field_dd(llc_x+1,llc_y+1,llc_z,2));

  f_interp = interpolate(y_offset, b_interp, t_interp);

  b_interp = interpolate(x_offset, _field_dd(llc_x,llc_y,llc_z+1,2),
                                   _field_dd(llc_x+1,llc_y,llc_z+1,2));
  t_interp = interpolate(x_offset, _field_dd(llc_x,llc_y+1,llc_z+1,2),
                                   _field_dd(llc_x+1,llc_y+1,llc_z+1,2));

  bk_interp = interpolate(y_offset, b_interp, t_interp);

  params._Fxx(0,2) = params._Fxx(2,0) = 
    interpolate(z_offset, f_interp, bk_interp);


  // second derivative in the yy direction
  b_interp = interpolate(x_offset, _field_dd(llc_x,llc_y,llc_z,3),
                                   _field_dd(llc_x+1,llc_y,llc_z,3));
  t_interp = interpolate(x_offset, _field_dd(llc_x,llc_y+1,llc_z,3),
                                   _field_dd(llc_x+1,llc_y+1,llc_z,3));

  f_interp = interpolate(y_offset, b_interp, t_interp);

  b_interp = interpolate(x_offset, _field_dd(llc_x,llc_y,llc_z+1,3),
                                   _field_dd(llc_x+1,llc_y,llc_z+1,3));
  t_interp = interpolate(x_offset, _field_dd(llc_x,llc_y+1,llc_z+1,3),
                                   _field_dd(llc_x+1,llc_y+1,llc_z+1,3));

  bk_interp = interpolate(y_offset, b_interp, t_interp);

  params._Fxx(1,1) = interpolate(z_offset, f_interp, bk_interp);


  // second derivative in the yz direction
  b_interp = interpolate(x_offset, _field_dd(llc_x,llc_y,llc_z,4),
                                   _field_dd(llc_x+1,llc_y,llc_z,4));
  t_interp = interpolate(x_offset, _field_dd(llc_x,llc_y+1,llc_z,4),
                                   _field_dd(llc_x+1,llc_y+1,llc_z,4));

  f_interp = interpolate(y_offset, b_interp, t_interp);

  b_interp = interpolate(x_offset, _field_dd(llc_x,llc_y,llc_z+1,4),
                                   _field_dd(llc_x+1,llc_y,llc_z+1,4));
  t_interp = interpolate(x_offset, _field_dd(llc_x,llc_y+1,llc_z+1,4),
                                   _field_dd(llc_x+1,llc_y+1,llc_z+1,4));

  bk_interp = interpolate(y_offset, b_interp, t_interp);

  params._Fxx(1,2) = params._Fxx(2,1) = 
    interpolate(z_offset, f_interp, bk_interp);


  // second derivative in the zz direction
  b_interp = interpolate(x_offset, _field_dd(llc_x,llc_y,llc_z,5),
                                   _field_dd(llc_x+1,llc_y,llc_z,5));
  t_interp = interpolate(x_offset, _field_dd(llc_x,llc_y+1,llc_z,5),
                                   _field_dd(llc_x+1,llc_y+1,llc_z,5));

  f_interp = interpolate(y_offset, b_interp, t_interp);

  b_interp = interpolate(x_offset, _field_dd(llc_x,llc_y,llc_z+1,5),
                                   _field_dd(llc_x+1,llc_y,llc_z+1,5));
  t_interp = interpolate(x_offset, _field_dd(llc_x,llc_y+1,llc_z+1,5),
                                   _field_dd(llc_x+1,llc_y+1,llc_z+1,5));

  bk_interp = interpolate(y_offset, b_interp, t_interp);

  params._Fxx(2,2) = interpolate(z_offset, f_interp, bk_interp);

  //params._F = -params._F;
  //params._Fx = -params._Fx;

#endif

  return true;
}

//------------------------------------------------------------------------
// Function    : readDataFile()
// Description : 
//------------------------------------------------------------------------
void DistanceField::readDistanceFieldFile( const char *df_file_base,
                                           int num_header_lines )
{
  cout << "Reading in the distance field files ..." << endl;
#ifdef TWO_D
#endif

#ifdef THREE_D

  _min_isovalue = 1.0e6;
  _max_isovalue = -1.0e6;

  // create the file names
  char df_file[150];     sprintf( df_file, "%s.vol", df_file_base );
  char df_dx_file[150];  sprintf( df_dx_file, "%s_dx.vol", df_file_base );
  char df_dy_file[150];  sprintf( df_dy_file, "%s_dy.vol", df_file_base );
  char df_dz_file[150];  sprintf( df_dz_file, "%s_dz.vol", df_file_base );
  char df_dxx_file[150]; sprintf( df_dxx_file, "%s_dxx.vol", df_file_base );
  char df_dyy_file[150]; sprintf( df_dyy_file, "%s_dyy.vol", df_file_base );
  char df_dzz_file[150]; sprintf( df_dzz_file, "%s_dzz.vol", df_file_base );
  char df_dxy_file[150]; sprintf( df_dxy_file, "%s_dxy.vol", df_file_base );
  char df_dxz_file[150]; sprintf( df_dxz_file, "%s_dxz.vol", df_file_base );
  char df_dyz_file[150]; sprintf( df_dyz_file, "%s_dyz.vol", df_file_base );

  // open the files
  FILE *in_file_field = fopen( df_file, "rb" );
  if ( !in_file_field )
  {
    cout << "DistanceField::Error reading field file." << endl;
    exit( 1 );
  }

  // first derivatives
  FILE *in_file_field_dx = fopen( df_dx_file, "rb" );
  if ( !in_file_field_dx )
  {
    cout << "DistanceField::Error reading field dx file." << endl;
    exit( 1 );
  }

  FILE *in_file_field_dy = fopen( df_dy_file, "rb" );
  if ( !in_file_field_dy )
  {
    cout << "DistanceField::Error reading field dy file." << endl;
    exit( 1 );
  }

  FILE *in_file_field_dz = fopen( df_dz_file, "rb" );
  if ( !in_file_field_dz )
  {
    cout << "DistanceField::Error reading field dz file." << endl;
    exit( 1 );
  }

  // second derivatives
  FILE *in_file_field_dxx = fopen( df_dxx_file, "rb" );
  if ( !in_file_field_dxx )
  {
    cout << "DistanceField::Error reading field dxx file." << endl;
    exit( 1 );
  }

  FILE *in_file_field_dyy = fopen( df_dyy_file, "rb" );
  if ( !in_file_field_dyy )
  {
    cout << "DistanceField::Error reading field dyy file." << endl;
    exit( 1 );
  }

  FILE *in_file_field_dzz = fopen( df_dzz_file, "rb" );
  if ( !in_file_field_dzz )
  {
    cout << "DistanceField::Error reading field dzz file." << endl;
    exit( 1 );
  }

  FILE *in_file_field_dxy = fopen( df_dxy_file, "rb" );
  if ( !in_file_field_dxy )
  {
    cout << "DistanceField::Error reading field dxy file." << endl;
    exit( 1 );
  }

  FILE *in_file_field_dxz = fopen( df_dxz_file, "rb" );
  if ( !in_file_field_dxz )
  {
    cout << "DistanceField::Error reading field dxz file." << endl;
    exit( 1 );
  }

  FILE *in_file_field_dyz = fopen( df_dyz_file, "rb" );
  if ( !in_file_field_dyz )
  {
    cout << "DistanceField::Error reading field dyz file." << endl;
    exit( 1 );
  }
  
  // read in and discard the headers
  char tmp[100];
  for ( int i = 0; i < num_header_lines; i++ )
  {
    fgets( tmp, 100, in_file_field );
    fgets( tmp, 100, in_file_field_dx );
    fgets( tmp, 100, in_file_field_dy );
    fgets( tmp, 100, in_file_field_dz );
    fgets( tmp, 100, in_file_field_dxx );
    fgets( tmp, 100, in_file_field_dyy );
    fgets( tmp, 100, in_file_field_dzz );
    fgets( tmp, 100, in_file_field_dxy );
    fgets( tmp, 100, in_file_field_dxz );
    fgets( tmp, 100, in_file_field_dyz );
  }

  //
  // for the second derivatives, here is how the matirx is stored
  // 
  // field_dd -> indices   | 0 1 2 |
  //                       |   3 4 |
  //                       |     5 |
  //

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

        // the first derivatives
        fread( &_field_d(i,j,k,0), sizeof(float), 1, in_file_field_dx );
        fread( &_field_d(i,j,k,1), sizeof(float), 1, in_file_field_dy );
        fread( &_field_d(i,j,k,2), sizeof(float), 1, in_file_field_dz );

        // the second derivatives
        fread( &_field_dd(i,j,k,0), sizeof(float), 1, in_file_field_dxx );
        fread( &_field_dd(i,j,k,1), sizeof(float), 1, in_file_field_dxy );
        fread( &_field_dd(i,j,k,2), sizeof(float), 1, in_file_field_dxz );
        fread( &_field_dd(i,j,k,3), sizeof(float), 1, in_file_field_dyy );
        fread( &_field_dd(i,j,k,4), sizeof(float), 1, in_file_field_dyz );
        fread( &_field_dd(i,j,k,5), sizeof(float), 1, in_file_field_dzz );
      }

  // close the files
  fclose( in_file_field );
  fclose( in_file_field_dx );
  fclose( in_file_field_dy );
  fclose( in_file_field_dz );
  fclose( in_file_field_dxx );
  fclose( in_file_field_dyy );
  fclose( in_file_field_dzz );
  fclose( in_file_field_dxy );
  fclose( in_file_field_dxz );
  fclose( in_file_field_dyz );
#endif

  cout << "              DONE" << endl;
}

