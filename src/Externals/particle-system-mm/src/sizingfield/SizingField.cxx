#include <cstdlib>
#include <iostream>
#include <LFS.h>
#include <SurfaceParameters.h>
#include <SizingField.h>
#include <NrrdIO.h>

using namespace std;

//------------------------------------------------------------------------
// Function    : Constructor and Destructor
// Description : 
//------------------------------------------------------------------------
SizingField::SizingField( const char *file_basename,
                          const char *volume_names, int main_indicator, bool useNrrdIO )
  : _useNrrdIO(useNrrdIO)
{  
  _base_filename = new char[250];
  sprintf( _base_filename, "%s", file_basename );

  _volume_files = new char[250];
  sprintf( _volume_files, "%s", volume_names );
  
  
  _redo_file = NULL;

  _main_indicator = main_indicator;

  _delta = 0.4f;
  _epsilon = 0.5f;

  _min_sf = 0.0f;
  _restart_files=true;
}

SizingField::SizingField( const char *file_basename,
                          const char *volume_names,
                          const char *redo_name,
                          bool useNrrdIO )
  : _useNrrdIO(useNrrdIO)
{  
  _base_filename = new char[250];
  sprintf( _base_filename, "%s", file_basename );

  _volume_files = new char[250];
  sprintf( _volume_files, "%s", volume_names );

  _redo_file = new char[250];
  sprintf( _redo_file, "%s", redo_name );

  _delta = 0.4f;
  _epsilon = 0.5f;

  _min_sf = 0.0f;
  _restart_files=true;
}

SizingField::~SizingField()
{
  delete [] _base_filename;
  delete [] _volume_files;
  if ( _redo_file ) delete [] _redo_file;
}

//------------------------------------------------------------------------
// Function    : generateInitialSizingField()
// Description : 
//------------------------------------------------------------------------
void SizingField::generateInitialSizingField()
{
  // create the IO scalar field, even if there is only one material
  _io_field = new IOScalarField( _volume_files,
                                 IOScalarField::INTERPOLATING,
                                 _main_indicator );


  // get the dimensions of the field
  vec<3> start, end;
  _io_field->domain( start, end );
  
  _xdim = _io_field->xdim();
  _ydim = _io_field->ydim();
  _zdim = _io_field->zdim();

  _scaling[0] = end(0) / (float)(_xdim-1);
  _scaling[1] = end(1) / (float)(_ydim-1);
  _scaling[2] = end(2) / (float)(_zdim-1);
  
  _centers[0] = _io_field->xcenter();
  _centers[1] = _io_field->ycenter();
  _centers[2] = _io_field->zcenter();

  cout << start << "  " << end << endl;
  cout << _xdim << " " << _ydim << " " << _zdim << endl;
  cout << _scaling[0] << " " << _scaling[1] << " " << _scaling[2] << endl;

  // create the binary crossing/not-crossing volume
  createCrossingVolume();

  // create the initial sizing field and initialize the maxk values
  _h0.resize( _xdim, _ydim, _zdim );
  
  //D Swenson  06/29/09  We aren't using this currently because the medial axis distance is sufficient.
  //initializeh0ToMaxk();

  // compare the maxk to the lfs
  compareMaxkToLFS();

  // scale the field by epsilon, and find the max value
  _max_size = 0;
  int cnt=0;
  cout << "Scaling the initial sizing field by epsilon..." << endl;
  for ( int k = 0; k < _zdim; k++ )
    for ( int j = 0; j < _ydim; j++ )
      for ( int i = 0; i < _xdim; i++ )
      {
        if ( _crossing(i,j,k) )
        {
         //cout<<cnt<<endl;
           _h0(i,j,k) *= _epsilon;

          // compare again the minimum sizing field value
          _h0(i,j,k)  = max( _min_sf, _h0(i,j,k) );
          //_h0(i,j,k)  = min( _max_size, _h0(i,j,k) );
          _max_size = max( _max_size, _h0(i,j,k) );
        cnt++;
        }
      }
  
  cout << "   DONE." << endl << endl;

  // writing out the h0 file
  cout << "Writing out the initial sizing field file..." << endl;

  if (_restart_files)
  {  
    if (_useNrrdIO)
    {
      char h0_filename[250];
      sprintf( h0_filename, "%s_sf_init.nrrd", _base_filename );
      cout << "Writing .nrrd..." << endl;
      outputNrrdFile(h0_filename,_h0);
    }
    else
    {
      char h0_filename[250];
      sprintf( h0_filename, "%s_sf_init.vol", _base_filename );
      cout << "Writing out .vol..." << endl;
      outputVolFile(h0_filename,_h0);
    }
  }

  // delete IO field
  delete _io_field;
  cout << "   DONE." << endl << endl;
}

//------------------------------------------------------------------------
// Function    : redoInitialSizingField()
// Description : 
//------------------------------------------------------------------------
void SizingField::redoInitialSizingField()
{

  // create the IO scalar field, even if there is only one material
  _io_field = new IOScalarField( _volume_files,
                                 IOScalarField::INTERPOLATING,
                                 _main_indicator );

  // get the dimensions of the field
  vec<3> start, end;
  _io_field->domain( start, end );
  _xdim = _io_field->xdim();
  _ydim = _io_field->ydim();
  _zdim = _io_field->zdim();

  _scaling[0] = end(0) / (float)(_xdim-1);
  _scaling[1] = end(1) / (float)(_ydim-1);
  _scaling[2] = end(2) / (float)(_zdim-1);

  // create the binary crossing/not-crossing volume
  createCrossingVolume();

  // create the initial sizing field and initialize the maxk values
  _h0.resize( _xdim, _ydim, _zdim );
  initializeh0ToMaxk();

  // delete IO field
  delete _io_field;

  // compare the maxk to the lfs
  compareMaxkToLFS();

  // make a redo volume
  ScalarField *redo = new ScalarField( _redo_file, ScalarField::LINEAR,
                                       1.0, 1.0, 1.0 );

  // scale the field by epsilon, and find the max value
  _max_size = 0.0;
  vec<3> pos;
  SurfacePointParams params;
  cout << "Scaling the initial sizing field by epsilon..." << endl;
  bool doredo;
  for ( int k = 0; k < _zdim; k++ )
    for ( int j = 0; j < _ydim; j++ )
      for ( int i = 0; i < _xdim; i++ )
      {
        doredo = false;
        // find out if this point is a redo point
        pos.set( _scaling[0]*(float)i,
                 _scaling[1]*(float)j,
                 _scaling[2]*(float)k );
        if ( redo->inBounds( pos ) )
        {
          redo->computeSurfacePointParams( pos, params, false,
                                           false );
          if ( params._F )
            doredo = true;         
        }
        if ( doredo )
          _h0(i,j,k)  = _min_sf/2.0f;
        else if ( _crossing(i,j,k) )
        {
          _h0(i,j,k) *= _epsilon;

          _h0(i,j,k) = max( _min_sf, _h0(i,j,k) );
          _max_size = max( _max_size, _h0(i,j,k) );
        }
      }
  delete redo;
  
  cout << "   DONE." << endl << endl;

  // writing out the h0 file
  cout << "Writing out the initial sizing field file..." << endl;

    
  if (_useNrrdIO)
  {
    char h0_filename[250];
    sprintf( h0_filename, "%s_sf_init.nrrd", _base_filename );
    cout << "Writing .nrrd..." << endl;
    outputNrrdFile(h0_filename,_h0);
  }
  else
  {
    char h0_filename[250];
    sprintf( h0_filename, "%s_sf_init.vol", _base_filename );
    cout << "Writing .vol..." << endl;
    outputVolFile(h0_filename,_h0);
  }

  cout << "   DONE." << endl << endl;
}


//------------------------------------------------------------------------
// Function    : generateSizingField()
// Description : 
//------------------------------------------------------------------------
void SizingField::generateSizingField()
{
  // smooth h0
  cout << "Gradient limiting the sizing field..." << endl;

  // create an update volume
  array3D<float> update( _xdim, _ydim, _zdim );

  // copy the initial sizing field into the smooth sizing field, making
  //   the outside part be only slightly larger valued than the crossings
  _h.resize( _xdim, _ydim, _zdim );
  for ( int k = 0; k < _zdim; k++ )
    for ( int j = 0; j < _ydim; j++ )
      for ( int i = 0; i < _xdim; i++ )
      {
        if ( _crossing(i,j,k) )
          _h(i,j,k) = _h0(i,j,k);
        else
          _h(i,j,k) = _max_size+_delta;
      }

  float 
    max_change, 
    max_relative_change,
    this_value,
    this_change, 
    grad_before_x, 
    grad_before_y, 
    grad_before_z, 
    grad_after_x, 
    grad_after_y, 
    grad_after_z, 
    grad_mag, 
    grad_before, 
    grad_after, 
    grad_x, 
    grad_y, 
    grad_z, 
    dt;
 
  // iteratively smooth the field
  bool done = false;
  size_t itercount = 0;
  while (!done)
  {
    max_change = max_relative_change = 0.0;
    for ( int k = 0; k < _zdim; k++ ) 
      for ( int j = 0; j < _ydim; j++ ) 
        for ( int i = 0; i < _xdim; i++ ) 
        {
          this_value = _h(i,j,k);

          if ( i < (_xdim - 1) ) 
            grad_before_x = (_h(i+1,j,k) - this_value) / _scaling[0];
          else 
            //            grad_before_x = this_value - _h(i-1,j,k);
            grad_before_x = 0.0f;
    			
          if ( i > 0 )
            grad_after_x = (this_value - _h(i-1,j,k)) / _scaling[0];
          else 
            //            grad_after_x = _h(i+1,j,k) - this_value;
            grad_after_x = 0.0f;

          // things are allways *decreasing* to preserve feature size.
          grad_before = min(grad_before_x, 0.0f);
          grad_after = max(grad_after_x, 0.0f);		       
          grad_x = grad_after*grad_after + grad_before*grad_before;

          if ( j < (_ydim - 1) ) 
            grad_before_y = (_h(i,j+1,k) - this_value) / _scaling[1];
          else 
            //            grad_before_y = this_value - _h(i,j-1,k);
            grad_before_y = 0.0f;

          if ( j > 0 )
            grad_after_y = (this_value - _h(i,j-1,k)) / _scaling[1];
          else
            //            grad_after_y = _h(i,j+1,k) - this_value;
            grad_after_y = 0.0f;

          grad_before = min(grad_before_y, 0.0f);
          grad_after = max(grad_after_y, 0.0f);		       
          grad_y = grad_after*grad_after + grad_before*grad_before;

          if ( k < (_zdim - 1) ) 
            grad_before_z = (_h(i,j,k+1) - this_value) / _scaling[2];
          else 
            //            grad_before_z = this_value - _h(i,j,k-1);
            grad_before_z = 0.0f;

          if ( k > 0 )
            grad_after_z = (this_value - _h(i,j,k-1)) / _scaling[2];
          else
            //            grad_after_z = _h(i,j,k+1) - this_value;
            grad_after_z = 0.0f;

          grad_before = min(grad_before_z, 0.0f);
          grad_after = max(grad_after_z, 0.0f);		       
          grad_z = grad_after*grad_after + grad_before*grad_before;

          // take the magnitude of the gradient
          grad_mag = sqrt(grad_x + grad_y + grad_z);

          // only want an update value if the gradient needs to be limited
          // ---> update is always <= 0

          update(i,j,k) = min(grad_mag, _delta) - grad_mag;

//           cerr << "grad_mag: " << grad_mag << endl;
//           cerr << "_delta: " << _delta << endl;
//           cerr << "min(grad_mag, _delta): " << min(grad_mag, _delta) << endl;
//           cerr << "update: " << update(i,j,k) << endl;
        
          // get the relative change
          this_change = update(i,j,k);
          max_change = max(fabs(this_change), max_change);

          this_change = fabs(this_change)/(_h(i,j,k) + static_cast<float>(EPSILON) );
          max_relative_change = max(this_change, max_relative_change);
        }
  
    cout << "max change is " << max_change << "   ";

    float min_scaling_tmp = min(_scaling[0], _scaling[1]);
    float min_scaling = min(min_scaling_tmp, _scaling[2]);

    dt = (min_scaling*1.0f/6.0f);
    for ( int k = 0; k < _zdim; k++ ) 
      for ( int j = 0; j < _ydim; j++ ) 
        for ( int i = 0; i < _xdim; i++ ) { 
           _h(i,j,k) += dt*update(i,j,k);
        }
    // stop when the maximum relative change gets small
    // probably should stop if too many iterations also
    if (max_change < (1.0e-4/min_scaling))
      done = true;
    cout << "Done with iteration number: " << itercount << endl;
    itercount++;
  }
  cout << "   DONE." << endl << endl;
  
  if (_useNrrdIO)
  {
      char h_filename[250];
      sprintf( h_filename, "%s_sf.nrrd", _base_filename );
      outputNrrdFile(h_filename,_h);
  }
  else
  {
      char h_filename[250];
      sprintf( h_filename, "%s_sf.vol", _base_filename );
      outputVolFile(h_filename,_h);
  }


  // writing out the h file
  cout << "Writing out the sizing field file..." << endl;

  cout << "   DONE." << endl << endl;
}

//------------------------------------------------------------------------
// Function    : generateConstantSizingField()
// Description : Generate a constant-valued sizing field.
//------------------------------------------------------------------------
void SizingField::generateConstantSizingField( float sf_value )
{
	// create the IO scalar field, even if there is only one material
	_io_field = new IOScalarField( _volume_files,
								  IOScalarField::INTERPOLATING,
								  _main_indicator );
	
	// get the dimensions of the field
	vec<3> start, end;
	_io_field->domain( start, end );
	
	_xdim = _io_field->xdim();
	_ydim = _io_field->ydim();
	_zdim = _io_field->zdim();
	
	_scaling[0] = end(0) / (float)(_xdim-1);
	_scaling[1] = end(1) / (float)(_ydim-1);
	_scaling[2] = end(2) / (float)(_zdim-1);
	
	_centers[0] = _io_field->xcenter();
	_centers[1] = _io_field->ycenter();
	_centers[2] = _io_field->zcenter();
	
	cout << start << "  " << end << endl;
	cout << _xdim << " " << _ydim << " " << _zdim << endl;
	cout << _scaling[0] << " " << _scaling[1] << " " << _scaling[2] << endl;
	
	// create the binary crossing/not-crossing volume
	createCrossingVolume();
	
	// smooth h0
	cout << "Setting up the constant sizing field..." << endl;
	
	// create an update volume
	array3D<float> update( _xdim, _ydim, _zdim );
	
	// copy the initial sizing field into the smooth sizing field, making
	//   the outside part be only slightly larger valued than the crossings
	_h0.resize(_xdim, _ydim, _zdim );
	_h.resize( _xdim, _ydim, _zdim );
	
	_max_size = sf_value;
	
	for ( int k = 0; k < _zdim; k++ )
		for ( int j = 0; j < _ydim; j++ )
			for ( int i = 0; i < _xdim; i++ )
			{
				_h0(i,j,k) = _h(i,j,k) = sf_value; //assign constant sizing field
			}
	
	// writing out the h file
	cout << "Writing out the sizing field file..." << endl;
	
	if (_useNrrdIO)
	{
		char h_filename[250];
		sprintf( h_filename, "%s_sf_init.nrrd", _base_filename );
		outputNrrdFile(h_filename,_h0);
		sprintf( h_filename, "%s_sf.nrrd", _base_filename );
		outputNrrdFile(h_filename,_h);
	}
	else
	{
		char h_filename[250];
		sprintf( h_filename, "%s_sf_init.vol", _base_filename );
		outputNrrdFile(h_filename,_h0);
		sprintf( h_filename, "%s_sf.vol", _base_filename );
		outputVolFile(h_filename,_h);
	}
	
	// delete IO field
	delete _io_field;
	cout << "   DONE." << endl << endl;
}

//------------------------------------------------------------------------
// Function    : createCrossingVolume()
// Description : 
//------------------------------------------------------------------------
void SizingField::createCrossingVolume()
{
  _crossing.resize( _xdim, _ydim, _zdim );
  
  // first, check and see if a crossing file already exists

  cout << "Generating the crossing volume..." << endl;

  vec<3> pos;
  float center_val, neigh_val;
  bool inside, neigh_inside;
  SurfacePointParams params;

  // go through the volume and store true to indicate if the grid node
  //   bounds the isosurface
  for ( int k = 0; k < _zdim; k++ )
  {
    //cout << k << endl;
    for ( int j = 0; j < _ydim; j++ )
      for ( int i = 0; i < _xdim; i++ )
      {
         // get the surface value at this grid node
        pos.set( _scaling[0]*(float)i,
                 _scaling[1]*(float)j,
                 _scaling[2]*(float)k );
        _io_field->computeSurfacePointParams( pos, params );
        center_val = params._F;

        // check if this point is actually on the surface
        if ( center_val == 0.0 )
        {
          _crossing(i,j,k) = true;
          continue;
        }

        // find out if we are inside or outside
        inside = (center_val < 0.0);

        _crossing(i,j,k) = false;

        // look at the x neighbors
        for ( int x = -1; x < 2; x += 2 )
        {
          // check that this x is inside
          if ( ((i+x) < 0) || ((i+x) >= _xdim) )
            continue;

          // find this neighbor's value
          pos.set( _scaling[0]*(float)(i+x),
                   _scaling[1]*(float)j,
                   _scaling[2]*(float)k );
          _io_field->computeSurfacePointParams( pos, params );
          neigh_val = params._F;

          // determine if the neighbor is inside or outside
          neigh_inside = (neigh_val < 0.0);

          // check if we have a crossing
          if ( (inside && !neigh_inside) || (!inside && neigh_inside) )
            _crossing(i,j,k) = true;
        }

        // look at the y neighbors
        for ( int y = -1; y < 2; y += 2 )
        {
          // check that this y is inside
          if ( ((j+y) < 0) || ((j+y) >= _ydim) )
            continue;

          // find this neighbor's value
          pos.set( _scaling[0]*(float)i,
                   _scaling[1]*(float)(j+y),
                   _scaling[2]*(float)k );
          _io_field->computeSurfacePointParams( pos, params );
          neigh_val = params._F;
 
          // determine if the neighbor is inside or outside
          neigh_inside = (neigh_val < 0.0);

          // check if we have a crossing
          if ( (inside && !neigh_inside) || (!inside && neigh_inside) )
            _crossing(i,j,k) = true;
        }

        // look at the z neighbors
        for ( int z = -1; z < 2; z += 2 )
        {
          // check that this z is inside
          if ( ((k+z) < 0) || ((k+z) >= _zdim) )
            continue;

          // find this neighbor's value
          pos.set( _scaling[0]*(float)i,
                   _scaling[1]*(float)j,
                   _scaling[2]*(float)(k+z) );
          _io_field->computeSurfacePointParams( pos, params );
          neigh_val = params._F;

          // determine if the neighbor is inside or outside
          neigh_inside = (neigh_val < 0.0);

          // check if we have a crossing
          if ( (inside && !neigh_inside) || (!inside && neigh_inside) )
            _crossing(i,j,k) = true;
        }
      }
  }

  cout << "    DONE." << endl << endl;

  // writing out the crossing file
  cout << "Writing out the crossing file..." << endl;


  if (_useNrrdIO)
  {
     char c_filename[250];
     sprintf( c_filename, "%s_crossing.nrrd", _base_filename ); 
     outputNrrdFile(c_filename,_crossing);
     cout << "Writing nrrd..." << endl;
  }
  else
  {
    char c_filename[250];
    sprintf( c_filename, "%s_crossing.vol", _base_filename );
    outputVolFile(c_filename,_crossing);
    cout << "Writing vol..." << endl;
  }


  cout << "   DONE." << endl << endl;
}

//------------------------------------------------------------------------
// Function    : initializeh0ToMaxk()
// Description : 
//------------------------------------------------------------------------
void SizingField::initializeh0ToMaxk()
{
  // check if we already have a maxk volume

//   FILE *maxk_in_file = fopen( maxk_filename, "rb" );
//   if ( maxk_in_file )
//   {
//     char yesorno;
//     do
//     {
//       cout << "maxk file exists; do you want to use it (y|n)?  ";
//       cin >> yesorno;
//     } while ( (yesorno != 'y') && (yesorno != 'n') );

//     if ( yesorno == 'y' )
//     {
//       cout << "Reading in existing maxk file..." << endl;
//       char tmp[100];
//       for ( int i = 0; i < 2; i++ )
//         fgets( tmp, 100, maxk_in_file );

//       // read in the data values
//       fread( _h0._array, sizeof(float), _xdim*_ydim*_zdim, maxk_in_file );

//       // close the files
//       fclose( maxk_in_file );
//       cout << "   DONE." << endl << endl;

//       return;
//     }

//     fclose( maxk_in_file );
//   }
  
  cout << "Computing the maxk..." << endl;

  // compute the maxk at each grid node that is a crossing and store
  //   the value in h0 --- otherwise, just store MAXFLOAT
  SurfacePointParams params;
  vec<3> pos;
  for ( int k = 0; k < _zdim; k++ )  
    for ( int j = 0; j < _ydim; j++ )
      for ( int i = 0; i < _xdim; i++ )
      {
        if ( _crossing(i,j,k) )
        {
          pos.set( _scaling[0]*(float)i,
                   _scaling[1]*(float)j,
                   _scaling[2]*(float)k );
          _io_field->computeSurfacePointParams( pos, params, true, true );

          // NOTE: we are actually computing the RADIUS of the max
          //   curvature, not the kappa value!!!!!
          _h0(i,j,k) = 1.0f/(_io_field->computeCurvature(params) + static_cast<float>(EPSILON) );
        }
        else
          _h0(i,j,k) = FLT_MAX;
      }
    
  cout << "   DONE." << endl << endl;

  // write out the maxk file
  cout << "Writing out the maxk file..." << endl;

  if (_restart_files)
  {
    if (_useNrrdIO)
    {
      char maxk_filename[250];
      sprintf( maxk_filename, "%s_maxk.nrrd", _base_filename );
      cout << "Writing .nrrd..." << endl;
      outputNrrdFile(maxk_filename,_h0);
    }
    else
    {
      char maxk_filename[250];
      sprintf( maxk_filename, "%s_maxk.vol", _base_filename );
      cout << "Writing .vol..." << endl;
      outputVolFile(maxk_filename,_h0);
    }
  }

  cout << "   DONE." << endl << endl;
}

//------------------------------------------------------------------------
// Function    : compareMaxkToLFS()
// Description : 
//------------------------------------------------------------------------
void SizingField::compareMaxkToLFS()
{
  // check if we already have a lfs volume
 //   FILE *lfs_in_file = fopen( lfs_filename, "rb" );
//   if ( lfs_in_file )
//   {
//     char yesorno;
//     do
//     {
//       cout << "lfs file exists; do you want to use it (y|n)?  ";
//       cin >> yesorno;
//     } while ( (yesorno != 'y') && (yesorno != 'n') );

//     if ( yesorno == 'y' )
//     {
//       cout << "Reading in existing lfs file..." << endl;
//       char tmp[100];
//       for ( int i = 0; i < 2; i++ )
//         fgets( tmp, 100, lfs_in_file );

//       // read in the data values
//       float lfs;
//       for ( int k = 0; k < _zdim; k++ )  
//         for ( int j = 0; j < _ydim; j++ )
//           for ( int i = 0; i < _xdim; i++ )
//           {
//             fread( &lfs, sizeof(float), 1, lfs_in_file );
//             //_h0(i,j,k) = min( _h0(i,j,k), lfs );
//             _h0(i,j,k) = lfs;
//           }

//       // close the files
//       fclose( lfs_in_file );
//       cout << "   DONE." << endl << endl;

//       return;
//     }

//     fclose( lfs_in_file );
//   }


  // create the LFS
  vec<3> start, end;
  _io_field->domain( start, end );
  char ma_pt_filename[250];
  sprintf( ma_pt_filename, "%s_ma.ptcl", _base_filename );
  LocalFeatureSize *lfs =
    new LocalFeatureSize( ma_pt_filename, start, end,
                          _xdim, _ydim, _zdim );
    
  // compute the lfs at each grid node that is a crossing and compare
  //   the value to the maxk value already stored in h0
  cerr << "FLT_MAX is : " << FLT_MAX << endl;
  array3D<float> lfs_vals(_xdim, _ydim, _zdim);
  vec<3> pos;
  for ( int k = 0; k < _zdim; k++ )  
    for ( int j = 0; j < _ydim; j++ )
      for ( int i = 0; i < _xdim; i++ )
      {
        if ( _crossing(i,j,k) )
        {
          pos.set( _scaling[0]*(float)i,
                   _scaling[1]*(float)j,
                   _scaling[2]*(float)k );
          lfs_vals(i,j,k) = lfs->lfs(pos);
          //_h0(i,j,k) = min( _h0(i,j,k), lfs_vals(i,j,k) );
          _h0(i,j,k) = lfs_vals(i,j,k);
        }
        else {
          lfs_vals(i,j,k) = FLT_MAX;
          _h0(i,j,k) = lfs_vals(i,j,k);
        }
      }

  // release the memory
  delete lfs;    
  cout << "   DONE." << endl << endl;

  // write out the lfs file
  cout << "Writing out the lfs file..." << endl;
 // FILE *lfs_out_file = fopen( lfs_filename, "w" );
//  if ( !lfs_out_file )
//  {
//    cout << "ERROR: output maxk file open failed\n";
//    return;
//  }
//
//  cout << "    vol width " << _xdim << ", height " <<
//    _ydim << ", depth " << _zdim << " ..... " << endl; 
//
//  // write out the volume size
//  fprintf( lfs_out_file, "4\n %d %d %d \n", _xdim, _ydim, _zdim );
//  fclose( lfs_out_file );
//
//  // now open in binary mode
//  if ( (lfs_out_file = fopen(lfs_filename, "ab")) == NULL )
//  {
//    cout << "ERROR: output maxk file open failed\n";
//    return;
//  }
//
//  fwrite( lfs_vals._array, sizeof(float), _xdim*_ydim*_zdim, lfs_out_file );
//  fclose( lfs_out_file );
  
  if (_restart_files)
  {
    if (_useNrrdIO)
    {
      char lfs_filename[250];
      sprintf( lfs_filename, "%s_lfs.nrrd", _base_filename );
      cout << "Writing out the lfs nrrd file..." << endl;
      outputNrrdFile(lfs_filename,lfs_vals);
    }
    else
    {
      char lfs_filename[250];
      sprintf( lfs_filename, "%s_lfs.vol", _base_filename );
      cout << "Writing out the lfs file..." << endl;
      outputVolFile(lfs_filename,lfs_vals);
    }
  }
  
  

  cout << "   DONE." << endl << endl;
}

// Utility functions to write output files.
// Nrrd functions are used by the BioMesh3D pipeline.


void SizingField::outputVolFile(char* h_filename,  array3D<float> & _h_pass)
{
 
  FILE *h_out_file = fopen( h_filename, "w" );
  if ( !h_out_file )
  {
    cout << "ERROR: output h file open failed\n";
    return;
  }

  cout << "    vol width " << _xdim << ", height " <<
    _ydim << ", depth " << _zdim << " ..... " << endl; 

  // write out the volume size
  fprintf( h_out_file, "4\n %d %d %d \n", _xdim, _ydim, _zdim );
  fclose( h_out_file );

  // now open in binary mode
  if ( (h_out_file = fopen(h_filename, "ab")) == NULL )
  {
    cerr << "ERROR: output h file open failed\n";
    return;
  }

  fwrite( _h_pass._array, sizeof(float), _xdim*_ydim*_zdim, h_out_file );
  fclose( h_out_file );
}

void SizingField::outputVolFile(char* h_filename,  array3D<bool> & _h_pass)
{
 
  FILE *h_out_file = fopen( h_filename, "w" );
  if ( !h_out_file )
  {
    cerr << "ERROR: output h file open failed\n";
    return;
  }

  cout << "    vol width " << _xdim << ", height " <<
    _ydim << ", depth " << _zdim << " ..... " << endl; 

  // write out the volume size
  fprintf( h_out_file, "4\n %d %d %d \n", _xdim, _ydim, _zdim );
  fclose( h_out_file );

  // now open in binary mode
  if ( (h_out_file = fopen(h_filename, "ab")) == NULL )
  {
    cout << "ERROR: output h file open failed\n";
    return;
  }

  fwrite( _h_pass._array, sizeof(float), _xdim*_ydim*_zdim, h_out_file );
  fclose( h_out_file );
}

void SizingField::outputNrrdFile( char* h_filename, array3D<float> & _h_pass)
{
  NrrdIO::write_nrrd(h_filename, _h_pass, _scaling[0], _scaling[1], _scaling[2], _centers[0], _centers[1], _centers[2]);
}

void SizingField::outputNrrdFile( char* h_filename, array3D<bool> & _h_pass)
{
  NrrdIO::write_nrrd(h_filename, _h_pass, _scaling[0], _scaling[1], _scaling[2], _centers[0], _centers[1], _centers[2]);
}



