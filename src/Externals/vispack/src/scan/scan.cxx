#include "scan.h"
#include "image/imagefile.h"
#include "image/image.h"
#include <string.h>
#include <limits.h>
#include <math.h>

#define TRYKEYWORDREQ(var, pref, key) { \
sprintf(scan_keyword, "%s%s", pref, keywords[key]); \
      if (VPF::set(var, pf[scan_keyword][0]) != VPF::VALID) { \
      cout << "Scan.C::myload -- missing required keyword : " << scan_keyword << endl; \
      exit(-1); } \
}
#define TRYKEYWORD(var, pref, key) { \
      sprintf(scan_keyword, "%s%s", pref, keywords[key]); \
      VPF::set(var, pf[scan_keyword][0]); \
}

//
// initialize static member variables here.....
//
char* Scan::keywords[6] = 
{
  "NUM_SCALES", "CONFIDENCE", "FILE", "NORMALS_SCALE",
  "RANGE_OFFSET", 
  "RANGE_SCALE"
};

// char* OrthoScan::keywords[2] = 
// {
//   "DELTA_X", "DELTA_Y"
// };

char *SphericalScan::keywords[NUMSPHEREPARAMS] = 
{
  "R0", "C0", "DELTA_THETA", "DELTA_PHI", 
};

char *SphericalAmplitudeRangeScan::keywords[2] = 
{
  "SPATIALUNCERTAINTY", "NOISELEVEL"
};

Scan::Scan(const VISImage<float> &range_map)
{
  //  	cout << "have entered scan create in" << endl;
  _range_map_orig = transformRange(range_map);
  _range_map_orig_dx = (_range_map_orig.dx()).setBorder(0.0f, 1);
  _range_map_orig_dy = (_range_map_orig.dy()).setBorder(0.0f, 1);
  _num_scales = 0;
  _scale_factor = NULL;
//	_dx_kernel = dx_kernel(1);
//	_dy_kernel = dy_kernel(1);
}

Scan::Scan(const VISImage<float> &range_map, float normals_scale)
{
  //	cout << "have entered scan create in" << endl;
  _range_map_orig = transformRange(range_map);
  VISImage<float> gauss_im = _range_map_orig.gaussDiffuse(normals_scale);
   // do you smooth the map or not
  //   _range_map_orig = gauss_im;
   // must fix this bug in dx and convolve!!!!!!!!!!!!!!!!!!!!!!!1
   _range_map_orig_dx = -1.0f*gauss_im.dx();
   int i, w = gauss_im.width(), h = gauss_im.height();
   for (i = 0; i < h; i++)
     {
       _range_map_orig_dx.poke(0, i) = (gauss_im.peek(1, i)  - gauss_im.peek(0, i))/2.0f;
       _range_map_orig_dx.poke(w-1, i) = (gauss_im.peek(w-1, i)  - gauss_im.peek(w - 2, i))/2.0f;
     }
    
   // must fix this bug in dx and convolve!!!!!!!!!!!!!!!!!!!!!!!1
   _range_map_orig_dy = -1.0f*gauss_im.dy();
   for (i = 0; i < w; i++)
     {
       _range_map_orig_dy.poke(i, 0) = (gauss_im.peek(i, 1)  - gauss_im.peek(i, 0))/2.0f;
       _range_map_orig_dy.poke(i, h-1) = (gauss_im.peek(i, h-1)  - gauss_im.peek(i, h - 2))/2.0f;
     }
  //	_dx_kernel = dx_kernel(1);
  //	_dy_kernel = dy_kernel(1);
}

void Scan::initRangeImages(const VISImage<float> &range_map,  float normals_scale)
{
  //  cout << "about to transform range" << endl;
  _range_map_orig = transformRange(range_map);
  VISImage<float> gauss_im = _range_map_orig.gaussDiffuse(normals_scale);
   // do you smooth the map or not
  //   _range_map_orig = gauss_im;
  // must fix this bug in dx and convolve!!!!!!!!!!!!!!!!!!!!!!!1
  _range_map_orig_dx = -1.0f*gauss_im.dx();
  int i, w = gauss_im.width(), h = gauss_im.height();
   for (i = 0; i < h; i++)
     {
       _range_map_orig_dx.poke(0, i) = (gauss_im.peek(1, i)  
					- gauss_im.peek(0, i))/2.0f;
       _range_map_orig_dx.poke(w-1, i) = (gauss_im.peek(w-1, i)  
					  - gauss_im.peek(w - 2, i))/2.0f;
     }
    
   // must fix this bug in dx and convolve!!!!!!!!!!!!!!!!!!!!!!!1
   _range_map_orig_dy = -1.0f*gauss_im.dy();
   for (i = 0; i < w; i++)
     {
       _range_map_orig_dy.poke(i, 0) = (gauss_im.peek(i, 1)  - gauss_im.peek(i, 0))/2.0f;
       _range_map_orig_dy.poke(i, h-1) = (gauss_im.peek(i, h-1)  - gauss_im.peek(i, h - 2))/2.0f;
     }
   createScales(_num_scales);
}

void Scan::readParams(VPF::ParameterFile pf, const char* prefix)
{
  // must read depth image, number of scales, confidience
  char scan_keyword[80];
  char filename[80];
  VISImageFile im_file;
  VISImage<float> range_map;

  TRYKEYWORD(_num_scales, prefix, NUMSCALES);
  TRYKEYWORD(filename, prefix, RANGEIMAGE);
  TRYKEYWORD(_r_offset, prefix, ROFFSET);
  TRYKEYWORD(_r_delta, prefix, RDELTA);
  TRYKEYWORD(filename, prefix, RANGEIMAGE);
  TRYKEYWORD(_normals_scale, prefix, NORMALSSCALE);
  if (!(range_map = VISImage<float>(im_file.read(filename))).isValid())
    cout << "Scan--readParams:  bad depth file: " << filename << endl;
  //  cout << "normals_scale: " << _normals_scale << endl;
  initRangeImages(range_map, _normals_scale);
  TRYKEYWORD(filename, prefix, CONFIDENCEIMAGE);
  VISImage<float> conf;
  if (!(conf = VISImage<float>(im_file.read(filename))).isValid())
    cout << "Scan--readParams:  bad conf file" << endl;
  else
    {
      //      cout << "about to set conf from read params on Scan" << endl;
      setConfidence(conf);
    }
}


void Scan::createScales(int num_scales)
{
    int k;
    _num_scales = num_scales;
    _range_maps.clear();
    _range_maps_dx.clear();
    _range_maps_dy.clear();
    _range_maps.appendItem(_range_map_orig);
    _range_maps_dx.appendItem(_range_map_orig_dx);
    _range_maps_dy.appendItem(_range_map_orig_dy);
    _scale_factor = new float[num_scales];
    _scale_factor[0] = 1.0f;
    char filename[80];
    VISImageFile im_file;
    for (k = 1; k < num_scales; k++)
	{
	    _range_maps.appendItem(subSampleAverage((_range_maps.peek(k - 1))));
	    _range_maps_dx.appendItem(subSampleAverage((_range_maps_dx.peek(k - 1))));
	    _range_maps_dy.appendItem(subSampleAverage((_range_maps_dy.peek(k - 1))));
	    _scale_factor[k] = 2.0f*_scale_factor[k - 1];
	    	    sprintf(filename, "range_dy_sub_%d.fit", k);
	    	    im_file.write(_range_maps_dy.peek(k), filename);
	    	    sprintf(filename, "range_sub_%d.fit", k);
	    	    im_file.write(_range_maps.peek(k), filename);
	}
    createScalesConf(_num_scales);
}

void Scan::createScalesConf(int num_scales)
{

  char filename[80];
  VISImageFile im_file;

  _conf_maps.clear();
  _conf_maps.appendItem(_conf_map);
  _conf_maps_dx.appendItem(_conf_map_dx);
  _conf_maps_dy.appendItem(_conf_map_dy);
  if (_conf_map.isValid())
    for (int k = 1; k < num_scales; k++)
      {
	_conf_maps.appendItem(subSampleAverage((_conf_maps.peek(k - 1))));
	_conf_maps_dx.appendItem((_conf_maps.peek(k)).dx());
	_conf_maps_dy.appendItem((_conf_maps.peek(k)).dy());
      }
  else
    for (int k = 1; k < num_scales; k++)
      {
	_conf_maps.appendItem(VISImage<float>());
	_conf_maps_dx.appendItem(VISImage<float>());
	_conf_maps_dy.appendItem(VISImage<float>());
      }
}


VISImage<float> subSampleAverage(VISImage<float> im)
{
    int w = im.width(), h = im.height();
    int new_w, new_h;

    if (!(im.isValid()))
      return VISImage<float>();

    VISImage<float> 
	r((new_w = w/2) + w%2, 
	  (new_h = h/2) + h%2),
	tmp(new_w + w%2, h);
    
    int i, j, this_x, this_y;
    
    for (j = 0; j < h; j++)
	{
	    tmp.poke(0, j) = 0.5f*(im.peek(0, j) 
				   + im.peek(1, j));
	    
	    for (i = 1; i < new_w; i++)
		{
		    this_x = 2*i;
		    tmp.poke(i, j) 
			= 0.25f*(im.peek(this_x - 1, j) + 
				 im.peek(this_x + 1, j)) + 
			0.5f*im.peek(this_x, j);
		}

	    if (w%2)
		tmp.poke(new_w, j) 
		    = 0.5f*(im.peek(w - 1, j) + im.peek(w - 2, j));
	}

    new_w = new_w + w%2;

    if (h < 2)
	{
	    return(tmp);
	}
    
    for (i = 0; i < new_w; i++)
	{
	    r.poke(i, 0) = 0.5f*(tmp.peek(i, 0) 
				   + tmp.peek(i, 1));
	    
	    for (j = 1; j < new_h; j++)
		{
		    this_y = 2*j;
		    r.poke(i, j) 
			= 0.25f*(tmp.peek(i, this_y - 1) + 
				 tmp.peek(i, this_y + 1)) + 
			0.5f*tmp.peek(i, this_y);
		}

	    if (h%2)
		r.poke(i, new_h) 
		    = 0.5f*(tmp.peek(i, h - 1) + tmp.peek(i, h - 2));
	}
    
    return(r);
}

VISImage<float> subSampleMin(VISImage<float> im)
{
  int w = im.width(), h = im.height();
  int new_w, new_h;

  VISImage<float> 
    r((new_w = w/2) + w%2, 
      (new_h = h/2) + h%2),
    tmp(new_w + w%2, h);
    
  int i, j, this_x, this_y;
    
  for (j = 0; j < h; j++)
    {
      tmp.poke(0, j) = VISmin(im.peek(0, j), im.peek(1, j));
	    
      for (i = 1; i < new_w; i++)
	{
	  this_x = 2*i;
	  tmp.poke(i, j) 
	    = VISmin(VISmin(im.peek(this_x - 1, j), im.peek(this_x + 1, j)), 
		     im.peek(this_x, j));
	}

      if (w%2)
	tmp.poke(new_w, j) 
	  = VISmin(im.peek(w - 1, j), im.peek(w - 2, j));
    }

  new_w = new_w + w%2;

  if (h < 2)
    {
      return(tmp);
    }
    
  for (i = 0; i < new_w; i++)
    {
      r.poke(i, 0) = VISmin(tmp.peek(i, 0), tmp.peek(i, 1));
	    
      for (j = 1; j < new_h; j++)
	{
	  this_y = 2*j;
	  r.poke(i, j) 
	    = VISmin(VISmin(tmp.peek(i, this_y - 1), tmp.peek(i, this_y + 1)),  tmp.peek(i, this_y));
	}

      if (h%2)
	r.poke(i, new_h) 
	  = VISmin(tmp.peek(i, h - 1), tmp.peek(i, h - 2));
    }
    
  return(r);
}


// *****************************************************
// *****************************************************
// ****************** TwoDOrthoScane *******************
// *****************************************************
// *****************************************************


//*********************************************************************
//*********************************************************************
//********************** 3D Spherical Scan ****************************
//*********************************************************************
//*********************************************************************


SphericalScan::SphericalScan(const VISImage<float> &range_image):Scan(range_image)
{
//    cout << "about to enter sphere create in" << endl;
  setDefaultParams();
  VISImage<float> conf(range_image.width(), range_image.height());
  conf = 1.0f;
  //  cout << "about to set conf from range const in sphere" << endl;
  setConfidence(conf);
}

SphericalScan::SphericalScan(const VISImage<float> &range_image, float normals_scale)
  :Scan(range_image, normals_scale)
{
//    cout << "about to enter sphere create in" << endl;
  setDefaultParams();
  VISImage<float> conf(range_image.width(), range_image.height());
  conf = 1.0f;
  //  cout << "about to set conf from range + scale const in sphere" << endl;
  setConfidence(conf);
}


VISVector TwoDOrthoScan::getSurfaceNormal(int i, int j) const
{
    VISVector r(4);
    float di = _range_map_orig_dx(i, j);
    float dj = _range_map_orig_dy(i, j);
    float range = _range_map_orig.peek(i, j);

    float norm = 1.0f/sqrt(1.0f + pow(di,2) * pow(dj,2));

    r.poke(0) = di;
    r.poke(1) = dj;
    r.poke(2) = -1.0f;
    r.poke(3) = 0.0f;

    return(r);
}

VISVector TwoDOrthoScan::getSurfaceNormal(int i, int j, int scale) const
{

    VISVector r(4);
    //    VISImage<float> range_map = _range_maps.peek(scale);
    float di = (_range_maps_dx.peek(scale)).peek(i, j);
    float dj = (_range_maps_dy.peek(scale)).peek(i, j);
    float this_scale = _scale_factor[scale];
    float range = (_range_maps.peek(scale)).interp(i, j);

    float norm = 1.0f/sqrt(1.0f + pow(di,2) * pow(dj,2));

    r.poke(0) = di;
    r.poke(1) = dj;
    r.poke(2) = -1.0f;
    r.poke(3) = 0.0f;

    return(r);
}



float TwoDOrthoScan::depth(const VISVector &p, int scale) const
{
  float x = p.peek(0)/_scale_factor[scale];
  VISImage<float> this_map =  _range_maps.peek(scale);
  if (!this_map.checkBounds(x, 0.0f))
	{
//	    cout << "check bounds bad " << p << "scale " 
//		 << scale << " x " << x << endl;
	    return(errorCode());
	}
    else
	return(this_map.interp(x, 0.0f));
}

VISVector TwoDOrthoScan::grad_depth(const VISVector &p, int scale) const
{
  float x = p.peek(0)/_scale_factor[scale];
  VISImage<float> this_map_dx =  _range_maps_dx.peek(scale);
  //    cout << "grad depth at " << x << " is " << 
  //	this_map_dx.interp(x, 0.0f) << endl;
  if (!this_map_dx.checkBounds(x, 0.0f))
    {
      //	    cout << "check bounds bad " << p << "scale " 
      //		 << scale << " x " << x << endl;
      return VISVector(0.0f, 0.0f, 0.0f);
    }
  else
    return VISVector(this_map_dx.interp(x, 0.0f), 0.0f, 0.0f);
}

float TwoDOrthoScan::depth(const VISVector &p) const
{
    float x = p.peek(0);
    if (!_range_map_orig.checkBounds(x, 0.0f))
	return(errorCode());
    else
	return(_range_map_orig.interp(x, 0.0f));
}

VISVector TwoDOrthoScan::lineOfSight(const VISVector &p) const
{
  return(VISVector(0.0f, -1.0f, 0.0f));
}

VISVector TwoDOrthoScan::lineOfSight(int i, int j) const
{
  return(VISVector(0.0f, -1.0f, 0.0f));
}

//float TwoDOrthoScan::distance(const VISVector &p)
//{
//    return(p.peek(0));
//}



VISVector TwoDOrthoScan::get3DPoint(int i, int j) const
{
    return(VISVector((float)i, 
		     _range_map_orig.peek(i, 0), 
		     1.0f));
}

VISVector TwoDOrthoScan::get3DPoint(int i, int j, int scale) const
{
    return(VISVector((float)i, (_range_maps.peek(scale))
		     .interp((float)i/_scale_factor[scale], 0), 
		     1.0f));
}


//*********************************************************************
//*********************************************************************
//********************** 3D Spherical Scan ****************************
//*********************************************************************
//*********************************************************************

void SphericalScan::setDefaultParams()
{
  setC0(_range_map_orig.width()/2.0f);
  setR0(_range_map_orig.height()/2.0f);
  setDeltaTheta(DEFAULT_FOV/_range_map_orig.width()); 
  setDeltaPhi(DEFAULT_FOV/_range_map_orig.height());
  //  cout << "about to set conf from set Def params" << endl;
  //  setConfidence(VISImage<float>());
}

void SphericalScan::readParams(VPF::ParameterFile pf, const char* prefix)
{
//       setR0(params.peek(SphereScan::R0)); 
//       setC0(params.peek(SphereScan::C0)); 
//       setDeltaPhi(params.peek(SphereScan::DPHI)); 
//       setDeltaTheta(params.peek(SphereScan::DTHETA)); 
//       setRangeDelta(params.peek(SphereScan::RDELTA)); 
//       setRangeOffset(params.peek(SphereScan::ROFFSET)); 
  char scan_keyword[80];
  // get the old params, this changes them, but does nothing if
  // the file contains to entry for that parameter
  VISVector params = getParams();
  int i;
  // all of the additional parameters are floats, so you can do this
  // must keep track of (6) parameters
  float float_tmp;
  for (i = 0; i < NUMSPHEREPARAMS; i++)
    {
      TRYKEYWORD(float_tmp, prefix, i);
      params.poke(i) = float_tmp;
    }
  setParams(params);
}

void SphericalScan::setParams(const VISVector &params)
{
  if (params.n() == NUMSPHEREPARAMS)
    {
      setR0(params.peek(R0)); 
      setC0(params.peek(C0)); 
      setDeltaPhi(params.peek(DPHI)); 
      setDeltaTheta(params.peek(DTHETA)); 
      //      setRangeDelta(params.peek(RDELTA)); 
      //      setRangeOffset(params.peek(ROFFSET)); 
    }
}

VISVector SphericalScan::getParams() const
{
  VISVector r(NUMSPHEREPARAMS);
  r.poke(R0) = _r0; 
  r.poke(C0) = _c0; 
  r.poke(DPHI) = _delta_phi; 
  r.poke(DTHETA) = _delta_theta; 
  //  r.poke(RDELTA) = _r_delta; 
  //  r.poke(ROFFSET) = _r_offset; 
  return(r);
}

VISVector SphericalScan::getSurfaceNormal(int i, int j) const
{
    VISVector r(4);
    //    VISImage<float> range_map = _range_map_orig;
    float di = _r_delta*(_range_map_orig_dx.peek(i, j));
    float dj = _r_delta*(_range_map_orig_dy.peek(i, j));
    //    cout << "getSurfaceNormal " << i << " " << j  << " "<< di  << " "<< dj << endl;
    float cos_theta, sin_theta, cos_phi, sin_phi, 
	dR_dPhi, dR_dTheta;

    float range = _range_map_orig.peek(i, j);

    float phi = (j - _r0)*_delta_phi,
      theta = (i - _c0)*_delta_theta;

    //    cout << "_delta " << _delta_phi << " " << _delta_theta  <<  endl;
    //    cout << "getSurfaceNormal " << phi << " " << theta  << " "<< range  <<  endl;

    cos_theta = cos(theta);
    sin_theta = sin(theta);
    cos_phi = cos(phi);
    sin_phi = sin(phi);


    dR_dTheta = di/_delta_theta;
    dR_dPhi = dj/_delta_phi;

    r.poke(0) = cos_theta*(range*sin_theta - dR_dTheta*cos_theta);
    r.poke(1) = cos_theta*sin_phi*(range*cos_theta + dR_dTheta*sin_theta)
	- dR_dPhi*cos_theta;
    r.poke(2) = dR_dPhi*sin_phi
	+ cos_phi*cos_theta*(range*cos_theta + dR_dTheta*sin_theta);
    r.poke(3) = 0.0f;

// should this be positive or negative???
    r /= -1.0f*r.norm();
    //    r /= r.norm();
//    cout << "normal length test " << r.dot(r) << endl;
    return(r);
}

VISVector SphericalScan::getSurfaceNormal(int i, int j, int scale) const
{

    VISVector r(4);
    float this_scale = _scale_factor[scale];

    float di = _r_delta*(_range_maps_dx.peek(scale)).interp(i, j);
    float dj = _r_delta*(_range_maps_dy.peek(scale)).interp(i, j);
    float cos_theta, sin_theta, cos_phi, sin_phi, 
	dR_dPhi, dR_dTheta, norm;

    float range = (_range_maps.peek(scale)).interp(i, j);

    float phi = (this_scale*j - _r0)*_delta_phi,
      theta = (this_scale*i - _c0)*_delta_theta;

    cos_theta = cos(theta);
    sin_theta = sin(theta);
    cos_phi = cos(phi);
    sin_phi = sin(phi);


    dR_dTheta = di/_delta_theta;
    dR_dPhi = dj/_delta_phi;

//    norm = 1.0f/sqrt(pow(dR_dTheta*cos_theta, 2) 
//		     + pow(range*cos_theta, 2) 
//		     + pow(dR_dPhi, 2));

    r.poke(0) = cos_theta*(range*sin_theta - dR_dTheta*cos_theta);
    r.poke(1) = cos_theta*sin_phi*(range*cos_theta + dR_dTheta*sin_theta)
      - dR_dPhi*cos_theta;
    r.poke(2) = dR_dPhi*sin_phi
      + cos_phi*cos_theta*(range*cos_theta + dR_dTheta*sin_theta);
    r.poke(3) = 0.0f;

//
// this is what you get if you have (x, y, z) = r*(sin_theta*cos_phi, sin_phi, 
// cos_phi*cos_theta)
//
//    r.poke(0) = cos_phi*sin_theta*(range*cos_phi + dR_dPhi*sin_phi)
//	- dR_dTheta*cos_theta;
//    r.poke(1) = cos_phi*(range*sin_phi - dR_dPhi*cos_phi);
//    r.poke(2) = dR_dTheta*sin_theta
//	+ cos_phi*cos_theta*(range*cos_phi + dR_dPhi*sin_phi);
//    r.poke(3) = 0.0f;

// should this be positive or negative???
    r /= -1.0f*r.norm();
//    r /= r.norm();

//    cout << "normal length test " << r.dot(r) << endl;

    return(r);

}

VISVector SphericalScan::myImageCoord(const VISVector &p) const
{
  float x = p.peek(0), y = p.peek(1), z = p.peek(2);
  float theta, phi;
  float r, c;

  if (z > 0)
    {
      theta = ::atan(x/sqrt(y*y + z*z));
      phi = ::atan(y/z);
      c = (theta/_delta_theta + _c0);
      r = (phi/_delta_phi + _r0);
      //		if ((r >= 0.0)&&(c >= 0.0)&&(c <= (_depth.width() - 1))
      //		    &&(r <= (_depth.height() - 1)))
      return(VISVector(c, r));
    }
  else
    {
      //	    cout << "imageCoord: ERROR " << p << endl;
      return VISVector();
    }
}


VISVector SphericalScan::imageCoord(const VISVector &p, int scale) const
{
    VISVector coord = imageCoord(p);    
    return(coord/_scale_factor[scale]);
}


float SphericalScan::depth(const VISVector &p, int scale) const
{
  VISVector coord = imageCoord(p);
  float c, r;
  VISImage<float> this_map;
    
  if (coord.isValid()&&
      ((this_map = _range_maps.peek(scale))
       .checkBounds(c = (coord.peek(0)/
			 _scale_factor[scale]),
		    r = (coord.peek(1)/
			 _scale_factor[scale]))))
    {
      //      cout << "in Scan::depth " << endl;
      //      cout << "point " << p << endl;
      //             (this_map).print();
      //      cout << "coord " << coord << endl;
      //      cout << "c " << c << " r " << r << endl;
      //             cout << "depth " << this_map.VISImage<float>::interp(c, r) << endl;
      return(_r_offset + _r_delta*this_map.VISImage<float>::interp(c, r));
    }
  else 		
    return(errorCode());
}

float SphericalScan::depth(const VISVector &p) const
{
  // geometry of the range map is defined here.....
  VISVector coord = imageCoord(p);
  float c, r;


  //  cout << "in Scan::depth " << endl;
  //  cout << "point " << p << endl;
  //             (this_map).print();
  //  cout << "coord " << coord << endl;
  //             cout << "depth " << this_map.VISImage<float>::interp(c, r) << endl;
  //  cout << "depth image coord " << coord << endl;
    
  if (coord.isValid()&&(_range_map_orig
			.checkBounds(c = coord.peek(0),r = coord.peek(1))))
    return(_r_offset + _r_delta*_range_map_orig.VISImage<float>::interp(c, r));
  else 		
    return(errorCode());
}


VISVector SphericalScan::grad_depth(const VISVector &p, int scale) const
{

  VISVector coord = myImageCoord(p);
  float c, r;
  VISImage<float> this_map, im_dx, im_dy;
  float dc, dr;
  float x = p.peek(0), y = p.peek(1), z = p.peek(2);

  if (!(coord.isValid()&&
	((im_dx  = _range_maps_dx.peek(scale))
	 .checkBounds(c = (coord.peek(0)/_scale_factor[scale]),
		      r = (coord.peek(1)/_scale_factor[scale])))))
    return VISVector(0.0f, 0.0f, 0.0f, 0.0f);
  im_dy = _range_maps_dy.peek(scale);

  dc = _r_delta*im_dx.interp(c, r);
  dr = _r_delta*im_dy.interp(c, r);

  //  cout << "dc " << dc << " and dr " << dr << endl;

  float dtheta_dx, dtheta_dy, dtheta_dz, dphi_dx, dphi_dy, dphi_dz;
  float mag1 = y*y + z*z;
  float mag2 = mag1 + x*x;
  float mag1_sqrt = sqrt(mag1);

  dtheta_dx = mag1/(mag2*mag1_sqrt);
  dtheta_dy = -x*y/(mag2*(mag1_sqrt));
  dtheta_dz = -x*z/(mag2*(mag1_sqrt));
  dphi_dx = 0.0f;
  dphi_dy = z/(mag1);
  dphi_dz = -y/(mag1);

  return VISVector
    (
     dtheta_dx*dc/_delta_theta, 
     dtheta_dy*dc/_delta_theta + dphi_dy*dr/_delta_phi,
     dtheta_dz*dc/_delta_theta + dphi_dz*dr/_delta_phi, 
     0.0f
     );
}

VISVector SphericalScan::grad_conf(const VISVector &p, int scale) const
{

  VISVector coord = myImageCoord(p);
  float c, r;
  VISImage<float> this_map, im_dx, im_dy;
  float dc, dr;
  float x = p.peek(0), y = p.peek(1), z = p.peek(2);

  if (!(coord.isValid()&&
	((im_dx  = _conf_maps_dx.peek(scale))
	 .checkBounds(c = (coord.peek(0)/_scale_factor[scale]),
		      r = (coord.peek(1)/_scale_factor[scale])))))
    return VISVector(0.0f, 0.0f, 0.0f, 0.0f);
  im_dy = _conf_maps_dy.peek(scale);

  dc = _r_delta*im_dx.interp(c, r);
  dr = _r_delta*im_dy.interp(c, r);

  //  cout << "dc " << dc << " and dr " << dr << endl;

  float dtheta_dx, dtheta_dy, dtheta_dz, dphi_dx, dphi_dy, dphi_dz;
  float mag1 = y*y + z*z;
  float mag2 = mag1 + x*x;
  float mag1_sqrt = sqrt(mag1);

  dtheta_dx = mag1/(mag2*mag1_sqrt);
  dtheta_dy = -x*y/(mag2*(mag1_sqrt));
  dtheta_dz = -x*z/(mag2*(mag1_sqrt));
  dphi_dx = 0.0f;
  dphi_dy = z/(mag1);
  dphi_dz = -y/(mag1);

  return VISVector
    (
     dtheta_dx*dc/_delta_theta, 
     dtheta_dy*dc/_delta_theta + dphi_dy*dr/_delta_phi,
     dtheta_dz*dc/_delta_theta + dphi_dz*dr/_delta_phi, 
     0.0f
     );
}


float SphericalScan::confidence(const VISVector &p) const
{
  VISVector coord = imageCoord(p);
  float c, r;

  //  cout << "image coord" << coord << endl;
  //  cout << "point " << p << endl;
  if (_conf_map.isValid())  
    if ((_conf_map.checkBounds(c = coord.peek(0),r = coord.peek(1))))
	return(_conf_map.VISImage<float>::interp(c, r));
    else 		
      return(0.0f);
  else
    {
      return(1.0f);
    }
}

float SphericalScan::confidence(const VISVector &p, int scale) const
{
  //  **************
  //    should we subsample confidence ????
  //  **************
  return(confidence(p));
  
//    VISVector coord = myImageCoord(p);
//    float c, r;
//    VISImage<float> this_map;



//    if (coord.isValid()&&((this_map = _conf_maps.peek(scale))
//   			.checkBounds(c = (coord.peek(0)/_scale_factor[scale]),
//   				     r = (coord.peek(1)/_scale_factor[scale]))))
//      {
//        if (this_map.VISImage<float>::interp(c, r) < 0.5)
//  	cout << this_map.VISImage<float>::interp(c, r)  << 
//  	  "got one less than 0.5 " << p << " " << r << " " << c << endl;
    
//        else
//  	cout << "got one greater ";
//        return(this_map.VISImage<float>::interp(c, r));
//      }
//    else 		
//      return(0.0f);

}

VISVector SphericalScan::lineOfSight(const VISVector &p) const
{
  float mag, z = p[2], y = p[1], x = p[0];
  mag = sqrt(x*x + y*y + z*z);

  //  cout << "sphere los " << p << endl;

  if (depth(p) == errorCode())
    return(VISVector(0.0f, 0.0f, 0.0f, 0.0f));

  if (mag > 0.0)
    mag = -1.0f/mag;
  return(VISVector(x*mag, y*mag, z*mag, 0.0f));
}

VISVector SphericalScan::lineOfSight(int i, int j) const
{
  float 
    phi = (j - _r0)*_delta_phi, 
    theta = (i - _c0)*_delta_theta;
  return(VISVector(-sin(theta), -sin(phi)*cos(theta), 
		   -cos(phi)*cos(theta), 0.0f)); 
}

VISVector SphericalScan::get3DPoint(int i, int j) const
{
  float 
    phi = (j - _r0)*_delta_phi, 
    theta = (i - _c0)*_delta_theta,
    r = _r_offset + _r_delta*_range_map_orig.peek(i, j);
  cout << "map " << _range_map_orig.peek(i, j) << " times " << _r_delta << " + " << "_r_offset " << _r_offset << 
    " is " << r << endl;
  cout << "point is " << endl << VISVector(r*sin(theta), r*sin(phi)*cos(theta), r*cos(phi)*cos(theta), 1.0f) << endl;
  return(VISVector(r*sin(theta), r*sin(phi)*cos(theta), 
		   r*cos(phi)*cos(theta), 1.0f)); 
}

VISVector SphericalScan::get3DPoint(int i, int j, int scale) const
// should this assume i, j are scaled or not?
// for now, yes
{
    float this_scale = _scale_factor[scale];
    float phi = (this_scale*j - _r0)*_delta_phi;
    float theta = (this_scale*i - _c0)*_delta_theta;    
//    cout << "this scale factor " << this_scale << endl;
    float r = _r_offset + _r_delta*(_range_maps.peek(scale)).peek(i, j);

    //    cout << "phi " << _delta_phi << " " << j << " " << _r0 << endl;
    //    cout << "theta " << _delta_theta << " " << i << " " << _c0 << endl;
    //    cout << "phi " << phi << "theta " << theta << "r " << r << endl;

    //  cout << "map " << _range_map_orig.peek(i, j) << " times " << _r_delta << " + " << "_r_offset " << _r_offset << 
    //    " is " << r << endl;
    //  cout << "point is " << endl << VISVector(r*sin(theta), r*sin(phi)*cos(theta), r*cos(phi)*cos(theta), 1.0f) << endl;

    return(VISVector(r*sin(theta), r*sin(phi)*cos(theta), 
		     r*cos(phi)*cos(theta), 1.0f)); 
}


VISArray< VISImage<float> > Scan
::get3DPoints(const VISMatrix &transform) const
{
  return get3DPoints(transform, 0);
}

VISArray< VISImage<float> > Scan::get3DPoints(const VISMatrix &transform, int scale) const
{
  VISImage<float> range = _range_maps.peek(scale);
  int w = range.width(), h = range.height();
  VISImage<float>
    x(w - 2*BORDER, h - 2*BORDER), 
    y(w - 2*BORDER, h - 2*BORDER), 
    z(w - 2*BORDER, h - 2*BORDER);
  

  VISVector p(4);
  p.poke(3) = 1.0f;
  int i, j;

  for (i = 0; i < (w - 2*BORDER); i++)
    for (j = 0; j < (h - 2*BORDER); j++)
      {
	//	    cout << "about to do get 3D point, scale = " << scale << endl;
	p.pokeROI(0, get3DPoint(i + BORDER, j + BORDER, scale));
	p = transform*p;
	x.poke(i, j) = p.peek(0);
	y.poke(i, j) = p.peek(1);
	z.poke(i, j) = p.peek(2);
      }

  VISArray< VISImage<float> > r;

  r.poke(0) = x;
  r.poke(1) = y;
  r.poke(2) = z;

  return(r);
}


//*********************************************************************
//*********************************************************************
//************************ 3D Ortho Scan ****************************
//*********************************************************************
//*********************************************************************

VISVector OrthoScan::getSurfaceNormal(int i, int j) const
{
    VISVector r(4);
    r.poke(0) = _range_map_orig_dx.peek(i, j);
    r.poke(1) = _range_map_orig_dy.peek(i, j);
    r.poke(2) = -1.0f;
    r.poke(3) = 0.0f;
    r /= r.norm();
    return(r);
}

VISVector OrthoScan::getSurfaceNormal(int i, int j, int scale) const
{

    VISVector r(4);
    float this_scale = _scale_factor[scale];
    r.poke(0) = (_range_maps_dx.peek(scale)).interp(i, j);
    r.poke(1) = (_range_maps_dy.peek(scale)).interp(i, j);
    r.poke(2) = -1.0f;
    r.poke(3) = 0.0f;
    r /= r.norm();

    return(r);

}

VISVector OrthoScan::myImageCoord(const VISVector &p) const
{
  float x = p.peek(0), y = p.peek(1), z = p.peek(2);
  float theta, phi;
  float r, c;

  if (z > 0)
    {
      c = x/_delta_x + _c0;
      r = y/_delta_y + _r0;
      return(VISVector(c, r));
    }
  else
    {
      return VISVector();
    }
}


VISVector OrthoScan::grad_depth(const VISVector &p, int scale) const
{
  VISVector coord = myImageCoord(p);
  float c, r;
  VISImage<float> this_map, im_dx, im_dy;
  float dc, dr;
  float x = p.peek(0), y = p.peek(1), z = p.peek(2);

  if (!(coord.isValid()&&
	((im_dx  = _range_maps_dx.peek(scale))
	 .checkBounds(c = (coord.peek(0)/_scale_factor[scale]),
		      r = (coord.peek(1)/_scale_factor[scale])))))
    {
      //  cout << "got invalid grad_depth" << p << coord <<endl;
      return VISVector(0.0f, 0.0f, 0.0f, 0.0f);
    }
  //  else
  //    cout << "got valid grad_depth" << p << coord <<endl;
  im_dy = _range_maps_dy.peek(scale);

  dc = im_dx.interp(c, r);
  dr = im_dy.interp(c, r);

  return VISVector
    (
     dc/_delta_x, 
     dr/_delta_y, 
     0.0f, 
     0.0f
     );
}


VISVector OrthoScan::get3DPoint(int i, int j) const
{
    float 
      y = (j - _r0)*_delta_y, 
      x = (i - _c0)*_delta_x, 
      z = _range_map_orig.peek(i, j);
    return(VISVector(x, y, z, 1.0f)); 
}

VISVector OrthoScan::get3DPoint(int i, int j, int scale) const
// should this assume i, j are scaled or not?
// for now, yes
{
    float this_scale = _scale_factor[scale];
    float 
      y = (this_scale*j - _r0)*_delta_y, 
      x = (this_scale*i - _c0)*_delta_x, 
      z = (_range_maps.peek(scale)).peek(i, j);
    return(VISVector(x, y, z, 1.0f)); 
}

//
// Need (only) for fitting models
//
// VISMatrix SphericalScan::dx_dparams(int c, int r) const
// {
//   // 4 coords (homo),  4 free parameters
//   VISMatrix ret(4, 4);

//   float 
//     phi = (r - _r0)*_delta_phi,
//     theta = (c - _c0)*_delta_theta,
//     sin_theta = sin(theta),
//     cos_theta = cos(theta),
//     sin_phi = sin(phi),
//     cos_phi = cos(phi),
//     range_raw, range;
  
//   range_raw = _range_map_orig.peek(c, r);
//   range = _r_offset + _r_delta*range_raw;

//     ret.pokeROI(0, Calib::ROFFSET - 7, 
//     VISVector(sin_theta, sin_phi*cos_theta, 
//   	    cos_phi*cos_theta, 0.0f));
//     //ret.pokeROI(0, Calib::ROFFSET - 7, VISVector(0.0f, 0.0f, 0.0f, 0.0f));
//   ret.pokeROI(0, Calib::RDELTA - 7, 
// 	      range_raw*VISVector(sin_theta, sin_phi*cos_theta, 
// 				  cos_phi*cos_theta, 0.0f));
//   ret.pokeROI(0, Calib::DTHETA - 7, ((c - _c0)*range)
// 	      *VISVector(cos_theta, -sin_phi*sin_theta, 
// 			 -cos_phi*sin_theta, 0.0f));
//   ret.pokeROI(0, Calib::DPHI - 7, ((r - _r0)*range)
// 	      *VISVector(sin_theta,  cos_phi*cos_theta, 
// 			 -sin_phi*cos_theta, 0.0f));
//   return(ret);
// }


// *****************************************************
// *****************************************************
// ****** SphericalAmplitudeRangeScan  ***********
// *****************************************************
// *****************************************************
void SphericalAmplitudeRangeScan::readParams(VPF::ParameterFile pf, const char* prefix)
{
  // must read in various parameters
  cout << "about to set conf from read params in sphAmpRangeScan" << endl;
  setConfidence(createConfidence());
}

VISImage<float> SphericalAmplitudeRangeScan::createConfidence() const
{
    VISImage<float> conf = _range_map_orig.createToSize();
    int w = _range_map_orig.width();
    int h = _range_map_orig.height();
    int r, c;
    //    float dr, dc;
//    float x, y, z;
    VISVector d_depth;
    VISVector p(3);
    float image_factor;
    
    float d;
    float factor;

    //where do dx and dy's get created?  RTW

// this is the OA coordinate system - see page 312, "Surfaces in Range Image 
// Understanding" by Besl 

    for (r = 0; r < h; r++)
	for (c = 0; c < w; c++) // here previously was c<h
	    {
	      p = depth(c, r)*lineOfSight(c, r);
	      d_depth = grad_depth(p, 0);
	      image_factor = d_depth.norm();
	      conf.at(c, r) = 1.0f/
		((1.0f + image_factor*image_factor)
		 *_spatial_uncertainty
		 *_spatial_uncertainty
		 + _noise_level*_noise_level 
// put this in to account for discretization of image
		 //		 + 0.5f*(dc*dc + dr*dr)
		 );
	      // << endl << flush;
	    }	    

    if (_amplitude.isValid())
      for (r = 0; r < h; r++)
	for (c = 0; c < w; c++) // here previously was c<h
	  {
	    conf.poke(c, r) *= amplConfMap(_amplitude.peek(c, r))*
	      coordConfMap(c, r, w, h);
	  }
    return(conf);
}

//
// this stuff is tuned for the perceptron data
// must fix and make generic
//
#define LOW_THRESH (150.0f)
#define HIGH_THRESH (2700.0f)
#define LOW_TRANS_WIDTH (1.0f)
#define HIGH_TRANS_WIDTH (100.0f)

float SphericalAmplitudeRangeScan::amplConfMap(float ampl) const
{

    float float_tmp;

    if ((ampl > (LOW_THRESH + LOW_TRANS_WIDTH/2.0f))&&
	(ampl < (HIGH_THRESH - HIGH_TRANS_WIDTH/2.0f)))
	{
	    return(1.0f);
	}

    if ((ampl < (LOW_THRESH - LOW_TRANS_WIDTH/2.0f))||
	(ampl > (HIGH_THRESH + HIGH_TRANS_WIDTH/2.0f)))
	{
	    return(0.0f);
	}
    
    if (ampl > (HIGH_THRESH - HIGH_TRANS_WIDTH/2.0f))
	{
	    float_tmp = (ampl - (HIGH_THRESH - HIGH_TRANS_WIDTH/2.0f))
		/HIGH_TRANS_WIDTH;
	    // else calculate the transition
	    return(0.5f*(cosFast(M_PI*float_tmp) + 1.0f));
	}
    else
	{
	    float_tmp = (ampl - (LOW_THRESH - LOW_TRANS_WIDTH/2.0f))
		/LOW_TRANS_WIDTH;
	    // else calculate the transition
	    return(0.5f*(cosFast(M_PI*float_tmp) + 1.0f));
	}
}

#define WINDOW_STD (0.25f)
float SphericalAmplitudeRangeScan::coordConfMap(unsigned u, unsigned v, 
						unsigned w, unsigned h) const
{
    float c_u = w/2.0f, c_v = h/2.0f;

    float x = (u - c_u)/(float)(w);
    float y = (v - c_v)/(float)(h);

    return(exp((-1.0f/WINDOW_STD)*(x*x + y*y)));
}
