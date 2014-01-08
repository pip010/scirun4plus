
#ifndef vis_scan_h
#define vis_scan_h
// -------------------------------------------------------------------------

//#include "util/classtypes.h"
//#include "ui/application.h"
//#include "ui/timer.h"

#include "util/array.h"
//#include "util/string.h"
//#include "util/geometry.h"
#include "util/mathutil.h"
#include "image/indexlist.h"
#include <limits.h>
#include "scanimage.h"
#include "mat/matrix.h"
#include <float.h>
#include "param/param.h"

// this is for the morphing
//class VISVolTransformInterp;

// -------------------------------------------------------------------------

// ****************************************************
// ****************************************************
// ******************** Scan **************************
// ****************************************************
// ****************************************************

#define SCAN_TYPE (0)
#define SPHERICAL_SCAN_TYPE (1)
#define ORTHO_SCAN_TYPE (2)

#define BORDER (3)

class Scan
{
 protected:
  // multi-resolutions
  VISArray< VISImage<float> > _range_maps;
  VISArray< VISImage<float> > _range_maps_dx;
  VISArray< VISImage<float> > _range_maps_dy;
  // full resolutions
  VISImage<float>  _range_map_orig;
  VISImage<float>  _range_map_orig_dx;
  VISImage<float>  _range_map_orig_dy;
  // conf stuff
  VISImage<float> _conf_map;
  VISImage<float> _conf_map_dx;
  VISImage<float> _conf_map_dy;
  VISArray< VISImage<float> > _conf_maps;
  VISArray< VISImage<float> > _conf_maps_dx;
  VISArray< VISImage<float> > _conf_maps_dy;
  //
  int _num_scales;
  float _normals_scale;
  // factor of pixel distance as you go up in scale
  float *_scale_factor;
  // used to calculate surface normals
  //    VISImage<float>  _dx_kernel, _dy_kernel;

  float _r_delta, _r_offset;

  void readParams(VPF::ParameterFile pf, const char* prefix = "");
  void initRangeImages(const VISImage<float> &range_map, float normals_scale);

  // this gives subclasses a chance to directly modify the range map
  // to compensate for problems in the measuring device.
  virtual VISImage<float>  transformRange(const VISImage<float> &im) const { return im; }

 public:
  virtual int type()
    { return(SCAN_TYPE); }

  const Scan& operator=(const Scan &other)
  {
    if (this == &other)
      return(*this);
    _range_maps = other._range_maps;
    _range_maps_dx = other._range_maps_dx;
    _range_maps_dy = other._range_maps_dy;
    _range_map_orig = other._range_map_orig;
    _range_map_orig_dx = other._range_map_orig_dx;
    _range_map_orig_dy = other._range_map_orig_dy;
    _num_scales = other._num_scales;

    _conf_map = other._conf_map;
    _conf_map_dx = other._conf_map_dx;
    _conf_map_dy = other._conf_map_dy;
    _conf_maps = other._conf_maps;
    _conf_maps_dx = other._conf_maps_dx;
    _conf_maps_dy = other._conf_maps_dy;

      setRangeDelta(other._r_offset);  
      setRangeOffset(other._r_delta); 

    //	cout << "num scales now is " << _num_scales << endl << flush;
    _scale_factor = new float[_num_scales];
    for (int i = 0; i < _num_scales; i++)
      cout << (_scale_factor[i] = other._scale_factor[i]) << " ";
    return(*this);
  }

  Scan(const Scan &other)
    {
      operator=(other);
    }

  Scan(const VISImage<float> &range_map);
  // you can smooth while taking the derivatives.
  Scan(const VISImage<float> &range_map, float normals_scale);
  Scan(VPF::ParameterFile pf, const char* prefix = "")
    {
      setDefaultParams();
      //      cout << "normals scale " << _normals_scale << endl;
      readParams(pf, prefix);
    }

  void loadParams(VPF::ParameterFile pf, const char* prefix = "")
    {      
      setDefaultParams();
      //      cout << "normals scale " << _normals_scale << endl;
      readParams(pf, prefix);
    }

  void setRangeOffset(float r_offset) {_r_offset = r_offset;}
  void setRangeDelta(float r_delta) {_r_delta = r_delta;} 

  Scan()
    {
      _range_map_orig_dx = _range_map_orig_dy 
	= _range_map_orig = _conf_map = _conf_map_dx 
	= _conf_map_dy = VISImage<float>();
      _scale_factor = NULL;
      setDefaultParams();
      //	_dx_kernel = dx_kernel(1);
      //	_dy_kernel = dy_kernel(1);
    }

  ~Scan()
    {
      if (_scale_factor != NULL)
	delete[] _scale_factor;
    }

  // you can tell it how man sub resolutions
  virtual void createScales(int num_scales);
  virtual void createScalesConf(int num_scales);
  // this sets a weighting factor for each range measurment
  virtual void setConfidence(const VISImage<float> &conf_map)
    {
      //      cout << "got set confidence" << endl;
      //      conf_map.print();
      if (conf_map.isValid())
	{
	  _conf_map = conf_map.setBorder(0.0f, 1);
	  _conf_map_dx = -1.0f*conf_map.dx();
	  _conf_map_dy = -1.0f*conf_map.dy();
	}
      else
	{
	  _conf_map = VISImage<float>();
	  _conf_map_dx = VISImage<float>();
	  _conf_map_dy = VISImage<float>();
	}
      createScalesConf(_num_scales);
    }

  // the value of the range map along the line of sight
  // through the point p.
  // if it's not in the viewing frustrum of scan, return "errorCode()"
  virtual float depth(const VISVector &p) const
    {return(errorCode());}

  // the value of the range map along the line of sight
  // through the point p, at a reduced scale (=1 -> original).
  //
  virtual float depth(const VISVector &p, int scale) const
    {return(errorCode());}

  //
  //  The change in the depth per unit change of the vector p
  // 
  //  only inportant for model fitting
  //
  virtual VISVector grad_depth(const VISVector &p, int scale) const
    {
      return VISVector();
    }

  virtual VISVector grad_conf(const VISVector &p, int scale) const
    {
      return VISVector();
    }

  // Gives change in a point position as a function of the parameters
  // 
  //  only inportant for model fitting
  //
  //virtual VISMatrix dx_dparams(int c, int r) const
  //    { return VISMatrix(); }
  //  virtual VISMatrix dx_dparams(int c, int r, int scale) const
  //    { return VISMatrix(); }

  // the "inverse of the variance" of the range measurements
  // taken along the line of sight through the point p
  virtual float confidence(const VISVector &p) const
    {return(errorCode());}

  // the "inverse of the variance" of the range measurements
  // taken along the line of sight through the point p, at 
  // a reduced scale (=1 -> original).
  virtual float confidence(const VISVector &p, int scale) const
    {return(errorCode());}

  //
  // line of sight passint through p --- in the form of direction cosines
  //
  virtual VISVector lineOfSight(const VISVector &p) const
    {
      return VISVector();
    }

  virtual VISVector lineOfSight(int i, int j) const
    {
      return VISVector();
    }

  //
  // a generic way of setting internal parameters
  //
  virtual void setParams(const VISVector &params) {;}
  virtual VISVector getParams() const {return VISVector();}

  // gives you the coordinates in the image associated with this point
  // projects through this point
  virtual VISVector imageCoord(const VISVector &p) const
    {
      return VISVector();
    }	

  virtual VISVector imageCoord(const VISVector &p, int scale) const
    {
      return VISVector();
    }

  virtual VISImage<float> depth() const
    {return(_range_map_orig);}

  virtual void depth(const VISImage<float> &d)
    { initRangeImages(d, _normals_scale);  }

  virtual VISImage<float> depth(int scale) const
    {
      //	cout << "yeah, we got def bitz" << endl;
      return(_range_maps.peek(scale));
    }

  // distance to the point along the line of sight    
  virtual float distance(const VISVector &p) const
    {
      return(-1.0f*(lineOfSight(p).dot(p)));
    }

  float errorCode() const {return(-FLT_MAX);}

  virtual VISVector get3DPoint(int i, int j) const
    {return VISVector();}
  virtual VISVector get3DPoint(int i, int j, int scale) const
    {return VISVector();}
  // 
  //  Only important for model fitting
  // 
  virtual VISVector getSurfaceNormal(int i, int j) const
    {return VISVector();}
  virtual VISVector getSurfaceNormal(int i, int j, int scale) const
    {return VISVector();}

  float scaleFactor(int scale_level) const
    {
      if ((scale_level >= 0)&&(scale_level < _num_scales))
	{
	  return(_scale_factor[scale_level]);
	}
    }

  int numScales() const {return(_num_scales);}

  VISArray< VISImage<float> > get3DPoints(const VISMatrix &transform) const;
  VISArray< VISImage<float> > get3DPoints(const VISMatrix &transform, int scale) const;

  void setDefaultParams() 
    {
      //      cout << "got scan set def params " << endl;
      _num_scales = 1; 
      _normals_scale = 0.0f;
      setConfidence(VISImage<float>()); 
    }

  // static variables
  // put the ones first that you will need for calibration, etc
  typedef enum {NUMSCALES = 0, CONFIDENCEIMAGE = 1, RANGEIMAGE = 2, 
		NORMALSSCALE = 3, ROFFSET = 4, RDELTA = 5} KeyWordNum;
  static char* keywords[6];
};



// ****************************************************
// ****************************************************
// ***************** TwoDOrthoScan  *******************
// ****************************************************
// ****************************************************


class TwoDOrthoScan: public Scan
{
  public:

    TwoDOrthoScan(const VISImage<float> &range_image):Scan(range_image)
    {
    }

    TwoDOrthoScan(const VISImage<float> &range_image, float scale)
	:Scan(range_image, scale)
    {
      ;
    }

    TwoDOrthoScan(VPF::ParameterFile pf, const char* prefix = ""):Scan(pf, prefix)
      {
	readParams(pf, prefix);
      }

// see base class for documentation of methods.    
    virtual float depth(const VISVector &p) const;
    virtual float depth(const VISVector &p, int scale) const;
    virtual float confidence(const VISVector &p) const 
    { return(1.0); }
    virtual float confidence(const VISVector &p, int scale) const 
    { return(1.0); }

    virtual VISVector grad_depth(const VISVector &p, int scale) const;

    virtual VISVector lineOfSight(const VISVector &p) const;
    virtual VISVector lineOfSight(int i, int j) const;

    virtual VISVector get3DPoint(int i, int j) const;
    virtual VISVector get3DPoint(int i, int j, int scale) const;

    virtual VISVector getSurfaceNormal(int i, int j) const;
    virtual VISVector getSurfaceNormal(int i, int j, int scale) const;

};


//*********************************************************************
//*********************************************************************
//********************** 3D Spherical Scan ****************************
//*********************************************************************
//*********************************************************************

#define DEFAULT_FOV ((45.0f)*(M_PI/180.0))

class SphericalScan: public Scan
{
protected:
  float _r0, _c0, _delta_phi, _delta_theta;
  // added these to do calibration

  // this one is not virtual, and is faster for methods of this object.
  VISVector myImageCoord(const VISVector &p) const;

  virtual void readParams(VPF::ParameterFile pf, const char* prefix = "");

public:

  virtual int type()
    { return(SPHERICAL_SCAN_TYPE); }

  const SphericalScan& operator=(const SphericalScan &other)
    {
      Scan::operator=(other);
      if (this == &other)
	return(*this);
      setR0(other._r0); setC0(other._c0); setDeltaPhi(other._delta_phi); 
      return(*this);
    }

  SphericalScan(const SphericalScan &other)
    {
      operator=(other);
    }

  SphericalScan(const VISImage<float> &range_image);
  SphericalScan(const VISImage<float> &range_image, float normals_scale);

  SphericalScan(VPF::ParameterFile pf, const char* prefix = ""):Scan(pf, prefix)
    {
      setDefaultParams();
      //      cout << "_normals scales " << _normals_scale << endl;
      readParams(pf, prefix);
    }

  void loadParams(VPF::ParameterFile pf, const char* prefix = "")
    {
      Scan::loadParams(pf, prefix);
      readParams(pf, prefix);
    }

  SphericalScan()
    {
      setDefaultParams();
      //      setR0(0.0f); setC0(0.0f); setDeltaPhi(0.0f); setDeltaTheta(0.0f); 
      //      setRangeDelta(1.0f);  setRangeOffset(0.0f); 
      //      setConfidence(VISImage<float>());
    }

  virtual void setParams(const VISVector &params);
  virtual VISVector getParams() const;

    void setR0(float r0) {_r0 = r0;}
    void setC0(float c0) {_c0 = c0;}
    void setDeltaPhi(float delta_phi) {_delta_phi = delta_phi;}
    void setDeltaTheta(float delta_theta) {_delta_theta = delta_theta;}


  virtual VISImage<float> depth() const
  {
    return(_r_offset + _r_delta*_range_map_orig);
  }

  virtual VISImage<float> depthRaw() const
  {
    return(_range_map_orig);
  }

  virtual VISImage<float> depth(int scale) const
    {
//	cout << "yeah, we got def bitz" << endl;
      return(_r_offset + _r_delta*_range_maps.peek(scale));
    }


  // see base class for documentation of methods.    

  virtual float depth(const VISVector &p) const;
  virtual float depth(const VISVector &p, int scale) const;
  virtual VISVector grad_depth(const VISVector &p, int scale) const;
  virtual VISVector grad_conf(const VISVector &p, int scale) const;

// Gives change in a point position as a function of the parameters
  //  virtual VISMatrix dx_dparams(int c, int r) const;

  virtual float confidence(const VISVector &p) const;
  virtual float confidence(const VISVector &p, int scale) const;

  virtual VISVector lineOfSight(const VISVector &p)const;
  virtual VISVector lineOfSight(int i, int j)const;

  virtual VISVector imageCoord(const VISVector &p) const
    {return(myImageCoord(p));}
  virtual VISVector imageCoord(const VISVector &p, int scale) const;

  virtual VISVector get3DPoint(int i, int j) const;
  virtual VISVector get3DPoint(int i, int j, int scale) const;

  virtual VISVector getSurfaceNormal(int i, int j) const;
  virtual VISVector getSurfaceNormal(int i, int j, int scale) const;
  //  virtual VISVector rawToPoint(const VISVector &v);
  //  virtual VISVector rawToNormal(const VISVector &v);
  void setDefaultParams();

  // static variables
typedef enum {R0 = 0, C0 = 1, DTHETA = 2, DPHI = 3} KeyWordNum;
  // put the ones first that you will need for calibration, etc
#define NUMSPHEREPARAMS (4) 
  static char *keywords[NUMSPHEREPARAMS];
  // range = ROFFSET + RDELTA*range_raw

};

class OrthoScan: public SphericalScan
{
protected:
  float _delta_x, _delta_y;

  // this one is not virtual, and is faster for methods of this object.
  VISVector myImageCoord(const VISVector &p) const;

  virtual void readParams(VPF::ParameterFile F, const char* prefix = "")
    {
      // no special params for this object
   }


public:

  virtual int type()
    { return(ORTHO_SCAN_TYPE); }

  const OrthoScan& operator=(const OrthoScan &other)
    {
      Scan::operator=(other);
      if (this == &other)
	return(*this);
      _conf_map = other._conf_map;
      setParams(other._c0,
		other._r0,
		other._delta_x,
		other._delta_y);
      return(*this);
    }

  OrthoScan(const OrthoScan &other)
    {
      operator=(other);
    }

  //  OrthoScan(VPF::ParameterFile F);

  OrthoScan(const VISImage<float> &range_image)
    :SphericalScan(range_image) 
    {
      setDefaultParams();
    }

  OrthoScan(const VISImage<float> &range_image, float normals_scale)
    :SphericalScan(range_image, normals_scale) 
    {
      setDefaultParams();
    }

  OrthoScan()
    {
      setParams(VISVector(0.0f, 0.0f, 0.0f, 0.0f));
    }

  virtual void setParams(const VISVector &params)
    {
      if (params.n() >= 4)
	{
	  _r0 = params.peek(0); _c0 =  params.peek(1);
	  _delta_x = params.peek(2); _delta_y = params.peek(3);
	}
      else
	cout << "WARN: Setparams --- bad vector size" << endl;
    }

  void setParams(float r0, float c0, float delta_x, 
		 float delta_y)
    {
	  _r0 = r0; _c0 =  c0;
	  _delta_x = delta_x; _delta_y = delta_y;
    }

  OrthoScan(VPF::ParameterFile pf, const char* prefix = ""):SphericalScan(pf, prefix)
    {
      setDefaultParams();
      readParams(pf, prefix);
      // piggy back off of sphere scan
      _delta_x = _delta_theta;
      _delta_y = _delta_phi;
    }

  void loadParams(VPF::ParameterFile pf, const char* prefix = "")
    {
      SphericalScan::loadParams(pf, prefix);
      readParams(pf, prefix);
      // piggy back off of sphere scan
      _delta_x = _delta_theta;
      _delta_y = _delta_phi;
    }

  void getParams(float &r0, float &c0, float &delta_x, 
		 float &delta_y) const
    {
      r0 = _r0; c0 = _c0; 
      delta_x = _delta_x; delta_y = _delta_y; 
    }

  //     void getParams(float &r0, float &c0, float &delta_x, 
  // 		   float &delta_y, int scale) const
  // 	{
  // 	  r0 = _r0/_scale_factor[scale]; c0 = _c0/_scale_factor[scale]; 
  // 	  delta_x = _delta_x*_scale_factor[scale]; 
  // 	  delta_y = _delta_y*_scale_factor[scale]; 
  // 	}

  // see base class for documentation of methods.    

  virtual VISVector grad_depth(const VISVector &p, int scale) const;
  virtual VISVector lineOfSight(const VISVector &p) const
    { return(VISVector(0.0f, 0.0f, -1.0f, 0.0f)); }
  virtual VISVector lineOfSight(int i, int j) const
    { return(VISVector(0.0f, 0.0f, -1.0f, 0.0f)); }
  virtual VISVector imageCoord(const VISVector &p) const
    {return(myImageCoord(p));}

  virtual VISVector get3DPoint(int i, int j) const;
  virtual VISVector get3DPoint(int i, int j, int scale) const;

  virtual VISVector getSurfaceNormal(int i, int j) const;
  virtual VISVector getSurfaceNormal(int i, int j, int scale) const;

  void setDefaultParams()
    {
      cout << "setDefaultParams in ortho scan" << endl;
      // nothing new to add
    }

  // static variables
  //  typedef enum {DX = 0, DY = 1} KeyWordNum;
  // put the ones first that you will need for calibration, etc
  //  static char *keywords[2];
  // range = ROFFSET + RDELTA*range_raw

};

class SphericalAmplitudeRangeScan: public SphericalScan
{
 protected:
  VISImage<float> _amplitude;
  float _noise_level, _spatial_uncertainty;
  
  // read amplitude and create confidence here
// also read noise_level and spatial_uncertainty
  void readParams(VPF::ParameterFile F, const char* prefix = "");

 public:
  virtual void loadParams(VPF::ParameterFile F, const char* prefix = "")
    {
      readParams(F, prefix);
      SphericalScan::loadParams(F, prefix);
    }

  virtual VISImage<float> createConfidence() const;

  SphericalAmplitudeRangeScan()
    :SphericalScan()
    {
      setDefaultParams();
    }

  SphericalAmplitudeRangeScan(VPF::ParameterFile pf, const char* prefix = "")
    :SphericalScan(pf, prefix)
    {
      setDefaultParams();
      readParams(pf, prefix);
    }

  void setDefaultParams()
    {
      _noise_level = 1.0f; _spatial_uncertainty = 0.0f;
      _amplitude = VISImage<float>();
    }

  virtual float amplConfMap(float ampl) const;
  virtual float coordConfMap(unsigned u, unsigned v, unsigned w, unsigned h) const;

  // static variables
  typedef enum {SPATIALUNCERTAINTY = 0, NOISE_LEVEL = 1} KeyWordNum;

  // put the ones first that you will need for calibration, etc
  static char *keywords[2];  // range = ROFFSET + RDELTA*range_raw

  const SphericalAmplitudeRangeScan& operator=(const SphericalAmplitudeRangeScan &other)
    {
      SphericalScan::operator=(other);
      if (this == &other)
	return(*this);
      _amplitude = other._amplitude;
      _noise_level = other._noise_level;
      _spatial_uncertainty = other._spatial_uncertainty;
      return(*this);
    }

  SphericalAmplitudeRangeScan(const SphericalAmplitudeRangeScan &other)
    {
      operator=(other);
    }
};

VISImage<float> subSampleAverage(VISImage<float> im);
VISImage<float> subSampleMin(VISImage<float> im);

#endif



