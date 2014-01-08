#ifndef __NavyScan_h
#define __NavyScan_h

#include <scan/scan.h>
#define NAVY_SCAN_TYPE (3)
#define DEF_ALPHA   (0.177235)
#define DEF_BETA    (0.0589219)
#define DEF_D_THETA (-0.00017079)
#define DEF_THETA_0 (1.20059)
#define DEF_D_PHI   (0.000141459)
#define DEF_PHI_0   (-0.365692)
float navy_scan_ripple_adj[7] =
  {
    0.999643009,
    1.000642652,
    0.999643009,
    1.000142581,
    0.99944332,
    0.999842778,
    1.000642652
  };

//*********************************************************************
//*********************************************************************
//********************** 3D NAVY   Scan ****************************
//*********************************************************************
//*********************************************************************

class NavyScan: public Scan
{
protected:
  float _phi_0, _theta_0, _delta_phi, _delta_theta, _alpha, _beta;
  // this one is not virtual, and is faster for methods of this object.
  VISVector myImageCoord(const VISVector &p) const;
  int getIJ(float &i, float &j, float X, float Y, float Z) const;
  virtual void readParams(VPF::ParameterFile pf, const char* prefix = "");
  virtual VISImage<float> transformRange(const VISImage<float> &im) const;

public:

  virtual int type()
    { return(NAVY_SCAN_TYPE); }

  const NavyScan& operator=(const NavyScan &other)
    {
      Scan::operator=(other);
      if (this == &other)
	return(*this);
      _conf_map = other._conf_map;
      setphi0(other._phi_0); settheta0(other._phi_0); 
      setDeltaPhi(other._delta_phi); 
      setDeltaTheta(other._delta_theta); 
      return(*this);
    }

  NavyScan(const NavyScan &other)
    {
      operator=(other);
    }

  NavyScan(const VISImage<float> &range_image);
  NavyScan(const VISImage<float> &range_image, float normals_scale);


  NavyScan()
    {
      setDefaultParams();
    }

  virtual void setParams(const VISVector &params);
  virtual VISVector getParams() const;

  void loadParams(VPF::ParameterFile pf, const char* prefix = "")
    {
      Scan::loadParams(pf, prefix);
      readParams(pf, prefix);
    }

  NavyScan(VPF::ParameterFile pf, const char* prefix = "")
    {
      Scan::loadParams(pf, prefix);
      readParams(pf, prefix);
    }

  virtual VISImage<float> depth() const
  {
    return(_range_map_orig);
  }

  virtual VISImage<float> depth(int scale) const
    {
      return(_range_maps.peek(scale));
    }

  // see base class for documentation of methods.    

  virtual float depth(const VISVector &p) const;
  virtual float depth(const VISVector &p, int scale) const;
  virtual float distance(const VISVector &p) const;
  virtual VISVector grad_depth(const VISVector &p, int scale) const;
  virtual VISVector grad_conf(const VISVector &p, int scale) const;

// Gives change in a point position as a function of the parameters
  virtual VISMatrix dx_dparams(int c, int r) const;

  virtual float confidence(const VISVector &p) const;
  virtual float confidence(const VISVector &p, int scale) const;

  virtual VISVector lineOfSight(const VISVector &p)const;

  virtual VISVector imageCoord(const VISVector &p) const {return(myImageCoord(p));}
  virtual VISVector imageCoord(const VISVector &p, int scale) const;

  virtual VISVector get3DPoint(int i, int j) const;
  virtual VISVector get3DPoint(int i, int j, int scale) const;

  virtual VISVector getSurfaceNormal(int i, int j) const;
  virtual VISVector getSurfaceNormal(int i, int j, int scale) const;

  //  virtual VISVector rawToPoint(const VISVector &v);
  //  virtual VISVector rawToNormal(const VISVector &v);


  virtual void setDefaultParams()
    {
      setphi0(DEF_PHI_0);
      settheta0(DEF_THETA_0);
      setDeltaPhi(DEF_D_PHI);
      setDeltaTheta(DEF_D_THETA);
      setalpha(DEF_ALPHA);
      setbeta(DEF_BETA);
      setRangeOffset(0.0f);
      setRangeDelta(1.0f);
    }

  void setphi0(float phi_0) {_phi_0 = phi_0;}
  void settheta0(float theta_0) {_theta_0 = theta_0;}
  void setDeltaPhi(float delta_phi) {_delta_phi = delta_phi;}
  void setDeltaTheta(float delta_theta) {_delta_theta = delta_theta;}
  void setalpha(float a) {_alpha = a;}
  void setbeta(float b)  {_beta  = b;}

  // static variables
  // put the ones first that you will need for calibration, etc
#define NUMNAVYPARAMS (6) 
  // range = ROFFSET + RDELTA*range_raw
typedef enum 
    {
      PHI0 = 0,  THETA0 = 1,  DTHETA = 2,  DPHI = 3,  ALPHA = 4,   BETA = 5
    }  KeyWordNum;
  static char *keywords[NUMNAVYPARAMS];
};

#endif
