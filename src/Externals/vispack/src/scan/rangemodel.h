/* *	sccsid "@(#) voxmodel.h     2.0 10/27/94" */

//
// rangemodel.h 
//
//

#ifndef __rangemodel_h
#define __rangemodel_h

#include "voxmodel/voxmodel.h"
#include "scan.h"

#include <limits.h>
// this is so that the volume is initially empty
#define INIT_VALUE (-FLT_MIN)

// template over the type of scan
template<class ScanType>
class RangeModel: public VoxModel
{
 protected:
  // must clean this up if you use array of ptrs
  VISArray<ScanType*> _scans;
  VISArray<VISMatrix> _transforms, _transforms_inv;
  float _force_window;
  int _num_scales;
  
  // sets relationship between vol units and image units
  // bigger means that object come out smaller
  // must establish relationship between volume coordinates
  // and world coordinates
  // Scaling and rigid transformation
  int _num_scans;

  float 
    _vol_units, 
    // sets the distance over which the range maps are active.
    _center_x, _center_y, _center_z,
    // the scale used to sample the depth map during initialization
    _init_scale;
    
  const ScanType & scan(unsigned i) const
    { return(*scan.peek(i));}

  const ScanType* scanRef(unsigned i)
    { return(scan.peek(i));}

  void setDefaults()
    {
      _num_scans = 0;
      _num_scales = 1;
      _vol_units = 1.0f;
      //
      _energy_weight = 0.0f;
      _normal_energy_weight = 0.0f;
      _grow_weight = 1.0f;          
      _const_grow_weight = 0.0f;
      _curve_weight = 0.0f;
      _gauss_curve_weight = _neg_curve_weight = _pos_curve_weight = 0.0f;
      _weighted_curve_weight = 0.0f;
      _special_curve_weight = 0.0f;
    }

  void initializeConstants();

  void myload(VPF::ParameterFile pf, const char* prefix = "");
  VISMatrix readTransformation(VPF::ParameterFile pf, const char* string) const;

 public:

  RangeModel(unsigned w, unsigned h, unsigned d)
    :VoxModel(VISVolume<float>(w, h, d))
    {
      setDefaults();
    }

  RangeModel()
    :VoxModel()
    {
      setDefaults();
    }

  ~RangeModel()
    {
      for (int i = 0; i < _num_scans; i++)
	delete _scans.poke(i);
    }

  RangeModel(VISVolume<float> vol)
    :VoxModel(vol)
    {
      setDefaults();
    }

  RangeModel(VPF::ParameterFile pf, const char* prefix = "")
    {
      setDefaults();
      loadParams(pf, prefix);
      initializeConstants();
    }

  int addScan(ScanType* scan, const VISMatrix &trans);

  void loadParams(VPF::ParameterFile pf, const char* prefix = "")
    {myload(pf, prefix);}
  
  void initializeScan();

  virtual float grow(float x, float y, float z, float n_x, 
		     float n_y, float n_z);

  virtual int numScans() const
    {return(_num_scans);}
  // wowaaa!  must fix this....

  void transformAllScans(const VISMatrix &T)
    {
      // which way should the inverse be applied...
      VISMatrix T_inv = T.inverseSVD();
      for (int i = 0; i < _num_scans; i++)
	{
	  (_transforms.poke(i)) = (_transforms.peek(i))*T;
	  (_transforms_inv.poke(i)) = T_inv*(_transforms_inv.peek(i));
	}
    }

  float g(VISVector point, int map_num, int scale) const;

  // static variables
  // put the ones first that you will need for calibration, etc
  typedef enum {NUMSCANS = 0, VOLWIDTH = 1, VOLHEIGHT = 2, VOLDEPTH = 3, 
		SCAN = 4, TRANSFORM = 5, VOLSCALE = 6, VOLUNITS = 7, 
		FORCEWINDOW = 8, NUMSCALES = 9, CURVEWEIGHT = 10} KeyWordNum;
  static char *keywords[11];
};

VISVolume<float> volumeResampleSpecial(const VISVolume<float> &vol_in, 
					float scale, 
					float x, float y, float z);

float partition(float);

#ifndef MANUAL_INSTANTIATION
#include "rangemodel.txx"
#endif

#endif
