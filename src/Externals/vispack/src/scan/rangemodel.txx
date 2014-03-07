#ifndef	RANGEMODEL_C
#define RANGEMODEL_C

#include "vol/volumefile.h"
#include "image/imagefile.h"
#include <float.h>

#define TRYKEYWORDREQ(var, pref, key) { \
sprintf(scan_keyword, "%s%s", pref, keywords[key]); \
      if (VPF::set(var, pf[scan_keyword][0]) != VPF::VALID) { \
      cout << "Rangemodel::myload -- missing required keyword : " << scan_keyword << endl; \
      exit(-1); } \
}

#define TRYKEYWORD(var, pref, key) { \
      sprintf(scan_keyword, "%s%s", pref, keywords[key]); \
      VPF::set(var, pf[scan_keyword][0]); \
}

//
// initialize static member variables here.....
//
template<class ScanType>
char* RangeModel<ScanType>::keywords[11] = 
    {
      "NUM_SCANS", "VOL_WIDTH", "VOL_HEIGHT", "VOL_DEPTH", 
      "SCAN", "TRANSFORM", "VOL_SCALE", "VOL_UNITS", "FORCE_WINDOW", "NUM_SCALES", 
      "CURVE_WEIGHT"
    };

template<class ScanType>
void RangeModel<ScanType>::initializeConstants()
{
  if (RangeModel<ScanType>::_values.isValid())// **mdm**
    {
      _center_x = ((int)RangeModel<ScanType>::_values.width() - 1)/2.0f;// **mdm**
      _center_y = ((int)RangeModel<ScanType>::_values.height() -1)/2.0f;// **mdm**
      _center_z = ((int)RangeModel<ScanType>::_values.depth() - 1)/2.0f;// **mdm**
    }
  else
    {
      _center_x = _center_y = _center_z = 0.0f;
    }
}

template<class ScanType>
void RangeModel<ScanType>::initializeScan()
{
    initializeConstants();
    VISVector p(4);
    VISVector the_normal, p3;
    int l;
    RangeModel<ScanType>::_values.setBorder(-1.0f, 3);// **mdm**
    const ScanType* the_scan; 

// for debugging transformations
//     float tmp_x = 79, tmp_y = 29, tmp_z = 65;
//     p.poke(0) = (tmp_x - _center_x)*_vol_units;
//     p.poke(1) = (tmp_y - _center_y)*_vol_units;
    //     p.poke(2) = (tmp_z - _center_z)*_vol_units;
    //    VISVector imageCoord;

//    cout << "vol pt \n" << tmp_x << " " 
//	 << tmp_y << " " 
//	 << tmp_z << endl;

//    cout << "vol units " << _vol_units << 
//	"c_x " << _center_x << 
//	"c_y " << _center_y << 
//	"c_z " << _center_z << endl;

//    cout << "trans vol pt \n" << p << endl;
    

//     for (l = 0; l < _num_scans; l++)	
// 	{
// 	  the_scan = _scans.peek(l);
// //	    cout << "map " << l << endl;
// 	  T = _transforms.peek(l);
// //	    cout << "transform " << T << endl;
// 	    p_tmp = T*p;
// //	    cout << "scan coords: " << p_tmp << endl;
// 	    imageCoord = the_scan->imageCoord(p_tmp);
// //	    cout << "image coords: " << imageCoord.a() << " " <<
// //		 imageCoord.b() << endl;
// 	}

    boolean got_hit;
    p.poke(3) = 1.0f;

    cout << "trans " << _transforms.peek(0) << endl;

    //    cout << "vol units " << _vol_units << endl;

    int w = RangeModel<ScanType>::_values.width(), h = RangeModel<ScanType>::_values.height(), d = RangeModel<ScanType>::_values.depth();// **mdm**
    for (int i = 3; i < (d - 3); i++)
      for (int j = 3; j < (h - 3); j++)
	for (int k = 3; k < (w - 3); k++)
		{
		    p.poke(0) = ((float)k - _center_x)*_vol_units;
		    p.poke(1) = ((float)j - _center_y)*_vol_units;
		    p.poke(2) = ((float)i - _center_z)*_vol_units;
		    
		    got_hit = FALSE;
		    RangeModel<ScanType>::_values.poke(k, j, i) = -FLT_MIN;// **mdm**
		    
		    for (l = 0; l < _num_scans; l++)	
		      {
			the_scan = _scans.peek(l);
			if (the_scan->depth(_transforms.peek(l)*p) != the_scan->errorCode())
			  {
			    RangeModel<ScanType>::_values.poke(k, j, i) += g(p, l, 0);// **mdm**
			    got_hit = TRUE;
			  }
			//else 
			//cout << "got no hit" << _transforms.peek(l)*p << endl;
		      }
		    if (!got_hit)
		      RangeModel<ScanType>::_values.poke(k, j, i) = -1.0f;// **mdm**
		}

        RangeModel<ScanType>::_values.print();// **mdm**
        VISVolumeFile vol_file;
        vol_file.writeVTK(RangeModel<ScanType>::_values, "values_init.vtk");// **mdm**

	//    VISImageFile im_file;
	//    im_file.write((_values.image()).becomeFlat(), "values_init.fit");

    construct_lists();

}


template<class ScanType>
float RangeModel<ScanType>::grow(float x, float y, float z, 
			float n_x, float n_y, float n_z)
{
    VISVector p(4), p_tmp;
    VISVector line_of_sight, normal(4);
    float cosine;

    p.poke(0) = (x - _center_x)*_vol_units;
    p.poke(1) = (y - _center_y)*_vol_units;
    p.poke(2) = (z - _center_z)*_vol_units;
    p.poke(3) = 1.0f;

    normal.poke(0) = n_x;
    normal.poke(1) = n_y;
    normal.poke(2) = n_z;
    normal.poke(3) = 0.0f;

    float grow_mag = 0.0f, weighting;;
    const ScanType* the_scan;
    int l;

    for (l = 0; l < _num_scans; l++)	
      {
	the_scan = _scans.peek(l);
	p_tmp = _transforms.peek(l) * p;
	line_of_sight = _transforms_inv.peek(l)*(the_scan)->lineOfSight(p_tmp);
	//
	// this depends on how you define the normals and LOS
	// but if the normal are "inward facing" it should probably be negative
	//
	cosine = -normal.dot(line_of_sight);
	//
	//	    if (cosine < 0.0f)
	// should be negative the old style range maps
	//
	if (cosine >= 0.0f)
	  {
	    weighting = partition(acosFast(cosine));
	    grow_mag += weighting*g(p, l, 0);
	  }
      }
    return(grow_mag);
}

template<class ScanType>
VISMatrix RangeModel<ScanType>::readTransformation(VPF::ParameterFile pf, const char* string) 
const
{
  int r, c, j;
  VISMatrix transform(4, 4);
  try
    {
      VPF::set(r, pf[string][0]);
      VPF::set(c, pf[string][1]);
      if ((c == 4)&&(r == 4))
	for (j = 0; j < 16; j++)
	  {
	    VPF::set(
		     transform.poke(j/4, j%4),
		     pf[string][j+2]);
	  }
      else 
	{
	  cout << "RangeModel:: readTransformation---did not get valid transform" 
	       << endl;
	  transform = VISIdentity(4);
	}
    }
  catch (VPF::Exception e)
    {
      // igore transform
      cout << "RangeModel:: readTransformation---did not get valid transform" 
	   << endl;
      transform = VISIdentity(4);
    }
  return(transform);
}

template<class ScanType>
void RangeModel<ScanType>::myload(VPF::ParameterFile pf, const char* prefix)
{
  // must read depth image, number of scales, confidience
  char scan_keyword[80];
  char filename[80];
  int w, h, d;
  float tmp;
  VISImageFile im_file;
  char this_prefix[80], trans_string[80];
  VISMatrix transformation;  

   TRYKEYWORDREQ(w, prefix, VOLWIDTH);
   TRYKEYWORDREQ(h, prefix, VOLHEIGHT);
   TRYKEYWORDREQ(d, prefix, VOLDEPTH);
      //
      RangeModel<ScanType>::_values = VISVolume<float>(w, h, d);// **mdm**
      //
   TRYKEYWORDREQ(_force_window, prefix, FORCEWINDOW);
   TRYKEYWORDREQ(_num_scans, prefix, NUMSCANS);
   //
  cout << "num_scans " << _num_scans << endl;
  //

  for (int i = 0; i < _num_scans; i++)
    {
      sprintf(trans_string, "%s%s_%d_%s", prefix, keywords[SCAN],
	      i, keywords[TRANSFORM]);
      sprintf(this_prefix, "%s%s_%d_", prefix, keywords[SCAN], i);
      addScan(new ScanType(pf, this_prefix), readTransformation(pf, trans_string));
      cout << "trans string " << trans_string  << readTransformation(pf, trans_string) << endl;
    }

  //
  tmp = -1.0f;
  //
   TRYKEYWORD(tmp, prefix, VOLSCALE);
   //
   if (tmp > 0.0f)
     _vol_units = tmp/w;
   //
   TRYKEYWORD(_vol_units, prefix, VOLUNITS);
   //
   cout << _vol_units << endl;
   //
   TRYKEYWORD(_num_scales, prefix, NUMSCALES);
   TRYKEYWORD(_curve_weight, prefix, CURVEWEIGHT);
      
  try
    {
      sprintf(scan_keyword, "%s%s", 
	      prefix, keywords[TRANSFORM]);
      // should this be inverted, or not!!!
      transformAllScans(readTransformation(pf, scan_keyword));
    }
  catch(VPF::Exception e)
    {
      // not required
    }


}

template <class ScanType>
float RangeModel<ScanType>::g(VISVector point, int map_num, int scale) const
{
  VISVector this_point;
  //    cout << point << flush;
  this_point = _transforms.peek(map_num)*point;
  float r = (_scans.peek(map_num))->depth(this_point, scale);
  float d = (_scans.peek(map_num))->distance(this_point);
  float c = (_scans.peek(map_num))->confidence(this_point, scale);
  //   cout << "point " << point;
  //   cout << "this_point " << this_point;
  //     cout << "r " << r << " d " << d << " c " << c << endl << flush;
     //     cout << "g " << c*(r - d) << endl << flush;
  //     cout << "g " << c*(r - d)*gaussFast((r - d)/_force_window) 
  //          << endl << flush;

  // this is for window cutoff
  //   d = d - r;
  //   if (fabs(d) > _force_window)
  //     return(0.0);
  //   else
  //     return(c*d);

  // this is for gaussian window
   float epsilon; 
   d = d - r; 
   if (d >= 0.0f)  
     epsilon = 0.0f; 
   else 
     epsilon = -FLT_MIN; 
   if (c > 0.0f) 
     return(c*d*gaussFast(d/_force_window) + epsilon);      
   else 
     return(0.0f);      
}

template <class ScanType>
int RangeModel<ScanType>::addScan(ScanType* scan, const VISMatrix &trans)
{
    _scans.appendItem(scan);
    _transforms.appendItem(trans);
    _transforms_inv.appendItem(trans.inverseSVD());
//    if (_num_scales > 1)
//    cout << "about to create scales " << _num_scales << endl << flush;
    scan->createScales(_num_scales);
//    cout << "created scales " << scan->numScales() << endl << flush;
    return(_scans.n());
}

#endif
