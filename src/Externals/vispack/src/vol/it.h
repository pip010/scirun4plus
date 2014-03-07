/*
****************************************************************
*  it.C - image test
*
* @(#)it.h 1.1 94/04/26 11:00:35
*
****************************************************************
*/

#ifndef IT_H
#define IT_H

// -------------------------------------------------------------------------

#include "util/classtypes.h"
#include "ui/application.h"
#include "ui/timer.h"
#include "image/imageRGBA.h"
#include "image/imagefile.h"
#include "image/imgmodels/matrix.h"
#include "util/array.h"
#include "util/string.h"
#include "localimage.h"
#include "volume.h"
#include <limits.h>

class GfxVolTransformInterp;


// -------------------------------------------------------------------------
// some constants needed to keep track of the the status of pixels.

#define ACTIVE_STATUS 0x01    // 1
#define INSIDE_STATUS 0x02    // 2 
#define OUTSIDE_STATUS 0x04   // 4 
#define CHANGING_STATUS 0x08   // 8
#define INSIDE_STATUS_2 0x016   // 22
#define OUTSIDE_STATUS_2 0x032  // 50

// -------------------------------------------------------------------------



class GfxFrameGrabber;
class GfxMenu;
class GfxImageViewer;
class GfxSquareFinder;
class GfxImIndexList;
class GfxVolIndexList;





class GfxRangeMap
{
  protected:
    LocalImage<float> _depth;
    GfxTransform _transform;
    GfxTransform _transform_inv;

  public:
    const GfxTransform& transform() const {return(_transform);}
    const GfxTransform& transform_inv() const {return(_transform_inv);}
    const GfxImage<float>& depth() const {return(_depth);}
    GfxRangeMap() {;}

    GfxRangeMap(const GfxTransform& t, 
		const GfxImage<float>& d) 
    {_transform = t; _depth = d; _transform_inv = inverse(t);}
    GfxRangeMap(const GfxRangeMap& rm):_transform(rm.transform()),
    _depth(rm.depth()),_transform_inv(rm.transform_inv()) {}
    float depth(const Point3& p)  const;
    
};


class GfxScannerMap: public GfxRangeMap
{
  protected:
    GfxImage<float> _confidence;
    float _x_scale, _r0, _c0, _y0, _fy, _z0;
    float _error_flag;
    GfxImage<float> _depth_dx, _depth_dy;

  public:

    GfxScannerMap() {;}
    GfxScannerMap(const GfxTransform& t, 
		  const GfxImage<float>& d, 
		  const GfxImage<float>& confidence, 
		  float r0, float c0, float y0, float z0, 
		  float x_scale, float fy, float derivative_scale);

    GfxScannerMap(const GfxScannerMap& other):GfxRangeMap(other.transform(), 
							  other.GfxRangeMap::
							  depth())
    {
	_r0 = other._r0; 
	_c0 = other._c0; 
	_x_scale = other._x_scale; 
	_y0 = other._y0; 
	_z0 = other._z0; 
	_fy = other._fy;
	_error_flag = other._error_flag;
	_confidence = other._confidence;
	_depth_dx = other._depth_dx;
	_depth_dy = other._depth_dy;
    }


    float error_flag() const {return(_error_flag);}
    float depth(const Point3& p) const;
    float depth(const Point3& p, GfxImage<float> kernel) const;
    Point3 d_depth(const Point3& p) const;
    Point3 normal(const Point3& p) const;
    float confidence(const Point3&) const;
};


class GfxScan: public Array<GfxScannerMap>
{
  public:
    GfxScan(char* filename);
    GfxScan() {;}
// sets relationship between vol units and image units
// bigger means that object come out smaller
    float _vol_scale, 
// sets the distance over which the range maps are active.
	_force_window, 
// these determine the position of the center of the volume relative to
// its size
	_center_x_ratio, _center_y_ratio, _center_z_ratio, 
	_init_scale;
    float volScale() {return(_vol_scale);}
};



// -------------------------------------------------------------------------

class IT : public GfxApplication
{
  GfxImageViewer*	_viewer;

  GfxMenu* build_menubar();		// Create menubar

  GfxTool* _select_tool;
  GfxTool* _rotate_tool;

    float _center_x, _center_y, _center_z;

// these determine the relationship between a unit in the volume and
// a unit in 3D space
    float _vol_units;

    GfxVolume<float> _values;
    GfxVolume<float> _distance;
    GfxVolume<float> _image, _image_dx, _image_dy, _image_dz;
    GfxVolume<float> _image_2, _image_2_dx, _image_2_dy, _image_2_dz;
    GfxScan _scan;

    GfxVolume<float> _distance_1;
    GfxVolume<float> _distance_1_dx, _distance_1_dy, _distance_1_dz;

    GfxVolume<float> _distance_2;
    GfxVolume<float> _distance_2_dx, _distance_2_dy, _distance_2_dz;

    GfxVolume<unsigned> _next_state;
    GfxVolume<int> _active_pixels;
    GfxVolume<float> _change;
    Array< GfxVolume<float> > *_derivatives, *_derivatives_first_order;
    
    int _active;
    int _iterations;
    int _i;
    float _energy_weight;
    float _energy_volume_weight;
    float _curve_weight;
    float _weighted_curve_weight;
    float _grow_weight;
    float _distance_grow_weight;
    float _texture_weight;
    float _positive_curve_weight;
    float _negative_curve_weight;

    int _display_iterations;
    int _save_iterations;
    
    GfxVolIndexList *_active_list, *_status_up_list, *_status_down_list,
    *_inside_list, *_outside_list;

    Array<GfxString> _filenames;
    
  public:
    IT();
    virtual ~IT();

    void init();
    void init_images();
    void object_panel();
    void start();
    void stop();
    void moments();
    void quit();
    void save();
    void load();
    void iterate();
    void fit_model();
    void calculate();
    void calculate_list();
    void construct_lists(float scale_factor = 1.0f);
    void read_range_data();
    void init_range_data();
    void init_range_data2();
    void init_range_data_union();
    float calculate_list3(float);
    float calculate_list3(float,
			  const GfxTransform&, 
			  float& max_change);
    float calculate_list(const GfxMatrix& affine_trans, 
			const GfxPoint& trans, float alpha);
    void calculate_list(const GfxMatrix& affine_trans, 
			const GfxPoint& trans);
    float calculate_list_scan(float max_dt);
    void volume_test();
    void volume_test2();
    void volume_test3();
    void double_size();
    void read_params(char* fname);
    void get_and_read_params();
    void range_sphere();
    void read_and_write_ply();
    void positive_smooth();

    GfxTimer* _timer;

    void init_morph();
    void save_morph();
    void init_crossing_rects();
    void do_morph();

    void surface_texture();
    GfxImage<float> _image_texture[2];

    void alpha_blend();

    GfxVolTransformInterp* _vol_trans;
    GfxTransform _transform;

};



// -------------------------------------------------------------------------

GfxVolume<int> compareDistance(const GfxVolume<float>& d1, 
			      const GfxVolume<float>& d2,
			      const GfxMatrix& a, 
			      const GfxPoint& b, 
			      float alpha);

GfxVolume<float> signedDistance(const GfxVolume<float>& vol_in);
GfxVolume<float> signedDistanceSq(const GfxVolume<float>& vol_in);


#define DIFFERENCE_FACTOR (2.0f/5.0f)
#define DIFFERENCE_FACTOR_2 (2.0f/5.0f)
#define CHANGE_FACTOR (1.0f/5.0f)


#endif
