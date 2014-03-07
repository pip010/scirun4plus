#ifndef voxmodel
#define voxmodel

#include "util/mathutil.h"
#include "mat/matrix.h"
#include "vol/volume.h"  
#include "vol/volindexlist.h"
#include "voxmodel/activeindex.h"
#include "param/param.h"

#define ACTIVE_STATUS 10
#define CHANGING_STATUS 1  // new active pixel 
// the rest of status flags are computed with the formula 
// ACTIVE_STATUS+10*level+1 for inside layers
// ACTIVE_STATUS+10*level-1 for outside layers 
#define INSIDE_STATUS_1 11
#define OUTSIDE_STATUS_1 9
#define INSIDE_STATUS_2 21
#define OUTSIDE_STATUS_2 19
#define INSIDE_STATUS_3 31
#define OUTSIDE_STATUS_3 29
#define INSIDE_STATUS_4 41
#define OUTSIDE_STATUS_4 39

#define CHANGE_FACTOR (0.5f)     // layers are separated at odd multiples of this number
//#define CURVE_DT (0.165f)
//#define WAVE_DT  (0.165f)        // 1/4 for 2D, 1/6 for 3D
#define CURVE_DT (0.12f)
#define WAVE_DT  (0.12f)        // 1/4 for 2D, 1/6 for 3D
#define MIN_NORM ((float)1.0e-6) // added to norm. of vectors to avoid divide by 0 

#define NUM_DERIVS (21)
#define DX (0)
#define DY (1)
#define DZ (2)
#define DPX (3)
#define DPY (4)
#define DPZ (5)
#define DMX (6)
#define DMY (7)
#define DMZ (8)
#define DYPX (9)
#define DZPX (10)
#define DXPY (11)
#define DZPY (12)
#define DXPZ (13)
#define DYPZ (14)
#define DYMX (15)
#define DZMX (16)
#define DXMY (17)
#define DZMY (18)
#define DXMZ (19)
#define DYMZ (20)
#define NSIZE (27)

// VoxModelData is the class to be used with the templated class ActiveIndex
class VoxModelData {
 public: 
  VoxModelData() {
    _value=_value_prev=0.0f;
    _Nphi[0]=_Nphi[1]=_Nphi[2]=0.0f;
    _active_pixel=0;
  }
  VoxModelData (float a, float b, byte c) {
    _value=a;
    _value_prev=b;
    _active_pixel=c;
    _Nphi[0]=_Nphi[1]=_Nphi[2]=0.0f;
  }
  VoxModelData (float a, byte b) {
    _value=_value_prev=a;
    _active_pixel=b;
    _Nphi[0]=_Nphi[1]=_Nphi[2]=0.0f;
  }
  
  float  _value, _value_prev, _Nphi[3];
  byte   _active_pixel;
  byte   _flag;

  void operator=(const VoxModelData &right) {
    _value=right._value;
    _value_prev=right._value_prev;
    _active_pixel=right._active_pixel;
    _Nphi[0]=right._Nphi[0];
    _Nphi[1]=right._Nphi[1];
    _Nphi[2]=right._Nphi[2];
  }
  void operator=(const VoxModelData *right) {
    _value=right->_value;
    _value_prev=right->_value_prev;
    _active_pixel=right->_active_pixel;
    _Nphi[0]=right->_Nphi[0];
    _Nphi[1]=right->_Nphi[1];
    _Nphi[2]=right->_Nphi[2];
  }
  static char *name () {
    return "VoxModelData";
  }
}; 

class LevelSetModel {
 protected: 
  // this is the volume that holds the indices into _data
  VISVolume <unsigned>       _index; 
  // the values and status flags are stored in here (see VoxModelData above)
  ActiveIndex <VoxModelData> *_data;
  // lists of pixels that are changing status because they have 
  // passed either above or below their active range
  VISVolIndexVISList *_status_up_list, *_status_down_list;
  // lists of pixels in the inside and outside layers
  // these lists should not be of type VolIndexValueVISList (inefficient memory) 
  VISVolIndexVISList *_inside_list, *_outside_list;
  // the list of active pixels. this data structure also has a 
  // float field to store changes to the active pixels
  VolIndexValueVISList *_active_list;
  int _width, _height, _depth; // size of volume
  unsigned _number_of_layers;  // the moving band is 2*_number_of_layers+1 pixels wide
  unsigned _border_width;      // active set can't move into border
  unsigned _INSIDE, _OUTSIDE;  // fixed indices (into _data) for INSIDE and OUTSIDE voxels
  float _INSIDE_VALUE;         // = _number_of_layers + 0.5
  float _wave_dt, _curve_dt;
  void dump                 (int, int, int, int = 1); // print local neighborhood
  // fetch a neighborhood of values from _data
  inline int  get_phi_neighborhood (int, int, int, float*);  
  
  // functions used by the constructors
  void allocate_lists  (); 
  void initialize_data ();
  void load_lists (VISVolume <float>&, VISVolume <byte>&);

  float distanceIterate (int num_extra_levels, 
			 VISVolIndexVISList *buff1, VISVolIndexVISList *buff2,
			 VISVolume<float> &val_exact, VISVolume<float> &val_prev,
			 VISVolume<boolean> &active_pixels);

 public:
  LevelSetModel (VISVolume <float>&, unsigned); 
  LevelSetModel (VISVolume <float>*, unsigned);
  LevelSetModel (unsigned);
  LevelSetModel ();            
  LevelSetModel (char *, int *); // load from previous save_internal_status
  ~LevelSetModel();

  void load_constants (VPF::ParameterFile);

  int get_number_of_layers () {
    return _number_of_layers;
  }

  int width  () {return _width;}
  int height () {return _height;}
  int depth  () {return _depth;}

  const VolIndexValueVISList* list() const {
    return(_active_list);
  }
  
  VolIndexValueVISList* listRef() {
    return(_active_list);
  }
  
  const VISVolIndexVISList* in_out_listRef (int m, int l) const {
    if (m>=_number_of_layers) return (VISVolIndexVISList*)NULL;
    if (l) return(&(_inside_list[m])); else return (&(_outside_list[m]));
  }

  VISVolIndexVISList* in_out_listRef (int m, int l) {
    if (m>=_number_of_layers) return (VISVolIndexVISList*)NULL;
    if (l) return(&(_inside_list[m])); else return (&(_outside_list[m]));
  }
  
  // saves internal status for use with constructor
  void save_internal_status (char*, int, boolean = FALSE);
  // checks the validity of internal state 
  boolean check_internal_status(boolean = FALSE);   
  
  VISVolume <float> values();  // constructs and returns volume of VoxModelData->_value
  VISVolume <byte>  status();  // constructs and returns volume of VoxModelData->_active_pixel
  VISVolume <float> Nphi(int); // constructs and returns volume of VoxModelData->_Nphi[i]
  VISVolume <float> calculateDistanceTransform (int); // exact distance transform
  VISVolume <float> calculateDistanceTransform (int, VISVolume<float> &); 
  VISVolume <float> calculateDistanceTransform (VISVolume<float> &); 
  
  // took out the following as they do not readily apply to 
  // the new dynamic memory storage scheme 
  // VISVolume<float>& valuesRef() {return(_values);}
  // void values(VISVolume<float>& new_values) 
  //  {_values = new_values;}

  virtual float constrain (float value, int x, int y, int z) { return(value); }
  
  // these functions replace construct_lists and construct_lists_9layers 
  // of old VoxModel. the number of layers is now a variable (_number_of_layers)
  void construct_lists (VISVolume<float>&);
  void construct_lists (VISVolume<float>*);
  
  float         update           (float dt);               
  virtual float calculate_change () = 0;      
  float         iterate          () {
    return(update(calculate_change()));} // the main processing function

  typedef enum {BORDERWIDTH = 0} KeyWordNum;
  static char *keywords[1];
};

class VoxModel: public LevelSetModel  {
 protected:
  // these are the weights for the various forcesthat control the model
  // The update is of the form
  // d \phi/ d t = A \nabla \phi \cdot F(x) + B \nabla \phi \cdot F(x, n)
  //  + C | \nabla \phi | G(x, n) + D | \nabla \phi | + E |\nabla \phi| \kappa
  //  + F | \nabla \phi | curvfit(\kappa)
  // 
  // where x is the position in 3D, n is the surface normal, and \kappa 
  // is the mean curvature of the level set.
  //
  float _energy_weight;         // A
  float _normal_energy_weight;  // B
  float _grow_weight;           // C
  float _const_grow_weight;     // D 
  float _curve_weight;          // E 
  float _curvature_grow_weight; // F
  // sub classes can construct arbitrary functions of the shape matrix
  float _special_curve_weight;
  // experimental weights 
  float _pos_curve_weight, _neg_curve_weight;
  float _weighted_curve_weight;
  float _gauss_curve_weight;

  inline float     calculate_curvature (const float *n);
  inline float     calculate_curvature (int x, int y, int z);
  inline VISMatrix calculate_curvature (const float* derivs, 
					float &curve_trace, float &curve_norm);
  
 public:
  VoxModel  (VISVolume<float> &vol_float, unsigned n) : LevelSetModel(vol_float,n) {
    set_default_parameters();
  }
  VoxModel  (VISVolume<float> *vol_float, unsigned n) : LevelSetModel(vol_float,n) {
    set_default_parameters();
  }
  VoxModel  (unsigned n) : LevelSetModel(n) {
    set_default_parameters();
  }
  VoxModel  () {
    set_default_parameters();
  }
  VoxModel (char *fname, int *k) : LevelSetModel (fname,k) {
    set_default_parameters();
  }
  ~VoxModel () {};
  
  void set_default_parameters ();
  void compute_normals_from_phi_on_active_set ();
  virtual float calculate_change();
  
  virtual VISVector     force (float x, float y, float z) {
    return(VISVector(0.0f, 0.0f, 0.0f));}

  virtual VISVector     normal_force (float x, float y, float z, 
				      float n_x, float n_y, float n_z) { 
    return(VISVector (0.0f, 0.0f, 0.0f)); }

  virtual float grow (float x, float y, float z, 
		      float n_x, float n_y, float n_z) {
    return(0.0f);}

  virtual float grow (float x, float y, float z) {
    return(0.0f);}

  virtual float curvature_grow (float, int x, int y, int z) {
    return(0.0f);}
  
  virtual float curvature_grow (float, float x, float y, float z) {
    return(0.0f);}

  virtual float specialCurve (const VISMatrix &curve, const VISVector &p, 
			      const VISVector &n) {
    return(0.0f); }

  float iterate () {
    compute_normals_from_phi_on_active_set ();
    return(update(calculate_change()));
  }
};

// functions to fetch a neighborhoods from volumes
void get_neighborhood (int, int, int, float*, const VISVolume<float>&);
void get_neighborhood (int, int, int, float*, VISVolume<float>*);

// compute only forward&backward differences
void         get_derivs_1sided    (const float*, float*); 
// compute all differences
void         get_derivs           (const float*, float*);

float VISMaxAbs            (float a, float b) {
  if (fabs(a) > fabs(b)) return(a); else  return(b);}

float distance_symmetric (VISVolume<float> &A, VISVolume<float> &B, float *d1, float *d2);

#endif
