//
// modelnode.h 
//


#ifndef iris_octmodel_h
#define iris_octmodel_h

#define UNKNOWN_STATUS 0x0  // 0
#define ACTIVE_STATUS 0x01    // 1
#define INSIDE_STATUS 0x02    // 2 
#define OUTSIDE_STATUS 0x04   // 4 
#define CHANGING_STATUS 0x08   // 8
#define INSIDE_STATUS_2 0x016   // 22
#define OUTSIDE_STATUS_2 0x032  // 50

#define BORDER_WIDTH 3

#define CURVE_DT (1.0f/6.0f)
// this should normally be 0.5
// #define WAVE_DT (0.5f)
#define WAVE_DT (1.0f/3.0f)

#include "oct/octmesh.h"
#include "vol/volindexlist.h"
#include "vol/volume.h"
#include "mat/matrix.h"
#include "util/geometry.h"
#include <string.h>


class VolHood;

//template <class T>
//class VISArray<T>;

float sphere(float, float, float);

class ModelNode
{
    float _value;
//    unsigned _x, _y, _z;
    byte _state;
//    float _change;

  public:

    ModelNode() {;}
    ModelNode(const ModelNode& other) 
    {
	_value = other._value;
//	_x = other._x;
//	_y = other._y;
//	_z = other._z;
	_state = other._state;
//	_change = other._change;
    }

    ModelNode(float v) 
    {
	_value = v;
	_state = UNKNOWN_STATUS;
    }

    ModelNode(float v, byte state) 
    {
	_value = v;
	_state = state;
    }

    boolean operator==(const ModelNode& other) const
    {
	return((_value == other._value)
	       &&(_state == other._state));
    }

    const ModelNode& operator+=(const ModelNode& other)
    {
	_value += other._value;
	_state = UNKNOWN_STATUS;
	return(*this);
    }

    ModelNode operator/(int d)
    {
	ModelNode r;
	r._value = _value/d;
	_state = UNKNOWN_STATUS;
	return(r);
    }

    boolean operator!=(const ModelNode& other) const
    {
	return((_value != other._value)
	       ||(_state != other._state));
    }
    
//    x(unsigned xx) {_x = x;}
//    y(unsigned yy) {_y = y;}
//    z(unsigned zz) {_z = z;}
//
//    unsigned x() {return(_x);}
//    unsigned y() {return(_y);}
//    unsigned z() {return(_z);}

    void value(float vv) {_value = vv;} 
    float value() const {return(_value);} 

    void state(byte ss) {_state = ss;}
    byte state() const {return(_state);}
    
//    change(float cc) {_change = cc;}
//    float change() {return(_change);}
};


typedef OctMeshLeaf<ModelNode> OctLeafNode; 

class OctNode: public OctMeshNode<ModelNode>
{
  public:
    
    const OctNode& assignVol(const VISVolume<float>& vol, float padding);

    OctNode(const ModelNode& v, unsigned levels):
	OctMeshNode<ModelNode>(v, levels)
    {;}

    OctNode():OctMeshNode<ModelNode>()
    {;}

    const OctNode& operator=(const VISVolume<float>& vol)
    { 
	return(assignVol(vol, 0.0f)); 
    }

    VolIndexValueVISList* makeActiveVISList();

    VISVolume<float> values() const;
    VISVolume<byte> state() const;
    
    // returns a 3x3 neighborhood of the float values 
    VolHood Neighborhood(unsigned, unsigned, unsigned);
    // returns a 3x3 neighborhood of the float values 
    // only 6-connected neighbors
    VolHood Neighborhood6(unsigned, unsigned, unsigned);
    // returns a 3x3 neighborhood of the float values 
    // only 18-connected neighbors
    VolHood Neighborhood18(unsigned, unsigned, unsigned);
    // returns a 3x3 neighborhood of the float values 
    // only 26-connected neighbors
    VolHood Neighborhood26(unsigned, unsigned, unsigned);

    void setBorder(ModelNode value, unsigned w);
};

#define DIFFERENCE_FACTOR (2.0f/5.0f)
#define DIFFERENCE_FACTOR_2 (2.0f/5.0f)
#define CHANGE_FACTOR (1.0f/5.0f)

class OctModel
{
  protected:
    OctNode* _values;
    VISVolIndexVISList *_status_up_list, 
    *_status_down_list;
    VolIndexValueVISList *_inside_list, *_outside_list;
    VolIndexValueVISList *_active_list;
    float _force_weight, 
	_normal_force_weight, 
	_curve_weight,
	_grow_weight, 
	_const_grow_weight, 
	_weighted_curve_weight, 
	_gauss_curve_weight,
	_pos_curve_weight, 
	_neg_curve_weight;
		   
    float _wave_dt, _curve_dt;

    void setup();
    
    
  public:
    OctModel(const VISVolume<float>& vol);
    OctModel(unsigned levels);
    OctModel() {setup();}

//    OctModel* resample(float scale);
//    OctModel* resample(unsigned w, unsigned h, unsigned d);
    
    void construct_lists(float scale_factor);
    void makeLayers(float scale_factor);
    virtual void construct_lists()  {construct_lists(DIFFERENCE_FACTOR);}
    virtual float calculate_change();
    virtual float calculate_change_curve();
    virtual float calculate_change_all();
    void iterate();
    float update(float dt);
    void updateLayers();
    void load();
    
    void rescale(float scale);

    const OctNode* values() const {return(_values);} 
    OctNode* valuesRef() const {return(_values);} 
    const VolIndexValueVISList* list() const {return(_active_list);}
    VolIndexValueVISList* listRef() {return(_active_list);}
    void values(OctNode* new_values) 
    {
	if (_values)
	    delete _values;
	_values = new_values;
    }

// these function control the movement of the surface
    virtual Point3 force(float x, float y, float z)
    { return(Point3(0.0, 0.0, 0.0)); }

    virtual Point3 normal_force(float x, float y, float z, 
				float n_x, float n_y, float n_z)
    { return(Point3(0.0, 0.0, 0.0)); }

    virtual float grow(float x, float y, float z, 
		       float n_x, float n_y, float n_z)
    { return(0.0); }

    virtual float grow(float x, float y, float z)
    { return(0.0); }


    void setParams(float force_weight, 
		   float normal_force_weight, 
		   float curve_weight,
		   float grow_weight, 
		   float const_grow_weight, 
		   float weighted_curve_weight, 
		   float gauss_curve_weight,
		   float pos_curve_weight, 
		   float neg_curve_weight
		   )
    {
	_force_weight = force_weight;
	_normal_force_weight = normal_force_weight;
	_curve_weight = curve_weight;
	_grow_weight = grow_weight;
	_const_grow_weight = const_grow_weight;
	_weighted_curve_weight = weighted_curve_weight;
	_gauss_curve_weight = gauss_curve_weight;
	_pos_curve_weight = pos_curve_weight;
	_neg_curve_weight = neg_curve_weight;
    }

    void loadParams(char* fname)
    {
	float energy_weight, 
	    normal_energy_weight, 
	    curve_weight, 
	    grow_weight, 
	    const_grow_weight, 
	    weighted_curve_weight, 
	    gauss_curve_weight,
	    pos_curve_weight, 
	    neg_curve_weight;

	FILE* config_file;
	char error_string[80];

	if ((config_file = fopen(fname,"r")) == NULL)
	    {
		strcpy(error_string, 
		       "voxmodel::read_params could not open file: ");
		strcat(error_string, fname);
		ERROR(error_string);
		setParams(-1.0, 0.0f, 0.0f, 0.0f, 0.0f,
			  0.0f, 0.0f, 0.0f, 0.0f);
	    }
	else
	    {
		fscanf(config_file, "%f", &energy_weight);
		while (fgetc(config_file) != '\n');

		fscanf(config_file, "%f", &normal_energy_weight);
		while (fgetc(config_file) != '\n');

		fscanf(config_file, "%f", &curve_weight);
		while (fgetc(config_file) != '\n');

		fscanf(config_file, "%f", &grow_weight);
		while (fgetc(config_file) != '\n');

		fscanf(config_file, "%f", &const_grow_weight);
		while (fgetc(config_file) != '\n');

		fscanf(config_file, "%f", &weighted_curve_weight);
		while (fgetc(config_file) != '\n');

		fscanf(config_file, "%f", &gauss_curve_weight);
		while (fgetc(config_file) != '\n');

		fscanf(config_file, "%f", &neg_curve_weight);
		while (fgetc(config_file) != '\n');

		fscanf(config_file, "%f", &pos_curve_weight);
		while (fgetc(config_file) != '\n');

		setParams(energy_weight, 
			  normal_energy_weight, 
			  curve_weight, 
			  grow_weight, 
			  const_grow_weight, 
			  weighted_curve_weight, 
			  gauss_curve_weight,
			  pos_curve_weight, 
			  neg_curve_weight);

		fclose(config_file);
	    }
	    

    }

};

#define D_X 0
#define D_Y 1
#define D_Z 2
#define D_XX 3
#define D_XY 4
#define D_YY 5
#define D_XZ 6
#define D_YZ 7
#define D_ZZ 8

#define D_XF 0
#define D_XB 1
#define D_YF 2
#define D_YB 3
#define D_ZF 4
#define D_ZB 5

// this is a 3x3 volume that defines a neighborhood
class VolHood: public VISVolume<float>
{

  public:
    
    VolHood():VISVolume<float>(3,3,3)
    {;}

    VolHood(float v):VISVolume<float>(3,3,3)
    {
	operator=(v);
    }

    VolHood(const VolHood& other):VISVolume<float>(other)
    {;}

    const VolHood& operator=(const VolHood& other)
    {
	VISVolume<float>::operator=(other);
	return(*this); 
    }

    const VolHood& operator=(const VISVolume<float>& other)
    {
	if ((other.width() == 3)
	    &&(other.height() == 3)
	    &&(other.depth() == 3))
	    VISVolume<float>::operator=(other);
	else
	    other.getROI(0, 0, 0, 3, 3, 3);
	return(*this); 
    }

    const VolHood& operator=(float v)
    {
	VISVolume<float>::operator=(v);
	return(*this); 
    }

    VISArray<float> derivatives();
    VISArray<float> derivsHalf();

};


class OctMultiModel: public OctModel
{
  protected:
    float _x_offset, _y_offset, _z_offset;
    float _scale;
//    int _levels;
    float (*_field_function)(float, float, float);

  public:
    OctMultiModel():OctModel()  
    {_field_function = NULL;}

    OctMultiModel(const VISVolume<float> vol):OctModel(vol)  
    {
	_field_function = NULL;
    }

    OctMultiModel(unsigned levels):OctModel()
    {
	_field_function = NULL;
	setLevels(levels);
    }

    OctMultiModel(float (*field_function)(float, float, float)):OctModel()  
    {
	_field_function = field_function;
    }

  //  this is the function around which the level-set will be constructed

    virtual float fieldFunction(float x, float y, float z)
    {
	if (_field_function)
	    return(_field_function(x, y, z));
	else return(0.0f);
    }

    float fieldFunction(Point3 p)
    {
//	if (x.size() >= 3)
	    return(fieldFunction(p.x(), p.y(), p.z()));
//	else 
//	    return(0.0f);
    }

    Point3 position(unsigned, unsigned, unsigned);
    Point3 position(float, float, float);

    void higherResInterp();
    void higherResFunc();
    void lowerRes();

    void offset_x(float x) {_x_offset = x;}
    void offset_y(float y) {_y_offset = y;}
    void offset_z(float z) {_z_offset = z;}
    void scale(float s) {_scale = s;}

    void initialize(unsigned levels);
    void initialize();
    void setLevels(unsigned levels);

    VolIndexValueVISList* makeActiveVISList(VISVolIndexVISList* list);
};



#endif
