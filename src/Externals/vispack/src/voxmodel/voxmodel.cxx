#include <string.h>
#include <float.h>
#include <math.h> 
#include <stdio.h> 
   
#include "mat/matrix.h"
#include "vol/volume.h"
#include "vol/volumefile.h"
#include "vol/volindexlist.h"
#include "util/mathutil.h" 
#include "voxmodel/activeindex.h"
#include "param/param.h"
#include "voxmodel/voxmodel.h"

#define TRYKEYWORDREQ(var, pref, key) { \
sprintf(scan_keyword, "%s%s", pref, keywords[key]); \
      if (VPF::set(var, pf[scan_keyword][0]) != VPF::VALID) { \
      cout << "NormModel3D -- missing required keyword : " << scan_keyword << endl; \
      exit(-1); } \
cout << "VoxModel keyword is " << scan_keyword; \
cout << "      --- variable is " << var << endl; \
}

#define TRYKEYWORD(var, pref, key) { \
      sprintf(scan_keyword, "%s%s", pref, keywords[key]); \
      VPF::set(var, pf[scan_keyword][0]); \
cout << "VoxModel keyword is " << scan_keyword; \
cout << "     --- variable is " << var << endl; \
}

char* LevelSetModel::keywords[1] = {
  "BORDER_WIDTH"
};

// if the input volume is an actual signed distance transform 
// (not just an inside/outside approx. to the distance transform)
// then use the next statement
// #define TRUE_DISTANCE_INPUT

// define the number and positions of 4-connected neighbors in 3D
#define NUM_N 6         
int Neighbors[6][3]= {{0,0,1},{0,0,-1},{0,1,0},{0,-1,0},{1,0,0},{-1,0,0}}; 

// positions of normals needed for curvature calculation at phi(x,y,z)
#define NUMC 7
int NeighC[NUMC][3]= {{0,0,1},{1,0,0},{1,0,1},{0,1,0},{0,1,1},{1,1,0},{1,1,1}}; 

void LevelSetModel::dump (int x, int y, int z, int m) {
  int i, j, k;
  VoxModelData *thisvoxel;
  
  cout<<"values ----------\n";
  for (k=-m;k<=m;k++) { for (j=-m;j<=m;j++) { for (i=-m;i<=m;i++) {
    thisvoxel=_data->pop(_index(x+i,y+j,z+k));
    cout<<thisvoxel->_value<<" ";
  } cout<<endl; } cout<<endl; }
  cout<<"active ----------\n";
  for (k=-m;k<=m;k++) { for (j=-m;j<=m;j++) { for (i=-m;i<=m;i++) {
    thisvoxel=_data->pop(_index(x+i,y+j,z+k));
    cout<<(int)(thisvoxel->_active_pixel)<<" ";
  } cout<<endl; } cout<<endl; }
  cout<<"=================\n";
}

LevelSetModel::LevelSetModel(VISVolume<float> &vol, unsigned n) { 
  // this is a full constructor
  _height = vol.height(); 
  _width  = vol.width(); 
  _depth  = vol.depth();
  _number_of_layers = n;  
  _border_width=n+1;
  _INSIDE_VALUE=(float)(CHANGE_FACTOR*(2*n+1));
  _index=VISVolume <unsigned> (vol.width(),vol.height(),vol.depth());
  allocate_lists();
  initialize_data();
  construct_lists(vol);
}

LevelSetModel::LevelSetModel(VISVolume<float> *vol, unsigned n) { 
  // this is a full constructor
  _height = vol->height(); 
  _width = vol->width(); 
  _depth = vol->depth();
  _number_of_layers=n;  
  _border_width=n+1;
  _INSIDE_VALUE=(float)(CHANGE_FACTOR*(2*n+1));
  _index=VISVolume <unsigned> (_width,_height,_depth);
  allocate_lists();
  initialize_data();
  construct_lists(vol);
}

LevelSetModel::LevelSetModel (unsigned n) {
  // this is a partial constructor -- the following need to be done later:
  // _index has to be allocated after defining the dimensions of the volume
  // _data has to be initialized and construct_lists need to be called
  _number_of_layers=n;  
  _border_width=n+1;
  _INSIDE_VALUE=(float)(CHANGE_FACTOR*(2*n+1));
  allocate_lists();
}

LevelSetModel::LevelSetModel () {   
  // default constructor ... all initializations have to be done later  
}

LevelSetModel::LevelSetModel (char *fname, int *cnt) {
  // this constructor is used for loading internal status saved from
  // another instance of LevelSetModel by calling save_internal_status
  char fname2[100], dummy[100];
  int n;
  VISVolumeFile volfile;
  cout<<"Loading LevelSetModel\n";
  sprintf(fname2,"%s.value",fname);VISVolume <float> vol=volfile.read_float(fname2); 
  sprintf(fname2,"%s.active",fname);VISVolume <byte> act=volfile.read(fname2); 
  sprintf(fname2,"%s.save",fname);
  FILE *loadfile=fopen(fname2,"r");
  fscanf(loadfile,"%s",&dummy);
  cout<<dummy;
  fscanf(loadfile,"%d %d",&n,cnt);
  fclose(loadfile);
  cout<<"Number of layers: "<<n<<endl;
  _height = vol.height(); 
  _width  = vol.width(); 
  _depth  = vol.depth();
  _number_of_layers = n;  
  _border_width=n+1;
  _INSIDE_VALUE=(float)(CHANGE_FACTOR*(2*n+1));
  _index=VISVolume <unsigned> (vol.width(),vol.height(),vol.depth());
  allocate_lists();
  initialize_data();
  load_lists(vol,act);
  cout<<"Done loading LevelSetModel\n";
}

void LevelSetModel::load_constants (VPF::ParameterFile pf) {
  char scan_keyword[80];
  TRYKEYWORD(_border_width,"", BORDERWIDTH);
}

void LevelSetModel::load_lists (VISVolume <float>& vol, VISVolume <byte>& act) {
  unsigned x, y, z;
  int i, k, l, level;
  VoxModelData thisvoxel;

  cout<<"Initializing lists...\n";
  thisvoxel._value=thisvoxel._value_prev=_INSIDE_VALUE;
  thisvoxel._active_pixel=0;
  _INSIDE=_data->push(thisvoxel);  // initialize fixed inside index
  thisvoxel._value=thisvoxel._value_prev=-_INSIDE_VALUE;
  _OUTSIDE=_data->push(thisvoxel); // initialize fixed outside index
  _active_list->clean();
  for (i=0; i<_number_of_layers; i++) {
    _inside_list[i].clean();
    _outside_list[i].clean();
  }
  for (x = 0; x < _width; x++) for (y = 0; y < _height; y++) for (z = 0; z < _depth; z++) 
    if (act(x,y,z)) {
      thisvoxel._active_pixel=act(x,y,z);
      thisvoxel._value=thisvoxel._value_prev=vol(x,y,z);
      _index.poke(x,y,z)=_data->push(thisvoxel);
      k=act(x,y,z)-ACTIVE_STATUS;
      if (k!=0) {
	level=(int)rint(0.1*(float)k); 
	l=10*level;
	if (k<l) _outside_list[level].appendItem(VISVolIndex(x,y,z));
	else _inside_list[level].appendItem(VISVolIndex(x,y,z));
      } else _active_list->appendItem(VolIndexValue(x, y, z,0.0f));;
    }
    else {
      if (vol(x,y,z)>0.0f) _index.poke(x,y,z)=_INSIDE;
      else  _index.poke(x,y,z)=_OUTSIDE;
    }
  //cout<<check_internal_status(1)<<endl;exit(-1);
}

void LevelSetModel::allocate_lists () {
  // allocates memory for lists -- needs number of layers to be defined
  if (_number_of_layers<1) {
    cout<<"number of layers = "<<_number_of_layers<<" is not valid.\n";
    exit(-1);
  }
  _active_list      = new VolIndexValueVISList();
  _status_up_list   = new VISVolIndexVISList[2*_number_of_layers+1];
  _status_down_list = new VISVolIndexVISList[2*_number_of_layers+1];
  _inside_list      = new VISVolIndexVISList[_number_of_layers];
  _outside_list     = new VISVolIndexVISList[_number_of_layers];
}

void LevelSetModel::initialize_data () {
  // initializes _data variable -- needs dimensions of the volume
  // and the number of layers to be defined
  int temp=(2*_number_of_layers+1)*(int)rint(power((_width+_height+_depth)/3.0,2));
    _data= new ActiveIndex <VoxModelData> (2*temp,temp);
}

LevelSetModel::~LevelSetModel() {
  delete _data;
  delete _active_list;
  delete [] _status_up_list;
  delete [] _status_down_list;
  delete [] _inside_list;
  delete [] _outside_list;
}

void LevelSetModel::save_internal_status (char *fname, int k, boolean flag) {
  char fname2[100];

  VISVolumeFile volfile;
  cout<<"Saving internal status.\n";
  sprintf(fname2,"%s.active",fname);volfile.write_byte(status(),fname2); 
  if (!check_internal_status(flag)) {
    cout<<"Invalid internal status. Aborting save.\n";
    exit(-1);
  }
  sprintf(fname2,"%s.value",fname);volfile.write_float(values(),fname2); 
  sprintf(fname2,"%s.save",fname);
  FILE *savefile=fopen(fname2,"w");
  fprintf(savefile,"LevelSetModel_save_file\n%d %d\n",_number_of_layers,k);
  fclose(savefile);
  if (flag) exit(-1);
}

boolean LevelSetModel::check_internal_status (boolean flag) {
  // recall that:
  // ACTIVE_STATUS+10*level+1 for inside layers
  // ACTIVE_STATUS+10*level-1 for outside layers 
  int k, l, level;
  float lolim, hilim, val;
  VoxModelData *thisvoxel;
  unsigned x, y, z;
  cout<<"Testing VoxModelData...\n";
  for (x = 0; x < _width; x++) for (y = 0; y < _height; y++) for (z = 0; z < _depth; z++){
    thisvoxel=_data->pop(_index(x,y,z));
    if (thisvoxel->_active_pixel) {
      k=(thisvoxel->_active_pixel)-ACTIVE_STATUS;
      val=thisvoxel->_value;
      if (k!=0) {
	level=(int)rint(0.1*(float)k); 
	l=10*level;
	level++;
	if (k<l) level=-level;
      } else level=0;
      if (level>0) {
	hilim=CHANGE_FACTOR*((float(2*level+1)));
	lolim=hilim-2.0*CHANGE_FACTOR;
	if ((val>hilim)||(val<=lolim)) {
	  cout<<"LevelSetModel: Invalid internal status!\n";
	  cout<<level<<" "<<val<<endl;
	  return FALSE;
	}
      } else if (level<0) {
	lolim=CHANGE_FACTOR*((float(2*level-1)));
	hilim=lolim+2.0*CHANGE_FACTOR;
	if ((val<lolim)||(val>=hilim)) {
	  cout<<"LevelSetModel: Invalid internal status!\n";
	  cout<<level<<" "<<val<<endl;
	  return FALSE;
	}
      } else {
	if (fabs(val)>CHANGE_FACTOR) {
	  cout<<"LevelSetModel: Invalid internal status!\n";
	  cout<<level<<" "<<val<<endl;
	  return FALSE;
	}
      }
    } 
  }
  cout<<"VoxModelData test: pass\n";
  if (flag) {
    // this is a through test that will actually destroy this instance of VoxModelData
    // call it only in a call to save. Note: VoxModelData._value is not affected by this call
    cout<<"Testing lists...\n";
    VISVolIndex index;
    VolIndexValue v_index;
    VISVolIndexVISList *list;
    int i, j;
    _active_list->reset();
    while (_active_list->valid()) {
      v_index = _active_list->itemAtCurrent();
      x=v_index.a();y=v_index.b();z=v_index.c();
      thisvoxel=_data->pop(_index(x,y,z));
      if (thisvoxel->_active_pixel==ACTIVE_STATUS) thisvoxel->_active_pixel=0;
      else {
	cout<<"LevelSetModel: VoxModelData->_active_pixel ("<<thisvoxel->_active_pixel<<
	  ") does not match _active_list voxel\n";
	return FALSE;
      }
      _active_list->stepForward();
    }
    for (j=-1;j<=1;j+=2) for (i=0;i<_number_of_layers;i++) {
      if (j>0) list=&(_inside_list[i]); else list=&(_outside_list[i]);
      list->reset();
      while (list->valid()) {
	index=list->itemAtCurrent();
	x=index.a();y=index.b();z=index.c();
	thisvoxel=_data->pop(_index(x,y,z));
	level=ACTIVE_STATUS+i*10+j;
	if ((thisvoxel->_active_pixel)==level) thisvoxel->_active_pixel=0;
	else {
	  cout<<"LevelSetModel: VoxModelData->_active_pixel ("<<thisvoxel->_active_pixel<<
	    ") does not match list "<<j<<" "<<i<<endl;
	  return FALSE;
	}
	list->stepForward();
      }
    }
    for (x = 0; x < _width; x++) for (y = 0; y < _height; y++) for (z = 0; z < _depth; z++) {
      thisvoxel=_data->pop(_index(x,y,z));
      if (thisvoxel->_active_pixel) {
	cout<<"LevelSetModel: voxel ("<<thisvoxel->_active_pixel<<") was not found in any lists\n";
	return FALSE;
      }
    }
  }
  cout<<"LevelSetModel test: pass\n";
  return TRUE;
}

VISVolume<float> LevelSetModel::values() {
  // this function creates a volume of values from _data
  // normally used for output purposes
  VISVolume<float> ret(_width,_height,_depth);
  unsigned x, y, z;
  VoxModelData *thisvoxel;
  for (x = 0; x < _width; x++) for (y = 0; y < _height; y++) for (z = 0; z < _depth; z++){
    thisvoxel=_data->pop(_index(x,y,z));
    ret.poke(x,y,z)=thisvoxel->_value;
  }
  return ret;
}
 
VISVolume<byte> LevelSetModel::status() {
  // this function creates a volume of status flags from _data
  // normally used for output (debugging) purposes
  VISVolume<byte> ret(_width,_height,_depth);
  unsigned x, y, z;
  VoxModelData *thisvoxel;
  for (x = 0; x < _width; x++) for (y = 0; y < _height; y++) for (z = 0; z < _depth; z++) {
    thisvoxel=_data->pop(_index(x,y,z));
    ret.poke(x,y,z)=thisvoxel->_active_pixel;
  }
  return ret;
}

VISVolume<float> LevelSetModel::Nphi(int i) {
  // this function creates a volume of _Nphi[i] from _data
  // normally used for output purposes
  VISVolume<float> ret(_width,_height,_depth);
  unsigned x, y, z;
  VoxModelData *thisvoxel;
  for (x = 0; x < _width; x++) for (y = 0; y < _height; y++) for (z = 0; z < _depth; z++){
    thisvoxel=_data->pop(_index(x,y,z));
    ret.poke(x,y,z)=thisvoxel->_Nphi[i];
  }
  return ret;
}

inline boolean check_layer (int m, int l, int status) {
  // this function is used only by LevelSetModel::update ()
  if (!status) return FALSE; 
  // Inside 1 and outside 1 should always use active set (including new actives) for updates
  if (m==0) return ((status==ACTIVE_STATUS)||(status==CHANGING_STATUS));
  else {
    if (m==1) // 2nd layer looks at layer 1 + actives
      return ((status==ACTIVE_STATUS+2*l-1)||(status==CHANGING_STATUS)||(status==ACTIVE_STATUS));
    else {
      // other layers look at all previous layers 
      int diff=status-10*(int)floor(0.1*(float)status);
      if (l) return ((diff==1)&&(status<(ACTIVE_STATUS+10*m+1)));
      else  return ((diff==9)&&(status<(ACTIVE_STATUS+10*m-1)));
    }
  }
}

inline boolean check_layer2 (int m, int l, int status) {
  return (status==(ACTIVE_STATUS+10*m+2*l-1));
}

void LevelSetModel::construct_lists (VISVolume<float> &values) {
  // constructs the lists -- can only be called after full initialization
  unsigned x, y, z, bbox[2][3];
  int k, i, j, dist_sign;
  float len, neighborhood[NSIZE], derivs[NUM_DERIVS], dist;
  VISVector d_F(2);
  VISVolIndex index;    
  VolIndexValue v_index;    
  VISVolIndexVISList *current[2], *nextlayer[2];
  VoxModelData thisvoxel;

  cout<<"Initializing lists...\n";
  VISVolume<byte> temp = VISVolume<byte>(values.zeroCrossings()); 
  thisvoxel._value=thisvoxel._value_prev=_INSIDE_VALUE;
  thisvoxel._active_pixel=0;
  _INSIDE=_data->push(thisvoxel);  // initialize fixed inside index
  thisvoxel._value=thisvoxel._value_prev=-_INSIDE_VALUE;
  _OUTSIDE=_data->push(thisvoxel); // initialize fixed outside index
#define CHECK_BBOX
#ifdef CHECK_BBOX
  bbox[0][0]=bbox[0][1]=bbox[0][2]=_width+_height+_depth;
  bbox[1][0]=bbox[1][1]=bbox[1][2]=0;
#endif
  _active_list->clean();
  for (x = 0; x < _width; x++) for (y = 0; y < _height; y++) for (z = 0; z < _depth; z++)
    if (temp.itemAt(x, y, z)||(values.itemAt(x, y, z) == 0.0)) {
#ifdef CHECK_BBOX
      if (x>bbox[1][0]) bbox[1][0]=x;
      if (x<bbox[0][0]) bbox[0][0]=x;
      if (y>bbox[1][1]) bbox[1][1]=y;
      if (y<bbox[0][1]) bbox[0][1]=y;
      if (z>bbox[1][2]) bbox[1][2]=z;
      if (z<bbox[0][2]) bbox[0][2]=z;
#endif
      temp.poke(x,y,z)=1;
#ifdef TRUE_DISTANCE_INPUT 
      // use following code only if input image is an exact signed distance transform
      cout<<"WARNING: Using true distance input.\n";
      if (fabs(values.at(x,y,z))>CHANGE_FACTOR) {
	cout<<"Active pixel input value out of range ("<<x<<","<<y<<","<<z<<"):"
	    <<values.at(x,y,z)<<"\n";
	exit(-1);
      }
      _active_list->appendItem(VolIndexValue(x, y, z,values(x,y,z))); 
#else
      // use this code when the input image is a general inside/outside function
      get_neighborhood(x, y, z, neighborhood, values);
      get_derivs_1sided(neighborhood, derivs); 
      d_F = VISVector(VISMaxAbs(derivs[DPX], derivs[DMX]), 
		      VISMaxAbs(derivs[DPY], derivs[DMY]),
		      VISMaxAbs(derivs[DPZ], derivs[DMZ]));
      len = MIN_NORM+sqrt(d_F(0)*d_F(0) + d_F(1)*d_F(1) + d_F(2)*d_F(2));
      dist=values.at(x,y,z)/len;
      dist=VISmin(VISmax(-CHANGE_FACTOR,dist),CHANGE_FACTOR);
      _active_list->appendItem(VolIndexValue(x, y, z, dist)); 
#endif
    }

#ifdef CHECK_BBOX
  cout<<"Active set bounding box ("<<bbox[0][0]<<","<<bbox[0][1]<<","<<
    bbox[0][2]<<") - ("<<bbox[1][0]<<","<<bbox[1][1]<<","<<bbox[1][2]<<")\n";
  cout<<"Volume size: "<<_width<<" "<<_height<<" "<<_depth<<endl;
#endif
  
  thisvoxel._active_pixel=ACTIVE_STATUS;
  _active_list->reset();
  // put active pixels on the list
  while (_active_list->valid()) {
    v_index = _active_list->itemAtCurrent();
    x=v_index.a();y=v_index.b();z=v_index.c();
    thisvoxel._value=thisvoxel._value_prev=v_index.value();
    _index.poke(x,y,z)=_data->push(thisvoxel);
    _active_list->stepForward();
  } 

  // the rest are inside or outside at this point
  for (x = 0; x < _width; x++) for (y = 0; y < _height; y++) for (z = 0; z < _depth; z++)
    if (!temp.itemAt(x, y, z)) {
      if (values.itemAt(x, y, z) > 0.0) _index.poke(x,y,z)=_INSIDE;
      else _index.poke(x,y,z)=_OUTSIDE;
    }
  
  // now fill in the rest of the band
  _inside_list[0].clean();_outside_list[0].clean();
  nextlayer[0]=&_outside_list[0];nextlayer[1]=&_inside_list[0];
  dist=1.0;
  for (j=0;j<2;j++) {
    _active_list->reset();
    dist_sign=j*2-1;  // j=1 is inside (_value positive)
    while (_active_list->valid()) {  
      v_index=_active_list->itemAtCurrent();
      x = v_index.a();y = v_index.b();z=v_index.c();
      if ((x < (_width - 1))&&(x > 0)&&
	  (y < (_height - 1))&&(y > 0)&&
	  (z < (_depth - 1))&&(z > 0))
	for (k = 0; k < NUM_N; k++)
	  if (temp.itemAt(x + Neighbors[k][0], y + Neighbors[k][1], z + Neighbors[k][2]) == 0) {	
	    if ((values.itemAt(x + Neighbors[k][0], y + Neighbors[k][1], z + Neighbors[k][2]) > 0.0f) == j) {
	      nextlayer[j]->appendItem(VISVolIndex(x + Neighbors[k][0], y + Neighbors[k][1], z + Neighbors[k][2]));
	      thisvoxel._value=thisvoxel._value_prev=dist*(float)dist_sign;
	      thisvoxel._active_pixel=ACTIVE_STATUS+dist_sign;
	      _index.poke(x + Neighbors[k][0], y + Neighbors[k][1], z + Neighbors[k][2])=_data->push(thisvoxel);
	      temp.poke(x + Neighbors[k][0], y + Neighbors[k][1], z + Neighbors[k][2])=1;
	    }
	  } 
      _active_list->stepForward();
    }
  }

  for (i=1;i<_number_of_layers;i++) {
    _inside_list[i].clean();_outside_list[i].clean();
    current[1]=&(_inside_list[i-1]);current[0]=&(_outside_list[i-1]);
    nextlayer[0]=&_outside_list[i];nextlayer[1]=&_inside_list[i];
    dist=(float)(i+1);
    for (j=0;j<2;j++) {
      current[j]->reset();
      dist_sign=j*2-1;  // j=1 is inside (_value positive)
      while (current[j]->valid()) {  
	index=current[j]->itemAtCurrent();
	x = index.a();y = index.b();z=index.c();
	if ((x < (_width - 1))&&(x > 0)&&
	    (y < (_height - 1))&&(y > 0)&&
	    (z < (_depth - 1))&&(z > 0))
	  for (k = 0; k < NUM_N; k++)
	    if (temp.itemAt(x + Neighbors[k][0], y + Neighbors[k][1], z + Neighbors[k][2]) == 0) {	
	      nextlayer[j]->appendItem(VISVolIndex(x + Neighbors[k][0], y + Neighbors[k][1], z + Neighbors[k][2]));
	      thisvoxel._value=thisvoxel._value_prev=dist*(float)dist_sign;
	      thisvoxel._active_pixel=ACTIVE_STATUS+i*10+dist_sign;
	      _index.poke(x + Neighbors[k][0], y + Neighbors[k][1], z + Neighbors[k][2])=_data->push(thisvoxel);
	      temp.poke(x + Neighbors[k][0], y + Neighbors[k][1], z + Neighbors[k][2])=1;
	    } 
	current[j]->stepForward();
      }
    }
  }

  for (i=0;i<2*_number_of_layers;i++) {
    _status_up_list[i].clean();
    _status_down_list[i].clean();
  }

  update(0.0); // since we don't set _values to 0 on active_list, we have to update other lists accordingly
}

void LevelSetModel::construct_lists (VISVolume<float> *values) {
  // constructs the lists -- can only be called after full initialization
  unsigned x, y, z, bbox[2][3];
  int k, i, j, dist_sign;
  float len, neighborhood[NSIZE], derivs[NUM_DERIVS], dist;
  VISVector d_F(2);
  VISVolIndex index;    
  VolIndexValue v_index;    
  VISVolIndexVISList *current[2], *nextlayer[2];
  VoxModelData thisvoxel;
  
  cout<<"Initializing lists...\n";
  VISVolume<byte> temp = VISVolume<byte>(values->zeroCrossings()); 
  thisvoxel._value=thisvoxel._value_prev=_INSIDE_VALUE;
  thisvoxel._active_pixel=0;
  _INSIDE=_data->push(thisvoxel);
  thisvoxel._value=thisvoxel._value_prev=-_INSIDE_VALUE;
  _OUTSIDE=_data->push(thisvoxel); 

#ifdef CHECK_BBOX
  bbox[0][0]=bbox[0][1]=bbox[0][2]=_width+_height+_depth;
  bbox[1][0]=bbox[1][1]=bbox[1][2]=0;
#endif
  
  _active_list->clean();
  for (x = 0; x < _width; x++) for (y = 0; y< _height; y++) for (z = 0; z < _depth; z++)
    if (temp.itemAt(x, y, z)||(values->itemAt(x, y, z) == 0.0)) {
#ifdef CHECK_BBOX
      if (x>bbox[1][0]) bbox[1][0]=x;
      if (x<bbox[0][0]) bbox[0][0]=x;
      if (y>bbox[1][1]) bbox[1][1]=y;
      if (y<bbox[0][1]) bbox[0][1]=y;
      if (z>bbox[1][2]) bbox[1][2]=z;
      if (z<bbox[0][2]) bbox[0][2]=z;
#endif
      temp.poke(x,y,z)=1;
#ifdef TRUE_DISTANCE_INPUT 
      // use following code only if input image is an exact signed distance transform
      cout<<"WARNING: Using true distance input.\n";
      if (fabs(values->at(x,y,z))>CHANGE_FACTOR) {
	cout<<"Active pixel input value out of range ("<<x<<","<<y<<","<<z<<"):"
	    <<values->at(x,y,z)<<"\n";
	exit(-1);
      }
      _active_list->appendItem(VolIndexValue(x, y, z,values->peek(x,y,z))); 
#else
      // use this code when the input image is a general inside/outside function
      get_neighborhood (x, y, z, neighborhood, values);
      get_derivs_1sided(neighborhood, derivs); 
      d_F = VISVector(VISMaxAbs(derivs[DPX], derivs[DMX]), 
		      VISMaxAbs(derivs[DPY], derivs[DMY]),
		      VISMaxAbs(derivs[DPZ], derivs[DMZ]));
      len = MIN_NORM+sqrt(d_F(0)*d_F(0) + d_F(1)*d_F(1) + d_F(2)*d_F(2));
      dist=values->at(x,y,z)/len;
      dist=VISmin(VISmax(-CHANGE_FACTOR,dist),CHANGE_FACTOR);
      _active_list->appendItem(VolIndexValue(x, y, z, dist)); 
#endif
    }

#ifdef CHECK_BBOX
  cout<<"Active set bounding box ("<<bbox[0][0]<<","<<bbox[0][1]<<","<<
    bbox[0][2]<<") - ("<<bbox[1][0]<<","<<bbox[1][1]<<","<<bbox[1][2]<<")\n";
  cout<<"Volume size: "<<_width<<" "<<_height<<" "<<_depth<<endl;
#endif 
  
  thisvoxel._active_pixel=ACTIVE_STATUS;
  _active_list->reset();
  while (_active_list->valid()) {
    v_index = _active_list->itemAtCurrent();
    x=v_index.a();y=v_index.b();z=v_index.c();
    thisvoxel._value=thisvoxel._value_prev=v_index.value();
    _index.poke(x,y,z)=_data->push(thisvoxel);
    _active_list->stepForward();
  } 

  for (x = 0; x < _width; x++) for (y = 0; y< _height; y++) for (z = 0; z < _depth; z++)
    if (!temp.itemAt(x, y, z)) {
      if (values->itemAt(x, y, z) > 0.0) _index.poke(x,y,z)=_INSIDE;
      else _index.poke(x,y,z)=_OUTSIDE;
    }

  _inside_list[0].clean();_outside_list[0].clean();
  nextlayer[0]=&_outside_list[0];nextlayer[1]=&_inside_list[0];
  dist=1.0;
  for (j=0;j<2;j++) {
    _active_list->reset();
    dist_sign=j*2-1;  // j=1 is inside (_value positive)
    while (_active_list->valid()) {  
      v_index=_active_list->itemAtCurrent();
      x = v_index.a();y = v_index.b();z=v_index.c();
      if ((x < (_width - 1))&&(x > 0)&&
	  (y < (_height - 1))&&(y > 0)&&
	  (z < (_depth - 1))&&(z > 0))
	for (k = 0; k < NUM_N; k++)
	  if (temp.itemAt(x + Neighbors[k][0], y + Neighbors[k][1], z + Neighbors[k][2]) == 0) {	
	    if ((values->itemAt(x + Neighbors[k][0], y + Neighbors[k][1], z + Neighbors[k][2]) > 0.0f) == j) {
	      nextlayer[j]->appendItem(VISVolIndex(x + Neighbors[k][0], y + Neighbors[k][1], z + Neighbors[k][2]));
	      thisvoxel._value=thisvoxel._value_prev=dist*(float)dist_sign;
	      thisvoxel._active_pixel=ACTIVE_STATUS+dist_sign;
	      _index.poke(x + Neighbors[k][0], y + Neighbors[k][1], z + Neighbors[k][2])=_data->push(thisvoxel);
	      temp.poke(x + Neighbors[k][0], y + Neighbors[k][1], z + Neighbors[k][2])=1;
	    }
	  } 
      _active_list->stepForward();
    }
  }
  
  for (i=1;i<_number_of_layers;i++) {
    _inside_list[i].clean();_outside_list[i].clean();
    current[1]=&(_inside_list[i-1]);current[0]=&(_outside_list[i-1]);
    nextlayer[0]=&_outside_list[i];nextlayer[1]=&_inside_list[i];
    dist=(float)(i+1);
    for (j=0;j<2;j++) {
      current[j]->reset();
      dist_sign=j*2-1;  // j=1 is inside (_value positive)
      while (current[j]->valid()) {  
	index=current[j]->itemAtCurrent();
	x = index.a();y = index.b();z=index.c();
	if ((x < (_width - 1))&&(x > 0)&&
	    (y < (_height - 1))&&(y > 0)&&
	    (z < (_depth - 1))&&(z > 0))
	  for (k = 0; k < NUM_N; k++)
	    if (temp.itemAt(x + Neighbors[k][0], y + Neighbors[k][1], z + Neighbors[k][2]) == 0) {	
	      nextlayer[j]->appendItem(VISVolIndex(x + Neighbors[k][0], y + Neighbors[k][1], z + Neighbors[k][2]));
	      thisvoxel._value=thisvoxel._value_prev=dist*(float)dist_sign;
	      thisvoxel._active_pixel=ACTIVE_STATUS+i*10+dist_sign;
	      _index.poke(x + Neighbors[k][0], y + Neighbors[k][1], z + Neighbors[k][2])=_data->push(thisvoxel);
	      temp.poke(x + Neighbors[k][0], y + Neighbors[k][1], z + Neighbors[k][2])=1;
	    } 
	current[j]->stepForward();
      }
    }
  }
  
  for (i=0;i<2*_number_of_layers;i++) {
    _status_up_list[i].clean();
    _status_down_list[i].clean();
  }
  
  update(0.0); // since we don't set _values to 0 on active_list, we have to update other lists accordingly
}

float LevelSetModel::update(float dt) {
  //dump(246,101,38,1);
  // this function is called after calculate_change to move the level set front
  int x, y, z, k, l, ii=0, total_active_neighbors, dist_sign;
  unsigned m;
  float scale, min_scale, diff1, diff2, new_value, total_change=0.0, change, this_value, neighbor_value;
  byte stat;
  VISVolIndex index;
  VolIndexValue v_index;
  VISVolIndexVISList *status_up[2], *status_down[2], *current[2];
  VoxModelData *thisvoxel, *nvoxel, newvoxel;
  
  _active_list->reset();
  while (_active_list->valid()) {
    v_index = _active_list->itemAtCurrent();
    x = v_index.a();y = v_index.b();z=v_index.c();
    thisvoxel=_data->pop(_index(x,y,z));
    change = _active_list->itemAtCurrent().value();
    if (fabs(dt*change)>=1.0) {
      cout<<x<<","<<y<<","<<z<<" : change is "<<dt*change<<endl;
      exit(-1);
    }
    new_value = thisvoxel->_value + dt*change;
    // constrain does not do anything currently 
    // new_value = constrain(new_value, x, y, z); 
    /*if ( (x<=(_number_of_layers+1)) || (x>=(_width-_number_of_layers-2)) ||
      (y<=(_number_of_layers+1)) || (y>=(_height-_number_of_layers-2)) ||
      (z<=(_number_of_layers+1)) || (z>=(_depth-_number_of_layers-2)) )*/
    if ( (x<=(_border_width+1)) || (x>=(_width-_border_width-1)) ||
	 (y<=(_border_width+1)) || (y>=(_height-_border_width-1)) ||
	 (z<=(_border_width+1)) || (z>=(_depth-_border_width-1)) )
      new_value=VISmin(new_value,CHANGE_FACTOR);
    thisvoxel->_value_prev = thisvoxel->_value;
    thisvoxel->_value = new_value;
    if (new_value > CHANGE_FACTOR) _status_up_list[0].appendItem(VISVolIndex(x, y, z));
    _active_list->stepForward();         
  }

  min_scale=1.0;
  _status_up_list[0].reset();
  while (_status_up_list[0].valid()) {
    index = _status_up_list[0].itemAtCurrent();
    x = index.a();y = index.b();z = index.c();
    thisvoxel=_data->pop(_index(x,y,z));
    for (k=0;k<NUM_N;k++) {
      nvoxel=_data->pop(_index(x + Neighbors[k][0], y + Neighbors[k][1], z + Neighbors[k][2]));
      if ((nvoxel->_active_pixel == ACTIVE_STATUS)&&(nvoxel->_value < (-1.0*CHANGE_FACTOR))) {
	ii++;
	diff1=thisvoxel->_value-nvoxel->_value;
	diff2=thisvoxel->_value_prev-nvoxel->_value_prev;
	if (diff2>=1.0) scale=min_scale=0.0;
	else {
	  scale=(1.0-diff2)/(diff1-diff2);
	  if (scale<min_scale) min_scale=scale;
	}
	thisvoxel->_value=(1.0f-scale)*thisvoxel->_value_prev+scale*thisvoxel->_value;
	nvoxel->_value=(1.0f-scale)*nvoxel->_value_prev+scale*nvoxel->_value;
	//cout<<"local scaling at "<<x<<" "<<y<<" "<<z<<endl;
	if (thisvoxel->_value<=CHANGE_FACTOR) break;
      }
    }
    _status_up_list[0].stepForward();
  }
  
  //if (ii) cout<<ii<<" local scalings. min scale: "<<min_scale<<endl;
  
  ii=0;
  _status_up_list[0].clean();
  _active_list->reset();
  while (_active_list->valid()) { 
    ii++;
    v_index = _active_list->itemAtCurrent();
    x = v_index.a();y = v_index.b();z = v_index.c();
    thisvoxel=_data->pop(_index(x,y,z));

    change=thisvoxel->_value-thisvoxel->_value_prev;
    total_change+=(change*change);
    new_value=thisvoxel->_value;
    if (new_value > CHANGE_FACTOR) { // moving inside 
      _active_list->removeCurrent();
      _status_up_list[0].appendItem(VISVolIndex(x, y, z));
      for (k=0;k<NUM_N;k++) {
	nvoxel=_data->pop(_index(x + Neighbors[k][0], y + Neighbors[k][1], z + Neighbors[k][2]));
	if (nvoxel->_active_pixel==OUTSIDE_STATUS_1) {
	  // this pixel has to be made active now!
	  nvoxel->_active_pixel=CHANGING_STATUS;
	  this_value=new_value-1.0f;
	  if (this_value>nvoxel->_value) nvoxel->_value=this_value; 
	}
      }
    } else if (new_value < (CHANGE_FACTOR*-1.0)) { // moving outside 
      _active_list->removeCurrent();
      _status_down_list[0].appendItem(VISVolIndex(x, y, z));
      for (k=0;k<NUM_N;k++) {
	nvoxel=_data->pop(_index(x + Neighbors[k][0], y + Neighbors[k][1], z + Neighbors[k][2]));
	if (nvoxel->_active_pixel==INSIDE_STATUS_1) {
	  // this pixel has to be made active now!
	  nvoxel->_active_pixel=CHANGING_STATUS;
	  this_value=new_value+1.0f;
	  if (this_value<nvoxel->_value) nvoxel->_value=this_value;
	}
      }
    } else _active_list->stepForward();
  }
  
  //cout<<"Update total change: "<<total_change<<endl;
  total_change = sqrt(total_change/(float)ii);
  
  for (m=0;m<_number_of_layers;m++) {
    current[0]=&_outside_list[m];
    current[1]=&_inside_list[m];
    status_up[0]=&_status_up_list[1+m];
    status_down[0]=&_status_down_list[1+m];
    status_up[1]=&_status_up_list[1+m+_number_of_layers];
    status_down[1]=&_status_down_list[1+m+_number_of_layers];
    for (l=0;l<2;l++) {
      dist_sign=2*l-1;
      current[l]->reset();
      while (current[l]->valid()) {
	index=current[l]->itemAtCurrent();
	x = index.a();y = index.b();z = index.c();
	if ((!x)||(x==(_width-1))||
	    (!y)||(y==(_height-1))||
	    (!z)||(z==(_depth-1))) {
	  cout<<"flowing out of boundaries\n";
	  cout<<m<<" "<<l<<endl;
	  cout<<x<<","<<y<<","<<z<<endl;
	  VISImageFile imfile;
	  imfile.write(((status()).image()).becomeFlat(), "layers.fit");
	  exit(-1);
	}
	
	thisvoxel=_data->pop(_index(x,y,z));
	if (thisvoxel->_active_pixel==CHANGING_STATUS) {
	  if (m) {
	    cout<<"m>0 for changing_status pixel!\n";
	    exit(-1);
	  }
	  current[l]->removeCurrent();
	  if (l) status_down[1]->appendItem(VISVolIndex(x, y, z)); 
	  else status_up[0]->appendItem(VISVolIndex(x, y, z));
	  continue; // go directly to next pixel on the list
	}
		
	total_active_neighbors = 0;
	if (l) neighbor_value = FLT_MAX; else neighbor_value = -FLT_MAX;
	for (k = 0; k < NUM_N; k++) {
	  nvoxel=_data->pop(_index(x + Neighbors[k][0], y + Neighbors[k][1], z + Neighbors[k][2]));
	  if (check_layer(m,l,nvoxel->_active_pixel)) {
	    this_value = nvoxel->_value;
	    if (l) {if (this_value < neighbor_value) neighbor_value = this_value; }
	    else if (this_value > neighbor_value) neighbor_value = this_value;
	    total_active_neighbors++;
	  }
	}
	
	if (total_active_neighbors == 0) {
	  // check to see when this happens
	  if (m==(_number_of_layers-1)) {
	    _data->remove(_index(x,y,z));
	    if (l) _index.poke(x,y,z)=_INSIDE; else _index.poke(x,y,z)=_OUTSIDE;
	    current[l]->removeCurrent();
	  } else {
	    current[l]->removeCurrent();
	    thisvoxel->_active_pixel=ACTIVE_STATUS+(m+1)*10+dist_sign;
	    if (l) _inside_list[m+1].appendItem(VISVolIndex(x, y, z)); 
	    else _outside_list[m+1].appendItem(VISVolIndex(x, y, z)); 
	  }
	} else thisvoxel->_value = neighbor_value + (float)dist_sign;
	current[l]->stepForward();
      }

      current[l]->reset();
      while (current[l]->valid()) {
	index=current[l]->itemAtCurrent();
	x = index.a();y = index.b();z = index.c();
	thisvoxel=_data->pop(_index(x,y,z));
	
	total_active_neighbors = 0;
	neighbor_value = thisvoxel->_value;
	for (k = 0; k < NUM_N; k++) {
	  nvoxel=_data->pop(_index(x + Neighbors[k][0], y + Neighbors[k][1], z + Neighbors[k][2]));
	  if (check_layer2(m,l,nvoxel->_active_pixel)) {
	    this_value = nvoxel->_value+ (float)dist_sign;
	    if (l) {if (this_value < neighbor_value) neighbor_value = this_value; }
	    else if (this_value > neighbor_value) neighbor_value = this_value;
	    total_active_neighbors++;
	  }
	}
	
	thisvoxel->_value = neighbor_value;
	if (l==0) {
	  // outside layer
	  if (m<(_number_of_layers-1)) {
	    if (thisvoxel->_value < (-CHANGE_FACTOR * (float)(2*m+3))) {
	      current[l]->removeCurrent();
	      status_down[l]->appendItem(VISVolIndex(x, y, z));
	    } else if (thisvoxel->_value >= (-CHANGE_FACTOR * (float)(2*m+1))) {
	      current[l]->removeCurrent();
	      status_up[l]->appendItem(VISVolIndex(x, y, z));
	    } else current[l]->stepForward();
	  } else {
	    if (thisvoxel->_value< -_INSIDE_VALUE) {
	      current[l]->removeCurrent();
	      _data->remove(_index(x,y,z));
	      _index.poke(x,y,z)=_OUTSIDE;
	    } else if (thisvoxel->_value >= (-CHANGE_FACTOR * (float)(2*m+1))) {
	      current[l]->removeCurrent();
	      status_up[l]->appendItem(VISVolIndex(x, y, z));
	      for (k=0; k < NUM_N; k++) {
		nvoxel=_data->pop(_index(x + Neighbors[k][0], y + Neighbors[k][1], z + Neighbors[k][2]));
		if ((nvoxel->_active_pixel==0)&&(nvoxel->_value<=-_INSIDE_VALUE)) { 
		  //newvoxel._value=newvoxel._value_prev=-_INSIDE_VALUE;
		  newvoxel._value=newvoxel._value_prev=thisvoxel->_value-1.0f;
		  newvoxel._active_pixel=CHANGING_STATUS;
		  _index.poke(x + Neighbors[k][0], y + Neighbors[k][1], z + Neighbors[k][2])=_data->push(newvoxel);
		  status_down[l]->appendItem(VISVolIndex(x+Neighbors[k][0],y+Neighbors[k][1],z+Neighbors[k][2]));
		}
	      }
	    } else current[l]->stepForward();
	  }
	} else {
	  // inside layer
	  if (m<(_number_of_layers-1)) {
	    if (thisvoxel->_value > (CHANGE_FACTOR * (float)(2*m+3))) {
	      current[l]->removeCurrent();
	      status_up[l]->appendItem(VISVolIndex(x, y, z));
	    } else if (thisvoxel->_value <= (CHANGE_FACTOR * (float)(2*m+1))) {
	      current[l]->removeCurrent();
	      status_down[l]->appendItem(VISVolIndex(x, y, z));
	    } else current[l]->stepForward();
	  } else {
	    if (thisvoxel->_value > _INSIDE_VALUE) {
	      current[l]->removeCurrent();
	      _data->remove(_index(x,y,z));
	      _index.poke(x,y,z)=_INSIDE;
	    } else if (thisvoxel->_value <= (CHANGE_FACTOR * (float)(2*m+1))) {
	      current[l]->removeCurrent();
	      status_down[l]->appendItem(VISVolIndex(x, y, z));
	      for (k=0; k < NUM_N; k++) {
		nvoxel=_data->pop(_index(x + Neighbors[k][0], y + Neighbors[k][1], z + Neighbors[k][2]));
		if ((nvoxel->_active_pixel==0)&&(nvoxel->_value>=_INSIDE_VALUE)) {
		  //newvoxel._value=newvoxel._value_prev=_INSIDE_VALUE;
		  newvoxel._value=newvoxel._value_prev=thisvoxel->_value+1.0f;
		  newvoxel._active_pixel=CHANGING_STATUS;
		  _index.poke(x + Neighbors[k][0], y + Neighbors[k][1], z + Neighbors[k][2])=_data->push(newvoxel);
		  status_up[l]->appendItem(VISVolIndex(x + Neighbors[k][0], y + Neighbors[k][1], z + Neighbors[k][2]));
		}
	      }
	    } else current[l]->stepForward();
	  }
	}
      }
    }
  }
  
  _status_up_list[0].reset();    
  while (_status_up_list[0].valid()) {
    index = _status_up_list[0].itemAtCurrent();
    x = index.a();y = index.b();z=index.c();
    thisvoxel=_data->pop(_index(x,y,z));
    thisvoxel->_active_pixel=ACTIVE_STATUS+1;
    _inside_list[0].appendItem(VISVolIndex(x, y, z));
    _status_up_list[0].removeCurrent();
    }
 
  _status_down_list[0].reset();
  while (_status_down_list[0].valid()) {
    index = _status_down_list[0].itemAtCurrent();
    x = index.a();y = index.b();z=index.c();
    thisvoxel=_data->pop(_index(x,y,z));
    thisvoxel->_active_pixel=ACTIVE_STATUS-1;
    _outside_list[0].appendItem(VISVolIndex(x, y, z));
    _status_down_list[0].removeCurrent();
  }

  for (l=0;l<2;l++) {
    for (m=0;m<_number_of_layers;m++) {
      status_up[0]=&_status_up_list[1+m];
      status_down[0]=&_status_down_list[1+m];
      status_up[1]=&_status_up_list[1+m+_number_of_layers];
      status_down[1]=&_status_down_list[1+m+_number_of_layers];
      
      status_up[l]->reset();
      while (status_up[l]->valid()) {
	index=status_up[l]->itemAtCurrent();
	x = index.a();y = index.b();z=index.c();
	thisvoxel=_data->pop(_index(x,y,z));
	if (l) {
	  // inside layer
	  k=VISmin(m+1,_number_of_layers-1);
	  stat=ACTIVE_STATUS+k*10+1;
	  _inside_list[k].appendItem(VISVolIndex(x, y, z));
	} else {
	  // outside layer
	  if (m) {
	    stat=ACTIVE_STATUS+(m-1)*10-1; 
	    _outside_list[m-1].appendItem(VISVolIndex(x, y, z));
	  } else {
	    stat=ACTIVE_STATUS;
	    _active_list->appendItem(VISVolIndex(x, y, z));
	  }
	}
	thisvoxel->_active_pixel = stat;
	status_up[l]->removeCurrent();
      }
      status_down[l]->reset();
      while (status_down[l]->valid()) {
	index=status_down[l]->itemAtCurrent();
	x = index.a();y = index.b();z=index.c();
	thisvoxel=_data->pop(_index(x,y,z));
	if (l) {
	  // inside layer
	  if (m) {
	    stat=ACTIVE_STATUS+(m-1)*10+1; 
	    _inside_list[m-1].appendItem(VISVolIndex(x, y, z));
	  } else {
	    stat=ACTIVE_STATUS;
	    _active_list->appendItem(VISVolIndex(x, y, z));
	  }
	} else {
	  // outside layer
	  k=VISmin(m+1,_number_of_layers-1);
	  stat=ACTIVE_STATUS+k*10-1;
	  _outside_list[k].appendItem(VISVolIndex(x, y, z));
	}
	thisvoxel->_active_pixel = stat;
	status_down[l]->removeCurrent();
      }
    }
  }
  //cout<<"after update:\n";dump(246,101,38,1);
  return (total_change);
}

inline int LevelSetModel::get_phi_neighborhood (int x, int y, int z, float* data) {
  int i, j, k, l = 0, ret=0;
  VoxModelData *thisvoxel;
  
  for (k = -1; k <= 1; k++) for (j = -1; j <= 1; j++) for (i = -1; i <= 1; i++) {
    thisvoxel=_data->pop(_index(VISmax(VISmin(x + i, _width - 1), 0), 
			 VISmax(VISmin(y + j, _height - 1), 0),
			 VISmax(VISmin(z + k, _depth - 1), 0)));
    data[l++]=thisvoxel->_value;
    if (!(thisvoxel->_active_pixel)) ret=-(abs(i)+abs(j)+abs(k));
  }        
  return ret;
} 

float LevelSetModel::distanceIterate (int num_extra_levels, 
				      VISVolIndexVISList *buff1, VISVolIndexVISList *buff2,
				      VISVolume<float> &val_exact, VISVolume<float> &val_prev,
				      VISVolume<boolean> &active_pixels) {
  int m, l, k, p, r, x, y, z, xn, yn, zn, cnt, nl, total_active_neighbors;
  VISVolIndexVISList *current[2];
  VISVolIndex index;
  float this_value, neighbor_value[6], sum, sq_sum, disc, total_change, change;
  
  val_prev=val_exact;

  total_change=0.0;
  cnt=0;
  nl=_number_of_layers+num_extra_levels; 
  for (m=0;m<nl;m++) {
    if (m<_number_of_layers) {
      current[0]=&_outside_list[m];
      current[1]=&_inside_list[m];
    } else {
      current[0]=&(buff1[m-_number_of_layers]);
      current[1]=&(buff2[m-_number_of_layers]);
    }
    for (l=0;l<2;l++) { // l=0 is outside , l=1 is inside
      current[l]->reset();
      while (current[l]->valid()) {
	index=current[l]->itemAtCurrent();
	x = index.a();y = index.b();z = index.c();
	total_active_neighbors = 0;
	
	for (k=0;k<6;k++) neighbor_value[k] = val_prev(x,y,z);
	for (k = 0; k < 6; k++) {
	  xn=x + Neighbors[k][0]; 
	  yn=y + Neighbors[k][1];
	  zn=z + Neighbors[k][2];
	  if (active_pixels(xn,yn,zn)!=0) {
	    this_value = val_prev(xn,yn,zn);
	    if (l) {
	      for (p=0;p<6;p++) if (this_value < neighbor_value[p]) {
		for (r=5;r>p;r--) neighbor_value[r]=neighbor_value[r-1];
		neighbor_value[p]=this_value;
		total_active_neighbors++;
		break;
	      }
	    } else {
	      for (p=0;p<6;p++) if (this_value > neighbor_value[p]) {
		for (r=5;r>p;r--) neighbor_value[r]=neighbor_value[r-1];
		neighbor_value[p]=this_value;
		total_active_neighbors++;
		break;
	      }
	    }
	  }
	}
	
	if ((m<_number_of_layers)&&(total_active_neighbors == 0)) {
	  cout<<"distanceIterate: no valid neighbours("<<x<<","<<y<<","<<z<<")\n";
	  cout<<"m: "<<m<<" l: "<<l<<endl;
	  int i2, j2, k2;
	  
	  for (k2=-1;k2<=1;k2++) { for (j2=-1;j2<=1;j2++) { for (i2=-1;i2<=1;i2++) {
	    cout<<((int)(active_pixels(x+i2,y+j2,z+k2)))<<" ";
	  } cout<<endl; } cout<<endl; }
	  cout<<"----------\n";
	  for (k2=-1;k2<=1;k2++) { for (j2=-1;j2<=1;j2++) { for (i2=-1;i2<=1;i2++) {
	    cout<<val_prev(x+i2,y+j2,z+k2)<<" ";
	  } cout<<endl; } cout<<endl; }
	  cout<<"----------\n";
	  for (k2=-1;k2<=1;k2++) { for (j2=-1;j2<=1;j2++) { for (i2=-1;i2<=1;i2++) {
	    cout<<val_exact(x+i2,y+j2,z+k2)<<" ";
	  } cout<<endl; } cout<<endl; }
	  
	  //dump(x,y,z,TRUE);
	  //exit(-1);
	} else {
	  sum=sq_sum=0.0;
	  for (p=0;p<total_active_neighbors;p++) {
	    sum+=neighbor_value[p];
	    sq_sum+=neighbor_value[p]*neighbor_value[p];
	  }
	  
	  for (p=total_active_neighbors;p>=1;p--) // use the largest number of consistent neighbors possible
	    if (p>1) { 
	      disc=sum*sum-(sq_sum-1.0)*(float)p;
	      if (disc>=0.0) { // found a solution
		val_exact.poke(x,y,z)=(2.0*sum+2.0*sqrt(disc)*(float)(2*l-1))/(2.0*(float)p);
		break;
	      } else { // use less neighbors
		sum-=neighbor_value[p-1];
		sq_sum-=(neighbor_value[p-1]*neighbor_value[p-1]);
	      }
	    }
	    else val_exact.poke(x,y,z)=neighbor_value[0]+(float)(2*l-1);
	  
	  change=val_prev(x,y,z)-val_exact(x,y,z);
	  total_change+=change*change;
	  cnt++;
	}
	current[l]->stepForward();
      }
    }
  }
  return sqrt(total_change/cnt);
}

VISVolume<float> LevelSetModel::calculateDistanceTransform (VISVolume<float> &GT) {
  // this function is extremely slow! 
  VISVolume<float> val_exact (_width,_height,_depth);
  int x, y, z, xa, ya, za;
  VolIndexValue index;
  VoxModelData *thisvoxel, *activevoxel;
  float min, max, val, dist;
  
  for (x=0;x<_width;x++) for (y=0;y<_height;y++) { cout<<x<<" "<<y<<endl; for (z=0;z<_depth;z++) {
    _active_list->reset();
    max=1.0e+6;
    min=-1.0e+6;
    thisvoxel=_data->pop(_index(x,y,z));
    while (_active_list->valid()) { 
      index=_active_list->itemAtCurrent();
      xa = index.a();ya = index.b(); za = index.c();
      dist=sqrt((float)((xa-x)*(xa-x)+(ya-y)*(ya-y)+(za-z)*(za-z)));
      activevoxel=_data->pop(_index(xa,ya,za));
      if (thisvoxel->_value>=0.0f) {
	val=activevoxel->_value+dist;
	if (val<max) max=val;
      } else {
	val=activevoxel->_value-dist;
	if (val>min) min=val;
      }
      _active_list->stepForward();
    }
    if (thisvoxel->_value>=0.0f) val_exact.poke(x,y,z)=max;
    else val_exact.poke(x,y,z)=min;
  }
  }
  int acnt=0;
  float diff=0.0f;
  for (x=0;x<_width;x++) for (y=0;y<_height;y++) for (z=0;z<_depth;z++) {
    acnt++;
    diff+=((GT(x,y,z)-val_exact(x,y,z))*(GT(x,y,z)-val_exact(x,y,z)));
  }
  cout<<"rms error: "<<sqrt(diff/(float)acnt)<<endl;
  return val_exact;
}

VISVolume<float> LevelSetModel::calculateDistanceTransform (int num_extra_levels, 
							    VISVolume<float> &GT) { 
  VISVolume<float> val_exact (_width,_height,_depth), val_prev (_width,_height,_depth);
  VISVolume<boolean> active_pixels (_width,_height,_depth);
  int x, y, z, xn, yn, zn, m, l, k, i;
  int max_it=(_number_of_layers+num_extra_levels+1)*3;
  VoxModelData *thisvoxel;
  VISVolIndexVISList *current[2], *inside_buffer, *outside_buffer;
  VISVolIndex index;
  float dist, change;

  for (x=0;x<_width;x++) for (y=0;y<_height;y++) for (z=0;z<_depth;z++) {
    thisvoxel=_data->pop(_index(x,y,z));
    val_exact.poke(x,y,z)=thisvoxel->_value;
    active_pixels.poke(x,y,z)=thisvoxel->_active_pixel;
  }

  if (num_extra_levels>0) {
    inside_buffer=new VISVolIndexVISList [num_extra_levels];
    outside_buffer=new VISVolIndexVISList [num_extra_levels];
    
    for (m=0;m<num_extra_levels;m++) {
      if (m) {
	current[1]=&(inside_buffer[m-1]);
	current[0]=&(outside_buffer[m-1]);
      } else {
	current[1]=&(_inside_list[_number_of_layers-1]);
	current[0]=&(_outside_list[_number_of_layers-1]);
      } 
      inside_buffer[m].clean();
      outside_buffer[m].clean();
      for (l=0;l<2;l++) {
	if (l) dist=FLT_MAX; else dist=-FLT_MAX;
	current[l]->reset();
	while (current[l]->valid()) { 
	  index=current[l]->itemAtCurrent();
	  x = index.a();y = index.b(); z=index.c();
	  if ((x < (_width - 1))&&(x > 0)&&
	      (y < (_height - 1))&&(y > 0)&&
	      (z < (_depth - 1))&&(z > 0))
	    for (k = 0; k < 6; k++) {
	      xn=x + Neighbors[k][0]; 
	      yn=y + Neighbors[k][1];
	      zn=z + Neighbors[k][2];
	      if (active_pixels(xn,yn,zn) == 0) {
		if (l) inside_buffer[m].appendItem(VISVolIndex(xn,yn,zn));
		else outside_buffer[m].appendItem(VISVolIndex(xn,yn,zn));
		
		active_pixels.poke(xn,yn,zn)=CHANGING_STATUS;
		val_exact.poke(xn,yn,zn)=dist;
	      }
	    }
	  current[l]->stepForward();
	}
      }
    }
  }
  change=1.0;
  i=0;
  int acnt;
  while((change>0.0001)&&(i<max_it)) {
    change=distanceIterate(num_extra_levels,outside_buffer,inside_buffer,val_exact,val_prev,active_pixels);
    cout<<"change: "<<change;
    acnt=0;
    float diff=0.0;
    for (x=0;x<_width;x++) for (y=0;y<_height;y++) for (z=0;z<_depth;z++) 
      if (active_pixels(x,y,z)) {
	acnt++;
	diff+=((GT(x,y,z)-val_exact(x,y,z))*(GT(x,y,z)-val_exact(x,y,z)));
      }
    cout<<" , rms error: "<<sqrt(diff/(float)acnt)<<endl;
    i++;
  }
  if (num_extra_levels>0) {
    for (m=0;m<num_extra_levels;m++) {
      inside_buffer[m].clean();
      outside_buffer[m].clean(); 
    }
      
    delete [] inside_buffer;
    delete [] outside_buffer;
  }
  return val_exact;
}

VISVolume<float> LevelSetModel::calculateDistanceTransform (int num_extra_levels) { 
  // warning: this function uses a lot of memory! 
  VISVolume<float> val_exact (_width,_height,_depth), val_prev (_width,_height,_depth);
  VISVolume<boolean> active_pixels (_width,_height,_depth);
  int x, y, z, xn, yn, zn, m, l, k, i;
  int max_it=(_number_of_layers+num_extra_levels+1)*3;
  VoxModelData *thisvoxel;
  VISVolIndexVISList *current[2], *inside_buffer, *outside_buffer;
  VISVolIndex index;
  float dist, change;

  for (x=0;x<_width;x++) for (y=0;y<_height;y++) for (z=0;z<_depth;z++) {
    thisvoxel=_data->pop(_index(x,y,z));
    val_exact.poke(x,y,z)=thisvoxel->_value;
    active_pixels.poke(x,y,z)=thisvoxel->_active_pixel;
  }

  if (num_extra_levels>0) {
    inside_buffer=new VISVolIndexVISList [num_extra_levels];
    outside_buffer=new VISVolIndexVISList [num_extra_levels];
    
    for (m=0;m<num_extra_levels;m++) {
      if (m) {
	current[1]=&(inside_buffer[m-1]);
	current[0]=&(outside_buffer[m-1]);
      } else {
	current[1]=&(_inside_list[_number_of_layers-1]);
	current[0]=&(_outside_list[_number_of_layers-1]);
      } 
      inside_buffer[m].clean();
      outside_buffer[m].clean();
      for (l=0;l<2;l++) {
	if (l) dist=FLT_MAX; else dist=-FLT_MAX;
	current[l]->reset();
	while (current[l]->valid()) { 
	  index=current[l]->itemAtCurrent();
	  x = index.a();y = index.b(); z=index.c();
	  if ((x < (_width - 1))&&(x > 0)&&
	      (y < (_height - 1))&&(y > 0)&&
	      (z < (_depth - 1))&&(z > 0))
	    for (k = 0; k < 6; k++) {
	      xn=x + Neighbors[k][0]; 
	      yn=y + Neighbors[k][1];
	      zn=z + Neighbors[k][2];
	      if (active_pixels(xn,yn,zn) == 0) {
		if (l) inside_buffer[m].appendItem(VISVolIndex(xn,yn,zn));
		else outside_buffer[m].appendItem(VISVolIndex(xn,yn,zn));
		
		active_pixels.poke(xn,yn,zn)=CHANGING_STATUS;
		val_exact.poke(xn,yn,zn)=dist;
	      }
	    }
	  current[l]->stepForward();
	}
      }
    }
  }
  change=1.0;
  i=0;
  while((change>0.0001)&&(i<max_it)) {
    change=distanceIterate(num_extra_levels,outside_buffer,inside_buffer,val_exact,val_prev,active_pixels);
    i++;
  }
  
  if (num_extra_levels>0) {
    for (m=0;m<num_extra_levels;m++) {
      inside_buffer[m].clean();
      outside_buffer[m].clean(); 
    }
      
    delete [] inside_buffer;
    delete [] outside_buffer;
  }
  return val_exact;
}

// -------------------- VOXMODEL -----------------------------

void VoxModel::set_default_parameters () {
  // the subclass should set some of these to non-zero values
  _energy_weight=0.0f;
  _normal_energy_weight=0.0f;
  _grow_weight=0.0f;
  _const_grow_weight=0.0f;
  _curve_weight=0.0f;
  _curvature_grow_weight=0.0f;
  _weighted_curve_weight=0.0f;
  _gauss_curve_weight=0.0f;
  _special_curve_weight=0.0f;
  _pos_curve_weight=0.0f;
  _neg_curve_weight=0.0f;

  _wave_dt=WAVE_DT;
  _curve_dt=CURVE_DT;
}

void VoxModel::compute_normals_from_phi_on_active_set () {
  VISVolIndex index;
  VolIndexValue v_index;
  VISVolIndexVISList *list;
  int i, j, x, y, z, x2, y2, z2, k, l;
  float phi[2][2][2], mag, N[3];
  VoxModelData *thisvoxel;
  
  _active_list->reset();
  while (_active_list->valid()) {
    v_index = _active_list->itemAtCurrent();
    x=v_index.a();y=v_index.b();z=v_index.c();
    thisvoxel=_data->pop(_index(x,y,z));
    thisvoxel->_flag=0;
    _active_list->stepForward();
  }
  for (j=-1;j<=1;j+=2) for (i=0;i<_number_of_layers;i++) {
    if (j>0) list=&(_inside_list[i]); else list=&(_outside_list[i]);
    list->reset();
    while (list->valid()) {
      index = list->itemAtCurrent();
      x=index.a();y=index.b();z=index.c();
      thisvoxel=_data->pop(_index(x,y,z));
      thisvoxel->_flag=0;
      list->stepForward();
    }
  }
  
  _active_list->reset();
  while (_active_list->valid()) {
    v_index = _active_list->itemAtCurrent();
    x=v_index.a();y=v_index.b();z=v_index.c();
    for (l=-1;l<NUMC;l++) {
      if (l>=0) {
	x2=x+NeighC[l][0];y2=y+NeighC[l][1];z2=z+NeighC[l][2];
      } else {
	x2=x;y2=y;z2=z;
      }
      thisvoxel=_data->pop(_index(x2,y2,z2));
      if (!(thisvoxel->_flag)) {
	thisvoxel->_flag=1;
	for (i=0;i<=1;i++) for (j=0;j<=1;j++) for (k=0;k<=1;k++) 
	  phi[i][j][k]=(_data->pop(_index(x2-i,y2-j,z2-k)))->_value;
	
	N[0]=0.25f*((phi[0][0][0]+phi[0][1][0]+phi[0][0][1]+phi[0][1][1])-
		    (phi[1][0][0]+phi[1][1][0]+phi[1][0][1]+phi[1][1][1]));
	N[1]=0.25f*((phi[0][0][0]+phi[1][0][0]+phi[0][0][1]+phi[1][0][1])-
		    (phi[0][1][0]+phi[1][1][0]+phi[0][1][1]+phi[1][1][1]));
	N[2]=0.25f*((phi[0][0][0]+phi[1][0][0]+phi[0][1][0]+phi[1][1][0])-
		    (phi[0][0][1]+phi[1][0][1]+phi[0][1][1]+phi[1][1][1]));
	mag=MIN_NORM+sqrt(N[0]*N[0]+N[1]*N[1]+N[2]*N[2]);
	(thisvoxel->_Nphi[0])=N[0]/mag;
	(thisvoxel->_Nphi[1])=N[1]/mag;
	(thisvoxel->_Nphi[2])=N[2]/mag; 
      }
    }
    _active_list->stepForward();
  }
}

inline float VoxModel::calculate_curvature (int x, int y, int z) {
  VoxModelData *pos[2][2][2];
  pos[0][0][0]=_data->pop(_index(x,y,z));
  pos[0][0][1]=_data->pop(_index(x,y,z+1));
  pos[0][1][0]=_data->pop(_index(x,y+1,z));
  pos[0][1][1]=_data->pop(_index(x,y+1,z+1));
  pos[1][0][0]=_data->pop(_index(x+1,y,z));
  pos[1][0][1]=_data->pop(_index(x+1,y,z+1));
  pos[1][1][0]=_data->pop(_index(x+1,y+1,z));
  pos[1][1][1]=_data->pop(_index(x+1,y+1,z+1));

  return 0.25f*((pos[1][0][0]->_Nphi[0]+pos[1][1][0]->_Nphi[0]+
		 pos[1][0][1]->_Nphi[0]+pos[1][1][1]->_Nphi[0])-
		(pos[0][0][0]->_Nphi[0]+pos[0][1][0]->_Nphi[0]+
		 pos[0][0][1]->_Nphi[0]+pos[0][1][1]->_Nphi[0])+
		(pos[0][1][0]->_Nphi[1]+pos[1][1][0]->_Nphi[1]+
		 pos[0][1][1]->_Nphi[1]+pos[1][1][1]->_Nphi[1])-
		(pos[0][0][0]->_Nphi[1]+pos[1][0][0]->_Nphi[1]+
		 pos[0][0][1]->_Nphi[1]+pos[1][0][1]->_Nphi[1])+
		(pos[0][0][1]->_Nphi[2]+pos[1][0][1]->_Nphi[2]+
		 pos[0][1][1]->_Nphi[2]+pos[1][1][1]->_Nphi[2])-
		(pos[0][0][0]->_Nphi[2]+pos[1][0][0]->_Nphi[2]+
		 pos[0][1][0]->_Nphi[2]+pos[1][1][0]->_Nphi[2]));
}

inline float VoxModel::calculate_curvature (const float *n) {
  float N[2][2][2][3], mag;
  int x, y, z, i;

  for (x=0;x<=1;x++) for (y=0;y<=1;y++) for (z=0;z<=1;z++) {
    N[x][y][z][0]=((n[x+1+y*3+z*9]+n[x+1+y*3+(z+1)*9]+n[x+1+(y+1)*3+z*9]+n[x+1+(y+1)*3+(z+1)*9])-
		   (n[x+y*3+z*9]+n[x+y*3+(z+1)*9]+n[x+(y+1)*3+z*9]+n[x+(y+1)*3+(z+1)*9]));
    N[x][y][z][1]=((n[x+(y+1)*3+z*9]+n[x+(y+1)*3+(z+1)*9]+n[x+1+(y+1)*3+z*9]+n[x+1+(y+1)*3+(z+1)*9])-
		   (n[x+y*3+z*9]+n[x+y*3+(z+1)*9]+n[x+1+y*3+z*9]+n[x+1+y*3+(z+1)*9]));
    N[x][y][z][2]=((n[x+y*3+(z+1)*9]+n[x+1+y*3+(z+1)*9]+n[x+(y+1)*3+(z+1)*9]+n[x+1+(y+1)*3+(z+1)*9])-
		   (n[x+y*3+z*9]+n[x+1+y*3+z*9]+n[x+(y+1)*3+z*9]+n[x+1+(y+1)*3+z*9]));
    mag=MIN_NORM+sqrt(N[x][y][z][0]*N[x][y][z][0]+
		      N[x][y][z][1]*N[x][y][z][1]+
		      N[x][y][z][2]*N[x][y][z][2]);
    N[x][y][z][0]/=mag;
    N[x][y][z][1]/=mag;
    N[x][y][z][2]/=mag;
  } 
  
  return 0.25f*((N[1][0][0][0]+N[1][1][0][0]+N[1][0][1][0]+N[1][1][1][0])-
		(N[0][0][0][0]+N[0][1][0][0]+N[0][0][1][0]+N[0][1][1][0])+
		(N[0][1][0][1]+N[1][1][0][1]+N[0][1][1][1]+N[1][1][1][1])-
		(N[0][0][0][1]+N[1][0][0][1]+N[0][0][1][1]+N[1][0][1][1])+
		(N[0][0][1][2]+N[1][0][1][2]+N[0][1][1][2]+N[1][1][1][2])-
		(N[0][0][0][2]+N[1][0][0][2]+N[0][1][0][2]+N[1][1][0][2]));
}

 inline VISMatrix VoxModel::calculate_curvature (const float *derivatives, 
						float &curve_trace, float &curve_norm) {
  VISMatrix curve(3, 3);
  VISVector Nf(3), Nb(3);

  // do the x derivative
  Nf.poke(0) = derivatives[DPX];
  Nf.poke(1) = 0.5f*(derivatives[DYPX] + derivatives[DY]);
  Nf.poke(2) = 0.5f*(derivatives[DZPX] + derivatives[DZ]);
  Nf = Nf/(Nf.norm()+MIN_NORM);
  Nb.poke(0) = derivatives[DMX];
  Nb.poke(1) = 0.5f*(derivatives[DYMX] + derivatives[DY]);
  Nb.poke(2) = 0.5f*(derivatives[DZMX] + derivatives[DZ]);
  Nb = Nb/(Nb.norm()+MIN_NORM);
  curve.pokeROI(0, 0, Nf - Nb);

  // do the y derivative
  Nf.poke(0) = 0.5f*(derivatives[DXPY] + derivatives[DX]);
  Nf.poke(1) = derivatives[DPY];
  Nf.poke(2) = 0.5f*(derivatives[DZPY] + derivatives[DZ]);
  Nf = Nf/(Nf.norm()+MIN_NORM);
  Nb.poke(0) = 0.5f*(derivatives[DXMY] + derivatives[DX]);
  Nb.poke(1) = derivatives[DMY];
  Nb.poke(2) = 0.5f*(derivatives[DZMY] + derivatives[DZ]);
  Nb = Nb/(Nb.norm()+MIN_NORM);
  curve.pokeROI(0, 1, Nf - Nb);

  // do the z derivative
  Nf.poke(0) = 0.5f*(derivatives[DXPZ] + derivatives[DX]);
  Nf.poke(1) = 0.5f*(derivatives[DYPZ] + derivatives[DY]);
  Nf.poke(2) = derivatives[DPZ];
  Nf = Nf/(Nf.norm()+MIN_NORM);
  Nb.poke(0) = 0.5f*(derivatives[DXMZ] + derivatives[DX]);
  Nb.poke(1) = 0.5f*(derivatives[DYMZ] + derivatives[DY]);
  Nb.poke(2) = derivatives[DMZ];
  Nb = Nb/(Nb.norm()+MIN_NORM);
  curve.pokeROI(0, 2, Nf - Nb);
  
  curve_trace = curve.peek(0, 0) + curve.peek(1, 1) + curve.peek(2, 2);
  curve_norm = 0.0f;
  for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++)
    curve_norm += pow(curve.peek(i, j), 2);
  return(curve);
}

float VoxModel::calculate_change()  {
  float max_total_force, max_total_grow, total_force_abs, total_change = 0.0f;
  float the_change, force_change, grow_change, curve_norm, curve_trace, total_curve, total_grow;
  float dt, grad, len, k1, k2, curve_gauss; 
  float grad_p, grad_m, grad_x, grad_y, grad_z;
  int x, y, z, k, ii; 
  VolIndexValue v_index; 
  float neighborhood[NSIZE], derivs[NUM_DERIVS];
  float tmp, rms_grow=0.0f, rms_curve=0.0f;
  int cnt=0;
  boolean calc_curve, fast_calc_curve, calc_force, calc_grow, calc_var_grow;
  VoxModelData *thisvoxel;
  VISVector p(3), the_force(3), pc(3), d_F(3);

  max_total_force = max_total_grow = 0.0f;

  calc_curve =   ((_curve_weight != 0.0f)||
		  (_weighted_curve_weight != 0.0f)||
		  (_gauss_curve_weight != 0.0f)||
		  (_pos_curve_weight != 0.0f)||
		  (_neg_curve_weight != 0.0f)||
		  (_curvature_grow_weight != 0.0f)||
		  (_special_curve_weight != 0.0f));
  fast_calc_curve=(((_curve_weight != 0.0f)||(_curvature_grow_weight != 0.0f))&&
		   (_weighted_curve_weight == 0.0f)&&
		   (_gauss_curve_weight == 0.0f)&&
		   (_pos_curve_weight == 0.0f)&&
		   (_neg_curve_weight == 0.0f)&&
		   (_special_curve_weight == 0.0f));
  calc_force =    (_energy_weight != 0.0f)||(_normal_energy_weight != 0.0f);
  calc_grow =     ((_grow_weight != 0.0f)||(_const_grow_weight != 0.0f));
  calc_var_grow = (_grow_weight != 0.0f);
  
  fast_calc_curve=0;
  //cout<<_grow_weight<<" "<<_curve_weight<<" "<<_curvature_grow_weight<<endl;
  //cout<<"calculation flags: "<<calc_curve<<" "<<fast_calc_curve<<" "
  //<<calc_force<<" "<<calc_grow<<" "<<calc_var_grow<<endl;
  _active_list->reset();
  while (_active_list->valid()) {
    cnt++;
    total_grow = total_curve = 0.0f;
    v_index = _active_list->itemAtCurrent();
    x = v_index.a();y = v_index.b();z = v_index.c();
    thisvoxel=_data->pop(_index(x,y,z));
    p = VISVector((float)x, (float)y, (float)z);
    
    if (get_phi_neighborhood(x, y, z, neighborhood)) {
      //cout<<get_phi_neighborhood(x, y, z, neighborhood);
      cout<<"Voxel evolved out of range ("<<x<<","<<y<<","<<z<<"): "<<_INSIDE_VALUE<<endl;
      dump(x,y,z,1);
      exit(-1);    }

    if (fast_calc_curve) get_derivs_1sided(neighborhood, derivs); 
    else get_derivs (neighborhood, derivs); 
    //get_derivs (neighborhood, derivs); 
    
    // do position correction 
    d_F = VISVector(VISMaxAbs(derivs[DPX], derivs[DMX]), 
		    VISMaxAbs(derivs[DPY], derivs[DMY]),
		    VISMaxAbs(derivs[DPZ], derivs[DMZ]));
    len = sqrt(d_F(0)*d_F(0) + d_F(1)*d_F(1)+ d_F(2)*d_F(2));
    if (len > 0.0) {
      d_F = d_F/len;
      pc = p - (d_F)*(thisvoxel->_value/(len));
    } else pc=p;
    
    if (calc_curve) {
      if (fast_calc_curve) {
	//if (FALSE) {
	curve_trace=calculate_curvature(x,y,z);
	if (_curve_weight>0.0) total_curve += (_curve_weight*curve_trace);
	if (_curvature_grow_weight>0.0) {	 
	  //total_curve += (_curvature_grow_weight*curvature_grow(curve_trace,x,y,z));
	  tmp=curvature_grow(curve_trace,x,y,z);
	  total_curve += (_curvature_grow_weight*tmp);
	  rms_curve+=tmp*tmp;
	}
      } else {
	total_curve += (_special_curve_weight*
			specialCurve(calculate_curvature(derivs, curve_trace, curve_norm),p,d_F));
	if (_curvature_grow_weight>0.0) {
	  //total_curve +=_curvature_grow_weight*curvature_grow(curve_trace,x,y,z);
	  curve_trace=0.2f*curve_trace+0.8f*calculate_curvature(x,y,z);
	  tmp=curvature_grow(curve_trace,x,y,z);
	  total_curve +=_curvature_grow_weight*tmp;
	  rms_curve+=tmp*tmp;
	}
	if (curve_norm > 0.0f) {
	  curve_gauss = (curve_trace*curve_trace - curve_norm)/2.0f;
	  total_curve += (_curve_weight*curve_trace + _gauss_curve_weight*curve_gauss);
	  if (_curve_weight>0.0f) rms_curve+=(curve_trace*curve_trace);
	  
	  //  L_1 norm weighting
	  if ((_weighted_curve_weight != 0.0)&&(curve_gauss > 0.0f))
	    total_curve += (2.0f*_weighted_curve_weight*curve_gauss/curve_trace);
	  if ((_pos_curve_weight != 0.0f)||(_neg_curve_weight != 0.0f)) {
	    k1 = (curve_trace + sqrt(-curve_trace*curve_trace + 2.0f*curve_norm))/2.0f;
	    k2 = (curve_trace - sqrt(-curve_trace*curve_trace + 2.0f*curve_norm))/2.0f;
	    total_curve += _pos_curve_weight*(VISmax(k1, 0.0f) + VISmax(k2, 0.0f))
	      + _neg_curve_weight*(VISmin(k1, 0.0f) + VISmin(k2, 0.0f));
	  }
	}
      }
    }
    //cout<<"t. curve: "<<total_curve;
    total_grow += total_curve;  // 1st addition to total_grow
    
    if (calc_force) {
      if (_energy_weight != 0.0f)
	the_force = _energy_weight*force(pc(0), pc(1), pc(2));
      else the_force = VISVector(0.0f, 0.0f, 0.0f);
      
      if (_normal_energy_weight != 0.0f)
	the_force += (_normal_energy_weight*normal_force(pc(0), pc(1), pc(2), 
							 d_F(0), d_F(1), d_F(2)));
      
      if ((the_force(0)) < 0.0) grad_x = derivs[DPX];
      else grad_x = derivs[DMX];
      
      if ((the_force(1)) < 0.0) grad_y = derivs[DPY];
      else grad_y = derivs[DMY];
      
      if ((the_force(2)) < 0.0) grad_z = derivs[DPZ];
      else grad_z = derivs[DMZ];
      
      force_change = -1.0f*(grad_x*the_force(0) + grad_y*the_force(1) + grad_z*the_force(2));
      total_force_abs = fabs(the_force(0)) + fabs(the_force(1)) + fabs(the_force(2));
      max_total_force = VISmax(max_total_force, total_force_abs);
    } else force_change = 0.0;
    
    // now compute the term for constant vection (and the weighted curv. terms)
    // these terms are weighted averages that ignore direction in order
    // to normalize the curvature terms
    
    if (calc_grow||calc_curve) {
      // constant grow
      total_grow += _const_grow_weight; 
      // general grow
      if (calc_var_grow) {
	//total_grow +=  _grow_weight * grow(pc(0), pc(1), pc(2),d_F(0), d_F(1), d_F(2)); 
	tmp=grow(pc(0), pc(1), pc(2),d_F(0), d_F(1), d_F(2)); 
	total_grow+=_grow_weight*tmp;
	rms_grow+=tmp*tmp;
	//cout<<" grow : "<<tmp<<" "<<total_grow<<endl;
      }
      grad = 0.0f;
      if (total_grow > 0.0f)
	for (ii = 0; ii < 3; ii++) {
	  grad_p = VISmax(derivs[DPX + ii], 0.0f);
	  grad_m = VISmin(derivs[DMX + ii], 0.0f);	      
	  grad += grad_p*grad_p + grad_m*grad_m;
	} else if (total_grow < 0.0f) {
	  for (ii = 0; ii < 3; ii++) {
	    grad_p = VISmin(derivs[DPX + ii], 0.0f);
	    grad_m = VISmax(derivs[DMX + ii], 0.0f);	      
	    grad += grad_p*grad_p + grad_m*grad_m;
	  }
	}
      if (grad>6.0) {
	//cout<<grad<<" grad too big "<<x<<" "<<y<<" " <<z<<endl;
	//dump(x,y,z); 
	grad=6.0f;
	//exit(-1);
      }
      grad = sqrt(grad);
      max_total_grow = VISmax(max_total_grow, (float)fabs(total_grow)); // max_total_grow
      grow_change = grad*total_grow;
    } else grow_change = 0.0f;
    
    the_change = force_change + grow_change;
    
    (_active_list->atCurrent()).value(the_change);
    _active_list->stepForward();
    
    total_change += (the_change*the_change);
  }
  
  max_total_force += max_total_grow;

  if (max_total_force != 0.0f) dt = _wave_dt/max_total_force; else dt = 0.0;

  float curve_weight_tmp = _curve_weight + _curvature_grow_weight + _weighted_curve_weight + 
    _neg_curve_weight + _pos_curve_weight + _gauss_curve_weight + _special_curve_weight; 
  
  if (curve_weight_tmp > 0.0f) dt = VISmin(dt, _curve_dt/curve_weight_tmp);

  cout<<" calc. change info: rms_curve = "<<sqrt(rms_curve/(float)cnt)
      <<"  rms_grow = "<<sqrt(rms_grow/(float)cnt)<<endl;
  
  //cout<<"Total change: "<<total_change<<endl;
  //cout<<"Max. total f: "<<max_total_force<<endl;
  //cout<<"dt          : "<<dt<<endl;
  return (dt);
  //return (_curve_dt,dt);
}

void get_neighborhood (int x, int y, int z, float* data, const VISVolume<float> &values) {
  int i, j, k, l = 0;
  
  for (k = -1; k <= 1; k++) for (j = -1; j <= 1; j++) for (i = -1; i <= 1; i++) 
    data[l++] = values.peek(VISmax(VISmin(x + i, (int)values.width() - 1), 0), 
			    VISmax(VISmin(y + j, (int)values.height() - 1), 0),
			    VISmax(VISmin(z + k, (int)values.depth() - 1), 0));
} 

void get_neighborhood (int x, int y, int z, float* data, VISVolume<float> *values) {
  int i, j, k, l = 0;
  
  for (k = -1; k <= 1; k++) for (j = -1; j <= 1; j++) for (i = -1; i <= 1; i++) 
    data[l++] = values->peek(VISmax(VISmin(x + i, (int)values->width() - 1), 0), 
			     VISmax(VISmin(y + j, (int)values->height() - 1), 0),
			     VISmax(VISmin(z + k, (int)values->depth() - 1), 0));
} 

void get_derivs (const float* n, float* data) {
  data[DX] = 0.5f*(n[14] - n[12]);
  data[DY] = 0.5f*(n[16] - n[10]);
  data[DZ] = 0.5f*(n[22] - n[4]);
  data[DPX] = n[14] - n[13];
  data[DPY] = n[16] - n[13];
  data[DPZ] = n[22] - n[13];  
  data[DMX] = n[13] - n[12];
  data[DMY] = n[13] - n[10];
  data[DMZ] = n[13] - n[4];
  data[DYPX] = 0.5f*(n[17] - n[11]);
  data[DZPX] = 0.5f*(n[23] - n[5]);
  data[DYMX] = 0.5f*(n[15] - n[9]);
  data[DZMX] = 0.5f*(n[21] - n[3]);
  data[DXPY] = 0.5f*(n[17] - n[15]);
  data[DZPY] = 0.5f*(n[25] - n[7]);
  data[DXMY] = 0.5f*(n[11] - n[9]);
  data[DZMY] = 0.5f*(n[19] - n[1]);
  data[DXPZ] = 0.5f*(n[23] - n[21]);
  data[DYPZ] = 0.5f*(n[25] - n[19]);
  data[DXMZ] = 0.5f*(n[5] - n[3]);
  data[DYMZ] = 0.5f*(n[7] - n[1]);
}

void get_derivs_1sided (const float* n, float* data) {
  data[DPX] = n[14] - n[13];
  data[DPY] = n[16] - n[13];
  data[DPZ] = n[22] - n[13];  
  data[DMX] = n[13] - n[12];
  data[DMY] = n[13] - n[10];
  data[DMZ] = n[13] - n[4];
}

float distance_symmetric (VISVolume<float> &A, VISVolume<float> &B, float *d1, float *d2) {
  int w=A.width(), h=A.height(), d=A.depth();
  int x, y, z, xn, yn, zn, i, cnt1, cnt2;
  float xp, yp, zp, val, dist1, dist2;
  int N4c[3][3]={{-1,0,0},{0,-1,0},{0,0,-1}};
  
  if ((w!=B.width())||(h!=B.height())||(d!=B.depth())) {
    cout<<"Error (distance_symmetric) : volume sizes must match\n";
    exit(-1);
  }

  cnt1=cnt2=0;
  dist1=dist2=0.0f;
  for (x=1;x<w;x++) for (y=1;y<h;y++) for (z=1;z<d;z++) {
    if (A(x,y,z)==0.0f) {
      dist1+=B(x,y,z)*B(x,y,z);
      cnt1++;
    } else {
      for (i=0;i<3;i++) {
	xn=x+N4c[i][0];
	yn=y+N4c[i][1];
	zn=z+N4c[i][2];
	if ((A(x,y,z)*A(xn,yn,zn))<0.0f) {
	  xp=((float)xn)-A(xn,yn,zn)*((float)(x-xn))/(A(x,y,z)-A(xn,yn,zn));
	  yp=((float)yn)-A(xn,yn,zn)*((float)(y-yn))/(A(x,y,z)-A(xn,yn,zn));
	  zp=((float)zn)-A(xn,yn,zn)*((float)(z-zn))/(A(x,y,z)-A(xn,yn,zn));
	  val=B.interp(xp,yp,zp);
	  cnt1++;
	  dist1+=val*val;
	}
      }
    }
    if (B(x,y,z)==0.0f) {
      dist2+=A(x,y,z)*A(x,y,z);
      cnt2++;
    } else  {
      for (i=0;i<3;i++) {
	if ((B(x,y,z)*B(xn,yn,zn))<0.0f) {
	  xp=((float)xn)-B(xn,yn,zn)*((float)(x-xn))/(B(x,y,z)-B(xn,yn,zn));
	  yp=((float)yn)-B(xn,yn,zn)*((float)(y-yn))/(B(x,y,z)-B(xn,yn,zn));
	  zp=((float)zn)-B(xn,yn,zn)*((float)(z-zn))/(B(x,y,z)-B(xn,yn,zn));
	  val=A.interp(xp,yp,zp);
	  cnt2++;
	  dist2+=val*val;
	}
      }
    }
  }
  //cout<<cnt1<<" pts. in A , dist: "<<sqrt(dist1/((float)cnt1))<<" , ";
  //cout<<cnt2<<" pts. in B , dist: "<<sqrt(dist2/((float)cnt2));
  val=sqrt((dist1+dist2)/((float)(cnt1+cnt2)));
  //cout<<" , dist : "<<val<<endl;
  (*d1)=sqrt(dist1/((float)cnt1));
  (*d2)=sqrt(dist2/((float)cnt2));

  return val;
}
