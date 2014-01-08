// ******************************************************************
// VISPack. Copyright (c) 1994-2000 Ross Whitaker rtw@utk.edu       *
// For conditions of distribution and use, see the file LICENSE.txt *
// accompanying this distribution.                                  *
// ******************************************************************

// $Id: volume.h,v 1.2 2005/12/01 03:32:51 miriah Exp $


// File:           vol.h
// Author:         Ross T. Whitaker
// Institution:    The University of Tennessee, Knoxville
// Contents:       Several classes all supporting a volume library.  The 'main'
//                 class normally used is the templated VISVolume class.
// Log of changes: July 14, 1999 -- Added comments


#ifdef DEBUG
#define DYNAMIC_CHECKING
#endif


#ifndef vis_volume_h
#define vis_volume_h

#ifdef DEBUG
class VISRep;
extern VISRep *__rep_tracker;
#endif


// *******************************************************************
// *******************************************************************

#include "vol.h"
#include "image/image.h"
#include "util/array.h"

template<class T>
class VISVolume;


template< class T=float >
class VISVolumeRep: public VISImRep<T>
{
  // this class holds the 'data' of a volume.  this class is used to implement
  // the copy-on-write strategy discussed in the document for this library
  friend class VISVolume<T>;
  //    template< class T1, class T2 >    
  //   friend void assignVolume(const VISVolume<T1>& from, VISVolume<T2>& to);
private:
  unsigned int _width, _height, _depth;
  // this is for being able to hand images back.
  // this image shares the buffers
  void initialize(unsigned int w, unsigned int h, unsigned int d,
		  T* buffer_in);
  void initialize(VISVolumeRep<T>*);
protected:
  unsigned int index(int x, int y, int z)
  {return(z*_height*_width + y*_width + x);}
public:
  ~VISVolumeRep(){}    
  VISVolumeRep(unsigned int w, unsigned int h, unsigned int d){
    initialize(w, h, d, NULL);}
  VISVolumeRep(unsigned int w, unsigned int h, unsigned int d, T* buf){
    initialize(w, h, d, buf);}
  VISVolumeRep(){
    initialize(0, 0, 0, NULL);}
  VISVolumeRep(VISVolumeRep<T>* r){
    initialize(r->_width, r->_height, r->_depth, NULL);
    VISImRep<T>::copy((VISImRep<T>*)r);}
  VISVolumeRep<T>* createToSize() const{
    VISVolumeRep<T>* new_rep 
      = new VISVolumeRep<T>(_width, _height, _depth);
    return(new_rep);}
  void insetRep(const VISVolumeRep<T> *in, 
		unsigned int x_pos, 
		unsigned int y_pos, 
		unsigned int z_pos);
  void getROI(const VISVolumeRep<T> *in, 
	      unsigned int x_pos, 
	      unsigned int y_pos, 
	      unsigned int z_pos);
  unsigned int width() const {return(_width);}
  unsigned int height() const {return(_height);}
  unsigned int depth() const {return(_depth);}
  int checkBounds(unsigned int x, unsigned int y, unsigned int z) const
  {return((x<_width)&&(y<_height)&&(z<_depth));}
  int checkBounds(int x, int y, int z) const
  {return((x<_width)&&(y<_height)&&(z<_depth)
	  &&(x >= 0)&&(y >= 0)&&(z >= 0));}
  int checkBounds(float x, float y, float z) const{
    return((x >= 0.0f)&&(x<=(float)(_width - 1))&&(y >= 0.0f)
	   &&(y<=(float)(_height - 1))
	   &&(z >= 0.0f)&&(z<=(float)(_depth - 1)));}
  T& at(unsigned int x, unsigned int y, unsigned int z)
  {
#ifdef DYNAMIC_CHECKING
    if (!checkBounds(x, y, z))
    {
      ERROR("VISVolume<>: channel(int ch) - channel out of range");
    }
#endif 
    return(VISVolumeRep<T>::_buffer[z*_height*_width + y*_width + x]);}// **mdm**
  const T& itemAt(unsigned int x, unsigned int y, unsigned int z) const
  {
#ifdef DYNAMIC_CHECKING
    if (!checkBounds(x, y, z))
    {
      ERROR("VISVolume<>: itemAt(int, int, int) - channel out of range");
    }
#endif 
    return(VISVolumeRep<T>::_buffer[z*_height*_width + y*_width + x]);}// **mdm**
  int compareSize(const VISVolumeRep<T>* other) const{
    return((_width==other->_width)
	   &&(_height == other->_height)
	   &&(_depth == other->_depth)
	   &&(VISImRep<T>::compareSize(other)));}
  void evaluate(T (*f)(unsigned int, unsigned int, unsigned int));
  void evaluate(T (*f)(unsigned int, unsigned int, unsigned int), 
		unsigned int, unsigned int, 
		unsigned int, unsigned int, 
		unsigned int, unsigned int);
  void evaluate(T (*f)(const T&), const VISVolumeRep<T>* other);

};


template<class T=float>
class VISVolume: public VISVol
{
  // IRISVOLUME
  // This class is the 'main' volume class, meaning this is the
  // class users will normaly use as an volume.  Most of the 
  // methods usually used by users are in this class, but some
  // methods exist up the inheritance tree.  (Notice that
  // VISVolume is inherited from VISVol.)


  //  private:
  //    VISVolumeRep<T>* _rep;
    
protected:
  void initialize(unsigned int, unsigned int, unsigned int);
  void initialize(const VISVolumeRep<T>* rep);
  void copy(const VISVolume<T>& a);
  void copy_on_write(VISRep*& r);


  //*****************************//
  //  FUNCTIONS USED BY THE USER //
  //*****************************//
public:
  //****Constructors and Destructors, the usuall stuff
  VISVolume(){initialize(0, 0, 0);}//create a volume of size zero
  VISVolume(unsigned int w, unsigned int h, unsigned int d){
    initialize(w, h, d);}//create a volume of size width, height, depth
  VISVolume(const VISVol& vol);//create a volume from an VISVol
  //create a volume from a multi-channel image.
  VISVolume(const VISImage<T>& image);
  //create a volume from another VISVolume
  VISVolume(const VISVolume<T>& volume);
  //create an image from a pointer to a buffer of type T of size
  //width x height x depth
  VISVolume(unsigned int w, unsigned int h, unsigned int d, T* buf);
  //create a volume from an VISVolumeRep
  VISVolume(VISVolumeRep<T>* rep);
  //destructor
  ~VISVolume();   
  //create an image to match size of another image
  VISVolume<T> createToSize() const{ 
    VISVolume<T> new_vol(width(), height(), depth()); return(new_vol);}
    
  //****Finding the size of the image
  // A volume has three indices width,height,depth
  unsigned int width() const {return(rep()->width());}
  unsigned int height() const {return(rep()->height());}
  unsigned int depth() const {return(rep()->depth());}
  boolean isValid() const{
    return((VISVol::isValid())&&(height() > 0)&&(width() > 0)
	   &&(depth() > 0));}
  void print() const;//prints out the height, width, and depth of a volume
  void print(const char* str) const{//prints out some string with the height
    printf("%s\n", str); print();}//width, depth, and type of the volume

  //****Checking bounds
  //these functions return 0 if 'in bounds' or 1 if 'out bounds'
  //x is width, y is height, and z is depth to be checked.
  //x, y, and z can be floats
  //compareSize compares size of one volume to another (other)
  int checkBounds(unsigned int x, unsigned int y, unsigned int z) const{
    return(rep()->checkBounds(x, y, z));}
  int checkBounds(int x, int y, int z) const{
    return(rep()->checkBounds(x, y, z));}
  int checkBounds(float x, float y, float z) const{
    return(rep()->checkBounds(x, y, z));}
  int compareSize(const VISVolume<T>& other) const{
    return(rep()->compareSize(other.rep()));}
  // returns a const pointer to the begining of the data of a volume
  const VISVolumeRep<T>* rep() const{
    return((const VISVolumeRep<T>*)_rep);}
  // returns a non-const pointer to the data of a volume
  // note: this routine copies the data (copy-on-write)
  VISVolumeRep<T>* repRef(){
    copy_on_write(_rep); return((VISVolumeRep<T>*)_rep);}

  //****Interacting with the data of the image
  // This group of functions includes:
  //    general functions for 'getting' and 'putting' values.
  //    functions that return pointers to the data of an volume
  //    functions that interpolate to return values 'between' voxels
  //    functions that 'get' and 'put' regions of interest
  //    functions that apply a user created function to each voxel of an image
  //
  // General functions for 'getting' voxel values from a volume
  // and 'putting' voxel values into a volume 
  // Make sure you use 'at' on the lhs and 'itemAt' on the rhs
  // NOTE: the pixels are refered to in a col,row manner
  // examples:    VISVolume m(5,5,5),n(5,5,5);
  //              m.at(2,2,2) = 5.0f; //x=2,y=2,z=2 is 5.0f
  //              m.at(2,1,2) = m.itemAt(2,2,2); //now x=2,y=1,z=2 is 5.0f
  T& at(unsigned int x, unsigned int y, unsigned int z){
    return(repRef()->at(x, y, z));}
  T& poke(unsigned int x, unsigned int y, unsigned int z){
    return(repRef()->at(x, y, z));}
    
  const T& itemAt(unsigned int x, unsigned int y, unsigned int z) const{
    return(rep()->itemAt(x, y, z));}
  const T& peek(unsigned int x, unsigned int y, unsigned int z) const{
    return(rep()->itemAt(x, y, z));}
  const T& operator()(unsigned int x, unsigned int y, unsigned int z) const{
    return(rep()->itemAt(x, y, z));}

  // interpolate between voxels to determine the value returned.
  // this is a linear interpolation
  // interpNoBounds does not use bounds checking.  These functions
  // are therefore faster than their counterparts, BUT you can
  // cause a core dump if you use these with numbers that are out of bounds
  // examples:    VISImage m(5,5),n(5,5);
  //              m.at(0,0) = n.interp(1.45,2.45);
  T interp(float x, float y, float z) const;
  // these functions allow the user to 'get' and 'put' a region of interest
  // this region of interest has a position (x_pos,y_pos,z_pos; its upper left
  // corner) and a size (w_roi,h_roi,d_roi).  'putting' the roi involves the
  // roi volume and an upper left corner position x_pos,y_pos,z_pos
  void putROI(const VISVolume<T>& volume_in, 
	      unsigned int x_pos, unsigned int y_pos, unsigned int z_pos);
  VISVolume<T> getROI(unsigned int x_pos, unsigned int y_pos, 
		      unsigned int z_pos, unsigned int w_roi, 
		      unsigned int h_roi, unsigned int d_roi) const;
  // evalutate evaluates a function on a volume
  // the user writes a function that has as input three unsigned int's:
  // x, y, and z -- the output is the type of the volume
  // The idea is that this function is writen to work on one voxel
  // of the volume and evaluate applies this function to each voxel
  // in the image.
  void evaluate(T (*f)(unsigned int, unsigned int, unsigned int), 
		unsigned int rect_ul_x, unsigned int rect_ul_y, 
		unsigned int rect_ul_z, unsigned int rect_width, 
		unsigned int rect_height, unsigned int rect_depth){
    repRef()->evaluate(f, rect_ul_x, rect_ul_y, rect_ul_z, 
		       rect_width, rect_height, rect_depth);}	
  void evaluate(T (*f)(unsigned int, unsigned int, unsigned int)){
    repRef()->evaluate(f);}
  VISVolume<T> evaluate(T (*f)(const T&)) const 
  {VISVolume<T> r = this->createToSize();
    r.repRef()->evaluate(f, this->repRef());
    return(r);
  } // **mdm**
    // this does a floating point printf on all of the data.  Be careful!
    //void VISVolume<T>::printData() const;
  void printData() const;

  //****Operators that have been overloaded (other overloaded operators
  // occur outside the class, search for 'NMOO' meaning non-member overloaded
  // operators)
  // every voxel in an image is assigned value with this operator=
  VISVolume<T>& operator=(T value); 
  // one volume equals another with this operator=
  VISVolume<T>& operator=(const VISVolume<T>& from){
    assign(from); return(*this);}
  // ask
#ifdef this_is_dangerous // because of automatic casting done by C++
  VISVolume<T>& operator=(const VISVol& from); 
#endif
  // add, subtract, multiply, divide each voxel in
  // an volume by the respective voxel in another volume
  VISVolume<T> operator+(const VISVolume<T>& volume) const;
  VISVolume<T> operator*(const VISVolume<T>& volume) const;
  VISVolume<T> operator-(const VISVolume<T>& volume) const;
  VISVolume<T> operator/(const VISVolume<T>& volume) const;
  // add, subtract, multiply, divide each voxel in
  // an volume by the respective voxel in another volume
  // and assign the result to the volume being operated on
  VISVolume<T>& operator*=(const VISVolume<T>& volume); 
  VISVolume<T>& operator+=(const VISVolume<T>& volume);
  VISVolume<T>& operator-=(const VISVolume<T>& volume);
  VISVolume<T>& operator/=(const VISVolume<T>& volume); 
  // add, subtract, multiply, divide each volume by
  // a value and assign the result to the volume being
  // operated on
  VISVolume<T>& operator*=(const T& value);
  VISVolume<T>& operator+=(const T& value);
  VISVolume<T>& operator-=(const T& value);
  VISVolume<T>& operator/=(const T& value);

  //****Math functions
  VISVolume<T> power(const int& power);// raise the volume to a power
  // operator^ also works
  VISVolume<T> sqrt() const;// returns square root of the volume
  //divide by volume, if division by zero is encountered return zeroCondition
  //ask
  VISVolume<T> div(const VISVolume<T>& volume, T value) const;
  VISVolume<T> abs() const;// absolute value of the volume
  T max() const;// returns max value in volume
  T min() const;// returns min value in volume
  // below min and max functions compare each pixel of two volumes and
  // return an volume containing the min or max values for each pixel
  VISVolume<T> min(const VISVolume<T> &other) const;
  VISVolume<T> max(const VISVolume<T> &other) const;
  VISVolume<T> min(const T &value) const;
  VISVolume<T> max(const T &value) const;
  float sum() const;// returns the sum of all pixels in the volume
  float average();// returns the average of all pixels in the volume
  VISVolume<T> scale(float value) const;// multiplies each pixel by value

  //****General Volume Processing functions
  //filtering functions
  VISVolume<T> gauss(float sigma) const;
  // convolution with any filter you choose
  // NOTE: filters can be created with a number of functions
  // search for 'Filter Creation' in this file
  VISVolume<T> convolve(const VISVolume<T>& kernel) const;
  //determines the inner product of a mask and an area in an volume.
  //the center of the mask is indicated by center_x, center_y, center_z
  //and the center of the mask is positioned at x_pos, y_pos, z_pos of the
  //volume and then the inner product is calculated.
  float maskFloat(const VISVolume<float>& mask, unsigned int center_x, 
		  unsigned int center_y, unsigned int center_z, 
		  unsigned int x_pos, unsigned int y_pos,
		  unsigned int z_pos) const;
  // does the same as maskFloat except an array of masks is applied
  // to the volume, and thus an array of floats is returned
  //VISArray<float>* VISVolume<T>::maskFloat(const VISArray<VISVolume<float> >& masks, unsigned int x_pos, unsigned int y_pos,unsigned int z_pos) const;
  VISArray<float>* maskFloat(const VISArray<VISVolume<float> >& masks, unsigned int x_pos, unsigned int y_pos,unsigned int z_pos) const;

  // taking derivatives
  VISVolume<T> dx() const;//derivative in x direction (along width)
  VISVolume<T> dy() const;//derivative in y direction (along height)
  VISVolume<T> dz() const;//derivative in z direction (along depth)
  VISVolume<T> dx(unsigned int order) const;// n-th order derivative
  VISVolume<T> dy(unsigned int order) const;// n-th order derivative
  VISVolume<T> dz(unsigned int order) const;// n-th order derivative
  // taking the n-th order derivative in x, y, and z directions
  // derivatives can be of different order in the x vs. y or z directions
  VISVolume<T> derivative(unsigned int order_x, unsigned int order_y, 
			  unsigned int order_z) const;
  // these half derivatives are useful in diffusion processes
  // need to be zeroed out along the border, by Zhong on 5/31/100
  VISVolume<T> dxHalfForward() const;
  VISVolume<T> dyHalfForward() const;
  VISVolume<T> dzHalfForward() const;
  VISVolume<T> dxHalfBack() const;
  VISVolume<T> dyHalfBack() const;
  VISVolume<T> dzHalfBack() const;
  // reduces the size of a volume by subsampling, taking ever n'th voxel
  // where n = scale an integer
  VISVolume<T> reduce(int scale) const;
  // resample scales the size of the volume by interpolating between
  // voxels.  the volume size can be changed by a scale.
  VISVolume<T> resample(float the_scale) const;
  // here the volume size is changed to w x h x d
  VISVolume<T> resample(unsigned w, unsigned h, unsigned d) const;
  // here the volume size is scaled by scale_x in width and scaled by
  // scale_y in height and scale_z in depth.  x, y, and z indicate a
  // position in the volume upon which the center of the new volume
  // will be placed
  VISVolume<T> resample(float scale_x, float scale_y, float scale_z,
			float x, float y, float z) const;
  //ask
  VISVolume<float> noise() const;
  // sets a border of voxels of "thickness" = w to "value" around the border
  void setBorder(T value, unsigned w);
  // sets a border of voxels of "thickness" = 1 to "value" around the border
  void setBorder(T value) { setBorder(value, 1);}


  // returns an image of the same size with values between 0.0 and 1.0 
  VISVolume<float> cityBlockDistTrans() const;

  VISVolume<T> zeroCrossings() const;
  VISVolume<T> floodFill(T thresh, unsigned x, unsigned y, unsigned z);

  // these return images with one channel per slice
  // right now can't work out difference between the "image" and 
  // "imageRef" cases --- Ross 5-28-97
  const VISImage<T> image() const;
  VISImage<T> imageRef();
  VISImage<T> image(unsigned int slice);
  VISImage<T> imageRef(unsigned int slice);


  //**** Support Functions (NOT USUALLY USED DIRECTLY BY THE USER)
  VISVolume<T> mult(T value) const;
  VISVolume<T> add(T value) const;
  VISVolume<T> div(T value) const;
  VISVolume<T> div_by(T value) const;
  VISVolume<T> sub_from(T value) const;
  VISVolume<T> sub(T value) const;
  VISVolume<T> multAssign(const VISVolume<T>& volume);
  VISVolume<T> addAssign(const VISVolume<T>& volume);
  VISVolume<T> divAssign(const VISVolume<T>& volume);
  VISVolume<T> subAssign(const VISVolume<T>& volume);
  VISVolume<T> multAssign(const T& value);
  VISVolume<T> divAssign(const T& value);
  VISVolume<T> addAssign(const T& value);
  VISVolume<T> subAssign(const T& value);
  void assign(const VISVol& from);		
  void assign(const VISVolume<T>& from);
  VISVolume<T> gt(T value) const;
  VISVolume<T> gteq(T value) const;
  VISVolume<T> lt(T value) const;
  VISVolume<T> lteq(T value) const;
  VISVolume<T> eq(T value) const;

  // used to assign volumes to each other
  // public but only to be used by friendly 
  // third parties.  (problems with templates&friends)
  void putRep(const VISVolumeRep<T>* r);

};  

//#define dIndex(a, b) (((((a) + (b) + 1)*((a) + (b)))/2) + (b))

unsigned dVolIndex(unsigned c, unsigned b, unsigned a);
/*{
    return((((a + b + c)*(2*(a + b + c)*(a + b + c) + 3*(a + b + c) + 1))/6
	    + (a + b + c)*(a + b + c + 1)/2)/2
	   + (a + b)*(a + b + 1)/2 + a);}
*/
//**** NMOO -- Non-Member Overloaded Operators
template< class T >
inline VISVolume<T> operator+(const VISVolume<T>& from, T value) {
    return(from.add(value));}
template< class T >
inline VISVolume<T> operator+(T value, const VISVolume<T>& from) {
    return(from.add(value));}
template< class T >
inline VISVolume<T> operator-(const VISVolume<T>& from, T value) {
    return(from.sub(value));}
template< class T >
inline VISVolume<T> operator-(T value, const VISVolume<T>& from) {
    return(from.sub_from(value));}
template< class T >
inline VISVolume<T> operator/(const VISVolume<T>& from, T value) {
    return(from.div(value));}
template< class T >
inline VISVolume<T> operator/(T value, const VISVolume<T>& from) {
    return(from.div_by(value));}
template< class T >
inline VISVolume<T> operator>(T value, const VISVolume<T>& from) {
    return(from.lt(value));}
template< class T >
inline VISVolume<T> operator>=(T value, const VISVolume<T>& from){
    return(from.lteq(value));}
template< class T >
inline VISVolume<T> operator<(T value, const VISVolume<T>& from) {
    return(from.gt(value));}
template< class T >
inline VISVolume<T> operator<=(T value, const VISVolume<T>& from){
    return(from.gteq(value));}
template< class T >
inline VISVolume<T> operator>(const VISVolume<T>& from, T value) {
    return(from.gt(value));}
template< class T >
inline VISVolume<T> operator>=(const VISVolume<T>& from, T value){
    return(from.gteq(value));}
template< class T >
inline VISVolume<T> operator<(const VISVolume<T>& from, T value) {
    return(from.lt(value));}
template< class T >
inline VISVolume<T> operator<=(const VISVolume<T>& from, T value){
    return(from.lteq(value));}
template< class T >
inline VISVolume<T> operator==(T value, const VISVolume<T>& from){
    return(from.eq(value));}
template< class T >
inline VISVolume<T> operator==(const VISVolume<T>& from, T value){
    return(from.eq(value));}

template< class T1, class T2 >
void copy(const VISVolume<T1>& from, VISVolume<T2>& to);
template< class T1, class T2 >
void assignVolume(const VISVolume<T1>& from, VISVolume<T2>& to);
template< class T1, class T2 >
int compareSize(const VISVolume<T1>& a, const VISVolume<T2>& b);
template< class T1, class T2 >
void copy(const VISVolumeRep<T1>* from, VISVolumeRep<T2>* to);
template< class T1, class T2 >
int compareSize(const VISVolumeRep<T1>* a, const VISVolumeRep<T2>* b);
template< class T >
VISVolume<T> operator*(const VISVolume<T>& from, T value);
template< class T >
VISVolume<T> operator*(T value, const VISVolume<T>& from);

VISArray< VISVolume<float> >* derivativeVolMasks(unsigned degree);
VISVolume<float> dz_kernel(int order);
VISVolume<float> gauss_slice_kernel(float sigma);
VISVolume<float> gauss_dz_kernel(int order, float sigma);


#ifndef MANUAL_INSTANTIATION
#include "volume.txx"
#endif

#endif





