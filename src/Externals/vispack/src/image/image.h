// ******************************************************************
// VISPack. Copyright (c) 1994-2000 Ross Whitaker rtw@utk.edu       *
// For conditions of distribution and use, see the file LICENSE.txt *
// accompanying this distribution.                                  *
// ******************************************************************

// $Id: image.h,v 1.2 2005/12/01 03:32:51 miriah Exp $

// File:           image.h
// Author:         Ross T. Whitaker
// Institution:    The University of Tennessee, Knoxville
// Contents:       Several classes all supporing an image library.  The 'main'
//                 class normally used is the templated VISImage class.
// Log of changes: June 28, 1999 -- Added comments

// Table of Contents for Comments (search strings):
//     IRISIMREP
//     IRISIMAGEREP
//     IRISIMAGE
//         FUNCTIONS USED BY THE USER
//             Constructors and Destructors, the usuall stuff
//             Finding the size of the image
//             Interacting with the data of the image
//             Checking bounds
//             Operators that have been overloaded
//             Math functions
//             General Image Processing functions
//             Specific Task Image Processing Functions
//             Support Functions (NOT USUALLY USED DIRECTLY BY THE USER)
//    Filter Creation
//    NMOO -- Non-Member Overloaded Operators

// #ifndef AUTO_IMAGE_CONVERSION
// #define AUTO_IMAGE_CONVERSION
// #endif


#ifdef DEBUG
#define DYNAMIC_CHECKING
#endif

#ifndef vis_image_h
#define vis_image_h

#define __memory_count 
#ifdef __memory_count
extern unsigned __MEMORY_SIZE;
#define VISImageMemory() printf("memory total is %d\n", __MEMORY_SIZE);
#endif


// this determines whether or not you expect to include the 
//hips library 
// #define USE_HIPS 

#ifdef DEBUG
class VISRep;
extern VISRep *__rep_tracker;
#endif


#define DEFAULT_NUM_CHANNELS (1)
#define DEFAULT_CHANNEL (0)

#include "image/rgba.h"

// *******************************************************************
// *******************************************************************

#include "im.h"

template< class T >
class VISArray;

template< class T >
class VISList;

template< class T=float >
class VISImRep: public VISRep
{
//IRISIMREP
  private:
    int _size;	
  protected:
    T* _buffer;
    boolean _free_buffer;
    void initialize(unsigned int size, T* buffers);
  public:
    ~VISImRep();
    const T* buffer() const {return((const T*)_buffer);}
    T* bufferRef()  {return(_buffer);}
    VISImRep<T>* createToSize() const   {
	VISImRep<T>* new_rep = new VISImRep<T>(_size);
	return(new_rep);}
    int compareSize(const VISImRep<T>* other) const{
	return(_size == other->_size);}
    int checkBounds(unsigned int i) const
	{return(i < _size);}
    VISImRep(unsigned int size){
	initialize(size, (T*)NULL);}
    VISImRep(VISImRep<T>* r){
	initialize(r->_size, (T*)NULL);
	copy(r);}
    VISImRep(){
	initialize(0, (T*)NULL);}
    VISImRep(unsigned int size, T* buffer){
	initialize(size, buffer);}
    void clear(T value);
    void copy(const VISImRep<T>* a);
    void copyBuffer(const T* buf);
    void add(const VISImRep<T>* a, const VISImRep<T>* b);
    void add(const VISImRep<T>* a, const T& value);
    void evaluate(T (*f)(T), const VISImRep<T>* a);
    void mult(const VISImRep<T>* a, const VISImRep<T>* b);  
    void mult(const VISImRep<T>* a, const T& value);
    void scale(const VISImRep<T>* a, float value);
    void div(const VISImRep<T>* a, const VISImRep<T>* b);  
    void div(const T& value, const VISImRep<T>* a);
    void div(const VISImRep<T>* a, const T& value);
// these two return "zeroCondition" when the denom is zero
    void div(const VISImRep<T>* a, const VISImRep<T>* b, T zeroCondition);
    void div(const T& value, const VISImRep<T>* a, T zeroCondition);
    void sub(const VISImRep<T>* a, const VISImRep<T>* b);    
    void sub(const T& value, const VISImRep<T>* a);
    void sub(const VISImRep<T>* a, const T& value);
    void power(const VISImRep<T>* other, int power1);
    void sqrt(const VISImRep<T>* a);
    void ln(const VISImRep<T>* a);
    void exp(const VISImRep<T>* a);
// uses a lookup table with linear interp
//    void expFast(const VISImRep<T>* a);
    void pos(const VISImRep<T>* a);
    void neg(const VISImRep<T>* a);
    void abs(const VISImRep<T>* a);
    void sign(const VISImRep<T>* a);
    void gt(const VISImRep<T>* other, T value);
    void gteq(const VISImRep<T>* other, T value);
    void lt(const VISImRep<T>* other, T value);
    void  lteq(const VISImRep<T>* other, T value);
    void eq(const VISImRep<T>* other, T value);
    T max() const;
    T min() const;
    float sum() const;
    void max(const VISImRep<T>* a, const VISImRep<T>* b);    
    void min(const VISImRep<T>* a, const VISImRep<T>* b);    
    void max(const VISImRep<T>* a, T value);    
    void min(const VISImRep<T>* a, T value);    
    int size() const{
	return(_size);}
    T& at(unsigned int i)
	{
#ifdef DYNAMIC_CHECKING
	    if (!checkBounds(i))
		ERROR("VISImRep: element - index out of bounds");
#endif
            return(_buffer[i]);}
    const T& itemAt(unsigned int i) const
	{
#ifdef DYNAMIC_CHECKING
	    if (!checkBounds(i))
		ERROR("VISImRep: element - index out of bounds");
#endif
		return(_buffer[i]);
	}
    const T& operator[](unsigned int i) const
	{
#ifdef DYNAMIC_CHECKING
	    if (!checkBounds(i))
		ERROR("VISImRep: element - index out of bounds");
#endif
		return(_buffer[i]);
	}

};


template<class T=float>
class VISImage;

template< class T=float >
class VISImageRep: public VISImRep<T>
{
//IRISIMAGEREP
// this class holds the 'data' of an image.  this class is used to implement
// the copy-on-write strategy discussed in the document for this library
    friend class VISImage<T>;
  private:
    unsigned int _width, _height;
    void initialize(unsigned int w, unsigned int h,
			    T* buffer_in);
    int checkBounds(unsigned int x, unsigned int y) const
      {return((x<_width)&&(y<_height));}
    int checkBounds(int x, int y) const
    {return((x<_width)&&(y<_height)&&(x >= 0)&&(y >= 0));}
    int checkBounds(float x, float y) const{
	return((x <= ((float)_width - 1.0f))&&(y <= ((float)_height - 1.0f))
	       &&(x >= 0.0f)&&(y >= 0.0f));}
  protected:
    unsigned int index(int x, int y)
	{return(y*_width + x);}
 public:
    ~VISImageRep(){}
    VISImageRep(unsigned int w, unsigned int h){
	initialize(w, h, NULL);}
    VISImageRep(unsigned int w, unsigned int h, T* buf){
	initialize(w, h, buf);}
    VISImageRep(){
	initialize(0, 0, NULL);}
  /*
    VISImageRep(VISImageRep<T>& r){
	initialize(r._width, r._height, NULL);
	copy(&r);}*/
    VISImageRep(const VISImageRep<T>& r){
	initialize(r._width, r._height, NULL);
	VISImRep<T>::copy(&r);} // **mdm** 
    VISImageRep<T>* createToSize() const   {
	    VISImageRep<T>* new_rep = new VISImageRep<T>(_width, _height);
	    return(new_rep);}
    unsigned int width() const {return(_width);}
    unsigned int height() const {return(_height);}
    T& at(unsigned int x, unsigned int y)
      {
#ifdef DYNAMIC_CHECKING
	  if (!checkBounds(x, y))
	      {
		  ERROR("VISImage<>: channel(int ch) - channel out of range");
	      }
#endif 
	    return(VISImageRep<T>::_buffer[y*_width + x]);// **mdm**
      }

    const T& itemAt(unsigned int x, unsigned int y) const
      {
#ifdef DYNAMIC_CHECKING
	if (!checkBounds(x, y))
	    {
		ERROR("VISImage<>: channel(int ch) - channel out of range");
	    }
#endif 
	    return(VISImageRep<T>::_buffer[y*_width + x]);// **mdm**
      }
    int compareSize(const VISImageRep<T>* other) const{
	    return((_width==other->_width)
		       &&(_height == other->_height)
		   &&(VISImRep<T>::compareSize(other)));}
    void evaluate(T (*f)(unsigned int, unsigned int));
    void evaluate(T (*f)(unsigned int, unsigned int), unsigned int,
		unsigned int, unsigned int, unsigned int);
    void insetRep(const VISImageRep<T> *in, unsigned int x_pos, 
		  unsigned int y_pos);
    void getROI(const VISImageRep<T> *in, unsigned int x_pos, 
		unsigned int y_pos);

};

template<class T>
class VISImage: public VISIm
{
// IRISIMAGE
// This class is the 'main' image class, meaning this is the
// class users will normally use as an image.  Most of the 
// methods usually used by users are in this class, but some
// methods exist up the inheritance tree.  (Notice that
// VISImage is inherited from VISIm.)

    friend class VISImageRGBA;
  //friend void assignImage(const VISIm& from, VISIm& to);

  private:

  protected:
    void initialize(unsigned int, unsigned int, unsigned int);
    void initialize(VISImageRep<T> const * const rep[],  unsigned int ch);
    void copy(const VISImage<T>& a);
    void copy_on_write(VISRep*& r, unsigned ch);
    // take all of the connect regions and put them on one list
    // used in the watershed routine
    void getLists(VISList<int> &these_regions, VISArray< VISList<int> > &connections, 
		  int starting_region);


//*****************************//
//  FUNCTIONS USED BY THE USER //
//*****************************//
  public:

//****Constructors and Destructors, the usuall stuff
    VISImage(){//create an image of size zero
	initialize(0, 0, DEFAULT_NUM_CHANNELS);}
    VISImage(unsigned int w, unsigned int h){//create an image of size w x h
	initialize(w, h, DEFAULT_NUM_CHANNELS);}
    VISImage(unsigned int w, unsigned int h, unsigned int ch){
	initialize(w, h, ch);}//create an image w x h with ch channels
    //create an image w x h with ch channels, using buf as the data
    VISImage(unsigned int w, unsigned int h, unsigned int ch, T** buf);
    //create an image 
    VISImage(VISImageRep<T> const * const rep[], unsigned int ch);
    //create an image that is a copy of another image
    VISImage(const VISImage<T>& image);
    VISImage(const VISIm& im);
    VISImage(const VISImage<T>& a, const VISImage<T>& b);//ask
    ~VISImage();//destructor
    //create an image to match size of another image
    VISImage<T> createToSize() const{ 
	VISImage<T> new_im(width(), height(), _channels); return(new_im);}

//****Finding the size of the image
    //in the base class -- height() and width()
    // put these in the base class - Ross 1-19-95
    //    unsigned int width() const {return(rep()->width());}
    //    unsigned int height() const {return(rep()->height());}
    unsigned int channels() const {return(_channels);}
    // printing out the size
    void print() const;// prints out the height, width, and type of an image
    void print(const char* str){// prints out some string with the height
	printf("%s\n", str); print();}// width and type or the image
    boolean isValid() const{// returns 1 if size is > 0, returns 0 otherwise
	return((VISIm::isValid())&&(height() > 0)&&(width() > 0));}

//****Interacting with the data of the image
    // This group of functions includes:
    //    general functions for 'getting' and 'putting' values
    //    functions to add channels and get data from channels of an image
    //    functions that return pointers to the data of an image
    //    functions that interpolate to return values 'between' pixels
    //    functions that 'get' and 'put' regions of interest
    //    functions that apply a user created function to each pixel of an image
    //
    // General functions for 'getting' pixel values from an image
    // and 'putting' pixel values into an image 
    // Make sure you use 'at' on the lhs and 'itemAt' on the rhs
    // NOTE: the pixels are refered to in a col,row manner
    // examples:    VISImage m(5,5),n(5,5);
    //              m.at(2,2) = 5.0f; //column=2,row=2 is 5.0f
    //              m.at(2,1) = m.itemAt(2,2); //now col=2,row=1 is 5.0f
    //              VISImage p(5,5,3); //a 5x5 image with 3 channels
    //              p.at(2,2,3) = m.itemAt(2,2); //now channel 3 of
    //                                           //p col=2,row=2 is 5.0f
    // NOTE (9-16-99): at and itemAt are old function names
    // use poke and peek instead or use operator() for peek
    const T& itemAt(unsigned int x, unsigned int y) const{
	return(rep()->itemAt(x, y));}
    const T& itemAt(unsigned int x, unsigned int y, unsigned int ch) const{
	return(rep(ch)->itemAt(x, y));}
    T& at(unsigned int x, unsigned int y){
	return(repRef()->at(x, y));}
    T& at(unsigned int x, unsigned int y, unsigned int ch){
	    return(repRef(ch)->at(x, y));}

    const T& peek(unsigned int x, unsigned int y) const{
	return(rep()->itemAt(x, y));}
    const T& peek(unsigned int x, unsigned int y, unsigned int ch) const{
	return(rep(ch)->itemAt(x, y));}
    const T& operator()(unsigned int x, unsigned int y) const{
	return(rep()->itemAt(x, y));}
    const T& operator()(unsigned int x, unsigned int y, unsigned int ch) const{
	return(rep(ch)->itemAt(x, y));}    
    T& poke(unsigned int x, unsigned int y){
	return(repRef()->at(x, y));}
    T& poke(unsigned int x, unsigned int y, unsigned int ch){
	    return(repRef(ch)->at(x, y));}

    // 'putting' and 'getting' channels in an image
    void putChannel(VISImage<T>& other, unsigned int ch);
    void putChannel(VISImage<T>& other);// place an image into a new channel
    VISImage<T> channel(unsigned int ch){// returns the image at channel ch
        #ifdef DYNAMIC_CHECKING
	if (!(ch < channels()))
	    ERROR("VISImage: channel - channel out of bounds");
        #endif
	VISImageRep<T> const* tmp_array[1];
	tmp_array[0] = rep(ch);
	VISImage<T> new_image(tmp_array, 1);
        //VISImage<T> new_image(repVISArray() + ch, 1);
	return(new_image);}
    // the following 'rep' functions return pointers to the data of an image
    // this allows the user to improve efficiency of code if desired.
    // returns a const pointer to the data of a choosen channel of an image
    const VISImageRep<T>* rep(unsigned int ch) const{
    #ifdef DYNAMIC_CHECKING
	if (!(ch < _channels)){
		WARN("VISImage<>: channel(int ch) - channel out of range");
		return(NULL);}
    #endif 
	return((const VISImageRep<T>*)_rep[ch]);}
    // returns a const pointer to the first channel of an image
    const VISImageRep<T>* rep() const{return(rep(DEFAULT_CHANNEL));}
    // returns a non-const pointer to the data of a choosen channel of an image
    // note: this routine copies the data (copy-on-write)
    VISImageRep<T>* repRef(unsigned int ch) {
        #ifdef DYNAMIC_CHECKING
	if (!(ch < _channels)){
		ERROR("VISImage<>: channel(int ch) - channel out of range");
		return(NULL);}
        #endif 
	copy_on_write(_rep[ch], ch);
	return((VISImageRep<T>*)rep(ch));}
    // returns a non-const pointer to the data of the first channel of an image
    // this routine copies the data (copy-on-write)
    VISImageRep<T>* repRef() {return(repRef(DEFAULT_CHANNEL));}
    // interpolate between pixels to determine the value returned
    // this is a linear interpolation
    // interpNoBounds does not use bounds checking.  These functions
    // are therefore faster than their counterparts, BUT you can
    // cause a core dump if you use these with numbers that are out of bounds
    // examples:    VISImage m(5,5),n(5,5);
    //              m.at(0,0) = n.interp(1.45,2.45);
    T interp(float x, float y, unsigned int ch) const;
    T interp(float x, float y) const{
	return(interp(x, y, DEFAULT_CHANNEL));}
    T interpNoBounds(float x, float y, unsigned int ch) const;
    T interpNoBounds(float x, float y) const{
	return(interpNoBounds(x, y, DEFAULT_CHANNEL));}
    // these functions allow the user to 'get' and 'put' a region of interest
    // this region of interest has a position (x_pos,y_pos; its upper left corner)
    // and a size (w_roi,h_roi).  'putting' the roi involves the roi image and
    // a upper left corner position x_pos,y_pos
    VISImage<T> getROI(unsigned int x_pos, unsigned int y_pos, 
		       unsigned int w_roi, unsigned int h_roi) const; 
    void putROI(const VISImage<T>& image_in, unsigned int x_pos, 
		unsigned int y_pos);
    // this does a floating point printf on all of the data.  Be careful!
    void printData() const;
    // evalutate evaluates a function on an image
    // the user writes a function that has as input two unsigned int's
    // width and height the output is the type of the image
    // The idea is that this function is writen to work on one pixel
    // of the image and evaluate applies this function to each pixel
    // in the image.
    const VISImage<T>& evaluate(T (*f)(unsigned int, unsigned int));
    VISImage<T> evaluate(T (*f)(T)) const;
    void evaluate(T (*f)(unsigned int, unsigned int), 
		  unsigned int rect_ul_x, unsigned int rect_ul_y, 
		  unsigned int rect_width, unsigned int rect_height){
	repRef(DEFAULT_CHANNEL)->evaluate(f, rect_ul_x, rect_ul_y, 
					  rect_width, rect_height);}
    const VISImage<T>& evaluate(T (*f)(unsigned int, unsigned int), 
				 unsigned int ch){
        #ifdef DYNAMIC_CHECKING
	if (ch > _channels)
	    ERROR("VISImage<>: channel(int ch) - channel out of range");
        #endif 
	repRef(ch)->evaluate(f);
	return(*this);}


//****Checking bounds
    //these functions return 0 if 'in bounds' or 1 if 'out bounds'
    //x is width and y is height to be checked.  ch is channel to be checked
    //x and y can be floats
    //compareSize compares size of one image to another (other)
    boolean checkBounds(unsigned x, unsigned y) const{
	return((channels() > 0)&&(rep(0)->checkBounds(x, y)));}
    boolean checkBounds(int x, int y) const{
	return((channels() > 0)&&(rep(0)->checkBounds(x, y)));}
    boolean checkBounds(float x, float y) const{
	return((channels() > 0)&&(rep(0)->checkBounds(x, y)));}
    boolean checkBounds(unsigned x, unsigned y, unsigned ch) const{
	return((channels() > ch)&&(rep(ch)->checkBounds(x, y)));}
    int compareSize(const VISImage<T>& other) const{
	    return((_channels == other._channels)
		   &&(rep()->compareSize(other.rep())));}

//****Operators that have been overloaded (other overloaded operators
    // occur outside the class, search for 'NMOO' meaning non-member overloaded
    // operators)
    // every pixel in an image is assigned value with this operator=
    // this was in line, but the SUNCC compiler hangs
    VISImage<T>& operator=(T value); 
    // one image equals another with this operator=
    VISImage<T>& operator=(const VISImage<T>& from){
	assign(from); return(*this);}
    // add, subtract, multiply, divide each pixel in
    // an image by the respective pixel in another image
    VISImage<T> operator+(const VISImage<T>& image) const;
    VISImage<T> operator*(const VISImage<T>& image) const;
    VISImage<T> operator-(const VISImage<T>& image) const;
    VISImage<T> operator/(const VISImage<T>& image) const;
    // exponential function  i.e. raise the image to the power of exponent
    VISImage<T> operator^(const int& exponent) const;
    // add, subtract, multiply, divide each pixel in
    // an image by the respective pixel in another image
    // and assign the result to the image being operated on
    VISImage<T>& operator*=(const VISImage<T>& image);
    VISImage<T>& operator+=(const VISImage<T>& image);
    VISImage<T>& operator-=(const VISImage<T>& image);
    VISImage<T>& operator/=(const VISImage<T>& image);
    // add, subtract, multiply, divide each pixel by
    // a value and assign the result to the image being
    // operated on
    VISImage<T>& operator*=(const T& value);
    VISImage<T>& operator+=(const T& value);
    VISImage<T>& operator-=(const T& value);
    VISImage<T>& operator/=(const T& value);

//****Math functions
    VISImage<T> power(const int& exponent) const;// raise the image to a power
                                                  // operator^ also works
    VISImage<T> sqrt() const;// returns square root of the image
    //divide by value, if division by zero is encountered return zeroCondition
    VISImage<T> div_by(T value, T zeroCondition) const;
    //divide by image, if division by zero is encountered return zeroCondition
    VISImage<T> div(const VISImage<T>& image, T zeroCondition) const;
    VISImage<T> ln() const;// returns natural log of the image
    VISImage<T> exp() const;// exponential of the image
    VISImage<T> pos() const;// replaces negative numbers with zeros
                             // does not change positive numbers
    VISImage<T> neg() const;// replaces positive numbers with zeros
                             // does not change negative numbers
    VISImage<T> abs() const;// absolute value of the image
    VISImage<T> sign() const;// replaces all zeros with +1
                              // replaces all positive numbers with +1
                              // replaces all negative numbers with -1
    T max() const;// returns max value in image
    T min() const;// returns min value in image
    // below min and max functions compare each pixel of two images and
    // return an image containing the min or max values for each pixel
    VISImage<T> min(const VISImage<T> &other) const;
    VISImage<T> max(const VISImage<T> &other) const;
    VISImage<T> min(T value) const;
    VISImage<T> max(T value) const;
    float sum() const;// returns the sum of all pixels in the image
    float average() const;// returns the average of all pixels in the image
    VISImage<T> scale(float value) const;// multiplies each pixel by value

//****General Image Processing functions
    // filtering functions
    VISImage<T> median(int window_size) const;// median filter
    VISImage<T> gauss(float sigma) const;// gaussian smoothing
    //
    // gaussian smoothing done with diffusion to enforce boundary conditions
  VISImage<T> gaussDiffuse(float sigma) const;
    // gaussian smoothing done with diffusion to enforce boundary conditions 
  // done only in the x direction
  VISImage<T> gaussDiffuseX(float sigma) const;
  // done only in the x direction
  VISImage<T> gaussDiffuseY(float sigma) const;
    
    // convolution with any filter you choose
    // NOTE: filters can be created with a number of functions
    // search for 'Filter Creation' in this file
    VISImage<T> convolve(const VISImage<T>& kernel) const;
    //determines the inner product of a mask and an area in an image.
    VISImage<T> mask(const VISImage<T>& kernel) const;
    //determines the inner product of a mask and an area in an image.
    //the center of the mask is indicated by center_x, center_y
    //and the center of the mask is positioned at x_pos, y_pos of the
    //image and then the inner product is calculated.
   float maskFloat(const VISImage<float>& mask, unsigned int center_x, 
	unsigned int center_y, unsigned int x_pos, unsigned int y_pos) const;
    // does the same as maskFloat except an array of masks is applied
    // to the image, and thus an array of floats is returned
    VISArray<float>* maskFloat(const VISArray< VISImage<float> >& masks, 
        unsigned int x_pos, unsigned int y_pos) const;
    //scales an image to the range 0 to 255
    VISImage<T> scaleToRGB();
    // taking derivatives
    VISImage<T> dx() const;//derivative in x direction (along rows)
    VISImage<T> dy() const;//derivative in y direction (along cols)
    VISImage<T> dx(unsigned int order) const;// n-th order derivative
    VISImage<T> dy(unsigned int order) const;// n-th order derivative
    // taking the n-th order derivative in x and y directions
    // derivatives can be of different order in the x vs. y direction
    VISImage<T> derivative(unsigned int order_x, unsigned int order_y) const;
    VISImage<T> derivative(unsigned int order_x, unsigned int order_y, 
			    float scale) const;
    // returns a multichannel image with derivatives up to degree.
    VISImage<T> derivatives(unsigned degree) const;
    VISImage<T> derivatives(unsigned degree, float scale) const;
    // these half derivatives are useful in diffusion processes
    // these things all "wrap", and thus need to be zeroed out along the border
    VISImage<T> dxHalfForward() const;
    VISImage<T> dxHalfBack() const;
    VISImage<T> dyHalfForward() const;
    VISImage<T> dyHalfBack() const;
    // retuns 0's for pixels closest to zero crossings along grid lines
    VISImage<T> zeroCrossings() const;
    // reduces the size of an image by subsampling, taking ever n'th pixel
    // where n = scale an integer
    VISImage<T> reduce(int scale) const;
    // does the same as reduce but averages all pixels around the n'th pixel
    // and places the average value in the returned image
    VISImage<T> reduceAverage(int scale) const;
    // resample scales the size of the image by interpolating between
    // pixels.  the image size can be changed by a scale.
    VISImage<T> resample(float the_scale) const;
    // here the image size is changed to wxh
    VISImage<T> resample(unsigned w, unsigned h) const;
    // here the image size is scaled by scale_x in width and scaled by
    // scale_y in height.  x and y indicate a position on the image upon which
    // the center of the new image grid will be placed
    VISImage<T> resample(float scale_x, float scale_y,
			  float x, float y) const;
    // shot noise is added to the image being operated on
    // percent is the percent of pixels effected by the noise
    // low and high are the limits of the uniformally distributed noise
    // sigma is the sigma value of the gaussian distributed noise
    // NOTE: noise can also be added to an image by adding images created by
    //       noiseUniform or noiseGauss to the image.  search for these
    //       for more information.
    VISImage<T> noiseShot(float percent, T low, T high) const;
    VISImage<T> noiseShot(float percent, float sigma) const;

    // sets a strip of pixels "thickness" to "value" around the border
    VISImage<T> setBorder(T value, unsigned thickness) const;
    // become flat creates one image from an image of multiple channels
    // by shrinking each channel and creating an image that is a
    // 'matrix' of the shrunken channel images.
    VISImage<T> becomeFlat() const;
    // these create bigger images with repeated copies of the smaller ones
    VISImage<T> repeatDown(unsigned num) const;
    VISImage<T> repeatAcross(unsigned num) const;
    // Distance transform - this distance transform returns an image containing
    // distance values.  the distance values indicate the euclidean distance
    // from the respective point in the input image to the closest non zero
    // point in the input image.
    VISImage<float> distanceTrans() const;
    // this distance transform has a max_distance if the distance is
    // greater than max_distance then a zero is placed in the returned image
    VISImage<float> distanceTrans(float max_distance) const;
    // city block distance transform
    VISImage<float> cityBlockDistTrans() const;
    // returns an image of 1's indicating local extrema and 0's elsewhere
    VISImage<int> extrema() const;

//****Specific Task Image Processing Functions
    // this returns one iteration of an anisotropic diffusion    
    VISImage<float> anisoDiffuse(float k) const;
    // this returns one iteration of an anisotropic diffusion    
    VISImage<float> anisoDiffuse(VISImage<float> image_dx, 
				  VISImage<float> image_dy) const;
    // this returns one iteration of an *isotropic* diffusion    
    VISImage<float> diffuse() const;
    // applies watershed algorithm to an image
    // thresh is the value of a threshold applied before performing the watershed
    //        this parameter can simply speed up the algorithm.
    //        set to zero and no pre-thresholding occurs
    // depth is the maximum depth of the watersheds created
    // size is the minimum size of a region.  regions above the size are merged
    // pad is the size of padding placed at the boarder of the image
    // NOTE: the image can be an approximation of the first derivative of
    // the image
    VISImage<int> watershed(float thresh, float depth, int pad){
	    watershed(thresh, depth, -1, pad);}
    VISImage<int> watershed(float thresh, float depth, int size, int pad);
    // returns an image of ones indicating where the flood fill occured
    // all other values are uneffected
    // x and y are the starting "seed" position
    // thresh is the value of the maximum level filled
    VISImage<T> floodFill(T thresh, unsigned x, unsigned y);
    // this one does the flood fill and changes a one value
    // (label_from) to another value (label_to)
    // e.g. floodFill(1,5,4,4) means that all the 1's connected to the
    // pixel at 4,4 by a path of 1's are replaced with 5's
    const VISImage<T>& floodFill(T label_from, T label_to, 
				  unsigned x, unsigned y);
    // returns a bounding box via passed parameters
    VISImage<T> floodFillBoundBox(T thresh, unsigned x, unsigned y,
				  unsigned &x_lo, unsigned &y_lo, 
				  unsigned &x_hi, unsigned &y_hi);
    // Canny edge detection
    VISImage<int> cannyEdges() const;
    VISImage<int> cannyEdges(T threshold) const;

//**** Support Functions (NOT USUALLY USED DIRECTLY BY THE USER)
    VISImage<T> mult(T value) const;
    VISImage<T> add(T value) const;
    VISImage<T> div(T value) const;
    VISImage<T> div_by(T value) const;
    VISImage<T> sub_from(T value) const;
    VISImage<T> sub(T value) const;
    VISImage<T> multAssign(const VISImage<T>& image);
    VISImage<T> addAssign(const VISImage<T>& image);
    VISImage<T> divAssign(const VISImage<T>& image);
    VISImage<T> subAssign(const VISImage<T>& image);
    VISImage<T> multAssign(const T& value);
    VISImage<T> divAssign(const T& value);
    VISImage<T> addAssign(const T& value);
    VISImage<T> subAssign(const T& value);
    VISImage<T> gt(T value) const;
    VISImage<T> gteq(T value) const;
    VISImage<T> lt(T value) const;
    VISImage<T> lteq(T value) const;
    VISImage<T> eq(T value) const;
    void assign(const VISIm& from);
    void assign(const VISImage<T>& from) { VISIm::assign (from);}
};  




unsigned dIndex(unsigned a, unsigned b) ;
//{return(((((a) + (b) + 1)*((a) + (b)))/2) + (b));}


//****Filter Creation
VISImage<float> gauss_col_kernel(float sigma);
VISImage<float> gauss_row_kernel(float sigma);
    // NOTE window sizes are given in #'s of stdeviations 
VISImage<float> gauss_kernel(float sigma, float window_size);
VISImage<float> gauss_kernel(float sigma);
VISImage<float> gauss_col_kernel(float sigma, float window_size);
VISImage<float> gauss_row_kernel(float sigma, float window_size);
VISImage<float> gauss_dx_kernel(int order, float sigma);
VISImage<float> gauss_dy_kernel(int order, float sigma);
VISImage<float> gauss_dx_kernel(int order, float sigma, float window_size);
VISImage<float> gauss_dy_kernel(int order, float sigma, float window_size);
VISImage<float> gauss_col_kernel_small(float sigma);
VISImage<float> gauss_row_kernel_small(float sigma);
    // Derivative filters
VISImage<float> dx_kernel(int order);
VISImage<float> dy_kernel(int order);
    // returns a set of masks or kernels to do derivatives
    // masks are in an array begining with degree=0
    // ask
VISArray< VISImage<float> >* derivativeMasks(unsigned degree);
VISArray< VISImage<float> >* derivativeMasks(unsigned degree, float scale);

    // fills an image with (pseudo)random noise between 0.0 and 1.0
    // returns an image of uniform noise
VISImage<float> noiseUniform(unsigned int width, unsigned int height); 
    // returns an image of gaussian noise
VISImage<float> noiseGauss(unsigned int width, unsigned int height,
			    float stdev);


//**** NMOO -- Non-Member Overloaded Operators
template< class T >
inline VISImage<T> operator+(const VISImage<T>& from, T value){
    return(from.add(value));}
template< class T >
inline VISImage<T> operator+(T value,const VISImage<T>& from){
    return(from.add(value));}
template< class T >
inline VISImage<T> operator*(const VISImage<T>& from, T value){
    return(from.mult(value));}
template< class T >
inline VISImage<T> operator*(T value, const VISImage<T>& from){
    return(from.mult(value));}
template< class T >
inline VISImage<T> operator-(const VISImage<T>& from, T value){
    return(from.sub(value));}
template< class T >
inline VISImage<T> operator-(T value, const VISImage<T>& from){
    return(from.sub_from(value));}
template< class T >
inline VISImage<T> operator/(const VISImage<T>& from, T value){
    return(from.div(value));}
template< class T >
inline VISImage<T> operator/(T value, const VISImage<T>& from){
    return(from.div_by(value));}
template< class T >
inline VISImage<T> operator>(T value, const VISImage<T>& from){
    return(from.lt(value));}
template< class T >
inline VISImage<T> operator>=(T value, const VISImage<T>& from){
    return(from.lteq(value));}
template< class T >
inline VISImage<T> operator<(T value, const VISImage<T>& from){
    return(from.gt(value));}
template< class T >
inline VISImage<T> operator<=(T value, const VISImage<T>& from){
    return(from.gteq(value));}
template< class T >
inline VISImage<T> operator>(const VISImage<T>& from, T value){
    return(from.gt(value));}
template< class T >
inline VISImage<T> operator>=(const VISImage<T>& from, T value){
    return(from.gteq(value));}
template< class T >
inline VISImage<T> operator<(const VISImage<T>& from, T value) {
    return(from.lt(value));}
template< class T >
inline VISImage<T> operator<=(const VISImage<T>& from, T value){
    return(from.lteq(value));}
template< class T >
inline VISImage<T> operator==(T value, const VISImage<T>& from){
    return(from.eq(value));}
template< class T >
inline VISImage<T> operator==(const VISImage<T>& from, T value){
    return(from.eq(value));}

//****Support Functions
template< class T1, class T2 >    
void copy(const VISImRep<T2>* a, VISImRep<T1>* b);
template< class T, class T2 >    
int compareSize(const VISImRep<T>* a, const VISImRep<T2>* b);

template< class T1, class T2 >
void copy(const VISImage<T1>& from, VISImage<T2>& to);

template< class T1, class T2 >
void assignImage(const VISImage<T1>& from, VISImage<T2>& to);

template< class T1, class T2 >
void assignImage(const VISImage<T1> *from, VISImage<T2> *to);

template< class T1, class T2 >
int compareSize(const VISImage<T1>& a, const VISImage<T2>& b);

template< class T1, class T2 >
void copy(const VISImageRep<T1>* from, VISImageRep<T2>* to);
    
template< class T1, class T2 >
int compareSize(const VISImageRep<T1>* a, 
			const VISImageRep<T2>* b);


//#ifdef DEPENDING
//#include "image/image.c"
//#endif

#ifndef MANUAL_INSTANTIATION
#include "image/image.txx"
#endif

#endif




