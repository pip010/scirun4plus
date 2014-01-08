// ******************************************************************
// VISPack. Copyright (c) 1994-2000 Ross Whitaker rtw@utk.edu       *
// For conditions of distribution and use, see the file LICENSE.txt *
// accompanying this distribution.                                  *
// ******************************************************************

// $Id: imageRGBA.h,v 1.1.1.1 2003/02/12 16:51:52 whitaker Exp $


// File:           imageRGBA.h
// Author:         Ross T. Whitaker
// Institution:    The University of Tennessee, Knoxville
// Contents:       Contains a class VISImageRGBA.  This class is inherited
//                 from VISImage<rgba>.  Therefore this class has all the
//                 functionality of an VISImage, but it also has some added
//                 specialized functionality for color images.
// Log of changes: June 30, 1999 -- Added comments



#ifndef iris_imageRGBA_h
#define iris_imageRGBA_h

#include "image/rgba.h"
#include "image/image.h"

// An VISImage<rgba> is an VISImage containing 4 channels of type byte.
// Channel 0 represents the alpha values, channel 1 contains the blue values,
// channel 2 contains the green values, and channel 3 contains the red values.
// VISColorChannel is used to address these channels.  An VISColorChannel
// can have four different values:  R, G, B, or A where R represents red
// (channel 3), G represents green (channel 2), B represents blue (channel 1),
// and A represents alpha (channel 0).

class VISImageRGBA: public VISImage<rgba>
{
    friend
	VISImageRGBA operator*(byte value, const VISImageRGBA& a);
    friend
	VISImageRGBA operator*(const VISImageRGBA& a, byte value);

    friend
	VISImageRGBA operator*(float value, const VISImageRGBA& a);
    friend
	VISImageRGBA operator*(const VISImageRGBA& a, float value);

  private:
  protected:
  public:
    VISImageRGBA(unsigned int w, unsigned int h)
	:VISImage<rgba>(w, h, 1){}
    VISImageRGBA(unsigned int w, unsigned int h, rgba *buf)
	:VISImage<rgba>(w, h, 1, &buf){}
    VISImageRGBA(){}
    VISImageRGBA(const VISImage<rgba>& image)
	:VISImage<rgba>(image){}
//****Interacting with the data of an VISImageRGBA
    // put a byte value into a color channel where the VISColorChannel is
    // R, G, B, or A for red, green, blue, and alpha respectively
    void put(unsigned int x, unsigned int y, VISColorChannel c, byte b){
	    ((repRef())->at(x, y)).at(c) = b;}
    // 'put' an VISImage<byte> into a color channel of the VISImageRGBA
    void putColor(const VISImage<byte>& im, VISColorChannel c){
	    rgba* buf = (repRef())->bufferRef();
	    const byte* buf_byte;
	    buf_byte = (im.rep())->buffer();
	    for (int i = 0; i < rep(0)->size(); i++)
		buf[i].at(c) = buf_byte[i];}
    // 'get' a channel of the VISImageRGBA
    VISImage<byte> channel(const VISColorChannel c) const{
	    VISImage<byte> new_im(width(), height(), 1);
	    
	    const rgba* buf = (rep())->buffer();
	    byte* buf_byte = (new_im.repRef())->bufferRef();
	    for (int i = 0; i < rep()->size(); i++)
		buf_byte[i] = (buf[i])[c];		
	    return(new_im);}
    // creates an VISImageRGBA having the same size as *this
    VISImageRGBA createToSize() const{ 
	    VISImageRGBA new_im(width(), height());
	    return(new_im);}
//****Overloaded Operators
    // Several operators that do not exist in the VISImage class
    // are added to the VISImageRGBA class here.
    // Several assignment operators to set an entire image to a value
    VISImageRGBA& operator=(byte value);
    VISImageRGBA& operator=(const rgba& value);
    // assign one image to another
    VISImageRGBA& operator=(const VISImageRGBA& from){
	assign(from); return(*this);}
    VISImageRGBA& operator=(const VISIm& from);
    // add an VISImageRGBA to another VISImageRGBA
    VISImageRGBA operator+(const VISImageRGBA& image) const{
	VISImageRGBA image_return(width(), height());
	(image_return.repRef())->add((this->rep()), 
				     (image.rep()));
	return(image_return);}
};


#endif











