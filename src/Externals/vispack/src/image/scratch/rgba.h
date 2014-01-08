// ******************************************************************
// VISPack. Copyright (c) 1994-2000 Ross Whitaker rtw@utk.edu       *
// For conditions of distribution and use, see the file LICENSE.txt *
// accompanying this distribution.                                  *
// ******************************************************************

// $Id: rgba.h,v 1.1.1.1 2003/02/12 16:51:53 whitaker Exp $

// File:           image.h
// Author:         Ross T. Whitaker
// Institution:    The University of Tennessee, Knoxville
// Contents:       Contains the rgba class used to create the RGBA type
//                 so that an VISImage<rgba> can be declared and therefore
//                 an VISImageRGBA.  See imageRGBA.h for more details about
//                 the VISImageRGBA class.
// Log of changes: June 30, 1999 -- Added comments

#ifndef VIS_RGBA_H
#define VIS_RGBA_H

typedef unsigned char byte;
typedef enum {R=3, G=2, B=1, A=0} VISColorChannel;

// An VISImage<rgba> is an VISImage containing 4 channels of type byte.
// Channel 0 represents the alpha values, channel 1 contains the blue values,
// channel 2 contains the green values, and channel 3 contains the red values.
// VISColorChannel is used to address these channels.  An VISColorChannel
// can have four different values:  R, G, B, or A where R represents red
// (channel 3), G represents green (channel 2), B represents blue (channel 1),
// and A represents alpha (channel 0).


// these are the weights used to make greyscale data
// these should be declared const float and not #defined
const float R_WEIGHT = 1.0/3.0;
const float G_WEIGHT = 1.0/3.0;
const float B_WEIGHT = 1.0/3.0;
const float A_WEIGHT = 0.0;

#include <math.h>

class rgba
{
  public:
    byte _color[4];

    rgba(byte r,byte g,byte b,byte a) {
	_color[R] = r;
	_color[G] = g;
	_color[B] = b;
	_color[A] = a;
    }
    rgba(byte r,byte g,byte b) {
	_color[R] = r;
	_color[G] = g;
	_color[B] = b;
	_color[A] = 255;
    }
    rgba(float f)
	{
	    _color[R] = _color[G] = _color[B] = (byte)f;
	    _color[A] = 255;
	}

    rgba(double f)
	{
	    _color[R] = _color[G] = _color[B] = (byte)f;
	    _color[A] = 255;
	}
  

    rgba(int i)
	{
	    _color[R] = _color[G] = _color[B] = (byte)i;
	    _color[A] = 255;
	}

    rgba(unsigned i)
	{
	    _color[R] = _color[G] = _color[B] = (byte)i;
	    _color[A] = 255;
	}
    
    
    const byte& operator[] (VISColorChannel c) const {return(_color[c]);}
    byte& at(VISColorChannel c) {return(_color[c]);}
    const byte& itemAt(VISColorChannel c)	const {return(_color[c]);}
    
    const byte& r() const {return(_color[R]);}
    const byte& g() const {return(_color[G]);}
    const byte& b() const {return(_color[B]);}
    const byte& a() const {return(_color[A]);}

    void r(const byte b) {_color[R] = b;}
    void g(const byte b) {_color[G] = b;}
    void b(const byte b) {_color[B] = b;}
    void a(const byte b) {_color[A] = b;}
	

    rgba& operator=(byte value)
	{
	    _color[R] = value;    
	    _color[G] = value;    
	    _color[B] = value;    
	    _color[A] = value;    
	    return(*this);
	}

    rgba& operator=(const rgba& other)
	{
	    _color[R] = other._color[R];
	    _color[G] = other._color[G];
	    _color[B] = other._color[B];
	    _color[A] = other._color[A];
	    return(*this);
	}


    rgba(const rgba& other)
	{
	    _color[R] = other._color[R];
	    _color[G] = other._color[G];
	    _color[B] = other._color[B];
	    _color[A] = other._color[A];
	}
    
    rgba operator+(const rgba& other) const
	{
	    rgba new_rgba;
	    float factor_this = (float)_color[A]/
		(float)(_color[A] + other._color[A]);
	    float factor_other = (float)other._color[A]/
		(float)(_color[A] + other._color[A]);

	    new_rgba._color[R] = (byte)rint(factor_this*_color[R]
		+ factor_other*other._color[R]);
	    new_rgba._color[G] = (byte)rint(factor_this*_color[G]
		+ factor_other*other._color[G]);
	    new_rgba._color[B] = (byte)rint(factor_this*_color[B]
		+ factor_other*other._color[B]);
	    new_rgba._color[A] = _color[A]+other._color[A];
	    return(new_rgba);
	}

    rgba& operator+=(const rgba& other) 
	{
	    float factor_this = (float)_color[A]/
		(float)(_color[A] + other._color[A]);
	    float factor_other = (float)other._color[A]/
		(float)(_color[A] + other._color[A]);

	    _color[R] = (byte)rint(factor_this*_color[R]
		+ factor_other*other._color[R]);
	    _color[G] = (byte)rint(factor_this*_color[G]
		+ factor_other*other._color[G]);
	    _color[B] = (byte)rint(factor_this*_color[B]
		+ factor_other*other._color[B]);
	    _color[A] = _color[A]+other._color[A];
	    return(*this);
	}

    int operator<(const rgba& other) const
	{
	    return((float)(*this) < (float)(other));
	}

    int operator>(const rgba& other) const
	{
	    return((float)(*this) > (float)(other));
	}

// I added the <= and >= operators. Are they needed??? -colas
// Yes, they are - Ross.  6-15-94.

    int operator<=(const rgba& other) const
	{
	    return((float)(*this) <= (float)(other));
	}

    int operator>=(const rgba& other) const
	{
	    return((float)(*this) >= (float)(other));
	}

    int operator!=(const rgba& other) const
	{
	    return((float)(*this) != (float)(other));
	}

    int operator==(const rgba& other) const
	{
	    return((float)(*this) == (float)(other));
	}

    
    rgba operator-(const rgba& other) const
	{
	    rgba new_rgba;
	    new_rgba._color[R] = _color[R]-other._color[R];
	    new_rgba._color[G] = _color[G]-other._color[G];
	    new_rgba._color[B] = _color[B]-other._color[B];
	    new_rgba._color[A] = _color[A]-other._color[A];
	    return(new_rgba);
	}

    rgba operator*(const rgba& other) const
	{
	    rgba new_rgba;
	    new_rgba._color[R] = _color[R]*other._color[R];
	    new_rgba._color[G] = _color[G]*other._color[G];
	    new_rgba._color[B] = _color[B]*other._color[B];
	    new_rgba._color[A] = _color[A]*other._color[A];
	    return(new_rgba);
	}


    rgba operator/(const rgba& other) const
	{
	    rgba new_rgba;
	    new_rgba._color[R] = _color[R]/other._color[R];
	    new_rgba._color[G] = _color[G]/other._color[G];
	    new_rgba._color[B] = _color[B]/other._color[B];
	    new_rgba._color[A] = _color[A]/other._color[A];
	    return(new_rgba);
	}



    operator int() const
	{
	    int int_out = (int)(R_WEIGHT*_color[R] + G_WEIGHT*_color[G]
				+ B_WEIGHT*_color[B] + A_WEIGHT*_color[A]);
	    return(int_out);
	}

    operator unsigned() const
	{
	    // Temporary solution - there is another (better?) way
	    unsigned int_out = (unsigned)(R_WEIGHT*_color[R] + G_WEIGHT*_color[G]
				+ B_WEIGHT*_color[B] + A_WEIGHT*_color[A]);
	    return(int_out);
	}

    operator byte() const
	{
	    byte byte_out;
	    byte_out = (byte)(R_WEIGHT*_color[R] + G_WEIGHT*_color[G]
			      + B_WEIGHT*_color[B] + A_WEIGHT*_color[A]);
	    return(byte_out);
	}

    operator float() const
	{
	    float float_out;
	    float_out = (float)(R_WEIGHT*_color[R] + G_WEIGHT*_color[G]
				+ B_WEIGHT*_color[B] + A_WEIGHT*_color[A]);
	    return(float_out);
	}

    rgba()
	{
	}

};

	rgba operator*(const rgba&, byte value);
	rgba operator*(byte value, const rgba&);
	rgba operator*(const rgba&, float value);
	rgba operator*(float value, const rgba&);

#endif


