// ******************************************************************
// VISPack. Copyright (c) 1994-2000 Ross Whitaker rtw@utk.edu       *
// For conditions of distribution and use, see the file LICENSE.txt *
// accompanying this distribution.                                  *
// ******************************************************************

// $Id: image_type.h,v 1.1.1.1 2003/02/12 16:51:52 whitaker Exp $

// File:           image_type.h
// Author:         Ross T. Whitaker
// Institution:    The University of Tennessee, Knoxville
// Contents:       Contains all the possible types for an image.  The possible
//                 types are:  byte, int, unsigned, float, and RGBA.
// Log of changes: June 30, 1999 -- Added comments


// Most of the types here are familiar, but the rgba type is not.
// An VISImage<rgba> is an VISImage containing 4 channels of type byte.
// One channel of each color (red, green, and blue) and one channel for
// an alpha value, which can be used to represent intensity.  For more
// information on this type see rgba.h.  Also, imageRGBA.h contains
// an VISImageRGBA class, which may be more desirable than using
// VISImage<rgba>'s.


#ifndef iris_image_type_h
#define iris_image_type_h

#include "rgba.h"

#define OTHER_MIN (0)
#define OTHER_MAX (255)

#define COLOR_MIN (0)
#define COLOR_MAX (UCHAR_MAX)

#define FLOAT_MIN (0.0)
#define FLOAT_MAX (1.0)

#include "image/im.h"
#include "limits.h"

template< class T >
class VISImageType
{
  public:
    VISImageType() {}
    VISIm::Type whatAmI() const { return(VISIm::OTHER); }
    float rangeMin() const { return(OTHER_MIN); }
    float rangeMax() const { return(OTHER_MAX); }
};

template< >
class VISImageType<byte>
{
  public:
    VISImageType() {}
    VISIm::Type whatAmI() const { return(VISIm::BYTE); }
    float rangeMin() const {return(0);}
    float rangeMax() const {return(UCHAR_MAX);}
};

template< >
class VISImageType<int>
{
  public:
    VISImageType() {}
    VISIm::Type whatAmI() const { return(VISIm::INT); }
    float rangeMin() const {return(INT_MIN);}
    float rangeMax() const {return(INT_MAX);}
};

template< >
class VISImageType<unsigned>
{
  public:
    VISImageType() {}
    VISIm::Type whatAmI() const { return(VISIm::UNSIGNED); }
    float rangeMin() const {return(0);}
    float rangeMax() const {return(UINT_MAX);}
};

template< >
class VISImageType<short>
{
  public:
    VISImageType() {}
    VISIm::Type whatAmI() const { return(VISIm::SHORT); }
    float rangeMin() const {return(SHRT_MIN);}
    float rangeMax() const {return(SHRT_MAX);}
};

template< >
class VISImageType<float>
{
  public:
    VISImageType() {}
    VISIm::Type whatAmI() const { return(VISIm::FLOAT); }
    float rangeMin() const {return(FLOAT_MIN);}
    float rangeMax() const {return(FLOAT_MAX);}
};

template< >
class VISImageType<rgba>
{
  public:
    VISImageType() {}
    VISIm::Type whatAmI() const { return(VISIm::RGBA); }
    float rangeMin() const {return(COLOR_MIN);}
    float rangeMax() const {return(COLOR_MAX);}
};


#endif






