#ifndef vis__scan_image_hdr
#define vis__scan_image_hdr


#include "image/image.h"
#include "util/array.h"

//#ifdef DEPENDING
//#include "scanimage.c"
//#endif

// NOTE: this is a prototype for developing new METHODS for templated image
// objects.  It is not sufficient for adding new MEMBERS.  For such work
// the constructors would need to be changed as well

template< class T >
class ScanImage: public VISImage<T>
{
  public:
    ScanImage(unsigned int w, unsigned int h, unsigned int ch)
	:VISImage<T>(w, h, ch){}
    ScanImage(unsigned int w, unsigned int h, unsigned int ch, T** buf)
	: VISImage<T>(w, h, ch, buf) {}
    ScanImage(unsigned int w, unsigned int h):VISImage<T>(w, h) {}
    ScanImage(): VISImage<T>() {}
    ScanImage(VISImageRep<T>* rep[], unsigned int ch)
	:VISImage<T>(rep, ch) {}
    ScanImage(const VISImage<T>& image):VISImage<T>(image){}
    ScanImage(VISImage<T>& a, VISImage<T>& b): VISImage<T>(a, b){}

    ScanImage<T>& operator=(const VISImage<T>& from)
	{
	    assign(from);
	    return(*this);
	}

    ScanImage<T>& operator=(const ScanImage<T>& from)
	{
	    VISImage<T>::assign(from);
	    return(*this);
	}

    ScanImage<T>& operator=(T value);

    ScanImage<T>& operator*=(const VISImage<T>& image)
	{
	    multAssign(image);
	    return(*this);
	}

    ScanImage<T>& operator*=(const T& value)
	{
	    multAssign(value);
	    return(*this);
	}

    ScanImage<T>& operator+=(const VISImage<T>& image)
	{
	    addAssign(image);
	    return(*this);
	}
    ScanImage<T>& operator+=(const T& value)
	{
	    addAssign(value);
	    return(*this);
	}
    ScanImage<T>& operator-=(const VISImage<T>& image)
	{
	    subAssign(image);
	    return(*this);
	}

    ScanImage<T>& operator-=(const T& value)
	{
	    subAssign(value);
	    return(*this);
	}

    ScanImage<T>& operator/=(const VISImage<T>& image)
	{
	    divAssign(image);
	    return(*this);
	}

    ScanImage<T>& operator/=(const T& value)
	{
	    divAssign(value);
	    return(*this);
	}
    
    ~ScanImage();

/* put new methods in here */

//    ScanImage<T> resample(float the_scale) const;
    ScanImage<T> resampleMin(float the_scale) const;
    ScanImage<T> resampleMax(float the_scale) const;
//    ScanImage<T> resample(unsigned w, unsigned h) const;
//    ScanImage<T> resample(float scale_x, float scale_y,
//			 float x, float y) const;

    ScanImage<T> resampleMin(float scale_x, float scale_y,
			 float x, float y) const;
    ScanImage<T> resampleMax(float scale_x, float scale_y,
			 float x, float y) const;
    ScanImage<T> reduceMin(int scale) const;
    ScanImage<T> reduceMax(int scale) const;
    
    // take the min and max neighbors over the interpolating footprint
    // i.e. 4 nearest pixels

    T interpMin(float x, float y) const;
    T interpMax(float x, float y) const;
    T interp(float x, float y, VISImage<float> mask) const;

};

#include "scanimage.txx"

#endif
