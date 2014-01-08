#ifndef gfx__local_image_hdr
#define gfx__local_image_hdr


#include "image/image.h"
#include "util/array.h"

#ifdef DEPENDING
#include "localimage.c"
#endif

// NOTE: this is a prototype for developing new METHODS for templated image
// objects.  It is not sufficient for adding new MEMBERS.  For such work
// the constructors would need to be changed as well

template< class T >
class LocalImage: public GfxImage<T>
{
  public:
    LocalImage(unsigned int w, unsigned int h, unsigned int ch)
	:GfxImage<T>(w, h, ch){}
    LocalImage(unsigned int w, unsigned int h, unsigned int ch, T** buf)
	: GfxImage<T>(w, h, ch, buf) {}
    LocalImage(unsigned int w, unsigned int h):GfxImage<T>(w, h) {}
    LocalImage(): GfxImage<T>() {}
    LocalImage(GfxImageRep<T>* rep[], unsigned int ch)
	:GfxImage<T>(rep, ch) {}
    LocalImage(const GfxImage<T>& image):GfxImage<T>(image){}
    LocalImage(GfxImage<T>& a, GfxImage<T>& b): GfxImage<T>(a, b){}

    LocalImage<T>& operator=(const GfxImage<T>& from)
	{
	    assign(from);
	    return(*this);
	}

    LocalImage<T>& operator=(const LocalImage<T>& from)
	{
	    GfxImage<T>::assign(from);
	    return(*this);
	}

    LocalImage<T>& operator=(T value);

    LocalImage<T>& operator*=(const GfxImage<T>& image)
	{
	    multAssign(image);
	    return(*this);
	}

    LocalImage<T>& operator*=(const T& value)
	{
	    multAssign(value);
	    return(*this);
	}

    LocalImage<T>& operator+=(const GfxImage<T>& image)
	{
	    addAssign(image);
	    return(*this);
	}
    LocalImage<T>& operator+=(const T& value)
	{
	    addAssign(value);
	    return(*this);
	}
    LocalImage<T>& operator-=(const GfxImage<T>& image)
	{
	    subAssign(image);
	    return(*this);
	}

    LocalImage<T>& operator-=(const T& value)
	{
	    subAssign(value);
	    return(*this);
	}

    LocalImage<T>& operator/=(const GfxImage<T>& image)
	{
	    divAssign(image);
	    return(*this);
	}

    LocalImage<T>& operator/=(const T& value)
	{
	    divAssign(value);
	    return(*this);
	}
    
    ~LocalImage();

/* put new methods in here */


    LocalImage<T> weightedCurvature() const;
    float maskFloat(const LocalImage<float>& mask, 
				   unsigned int center_x, 
				   unsigned int center_y, 
				   unsigned int x_pos, 
				   unsigned int y_pos
				   ) const;

    Array<float>* maskFloat(const Array< LocalImage<float> >& 
					   masks, 
					   unsigned int x_pos, 
					   unsigned int y_pos
					   ) const;

    void printData() const;

    LocalImage<T> resample(float the_scale) const;
    LocalImage<T> resample(unsigned w, unsigned h) const;
    LocalImage<T> resample(float scale_x, float scale_y,
			 float x, float y) const;

    // take the min and max neighbors over the interpolating footprint
    // i.e. 4 nearest pixels

    T interpMin(float x, float y) const;
    T interpMax(float x, float y) const;
    T LocalImage<T>::interp(float x, float y, GfxImage<float> mask) const;
    
};


#endif
