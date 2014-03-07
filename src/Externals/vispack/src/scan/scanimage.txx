//#include "scanimage.h"
#include <math.h>
#include <limits.h>
#include "util/mathutil.h"
#include "image/image.h"

/* new methods can be defined here */

template< class T >
ScanImage<T>::~ScanImage()
{
}


template< class T >
ScanImage<T>& ScanImage<T>::operator=(T value)
{
    for (int i = 0; i < VISImage<T>::channels(); i++)// **mdm**
	{
	    ScanImage<T>::repRef(i)->clear(value);// **mdm**
	}
    return(*this);
}


template< class T >
T ScanImage<T>::interpMin(float x, float y) const
{
    unsigned int x_lo, x_hi, y_lo, y_hi;
    x_lo = (unsigned int)floor(x);
    y_lo = (unsigned int)floor(y);
    x_hi = (unsigned int)ceil(x);
    y_hi = (unsigned int)ceil(y); 

    if ((x < 0.0)||(y < 0.0)||(x > ScanImage<T>::width())||(y >
	ScanImage<T>::height())) // **mdm**
	{
	    printf("interp out of bounds %f %f \n", x, y);
	    return((T)0);
	}
	   
    return(VISmin(VISmin(ScanImage<T>::itemAt(x_lo, y_lo), ScanImage<T>::itemAt(x_lo, y_hi)), 
		 VISmin(ScanImage<T>::itemAt(x_hi, y_lo), ScanImage<T>::itemAt(x_hi, y_hi))));// **mdm**
}

template< class T >
T ScanImage<T>::interp(float x, float y, VISImage<float> mask) const
{
    unsigned int x_lo, x_hi, y_lo, y_hi;
    x_lo = (unsigned int)floor(x);
    y_lo = (unsigned int)floor(y);
    x_hi = (unsigned int)ceil(x);
    y_hi = (unsigned int)ceil(y); 

    float x_in = x - (float)x_lo;
    float y_in = y - (float)y_lo;


    if ((x < 0.0)||(y < 0.0)||(x > ScanImage<T>::width())||(y > ScanImage<T>::height()))// **mdm**
	{
	    printf("interp out of bounds %f %f \n", x, y);
	    return((T)0);
	}
    float xl_yl, xh_yh, xl_yh, xh_yl;

    xl_yl = VISImage<T>::maskFloat(mask, mask.width()/2, 
				   mask.height()/2, x_lo, y_lo);

    xh_yl = VISImage<T>::maskFloat(mask, mask.width()/2, 
				   mask.height()/2, x_hi, y_lo);

    xh_yh = VISImage<T>::maskFloat(mask,mask.width()/2, 
				   mask.height()/2,  x_hi,  y_hi);

    xl_yh = VISImage<T>::maskFloat(mask, mask.width()/2, 
				   mask.height()/2, x_lo, y_hi);

    return((T)(((float)1.0 - y_in)*
	       (x_in*xh_yl 
		+ (1.0f - x_in)*xl_yl)
	       + y_in*(x_in*xh_yh 
		       + (1.0f - x_in)*xl_yh)));
}

template< class T > 
T ScanImage<T>::interpMax(float x, float y) const
{

    unsigned int x_lo, x_hi, y_lo, y_hi;
    x_lo = (unsigned int)floor(x);
    y_lo = (unsigned int)floor(y);
    x_hi = (unsigned int)ceil(x);
    y_hi = (unsigned int)ceil(y); 

    if ((x < 0.0)||(y < 0.0)||(x > ScanImage<T>::width())||(y > ScanImage<T>::height()))// **mdm**
	{
	    printf("interp out of bounds %f %f \n", x, y);
	    return((T)0);
	}
	   
    return(VISmax(VISmax(ScanImage<T>::itemAt(x_lo, y_lo), ScanImage<T>::itemAt(x_lo, y_hi)), 
		 VISmax(ScanImage<T>::itemAt(x_hi, y_lo), ScanImage<T>::itemAt(x_hi, y_hi))));// **mdm**

}


template< class T >
ScanImage<T> ScanImage<T>::resampleMax(float the_scale) const
{
    return(resampleMax(the_scale, the_scale, ScanImage<T>::width()/2.0 - 0.5, 
		    ScanImage<T>::height()/2.0 - 0.5));// **mdm**
}

template< class T >
ScanImage<T> ScanImage<T>::resampleMin(float the_scale) const
{
    return(resampleMin(the_scale, the_scale, ScanImage<T>::width()/2.0 - 0.5, 
		    ScanImage<T>::height()/2.0 - 0.5));// **mdm**
}


// this returns a anistropically resampled image scaled according to the size
// indicated by "scale_x" and "scale_y".  It rescales around the point
// given by "x" and "y".

template< class T >
ScanImage<T> ScanImage<T>::resampleMin(float scale_x, float scale_y, 
				  float x, float y) const
{
    unsigned new_w, new_h;
    float new_x, new_y;
    float tmp_x, tmp_y;
    unsigned w = ScanImage<T>::width(), h = ScanImage<T>::height();// **mdm**
    int i, j;
    
    new_w = (unsigned)w*scale_x; //why float
    new_h = (unsigned)h*scale_y;

    new_x = (x + 0.5)*scale_x - 0.5;
    new_y = (y + 0.5)*scale_y - 0.5;

    ScanImage<T> image_return(new_w, new_h);
    for (j = 0; j < new_h; j++)
	for (i = 0; i < new_w; i++)
		{
		    tmp_x = ((i - new_x)/scale_x + x); 
		    tmp_y = ((j - new_y)/scale_y + y); 
		    
		    tmp_x = MAX(tmp_x, 0);
		    tmp_x = MIN(tmp_x, w - 1);

		    tmp_y = MAX(tmp_y, 0);
		    tmp_y = MIN(tmp_y, h - 1);

		    image_return.at(i, j) 
			= interpMin(tmp_x, tmp_y);
		}
    return(image_return);
}

template< class T >
ScanImage<T> ScanImage<T>::resampleMax(float scale_x, float scale_y, 
				  float x, float y) const
{
    unsigned new_w, new_h;
    float new_x, new_y;
    float tmp_x, tmp_y;
    unsigned w = ScanImage<T>::width(), h = ScanImage<T>::height();// **mdm**
    int i, j;
    
    new_w = (unsigned)w*scale_x; //why float
    new_h = (unsigned)h*scale_y;


    new_x = (x + 0.5)*scale_x - 0.5;
    new_y = (y + 0.5)*scale_y - 0.5;

    ScanImage<T> image_return(new_w, new_h);
    for (j = 0; j < new_h; j++)
	for (i = 0; i < new_w; i++)
		{
		    tmp_x = ((i - new_x)/scale_x + x); 
		    tmp_y = ((j - new_y)/scale_y + y); 
		    
		    tmp_x = MAX(tmp_x, 0);
		    tmp_x = MIN(tmp_x, w - 1);

		    tmp_y = MAX(tmp_y, 0);
		    tmp_y = MIN(tmp_y, h - 1);

		    image_return.at(i, j) 
			= interpMax(tmp_x, tmp_y);
		}
    return(image_return);
}

template < class T >
ScanImage<T> ScanImage<T>::reduceMin(int scale) const
{
    ScanImage<T> r = ScanImage<T>(ScanImage<T>::width()/scale, ScanImage<T>::height()/scale, ScanImage<T>::channels());// **mdm**
    int i, j, k, l, m;
    int x_mod, y_mod;
    int x_offset, y_offset;
    const T* buf_in;
    T* buf_out;

    if (scale <= 0)
	return(*this);

    int w = ScanImage<T>::width();// **mdm**
    int h = ScanImage<T>::height();// **mdm**
    int final_w = r.width();

    x_mod = w%scale;
    y_mod = h%scale;

    x_offset = (scale - 1)/2 + (int)ceil((double)x_mod/2.0);
    y_offset = (scale - 1)/2 + (int)ceil((double)y_mod/2.0);

#ifdef DEBUG
    printf("Image::reduce x_mod %d x_offset %d, y_mod %d y_offset %d\n", 
	   x_mod, x_offset, y_mod, y_offset);
#endif
    
    T the_min;
//    int position = 0;
    for (i = 0; i < ScanImage<T>::channels(); i++)// **mdm**
	{
	    buf_in = ScanImage<T>::rep(i)->buffer() + w*(y_offset)// **mdm** 
		+ x_offset;
//	    position = w*(y_offset)+ x_offset;
	    buf_out = r.repRef(i)->bufferRef();
	    for (j = 0; j < r.height(); j++)
		{
		    for (k = 0; k < final_w; k++)
			{
//
// find the min of the neigborhood
//
			    the_min = *buf_in;
			    for (l = 0; l < scale; l++)
				for (m = 0; m < scale; m++)
				    the_min = VISmin(*(buf_in + m + 
						      l*w), the_min);
			    (*buf_out++) = the_min;
			    buf_in += scale;
//			    position += scale;
			}
		    buf_in += x_mod + (scale - 1)*w;
//		    position += x_mod + (scale - 1)*w;
		}
	}
    return(r);
}

template < class T >
ScanImage<T> ScanImage<T>::reduceMax(int scale) const
{
    ScanImage<T> r = ScanImage<T>(ScanImage<T>::width()/scale, ScanImage<T>::height()/scale, ScanImage<T>::channels());// **mdm**
    int i, j, k, l, m;
    int x_mod, y_mod;
    int x_offset, y_offset;
    const T* buf_in;
    T* buf_out;

    if (scale <= 0)
	return(*this);

    int w = ScanImage<T>::width();// **mdm**
    int h = ScanImage<T>::height();// **mdm**
    int final_w = r.width();

    x_mod = w%scale;
    y_mod = h%scale;

    x_offset = (scale - 1)/2 + (int)ceil((double)x_mod/2.0);
    y_offset = (scale - 1)/2 + (int)ceil((double)y_mod/2.0);

#ifdef DEBUG
    printf("Image::reduce x_mod %d x_offset %d, y_mod %d y_offset %d\n", 
	   x_mod, x_offset, y_mod, y_offset);
#endif
    
    T the_max;
//    int position = 0;
    for (i = 0; i < ScanImage<T>::channels(); i++)// **mdm**
	{
	    buf_in = ScanImage<T>::rep(i)->buffer() + w*(y_offset)// **mdm** 
		+ x_offset;
//	    position = w*(y_offset)+ x_offset;
	    buf_out = r.repRef(i)->bufferRef();
	    for (j = 0; j < r.height(); j++)
		{
		    for (k = 0; k < final_w; k++)
			{
//
// find the max of the neigborhood
//
			    the_max = *buf_in;
			    for (l = 0; l < scale; l++)
				for (m = 0; m < scale; m++)
				    the_max = VISmax(*(buf_in + m + 
						      l*w), the_max);
			    (*buf_out++) = the_max;
			    buf_in += scale;
//			    position += scale;
			}
		    buf_in += x_mod + (scale - 1)*w;
//		    position += x_mod + (scale - 1)*w;
		}
	}
    return(r);
}



// *************************************************
// ****************  END ***************************
// *************************************************

#ifdef not_for_now
// these things should all be in the base class

// **************  where does this go ????? ************

// void ScanImage<float>::printData() const
//{
//    for (int i = 0; i < height(); i++)
//	{
//	for (int j = 0; j < width(); j++)
//	    printf("%3.2f ", itemAt(j, i));
//	printf("\n");
//	}
//}



template< class T >
ScanImage<T> ScanImage<T>::resample(unsigned w, unsigned h) const
{
    float scale_w, scale_h;
    int i, j;
    
    scale_w = (float)w/(float)width();
    scale_h = (float)h/(float)height();

    return(resample(scale_w, scale_h, 
		    (float)width()/2.0 - 0.5, 
		    (float)height()/2.0 - 0.5
		    ));
}

template< class T >
ScanImage<T> ScanImage<T>::resample(float the_scale) const
{
    return(resample(the_scale, the_scale, width()/2.0 - 0.5, 
		    height()/2.0 - 0.5);
}

template< class T >
ScanImage<T> ScanImage<T>::resample(float scale_x, float scale_y, 
				  float x, float y) const
{
    unsigned new_w, new_h;
    float float_w, float_h;
    float new_x, new_y;
    float tmp_x, tmp_y;
    unsigned w = width(), h = height();
    int i, j;
    
    new_w = float_w = (float)w*scale_x;
    new_h = float_h = (float)h*scale_y;


    new_x = (x + 0.5)*scale_x - 0.5;
    new_y = (y + 0.5)*scale_y - 0.5;

    ScanImage<T> image_return(new_w, new_h);
    for (j = 0; j < new_h; j++)
	for (i = 0; i < new_w; i++)
		{
		    tmp_x = ((i - new_x)/scale_x + x); 
		    tmp_y = ((j - new_y)/scale_y + y); 
		    
		    tmp_x = MAX(tmp_x, 0);
		    tmp_x = MIN(tmp_x, w - 1);

		    tmp_y = MAX(tmp_y, 0);
		    tmp_y = MIN(tmp_y, h - 1);

		    image_return.at(i, j) 
			= interp(tmp_x, tmp_y);
		}
    return(image_return);
}



#endif
