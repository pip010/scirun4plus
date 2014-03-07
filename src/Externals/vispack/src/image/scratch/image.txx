#ifndef	IMAGE_C
#define IMAGE_C

// ******************************************************************
// VISPack. Copyright (c) 1994-2000 Ross Whitaker rtw@utk.edu       *
// For conditions of distribution and use, see the file LICENSE.txt *
// accompanying this distribution.                                  *
// ******************************************************************

// $Id: image.txx,v 1.1.1.1 2003/02/12 16:51:52 whitaker Exp $

//
// *******************

//VISImage

// *******************
//
//#include "image/image.h"

#include <math.h>

#include "image_type.h"
#include "util/mathutil.h"
#include "util/array.h"
//#include "image/rgba.h"
#include "image/indexlist.h"
#include "image/imagefile.h"

template< class T >
VISImage<T>::~VISImage()
{
// let the base class do this
/**
    if (_rep)
	{
	    for (int i = 0; i < _channels; i++)
		unref(_rep[i]);
	    delete _rep;
	    _rep = NULL;
	}
***/
}

template< class T>
VISImage<T> VISImage<T>::noiseUniform(float level)
{return((*this) + VISImage<T>(level*noiseUniform(_width, _height)));}

template< class T>
VISImage<T> VISImage<T>::noiseGauss(float stdev)
{return((*this) + VISImage<T>(::noiseGauss(width(), height(), stdev)));}

template< class T>
VISImage<T>::VISImage(VISImageRep<T> const* const the_rep[], 
		      unsigned int ch)
{
    initialize(the_rep, ch);
}
 
template< class T >
VISImage<T>::VISImage(unsigned int w, unsigned int h, 
		      unsigned int ch, T** buf)
{
    VISImageType<T> type_tmp;
    _type = type_tmp.whatAmI();

// this is now done in the parent constructor
//    _channels = ch;
//    _rep = new VISImRep<T>*[_channels];

    VISIm::initialize(ch);
    for (int i = 0; i < _channels; i++)
	{
	    putRep(new VISImageRep<T>(w, h, buf[i]), i);
	}
}

template< class T >
VISImage<T>::VISImage(const VISImage<T>& image)
{
    VISIm::initialize(image.repVISArray(), image.channels());
    VISImageType<T> type_tmp;
    _type = type_tmp.whatAmI();
}

template< class T >
VISImage<T>::VISImage(const VISIm& im)
{
    VISIm::initialize(0);
    VISImageType<T> type_tmp;
    _type = type_tmp.whatAmI();
    assign(im);
}

template<class T>
VISImage<T> VISImage<T>::addAssign(const T& value)
{
    for (int j = 0; j < _channels; j++)
	(this->repRef(j))->add(this->rep(j), value);
    return(*this);
}

template<class T>
VISImage<T> VISImage<T>::subAssign(const T& value)
{
    for (int j = 0; j < _channels; j++)
	(this->repRef(j))->sub(this->rep(j), value);
    return(*this);
}


template<class T>
VISImage<T> VISImage<T>::multAssign(const T& value)
{
    for (int j = 0; j < _channels; j++)
	(this->repRef(j))->mult(this->rep(j), value);
    return(*this);
}

template<class T>
VISImage<T> VISImage<T>::multAssign(const VISImage<T>& image)
{
    for (int j = 0; j < _channels; j++)
	(this->repRef(j))->mult(this->rep(j), image.rep(j));
    return(*this);
}

template<class T>
VISImage<T> VISImage<T>::divAssign(const T& value)
{
    for (int j = 0; j < _channels; j++)
	(this->repRef(j))->div(this->rep(j), value);
    return(*this);
}

template<class T>
VISImage<T> VISImage<T>::divAssign(const VISImage<T>& image)
{
    for (int j = 0; j < _channels; j++)
	(this->repRef(j))->div(this->rep(j), image.rep(j));
    return(*this);
}

template<class T>
VISImage<T>  VISImage<T>::addAssign(const VISImage<T>& image)
{
    for (int j = 0; j < _channels; j++)
	(this->repRef(j))->add(this->rep(j), image.rep(j));
    return(*this);
}


template<class T>
VISImage<T>  VISImage<T>::subAssign(const VISImage<T>& image)
{
    for (int j = 0; j < _channels; j++)
	(this->repRef(j))->sub(this->rep(j), image.rep(j));
    return(*this);
}

template< class T >
void VISImage<T>::assign(const VISIm& from)
{
    if (!(&from == this))
	{
	    rechannel(from.channels());
	    if (_type == from.type())
		{
		    for (int i = 0; i < _channels; i++)
			ref(_rep[i] = (VISImageRep<float>*)from.rep(i));
		}
	    else
		{
		    VISImageRep<float> *rep_float;
		    VISImageRep<byte> *rep_byte;
		    VISImageRep<unsigned> *rep_unsigned;
		    VISImageRep<int> *rep_int;
		    VISImageRep<short> *rep_short;
		    VISImageRep<rgba> *rep_rgba;
		    int i;

		switch(from.type())
		    {
		      case FLOAT:
			for (i = 0; i < _channels; i++)
			    {
				rep_float = (VISImageRep<float>*)from.rep(i);
				ref(_rep[i] = new VISImageRep<T>
				    (rep_float->width(), rep_float->height()));
				::copy(rep_float, repRef(i));
			    }
			break;
		      case BYTE:
			for (i = 0; i < _channels; i++)
			    {
				rep_byte = (VISImageRep<byte>*)from.rep(i);
				ref(_rep[i] = new VISImageRep<T>
				    (rep_byte->width(), rep_byte->height()));
				::copy(rep_byte, repRef(i));
			    }
			break;
		      case UNSIGNED:
			for (i = 0; i < _channels; i++)
			    {
				rep_unsigned 
				    = (VISImageRep<unsigned>*)from.rep(i);
				ref(_rep[i] = new VISImageRep<T>
				    (rep_unsigned->width(), 
				     rep_unsigned->height()));
				::copy(rep_unsigned, repRef(i));
			    }
			break;
		      case INT:
			for (i = 0; i < _channels; i++)
			    {
				rep_int = (VISImageRep<int>*)from.rep(i);
				ref(_rep[i] = new VISImageRep<T>
				    (rep_int->width(), rep_int->height()));
				::copy(rep_int, repRef(i));
			    }
			break;
		    case SHORT:
			for (i = 0; i < _channels; i++)
			    {
				rep_short = (VISImageRep<short>*)from.rep(i);
				ref(_rep[i] = new VISImageRep<T>
				    (rep_short->width(), rep_short->height()));
				::copy(rep_short, repRef(i));
			    }
			break;
		      case RGBA:
			for (i = 0; i < _channels; i++)
			    {
				rep_rgba = (VISImageRep<rgba>*)from.rep(i);
				ref(_rep[i] = new VISImageRep<T>
				    (rep_rgba->width(), rep_rgba->height()));
				::copy(rep_rgba, repRef(i));
			    }
			break;
		      case OTHER:
			for (i = 0; i < _channels; i++)
			    {
				ref(_rep[i] = new VISImageRep<T>(0,0));
			    }
			WARN("VISImage<T>::assign - auto conversion from type other not supported\n");
			break;
		      case NONE:
			for (i = 0; i < _channels; i++)
			    {
				ref(_rep[i] = new VISImageRep<T>(0,0));
			    }
			WARN("VISImage<T>::assign - auto conversion from type none not supported\n");
			break;
		    }
		}
	}
}

  
template< class T >
void VISImage<T>::copy(const VISImage<T>& a)
{
  if (compareSize(a))
      for (int i = 0; i < _channels; i++)
	  repRef(i)->copy(a.rep(i));
  else
      ERROR("VISImage<>: copy - image sizes not compatible");
}


#ifdef not_for_now    // see comment in header file
template<class T>
VISImage<T>& VISImage<T>::operator=(const VISIm& from)
{
    printf("got type assign \n");
    if (from.type() == type())
	{
	    printf("got equal type assign\n");
	    assign(*((VISImage<T>*)(&from)));
	}
    else
#ifdef AUTO_IMAGE_CONVERSION
	switch(from.type())
	    {
	      case VISIm::BYTE:
		assignImage(*((VISImage<byte>*)(&from)), *this);
		break;
	      case VISIm::INT:
		assignImage((*(VISImage<int>*)(&from)), *this);
		break;
	      case VISIm::FLOAT:
		assignImage((*(VISImage<float>*)(&from)), *this);
		break;
	      default:
		ERROR("unrecognized image type conversion\n");
	    }
#else
    {
	printf("got unequal type assign\n");
	ERROR("automatic image conversion must be set as compiler option");
    }
#endif
    return(*this);
}
#endif


template< class T >
VISImage<T> VISImage<T>::operator+(const VISImage<T>& image) const
{
    VISImage<T> image_return = image.createToSize();
    
    if (compareSize(image))
	{
	    for (int i = 0; i < _channels; i++)
		image_return.repRef(i)->add((this->rep(i)), 
					 (image.rep(i)));
	}
    else 
	{
	    ERROR("VISIm:operator +; image size mismatch");
	}
    return(image_return);
}

template< class T >
VISImage<T> VISImage<T>::min(const VISImage<T>& image) const
{
    
    VISImage<T> image_return = image.createToSize();
    
    if (compareSize(image))
	{
	    for (int i = 0; i < _channels; i++)
		image_return.repRef(i)->min((this->rep(i)), 
					    (image.rep(i)));
	}
    else 
	{
	    ERROR("VISIm:min; image size mismatch");
	}
    return(image_return);

}

template< class T >
VISImage<T> VISImage<T>::max(const VISImage<T>& image) const
{
    
    VISImage<T> image_return = image.createToSize();
    
    if (compareSize(image))
	{
	    for (int i = 0; i < _channels; i++)
		image_return.repRef(i)->max((this->rep(i)), 
					    (image.rep(i)));
	}
    else 
	{
	    ERROR("VISIm:max; image size mismatch");
	}
    return(image_return);

}


template< class T >
VISImage<T> VISImage<T>::operator*(const VISImage<T>& image) const
{
    VISImage<T> image_return = image.createToSize();
    
    if (compareSize(image))
	{
	    for (int i = 0; i < _channels; i++)
		image_return.repRef(i)->mult((this->rep(i)), 
					  (image.rep(i)));
	}
    else 
	{
	    ERROR("VISIm:operator*; image size mismatch");
	}
    return(image_return);
}


template< class T >
VISImage<T> VISImage<T>::operator-(const VISImage<T>& image) const
{
    VISImage<T> image_return = image.createToSize();
    
    if (compareSize(image))
	{
	    for (int i = 0; i < _channels; i++)
		image_return.repRef(i)->sub((this->rep(i)), 
					 (image.rep(i)));
	}
    else 
	{
	    ERROR("VISIm:operator -; image size mismatch");
	}
    return(image_return);
}


template< class T >
VISImage<T> VISImage<T>::operator/(const VISImage<T>& image) const
{
    VISImage<T> image_return = image.createToSize();
    
    if (compareSize(image))
	{
	    for (int i = 0; i < _channels; i++)
		image_return.repRef(i)->div((this->rep(i)), 
					    (image.rep(i)));
	}
    else 
	{
	    ERROR("VISIm:operator /; image size mismatch");
	}
    return(image_return);
}

template< class T >
VISImage<T> VISImage<T>::div(const VISImage<T>& image, 
			       T zeroCondition) const
{
    VISImage<T> image_return = image.createToSize();
    
    if (compareSize(image))
	{
	    for (int i = 0; i < _channels; i++)
		image_return.repRef(i)->div((this->rep(i)), 
					    (image.rep(i)), zeroCondition);
	}
    else 
	{
	    ERROR("VISIm:operator /; image size mismatch");
	}
    return(image_return);
}


template< class T >
VISImage<T>::VISImage(const VISImage<T>& a, const VISImage<T>& b)
{	
    if (!(a.rep(DEFAULT_CHANNEL)->compareSize(b.rep(DEFAULT_CHANNEL))))
/**	ERROR("VISImage combine-constructor \n"); **/
	printf("VISImage combine-constructor \n");
    else		
	{
	    VISImageType<T> type_tmp;
	    _type = type_tmp.whatAmI();
//
//	    _channels = a._channels + b._channels;
//	    _rep = new VISImRep<T>*[_channels];
//
	    VISIm::initialize(a.channels() + b.channels());
	    for (int i = 0; i < a._channels; i++)
		putRep((VISImageRep<T>*)a.rep(i), i);
	    for (int i = 0; i < b._channels; i++)
		putRep((VISImageRep<T>*)b.rep(i), i + a._channels);
	}
}


template< class T >
void VISImage<T>::initialize(unsigned int w, unsigned int h, unsigned int ch)
{
    VISIm::initialize(ch);
    VISImageType<T> type_tmp;
    _type = type_tmp.whatAmI();
//
//    _channels = ch;
//    _rep = new VISImageRep<T>*[_channels];
//
    for (int i = 0; i < _channels; i++)
	{
	    putRep(new VISImageRep<T>(w, h), i);
	}
}



template< class T >
void VISImage<T>::initialize(VISImageRep<T> const * const the_rep[],  unsigned int ch)
{
    VISIm::initialize(ch);
    VISImageType<T> type_tmp;
    _type = type_tmp.whatAmI();

//
//    
//
//    if (_channels > 0)
//	
//    else 
//	_rep = NULL;
//
    for (int i = 0; i < channels(); i++)
	{
	    putRep(the_rep[i], i);
	}
}

template< class T >
T VISImage<T>::max() const
{
    T tmp;
    T max_tmp;

    if (_channels > 0)
	max_tmp = rep(0)->max();

    for (int i = 1; i < _channels; i++)
	{
	    tmp = rep(i)->max();
	    if (tmp > max_tmp)
		max_tmp = tmp;
	}
    return(max_tmp);
}

template< class T >
T VISImage<T>::min() const
{
    T tmp;
    T min_tmp;

    if (_channels > 0)
	min_tmp = rep(0)->min();

    for (int i = 1; i < _channels; i++)
	{
	    tmp = rep(i)->min();
	    if (tmp < min_tmp)
		min_tmp = tmp;
	}
    return(min_tmp);
}

template< class T >
float VISImage<T>::sum() const
{
    float total = 0.0;

    for (int i = 0; i < _channels; i++)
	{
	    total += rep(i)->sum();
	}
    return(total);
}

template< class T >
VISImage<T> VISImage<T>::scaleToRGB()
{
   if( type() == VISIm::RGBA ) return *this;
   register T ImgMax = max(), ImgMin = min();
   register float scale = 255.0/(ImgMax-ImgMin);
   VISImage<T> image_return = createToSize();

   for (int c = 0; c < _channels; c++) {
       register T *ip = (T*) this->rep(c)->buffer();
       register T *op = image_return.repRef(c)->bufferRef();
       for( register int i=width()*height(); i; --i, ip++, op++ ){
	   *op = (T)((float)(*ip-ImgMin))*scale;
       }
   }

   return image_return;
}

//
// *******************

//VISRep

// *******************
//

template< class T >
void VISImageRep<T>::initialize(unsigned int w, unsigned int h,
				T* buffer_in)
{
    VISImRep<T>::initialize(w*h, buffer_in);
    _width = w;
    _height = h;
}

template< class T >
void VISImageRep<T>::evaluate(T (*f)(unsigned int, unsigned int))
{
    for (int i = 0; i < _width; i++)
	for (int j = 0; j < _height; j++)
	    _buffer[index(i, j)] = f(i, j);
}

template< class T >
void VISImageRep<T>::evaluate(T (*f)(unsigned int, unsigned int), 
	      unsigned int rect_ul_x, 
	      unsigned int rect_ul_y, 
	      unsigned int rect_width, 
	      unsigned int rect_height)
    {
#ifdef DYNAMIC_CHECKING
	if (!(((rect_ul_x + rect_width) < _width)
	      &&((rect_ul_y + rect_height) < _height)))
	    {
		ERROR("Image evaluate patch out of bounds\n");
	    }
	else
#endif
/* fix _height in image.c */
	    for (int i = 0; i < rect_height; i++)
		for (int j = 0; j < rect_width; j++)
		    _buffer[index((j + rect_ul_x), (i + rect_ul_y))] = f(j, i);
	
    }


	

template< class T >
VISImage<T> VISImage<T>::operator^(const int& exponent) const
{
    int j;
    VISImage<T> im_return = createToSize();
    for (j = 0; j < _channels; j++)
	{
	    (im_return.repRef(j))->power(this->rep(j), exponent);
	}
    return(im_return);
}

template< class T >
VISImage<T> VISImage<T>::power(const int& exponent) const
{
    int j;
    VISImage<T> im_return = createToSize();
    for (j = 0; j < _channels; j++)
	{
	    (im_return.repRef(j))->power(this->rep(j), exponent);
	}
    return(im_return);
}



template<class T>
VISImage<T>& VISImage<T>::operator*=(const VISImage<T>& image)
{
    multAssign(image);
    return(*this);
}

template<class T>
VISImage<T>& VISImage<T>::operator*=(const T& value)
{
    multAssign(value);
    return(*this);
}


template<class T>
VISImage<T>& VISImage<T>::operator+=(const VISImage<T>& image)
{
    addAssign(image);
    return(*this);
}

template<class T>
VISImage<T>& VISImage<T>::operator+=(const T& value)
{
    addAssign(value);
    return(*this);
}


template<class T>
VISImage<T>& VISImage<T>::operator-=(const VISImage<T>& image)
{
    subAssign(image);
    return(*this);
}

template<class T>
VISImage<T>& VISImage<T>::operator-=(const T& value)
{
    subAssign(value);
    return(*this);
}


template<class T>
VISImage<T>& VISImage<T>::operator/=(const VISImage<T>& image)
{
    divAssign(image);
    return(*this);
}


template<class T>
VISImage<T>& VISImage<T>::operator/=(const T& value)
{
    multAssign((T)1 / value);
    return(*this);
}


template< class T >
VISImage<T>  VISImage<T>::add(T value) const
{
    VISImage<T> image_return = createToSize();
    for (int i = 0; i < image_return._channels; i++)
	(image_return.repRef(i))->add((rep(i)), 
				    value);
    return(image_return);
}

template< class T >
VISImage<T>  VISImage<T>::sub(T value) const
{
    VISImage<T> image_return = createToSize();
    for (int i = 0; i < image_return._channels; i++)
	(image_return.repRef(i))->sub((rep(i)), 
				    value);
    return(image_return);
}

template< class T >
VISImage<T>  VISImage<T>::sub_from(T value) const
{
    VISImage<T> image_return = createToSize();
    for (int i = 0; i < image_return._channels; i++)
	(image_return.repRef(i))->sub(value, rep(i));
    return(image_return);
}

template< class T >
VISImage<T> VISImage<T>::mult(T value) const
{
    VISImage<T> image_return = createToSize();
    for (int i = 0; i < image_return._channels; i++)
	(image_return.repRef(i))->mult((rep(i)), 
				    value);
    return(image_return);
}

template< class T >
VISImage<T>  VISImage<T>::div(T value) const
{
    VISImage<T> image_return = createToSize();
    for (int i = 0; i < image_return._channels; i++)
	(image_return.repRef(i))->div((rep(i)), 
				      value);
    return(image_return);
}

template< class T >
VISImage<T>  VISImage<T>::div_by(T value) const
{
    VISImage<T> image_return = createToSize();
    for (int i = 0; i < image_return._channels; i++)
	(image_return.repRef(i))->div(value, (rep(i)));
    return(image_return);
}

template< class T >
VISImage<T>  VISImage<T>::div_by(T value, T zeroCondition) const
{
    VISImage<T> image_return = createToSize();
    for (int i = 0; i < image_return._channels; i++)
	(image_return.repRef(i))->div(value, (rep(i)), zeroCondition);
    return(image_return);
}


template< class T >
VISImage<T>  VISImage<T>::gt(T value) const
{
    VISImage<T> image_return = createToSize();
    for (int i = 0; i < image_return._channels; i++)
	(image_return.repRef(i))->gt((rep(i)), 
				    value);
    return(image_return);
}

template< class T >
VISImage<T>  VISImage<T>::gteq(T value) const
{
    VISImage<T> image_return = createToSize();
    for (int i = 0; i < image_return._channels; i++)
	(image_return.repRef(i))->gteq((rep(i)), 
				    value);
    return(image_return);
}

template< class T >
VISImage<T>  VISImage<T>::eq(T value) const
{
    VISImage<T> image_return = createToSize();
    for (int i = 0; i < image_return._channels; i++)
	(image_return.repRef(i))->eq((rep(i)), 
				    value);
    return(image_return);
}

template< class T >
VISImage<T>  VISImage<T>::lt(T value) const
{
    VISImage<T> image_return = createToSize();
    for (int i = 0; i < image_return._channels; i++)
	(image_return.repRef(i))->lt((rep(i)), 
				    value);
    return(image_return);
}

template< class T >
VISImage<T>  VISImage<T>::lteq(T value) const
{
    VISImage<T> image_return = createToSize();
    for (int i = 0; i < image_return._channels; i++)
	(image_return.repRef(i))->lteq((rep(i)), 
				    value);
    return(image_return);
}


template<class T>
VISImage<T> VISImage<T>::dx(unsigned int order) const
{
  VISImage<T> local_dx_kernel;
  VISImage<float> float_dx_kernel = dx_kernel(order);
  assignImage((float)2.0*float_dx_kernel, local_dx_kernel);
  return(convolve(local_dx_kernel).divAssign((T)2)); 
}

template<class T>
VISImage<T> VISImage<T>::dxHalfForward() const
{
    VISImage<T> r = this->createToSize();
    const T* buff_in[2];
    T *buff_out;
    int h = r.height();
    int w = r.width();
    for (int i = 0; i < _channels; i++)
	{
	    buff_in[0] = rep(i)->buffer();
	    buff_in[1] = buff_in[0] + 1;
	    buff_out = (r.repRef(i))->bufferRef();
	    for (int j = 0; j < h; j++)
		{
		    for (int k = 0; k < w - 1; k++)
			*(buff_out++) = *(buff_in[1]++) - *(buff_in[0]++);
		    *(buff_out++) = (T)0;
		    buff_in[1]++; buff_in[0]++;
		}
	}
    return(r);
}

template<class T>
VISImage<T> VISImage<T>::dyHalfForward() const
{
    VISImage<T> r = this->createToSize();
    const T* buff_in[2];
    T *buff_out;
    int h = r.height();
    int w = r.width();
    for (int i = 0; i < _channels; i++)
	{
	    buff_in[0] = rep(i)->buffer();
	    buff_in[1] = buff_in[0] + w;
	    buff_out = (r.repRef(i))->bufferRef();
	    for (int j = 0; j < h - 1; j++)
		{
		    for (int k = 0; k < w; k++)
			*(buff_out++) = *(buff_in[1]++) - *(buff_in[0]++);
		}
	    for (int k = 0; k < w; k++)
		{
		*(buff_out++) = (T)0;
		buff_in[1]++; buff_in[0]++;
		}
	}
    return(r);
}

template<class T>
VISImage<T> VISImage<T>::dxHalfBack() const
{
    VISImage<T> r = this->createToSize();
    const T* buff_in[2];
    T *buff_out;
    int h = r.height();
    int w = r.width();
    for (int i = 0; i < _channels; i++)
	{
	    buff_in[0] = rep(i)->buffer() - 1;
	    buff_in[1] = buff_in[0] + 1;
	    buff_out = (r.repRef(i))->bufferRef();
	    for (int j = 0; j < h; j++)
		{
		    *(buff_out++) = (T)0;
		    buff_in[1]++; buff_in[0]++;
		    for (int k = 1; k < w; k++)
			*(buff_out++) = *(buff_in[1]++) - *(buff_in[0]++);
		}
	}
    return(r);
}

template<class T>
VISImage<T> VISImage<T>::dyHalfBack() const
{
    VISImage<T> r = this->createToSize();
    const T* buff_in[2];
    T *buff_out;
    int h = r.height();
    int w = r.width();
    for (int i = 0; i < _channels; i++)
	{
	    buff_in[0] = rep(i)->buffer() - w;
	    buff_in[1] = buff_in[0] + w;
	    buff_out = (r.repRef(i))->bufferRef();
	    for (int k = 0; k < w; k++)
		{
		    *(buff_out++) = (T)0;
		    buff_in[1]++; buff_in[0]++;
		}
	    for (int j = 0; j < h - 1; j++)
		{
		    for (int k = 0; k < w; k++)
			*(buff_out++) = *(buff_in[1]++) - *(buff_in[0]++);
		}
	}
    return(r);
}



template<class T>
VISImage<T> VISImage<T>::dy(unsigned int order) const
{
    VISImage<T> local_dy_kernel;
    assignImage((float)2*dy_kernel(order), local_dy_kernel); 
    return((convolve(local_dy_kernel)).divAssign((T)2)); 
}


template<class T>
VISImage<T> VISImage<T>::derivative(unsigned int order_x, unsigned int order_y) const
{
    VISImage<T> local_dx_kernel;
    assignImage((float)2*dx_kernel(order_x), local_dx_kernel);

    VISImage<T> local_dy_kernel;
    assignImage((float)2*dy_kernel(order_y), local_dy_kernel);

    return(((convolve(local_dx_kernel)).convolve(local_dy_kernel))
	   .divAssign((T)4));
}


template<class T>
VISImage<T> VISImage<T>::derivative(unsigned int order_x, unsigned int order_y, float scale) const
{
    VISImage<T> local_dx_kernel;
    assignImage(gauss_dx_kernel(order_x, scale), local_dx_kernel);

    VISImage<T> local_dy_kernel;
    assignImage(gauss_dy_kernel(order_y, scale), local_dy_kernel);

// debugging
//    printf("local dx\n");
//    local_dx_kernel.printData();
//    printf("\n");
//
//    printf("local dy\n");
//    local_dy_kernel.printData();
//    printf("\n");

    return(((convolve(local_dx_kernel)).convolve(local_dy_kernel)));
}

// does Gaussian blurring using diffusion to handle boundary conditions
//
template<class T>
VISImage<T> VISImage<T>::gaussDiffuse(float sigma) const
{
  VISImage<T> r(*this);
  float t = sigma*sigma/2.0f;
  int i = (int)ceil(::VISmax(4.0f, t/0.2f)), j;
  float dt = t/i;
  for (j = 0; j < i; j++)
    r += dt*r.diffuse();
  return(r);
}

//
// this returns one iteration of an anisotropic diffusion    
// it computes the derivatives at smallest scale, and weights
// the factor "k" by the rms of the grad mag
//

template<class T>
VISImage<float> VISImage<T>::diffuse() const
{
    VISImage<float> im(*this);
    VISImage<float> dxx, dyy;
    
    dxx = im.dx(2);
    dyy = im.dy(2);

    int i, h = im.height(), w = im.width();
      
    for (i = 0; i < h; i++)
      {
	dxx.poke(0, i) = im.peek(1, i) - im.peek(0, i);
	dxx.poke(w - 1, i) = im.peek(w - 2, i) - im.peek(w - 1, i);
      }

    for (i = 0; i < w; i++)
      {
	dyy.poke(i, 0) = im.peek(i, 1) - im.peek(i, 0);
	dyy.poke(i, h - 1) = im.peek(i, h - 2) - im.peek(i, h - 1);
      }

    return(dxx + dyy);
}

//
// this returns one iteration of an anisotropic diffusion    
// it computes the derivatives at smallest scale, and weights
// the factor "k" by the rms of the grad mag
//

template<class T>
VISImage<float> VISImage<T>::anisoDiffuse(float k) const
{
    VISImage<float> im(*this);
    
    VISImage<float> f_x, f_y, 
	dx_half, dy_half, 
	dx, dy, 
	x_kernel, y_kernel;

    f_x = im.dx();
    f_y = im.dy();
    dx_half = im.dxHalfForward();
    dy_half = im.dyHalfForward();

    float grad_mag_average = (f_x.power(2) + f_y.power(2)).average();
    float k_adj = k*k*grad_mag_average;

    x_kernel = VISImage<T>(3, 1);
    x_kernel = (T)0.5f;
    x_kernel.at(0, 0) = 0.0f;

    y_kernel = VISImage<T>(1, 3);
    y_kernel = (T)0.5f;
    y_kernel.at(0, 0) = 0.0f;

    f_x = f_x.convolve(x_kernel);
    f_y = f_y.convolve(y_kernel);
	    
    dx = dx_half*(((dx_half.power(2) + f_y.power(2))
		   /(-1.0f*k_adj)).exp());
    dy = dy_half*(((dy_half.power(2) + f_x.power(2))
		   /(-1.0f*k_adj)).exp());

    dx = dx.dxHalfBack();
    dy = dy.dyHalfBack();

    return(dx + dy);
}

//
// this returns one iteration of an anisotropic diffusion    
// it computes the derivatives at smallest scale, and weights
// the factor "k" by the rms of the grad mag
//

template<class T>
VISImage<float> VISImage<T>::anisoDiffuse(VISImage<float> image_dx, 
					    VISImage<float> image_dy) const
{
    VISImage<float> im(*this);

    if ((!compareSize(image_dx))&&(!compareSize(image_dx)))
	{
	    im = 0.0f;
	    return(im);
	}
    
    VISImage<float> f_x, f_y, 
	dx_half, dy_half, 
	dx, dy, 
	x_kernel, y_kernel;

    f_x = im.dx();
    f_y = im.dy();
    dx_half = im.dxHalfForward();
    dy_half = im.dyHalfForward();

    dx = dx_half*image_dx;
    dy = dy_half*image_dy;

    dx = dx.dxHalfBack();
    dy = dy.dyHalfBack();

    // set the boundary conditions
    int w = width(), h = height();
    int i;
    
    for (i = 0; i < h; i++)
	{
	    dx.at(0, i) = 0.0f;
	    dx.at(w - 1, i) = 0.0f;
	}

    for (i = 0; i < w; i++)
	{
	    dy.at(i, 0) = 0.0f;
	    dy.at(i, h - 1) = 0.0f;
	}

    return(dx + dy);
}

//
// this returns the median within the window given
//

template<class T>
VISImage<T> VISImage<T>::median(int window_size) const
{
    if (window_size%2 == 0)
	window_size++;

    VISImage<T> median = createToSize();

    int i, j, k, l;
    int index_x, index_y;
    int half_window = window_size/2;
    VISArray<T> array;

    for (i = 0; i < height(); i++)
	for (j = 0; j < width(); j++)
	    {
		for (k = -half_window; k < half_window + 1; k++)
		    for (l = -half_window; l < half_window + 1; l++)
		    {
			index_x = j + l;
			index_y = i + k;
			if (checkBounds((unsigned)index_x, (unsigned)index_y))
			    array.appendItem(itemAt(index_x, index_y));
		    }
		array.sort();
		median.at(j, i) = array.itemAt(array.n()/2);
		array.clear();
	    }
    return(median);
}

template<class T>
VISImage<T> VISImage<T>::gauss(float sigma) const
{

// if sigma <= 0.0 this is a signal for don't blur
    if (sigma <= 0.0) return(*this);

// this needs to be optimized for the float case
    if (type() != VISIm::FLOAT)
	{
	    float scale_x;
	    float scale_y;
	    VISImageType<T> type;
	    VISImage<T> local_col_kernel;
	    VISImage<float> col_kernel;
	    col_kernel = gauss_col_kernel(sigma);
	    scale_y = type.rangeMax()/(col_kernel.max()*
				       col_kernel.height());
	    assignImage(scale_y*col_kernel, local_col_kernel);

	    VISImage<T> local_row_kernel;
	    VISImage<float> row_kernel;
	    row_kernel = gauss_row_kernel(sigma);
	    scale_x = type.rangeMax()/
		(row_kernel.max()*col_kernel.width());
	    assignImage(scale_x*row_kernel, local_row_kernel);
	    
	    return(
		((convolve(local_col_kernel))
		 .convolve(local_row_kernel)).scale(1.0/(scale_x*scale_y))
		);
	}
    else
	{
	    VISImage<T> local_col_kernel;
	    VISImage<float> col_kernel;
	    col_kernel = gauss_col_kernel(sigma);
	    assignImage(col_kernel, local_col_kernel);

	    VISImage<T> local_row_kernel;
	    VISImage<float> row_kernel;
	    row_kernel = gauss_row_kernel(sigma);
	    assignImage(row_kernel, local_row_kernel);
	    
	    return(((convolve(local_col_kernel))
		 .convolve(local_row_kernel)));
	}
}


template< class T >
VISImage<T> VISImage<T>::sqrt() const
{
    VISImage<T> r = this->createToSize();
    for (int i = 0; i < _channels; i++)    
	r.repRef(i)->sqrt(rep(i));
    return(r);
}

template< class T >
VISImage<T> VISImage<T>::exp() const
{
    VISImage<T> r = this->createToSize();
    for (int i = 0; i < _channels; i++)    
	r.repRef(i)->exp(rep(i));
    return(r);
}


template< class T >
VISImage<T> VISImage<T>::ln() const
{
    VISImage<T> r = this->createToSize();
    for (int i = 0; i < _channels; i++)    
	r.repRef(i)->ln(rep(i));
    return(r);
}

template< class T >
VISImage<T> VISImage<T>::abs() const
{
    VISImage<T> r = this->createToSize();
    for (int i = 0; i < _channels; i++)    
	r.repRef(i)->abs(rep(i));
    return(r);
}

template< class T >
VISImage<T> VISImage<T>::sign() const
{
    VISImage<T> r = this->createToSize();
    for (int i = 0; i < _channels; i++)    
	r.repRef(i)->sign(rep(i));
    return(r);
}

template< class T >
VISImage<T> VISImage<T>::pos() const
{
    VISImage<T> r = this->createToSize();
    for (int i = 0; i < _channels; i++)    
	r.repRef(i)->pos(rep(i));
    return(r);
}


template< class T >
VISImage<T> VISImage<T>::neg() const
{
    VISImage<T> r = this->createToSize();
    for (int i = 0; i < _channels; i++)    
	r.repRef(i)->neg(rep(i));
    return(r);
}

template< class T >
VISImage<T> VISImage<T>::scale(float value) const
{
    VISImage<T> r = this->createToSize();
    for (int i = 0; i < _channels; i++)    
	r.repRef(i)->scale(rep(i), value);
    return(r);
}

template< class T>
void VISImage<T>::copy_on_write(VISRep*& r, unsigned ch)
{
    if (r->ref_count() >  1)
	putRep(new VISImageRep<T>(*((VISImageRep<T>*)r)), ch);
}

template< class T >
VISImage<T> VISImage<T>::dx() const 
{return(dx(1));}

template< class T >
VISImage<T> VISImage<T>::dy() const {return(dy(1));}

template< class T >
T VISImage<T>::interp(float x, float y, unsigned int ch) const
{
    unsigned int x_lo, x_hi, y_lo, y_hi;
    float x_in, y_in;
    x_lo = (unsigned int)floor(x);
    y_lo = (unsigned int)floor(y);
    x_hi = (unsigned int)ceil(x);
    y_hi = (unsigned int)ceil(y); 

    x_in = x - (float)x_lo;
    y_in = y - (float)y_lo;
//    printf("y_in %3.2f x_in %3.2f items %3.2f %3.2f %3.2f %3.2f\n", 
//	   y_in, x_in,
//	   rep(ch)->itemAt(x_hi, y_lo), 
//	   rep(ch)->itemAt(x_lo, y_lo), 
//	   rep(ch)->itemAt(x_hi, y_hi),
//	   rep(ch)->itemAt(x_lo, y_hi));

    if ((x < 0.0)||(y < 0.0)||(x > width())||(y > height()))
	{
	    printf("interp out of bounds %f %f \n", x, y);
	    return((T)0);
	}
	   
    return((T)(((float)1.0 - y_in)*
	       (x_in*rep(ch)->itemAt(x_hi, y_lo) 
		+ ((float)1.0 - x_in)*rep(ch)->itemAt(x_lo, y_lo))
	       + y_in*(x_in*rep(ch)->itemAt(x_hi, y_hi) 
		       + ((float)1.0 - x_in)*rep(ch)->itemAt(x_lo, y_hi))));
}


template< class T >
T VISImage<T>::interpNoBounds(float x, float y, unsigned int ch) const
{
    unsigned int x_lo, x_hi, y_lo, y_hi;
    float x_in, y_in;
    x_lo = (unsigned int)floor(x);
    y_lo = (unsigned int)floor(y);
    x_hi = (unsigned int)ceil(x);
    y_hi = (unsigned int)ceil(y); 

    x_in = x - (float)x_lo;
    y_in = y - (float)y_lo;
//    printf("y_in %3.2f x_in %3.2f items %3.2f %3.2f %3.2f %3.2f\n", 
//	   y_in, x_in,
//	   rep(ch)->itemAt(x_hi, y_lo), 
//	   rep(ch)->itemAt(x_lo, y_lo), 
//	   rep(ch)->itemAt(x_hi, y_hi),
//	   rep(ch)->itemAt(x_lo, y_hi));

    return((T)(((float)1.0 - y_in)*
	       (x_in*rep(ch)->itemAt(x_hi, y_lo) 
		+ ((float)1.0 - x_in)*rep(ch)->itemAt(x_lo, y_lo))
	       + y_in*(x_in*rep(ch)->itemAt(x_hi, y_hi) 
		       + ((float)1.0 - x_in)*rep(ch)->itemAt(x_lo, y_hi))));
}


template < class T > VISImRep<T>::~VISImRep()
{
    if ((_buffer != NULL)&&(_free_buffer))
	{
	    delete[] _buffer;
#ifdef __memory_count
	    __MEMORY_SIZE -= sizeof(T)*_size;
#endif
	}
    _buffer = NULL;
}

template < class T > 
void VISImRep<T>::initialize(unsigned int size, T* buffer_in)
{
    _size = size;
    if (buffer_in == NULL)
	if (_size > 0)
	    {
		_buffer = new T[_size];
#ifdef __memory_count
		__MEMORY_SIZE += sizeof(T)*_size;
#endif
		_free_buffer = TRUE;
	    }
	else 
	    {
		_buffer = NULL;
	    }
    else
	{
	    _buffer = buffer_in;
	    _free_buffer = FALSE;
	}
}

template < class T > 
void VISImRep<T>::add(const VISImRep<T>* a, const VISImRep<T>* b)
{
 if (compareSize(a)&&compareSize(b))
     for (int i = 0; i < _size; i++)
	  _buffer[i] = a->_buffer[i] + b->_buffer[i];
 else
     ERROR("VISImRep: add - image sizes not compatible");

}


template< class T >
void VISImRep<T>::sqrt(const VISImRep<T>* a)
{
 if (compareSize(a))
    for (int i = 0; i < _size; i++)
	_buffer[i] = (T)(::sqrt((float)a->_buffer[i]));
 else
     ERROR("VISImRep: sqrt - image sizes not compatible");
}


template< class T >
void VISImRep<T>::exp(const VISImRep<T>* a)
{
 if (compareSize(a))
    for (int i = 0; i < _size; i++)
	_buffer[i] = (T)(::exp((float)a->_buffer[i]));
 else
     ERROR("VISImRep: exp - image sizes not compatible");
}

template< class T >
void VISImRep<T>::sign(const VISImRep<T>* a)
{
    if (compareSize(a))
	for (int i = 0; i < _size; i++)
	    {
		if (a->_buffer[i] != (T)0)
		    _buffer[i] = a->_buffer[i]/(::tabs(a->_buffer[i]));
		else
		    _buffer[i] = 1.0f;
	    }
    else
	ERROR("VISImRep: exp - image sizes not compatible");
}

template< class T >
void VISImRep<T>::ln(const VISImRep<T>* a)
{
 if (compareSize(a))
    for (int i = 0; i < _size; i++)
	_buffer[i] = (T)(::log((double)a->_buffer[i]));
 else
     ERROR("VISImRep: ln - image sizes not compatible");
}

template< class T >
void VISImRep<T>::pos(const VISImRep<T>* a)
{
 if (compareSize(a))
    for (int i = 0; i < _size; i++)
	_buffer[i] = ::VISmax(a->_buffer[i], (T)0);
 else
     ERROR("VISImRep: pos - image sizes not compatible");
}

template< class T >
void VISImRep<T>::abs(const VISImRep<T>* a)
{
 if (compareSize(a))
    for (int i = 0; i < _size; i++)
	_buffer[i] = ::tabs(a->_buffer[i]);
 else
     ERROR("VISImRep: abs - image sizes not compatible");
}

template< class T >
void VISImRep<T>::neg(const VISImRep<T>* a)
{
 if (compareSize(a))
    for (int i = 0; i < _size; i++)
	_buffer[i] = min(a->_buffer[i], (T)0);
 else
     ERROR("VISImRep: neg - image sizes not compatible");
}

template < class T > 
void VISImRep<T>::scale(const VISImRep<T>* a, float value)
{
    if (compareSize(a))
	for (int i = 0; i < _size; i++)
	    _buffer[i] = (T)((float)a->_buffer[i]*value);
    else
	ERROR("VISImRep: scale - image sizes not compatible");
}



template< class T >
void VISImRep<T>::power(const VISImRep<T>* other, int power1)
{
    int i;
    if (power1 <= 0)
	clear((T)1);
    else
	{
	    copy(other);
// new code
	    for (i = 1; i < power1; i++)
		mult(this, other);
/* Updated on April 18, 2000
	    for (i = 1; i < power1; i *= 2)
		mult(this, this);
	    for (; i < power1; i++)		    
		mult(this, other);
*/
	}
}


template < class T > 
void VISImRep<T>::add(const VISImRep<T>* a, const T& value)
{
     for (int i = 0; i < _size; i++)
	  _buffer[i] = a->_buffer[i] + value;
}

template < class T > 
void VISImRep<T>::evaluate(T (*f)(T), const VISImRep<T>* a)
{
     for (int i = 0; i < _size; i++)
	 _buffer[i] = f(a->_buffer[i]);
}


template < class T > 
void VISImRep<T>::gt(const VISImRep<T>* a, T value)
{
     for (int i = 0; i < _size; i++)
	 _buffer[i] = (T)(a->_buffer[i] > value);
}

template < class T > 
void VISImRep<T>::lt(const VISImRep<T>* a, T value)
{
     for (int i = 0; i < _size; i++)
	 _buffer[i] = (T)(a->_buffer[i] < value);
}

template < class T > 
void VISImRep<T>::gteq(const VISImRep<T>* a, T value)
{
     for (int i = 0; i < _size; i++)
	 _buffer[i] = (T)(a->_buffer[i] >= value);
}

template < class T > 
void VISImRep<T>::eq(const VISImRep<T>* a, T value)
{
     for (int i = 0; i < _size; i++)
	 _buffer[i] = (T)(a->_buffer[i] == value);
}

template < class T > 
void VISImRep<T>::lteq(const VISImRep<T>* a, T value)
{
     for (int i = 0; i < _size; i++)
	 _buffer[i] = (T)(a->_buffer[i] <= value);
}


template < class T > 
void VISImRep<T>::mult(const VISImRep<T>* a, const VISImRep<T>* b)
{
 if (compareSize(a)&&compareSize(b))
     for (int i = 0; i < _size; i++)
	  _buffer[i] = a->_buffer[i] * b->_buffer[i];
 else
     ERROR("VISImRep: mult - image sizes not compatible");
}


template < class T > 
void VISImRep<T>::mult(const VISImRep<T>* a, const T& value)
{
     for (int i = 0; i < _size; i++)
	  _buffer[i] = a->_buffer[i]*value;
}

template < class T > 
void VISImRep<T>::div(const VISImRep<T>* a, const VISImRep<T>* b)
{
 if (compareSize(a)&&compareSize(b))
     for (int i = 0; i < _size; i++)
	 _buffer[i] = a->_buffer[i]/b->_buffer[i];
 else
     ERROR("VISImRep: div - image sizes not compatible");
}

template < class T > 
void VISImRep<T>::div(const VISImRep<T>* a, const VISImRep<T>* b, 
		       T zeroCondition)
{
 if (compareSize(a)&&compareSize(b))
     for (int i = 0; i < _size; i++)
	 if (b->_buffer[i] != (T)0)
	     _buffer[i] = a->_buffer[i]/b->_buffer[i];
	else 
	    _buffer[i] = zeroCondition;
 else
     ERROR("VISImRep: div - image sizes not compatible");
}

template < class T > 
void VISImRep<T>::sub(const VISImRep<T>* a, const VISImRep<T>* b)
{
 if (compareSize(a)&&compareSize(b))
     for (int i = 0; i < _size; i++)
	  _buffer[i] = a->_buffer[i] - b->_buffer[i];
 else
     ERROR("VISImRep: sub - image sizes not compatible");

}

template < class T > 
void VISImRep<T>::sub(const VISImRep<T>* a, const T& value)
{
     for (int i = 0; i < _size; i++)
	  _buffer[i] = a->_buffer[i] - value;
}


template < class T > 
void VISImRep<T>::div(const T& value, const VISImRep<T>* a)
{
     for (int i = 0; i < _size; i++)
	  _buffer[i] = value/a->_buffer[i];
}

template < class T > 
void VISImRep<T>::div(const T& value, const VISImRep<T>* a, T zeroCondition)
{
    for (int i = 0; i < _size; i++)
	if (a->_buffer[i] != (T)0)
	    _buffer[i] = value/a->_buffer[i];
	else 
	    _buffer[i] = zeroCondition;
}


template < class T > 
void VISImRep<T>::div(const VISImRep<T>* a, const T& value)
{
     for (int i = 0; i < _size; i++)
	 _buffer[i] = a->_buffer[i]/value;
}


template < class T > 
void VISImRep<T>::sub(const T& value, const VISImRep<T>* a)
{
     for (int i = 0; i < _size; i++)
	  _buffer[i] = value - a->_buffer[i];
}


template < class T > 
void VISImRep<T>::clear(T value)
{
 for (int i = 0; i < _size; i++)
     _buffer[i] = value;
}


template < class T >
void VISImRep<T>::copy(const VISImRep<T>* a)
{
 if (compareSize(a))
     for (int i = 0; i < _size; i++)
	  _buffer[i] = (T)a->_buffer[i];
 else
     ERROR("VISImRep: clear - image sizes not compatible");
}


template < class T >
void VISImRep<T>::copyBuffer(const T* buf)
{
    for (int i = 0; i < _size; i++)
	_buffer[i] = buf[i];
}

template< class T >
void VISImageRep<T>::insetRep(const VISImageRep<T> *in, unsigned int x_pos, 
	      unsigned int y_pos)
{
    if (!(((x_pos + in->width()) <= width())
	  &&((y_pos + in->height()) <= height())))
	{
	    WARN("Image rep inset destination out of bounds\n");
	}
    else
	{	    
	    const T* buf_in = in->buffer();
	    T* buf_out = bufferRef();
	    const T* buf_in_tmp;
	    T* buf_out_tmp;
	    for (int i = y_pos; i < (y_pos + in->height()); i++)
		{
		    buf_out_tmp = buf_out + i*width() + x_pos;
		    buf_in_tmp = buf_in + (i - y_pos)*in->width();
		    for (int j = 0; j < in->width(); j++)
			*(buf_out_tmp++) = *(buf_in_tmp++);
		}
	}
}

template< class T >
void VISImageRep<T>::getROI(const VISImageRep<T> *in, unsigned int x_pos, 
			    unsigned int y_pos)
{
    int w = width(), h = height();
    int in_w = in->width(), 
	in_h = in->height(); 	
    if (!(((x_pos + w) <= in_w)
	  &&((y_pos + h) <= in_h)))
	{
	    WARN("Image rep inset destination out of bounds\n");
	}
    else
	{	    
	    const T* buf_in = in->buffer();
	    T* buf_out = bufferRef();
	    const T* buf_in_tmp;
	    T* buf_out_tmp;
	    for (int i = y_pos; i < (y_pos + h); i++)
		{
		    buf_out_tmp = buf_out + (i - y_pos)*w;
		    buf_in_tmp = buf_in + i*in_w + x_pos;
		    for (int j = 0; j < w; j++)
			*(buf_out_tmp++) = *(buf_in_tmp++);
		}
	}
}


/***************************************************/

template <class T>
T VISImRep<T>::max() const
{
    T max_tmp, t_tmp;
    if (_size > 0)
	max_tmp = _buffer[0];
    else max_tmp = (T)0;
    T* buf = _buffer;
    
    for (int i = 0; i < _size; i++)
	{
	    t_tmp = *(buf++);
	    if (t_tmp > max_tmp)
		max_tmp = t_tmp;
	}

    return(max_tmp);
}

template <class T>
T VISImRep<T>::min() const
{
    T min_tmp, t_tmp;
    if (_size > 0)
	min_tmp = _buffer[0];
    else min_tmp = (T)0;
    T* buf = _buffer;
    
    for (int i = 0; i < _size; i++)
	{
	    t_tmp = *(buf++);
	    if (t_tmp < min_tmp)
		min_tmp = t_tmp;
	}

    return(min_tmp);

}

template <class T>
void VISImRep<T>::max(const VISImRep<T>* a, const VISImRep<T>* b)
{
    const T* buf_a = a->buffer();
    const T* buf_b = b->buffer();
    T* buf_this = _buffer;

    for (int i = 0; i < _size; i++)
	*(buf_this++) = ::VISmax(*(buf_a++), *(buf_b++));
}

template <class T>
void VISImRep<T>::min(const VISImRep<T>* a, const VISImRep<T>* b)
{
    const T* buf_a = a->buffer();
    const T* buf_b = b->buffer();
    T* buf_this = _buffer;

    for (int i = 0; i < _size; i++)
	*(buf_this++) = ::VISmin(*(buf_a++), *(buf_b++));
}

template <class T>
float VISImRep<T>::sum() const
{
    T* buf = _buffer;
    float total = 0.0;
    for (int i = 0; i < _size; i++)
	total += (float)(*(buf++));
    return(total);
}

template < class T >
void VISImage<T>::putChannel(VISImage<T>& other, unsigned int ch)
{
    putRep(other._rep[DEFAULT_CHANNEL], ch);
}

template < class T >
void VISImage<T>::putChannel(VISImage<T>& other)
{
    putRep(other._rep[DEFAULT_CHANNEL], DEFAULT_CHANNEL);
}

// the constructor VISImage<T>(const VISIm&) should take care of these
#ifdef not_for_now
template < class T >
    VISImage<T>::operator VISImage<byte>() const
	{
	    if (type() == BYTE)
		return(*this);
	    else
		{
		    VISImage<byte> r;
		    assignImage(*this, r);
		    return(r);
		}
	}


template < class T >
    VISImage<T>::operator VISImage<unsigned>() const
	{
	    if (type() == UNSIGNED)
		return(*this);
	    else
		{
		    VISImage<unsigned> r;
		    assignImage(*this, r);
		    return(r);
		}
	}

template < class T >
    VISImage<T>::operator VISImage<int>() const
	{
	    if (type() == INT)
		return(*this);
	    else
		{
		    VISImage<int> r;
		    assignImage(*this, r);
		    return(r);
		}
	}

template < class T >
    VISImage<T>::operator VISImage<float>() const
	{
	    if (type() == FLOAT)
		return(*this);
	    else
		{
		    VISImage<float> r;
		    assignImage(*this, r);
		    return(r);
		}
	}

template < class T >
    VISImage<T>::operator VISImage<rgba>() const
	{
	    if (type() == RGBA)
		return(*this);
	    else
		{
		    VISImage<rgba> r;
		    assignImage(*this, r);
		    return(r);
		}
	}
#endif

template < class T >
void VISImage<T>::print() const
{
    printf("VISImage<T> height %d width %d and type %d\n",
	   height(), width(), type());
}

template < class T >
VISImage<T>& VISImage<T>::operator=(T value)
{
    for (int i = 0; i < _channels; i++)
	{
	    repRef(i)->clear(value);
	}
    return(*this);
}

template < class T >
const VISImage<T>& VISImage<T>::evaluate(T (*f)(unsigned int, unsigned int))
{
    for (int i = 0; i < _channels; i++)
	{
	    repRef(i)->evaluate(f);
	}
    return(*this);
}

template < class T >
VISImage<T> VISImage<T>::evaluate(T (*f)(T)) const
{
    VISImage<T> r = createToSize();
    for (int i = 0; i < _channels; i++)
	{
	    r.repRef(i)->VISImRep<T>
		::evaluate(f, (const VISImRep<T>*)(rep(i)));
	}
    return(r);
}

#ifdef USE_HIPS
#include "image.hips.c"
#endif
#ifndef USE_HIPS
template<class T>
VISImage<T> VISImage<T>::mask(const VISImage<T>& kernel) const
{

    unsigned int w = kernel.width();
    unsigned int h = kernel.height();
    unsigned int x_center = w/2;
    unsigned int y_center =  h/2; 

    const T *buffer_in, *buf_in;
    T *buffer_out, *buf_out;
    const T *buf_k;
    T k_value;
    int h_offset, w_offset, w_in;
    w_in = width();

    int i, j, k, l, m;

    if (kernel.channels() != 1)
	{
	    ERROR("VISImage convolve: can only do one channel kernels\n");
	}

    VISImage<T> im_return = this->createToSize();
    im_return = (T)0.0;


    for (i = 0; i < _channels; i++)
	{
	    buf_k = kernel.rep()->buffer();
	    buffer_out = im_return.repRef(i)->bufferRef();
	    buffer_in = rep(i)->buffer();
	    for (j = 0; j < h; j++)
		for (k = 0; k < w; k++)
		    {
			h_offset = j - y_center;
			w_offset = k - x_center;
			k_value = *buf_k++;
			buf_out = buffer_out;
			buf_in = buffer_in + h_offset*w_in + w_offset;

			if (h_offset < 0)
			    {
				buf_out += -(h_offset*w_in);
				buf_in += -(h_offset*w_in);
			    }
			for (l = ::VISmax(0, h_offset);
			     l < ::VISmin(height(), 
				     height() + h_offset); l++)
			    {
				if (w_offset < 0)
				    {
					buf_out += -(w_offset);
					buf_in += -(w_offset);
				    }
				for (m = ::VISmax(0, w_offset); 
				     m < ::VISmin(width(), 
					     width() + w_offset); m++)
				    *buf_out++ += k_value*(*buf_in++);
				if (w_offset > 0)
				    {
					buf_out += (w_offset);
					buf_in += (w_offset);
				    }
			    }
		    }
	}
    return(im_return);
}
#endif


template<class T>
VISImage<T> VISImage<T>::convolve(const VISImage<T>& kernel) const
{

    unsigned int w = kernel.width();
    unsigned int h = kernel.height();
    unsigned int x_center = w/2;
    unsigned int y_center =  h/2; 

    const T *buffer_in, *buf_in;
    T *buffer_out, *buf_out;
    const T *buf_k;
    T k_value;
    int h_offset, w_offset, w_in;
    w_in = width();

    int i, j, k, l, m;

    if (kernel.channels() != 1)
	{
	    ERROR("VISImage convolve: can only do one channel kernels\n");
	}

    VISImage<T> im_return = this->createToSize();
    im_return = (T)0.0;


    for (i = 0; i < _channels; i++)
	{
	    buf_k = kernel.rep()->buffer() + w*h;
	    buffer_out = im_return.repRef(i)->bufferRef();
	    buffer_in = rep(i)->buffer();
	    for (j = 0; j < h; j++)
		for (k = 0; k < w; k++)
		    {
//			h_offset = j - y_center;
//			w_offset = k - x_center;
// this has to reflect the actual definition of convolution
// c(x)  = /int f(x) g(x - y) dy
			h_offset = y_center - j;
			w_offset = x_center - k;
			k_value = *(--buf_k);
			buf_out = buffer_out;
			buf_in = buffer_in + h_offset*w_in + w_offset;

			if (h_offset < 0)
			    {
				buf_out += -(h_offset*w_in);
				buf_in += -(h_offset*w_in);
			    }
			for (l = ::VISmax(0, h_offset);
			     l < ::VISmin(height(), 
				     height() + h_offset); l++)
			    {
				if (w_offset < 0)
				    {
					buf_out += -(w_offset);
					buf_in += -(w_offset);
				    }
				for (m = ::VISmax(0, w_offset); 
				     m < ::VISmin(width(), 
					     width() + w_offset); m++)
				    *buf_out++ += k_value*(*buf_in++);
				if (w_offset > 0)
				    {
					buf_out += (w_offset);
					buf_in += (w_offset);
				    }
			    }
		    }
	}
    return(im_return);
}




template< class T >
VISImage<T> VISImage<T>::noiseShot(float percent, T low, T high) const
{
    VISImage<T> r = *this;
    float* buf;
    int size = r.height()*r.width();
    for (int j = 0; j < channels(); j++)
	{
	    buf = r.repRef(j)->bufferRef();
	    for (int i = 0; i < size; i++)
		{
		    if (rand1() < percent)
			*buf = (T)(rand1()*(float)(high - low) + (float)low);
		    buf++;
		}
	}
    return(r);
}


// shot noise where the false shots are taken from a gaussian deviate
template< class T >
VISImage<T> VISImage<T>::noiseShot(float percent, float sigma) const
{
    long value;
    VISImage<T> r = *this;
    float* buf;
    int size = r.height()*r.width();
    for (int j = 0; j < channels(); j++)
	{
	    buf = r.repRef(j)->bufferRef();
	    for (int i = 0; i < size; i++)
		{
		    if (rand1() < percent)
			*buf += (T)(gasdev()*sigma);
		    buf++;
		}
	}
    return(r);
}


template< class T >
VISImage<T> VISImage<T>::getROI(unsigned int x_pos, 
		       unsigned int y_pos, 
		       unsigned int w_roi, 
		       unsigned int h_roi)
{
    VISImage<T> r(w_roi, h_roi, channels());
    for (int i = 0; i < channels(); i++)
	(r.repRef(i))->getROI(rep(i), x_pos, y_pos);
    return(r);
}

template< class T >
void VISImage<T>::putROI(const VISImage<T>& image_in, 
			 unsigned int x_pos, 
			 unsigned int y_pos)
{
    if (channels() == image_in.channels())
	for (int i = 0; i < channels(); i++)
	    repRef(i)->insetRep(image_in.rep(i), x_pos, y_pos);
}


template< class T >
VISImage<T> VISImage<T>::becomeFlat() const
{

    int new_height, new_width, num_col_tiles, num_row_tiles;
    float sqrt_depth;
    sqrt_depth = ::sqrt((double)channels());
    VISImage<T> image_tmp = VISImage<T>(width(), height());
    image_tmp = (T)0;
    if (channels() == 1) return(*this);

    num_col_tiles = (int)ceil(sqrt_depth);
    num_row_tiles = (int)ceil((float)channels()/
			      (float)num_col_tiles);

    new_width = num_col_tiles*width();
    new_height = num_row_tiles*height();
    
    VISImage<T> image_out(new_width, new_height, 1);

    for (int j = 0; j < num_row_tiles; j++)
	for (int i = 0; i < num_col_tiles; i++)
	{
	    if ((j*num_col_tiles + i) < channels())
		(image_out.repRef())->insetRep(rep(j*num_col_tiles + i),
			 i*width(),
			 j*height());
	    else
		(image_out.repRef())->insetRep(image_tmp.rep(),
					      i*width(),
					      j*height());
	}
    return(image_out);
} 


template< class T >
float VISImage<T>::maskFloat(const VISImage<float>& mask, 
			       unsigned int center_x, 
			       unsigned int center_y, 
			       unsigned int x_pos, 
			       unsigned int y_pos
			       ) const
{

    int i, j;
    float total = 0.0;
    const float *buf_m = mask.rep()->buffer() + (mask.width()*mask.height() - 1);

    int y_lo  = (y_pos - ((mask.height() - 1) - center_y));
    int y_hi  = (y_pos + ((mask.height() - 1) - center_y));

    int x_lo  = (x_pos - ((mask.width() - 1) - center_x));
    int x_hi  = (x_pos + ((mask.width() - 1) - center_x));
    
    int x_pos_tmp, y_pos_tmp;
    
    if ((x_lo >= 0)&&(x_hi < width())&&(y_lo >= 0)&&(y_hi < height()))
	{
//	    const T *buf_this = rep()->buffer() + (y_pos - center_y)*height() 
//		+ (x_pos - center_x);
	    const T *buf_this = rep()->buffer() + (y_lo)*width()
                + (x_lo);

	    for (i = 0; i < mask.height(); i++)
		{
		    for (j = 0; j < mask.width(); j++)	
			{    
//			    printf("at %d %d mask %f, image %f, total %f\n",
//                                   j, i, 
//                                   *buf_this,
//                                   *buf_m,
//                                   total);
			    total += *(buf_this++)*(*(buf_m--));
			}
		    buf_this += width() - mask.width();
		}
	}
    else
	{
	    for (i = 0; i < mask.height(); i++)
		{
		    for (j = 0; j < mask.width(); j++)	
			{
			    x_pos_tmp = x_pos - center_x + j;
			    y_pos_tmp = y_pos - center_y + i;
			    if ((x_pos_tmp >= 0)&&(x_pos_tmp < width())&&
				(y_pos_tmp >= 0)&&(y_pos_tmp < height()))
				{
				    total += itemAt(x_pos_tmp, y_pos_tmp)
					*mask.itemAt((mask.height() - 1) - j, 
						     (mask.width() - 1) - i);
//				    printf("at %d %d mask %f, image %f, total %f\n",
//					   x_pos_tmp, y_pos_tmp,
//					   itemAt(x_pos_tmp, y_pos_tmp), 
//					   mask.itemAt(j, i), 
//					   total);
				}
			}
		}
	}
    return(total);
}



template< class T >
VISArray<float>* VISImage<T>::maskFloat(const VISArray<VISImage<float> >& masks, 
			       unsigned int x_pos, 
			       unsigned int y_pos
			       ) const
{
    VISArray<float>* r = new VISArray<float>();
    for (int i = 0; i < masks.n(); i++)
	r->appendItem(maskFloat(masks.itemAt(i), (masks.itemAt(i)).width()/2, 
				(masks.itemAt(i)).height()/2, x_pos, y_pos));
    return(r);
}


// changed this to do four-connected neighborhood - Ross 24.7.95

template< class T >
VISImage<T> VISImage<T>::zeroCrossings() const
{
    VISImage<T> r = VISImage<T>(this->createToSize());
    r = (T)0;
    const T* buffers[2];
    T* buffers_out[2];
    const T* center;
    T* center_out;
    T this_one, that, abs_this_one, abs_that;
    const T* buff_in = rep()->buffer();
    T* buff_out = r.repRef()->bufferRef();
    unsigned int w = width(), h = height();
    int i, j, k;
 
    /* do first row here specially */
    center = buff_in + 1;
    center_out = buff_out + 1;
    buffers[0] = buff_in;
    buffers_out[0] = buff_out;
    for (j = 1; j < w; j++)
	{
	    this_one = *center;
	    for (k = 0; k < 1; k++)
		{	
		    that = *(buffers[k]++);
				/* use exclusive or here */	
		    if (((this_one < (T)0)&&(that > (T)0))||
			((this_one > (T)0)&&(that < (T)0))||
			((this_one == (T)0)&&(that != (T)0))||
			    ((this_one != (T)0)&&(that == (T)0)))
			{
			    abs_this_one 
				= ::tabs(this_one);
			    abs_that 
				= ::tabs(that);
			    if (abs_this_one < abs_that)
				*center_out = (T)1;
			    else
				*buffers_out[k] = (T)1;
			}
		    buffers_out[k]++;
		}
	    center++;
	    center_out++;
	}

    center = buff_in + w + 1;
    buffers[0] = buff_in + 1;
    buffers[1] = buff_in + w;

    center_out = buff_out + w + 1;
    buffers_out[0] = buff_out + 1;
    buffers_out[1] = buff_out + w;

    for (i = 1; i < h; i++)
	{
	    that = *(buffers[0] - 1);
	    this_one = *buffers[1];
	    if (((this_one < (T)0)&&(that > (T)0))
		||((this_one > (T)0)&&(that < (T)0))
		||((this_one == (T)0)&&(that != (T)0))
		||((this_one != (T)0)&&(that == (T)0))
		)
		{
		    abs_this_one = ::tabs(this_one);
		    abs_that = ::tabs(that);

		    if (abs_this_one <= abs_that)
			*(buffers_out[0] - 1) = (T)1;
		    else
			*buffers_out[1] = (T)1;
		}
	    for (j = 1; j < w; j++)
		{
		    this_one = *center;
		    abs_this_one = ::tabs(this_one);
		    
		    for (k = 0; k < 2; k++)
			{
			    that = *(buffers[k]++);
			    if (((this_one < (T)0)&&(that > (T)0))
				||((this_one > (T)0)&&(that < (T)0))
				||((this_one == (T)0)&&(that != (T)0))
				||((this_one != (T)0)&&(that == (T)0))
				)
				{
				    abs_that = ::tabs(that);
				    if (abs_this_one <= abs_that)
					*center_out = (T)1;
				    else
					*buffers_out[k] = (T)1;
				}
			    (buffers_out[k]++);
			}
		    center++;
		    center_out++;
		}
	    for (k = 0; k < 2; k++)		
	   	{
		    buffers_out[k]++;
		    buffers[k]++;
		}
	    center++;
	    center_out++;
	}

    return(r);
}


template< class T >
VISImage<float> VISImage<T>::distanceTrans() const
{return(distanceTrans((float)(width()+height())));}

template< class T >
VISImage<float> VISImage<T>::distanceTrans(float max_distance) const
{
    int num_neighbors = 4;
    float new_x_dist, new_y_dist, new_dist;
    int w = width(), h = height();
    VISImage<int> x_pos, y_pos, visited;
    x_pos = VISImage<int>(w, h);
    y_pos = VISImage<int>(w, h);
    visited = VISImage<int>(w, h);

    max_distance = max_distance*max_distance; 

    VISImage<float> distance;
    distance = VISImage<float>(w, h);
//    distance = (float)-1.0;
    int N[4][2];
    
    N[0][0] = 0;
    N[0][1] = 1;
    N[1][0] = 0;
    N[1][1] = -1;
    N[2][0] = -1;
    N[2][1] = 0;
    N[3][0] = 1;
    N[3][1] = 0;
    
    int i, j, k;
    int x, y;
    int changes;

    for (i = 0; i < h; i++)
	for (j = 0; j < w; j++)
	    {
		if ((int)itemAt(j, i))
		    {
			distance.at(j, i) = (float)0.0;
			x_pos.at(j, i) = j;
			y_pos.at(j, i) = i;
			visited.at(j, i)  = 1;
		    }
		else
		    {
			distance.at(j, i) = (float)(w + h + 1);
			x_pos.at(j, i) = 0;
			y_pos.at(j, i) = 0;
			visited.at(j, i)  = 0;
		    }
	    }

    int count = 0;
    changes = 1;
    while (changes)
	{
	    changes = 0;
	    for (i = 0; i < (h); i++)
		for (j = 0; j < (w); j++)
		    {
			for (k = 0; k < num_neighbors; k++)
			    {
				x = j + N[k][0];
				y = i + N[k][1];
				if (((x >= 0)&&(x < w)&&(y >= 0)&&(y < h))
				    &&(visited.itemAt(x, y)))
				    {
					new_x_dist 
					    = (float)x_pos.at(x, y) - (float)j;
					new_y_dist 
					    = (float)y_pos.at(x, y) - (float)i;
					new_dist = new_x_dist*new_x_dist + 
					    new_y_dist*new_y_dist;
					if (new_dist <= max_distance)
					    {
					if (!visited.at(j, i))
					    {
						distance.at(j , i) = new_dist;
						x_pos.at(j, i) = x_pos.itemAt(x, y);
						y_pos.at(j, i) = y_pos.itemAt(x, y);
						visited.at(j, i) = 1;
						changes = 1;
					    }
					else if (new_dist < distance.at(j, i))
					    {
						distance.at(j , i) = new_dist;
						x_pos.at(j, i) = x_pos.itemAt(x, y);
						y_pos.at(j, i) = y_pos.itemAt(x, y);
						changes = 1;
					    }
					    }
				    }
			    }
		    }
	}
    return(distance.sqrt());
}


template< class T >
VISImage<float> VISImage<T>::cityBlockDistTrans() const
{
    float new_x_dist, new_y_dist, new_dist;
    int w = width(), h = height();
    VISImage<int> x_pos, y_pos, visited;
    x_pos = VISImage<int>(w, h);
    y_pos = VISImage<int>(w, h);
    visited = VISImage<int>(w, h);
    VISImage<float> distance = VISImage<float>(w, h);

    int i, j;

    for (i = 0; i < h; i++)
	for (j = 0; j < w; j++)
	    {
		if ((int)itemAt(j, i))
		    {
			distance.at(j, i) = (float)0.0;
			x_pos.at(j, i) = j;
			y_pos.at(j, i) = i;
			visited.at(j, i)  = 1;
		    }
		else
		    {
			distance.at(j, i) = (float)((w + h + 1)*(w + h + 1));
			x_pos.at(j, i) = 0;
			y_pos.at(j, i) = 0;
			visited.at(j, i)  = 0;
		    }
	    }

    for (i = 0; i < (h); i++)
	for (j = 1; j < (w); j++)
	    {
		if (visited.itemAt(j - 1, i))
		    {
			new_x_dist 
			    = (float)x_pos.at(j - 1, i) - (float)(j);
			new_y_dist 
			    = (float)y_pos.at(j - 1, i) - (float)(i);
			new_dist = new_x_dist*new_x_dist + 
			    new_y_dist*new_y_dist;
			if (new_dist < distance.at(j, i))
			    {
				distance.at(j , i) = new_dist;
				x_pos.at(j, i) = x_pos.itemAt(j - 1, i);
				y_pos.at(j, i) = y_pos.itemAt(j - 1, i);
				visited.at(j, i) = 1;
			    }
		    }
	    }

    for (i = 0; i < (h); i++)
	for (j = (w - 2); j >= 0; j--)
	    {
		if (visited.itemAt(j + 1, i))
		    {
			new_x_dist 
			    = (float)x_pos.at(j + 1, i) - (float)(j);
			new_y_dist 
			    = (float)y_pos.at(j + 1, i) - (float)(i);
			new_dist = new_x_dist*new_x_dist + 
			    new_y_dist*new_y_dist;
			
			if (new_dist < distance.at(j, i))
			    {
				distance.at(j , i) = new_dist;
				x_pos.at(j, i) = x_pos.itemAt(j + 1, i);
				y_pos.at(j, i) = y_pos.itemAt(j + 1, i);
				visited.at(j, i) = 1;
			    }
		    }
	    }

    for (j = 0; j < w; j++)
	for (i = 1; i < (h); i++)
	    {
		if (visited.itemAt(j, i - 1))
		    {
			new_x_dist 
			    = (float)x_pos.at(j, i - 1) - (float)(j);
			new_y_dist 
			    = (float)y_pos.at(j, i - 1) - (float)(i);
			new_dist = new_x_dist*new_x_dist + 
			    new_y_dist*new_y_dist;
			
			if (new_dist < distance.at(j, i))
			    {
				distance.at(j , i) = new_dist;
				x_pos.at(j, i) = x_pos.itemAt(j, i - 1);
				y_pos.at(j, i) = y_pos.itemAt(j, i - 1);
				visited.at(j, i) = 1;
			    }
		    }
	    }

    for (j = 0; j < w; j++)
	for (i = h - 2; i >= 0; i--)
	    {
		if (visited.itemAt(j, i + 1))
		    {
			new_x_dist 
			    = (float)x_pos.at(j, i + 1) - (float)(j);
			new_y_dist 
			    = (float)y_pos.at(j, i + 1) - (float)(i);
			new_dist = new_x_dist*new_x_dist + 
			    new_y_dist*new_y_dist;
			
			if (new_dist < distance.at(j, i))
			    {
				distance.at(j , i) = new_dist;
				x_pos.at(j, i) = x_pos.itemAt(j, i + 1);
				y_pos.at(j, i) = y_pos.itemAt(j, i + 1);
				visited.at(j, i) = 1;
			    }
		    }
	    }
    return(distance.sqrt());
}




template< class T >
VISImage<T>  VISImage<T>::reduce(int scale) const
{
    VISImage<T> r = VISImage<T>(width()/scale, height()/scale, channels());
    int i, j, k;
    int x_mod, y_mod;
    int x_offset, y_offset;
    const T* buf_in;
    T* buf_out;


    if (scale <= 0)
	return(*this);

    x_mod = width()%scale;
    y_mod = height()%scale;

    x_offset = (scale - 1)/2 + (int)ceil((double)x_mod/2.0);
    y_offset = (scale - 1)/2 + (int)ceil((double)y_mod/2.0);

#ifdef DEBUG
    printf("Image::reduce x_mod %d x_offset %d, y_mod %d y_offset %d\n", 
	   x_mod, x_offset, y_mod, y_offset);
#endif
    
    int position = 0;
    for (i = 0; i < channels(); i++)
	{
	    buf_in = rep(i)->buffer() + width()*(y_offset) 
		+ x_offset;
	    position = width()*(y_offset)+ x_offset;
	    buf_out = r.repRef(i)->bufferRef();
	    for (j = 0; j < r.height(); j++)
		{
		    for (k = 0; k < r.width(); k++)
			{
			    (*buf_out++) = *buf_in;
			    buf_in += scale;
			    position += scale;
			}
		    buf_in += x_mod + (scale - 1)*width();
		    position += x_mod + (scale - 1)*width();
		}
	}
    return(r);
}

// this averages blocks of "scale x scale" pixels and returns 
// the appropriately smaller image
// Note: scale must be greater than zero

template< class T >
VISImage<T>  VISImage<T>::reduceAverage(int scale) const
{
    VISImage<T> r = VISImage<T>(width()/scale, height()/scale, channels());
    r = 0.0f;
    int i, j, k, l, m;
    int x_mod, y_mod;
    int x_offset, y_offset;

    if (scale <= 0)
	return(*this);

    x_mod = width()%scale;
    y_mod = height()%scale;


//
// try to keep things centered
//
    x_offset = (int)ceil((double)x_mod/2.0);
    y_offset = (int)ceil((double)y_mod/2.0);

    int w = r.width(), h = r.height();
    int x, y;
    
    int position = 0;
    for (i = 0; i < channels(); i++)
	{
	    for (j = 0; j < h; j++)
		{
		    y = j*scale + y_offset;
		    for (k = 0; k < w; k++)
			{
			    x = k*scale + x_offset;
			    for (l = 0; l < scale; l++)
				for (m = 0; m < scale; m++)
				    r.at(k, j, i) 
					+= itemAt(x + l, y + m, i);
			}
		}
	}
    return(r/(scale*scale));
}


// this returns a resampled image scaled according to the size
// indicated by "the_scale".  It rescales around the middle of
// the image.

template< class T >
VISImage<T> VISImage<T>::resample(float the_scale) const
{
    return(resample(the_scale, the_scale, width()/2.0 - 0.5, 
		    height()/2.0 - 0.5));
}

// this returns a anistropically resampled image scaled according to the size
// indicated by "scale_x" and "scale_y".  It rescales around the point
// given by "x" and "y".

template< class T >
VISImage<T> VISImage<T>::resample(float scale_x, float scale_y, 
				  float x, float y) const
{
    unsigned new_w, new_h;
    float new_x, new_y;
    float tmp_x, tmp_y;
    unsigned w = width(), h = height();
    int i, j;
    
    new_w  = (float)w*scale_x;
    new_h  = (float)h*scale_y;


    new_x = (x + 0.5)*scale_x - 0.5;
    new_y = (y + 0.5)*scale_y - 0.5;

    VISImage<T> image_return(new_w, new_h);
    for (j = 0; j < new_h; j++)
	for (i = 0; i < new_w; i++)
		{
		    tmp_x = ((i - new_x)/scale_x + x); 
		    tmp_y = ((j - new_y)/scale_y + y); 
		    
		    tmp_x = ::VISmax(tmp_x, 0.0f);
		    tmp_x = ::VISmin(tmp_x, (float)(w - 1));

		    tmp_y = ::VISmax(tmp_y, 0.0f);
		    tmp_y = ::VISmin(tmp_y, (float)(h - 1));

		    image_return.at(i, j) 
			= interp(tmp_x, tmp_y);
		}
    return(image_return);
}

// this returns a anistropically resampled image of the indicated
// size that has been rescaled about the middle


template< class T >
VISImage<T> VISImage<T>::resample(unsigned w, unsigned h) const
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
VISImage<int>  VISImage<T>::cannyEdges() const
{
    unsigned int w = width();
    unsigned int h = height();
    VISImage<T> derivatives[dIndex(3,0)];
    VISImage<int> r;
    VISImage<T> grad_mag;
    VISImage<T> curve;
    VISImage<T> zero_crossings, curvature_pos, grad_mag_pos;    
    int i, j;
    for (i = 0; i <= 2; i++)
	for (j = 0; j <= (2 - i); j++)
		derivatives[dIndex(j, i)] = derivative(j, i);

    grad_mag = (((derivatives[dIndex(1, 0)])*(derivatives[dIndex(1, 0)]))
		+ ((derivatives[dIndex(0, 1)])
		   *(derivatives[dIndex(0, 1)]))).sqrt();
    derivatives[dIndex(1, 0)] = derivatives[dIndex(1, 0)].div(grad_mag, 0.0f);
    derivatives[dIndex(0, 1)] = derivatives[dIndex(0, 1)].div(grad_mag, 0.0f);

    curve = ((derivatives[dIndex(1, 0)]^2)*derivatives[dIndex(2, 0)] 
	     + (derivatives[dIndex(0, 1)]^2)*derivatives[dIndex(0, 2)] 
	     + 2.0f*derivatives[dIndex(0, 1)]*derivatives[dIndex(1, 0)]
	     *derivatives[dIndex(1, 1)]);
    
    zero_crossings = curve.zeroCrossings();
    curvature_pos = (((curve.dx()*derivatives[dIndex(1, 0)]) + 
			   (curve.dy()*derivatives[dIndex(0, 1)]))
// changed this to <= to handle special cases Ross 2-5-97
			   <= (T)0.0);
    grad_mag_pos = (grad_mag > (T)0.0);

    r = (VISImage<int>)(zero_crossings*curvature_pos*grad_mag_pos);
//    r = (VISImage<int>)((curve.zeroCrossings())
//			*(((curve.dx()*derivatives[dIndex(1, 0)]) + 
//			   (curve.dy()*derivatives[dIndex(0, 1)]))
//			   < (T)0.0)*(grad_mag > (T)0.0));
    

    return(r);
}


template< class T >
VISImage<int>  VISImage<T>::cannyEdges(T threshold) const
{
    unsigned int w = width();
    unsigned int h = height();
    VISImage<T> derivatives[dIndex(3,0)];
    VISImage<int> r;
    VISImage<T> grad_mag, grad_mag_tmp;
    VISImage<T> curve;
    VISImage<T> zero_crossings, curvature_pos, grad_mag_pos;    
    int i, j;
    for (i = 0; i <= 2; i++)
	for (j = 0; j <= (2 - i); j++)
		derivatives[dIndex(j, i)] = derivative(j, i);

    grad_mag = (((derivatives[dIndex(1, 0)])*(derivatives[dIndex(1, 0)]))
		+ ((derivatives[dIndex(0, 1)])
		   *(derivatives[dIndex(0, 1)]))).sqrt();

// NOTE: if you don't check for zeros, the fp process really drags.
// This might be hardware dependent.
grad_mag_tmp = VISImage<T>(grad_mag.width(), grad_mag.height());
  for (int m = 0; m < grad_mag_tmp.height(); m++)
     for (int n = 0; n < grad_mag_tmp.width(); n++)
   if (grad_mag.itemAt(n, m) == (T)0.0)
      grad_mag_tmp.at(n, m) = (T)1.0;
   else
      grad_mag_tmp.at(n, m) = (T)1.0;
    
    derivatives[dIndex(1, 0)] /= grad_mag_tmp;
    derivatives[dIndex(0, 1)] /= grad_mag_tmp;
    curve = ((derivatives[dIndex(1, 0)]^2)*derivatives[dIndex(2, 0)] 
	     + (derivatives[dIndex(0, 1)]^2)*derivatives[dIndex(0, 2)] 
	     + derivatives[dIndex(0, 1)]*derivatives[dIndex(1, 0)]
	     *derivatives[dIndex(1, 1)]);
    zero_crossings = curve.zeroCrossings();
    curvature_pos = (((curve.dx()*derivatives[dIndex(1, 0)]) + 
			   (curve.dy()*derivatives[dIndex(0, 1)]))
			   < (T)0.0);
    grad_mag_pos = (grad_mag > threshold);
    r = (VISImage<int>)(zero_crossings*curvature_pos*grad_mag_pos);


    return(r);
}

template< class T >
VISImage<T> VISImage<T>::floodFill(T thresh, unsigned x, unsigned y)
{
    VISImIndexList list;
    VISImage<T> r = createToSize();

    // these are the neighbors
    int N[4][2];
    N[0][0] = 1;
    N[0][1] = 0;
    N[1][0] = -1;
    N[1][1] = 0;
    N[2][0] = 0;
    N[2][1] = 1;
    N[3][0] = 0;
    N[3][1] = -1;

    r = (T)0;
    unsigned at_x, at_y, next_x, next_y;
    int i;

    if (itemAt(x, y) > thresh)
	return(r);

//  if it passes the threshold put the starting pixel on the list
    
    list.appendItem(VISImIndex(x, y));
//  mark that pixel as "in"
    r.at(x, y) = (T)1.0;
    list.reset();
    
//    int n = 0;
    while (list.valid())
	{
	    at_x = list.atCurrent().a();
	    at_y = list.atCurrent().b();
//	    printf("flood fill iteration %d x %d y %d z %d\n", n++, at_x, 
//		   at_y );
	    
// look at your neighbors and if they are: 1) in bounds, 2) more than the
// threshod and 3) not yet visited, put them on the list and mark them as
// visited.     

	    for (i = 0; i < 4; i++)
		{
		    next_x = at_x + N[i][0];
		    next_y = at_y + N[i][1];
		    if (checkBounds((unsigned)next_x, (unsigned)next_y))
			{
			    if ((r.itemAt(next_x, next_y) == (T)0.0)
				&&(itemAt(next_x, next_y) <= thresh))
			    {
				    r.at(next_x, next_y) = (T)1.0;
				    list.appendItem(VISImIndex(next_x, 
								next_y));
//				    printf("append item %d x %d y %d z %d\n", 
//					  next_x, next_y);
				}
			}
		}
//
// remove the guy whose neighbors you have just visited
// when the list is empty, you are done.
//
	    list.removeCurrent();
	}
    return(r);
}

// this one replaces one value with another
template< class T >
const VISImage<T>& VISImage<T>::floodFill(T label_from, T label_to, 
				     unsigned x, unsigned y)
{
    VISImIndexList list;

    // these are the neighbors
    int N[4][2];
    N[0][0] = 1;
    N[0][1] = 0;
    N[1][0] = -1;
    N[1][1] = 0;
    N[2][0] = 0;
    N[2][1] = 1;
    N[3][0] = 0;
    N[3][1] = -1;

    unsigned at_x, at_y, next_x, next_y;
    int i;

    if ((itemAt(x, y) != label_from)
	||(label_from == label_to))
	return(*this);

//  if it passes the threshold put the starting pixel on the list
    
    list.appendItem(VISImIndex(x, y));
//  mark that pixel as "in"
    //cout << "inside ff";exit(1)
    at(x, y) = label_to;
    list.reset();
    
//    int n = 0;
    while (list.valid())
	{
	    at_x = list.atCurrent().a();
	    at_y = list.atCurrent().b();
//	    printf("flood fill iteration %d x %d y %d\n", n++, at_x, 
//		   at_y );
	    
// look at your neighbors and if they are: 1) in bounds, 2) more than the
// threshod and 3) not yet visited, put them on the list and mark them as
// visited.     

	    for (i = 0; i < 4; i++)
		{
		    next_x = at_x + N[i][0];
		    next_y = at_y + N[i][1];
		    if (checkBounds((unsigned)next_x, (unsigned)next_y))
			{
			    if (itemAt(next_x, next_y) == label_from)
				{
				    at(next_x, next_y) = label_to;
				    list.appendItem(VISImIndex(next_x, 
								next_y));
//				    printf("append item %d x %d y %d z %d\n", 
//					   next_x, next_y);
				}
			}
		}
//
// remove the guy whose neighbors you have just visited
// when the list is empty, you are done.
//
	    list.removeCurrent();
	}
    return(*this);
}


template< class T >
VISImage<T> VISImage<T>::floodFillBoundBox(T thresh, unsigned x, unsigned y,
					   unsigned &x_lo, unsigned &y_lo, 
					   unsigned &x_hi, unsigned &y_hi)
{
    VISImIndexList list;
    VISImage<T> r = createToSize();

    // these are the neighbors
    int N[4][2];
    N[0][0] = 1;
    N[0][1] = 0;
    N[1][0] = -1;
    N[1][1] = 0;
    N[2][0] = 0;
    N[2][1] = 1;
    N[3][0] = 0;
    N[3][1] = -1;

    r = (T)0;
    unsigned at_x, at_y, next_x, next_y;
    int i;

    x_hi = x_lo = x;
    y_hi = y_lo = y;
    

//    printf("x lo %d hi %d y lo %d hi %d\n",
//	   x_lo, x_hi, y_lo, y_hi);

//    printf("first new item at x %d y %d with at %f\n", 
//	   x, y, 
//	   itemAt(x, y));

    if (itemAt(x, y) > thresh)
	return(r);


//  if it passes the threshold put the starting pixel on the list
    
    list.appendItem(VISImIndex(x, y));
    x_hi = x_lo + 1;
    y_hi = y_lo + 1;

//  mark that pixel as "in"
    r.at(x, y) = (T)1.0;
    list.reset();
    
    while (list.valid())
	{
	    at_x = list.atCurrent().a();
	    at_y = list.atCurrent().b();
//	    printf("flood fill iteration %d x %d y %d z %d\n", n++, at_x, 
//		   at_y );
	    
// look at your neighbors and if they are: 1) in bounds, 2) more than the
// threshod and 3) not yet visited, put them on the list and mark them as
// visited.     
	    for (i = 0; i < 4; i++)
		{
		    next_x = at_x + N[i][0];
		    next_y = at_y + N[i][1];
		    if (checkBounds((unsigned)next_x, (unsigned)next_y))
			{
			    if ((r.itemAt(next_x, next_y) == (T)0.0)
				&&(itemAt(next_x, next_y) <= thresh))
			    {
				    r.at(next_x, next_y) = (T)1.0;
				    list.appendItem(VISImIndex(next_x, next_y));
				    x_lo = ::VISmin(x_lo, next_x);
				    x_hi = ::VISmax(x_hi, (next_x + 1));
				    y_lo = ::VISmin(y_lo, next_y);
				    y_hi = ::VISmax(y_hi, (next_y + 1));
			    }
			}
		}
//
// remove the guy whose neighbors you have just visited
// when the list is empty, you are done.
//
	    list.removeCurrent();
	}
    return(r);
}



// this routine calculates derivative up to order "degree" and
// then puts them into a multi channel image 
// the channel of a particular deritivative is 
// given by dIndex(dx, dy).
// If you give it a multichannel image, then it's up to you 
// to figure out which image is which.

template< class T >
VISImage<T> VISImage<T>::derivatives(unsigned degree) const
{
    VISImage<T> r;
    VISImage<T> *derivs;
    derivs = new VISImage<T>[dIndex(degree + 1,0)];
    int i, j;
    
    for (i = 0; i <= degree; i++)
	for (j = 0; j <= (degree - i); j++)
	    derivs[dIndex(j, i)] = derivative(j, i);
    
    r = derivs[0];
    for (i = 1; i < dIndex(degree + 1, 0); i++)
	r = VISImage<T>(r, derivs[i]);
	
    delete[] derivs;
    return(r);
}

template< class T >
VISImage<T> VISImage<T>::derivatives(unsigned degree, float scale) const
{
    VISImage<T> r;
    VISImage<T> *derivs;
    derivs = new VISImage<T>[dIndex(degree + 1,0)];
    int i, j;
    
    for (i = 0; i <= degree; i++)
	for (j = 0; j <= (degree - i); j++)
	    derivs[dIndex(j, i)] = derivative(j, i, scale);
    
    r = derivs[0];
    for (i = 1; i < dIndex(degree + 1, 0); i++)
	r = VISImage<T>(r, derivs[i]);
	
    delete[] derivs;
    return(r);
}


template< class T >
float VISImage<T>::average() const
{
    float total = 0.0;
    for (int i = 0; i < channels(); i++)
	{
	    total += rep(i)->sum();
	}
    return(total/(float)(height()*width()*channels()));
}



template< class T >
void VISImage<T>::printData() const
{
    for (int i = 0; i < height(); i++)
	{
	for (int j = 0; j < width(); j++)
	    printf("%3.2f ", (float)itemAt(j, i));
	printf("\n");
	}
}

template< class T >
VISImage<T> VISImage<T>::setBorder(T value, unsigned thickness) const
{
    unsigned i, j;
    unsigned t = ::VISmin(thickness, height());
    VISImage<T> r(*this);
    for (i = 0; i < t; i++)
	{
	    for (j = 0; j < width(); j++)
		{
		    r.at(j, i) = value;
		    r.at(j, (height() - 1) - i) = value;
		}
	}

    t = ::VISmin(thickness, width());
    for (i = 0; i < t; i++)
	{
	    for (j = 0; j < height(); j++)
		{
		    r.at(i, j) = value;
		    r.at((width() - 1) - i, j) = value;
		}
	}
    return(r);
}


/* these are non-member functions */


template< class T1, class T2 >
void copy(const VISImage<T1>& from, VISImage<T2>& to)
{
  if (compareSize(to, from))
      for (int i = 0; i < from.channels(); i++)
	  copy((from.rep(i)), (to.repRef(i)));
  else
      ERROR("VISImage<>: copy - image sizes not compatible");
}




template< class T1, class T2 >
void assignImage(const VISImage<T1>& from, VISImage<T2>& to)
{
    if (from.type() != to.type())
    {
	VISImageRep<T2>* rep_tmp;
	to.rechannel(from.channels());
	for (int i = 0; i < to.channels(); i++)
	{
	    rep_tmp = new VISImageRep<T2>(from.rep(i)->width(),
					  from.rep(i)->height());
	    copy(from.rep(i), rep_tmp);
	    to.putRep(rep_tmp, i);
	}
    }
    else
	to.assign(from);
}

template< class T1, class T2 >
void assignImage(const VISImage<T1>* from, VISImage<T2>* to)
{
    if (from->type() != to->type())
    {
	VISImageRep<T2>* rep_tmp;
	to->rechannel(from->channels());
	for (int i = 0; i < to->channels(); i++)
	    {
		rep_tmp = new VISImageRep<T2>(from->rep(i)->width(), 
					      from->rep(i)->height());
		copy(from->rep(i), rep_tmp);
		to->putRep(rep_tmp, i);
	    }
    }
    else
	to->assign(*from);
}

template< class T1, class T2 >
int compareSize(const VISImage<T1>& a, const VISImage<T2>& b)
	{
	    return((a.channels() == b.channels())
		   &&compareSize(a.rep(DEFAULT_CHANNEL), 
				 b.rep(DEFAULT_CHANNEL)));
	}

template< class T1, class T2 >
int compareSize(const VISImageRep<T1>* a, const VISImageRep<T2>* b)
	{
		return((a->width() == b->width())
		       &&(a->height() == b->height()));
	}
	

template< class T1, class T2 >
void copy(const VISImageRep<T1>* a, VISImageRep<T2>* b)
{
    if (!(compareSize(a, b)))
	printf("ERROR copy\n");
    else
	copy((VISImRep<T1>*)a, (VISImRep<T2>*) b);
}

template< class T, class T2 >    
int compareSize(const VISImRep<T>* a, const VISImRep<T2>* b)
{
    return(a->size() == b->size());
}

template< class T, class T2 >    
void copy(const VISImRep<T2>* a, VISImRep<T>* b)
{
    const T2 *buf_from;
    T *buf_to;
    buf_to = b->bufferRef();
    buf_from = a->buffer();
    if (compareSize(a, b))
	for (int i = 0; i < a->size(); i++)
	    buf_to[i] = (T)buf_from[i];
    else
	ERROR("copy(const VISImRep<T>& a, VISImRep<T2>& b) - image sizes not compatible");
}



template< class T >
VISImage<int> VISImage<T>::watershed(float thresh, 
				       float depth, int size,
				       int pad)
{
    int 		h = height(), w = width(),x,y;
    int                 label_max, max_array_size, min_array_size,true_label;
    int 	 	go, reg_min,reg_max;
    int			i, j, k, index,tmp_size;		    
    float		diff,min_value,max_val,real_label;
    int 		N[24][2], x_pos, y_pos, flat, u, v;
    int 	  	new_label,count,min_found,min_x,min_y;
    float		test;
    VISArray<int>          max_array,order;
    VISArray<VISImIndex> 	array,min_border,change,tmp_set, boundary_low;
    VISImIndex		coord, pixel_loc,current,location;
    VISImage<int>      label(w,h), bin(w,h),marker(w,h),final_bin(w,h);
    VISImageFile	label_file;
    T                  data_max, data_min;
    

// ---------------------------
// Begin Watershed Alogorithm.
// ---------------------------

   data_max = max();
   data_min = min();

// ----------------------------------------------------
// Threshold low and high values to a reasonable value.
// ----------------------------------------------------

   label = -1;

   for (i = pad; i < (w - pad); i++)
   {
    for (j = pad; j < (h - pad); j++)
    {
	if (itemAt(i,j) <= thresh)
	{
		at(i,j) = data_min;
	}
      }
    }

   *this = setBorder(data_max + (T)1, pad);
   label = label.setBorder(-2, pad);

// -----------------------------------
// Establish the eight-neighbor array.
// -----------------------------------


  N[8][0] = -2;    N[9][0] = -1;   N[10][0] = 0;    N[11][0] = 1;    N[12][0] = 2;
  N[8][1] = -2;    N[9][1] = -2;   N[10][1] = -2;   N[11][1] = -2;   N[12][1] = -2;
                 
  N[23][0] = -2;   N[4][0] = -1;   N[1][0] = 0;     N[5][0] = 1;     N[13][0] = 2;
  N[23][1] = -1;   N[4][1] = -1;   N[1][1] = -1;    N[5][1] = -1;    N[13][1] = -1;

  N[22][0] = -2;   N[0][0] = -1;                    N[2][0] = 1;     N[14][0] = 2;
  N[22][1] = 0;    N[0][1] = 0;                     N[2][1] = 0;     N[14][1] = 0;
 
  N[21][0] = -2;   N[7][0] = -1;   N[3][0] = 0;     N[6][0] = 1;     N[15][0] = 2;
  N[21][1] = 1;    N[7][1] = 1;    N[3][1] = 1;     N[6][1] = 1;     N[15][1] = -1;

  N[20][0] = -2;   N[19][0] = -1;  N[18][0] = 0;    N[17][0] = 1;    N[16][0] = 2;
  N[20][1] = 2;    N[19][1] = 2;   N[18][1] = 2;    N[17][1] = 2;    N[16][1] = 2;

// ###########################################################
       
// -----------------------------------------
// Loop through the 4 neighbors and put all
// neighboring pixels which are equal to the
// center pixel on a list.
// ----------------------------------------- 

	new_label = 0;
        marker = 0;

 for (i = pad; i < (w - pad); i++)
 {
  for (j = pad; j < (h - pad); j++)
  {
	flat = 0;
	go = 1;
	reg_min = 0;
	
// -------------------------------------------
// Check to see if the pixel has already been
// visited and labeled.  If it has, set a flag
// so that the next pixel will be visited 
// directly with no operations done on this 
// pixel.
// -------------------------------------------

	if (label.itemAt(i,j) > 0)
	{
		go = 0;
	}

// ----------------------------------------------------------
// Check to see if the pixel you are at is a member of a
// flat region.  If it is not, then proceed and do not label
// that pixel.  If it is, and has a neighboring pixel with
// a label, then label the current pixel with its neighbors
// label and exit.  If the pixel is a member of a flat 
// region, but none of its neighbors are labelled, a new
// flat region has been found.  Note that and go on to label
// a new flat region with a new label.
// ----------------------------------------------------------

        else
	{
  	 for (k = 0; k < 4; k++)
	 {
		x_pos = i + N[k][0];
		y_pos = j + N[k][1];

		if (itemAt(x_pos,y_pos) == itemAt(i,j))
		{ 
		    flat = 1;
		}

		if (itemAt(x_pos,y_pos) > itemAt(i,j))
		{
			    reg_min++;
		}			
	 }
	}

// ----------------------------------------------------------
// If you are at a single pixel regional minimum, 
// increment the label counter, set the marker, and continue.
// ----------------------------------------------------------
        
	if (reg_min == 4)
	{
		new_label++;
		label.at(i,j) = new_label;
		marker.at(i,j) = 1;
	}

//-----------------------------------------------------
// If you are at a new flat region, increment the label
// counter and continue.
// ----------------------------------------------------

	if (flat == 1)
	{
		new_label++;
		label.at(i,j) = new_label;
		array.appendItem(VISImIndex(i,j));
		tmp_set.appendItem(VISImIndex(i,j));
		marker.at(i,j) = 1;
		min_value = itemAt(i,j);
		max_val = min_value;
		reg_max = 1;
		min_found = 0;
	}

// --------------------------------------------
// If the pixel is a member of a flat region
// and has not been labeled yet, locate all 
// the pixels in the flat region, and mark them
// as such.  Also check to see if the flat
// region is a regional minimum, regional max-
// imum, or neither and mark them as such.
// --------------------------------------------

	while ((go == 1) && (flat == 1))
	{
	 coord = array.itemAt(array.n() - 1);
	 array.removeItemAt(array.n() - 1);

	 for (k = 0; k < 4; k++)
	 { 	
		x_pos = coord.a() + N[k][0];
		y_pos = coord.b() + N[k][1];

		if ((itemAt(x_pos,y_pos) == itemAt(i,j)) && 
		    (label.itemAt(x_pos,y_pos) == -1))
		{
			label.at(x_pos,y_pos) = new_label;
			array.appendItem(VISImIndex(x_pos,y_pos));
			tmp_set.appendItem(VISImIndex(x_pos,y_pos));
			marker.at(x_pos,y_pos) = 1;
		}
		
		if ((itemAt(x_pos,y_pos) > max_val) && 
		    (itemAt(x_pos,y_pos) != (float)(data_max + 1)))
		{
		    reg_max = 0;
		}

		if (itemAt(x_pos,y_pos) < min_value)
		{
		   min_value = itemAt(x_pos,y_pos);
		   min_found = 1;
		   min_x = x_pos;
		   min_y = y_pos;
		   
		}
	 }
	 
// -------------------------------------------
// Once all members of a flat region have been
// labeled, set a marker for those associated
// which areas of regional maxima.
// -------------------------------------------

	 if (array.n() == 0)
	 {
	     if (reg_max == 1)
	     {
		 tmp_size = tmp_set.n();
		 for (k = 0; k < tmp_size; k++)
		 {
		     coord = tmp_set.itemAt(k);
		     marker.at(coord.a(),coord.b()) = 2;
		 }
		 
		 tmp_set.clear();
	     }
	     
// -------------------------------------------
// Mark an array assigning a minimum location
// to those members of a flat region which are
// not basins.
// -------------------------------------------

	      if ((min_found == 1) && (reg_max != 1))
	      {
		     change.insertItemAt(VISImIndex(min_x,min_y),new_label);
		     tmp_set.clear();
	      }
         
              array.clear();
              tmp_set.clear();
	      go = 0;
	 }
	}
  }
 }

//    printf("Flat regions and regional minima marked.\n");

// if (create == 1)
// {
//    remove("data.fit");
//    remove("initial_label.fit");
//    remove("marker.fit");
//    label_file.write_fits(label,"initial_label.fit");
//    label_file.write(*this,"data.fit");
//    label_file.write_fits(marker,"marker.fit");
// }

// ###########################################################
// ###########################################################

// -------------------------------------------------
// Follow the pixel to its closest minimum neighbor
// and label it with the same label as that minimum
// neighbor.
// -------------------------------------------------
  
    int view, button;

    for (i = pad; i < (w - pad); i++)
	{
	    for (j = pad; j < (h - pad); j++)
		{ 
		    view = 1;
		    u = i;
		    v = j;
		    current.a(u);
		    current.b(v);
		    coord.a(0);
		    coord.b(0);
	
		    if (marker.itemAt(i,j) == 0)
			{
			    while (label.itemAt(current.a(),current.b()) <= 0)
				{
				    min_value = itemAt(u,v);
				    for (k = 0; k < 4; k++)
					{
					    x_pos = u + N[k][0];
					    y_pos = v + N[k][1];

					    if (itemAt(x_pos,y_pos) < min_value)
						{
						    coord.a(x_pos);
						    coord.b(y_pos);
						    min_value 
							= itemAt(coord.a(),
								 coord.b());
						}
		
					    if (k == 3)
						{
						    array.appendItem
							(VISImIndex
							 (current.a(),
							  current.b()));
						    current = coord;
						    u = coord.a();
						    v = coord.b();
						}
					}
				}

			    for (k = 0; k < array.n(); k++)
				{
				    location = array.itemAt(k);
				    label.at(location.a(),location.b()) =
					label.itemAt(current.a(),current.b());
				}
			}

		    array.clear();
		}
 }

//    printf("Pixels tracked down to minimum.\n");

// if (create == 1)
// {
//    remove("label.fit");
//    label_file.write_fits(label.setBorder(0,pad),"label.fit");
// }

// ###########################################################
// ###########################################################

// ------------------------------------------------------
// Handle the flat regions which are not regional minima.
// ------------------------------------------------------

// this is very inefficient --- change it for speed RTW 2-15-97
// 
// min_array_size = change.n();
//
// for (k = 0; k < min_array_size; k++)
// {
//     coord = change.itemAt(k);
//     true_label = label.itemAt(coord.a(),coord.b());
//     
//     if ((coord.a() != 0) && (coord.b() != 0))
//     {
//     for (i = pad; i < (w - pad); i++)
//     {
//      for (j = pad; j < (h - pad); j++)
//      { 
//        if (label.itemAt(i,j) == k)
//        { 
//	   label.at(i,j) = true_label;
//        }
//      }
//     }
//     }
// }

 min_array_size = change.n();
 for (i = pad; i < (w - pad); i++)
     {
	 for (j = pad; j < (h - pad); j++)
	     { 
		 if (((k = label.itemAt(i,j)) >= 0)&&
		     (k < min_array_size))
		     { 
			 coord = change.itemAt(k);
			 if ((coord.a() != 0) && (coord.b() != 0))
			     {
				 true_label 
				     = label.itemAt(coord.a(),coord.b());
				 label.at(i,j) = true_label;
			     }
		     }
	     }
     }
 
  change.clear();

//  printf("Flat regions and those associated have been labeled.\n");

// if (create == 1)
// {
//  remove("label_II.fit");
//  label_file.write_fits(label.setBorder(0,pad),"label_II.fit");
// }

// ###########################################################
// ###########################################################

// ------------------------------------------------------
// Handle the flat regions which are regional maxima.
// ------------------------------------------------------

  for (i = pad; i < (w - pad); i++)
  {
   for (j = pad; j < (h - pad); j++)
   { 
    if (marker.itemAt(i,j) == 2)
    {
	 min_value = itemAt(i,j);

	 for (k = 0; k < 24; k++)
	 {
	   x_pos = i + N[k][0];
	   y_pos = j + N[k][1];

           if ((x_pos >= 0) && (y_pos >= 0) && (x_pos <= (w-1)) && (y_pos <= (h-1)))
	   { 
	   
	    if (itemAt(x_pos,y_pos) < min_value)
	    {
	       min_value = itemAt(x_pos,y_pos);
	       coord.a(x_pos);
	       coord.b(y_pos);
	    }
	   }
	 }

	 label.at(i,j) = label.itemAt(coord.a(),coord.b());

    }
   }
  }

//  printf("Flat regions (maxima) and those associated have been labeled.\n");

// if (create == 1)
// {
//  remove("label_III.fit");
//  label_file.write_fits(label.setBorder(0,pad),"label_III.fit");
// }

// ###########################################################
// ###########################################################

// ---------------------------------
// Locate the lowest pixel value in
// any given region.
// ---------------------------------
    

    float ws_depth;
    boolean done_consolidate = FALSE;
    boolean done;
    int number_consolidations;
// shouldn't this be kept from something above rather than recalculated?
    label_max = label.max();
    VISArray<int> size_array;
    VISArray<float> min_array(label_max + 1);
    VISArray<int> change_array(label_max + 1);
    VISArray< VISList<int> > connections(label_max + 1);
    int num_consolidations = 0;
    
    while (!done_consolidate)
	{
	    for (i = 0; i <= label_max ; i++)
		{
		    min_array.at(i) = (float)(data_max + 1.0f);
		    change_array.at(i) = -1;
		    size_array.at(i) = 0;
		    boundary_low.at(i) = VISImIndex(-1, -1);
		}

	    for (i = pad; i < (w - pad); i++)
		{
		    for (j = pad; j < (h - pad); j++)
			{ 
			    index = label.itemAt(i,j);
			    if (index > 0)
				{
				    if (itemAt(i,j) < min_array.itemAt(index))
					{
					    min_array.at(index) = itemAt(i,j);
					}
				    size_array.at(index)++;
				}
			}
		}

//	    printf("Lowest pixel in a region located.\n");

// ###########################################################
// ###########################################################

// ----------------------------------
// Loop through the image and locate
// all the boundary pixels of each 
// watershed region, and mark it with
// a value of -1.
// ----------------------------------
// Also, keep track of low boundary point
// 
// 
//
// ----------------------------------

	    float this_value;
	    VISImIndex low_index;
//	    bin = 0;
	    int this_label;

	    for (i = pad; i < (w - pad); i++)
		{
		    for (j = pad; j < (h - pad); j++)
			{ 
			    count = 0;
			    this_value = itemAt(i, j);
			    for (k = 0; k < 4; k++)
				{
				    x_pos = i + N[k][0];
				    y_pos = j + N[k][1];
	
				    if (itemAt(x_pos, y_pos) < this_value)
					{
					    if ((((this_label 
						  = label.itemAt(x_pos,y_pos))
						 != label.itemAt(i,j)))
						&&(this_label > 0))
						{
						    count = 1;
						    low_index = 
							boundary_low
							.itemAt(this_label);
						    if 
							(checkBounds(
							    (unsigned)
							    low_index.a(), 
							    (unsigned)
							    low_index.b()))
							{
							    if (itemAt
								(low_index.a(),
								 low_index.b())
								> this_value)
								boundary_low.at
								    (this_label)
								    = VISImIndex(i,
										  j);
							}
						    else
							boundary_low
							    .at(this_label)
							    = VISImIndex(i, j);
						}
					} // end debug
				}
			    if (count == 1)
				{
				    bin.at(i,j) = -1;
				}
			    else
				{
				    bin.at(i,j) = 0;
				}
			}
		}

	    bin = bin.setBorder(-2,pad);
	    if (num_consolidations == 0)
		label_file.write_fits(bin,"initial_edges.fit");
//	    printf("bin min %d max %d\n", bin.min(), bin.max());

//	    printf("First pass to located edges of regions.\n");

// if (create == 1)
// {
//   remove("binary_original.fit");
//   label_file.write_fits(bin, "binary_original.fit");
// }

// *************************
// Boundary Pixels have been
// properly marked.
// *************************

// ###########################################################

// --------------------------
// Loop through the image and
// locate shallow watersheds.
// If they are too shallow, 
// puncture the boundary.
// --------------------------

// Change this to punch through only the lowest - Ross 9-9-97
// for (i = pad; i < (w - pad); i++)
// {
//  for (j = pad; j < (h - pad); j++)
//  {
//	if ((bin.itemAt(i,j) == -1))
//	{
//	 index = label.itemAt(i,j);
//
//	 if ((itemAt(i,j) - min_array.itemAt(index)) < depth)	     
//	 {
//	     bin.at(i,j) = 0;
//	 }
//	}
//  }
// } 
 

	    done_consolidate = TRUE;
	    number_consolidations = 0;
	    int x_next, y_next, that_label;
 
	    for (i = 0; i <= label_max; i++)
		{
		    low_index = boundary_low.itemAt(i);
		    if (checkBounds((unsigned)low_index.a(), (unsigned)low_index.b()))
			{
			    ws_depth = itemAt(low_index.a(), low_index.b()) 
				- min_array.itemAt(i);
			    if ((ws_depth < depth)||
				(size_array.itemAt(i) < size))
				{
// if a region is going to go away you must
// 1) get rid of the boundary
// 2) relabel the regions that are adjacent (i.e. , they 
//   have an adjacent nonbin pixel.
				    bin.at(low_index.a(), low_index.b()) = 0;
				    for (k = 0; k < 4; k++)
					{
					    x_next = low_index.a() + N[k][0];
					    y_next = low_index.b() + N[k][1];
					    if ((bin.itemAt(x_next, y_next) 
						 == 0)&&
						((that_label = 
						  label.itemAt(x_next, y_next))
						 != i)&&
						(that_label >= 0))
						{
						    number_consolidations++;
						    (connections.at(i))
							.appendItem
							(that_label);
						    (connections.at
						     (that_label))
							.appendItem(i);
//						    printf("got connection %d and %d x %d y %d min %f and depth %f\n", 
//							   i, that_label, 
//							   x_next, y_next, 
//							   min_array.itemAt(i),
//							   ws_depth);
						}
					}
				}
			}
		}

//	    printf("Second pass to located edges of regions.\n");
	    bin = bin.setBorder(-2,pad);

// if (create == 1)
// {
//   remove("binary.fit");
//   label_file.write_fits(bin, "binary.fit");
// }

//	    printf("got one consolidation pass %d regions %d\n", 
//		   number_consolidations, 
//		   label_max);

	    if (number_consolidations > 0)
		done_consolidate = FALSE;


//	    now convert to the new labels
//          first set up array of new values
	    int num_regions = 0;
	    VISList<int> these_regions;
	    Link<int> *this_link;
	    
	    for (i = 0; i <= label_max; i++)
		{
		    if (change_array.itemAt(i) == -1)
			{
			    change_array.at(i) 
				= num_regions;

			    these_regions.deleteItems();
			    getLists(these_regions, connections, i);

// now set up the change array
			    this_link = these_regions.head();
			    while (this_link != NULL)
				{
				    change_array.at(this_link->data()) 
					= num_regions;
//				    printf("change region %d to %d with lablel %d, num reg %d\n", 
//					   this_link->data(), i, 
//					   change_array.itemAt
//					   (this_link->data()), 
//					   num_regions);
				    this_link = this_link->next();
				}
			    num_regions++;
			}
		}

	    
//	    for (i = 0; i < change_array.n(); i++)
//		printf("change %d to %d \n", i, change_array.itemAt(i));
	    
	    VISImage<int> labels_new(w, h);
	    int old_label;

	    for (i = 0; i < w; i++)
		{
		    for (j = 0; j < h; j++)
			{
			    if ((old_label = label.itemAt(i, j)) >= 0)
				if (old_label <= label_max)
				    labels_new.at(i, j) 
					= change_array.itemAt(old_label);
				else
				    printf("got label too big %d at %d %d\n", 
					   old_label, i, j);
			    else
				labels_new.at(i, j)  = old_label;
			}
		}

	    if (num_regions >= label_max)
		done_consolidate = TRUE;
	    else
		label_max = num_regions - 1;
	    
	    label = labels_new;
	    num_consolidations++;
	    
	} // while not done_consolidate 

    label_file.write_fits(bin,"final_edges.fit");

#ifdef needs_to_be_rewritten

// ###########################################################
 
// ---------------
// Final Labelling
// ---------------

	    new_label = 0;

	    for (i = pad; i < (w - pad); i++)
		{
		    for (j = pad; j < (h - pad); j++)
			{
			    go = 1;

			    if (bin.itemAt(i,j) == 0)
				{
				    new_label++;
				    array.appendItem(VISImIndex(i,j));
				}
		
			    if (bin.itemAt(i,j) == -1)
				{
				    go = 0;
				}
			    
			    if ((bin.itemAt(i,j) == 0) && (go == 1))
				{
				    while (go == 1)
					{
					    coord = array.itemAt(array.n() - 1);
//					    bin.at(coord.a(), coord.b()) = new_label;
					    array.removeItemAt(array.n() - 1);
					    for (k = 0; k < 4; k++)
						{ 	
						    x_pos = coord.a() + N[k][0];
						    y_pos = coord.b() + N[k][1];

						    if (bin.itemAt(x_pos,y_pos) == 0)
							{
							    bin.at(x_pos,y_pos) 
								= new_label;
							    array.appendItem
								(VISImIndex(x_pos,
									     y_pos));
							}
						}
					    if (array.n() == 0)
						{
						    go = 0;
						}
					}
				}
			    array.clear();
			}
		}


//	    printf("Done relabeling\n");
// ###########################################################
// -----------------
// Remove appendages
// -----------------
	    int done = 0, mark = 1;
	    int rand;
	    int kk;
	    
	    while (!done)
		{
		    done = 1;
		    for (i = pad; i < (w - pad); i++)
			{
			    for (j = pad; j < (h - pad); j++)
				{
				    if (bin.itemAt(i,j) < 0)
					{
					    done = 0;
					    rand 
						= (int)rint
						(4.0f*rand1() - 0.5);
					    for (k = 0; k < 4; k++)
						{
						    kk = (k+rand)%4;
						    x_pos = i + N[kk][0];
						    y_pos = j + N[kk][1];
						    
						    if (final_bin
							.itemAt(x_pos,y_pos) 
							> 0)
							{
							    bin.at(i,j) = bin
								.itemAt(x_pos,
									y_pos);
							}
						}
					}
				}
			}
		}

	    label_file.write_fits(bin,"final_final_edges.fit");

	    label = bin;
	    label_max = label.max();
#endif




    return(label);
}

template <class T>
void VISImage<T>::getLists(VISList<int> &these_regions, 
			   VISArray< VISList<int> > &connections, 
			   int starting_region)
{
    Link<int> *this_link;
//    , *next_link;
    int next_region;
    these_regions.appendItem(starting_region);
    this_link = (connections.at(starting_region)).head();
    while (this_link != NULL)
	{
	    next_region = this_link->data();
//	    next_link = this_link->next();
	    (connections.at(starting_region)).removeItem(this_link);
	    getLists(these_regions, connections, next_region);
	    this_link = (connections.at(starting_region)).head();
	}
}

// these create bigger images with repeated copies of the smaller ones
template <class T>
VISImage<T> VISImage<T>::repeatDown(unsigned num) const
{
  
  VISImage<T> im_r = VISImage<T>();
    if (num <= 0)
	return(im_r);
    int w = width(), h = height();
    VISImage<T> r(w, num*h);
    T *buf_out = (r.repRef())->bufferRef();
    int i, j, k;

    const T *buf_in; 

    for (k = 0; k < num; k++)
	{
	    buf_in = (rep())->buffer();
	    for (i = 0; i < w*h; i++)
		*(buf_out++) = *buf_in++;
	}
    return(r);
}

template <class T>
VISImage<T> VISImage<T>::repeatAcross(unsigned num) const
{  
  VISImage<T> im_r = VISImage<T>();
    if (num <= 0)
	return(im_r);
    int w = width(), h = height();
    VISImage<T> r(num*w, h);
    T *buf_out;
    int i, j, k;
    buf_out = (r.repRef())->bufferRef();
    const T *buf_in;
    

    for (j = 0; j < h; j++)
	for (k = 0; k < num; k++)
	    {
		buf_in = (rep())->buffer() + j*w;
		for (i = 0; i < w; i++)
		    *(buf_out++) = *buf_in++;
	    }

    return(r);
}

template< class T >
VISImage<int> VISImage<T>::extrema() const
{
    unsigned int w = width();
    unsigned int h = height();
    VISImage<int> r = VISImage<int>(w, h);
    VISImage<T> derivatives[dIndex(3,0)];
    int i, j;
    for (i = 0; i <= 2; i++)
	for (j = 0; j <= (2 - i); j++)
	    derivatives[dIndex(j, i)] = derivative(j, i);
    
    r = (VISImage<int>)VISImage<T>(
	(derivatives[dIndex(1, 0)].zeroCrossings()*
	 derivatives[dIndex(0, 1)].zeroCrossings())*
	((derivatives[dIndex(2, 0)]
	   *derivatives[dIndex(0, 2)]
	   - ((derivatives[dIndex(1, 1)])^2))
	  >= (T)0));
    return(r);
}


#endif /* IMAGE_C */
