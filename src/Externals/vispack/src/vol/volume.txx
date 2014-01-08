// ******************************************************************
// VISPack. Copyright (c) 1994-2000 Ross Whitaker rtw@utk.edu       *
// For conditions of distribution and use, see the file LICENSE.txt *
// accompanying this distribution.                                  *
// ******************************************************************

// $Id: volume.txx,v 1.6 2005/12/01 03:32:51 miriah Exp $

/* *	sccsid "@(#) image.C     2.0 2/16/94" */

//#define USE_HIPS

// #include "image.h" 
// #include "im.h" 

//
// *******************

//VISImage

// *******************
//

#ifndef iris_volume_c
#define iris_volume_c

#include <math.h>

#include "image/image_type.h"
#include <limits.h>
//#include "vol/volume.h"
#include "vol/volindexlist.h"
#include "util/mathutil.h"

template< class T >
VISVolume<T>::~VISVolume()
{
  // this is done in the base class    
  // unref(_rep);
}

template< class T>
VISVolume<T>::VISVolume(VISVolumeRep<T>* rep)
{
  initialize(rep);
}
 
template< class T >
VISVolume<T>::VISVolume(unsigned int w, unsigned int h, unsigned int d, 
			T* buf)
{
  VISImageType<T> type_tmp;
  _type = type_tmp.whatAmI();
  ref(_rep = (VISRep*)(new VISVolumeRep<T>(w, h, d, buf)));
}

template< class T >
VISVolume<T>::VISVolume(const VISVolume<T>& volume)
{
  initialize(volume.rep());
}


template< class T>
void VISVolume<T>::copy_on_write(VISRep*& r)
{
  VISRep* r_tmp = r;
  if (r->ref_count() >  1)
  {
    ref(r = (VISRep*)(new VISVolumeRep<T>(((VISVolumeRep<T>*)r))));
    unref(r_tmp);
  }
}

template< class T >
void VISVolumeRep<T>::insetRep(const VISVolumeRep<T> *in, 
			       unsigned int x_pos, 
			       unsigned int y_pos, 
			       unsigned int z_pos)
{
  if (!(((x_pos + in->width()) <= width())
	&&((y_pos + in->height()) <= height())
	&&((z_pos + in->depth()) <= depth())))
  {
    WARN("Image rep inset destination out of bounds\n");
  }
  else
  {	    
    const T* buf_in = in->buffer();
    T* buf_out = VISVolumeRep<T>::bufferRef();// **mdm**
    const T* buf_in_tmp;
    T* buf_out_tmp;
    for (int j = z_pos; j < (z_pos + in->depth()); j++)
      for (int i = y_pos; i < (y_pos + in->height()); i++)
      {
	buf_out_tmp = buf_out + (j*height() + i)*width() 
	  + x_pos;
	buf_in_tmp = buf_in + ((j - z_pos)*in->height() 
			       + (i - y_pos))*in->width();
	for (int j = 0; j < in->width(); j++)
	  *(buf_out_tmp++) = *(buf_in_tmp++);
      }
  }
}

template< class T >
void VISVolumeRep<T>::getROI(const VISVolumeRep<T> *in, 
			     unsigned int x_pos, 
			     unsigned int y_pos, 
			     unsigned int z_pos)
{
  if (!(((x_pos + width()) <= in->width())
	&&((y_pos + height()) <= in->height())
	&&((z_pos + depth()) <= in->depth())))
  {
    WARN("Image rep inset destination out of bounds\n");
  }
  else
  {	    
    const T* buf_in = in->buffer();
    T* buf_out = VISVolumeRep<T>::bufferRef();// **mdm**
    const T* buf_in_tmp;
    T* buf_out_tmp;
    for (int j = z_pos; j < (z_pos + depth()); j++)
      for (int i = y_pos; i < (y_pos + height()); i++)
      {
	buf_out_tmp = buf_out + ((j - z_pos)*height()
				 + (i - y_pos))*width();
	buf_in_tmp = buf_in + (j*in->height() + i)
	  *in->width() + x_pos;
	for (int j = 0; j < width(); j++)
	  *(buf_out_tmp++) = *(buf_in_tmp++);
      }
  }
}

template< class T >
void VISVolume<T>::assign(const VISVolume<T>& from)
{
  //  cout << "about to print" << endl;
  //  from.print();
  if (!(&from == this))
  {
    unref(_rep);
    ref(_rep = (VISRep*)(from.rep()));
  }
}

template< class T >
void VISVolume<T>::assign(const VISVol& from)
{
  if (!(&from == this))
  {
    if (_type == from.type())
    {
      unref(_rep);
      ref(_rep = (VISRep*)from.rep());
    }
    else
    {
      VISVolumeRep<float> *rep_float;
      VISVolumeRep<byte> *rep_byte;
      VISVolumeRep<unsigned> *rep_unsigned;
      VISVolumeRep<int> *rep_int;
      VISVolumeRep<rgba> *rep_rgba;
      VISVolumeRep<short> *rep_short;
      VISVolumeRep<unsigned short> *rep_ushort;
		    
      unref(_rep);

      switch(from.type())
      {
      case VISIm::FLOAT:
	rep_float = (VISVolumeRep<float>*)from.rep();
	ref(_rep = new VISVolumeRep<T>
	    (rep_float->width(), rep_float->height(), 
	     rep_float->depth()));
	::copy(rep_float, repRef());
	break;
      case VISIm::BYTE:
	rep_byte = (VISVolumeRep<byte>*)from.rep();
	ref(_rep = new VISVolumeRep<T>
	    (rep_byte->width(), rep_byte->height(), 
	     rep_byte->depth()));
	::copy(rep_byte, repRef());
	break;
      case VISIm::SHORT:
	rep_short = (VISVolumeRep<short>*)from.rep();
	ref(_rep = new VISVolumeRep<T>
	    (rep_short->width(), rep_short->height(), 
	     rep_short->depth()));
	::copy(rep_short, repRef());
	break;
      case VISIm::USHORT:
	rep_ushort = (VISVolumeRep<unsigned short>*)from.rep();
	ref(_rep = new VISVolumeRep<T>
	    (rep_ushort->width(), rep_ushort->height(), 
	     rep_ushort->depth()));
	::copy(rep_ushort, repRef());
	break;
      case VISIm::UNSIGNED:
	rep_unsigned 
	  = (VISVolumeRep<unsigned>*)from.rep();
	ref(_rep = new VISVolumeRep<T>
	    (rep_unsigned->width(), 
	     rep_unsigned->height(), 
	     rep_unsigned->depth()));
	::copy(rep_unsigned, repRef());
	break;
      case VISIm::INT:
	rep_int = (VISVolumeRep<int>*)from.rep();
	ref(_rep = new VISVolumeRep<T>
	    (rep_int->width(), rep_int->height(), 
	     rep_int->depth()));
	::copy(rep_int, repRef());
	break;
      case VISIm::RGBA:
	rep_rgba = (VISVolumeRep<rgba>*)from.rep();
	ref(_rep = new VISVolumeRep<T>
	    (rep_rgba->width(), rep_rgba->height(), 
	     rep_rgba->depth()));
	::copy(rep_rgba, repRef());
	break;
      case VISIm::OTHER:
	ref(_rep = new VISVolumeRep<T>(0,0,0));
	WARN("VISVolume<T>::assign - auto conversion from type other not supported\n");
	break;
      case VISIm::NONE:
	ref(_rep = new VISVolumeRep<T>(0,0,0));		
	WARN("VISVolume<T>::assign - auto conversion from type none not supported\n");
	break;
      }
    }
  }
}



template<class T>
VISVolume<T> VISVolume<T>::addAssign(const T& value)
{
  (this->repRef())->add(this->rep(), value);
  return(*this);
}

template<class T>
VISVolume<T> VISVolume<T>::subAssign(const T& value)
{
  (this->repRef())->sub(this->rep(), value);
  return(*this);
}


template<class T>
VISVolume<T> VISVolume<T>::multAssign(const T& value)
{
  (this->repRef())->mult(this->rep(), value);
  return(*this);
}

template<class T>
VISVolume<T> VISVolume<T>::multAssign(const VISVolume<T>& volume)
{
  (this->repRef())->mult(this->rep(), volume.rep());
  return(*this);
}

template<class T>
VISVolume<T> VISVolume<T>::divAssign(const T& value)
{
  (this->repRef())->div(this->rep(), value);
  return(*this);
}

template<class T>
VISVolume<T> VISVolume<T>::divAssign(const VISVolume<T>& volume)
{
  (this->repRef())->div(this->rep(), volume.rep());
  return(*this);
}

template<class T>
VISVolume<T>  VISVolume<T>::addAssign(const VISVolume<T>& volume)
{
  (this->repRef())->add(this->rep(), volume.rep());
  return(*this);
}


template<class T>
VISVolume<T>  VISVolume<T>::subAssign(const VISVolume<T>& volume)
{
  (this->repRef())->sub(this->rep(), volume.rep());
  return(*this);
}

  
template< class T >
void VISVolume<T>::copy(const VISVolume<T>& a)
{
  if (compareSize(a))
    repRef()->copy(a.rep());
  else
    ERROR("VISVolume<>: copy - volume sizes not compatible");
}


#ifdef this_is_dangerous
// Ross 21.3.95
template<class T>
VISVolume<T>& VISVolume<T>::operator=(const VISVol& from)
{
  if (from.type() == type())
  {
    unref(_rep);
    ref(_rep = (VISRep*)from.rep());
  }
  else
#ifdef AUTO_VOLUME_CONVERSION
    switch(from.type())
    {
    case VISIm::BYTE:
      assignVolume(*((VISVolume<byte>*)(&from)), *this);
      break;
    case VISIm::INT:
      assignVolume((*(VISVolume<int>*)(&from)), *this);
      break;
    case VISIm::FLOAT:
      assignVolume((*(VISVolume<float>*)(&from)), *this);
      break;
    default:
      ERROR("unrecognized volume type conversion\n");
    }
#else
  {
    printf("got unequal type assign\n");
    ERROR("automatic volume conversion must be set as compiler option");
  }
#endif
  return(*this);
}
#endif

template<class T>
VISVolume<T>::VISVolume(const VISVol& from)
{
  VISImageType<T> type_tmp;
  _type = type_tmp.whatAmI();
  assign(from);
}


template<class T>
void VISVolume<T>::setBorder(T value, unsigned w)
{
  unsigned width_tmp;
  int i, j, k;

  width_tmp = MIN(w, depth());
  for (i = 0; i < width_tmp; i++)
    for (j = 0; j < height(); j++)	
      for (k = 0; k < width(); k++)	
      {
	at(k, j, i) = value;
	at(k, j, (depth() - 1) - i) = value;
      }

  width_tmp = MIN(w, height());
  for (i = 0; i < depth(); i++)
    for (j = 0; j < width_tmp; j++)	
      for (k = 0; k < width(); k++)	
      {
	at(k, j, i) = value;
	at(k, (height() - 1) - j, i) = value;
      }

  width_tmp = MIN(w, width());
  for (i = 0; i < depth(); i++)
    for (j = 0; j < height(); j++)	
      for (k = 0; k < width_tmp; k++)	
      {
	at(k, j, i) = value;
	at((width() - 1) - k, j, i) = value;
      }

}

template< class T >
VISVolume<T> VISVolume<T>::operator+(const VISVolume<T>& volume) const
{
  VISVolume<T> volume_return = volume.createToSize();
    
  if (compareSize(volume))
  {
    volume_return.repRef()->add((this->rep()), 
				(volume.rep()));
  }
  else 
  {
    ERROR("VISIm:operator +; volume size mismatch");
    print("vol 1");
    volume.print("vol 2");
  }
  return(volume_return);
}


template< class T >
VISVolume<T> VISVolume<T>::operator*(const VISVolume<T>& volume) const
{
  VISVolume<T> volume_return = volume.createToSize();
    
  if (compareSize(volume))
  {
    volume_return.repRef()->mult((this->rep()), 
				 (volume.rep()));
  }
  else 
  {
    ERROR("VISIm:operator*; volume size mismatch");
  }
  return(volume_return);
}


template< class T >
VISVolume<T> VISVolume<T>::operator-(const VISVolume<T>& volume) const
{
  VISVolume<T> volume_return = volume.createToSize();
    
  if (compareSize(volume))
  {
    volume_return.repRef()->sub((this->rep()), 
				(volume.rep()));
  }
  else 
  {
    ERROR("VISIm:operator -; volume size mismatch");
  }
  return(volume_return);
}


template< class T >
VISVolume<T> VISVolume<T>::operator/(const VISVolume<T>& volume) const
{
  VISVolume<T> volume_return = volume.createToSize();
    
  if (compareSize(volume))
  {
    volume_return.repRef()->div((this->rep()), 
				(volume.rep()));
  }
  else 
  {
    ERROR("VISIm:operator /; volume size mismatch");
  }
  return(volume_return);
}

template< class T >
VISVolume<T> VISVolume<T>::div(const VISVolume<T>& volume, 
			       T value) const
{
  VISVolume<T> volume_return = volume.createToSize();
    
  if (compareSize(volume))
  {
    volume_return.repRef()->div((this->rep()), 
				(volume.rep()), value);
  }
  else 
  {
    ERROR("VISIm:operator /; volume size mismatch");
  }
  return(volume_return);
}



template< class T >
void VISVolume<T>::initialize(unsigned int w, unsigned int h, unsigned int d)
{
  VISImageType<T> type_tmp;
  _type = type_tmp.whatAmI();
  ref(_rep = (VISRep*)(new VISVolumeRep<T>(w, h, d)));
}



template< class T >
void VISVolume<T>::initialize(const VISVolumeRep<T>* rep)
{
  VISImageType<T> type_tmp;
  _type = type_tmp.whatAmI();
  ref(_rep = (VISRep*)rep);
}

template< class T >
T VISVolume<T>::max() const
{
  return(rep()->max());
}

template< class T >
T VISVolume<T>::min() const
{
  return(rep()->min());
}

template< class T >
float VISVolume<T>::sum() const
{
  return(rep()->sum());
}

template< class T >
VISVolume<T> VISVolume<T>::max(const T &value) const
{
  VISVolume<T> volume_return = this->createToSize();
  volume_return.repRef()->max((this->rep()), value);
  return(volume_return);
}

template< class T >
VISVolume<T> VISVolume<T>::min(const T &value) const
{
  VISVolume<T> volume_return = this->createToSize();
  volume_return.repRef()->min((this->rep()), value);
  return(volume_return);
}

template< class T >
VISVolume<T> VISVolume<T>::min(const VISVolume<T>& volume) const
{
    
  VISVolume<T> volume_return = volume.createToSize();
    
  if (compareSize(volume))
  {
    volume_return.repRef()->min((this->rep()), 
				(volume.rep()));
  }
  else 
  {
    ERROR("VISVolume:min; volume size mismatch");
  }
  return(volume_return);
}

template< class T >
VISVolume<T> VISVolume<T>::max(const VISVolume<T>& volume) const
{
  VISVolume<T> volume_return = volume.createToSize();
  if (compareSize(volume))
  {
    volume_return.repRef()->max((this->rep()), 
				(volume.rep()));
  }
  else 
  {
    ERROR("VISVolume:max; volume size mismatch");
  }
  return(volume_return);
}



//
// *******************

//VISRep

// *******************
//

template< class T >
void VISVolumeRep<T>::initialize(unsigned int w, unsigned int h, 
				 unsigned int d, T* buffer_in)
{
  VISImRep<T>::initialize(w*h*d, buffer_in);
  _width = w;
  _height = h;
  _depth = d;
}

template< class T >
void VISVolumeRep<T>::evaluate(T (*f)(unsigned int, unsigned int, 
				      unsigned int))
{
  for (int h = 0; h < _depth; h++)
    for (int i = 0; i < _height; i++)
      for (int j = 0; j < _width; j++)
	//		VISVolumeRep<T>::_buffer[index(j, i, h)] = f(j, i, h);// **mdm**
	VISVolumeRep<T>::_buffer[index(j, i, h)] = f(j, i, h);// **mdm**
}

template< class T >
void VISVolumeRep<T>::evaluate(T (*f)(const T&), const VISVolumeRep<T> *other)
{
  for (int h = 0; h < _depth; h++)
    for (int i = 0; i < _height; i++)
      for (int j = 0; j < _width; j++)
	VISVolumeRep<T>::_buffer[index(j, i, h)] = f(other->_buffer[index(j, i, h)]);// **mdm**
}

template< class T >
void VISVolumeRep<T>::evaluate(T (*f)(unsigned int, unsigned int, 
				      unsigned int), 
			       unsigned int rect_ul_x, 
			       unsigned int rect_ul_y, 
			       unsigned int rect_ul_z, 
			       unsigned int rect_width, 
			       unsigned int rect_height, 
			       unsigned int rect_depth)
{
#ifdef DYNAMIC_CHECKING
  if (!(((rect_ul_x + rect_width) < _width)
	&&((rect_ul_y + rect_height) < _height)
	&&((rect_ul_z + rect_depth) < _depth)))
  {
    ERROR("Volume evaluate patch out of bounds\n");
  }
  else
#endif
    /* fix _height in volume.c */
    for (int h = 0; h < rect_height; h++)
      for (int i = 0; i < rect_height; i++)
	for (int j = 0; j < rect_width; j++)
	  VISVolumeRep<T>::_buffer[index((j + rect_ul_x), (i + rect_ul_y), 
					 (h + rect_ul_z))] = f(j, i, h);// **mdm**
	
}


	

template< class T >
VISVolume<T> VISVolume<T>::power(const int& power)
{
  VISVolume<T> im_return = createToSize();
  (im_return.repRef())->power(this->rep(), power);
  return(im_return);
}



template<class T>
VISVolume<T>& VISVolume<T>::operator*=(const VISVolume<T>& volume)
{
  multAssign(volume);
  return(*this);
}

template<class T>
VISVolume<T>& VISVolume<T>::operator*=(const T& value)
{
  multAssign(value);
  return(*this);
}


template<class T>
VISVolume<T>& VISVolume<T>::operator+=(const VISVolume<T>& volume)
{
  addAssign(volume);
  return(*this);
}

template<class T>
VISVolume<T>& VISVolume<T>::operator+=(const T& value)
{
  addAssign(value);
  return(*this);
}


template<class T>
VISVolume<T>& VISVolume<T>::operator-=(const VISVolume<T>& volume)
{
  subAssign(volume);
  return(*this);
}

template<class T>
VISVolume<T>& VISVolume<T>::operator-=(const T& value)
{
  subAssign(value);
  return(*this);
}


template<class T>
VISVolume<T>& VISVolume<T>::operator/=(const VISVolume<T>& volume)
{
  divAssign(volume);
  return(*this);
}

template<class T>
VISVolume<T>& VISVolume<T>::operator/=(const T& value)
{
  multAssign((T)1 / value);
  return(*this);
}


template< class T >
VISVolume<T>  VISVolume<T>::add(T value) const
{
  VISVolume<T> volume_return = createToSize();
  (volume_return.repRef())->add((rep()), 
				value);
  return(volume_return);
}

template< class T >
VISVolume<T>  VISVolume<T>::sub(T value) const
{
  VISVolume<T> volume_return = createToSize();
  (volume_return.repRef())->sub((rep()), 
				value);
  return(volume_return);
}

template< class T >
VISVolume<T>  VISVolume<T>::sub_from(T value) const
{
  VISVolume<T> volume_return = createToSize();
  (volume_return.repRef())->sub(value, rep());
  return(volume_return);
}

template< class T >
VISVolume<T> VISVolume<T>::mult(T value) const
{
  VISVolume<T> volume_return = createToSize();
  (volume_return.repRef())->mult((rep()), value);
  return(volume_return);
}

template< class T >
VISVolume<T>  VISVolume<T>::div(T value) const
{
  VISVolume<T> volume_return = createToSize();
  (volume_return.repRef())->div((rep()), 
				value);
  return(volume_return);
}

template< class T >
VISVolume<T>  VISVolume<T>::div_by(T value) const
{
  VISVolume<T> volume_return = createToSize();
  (volume_return.repRef())->div(value, (rep()));
  return(volume_return);
}


template< class T >
VISVolume<T>  VISVolume<T>::gt(T value) const
{
  VISVolume<T> volume_return = createToSize();
  (volume_return.repRef())->gt((rep()), 
			       value);
  return(volume_return);
}

template< class T >
VISVolume<T>  VISVolume<T>::gteq(T value) const
{
  VISVolume<T> volume_return = createToSize();
  (volume_return.repRef())->gteq((rep()), 
				 value);
  return(volume_return);
}

template< class T >
VISVolume<T>  VISVolume<T>::lt(T value) const
{
  VISVolume<T> volume_return = createToSize();
  (volume_return.repRef())->lt((rep()), value);
  return(volume_return);
}

template< class T >
VISVolume<T>  VISVolume<T>::lteq(T value) const
{
  VISVolume<T> volume_return = createToSize();
  (volume_return.repRef())->lteq(rep(), value);
  return(volume_return);
}

template< class T >
VISVolume<T>  VISVolume<T>::eq(T value) const
{
  VISVolume<T> volume_return = createToSize();
  (volume_return.repRef())->eq(rep(), value);
  return(volume_return);
}


template<class T>
VISVolume<T> VISVolume<T>::gauss(float sigma) const
{
  if (sigma <= 0.0)
  {
    return(*this);
  }
  /* this needs to be optimized for the float case */

  if (type() != VISIm::FLOAT)
  {
    // The idea here is to scale the kernels up to fit into the range
    // indicated by the type.

    float scale_x;
    float scale_y;
    float scale_z;
    VISImageType<T> type;
    VISVolume<T> local_col_kernel;
    VISVolume<float> col_kernel;
    col_kernel = gauss_col_kernel(sigma);
    scale_y = type.rangeMax()/(col_kernel.max()*
			       col_kernel.height());
    assignVolume(scale_y*col_kernel, local_col_kernel);

    VISVolume<T> local_row_kernel;
    VISVolume<float> row_kernel;
    row_kernel = gauss_row_kernel(sigma);
    scale_x = type.rangeMax()/
      (row_kernel.max()*col_kernel.width());
    assignVolume(scale_x*row_kernel, local_row_kernel);

    VISVolume<T> local_slice_kernel;
    VISVolume<float> slice_kernel;
    row_kernel = gauss_slice_kernel(sigma);
    scale_z = type.rangeMax()/
      (slice_kernel.max()*slice_kernel.depth());
    assignVolume(scale_z*slice_kernel, local_slice_kernel);
	    
    return(
	   (((convolve(local_col_kernel))
	     .convolve(local_row_kernel))
	    .convolve(local_slice_kernel))
	   .scale(1.0/(scale_x*scale_y*scale_z))
	   );
  }
  else
  {

    // these funny assignments are done so that the thing doesn't choke
    // on the template stuff.  This is compiler dependent
    VISVolume<T> local_col_kernel;
    VISVolume<float> col_kernel;
    col_kernel = gauss_col_kernel(sigma);
    assignVolume(col_kernel, local_col_kernel);

    VISVolume<T> local_row_kernel;
    VISVolume<float> row_kernel;
    row_kernel = gauss_row_kernel(sigma);
    assignVolume(row_kernel, local_row_kernel);

    VISVolume<T> local_slice_kernel;
    VISVolume<float> slice_kernel;
    slice_kernel = gauss_slice_kernel(sigma);
    assignVolume(slice_kernel, local_slice_kernel);

    return(
	   (((convolve(local_col_kernel))
	     .convolve(local_row_kernel))
	    .convolve(local_slice_kernel))
	   );
  }
}

template<class T>
VISVolume<T> VISVolume<T>::derivative(unsigned int order_x, 
				      unsigned int order_y, 
				      unsigned int order_z) const
{
  VISVolume<T> local_dx_kernel;
  VISVolume<T> tmp_dx_kernel = 2.0f*dx_kernel(order_x);// **mdm**
  assignVolume(tmp_dx_kernel, local_dx_kernel);// **mdm **
  //assignVolume(2.0f*dx_kernel(order_x), local_dx_kernel);// **mdm**

  VISVolume<T> local_dy_kernel;
  VISVolume<T> tmp_dy_kernel = 2.0f*dy_kernel(order_y);// **mdm**
  assignVolume(tmp_dy_kernel, local_dy_kernel);// **mdm**
  //assignVolume(2.0f*dy_kernel(order_y), local_dy_kernel);// **mdm**

  VISVolume<T> local_dz_kernel;
  assignVolume(2.0f*dz_kernel(order_z), local_dz_kernel);

  return((((convolve(local_dx_kernel)).convolve(local_dy_kernel))
	  .convolve(local_dz_kernel))
	 .divAssign((T)8));
}

template< class T >
VISVolume<T> VISVolume<T>::sqrt() const
{
  VISVolume<T> r = this->createToSize();
  r.repRef()->sqrt(rep());
  return(r);
}

template< class T >
VISVolume<T> VISVolume<T>::abs() const
{
  VISVolume<T> r = this->createToSize();
  r.repRef()->abs(rep());
  return(r);
}

template< class T >
VISVolume<T> VISVolume<T>::scale(float value) const
{
  VISVolume<T> r = this->createToSize();
  r.repRef()->scale(rep(), value);
  return(r);
}

template<class T>
VISVolume<T> VISVolume<T>::dx(unsigned int order) const
{
  VISVolume<T> local_dx_kernel;
  VISVolume<float> float_dx_kernel = dx_kernel(order);
  assignVolume((float)2.0*float_dx_kernel, local_dx_kernel);
  return(convolve(local_dx_kernel).divAssign((T)2)); 
}


template<class T>
VISVolume<T> VISVolume<T>::dy(unsigned int order) const
{
  VISVolume<T> local_dy_kernel;
  assignVolume(2.0f*VISVolume(dy_kernel(order)), local_dy_kernel); 
  return((convolve(local_dy_kernel)).divAssign((T)2)); 
}

template<class T>
VISVolume<T> VISVolume<T>::dz(unsigned int order) const
{
  VISVolume<T> local_dz_kernel;
  VISVolume<float> float_dz_kernel = dz_kernel(order);
  assignVolume(2.0f*float_dz_kernel, local_dz_kernel);
  return(convolve(local_dz_kernel).divAssign((T)2)); 
}


template< class T >
VISVolume<T> VISVolume<T>::dx() const 
{
  return(dx(1));
}

template< class T >
VISVolume<T> VISVolume<T>::dy() const 
{
  return(dy(1));
}

template< class T >
VISVolume<T> VISVolume<T>::dz() const 
{
  return(dz(1));
}

// these half derivatives are useful in diffusion processes
// need to be zeroed out along the border, by Zhong on 5/31/100
template< class T >
VISVolume<T> VISVolume<T>::dxHalfForward() const
{
  VISVolume<T> r = this->createToSize();
  int h = r.height();
  int w = r.width();
  int d = r.depth();

  for (int k = 0; k < d; k++)
  {
    for (int j = 0; j < h; j++) {
      for (int i = 0; i < w - 1; i++)
	r.at(i, j, k) = this->itemAt(i + 1, j, k) 
	  - this->itemAt(i, j, k);
      r.at(w - 1, j, k) = (T) 0;
    }
  }
  return(r);
}

template< class T >
VISVolume<T> VISVolume<T>::dyHalfForward() const
{
  VISVolume<T> r = this->createToSize();
  int h = r.height();
  int w = r.width();
  int d = r.depth();

  for (int k = 0; k < d; k++)
  {
    for (int i = 0; i < w; i++) {
      for (int j = 0; j < h - 1; j++)
	r.at(i, j, k) = this->itemAt(i, j + 1, k) 
	  - this->itemAt(i, j, k);
      r.at(i, h - 1, k) = (T) 0;
    }
  }
  return(r);
}

template< class T >
VISVolume<T> VISVolume<T>::dzHalfForward() const
{
  VISVolume<T> r = this->createToSize();
  int h = r.height();
  int w = r.width();
  int d = r.depth();

  for (int j = 0; j < h; j++)
  {
    for (int i = 0; i < w; i++) {
      for (int k = 0; k < d - 1; k++)
	r.at(i, j, k) = this->itemAt(i, j, k + 1) 
	  - this->itemAt(i, j, k);
      r.at(i, j, d - 1) = (T) 0;
    }

  }
  return(r);
}

template< class T >
VISVolume<T> VISVolume<T>::dxHalfBack() const
{
  VISVolume<T> r = this->createToSize();
  int h = r.height();
  int w = r.width();
  int d = r.depth();

  for (int k = 0; k < d; k++)
  {
    for (int j = 0; j < h; j++) {
      r.at(0, j, k) = (T) 0;
      for (int i = 1; i < w; i++)
	r.at(i, j, k) = this->itemAt(i, j, k) 
	  - this->itemAt(i - 1, j, k);
    }
  }
  return(r);
}

template< class T >
VISVolume<T> VISVolume<T>::dyHalfBack() const
{
  VISVolume<T> r = this->createToSize();
  int h = r.height();
  int w = r.width();
  int d = r.depth();

  for (int k = 0; k < d; k++)
  {
    for (int i = 0; i < w; i++) {
      r.at(i, 0, k) = (T) 0;	
      for (int j = 1; j < h; j++)
	r.at(i, j, k) = this->itemAt(i, j, k) 
	  - this->itemAt(i, j - 1, k);
    }
  }
  return(r);
}

template< class T >
VISVolume<T> VISVolume<T>::dzHalfBack() const
{
  VISVolume<T> r = this->createToSize();
  int h = r.height();
  int w = r.width();
  int d = r.depth();

  for (int j = 0; j < h; j++)
  {
    for (int i = 0; i < w; i++) {
      r.at(i, j, 0) = (T) 0;
      for (int k = 1; k < d; k++)
	r.at(i, j, k) = this->itemAt(i, j, k) 
	  - this->itemAt(i, j, k - 1);
    }

  }
  return(r);
}


template< class T >
T VISVolume<T>::interp(float x, float y, float z) const
{
  unsigned int x_lo, x_hi, y_lo, y_hi, z_lo, z_hi;
  float x_in, y_in, z_in;
  x_lo = (unsigned int)floor(x);
  y_lo = (unsigned int)floor(y);
  z_lo = (unsigned int)floor(z);
  x_hi = (unsigned int)ceil(x);
  y_hi = (unsigned int)ceil(y); 
  z_hi = (unsigned int)ceil(z); 

  //    cout << "z " << z_lo << " " << z_hi << endl;
  //    cout << "y " << y_lo << " " << y_hi << endl;
  //    cout << "x " << z_lo << " " << x_hi << endl;

  x_in = x - (float)x_lo;
  y_in = y - (float)y_lo;
  z_in = z - (float)z_lo;

  if ((x < 0.0)||(y < 0.0)||(x > (width() - 1))||(y > (height() - 1))
      ||(z < 0.0)||(z > (depth() - 1)))
  {
    printf("interp out of bounds %f %f %f\n", x, y, z);
    return((T)0);
  }

  T z_lo_val, z_hi_val;

  //    cout << rep()->itemAt(x_hi, y_lo, z_lo) << endl;
  //    cout << rep()->itemAt(x_lo, y_i, z_lo) << endl;

  z_lo_val = 
    (1.0f - y_in)*(x_in*rep()->itemAt(x_hi, y_lo, z_lo) 
		   + (1.0f - x_in)*rep()->itemAt(x_lo, y_lo, z_lo))
    + y_in*(x_in*rep()->itemAt(x_hi, y_hi, z_lo)
	    + (1.0f - x_in)*rep()->itemAt(x_lo, y_hi, z_lo));
    
  z_hi_val = 
    ((1.0f - y_in)*(x_in*rep()->itemAt(x_hi, y_lo, z_hi) 
		    + (1.0f - x_in)*rep()->itemAt(x_lo, y_lo, z_hi))
     + y_in*(x_in*rep()->itemAt(x_hi, y_hi, z_hi) 
	     + (1.0f - x_in)*rep()->itemAt(x_lo, y_hi, z_hi)));
    
  return((1.0f - z_in)*z_lo_val + z_in*z_hi_val);
}

template< class T >
void VISVolume<T>::putRep(const VISVolumeRep<T>* r)
{
  unref(_rep);
  ref(_rep = (VISRep*)r);
} 



template < class T >
void VISVolume<T>::print() const
{
  printf("VISVolume<T> height %d width %d depth %d and type %d\n",
	 height(), width(), depth(), type());
}

template < class T >
VISVolume<T>& VISVolume<T>::operator=(T value)
{
  repRef()->clear(value);
  return(*this);
}


template<class T>
VISVolume<T> VISVolume<T>::convolve(const VISVolume<T>& kernel) const
{

  unsigned int w = kernel.width();
  unsigned int h = kernel.height();
  unsigned int d = kernel.depth();
  unsigned int x_center = w/2;
  unsigned int y_center =  h/2; 
  unsigned int z_center =  d/2; 

  //    printf("w %d h %d  d %d\n", w, h, d);

  const T *buffer_in, *buf_in;
  T *buffer_out, *buf_out;
  const T *buf_k;
  T k_value;
  int h_offset, w_offset, d_offset;
  // , w_in, h_in;
  //    w_in = width();
  //    h_in = height();

  int j, k, l, m, n, p;

  int this_d = depth(), this_h = height(), this_w = width();

  VISVolume<T> im_return = this->createToSize();
  im_return = (T)0.0;


  buf_k = kernel.rep()->buffer() + (w*h*d - 1);
  buffer_out = im_return.repRef()->bufferRef();
  buffer_in = rep()->buffer();

  // the outer loop takes place over the kernel
  for (n = 0; n < d; n++)
    for (j = 0; j < h; j++)
      for (k = 0; k < w; k++)
      {
	//		    cout << "kernel " << n << " " << 
	//			j << " " << k << endl;
	h_offset = y_center - j;
	w_offset = x_center - k;
	d_offset = z_center - n;
	k_value = *(buf_k--);
	buf_out = buffer_out;
	buf_in = buffer_in + d_offset*this_w*this_h 
	  + h_offset*this_w + w_offset;

	//		    cout << "offset " << d_offset*this_w*this_h 
	//			+ h_offset*this_w + w_offset << endl;
	//		    cout << "k_value " << k_value << endl;
		
	if (d_offset < 0)
	{
	  buf_out += -(d_offset*this_w*this_h);
	  buf_in += -(d_offset*this_w*this_h);
	}
	for (p = ::VISmax(0, d_offset);
	     p < ::VISmin(this_d, 
			  this_d + d_offset); p++)

	{
	  if (h_offset < 0)
	  {
	    buf_out += -(h_offset*this_w);
	    buf_in += -(h_offset*this_w);
	  }
		
	  for (l = ::VISmax(0, h_offset);
	       l < ::VISmin(this_h, 
			    this_h + h_offset); l++)
	  {
	    if (w_offset < 0)
	    {
	      buf_out += -(w_offset);
	      buf_in += -(w_offset);
	    }
	    for (m = ::VISmax(0, w_offset); 
		 m < ::VISmin(this_w, 
			      this_w + w_offset); m++)
	      *buf_out++ += k_value*(*buf_in++);
	    if (w_offset > 0)
	    {
	      buf_out += (w_offset);
	      buf_in += (w_offset);
	    }
	  }
	}
	if (h_offset > 0)
	{
	  buf_out += (h_offset*this_w);
	  buf_in += (h_offset*this_w);
	}
      }
  return(im_return);
}

template< class T >
float VISVolume<T>::maskFloat(const VISVolume<float>& mask, 
			      unsigned int center_x, 
			      unsigned int center_y, 
			      unsigned int center_z, 
			      unsigned int x_pos, 
			      unsigned int y_pos,
			      unsigned int z_pos
			      ) const
{
  int i, j, k;
  int w = width(), h = height(), d = depth();
  int m_w = mask.width(), m_h = mask.height(), m_d = mask.depth();
  float total = 0.0;
  const float *buf_m = mask.rep()->buffer() + m_w*m_h*m_d - 1;
  int y_lo  = (y_pos + (center_y - m_h));
  int y_hi  = (y_pos - (center_y - m_h));

  int x_lo  = (x_pos + (center_x - m_w));
  int x_hi  = (x_pos - (center_x - m_w));

  int z_lo  = (z_pos + (center_z - m_d));
  int z_hi  = (z_pos - (center_z - m_d));

  int x_pos_tmp, y_pos_tmp, z_pos_tmp;
    
    
  if ((x_lo >= 0)&&(x_hi < w)&&(y_lo >= 0)
      &&(y_hi < h)
      &&(z_lo >= 0)&&(z_hi < d))
  {
    const T *buf_this = rep()->buffer() 
      + (z_pos - center_z)*w*h 
      // replaced "*height()" with "*width()", which is correct
      + (y_pos - center_y)*w
      + (x_pos - center_x);
	    
    for (k = 0; k < m_d; k++)
    {
      for (i = 0; i < m_h; i++)
      {
	for (j = 0; j < m_w; j++)	
	{
	  total += *(buf_this++)*(*(buf_m--));
	}
	buf_this += (w - m_w);
      }
      buf_this += (h - m_h)*w;
    }
  }
  else
  {
    for (k = m_d - 1; k >= 0; k--)
      for (i = m_h - 1; i >= 0; i--)
	for (j = m_w - 1; j >= 0; j--)
	{
	  x_pos_tmp = x_pos - center_x + j;
	  y_pos_tmp = y_pos - center_y + i;
	  z_pos_tmp = z_pos - center_z + k;
	  if ((x_pos_tmp >= 0)&&(x_pos_tmp < w)
	      &&(y_pos_tmp >= 0)&&(y_pos_tmp < h)
	      &&(z_pos_tmp >= 0)&&(z_pos_tmp < d))
	  {
	    total += itemAt(x_pos_tmp, 
			    y_pos_tmp, 
			    z_pos_tmp)
	      *mask.itemAt(j, i, k);
	  }
	}
  }
  return(total);
}

template< class T >
VISVolume<T> VISVolume<T>::resample(float the_scale) const
{
  return(resample(the_scale, the_scale, the_scale, width()/2.0 - 0.5, 
		  height()/2.0 - 0.5, depth()/2.0 -0.5));
}

template< class T >
VISVolume<T> VISVolume<T>::resample(float scale_x, float scale_y, 
				    float scale_z, float x, float y, 
				    float z) const
{
  unsigned new_w, new_h, new_d;
  float new_x, new_y, new_z;    
  float tmp_x, tmp_y, tmp_z;    
  unsigned w = width(), h = height(), d = depth();
  int i, j, k;
    
  new_w  = (float)w*(float)scale_x;
  new_h  = (float)h*(float)scale_y;
  new_d  = (float)d*(float)scale_z;

  //    new_x = x*scale_x + (float_w - (float)new_w)/2.0;
  //    new_y = y*scale_y + (float_h - (float)new_h)/2.0;
  //    new_z = z*scale_z + (float_d - (float)new_d)/2.0;

  new_x = (x + 0.5)*scale_x - 0.5;
  new_y = (y + 0.5)*scale_y - 0.5;
  new_z = (z + 0.5)*scale_z - 0.5;

  VISVolume<T> volume_return(new_w, new_h, new_d);
  for (k = 0; k < new_d; k++)
    for (j = 0; j < new_h; j++)
      for (i = 0; i < new_w; i++)
      {
	tmp_x = ((i - new_x)/scale_x + x); 
	tmp_y = ((j - new_y)/scale_y + y); 
	tmp_z = ((k - new_z)/scale_z + z); 		    
		    
	tmp_x = MAX(tmp_x, 0);
	tmp_x = MIN(tmp_x, w - 1);

	tmp_y = MAX(tmp_y, 0);
	tmp_y = MIN(tmp_y, h - 1);

	tmp_z = MAX(tmp_z, 0);
	tmp_z = MIN(tmp_z, d - 1);

	volume_return.at(i, j, k) 
	  = interp(tmp_x, tmp_y, tmp_z);
      }
  return(volume_return);
}

template< class T >
VISVolume<T> VISVolume<T>::resample(unsigned w, unsigned h, unsigned d) const
{
  float scale_w, scale_h, scale_d;

  scale_w = (float)w/(float)width();
  scale_h = (float)h/(float)height();
  scale_d = (float)d/(float)depth();

  return(resample(scale_w, scale_h, scale_d, 
		  (float)width()/2.0 - 0.5, 
		  (float)height()/2.0 - 0.5, 
		  (float)depth()/2.0 - 0.5
		  ));
}

template< class T >
VISVolume<float> VISVolume<T>::noise() const
{
  VISVolume<float> r(width(), height(), depth());
  float* buf;
  int size = r.height()*r.width()*r.depth();
  buf = r.repRef()->bufferRef();
  for (int i = 0; i < size; i++)
  {
    *buf++ = rand1();
  }
  return(r);
}





template< class T >
VISVolume<T>  VISVolume<T>::reduce(int scale) const
{
  int new_w, new_h, new_d, w, h, d;
  VISVolume<T> r = VISVolume<T>(new_w = (w = width())/scale, 
				new_h = (h = height())/scale, 
				new_d = (d = depth())/scale);
  int i, j, k;
  int x_mod, y_mod, z_mod;
  int x_offset, y_offset, z_offset;
  const T* buf_in;
  T* buf_out;


  if (scale <= 0)
    return(*this);

  x_mod = width()%scale;
  y_mod = height()%scale;
  z_mod = depth()%scale;

  x_offset = (scale - 1)/2 + (int)ceil((double)x_mod/2.0);
  y_offset = (scale - 1)/2 + (int)ceil((double)y_mod/2.0);
  z_offset = (scale - 1)/2 + (int)ceil((double)z_mod/2.0);

#ifdef DEBUG
  printf("Volume::reduce x_mod %d x_offset %d, y_mod %d y_offset %d, z_mod %d z_offset %d\n", 
	 x_mod, x_offset, y_mod, y_offset, z_mod, z_offset);
#endif
    
  int position = 0;
    
  buf_in = rep()->buffer() + height()*width()*(z_offset) 
    + width()*(y_offset) + x_offset;

  //    position = height()*width()*(z_offset) + width()*(y_offset) + x_offset;
  buf_out = r.repRef()->bufferRef();
  for (i = 0; i < new_d; i++)
  {
    for (j = 0; j < new_h; j++)
    {
      for (k = 0; k < new_w; k++)
      {
	//		printf(" the 3d position %d\n", position);
	(*buf_out++) = *buf_in;
	buf_in += scale;
	//		position += scale;
      }
      buf_in += x_mod + (scale - 1)*w;
      //	    position += x_mod + (scale - 1)*width();
    }
    buf_in += y_mod*w + (scale - 1)*w*h;
    //	position += y_mod + (scale - 1)*width()*height();
  }
  return(r);
}


template< class T >
float VISVolume<T>::average()
{
  return(rep()->sum()/(float)(height()*width()*depth()));
}



/* these are non-member functions */


template< class T1, class T2 >
void copy(const VISVolume<T1>& from, VISVolume<T2>& to)
{
  if (compareSize(to, from))
    copy((from.rep()), (to.repRef()));
  else
    ERROR("VISVolume<>: copy - volume sizes not compatible");
}




template< class T1, class T2 >
void assignVolume(const VISVolume<T1>& from, VISVolume<T2>& to)
{
  if (from.type() != to.type())
  {
    VISVolumeRep<T2>* rep_tmp;
    rep_tmp = new VISVolumeRep<T2>(from.rep()->width(), 
				   from.rep()->height(), 
				   from.rep()->depth());
    copy(from.rep(), rep_tmp);
    to.putRep(rep_tmp);
  }
  else
    to.assign(from);
}

template< class T1, class T2 >
int compareSize(const VISVolume<T1>& a, const VISVolume<T2>& b)
{
  compareSize(a.rep(), b.rep());
}

template< class T1, class T2 >
int compareSize(const VISVolumeRep<T1>* a, const VISVolumeRep<T2>* b)
{
  return((a->width() == b->width())
	 &&(a->height() == b->height())
	 &&(a->depth() == b->depth()));
}
	

template< class T1, class T2 >
void copy(const VISVolumeRep<T1>* a, VISVolumeRep<T2>* b)
{
  if (!(compareSize(a, b)))
    printf("ERROR copy\n");
  else
    copy((VISImRep<T1>*)a, (VISImRep<T2>*) b);
}

template< class T >
const VISImage<T> VISVolume<T>::image() const
{
  const T *buf;
  T **buf_array;
  VISImage<T> r;
  buf = rep()->buffer();
  buf_array = new T*[depth()];

  for (int i = 0; i < depth(); i++)
    buf_array[i] = (T*)buf + i*height()*width();
  r = VISImage<T>(width(), height(), depth(), buf_array);
    
  delete buf_array;
  return(r);
}

template<class T>
VISImage<T> VISVolume<T>::imageRef()
{
  T *buf;
  T **buf_array;
  VISImage<T> r;
  buf = rep()->buffer();
  buf_array = new T*[depth()];

  for (int i = 0; i < depth(); i++)
    buf_array[i] = buf + i*height()*width();
  r = VISImage<T>(width(), height(), depth(), buf_array);
    
  delete buf_array;
  return(r);
}

template<class T>
VISImage<T>  VISVolume<T>::image(unsigned int slice) 
{
  const T *buf;
  T *image_buf;
  int w, h;
  int i;
  VISImage<T> r(w = width(), h = height());
  buf = rep()->buffer();
  image_buf = r.repRef()->bufferRef();

  for (i = 0; i < w*h; i++)
    *image_buf++ = *(buf + slice*h*w + i);
  //    r = VISImage<T>(width(), height(), 1, buf_array);
  return(r);
}

template<class T>
VISImage<T> VISVolume<T>::imageRef(unsigned int slice)
{
  T *buf;
  T *buf_array[1];
  VISImage<T> r;
  buf = repRef()->bufferRef();

  buf_array[0] = buf + slice*height()*width();
  r = VISImage<T>(width(), height(), 1, buf_array);

  return(r);
}

template< class T >
VISVolume<T>::VISVolume(const VISImage<T>& image)
{
  initialize(image.width(), image.height(), image.channels());

  const T *buf_in;
  T *buf_out;
  buf_out = repRef()->bufferRef();
    
  for (int i = 0; i < depth(); i++)
  {
    buf_in = image.rep(i)->buffer();
    for (int j = 0; j < height(); j++)
      for (int k = 0; k < width(); k++)
	*buf_out++ = *buf_in++;
  }
}


template< class T >
VISVolume<T> operator*(const VISVolume<T>& from, T value) 
{
  return(from.mult(value));
}

template< class T >
VISVolume<T> operator*(T value, const VISVolume<T>& from) 
{
  return(from.mult(value));
}



template< class T >
VISVolume<float> VISVolume<T>::cityBlockDistTrans() const
{
  float new_x_dist, new_y_dist, new_z_dist, new_dist;
  int w = width(), h = height(), d = depth();
  VISVolume<short> x_pos, y_pos, z_pos, visited;
  x_pos = VISVolume<short>(w, h, d);
  y_pos = VISVolume<short>(w, h, d);
  z_pos = VISVolume<short>(w, h, d);
  visited = VISVolume<short>(w, h, d);
  VISVolume<float> distance = VISVolume<float>(w, h, d);

  int i, j, k;

  for (k = 0; k < d; k++)
    for (i = 0; i < h; i++)
      for (j = 0; j < w; j++)
      {
	if ((int)itemAt(j, i, k))
	{
	  distance.at(j, i, k) = (float)0.0;
	  x_pos.at(j, i, k) = j;
	  y_pos.at(j, i, k) = i;
	  z_pos.at(j, i, k) = k;
	  visited.at(j, i, k)  = 1;
	}
	else
	{
	  distance.at(j, i, k) 
	    = (float)((w + h + d + 1)*(w + h + d + 1));
	  x_pos.at(j, i, k) = 0;
	  y_pos.at(j, i, k) = 0;
	  z_pos.at(j, i, k) = 0;
	  visited.at(j, i, k)  = 0;
	}
      }

  for (k = 0; k < d; k++ )
    for (i = 0; i < (h); i++)
      for (j = 1; j < (w); j++)
      {
	if (visited.itemAt(j - 1, i, k))
	{
	  new_x_dist 
	    = (float)x_pos.at(j - 1, i, k) - (float)(j);
	  new_y_dist 
	    = (float)y_pos.at(j - 1, i, k) - (float)(i);
	  new_z_dist 
	    = (float)z_pos.at(j - 1, i, k) - (float)(k);
	  new_dist = new_x_dist*new_x_dist + 
	    new_y_dist*new_y_dist + new_z_dist*new_z_dist;
	  if (new_dist < distance.at(j, i, k))
	  {
	    distance.at(j , i, k) = new_dist;
	    x_pos.at(j, i, k) 
	      = x_pos.itemAt(j - 1, i, k);
	    y_pos.at(j, i, k) 
	      = y_pos.itemAt(j - 1, i, k);
	    z_pos.at(j, i, k) 
	      = z_pos.itemAt(j - 1, i, k);
	    visited.at(j, i, k) = 1;
	  }
	}
      }

  for (k = 0; k < d; k++ )
    for (i = 0; i < (h); i++)
      for (j = (w - 2); j >= 0; j--)
      {
	if (visited.itemAt(j + 1, i, k))
	{
	  new_x_dist 
	    = (float)x_pos.at(j + 1, i, k) - (float)(j);
	  new_y_dist 
	    = (float)y_pos.at(j + 1, i, k) - (float)(i);
	  new_z_dist 
	    = (float)z_pos.at(j + 1, i, k) - (float)(k);
	  new_dist = new_x_dist*new_x_dist + 
	    new_y_dist*new_y_dist +
	    new_z_dist*new_z_dist;
			
	  if (new_dist < distance.at(j, i, k))
	  {
	    distance.at(j , i, k) = new_dist;
	    x_pos.at(j, i, k) 
	      = x_pos.itemAt(j + 1, i, k);
	    y_pos.at(j, i, k) 
	      = y_pos.itemAt(j + 1, i, k);
	    z_pos.at(j, i, k) 
	      = z_pos.itemAt(j + 1, i, k);
	    visited.at(j, i, k) = 1;
	  }
	}
      }


  for (k = 0; k < d; k++ )
    for (j = 0; j < w; j++)
      for (i = 1; i < (h); i++)
      {
	if (visited.itemAt(j, i - 1, k))
	{
	  new_x_dist 
	    = (float)x_pos.at(j, i - 1, k) - (float)(j);
	  new_y_dist 
	    = (float)y_pos.at(j, i - 1, k) - (float)(i);
	  new_z_dist 
	    = (float)z_pos.at(j, i - 1, k) - (float)(k);

	  new_dist = new_x_dist*new_x_dist + 
	    new_y_dist*new_y_dist +
	    new_z_dist*new_z_dist;
			
	  if (new_dist < distance.at(j, i, k))
	  {
	    distance.at(j , i, k) = new_dist;
	    x_pos.at(j, i, k) 
	      = x_pos.itemAt(j, i - 1, k);
	    y_pos.at(j, i, k) 
	      = y_pos.itemAt(j, i - 1, k);
	    z_pos.at(j, i, k) 
	      = z_pos.itemAt(j, i - 1, k);
	    visited.at(j, i, k) = 1;
	  }
	}
      }

  for (k = 0; k < d; k++ )
    for (j = 0; j < w; j++)
      for (i = h - 2; i >= 0; i--)
      {
	if (visited.itemAt(j, i + 1, k))
	{
	  new_x_dist 
	    = (float)x_pos.at(j, i + 1, k) - (float)(j);
	  new_y_dist 
	    = (float)y_pos.at(j, i + 1, k) - (float)(i);
	  new_z_dist 
	    = (float)z_pos.at(j, i + 1, k) - (float)(k);
	  new_dist = new_x_dist*new_x_dist + 
	    new_y_dist*new_y_dist + 
	    new_z_dist*new_z_dist;
			
	  if (new_dist < distance.at(j, i, k))
	  {
	    distance.at(j, i, k) = new_dist;
	    x_pos.at(j, i, k) 
	      = x_pos.itemAt(j, i + 1, k);
	    y_pos.at(j, i, k) 
	      = y_pos.itemAt(j, i + 1, k);
	    z_pos.at(j, i, k) 
	      = z_pos.itemAt(j, i + 1, k);
	    visited.at(j, i, k) = 1;
	  }
	}
      }


  for (i = 0; i < h; i++ )
    for (j = 0; j < w; j++)
      for (k = 1; k < (d); k++)
      {
	if (visited.itemAt(j, i, k - 1))
	{
	  new_x_dist 
	    = (float)x_pos.at(j, i, k - 1) - (float)(j);
	  new_y_dist 
	    = (float)y_pos.at(j, i, k - 1) - (float)(i);
	  new_z_dist 
	    = (float)z_pos.at(j, i, k - 1) - (float)(k);

	  new_dist = new_x_dist*new_x_dist + 
	    new_y_dist*new_y_dist +
	    new_z_dist*new_z_dist;
			
	  if (new_dist < distance.at(j, i, k))
	  {
	    distance.at(j , i, k) = new_dist;
	    x_pos.at(j, i, k) 
	      = x_pos.itemAt(j, i, k - 1);
	    y_pos.at(j, i, k) 
	      = y_pos.itemAt(j, i, k - 1);
	    z_pos.at(j, i, k) 
	      = z_pos.itemAt(j, i, k - 1);
	    visited.at(j, i, k) = 1;
	  }
	}
      }

  for (i = 0; i < h; i++ )
    for (j = 0; j < w; j++)
      for (k = d - 2; k >= 0; k--)
      {
	if (visited.itemAt(j, i, k + 1))
	{
	  new_x_dist 
	    = (float)x_pos.at(j, i, k + 1) - (float)(j);
	  new_y_dist 
	    = (float)y_pos.at(j, i, k + 1) - (float)(i);
	  new_z_dist 
	    = (float)z_pos.at(j, i, k + 1) - (float)(k);
	  new_dist = new_x_dist*new_x_dist + 
	    new_y_dist*new_y_dist + 
	    new_z_dist*new_z_dist;
			
	  if (new_dist < distance.at(j, i, k))
	  {
	    distance.at(j, i, k) = new_dist;
	    x_pos.at(j, i, k) 
	      = x_pos.itemAt(j, i, k + 1);
	    y_pos.at(j, i, k) 
	      = y_pos.itemAt(j, i, k + 1);
	    z_pos.at(j, i, k) 
	      = z_pos.itemAt(j, i, k + 1);
	    visited.at(j, i, k) = 1;
	  }
	}
      }

  return(distance.sqrt());
}


template< class T >
VISVolume<T> VISVolume<T>::zeroCrossings() const
{
  VISVolume<T> r = VISVolume<T>(this->createToSize());
  r = (T)0;
  const int num_buffers = 9;
  const T* buffers[num_buffers];
  T* buffers_out[num_buffers];
  const T* center;
  T* center_out;
  T this_one, that, abs_this_one;
  const T* buff_in = rep()->buffer();
  T* buff_out = r.repRef()->bufferRef();
  unsigned int w = width(), h = height(), d = depth();
  int i, j, k, l;
 

  // first solve for the slices and then do the z direction separately

  /* do first row here specially */

  for (l = 0; l < d; l++)
  {
    center = buff_in + 1 + l*w*h;
    center_out = buff_out + 1 + l*w*h;
    buffers[0] = buff_in + l*w*h;
    //		buffers[1] = buff_in + w + l*w*h;
    buffers_out[0] = buff_out + l*w*h;
    //		buffers_out[1] = buff_out + w + l*w*h;
    for (j = 1; j < w; j++)
    {
      this_one = *center;
      abs_this_one = tabs(this_one);
      for (k = 0; k < 1; k++)
      {	
	that = *(buffers[k]++);
	/* use exclusive or here */	
	if (((this_one < (T)0)&&(that > (T)0))||
	    ((this_one > (T)0)&&(that < (T)0))||
	    ((this_one == (T)0)&&(that != (T)0))||
	    ((this_one != (T)0)&&(that == (T)0)))
	{ 
	  that = tabs(that);
	  if (abs_this_one <= that)
	    *center_out = (T)1;
	  else
	    *buffers_out[k] = (T)1;

	} 
	buffers_out[k]++;
      } // k loop
      center++;
      center_out++;
    } // d  loop

    center = buff_in + w + 1 + l*w*h;
    buffers[0] = buff_in + 1 + l*w*h;
    buffers[1] = buff_in + w + l*w*h;
    //		buffers[2] = buff_in + l*w*h;
    //		buffers[3] = buff_in + 2*w + l*w*h;

    center_out = buff_out + w + 1 + l*w*h;
    buffers_out[0] = buff_out + 1 + l*w*h;
    buffers_out[1] = buff_out + w + l*w*h;
    //		buffers_out[2] = buff_out + l*w*h;
    //		buffers_out[3] = buff_out + 2*w + l*w*h;

    for (i = 1; i < h - 1; i++)
    {
      that = *(buffers[0] - 1);
      this_one = *buffers[1];
      if (((this_one < (T)0)&&(that > (T)0))||
	  ((this_one > (T)0)&&(that < (T)0))||
	  ((this_one == (T)0)&&(that != (T)0))||
	  ((this_one != (T)0)&&(that == (T)0))
	  )
      {
	abs_this_one = tabs(this_one);
	that = tabs(that);

	if (abs_this_one <= that)
	  *(buffers_out[0] - 1) = (T)1;
	else
	  *buffers_out[1] = (T)1;
      }
      for (j = 1; j < w; j++)
      {
	this_one = *center;
	abs_this_one 
	  = tabs(this_one);

	for (k = 0; k < 2; k++)
	{
	  that = *(buffers[k]++);
	  if (((this_one < (T)0)
	       &&(that > (T)0))||
	      ((this_one > (T)0)
	       &&(that < (T)0))||
	      ((this_one == (T)0)
	       &&(that != (T)0))||
	      ((this_one != (T)0)
	       &&(that == (T)0)))
	  {
	    that = tabs(that);
	    if (abs_this_one < that)
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
  }

  // now do zero crossings in the z-direction


  for (l = 0; l < d - 1; l++)
  {
    center = buff_in + (l + 1)*w*h;
    buffers[0] = buff_in + l*w*h;

    center_out = buff_out + (l + 1)*w*h;
    buffers_out[0] = buff_out + l*w*h;

    for (i = 0; i < h; i++)
    {
      for (j = 0; j < w; j++)
      {
	this_one = *center;
	abs_this_one = tabs(this_one);

	for (k = 0; k < 1; k++)
	{
	  that = *(buffers[k]++);
	  if (((this_one < (T)0)
	       &&(that > (T)0))||
	      ((this_one > (T)0)
	       &&(that < (T)0))||
	      ((this_one == (T)0)
	       &&(that != (T)0))||
	      ((this_one != (T)0)
	       &&(that == (T)0)))
	  {
	    that 
	      = tabs(that);
	    if (abs_this_one < that)
	      *center_out = (T)1;
	    else
	      *buffers_out[k] = (T)1;
	  }
	  (buffers_out[k]++);
	}
	center++;
	center_out++;
      }
    }
  }
  return(r);
}



template< class T >
void VISVolume<T>::printData() const
{
  for (int k = 0; k < depth(); k++)
  {
    for (int i = 0; i < height(); i++)
    {
      for (int j = 0; j < width(); j++)
	printf("%3.2f ", (float)itemAt(j, i, k));
      printf("\n");
    }
    printf("\n");
  }
}


template< class T >
VISArray<float>* VISVolume<T>::maskFloat(const VISArray<VISVolume<float> >& masks, 
					 unsigned int x_pos, 
					 unsigned int y_pos,
					 unsigned int z_pos
					 ) const
{
  VISArray<float>* r = new VISArray<float>();
  for (int i = 0; i < masks.n(); i++)
    r->appendItem(maskFloat(masks.itemAt(i), (masks.itemAt(i)).width()/2, 
			    (masks.itemAt(i)).height()/2, 
			    (masks.itemAt(i)).depth()/2, 
			    x_pos, y_pos, z_pos));
  return(r);
}

template< class T >
VISVolume<T> VISVolume<T>::floodFill(T thresh, unsigned x, unsigned y, 
				     unsigned z)
{
  VISVolIndexVISList list;
  VISVolume<T> r = createToSize();

  int N[6][3];
  N[0][0] = 1;
  N[0][1] = 0;
  N[0][2] = 0;
  N[1][0] = -1;
  N[1][1] = 0;
  N[1][2] = 0;
  N[2][0] = 0;
  N[2][1] = 1;
  N[2][2] = 0;
  N[3][0] = 0;
  N[3][1] = -1;
  N[3][2] = 0;
  N[4][0] = 0;
  N[4][1] = 0;
  N[4][2] = 1;
  N[5][0] = 0;
  N[5][1] = 0;
  N[5][2] = -1;

  r = (T)0;
  unsigned at_x, at_y, at_z, next_x, next_y, next_z;
  int i;

  if (itemAt(x, y, z) > thresh)
    return(r);

    
  list.appendItem(VISVolIndex(x, y, z));
  r.at(x, y, z) = (T)1.0;
  list.reset();
    
  int n = 0;
  while (list.valid())
  {
    at_x = list.atCurrent().a();
    at_y = list.atCurrent().b();
    at_z = list.atCurrent().c();
    //	    printf("flood fill iteration %d x %d y %d z %d\n", n++, at_x, 
    //		   at_y , at_z);
	    
	    
    for (i = 0; i < 6; i++)
    {
      next_x = at_x + N[i][0];
      next_y = at_y + N[i][1];
      next_z = at_z + N[i][2];
      if (checkBounds(next_x, next_y, next_z))
      {
	if ((r.itemAt(next_x, next_y, next_z) == (T)0.0)
	    &&(itemAt(next_x, next_y, next_z) <= thresh))
	{
	  r.at(next_x, next_y, next_z) = (T)1.0;
	  list.appendItem(VISVolIndex(next_x, 
				      next_y, 
				      next_z));
	  //				    printf("append item %d x %d y %d z %d\n", 
	  //					  next_x, next_y , next_z);
	}
      }
    }
    list.removeCurrent();
  }
  return(r);
}

template< class T >
VISVolume<T> VISVolume<T>::getROI(unsigned int x_pos, 
				  unsigned int y_pos,
				  unsigned int z_pos,
				  unsigned int w_roi, 
				  unsigned int h_roi, 
				  unsigned int d_roi) const
{
  VISVolume<T> r(w_roi, h_roi, d_roi);
  (r.repRef())->getROI(rep(), x_pos, y_pos, z_pos);
  return(r);
}


template< class T >
void VISVolume<T>::putROI(const VISVolume<T>& volume_in, 
			  unsigned int x_pos, 
			  unsigned int y_pos,
			  unsigned int z_pos)
{
  repRef()->insetRep(volume_in.rep(), x_pos, y_pos, z_pos);
}




#endif
