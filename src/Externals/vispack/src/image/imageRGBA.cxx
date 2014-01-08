// ******************************************************************
// VISPack. Copyright (c) 1994-2000 Ross Whitaker rtw@utk.edu       *
// For conditions of distribution and use, see the file LICENSE.txt *
// accompanying this distribution.                                  *
// ******************************************************************

// $Id: imageRGBA.cxx,v 1.1.1.1 2003/02/12 16:51:52 whitaker Exp $

#include "image/imageRGBA.h"

VISImageRGBA& VISImageRGBA::operator=(byte value)
{
    rgba* buf = (repRef())->bufferRef();
    for (int i = 0; i < rep(0)->size(); i++)
	{
	    *(buf++) = value;
	}
    return(*this);
}

VISImageRGBA& VISImageRGBA::operator=(const rgba& value)
{
    rgba* buf = (repRef())->bufferRef();
    for (int i = 0; i < rep(0)->size(); i++)
	{
	    *(buf++) = value;
	}
    return(*this);
}

VISImageRGBA& VISImageRGBA::operator=(const VISIm& from)
{
    if (from.type() == type())
	assign(*((VISImageRGBA*)(&from)));
    else
#ifdef AUTO_IMAGE_CONVERSION
	switch(from.type())
	    {
	      case BYTE:
		assignImage((*(const VISImage<byte>*)(&from)), *this);
		break;
	      case INT:
		assignImage((*(const VISImage<int>*)(&from)), 
			 *((VISImage<rgba>*)(this)));
		break;
	      case FLOAT:
		assignImage((*(const VISImage<float>*)(&from)), 
		       *((VISImage<rgba>*)(this)));
		break;
	      default:
		WARN("unrecognized image type conversion");
	    }
#else
    WARN("automatic image conversion must be set as compiler option");
#endif
    return(*this);
}


/* these are friend functions */

VISImageRGBA operator*(byte value, VISImageRGBA& a)
{
    VISImageRGBA image_return = a.createToSize();
    rgba* buf = (image_return.repRef(0))->bufferRef();
    const rgba* buf_from = (a.rep(0))->buffer();
    for (int i = 0; i < image_return.rep(0)->size(); i++)
	{
	    *(buf++) = *(buf_from++)*value;
	}
    return(image_return);
}

VISImageRGBA operator*(VISImageRGBA& a, byte value)
{
    VISImageRGBA image_return = a.createToSize();
    rgba* buf = (image_return.repRef())->bufferRef();
    const rgba* buf_from = (a.rep(0))->buffer();
    for (int i = 0; i < image_return.rep(0)->size(); i++)
	{
	    *(buf++) = *(buf_from++)*value;
	}
    return(image_return);
}


VISImageRGBA operator*(float value, const VISImageRGBA& a)
{
    return(a.scale(value));
}

VISImageRGBA operator*(const VISImageRGBA& a, float value)
{
    return(a.scale(value));
}
