/* *	sccsid "@(#) image.C     2.0 2/16/94" */

#include "image/im.h"

VISRep* __rep_tracker = NULL;

const VISIm nullImage;

void VISIm::unref(VISRep*& r)
{
    if (r)
	{
	    r->deref();
	    if (r->ref_count() < 1)
		{
		    delete r;
		    r = NULL;
		}
	}
}


void VISIm::ref(const VISRep* r)
{
    VISRep* r_ref = (VISRep*)r;
    r_ref->ref();
}


void VISIm::assign(const VISIm& from)
{
    if (!(&from == this))
	{
	    _type = from._type;
	    rechannel(from.channels());
	    for (int i = 0; i < _channels; i++)
		ref(_rep[i] = from._rep[i]);
	}
}

//
// tells whether or not all of the reps of two images are
// thus same (i.e. sufficient but not necessary condition
// for all of the data being same).
//

jpegBoolean VISIm::repsEqual(const VISIm& im) const
{
    jpegBoolean equal = TRUE;
    if (channels() == im.channels())
	{
	    for (int i = 0; (i < channels())&&equal; i++)
		if (rep(i) != im.rep(i))
		    equal = FALSE;
	}
    return(equal);
}


VISIm::~VISIm()
{
    if (_rep)
	{
	    for (int i = 0; i < _channels; i++)
		unref(_rep[i]);
	    delete[] _rep;
	    _rep = NULL;
	}
}


void VISIm::rechannel(unsigned int ch)
{
    int i;
    
    for (i = 0; i < _channels; i++)
	{
	    unref(_rep[i]);
	}

    if (ch != _channels)
	{
	    if (_rep)
		delete[] _rep;
	    _channels = ch;
	    if (_channels > 0)
		_rep = new VISRep*[_channels];
	    else 
		_rep = NULL;
	}

    for (i = 0; i < _channels; i++)
	_rep[i] = NULL;
}

void VISIm::initialize(VISRep** the_rep,  unsigned int ch)
{
    _channels = ch;
    _rep = new VISRep*[ch];
    for (int i = 0; i < channels(); i++)
	{
	    _rep[i] = NULL;
	    putRep(the_rep[i], i);
	}
}

void VISIm::initialize(unsigned int ch)
{
    _channels = ch;
    _rep = new VISRep*[_channels];
    for (int i = 0; i < channels(); i++)
	{
	    _rep[i] = NULL;
	}
}

void VISIm::putRep(const VISRep* r)
{
    unref(_rep[0]);
    ref(_rep[0] = (VISRep*)r);
} 

void VISIm::putRep(const VISRep* r, unsigned int ch)
{
    unref(_rep[ch]);
    ref(_rep[ch] = (VISRep*)r);
}




