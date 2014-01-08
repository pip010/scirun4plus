// ******************************************************************
// VISPack. Copyright (c) 1994-2000 Ross Whitaker rtw@utk.edu       *
// For conditions of distribution and use, see the file LICENSE.txt *
// accompanying this distribution.                                  *
// ******************************************************************

// $Id: vol.cxx,v 1.1.1.1 2003/02/12 16:51:54 whitaker Exp $

/* *	sccsid "@(#) image.C     2.0 2/16/94" */

#include "vol/vol.h"

void VISVol::unref(VISRep*& r)
{
    if (r)
	{
//	    printf("got an unref on %x with count %d\n", 
//		   (unsigned)r, r->ref_count());
	    r->deref();
	    if (r->ref_count() < 1)
		{
		    delete r;
		    r = NULL;
		}
	}
}

void VISVol::ref(const VISRep* r)
{
    if (r)
	{
//	    printf("got a ref on %x with count %d\n", 
//		   (unsigned)r, r->ref_count());
	    VISRep* r_ref = (VISRep*)r;
	    r_ref->ref();
	}
}

VISVol::~VISVol()
{
unref(_rep);
}
