// ******************************************************************
// VISPack. Copyright (c) 1994-2000 Ross Whitaker rtw@utk.edu       *
// For conditions of distribution and use, see the file LICENSE.txt *
// accompanying this distribution.                                  *
// ******************************************************************

// $Id: vol.h,v 1.1.1.1 2003/02/12 16:51:54 whitaker Exp $



// File:           vol.h
// Author:         Ross T. Whitaker
// Institution:    The University of Tennessee, Knoxville
// Contents:       VISVol:  A class that is a base class for VISVolume
//                 (see volume.h) 
// Log of changes: July 14, 1999 -- Added comments

#ifndef iris_vol_h
#define iris_vol_h


#include <stdio.h>
#include "image/im.h"




/************************************************************/
class VISVol
{
  protected:
    VISIm::Type _type;
    VISRep* _rep;

  public:

    VISVol()  {_type = VISIm::NONE; _rep = NULL;}
    VISIm::Type type() const {return(_type);}
    unsigned int width() const 
    {
	if (_rep == NULL)
	    return 0;
	else
	    return(_rep->width());	    
    }
    unsigned int height() const 
    {
	if (_rep == NULL)
	    return 0;
	else
	    return(_rep->height());	    
    }
    unsigned int depth() const 
    {
	if (_rep == NULL)
	    return 0;
	else
	    return(_rep->depth());	    
    }

    void unref(VISRep*& r);
    void ref(const VISRep* r);

//   Volplement a deep copy later on    
    VISVol operator=(const VISVol& from)
	{
	    assign(from);
	    return(*this);
	}

    VISVol(const VISVol& from)
	{
	    _rep = NULL;
	    assign(from);
	}

    void assign(const VISVol& from)
	{
//	    printf("got vol assign\n");
	    if (&from != this)
		{
		    _type = from._type;
		    unref(_rep);
		    ref(_rep = from._rep);
		}
	}

    const VISRep* rep() const {return(_rep);}

    virtual ~VISVol();

    virtual boolean isValid() const
    {
	return(_type != VISIm::NONE);
    }

    virtual void print() const
    {
	printf("VISVol type %d\n", type());
    }
};


/*******************************************************************/
// these are methods which stand on their own 

#endif









