// ******************************************************************
// VISPack. Copyright (c) 1994-2000 Ross Whitaker rtw@utk.edu       *
// For conditions of distribution and use, see the file LICENSE.txt *
// accompanying this distribution.                                  *
// ******************************************************************

// $Id: im.h,v 1.2 2003/02/26 03:36:44 whitaker Exp $

// File:           im.h
// Author:         Ross T. Whitaker
// Institution:    The University of Tennessee, Knoxville
// Contents:       Contains base classes for the VISImage and VISVol classes.
// Log of changes: June 30, 1999 -- Added comments

#ifndef iris_im_h
#define iris_im_h


#include <stdio.h>

#ifdef DEBUG
class VISRep;
extern VISRep *__rep_tracker;
#endif

typedef signed char schar;

#ifdef MAX
#undef MAX
#endif
#define MAX(a, b)      ((a)>(b)?(a):(b))

#ifdef MIN
#undef MIN
#endif
#define MIN(a, b)      ((a)<(b)?(a):(b))

#ifdef ABS
#undef ABS
#endif
#define ABS(a)      ((a)<(0)?(-a):(a))


#include "util/defs.h"


//
//
inline void ERROR(const char* a)
{
    printf("Image Lib Error: %s \n", a);
}


//
//
inline void WARN(const char* a)
{
    printf("Image Lib Warning: %s \n", a);
}


// This class is the abstract buffer for both VISImage and VISVol classes
class VISRep
{
  protected:
    int _ref_count;
  public:
    virtual ~VISRep(){}
    VISRep() {_ref_count = 0;}
    void deref() 
    {
#ifdef DEBUG 
	if (this == __rep_tracker)
	    printf("deref on rep %d with ref count %d\n", 
		   (int)this, _ref_count);
#endif
	_ref_count--;
    }
    void ref() 
    {
#ifdef DEBUG 
	if (this == __rep_tracker)
	    printf("ref on rep %d with ref count %d\n", 
		   (int)this, _ref_count);
#endif
	_ref_count++;
    }

    int ref_count() const {return(_ref_count);}
    virtual unsigned int width() const {return(0);}
    virtual unsigned int height() const {return(0);}

// this is so the same thing can be used for volumes
    virtual unsigned int depth() const {return(0);}
};


/************************************************************/
class VISIm
{
  public:
    typedef enum {OTHER, SCHAR, BYTE, FLOAT, INT, RGBA, UNSIGNED, SHORT, USHORT, NONE} Type;


//    friend
//	void assignImage(const VISIm& from, VISIm& to);
    
  protected:
    Type _type;
    VISRep** _rep;
    unsigned int _channels;

    VISIm(unsigned int ch)
    {
	_rep = NULL;
	initialize(ch);
    }

  public:

    VISIm()  {_type = NONE; _rep = NULL; _channels = 0;}
    Type type() const {return(_type);}
    unsigned int channels() const {return(_channels);}
    unsigned int width() const 
    {
	return((_channels > 0)?(rep(0)->width()):0);
    }
    unsigned int height() const 
    {
	return((_channels > 0)?(rep(0)->height()):0);
    }

    VISIm& operator=(const VISIm& from)
	{
	    assign(from);
	    return(*this);
	}

    VISIm(const VISIm& from)
    {
	_type = NONE; _rep = NULL; _channels = 0;
	assign(from);
    }
    
//    VISIm(VISIm& from)
//    {
//	_type = NONE; _rep = NULL; _channels = 0;	
//	assign(from);
//    }

    void assign(const VISIm& from);

    virtual ~VISIm();

    boolean isValid() const
    {
	return((_type != NONE)&&(width()>0)&&(height()>0));
    }

//
// tells whether or not all of the reps of two images are
// thus same (i.e. sufficient but not necessary condition
// for all of the data being same).
//

    boolean repsEqual(const VISIm& im) const;
    boolean operator==(const VISIm& im) const { return(repsEqual(im)); }

    const VISRep* rep(unsigned int ch) const
    {
#ifdef DYNAMIC_CHECKING
	if (!(ch < _channels))
	    {
		WARN("VISImage<>: channel(int ch) - channel out of range");
		return(NULL);
	    }
#endif 
	return(_rep[ch]);
    }

//
// for now this doesn't make sense, should this be available ?
//    VISRep* repRef(unsigned int ch);

    virtual void print() const
    {
	printf("VISIm type %d\n", type());
    }

    void rechannel(unsigned int);
    void initialize(VISRep** the_rep,  unsigned int ch);
    void initialize(unsigned int ch);
    void unref(VISRep*& r);
    void ref(const VISRep* r);
    void putRep(const VISRep* r);
    void putRep(const VISRep* r, unsigned int ch);
    VISRep** repVISArray() const {return(_rep);}
};

extern const VISIm nullImage;


/*******************************************************************/
// these are methods which stand on their own 

#endif






