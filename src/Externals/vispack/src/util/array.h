// ******************************************************************
// VISPack. Copyright (c) 1994-2000 Ross Whitaker rtw@utk.edu       *
// For conditions of distribution and use, see the file LICENSE.txt *
// accompanying this distribution.                                  *
// ******************************************************************

// $Id: array.h,v 1.1.1.1 2003/02/12 16:51:54 whitaker Exp $

// File:           array.h
// Author:         Ross T. Whitaker
// Institution:    The University of Tennessee, Knoxville
// Contents:       Several classes all supporing an array library.  The 'main'
//                 class normally used is the templated VISArray class.
// Log of changes: June 28, 1999 -- Added comments


#ifndef	util_array_h
#define	util_array_h

#include <iostream>
using namespace std;
#include "util/defs.h"
//#include "util/handler.h"
//#include "util/resource.h"
#include "util/utilExports.h"


// ---------------------------------------------------------------------------

class GfxTransform;
class GfxResource;

// ---------------------------------------------------------------------------
// Put things in here you want all arrays to understand

class VISArrayBase
{
  public:
    virtual ~VISArrayBase() {}
    virtual unsigned	n() const = 0;
};

// ---------------------------------------------------------------------------
// This is a general VISArray class -- it can grow if it wants to

template <class T> class VISArray : public VISArrayBase
{
// VISArray
// This class is the 'main' array class, meaning this is the
// class users will normally use as an array.

  public:
//****Constructors and Destructors, the usuall stuff
    VISArray(){_size = 0; _n = 0;_array = NULL; }//create an empty array
    // create an array that is a copy of another array
    VISArray(const VISArray& s){ _size = 0; _n = 0; _array = NULL; operator=(s); }
    // create an array of a particular size
    VISArray(unsigned size0) { _size = size0; _n = 0; _array = new T[_size]; }
    // These allow one to create arrays with up to four elements by default
    VISArray(T x0,T x1) {
	_size = 2; _n = 2; _array = new T[_size]; at(0) = x0; at(1) = x1; }
    VISArray(T x0,T x1,T x2) {
	_size = 3; _n = 3; _array = new T[_size]; at(0) = x0; at(1) = x1;
	at(2) = x2; }
    VISArray(T x0,T x1,T x2,T x3) {
	_size = 4; _n = 4; _array = new T[_size]; at(0) = x0; at(1) = x1;
	at(2) = x2; at(3) = x3; }
    virtual ~VISArray();// destructor

//****Finding the size of an array
    virtual unsigned	n() const;// return the number of elements in an array
    // WARNING:  Size does NOT return the number of elements in an array
    // Do NOT use size() in a loop to loop through the elements of the array
    // it returns the size of the array.  The size of an array begins
    // as 8 and is doubled when new space is needed for elements.  For example:
    // VISArray<float> A;      //size()=0   n()=0
    // A.appenditem(0.0f);  //size()=8   n()=1
    virtual unsigned	size() const;

//****Interacting with the elements of an array
    virtual void clearItem(unsigned);// clear item at index
    virtual void clearItems(unsigned i0,unsigned i1);//clear items between indicies
    virtual void clear();//clear all items in array and set size=0
    virtual void deleteItem(unsigned);
    virtual void deleteItems();
    //returns the item at index idx.  Use on the rhs of an operator=
    virtual const T& peek(unsigned idx) const;
    virtual const T& operator[](unsigned idx) const;//same as peek

    //returns a reference to the item at index idx.  This reference can
    //be used to change the item.  Use on the lhs of an operator=
    virtual T& poke(unsigned idx);
    //returns a reference to the item at index idx.  This is a const
    //reference so that the data can not be changed.
    virtual T& refAt(unsigned idx) const;
    //adds an item onto the end of the array
    virtual void appendItem(const T& t);
    //adds an item onto the begining of the array
    virtual void prependItem(const T& t);
    //inserts an item at index idx
    virtual void insertItemAt(const T& t,unsigned idx);
    //replaces an item in an array with a new item at index idxx
    virtual void replaceItemAtWith(unsigned idx, const T& t);
    //removes item from array at index idx
    virtual void removeItemAt(unsigned idx);
    //reverses the order of the array
    virtual void reverse();
    virtual void reverse(unsigned,unsigned=0);
    virtual void rowReverse(unsigned rows,unsigned cols);
    //this does a bubble sort (descending)
    void sort();

    virtual VISArray<T>* copy() const { return new VISArray<T>(*this); }

    // Assignment operators
     VISArray& operator=(const VISArray& s);
     virtual void sizeTo(unsigned new_size);
    //functions not used by the user
    virtual const T& itemAt(unsigned idx) const;
    virtual T& at(unsigned idx);

  protected:
    unsigned	_size;
    unsigned	_n;
    T*		_array;
    
     virtual void copyItem(T& dst,const T& src);
    
     virtual void grow();
     virtual void growTo(unsigned new_size);

    VISArray(const VISArray* s) { _size = 0; _n = 0; _array = NULL; operator=(s); }
     VISArray& operator=(const VISArray* s);
};



//#ifdef	DEPENDING
//#  include "util/array.c"
//#else
//#  ifdef  INLINE_TEMPLATES
//#    include "util/array.c"
//#  endif
//#endif

#ifndef MANUAL_INSTANTIATION
#include "array.txx"
#endif



#endif	/* util_array_h */
