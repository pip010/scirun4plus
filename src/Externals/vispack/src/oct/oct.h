// this is the foundation for an nth (quad or oct) tree object
// it has a the ability to keep track of location withing a tree 
// and return an object at a specific location

// Ross Whitaker 4-5-96
//

#ifndef iris_oct_h
#define iris_oct_h


#include <stdio.h>
#include "util/defs.h"
#include "vol/volume.h"


#define __MAX_LEVELS (17)

static unsigned __LEVEL_FACTOR[17]
 = 
{
1,
1,
2,
4,
8,
16,
32,
64,
128,
256,
512,
1024,
2048,
4096,
8192,
16384,
32768
};



//************************************************************
// 
template <class T>
class VISOct
{
  protected:
    VISOct* _children[8];
    VISOct* _parent; 
    T _value;
    boolean _leaf;

    virtual void createChildren();

// you have to know where you are in the tree in order to do indexing
    int _level;
//    int _level_factor;

    void level(int l); 

    unsigned height() const {return(size());}
    unsigned width() const {return(size());}
    unsigned depth() const {return(size());}

    unsigned index(unsigned i) const
    {
	return(i/__LEVEL_FACTOR[_level]);
    }

    unsigned index(unsigned x, unsigned y, unsigned z) const
    {
	return(index(x) + 2*(index(y) + 2*index(z)));
    }

    unsigned remainder(unsigned i) const
    {
	return(i%__LEVEL_FACTOR[_level]);
    }

  public:
    VISOct(const T& v, VISOct<T>* my_parent, int l);
    VISOct();
    VISOct(int l);
    VISOct(const T& v, int l);

//
//    VISOct(VISOct* _children[8]);
//    VISOct(VISOct* _children[8], int l);
//

    int level() const {return(_level);}
    boolean isLeaf() const {return(_leaf);}

    const VISOct* child(unsigned x, unsigned y, unsigned z) const
    {
	return(child(index(x, y, z)));
    }

    const VISOct* parent() {return _parent;}
    
    VISOct* childRef(unsigned x, unsigned y, unsigned z) const
    {
	return(_children[index(x, y, z)]);
    }

    const T& value(unsigned x, unsigned y, unsigned z) const;
    void value(const T& v) {_value = v;}

    const VISOct* leaf(unsigned x, unsigned y, unsigned z) const;
    const VISOct* leaf(unsigned x, unsigned y, unsigned z, unsigned l) const;

    VISOct* leafRef(unsigned x, unsigned y, unsigned z);
    VISOct* leafRef(unsigned x, unsigned y, unsigned z, unsigned l);
    VISOct* leafRefCreate(unsigned x, unsigned y, unsigned z, unsigned l);

    const VISOct* child(unsigned num) const
    {
#ifdef BOUNDS_CHECKING_OCT
	if (num >= 8)
	    return(NULL);
	else
#endif
	    return(_children[num]);
    }


    VISOct* childRef(unsigned num) 
    {
#ifdef BOUNDS_CHECKING_OCT
	if (num >= 8)
	    return(NULL);
	else
#endif
	    return(_children[num]);
    }

    const T& value() const
    { 
	    return(_value);
    }

    T& valueRef()   
    { 
	    return(_value);
    }

    virtual ~VISOct();

    const VISOct& operator=(const VISOct& from)
	{
	    assign(&from);
	    return(*this);
	}

    const VISOct& operator=(const T& v)
	{
	    setValue(v, 0, 0, 0, level());
	    return(*this);
	}

    const VISOct& assignVol(const VISVolume<T>& vol, T padding);
    const VISOct& operator=(const VISVolume<T>& vol)
    { 
	return(assignVol(vol, (T)0)); 
    }


    VISOct(const VISOct& from)
	{
	    assign(&from);
	}

    VISOct(const VISOct* from)
	{
	    if (from != NULL)
		assign(from);
	}

//
// These cause problems with operator= and they don't assign levels
//
//    VISOct(const T* val_in)
//	{
//	    for (int i = 0; i < 8; i++)
//		_children[i] = NULL;
//
//	    _value = T(*val_in);
//	}
//
//    VISOct(const T& val_in)
//	{
//	    for (int i = 0; i < 8; i++)
//		_children[i] = NULL;
//
//	    _value = T(val_in);
//	}

    void assign(const VISOct* from)
	{
	    deleteChildren();
	    for (int i = 0; i < 8; i++)
		if (from != NULL)
		    _children[i] = new VISOct(from->_children[i]);
		else
		    _children[i] = NULL;

	    _value = from->_value;
	    _level = from->_level;
	}


// allows a node to consolidate it's children
    boolean childChanged();
// propogates a change up the tree
    void resolveChange();

    virtual void deleteChildren();
    T deleteChildrenAverage();
    T average();

// consolidates all the nodes below in the tree.
    boolean prune();

    void setValue(const T& v, unsigned x, unsigned y, unsigned z);
    void setValue(const T& v, unsigned x, unsigned y, 
		       unsigned z, int l);

// gets the value associate with a certain level
    const T& value(unsigned x, unsigned y, 
		   unsigned z, int l) const;
    
    unsigned size() const {return(2*__LEVEL_FACTOR[_level]);}
    void printData() const;
    void printDataLevels() const;

    virtual int promote(int d);

    virtual void print() const
    {
	printf("Oct:print()\n");
    }
};



// ************************************************************/
// this is the version that works with pointers to values
// Rules of the OctPtr tree:
//    The Oct works with pointers to values and other nodes
//    It is assumed that any pointer given to Oct is cleaned up
//    by Oct.
// 
template <class T>
class VISOctPtr
{
  protected:
//
// changed this to make children only if there is no value
//
    VISOctPtr **_children;
    VISOctPtr* _parent; 
    T* _value;

    virtual void createChildren();

// you have to know where you are in the tree in order to do indexing
    int _level;
//    int _level_factor;

    void level(int l); 

    unsigned index(unsigned i) const
    {
	return(i/__LEVEL_FACTOR[_level]);
    }

    unsigned index(unsigned x, unsigned y, unsigned z) const
    {
	return(index(x) + 2*(index(y) + 2*index(z)));
    }

    unsigned remainder(unsigned i) const
    {
	return(i%__LEVEL_FACTOR[_level]);
    }

  public:
    VISOctPtr();
    VISOctPtr(int l);
    VISOctPtr(const T& v, int l);
    VISOctPtr(T* v, int l);
    VISOctPtr(const T& val_in, VISOctPtr* my_parent, int l);
    VISOctPtr(T* val_in, const VISOctPtr* my_parent, int l);
//    VISOctPtr(VISOctPtr* _children[8]);
//    VISOctPtr(VISOctPtr* _children[8], int l);

    boolean isLeaf() const {return(_value != NULL);}
    int level() const {return(_level);}
    unsigned height() const {return(size());}
    unsigned width() const {return(size());}
    unsigned depth() const {return(size());}


    const VISOctPtr* child(unsigned x, unsigned y, unsigned z) const
    {
	return(child(index(x, y, z)));
    }

    const VISOctPtr* parent() {return _parent;}
    
    VISOctPtr* childRef(unsigned x, unsigned y, unsigned z) const
    {
	return(_children[index(x, y, z)]);
    }

    const T& value(unsigned x, unsigned y, unsigned z) const;
    void value(const T& v) {if (_value) *_value = v; else _value = new T(v);}

    const VISOctPtr* leaf(unsigned x, unsigned y, unsigned z) const;
    const VISOctPtr* leaf(unsigned x, unsigned y, unsigned z, 
			   unsigned l) const;

    VISOctPtr* leafRef(unsigned x, unsigned y, unsigned z);    
    VISOctPtr* leafRef(unsigned x, unsigned y, unsigned z, unsigned l);
    VISOctPtr* leafRefCreate(unsigned x, unsigned y, unsigned z, unsigned l);

    const VISOctPtr* child(unsigned num) const
    {
#ifdef BOUNDS_CHECKING_OCT
	if (num >= 8)
	    return(NULL);
	else
#endif
	    if (_children)
		return(_children[num]);
	else
	    return(NULL);
    }


    VISOctPtr* childRef(unsigned num) 
    {
#ifdef BOUNDS_CHECKING_OCT
	if (num >= 8)
	    return(NULL);
	else
#endif
	    if (_children)
		return(_children[num]);
	else
	    return(NULL);
    }

    const T* value() const
    { 
	    return(_value);
    }

    T* valueRef()   
    { 
	    return(_value);
    }

    virtual ~VISOctPtr();

    VISOctPtr operator=(const VISOctPtr& from)
	{
	    assign(&from);
	    return(*this);
	}

    const VISOctPtr& assignVol(const VISVolume<T>& vol, T padding);

    const VISOctPtr& operator=(const VISVolume<T>& vol)
    { 
	return(assignVol(vol, (T)0)); 
    }

    const VISOctPtr& operator=(const T& v)
	{
	    setValue(v, 0, 0, 0, level());
	    return(*this);
	}


    VISVolume<T> volume() const;


    VISOctPtr(const VISOctPtr& from)
	{
	    assign(&from);
	}

    VISOctPtr(const VISOctPtr* from)
	{
	    if (from != NULL)
		assign(from);
	}

    VISOctPtr(const T* val_in)
	{
//	    for (int i = 0; i < 8; i++)
//		_children[i] = NULL;

	    _children = NULL;
	    _value = val_in;
	}

    VISOctPtr(const T& val_in)
	{
//	    for (int i = 0; i < 8; i++)
//		_children[i] = NULL;

	    _children = NULL;
	    _value = new T(val_in);
	}

    void assign(const VISOctPtr* from)
	{
//	    for (int i = 0; i < 8; i++)
//		if (from != NULL)
//		    _children[i] = new VISOctPtr(from->_children[i]);
//		else
//		    _children[i] = NULL;

	    if (from->_value != NULL)
		_value = new T(*(from->_value));
	    else
		{
		    _value = NULL;
		    _children = new VISOctPtr**[8];
		    for (int i = 0; i < 8; i++)
			_children[i] = new VISOctPtr(from->_children[i]);
		}
	}

// allows a node to consolidate it's children
    boolean childChanged();
// propogates a change up the tree
    void resolveChange();

    virtual void deleteChildren();
    T deleteChildrenAverage();
    T average();

// consolidates all the nodes below in the tree.
    boolean prune();

    void setValue(const T& v, unsigned x, unsigned y, unsigned z);
    void setValue(const T& v, unsigned x, unsigned y, 
		       unsigned z, int l);

    const T* value(unsigned x, unsigned y, 
			   unsigned z, int l) const;

    unsigned size() const 
    {
	if (_level != 0)
	    return(2*__LEVEL_FACTOR[_level]);
	else
	    return(1);
    }
    
    void printData() const;
    void printDataLevels() const;

    virtual int promote(int d);

    virtual void print() const
    {
	printf("Oct:print()\n");
    }
};

#include "oct.txx"

/*******************************************************************/
// these are methods which stand on their own 



#endif








