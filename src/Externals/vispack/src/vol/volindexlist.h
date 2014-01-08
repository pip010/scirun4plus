// ******************************************************************
// VISPack. Copyright (c) 1994-2000 Ross Whitaker rtw@utk.edu       *
// For conditions of distribution and use, see the file LICENSE.txt *
// accompanying this distribution.                                  *
// ******************************************************************

// $Id: volindexlist.h,v 1.1.1.1 2003/02/12 16:51:54 whitaker Exp $

/*
****************************************************************
*  indexlist.h
*
* @(#)indexlist.h 1.1 94/04/26 11:00:35
*
****************************************************************
*/

#ifndef vis_volindexlist_h
#define vis_volindexlist_h

// -------------------------------------------------------------------------

#include <util/list.h>

// ---------------------------------------------------------------------------
// triple of anything
// ---------------------------------------------------------------------------

template < class T >
class VISTriple
{
  public:
    VISTriple()	{ a(0); b(0); c(0);}
    VISTriple(T aa, T bb, T cc) { a(aa); b(bb); c(cc);}
    VISTriple(const T *p) 	{ a(p[0]);  b(p[1]); c(p[2]); }
    VISTriple(const VISTriple<T>& p)	{ operator=(p); }
    
    // provide access to indivual components.;
    // There method will need to be made virtual for proper sub-classing;
    T	a()	const	{ return _d[0]; }
    T	b()	const	{ return _d[1]; }
    T	c()	const	{ return _d[2]; }
    
    void	a(T a)	{ _d[0] = a; }
    void	b(T b)	{ _d[1] = b; }
    void	c(T c)	{ _d[2] = c; }
    
    // Allow indexing on point objects
    // Range checking is performed when debugging is enabled
    T&	at(const unsigned int i) {
#ifdef DEBUG
	if (i > 2)
	    ERROR("VISTriple:Range check failure");
#endif
	return _d[i];
    }
    const T&	operator[](const unsigned i) const {
#ifdef DEBUG
	if (i > 2)
	    ERROR("VISTriple<T>:Range check failure");
#endif
	return _d[i];
    }
    
    boolean operator==(const VISTriple<T>& p) const {
	return (p.a() == a() && p.b() == b() && p.c() == c());
    }
    boolean operator!=(const VISTriple<T>& p) const	
    { return !operator==(p); }

    VISTriple<T>& operator*=(T s) 
    { 
	_d[0] *= s;  _d[1] *= s; _d[2] *= s;
	return *this; 
    }
    VISTriple<T>& operator/=(T s) 
    { 
	_d[0] /= s;  _d[1] /= s; _d[2] /= s;
	return *this; 
    }

    VISTriple<T>& operator+=(const VISTriple<T>& p)	
    { 
	_d[0] += p.a(); _d[1] += p.b(); _d[2] += p.c();
	return *this; 
    }
    
    VISTriple<T>& operator-=(const VISTriple<T>& p)	
    { 
	_d[0] -= p.a(); _d[1] -= p.b(); _d[2] -= p.c();
	return *this; 
    }
    // misc

    VISTriple<T>& operator=(const VISTriple<T>& p) 
    { 
	a(p.a()); b(p.b());  c(p.c());
	return *this; 
    }

    VISTriple<T>& operator=(const T* p)
    { 
	a(p[0]); b(p[1]); c(p[2]);
	return *this; 
    }

    virtual void  print(ostream& ostr = cout) const
    {
	// { X, Y };
	ostr << "[ " <<  (float)a() << ", " << (float)b()
	     << ", " << (float)c()
	     << " ]";
    }

  protected:
    T  _d[3];
};



class VISVolIndex: public VISTriple<unsigned>
{
  public:
    VISVolIndex():VISTriple<unsigned>() {}
    VISVolIndex(unsigned aa, unsigned bb, unsigned cc)
	:VISTriple<unsigned>(aa, bb, cc)
    {}
    VISVolIndex(const unsigned *p)
    	:VISTriple<unsigned>(p)
    {}
    
    VISVolIndex(const VISTriple<unsigned>& p)	
    { VISTriple<unsigned>::operator=(p);}
    
    VISVolIndex& operator=(const VISTriple<unsigned>& p) 
    {
	VISTriple<unsigned>::operator=(p);
	return(*this);
    }

    VISVolIndex& operator=(const VISVolIndex& p) 
    {
	VISTriple<unsigned>::operator=(p);
	return(*this);
    }


    VISVolIndex& operator=(const unsigned* p)
    {
	VISTriple<unsigned>::operator=(p);
	return(*this);
    }
    
};


class VISVolIndexVISList: public VISList<VISVolIndex>
{

  public:
    VISVolIndexVISList()
    {
	_current_element = NULL;
    }
    
    void clean();
    boolean valid() {return(!(_current_element == NULL));}
    boolean reset() {return(!((_current_element=head())== NULL));}
    boolean stepBack();
    boolean stepForward();

    VISVolIndex& atCurrent() const
    {
	return(_current_element->data());
    }

    VISVolIndex itemAtCurrent() const
    {
	return(_current_element->data());
    }

    void atCurrent(const VISVolIndex& other)
    {
	_current_element->data(other);
    }


    void copy(const VISVolIndexVISList& other);
    
    boolean removeCurrent();
    
  protected:
    Link<VISVolIndex> *_current_element;
    
};



class VolIndexValue: public VISTriple<unsigned>
{
    float _value;

  public:
    VolIndexValue():VISTriple<unsigned>() {}

    VolIndexValue(unsigned aa, unsigned bb, unsigned cc, float vv)
	:VISTriple<unsigned>(aa, bb, cc)
    {_value = vv;}

    VolIndexValue(unsigned aa, unsigned bb, unsigned cc)
	:VISTriple<unsigned>(aa, bb, cc)
    {;}

    VolIndexValue(const unsigned *p, float vv)
    	:VISTriple<unsigned>(p)
    {_value = vv;}
    
    VolIndexValue(const VISTriple<unsigned>& p, float vv)	
    { VISTriple<unsigned>::operator=(p);
    _value = vv;
    }

    VolIndexValue(const VISTriple<unsigned>& p)	
    { 
	VISTriple<unsigned>::operator=(p);
	_value = 0.0;
    }

    VolIndexValue(const VolIndexValue& other)	
    { 
	VolIndexValue::operator=(other);
    }
    
    VolIndexValue& operator=(const VolIndexValue& p) 
    {
	VISTriple<unsigned>::operator=(p);
	_value = p.value();
	return(*this);
    }

    float value() const {return(_value);}
    void value(float v) {_value = v;}
};


class VolIndexValueVISList: public VISList<VolIndexValue>
{

  public:
    VolIndexValueVISList()
    {
	_current_element = NULL;
    }
    
    void clean();
    boolean valid() {return(!(_current_element == NULL));}
    boolean reset() {return(!((_current_element=head())== NULL));}
    boolean stepBack();
    boolean stepForward();

    VolIndexValue& atCurrent() 
    {
	return(_current_element->data());
    }

    const VolIndexValue& itemAtCurrent() const
    {
	return(_current_element->data());
    }

    void atCurrent(const VolIndexValue& other)
    {
	_current_element->data(other);
    }


    void copy(const VolIndexValueVISList& other);
    
    boolean removeCurrent();
    
  protected:
    Link<VolIndexValue> *_current_element;
    
};



#endif
