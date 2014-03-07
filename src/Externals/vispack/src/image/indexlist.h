// ******************************************************************
// VISPack. Copyright (c) 1994-2000 Ross Whitaker rtw@utk.edu       *
// For conditions of distribution and use, see the file LICENSE.txt *
// accompanying this distribution.                                  *
// ******************************************************************

// $Id: indexlist.h,v 1.1.1.1 2003/02/12 16:51:52 whitaker Exp $

/*
****************************************************************
*  indexlist.h
*
* @(#)indexlist.h 1.1 94/04/26 11:00:35
*
****************************************************************
*/

#ifndef INDEXLIST_H
#define INDEXLIST_H

// -------------------------------------------------------------------------

#include <util/list.h>

// ---------------------------------------------------------------------------
// Pair of anything
// ---------------------------------------------------------------------------


template < class T >
class VISPair
{
  public:
    VISPair()	{ a(0); b(0); }
    VISPair(T aa, T bb) { a(aa); b(bb); }
    VISPair(const T *p) 	{ a(p[0]);  b(p[1]); }
    VISPair(const VISPair<T>& p)	{ operator=(p); }
    
    // provide access to indivual components.;
    // There method will need to be made virtual for proper sub-classing;
    T	a()	const	{ return _d[0]; }
    T	b()	const	{ return _d[1]; }
    
    void	a(T a)	{ _d[0] = a; }
    void	b(T b)	{ _d[1] = b; }
    
    // Allow indexing on point objects
    // Range checking is performed when debugging is enabled
    T&	at(const unsigned int i) {
#ifdef DEBUG
	if (i > 1)
	    ERROR("VISPair<T>:Range check failure");
#endif
	return _d[i];
    }
    const T&	operator[](const unsigned i) const {
#ifdef DEBUG
	if (i > 1)
	    ERROR("VISPair<T>:Range check failure");
#endif
	return _d[i];
    }
    
    boolean operator==(const VISPair<T>& p) const {
	return (p.a() == a() && p.b() == b());
    }
    boolean operator!=(const VISPair<T>& p) const	{ return !operator==(p); }

    VISPair<T>& operator*=(T s) { _d[0] *= s;  _d[1] *= s; return *this; }
    VISPair<T>& operator/=(T s) { _d[0] /= s;  _d[1] /= s; return *this; }
    VISPair<T>& operator+=(const VISPair<T>& p)	{ _d[0] += p.a(); _d[1] += p.b();
					  return *this; }
    VISPair<T>& operator-=(const VISPair<T>& p)	{ _d[0] -= p.a(); _d[1] -= p.b();
					  return *this; }
    // misc

    VISPair<T>& operator=(const VISPair<T>& p) { a(p.a()); b(p.b()); return *this; }
    VISPair<T>& operator=(const T* p){ a(p[0]); b(p[1]); return *this; }

    virtual void  print(ostream& ostr = cout) const
    {
	// { X, Y };
	ostr << "[ " << (float)a()
	     << ", " << (float)b()
	     << " ]";
    }

  protected:
    T  _d[2];
};




class VISImIndex: public VISPair<unsigned>
{
  public:
    VISImIndex() { a(0); b(0); }
    VISImIndex(unsigned aa, unsigned bb) { a(aa); b(bb); }
    VISImIndex(const unsigned *p) 	{ a(p[0]);  b(p[1]); }

    VISImIndex(const VISPair<unsigned>& p)	
    { VISPair<unsigned>::operator=(p);} 

    VISImIndex(const VISImIndex& p)	
    { VISPair<unsigned>::operator=(p);}
   
    VISImIndex& operator=(const VISPair<unsigned>& p) 
    {
	VISPair<unsigned>::operator=(p);
	return(*this);
    }

    VISImIndex& operator=(const VISImIndex& p) 
    {
	VISPair<unsigned>::operator=(p);
	return(*this);
    }

    VISImIndex& operator=(const unsigned* p)
    {
	VISPair<unsigned>::operator=(p);
	return(*this);
    }

    
};

class VISImIndexList: public VISList<VISImIndex>
{
  public:
    VISImIndexList()
    {
	_current_element = NULL;
    }
    
    void clean();
    boolean valid() {return(!(_current_element == NULL));}
    boolean reset() {return(!((_current_element=head())== NULL));}
    boolean stepBack();
    boolean stepForward();

    VISImIndex& atCurrent() const
    {
	return(_current_element->data());
    }

    VISImIndex itemAtCurrent() const
    {
	return(_current_element->data());
    }

    void atCurrent(const VISImIndex& other)
    {
	_current_element->data(other);
    }
    
    boolean removeCurrent();
    
  protected:
    Link<VISImIndex> *_current_element;
    
};





#endif
