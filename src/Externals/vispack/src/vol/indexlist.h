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
class IrisPair
{
  public:
    IrisPair()	{ a(0); b(0); }
    IrisPair(T aa, T bb) { a(aa); b(bb); }
    IrisPair(const T *p) 	{ a(p[0]);  b(p[1]); }
    IrisPair(const IrisPair<T>& p)	{ operator=(p); }
    
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
	    ERROR("IrisPair<T>:Range check failure");
#endif
	return _d[i];
    }
    const T&	operator[](const unsigned i) const {
#ifdef DEBUG
	if (i > 1)
	    ERROR("IrisPair<T>:Range check failure");
#endif
	return _d[i];
    }
    
    boolean operator==(const IrisPair<T>& p) const {
	return (p.a() == a() && p.b() == b());
    }
    boolean operator!=(const IrisPair<T>& p) const	{ return !operator==(p); }

    IrisPair<T>& operator*=(T s) { _d[0] *= s;  _d[1] *= s; return *this; }
    IrisPair<T>& operator/=(T s) { _d[0] /= s;  _d[1] /= s; return *this; }
    IrisPair<T>& operator+=(const IrisPair<T>& p)	{ _d[0] += p.a(); _d[1] += p.b();
					  return *this; }
    IrisPair<T>& operator-=(const IrisPair<T>& p)	{ _d[0] -= p.a(); _d[1] -= p.b();
					  return *this; }
    // misc

    IrisPair<T>& operator=(const IrisPair<T>& p) { a(p.a()); b(p.b()); return *this; }
    IrisPair<T>& operator=(const T* p){ a(p[0]); b(p[1]); return *this; }

    virtual void  print(ostream& ostr = cout) const
    {
	// { X, Y };
	ostr << "[ " << form("%3.2f", (float)a())
	     << ", " << form("%3.2f", (float)b())
	     << " ]";
    }

  protected:
    T  _d[2];
};



template < class T >
class IrisTriple
{
  public:
    IrisTriple()	{ a(0); b(0); c(0);}
    IrisTriple(T aa, T bb, T cc) { a(aa); b(bb); c(cc);}
    IrisTriple(const T *p) 	{ a(p[0]);  b(p[1]); c(p[2]); }
    IrisTriple(const IrisTriple<T>& p)	{ operator=(p); }
    
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
	    ERROR("IrisTriple:Range check failure");
#endif
	return _d[i];
    }
    const T&	operator[](const unsigned i) const {
#ifdef DEBUG
	if (i > 2)
	    ERROR("IrisTriple<T>:Range check failure");
#endif
	return _d[i];
    }
    
    boolean operator==(const IrisTriple<T>& p) const {
	return (p.a() == a() && p.b() == b() && p.c() == c());
    }
    boolean operator!=(const IrisTriple<T>& p) const	
    { return !operator==(p); }

    IrisTriple<T>& operator*=(T s) 
    { 
	_d[0] *= s;  _d[1] *= s; _d[2] *= s;
	return *this; 
    }
    IrisTriple<T>& operator/=(T s) 
    { 
	_d[0] /= s;  _d[1] /= s; _d[2] /= s;
	return *this; 
    }

    IrisTriple<T>& operator+=(const IrisTriple<T>& p)	
    { 
	_d[0] += p.a(); _d[1] += p.b(); _d[2] += p.c();
	return *this; 
    }
    
    IrisTriple<T>& operator-=(const IrisTriple<T>& p)	
    { 
	_d[0] -= p.a(); _d[1] -= p.b(); _d[2] -= p.c();
	return *this; 
    }
    // misc

    IrisTriple<T>& operator=(const IrisTriple<T>& p) 
    { 
	a(p.a()); b(p.b());  c(p.c());
	return *this; 
    }

    IrisTriple<T>& operator=(const T* p)
    { 
	a(p[0]); b(p[1]); c(p[2]);
	return *this; 
    }

    virtual void  print(ostream& ostr = cout) const
    {
	// { X, Y };
	ostr << "[ " << form("%3.2f", (float)a())
	     << ", " << form("%3.2f", (float)b())
	     << ", " << form("%3.2f", (float)c())
	     << " ]";
    }

  protected:
    T  _d[3];
};



class IrisImIndex: public IrisPair<unsigned>
{
  public:
    IrisImIndex() { a(0); b(0); }
    IrisImIndex(unsigned aa, unsigned bb) { a(aa); b(bb); }
    IrisImIndex(const unsigned *p) 	{ a(p[0]);  b(p[1]); }
    IrisImIndex(const IrisPair<unsigned>& p)	
    { IrisPair<unsigned>::operator=(p);}
    
    IrisImIndex& operator=(const IrisPair<unsigned>& p) 
    {
	IrisPair<unsigned>::operator=(p);
	return(*this);
    }
    IrisImIndex& operator=(const unsigned* p)
    {
	IrisPair<unsigned>::operator=(p);
	return(*this);
    }

    
};


class IrisVolIndex: public IrisTriple<unsigned>
{
  public:
    IrisVolIndex():IrisTriple<unsigned>() {}
    IrisVolIndex(unsigned aa, unsigned bb, unsigned cc)
	:IrisTriple<unsigned>(aa, bb, cc)
    {}
    IrisVolIndex(const unsigned *p)
    	:IrisTriple<unsigned>(p)
    {}
    
    IrisVolIndex(const IrisTriple<unsigned>& p)	
    { IrisTriple<unsigned>::operator=(p);}
    
    IrisVolIndex& operator=(const IrisTriple<unsigned>& p) 
    {
	IrisTriple<unsigned>::operator=(p);
	return(*this);
    }

    IrisVolIndex& operator=(const IrisVolIndex& p) 
    {
	IrisTriple<unsigned>::operator=(p);
	return(*this);
    }


    IrisVolIndex& operator=(const unsigned* p)
    {
	IrisTriple<unsigned>::operator=(p);
	return(*this);
    }
    
};


class IrisImIndexList: public List<IrisImIndex>
{

  public:
    IrisImIndexList()
    {
	_current_element = NULL;
    }
    
    void clean();
    boolean valid() {return(!(_current_element == NULL));}
    boolean reset() {return(!((_current_element=head())== NULL));}
    boolean stepBack();
    boolean stepForward();

    IrisImIndex atCurrent() const
    {
	return(_current_element->data());
    }

    void atCurrent(const IrisImIndex& other)
    {
	_current_element->data(other);
    }
    
    boolean removeCurrent();
    
  protected:
    Link<IrisImIndex> *_current_element;
    
};


class IrisVolIndexList: public List<IrisVolIndex>
{

  public:
    IrisVolIndexList()
    {
	_current_element = NULL;
    }
    
    void clean();
    boolean valid() {return(!(_current_element == NULL));}
    boolean reset() {return(!((_current_element=head())== NULL));}
    boolean stepBack();
    boolean stepForward();

    IrisVolIndex atCurrent() const
    {
	return(_current_element->data());
    }

    void atCurrent(const IrisVolIndex& other)
    {
	_current_element->data(other);
    }


    void copy(const IrisVolIndexList& other);
    
    boolean removeCurrent();
    
  protected:
    Link<IrisVolIndex> *_current_element;
    
};




#endif
