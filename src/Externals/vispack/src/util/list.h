// ******************************************************************
// VISPack. Copyright (c) 1994-2000 Ross Whitaker rtw@utk.edu       *
// For conditions of distribution and use, see the file LICENSE.txt *
// accompanying this distribution.                                  *
// ******************************************************************

// $Id: list.h,v 1.2 2005/12/01 03:32:51 miriah Exp $

// sccsid  "@(#)	list.h	2.4 25 Aug 1994"

#ifndef	util_list_h
#define	util_list_h

#include <iostream>
using namespace std;
#include "util/defs.h"
#include "util/utilExports.h"

// ---------------------------------------------------------------------------
// Put things in here you want all list links to understand

class util_SHARE LinkBase
{
  protected:
    LinkBase	*_next,*_prev;
  public:
    LinkBase()		{ _next = NULL; _prev = NULL; }
    virtual ~LinkBase()	{ deleteItem(); }

    LinkBase *next() const	{ return _next; }
    LinkBase *prev() const	{ return _prev; }

    virtual void clearItem() {}
    virtual void deleteItem() {}

    void append(LinkBase*);
    void prepend(LinkBase*);
    void removeSelf();
};

// ---------------------------------------------------------------------------
// Put things in here you want all lists to understand

class util_SHARE VISListBase
{
  protected:
    unsigned	_n;
    LinkBase	*_head,*_tail;
  public:
    VISListBase()	{ _n = 0; _head = NULL; _tail = NULL; }
    virtual ~VISListBase()	{ deleteItems(); }
    unsigned n() const	{ return _n; }

    LinkBase *appendItem(LinkBase*);
    LinkBase *prependItem(LinkBase*);
    LinkBase *insertItemAt(LinkBase* link,LinkBase* at) { return insertItemAfter(link,at); }
    LinkBase *insertItemAfter(LinkBase*,LinkBase*);
    LinkBase *insertItemBefore(LinkBase*,LinkBase*);
    void removeItem(LinkBase*);

    LinkBase *head() const { return _head; }
    LinkBase *tail() const { return _tail; }
    
    virtual void clear();
    virtual void deleteItems();
};

// ---------------------------------------------------------------------------
// This is a general Link class

template <class T> class Link : public LinkBase
{
  protected:
    T		_data;
  public:
    Link()			{ }
    Link(T data)		{ _data = data; }
    Link(const Link& link)	{ _data = link._data(); }
    virtual ~Link()		{ deleteItem(); }
    
    Link<T> *next() const	{ return (Link<T>*)LinkBase::next(); }
    Link<T> *prev() const	{ return (Link<T>*)LinkBase::prev(); }

    virtual void clearItem() {}
    virtual void deleteItem() {}

    virtual void append(Link<T>* link)	{ LinkBase::append((LinkBase*)link); }
    virtual void prepend(Link<T>* link)	{ LinkBase::prepend((LinkBase*)link); }

    virtual T& data() { return _data; }
    virtual void data(T data)	{ _data = data; }
};

// ---------------------------------------------------------------------------
// This is a general Link class

template <class T> class PtrLink : public Link<T>
{
  public:
    PtrLink()			{ PtrLink<T>::_data = NULL; }// **mdm**
    PtrLink(T data) : Link<T>(data) {}
    PtrLink(const PtrLink& link) { PtrLink<T>::_data = link.data(); }  // **mdm**
    ~PtrLink()			{ deleteItem(); }

    PtrLink<T> *next() const	{ return (PtrLink<T>*)LinkBase::next(); }
    PtrLink<T> *prev() const	{ return (PtrLink<T>*)LinkBase::prev(); }
    
    void clearItem()	{ PtrLink<T>::_data = NULL; }  // **mdm**
    void deleteItem()	{ delete PtrLink<T>::_data; PtrLink<T>::_data = NULL; }// **mdm**
};

// ---------------------------------------------------------------------------
// This is a general linked-list class

template <class T> class VISList : public VISListBase
{
  protected:
  public:
    Link<T> *head() const { return (Link<T>*)VISListBase::head(); }
    Link<T> *tail() const { return (Link<T>*)VISListBase::tail(); }

    Link<T> *appendItem(T data) {
	return (Link<T>*)VISListBase::appendItem(new Link<T>(data)); }
    Link<T> *prependItem(T data) {
	return (Link<T>*)VISListBase::prependItem(new Link<T>(data)); }
    Link<T> *insertItemAt(T data,Link<T> *link) {
	return (Link<T>*)VISListBase::insertItemAt(new Link<T>(data),link); }
    Link<T> *insertItemBefore(T data,Link<T> *link) {
	return (Link<T>*)VISListBase::insertItemBefore(new Link<T>(data),link); }
    Link<T> *insertItemAfter(T data,Link<T> *link) {
	return (Link<T>*)VISListBase::insertItemAfter(new Link<T>(data),link); }
    void removeItem(Link<T> *link)	{ VISListBase::removeItem(link); }
};

// ---------------------------------------------------------------------------
// This is a general linked-list class for pointer objects

template <class T> class PtrVISList : public VISList<T>
{
  protected:
  public:
    virtual ~PtrVISList()	{ PtrVISList<T>::deleteItems(); }// **mdm**
    PtrLink<T> *head() const { return (PtrLink<T>*)VISListBase::head(); }
    PtrLink<T> *tail() const { return (PtrLink<T>*)VISListBase::tail(); }
    
    PtrLink<T> *appendItem(T data) {
	return (PtrLink<T>*)VISListBase::appendItem(new PtrLink<T>(data)); }
    PtrLink<T> *prependItem(T data) {
	return (PtrLink<T>*)VISListBase::prependItem(new PtrLink<T>(data)); }
    PtrLink<T> *insertItemAt(T data,PtrLink<T> *link) {
	return (PtrLink<T>*)VISListBase::insertItemAt(new PtrLink<T>(data),link); }
    PtrLink<T> *insertItemBefore(T data,PtrLink<T> *link) {
	return (PtrLink<T>*)VISListBase::insertItemBefore(new PtrLink<T>(data),link); }
    PtrLink<T> *insertItemAfter(T data,PtrLink<T> *link) {
	return (PtrLink<T>*)VISListBase::insertItemAfter(new PtrLink<T>(data),link); }
    void removeItem(PtrLink<T> *link)	{ VISListBase::removeItem(link); }
};

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

#endif	/* util_list_h */
