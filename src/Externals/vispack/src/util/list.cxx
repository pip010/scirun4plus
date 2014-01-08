// ******************************************************************
// VISPack. Copyright (c) 1994-2000 Ross Whitaker rtw@utk.edu       *
// For conditions of distribution and use, see the file LICENSE.txt *
// accompanying this distribution.                                  *
// ******************************************************************

// $Id: list.cxx,v 1.1.1.1 2003/02/12 16:51:54 whitaker Exp $

// sccsid  "@(#)	list.C	2.4 25 Aug 1994"

#include "util/list.h"

// ---------------------------------------------------------------------------

void LinkBase::append(LinkBase *link)
{
    if (link == NULL)
        return;
    link->_next = _next;
    if (_next != NULL)
        _next->_prev = link;
    _next = link;
    link->_prev = this;
}

// ---------------------------------------------------------------------------

void LinkBase::prepend(LinkBase *link)
{
    if (link == NULL)
        return;
    link->_next = this;
    if (_prev != NULL)
        _prev->_next = link;
    link->_prev = _prev;
    _prev = link;
}

// ---------------------------------------------------------------------------

void LinkBase::removeSelf()
{
    if (_prev != NULL) {
        _prev->_next = _next;
        
    }
    if (_next != NULL) {
        _next->_prev = _prev;
        
    }

    // here changed on 6/20/01 for pixmodel
    // for the need of parallel processing

    _prev = NULL; 
    _next = NULL;
}

// ---------------------------------------------------------------------------

LinkBase *VISListBase::appendItem(LinkBase* link)
{
    if (_tail == NULL) {
        _head = link;
        _tail = link;
    }
    else {
        _tail->append(link);
        _tail = link;
    }
    _n++;

    return link;
}

// ---------------------------------------------------------------------------

LinkBase *VISListBase::prependItem(LinkBase* link)
{
    if (_head == NULL) {
        _head = link;
        _tail = link;
    }
    else {
        _head->prepend(link);
        _head = link;
    }
    _n++;

    return link;
}

// ---------------------------------------------------------------------------

LinkBase *VISListBase::insertItemAfter(LinkBase *link,LinkBase *at)
{
    if (_tail == NULL) {
	_head = link;
	_tail = link;
    }
    else if (at == _tail) {
	_tail->append(link);
	_tail = link;
    }
    else if (at == NULL) {	// append
	_tail->append(link);
	_tail = link;
    }
    else {
	at->append(link);
    }
    _n++;

    return link;
}

// ---------------------------------------------------------------------------

LinkBase *VISListBase::insertItemBefore(LinkBase *link,LinkBase *at)
{
    if (_tail == NULL) {
	_head = link;
	_tail = link;
    }
    else if (at == _head) {
	_head->prepend(link);
	_head = link;
    }
    else if (at == NULL) {	// prepend
	_head->prepend(link);
	_head = link;
    }
    else {
	at->prepend(link);
    }
    _n++;

    return link;
}

// ---------------------------------------------------------------------------

void VISListBase::removeItem(LinkBase* link)
{
    if (link == _head) {
        _head = link->next();
    }
    if (link == _tail) {
        _tail = link->prev();
    }
    link->removeSelf();
    
    delete link;
    _n--;
}

// ---------------------------------------------------------------------------
// Clear before deleting

void VISListBase::clear()
{
    LinkBase *p,*q;

    for (p=head(); p!=NULL; q=p, p=p->next(), q->clearItem(), delete q);
    _head = NULL;
    _tail = NULL;
    _n = 0;
}

// ---------------------------------------------------------------------------
// delete(q) does a deleteItem()

void VISListBase::deleteItems()
{
    LinkBase *p,*q;

    for (p=head(); p!=NULL; q=p, p=p->next(), delete q);
    _head = NULL;
    _tail = NULL;
    _n = 0;
}

// ---------------------------------------------------------------------------
