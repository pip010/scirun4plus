// ******************************************************************
// VISPack. Copyright (c) 1994-2000 Ross Whitaker rtw@utk.edu       *
// For conditions of distribution and use, see the file LICENSE.txt *
// accompanying this distribution.                                  *
// ******************************************************************

// $Id: volindexlist.cxx,v 1.1.1.1 2003/02/12 16:51:54 whitaker Exp $

#include "vol/volindexlist.h"

void VISVolIndexVISList::clean()
{
    reset();
    while(valid())
	removeCurrent();
}


boolean VISVolIndexVISList::removeCurrent()
{
    Link<VISVolIndex> *element_tmp;
    if (valid())
	{
	    element_tmp = _current_element->next();
	    removeItem(_current_element);
	    _current_element = element_tmp;
	    return(TRUE);
	}
    else
	return(FALSE);
}
    
boolean VISVolIndexVISList::stepForward() 
{
	return((valid())&&
	       ((_current_element=_current_element->next())||TRUE));
    }

boolean  VISVolIndexVISList::stepBack()
{
    return((valid())&&
	   ((_current_element=_current_element->prev())||TRUE));
}


void VISVolIndexVISList::copy(const VISVolIndexVISList& other)
{
    clean();
    VISVolIndexVISList* list_tmp = (VISVolIndexVISList*)&other;
    Link<VISVolIndex> *current_element;

    current_element = list_tmp->head();
    while (current_element != NULL)
	{
	    this->appendItem(current_element->data());
	    current_element = current_element->next();
	}
}


void VolIndexValueVISList::clean()
{
    reset();
    while(valid())
	removeCurrent();
}


boolean VolIndexValueVISList::removeCurrent()
{
    Link<VolIndexValue> *element_tmp;
    if (valid())
	{
	    element_tmp = _current_element->next();
	    removeItem(_current_element);
	    _current_element = element_tmp;
	    return(TRUE);
	}
    else
	return(FALSE);
}
    
boolean VolIndexValueVISList::stepForward() 
{
    return((valid())&&
	   ((_current_element=_current_element->next())||TRUE));
}

boolean  VolIndexValueVISList::stepBack()
{
    return((valid())&&
	   ((_current_element=_current_element->prev())||TRUE));
}


void VolIndexValueVISList::copy(const VolIndexValueVISList& other)
{
    clean();
    VolIndexValueVISList* list_tmp = (VolIndexValueVISList*)&other;
    Link<VolIndexValue> *current_element;

    current_element = list_tmp->head();
    while (current_element != NULL)
	{
	    this->appendItem(current_element->data());
	    current_element = current_element->next();
	}
}
