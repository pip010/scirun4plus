// ******************************************************************
// VISPack. Copyright (c) 1994-2000 Ross Whitaker rtw@utk.edu       *
// For conditions of distribution and use, see the file LICENSE.txt *
// accompanying this distribution.                                  *
// ******************************************************************

// $Id: indexlist.cxx,v 1.1.1.1 2003/02/12 16:51:52 whitaker Exp $

#include "image/indexlist.h"

void VISImIndexList::clean()
{
    reset();
    while(valid())
	removeCurrent();
}


boolean VISImIndexList::removeCurrent()
{
    Link<VISImIndex> *element_tmp;
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
    
boolean VISImIndexList::stepForward() 
{
    return((valid())&&
	   ((_current_element=_current_element->next())||TRUE));
}

boolean  VISImIndexList::stepBack()
{
    return((valid())&&
	   ((_current_element=_current_element->prev())||TRUE));
}



