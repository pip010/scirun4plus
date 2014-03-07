#include "image/indexlist.h"

void VISImIndexVISList::clean()
{
    reset();
    while(valid())
	removeCurrent();
}


jpegBoolean VISImIndexVISList::removeCurrent()
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
    
jpegBoolean VISImIndexVISList::stepForward() 
{
    return((valid())&&
	   ((_current_element=_current_element->next())||TRUE));
}

jpegBoolean  VISImIndexVISList::stepBack()
{
    return((valid())&&
	   ((_current_element=_current_element->prev())||TRUE));
}



