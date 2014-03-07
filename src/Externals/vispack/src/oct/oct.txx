
#include "util/defs.h"

#include "util/mathutil.h"

#ifndef MANUAL_INSTANTIATION
#include "oct/oct.h"
#endif

// *********************************************************************
//        VISOct Class Implementation
// *********************************************************************

template <class T>
void VISOct<T>::level(int l)
    {
	_level = l; 
//	_level_factor = 1;
//	for (int i = 1; i < l; i++)
//	    _level_factor *= 2;
    }

template <class T>
boolean VISOct<T>::childChanged()
{
    VISOct<T>* node_ptr;
//    T val;

    boolean consolidate = TRUE;
    if ((node_ptr = _children[0]) != NULL)
	{
	    T v = (node_ptr->_value);
	    for (int i = 1; (i < 8)&&consolidate; i++)
		if ((!_children[i]->_leaf)||
		    ((_children[i]->_value) != v))
		    consolidate = FALSE;

	    if (consolidate)
		{
		    _value = v;
		    for (i = 0; i < 8; i++)
			{
			    delete _children[i];
			    _children[i] = NULL;
			}
		    _leaf = TRUE;
		}
	    return(consolidate);
	}
    else return(FALSE);
}

template <class T>
void VISOct<T>::resolveChange()
{
//    VISOct<T>* node_ptr;
//    T* value_ptr;
    
    if (_parent)
	{
	    if (_parent->childChanged())
		_parent->resolveChange();
	}
}

template <class T>
VISOct<T>::~VISOct()
{
    deleteChildren();
}

template <class T>
void VISOct<T>::deleteChildren()
{
    if (!isLeaf())
	for (int i = 0; i < 8; i++)
	    if (_children[i] != NULL)
		{
		    delete _children[i];
		    _children[i] = NULL;
		}
}

template <class T>
T VISOct<T>::deleteChildrenAverage()
{
    int number;
    T total;

    if (_leaf) return(_value);

    for (int i = 0; (i < 8)&&(_children[i] == NULL); i++);
    total = _children[i++]->deleteChildrenAverage();
    number = 1;
    for (; i < 8; i++)
	if (_children[i] != NULL)
	    {
		number++;
		total += _children[i]->deleteChildrenAverage();
		delete _children[i];
		_children[i] = NULL;
	    }
    return(total/number);
}

// this is a really neat recursive (slow!!!) algorithm
template <class T>
boolean VISOct<T>::prune()
{
    boolean children_pruned = FALSE;
    if (!_leaf)
	{
	    for (int i = 0; i < 8; i++)
		if (_children[i]->prune())
		    children_pruned = TRUE;
	    if (children_pruned)
		return(childChanged());
	}
    else
	return(TRUE);
}

template <class T>
void VISOct<T>::createChildren()
{
    for (int i = 0; i < 8; i++)
	{
	    _children[i] = new VISOct<T>(_value, this, (_level - 1));
	}
    _leaf = FALSE;
}

template <class T>
VISOct<T>::VISOct(const T& v, VISOct<T>* my_parent, int l)
{
    _parent = my_parent;
    _value = v;
    _leaf = TRUE;
    level(l);
}

template <class T>
VISOct<T>::VISOct(int l)
{
    _parent = NULL;
    _leaf = TRUE;
    level(l);
}

template <class T>
VISOct<T>::VISOct()
{
    _parent = NULL;
    _leaf = TRUE;    
    level(0);
}

template <class T>
VISOct<T>::VISOct(const T& v,  int l)
{
    _parent = NULL;
    _value = v;
    _leaf = TRUE;
    level(l);    
}


template <class T>
void VISOct<T>::setValue(const T& v, unsigned x, unsigned y, unsigned z)
{
//    unsigned int x_current, y_current, z_current;
    if ((isLeaf())&&(v != _value))
	{
	    if (level() > 0)
		{
		    createChildren();
		    childRef(x, y, z)->setValue(v, remainder(x),
						remainder(y),
						remainder(z));
		}
	    else
		{
		    value(v);
		    resolveChange();
		}
	}
    else
	childRef(x, y, z)->setValue(v, remainder(x),
				    remainder(y),
				    remainder(z));
}

template <class T>
void VISOct<T>::setValue(const T& v, unsigned x, unsigned y, unsigned z, int l)
{
    unsigned int x_current, y_current, z_current;
    if (isLeaf())
	{
	    if (v != _value)
		{
		    if (level() > l)
			{
			    createChildren();
			    childRef(x, y, z)->setValue(v, remainder(x),
							remainder(y),
							remainder(z), l);
			}
		    else
			{
			    deleteChildren();
			    value(v);
			    resolveChange();
			}
		}
	}
    else
    {
	childRef(x, y, z)->setValue(v, remainder(x),
				      remainder(y),
				      remainder(z), l);
    }
}


template <class T>
void VISOct<T>::printData() const
{
    int i, j, k;
    for (i = 0; i < size(); i++)
	{
	    for (j = 0; j < size(); j++)
		{
		    for (k = 0; k < size(); k++)
			printf("%3.2f ", (float)value(k, j, i));
		    printf("\n");
		}
	    printf("\n");
	}
}

template <class T>
void VISOct<T>::printDataLevels() const
{
    int i, j, k, l, m;
    int blocks;
    const VISOct<T>* node_ptr;
    for (l = level(); l >= 0; l--)
	{
	    printf("Level %d:\n", l);
	    
	    blocks = 1; 
	    for (m = 0; m < l; m++)
		blocks *= 2;

	    for (i = 0; i < size(); i += blocks)
		{
		    for (j = 0; j < size(); j += blocks)
			{
			    for (k = 0; k < size(); k += blocks)
				{
// get the leaf or node of depth l
				    if (!((node_ptr 
					   = leaf(k, j, i, l))->isLeaf()))
					printf("X.XX ");
				    else
					printf("%3.2f ", 
					       (float)(node_ptr->value()));
				}
			    printf("\n");
			}
		    printf("\n");
		}
	    printf("\n");
	}
}


template <class T>
VISOct<T>* VISOct<T>::leafRef(unsigned x, unsigned y, unsigned z, unsigned l)
    {
	if ((_leaf)||(_level <= l))
	    return(this);
	else
// you have to check for levels in here
	    return(childRef(index(x, y, z))->leafRef(
		remainder(x), remainder(y), remainder(z), l
		   ));
    }

template <class T>
VISOct<T>* VISOct<T>::leafRefCreate(unsigned x, unsigned y, 
				      unsigned z, unsigned l)
    {
	if (_level <= l)
	    return(this);
	else
	    {
		if (_leaf)
		    createChildren();
		return(childRef(index(x, y, z))->leafRef(
		    remainder(x), remainder(y), remainder(z), l
		    ));
	    }
    }

template <class T>
const VISOct<T>* VISOct<T>::leaf(unsigned x, unsigned y, unsigned z, 
				   unsigned l) const
    {
	if ((_leaf)||(_level <= l))
	    return(this);
	else
// you have to check for levels in here
	    return(child(index(x, y, z))->leaf(
		remainder(x), remainder(y), remainder(z), l 
		   ));
    }

template <class T>
VISOct<T>* VISOct<T>::leafRef(unsigned x, unsigned y, unsigned z)
{
	if (_leaf)
	    return(this);
	else
// you have to check for levels in here
	    return(childRef(index(x, y, z))->leafRef(
		remainder(x), remainder(y), remainder(z) 
		   ));
    }

template <class T>
const VISOct<T>* VISOct<T>::leaf(unsigned x, unsigned y, unsigned z) const
{
	if (_leaf)
	    return(this);
	else
// you have to check for levels in here
	    return(child(index(x, y, z))->leaf(
		remainder(x), remainder(y), remainder(z) 
		   ));
}

template <class T>
const T& VISOct<T>::value(unsigned x, unsigned y, unsigned z) const
{
    if (isLeaf())
	return(_value);
    else
// you have to check for levels in here
	return(child(index(x, y, z))->value(
	    remainder(x), remainder(y), remainder(z) 
	       ));
}


template <class T>
const T& VISOct<T>::value(unsigned x, unsigned y, unsigned z, int l) const
{
    if ((!_leaf)&&(_level > l))
	return(child(index(x, y, z))->value(
	    remainder(x), remainder(y), remainder(z), l 
	       ));
    else
	return(_value);

}

// this moves everything up (or down for (d < 0)) a level
// and removes anything that falls below the 0 level.
template <class T>
int VISOct<T>::promote(int d)
{
    level(level() + d);

    if (!_leaf)
	{
	    if (level() <= 0)
		{
		    value(deleteChildrenAverage());
		}
	    else
		for (int i = 0; i < 8; i++)
		    {
			_children[i]->promote(d);
		    }
	}
    return(0);
}



// *********************************************************************
//        VISOctPtr Class Implementation
// *********************************************************************


template <class T>
void VISOctPtr<T>::level(int l)
    {
	_level = l; 
//	_level_factor = 1;
//
// level 0 must maintain _level_factor 
// one or else you run into problems with the mod
// and the division
//
//	for (int i = 1; i < l; i++)
//	    _level_factor *= 2;
    }

template <class T>
boolean VISOctPtr<T>::childChanged()
{
//    VISOctPtr<T>* node_ptr;
    T* value_ptr;
    int i;

    for (i = 0; ((i < 8)&&(_children[i]->_value != NULL)); i++);
    if (i < 8)
	{
	    T v = *(_children[i]->_value);
	    boolean consolidate = TRUE;
	    for (i = i + 1; (i < 8)&&consolidate; i++)
		if ((value_ptr =_children[i]->_value) != NULL)
		    {
			if (*value_ptr != v)
			    consolidate = FALSE;
		    }
		else
		    consolidate = FALSE;

	    if (consolidate)
		{
		    _value = new T(v);
		    for (i = 0; i < 8; i++)
			{
			    delete _children[i];
			    _children[i] = NULL;
			}
		}
	    return(consolidate);
	}
    else return(FALSE);
}

template <class T>
void VISOctPtr<T>::resolveChange()
{
//    VISOctPtr<T>* node_ptr;
//    T* value_ptr;
    
    if (_parent)
	{
	    if (_parent->childChanged())
		_parent->resolveChange();
	}
}

template <class T>
VISOctPtr<T>::~VISOctPtr()
{
    if (_value != NULL)
	delete _value;
    else
	deleteChildren();
}

template <class T>
void VISOctPtr<T>::deleteChildren()
{
    if (!isLeaf())
	for (int i = 0; i < 8; i++)
	    if (_children[i] != NULL)
		{
		    delete _children[i];
//		    _children[i] = NULL;
		}
    delete _children;
    _children = NULL;
}

template <class T>
T VISOctPtr<T>::deleteChildrenAverage()
{
    int number;
    T total;

    if (_value) return(*_value);

    for (int i = 0; (i < 8)&&(_children[i] == NULL); i++);
    total = _children[i++]->deleteChildrenAverage();
    number = 1;
    for (; i < 8; i++)
	if (_children[i] != NULL)
	    {
		number++;
		total += _children[i]->deleteChildrenAverage();
		delete _children[i];
		_children[i] = NULL;
	    }
    delete _children;
    _children = NULL;
    return(total/number);
}

// this is a really neat recursive (slow!!!) algorithm
template <class T>
boolean VISOctPtr<T>::prune()
{
    boolean children_pruned = FALSE;
    if (_value == NULL)
	{
	    for (int i = 0; i < 8; i++)
		if (_children[i]->prune())
		    children_pruned = TRUE;
	    if (children_pruned)
		return(childChanged());
	}
    else
	return(TRUE);
}

template <class T>
void VISOctPtr<T>::createChildren()
{
    _children = new VISOctPtr<T>*[8];
    for (int i = 0; i < 8; i++)
	{
	    _children[i] = new VISOctPtr<T>(*_value, this, (_level - 1));
	}
    delete _value;
    _value = NULL;
}

template <class T>
VISOctPtr<T>::VISOctPtr(const T& val_in, VISOctPtr<T>* my_parent, int l)
{
    _parent = my_parent;
    _value = new T(val_in);
    level(l);
    _children = NULL;
}

template <class T>
VISOctPtr<T>::VISOctPtr(int l)
{
    _parent = NULL;
    _value = new T;
    level(l);
    _children = NULL;
}

template <class T>
VISOctPtr<T>::VISOctPtr()
{
    _parent = NULL;
    _value = new T;
    level(0);
    _children = NULL;
}

template <class T>
VISOctPtr<T>::VISOctPtr(const T& val_in,  int l)
{
    _parent = NULL;
    _value = new T(val_in);
    level(l);    
    _children = NULL;
}

template <class T>
VISOctPtr<T>::VISOctPtr(T* val_in,  int l)
{
    _parent = NULL;
    _value = val_in;
    level(l);    
    _children = NULL;
}

template <class T>
void VISOctPtr<T>::setValue(const T& v, unsigned x, unsigned y, unsigned z)
{
//    unsigned int x_current, y_current, z_current;
    if (_value != NULL)
	{
	    if (v != *_value)
		{
		    if (level() > 0)
			{
			    createChildren();
			    childRef(x, y, z)->setValue(v, 
							remainder(x),
							remainder(y),
							remainder(z));
			}
		    else
			{
			    value(v);
			    resolveChange();
			}
		}
	}
    else
	childRef(x, y, z)->setValue(v, remainder(x),
				    remainder(y),
				    remainder(z));
}

template <class T>
void VISOctPtr<T>::setValue(const T& v, unsigned x, unsigned y, 
			     unsigned z, int l)
{
//    unsigned int x_current, y_current, z_current;
    if (_value != NULL)
	{
	if ((v != *_value))
	    {
		if (level() > l)
		    {
			createChildren();
			childRef(x, y, z)->setValue(v, remainder(x),
						    remainder(y),
						    remainder(z), l);
		    }
		else
		    {
			deleteChildren();
			value(v);
			resolveChange();
		    }
	    }
	}
    else
	{
	    childRef(x, y, z)->setValue(v, remainder(x),
					remainder(y),
					remainder(z), l);
	}
}


template <class T>
void VISOctPtr<T>::printData() const
{
    int i, j, k;
    for (i = 0; i < size(); i++)
	{
	    for (j = 0; j < size(); j++)
		{
		    for (k = 0; k < size(); k++)
			printf("%3.2f ", (float)value(k, j, i));
		    printf("\n");
		}
	    printf("\n");
	}
}

template <class T>
void VISOctPtr<T>::printDataLevels() const
{
    int i, j, k, l, m;
    const T *v_ptr;
    int blocks;
    for (l = level(); l >= 0; l--)
	{
	    printf("Level %d:\n", l);
	    
	    blocks = 1; 
	    for (m = 0; m < l; m++)
		blocks *= 2;

	    for (i = 0; i < size(); i += blocks)
		{
		    for (j = 0; j < size(); j += blocks)
			{
			    for (k = 0; k < size(); k += blocks)
				{
				    if ((v_ptr = value(k, j, i, l))
					== NULL)
					printf("X.XX ");
				    else
					printf("%3.2f ", (float)(*v_ptr));
				}
			    printf("\n");
			}
		    printf("\n");
		}
	    printf("\n");
	}
}



template <class T>
const VISOctPtr<T>* VISOctPtr<T>::leaf(unsigned x, unsigned y, unsigned z) 
    const
{
	if (_value != NULL)
	    return(this);
	else
// you have to check for levels in here
	    return(child(index(x, y, z))->leaf(
		remainder(x), remainder(y), remainder(z) 
		   ));
}

template <class T>
const T& VISOctPtr<T>::value(unsigned x, unsigned y, unsigned z) const
{
    if (_value != NULL)
	return(*_value);
    else
// you have to check for levels in here
	return(child(index(x, y, z))->value(
	    remainder(x), remainder(y), remainder(z) 
	       ));
}


template <class T>
const T* VISOctPtr<T>::value(unsigned x, unsigned y, unsigned z, int l) const
{
    if ((_value == NULL)&&(_level > l))
	return(child(index(x, y, z))->value(
	    remainder(x), remainder(y), remainder(z), l 
	       ));
    else
	return(_value);

}

// this moves everything up (or down for (d < 0)) a level
// and removes anything that falls below the 0 level.

template <class T>
int VISOctPtr<T>::promote(int d)
{
    level(level() + d);

    if (!_value)
	{
	    if (level() <= 0)
		{
		    value(deleteChildrenAverage());
		}
	    else
		for (int i = 0; i < 8; i++)
		    {
			_children[i]->promote(d);
		    }
	}
    return(0);
}



template <class T>
VISOctPtr<T>* VISOctPtr<T>::leafRef(unsigned x, unsigned y, 
				      unsigned z, unsigned l) 
{
    if ((_value)||(_level <= l))
	return(this);
    else
// you have to check for levels in here
	return(childRef(index(x, y, z))->leafRef(
	    remainder(x), remainder(y), remainder(z), l
	    ));
}

template <class T>
VISOctPtr<T>* VISOctPtr<T>::leafRef(unsigned x, unsigned y, 
				      unsigned z) 
{
    if (_value)
	return(this);
    else
// you have to check for levels in here
	return(childRef(index(x, y, z))->leafRef(
	    remainder(x), remainder(y), remainder(z)));
}

template <class T>
VISOctPtr<T>* VISOctPtr<T>::leafRefCreate(unsigned x, unsigned y, 
					    unsigned z, unsigned l) 
    {
	if (_level <= l)
	    return(this);
	else
	    {
		if (_value)
		    createChildren();
		return(childRef(index(x, y, z))->leafRefCreate(
		    remainder(x), remainder(y), remainder(z), l
		    ));
	    }
    }

template <class T>
const VISOctPtr<T>* VISOctPtr<T>::leaf(unsigned x, unsigned y, unsigned z, 
				   unsigned l) const
    {
	if ((_leaf)||(_level <= l))
	    return(this);
	else
// you have to check for levels in here
	    return(child(index(x, y, z))->leaf(
		remainder(x), remainder(y), remainder(z), l
		   ));
    }


template <class T>
const VISOct<T>& VISOct<T>::assignVol(const VISVolume<T>& vol, T padding)
{
    unsigned max_size;
    max_size = max(vol.depth(), max(vol.height(), vol.width()));

    int levels = ceil(log((float)max_size)/log(2.0));

    deleteChildren();
    _leaf = TRUE;
    level(levels);

    setValue(padding, 0, 0, 0, levels);
    
    int i, j, k;

    for (k = 0; k < vol.depth(); k++)
	for (j = 0; j < vol.depth(); j++)
	    for (i = 0; i < vol.depth(); i++)
		value(i, j, k) = vol.itemAt(i, j, k);

    return(*this);

}


template <class T>
const VISOctPtr<T>& VISOctPtr<T>::assignVol(const VISVolume<T>& vol, 
					      T padding)
{
    unsigned max_size;
    max_size = max(vol.depth(), max(vol.height(), vol.width()));

    int levels = (int)ceil(log((float)max_size)/log(2.0));

    deleteChildren();
    _value = new T;
    level(levels);

    setValue(padding, 0, 0, 0, levels);
    
    int i, j, k;

    for (k = 0; k < vol.depth(); k++)
	for (j = 0; j < vol.depth(); j++)
	    for (i = 0; i < vol.depth(); i++)
		setValue(vol.itemAt(i, j, k), i, j, k);
    
    return(*this);

}

template <class T>
VISVolume<T> VISOctPtr<T>::volume() const
{

    int w, h, d;
    int i, j, k;
    VISVolume<T> r(w = width(), h = height(), d = depth());
    
    for (k = 0; k < d; k++)
	for (j = 0; j < h; j++)
	    for (i = 0; i < w; i++)
		{
		    r.at(i, j, k) = value(i, j, k);
		}
    return(r);
}

