
#include "oct/octmesh.h"

//
// give back a neighbor (if links exist) for a relative distance x, y, z
//

template <class T>
OctMeshLeaf<T>* OctMeshLeaf<T>::neighbor(int x, int y, int z)
{
    unsigned i;
    if ((x == 0)&&(y == 0)&&(z == 0))
	return(this);
    else
 	if (_neighbors[i = nIndex(x, y, z)])
	    return(_neighbors[i]->neighbor
		   (
		       nRemainderX(x, y, z),
		   nRemainderY(x, y, z),
		   nRemainderZ(x, y, z)
		   ));
	else
	    return(NULL);
}

template <class T>
OctMeshLeaf<T>* OctMeshLeaf<T>::neighbor(int i)
{
    return(_neighbors[i]);
}

// make it in line in the .h
//template <class T>
//OctMeshNode<T>* OctMeshNode<T>::childRef(unsigned i)
//{
//    return((OctMeshNode<T>*)_children[i]);
//}


//template <class T>
//OctMeshNode<T>* OctMeshNode<T>::leafRef(unsigned x, 
//						unsigned y, 
//						unsigned z)
//{
//    return((OctMeshNode<T>*)(VISOctPtr<T>::leafRef(x, y, z)));
//}
//
//template <class T>
//OctMeshNode<T>* OctMeshNode<T>::leafRefCreate(unsigned x, 
//			      unsigned y, 
//			      unsigned z, 
//			      int l)
//{
//    return((OctMeshNode<T>*)(VISOctPtr<T>::leafRefCreate(x, y, z, l)));
//}


template <class T>
void OctMeshNode<T>::createChildren()
{
    _children = new VISOctPtr<T>*[8];
    if (_level > 1)
	{
	    for (int i = 0; i < 8; i++)
		{
		  _children[i] 
		      = new OctMeshNode<T>(*_value, this, (_level - 1), 
					       i);
		}
	    delete _value;
	    _value = NULL;
	}
    else
	{
	    for (int i = 0; i < 8; i++)
		{
		    _children[i] = new OctMeshLeaf<T>
			(*_value, this, (_level - 1), i);
		}
	    delete _value;
	    _value = NULL;

// now connect it to it's neighbors in each direction.
// this is important because it involves (possibly) the whole tree
	    for (i = 0; i < 8; i++)
		{
		    for (int j = 0; j < 6; j++)
			((OctMeshLeaf<T>*)childRef(i))
			    ->connect(getNeighborLeaf(0, 0, 0, i, (Direct)j), 
				      (Direct)j);
		}
	}
}

template <class T>
void OctMeshLeaf<T>::createChildren()
{
    printf("OctMeshLeaf<T>::createChildren() - cannot create children for \
for OctMeshLeaf");
}

template <class T>
OctMeshNode<T>::OctMeshNode(const T& val_in, 
				    OctMeshNode<T>* my_parent, 
				    int l, int childnum)
    :VISOctPtr<T>(val_in,  (VISOctPtr<T>*)my_parent, l)
{
    _child_num = childnum;
}

template <class T>
OctMeshNode<T>::OctMeshNode(T* val_in, OctMeshNode<T>* my_parent, 
				    int l, int childnum)
    :VISOctPtr<T>(val_in,  (VISOctPtr<T>*)my_parent, l)
{
    _child_num = childnum;
}
;

template <class T>
OctMeshNode<T>::OctMeshNode(const T& val_in, 
				    OctMeshNode<T>* my_parent, 
				int l)
    :VISOctPtr<T>(val_in, my_parent, l)
{
    _child_num = -1;
}

template <class T>
OctMeshNode<T>::OctMeshNode():VISOctPtr<T>()
{
    _child_num = -1;
}

template <class T>
OctMeshLeaf<T>* OctMeshNode<T>::getNeighborCreate(int x, int y, 
						    int z, unsigned child_num, 
						    Direct direct)
{
//    int neighbor_child_num;
    OctMeshNode<T>* the_leaf;
    unsigned level_factor;

// the ChildRoute table returns a positive or negative number depending on 
// whether or not the neighbor the child is looking for is in your subtree.
// 
// If it's not in your subtree, then adjust the index to reflect your level
// and pass the request to the next level up
// If it is in your subtree, then return the leaf associated with that subtree

    // you need the new index when given a child number and a direction
//    if ((neighbor_child_num = 
    if (
	(ChildRoute[direct][child_num]) < 0)
	{
	    if (parent() != NULL)
		return(parentRef()->getNeighborCreate
		       (ChildPos[X_AXIS][child_num]*
			(level_factor 
			 = __LEVEL_FACTOR[_level]) + x, 
			ChildPos[Y_AXIS][child_num]*level_factor + y, 
			ChildPos[Z_AXIS][child_num]*level_factor + z,
			childNum(), direct));
	    else 
		return(NULL);
	}
    else 
	{
// when you go to look for a the neigbor leaf you have to add the 
// "DirectAdjust" term, which accounts for the fact you are looking for 
// a neighbor not the original child.
	    the_leaf = leafRefCreate
		(ChildPos[X_AXIS][child_num]*(level_factor 
					      = __LEVEL_FACTOR[_level]) + x
		 + DirectAdjust[X_AXIS][direct], 
		 ChildPos[Y_AXIS][child_num]*level_factor + y
		 + DirectAdjust[Y_AXIS][direct], 
		 ChildPos[Z_AXIS][child_num]*level_factor + z
		 + DirectAdjust[Z_AXIS][direct], 0);

// the assumption is that all level 0 leaves are meshleaf's and 
// all other levels are not
//	    if (the_leaf->level() == 0)
		return((OctMeshLeaf<T>*)the_leaf);
//	    else
// this means that the attempt to find a neighbor failed, i.e. there is none.
//		return(NULL);
	}
}



template <class T>
OctMeshLeaf<T>* OctMeshNode<T>::getNeighborLeaf(int x, int y, 
						    int z, unsigned child_num, 
						    Direct direct)
{
//    int neighbor_child_num;
    OctMeshNode<T>* the_leaf;
    unsigned level_factor;

// the ChildRoute table returns a positive or negative number depending on 
// whether or not the neighbor the child is looking for is in your subtree.
// 
// If it's not in your subtree, then adjust the index to reflect your level
// and pass the request to the next level up
// If it is in your subtree, then return the leaf associated with that subtree

    // you need the new index when given a child number and a direction
//    if ((neighbor_child_num = 
    if (
	(ChildRoute[direct][child_num]) < 0)
	{
	    if (parent() != NULL)
		return(parentRef()->getNeighborLeaf
		       (ChildPos[X_AXIS][child_num]*
			(level_factor = __LEVEL_FACTOR[_level])	+ x, 
			ChildPos[Y_AXIS][child_num]*level_factor + y, 
			ChildPos[Z_AXIS][child_num]*level_factor + z,
			childNum(), direct));
	    else 
		return(NULL);
	}
    else 
	{
// when you go to look for a the neigbor leaf you have to add the 
// "DirectAdjust" term, which accounts for the fact you are looking for 
// a neighbor not the original child.
	    the_leaf = leafRef
		(ChildPos[X_AXIS][child_num]*
		 (level_factor = __LEVEL_FACTOR[_level])+ x
		 + DirectAdjust[X_AXIS][direct], 
		 ChildPos[Y_AXIS][child_num]*level_factor + y
		 + DirectAdjust[Y_AXIS][direct], 
		 ChildPos[Z_AXIS][child_num]*level_factor + z
		 + DirectAdjust[Z_AXIS][direct]);

// the assumption is that all level 0 leaves are meshleaf's and 
// all other levels are not
	    if (the_leaf->level() == 0)
		return((OctMeshLeaf<T>*)the_leaf);
	    else
// this means that the attempt to find a neighbor failed, i.e. there is none.
		return(NULL);
	}
}

// if it can't find a leaf, this one will return a node
template <class T>
OctMeshNode<T>* OctMeshNode<T>::getNeighborNode(int x, int y, 
						    int z, unsigned child_num, 
						    Direct direct)
{
//    int neighbor_child_num;
    OctMeshNode<T>* the_leaf;
    unsigned level_factor;

// the ChildRoute table returns a positive or negative number depending on 
// whether or not the neighbor the child is looking for is in your subtree.
// 
// If it's not in your subtree, then adjust the index to reflect your level
// and pass the request to the next level up
// If it is in your subtree, then return the leaf associated with that subtree

    // you need the new index when given a child number and a direction
//    if ((neighbor_child_num = 
    if ((ChildRoute[direct][child_num]) < 0)
	{
	    if (parent() != NULL)
		return(parentRef()->getNeighborNode
		       (ChildPos[X_AXIS][child_num]*
			(level_factor = __LEVEL_FACTOR[_level]) + x,
			ChildPos[Y_AXIS][child_num]*level_factor + y, 
			ChildPos[Z_AXIS][child_num]*level_factor + z,
			childNum(), direct));
	    else 
// this means you hav probably moved off of the grid
		return(NULL);
	}
    else 
	{
// when you go to look for a the neigbor leaf you have to add the 
// "DirectAdjust" term, which accounts for the fact you are looking for 
// a neighbor not the original child.
	    the_leaf = leafRef(ChildPos[X_AXIS][child_num]*
			       (level_factor = __LEVEL_FACTOR[_level]) + x
			       + DirectAdjust[X_AXIS][direct], 
			       ChildPos[Y_AXIS][child_num]*level_factor + y
			       + DirectAdjust[Y_AXIS][direct], 
			       ChildPos[Z_AXIS][child_num]*level_factor + z
			       + DirectAdjust[Z_AXIS][direct]);
// the assumption is that all level 0 leaves are meshleaf's and 
// all other levels are not
	    return(the_leaf);
	}
}

// removes the links between you and your neighbors
template <class T>
void OctMeshLeaf<T>::deleteNlinks()
{    
    OctMeshLeaf<T>* other_leaf;
    for (int j = 0; j < 6; j++)
	{
	    if ((other_leaf = _neighbors[j]))
		{
		    _neighbors[j] = NULL;
		    other_leaf->_neighbors[OppositeDirection[j]] = NULL;
		}
	}
}

// this moves everything up (or down for (d < 0)) a level
// and removes anything that falls below the 0 level.
//
// In the OctMesh case you have to deal specially with nodes that come
// and go from the zero leaf level
// 
template <class T>
int OctMeshNode<T>::promote(int d)
{
    int stat;
    level(level() + d);
    if (d > 0)
	{
	    if (!_value)
		{
		    for (int i = 0; i < 8; i++)
			{
			    stat = childRef(i)->promote(d);
			    if (stat < 0)
				downgradeChild(i);
			    else if (stat > 0)
// this should only be returned by a leaf node
				upgradeChild(i);
			}
		}
	    return(0);
	}
    else
	{
	    if (level() <= 0)
		{
		    value(deleteChildrenAverage());
		    return(-1);
		}
// what should it return in this case
	    return(0);
	}
}

template <class T>
int OctMeshLeaf<T>::promote(int d)
{
    int stat;
    level(level() + d);
    if (d > 0)
	{
	    if (!_value)
		{
		    for (int i = 0; i < 8; i++)
			{
			    stat = childRef(i)->promote(d);
			    if (stat < 0)
				downgradeChild(i);
			    else if (stat > 0)
// this should only be returned by a leaf node
				upgradeChild(i);
			}
		}
	    return(1);
	}
    else
	{
	    if (level() <= 0)
		{
		    deleteChildren();
		    return(0);
		}
	    else return(-1);
	}
}

//
// these next two upgrade and downgrade children from leaves to nodes 
// used in the promotion process to keep the hiearchy clean and neat, 
// with leaves only at the zero level in the hiearchy
//
template <class T>
void OctMeshNode<T>::upgradeChild(unsigned child_num)
{
 
// care must be taken to see that this is only called after 
// a check has been made to insure that this child is a leaf

    OctMeshLeaf<T>* child_tmp 
	= (OctMeshLeaf<T>*)childRef(child_num);
    _children[child_num] = child_tmp->becomeNode();
    delete child_tmp;
}

template <class T>
void OctMeshNode<T>::downgradeChild(unsigned child_num)
{
// if the child is not level 1 then you really can't do this
    if ((child(child_num)->level()) == 0)
	{
	    OctMeshNode<T>* child_tmp 
		= childRef(child_num);
	    _children[child_num] = new OctMeshLeaf<T>(child_tmp);
	    delete child_tmp;
	    for (int j = 0; j < 6; j++)
		((OctMeshLeaf<T>*)childRef(child_num))
//
// fix this
//
		    ->connect(getNeighborLeaf(0, 0, 0, 
					      child_num, (Direct)j), 
			      (Direct)j);
	}
    else
	printf("Tried to downgrade none-leaf node");
}


template <class T>
OctMeshLeaf<T>::~OctMeshLeaf()
{
    if (_value != NULL)
	{
	    delete _value;
	    deleteNlinks();
	}
    else
	deleteChildren();
}

//
// take out all of the checks for parents, because the case that indicates 
// the position takes care of that (by checking the routing table to see
// if it's necessary to go up the tree)
//

// the connect will work it's way up the tree and then back down until 
// it finds the right leaf/node to connect to
template< class T >
OctMeshNode<T>* OctMeshLeaf<T>::getNeighborNode(Direct direct)
{
    if (_neighbors[direct] != NULL)
	return(_neighbors[direct]);
    else
	if (_parent)
	    return(((OctMeshNode<T>*)_parent)->
		   getNeighborNode(0, 0, 0, _child_num, direct));
//		   getNeighborNode(0, 0, 0, 0, direct));
//		   getNeighborNode(direct));
	else
	    return(NULL);
}

// the connect will work it's way up the tree and then back down until 
// it finds the right leaf/node to connect to
template< class T >
OctMeshLeaf<T>* OctMeshLeaf<T>::getNeighborCreate(Direct direct)
{
    if (_neighbors[direct] != NULL)
	return(_neighbors[direct]);
    else
	if (_parent)
	    return(((OctMeshNode<T>*)_parent)->
		   getNeighborCreate(0, 0, 0, _child_num, direct));
//		   getNeighborNode(0, 0, 0, _child_num, direct));
//		   getNeighborCreate(0, 0, 0, 0, direct));
	else
	    return(NULL);
}

// the connect will work it's way up the tree and then back down until 
// it finds the right leaf/node to connect to
template< class T >
OctMeshLeaf<T>* OctMeshLeaf<T>::getNeighborLeaf(Direct direct)
{
    if (_neighbors[direct] != NULL)
	return(_neighbors[direct]);
    else
	if (_parent)
	    return(((OctMeshNode<T>*)_parent)->
		   getNeighborLeaf(0, 0, 0, _child_num, direct));
//		   getNeighborCreate(0, 0, 0, 0, direct));
	else
	    return(NULL);
}

// the connect will work it's way up the tree and then back down until 
// it finds the right leaf/node to connect to
template< class T >
OctMeshNode<T>* OctMeshNode<T>::getNeighborNode(Direct direct)
{
    if (_parent)
	return(((OctMeshNode<T>*)_parent)
	       ->getNeighborNode(0, 0, 0, _child_num, direct));
    else
	return(NULL);
}


// the connect will work it's way up the tree and then back down until 
// it finds the right leaf/node to connect to
template< class T >
OctMeshLeaf<T>* OctMeshNode<T>::getNeighborCreate(Direct direct)
{
    if (_parent)
	return(((OctMeshNode<T>*)_parent)
	       ->getNeighborCreate(0, 0, 0, _child_num, direct));
    else
	return(NULL);
}

// the connect will work it's way up the tree and then back down until 
// it finds the right leaf/node to connect to
template< class T >
OctMeshLeaf<T>* OctMeshNode<T>::getNeighborLeaf(Direct direct)
{
    if (_parent)
	return(((OctMeshNode<T>*)_parent)
	       ->getNeighborCreate(0, 0, 0, _child_num, direct));
    else
	return(NULL);
}

// this should probably be in line
template <class T>
void OctMeshLeaf<T>::connect(OctMeshLeaf<T>* other_leaf, Direct direct)
{
    if (other_leaf)
	{
	    _neighbors[direct] = other_leaf;
	    other_leaf->_neighbors[OppositeDirection[direct]] = this;
	}
}



// *********************************************************************
// *********************************************************************
// *************************** OCTMESHVOL ******************************
// *********************************************************************
// *********************************************************************


//#ifdef work_in_progress

template <class T>
void OctMeshVol<T>::createChildren()
{
    int size = _w*_h*_d;
    _children = new VISOctPtr<T>*[size];
    if (_level > 1)
	{
	    for (int i = 0; i < size; i++)
		{
		  _children[i] 
		      = new OctMeshNode<T>(*_value, this, (_level - 1), 
					       i);
		}
	    delete _value;
	    _value = NULL;
	}
    else
	{
	    for (int i = 0; i < size; i++)
		{
		    _children[i] = new OctMeshLeaf<T>
			(*_value, this, (_level - 1), i);
		}
	    delete _value;
	    _value = NULL;

// now connect it to it's neighbors in each direction.
// this is important because it involves (possibly) the whole tree
	    for (i = 0; i < size; i++)
		{
		    for (int j = 0; j < 6; j++)
			((OctMeshLeaf<T>*)childRef(i))
			    ->connect(getNeighborLeaf(0, 0, 0, i, (Direct)j), 
				      (Direct)j);
		}
	}
}

template <class T>
void OctMeshVol<T>::deleteChildren()
{
    int size = _w*_h*_d;
    if (_children != NULL)
	{
	    for (int i = 0; i < size; i++)
		if (_children[i] != NULL)
		    {
			delete _children[i];
//    _children[i] = NULL;
		    }
	    delete _children;
	    _children = NULL;
	}
}


template <class T>
OctMeshLeaf<T>* OctMeshVol<T>::getNeighborCreate(int x, int y, 
						 int z, unsigned child_num, 
						 Direct direct)
{
//    int neighbor_child_num;
    OctMeshNode<T>* the_leaf;
    unsigned level_factor;

//    We assume you are as high up as you can go therefore we can just ask 
// for the node directly

    the_leaf = 
	leafRefCreate(childPosX(child_num)
		      *(level_factor 
			= __LEVEL_FACTOR[_level]) + x
		      + DirectAdjust[X_AXIS][direct], 
		      childPosY(child_num)*level_factor + y
		      + DirectAdjust[Y_AXIS][direct], 
		      childPosZ(child_num)*level_factor + z
		      + DirectAdjust[Z_AXIS][direct], 0);
    if (the_leaf == NULL)
	return((OctMeshLeaf<T>*)(&_outside));
    else return((OctMeshLeaf<T>*)the_leaf);


}


template <class T>
OctMeshNode<T>* OctMeshVol<T>::leafRefCreate(unsigned x, unsigned y, 
					     unsigned z, int l)
{
    int level_factor;
    if (checkBounds(x, y, z))
	{
	    _children[index(x/(level_factor = __LEVEL_FACTOR[_level]), 
			    y/level_factor, 
			    z/level_factor)]->leafRefCreate(
				remainder(x), remainder(y), remainder(z), l
				);
	}
    else return(&_outside);
}

template <class T>
OctMeshNode<T>* OctMeshVol<T>::leafRef(unsigned x, unsigned y, unsigned z)
{
    int level_factor;
    if (checkBounds(x, y, z))
	{
	    _children[index(x/(level_factor = __LEVEL_FACTOR[_level]), 
			    y/level_factor, 
			    z/level_factor)]->leafRef(
				remainder(x), remainder(y), remainder(z));
	}
    else return(&_outside);
}

template <class T>
const OctMeshNode<T>* OctMeshVol<T>::leaf(unsigned x, unsigned y, unsigned z) 
const
{
    int level_factor
    if checkBounds(x, y, z)
	{
	    _children[index(x/(level_factor = __LEVEL_FACTOR[_level]), 
			    y/level_factor, 
			    z/level_factor)]->leaf(
				remainder(x), remainder(y), remainder(z));
	}
    else return(&_outside);
}


template <class T>
OctMeshLeaf<T>* OctMeshVol<T>::getNeighborLeaf(int x, int y, 
					       int z, unsigned child_num, 
					       Direct direct)
{
//    int neighbor_child_num;
    OctMeshNode<T>* the_leaf;
    unsigned level_factor;

    the_leaf = leafRef(childPosX(child_num)
		       *(level_factor 
			 = __LEVEL_FACTOR[_level]) + x
		       + DirectAdjust[X_AXIS][direct], 
		       childPosY(child_num)*level_factor + y
		       + DirectAdjust[Y_AXIS][direct], 
		       childPosZ(child_num)*level_factor + z
		       + DirectAdjust[Z_AXIS][direct]);

    if (the_leaf == NULL)
	return((OctMeshLeaf<T>*)(&_outside));
    else if (the_leaf->level() == 0)
	return((OctMeshLeaf<T>*)the_leaf);
    else
// this means that the attempt to find a neighbor failed, i.e. there is none.
	return(NULL);
}



// if it can't find a leaf, this one will return a node
template <class T>
OctMeshNode<T>* OctMeshVol<T>::getNeighborNode(int x, int y, 
						int z, unsigned child_num, 
						Direct direct)
{
//    int neighbor_child_num;
    OctMeshNode<T>* the_leaf;
    unsigned level_factor;
    the_leaf = leafRef(childPosX(child_num)
		       *(level_factor 
			 = __LEVEL_FACTOR[_level]) + x
		       + DirectAdjust[X_AXIS][direct], 
		       childPosY(child_num)*level_factor + y
		       + DirectAdjust[Y_AXIS][direct], 
		       childPosZ(child_num)*level_factor + z
		       + DirectAdjust[Z_AXIS][direct]);
// the assumption is that all level 0 leaves are meshleaf's and 
// all other levels are not
    if (the_leaf == NULL)
	return((OctMeshLeaf<T>*)(&_outside));
    else return(the_leaf);
}


template <class T>
OctMeshVol<T>::OctMeshVol(const T& val_in, int l, 
			  unsigned x, unsigned y, unsigned z)
    :OctMeshNode<T>(val_in, l)
{
    _w = x; _h = y; _d = z;
    deleteChildren();
    createChildren();
}

template <class T>
OctMeshVol<T>::OctMeshVol(const T& val_in, int l, 
			  unsigned x, unsigned y, unsigned z, 
			  T outside_value)
    :OctMeshNode<T>(val_in, l)
{
    _w = x; _h = y; _d = z;
    OctMeshNode<T>::deleteChildren();
    createChildren();
    setOutside(outside_value);
}


template <class T>
void OctMeshVol<T>::setValue(const T& v, unsigned x, unsigned y, unsigned z)
{
    childRef(index(index(x), 
		   index(y), 
		   index(z)))
	->setValue(v, remainder(x), remainder(y), remainder(z));
}

template <class T>
void OctMeshVol<T>::setValue(const T& v, unsigned x, unsigned y, 
		       unsigned z, int l)
{
    childRef(index(index(x), index(y), index(z)))
		      ->setValue(v, remainder(x), 
				 remainder(y), 
				 remainder(z), l);
}

template <class T>
const T* OctMeshVol<T>::value(unsigned x, unsigned y, 
			   unsigned z, int l) const
{
    return(
	childRef(index(index(x), index(y), index(z)))
		       ->value(remainder(x), 
				  remainder(y), 
				  remainder(z), l)
	);
}

template <class T>
const T& OctMeshVol<T>::value(unsigned x, unsigned y, 
			   unsigned z) const
{
    if (checkBounds(x, y, z))
	return(child(index(index(x), index(y), index(z)))
	       ->value(remainder(x), 
		       remainder(y), 
		       remainder(z)));
    else
	return(*_outside.value());
}

template <class T>
int OctMeshVol<T>::promote(int d)
{
    int i, 
// size of the array of children
	size, 
// this tells whether or note you need to upgrade or downgrade a child
// from leaf to node or vise versa
	stat;

    size = _w*_h*_d;

    if ((_level + d) < 1)
	ERROR("OctMeshVol::promote() cannot promite top to less than 1");

    _level += d;

    for (i = 0; i < size; i++)
	{
// must check to see of children need upgrading	    
	    stat = _children[i]->promote(d);
			    if (stat < 0)
				downgradeChild(i);
			    else if (stat > 0)
// this should only be returned by a leaf node
				upgradeChild(i);
	}
}

//#endif
