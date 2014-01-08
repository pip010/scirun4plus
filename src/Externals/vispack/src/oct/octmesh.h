// this is the foundation for an nth (quad or oct) tree object
// it has a the ability to keep track of location withing a tree 
// and return an object at a specific location

// Ross Whitaker 5-3-96
//

#ifndef iris_octmesh_h
#define iris_octmesh_h

#include <stdio.h>
#include "util/defs.h"
#include "oct/oct.h"

template <class T>
class OctMeshLeaf;

typedef enum {X_POS, X_NEG, Y_POS, Y_NEG, 
	      Z_POS, Z_NEG, NONE, X_DIR, 
	      Y_DIR, Z_DIR} Direct;

template <class T>
class OctMeshNode: public VISOctPtr<T>
{
    
    friend OctMeshLeaf<T>;

  protected:
    int _child_num;

    virtual void createChildren();
    void upgradeChild(unsigned child_num);
    void downgradeChild(unsigned child_num);

//
// the connect will work it's way up the tree and then back down until 
// it finds the right leaf to connect to
//
// this is mostly meant to be a request that a child would make of a parent
//
    virtual OctMeshLeaf<T>* getNeighborCreate(int x, int y, int z, 
				    unsigned child_num, 
				    Direct direct);


    virtual OctMeshLeaf<T>* getNeighborLeaf(int x, int y, int z, 
				    unsigned child_num, 
				    Direct direct);

//
// the connect will work it's way up the tree and then back down until 
// it finds the right leaf/node to connect to
//
// this is mostly meant to be a request that a child would make of a parent
//
    virtual OctMeshNode<T>* getNeighborNode(int x, int y, int z, 
					unsigned child_num, 
					Direct direct);


  public:

    OctMeshNode(const T& val_in, OctMeshNode<T>* my_parent, 
		    int l, int childnum);
    OctMeshNode(T* val_in, OctMeshNode<T>* my_parent, 
		    int l, int childnum);
    OctMeshNode(const T& val_in, OctMeshNode<T>* my_parent, 
		    int l);
    OctMeshNode(const T& val_in, int l):VISOctPtr<T>(val_in, l)
	{
	    _child_num = -1;
	}
    OctMeshNode(T* val_in, int l):VISOctPtr<T>(val_in, l)
	{
	    _child_num = -1;
	}
    OctMeshNode();

    boolean isLeaf() const {return(_level == 0);}

    virtual OctMeshNode<T>* getNeighborNode(Direct direct);
    virtual OctMeshLeaf<T>* getNeighborCreate(Direct direct);
    virtual OctMeshLeaf<T>* getNeighborLeaf(Direct direct);

//    OctMeshNode(const T& val_in, OctMeshNode<T>* my_parent, int l);
//    OctMeshNode(T* val_in, const OctMeshNode<T>* my_parent, int l);

    unsigned childNum() {return(_child_num);}
    const OctMeshNode<T>* parent() {return((OctMeshNode<T>*)_parent);}
    OctMeshNode<T>* parentRef() {return((OctMeshNode<T>*)_parent);}

    const OctMeshNode<T>* child(unsigned i) const
    {
	return((OctMeshNode<T>*)_children[i]);
    }

    OctMeshNode<T>* childRef(unsigned i)
    {
	return((OctMeshNode<T>*)_children[i]);
    }

    virtual int promote(int d);



    const OctMeshNode<T>* leaf(unsigned x, 
				   unsigned y, 
				   unsigned z) const
    { return((OctMeshNode<T>*)VISOctPtr<T>::leaf(x, y, z));}

    OctMeshNode<T>* leafRef(unsigned x, 
			    unsigned y, 
			    unsigned z)
    {
	return((OctMeshNode<T>*)(VISOctPtr<T>::leafRef(x, y, z)));
    }

    OctMeshNode<T>* leafRefCreate(unsigned x, 
				  unsigned y, 
				  unsigned z, 
				  int l)
    {
	return((OctMeshNode<T>*)(VISOctPtr<T>::leafRefCreate(x, y, z, l)));
    }


};

// ******************************************************************
// ******************************************************************
// *************************** OctMeshVol ***************************
// ******************************************************************
// ******************************************************************

// #ifdef work_in_progress

template <class T>
class OctMeshVol: public OctMeshNode<T>
{
  protected:
//    VISOctPtr **_children_array;
    virtual void createChildren();
    unsigned _w, _h, _d;
//    void upgradeChild(unsigned child_num);
//    void downgradeChild(unsigned child_num);
    OctMeshLeaf<T> _outside;
    virtual void deleteChildren();

//
// the connect will work it's way up the tree and then back down until 
// it finds the right leaf to connect to
//
// this is mostly meant to be a request that a child would make of a parent
//
    virtual OctMeshLeaf<T>* getNeighborCreate(int x, int y, int z, 
					      unsigned child_num, 
					      Direct direct);
    virtual OctMeshLeaf<T>* getNeighborLeaf(int x, int y, int z, 
					    unsigned child_num, 
					    Direct direct);
//
// the connect will work it's way up the tree and then back down until 
// it finds the right leaf/node to connect to
//
// this is mostly meant to be a request that a child would make of a parent
//
    virtual OctMeshNode<T>* getNeighborNode(int x, int y, int z, 
					    unsigned child_num, 
					    Direct direct);

    unsigned index(unsigned x, unsigned y, unsigned z) const
    {
	return(x + _w*(y + _h*z));
    }

    unsigned index(unsigned i) const
    {
	return(i/__LEVEL_FACTOR[_level]);
    }


  public:

    OctMeshVol(const T& val_in, int l, unsigned x, unsigned y, unsigned z);
    OctMeshVol(const T& val_in, int l, unsigned x, unsigned y, unsigned z, 
	       T outside_value);

    OctMeshVol():OctMeshNode<T>()
	{
	    _w = _h = _d = 0;
	    VISOctPtr<T>::deleteChildren();
//	    createChildren();
	}

    boolean isLeaf() const {return(0);}

    void setOutside(T value)
	{
	    _outside.value(value);
	}

    unsigned size() {return(_w*_h*_d);}

    unsigned childNum() {return(-1);}
    const OctMeshNode<T>* parent() {return(NULL);}
    OctMeshNode<T>* parentRef() {return(NULL);}

    unsigned childPosX(unsigned i)
    {
	return((i%(_w*_h))%_w);
    }

    unsigned childPosY(unsigned i)
    {
	return((i%(_w*_h))/_w);
    }

    unsigned childPosZ(unsigned i)
    {
	return(i/(_w*_h));
    }


    virtual int promote(int d);

    boolean checkBounds(unsigned x, unsigned y, unsigned z) const
    {
	return((x < _w*__LEVEL_FACTOR[_level])&&(y < _h*__LEVEL_FACTOR[_level])
	       &&(z < _d*__LEVEL_FACTOR[_level]));
    }

    const OctMeshNode<T>* leaf(unsigned x, 
			       unsigned y, 
			       unsigned z) const;

    OctMeshNode<T>* leafRef(unsigned x, 
			    unsigned y, 
			    unsigned z);

    OctMeshNode<T>* leafRefCreate(unsigned x, 
				  unsigned y, 
				  unsigned z, 
				  int l);

    void setValue(const T& v, unsigned x, unsigned y, unsigned z);
    void setValue(const T& v, unsigned x, unsigned y, 
		       unsigned z, int l);

    const T* value(unsigned x, unsigned y, 
		   unsigned z, int l) const;
    const T& value(unsigned x, unsigned y, 
		   unsigned z) const;

    unsigned width() const 
    {
	if (_level != 0)
	    return(_w*__LEVEL_FACTOR[_level]);
	else
	    return(1);
    }

    unsigned height() const 
    {
	if (_level != 0)
	    return(_h*__LEVEL_FACTOR[_level]);
	else
	    return(1);
    }

    unsigned depth() const 
    {
	if (_level != 0)
	    return(_d*__LEVEL_FACTOR[_level]);
	else
	    return(1);
    }


};

//#endif // work in progress


const Direct OppositeDirection[6] 
= {X_NEG, 
   X_POS, 
   Y_NEG, 
   Y_POS, 
   Z_NEG, 
   Z_POS};

// give it the direction and the child_num, and it tells you where to direct
// the routing of this connection request.  -1 means route it above;
    const int ChildRoute[6][8] = 
    {
// direct X
	{1, -1, 3, -1, 5, -1, 7, -1},
// direct -X
	{-1, 0, -1, 2, -1, 4, -1, 6},
// direct Y
	{2, 3, -1, -1, 6, 7, -1, -1},
// direct -Y
	{-1, -1, 0, 1, -1, -1, 4, 5},
// direct Z
	{4, 5, 6, 7, -1, -1, -1, -1},
// direct -Y
	{-1, -1, -1, -1, 0, 1, 2, 3}
    };

#define X_AXIS 0
#define Y_AXIS 1
#define Z_AXIS 2

    const int ChildPos[3][8] = 
    {
//  X_AXIS
	{0, 1, 0, 1, 0, 1, 0, 1},
//  Y_AXIS
	{0, 0, 1, 1, 0, 0, 1, 1},
//  Z_AXIS
	{0, 0, 0, 0, 1, 1, 1, 1},
    };

    const int DirectAdjust[3][6] =
    {
//  X_AXIS
	{1, -1, 0, 0, 0, 0}, 
//  Y_AXIS
	{0, 0, 1, -1, 0, 0}, 
//  Z_AXIS
	{0, 0, 0, 0, 1, -1}, 
    };



template <class T>
class OctMeshLeaf: public OctMeshNode<T>
{

  protected:

    virtual void createChildren();

    OctMeshLeaf* _neighbors[6];
    int nIndex(int x, int y, int z)
// this complex expression actually works
// the table looks like
// 0 : (x < 0)
// 1 : (x > 0)
// 2 : (y < 0)&&(x == 0)
// 3 : (y > 0)&&(x == 0)
// 4 : (z < 0)&&(y == 0)&&(x == 0)
// 5 : (z > 0)&&(y == 0)&&(x == 0)
// 6 : (z == 0)&&(y == 0)&&(x == 0)
    {return((((2*(z == 0) + (z < 0) + 2)*(y == 0) + (y < 0) + 2)*(x == 0)) + 
	    (x < 0));}

    int nRemainderX(int x, int y, int z)
    {
	return(x + (-1*(x > 0) + (x < 0)));
    }

    int nRemainderY(int x, int y, int z)
    {
	return(y + (x == 0)*(-1*(y > 0) + (y < 0)));
    }

    int nRemainderZ(int x, int y, int z)
    {
	return(z + ((x == 0)&&(y == 0))*(-(z > 0) + (z < 0)));
    }


// removes the links between you and your neighbors
    void deleteNlinks();
    
  public:

    OctMeshLeaf(const T& val_in, OctMeshNode<T>* my_parent, 
		    int l, int childnum)
	:OctMeshNode<T>(val_in, my_parent, l, childnum)
    {
	for (int j = 0; j < 6; j++)
	    _neighbors[j] = NULL;
    }

    OctMeshLeaf(OctMeshNode<T>* my_clone)
	:OctMeshNode<T>(*(my_clone->value()), my_clone->parentRef(), 
			    my_clone->level(), my_clone->childNum())
    {
	for (int j = 0; j < 6; j++)
	    _neighbors[j] = NULL;
    }

    OctMeshNode<T>* becomeNode()
    {
	return(new OctMeshNode<T>(*(valueRef()), 
				      parentRef(), 
				      level(), 
				      childNum()));
    }

    OctMeshLeaf():OctMeshNode<T>()
	{
	    for (int j = 0; j < 6; j++)
		_neighbors[j] = NULL;
	}

//
// note: the current scheme for neighbors is not going to route around NULL
// pointers in the grid. There are schemes that could do this, but they 
// are slower (see below)
//
    OctMeshLeaf<T>* neighbor(int x, int y, int z);
    OctMeshLeaf<T>* neighbor(int i);

// 
// the connect will find a neighbor by either
// 1) connecting through the mesh or
// 2) working it's way up the tree
//
    virtual OctMeshNode<T>* getNeighborNode(Direct direct);
    virtual OctMeshLeaf<T>* getNeighborCreate(Direct direct);
    virtual OctMeshLeaf<T>* getNeighborLeaf(Direct direct);

    virtual ~OctMeshLeaf();
    virtual int promote(int d);
//    virtual int promoteInterp();

// connect this leaf to the other along the direction "direct"
    void connect(OctMeshLeaf<T>* other_leaf, Direct direct);
    
};

#endif
