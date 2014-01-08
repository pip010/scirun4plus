// ******************************************************************
// VISPack. Copyright (c) 1994-2000 Ross Whitaker rtw@utk.edu       *
// For conditions of distribution and use, see the file LICENSE.txt *
// accompanying this distribution.                                  *
// ******************************************************************


// Antoine Pardigon. antoine.pardigon@m4x.org
//
// This file declares the node and octree classes.


//-----------------------------------------------------------------
// node,leafnode,treenode ( the leaf-node is templated) and octree
//-----------------------------------------------------------------
#ifndef __octree
#define __octree


#include <new.h>
#include <ioiostream>
#include <stdexcept>



const int INSIDE=1;
const int OUTSIDE=-1;
const int NO_ACTIVE=-1;

typedef unsigned int uint;
 
/*---------------------------------------*/
/* class Node. Defines an octree-node    */
/*---------------------------------------*/

class Node{
 public:  
  signed char index;   
                
  signed char n_region;
  signed char n_region_sons[8]; // sum of the numbers of 
                                //outside/inside nodes below
  Node *parent;             // pointer to the node above
  //constructors
  Node(){}
  Node(short i,short n,short n_sons, Node *p){
    index=i;
    n_region=n;
    parent=p;
    for(int j=0;j < 8;j++) {      
      n_region_sons[j]=n_sons;
    }
  }

//destructors
  ~Node() { }
 
};


class treeNode : public Node {

 public:
  Node *pcells[8];        // pointers to the cells below
  treeNode(){}
  treeNode(short i,short n,short n_sons, Node *p){    
    index=i;
    n_region=n;
    parent=p;
    for(int j=0;j < 8;j++) {
      pcells[j]=0;
      n_region_sons[j]=n_sons;
    }
  }
  

//destructors
  ~treeNode() {}

};

template<class T>
class leafNode : public Node {
// protected:
 public:
  T *values;
    //constructors
  leafNode(){}
  leafNode(short i,short n,short n_sons, Node *p){
    index=i;
    n_region=n;
    parent=p;
    values=0;
    for(int j=0;j < 8;j++) {      
      n_region_sons[j]=n_sons;
    }
  }

  //destructors
  ~leafNode() { }

};




/*-------------------------------------------------*/
/* CLASS OCTREE : Defines an Octree and operations */
/*-------------------------------------------------*/

template<class T> 
class Octree : public Node {
public:
  Octree() {};
  Octree(int n,T zero );           // initializes an octree of n-levels
  void _init(int n, T zero);       // with T as inside value.
  ~Octree();                          
  

  
  T getValue(uint x,uint y,uint z);     //gets value at x,y,z

  T& setValue(uint x, uint y, uint z);  //returns reference at data xyz
  T& poke(uint x, uint y, uint z);      // same as setValue.

  void setOctree(T v,uint x,uint y,uint z);  //puts data v at leaf xyz

  void setToRegion(int sgn ,uint x,uint y,uint z);//sets cell xyz to +/- inside
                                                  //updates tree
  void setToRegion(int sgn ,uint x,uint y,uint z,uint l);

  void setBorder(T v, int wd);          //sets border of volume to value v.

  Octree* getSubTree(uint x, uint y, uint z, uint l); // returns subtree of level 
                                                      // l containing cell xyz.
  void insertTree(Octree *insT, uint x, uint y, uint z);//inserts tree at position xyz
  
  Octree* copyTree();                                 //deep copy of the tree
  Octree& operator=(const Octree& t);

  const T& operator()(uint x, uint y, uint z){       
    return getValue(x,y,z);
  }

  int checkBounds(int i, int j, int k){             //checks if asked point is in tree 
    if((i<0)||(i>=_width)||(j<0)||(j>=_width)||(k<0)||(k>=_width))
      return 0;
    else return 1;
  }
  int width() { return _width; }                     // returns width
  int level() { return lvl;}                         //returns level
  int getSize(){return treenode_number*treenode_size+ //returns size of tree(in bytes)
		  leafnode_number*leafnode_size+
		  leaf_number*leaf_size;}
  int getTreeSize(){                                   //returns size of tree structure
    return treenode_number*treenode_size+leafnode_number*leafnode_size;
  }
  int getDataSize(){return leaf_number*leaf_size;}      //returns size of data only

  void sizeInfo(){ cerr<<"size of treenode :"<<treenode_size  //outputs sizes of treenode, leafnode and data to err.
		       <<"\nsize of leafnode :"<<leafnode_size
		       <<"\nsize of leaf :"<<leaf_size<<endl;
  }
  
  void clean();
  
 
  T getInside(){return INSIDE_S;} //returns the INSIDE value
 

  void writeTree(char* treeName){ writeTree(treeName,6);}  //writes tree to file
  void writeTree(char* treeName, int prec);                // with asked precision 
  void readTree(char* treeName);                           //reads tree from file
  
  int getNextLeaf(int &px, int &py, int &pz);  //gets coordinates of next leaf in the octree

protected:
// public:
  T INSIDE_S;                        // inside value of data
  treeNode *root;                    // root of the tree

  Node *current, *currentNext;       // last visited nodes 
                                     //(currentNext for getNextLeaf method)

  uint mask1,mask2;                

  uint cur_X,cur_Y,cur_Z,            //current point for getV, setV, ...
    cur_XNL, cur_YNL, cur_ZNL;       //current point for getNextLeaf;

  uint _wi2,_width;                  //width, half-width of the octree
  int lNextLeaf;

  int treenode_number,leafnode_number, leaf_number;
  int treenode_size,leafnode_size,leaf_size;

  int lvl;                            // number of levels of the tree.
  
  Node* cpy(Node* pN, uint m1, Node* prt);// copy of a cell of the tree
                                          //(and of the cells below)
  void deleteTree(Node *pN,int m1);         
  void nodeCounter(Node *pN,uint m1);
  void nodeCountReset();

  void readNode(Node* par, ifstream &file, uint m1);  
  void writeNode( Node* par,ofstream &file, uint m1);


};





#ifndef MANUAL_INSTANTIATION
#include "octree.c"
#endif

#endif /* __octree */ 
