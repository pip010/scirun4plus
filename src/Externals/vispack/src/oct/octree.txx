// ******************************************************************
// VISPack. Copyright (c) 1994-2000 Ross Whitaker rtw@utk.edu       *
// For conditions of distribution and use, see the file LICENSE.txt *
// accompanying this distribution.                                  *
// ******************************************************************


#include <math.h>


// version 04 (node+treeNode + leafNode)



/*--------------------------------------------------------*/
/* Octree methods                                         */
/*--------------------------------------------------------*/

/*--------------------------------*/
/* void Octree::_init(int n)      */
/* initializes an n-levels Octree */
/*--------------------------------*/ 
template<class T>
void Octree<T>::_init(int n, T zero){
  INSIDE_S=zero;
  lvl=n;
  lNextLeaf=n;
  //_level=n;
  root=new treeNode(0,-8,-8,0); 
  current=root;
  currentNext=root;
  treenode_number=1;leaf_number=0;leafnode_number=0;
  treenode_size=sizeof(treeNode);leaf_size=sizeof(T);
  leafnode_size=sizeof(leafNode<T>);

  mask1= 1<<(lvl-1);
  mask2= 1<<lvl;

  _width=mask2;_wi2=mask1;
  cur_X=0;cur_Y=0;cur_Z=0;
  cur_XNL=cur_YNL=cur_ZNL=0;
}
/*---------------------------------------------------*/
/*-- Octree<T>::Octree(int n, T zero)                */
/*-- constructor ------------------------------------*/
template<class T>
Octree<T>::Octree(int n, T zero){
  //_level=n;
  INSIDE_S=zero;
  lvl=n;lNextLeaf=n;
  root=new treeNode(0,-8,-8,0);
  treenode_number=1;leaf_number=0;leafnode_number=0;
  treenode_size=sizeof(treeNode);leaf_size=sizeof(T);
  leafnode_size=sizeof(leafNode<T>);
  current=root;
  currentNext=root;
  mask1= 1<<(lvl-1);  
  mask2= 1<<lvl;
  _width=mask2;_wi2=mask1;
  cur_X=0;cur_Y=0;cur_Z=0;
  cur_XNL=cur_YNL=cur_ZNL=0;
}
/*-----------------------------------------------------*/
/* void Octree::setOctree(T v,uint x,uint y,uint z)    */
/* set the leaf xyz at value v-------------------------*/  
/*-----------------------------------------------------*/  
template<class T>
void Octree<T>::setOctree(T v,uint x,uint y,uint z){
  int i;
 
  Node *par=current;
  treeNode *parT;leafNode<T> *parL;

  //find common parent of current- and asked-leaf
  while(((x & mask2) != (cur_X & mask2)) ||
	((y & mask2) != (cur_Y & mask2)) ||
	((z & mask2) != (cur_Z & mask2))|| (mask1==1))
    {
      par=par->parent;
      mask1 = mask1<<1; mask2 -= mask1;
    }
  //cerr<<"a";
  //gets to the asked leaf.if no way exists, creates it
  while(1){
    i=0;
    if (z & mask1) i+=4;
    if (y & mask1) i+=2;
    if (x & mask1) i+=1;
    if ((mask1 & 1)==0) {
      parT=static_cast<treeNode*>(par);
      if (parT->pcells[i]==0){ 
	if((parT->n_region_sons[i])== 8 ) 
	  {	   
	    if(mask1 & 2) {
	      parT->pcells[i]=new leafNode<T>(i,8,8,par);
	      leafnode_number++;
	    }else {
	      parT->pcells[i]=new treeNode(i,8,8,par);
	      treenode_number++;
	    }	    
	    par->n_region--;
	    if(mask1<_wi2) {	      
	      parT->parent->n_region_sons[par->index]--;
	    }
	  } else if ((parT->n_region_sons[i])==-8 ) 
	    {
	      if(mask1 & 2) {
		parT->pcells[i]=new leafNode<T>(i,-8,-8,par);
		leafnode_number++;
	      }else {
		parT->pcells[i]=new treeNode(i,-8,-8,par);
		treenode_number++;
	      }	    	     
	      parT->n_region++;
	      if(mask1<_wi2) {	      
		parT->parent->n_region_sons[par->index]++;
	      }	      
	    } 
	else {
	  cerr<<"Bug setOct\n";exit(0);
	  parT->pcells[i]=new Node(i,0,0,par);	  	  
	}	
 	par=parT->pcells[i];
      }
      else{
 	par=parT->pcells[i];
      }
      mask2+=mask1;mask1=mask1>>1;
    }
    else break;
  }
  
  parL=static_cast<leafNode<T>*>(par);
  if (parL->values==0){
    parL->values=new T[8];
    leaf_number+=8;
    for(int j=0;j<8;j++) 
      parL->values[j]=INSIDE_S*((int)(parL->n_region_sons[j])/8);
  }
  //cerr<<"c";
//sets the leafs value.
  cur_X=x; cur_Y=y ; cur_Z=z;
  parL->values[i]=v;
  //cerr<<"b";
  parL->n_region-=(int)(par->n_region_sons[i]/8);
  parL->parent->n_region_sons[par->index]-=(int)(par->n_region_sons[i]/8);
  parL->n_region_sons[i]=0;
  //cerr<<"a"; 
  current=par;

 
} 

//------------------------------------------------------------------
//---  T& Octree<T>::poke(uint x, uint y, uint z)                ---
//---  returns a ref to the T value so that you can modify it .  ---
//------------------------------------------------------------------ 


template<class T>
T& Octree<T>::poke(uint x, uint y, uint z){
  int i; 
  Node *par=current;
  treeNode *parT;leafNode<T> *parL;

  //find common parent of current- and asked-leaf
  while(((x & mask2) != (cur_X & mask2)) ||
	((y & mask2) != (cur_Y & mask2)) ||
	((z & mask2) != (cur_Z & mask2))|| (mask1==1))
    {
      par=par->parent;
      mask1 = mask1<<1; mask2 -= mask1;
    }  
  //gets to the asked leaf.if no way exists, creates it
  while(1){
    i=0;
    if (z & mask1) i+=4;
    if (y & mask1) i+=2;
    if (x & mask1) i+=1;
    if ((mask1 & 1)==0) {
      parT=static_cast<treeNode*>(par);
      if (parT->pcells[i]==0){ 
	if((parT->n_region_sons[i])== 8 ) 
	  {	   
	    if(mask1 & 2) {
	      parT->pcells[i]=new leafNode<T>(i,8,8,par);
	      leafnode_number++;
	    } else {
	      parT->pcells[i]=new treeNode(i,8,8,par);
	      treenode_number++;
	    }
	    par->n_region--;
	    if(mask1<_wi2) {	      
	      parT->parent->n_region_sons[par->index]--;
	    }
	  } else if ((parT->n_region_sons[i])==-8 ) 
	    {
	      if(mask1 & 2){
		parT->pcells[i]=new leafNode<T>(i,-8,-8,par);
		leafnode_number++;
	      }else {
		parT->pcells[i]=new treeNode(i,-8,-8,par);
		treenode_number++;
	      }	    	     
	      parT->n_region++;
	      if(mask1<_wi2) {	      
		parT->parent->n_region_sons[par->index]++;
	      }	      
	    } 
	else {
	  parT->pcells[i]=new Node(i,0,0,par);
	  cerr<<"Bug setOct\n";exit(0);	  
	}	
 	par=parT->pcells[i];
      }
      else{
 	par=parT->pcells[i];
      }
      mask2+=mask1;mask1=mask1>>1;
    }
    else break;
  }
  parL=static_cast<leafNode<T>*>(par);
  if (parL->values==0){
    parL->values=new T[8];
    leaf_number+=8;
    for(int j=0;j<8;j++) 
      parL->values[j]=INSIDE_S*((int)(parL->n_region_sons[j])/8);
  }
 
//return ref to the leafs value.
  cur_X=x; cur_Y=y ; cur_Z=z;
  parL->n_region-=(int)(par->n_region_sons[i]/8);
  parL->parent->n_region_sons[par->index]-=(int)(par->n_region_sons[i]/8);
  parL->n_region_sons[i]=0;
  current=par;
  return parL->values[i];

 
}

/*-------------------------------------------------------*/
/* void Octree::setToRegion(short v,uint x,uint y,uint z)*/
/* sets a point to "inside/outside"->0                   */
/*-------------------------------------------------------*/
template<class T>
void Octree<T>::setToRegion(const int sgn,uint x,uint y,uint z){
  int i;
  
  Node *par=current;
  treeNode *parT;
  leafNode<T> *parL;
 
  // first find the nearest parent of current- and asked-position
  while(((x & mask2) != (cur_X & mask2)) ||
	((y & mask2) != (cur_Y & mask2)) ||
	((z & mask2) != (cur_Z & mask2))||(mask1==1))
    {
      par=par->parent;
      mask1 = mask1<<1; mask2 -= mask1;
    }
  cur_X=x;cur_Y=y;cur_Z=z;
 
// //find the asked leaf
  while(1){
    i=0;
    if (z & mask1) i+=4;
    if (y & mask1) i+=2;
    if (x & mask1) i+=1;
    if ((mask1 & 1)==0) {
      parT=static_cast<treeNode*>(par);
      if (parT->pcells[i]==0) {
 
	if (parT->n_region_sons[i]== 8 *sgn)
	  {
	    current=par;
	    return;
	  }else if (parT->n_region_sons[i] == -8 * sgn){
	    if(mask1==2){
	      parT->pcells[i]=new leafNode<T>(i,-8*sgn,-8*sgn,par);
	      leafnode_number++;
	    }else{
	      parT->pcells[i]=new treeNode(i,-8*sgn,-8*sgn,par);
	      treenode_number++;
	    }	  
	    parT->n_region+=sgn;
	    if(mask1<_wi2) 
	      parT->parent->n_region_sons[par->index]+=sgn;	  
	  } else {
	    cerr<<"Bug";exit(0);
	    parT->pcells[i]=new Node(i,0,0,par);
	    //node_number++;
	  }
	mask2+=mask1;mask1=mask1>>1;
 	par=parT->pcells[i];
      }
      else{	
	mask2+=mask1;mask1=mask1>>1;
   	par=parT->pcells[i];
      }
    }
    else break;
  }
  //cerr<<".";
  par->n_region+=sgn-(int)par->n_region_sons[i]/8;
  par->parent->n_region_sons[par->index]+=
    sgn-(int)par->n_region_sons[i]/8;
  par->n_region_sons[i]=8*sgn;
  parL=static_cast<leafNode<T> *>(par);
  if(parL->values!=0) parL->values[i]=INSIDE_S*(sgn);

  while((abs(par->n_region)==8)&&(mask1<_wi2)){
 
    i=par->index;
    if((mask1==1)&&(parL->values!=0)){
      delete[] parL->values;
      parL->values=0;
      leaf_number-=8;
    }
    par->parent->n_region_sons[i]=8*sgn;
    par=par->parent;parT=static_cast<treeNode*>(par);
    delete parT->pcells[i];

    parT->pcells[i]=0;
    if(mask1==1) leafnode_number--; else treenode_number--;
    mask1=mask1<<1;
    mask2-=mask1;
    par->n_region+=sgn;
  
    if(mask1<_wi2) par->parent->n_region_sons[par->index]+=sgn;
  }
  current=par; 
 
}
//------------------------------------------------------------------
//----    this ones for setting whole tree branches to regions    --
//----                                   purpose:faster init)     --
//------------------------------------------------------------------
template<class T>
void Octree<T>::setToRegion(const int sgn,uint x,uint y,uint z,uint l){
  int i;
  
  Node *par=current;
  treeNode *parT;
  leafNode<T> *parL;
  if((1<<l)>_wi2){ cerr<<" A VOIR \n"<<_wi2<<" "<<(1<<l)<<endl;exit(0);}
  else if(l==0) {setToRegion(sgn,x,y,z);return;}
  // first find the nearest parent of current- and asked-position
  while(((x & mask2) != (cur_X & mask2)) ||
	((y & mask2) != (cur_Y & mask2)) ||
	((z & mask2) != (cur_Z & mask2))||(mask1==1)||(mask1<(1<<l)))
    {
      par=par->parent;
      mask1 = mask1<<1; mask2 -= mask1;
    }
  cur_X=x;cur_Y=y;cur_Z=z;

 //find the asked leaf
  while(1){
    i=0;
    if (z & mask1) i+=4;
    if (y & mask1) i+=2;
    if (x & mask1) i+=1;
    if ((mask1 & (1<<l))==0) {
      parT=static_cast<treeNode*>(par);
      if (parT->pcells[i]==0) { 
	if (parT->n_region_sons[i]== 8 *sgn)
	  {
	    current=par;
	    return;
	  }else if (parT->n_region_sons[i] == -8 * sgn){
	    if(mask1==2){
	      parT->pcells[i]=new leafNode<T>(i,-8*sgn,-8*sgn,par);
	      leafnode_number++;
	    }else{
	      parT->pcells[i]=new treeNode(i,-8*sgn,-8*sgn,par);
	      treenode_number++;
	    }	  
	    parT->n_region+=sgn;
	    if(mask1<_wi2) 
	      parT->parent->n_region_sons[par->index]+=sgn;	  
	  } else {
	    cerr<<"Bug";exit(0);
	    parT->pcells[i]=new Node(i,0,0,par);
	    //node_number++;
	  }
	mask2+=mask1;mask1=mask1>>1;
 	par=parT->pcells[i];
      }
      else{	
	mask2+=mask1;mask1=mask1>>1;
   	par=parT->pcells[i];
      }
    }
    else break;
  }

  par->n_region+=sgn-(int)par->n_region_sons[i]/8;

  if(mask1<_wi2)
    par->parent->n_region_sons[par->index]+=  
      sgn-(int)par->n_region_sons[i]/8;

  par->n_region_sons[i]=8*sgn;

  if(mask1==1){
    parL=static_cast<leafNode<T> *>(par);
    if(parL->values!=0) parL->values[i]=INSIDE_S*(sgn);
  }

  if(mask1==(1<<l)){  
    parT=static_cast<treeNode *>(par);
    if(parT->pcells[i]!=0) {
      deleteTree(parT->pcells[i], mask1>>1);
      parT->pcells[i]=0;
    }
   
  } 

  while((abs(par->n_region)==8)&&(mask1<_wi2)){
    //cerr<<mask1<<" ";
    i=par->index;
    if((mask1==1)&&(parL->values!=0)){
      delete[] parL->values;
      parL->values=0;
      leaf_number-=8;
      leafnode_number--;
    }else{treenode_number--;}


    par->parent->n_region_sons[i]=8*sgn;
    par=par->parent;parT=static_cast<treeNode*>(par);
    //cerr<<par<<" ";
    delete parT->pcells[i];parT->pcells[i]=0;
    // if(mask1==2) leafnode_number--; else treenode_number--;
    mask1=mask1<<1;
    mask2-=mask1;
    par->n_region+=sgn;
    //cerr<<"T "<<par->parent;;
    if(mask1<_wi2) par->parent->n_region_sons[par->index]+=sgn;
    //cerr<<mask1<<" ";
  }
  current=par; 
 
}
 
/*------------------------------------------------*/
/* T Octree::getValue(uint x,uint y,uint z)       */
/* gets the value of the position x,y,z in Octree */
/*------------------------------------------------*/
template<class T>
T Octree<T>::getValue(uint x,uint y,uint z){
  T v;
  Node *par=current;
  leafNode<T>* parL;
  treeNode* parT;
  if((x >= _width)|| (y >= _width) || (z >= _width)) return INSIDE_S*-1;
  uint i;
//  x<<=1;y<<=1;z<<=1;
  //find common parent of current- and asked-leaf
  while(((x & mask2) != (cur_X & mask2)) ||
	((y & mask2) != (cur_Y & mask2)) ||
	((z & mask2) != (cur_Z & mask2)))
    {
      par=par->parent;
      mask1 = mask1<<1; mask2 -= mask1;
    }
  cur_X=x;cur_Y=y;cur_Z=z;
  //find the asked leaf.return +/- INSIDE_S if it doesn't exist
  while(1){
    i=0;
    if (z & mask1) i+=4;
    if (y & mask1) i+=2;
    if (x & mask1) i+=1;
    //cerr<<"MIgV"<<mask1<<" "<<i<<" ";
    if ((mask1 & 1)==0) {
      parT=static_cast<treeNode*>(par);
      if (parT->pcells[i]==0) {
	v=INSIDE_S*(par->n_region_sons[i]/8);
	current=par;
	return v;
      }
      else{
	mask2+=mask1;mask1=mask1>>1;
   	par=parT->pcells[i];  
      }
    }
    else break;
  }
  parL=static_cast<leafNode<T> *>(par);
  if(parL->values!=0) v= parL->values[i];
  else v= INSIDE_S*(parL->n_region_sons[i]/8);
  current=par;
  //mask1 = mask1<<1; mask2 -= mask1;
  return v;
  
}

//-------------------------------------------------------------------
//------- Octree copyTree();                                       --
//------- makes a copy of a tree.                                  --
//-------------------------------------------------------------------
template <class T>
Octree<T>* Octree<T>::copyTree(){
  Octree<T> *cpTree=new Octree<T>(lvl , INSIDE_S);
  cpTree->treenode_number=treenode_number;
  cpTree->leafnode_number=leafnode_number;
  cpTree->leaf_number=leaf_number;
  Node *tmp=cpy(root,_wi2,0); 
  cpTree->root=static_cast<treeNode* >(tmp);
  cpTree->current=cpTree->root;

  return cpTree;
}

//--------------------------------------------------------
//-submethod  of copyTree                               --
//--------------------------------------------------------
template<class T>
Node* Octree<T>::cpy(Node* pN, uint m1, Node* prt){
  if(m1>1){
    treeNode* b=new treeNode(pN->index, pN->n_region, 0, prt);   
    treeNode* par= static_cast<treeNode* >(pN);
    if(m1 & 2)
      for(int j=0;j<8;j++){
	b->n_region_sons[j]=pN->n_region_sons[j];
	if(par->pcells[j]!=0){
	 
	  b->pcells[j]=cpy(par->pcells[j],m1>>1,b);
	}else b->pcells[j]=0;
      }
    else
      for(int j=0;j<8;j++){
	b->n_region_sons[j]=pN->n_region_sons[j];
	if(par->pcells[j]!=0){
	  b->pcells[j]=cpy(par->pcells[j],m1>>1, b);
	} else b->pcells[j]=0;
      }  
    
    return b;
  }else{
    leafNode<T>* b=new leafNode<T>(pN->index, pN->n_region, 0, prt);
    leafNode<T>* par=static_cast<leafNode<T>* >(pN);
    if(par->values!=0){
      b->values=new T[8];
      for(int j=0;j<8;j++){
	b->n_region_sons[j]=par->n_region_sons[j];
	b->values[j]=par->values[j];  
      }
    }else { for(int j=0;j<8;j++) b->n_region_sons[j]=par->n_region_sons[j]; }
    return b;
  }
}

//-----------------------------------------------------------------
//---- Octree<T>& Octree<T>::operator=(const Octree<T>& t)       --
//----          copy operator                                    --
//-----------------------------------------------------------------
template<class T>
Octree<T>& Octree<T>::operator=(const Octree<T>& t)
{
    if (this != &t){     
      delete this->root;   
      this->_init(t.lvl, t.INSIDE_S);
      this->treenode_number=t.treenode_number;
      this->leafnode_number=t.leafnode_number;
      this->leaf_number=t.leaf_number; 
      //this->root=t.root->cpy<T>(t._wi2,0);   
      this->root=cpy(t.root,t._wi2,0);  
      this->current=this->root;
    } 
    return *this;
}

//-----------------------------------------------------------------------
//--  Octree<T>* Octree<T>::getSubTree(uint x, uint y, uint z, uint l) --
//--       returns subtree.                                            --
//-----------------------------------------------------------------------
template<class T>
Octree<T>* Octree<T>::getSubTree(uint x, uint y, uint z, uint l){
  Octree<T> *subTree=new Octree<T>(l,INSIDE_S);
  if((l<2)||(l>lvl)) {cerr<<"too small .exiting.\n";exit(0);}
  int i; 
  Node *par=current;
  treeNode *parT;
  //find common parent of current- and asked-leaf
  while(((x & mask2) != (cur_X & mask2)) ||
	((y & mask2) != (cur_Y & mask2)) ||
	((z & mask2) != (cur_Z & mask2))|| (mask1<(1<<(l-1))))
    {
      par=par->parent;
      mask1 = mask1<<1; mask2 -= mask1;
    }
 
  cur_X=x;cur_Y=y;cur_Z=z;
  //find the asked leaf.
  while(1){
    i=0;
    if (z & mask1) i+=4;
    if (y & mask1) i+=2;
    if (x & mask1) i+=1;  
    if ((mask1 & (1<<(l-1)))==0) {
      parT=static_cast<treeNode*>(par);
      if (parT->pcells[i]==0) {
	current=par;
	if(par->n_region_sons[i]<0) return subTree;
	else{
	  for(int j=0;j<8;j++) subTree->root->n_region_sons[j]=8;
	  subTree->root->n_region=8;
	  return subTree;
	}	
      }
      else{
	mask2+=mask1;mask1=mask1>>1;
   	par=parT->pcells[i];  
      }
    }
    else break;
  }
  current=par;
  cerr<<"test sub";cerr<<mask1<<endl;
  Node *b;
//  subTree->root=parT->cpy<T>(mask1,0);
//  subTree->root=cpy(parT,mask1,0);
  b=cpy(par,mask1,0);
  subTree->root=static_cast<treeNode *>(b); 
  subTree->current=subTree->root;
  subTree->nodeCountReset();
  return subTree;
}


//------------------------------------------------------------------------
//-------   void Octree<T>::deleteTree(treeNode *pT,int l)              --
//-------            deletes the tree                                   --
//------------------------------------------------------------------------
template<class T>
void Octree<T>::deleteTree(Node *pN,int m1){
  treeNode *pT;
  if(pN!=0){
    pT=static_cast<treeNode *>(pN);
    if(m1>2){
      treeNode *ptmp;
      for(int j=0; j<8; j++){
	if(pT->pcells[j]!=0) ptmp=static_cast<treeNode* >(pT->pcells[j]);
	else ptmp=0;
	deleteTree(ptmp,m1>>1);
      }
      delete pT;treenode_number--;
      pT=0;
    }else {
      for(int j=0;j<8;j++){
	if(pT->pcells[j]!=0) {
	  leafNode<T> *ptmp=static_cast<leafNode<T>* >(pT->pcells[j]);
	  if(ptmp->values!=0) {
	    delete[] ptmp->values;
	    leaf_number-=8;
	  }
	  delete ptmp;	ptmp=0;leafnode_number--;
	}
      }
      delete pT;pT=0;treenode_number--;
    }
    pN=0;
  }
}
//----------------------------------------------------------------
//----  Octree<T>::~Octree()                                    --
//----         destructor                                       --
//----------------------------------------------------------------
template<class T>
Octree<T>::~Octree(){
  deleteTree(root, _wi2);
  root=0;current=0;
}


//-----------------------------------------------------------------------
//-------  void Octree<T>::setBorder(T v, int wd)                      --
//-------         sets the border of the cube to value                 --
//-----------------------------------------------------------------------
template<class T>
void Octree<T>::setBorder(T v, int wd){
  if (v==INSIDE_S){
    for(int i=0;i<wd;i++)
      for(int j=0;j<_width;j++)
	for(int k=0;k<_width;k++){
	  setToRegion(1,i,j,k);
	  setToRegion(1,_width-i-1,j,k);
	  setToRegion(1,j,k,i);
	  setToRegion(1,j,k,_width-i-1);
	  setToRegion(1,k,i,j);
	  setToRegion(1,k,_width-i-1,j);
	}
  }else if (v==INSIDE_S*-1){
    for(int i=0;i<wd;i++)
      for(int j=0;j<_width;j++)
	for(int k=0;k<_width;k++){
	  setToRegion(-1,i,j,k);
	  setToRegion(-1,_width-i-1,j,k);
	  setToRegion(-1,j,k,i);
	  setToRegion(-1,j,k,_width-1-i);
	  setToRegion(-1,k,i,j);
	  setToRegion(-1,k,_width-i-1,j);
	}
  }else{
     for(int i=0;i<wd;i++)
      for(int j=0;j<_width;j++)
	for(int k=0;k<_width;k++){
	  setOctree(v,_width-i-1,j,k);
	  setOctree(v,i,j,k);
	  setOctree(v,j,k,_width-i-1);
	  setOctree(v,j,k,i);
	  setOctree(v,k,_width-i-1,j);
	  setOctree(v,k,i,j);
	}
  }
}


//-----------------------------------------------------------------


template<class T>
void Octree<T>::clean(){
  cout<<_width<<endl;
  for(int k=0;k<_width;k++){
    //cerr<<".";
    for(int j=0;j<_width;j++)
      for(int i=0;i<_width;i++){
	setToRegion(-1, i,j,k);
      }
  }
}

//-----------------------------------------------
//--- counts the number of nodes and leaves.   --
//---   just for use in getSubTree             --
//-----------------------------------------------
template<class T>
void Octree<T>::nodeCountReset(){
  treenode_number=0;
  leafnode_number=0;
  leaf_number=0;
  nodeCounter(root, _wi2);
}

template<class T>
void Octree<T>::nodeCounter(Node *pN,uint m1){
  if(pN!=0)    
    if(m1>1){
      treenode_number++;
      treeNode *pT=static_cast<treeNode *>(pN);
      for(int j=0; j<8;j++) nodeCounter(pT->pcells[j],m1>>1);
    }else{
      leafnode_number++;
      leafNode<T> *pL=static_cast<leafNode<T> *>(pN);
      if(pL->values!=0)
	leaf_number+=8;
    }
}

//-------------------------------------------------------------------------
//-- void Octree<T>::insertTree(Octree<T>* insT, uint x, uint y, uint z) --
//--       inserts the insT tree in the Octree at cell containing xyz    --
//-------------------------------------------------------------------------
template<class T>
void Octree<T>::insertTree(Octree<T>* insT, uint x, uint y, uint z){
  if((insT->_width >= _width)||(insT->lvl<2)) {
    
    cerr<<"Size problem in tree insertion.exiting\n";
    cerr<<insT->_width<<" "<<_width<<" "<<insT->lvl<<endl;
    exit(0);
  } 
  int l=insT->lvl;
  int i; 
  Node *par=current;
  treeNode *parT;
  //find common parent of current- and asked-leaf
  while(((x & mask2) != (cur_X & mask2)) ||
	((y & mask2) != (cur_Y & mask2)) ||
	((z & mask2) != (cur_Z & mask2))|| (mask1<(1<<(l))))
    {
      par=par->parent;
      mask1 = mask1<<1; mask2 -= mask1;
    }
  cur_X=x;cur_Y=y;cur_Z=z;
  while(1){
    i=0;
    if (z & mask1) i+=4;
    if (y & mask1) i+=2;
    if (x & mask1) i+=1;  
    if ((mask1 & ( 1<<l ) )==0) {
      parT=static_cast<treeNode*>(par);
      if (parT->pcells[i]==0) {
	current=par;
	if((parT->n_region_sons[i])== 8 ) 
	  {	   
	    if(mask1 & 2) {
	      parT->pcells[i]=new leafNode<T>(i,8,8,par);
	      leafnode_number++;
	    }else {
	      parT->pcells[i]=new treeNode(i,8,8,par);
	      treenode_number++;
	    }	    
	    par->n_region--;
	    if(mask1<_wi2) {	      
	      parT->parent->n_region_sons[par->index]--;
	    }
	  } else if ((parT->n_region_sons[i])==-8 ) 
	    {
	      if(mask1 & 2) {
		parT->pcells[i]=new leafNode<T>(i,-8,-8,par);
		leafnode_number++;
	      }else {
		parT->pcells[i]=new treeNode(i,-8,-8,par);
		treenode_number++;
	      }	    	     
	      parT->n_region++;
	      if(mask1<_wi2) {	      
		parT->parent->n_region_sons[par->index]++;
	      }	      
	    } 
	else {
	  cerr<<"Bug setOct\n";exit(0);		  
	}	
 	par=parT->pcells[i];
      }
      else{
 	par=parT->pcells[i];
      }
      mask2+=mask1;mask1=mask1>>1;
    }
    else break;
  }
  parT=static_cast<treeNode *>(par);
  if(parT->pcells[i]!=0) deleteTree(parT->pcells[i], mask1>>1);
/*   parT->pcells[i]=insT->root->cpy<T>(mask1>>1,par); */
  parT->pcells[i]=cpy(insT->root,mask1>>1,par);
  parT->n_region_sons[i]=insT->root->n_region;
  parT->pcells[i]->index=i;
  

//sets the leafs value.
  cur_X=x; cur_Y=y ; cur_Z=z;

  par->n_region-=(int)(par->n_region_sons[i]/8);
  if(!(mask1 & _wi2))
    par->parent->n_region_sons[par->index]-=(int)(par->n_region_sons[i]/8);
  par->n_region_sons[i]=0; 
  current=par;
  nodeCountReset();

 
}


template<class T>
void Octree<T>::writeTree(char* treeName, int prec){
  ofstream file;

  file.open(treeName);
  file.setf(ios::scientific, ios::floatfield);
  file.precision(prec);

  file<<"VISPack octree file\n"
      <<"level "<<lvl<<endl;
  //file<<root->n_region<<" "<<root->index<<" , ";
  file<<root->n_region<<" , ";
  for(int j=0;j<8; j++) file<<root->n_region_sons[j]<<" ";
  file<<endl;
  for(int j=0; j<8; j++){
    if(root->pcells[j]!=0){
      file<<j<<" { ";cerr<<"<";
      writeNode(root->pcells[j], file, _wi2>>1);
      //file<<" } ";
      cerr<<">";
      file<<" } ";
    }
  }
  file<<" }\n";
  
  file.close();
}

template<class T>
void Octree<T>::writeNode(Node* par, ofstream &file, uint m1){

  if(m1>1){
    
    treeNode* pT=static_cast<treeNode *>(par);
    //file<<pT->n_region<<" "<<pT->index<<" , ";
    file<<pT->n_region<<" , ";
    for(int j=0;j<8; j++) file<<pT->n_region_sons[j]<<" ";
    file<<endl;
    for(int j=0; j<8; j++){
      if(pT->pcells[j]!=0){
	file<<j<<" { ";
	 writeNode(pT->pcells[j], file, m1>>1);
	
	file<<" } ";
      }
    }
  }else{
    leafNode<T>* pL=static_cast<leafNode<T> *>(par);
    //file<<pL->n_region<<" "<<pL->index<<" , ";
    file<<pL->n_region<<" , ";
    for(int j=0;j<8; j++) file<<pL->n_region_sons[j]<<" ";
    file<<endl;
    if(pL->values!=0)
      for(int j=0;j<8;j++) file<<pL->values[j]<<" ";
  }
  file<<endl;
}

template<class T>
void Octree<T>::readTree(char* treeName){
  char s[80];
  ifstream file;
  int a;char c;
  file.open(treeName);
  file.getline(s, 50);
  if(strcmp(s,"VISPack octree file")!=0){ cerr<<"Bad File Type\n";exit(0);}
  file>>s>>a; cerr<<"level : "<<a<<endl;
  _init(a, INSIDE_S);
  file>>a;root->n_region=a; //cerr<<"n_reg "<<(int)root->n_region<<" ";
  //file>>a; 
  root->index=0; //cerr<<"n_reg "<<(int)root->index<<" ";
  file>>c; //cerr<<"n_reg "<<c<<" ";
  for(int j=0; j<8; j++) {
    file>>a;root->n_region_sons[j]=a; 
    cerr<<" "<<(int)root->n_region_sons[j]; 
  }
  file>>c; 
  
  while(c!= '}'){
    file.putback(c); file>>a;//cerr<<"PP "<<a<<" "<<_wi2<<endl;
    root->pcells[a]=new treeNode(a, 0,0,root);
    readNode(root->pcells[a], file, _wi2>>1);cerr<<">";
    file>>c;
  }
  file.close();
  nodeCountReset();
}

template<class T>
void Octree<T>::readNode(Node* par, ifstream &file, uint m1){  
  char c; int a;
  file>>c;
  // cerr<<"test "<<c<<" "<<m1<<endl;
  if(c!='{') {cerr<<"bug\n";exit(0);}
  if(m1>1){
    treeNode* pT=static_cast<treeNode *>(par);
    //cerr<<"=+";
    file>>a;pT->n_region=a;
    //file>>a>>c;
    file>>c;
    for(int j=0;j<8;j++) {
      file>>a;
      pT->n_region_sons[j]=a; 
    }
    file>>c;
    if(m1>2){ 
      while(c!='}'){
	file.putback(c);file>>a;
	pT->pcells[a]=new treeNode(a,0,0,par);
	readNode(pT->pcells[a], file, m1>>1);
	
	file>>c;
      }
    }else{
      while(c!='}'){
	file.putback(c);file>>a;
	pT->pcells[a]=new leafNode<T>(a,0,0,par);
	readNode(pT->pcells[a], file, m1>>1);
	file>>c;
      }
    }
  
  }else{
    //cerr<<"leaf ";
    leafNode<T> *pL=static_cast<leafNode<T>* >(par);
    file>>a;par->n_region=a;
    //file>>a>>c;
    file>>c;
    for(int j=0;j<8;j++) {
      file>>a;
      pL->n_region_sons[j]=a; 
    }
    file>>c;
    if(c!='}'){
      pL->values=new T[8];
      file.putback(c);
      for(int j=0;j<8;j++) file>>pL->values[j];
      file>>c;
    }
  }
}

template<class T>
/* int Octree<T>::getNextLeaf(int *px, int *py, int *pz){ */
int Octree<T>::getNextLeaf(int &px, int &py, int &pz){
  int valid=0;
  uint x=px, y=py, z=pz;
  x=cur_XNL;y=cur_YNL;z=cur_ZNL;
  Node *par=currentNext;
  leafNode<T>* parL;
  treeNode* parT;
  uint i=0;
  //find common parent of current- and asked-leaf
  if(lNextLeaf==1){
    i=0;
    if(1&x) i+=1;
    if(1&y) i+=2;
    if(1&z) i+=4;
    i+=1;
    x=x>>1;y=y>>1;z=z>>1;
  }
  int do_it=1;
  while(do_it){
    //cerr<<"1";
    // ascend tree until i<8
    while((i>7)&&(lNextLeaf<lvl)){
      //cerr<<"getup ";
      i=par->index+1;
      par=par->parent;if(par==0){cerr<<"bug parent\n";exit(0);}
      //mask1=mask1<<1;mask2-=mask1;
      lNextLeaf++;
      x=x>>1;y=y>>1;z=z>>1;
    }
    if((lNextLeaf==lvl)&&(i==8)){valid=0; break;}
    // descend tree or move right
    //cerr<<"par i mask1 "<<par<<" "<<i<<" "<<mask1<<endl;
    while(i<=7){
      //cerr<<" godown "<<lNextLeaf<<" "<<i<<" ";
      if(lNextLeaf>1){
	parT=static_cast<treeNode* >(par);
	if(parT->pcells[i]==0){    //empty son. move right
	  i=i+1;
	}else{                     //there is a son; update coords
	  x=(x<<1)|(i&1);
	  y=(y<<1)|((i>>1)&1);
	  z=(z<<1)|((i>>2)&1);
	  par=parT->pcells[i];     // treeNode : descend;
	  i=0;
	  //mask2+=mask1;mask1=mask1>>1;
	  lNextLeaf--;
	}
      }
      else{                  // leafNode; if there are values, return xyz;
	x=(x<<1)|(i&1);
	y=(y<<1)|((i>>1)&1);
	z=(z<<1)|((i>>2)&1);
	parL=static_cast<leafNode<T>* >(par);
	if(parL->values!=0){
	  px=x;py=y;pz=z;
	  do_it=0;valid=1;break;
	}else{
	  i=8;
	  x=x>>1;y=y>>1;z=z>>1;
	}
      }      
    }    
  }
  cur_XNL=x;cur_YNL=y;cur_ZNL=z;
  currentNext=par;
  if(valid) {px=x;py=y;pz=z;}
  else{
    currentNext=root;
    lNextLeaf=lvl;
    cur_XNL=cur_YNL=cur_ZNL=0;
  }
  return valid;
}
	  
