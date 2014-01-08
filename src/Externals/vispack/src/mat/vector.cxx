// ******************************************************************
// VISPack. Copyright (c) 1994-2000 Ross Whitaker rtw@utk.edu       *
// For conditions of distribution and use, see the file LICENSE.txt *
// accompanying this distribution.                                  *
// ******************************************************************

// $Id: vector.cxx,v 1.1.1.1 2003/02/12 16:51:53 whitaker Exp $

#include <math.h>
#include "mat/matrix.h"
#include <fstream>
#include <iostream>

VISVector::VISVector(const VISMatrix& other)
{
  if (other.c() == 1)
    // this doesn't work
    //	VISMatrix(other);
    //
	{
	  _cols = other._cols; 
	  _rows = other._rows;
	  _matrixdata = other._matrixdata;
	  refCntInc();
	}
  else
	{
	  std::cout << "Attempted to create a vector from a"
    << "matrix that has _cols > 1" << std::endl;
	  exit(0);
	  //std::cout << (*((VISVector*)(NULL))) << flush;
	}
}

VISVector::VISVector(char* filename){
  std::ifstream file(filename);
  int r;
  file >> r;
  initialize(r,1);
  int i;
  for (i = 0; i < _rows; i++) {
    if (!file)
	    std::cerr << "error reading matrix: ran out of data" << std::endl;
    file >> poke(i);
  }
}

VISVector& VISVector::operator=(const VISMatrix& other)
{
  if (&other == this)
    return(*this);
  if (other.c() > 1)
	{
	  std::cout <<
    "You have attempted to assign a matrix(cols>1) to a vector." 
    << std::endl;
	  exit(0);
	  return *this;
	}
  else
	{
	  VISMatrix::operator=(other);
	  return *this;
	}
}


VISVector VISVector::cross(const VISVector& other) const
{
  if ((r()==3)&&(other.r()==3)){
    VISVector ans(3);
    ans.poke(0)=peek(1)*other.peek(2)-peek(2)*other.peek(1);
    ans.poke(1)=-1*(peek(0)*other.peek(2)-other.peek(0)*peek(2));
    ans.poke(2)=peek(0)*other.peek(1)-peek(1)*other.peek(0);
    //std::cout << "in cross " << ans << flush;
    //std::cout << ans.query(0) << " " << ans.query(1) << " "
    //<< ans.query(2) << std::endl;
    return ans;
    //	return *this;
  }
  else{
    std::cout << "VISVector VISVector::cross const(const VISVector& other)" << std::endl;
    std::cout << "Attempting cross product of vectors of sizes"
    << " other than 3." << std::endl;
    exit(0);
    return 0;
  }
}

datatype VISVector::dot(const VISVector& other) const
{
  datatype sum = 0.0;
  if (r() != other.r())
	{
	  std::cout << "You have attempted to dot two vectors of"
    << "unequal dimensions." << std::endl;
	  exit(0);
	  return 0;
	}
  else
	{
	  for (int i=0;i<other.r();i++)
      sum=sum+peek(i)*other.peek(i);
	  return sum;
	}
}

VISMatrix VISVector::exterior(const VISVector& other) const
{
  int rows = r();
  VISMatrix ret(rows, rows);
  
  if (rows !=other.r())
	{
	  std::cout << "You have attempted to exterior two vectors of"
    << "unequal dimensions." << std::endl;
	  exit(0);
	  return 0;
	}
  else
	{
	  for (int i=0; i < rows; i++)
	    for (int j=0; j< rows; j++)
	      ret.poke(i, j) = peek(j)*other.peek(i);
	}
  return(ret);
}

datatype VISVector::mean() const
{
  int rows = r();
  datatype tot = 0.0;
  for (int r = 0; r < rows; r++)
    tot = tot + peek(r);
  return ( tot / static_cast<datatype>(rows) );
}

VISMatrix VISVector::mat(int cols) const{
  //create a matrix from a vector by repeating
  //the vector across all the columns of the matrix
  //cols = number of columns in the matrix
  int num = n();
  VISMatrix ans(num,cols);
  for (int row = 0; row < num; row++)
    for (int col = 0; col < cols; col++)
	    ans.poke(row,col) = peek(row);
  return ans;
}

VISMatrix VISVector::rowMatrix(int rows) const{
  //create a matrix from a vector by repeating
  //the vector down all the rows of the matrix
  //rows = number of rows in the matrix
  int num=n();
  VISMatrix ans(rows,num);
  for (int row = 0; row < rows; row++)
    for (int col = 0; col < num; col++)
	    ans.poke(row,col) = peek(col);
  return ans;
}

VISVector& VISVector::pokeROI(int row, const VISVector& roi)
{
  int n = roi.n();
  for (int i = 0; i < n; i++)
    poke(i + row) = roi.peek(i);
  
  return(*this);
}

// this prints out the matrix is a nice way (as test) to an std::ostream object
void VISVector::print(std::ostream *os) const
{
#ifdef SGI
  *os << setw(3) << _rows << std::endl;
#else //ifdef LINUX
  *os << _rows <<" " <<_cols<<std::endl;
#endif
  int i;
  for (i = 0; i < _rows; i++){
#ifdef SGI
    *os << setw(7) << setprecision(7)  << peek(i) << std::endl;
#else //ifdef LINUX
    *os<< peek(i)<<" ";
#endif
  }
}

void VISVector::read(std::istream *is)
{
  //initialize this
  int r;
  *is >> r;
  initialize(r,1);
  int i;
  for (i = 0; i < _rows; i++) {
    if (!is)
	    std::cerr << "error reading matrix: ran out of data";
    *is >> poke(i);
  }
}

datatype* VISVector::createLAPK() const{
  datatype* lapk;
  int elements;
  elements = n();
  lapk = new datatype [elements];
  for (int i=0; i<elements; i++){
    lapk[i] = peek(i);
  }
  return lapk;
}//datatype* VISVector::createLAPK()

void VISVector::becomeLAPK(datatype* lapk, int n){
  if (refCnt() == 1){
    delete _matrixdata;
    refCntDec();
  }
  else if (refCnt() > 1){
    refCntDec();
  }
  initialize(n,1);
  for (int i=0; i<_rows; i++){
    poke(i) = lapk[i];
  }
}//void VISVector::becomeLAPK(datatype* lapk, int n)

VISMatrix VISVector::concat(const VISMatrix& second)
{
  //concatinate cols i.e. a=[1   b=[ 7  8  9
  //                         4 ]    10 11 12]
  // c=concatcols(a,b)   c=[1  7  8  9
  //                        4 10 11 12]
  if ((n()==0) &&  ((second.c()==0)&&(second.r()==0)))
	{
	  std::cout << "VISMatrix concat(const VISVector& "
    << "first,const VISMatrix& second)" << std::endl;
	  std::cout << "You have attempted to col concatinate a matrix "
    << "and a vector, both of which have zero size." << std::endl;
	  exit(0);
	  return 0;
  }
  else if (n()==0)
	{
	  VISMatrix ans;
	  ans=second;
	  return ans;
  }
  else if ((second.c()==0)&&(second.r()==0))
	{
	  VISMatrix ans;
	  ans=mat();
	  return ans;
  }	
  else if (n()!=second.r())
	{
	  std::cout << "VISMatrix concat(const VISVector& "
    << "first,const VISMatrix& second)" << std::endl;
	  std::cout << "You have attempted to column concatinate a matrix "
    << "and a vector that have different numbers of rows/elements."
    << std::endl;
	  std::cout << "first.n()  " << n() << "     "
    << "second.r()  " << second.r() << std::endl;
	  exit(0);
	  return 0;
  }
  else
	{
	  int numrows = n();
	  int numcols = second.c()+1;
	  int firstcols = 1;
	  VISMatrix ans(numrows,numcols);
	  for (int row=0;row<numrows;row++)
	    for (int col=0;col<numcols;col++)
        if(col<firstcols)
          ans.poke(row,col) = peek(row);
        else
          ans.poke(row,col) = second[row][col-firstcols];
	  return ans;
  }
}

VISMatrix VISVector::concat(const VISVector& second)
{
  //concatinate cols i.e. a=[1  b=[ 7 
  //                         4]    10]
  // c=concatcols(a,b)   c=[1  7 
  //                        4 10]
  if ((second.n()==0) && (n()==0))
	{
	  std::cout << "VISMatrix concat(const VISVector& "
    << "first,const VISVector& second)" << std::endl;
	  std::cout << "You have attempted to col concatinate vectors, both "
    << "of which have zero size." << std::endl;
	  exit(0);
	  return 0;
  }
  else if (n()==0)
	{
	  VISMatrix ans;
	  ans=second.mat();
	  return ans;
  }
  else if (second.n()==0)
	{
	  VISMatrix ans;
	  ans=mat();
	  return ans;
  }	
  else if (n()!=second.n())
	{
	  std::cout << "VISMatrix concat(const VISVector& "
    << "first,const VISVector& second)" << std::endl;
	  std::cout << "You have attempted to column concatinate a vector "
    << "and a vector that have different numbers of elements."
    << std::endl;
	  std::cout << "n()  " << n() << "     "
    << "second.n()  " << second.n() << std::endl;
	  exit(0);
	  return 0;
  }
  else
	{
	  int numrows = n();
	  int numcols = 2;
	  int firstcols = 1;
	  VISMatrix ans(numrows,numcols);
	  for (int row=0;row<numrows;row++)
	    for (int col=0;col<numcols;col++)
        if(col<firstcols)
          ans.poke(row,col) = peek(row);
        else
          ans.poke(row,col) = second[row];
	  return ans;
  }
}
VISVector VISVector::concatRow(const VISVector& second)
{
  //concatinate rows i.e. a=[1  b=[ 7
  //                         4]    10]
  // c=concatRows(a,b)   c=[1
  //                        4 
  //                        7
  //                       10]
  if ((n()==0) && (second.n()==0))
	{
	  std::cout << "VISMatrix concatRows(const VISMatrix& "
    << "first,const VISVector& second)" << std::endl;
	  std::cout << "You have attempted to row concatinate two vectors"
    << ", both of which are zero size "
    << std::endl;
	  exit(0);
	  return 0;
  }
  else if (n()==0)
	{
	  VISVector ans;
	  ans=second;
	  return ans;
  }
  else if (second.n()==0)
	{
	  VISVector ans;
	  ans=*this;
	  return ans;
  }	
  else
	{
	  int numrows = n()+second.n();
	  int firstrows = n();
	  VISVector ans(numrows);
	  for (int row=0;row<numrows;row++)
	    if(row<firstrows)
        ans.poke(row) = peek(row);
	    else
        ans.poke(row) = second[row-firstrows];
	  return ans;
  }
}



/*HVISMatrix quaternion(VISVector q){
 HVISMatrix r;
 VISVector n = q/q.norm();
 datatype x = n[0];
 datatype y = n[1];
 datatype z = n[2];
 datatype w = n[3];
 datatype xx = 2*x*x;
 datatype yy = 2*y*y;
 datatype zz = 2*z*z;
 //GeomReal ww = 2*w*w;
 datatype xy = 2*x*y;
 datatype xz = 2*x*z;
 datatype xw = 2*x*w;
 datatype yz = 2*y*z;
 datatype yw = 2*y*w;
 datatype zw = 2*z*w;
 
 r = i(3);
 r.poke(0,0) = 1.0 - yy - zz;
 r.poke(0,1) = xy + zw;
 r.poke(0,2) = xz - yw;
 r.poke(1,0) = xy - zw;
 r.poke(1,1) = 1.0 - xx - zz;
 r.poke(1,2) = yz + xw;
 r.poke(2,0) = xz + yw;
 r.poke(2,1) = yz - xw;
 r.poke(2,2) = 1.0 - xx - yy;
 return r;
 }*/


