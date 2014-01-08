// ******************************************************************
// VISPack. Copyright (c) 1994-2000 Ross Whitaker rtw@utk.edu       *
// For conditions of distribution and use, see the file LICENSE.txt *
// accompanying this distribution.                                  *
// ******************************************************************

// $Id: matrix.h,v 1.4 2003/02/25 22:37:10 whitaker Exp $

// File:           matrix.h
// Author:         Samuel G. Burgiss & Zhong Tang 
// Institution:    The University of Tennessee, Knoxville
// Contents:       A VISMatrix class, a VISVector class
//                 that contains functions specifically for computer vision
//                 applications.
// Log of changes: Sept 7, 1999 -- Call LAPACK functions instead: LU, SVD, QR

#ifndef __matrix
#define __matrix

#include <stdlib.h>
#include <iostream>
#include <math.h>

#ifdef HAVE_UNISTD_H
#  include <unistd.h>
#endif

#include "matrixExports.h"

//#define VISMATRIXDOUBLE

//extern int MATRIX_TOTAL_SIZE;

// Through this new matrix class you can enjoy a MATLAB
// like interface and functionality in the comfort of your
// very own C++ programs.  You can add, subtract, multiply,
// take the inverse, etc.  All through a very simple and
// intuitive interface.  All the functions you need to start
// using the powers of linearalgebra are listed (and labled!)
// below in the VISMatrix class.
//
// The biggest difference is this class and MATLAB is that
// this class uses vectors.
//
// Be sure and check out the non-member functions listed below
// the class definition for VISMatrix

// IMPORTANT  fix all this SH@#)($*@

//#define VISMATRIXDOUBLE
#ifdef VISMATRIXDOUBLE
typedef double datatype;
#else
typedef float datatype;
#endif

class matrix_SHARE VISMatrixData{
    friend class VISMatrix;
    friend class VISVector;
  protected:
    datatype* _data;
    int _refcnt;
    int _size;
  public:
    
    VISMatrixData() : _data(NULL), _refcnt(0), _size(0) {}
    VISMatrixData(int size) : _data(new datatype[size]), _refcnt(0), _size(size) {
//        MATRIX_TOTAL_SIZE += size;
}
    VISMatrixData(const datatype* data,int size);
    
    ~VISMatrixData()
	{
	    delete [] _data;
//	    MATRIX_TOTAL_SIZE -= _size;
	}
};

class VISVector;


class matrix_SHARE VISMatrix{

    friend class VISVector;

  protected:
    int _cols;
    int _rows;
    VISMatrixData* _matrixdata;

    void initialize(int r,int c);
    //Below are functions used to manipulate the VISMatrixData object
    //from within the VISMatrix object.
    VISMatrixData& matrixDataRef(){return(*(_matrixdata));}
    const VISMatrixData* matrixDataPtr() const {return((_matrixdata));} 
    datatype& dataRef(){return dataRef(0);}
    datatype& dataRef(int i){return(*((_matrixdata->_data)+i));}
    void refCntInc(){_matrixdata->_refcnt++;}
    void refCntDec(){_matrixdata->_refcnt--;}

////Functions to convert to and from Numerical Recipies(NR) matrix.
    datatype** createNR() const;
    void becomeNR(datatype** nr, int r, int c);

  
  public:

////Functions to convert to and from LAPACK matrix data
    datatype* createLAPK() const;
    void becomeLAPK(datatype*, int, int);

    const datatype* dataBuf() const {return(_matrixdata->_data);}
    datatype* dataBufRef() {return(_matrixdata->_data);}

    // for debugging
    int refCnt() const{return _matrixdata->_refcnt;}


    //*****************************//
    //  FUNCTIONS USED BY THE USER //
    //*****************************//
//****Constructors and Destructors, the usual stuff
    VISMatrix();
    VISMatrix(const VISMatrix& other) : _cols(other._cols),_rows(other._rows),
      _matrixdata(other._matrixdata){refCntInc();}
    VISMatrix(const VISVector& other);

    VISMatrix(int row, int col);
    VISMatrix(char* filename);//initialize a matrix from a filename
    
    ~VISMatrix();
//****Finding the size of the matrix
    int r() const { return(_rows); }//return # of rows
    int c() const { return(_cols); }//return # of cols
    int isValid() const
	{return ((_rows > 0)&&(_cols > 0)); }
 // returns true of matrices are of same size
  int compareSize(const VISMatrix &other) const
  {return ((other._cols == _cols)&&(other._rows == _rows));}
//****Interacting with the data of the matrix
    //(if you liked BASIC you'll love peek and poke!)
    //The great thing with this set of functions is that your rhs
    //can use the shorter notation of M[1][2]+M[2][3]+ ...
    //while the lhs has the less convenient M.poke(2,2)
    //examples:    VISMatrix m(5,5),n(5,5);
    //             m.poke(2,2) = n[3][3];
    //             m.poke(2,2) = n.peek(3,3);
    const datatype* operator[](int r) const;// bounds checking only on rows
    datatype operator()(int r,int c) const{// no bounds checking!
      return (data(r + c*_rows));}
    datatype peek(int r,int c) const{// no bounds checking!
      return (data(r + c*_rows));}
    datatype& poke(int r,int c);
//****Operators that have been overloaded
  VISMatrix& operator=(const VISMatrix& other);
  VISMatrix& operator=(const VISVector& other);
  VISMatrix& operator=(datatype s);
  VISMatrix operator+(const VISMatrix& other) const;
  VISMatrix operator+(datatype s) const ;
  VISMatrix operator-(const VISMatrix& other) const;
  VISMatrix operator-(datatype s) const;
  VISMatrix operator-() const;
  VISMatrix operator*(const VISMatrix& other) const;
  VISVector operator*(const VISVector& other) const;
  VISMatrix operator*(datatype s) const;
  VISMatrix operator/(datatype s) const;
  VISMatrix& operator+=(const VISMatrix& other)
  {
	return (operator=(operator+(other)));
  }
  VISMatrix& operator-=(const VISMatrix& other)
  {
	return (operator=(operator-(other)));
  }
  VISMatrix& operator*=(const VISMatrix& other)
  {
	return (operator=(operator*(other)));
  }
  VISMatrix& operator+=(datatype s)
  {
	return (operator=(operator+(s)));
  }
  VISMatrix& operator-=(datatype s)
  {
	return (operator=(operator-(s)));
  }
  VISMatrix& operator*=(datatype s)
  {
	return (operator=(operator*(s)));
  }
  VISMatrix& operator/=(datatype s)
  {
	return (operator=(operator/(s)));
  }
  VISMatrix operator==(datatype s) const;//returns a matrix of 1/0's for equal
  //to s and not equal to s respectively
  VISMatrix operator==(const VISMatrix& other) const;//returns a matrix of 1/0's
  //for equal to and not equal elements of this and other


//****General functions
  VISMatrix t() const{return (~(*this));} //transpose
  datatype max(int &r, int &c) const;//return max value in matrix
  datatype max() const{int r,c;
  return (max(r,c));}//return max value in matrix
  datatype min(int &r, int &c) const;//return min value in matrix
  datatype min() const{int r,c;
  return (min(r,c));}//return min value in matrix
  VISVector mean() const;//returns mean of each row in a matrix
  VISVector meanOfCols() const;//returns mean of each col in a matrix
  datatype normOfCol(int col = 0) const;
  datatype norm() const;//returns l2 norma of matrix
  datatype norm(datatype e) const;//returns l2 norma of matrix with epsilon sq added
  datatype normOfRow(int row = 0) const;
  datatype sum() const;

  //VISVector normOfRow() const;//returns norms of rows
  VISMatrix dot(const VISMatrix &other) const;
  VISMatrix inv() const{return inverseSVD();} //default inverse is svd(from nr)
  VISMatrix cov() const;//return the cov of matrix
  //***concatination
  //concatinate cols i.e. a=[1 2 3  b=[ 7  8  9
  //                         4 5 6]    10 11 12]
  // c=a.concat(b)   c=[1 2 3  7  8  9
  //                    4 5 6 10 11 12]
  //use concatRow if you want to concatinate rows
  //concatinate rows i.e. a=[1 2 3  b=[ 7  8  9
  //                         4 5 6]    10 11 12]
  // c=a.concatRows(b)   c=[1  2  3
  //                        4  5  6
  //                        7  8  9
  //                       10 11 12]
  VISMatrix concat(const VISMatrix& second);
  VISMatrix concat(const VISVector& second);
  VISMatrix concatRow(const VISMatrix& second);
  //Functions implemented by calling LAPACK
  datatype det() const;
  datatype trace() const;
  void svd(VISMatrix &U, VISMatrix &W, VISMatrix &V) const;
  void QR(VISMatrix &Q, VISMatrix &R) const;
  VISMatrix inverseSVD() const;
  VISMatrix inverseGJ() const;

  //miscellaneous functions
  VISVector vec(int col = 0) const;//makes a vector from one col of a matrix
  //VISVector vecfromcol(int col) const{return (vec(col));}//same as vec
  VISVector vecFromRow(int row = 0) const; //returns a column vector
  VISVector vecFromDiag() const;//returns a vector from the diagonal
  VISMatrix col(int col) const; //returns a col matrix
  VISMatrix row(int row) const; //returns a row matrix
  //returns a submatrix of a matrix
  VISMatrix peekROI(int startrow,int endrow,int startcol,int endcol) const;
  //returns a matrix containing a submatrix; row,col are upper left corner
  VISMatrix&  pokeROI(int row,int col, const VISMatrix& roi);
  VISMatrix&  pokeROI(int row,int col, const VISVector& roi);
  //VISMatrix append_col(VISMatrix col_matrix, int col);
  //VISMatrix append_row(VISMatrix row_matrix, int row);
  //********************************************//
  //  FUNCTIONS USUALLY USED BY OTHER FUNCTIONS //
  //********************************************//
  //Functions used to get data from underlying VISMatrixData
  //object.  (These functions are used mostly by the
  //VISMatrix object but can be used by the user if he/she
  //needs a pointer to the data block.
  const datatype* dataPtr(int i = 0) const{
	return((_matrixdata->_data)+i);}
  datatype data(int i = 0) const{
	return(*((_matrixdata->_data)+i));}
  //Functions for fstream stuff
  void print(std::ostream *os) const;
  void read(std::istream *os);
  //Left in for backward compatibility
  VISMatrix operator~() const; //transpose, useage: ~M
  VISMatrix operator!() const{ 
	if(_cols!=_rows){
	  std::cout << "VISMatrix operator!(): matrix must be square" << std::endl;
	  exit(0);}
	return inverseSVD();} //inverse useage: !M
  
};

inline VISMatrix operator*(datatype s,const VISMatrix m){
  return (m.operator*(s));}
inline VISMatrix operator+(datatype s,const VISMatrix m){
  return (m.operator+(s));}
inline VISMatrix operator-(datatype s,const VISMatrix m){
  return (m.operator-(s));}
inline VISMatrix operator/(datatype s,const VISMatrix m){
  return (m.operator/(s));}


void printNR(datatype** nr, int r, int c);
void destroyNR(datatype** nr, int r, int c);
//corresponding functions with LAPACK
void printLAPK(datatype* lapk, int r, int c);
void destroyLAPK(datatype* lapk);

//solves a system of linear equation: Ax=b
matrix_SHARE VISVector solvLinearEqn(const VISMatrix& A,const VISVector& b);
matrix_SHARE VISMatrix VISIdentity(int size); //return an identity matrix 
matrix_SHARE VISMatrix homogeneousTransform(VISMatrix &mother, VISVector &vother);

inline std::istream& operator>> (std::istream &is, VISMatrix& m)
{
  m.read(&is);
  return is;
}
inline std::ostream& operator<<(std::ostream& os, const VISMatrix& m)
{
    m.print(&os);
    return os;
}

//VISVector Class
class matrix_SHARE VISVector : public VISMatrix
{

  friend class VISMatrix;

protected:

  //Replaced by ZT for use of Lapack
  ////Functions to convert to and from NR data
  datatype* createNR() const;
  void becomeNR(datatype* nr, int n);


public:

  ////Functions to convert to and from LAPACK data
  datatype* createLAPK() const;
  void becomeLAPK(datatype* lapk, int n);

  int refCnt() const{return VISMatrix::refCnt();}
    
  //Constructors and Destructors


  VISVector() : VISMatrix(){}
  VISVector(const VISVector& other) : VISMatrix(other){}
  VISVector(const VISMatrix& other);
  VISVector(char* filename);
  VISVector(int r) : VISMatrix(r,1){}
  VISVector(datatype x,datatype y) : VISMatrix(2,1){
    poke(0)=x; poke(1)=y;}//2D point
  VISVector(datatype x,datatype y,datatype z) : VISMatrix(3,1){
    poke(0)=x; poke(1)=y; poke(2)=z;}//3D point
  VISVector(datatype x,datatype y,datatype z,datatype h) : VISMatrix(4,1){
    poke(0)=x; poke(1)=y; poke(2)=z; poke(3)=h;}//homogeneous 3D pt
  ~VISVector(){};
  //Interacting with the data of the vector
  datatype peek(int r) const{return VISMatrix::peek(r,0);}
  datatype operator()(int r) const{return VISMatrix::peek(r,0);}
  datatype operator[](int r) const{return VISMatrix::peek(r,0);}
  datatype& poke(int r){return VISMatrix::poke(r,0);}
  VISVector peekROI(int startelement,int endelement) const{
    VISVector v;
    v=VISMatrix::peekROI(startelement,endelement,0,0);
    return v;
  }

  VISVector& pokeROI(int row, const VISVector& roi);

  //Functions from matrix class
  //Finding the size of the matrix
  int n() const{return r();}
  int isValid() const {return((_rows > 0));}
  //Operators that have been overloaded
  VISVector& operator=(const VISMatrix& other);
  VISVector& operator=(const VISVector& other){
    VISMatrix::operator=(other); return *this;}
  VISVector& operator=(datatype s){
    VISMatrix::operator=(s); return *this;}
  VISVector operator+(const VISVector& other) const{
    VISVector v; v=VISMatrix::operator+(other);
    return (v);}
  VISVector operator-(const VISVector& other) const{
    VISVector v; v=VISMatrix::operator-(other);
    return (v);}
  VISVector operator+(datatype s) const{
    VISVector v; v=VISMatrix::operator+(s);
    return (v);}    
  VISVector operator-(datatype s) const{
    VISVector v; v=VISMatrix::operator-(s);
    return (v);}
  VISVector operator*(datatype s) const{
    VISVector v; v=VISMatrix::operator*(s);
    return (v);}    
  VISVector operator/(datatype s) const{
    VISVector v; v=VISMatrix::operator/(s);
    return (v);}    
  VISVector operator-() const{
    VISVector v; v=VISMatrix::operator-();
    return (v);}
  VISVector& operator+=(const VISVector& other){
    return (operator=(operator+(other)));}
  VISVector& operator-=(const VISVector& other){
    return (operator=(operator-(other)));}
  VISVector& operator+=(datatype s){
    return (operator=(operator+(s)));}
  VISVector& operator-=(datatype s){
    return (operator=(operator-(s)));}
  VISVector& operator*=(datatype s){
    return (operator=(operator*(s)));}
  VISVector& operator/=(datatype s){
    return (operator=(operator/(s)));}
  VISVector operator==(datatype s) const{//returns a vector of 1/0's for equal
    VISVector v; v = VISMatrix::operator==(s);
    return (v);}//to s and not equal to s respectively
  VISVector operator==(const VISVector& other) const{
    VISVector v; v = VISMatrix::operator==(other);
    return (v);}
  VISMatrix operator~() const{
    VISMatrix m; m = rowMatrix();
    return (m);}//returns a row vector in the form of a matrix
  VISMatrix t() const{
    VISMatrix m; m=rowMatrix();
    return (m);}


  VISVector augment() const {
    VISVector v(_rows + 1); 
    for (int i = 0; i < _rows; i++)
      *(((v._matrixdata)->_data)+i) = *(((_matrixdata)->_data)+i);
    *(((v._matrixdata)->_data)+_rows) = 1.0f;
    return (v);
  }




  datatype max(int &n) const{int c;//return max value in vector
  return (VISMatrix::max(n,c));}//if there is a tie n will point to ONE of them
  datatype max() const{int r,c;
  return (VISMatrix::max(r,c));}//return max value in vector
  datatype norm() const {return (VISMatrix::norm());}
  datatype norm(datatype e) const {return (VISMatrix::norm(e));}
  //returns l2 norma of matrix with epsilon sq added    
  datatype min() const {int r,c; return (VISMatrix::min(r,c));}
  datatype mean() const;
  //Functions for vectors
  VISVector cross(const VISVector& other) const;
  datatype dot(const VISVector& other) const;

  // returns the exterior product of two vectors
  VISMatrix exterior(const VISVector &v) const;

  VISMatrix mat(int cols = 1) const;//create a colmatrix from the vector
  VISMatrix rowMatrix(int rows = 1) const;
  //***concatination
  //the default(concat) concatination is to concatinate columns
  //of second into columns of first
  //concatinate cols i.e. a=[1 2 3  b=[ 7  8  9
  //                         4 5 6]    10 11 12]
  // c=concat(a,b)   c=[1 2 3  7  8  9
  //                    4 5 6 10 11 12]
  //use concatRow if you want to concatinate rows
  //concatinate rows i.e. a=[1 2 3  b=[ 7  8  9
  //                         4 5 6]    10 11 12]
  // c=concatRows(a,b)   c=[1  2  3
  //                        4  5  6
  //                        7  8  9
  //                       10 11 12]
  //VISMatrix VISVector::concat(const VISMatrix& second);
  //VISMatrix VISVector::concat(const VISVector& second);
  //VISVector VISVector::concatRow(const VISVector& second);
  VISMatrix concat(const VISMatrix& second);
  VISMatrix concat(const VISVector& second);
  VISVector concatRow(const VISVector& second);

  //Functions for fstream stuff
  void print(std::ostream *os) const;
  void read(std::istream *os);
};

VISMatrix addToCols(const VISMatrix &mat, const VISVector &vec);

inline VISVector operator*(datatype s,const VISVector v){ return (v.operator*(s)); }
inline VISVector operator+(datatype s,const VISVector v){ return (v.operator+(s)); }
inline VISVector operator-(datatype s,const VISVector v){ return (v.operator-(s)); }
inline VISVector operator/(datatype s,const VISVector v){ return (v.operator/(s)); }

inline VISVector makept(datatype x,datatype y,datatype z)
{
  VISVector ans(x,y,z);
  return ans;
}
inline std::istream& operator>> (std::istream &is, VISVector& v)
{
  v.read(&is);
  return is;
}
inline std::ostream& operator<<(std::ostream& os, const VISVector& v)
{
  v.print(&os);
  return os;
}



#endif





