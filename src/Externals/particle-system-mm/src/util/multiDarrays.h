#ifndef __MULTI_D_ARRAYS_H__
#define __MULTI_D_ARRAYS_H__

#include <cstdlib>
#include <iostream>

template <class T> 
class array1D 
{
public:
  //--------------------------------------------------------------------
  // constructors and assignment
  array1D() 
  { _idim=0; _array = NULL; };
  array1D(int idim)
  { _array=NULL; resize( idim ); };
  array1D(const array1D &rhs) 
  { _array=NULL; *this = rhs; };
  ~array1D() 
  { if (_array) delete [] _array; };

  array1D& operator = (const array1D &rhs);
  array1D& operator = (float rhs);

  void resize(int idim)
  { if ( _array ) delete [] _array; 
    _idim=idim;
    _array = new T[_idim]; }; 
  void size(int &idim) const
  { idim = _idim; };

  //--------------------------------------------------------------------
  // array indexing
  inline T& operator () (int i) 
  { if ((i>=_idim)||(i<0)) std::cout<<"array1D: overbounds"<<std::endl;
    return _array[i]; };
  inline const T& operator () (int i) const 
  { if ((i>=_idim)||(i<0)) std::cout<<"array1D: overbounds"<<std::endl;
    return _array[i]; };

  inline T& operator [] (int i) 
  { if ((i>=_idim)||(i<0)) std::cout<<"array1D: overbounds"<<std::endl;
    return _array[i]; };
  inline const T& operator [] (int i) const 
  { if ((i>=_idim)||(i<0)) std::cout<<"array1D: overbounds"<<std::endl;
    return _array[i]; };

  inline bool checkBounds(int i) const
  { if ((i>=_idim)||(i<0)) return false;
    return true; };

  T  *_array;
  int _idim;
};

template <class T> 
class array2D 
{
public:
  //--------------------------------------------------------------------
  // constructors and assignment
  array2D() 
  { _idim=_jdim=_ijdim=0; _array = NULL; };
  array2D(int idim, int jdim)
  { _array=NULL; resize( idim, jdim ); };
  array2D(const array2D &rhs) 
  { _array=NULL; *this = rhs; };
  ~array2D() 
  { if (_array) delete [] _array; };

  array2D& operator = (const array2D &rhs);
  array2D& operator = (float rhs);

  void resize(int idim, int jdim)
  { if ( _array ) delete [] _array; 
    _idim=idim; _jdim=jdim; 
    _ijdim=idim*jdim; 
    _array = new T[_ijdim]; }; 
  void size(int &idim, int &jdim) const
  { idim = _idim; jdim = _jdim; };

  //--------------------------------------------------------------------
  // array indexing
  inline T& operator () (int i, int j) 
  { if ((i>=_idim)||(j>=_jdim)||(i<0)||(j<0))
      std::cout<<"array2D: overbounds"<<std::endl;
    return _array[i+j*_idim]; };
  inline const T& operator () (int i, int j) const 
  { if ((i>=_idim)||(j>=_jdim)||(i<0)||(j<0))
      std::cout<<"array2D: overbounds"<<std::endl;
    return _array[i+j*_idim]; };

  inline bool checkBounds(int i, int j) const
  { if ((i>=_idim)||(j>=_jdim)||(i<0)||(j<0)) return false;
    return true; };

  T  *_array;
  int _idim, _jdim, _ijdim;
};

template <class T> 
class array3D 
{
public:
  //--------------------------------------------------------------------
  // constructors and assignment
  array3D() 
  { _idim=_jdim=_kdim=_ijdim=_ijkdim=0; _array = NULL; };
  array3D(int idim, int jdim, int kdim)
  { _array=NULL; resize( idim, jdim, kdim ); };
  array3D(const array3D &rhs) 
  { _array=NULL; *this = rhs; };
  ~array3D() 
  { if (_array) delete [] _array; };

  array3D& operator = (const array3D &rhs);
  array3D& operator = (float rhs);

  void resize(int idim, int jdim, int kdim)
  { if ( _array ) delete [] _array; 
    _idim=idim; _jdim=jdim; _kdim=kdim; 
    _ijdim=idim*jdim; _ijkdim=_ijdim*kdim;
    _array = new T[_ijkdim]; }; 
  void size(int &idim, int &jdim, int &kdim) const
  { idim = _idim; jdim = _jdim; kdim = _kdim; };

  //--------------------------------------------------------------------
  // array indexing
  inline T& operator () (int i, int j, int k) 
  { if ((i>=_idim)||(j>=_jdim)||(k>=_kdim)||(i<0)||(j<0)||(k<0))
      std::cout<<"array3D: overbounds"<<std::endl;
    return _array[i+j*_idim+k*_ijdim]; };
  inline const T& operator () (int i, int j, int k) const 
  { if ((i>=_idim)||(j>=_jdim)||(k>=_kdim)||(i<0)||(j<0)||(k<0))
      std::cout<<"array3D: overbounds"<<std::endl;
    return _array[i+j*_idim+k*_ijdim]; };

  inline bool checkBounds(int i, int j, int k) const
  { if ((i>=_idim)||(j>=_jdim)||(k>=_kdim)||(i<0)||(j<0)||(k<0))
      return false;
    return true; };

  T  *_array;
  int _idim, _jdim, _kdim, _ijdim, _ijkdim;
};


template <class T> 
class array4D 
{
public:
  //--------------------------------------------------------------------
  // constructors and assignment
  array4D() 
  { _idim=_jdim=_kdim=_ldim=_ijdim=_ijkdim=_ijkldim=0; _array = NULL; };
  array4D(int idim, int jdim, int kdim, int ldim)
  { _array=NULL; resize( idim, jdim, kdim, ldim ); };
  array4D(const array4D &rhs) 
  { _array=NULL; *this = rhs; };
  ~array4D() 
  { if (_array) delete [] _array; };

  array4D& operator = (const array4D &rhs);
  array4D& operator = (float rhs);

  void resize(int idim, int jdim, int kdim, int ldim)
  { if ( _array ) delete [] _array; 
    _idim=idim; _jdim=jdim; _kdim=kdim; _ldim=ldim; 
    _ijdim=idim*jdim; _ijkdim=_ijdim*kdim; _ijkldim=_ijkdim*ldim;
    _array = new T[_ijkldim]; }; 
  void size(int &idim, int &jdim, int &kdim, int &ldim) const
  { idim = _idim; jdim = _jdim; kdim = _kdim; ldim = _ldim; };

  //--------------------------------------------------------------------
  // array indexing
  inline T& operator () (int i, int j, int k, int l) 
  { if ((i>=_idim)||(j>=_jdim)||(k>=_kdim)||(l>=_ldim)||
        (i<0)||(j<0)||(k<0)||(l<0))
      std::cout<<"array4D: overbounds"<<std::endl;
    return _array[i+j*_idim+k*_ijdim+l*_ijkdim]; };
  inline const T& operator () (int i, int j, int k, int l) const 
  { if ((i>=_idim)||(j>=_jdim)||(k>=_kdim)||(l>=_ldim)||
        (i<0)||(j<0)||(k<0)||(l<0))
      std::cout<<"array4D: overbounds"<<std::endl;
    return _array[i+j*_idim+k*_ijdim+l*_ijkdim]; };

  inline bool checkBounds(int i, int j, int k, int l) const
  { if ((i>=_idim)||(j>=_jdim)||(k>=_kdim)||(l>=_ldim)||
        (i<0)||(j<0)||(k<0)||(l<0)) return false;
    return true; };

  T  *_array;
  int _idim, _jdim, _kdim, _ldim, _ijdim, _ijkdim, _ijkldim;
};

#include <multiDarrays.T>

#endif // __MULTI_D_ARRAYS_H__
