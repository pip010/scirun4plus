#ifndef __MULTI_D_ARRAYS_H__
#define __MULTI_D_ARRAYS_H__

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
  { *this = rhs; };
  ~array3D() 
  { if (_array) delete [] _array; };

  array3D& operator = (const array3D &rhs);
  array3D& operator = (float rhs);

  void resize(int idim, int jdim, int kdim)
  { if ( _array ) delete _array; 
    _idim=idim; _jdim=jdim; _kdim=kdim; 
    _ijdim=idim*jdim; _ijkdim=_ijdim*kdim;
    _array = new T[_ijkdim]; }; 
  void size(int &idim, int &jdim, int &kdim)
  { idim = _idim; jdim = _jdim; kdim = _kdim; };

  //--------------------------------------------------------------------
  // array indexing
  inline T& operator () (int i, int j, int k) 
  { return _array[i+j*_idim+k*_ijdim]; };
  inline const T& operator () (int i, int j, int k) const 
  { return _array[i+j*_idim+k*_ijdim]; };

  float  *_array;
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
  { *this = rhs; };
  ~array4D() 
  { if (_array) delete [] _array; };

  array4D& operator = (const array4D &rhs);
  array4D& operator = (float rhs);

  void resize(int idim, int jdim, int kdim, int ldim)
  { if ( _array ) delete _array; 
    _idim=idim; _jdim=jdim; _kdim=kdim; _ldim=ldim; 
    _ijdim=idim*jdim; _ijkdim=_ijdim*kdim; _ijkldim=_ijkdim*ldim;
    _array = new T[_ijkldim]; }; 
  void size(int &idim, int &jdim, int &kdim, int &ldim)
  { idim = _idim; jdim = _jdim; kdim = _kdim; ldim = _ldim; };

  //--------------------------------------------------------------------
  // array indexing
  inline T& operator () (int i, int j, int k, int l) 
  { return _array[i+j*_idim+k*_ijdim+l*_ijkdim]; };
  inline const T& operator () (int i, int j, int k, int l) const 
  { return _array[i+j*_idim+k*_ijdim+l*_ijkdim]; };

  float  *_array;
  int _idim, _jdim, _kdim, _ldim, _ijdim, _ijkdim, _ijkldim;
};

#include "multiDarrays.T"

#endif // __MULTI_D_ARRAYS_H__
