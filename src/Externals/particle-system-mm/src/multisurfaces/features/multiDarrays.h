#ifndef __MULTI_D_ARRAYS_H__
#define __MULTI_D_ARRAYS_H__

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
  { *this = rhs; };
  ~array2D() 
  { if (_array) delete [] _array; };

  array2D& operator = (const array2D &rhs);
  array2D& operator = (float rhs);

  void resize(int idim, int jdim)
  { if ( _array ) delete _array; 
    _idim=idim; _jdim=jdim; 
    _ijdim=idim*jdim; 
    _array = new T[_ijdim]; }; 
  void size(int &idim, int &jdim)
  { idim = _idim; jdim = _jdim; };

  //--------------------------------------------------------------------
  // array indexing
  inline T& operator () (int i, int j) 
  { return _array[i+j*_idim]; };
  inline const T& operator () (int i, int j) const 
  { return _array[i+j*_idim]; };

  float  *_array;
  int _idim, _jdim, _ijdim;
};

//P.Petrov 2014
//changed to multidim arary to avoid index computation every time
//also why using float hard tyed when defined as templated class by type T?
template <class T> 
class array3D 
{
public:
  //--------------------------------------------------------------------
  // constructors and assignment
  array3D() 
  { 
      _idim=_jdim=_kdim=_ijdim=_ijkdim=0; 
      _array = 0; 
  };
    
  array3D(int idim, int jdim, int kdim)
  { 
      _array=0; resize( idim, jdim, kdim ); 
  };

  ~array3D() 
  { 
      if ( _array )
      {
           // De-Allocate memory to prevent memory leak
          for (int i = 0; i < _idim; ++i) {
            for (int j = 0; j < _jdim; ++j)
              delete [] _array[i][j];
        
            delete [] _array[i];
          }
          delete [] _array;
      }
  };

  array3D& operator = (const array3D &rhs);
  array3D& operator = (float rhs);

  void resize(int idim, int jdim, int kdim)
  { 
      if ( _array )
      {
           // De-Allocate memory to prevent memory leak
          for (int i = 0; i < idim; ++i) {
            for (int j = 0; j < jdim; ++j)
              delete [] _array[i][j];
        
            delete [] _array[i];
          }
          delete [] _array;
      }
      
    _idim=idim; _jdim=jdim; _kdim=kdim; 
    _ijdim=idim*jdim; _ijkdim=_ijdim*kdim;
    
    //_array = new T[_ijkdim];
        
      // Allocate memory
      _array = new T**[idim];
      for (int i = 0; i < idim; ++i) {
        _array[i] = new T*[jdim];   
        for (int j = 0; j < jdim; ++j)
          _array[i][j] = new T[kdim];
      }
  }; 
  
  void size(int &idim, int &jdim, int &kdim)
  { idim = _idim; jdim = _jdim; kdim = _kdim; };

  //--------------------------------------------------------------------
  // array indexing
  inline T& operator () (int i, int j, int k) 
  { return _array[i][j][k]; };
  inline const T& operator () (int i, int j, int k) const 
  { return _array[i][j][k]; };

  T  ***_array;
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
    _array = new T[_ijkldim]; }
 
  void size(int &idim, int &jdim, int &kdim, int &ldim) { 
    idim = _idim; jdim = _jdim; kdim = _kdim; ldim = _ldim;
  }

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
