//#include <cmath>
#include <cstdlib>
#include <iostream>

#ifdef _WIN32
//#pragma warning( disable: 4244 )
#endif


//-------------------------------------------------------------------------
// function   : = operators
// description: 
//-------------------------------------------------------------------------
template <class T>
array3D<T>& array3D<T>::operator = ( const array3D &rhs )
{
  if ( !(_idim == rhs._idim) || 
       !(_jdim == rhs._jdim) || 
       !(_kdim == rhs._kdim) )
    resize( rhs._idim, rhs._jdim, rhs._kdim );
    
  for ( int i = 0; i < _ijkdim; i++ )
    _array[i] = rhs._array[i];
    
  return *this;
}

template <class T>
array3D<T>& array3D<T>::operator = ( float rhs )
{
  for ( int i = 0; i < _ijkdim; i++ )
    _array[i] = rhs;
    
  return *this;
}


/************************************************************************/

template <class T>
array4D<T>& array4D<T>::operator = ( const array4D &rhs )
{
  if ( !(_idim == rhs._idim) || 
       !(_jdim == rhs._jdim) || 
       !(_kdim == rhs._kdim) ||
       !(_ldim == rhs._ldim) )
    resize( rhs._idim, rhs._jdim, rhs._kdim, rhs._ldim );
    
  for ( int i = 0; i < _ijkldim; i++ )
    _array[i] = rhs._array[i];
    
  return *this;
}

template <class T>
array4D<T>& array4D<T>::operator = ( float rhs )
{
  for ( int i = 0; i < _ijkldim; i++ )
    _array[i] = rhs;
    
  return *this;
}


