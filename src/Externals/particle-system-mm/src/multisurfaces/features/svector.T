#include <cstdlib>
#include <iostream>

using namespace custom_class;
using namespace std;

//-------------------------------------------------------------------------
// function   : front()
// description: returns the first element of the vector
//-------------------------------------------------------------------------
template <class T> T& svector<T>::front() 
{
  if ( _vector.empty() )
  {
    cout << "SVECTOR IS EMPTY: cannot return front()\n";
    exit( -1 );
  }

  return _vector.front();
}

//-------------------------------------------------------------------------
// function   : back()
// description: returns the last element of the vector
//-------------------------------------------------------------------------
template <class T> T& svector<T>::back() 
{
  if ( _vector.empty() )
  {
    cout << "SVECTOR IS EMPTY: cannot return back()\n";
    exit( -1 );
  }

  return _vector.back();
}

//-------------------------------------------------------------------------
// function   : pop_back()
// description: erase the last element in the vector
//-------------------------------------------------------------------------
template <class T> void svector<T>::pop_back()
{
  if ( _vector.empty() )
  {
    cout << "SVECTOR IS EMPTY: cannot pop_back()\n";
    exit( -1 );
  }

  _vector.pop_back();
}

//-------------------------------------------------------------------------
// function   : push_back()
// description: adds an element to the end of the vector
//-------------------------------------------------------------------------
template <class T> void svector<T>::push_back( const T& x )
{
  _vector.push_back( x );
}

//-------------------------------------------------------------------------
// function   : operator []
// description: returns a ref to the vector element at position i
//-------------------------------------------------------------------------
template <class T> T& svector<T>::operator [] ( unsigned int i )
{
  if ( i >= _vector.size() )
  {
    cout << "SVECTOR OUT OF BOUNDS ERROR: " << i << " >= " <<
      _vector.size() << endl;
    abort();
    exit( -1 );
  }

  return _vector[i];
}

 template <class T> const T& svector<T>::operator [] ( unsigned int i )
  const
{
  if ( i >= _vector.size() )
  {
    cout << "SVECTOR OUT OF BOUNDS ERROR: " << i << " >= " <<
      _vector.size() << endl;
    exit( -1 );
  }

  return _vector[i];
}

//-------------------------------------------------------------------------
// function   : at()
// description: returns a ref to the vector element at position i
//-------------------------------------------------------------------------
template <class T> T& svector<T>::at( unsigned int i )
{
  if ( i >= _vector.size() )
  {
    cout << "SVECTOR OUT OF BOUNDS ERROR: " << i << " >= " <<
      _vector.size() << endl;
    exit( -1 );
  }

  return _vector.at( i );
}

template <class T> const T& svector<T>::at( unsigned int i ) const
{
  if ( i >= _vector.size() )
  {
    cout << "SVECTOR OUT OF BOUNDS ERROR: " << i << " >= " <<
      _vector.size() << endl;
    exit( -1 );
  }

  return _vector.at( i );
}

//-------------------------------------------------------------------------
// function   : operator =
// description: erases the elements from this and copies that
//-------------------------------------------------------------------------
template <class T> svector<T>& svector<T>::operator = ( const
                                                               svector<T>&
                                                               that )
{
  _vector = that._vector;
  return *this;
}


