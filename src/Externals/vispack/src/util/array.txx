// ******************************************************************
// VISPack. Copyright (c) 1994-2000 Ross Whitaker rtw@utk.edu       *
// For conditions of distribution and use, see the file LICENSE.txt *
// accompanying this distribution.                                  *
// ******************************************************************

// $Id: array.txx,v 1.1.1.1 2003/02/12 16:51:54 whitaker Exp $

template<class T>
 void VISArray<T>::copyItem(T& dst,const T& src) {
    dst = src;
}

// ---------------------------------------------------------------------------

template<class T>
 void VISArray<T>::sizeTo(unsigned new_size) {
    unsigned i;
    if (new_size > _size) {
	T* new_array = new T[new_size];
	for (i=0; i<_n; ++i) { new_array[i] = _array[i]; }
	unsigned old_size = _size;
	_size = new_size;
	
	delete[] _array;
	_array = new_array;
	for (i=old_size; i<_size; ++i) clearItem(i);
    }
    else if (new_size < _n) {
	for (i=new_size; i<_n; ++i)	    deleteItem(i);
	_n = new_size;
    }
}

template<class T>
 void VISArray<T>::grow() {
    if (_size == 0)
	sizeTo(8);
    else
	sizeTo(_size*2);
}

// This should be used instead of sizeTo -- its more efficient
// in terms of doubling sizes.  This doubles the size until
// it is greater than or equal to new_size.

template<class T>
 void VISArray<T>::growTo(unsigned new_size) {
    if (new_size > _size) {
	int sz = _size;
	if (sz == 0) sz++;		// Watch out for zero!
	while (sz < new_size)
	    sz *= 2;
	sizeTo(sz);
    }
}

// ---------------------------------------------------------------------------

template<class T>
 VISArray<T>::~VISArray() {
    delete[] _array;
}

template<class T>
 void VISArray<T>::clearItem(unsigned) {}

template<class T>
 void VISArray<T>::deleteItem(unsigned) {}

template<class T>
 unsigned VISArray<T>::size() const	{ return _size; }

template<class T>
 unsigned VISArray<T>::n() const	{ return _n; }

template<class T>
 const T& VISArray<T>::peek(unsigned idx) const {
    return refAt(idx);
}

template<class T>
 const T& VISArray<T>::itemAt(unsigned idx) const {
    return refAt(idx);
}


template<class T>
 T& VISArray<T>::refAt(unsigned idx) const {
#ifdef DEBUG
    if (idx >= _n) {
	GfxHandler::panic ("Range check failure",
			   "VISArray",
			   "refAt()");
    }
#endif
    return _array[idx];
}

template<class T>
 void VISArray<T>::clearItems(unsigned i0,unsigned i1) {
    for (unsigned i=i0; i<i1; ++i) clearItem(i);
}

template<class T>
 void VISArray<T>::clear() {
    for (unsigned i=0; i<_n; ++i) clearItem(i);
    _n = 0;
}

template<class T>
 void VISArray<T>::deleteItems() {
    for (unsigned i=0; i<_n; ++i) deleteItem(i);
    _n = 0;
}

template<class T>
 const T& VISArray<T>::operator[](unsigned idx) const {
    return peek(idx);
}

template<class T>
 T& VISArray<T>::at(unsigned idx) {
    if (_size == 0 || idx >= _size) growTo(idx + 1);
    if (idx >= _n) {
	clearItems(_n,idx+1);
	//for (unsigned i=_n; i<=idx; ++i) clearItem(i);  // Doesn't work!
	_n = idx + 1;
    }
    return _array[idx];
}

template<class T>
T& VISArray<T>::poke(unsigned idx) {
    if (_size == 0 || idx >= _size) growTo(idx + 1);
    if (idx >= _n) {
	clearItems(_n,idx+1);
	//for (unsigned i=_n; i<=idx; ++i) clearItem(i);  // Doesn't work!
	_n = idx + 1;
    }
    return _array[idx];
}
    
template<class T>
 void VISArray<T>::appendItem(const T& t) {
    insertItemAt(t,_n);
}

template<class T>
 void VISArray<T>::prependItem(const T& t)
{
    insertItemAt(t,0);
}

template<class T>
 void VISArray<T>::insertItemAt(const T& t,unsigned idx) {
    if (_size == 0 || idx >= _size) growTo(idx+1);
    if (_n >= _size) grow();
    if (idx < _n) _n++;           // Extended the range
    else {
	clearItems(_n,idx);
	_n = idx + 1;
    }
    if (_n > 0) for (unsigned i=_n-1; i>idx; --i) _array[i] = _array[i-1];
    copyItem(_array[idx],t);
}

template<class T>
 void VISArray<T>::replaceItemAtWith(unsigned idx, const T& t)
{
    if (_size == 0 || idx >= _size) growTo(idx+1);
    if (_n >= _size) grow();
    deleteItem(idx);
    copyItem(_array[idx],t);
    if (idx >= _n) {
	clearItems(_n,idx);
	_n = idx + 1;
    }
}

template<class T>
 void VISArray<T>::removeItemAt(unsigned idx) {
#ifdef DEBUG
    if (idx >= _n) {
	GfxHandler::panic("Range exceeded",
			  "GrowVISArray",
			  "removeItemAt()");
    }
    else
#endif
    {
	deleteItem(idx);
	for (unsigned i=idx; i<_n-1; ++i) _array[i] = _array[i+1];
	--_n;
    }
}

template<class T>
 void VISArray<T>::reverse() {
    unsigned int	stop = _n/2, top = _n - 1;
    for (unsigned int i=0; i<stop; ++i) {
	T	temp = _array[i];
	_array[i] = _array[top - i];
	_array[top - i] = temp;
    }
}

template<class T>
void VISArray<T>::sort() 
{
    int i;
//    unsigned int stop = _n;
    boolean changed = TRUE;
    T tmp; 
    while (changed)
	{
	    changed = FALSE;
	    for (i=0; i < (int)_n - 1; i++) 
		{
		    if ((_array[i]) < (_array[i + 1]))
			{
			    changed = TRUE;
			    tmp = _array[i];
			    _array[i] = _array[i + 1];
			    _array[i + 1] = tmp;
			}
		}
	}
}


template<class T>
 void VISArray<T>::reverse(unsigned start,unsigned end) {
    if (end <= start)	// Fix end to end of array
	end = _n;
    if (start >= _n)	// Cant do anything here
	return;
    unsigned int	stop = (start+end)/2, top = end + start - 1;
    for (unsigned int i=start; i<stop; ++i) {
	T	temp = _array[i];
	_array[i] = _array[top - i];
	_array[top - i] = temp;
    }
}

template<class T>
 void VISArray<T>::rowReverse(unsigned rows,unsigned cols) {
    unsigned int	rstop = rows/2, top = rows - 1;
    for (unsigned i=0; i<rstop; ++i) {
	for (unsigned j=0; j<cols; ++j) {
	    unsigned idx1 = i*cols+j;
	    unsigned idx2 = (top-i)*cols+j;
	    
	    T temp = _array[idx1];
	    _array[idx1] = _array[idx2];
	    _array[idx2] = temp;
	}
    }
}

// Assignment operators
template<class T>
 VISArray<T>& VISArray<T>::operator=(const VISArray<T>& s) {
    if (&s == this) return *this;	// Check for assignment to self
    deleteItems();		// Clear out the array - delete what is there
    growTo(s._size);		// This wont copy cuz n=0
    for (unsigned i=0; i<s.n(); ++i) copyItem(_array[i],s.peek(i));
    _n = s.n();			// Need to set n to the value
    return *this;
}

template<class T>
 VISArray<T>& VISArray<T>::operator=(const VISArray<T>* s) {
    if (s != this) {
	if (s != NULL) operator=(*s);
	else { 
	    exit(-1);
	}
    }
    return *this;
}


