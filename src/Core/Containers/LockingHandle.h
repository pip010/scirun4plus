/*
   For more information, please see: http://software.sci.utah.edu

   The MIT License

   Copyright (c) 2009 Scientific Computing and Imaging Institute,
   University of Utah.

   
   Permission is hereby granted, free of charge, to any person obtaining a
   copy of this software and associated documentation files (the "Software"),
   to deal in the Software without restriction, including without limitation
   the rights to use, copy, modify, merge, publish, distribute, sublicense,
   and/or sell copies of the Software, and to permit persons to whom the
   Software is furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included
   in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
   OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
   THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
   DEALINGS IN THE SOFTWARE.
*/



/*
 *  LockingHandle.h: Smart Pointers
 *
 *  Written by:
 *   Steven G. Parker
 *   Department of Computer Science
 *   University of Utah
 *   April 1994
 *
 */

#ifndef SCI_Containers_LockingHandle_h
#define SCI_Containers_LockingHandle_h 1

#ifndef SCI_NOPERSISTENT
#include <Core/Persistent/Persistent.h>
#endif

#include <boost/type_traits.hpp>
#include <boost/static_assert.hpp>
#include <typeinfo>
#include <Core/Thread/UsedWithLockingHandle.h>

namespace SCIRun {


template<class T>
class LockingHandle;

#ifndef SCI_NOPERSISTENT
template<class T>
void Pio(Piostream& stream, LockingHandle<T>& data);
#endif

template<class T>
class LockingHandle {
  T* rep;
public:

  typedef T   value_type;
  typedef T * pointer_type;

  LockingHandle();
  LockingHandle(T*);
  LockingHandle(const LockingHandle<T>&);
  template <class U>
  LockingHandle(const LockingHandle<U>&);
  LockingHandle<T>& operator=(const LockingHandle<T>&);
  LockingHandle<T>& operator=(T*);
  bool operator==(const LockingHandle<T>&) const;
  bool operator!=(const LockingHandle<T>&) const;

  ~LockingHandle();

  void detach();

  inline T* operator->() const 
  { 
    if(rep == 0)
    {
      ASSERTFAIL("Empty LockingHandle of type " + std::string((typeid(*this)).name())); 
    }
    return (rep); 
  }
  
  inline T* get_rep() const 
  { 
    return (rep); 
  }

  inline T& operator*() const { return *get_rep(); }

#ifndef SCI_NOPERSISTENT
  friend void Pio <> (Piostream& stream,LockingHandle<T>& data);
#endif
};

template<class T>
LockingHandle<T>::LockingHandle()
: rep(0)
{
}

template<class T>
LockingHandle<T>::LockingHandle(T* rep)
: rep(rep)
{
  if(rep)
  {
    rep->lock.lock();
    rep->ref_cnt++;
    rep->lock.unlock();
  }
}

template<class T>
LockingHandle<T>::LockingHandle(const LockingHandle<T>& copy)
: rep(copy.rep)
{
  if(rep)
  {
    rep->lock.lock();
    rep->ref_cnt++;
    rep->lock.unlock();
  }
}

template<class T>
template<class U>
LockingHandle<T>::LockingHandle(const LockingHandle<U>& copy)
: rep(copy.get_rep())
{
  if(rep)
  {
    rep->lock.lock();
    rep->ref_cnt++;
    rep->lock.unlock();
  }
}

template<class T>
LockingHandle<T>& LockingHandle<T>::operator=(const LockingHandle<T>& copy)
{
  if(rep != copy.rep)
  {
    if(rep)
    {
	    rep->lock.lock();
	    if(--rep->ref_cnt==0)
      {
        rep->lock.unlock();
        delete rep;
	    } 
      else 
      {
        rep->lock.unlock();
	    }
    }
    if(copy.rep)
    {
	    copy.rep->lock.lock();
	    rep=copy.rep;
	    rep->ref_cnt++;
	    copy.rep->lock.unlock();
    } 
    else 
    {
	    rep=copy.rep;
    }
  }
  return *this;
}

template<class T>
LockingHandle<T>& LockingHandle<T>::operator=(T* crep)
{
  if(rep)
  {
    rep->lock.lock();
    if(--rep->ref_cnt==0)
    {
      rep->lock.unlock();
      delete rep;
    } 
    else 
    {
      rep->lock.unlock();
    }
  }
  if(crep)
  {
    crep->lock.lock();
    rep=crep;
    rep->ref_cnt++;
    crep->lock.unlock();
  } 
  else 
  {
    rep=crep;
  }
  return *this;
}

template<class T>
bool LockingHandle<T>::operator==(const LockingHandle<T>& crep) const
{
  return (get_rep() == crep.get_rep());
}

template<class T>
bool LockingHandle<T>::operator!=(const LockingHandle<T>& crep) const
{
  return (get_rep() != crep.get_rep());
}

template<class T>
LockingHandle<T>::~LockingHandle()
{
  typedef boost::is_base_of<UsedWithLockingHandleBase, T> base_of_type;
  BOOST_STATIC_ASSERT(base_of_type::value);
  if(rep)
  {
    rep->lock.lock();
    if(--rep->ref_cnt==0)
    {
	    rep->lock.unlock();
	    delete rep;
    } 
    else 
    {
	    rep->lock.unlock();
    }
  }
}

template<class T>
void LockingHandle<T>::detach()
{
  if(rep == 0) 
  {
    ASSERTFAIL("Empty LockingHandle of type " + std::string((typeid(*this)).name()) + " (detach handle)");
  }
  rep->lock.lock();
  if(rep->ref_cnt==1)
  {
    rep->lock.unlock();
    rep->generation = rep->compute_new_generation();
    return; // We have the only copy
  }
  T* oldrep=rep;
  rep=oldrep->clone();
  oldrep->ref_cnt--;
  oldrep->lock.unlock();
  rep->ref_cnt++;
}



#ifndef SCI_NOPERSISTENT
template<class T>
void Pio(Piostream& stream, LockingHandle<T>& data)
{
  if (stream.reading())
  {
    if (data.rep)
    {
      data.rep->lock.lock();
      if(--data.rep->ref_cnt==0)
      {
        data.rep->lock.unlock();
        delete data.rep;
      } 
      else 
      {
        data.rep->lock.unlock();
      }    
    }
  }
  stream.begin_cheap_delim();
  Persistent* trep=data.rep;
  stream.io(trep, T::type_id);
  if(stream.reading())
  {
    data.rep=static_cast<T*>(trep);
    if(data.rep)
      data.rep->ref_cnt++;
  }
  stream.end_cheap_delim();
}
#endif

} // End namespace SCIRun


#endif
