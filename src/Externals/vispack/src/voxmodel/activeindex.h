#ifndef ACTIVEINDEX
#define ACTIVEINDEX

#include <iostream>// **mdm**

template <class T> class ActiveIndex
{
 public: 
  ActiveIndex(unsigned, unsigned);
  ~ActiveIndex();

  T* pop (unsigned index) {
    if (index<max_accessed) return &(data[index]);
    else {
      std::cout<<"Attempt to pop invalid index ("<<T::name()<<"):"<<index<<" / "<<max_accessed<<" \n";// **mdm**
      exit(-1);
    }
  } 
  void remove (unsigned index) 
    {
      if (index<n_fixed) std::cout<<"Attempt to remove fixed element!\n";// **mdm**
      else unused[available++]=index;
    }
  void edit (unsigned, T&); 
  unsigned push (T&);
  void info ();

 private:
  T *data;
  unsigned *unused;        // list of unused data slots
  unsigned available;      // number of available slots in unused list
  unsigned size_increment; 
  unsigned current_size;
  unsigned max_accessed;
  unsigned n_fixed;        // the first n_fixed elements of data are wrote once and after that are readonly

  void expand ();
};

#ifndef MANUAL_INSTANTIATION
#include "voxmodel/activeindex.txx"
#endif

#endif





