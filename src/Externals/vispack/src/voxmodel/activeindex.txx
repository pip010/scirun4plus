#include "voxmodel/activeindex.h"

template <class T> 
ActiveIndex< T >::~ActiveIndex()
{
  delete [] data;
  delete [] unused;
}
 
template <class T>
ActiveIndex< T >::ActiveIndex (unsigned initial_size, unsigned inc_size)
{
  unsigned i;
  data = new T[initial_size];
  unused = new unsigned[initial_size];
  if ((data==0)||(unused==0)) {
    std::cout<<"--- "<<T::name()<<" information ---\n";// **mdm**
    std::cout<<"Out of memory!\n";// **mdm**
    exit(-1);
  }
  for (i=0;i<initial_size;i++) unused[initial_size-1-i]=i;
  std::cout<<"--- "<<T::name()<<" information ---\n";// **mdm**
  std::cout<<"Initial memory size is : "<<initial_size<<std::endl;// **mdm**
  available=initial_size;
  size_increment=inc_size;
  current_size=initial_size;
  max_accessed=0;
  n_fixed=0;
}

template <class T>
void ActiveIndex< T >::expand ()
{
  T *swap;
  unsigned *swap2, i;
  
  std::cout << "Expanding memory allocation ("<<T::name()<<") ...";// **mdm**
  swap = new T[current_size+size_increment];
  swap2 = new unsigned[current_size+size_increment];

  if ((swap==0)||(swap2==0)) {
    std::cout<<"Out of memory!\n";// **mdm**
    exit(-1);
  }

  for (i=0;i<current_size;i++) swap[i]=data[i]; // = operator must be defined for class T

  delete [] data;
  delete [] unused;
  data=swap;
  unused=swap2;

  for (i=0;i<size_increment;i++) unused[size_increment-1-i]=current_size++;
  available=size_increment;
  std::cout<<"done\n";// **mdm**
  info();
}

template <class T>
void ActiveIndex< T >::edit (unsigned index, T &new_data) { 
  if (index<n_fixed) std::cout<<"Attempt to edit fixed entry! ("<<T::name()<<")\n";// **mdm**
  else {
    if (index<current_size) data[index]=new_data; // = operator must be defined for class T
    else {
      std::cout<<"ActiveIndex error: Attempt to edit non existing data entry ("<<T::name()<<")\n";// **mdm**
      exit(-1);
    }
  }
}

template <class T>
unsigned ActiveIndex< T >::push (T &new_data) {
  if (available==0) expand();
  unsigned i=unused[--available];
  data[i]=new_data;             // = operator must be defined for class T
  if ((i+1)>max_accessed) max_accessed=i+1;
  return i;
} 

template <class T>
void ActiveIndex< T >::info ()
{
  std::cout<<"--- "<<T::name()<<" information ---\n";// **mdm**
  // T must have static name() function to return its class name
  std::cout<<"Data cells allocated             : "<<current_size<<std::endl;// **mdm**
  std::cout<<"Max. data cell used              : "<<max_accessed<<std::endl;// **mdm**
  std::cout<<"Current number of available cells: "<<available<<std::endl;// **mdm**
  std::cout<<"Next available slot              : "<<unused[available-1]<<std::endl;// **mdm**
}
