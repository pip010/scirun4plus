// smart vector template

#ifndef __SMART_VECTOR_H__
#define __SMART_VECTOR_H__

#ifdef _WIN32
#pragma warning (disable : 4267)
#endif

#include <vector>

namespace custom_class
{
  template <class T> class svector
  {
  public:
    svector() {};
    svector(unsigned int size) { _vector.resize(size); };
    svector(const svector<T>& that) { *this = that; };
    ~svector() {};

    unsigned int size() const { return _vector.size(); };
    void resize(unsigned int size) { _vector.resize(size); };
    bool empty() const { return _vector.empty(); };
    void clear() { _vector.clear(); };
    T& front();
    T& back();
    void pop_back();
    void push_back(const T& x);

    T& operator [] (unsigned int i);
    const T& operator [] (unsigned int i) const;
    T& at(unsigned int i);
    const T& at(unsigned int i) const;
    svector<T>& operator = (const svector<T>& that);

    std::vector<T> _vector;  
  };
}

#include <svector.T>


#endif // __SMART_VECTOR_H__
