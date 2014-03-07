/*
   For more information, please see: http://software.sci.utah.edu

   The MIT License

   Copyright (c) 2009 Scientific Computing and Imaging Institute,
   University of Utah.

   License for the specific language governing rights and limitations under
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


#ifndef DATAFLOW_NETWORK_PORT_MANAGER_H
#define DATAFLOW_NETWORK_PORT_MANAGER_H 1

#include <Core/Containers/LockingHandle.h>
#include <string>
#include <vector>
#include <map>

#include <Dataflow/Network/share.h>

namespace SCIRun {

class Module;

typedef std::multimap<std::string,int> port_map_type;
typedef std::pair<port_map_type::iterator,
		  port_map_type::iterator> port_range_type;


// PortManager class has been altered to be thread safe
// so the module code as well as the network editor can
// update the number of ports.
// The PortManager will use Handle instead of pointers 
// as well so we can delete a port on the network
// while executing and not crash the network.

template<class T, class MAKER> 
class PortManager {
private:
  port_map_type    namemap_;
  std::vector<T>   ports_;
  RecursiveMutex   lock_;   // Thread safety
  MAKER            port_maker_;
  Module*          module_;
  std::string      lastportname_;
  
public:
  PortManager();
  int size();
  void add(const T &item);
  void remove(int item);
  T operator[](int);
  port_range_type operator[](const std::string&);
  T get_port(int);
  std::vector<T> get_port_range(const std::string& item);  // atomic version of getting all handles with a certain name

  void lock() { lock_.lock(); }
  void unlock() { lock_.unlock(); }
  void set_dynamic_maker(MAKER maker) { port_maker_ = maker; }
  void set_module(Module* mod) { module_ = mod; }
  void set_lastportname(const std::string& name) { lastportname_ = name; }
};

template<class T,class MAKER>
PortManager<T,MAKER>::PortManager() :
  lock_("port manager lock"),
  port_maker_(0),
  module_(0)
{
}

template<class T,class MAKER>
int
PortManager<T,MAKER>::size()
{ 
  lock_.lock();
  size_t s = ports_.size();
  lock_.unlock();
  return (static_cast<int>(s)); 
}

template<class T,class MAKER>
T
PortManager<T,MAKER>::get_port(int item)
{
  lock_.lock();
  T handle(0);
  if (item < static_cast<int>(ports_.size())) 
  {
    handle = ports_[item];
  }
  lock_.unlock();
  return(handle);
}

template<class T,class MAKER>
void
PortManager<T,MAKER>::add(const T &item)
{ 
  lock_.lock();
  namemap_.insert(std::pair<std::string, int>(item->get_portname(), ports_.size())); 
  ports_.push_back(item);
  lock_.unlock();
}

template<class T,class MAKER>
void
PortManager<T,MAKER>::remove(int item)
{
  lock_.lock();
 
  if (static_cast<int>(ports_.size()) <= item)
  {
    lock_.unlock();
    throw "PortManager tried to remove a port that does not exist";
  }

  std::string name = ports_[item]->get_portname();
  port_map_type::iterator erase_me;

  port_range_type p = namemap_.equal_range(name);
  for (port_map_type::iterator i=p.first;i!=p.second;i++)
    if ((*i).second>item)
      (*i).second--;
    else if ((*i).second==item)
      erase_me = i;
    
  ports_.erase(ports_.begin() + item);
  namemap_.erase(erase_me);
  lock_.unlock();  
}

template<class T,class MAKER>
T
PortManager<T,MAKER>::operator[](int item)
{
  lock_.lock();

  // This is patch to load dynamic ports properly
  // With subnetworks the order in which ports are allocated is improper.  
  if (port_maker_)
  { // it is dynamic we need to have an empty port next to the current one
    if (ports_.size() <= static_cast<size_t>(item))
    {
      while (ports_.size() <= static_cast<size_t>(item))
      {
        ports_.push_back(port_maker_(module_,lastportname_));
      }
    }
  }
  else
  {
    if (ports_.size() <= static_cast<size_t>(item))
    {
      lock_.unlock();
      throw "PortManager tried to access a port that does not exist";
    }
  }

  T t = ports_[item];
  lock_.unlock();
  return (t);

}

template<class T,class MAKER>
port_range_type
PortManager<T,MAKER>::operator[](const std::string& item)
{
  lock_.lock();
  port_range_type prt = static_cast<port_range_type>(namemap_.equal_range(item));
  lock_.unlock();
  return (prt);
}

template<class T, class MAKER>
std::vector<T>
PortManager<T,MAKER>::get_port_range(const std::string& item)
{
  std::vector<T> ports;
  lock_.lock();
  port_range_type range = static_cast<port_range_type>(namemap_.equal_range(item));
  port_map_type::iterator pi = range.first;
  while (pi != range.second)
  {
    ports.push_back(ports_[pi->second]);
    ++pi;
  }
  lock_.unlock();
  return (ports);
}

}

#endif
