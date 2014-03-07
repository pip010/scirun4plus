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


#ifndef DATAFLOW_NETWORK_PORTS_SIMPLEPORT_H
#define DATAFLOW_NETWORK_PORTS_SIMPLEPORT_H 1

#include <Dataflow/Network/Connection.h>
#include <Dataflow/Network/Ports/Port.h>
#include <Core/Datatypes/PropertyManager.h>
#include <Core/Thread/Mailbox.h>
#include <Core/Util/Environment.h>
#include <Core/Containers/LockingHandle.h>
#include <Dataflow/Network/share.h>

namespace SCIRun {

class Module;

struct SimplePortComm
{
  SimplePortComm();
  SimplePortComm(const DatatypeHandle&);

  DatatypeHandle data_;
  bool have_data_;
};


class SCISHARE SimpleIPort : public IPort
{
public:
  SimpleIPort(Module*, const std::string& name, 
              const std::string& port_type, 
              const std::string& port_color );
  virtual ~SimpleIPort();

  Mailbox<SimplePortComm*> mailbox;

  virtual void reset();
  virtual void finish();

  template<class HandleType>
  int get(HandleType& data)
  {
    bool deactivated;

    lock.lock();
    deactivated = deactivated_;
    lock.unlock();
    if (deactivated) { return 0; }
    
    if (nconnections() == 0) { return 0; }
    turn_on_light();

    // Wait for the data.
    SimplePortComm* comm = mailbox.receive();
    got_something_ = true;
    if (comm->have_data_)
    {
      data = SCI_DATATYPE_CAST<typename HandleType::pointer_type>(comm->data_.get_rep());
      last_generation_ = generation_;
      if (data.get_rep()) 
        generation_ = data->generation;
      else
        generation_ = -1;

      delete comm;
      turn_off_light();
      return 1;
    }
    else
    {
      delete comm;
      turn_off_light();
      return 0;
    }
  }
  
  bool changed();
  
  //! when the port is deleted from a network it is
  //! deactivated. As another thread may keep a handle
  //! to the port and may be waiting for data we need
  //! to alter the state of the port. In this state
  //! the port will not receive any new data
  void deactivate();

private:
  bool got_something_;
  bool deactivated_;

  int generation_;
  int last_generation_;
};


class SCISHARE SimpleOPort : public OPort
{
public:
  SimpleOPort(Module*, 
              const std::string& port_name,
              const std::string& port_type, 
              const std::string& port_color);
  virtual ~SimpleOPort();

  virtual void reset();
  virtual void finish();
  virtual void do_not_send();
  virtual void detach(Connection* conn, bool blocked);

  template<class T> 
  void send(T& data )
  {
    if (data.get_rep())
    {
      PropertyManager *pm = dynamic_cast<PropertyManager *>(data.get_rep());
      if (pm && !pm->is_frozen())
        pm->freeze();
    }
    do_send(DatatypeHandle(data.get_rep()), SEND_NORMAL, DEREF_NEVER);
  }

  template<class T>
  void send_intermediate(T& data)
  {
    if (data.get_rep())
    {
      PropertyManager *pm = dynamic_cast<PropertyManager *>(data.get_rep());
      if (pm && !pm->is_frozen()) pm->freeze();
    }
    do_send(DatatypeHandle(data.get_rep()), SEND_INTERMEDIATE, DEREF_NEVER);
  }
  
  template<class T> void send_and_dereference(T& data, bool save_when_caching)
  {
    if (data.get_rep())
    {
      PropertyManager *pm = dynamic_cast<PropertyManager *>(data.get_rep());
      if (pm && !pm->is_frozen())
        pm->freeze();
    }
    do_send(DatatypeHandle(data.get_rep()), SEND_NORMAL, save_when_caching?DEREF_ALWAYS:DEREF_NOCACHE);
  }

  virtual bool cache_flag_supported();
  virtual bool get_cache();
  virtual void set_cache(bool cache);

  virtual bool have_data();
  virtual void resend(ConnectionHandle conn);

private:
  enum SendType {SEND_NORMAL=0, SEND_INTERMEDIATE=1};
  enum DerefType {DEREF_NEVER=0, DEREF_ALWAYS=1, DEREF_NOCACHE=2};

  void do_send(DatatypeHandle, SendType type, DerefType deref);

  bool           cache_;
  bool           sent_something_;
  DatatypeHandle handle_;
};

typedef LockingHandle<SimpleIPort> SimpleIPortHandle;
typedef LockingHandle<SimpleOPort> SimpleOPortHandle;

} // End namespace SCIRun


#endif
