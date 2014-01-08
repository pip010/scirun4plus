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

#include <Dataflow/Network/Ports/SimplePort.h>
#include <Dataflow/Network/Module.h>

namespace SCIRun {

SimpleIPort::SimpleIPort(Module* module,
			    const std::string& port_name,
          const std::string& port_type,
          const std::string& port_color)
  : IPort(module, port_type, port_name, port_color),
    mailbox("Port mailbox (SimpleIPort)", 2),
    deactivated_(false),
    generation_( -1 ),
    last_generation_( -1 )
{
  DEBUG_CONSTRUCTOR("SimpleIPort") 
}


SimpleIPort::~SimpleIPort()
{
  DEBUG_DESTRUCTOR("SimpleIPort") 
}


SimpleOPort::SimpleOPort(Module* module,
			    const std::string& port_name,
          const std::string& port_type,
          const std::string& port_color)
  : OPort(module, port_type, port_name, port_color),
    cache_(true),
    sent_something_(true),
    handle_(0)
{
  DEBUG_CONSTRUCTOR("SimpleOPort")
  if (sci_getenv_p("SCIRUN_NO_PORT_CACHING"))
  {
    issue_no_port_caching_warning();
    cache_ = false;
  }
}


SimpleOPort::~SimpleOPort()
{
  DEBUG_DESTRUCTOR("SimpleOPort")
}


void SimpleIPort::reset()
{
  got_something_ = false;
}


void
SimpleIPort::finish()
{
  if (!got_something_ && nconnections() > 0)
  {
    turn_on_light(Finishing);
    SimplePortComm* msg = mailbox.receive();
    delete msg;
    turn_off_light();
  }
  got_something_ = true;
}


void
SimpleOPort::reset()
{
  sent_something_ = false;
}


void
SimpleOPort::do_not_send()
{
  sent_something_ = true;
}


void
SimpleOPort::finish()
{
  if (!sent_something_ && nconnections() > 0)
  {
    // Tell them that we didn't send anything.
    turn_on_light(Finishing);
    
    for (int i = 0; i < nconnections(); i++)
    {
      ConnectionHandle conn = connections[i].get_rep();
      SimplePortComm* msg = new SimplePortComm(handle_);
      dynamic_cast< SimpleIPort* >(conn->iport.get_rep())->mailbox.send(msg);
    }
    turn_off_light();
  }
  sent_something_ = true;
}


void
SimpleOPort::detach(Connection* conn, bool blocked)
{
  if (!sent_something_)
  {
    SimplePortComm* msg = new SimplePortComm(0);
    dynamic_cast<SimpleIPort*>(conn->iport.get_rep())->mailbox.send(msg);
  }
  //sent_something_ = true;  // Only sent something on the one port.
  OPort::detach(conn, blocked);
}


void
SimpleOPort::do_send(DatatypeHandle data, SendType type, DerefType deref)
{
  handle_ = cache_ ? data : 0;

  if (nconnections() == 0) { return; }

  // Change oport state and colors on screen.
  turn_on_light();

  if( type == SEND_INTERMEDIATE )
  {
    module->request_multisend(this);
  }
  for (int i = 0; i < nconnections(); i++)
  {
    // Add the new message.
    ConnectionHandle conn = connections[i];
    SimplePortComm* msg = new SimplePortComm(data);
    if (i == nconnections()-1 &&
        (deref == DEREF_ALWAYS ||
         (deref == DEREF_NOCACHE && !handle_.get_rep())))
    {
      data = 0;
    }
    dynamic_cast<SimpleIPort*>(conn->iport.get_rep())->mailbox.send(msg);
  }
  
  if (nconnections() == 0 &&
      (deref == DEREF_ALWAYS ||
       (deref == DEREF_NOCACHE && !handle_.get_rep())))
  {
    data = 0;
  }

  sent_something_ = true;

  // Change oport state and colors on screen.
  turn_off_light();
}


void
SimpleIPort::deactivate()
{

  lock.lock();
  deactivated_ = true;
  lock.unlock();
 
  // free any waiting thread by sending
  // a dummy message.
  
  SimplePortComm* msg = new SimplePortComm(); 
  mailbox.send(msg);
}


bool
SimpleIPort::changed()
{
  return (generation_ != last_generation_);
}


bool
SimpleOPort::have_data()
{
  return (handle_.get_rep()!=0);
}


void
SimpleOPort::resend(ConnectionHandle conn)
{
  turn_on_light();
  for (int i = 0; i < nconnections(); i++)
  {
    ConnectionHandle c = connections[i];
    if (c == conn)
    {
      SimplePortComm* msg = new SimplePortComm(handle_);
      dynamic_cast<SimpleIPort*>(c->iport.get_rep())->mailbox.send(msg);
    }
  }
  turn_off_light();
}


SimplePortComm::SimplePortComm() :
  have_data_(false)
{
}


void 
SimpleOPort::set_cache(bool cache)
{
  cache_ = cache;
  if ( !cache )
  {
    handle_ = 0;
  }
}


bool 
SimpleOPort::cache_flag_supported() 
{ 
  return true; 
}

bool 
SimpleOPort::get_cache() 
{ 
  return cache_; 
}


SimplePortComm::SimplePortComm(const DatatypeHandle& data) :
  data_(data),
  have_data_(true)
{
}


} // End namespace SCIRun
