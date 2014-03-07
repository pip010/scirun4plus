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
 *  GeometryPort.cc: Handle to the Geometry Data type
 *
 *  Written by:
 *   Steven G. Parker
 *   Department of Computer Science
 *   University of Utah
 *   April 1994
 *
 */

#include <Dataflow/Network/Ports/GeometryPort.h>
#include <Dataflow/Network/Ports/GeometryComm.h>

#include <Dataflow/Network/Connection.h>
#include <Dataflow/Network/Module.h>
#include <Dataflow/Network/Ports/Port.h>
#include <Core/Geom/View.h>

#include <Core/Thread/FutureValue.h>
#include <Core/Util/Assert.h>
#include <Core/Util/MemoryUtil.h>

#include <iostream>
using std::cerr;

#undef SCISHARE
#if defined(_WIN32) && !defined(BUILD_SCIRUN_STATIC)
#  define SCISHARE __declspec(dllexport)
#else
#  define SCISHARE
#endif

namespace SCIRun {

extern "C" {
  SCISHARE IPort* make_GeometryIPort(Module* module, const std::string& name) 
  {
    return new GeometryIPort(module,name);
  }
  SCISHARE OPort* make_GeometryOPort(Module* module, const std::string& name) 
  {
    return new GeometryOPort(module,name);
  }
}


static std::string Geometry_type("Geometry");
static std::string Geometry_color("magenta3");

//! Create the geometry input port
GeometryIPort::GeometryIPort(Module* module, const std::string& portname)
  : IPort(module, Geometry_type, portname, Geometry_color)
{
}

//! Destructor
GeometryIPort::~GeometryIPort()
{
}


//! Create the geometry output port
GeometryOPort::GeometryOPort(Module* module, const std::string& portname)
  : OPort(module, Geometry_type, portname, Geometry_color),
    serial_(1), lserial_(4),
    dirty_(false),
    reply_mailbox_("GeometryOPort Response mailbox", 0)
{
}


//! Destructor of output port class
GeometryOPort::~GeometryOPort()
{
  //! Clean up all the messages stored in the port
  //! messages are the data cached in the port
  delete_all_items(saved_msgs_);
  saved_msgs_.clear();
}


void
GeometryOPort::reset()
{
  dirty_ = false;
}


//! Send a flush message to every module connected to this one
void
GeometryOPort::flush()
{
  //! Send a message to every connected module
  for (size_t i=0; i<modules_.size(); i++)
  {
    GeometryComm* msg = new GeometryComm(MessageTypes::GeometryFlush,portid_[i]);
    modules_[i]->mailbox_.send(msg);
  }
  //! Store the message for modules that will be connected in the future
  GeometryComm* msg = new GeometryComm(MessageTypes::GeometryFlush, 0);
  save_msg(msg);
}
  

//! Flush all the output.
//! This function is trigger at the end of do_execute in the module code
void
GeometryOPort::finish()
{
  if (dirty_)
  {
    turn_on_light(Finishing);
    flush();
    turn_off_light();
  }
}


void
GeometryOPort::synchronize()
{
  for (size_t i=0; i<modules_.size(); i++)
  {
    GeometryComm* msg =
      new GeometryComm(MessageTypes::GeometrySynchronize, portid_[i]);
    modules_[i]->mailbox_.send(msg);
  }
}


GeomID
GeometryOPort::addObj(GeomHandle obj, const std::string& name, CrowdMonitor* lock)
{
  turn_on_light();
  GeomID id = serial_++;
  for (unsigned int i = 0; i < modules_.size(); i++)
  {
    GeometryComm* msg = new GeometryComm(portid_[i], id, obj, name, lock);
    modules_[i]->mailbox_.send(msg);
  }
  
  GeometryComm* msg = new GeometryComm(0, id, obj, name, lock);
  save_msg(msg);
  dirty_ = true;
  turn_off_light();
  return id;
}

LightID
GeometryOPort::addLight(LightHandle obj, 
			const std::string& name, CrowdMonitor* lock)
{
  turn_on_light();
  LightID id = lserial_++;

  for (unsigned int i = 0; i < modules_.size(); i++)
  {
    GeometryComm* msg = new GeometryComm(portid_[i], id, obj, name, lock);
    modules_[i]->mailbox_.send(msg);
  }
  
  GeometryComm* msg = new GeometryComm(0, id, obj, name, lock);
  save_msg(msg);
  dirty_ = true;
  turn_off_light();
  return id;
}



bool
GeometryOPort::direct_forward(GeometryComm* msg)
{
  if (modules_.size() == 0) { return false; }

  // Send msg to last port directly, but copy to all of the prior ones.
  // Note that msg is not deleted if this function returns false.
  unsigned int i;
  for (i = 0; i < modules_.size()-1; i++)
  {
    GeometryComm *cpy = new GeometryComm(*msg);
    modules_[i]->mailbox_.send(cpy);
  }
  modules_[i]->mailbox_.send(msg);

  return true;
}



void
GeometryOPort::forward(GeometryComm* msg0)
{
  for (unsigned int i = 0; i < modules_.size(); i++)
  {
    GeometryComm *msg = new GeometryComm(*msg0);
    msg->portno = portid_[i];
    modules_[i]->mailbox_.send(msg);
  }
  save_msg(msg0);
}



void
GeometryOPort::delObj(GeomID id, int /*del*/)
{
  turn_on_light();

  for (unsigned int i=0; i < modules_.size(); i++)
  {
    GeometryComm* msg = new GeometryComm(portid_[i], id);
    modules_[i]->mailbox_.send(msg);
  }
  GeometryComm* msg = new GeometryComm(0, id);
  save_msg(msg);
  dirty_ = true;
  turn_off_light();
}

void
GeometryOPort::delLight(LightID id, int /*del*/)
{
  turn_on_light();

  for (unsigned int i=0; i < modules_.size(); i++)
  {
    GeometryComm* msg = new GeometryComm(portid_[i], id);
    modules_[i]->mailbox_.send(msg);
  }
  GeometryComm* msg = new GeometryComm(0, id);
  save_msg(msg);
  dirty_ = true;
  turn_off_light();
}



void
GeometryOPort::delAll()
{
  turn_on_light();
  for (unsigned int i = 0; i < modules_.size(); i++)
  {
    GeometryComm* msg = new GeometryComm(MessageTypes::GeometryDelAll,
					    portid_[i]);
    modules_[i]->mailbox_.send(msg);
  }
  GeometryComm* msg = new GeometryComm(MessageTypes::GeometryDelAll, 0);
  save_msg(msg);
  dirty_ = true;
  turn_off_light();
}



void
GeometryOPort::flushViews()
{
  turn_on_light();
  for (unsigned int i = 0; i < modules_.size(); i++)
  {
    GeometryComm* msg = new GeometryComm(MessageTypes::GeometryFlushViews,
					    portid_[i], (Semaphore*)0);
    modules_[i]->mailbox_.send(msg);
  }
  GeometryComm* msg = new GeometryComm(MessageTypes::GeometryFlushViews,
					  0, (Semaphore*)0);
  save_msg(msg);
  dirty_ = false;
  turn_off_light();
}

bool
GeometryOPort::get_view_bounds(BBox &bbox)
{

  GeometryData *data;
  data = getData(0, 0, GEOM_VIEW_BOUNDS);
  if (data) {
    bbox = data->view_bounds_;
    return true;
  }
  return false;
}


void
GeometryOPort::flushViewsAndWait()
{
  turn_on_light();
  Semaphore waiter("flushViewsAndWait wait semaphore", 0);
  for (unsigned int i = 0; i < modules_.size(); i++)
  {
    GeometryComm* msg = new GeometryComm(MessageTypes::GeometryFlushViews,
					    portid_[i], &waiter);
    modules_[i]->mailbox_.send(msg);
  }
  GeometryComm* msg = new GeometryComm(MessageTypes::GeometryFlushViews,
					  0, &waiter);
  save_msg(msg); // TODO:  Should a synchronized primitive be queued?

  // Wait on everyone.
  for (unsigned int i = 0; i < modules_.size(); i++)
  {
    waiter.down();
  }
  dirty_ = false;
  turn_off_light();
}


void
GeometryOPort::save_msg(GeometryComm* msg)
{
  //! Update the cached messages. These are used when a new viewer is connected
  //! to this port
  switch(msg->type)
  {
  case MessageTypes::GeometryDelObj:
    {
      //! Delete the object from the message queue.
      //! The object is identified with the serial number
      //! This message is never put on the stack
      //! instead the message that adds it is deleted.
      std::list<GeometryComm *>::iterator itr0 = saved_msgs_.begin();
      while (itr0 != saved_msgs_.end())
      {
        std::list<GeometryComm *>::iterator itr  = itr0;
        ++itr0;
        
        //! check whether this message added the object
        if ((*itr)->type == MessageTypes::GeometryAddObj &&
            (*itr)->serial == msg->serial)
        {
          delete *itr;
          saved_msgs_.erase(itr);
        }
      }
      delete msg;
      return;
    }
    break;
    
  case MessageTypes::GeometryDelLight:
    {
      //! Delete the object from the message queue.
      //! The object is identified with the serial number
      //! This message is never put on the stack
      //! instead the message that adds it is deleted.
      std::list<GeometryComm *>::iterator itr0 = saved_msgs_.begin();
      while (itr0 != saved_msgs_.end())
      {
        std::list<GeometryComm *>::iterator itr  = itr0;
        ++itr0;
        
        //! check whether this message added the object
        if ((*itr)->type == MessageTypes::GeometryAddLight &&
            (*itr)->lserial == msg->lserial)
        {
          delete *itr;
          saved_msgs_.erase(itr);
        }
      }
      delete msg;
      return;
    }
    break;

  case MessageTypes::GeometryDelAll:
    {
      //! Delete all AddObj messages from the queue. This way when
      //! new viewer are added no objects are created.
      std::list<GeometryComm *>::iterator itr0 = saved_msgs_.begin();
      while (itr0 != saved_msgs_.end())
      {
        std::list<GeometryComm *>::iterator itr  = itr0;
        ++itr0;
        
        //! Check whether it is a AddObj message
        if ((*itr)->type == MessageTypes::GeometryAddObj)
        {
          delete *itr;
          saved_msgs_.erase(itr);
        }
      }
      delete msg;
      return;
    }
    break;

  case MessageTypes::GeometryFlush:
    {
      //! Delete an older flush messages on the message stack and replace it
      //! with this new one.
      std::list<GeometryComm *>::iterator itr0 = saved_msgs_.begin();
      while (itr0 != saved_msgs_.end())
      {
        std::list<GeometryComm *>::iterator itr  = itr0;
        ++itr0;
        
        //! If it is a flush message delete it
        if ((*itr)->type == MessageTypes::GeometryFlush)
        {
          delete *itr;
          saved_msgs_.erase(itr);
        }
      }
      
      //! add the flush message to the list of pending messages
      saved_msgs_.push_back(msg);
      return;
    }
    break;

  case MessageTypes::GeometryFlushViews:
    {
      //! Delete an older flush messages on the message stack and replace it
      //! with this new one.
      std::list<GeometryComm *>::iterator itr0 = saved_msgs_.begin();
      while (itr0 != saved_msgs_.end())
      {
        std::list<GeometryComm *>::iterator itr  = itr0;
        ++itr0;
        
        //! If it is a flushviews message delete it
        if ((*itr)->type == MessageTypes::GeometryFlushViews)
        {
          delete *itr;
          saved_msgs_.erase(itr);
        }
      }

      //! add the flushviews message to the list of pending messages
      saved_msgs_.push_back(msg);
      return;
    }
    break;

  case MessageTypes::GeometrySetView:
    {
      //! Delete an older setview messages on the message stack and replace it
      //! with this new one.
      std::list<GeometryComm *>::iterator itr0 = saved_msgs_.begin();
      while (itr0 != saved_msgs_.end())
      {
        std::list<GeometryComm *>::iterator itr  = itr0;
        ++itr0;
        
        //! If it is a setview message delete it
        if ((*itr)->type == MessageTypes::GeometrySetView)
        {
          delete *itr;
          saved_msgs_.erase(itr);
        }
      }

      //! add the setview message to the list of pending messages
      saved_msgs_.push_back(msg);
    }
    break;

  default:
    //! any other message (addobj) is added to the queue
    saved_msgs_.push_back(msg);
  }
}


void
GeometryOPort::attach(Connection* c)
{
  //! store connection handle locally in the port
  OPort::attach(c);

  //! what is the new number for this one.
  int which = modules_.size();

  // Set up the modules_ and portid_ variables.
  turn_on_light(Resetting);
  
  //! Get the module handle
  ModuleHandle mod = c->iport->get_module();
  modules_.push_back(mod);
    
  //! Send the registration message.
  //! This piece of code sends a pointer to a mailbox
  //! to receive the portid
  modules_[which]->mailbox_.send(new GeometryComm(&reply_mailbox_));
  GeomReply reply = reply_mailbox_.receive();

  portid_.push_back(reply.portid);
  
  // Forward all of the queued up messages.
  turn_on_light();
  
  // Send all the geometry that is cached in this port to the
  // receiving port.
  std::list<GeometryComm *>::iterator itr = saved_msgs_.begin();
  while (itr != saved_msgs_.end())
  {
    GeometryComm *msg = new GeometryComm(**itr);
    msg->portno = portid_[which];
    modules_[which]->mailbox_.send(msg);
    ++itr;
  }
  
  //! done
  turn_off_light();
}


void
GeometryOPort::detach(Connection* c, bool blocked)
{
  // Determine which connection gets it.
  unsigned int i;
  for (i = 0; i < connections.size(); i++)
  {
    if (connections[i].get_rep() == c)
    {
      break;
    }
  }
  
  if (i<connections.size())
  {
    // Let the Viewer know that the port is shutting down.
    GeometryComm *msg =
      new GeometryComm(MessageTypes::GeometryDetach, portid_[i]);
    
    modules_[i]->mailbox_.send(msg);

    //! Clean up the modules_ and portid_ vectors.
    //! Release the handle on the port so the module is not kept in existence
    //! for this connection.
    modules_.erase(modules_.begin() + i);
    portid_.erase(portid_.begin() + i);
  }

  OPort::detach(c, blocked);
}


bool
GeometryOPort::have_data()
{
  //! check whether there are saved messages
  return saved_msgs_.size();
}


//! this one should not do anything
void
GeometryOPort::resend(ConnectionHandle)
{
}

GeometryData *
GeometryOPort::getData(int which_viewer, int which_viewwindow, int datamask)
{
  //! Check whether we are in range
  if (which_viewer >= (int)modules_.size() || which_viewer < 0) return 0;

  FutureValue<GeometryData*> reply("Geometry getData reply");
  GeometryComm *msg = new GeometryComm(MessageTypes::GeometryGetData,
					  portid_[which_viewer],
					  &reply, which_viewwindow,
					  datamask);
  //! Send out the request for data
  modules_[which_viewer]->mailbox_.send(msg);
  
  //! Wait for other thread to respond
  return reply.receive();
}


void
GeometryOPort::setView(int which_viewer, int which_viewwindow, View view)
{
  //! Check whether we are in range
  if (which_viewer >= (int)modules_.size() || which_viewer < 0) return;

  GeometryComm* msg = new GeometryComm(MessageTypes::GeometrySetView,
					  portid_[which_viewer],
					  which_viewwindow,
					  view);
  modules_[which_viewer]->mailbox_.send(msg);
}




//! The memory for reply is not owned by this.
GeometryComm::GeometryComm(Mailbox<GeomReply> *reply)
  : MessageBase(MessageTypes::GeometryInit),
    reply(reply)
{
}


GeometryComm::GeometryComm(int portno, GeomID serial, GeomHandle obj,
			   const std::string& name, CrowdMonitor* lock)
  : MessageBase(MessageTypes::GeometryAddObj),
    reply(0),
    portno(portno),
    serial(serial),
    obj(obj),
    name(name),
    lock(lock)
{
}
GeometryComm::GeometryComm(int portno, LightID serial, LightHandle light,
			   const std::string& name, CrowdMonitor* lock)
  : MessageBase(MessageTypes::GeometryAddLight),
    reply(0),
    portno(portno),
    lserial(serial),
    light(light),
    name(name),
    lock(lock)
{
}


GeometryComm::GeometryComm(int portno, GeomID serial)
  : MessageBase(MessageTypes::GeometryDelObj),
    reply(0),
    portno(portno),
    serial(serial)
{
}

GeometryComm::GeometryComm(int portno, LightID serial)
  : MessageBase(MessageTypes::GeometryDelLight),
    reply(0),
    portno(portno),
    lserial(serial)
{
}


GeometryComm::GeometryComm(MessageTypes::MessageType type,
			   int portno, Semaphore* wait)
  : MessageBase(type),
    reply(0),
    portno(portno),
    wait(wait)
{
}


GeometryComm::GeometryComm(MessageTypes::MessageType type, int portno)
  : MessageBase(type),
    reply(0),
    portno(portno)
{
}


GeometryComm::GeometryComm(MessageTypes::MessageType type, int portno,
			   FutureValue<GeometryData*>* datareply,
			   int which_viewwindow, int datamask)
  : MessageBase(type),
    reply(0),
    portno(portno),
    which_viewwindow(which_viewwindow),
    datamask(datamask),
    datareply(datareply)
{
}


GeometryComm::GeometryComm(MessageTypes::MessageType type, int portno,
			   int which_viewwindow, View view)
  : MessageBase(type),
    reply(0),
    portno(portno),
    view(view),
    which_viewwindow(which_viewwindow)
{
}


GeometryComm::GeometryComm(MessageTypes::MessageType type, int portno,
			   FutureValue<int>* nreply)
  : MessageBase(type),
    reply(0),
    portno(portno),
    nreply(nreply)
{
}


GeometryComm::GeometryComm(const GeometryComm &copy)
  : MessageBase(copy),
    reply(copy.reply),
    portno(copy.portno),
    serial(copy.serial),
    obj(copy.obj),
    lserial(copy.lserial),
    light(copy.light),
    name(copy.name),
    lock(0),
    wait(0),
    view(copy.view),
    next(0),
    which_viewwindow(copy.which_viewwindow),
    datamask(copy.datamask),
    datareply(0),
    nreply(0)
{
}


GeometryComm::~GeometryComm()
{
}



GeomReply::GeomReply()
{
}


GeomReply::GeomReply(int portid)
  : portid(portid)
{
}



GeometryData::GeometryData()
{
  view=0;
  colorbuffer=0;
  depthbuffer=0;
}

GeometryData::~GeometryData()
{
}


void
GeometryData::Print()
{
  cerr << "GEOMETRY data review\n\n";
  cerr << "X resolution: " << xres << " Y resolution: " << yres << std::endl;
  cerr << "Clipping planes.  Near = " << znear << " Far = " << zfar << std::endl;

  if ( depthbuffer == NULL )
    cerr << "depthbuffer has nothing\n";

  if ( colorbuffer == NULL )
    cerr << "colorbuffer has nothing\n";

  if ( view == NULL )
    cerr << "view has nothing\n";

  cerr << std::endl;
}

} // End namespace SCIRun

