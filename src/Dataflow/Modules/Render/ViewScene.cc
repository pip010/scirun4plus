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
 *  ViewScene.cc:  The Geometry ViewScene
 *
 *  Written by:
 *   Steven G. Parker
 *   Department of Computer Science
 *   University of Utah
 *   April 1994
 *
 */

#include <Dataflow/Modules/Render/ViewScene.h>
#include <Dataflow/Modules/Render/ViewWindow.h>
#include <Dataflow/Comm/MessageTypes.h>
#include <Dataflow/Network/Connection.h>
#include <Dataflow/Network/ModuleHelper.h>
#include <Dataflow/Network/Scheduler.h>
#include <Dataflow/Network/NetworkIO.h>
#include <Dataflow/Network/Network.h>
#include <Dataflow/Modules/Render/ViewGeom.h>
#include <Dataflow/Modules/Render/OpenGL.h>
#include <Core/Util/StringUtil.h>
#include <Core/Geom/DirectionalLight.h>
#include <Core/Geom/GeomObj.h>
#include <Core/Geom/HeadLight.h>
#include <Core/Geom/GeomViewerItem.h>

#include <Core/Math/MiscMath.h>
#include <Core/Thread/CleanupManager.h>
#include <Core/Thread/Guard.h>
#include <Core/Thread/FutureValue.h>
#include <Core/Thread/Time.h>
#include <Core/Util/Environment.h>

#include <iostream>

using std::cerr;
using std::endl;
using std::ostream;

namespace SCIRun {

//----------------------------------------------------------------------
DECLARE_MAKER(ViewScene)
//----------------------------------------------------------------------

ViewScene::ViewScene(GuiContext* ctx)
  : Module("ViewScene", ctx, ViewerSpecial,"Render","SCIRun"),
    geomlock_("ViewScene geometry lock"), 
    view_window_lock_("ViewWindowLock"),
    max_portno_(0),
    stop_rendering_(false),
    synchronized_debt_(0),
    delete_check_autoview_on_load_(false)
{
  DEBUG_CONSTRUCTOR("ViewScene")
  std::map<LightID, int> li;
  // Add a headlight
  lighting_.lights.add(new HeadLight("Headlight", Color(1,1,1)));
  li[0] = 0;
  for(int i = 1; i < 4; i++)
  { 
    // only set up 3 more lights
    char l[8];
    sprintf(l,"Light%d",i);
    lighting_.lights.add
      (new DirectionalLight(l, Vector(0,0,1), Color(1,1,1), false, false));
    li[i] = i;
  }
  pli_[0] = li;

  default_material_ =
    new Material(Color(.1,.1,.1), Color(.6,0,0), Color(.7,.7,.7), 50);
  have_own_dispatch_=true;

  // Create port 0 - we use this for global objects such as cameras,
  // light source icons, etc.
  ports_.addObj(new GeomViewerPort(0),0);
  max_portno_ = 1;
}

//----------------------------------------------------------------------
ViewScene::~ViewScene()
{
  view_window_lock_.lock();
  for(size_t i=0;i<view_window_.size();i++)
  {
    delete view_window_[i];
  }
  view_window_.clear();
  view_window_lock_.unlock();

// TO DO:
// Callbacks cannot be deleted to be removed in a callback call
// The callback system should be revised at some point.
/*
  if (sci_getenv("SCI_REGRESSION_TESTING"))
  {
    sched_->remove_callback(regression_save_image_callback, this);
  }
  if (sci_getenv("SCIRUN_SAVE_IMAGE_FILE"))
  {
    sched_->remove_callback(save_image_callback, this);  
  }  
*/

  DEBUG_DESTRUCTOR("ViewScene")
}

void
ViewScene::prepare_cleanup()
{
  stop_rendering_ = true;

  //! stop spinning windows.
  for(size_t i = 0; i < view_window_.size(); i++) 
  {
    view_window_[i]->gui_inertia_mode_.set(0);
  }      

  //! end all communication with gui
  for(size_t i = 0; i < view_window_.size(); i++) 
  {
    view_window_[i]->get_ctx()->inactivate();
  }    
}

//----------------------------------------------------------------------
void
ViewScene::do_execute()
{
  for(;;)
  {
    if(!stop_rendering_ && mailbox_.numItems() == 0)
    {
      // See if anything needs to be redrawn.
      int did_some=1;
      while(did_some)
      {
        did_some=0;
        view_window_lock_.lock();
        for(size_t i=0;i<view_window_.size();i++)
        {
          if (view_window_[i]) 
          {
            if (view_window_[i]->need_redraw_)
            {
              did_some++;
              ViewWindow* view = view_window_[i];
              view_window_lock_.unlock();
              view->redraw_if_needed();
              view_window_lock_.lock();
            }
          }
        }
        view_window_lock_.unlock();
      }
    }
    if (process_event() == 86)  
    {
      // Doesn't this get handled in the destructor of the viewwindow?
      view_window_lock_.lock();
      for(size_t i=0;i<view_window_.size();i++)
      {
      	ViewWindow* r=view_window_[i];
        view_window_lock_.unlock();
        if (r && r->renderer_)
        {
          r->renderer_->kill_helper();
          r->viewer_ = 0;
      	}
        view_window_lock_.lock();
      }
      view_window_lock_.unlock();
      return;
    }
  }
}


//----------------------------------------------------------------------
int 
ViewScene::process_event()
{
  MessageBase* msg=mailbox_.receive();
  GeometryComm* gmsg=(GeometryComm*)msg;
  
  switch(msg->type)
  {
    case MessageTypes::GoAway:
      return 86;

    case MessageTypes::GoAwayWarn:
      //stop spinning windows.
      break;
    case MessageTypes::ViewWindowKill:
      // only sent in the EXPERIMENTAL TCL THREAD, otherwise, as the TCL thread destroyss the ViewWindow,
      // and the OpenGL helper is drawing, they deadlock with each other
      {
        ViewSceneMessage* rmsg=(ViewSceneMessage*)msg;
        delete_viewwindow(rmsg->rid);
        break;
      }
    case MessageTypes::ExecuteModule:
      if (synchronized_debt_ < 0)
      {
        synchronized_debt_++;
        Scheduler_Module_Message *smmsg = (Scheduler_Module_Message *)msg;
        sched_->report_execution_finished(smmsg->serial);
      }
      else
      {
        Scheduler_Module_Message *smmsg = (Scheduler_Module_Message *)msg;
        synchronized_serials_.push_back(smmsg->serial);
      }
      break;

    case MessageTypes::SynchronizeModule:
      // We (mostly) ignore these messages.
      sched_->report_execution_finished(msg);
      break;

    case MessageTypes::ViewWindowRedraw:
      {
        ViewSceneMessage* rmsg=(ViewSceneMessage*)msg;
        ViewWindow* r=0;
        size_t i;
        view_window_lock_.lock();
        for(i=0;i<view_window_.size();i++)
        {
          r=view_window_[i];
          if(r->id_ == rmsg->rid) break;
        }
        if(i==view_window_.size())
        {
          warning("ViewWindow not found for redraw! (get_id()=" + rmsg->rid + ").");
        }
        else if(rmsg->nframes == 0)
        {
          // Normal redraw (lazy)
          r->need_redraw_=1;
        }
        else
        {
          // Do animation.
        view_window_lock_.unlock();
          r->redraw(rmsg->tbeg, rmsg->tend, rmsg->nframes,
              rmsg->framerate);
        view_window_lock_.lock();

        }
        view_window_lock_.unlock();   
      }
      break;
   
    case MessageTypes::ViewWindowEditLight:
      {
        ViewSceneMessage* rmsg=(ViewSceneMessage*)msg;
        view_window_lock_.lock();
        for(size_t i=0;i<view_window_.size();i++)
        {
          ViewWindow* r=view_window_[i];
          if(r->id_ == rmsg->rid)
          {
            ((lighting_.lights)[rmsg->index])->on = rmsg->on;
            if( rmsg->on )
            {
              if(DirectionalLight *dl = dynamic_cast<DirectionalLight *>
                 (((lighting_.lights)[rmsg->index]).get_rep())) 
              {
                dl->move( rmsg->dir );
                dl->setColor( rmsg->color );
              } 
              else if( HeadLight *hl = dynamic_cast<HeadLight *>
                   (((lighting_.lights)[rmsg->index]).get_rep())) 
              {
                hl->setColor( rmsg->color );
              }
            }
            r->need_redraw_ = 1;
            break;
          }
        }
        view_window_lock_.unlock();
      }
      break;
    case MessageTypes::ViewWindowUpdateClipFrame:
      {
        ViewSceneMessage* rmsg =(ViewSceneMessage*)msg;
        view_window_lock_.lock();

        for(size_t i=0;i<view_window_.size();i++)
        {
          ViewWindow* r=view_window_[i];
          if(r->id_ == rmsg->rid)
          {
            double scale = rmsg->framerate; //framerate used to store scale
            double width = rmsg->tbeg; // tbeg stores width
            double height = rmsg->tend; // tend stores height
            r->viewwindow_clip_frames_[rmsg->index]->Set(rmsg->loc, rmsg->dir,
                                               width, height, scale);
            break;
          }
        }
        view_window_lock_.unlock();
      }
      break;
    case MessageTypes::ViewWindowDumpImage:
      {
        ViewSceneMessage* rmsg=(ViewSceneMessage*)msg;
        view_window_lock_.lock();
        
        for(size_t i=0;i<view_window_.size();i++)
        {
          ViewWindow* r=view_window_[i];
          if(r->id_ == rmsg->rid)
          {
            view_window_lock_.unlock();
            r->renderer_->saveImage(rmsg->filename, rmsg->format, 
                  rmsg->resx,rmsg->resy);
            view_window_lock_.lock();
            break;
          }
        }
        view_window_lock_.unlock();
      }
      break;

    case MessageTypes::ViewWindowMouse:
      {
        ViewWindowMouseMessage* rmsg=(ViewWindowMouseMessage*)msg;
        float NX = 1.0, NY = 1.0;
        bool tracking = false;
        for(size_t i=0;i<view_window_.size();i++)
        {
          ViewWindow* r=view_window_[i];
          if(r->id_ == rmsg->rid)
          {
            (r->*(rmsg->handler))(rmsg->action, rmsg->x, rmsg->y, 
                rmsg->state, rmsg->btn, rmsg->time);
            
            r->gui_lock_view_window_.reset();
            if ( r->gui_lock_view_window_.get())
            {
              tracking = true;
              r->NormalizeMouseXY(rmsg->x, rmsg->y, &NX, &NY);
            }
          }
        }
        if (tracking) 
        {
          for(unsigned int i=0;i<view_window_.size();i++)
          {
            ViewWindow* r=view_window_[i];
            r->gui_lock_view_window_.reset();
            if(r->id_ != rmsg->rid && r->gui_lock_view_window_.get())
            {
              int xx,yy;
              r->UnNormalizeMouseXY(NX,NY,&xx,&yy);
              (r->*(rmsg->handler))(rmsg->action, xx, yy, 
                  rmsg->state, rmsg->btn, rmsg->time);
            }
          }
        }
      }
      break;

    case MessageTypes::GeometryAddLight:
      {
        // Add a light to the light list
        lighting_.lights.add(gmsg->light);

        // Now associate the port, and LightID to the index
        std::map<LightID, int> li;
        std::map<int, std::map<LightID, int> >::iterator it = 
        pli_.find( gmsg->portno );
        if( it == pli_.end() )
        {
          li[gmsg->lserial] =  lighting_.lights.size() - 1;
          pli_[ gmsg->portno ] = li;
        } 
        else 
        {
          ((*it).second)[gmsg->lserial] = lighting_.lights.size() - 1;
        }
      }
      break;

    case MessageTypes::GeometryDelLight:
      {
        std::map<LightID, int>::iterator li;
        std::map<int, std::map<LightID, int> >::iterator it = pli_.find( gmsg->portno );
        if( it == pli_.end() )
        {
          error("Error while deleting a light: no data base for port number " +
          to_string(gmsg->portno));
        } 
        else 
        {
          li = ((*it).second).find(gmsg->lserial);
          if( li == (*it).second.end() )
          {
            error("Error while deleting a light: no light with get_id() " +
            to_string(gmsg->lserial) + "in database for port number" +
            to_string(gmsg->portno));
          } 
          else 
          {
            int idx = (*li).second;
            int i;
            for(i = 0; i < lighting_.lights.size(); i++)
            {
              if( i == idx )
              {
                lighting_.lights[i] = 0;
                lighting_.lights.remove(i);
                ((*it).second).erase( li );
                break;
              }
              if( i == lighting_.lights.size() )
              {
                error("Error deleting light, light not in database.(lserial=" +
                  to_string(gmsg->lserial));
              }
            }
          }
        }
      }
      break;
    
    case MessageTypes::GeometryInit:
      geomlock_.writeLock();
      initPort(gmsg->reply);
      geomlock_.writeUnlock();
      break;

    case MessageTypes::GeometryDetach:
      geomlock_.writeLock();
      detachPort(gmsg->portno);
      geomlock_.writeUnlock();
      break;

    case MessageTypes::GeometryAddObj:
    case MessageTypes::GeometryDelObj:
    case MessageTypes::GeometryDelAll:
      append_port_msg(gmsg);
      msg=0; // Don't delete it yet.
      break;

    case MessageTypes::GeometryFlush:
      geomlock_.writeLock();
      flushPort(gmsg->portno);
      geomlock_.writeUnlock();
      flushViews();
      for(size_t i=0;i<view_window_.size();i++)
        view_window_[i]->redraw_if_needed();
      break;

    case MessageTypes::GeometrySynchronize:
      // Port finish message.  Synchronize viewer.
      //geomlock_.writeLock();
      finishPort(gmsg->portno);
      //geomlock_.writeUnlock();
      break;

    case MessageTypes::GeometryFlushViews:
      geomlock_.writeLock();
      flushPort(gmsg->portno);
      geomlock_.writeUnlock();
      flushViews();
      if(gmsg->wait)
      {
        // Synchronized redraw - do it now and signal them.
        for(size_t i=0;i<view_window_.size();i++)
          view_window_[i]->redraw_if_needed();
        gmsg->wait->up();
      }
      break;

    case MessageTypes::GeometryGetNViewWindows:
      gmsg->nreply->send(view_window_.size());
      break;

    case MessageTypes::GeometryGetData:
      if((unsigned int)(gmsg->which_viewwindow) >= view_window_.size())
      {
        gmsg->datareply->send(0);
      } 
      else 
      {
        view_window_[gmsg->which_viewwindow]->
          getData(gmsg->datamask, gmsg->datareply);
      }
      break;

    case MessageTypes::GeometrySetView:
      if((unsigned int)(gmsg->which_viewwindow) < view_window_.size())
      {
        view_window_[gmsg->which_viewwindow]->setView(gmsg->view);
      }
      break;


    default:
      error("Illegal Message type: " + to_string(msg->type));
      break;
  }

  if(msg) { delete msg; }

  return 1;
}

//----------------------------------------------------------------------
void
ViewScene::initPort(Mailbox<GeomReply>* reply)
{
  int portid=max_portno_++;
  portno_map_.push_back(portid);
  synchronized_map_.push_back(0);
  ports_.addObj(new GeomViewerPort(portid), portid);   // Create the port
  reply->send(GeomReply(portid));
}

//----------------------------------------------------------------------
int
ViewScene::real_portno(int portno)
{
  for (unsigned int i=0; i < portno_map_.size(); i++)
  {
    if (portno == portno_map_[i])
    {
      return i + 1;
    }
  }
  ASSERTFAIL("PORTNO NOT FOUND");

  return (0);
}

//----------------------------------------------------------------------
void
ViewScene::delete_patch_portnos(int portid)
{
  int found = -1;
  for (unsigned int i=0; i < portno_map_.size(); i++)
  {
    if (found >= 0)
    {
      GeomViewerPort *pi;
      if(!(pi = ((GeomViewerPort*)ports_.getObj(portno_map_[i]).get_rep())))
      {
        warning("Geometry message sent to bad port!");
        continue;
      }
      GeomIndexedGroup::IterIntGeomObj iter = pi->getIter();
      for (; iter.first != iter.second; iter.first++)
      {
        GeomViewerItem* si = 
          dynamic_cast<GeomViewerItem*>((*iter.first).second.get_rep());
        if (si)
        {
          const std::string::size_type loc = si->getString().find_last_of('(');
          std::string newname =
            si->getString().substr(0, loc+1) + to_string(i) + ")";
          std::string cached_name = si->getString();
          // Do a rename here.
          
          view_window_lock_.lock();
          for (size_t j = 0; j < view_window_.size(); j++)
          {
            // itemRenamed will set si->name_ = newname
            // so we need to reset it for other windows
            si->getString() = cached_name;
            view_window_[j]->itemRenamed(si, newname);
          }
          view_window_lock_.unlock();
          si->getString() = newname;
        }
      }
    }
    else if (portid == portno_map_[i])
    {
      found = i;
    }
  }

  if (found >= 0)
  {
    portno_map_.erase(portno_map_.begin() + found);
    synchronized_map_.erase(synchronized_map_.begin() + found);
  }

}

//----------------------------------------------------------------------
void
ViewScene::detachPort(int portid)
{
  GeomViewerPort* pi;
  if(!(pi = ((GeomViewerPort*)ports_.getObj(portid).get_rep())))
  {
    warning("Geometry message sent to bad port!");
    return;
  }
  delAll(pi);
  ports_.delObj(portid);
  delete_patch_portnos(portid);
  flushViews();
}

//----------------------------------------------------------------------
void
ViewScene::flushViews()
{
  for (size_t i=0; i<view_window_.size(); i++)
  {
    view_window_[i]->need_redraw_ = 1;
  }
}


//----------------------------------------------------------------------
void
ViewScene::addObj(GeomViewerPort* port, int serial, GeomHandle obj,
		    const std::string& name, CrowdMonitor* lock)
{
  std::string pname(name + " ("+to_string(real_portno(port->portno))+")");
  GeomViewerItem* si = new GeomViewerItem(obj, pname, lock);
  port->addObj(si,serial);
  for (size_t i=0; i<view_window_.size(); i++)
  {
    view_window_[i]->itemAdded(si);
  }
}

//----------------------------------------------------------------------
void
ViewScene::delObj(GeomViewerPort* port, int serial)
{
  GeomViewerItem* si;
  if ((si = ((GeomViewerItem*)port->getObj(serial).get_rep())))
  {
    for (size_t i=0; i<view_window_.size(); i++)
    {
      view_window_[i]->itemDeleted(si);
    }
    port->delObj(serial);
  }
  else
  {
    error("Error deleting object, object not in database.(serial=" +
	  to_string(serial) + ")" );
  }
}

//----------------------------------------------------------------------
void
ViewScene::delAll(GeomViewerPort* port)
{
  GeomIndexedGroup::IterIntGeomObj iter = port->getIter();
  if (!stop_rendering_) 
  {
    for ( ; iter.first != iter.second; iter.first++)
    {
      GeomViewerItem* si = (GeomViewerItem*)((*iter.first).second.get_rep());
      for (unsigned int i=0; i<view_window_.size(); i++)
      {
        view_window_[i]->itemDeleted(si);
      }
    }
  }
  port->delAll();
}

/// Function object which matches a ViewWindow on its (string) ID.
struct vw_id_match :
    public std::binary_function<ViewWindow*, std::string, bool>
{
  bool operator()(const ViewWindow *vw, const std::string &s) const 
  {
    return vw->id_ == s;
  }
};

//----------------------------------------------------------------------
void
ViewScene::delete_viewwindow(const std::string &id)
{
  std::vector<ViewWindow*>::iterator iter;
  Guard vw_lock(&view_window_lock_);

  iter = std::find_if(view_window_.begin(), view_window_.end(),
                      std::bind2nd(vw_id_match(), id));

  if(iter == view_window_.end()) 
  {
    std::cerr << "ERROR in delete_viewwindow, cannot find ID: " << id
              << std::endl;
    return;
  }
  delete *iter;
  view_window_.erase(iter);
}

//----------------------------------------------------------------------
void
ViewScene::tcl_command(GuiArgs& args, void* userdata)
{
  if(args.count() < 2)
  {
    args.error("ViewScene needs a minor command");
    return;
  }
  if(args[1] == "addviewwindow")
  {
    if(args.count() != 3)
    {
      args.error(args[1]+" must have a RID");
      return;
    }
    ViewWindow* r=new ViewWindow(this,  
				    TCLInterface::createContext(args[2]));

    view_window_lock_.lock();
    view_window_.push_back(r);
    view_window_lock_.unlock();

  } 
  else if (args[1] == "deleteviewwindow") 
  {
    if(args.count() != 3)
    {
      args.error(args[1]+" must have a RID");
      return;
    }

    ViewSceneMessage *msg = new ViewSceneMessage
      (MessageTypes::ViewWindowKill,args[2],"");
    mailbox_.send(msg);

  } 
  else 
  {
    Module::tcl_command(args, userdata);
  }

}
//----------------------------------------------------------------------
void
ViewScene::execute()
{
  // Never gets called.
  ASSERTFAIL("ViewScene::execute() should not ever be called.");
}

//----------------------------------------------------------------------
ViewSceneMessage::ViewSceneMessage(const std::string& rid)
  : MessageBase(MessageTypes::ViewWindowRedraw), rid(rid), nframes(0)
{
}

//----------------------------------------------------------------------
ViewSceneMessage::ViewSceneMessage(const std::string& rid, 
           double tbeg, double tend,
			     int nframes, double framerate)
  : MessageBase(MessageTypes::ViewWindowRedraw),
    rid(rid),
    tbeg(tbeg),
    tend(tend),
    nframes(nframes),
    framerate(framerate)
{
  if (nframes <= 0)
  {
    nframes = 1;
  }
}

//----------------------------------------------------------------------
ViewSceneMessage::ViewSceneMessage(MessageTypes::MessageType type,
			     const std::string& rid, const std::string& filename)
  : MessageBase(type), rid(rid), filename(filename)
{
}

//----------------------------------------------------------------------
ViewSceneMessage::ViewSceneMessage(MessageTypes::MessageType type,
			     const std::string& rid,
			     const std::string& filename,
			     const std::string& format,
			     const std::string& resx_string,
			     const std::string& resy_string)
  : MessageBase(type), rid(rid), filename(filename), format(format)
{
  from_string(resx_string,resx);
  from_string(resy_string,resy);
}

ViewSceneMessage::ViewSceneMessage(MessageTypes::MessageType type,
			     const std::string& rid, int lightNo, 
			     bool on, const Vector& lightDir,
			     const Color& lightColor)
  : MessageBase(type), rid(rid), dir(lightDir), color(lightColor), 
    index(lightNo), on(on)
{}


//----------------------------------------------------------------------
ViewSceneMessage::ViewSceneMessage(MessageTypes::MessageType type,
                             const std::string& rid, int clip_no,
                             const Point& center, const Vector& normal,
                             double width, double height, double scale)
  : MessageBase(type), rid(rid),
    tbeg(width), // tbeg stores width
    tend(height), // tend store height
    framerate( scale ), //use framerate to store scale
    loc(center), dir( normal ), 
    index( clip_no )
{}
//----------------------------------------------------------------------

ViewSceneMessage::~ViewSceneMessage()
{
}


//----------------------------------------------------------------------
void
ViewScene::append_port_msg(GeometryComm* gmsg)
{
  // Look up the right port number.
  GeomViewerPort *pi;
  if (!(pi = ((GeomViewerPort*)ports_.getObj(gmsg->portno).get_rep())))
  {
    warning("Geometry message sent to bad port!!!: "+to_string(gmsg->portno));
    return;
  }
  
  // Queue up the messages until the flush.
  if(pi->msg_tail)
  {
    pi->msg_tail->next=gmsg;
    pi->msg_tail=gmsg;
  } 
  else
  {
    pi->msg_head=pi->msg_tail=gmsg;
  }
  gmsg->next=0;
}

//----------------------------------------------------------------------
void
ViewScene::flushPort(int portid)
{
  // Look up the right port number.
  GeomViewerPort* pi;
  if(!(pi = ((GeomViewerPort*)ports_.getObj(portid).get_rep())))
  {
    warning("Geometry message sent to bad port!");
    return;
  }
  
  GeometryComm* gmsg=pi->msg_head;
  while(gmsg)
  {
    switch(gmsg->type)
    {
    case MessageTypes::GeometryAddObj:
      addObj(pi, gmsg->serial, gmsg->obj, gmsg->name, gmsg->lock);
      break;
    case MessageTypes::GeometryDelObj:
      delObj(pi, gmsg->serial);
      break;
    case MessageTypes::GeometryDelAll:
      delAll(pi);
      break;
    default:
      error("How did this message get in here???");
    }
    GeometryComm* next=gmsg->next;
    delete gmsg;
    gmsg=next;
  }
  if(pi->msg_head)
  {
    flushViews();
    pi->msg_head=pi->msg_tail=0;
  }
}


//----------------------------------------------------------------------
void
ViewScene::finishPort(int portid)
{
  synchronized_map_[real_portno(portid)-1]++;

  bool all = true;
  for (size_t i=0; i < synchronized_map_.size(); i++)
  {
    if (synchronized_map_[i] == 0)
    {
      all = false;
      break;
    }
  }
  
  if (all)
  {
    flushViews();
    for(size_t i=0;i<view_window_.size();i++)
      view_window_[i]->redraw_if_needed();


    // Clear the entries from the map.
    for (unsigned int i = 0; i < synchronized_map_.size(); i++)
    {
      synchronized_map_[i]--;
    }

    // Push the serial number back to the scheduler.
    if (synchronized_serials_.size())
    {
      unsigned int serial = synchronized_serials_.front();
      synchronized_serials_.pop_front();

      sched_->report_execution_finished(serial);

      // This turns on synchronous movie making.  It's very useful for
      // making movies that are driven by an event loop (such as
      // send-intermediate).  The movie frames are not taken during
      // user interaction but only on execute.  This is only
      // synchronous in one direction.  If a module executes too fast
      // then the latest one is used.  That is, it waits for the slow
      // modules but does not throttle the fast ones.
      for(size_t window = 0; window < view_window_.size(); window++) 
      {
        view_window_[window]->maybeSaveMovieFrame();
      }

    }
    else
    {
      // Serial number hasn't arrived yet, defer until we get it.
      synchronized_debt_--;
    }
  }
}


void
ViewScene::set_context(Network* network)
{
  if (sci_getenv("SCI_REGRESSION_TESTING"))
  {
    network->get_scheduler()->add_callback(regression_save_image_callback, this, -1);
  }
  if (sci_getenv("SCIRUN_SAVE_IMAGE_FILE"))
  {
    network->get_scheduler()->add_callback(save_image_callback, this, -1);  
  }

  Module::set_context(network);
}

bool
ViewScene::regression_save_image_callback(void *voidstuff)
{
  ViewScene *viewer = (ViewScene *)voidstuff;
  for (unsigned int i = 0; i < viewer->view_window_.size(); i++)
  {
    std::string name = std::string("snapshot") + to_string(i) + ".png";
    
    if (sci_getenv("SCIRUN_NETFILE") != std::string(""))
    {
      std::string netfile = sci_getenv("SCIRUN_NETFILE");
      std::string::size_type lastslash = netfile.find_last_of("/");
      std::string::size_type lastdot   = netfile.find_last_of(".");
      
      if (lastslash == std::string::npos)
      {
        name = netfile.substr(0,lastdot) + "-window" + to_string(i) + ".png";
      }
      else
      {          
        if (lastdot > lastslash)
        {
          name = netfile.substr(lastslash+1,lastdot-lastslash-1) + "-window" + to_string(i) + ".png";
        }
      }
      
      if (sci_getenv_p("SCI_REGRESSION_IMAGE_DIR"))
      {
        std::string path = sci_getenv("SCI_REGRESSION_IMAGE_DIR");
        if (path[path.size()-1] != '/') path += '/';
        name = path + name;
      }
    }

    std::cerr << "saving "<<name << std::endl;
    
    viewer->view_window_[i]->redraw_if_needed();
    // Make sure that the 640x480 here matches up with ViewWindow.cc defaults.
    viewer->view_window_[i]->renderer_->saveImage(name, "png", 640, 480);
    viewer->view_window_[i]->redraw(); // flushes saveImage.

    std::cerr << "saved "<<name << std::endl;


  }
  return true;
}



bool
ViewScene::save_image_callback(void *voidstuff)
{
  ViewScene *viewer = (ViewScene *)voidstuff;
  for (size_t i = 0; i < viewer->view_window_.size(); i++)
  {
    std::string name = std::string("snapshot") + to_string(i) + ".png";
    
    std::string idstr = viewer->get_id();
    size_t uloc = idstr.find_last_of("_");
    std::string viewscenenum = idstr.substr(uloc+1);
    
    if (sci_getenv("SCIRUN_SAVE_IMAGE_FILE") != std::string(""))
    {
      std::string netfile = sci_getenv("SCIRUN_SAVE_IMAGE_FILE");
      std::string::size_type lastslash = netfile.find_last_of("/");
      std::string::size_type lastdot   = netfile.find_last_of(".");
      if (lastdot == std::string::npos)
      {
        name = netfile + "-window" + viewscenenum +"-" +to_string(i) + ".png";      
      }
      else
      {      
        if (lastslash == std::string::npos)
        {
          name = netfile.substr(0,lastdot) + "-window" + viewscenenum +"-" +to_string(i) + ".png";                
        }
        else
        {          
          if (lastdot > lastslash)
          {
            name = netfile.substr(lastslash+1,lastdot) + "-window"  + viewscenenum +"-" + to_string(i) + ".png";
          }
          else
          {
            name = netfile.substr(lastslash+1) + "-window"  + viewscenenum +"-" + to_string(i) + ".png";        
          }
        }
      }
    }

    
    viewer->view_window_[i]->redraw_if_needed();
    // Make sure that the 640x480 here matches up with ViewWindow.cc defaults.

    std::cerr << "Saving "<<name<<"..." << std::endl;
    viewer->view_window_[i]->renderer_->saveImage(name, "png");
    viewer->view_window_[i]->redraw(); // flushes saveImage.

  }
  return true;
}


} // End namespace SCIRun
