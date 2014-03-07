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
 *  SynchronizeGeometry.cc:
 *
 *  Written by:
 *   Kai Li
 *   Jan, 2003
 *
 */
#include <Dataflow/Network/Module.h>
#include <Dataflow/Network/Scheduler.h>

#include <Dataflow/Network/Ports/GeometryPort.h>
#include <Dataflow/Network/Ports/GeometryComm.h>
#include <Dataflow/Comm/MessageTypes.h>
#include <Dataflow/Comm/MessageBase.h>
#include <Core/Util/StringUtil.h>

#include <vector>

namespace SCIRun {

class SynchronizeGeometry : public Module {
public:
  SynchronizeGeometry(GuiContext*);

  virtual ~SynchronizeGeometry();

  virtual void execute();
  virtual void do_execute();

private:
  std::vector<GeometryComm*> msg_heads_;
  std::vector<GeometryComm*> msg_tails_;
  std::vector<int>           id_;
  std::vector<int>           synchronized_map_;  
  std::vector<std::map<GeomID, GeomID, std::less<GeomID> > > geom_ids_;

  int max_id_;

  GeometryOPortHandle ogeom_;

  GuiInt gui_enforce_;

  int process_event(MessageBase* message);
  void forward_saved_msg();
  void flush_all_msgs();
  void append_msg(GeometryComm* gmsg);
  bool flush_port(int portno, int count);

  virtual void widget_moved(bool, BaseWidget*) { std::cout << "SynchronizeGeometry::widget_moved called!\n"; }
};



DECLARE_MAKER(SynchronizeGeometry)
SynchronizeGeometry::SynchronizeGeometry(GuiContext* ctx)
  : Module("SynchronizeGeometry", ctx, Filter, "Render", "SCIRun"),
    max_id_(0),
    gui_enforce_(get_ctx()->subVar("enforce"), 1)
{
  //! make sure do_execute is not done.
  have_own_dispatch_ = true;
}



SynchronizeGeometry::~SynchronizeGeometry()
{
}



void
SynchronizeGeometry::execute()
{
}



void
SynchronizeGeometry::do_execute()
{
  get_oport_handle("Output Geometry",ogeom_);
  for (;;)
  {
    MessageBase *msg = mailbox_.receive();
    if (process_event(msg) == 86)
    {
      return;
    }
  }
}



int
SynchronizeGeometry::process_event(MessageBase* msg)
{
  GeometryComm* gmsg = (GeometryComm*)msg;

  switch (msg->type)
  {
  case MessageTypes::GoAway:
    return 86;

  case MessageTypes::GeometryInit:
    {
      max_id_++;
      gmsg->reply->send(GeomReply(max_id_));
      
      id_.push_back(max_id_);
      synchronized_map_.push_back(0);
      msg_heads_.push_back(0);
      msg_tails_.push_back(0);
      geom_ids_.push_back(std::map<GeomID, GeomID>());
    }
    break;

  case MessageTypes::GeometryDetach:
    {
      // Remove all of the gmsg->portno objects
      
      int idx = -1;
      for (size_t p = 0; p < id_.size(); p++) 
        if (id_[p] == gmsg->portno) { idx = p; break; } 
      if (idx == -1) break;
      

      std::map<GeomID, GeomID>::iterator itr;
      
      itr = geom_ids_[idx].begin();
      int counter = 0;
      while (itr != geom_ids_[idx].end())
      {
        ogeom_->delObj((*itr).second);
        ++itr;
        counter++;
      }
      geom_ids_[idx].clear();      
      if (counter) { ogeom_->flush(); }

      // Maybe still connected geometries are now in sync?
      gui_enforce_.reset();
      if (gui_enforce_.get())
      { 
        forward_saved_msg();
      }
      else
      {
        flush_all_msgs();
      }

      // Fix the portnos.
      id_.erase(id_.begin()+idx);
      synchronized_map_.erase(synchronized_map_.begin()+idx);
      msg_heads_.erase(msg_heads_.begin()+idx);
      msg_tails_.erase(msg_tails_.begin()+idx);
      geom_ids_.erase(geom_ids_.begin()+idx);
      

      // Check whether we were waiting for this one
      bool all = true;
      for (unsigned int i=0; i < synchronized_map_.size(); i++)
      {
        if (synchronized_map_[i] < 0) continue;
        if (synchronized_map_[i] == 0)
        {
          all = false;
          break;
        }
      }
    
      if (all)
      {
        for (unsigned int i = 0; i < synchronized_map_.size(); i++)
        {
          if (synchronized_map_[i] < 0) continue;
          synchronized_map_[i]--;
        }
        ogeom_->synchronize();
      }          
    }
    break;

  case MessageTypes::GeometrySynchronize:  
    {
      int idx = -1;
      for (size_t p = 0; p < id_.size(); p++) 
        if (id_[p] == gmsg->portno) { idx = static_cast<int>(p); break; } 
      if (idx == -1) break;
      
      synchronized_map_[idx]++;
    
      bool all = true;
      for (unsigned int i=0; i < synchronized_map_.size(); i++)
      {
        if (synchronized_map_[i] < 0) continue;
        if (synchronized_map_[i] == 0)
        {
          all = false;
          break;
        }
      }
    
      if (all)
      {
        for (unsigned int i = 0; i < synchronized_map_.size(); i++)
        {
          if (synchronized_map_[i] < 0) continue;
          synchronized_map_[i]--;
        }

        ogeom_->synchronize();
      }
    }
    break;
    
  case MessageTypes::GeometryDelAll:
    {
      int idx = -1;
      for (size_t p = 0; p < id_.size(); p++) 
        if (id_[p] == gmsg->portno) { idx = static_cast<int>(p); break; } 
      if (idx == -1) break;
    
      gui_enforce_.reset();
      if (gui_enforce_.get())
      {
        append_msg(gmsg);
        msg = 0;
      }
      else
      {
        flush_all_msgs();
        std::map<GeomID, GeomID>::iterator itr;
        itr = geom_ids_[idx].begin();
        while (itr != geom_ids_[idx].end())
        {
          ogeom_->delObj((*itr).second);
          ++itr;
        }
        geom_ids_[idx].clear();
      }
    }
    break;

  case MessageTypes::GeometryDelObj:
    {
      int idx = -1;
      for (size_t p = 0; p < id_.size(); p++) 
        if (id_[p] == gmsg->portno) { idx = static_cast<int>(p); break; } 
      if (idx == -1) break;

      gui_enforce_.reset();
      if (gui_enforce_.get())
      {
        append_msg(gmsg);
        msg = 0;
      }
      else
      {
        flush_all_msgs();
        // TODO: verify get_id() found in map.
        ogeom_->delObj(geom_ids_[idx][gmsg->serial]);
        geom_ids_[idx].erase(gmsg->serial);
      }
    }
    break;

  case MessageTypes::GeometryAddObj:
    {
      int idx = -1;
      for (size_t p = 0; p < id_.size(); p++) 
        if (id_[p] == gmsg->portno) { idx = static_cast<int>(p); break; } 
      if (idx == -1) break;
    
      gui_enforce_.reset();
      if (gui_enforce_.get())
      {
        append_msg(gmsg);
        msg = 0;
      }
      else
      {
        flush_all_msgs();
        const std::string pnum = to_string(idx);
        const std::string newname =  gmsg->name + " (" + pnum + ")";
        const GeomID id = ogeom_->addObj(gmsg->obj, newname, gmsg->lock);
        geom_ids_[idx][gmsg->serial] = id;
      }
    }
    break;

  case MessageTypes::GeometryFlush:
    {
      gui_enforce_.reset();
      if (gui_enforce_.get())
      {
        append_msg(gmsg);
        forward_saved_msg();
        msg = 0;
      }
      else
      {
        flush_all_msgs();
        ogeom_->flush();
      }
    }
    break;

  case MessageTypes::GeometryFlushViews:
    {
      gui_enforce_.reset();
      if (gui_enforce_.get())
      {
        append_msg(gmsg);
        forward_saved_msg();
        msg = 0;
      }
      else
      {
        flush_all_msgs();
        ogeom_->flushViews();
      }
    }
    break;

  case MessageTypes::ExecuteModule:
    {
      //! only forward messages if we were waiting for them
      gui_enforce_.reset();
      if (!gui_enforce_.get())
      {
        flush_all_msgs();
      }
      sched_->report_execution_finished(msg);
    }
    break;
  
  case MessageTypes::SynchronizeModule:
    // We (mostly) ignore these messages.
    sched_->report_execution_finished(msg);
    break;

  default:
    break;
  }

  if (msg)
  {
    delete msg;
  }

  return 0;
}



void
SynchronizeGeometry::append_msg(GeometryComm* gmsg)
{
  int idx = -1;
  for (size_t p = 0; p < id_.size(); p++) 
    if (id_[p] == gmsg->portno) { idx = static_cast<int>(p); break; } 
  if (idx == -1) return;

  gmsg->next = 0;
  if (msg_heads_[idx])
  {
    msg_tails_[idx]->next = gmsg;
    msg_tails_[idx] = gmsg;
  }
  else
  {
    msg_heads_[idx] = msg_tails_[idx] = gmsg;
  }
}


      
void
SynchronizeGeometry::forward_saved_msg()
{
  std::ostringstream str;
  str << " Checking " << id_.size() << " ports.";
  remark( str.str() );

  int num_flush, valid;

  num_flush = 0;
  valid = 0;
  
  for (size_t i = 0; i < id_.size() ; i++)
  {
    valid++;
    GeometryComm *tmp_gmsg = msg_heads_[i];
    
    while (tmp_gmsg)
    {
      if (tmp_gmsg->type == MessageTypes::GeometryFlush ||
          tmp_gmsg->type == MessageTypes::GeometryFlushViews)
      {
        num_flush++;

        std::ostringstream str;
        str << "  port " << i << " is ready.";
        remark( str.str() );
        break;
      }
      tmp_gmsg = tmp_gmsg->next;
    }
  }

  if (num_flush == valid)
  {
    remark( " All were ready, flushing." );
    bool some = false;
    for (size_t i = 0; i < id_.size(); i++)
    {
      some |= flush_port(static_cast<int>(i), 1);
    }
    
    if (some) { ogeom_->flush(); }

    update_progress(1.0);
    update_state(Completed);
  }
  else
  {
    update_progress(num_flush/double(id_.size()));
  }
}



bool
SynchronizeGeometry::flush_port(int portno, int count)
{
  GeometryComm *gmsg = msg_heads_[portno];
  const bool some = gmsg;
 
  while (gmsg)
  {
    int idx = -1;
    for (size_t p = 0; p < id_.size(); p++) 
      if (id_[p] == gmsg->portno) { idx = static_cast<int>(p); break; } 
    if (idx == -1) break;  
  
    // Process messages here.
    // GeometryDelAll
    if (gmsg->type == MessageTypes::GeometryDelAll)
    {
      std::map<GeomID, GeomID>::iterator itr = geom_ids_[idx].begin();
      while (itr != geom_ids_[idx].end())
      {
        ogeom_->delObj((*itr).second);
        ++itr;
      }
      geom_ids_[idx].clear();
    }

    // GeometryDelObj
    else if (gmsg->type == MessageTypes::GeometryDelObj)
    {
      // TODO: verify get_id() found in map.
      ogeom_->delObj(geom_ids_[idx][gmsg->serial]);
      geom_ids_[idx].erase(gmsg->serial);
    }

    // GeometryAddObj
    else if (gmsg->type == MessageTypes::GeometryAddObj)
    {
      const std::string pnum = to_string(idx);
      const std::string newname =  gmsg->name + " (" + pnum + ")";
      const GeomID id = ogeom_->addObj(gmsg->obj, newname);
      geom_ids_[idx][gmsg->serial] = id;
    }

    // Eat up the flushes.
    else if (gmsg->type == MessageTypes::GeometryFlush ||
	     gmsg->type == MessageTypes::GeometryFlushViews)
    {
      count--;
      if (count == 0)
      {
        msg_heads_[portno] = gmsg->next;
        if (gmsg->next == 0)
        {
          msg_tails_[portno] = 0;
        }
        delete gmsg;
        break;
      }
    }
    else
    {
      // Unprocessed message.
    }

    GeometryComm *next = gmsg->next;
    delete gmsg;
    gmsg = next;
  }

  if (gmsg == 0)
  {
    msg_heads_[portno] = 0;
    msg_tails_[portno] = 0;
  }

  return some;
}



void
SynchronizeGeometry::flush_all_msgs()
{
  bool some = false;
  for (size_t i = 0; i < id_.size(); i++)
  {
    some |= flush_port(static_cast<int>(i), -1);
  }
 
  if (some)
  {
    ogeom_->flush();
  }

  update_progress(1.0);
  update_state(Completed);
}


} // End namespace SCIRun
