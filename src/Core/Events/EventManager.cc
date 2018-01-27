//
//  For more information, please see: http://software.sci.utah.edu
//
//  The MIT License
//
//  Copyright (c) 2009 Scientific Computing and Imaging Institute,
//  University of Utah.
//
//
//  Permission is hereby granted, free of charge, to any person obtaining a
//  copy of this software and associated documentation files (the "Software"),
//  to deal in the Software without restriction, including without limitation
//  the rights to use, copy, modify, merge, publish, distribute, sublicense,
//  and/or sell copies of the Software, and to permit persons to whom the
//  Software is furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included
//  in all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.
//
//    File   : EventManager.cc
//    Author : Martin Cole
//    Date   : Thu May 25 15:53:19 2006

#include <Core/Events/EventManager.h>
#include <Core/Util/Environment.h>
#include <Core/Util/Timer.h>
#include <iostream>

#include <sci_defs/x11_defs.h>

#if defined(__APPLE__) && !defined(HAVE_X11)
namespace Carbon {
#  include <Carbon/Carbon.h>
}
#endif


namespace SCIRun {



EventManager::id_tm_map_t EventManager::mboxes_;
Mutex EventManager::mboxes_lock_("EventManager mboxes_ lock");
Mailbox<event_handle_t> EventManager::mailbox_("EventManager", 512);
PiostreamPtr EventManager::trailfile_stream_;


EventManager::EventManager() :
  tm_("EventManager tools")
{
  if (sci_getenv_p("SCI_DEBUG")) {
    std::cerr << "EventManager has been created." << std::endl;
  }
}


EventManager::~EventManager()
{
  mailbox_.send(new QuitEvent());
  mboxes_lock_.lock();
  bool empty = mboxes_.empty();
  mboxes_lock_.unlock();
  if (!empty) {
    TimeThrottle timer;
    timer.start();
    const double start = timer.time();
    double t = start;
    do {
      mboxes_lock_.lock();
      empty = mboxes_.empty();
      mboxes_lock_.unlock();
      t = timer.time();
      if (!empty) {
        timer.wait_for_time(t + 1/20.0);
      }
    } while (!empty && (t < start + 5.0));
    timer.stop();

    if (!empty && sci_getenv_p("SCI_DEBUG")) {
      std::cerr << "EventManager quit after waiting " << t - start
           << " seconds for the following mailboxes to unregister:\n";

      mboxes_lock_.lock();
      for (id_tm_map_t::iterator it = mboxes_.begin();
           it != mboxes_.end(); ++it)
        {
          std::cerr << it->first << "\n";
        }
      mboxes_lock_.unlock();
    }
  }

  if (sci_getenv_p("SCI_DEBUG")) {
    std::cerr << "EventManager has been destroyed." << std::endl;
  }
  stop_trail_file();
}


void
EventManager::add_event(event_handle_t event)
{
  if (event->is_trail_enabled() && trailfile_is_playing()) {
    //    mailbox_.send(new TrailfilePauseEvent());
    //    mailbox_.send(event);
    return;
  }

  mailbox_.send(event);
}


bool
EventManager::open_trail_file(const std::string &filename, bool record)
{
  if (trailfile_stream_) {
    return false;
  }

  if (record) {
    trailfile_stream_ = auto_ostream(filename, "Text", 0);
  } else {
    trailfile_stream_ = auto_istream(filename, 0);
  }


  if (trailfile_stream_ && trailfile_stream_->error()) {
    trailfile_stream_.reset();
    return false;
  }

  if (record) {
    trailfile_stream_->disable_pointer_hashing();
  }

  return trailfile_stream_ ? true : false;
}


void
EventManager::stop_trail_file()
{
  if (!trailfile_stream_) return;
  trailfile_stream_.reset();
}


bool
EventManager::trailfile_is_playing()
{
  return trailfile_stream_ && trailfile_stream_->reading();
}


bool
EventManager::trailfile_is_recording()
{
  return trailfile_stream_ && trailfile_stream_->writing();
}


EventManager::event_mailbox_t*
EventManager::register_event_messages(const std::string& id)
{
  mboxes_lock_.lock();
  id_tm_map_t::iterator it = mboxes_.find(id);

  if (sci_getenv_p("SCI_DEBUG")) {
    std::cerr << "Mailbox id \"" << id << "\" registered\n";
  }
  std::string mbox_name = "event_mailbox_" + id;
  event_mailbox_t* mbox = new event_mailbox_t(mbox_name.c_str(), 8024);
  mboxes_.insert(make_pair(id, mbox));
  mboxes_lock_.unlock();
  return mbox;
}


void
EventManager::unregister_event_messages(const std::string& id)
{
  mboxes_lock_.lock();
  id_tm_map_t::iterator it = mboxes_.find(id);
  if (it == mboxes_.end()) {
    std::cerr << "Warning: maiboxl id """ << id
         << """ is not registered." << std::endl;
    std::cerr << "Valid ids are: ";

    for (id_tm_map_t::iterator it = mboxes_.begin();
         it != mboxes_.end(); ++it)
      std::cerr << it->first << ", ";
    std::cerr << std::endl;
    mboxes_lock_.unlock();
    return;
  }
  if (sci_getenv_p("SCI_DEBUG")) {
    std::cerr << "Mailbox id \"" << id << "\" un-registered\n";
  }
  delete it->second;
  mboxes_.erase(it);
  mboxes_lock_.unlock();
}


void
EventManager::unregister_mailbox(event_mailbox_t *mailbox)
{
  mboxes_lock_.lock();
  for (id_tm_map_t::iterator it = mboxes_.begin(); it != mboxes_.end(); ++it)
  {
    if (it->second == mailbox) {
      if (sci_getenv_p("SCI_DEBUG")) {
        std::cerr << "Mailbox id \"" << it->first << "\" un-registered\n";
      }
      delete it->second;
      mboxes_.erase(it);
      mboxes_lock_.unlock();
      return;
    }
  }

  std::cerr << "Cannot find mailbox to unregister!\n";
  mboxes_lock_.unlock();
}


void
EventManager::play_trail()
{
  ASSERT(trailfile_stream_ && trailfile_stream_->reading());

  event_handle_t event;
  signed long event_time = 0;
  signed long last_event_time = 0;
  double last_timer_time = 0;
  const double millisecond = 1.0 / 1000.0;
  TimeThrottle timer;
  timer.start();
  int l = 0;
  while (trailfile_stream_ && !trailfile_stream_->eof()) {
    ++l;
    event = 0;
    Pio(*trailfile_stream_, event);

    if (!event.get_rep()) {
      stop_trail_file();
      std::cerr << "Error in trailfile on Line:  " << l << std::endl;
      continue;
    }
    if ((event_time = event->get_time())) {

      if (last_event_time) {
        double diff = (event_time-last_event_time) * millisecond;
        if (diff > 10) {
          std::cerr << "Event Time: " << event_time << std::endl;
          std::cerr << "Last Event Time: " << last_event_time << std::endl;
          std::cerr << "diff: " << diff << std::endl;

          diff = 10;
        }
        if (diff < 0) {
          diff = 0;
        }

        timer.wait_for_time(last_timer_time + diff);
      }
      last_event_time = event_time;
      last_timer_time = timer.time();
    }

    mailbox_.send(event);
  }
  timer.stop();
}


void
EventManager::do_trails()
{
  std::string default_trailfile = std::string(sci_getenv("SCIRUN_TMP_DIR"))+"/"+
                           sci_getenv("EXECUTABLE_NAME")+".trail";
  std::string default_trailmode("");

  const char *trailfile = sci_getenv("SCIRUN_TRAIL_FILE");
  const char *trailmode = sci_getenv("SCIRUN_TRAIL_MODE");

  if (!trailfile) {
    trailfile = default_trailfile.c_str();
  }

  if (!trailmode) {
    trailmode = default_trailmode.c_str();
  }

  if (trailmode && trailfile) {
    if (trailmode[0] == 'R') {
      if (EventManager::open_trail_file(trailfile, true)) {
        std::cerr << "Recording trail file: ";
      } else {
        std::cerr << "ERROR recording trail file ";
      }
      std::cerr << trailfile << std::endl;
    }

    if (trailmode[0] == 'P') {
      if (EventManager::open_trail_file(trailfile, false)) {
        std::cerr << "Playing trail file: " << trailfile << std::endl;
        EventManager::play_trail();
        std::cerr << "Trail file completed.\n";
      } else {
        std::cerr << "ERROR playing trail file " << trailfile << std::endl;
      }
    }
  }
}



void
EventManager::run()
{
  bool debug_p = false;
  if (sci_getenv_p("SCI_DEBUG")) 
  {
    debug_p = true;
  }
  bool done = false;
  event_handle_t event;
  do 
  {
    event = tm_.propagate_event(mailbox_.receive());

    static signed long last_event_time = 0;
    static double last_timer_time = 0;
    static TimeThrottle timer;
    if (timer.state() == Timer::Stopped) 
    {
      timer.start();
    }

    if (event->get_time()) 
    {
      last_event_time = event->get_time();
    } 
    else if (last_event_time) 
    {
      last_event_time +=
        (signed long)((timer.time() - last_timer_time) * 100.0);
      event->set_time(last_event_time);
    }
    last_timer_time = timer.time();

    if (event->is_trail_enabled() && trailfile_is_recording()) 
    {
      Pio(*trailfile_stream_, event);
    }

    if (dynamic_cast<QuitEvent*>(event.get_rep()) != 0 &&
        event->get_target().empty()) 
    {
      done = true;
    }

    mboxes_lock_.lock();
    const std::string target = event->get_target();
    if (debug_p) 
    {
      std::cerr << "Event has target id '" << target << "'" << std::endl;
    }
    
    // If the event has a specific target mailbox
    if (!target.empty()) 
    {
      typedef std::pair<id_tm_map_t::iterator, id_tm_map_t::iterator> mbox_range_t;
      mbox_range_t range = mboxes_.equal_range(target);
      if (range.first == range.second) 
      {
        if (debug_p) 
        {
          std::cerr << "Event target mailbox id """ << target
               << """ not found." << std::endl;
        }
      } 
      else 
      {
        while(range.first != range.second)
        {
          if (debug_p) 
          {
            std::cerr << "Event target mailbox id " << target
                 << " found, sending." << std::endl;
            std::cerr << range.first->first << ": mailbox numItems: "
                 << range.first->second->numItems() << std::endl;
          }
          range.first->second->send(event);
          ++range.first;
        }
      }
    } 
    else 
    {
      if (debug_p) 
      {
        std::cerr << "Event target mailbox id empty, broadcasting." << std::endl;
      }

      id_tm_map_t::iterator it = mboxes_.begin();
      for (;it != mboxes_.end(); ++it) 
      {
        if (sci_getenv_p("SCI_DEBUG")) 
        {
          std::cerr << it->first << " size: " << it->second->numItems() << std::endl;
        }
        it->second->send(event);
      }
    }
    mboxes_lock_.unlock();
    if (timer.state() != Timer::Stopped) 
    {
      timer.stop();
    }
  } 
  while (!done);
  
#if defined(__APPLE__) && !defined(HAVE_X11)
  Carbon::QuitApplicationEventLoop();
#endif
}

} // namespace SCIRun
