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
//    File   : OpenGLBase.cc
//    Author : Martin Cole
//    Date   : Sat May 27 08:51:31 2006

#include <sci_gl.h>
#include <sci_values.h>

#include <sci_defs/bits_defs.h>

#include <Core/Geom/OpenGLContext.h>
#include <Core/Events/OpenGLBase.h>
#include <Core/Events/BaseTool.h>

#ifdef _WIN32
#  include <windows.h>
#  include <winbase.h>
#  include <Core/Thread/Time.h>
#  if !defined(BUILD_SCIRUN_STATIC)
#    define SCISHARE __declspec(dllimport)
#  else
#    define SCISHARE
#  endif
#else
#  define SCISHARE
#  include <sys/time.h>
#endif

namespace SCIRun {

OpenGLBase::OpenGLBase(OpenGLContext *oglc, const std::string& id) :
  id_(id),
  xres_(0),
  yres_(0),
  interval_(30),
  gl_context_(oglc),
  bgcolor_(Color(.0, .0, .0)),
  dead_(false),
  tm_("OpenGLBase tool manager"),
  events_(0)
{
  if (gl_context_) {
    events_ = EventManager::register_event_messages(id);
  }
}

OpenGLBase::~OpenGLBase()
{
  delete gl_context_;
  gl_context_ = 0;
}

void
OpenGLBase::update()
{
  // process the cached events since last draw.
  event_handle_t ev;
  while (events_ && events_->tryReceive(ev))
  {
    // Tools will set up the appropriate rendering state.
    QuitEvent *qe = dynamic_cast<QuitEvent*>(ev.get_rep());
    if (qe) 
    {
      dead_ = true;
      
      // Signal requesting thread that we are done
      qe->signal_done();
      
      // this is the terminate signal, so return.
      return;
    }
    tm_.propagate_event(ev);
  }
  redraw_frame();
  
  // replies should have been sent via events to the EventManager.
}


void
OpenGLBase::run()
{
  TimeThrottle throttle;
  throttle.start();
  // the rate at which we refresh the scene.
  const double inc = 1. / double(interval_);
  double t = throttle.time();
  
  while (!dead_) 
  {
    t = throttle.time();
    throttle.wait_for_time(t + inc);
    update();
  } // end while(!dead_)
  throttle.stop();
}


int
OpenGLBase::width() const
{
  ASSERT(gl_context_);
  return gl_context_->width();
}


int
OpenGLBase::height() const
{
  ASSERT(gl_context_);
  return gl_context_->height();
}


} // End namespace SCIRun
