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
//    File   : OpenGLBase.h
//    Author : Martin Cole
//    Date   : Tue Sep 18, 2007
//
//! Simple base class for any threaded opengl window.


#if !defined(OpenGLBase_h)
#define OpenGLBase_h

#include <Core/Events/EventManager.h>
#include <Core/Events/ToolManager.h>
#include <Core/Thread/Runnable.h>
#include <Core/Thread/Thread.h>
#include <Core/Util/Timer.h>
#include <Core/Datatypes/Color.h>
#include <Core/Events/share.h>


namespace SCIRun {

class OpenGLContext;

class SCISHARE OpenGLBase : public Runnable
{
public:

  //! The id is the name the mailbox gets. setting event targets to the
  //! name constructed here, uniquely sends events to this object.
  OpenGLBase(OpenGLContext*, const std::string& id = "OpenGLBase");
  virtual ~OpenGLBase();

  // Calls update on a timer. This is the top level rendering loop.
  virtual void          run();
  // process events, then call redraw_frame.
  virtual void          update();

  ToolManager&        get_tm() { return tm_; }
  const Color&        bgcolor() { return bgcolor_; }

  virtual int         width() const;
  virtual int         height() const;
  virtual void        set_timer_interval(int i) { interval_ = i; }
  virtual int         get_timer_interval() const { return interval_; }

  //! the following enum defines the stack priorities for the tool manager,
  //! and classifies the priority ranges. For now leave room for 100 stacks
  //! in each range.  See init_tool_manager() for use of these ranges.
  enum {
    EVENT_MODIFIERS_E = 0, // Tools to modify incoming events.
    TOOL_MODIFIERS_E = 100, // Tools to manipulate the set of active tools.
    DATA_TOOLS_E = 200, // Tools that handle data,
    SELECTION_TOOL_E = 298,  // Tool that gets pushed used and popped...
    ACTIVE_TOOL_E = 299,  // Tool that gets pushed used and popped...
    VIEWER_TOOLS_E = 300,  // Tools to manipulate the current view.
                           // always on the stack (so last)
    LAST_CHANCE_E = 500
  };

  std::string id() { return id_; }
  void   set_dead() { dead_ = true; } 

protected:
  // Main draw function, called after events have been handled.
  virtual void                     redraw_frame() = 0;

  std::string                      id_;
  int                              xres_;
  int                              yres_;
  int                              interval_;
  OpenGLContext                   *gl_context_;
  Color                            bgcolor_;
  bool                             dead_;
  ToolManager                      tm_;
  EventManager::event_mailbox_t   *events_;
};

} // End namespace SCIRun

#endif // OpenGLBase_h
