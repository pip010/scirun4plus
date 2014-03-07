//
//  For more information, please see: http://software.sci.utah.edu
//
//  The MIT License
//
//  Copyright (c) 2009 Scientific Computing and Imaging Institute,
//  University of Utah.
//
//  License for the specific language governing rights and limitations under
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
//    File   : WindowRedraw.cc
//    Author : McKay Davis
//    Date   : Thu Feb 22 12:22:37 2007

#include <Core/Skinner/WindowRedraw.h>
#include <Core/Skinner/Window.h>
#include <Core/Skinner/FilterRedrawEventsTool.h>
#include <Core/Util/Environment.h>
#include <iostream>


using std::cerr;
using std::endl;

namespace SCIRun {
namespace Skinner {


WindowRedraw::WindowRedraw(GLWindow *window, double hertz) :
  ThrottledRunnable(hertz),
  window_(window),
  mailbox_(EventManager::register_event_messages(window_->get_id())),
  redraw_event_(new WindowEvent(WindowEvent::REDRAW_E))
{
}


WindowRedraw::~WindowRedraw()
{
  EventManager::unregister_mailbox(mailbox_);
}


bool
WindowRedraw::iterate()
{
  try {
    event_handle_t event = 0;
    bool redraw = false;
    while (mailbox_->tryReceive(event)) {
      if (!redraw && event->is_redraw_event()) {
        redraw = true;
      }
    }

    if (redraw) {
      window_->lock_.lock();
      window_->process_event(redraw_event_);
      window_->lock_.unlock();
    }
  } catch (...) {
    cerr << "Ignoring ERROR in redraw: " << window_->get_id() << endl;
  }
  return true;
}


}
}

