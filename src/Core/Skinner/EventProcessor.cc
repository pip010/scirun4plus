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
//    File   : EventProcessor.cc
//    Author : McKay Davis
//    Date   : Thu Feb 22 12:35:51 2007

#include <Core/Events/EventManager.h>
#include <Core/Skinner/EventProcessor.h>
#include <Core/Skinner/Window.h>
#include <iostream>

namespace SCIRun {
namespace Skinner {


EventProcessor::EventProcessor(GLWindow *window, double hertz) :
  ThrottledRunnable(hertz),
  window_(window),
  window_id_(window->get_id()),
  mailbox_(SCIRun::EventManager::register_event_messages(window_id_))
{
}


EventProcessor::~EventProcessor()
{
  delete window_;
  SCIRun::EventManager::unregister_mailbox(mailbox_);
}


bool
EventProcessor::iterate()
{
  event_handle_t event = 0;
  window_->lock_.lock();
  while (mailbox_->tryReceive(event)) {
    if (event->is_redraw_event()) {
      continue;
    }

    try {
      window_->process_event(event);
    } catch (...) {
      cerr << "EventProcessor error passing event to window id: "
           << window_id_ << endl;
    }

    // Quit this throttled runnable if QuitEvent type received
    if (event->is_quit_event()) {
      window_->lock_.unlock();
      return false;
    }
  }
  window_->lock_.unlock();
  return true;
}


}
}

