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
//    File   : WindowRedraw.h
//    Author : McKay Davis
//    Date   : Thu Feb 22 12:18:23 2007

#ifndef SKINNER_WINDOW_REDRAW_H
#define SKINNER_WINDOW_REDRAW_H

#include <Core/Events/EventManager.h>
#include <Core/Util/ThrottledRunnable.h>

namespace SCIRun {
namespace Skinner {


class GLWindow;
class WindowRedraw : public ThrottledRunnable {
public:
  WindowRedraw(GLWindow *window, double hertz);
  virtual ~WindowRedraw();

private:
  virtual bool                            iterate();

  GLWindow *                              window_;
  SCIRun::EventManager::event_mailbox_t * mailbox_;
  event_handle_t                          redraw_event_;
};


}
}

#endif
