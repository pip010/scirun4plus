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
//    File   : ViewSubRegion.cc
//    Author : McKay Davis
//    Date   : Sat Aug 12 12:57:50 2006


#include <Core/Skinner/ViewSubRegion.h>
#include <Core/Skinner/Variables.h>
#include <Core/Events/EventManager.h>
#include <Core/Math/MinMax.h>
#include <Core/Geometry/Point.h>
#include <sci_gl.h>

namespace SCIRun {
namespace Skinner {


ViewSubRegion::ViewSubRegion(Variables *vars) :
  Parent(vars),
  x_(vars, "ViewSubRegion::x"),
  y_(vars, "ViewSubRegion::y"),
  width_(vars, "ViewSubRegion::width"),
  height_(vars, "ViewSubRegion::height"),
  invisible_width_(vars, "ViewSubRegion::invisible_width"),
  invisible_height_(vars, "ViewSubRegion::invisible_height"),
  cache_region_()
{
  REGISTER_CATCHER_TARGET(ViewSubRegion::do_PointerEvent);
}


ViewSubRegion::~ViewSubRegion()
{
}


BaseTool::propagation_state_e
ViewSubRegion::process_event(event_handle_t &event)
{
  cache_region_ = get_region();
  RectRegion virtual_region = cache_region_;
  if (x_.exists() && y_.exists() && width_.exists() && height_.exists())
  {
    virtual_region = RectRegion(cache_region_.x1() - x_, cache_region_.y2() + y_ - height_,
                                cache_region_.x1() - x_ + width_, cache_region_.y2() + y_);
    if (invisible_width_.exists())
      invisible_width_ = Max(width_ - cache_region_.width(),0.);
    if (invisible_height_.exists())
      invisible_height_ = Max(height_ - cache_region_.height(),0.);
  }

  WindowEvent *window = dynamic_cast<WindowEvent *>(event.get_rep());
  bool redraw = window && window->get_window_state() == WindowEvent::REDRAW_E;
  if (redraw) {
    RectRegion view_region = get_region();

    bool scissor_on = glIsEnabled(GL_SCISSOR_TEST);
    GLint scissor[4];
    if (scissor_on) {
      glGetIntegerv(GL_SCISSOR_BOX, scissor);

      double x1 = Max(view_region.x1(), (double)scissor[0]);
      double x2 = Min(view_region.x2(), (double)scissor[0]+(double)scissor[2]);

      double y1 = Max(view_region.y1(), (double)scissor[1]);
      double y2 = Min(view_region.y2(), (double)scissor[1]+(double)scissor[3]);

      if (x1 >= x2 || y1 >= y2) {
        return CONTINUE_E;
      }

      view_region = RectRegion(x1, y1, x2, y2);

    } else {
      glEnable(GL_SCISSOR_TEST);
    }

    glScissor((GLint)view_region.x1(), (GLint)view_region.y1(),
              (GLint)view_region.width(), (GLint)view_region.height());
    set_region(virtual_region);
    Parent::process_event(event);
    set_region(cache_region_);
    if (scissor_on) {
      glScissor(scissor[0], scissor[1], scissor[2], scissor[3]);
    } else {
      glDisable(GL_SCISSOR_TEST);
    }
  } else {
    set_region(virtual_region);
    if (event->is_pointer_event()) {
      PointerEvent *pe = dynamic_cast<PointerEvent *>(event.get_rep());
      bool rel = pe->get_pointer_state() & PointerEvent::BUTTON_RELEASE_E;
      ASSERT(pe);
      if (!rel && !cache_region_.inside(pe->get_x(), pe->get_y())) {
        return STOP_E;
      }
    }
    Parent::process_event(event);
    set_region(cache_region_);
  }
  return CONTINUE_E;
}


BaseTool::propagation_state_e
ViewSubRegion::do_PointerEvent(event_handle_t &event)
{
  PointerSignal *signal = dynamic_cast<PointerSignal *>(event.get_rep());
  PointerEvent *pointer = signal->get_pointer_event();

  bool pressed =
    pointer->get_pointer_state() & PointerEvent::BUTTON_PRESS_E;
  bool button4 = pointer->get_pointer_state() & PointerEvent::BUTTON_4_E;
  bool button5 = pointer->get_pointer_state() & PointerEvent::BUTTON_5_E;

  if (pressed && y_.exists() && get_region().inside(pointer->get_x(), pointer->get_y())) {
    if (button5) {
      y_ = Min(y_ + cache_region_.height()/2.0, height_());
      EventManager::add_event(new WindowEvent(WindowEvent::REDRAW_E));
    } else if (button4) {
      y_ = Max(y_ - cache_region_.height()/2.0, 0.0);
      EventManager::add_event(new WindowEvent(WindowEvent::REDRAW_E));
    }
  }
  return CONTINUE_E;
}


}
}
