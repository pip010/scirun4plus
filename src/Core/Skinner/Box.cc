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
//    File   : Box.cc
//    Author : McKay Davis
//    Date   : Tue Jun 27 12:59:23 2006

#include <Core/Skinner/Variables.h>
#include <Core/Skinner/Box.h>
#include <Core/Events/EventManager.h>
#include <Core/Events/keysyms.h>
#include <Core/Util/StringUtil.h>
#include <sci_gl.h>


namespace SCIRun {
namespace Skinner {

Box::Box(Variables *variables) :
  Parent(variables),
  color_(variables, "color", Color(1.0, 0.0, 0.0, 1.0)),
  inside_(false),
  pressed_inside_(false)
{
  REGISTER_CATCHER_TARGET(Box::redraw);
  REGISTER_CATCHER_TARGET(Box::do_PointerEvent);
  REGISTER_CATCHER_TARGET(Box::do_KeyEvent);
}


Box::~Box()
{
}


BaseTool::propagation_state_e
Box::redraw(event_handle_t &)
{
  glColor4dv(&(color_().r));
  glBegin(GL_QUADS);

  const RectRegion &reg = get_region();
  glVertex2d(reg[0], reg[1]);
  glVertex2d(reg[2], reg[1]);
  glVertex2d(reg[2], reg[3]);
  glVertex2d(reg[0], reg[3]);
  glEnd();
  return CONTINUE_E;
}


BaseTool::propagation_state_e
Box::do_PointerEvent(event_handle_t &event)
{
  ASSERT(dynamic_cast<PointerSignal *>(event.get_rep()));
  PointerSignal *ps = (PointerSignal *)(event.get_rep());
  PointerEvent *pointer = ps->get_pointer_event();
  Signal *signal = 0;
  event_handle_t result = 0;
  ASSERT(pointer);

  inside_ = get_region().inside(pointer->get_x(), pointer->get_y());
  bool pressed =
    pointer->get_pointer_state() & PointerEvent::BUTTON_PRESS_E;


  if (inside_ && pressed &&
      (pointer->get_pointer_state() & PointerEvent::BUTTON_1_E))
  {
    pressed_inside_ = true;
    result = throw_signal("button_1_clicked");
  }

  if (//pressed_inside_ &&
      (pointer->get_pointer_state() & PointerEvent::BUTTON_RELEASE_E) &&
      (pointer->get_pointer_state() & PointerEvent::BUTTON_1_E)) {
    pressed_inside_ = false;
    result = throw_signal("button_1_released");
  }

  if (result.get_rep()) {
    signal = dynamic_cast<Signal *>(result.get_rep());
    if (signal) {
      return signal->get_signal_result();
    }
  }

  return CONTINUE_E;
}


BaseTool::propagation_state_e
Box::do_KeyEvent(event_handle_t &event)
{
  ASSERT(dynamic_cast<KeySignal *>(event.get_rep()));
  KeySignal *ks = (KeySignal *)(event.get_rep());
  KeyEvent *key = ks->get_key_event();
  Signal *signal = 0;
  event_handle_t result = 0;
  ASSERT(key);

  bool pressed =
    key->get_key_state() & KeyEvent::KEY_PRESS_E;

  if (inside_ && pressed &&
      key->get_keyval() == SCIRun_1)
  {
    result = throw_signal("key_1_pressed");
  }

  if (inside_ && pressed &&
      key->get_keyval() == SCIRun_2)
  {
    result = throw_signal("key_2_pressed");
  }

  if (inside_ && pressed &&
      key->get_keyval() == SCIRun_3)
  {
    result = throw_signal("key_3_pressed");
  }

  if (inside_ && pressed &&
      key->get_keyval() == SCIRun_4)
  {
    result = throw_signal("key_4_pressed");
  }

  if (result.get_rep()) {
    signal = dynamic_cast<Signal *>(result.get_rep());
    if (signal) {
      return signal->get_signal_result();
    }
  }

  return CONTINUE_E;
}


int
Box::get_signal_id(const string &signalname) const
{
  if (signalname == "button_1_clicked") return 1;
  if (signalname == "button_1_released") return 2;
  if (signalname == "button_2_clicked") return 3;
  if (signalname == "key_1_pressed") return 4;
  if (signalname == "key_2_pressed") return 5;
  if (signalname == "key_3_pressed") return 6;
  if (signalname == "key_4_pressed") return 7;
  return 0;
}


}
}
