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
//    File   : Slider.cc
//    Author : McKay Davis
//    Date   : Tue Feb  6 17:22:56 2007

#include <Core/Skinner/Slider.h>
#include <Core/Skinner/SliderButton.h>
#include <Core/Skinner/Variables.h>
#include <Core/Events/EventManager.h>
#include <Core/Math/MiscMath.h>

namespace SCIRun {
namespace Skinner {


Slider::Slider(Variables *vars) :
  Parent(vars),
  button_(0),
  min_(vars,"Slider::min", 0.0),
  max_(vars,"Slider::max", 1.0),
  resolution_(vars,"Slider::resolution", 0.1),
  value_(vars,"Slider::value", 0.5),
  vertical_(vars,"Slider::vertical", false),
  sx_(0),
  sy_(0),
  drag_(0),
  drag_value_(0.0)
{
  REGISTER_CATCHER_TARGET(Slider::SliderButton_Maker);
  REGISTER_CATCHER_TARGET(Slider::do_PointerEvent);
}


Slider::~Slider()
{
}


BaseTool::propagation_state_e
Slider::SliderButton_Maker(event_handle_t &maker_signal)
{
  button_ =
    construct_child_from_maker_signal<SliderButton>(maker_signal);
  button_->slider_ = this;
  return STOP_E;
}


bool
Slider::update_value()
{
  double value = value_;
  value /= resolution_;
  value = int(value > 0.0 ? value + 0.5 : value - 0.5);
  value *= resolution_;
  value = Clamp(value, min_, max_);
  if (value_ == value) return false;
  value_ = value;
  return true;
}


BaseTool::propagation_state_e
Slider::do_PointerEvent(event_handle_t &event)
{
  PointerSignal *signal = dynamic_cast<PointerSignal *>(event.get_rep());
  PointerEvent *pointer = signal->get_pointer_event();

  bool pressed =
    pointer->get_pointer_state() & PointerEvent::BUTTON_PRESS_E;
  bool button1 = pointer->get_pointer_state() & PointerEvent::BUTTON_1_E;
  bool button4 = pointer->get_pointer_state() & PointerEvent::BUTTON_4_E;
  bool button5 = pointer->get_pointer_state() & PointerEvent::BUTTON_5_E;
  int x = pointer->get_x();
  int y = pointer->get_y();

  if (pressed && get_region().inside(x, y)) {
    if (button4) {
      value_ = value_ + 10*resolution_*(vertical_()?-1:1);
      throw_signal("Slider::value_changed");
      EventManager::add_event(new WindowEvent(WindowEvent::REDRAW_E));
    } else if (button5) {
      value_ = value_ + 10*resolution_*(vertical_()?1:-1);
      throw_signal("Slider::value_changed");
      EventManager::add_event(new WindowEvent(WindowEvent::REDRAW_E));
    } else if (!vertical_ && x < get_subregion().x1()) {
      value_ = value_ - resolution_;
      throw_signal("Slider::value_changed");
      EventManager::add_event(new WindowEvent(WindowEvent::REDRAW_E));
    } else if (!vertical_ && x > get_subregion().x2()) {
      value_ = value_ + resolution_;
      throw_signal("Slider::value_changed");
      EventManager::add_event(new WindowEvent(WindowEvent::REDRAW_E));
    } else if (vertical_ && y < get_subregion().y1()) {
      value_ = value_ + resolution_;
      throw_signal("Slider::value_changed");
      EventManager::add_event(new WindowEvent(WindowEvent::REDRAW_E));
    } else if (vertical_ && y > get_subregion().y2()) {
      value_ = value_ - resolution_;
      throw_signal("Slider::value_changed");
      EventManager::add_event(new WindowEvent(WindowEvent::REDRAW_E));
    } else {
      drag_ = true;
      drag_value_ = value_;
      sx_ = x;
      sy_ = y;
    }
  }

  if (drag_ && !button1) {
    drag_ = false;
    throw_signal("Slider::value_changed");
  }

  if (drag_ && button1) {
    if (!vertical_) {
      double dx = x - sx_;
      double width = get_region().width() - button_->width_;
      double interp = dx / width;
      value_ = drag_value_ + interp*(max_-min_);
      if (update_value()) {
        throw_signal("Slider::value_changed");
        EventManager::add_event(new WindowEvent(WindowEvent::REDRAW_E));
      }
    } else {
      double dy = y - sy_;
      double height = get_region().height() - button_->width_;
      double interp = - dy / height;
      value_ = drag_value_ + interp*(max_-min_);
      if (update_value()) {
        throw_signal("Slider::value_changed");
        EventManager::add_event(new WindowEvent(WindowEvent::REDRAW_E));
      }
    }

  }

  return CONTINUE_E;
}


RectRegion
Slider::get_subregion()
{
  if (!vertical_) {
    double width = get_region().width() - button_->width_;
    double interp = (value_ - min_) / (max_ - min_);
    const RectRegion &reg = get_region();
    double x1 = Floor(reg.x1() + interp * width);
    double x2 = Ceil(x1 + button_->width_);

    return RectRegion(x1, reg.y1(), x2, reg.y2());
  } else {
    double height = get_region().height() - button_->width_;
    double interp = (value_ - min_) / (max_ - min_);
    const RectRegion &reg = get_region();
    double y2 = Floor(reg.y2() - interp * height);
    double y1 = Ceil(y2 - button_->width_);

    return RectRegion(reg.x1(), y1, reg.x2(), y2);
  }
}


BaseTool::propagation_state_e
Slider::process_event(event_handle_t &event)
{
  update_value();
  return Parent::process_event(event);
}


int
Slider::get_signal_id(const string &signalname) const
{
  if (signalname == "Slider::value_changing") return 1;
  if (signalname == "Slider::value_changed") return 2;
  return 0;
}


}
}
