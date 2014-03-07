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
//    File   : FocusRegion.cc
//    Author : McKay Davis
//    Date   : Tue Aug 15 22:11:01 2006

#include <Core/Skinner/FocusRegion.h>
#include <Core/Skinner/Variables.h>
#include <Core/Events/BaseEvent.h>
#include <Core/Geometry/BBox.h>
#include <Core/Math/MinMax.h>
#include <Core/Math/MiscMath.h>
#include <Core/Events/keysyms.h>

namespace SCIRun {
namespace Skinner {


// Cache the last pointer position so that on a visibility event
// skinner can tell if this item has focus or not.  Note that we
// only cache pointer events going through this item so they may not
// be accurate.  May need to do this deeper in the event tree so
// that we can get a more accurate position.
static int last_x_ = -1;
static int last_y_ = -1;


FocusRegion::FocusRegion(Variables *vars) :
  Parent(vars),
  focus_(vars,"focus",false),
  reset_(0),
  x_(-1),
  y_(-1),
  vis_cache_(false)
{
  REGISTER_CATCHER_TARGET(FocusRegion::do_PointerEvent);
  REGISTER_CATCHER_TARGET(FocusRegion::do_KeyEvent);
  REGISTER_CATCHER_TARGET_BY_NAME
    (Skinner::FocusRegion_Maker, FocusRegion::FocusRegion_Maker);
  REGISTER_CATCHER_TARGET(FocusRegion::reset_focus);
}


FocusRegion::~FocusRegion()
{
}


BaseTool::propagation_state_e
FocusRegion::FocusRegion_Maker(event_handle_t &event)
{
  FocusRegion *focus_region =
    construct_child_from_maker_signal<FocusRegion>(event);
  ASSERT(focus_region);
  focus_regions_.push_back(focus_region);
  return STOP_E;
}


BaseTool::propagation_state_e
FocusRegion::do_PointerEvent(event_handle_t &event)
{
  ASSERT(dynamic_cast<PointerSignal *>(event.get_rep()));
  PointerSignal *signal = (PointerSignal *)(event.get_rep());

  PointerEvent *pointer = signal->get_pointer_event();

  unsigned int buttons = (PointerEvent::BUTTON_1_E |
                          PointerEvent::BUTTON_2_E |
                          PointerEvent::BUTTON_3_E |
                          PointerEvent::BUTTON_4_E |
                          PointerEvent::BUTTON_5_E);

  if (pointer->get_pointer_state() & buttons) {
    return CONTINUE_E;
  }

  x_ = pointer->get_x();
  y_ = pointer->get_y();
  last_x_ = x_;
  last_y_ = y_;

  return do_Focus(event);
}


BaseTool::propagation_state_e
FocusRegion::do_Focus(event_handle_t &)
{
  if (get_region().inside(x_, y_)) {
    set_focus(true);
    return CONTINUE_E;
  }
  set_focus(false);
  return STOP_E;
}



BaseTool::propagation_state_e
FocusRegion::do_KeyEvent(event_handle_t &event)
{
  if (!focus_()) return STOP_E;

  KeySignal *signal = dynamic_cast<KeySignal *>(event.get_rep());
  ASSERT(signal);

  KeyEvent *key = signal->get_key_event();
  if (key->get_key_state() & KeyEvent::KEY_PRESS_E &&
      key->get_keyval() == SCIRun_Tab &&
      focus_regions_.size() > 1) {

    bool backward = key->get_modifiers() & EventModifiers::SHIFT_E;
    FocusRegions_t subregions;
    FocusRegions_t::size_type num = focus_regions_.size();
    for (unsigned int i = 0; i < num; ++i) {
      if (focus_regions_[i]->focus_ &&
          focus_regions_[i]->focus_regions_.empty())
      {
        Var<bool> notab(focus_regions_[i]->get_vars(), "ignore_tab",0);

        if (!backward && notab()) {
          break;
        }
        focus_regions_[i]->set_focus(false);
        int next = 0;
        do {
          next = (int(i) + (backward ? ( int(num) - 1) : 1)) % num;
        } while (!focus_regions_[next]->get_vars()->get_bool("visible"));
        focus_regions_[next]->set_focus(true);
        break;
      }
    }
  }
  return CONTINUE_E;
}


void
FocusRegion::set_focus(bool focus)
{
  if (focus_ != focus) {
    if (focus) {
      throw_signal("FocusRegion::focus");
    }
    else {
      FocusRegions_t::size_type num = focus_regions_.size();
      for (unsigned int i = 0; i < num; ++i) {
        focus_regions_[i]->set_focus(focus);
      }

      throw_signal("FocusRegion::unfocus");
    }
    focus_ = focus;
  }
}


BaseTool::propagation_state_e
FocusRegion::reset_focus(event_handle_t &)
{
  set_focus(false);
  reset_ = true;

  return CONTINUE_E;
}


BaseTool::propagation_state_e
FocusRegion::process_event(event_handle_t &event)
{
  if (vis_cache_ != visible_) {
    reset_focus(event);
    if (visible_)
    {
      x_ = last_x_;
      y_ = last_y_;
    }
    else
    {
      x_ = -1;
      y_ = -1;
    }
  }
  vis_cache_ = visible_;

  const bool isa_pointer = event->is_pointer_event();
  const bool isa_key     = event->is_key_event();
  const bool isa_win     = event->is_window_event();

  if (reset_ && (isa_pointer || isa_key || isa_win)) {
    do_Focus(event);
    reset_ = false;
  }

  return Parent::process_event(event);
}


int
FocusRegion::get_signal_id(const string &id) const
{
  if (id == "FocusRegion::focus") return 1;
  if (id == "FocusRegion::unfocus") return 2;
  return 0;
}


}
}
