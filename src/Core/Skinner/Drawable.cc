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
//    File   : Drawable.cc
//    Author : McKay Davis
//    Date   : Tue Jun 27 13:04:57 2006

#include <Core/Skinner/Drawable.h>
#include <Core/Skinner/Variables.h>

namespace SCIRun {
namespace Skinner {


Drawable::Drawable(Variables *variables) :
  BaseTool(variables ? variables->get_id() : ""),
  SignalCatcher(),
  SignalThrower(),
  parent_(0),
  visible_(variables, "visible",1),
  visible_cache_(visible_()),
  class_(variables, "name", ""),
  redraw_callback_(this, class_()+"::redraw"),
  pointer_callback_(new PointerSignal(class_()+"::do_PointerEvent", this)),
  key_callback_(new KeySignal(class_()+"::do_KeyEvent", this)),
  region_(),
  variables_(variables)
{
  REGISTER_CATCHER_TARGET(Drawable::bridge);
}


Drawable::~Drawable()
{
  if (variables_) delete variables_;
}


void
Drawable::set_parent(Drawable *parent)
{
  parent_ = parent;
}


MinMax
Drawable::get_minmax(unsigned int)
{
  return SPRING_MINMAX;
}


BaseTool::propagation_state_e
Drawable::process_event(event_handle_t &event)
{
  const bool isa_pointer = event->is_pointer_event();
  const bool isa_key     = event->is_key_event();
  const bool isa_win     = event->is_window_event();

  if (!visible_) {
    return (isa_pointer || isa_key || isa_win) ? STOP_E : CONTINUE_E;
  }

  if (event->is_redraw_event()) {
    return redraw_callback_();
  } else if (isa_pointer) {
    PointerSignal *ps = (PointerSignal *)pointer_callback_.get_signal();
    ps->set_pointer_event((PointerEvent *)event.get_rep());
    return pointer_callback_();
  } else if (isa_key) {
    KeySignal *ps = (KeySignal *)key_callback_.get_signal();
    ps->set_key_event((KeyEvent *)event.get_rep());
    return key_callback_();
  }

  return CONTINUE_E;
}


void
Drawable::get_modified_event(event_handle_t &)
{
}

string
Drawable::get_id() const
{
  return get_name();
}


int
Drawable::get_signal_id(const string &name) const
{
  if (name == "redraw") return 1;
  return 0;
}


const RectRegion &
Drawable::get_region() const
{
  return region_;
}

void
Drawable::set_region(const RectRegion &region)
{
  region_ = region;
}


Variables *
Drawable::get_vars()
{
  return variables_;
}


event_handle_t
Drawable::throw_signal(const string &signalname)
{
  event_handle_t signal = new Signal(signalname, this, variables_);
  return SignalThrower::throw_signal(callbacks_, signal);
}



event_handle_t
Drawable::throw_signal_extended(event_handle_t &event)
{
  ASSERT(dynamic_cast<Signal *>(event.get_rep()));
  Signal *signal = (Signal *)event.get_rep();
  Drawable *object = this;
  const string name = signal->get_signal_name();
  while (object) {
    CallbacksMap_t::iterator citer = object->callbacks_.find(name);
    if (citer != object->callbacks_.end()) {
      event = SignalThrower::throw_signal(object->callbacks_, event);
      ASSERT(dynamic_cast<Signal *>(event.get_rep()));
      signal = (Signal *)event.get_rep();
      if (signal->get_signal_result() == BaseTool::STOP_E) {
        break;
      }
    }
    object = object->parent_;
  }

  return event;
}


bool
Drawable::add_callback(const string &name,
                       const string &target,
                       const Variables *vars)
{
  Drawable *object = this;
  const bool found = get_signal_id(name);
  Callbacks_t &callbacks = found ? callbacks_[name] : slots_[name];
  int count = 0;
  while (object) {
    CallbacksMap_t::iterator citer = object->slots_.find(target);
    CallbacksMap_t::iterator cend = object->slots_.end();
    object = object->parent_;
    if (citer == cend) continue;

    Callbacks_t::iterator siter = citer->second.begin();
    Callbacks_t::iterator slast = citer->second.end();
    for (; siter != slast; ++siter) {
      Callback *callback = (*siter)->clone();
      callback->targetname_ = name + " -> " + callback->targetname_;
      ASSERT(callback->variables_);
      callback->variables_->merge(vars);
      callback->variables_->set_parent(this->get_vars());
      callbacks.push_back(callback);
      count++;
    }
  }
#if 0
  // Debugging print, detect disconnected/bad targets in skin files.
  if (count == 0)
  {
    cerr << "Skinner bad target '" << target << "', signal '" << name << "'\n";
  }
#endif
  return false;
}


bool
Drawable::make_bridge(const string &name)
{
  if (bridges_.find(name) != bridges_.end()) {
    return false;
  }

  bridges_[name] = 0;
  Callback *cback = make_callback<Drawable>(name, this, &Drawable::bridge);
  cback->variables_->insert("bridge", name);
  return true;
}


bool
Drawable::use_bridge(const string &name,
                     Drawable *target,
                     const string &targetname)
{
  Drawable *object = this;
  while (object) {
    Bridges_t::iterator biter = object->bridges_.find(name);
    if (biter != object->bridges_.end()) {
      biter->second = &target->slots_[targetname];
      return true;
    }
    object = object->parent_;
  }

  return false;
}


BaseTool::propagation_state_e
Drawable::bridge(event_handle_t &event)
{
  ASSERT(dynamic_cast<Signal *>(event.get_rep()));
  Signal *signal = (Signal *)(event.get_rep());
  Variables *vars = signal->get_vars();
  Var<string> name(vars, "bridge", "");

  Bridges_t::iterator biter = bridges_.find(name());
  if (biter != bridges_.end()) {
    event = SignalThrower::throw_signal(*biter->second, event);
    signal = (Signal *)(event.get_rep());
    return signal->get_signal_result();
  }

  return BaseTool::CONTINUE_E;
}


}
}

