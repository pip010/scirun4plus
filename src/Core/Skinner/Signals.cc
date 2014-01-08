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
//    File   : Signals.cc
//    Author : McKay Davis
//    Date   : Thu Jun 29 18:14:47 2006


#include <Core/Skinner/Signals.h>
#include <Core/Skinner/Variables.h>
#include <Core/Util/Environment.h>
#include <Core/Thread/Thread.h>
#include <Core/Thread/Runnable.h>
#include <Core/Skinner/Drawable.h>
#include <iostream>

using std::cerr;
using std::endl;

namespace SCIRun {
namespace Skinner {


Persistent *make_Signal() { return new Signal("",0,0); }
PersistentTypeID Signal::type_id("Signal", "BaseEvent", &make_Signal);


Signal::Signal(const string &name,
               SignalThrower *thrower,
               Variables *vars) :
  BaseEvent("", 0),
  signal_name_(name),
  variables_(vars),
  thrower_(thrower),
  result_(BaseTool::STOP_E)
{
}


Signal::Signal(const Signal &copy) :
  BaseEvent("", 0),
  signal_name_(copy.signal_name_),
  variables_(copy.variables_),
  thrower_(copy.thrower_),
  result_(copy.result_)
{
}


Signal::~Signal()
{
}


string
Signal::get_signal_name()
{
  return signal_name_;
}


void
Signal::set_signal_name(const string &name)
{
  signal_name_ = name;
}


Variables *
Signal::get_vars()
{
  return variables_;
}


void
Signal::set_vars(Variables *vars)
{
  variables_ = vars;
}


SignalThrower *
Signal::get_signal_thrower()
{
  return thrower_;
}


void
Signal::set_signal_thrower(SignalThrower *t)
{
  thrower_ = t;
}


result_t
Signal::get_signal_result()
{
  return result_;
}


void
Signal::set_signal_result(result_t result)
{
  result_ = result;
}


void
Signal::io(Piostream &)
{
}


Callback::Callback() :
  targetname_(),
  variables_(new Variables()),
  threaded_(variables_, "threaded", 0)
{
}


Callback::Callback(const Callback &copy) :
  targetname_(copy.targetname_),
  variables_(new Variables(*copy.variables_)),
  threaded_(variables_, "threaded", 0)
{
}


Callback::~Callback()
{
  ASSERT(variables_);
  delete variables_;
}


SignalCatcher::SignalCatcher() :
  slots_()
{
}


SignalCatcher::~SignalCatcher()
{
  CallbacksMap_t::iterator cbeg = slots_.begin();
  CallbacksMap_t::iterator cend = slots_.end();
  for (;cbeg != cend; ++cbeg) {
    Callbacks_t::iterator citer = cbeg->second.begin();
    Callbacks_t::iterator clast = cbeg->second.end();
    for (;citer != clast; ++citer) {
      delete *citer;
    }
  }
}


SignalThrower::SignalThrower() :
  callbacks_()
{
}


SignalThrower::~SignalThrower()
{
}


void
SignalThrower::register_callback(Callback* callback)
{
  callbacks_[callback->targetname_].push_back(callback);
}


event_handle_t
SignalThrower::throw_signal(event_handle_t &signal)
{
  return throw_signal(callbacks_, signal);
}



class ThreadedCallback : public Runnable {
public:
  ThreadedCallback(Callback *callback,
                   const event_handle_t &signal) :
    callback_(callback), signal_(0), state_(BaseTool::STOP_E)
  {
    signal_ = signal;
    signal_.detach();
  }

  virtual ~ThreadedCallback() { signal_ = 0; }

  virtual void run() {
    ASSERT(callback_);
    ASSERT(dynamic_cast<Signal *>(signal_.get_rep()));
    state_ = (*callback_)(signal_);
  }

  BaseTool::propagation_state_e     get_state() { return state_; }
private:
  Callback *    callback_;
  event_handle_t                            signal_;
  BaseTool::propagation_state_e             state_;
};


event_handle_t
SignalThrower::throw_signal(CallbacksMap_t &catchers,
                            event_handle_t &signalh)
{
  ASSERT(dynamic_cast<Signal *>(signalh.get_rep()));
  Signal *signal = (Signal *)(signalh.get_rep());
  CallbacksMap_t::iterator iter = catchers.find(signal->get_signal_name());

  if (iter != catchers.end()) {
    return throw_signal(iter->second, signalh);
  }

  signal->set_signal_result(BaseTool::CONTINUE_E);
  return signalh;
}


event_handle_t
SignalThrower::throw_signal(Callbacks_t &callbacks,
                            event_handle_t &signalh)
{
  ASSERT(dynamic_cast<Signal *>(signalh.get_rep()));
  Signal *signal = (Signal *)(signalh.get_rep());
  result_t state = BaseTool::CONTINUE_E;
  Callbacks_t::iterator riter = callbacks.begin();
  Callbacks_t::iterator rend = callbacks.end();
  for(;riter != rend; ++riter) {
    Callback* callback = *riter;
    signal->set_vars(callback->variables_);
    if (!callback->threaded_) {
      state = (*callback)(signalh);
    } else {
      ThreadedCallback *thread = new ThreadedCallback(callback,signalh);
      (new Thread(thread, signal->get_signal_name().c_str()))->detach();
    }
    if (state == BaseTool::STOP_E) break;
  }
  signal->set_signal_result(state);
  return signalh;
}


SignalCallback::SignalCallback(Drawable *thrower, const string &signalname)
  : signal_(new Signal(signalname, thrower, thrower->get_vars())),
    signalh_(signal_),
    valid_(false),
    init_(false)
{
}


SignalCallback::SignalCallback(Signal *signal)
  : signal_(signal),
    signalh_(signal),
    valid_(false),
    init_(false)
{
}


SignalCallback::~SignalCallback()
{
}


BaseTool::propagation_state_e
SignalCallback::operator()()
{
  if (!init_) {
    CallbacksMap_t::iterator aciter =
      signal_->get_signal_thrower()->callbacks_.find
      (signal_->get_signal_name());

    if (aciter != signal_->get_signal_thrower()->callbacks_.end()) {
      cbegin_ = aciter->second.begin();
      cend_ = aciter->second.end();
      valid_ = true;
    }
    init_ = true;
  }

  BaseTool::propagation_state_e state = BaseTool::CONTINUE_E;
  if (!valid_) return state;

  Callbacks_t::iterator iter = cbegin_;
  for(;iter != cend_; ++iter) {
    Callback* callback = *iter;
    Variables *vars = callback->variables_;
    signal_->set_vars(vars);
    if (!callback->threaded_) {
      state = (*callback)(signalh_);
    } else {
      ThreadedCallback *runner = new ThreadedCallback(callback, signalh_);
      Thread *thread=new Thread(runner,signal_->get_signal_name().c_str());
      thread->detach();
    }
    if (state == BaseTool::STOP_E) break;
  }

  signal_->set_signal_result(state);
  return state;
}


}
}

