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
//    File   : Signals.h
//    Author : McKay Davis
//    Date   : Thu Jun 29 18:13:09 2006

#ifndef SKINNER_SIGNALS_H
#define SKINNER_SIGNALS_H

#include <Core/Events/BaseEvent.h>
#include <Core/Events/BaseTool.h>

#include <Core/Util/Environment.h>
#include <Core/Skinner/Variables.h>
#include <string>
#include <deque>
#include <vector>
#include <list>

#define REGISTER_CATCHER_TARGET(catcher_target_function_name) \
  this->register_callback \
    (this->make_callback(#catcher_target_function_name, this, \
       &catcher_target_function_name));

#define REGISTER_CATCHER_TARGET_BY_NAME(target_name, function_name) \
  this->register_callback \
    (this->make_callback(#target_name, this, &function_name));


#include <Core/Skinner/share.h>

namespace SCIRun {

class PointerEvent;

namespace Skinner {

class Signal;
class SignalThrower;
class SignalCatcher;
class Variables;
class Drawable;
class Callback;

typedef BaseTool::propagation_state_e result_t;

class SCISHARE Signal : public SCIRun::BaseEvent
{
public:
  Signal(const std::string &name, SignalThrower *, Variables *);
  Signal(const Signal &copy);
  virtual ~Signal();

  std::string               get_signal_name();
  void                      set_signal_name(const std::string &);

  Variables *               get_vars();
  void                      set_vars(Variables *);

  SignalThrower *           get_signal_thrower();
  void                      set_signal_thrower(SignalThrower *t);

  result_t                  get_signal_result();
  void                      set_signal_result(result_t);

  static result_t           perform_callback(Signal *, Callback*);

  virtual Signal *          clone() { return new Signal(*this); }
  virtual void              io(Piostream&);
  static PersistentTypeID   type_id;

private:
  std::string               signal_name_;
  Variables *               variables_;
  SignalThrower *           thrower_;
  result_t                  result_;
};


struct SCISHARE Callback {
  Callback();
  Callback(const Callback &copy);
  virtual ~Callback();
  virtual result_t          operator()(event_handle_t &) = 0;
  virtual Callback*         clone() = 0;
  std::string               targetname_;
  Variables*                variables_;
  Var<bool>                 threaded_;
};


template <class T>
struct CallbackT : public Callback {
  typedef result_t          (T::*FuncPtr)(event_handle_t &);

  CallbackT() : Callback(), catcher_(0), function_(0) {};
  CallbackT(const CallbackT &copy) :
    Callback(copy),
    catcher_(copy.catcher_),
    function_(copy.function_)
  {
    // empty
  }

  virtual ~CallbackT() {}

  void
  print_error(event_handle_t &e, const std::string &msg) {
    Signal *sig = dynamic_cast<Signal *>(e.get_rep());
    BaseTool *from = sig?dynamic_cast<BaseTool *>(sig->get_signal_thrower()):0;
    BaseTool *to = dynamic_cast<BaseTool *>(catcher_);
    cerr << "ERROR in callback:" << msg
         << "\n\tCallback ID:" << targetname_
         << "\n\tthrown from: " << (from ? from->get_name() : "???")
         << "\n\tthrown to  : " << (to ? to->get_name() : "???") << endl;
  }

  virtual result_t          operator()(event_handle_t &signal) {
    try {
      ASSERT(catcher_);
      ASSERT(function_);
      ASSERT(signal.get_rep());
      return (catcher_->*function_)(signal);
    } catch(const Exception& e){
      std::string err = std::string("Exception: ")+ e.message() +"\n";
      err = err + e.stackTrace()+"\n\n";
      print_error(signal, err);
    } catch (const std::string &err) {
      print_error(signal, err);
    } catch (const char * &err) {
      print_error(signal, err);
    } catch (...) {
      print_error(signal, " unknown problem");

    }
    return BaseTool::STOP_E;

  }

  virtual Callback*  clone() {
    return new CallbackT<T>(*this);
  }

  T*                        catcher_;
  FuncPtr                   function_;
};


typedef std::list<Callback*>           Callbacks_t;
typedef std::map<std::string, Callbacks_t>  CallbacksMap_t;


class SCISHARE SignalCatcher {
public:
  SignalCatcher();
  ~SignalCatcher();
  typedef result_t          CatcherFunction_t(event_handle_t &);
  CallbacksMap_t            slots_;

protected:
  template <class T>
  Callback* make_callback(const string &targetname,
                          T* caller,
                          typename CallbackT<T>::FuncPtr funptr)
  {
    ASSERT(funptr);
    CallbackT<T>* callback = new CallbackT<T>;
    callback->catcher_ = caller;
    callback->function_ = funptr;
    callback->targetname_ = targetname;

    slots_[targetname].push_back(callback);
    return callback;
  }
};


class SCISHARE SignalThrower {
public:
  SignalThrower();
  virtual ~SignalThrower();

  event_handle_t            throw_signal(event_handle_t &);
  static event_handle_t     throw_signal(CallbacksMap_t &,
                                         event_handle_t &);
  static event_handle_t     throw_signal(Callbacks_t &,
                                         event_handle_t &);

  void                      register_callback(Callback *);
  virtual int               get_signal_id(const std::string &) const =0;

  CallbacksMap_t            callbacks_;
};


class SCISHARE SignalCallback {
public:
  SignalCallback(Drawable *, const std::string &signalname);
  SignalCallback(Signal *signal);
  virtual ~SignalCallback();
  result_t                  operator()();
  Signal *                  get_signal() { return signal_; }

private:
  Signal *                  signal_;
  event_handle_t            signalh_;
  Callbacks_t::iterator     cbegin_;
  Callbacks_t::iterator     cend_;
  bool                      valid_;
  bool                      init_;
};


class SCISHARE MakerSignal : public Signal
{
  Variables * variables_;
public:
  MakerSignal(const std::string &name, Variables *vars) :
    Signal(name, 0, vars), variables_(vars) {}
  Variables *       get_vars() { return variables_; }
};


class SCISHARE PointerSignal : public Signal
{
  PointerEvent *pointer_;
public:
  PointerSignal(const std::string &name,
                SignalThrower *thrower,
                PointerEvent *pointer = 0) :
    Signal(name, thrower, 0), pointer_(pointer) {}
  PointerEvent *    get_pointer_event() { return pointer_; }
  void              set_pointer_event(PointerEvent *p) { pointer_ = p; }
};


class SCISHARE KeySignal : public Signal
{
  KeyEvent *key_;
public:
  KeySignal(const std:: &name,
            SignalThrower *thrower,
            KeyEvent *key = 0) :
    Signal(name, thrower, 0), key_(key) {}
  KeyEvent *        get_key_event() { return key_; }
  void              set_key_event(KeyEvent *k) { key_ = k; }
};


}
}

#endif
