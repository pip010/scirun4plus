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
//    File   : Arithmetic.cc
//    Author : McKay Davis
//    Date   : Tue Aug  8 20:02:54 2006

#include <Core/Skinner/Variables.h>
#include <Core/Skinner/Arithmetic.h>
#include <Core/Events/EventManager.h>
#include <Core/Util/StringUtil.h>


namespace SCIRun {
namespace Skinner {

Arithmetic::Arithmetic(Variables *variables) :
  Parent(variables)
{
  REGISTER_CATCHER_TARGET(Arithmetic::increment);
  REGISTER_CATCHER_TARGET(Arithmetic::decrement);
  REGISTER_CATCHER_TARGET(Arithmetic::add);
  REGISTER_CATCHER_TARGET(Arithmetic::subtract);
  REGISTER_CATCHER_TARGET(Arithmetic::continue_if);
  REGISTER_CATCHER_TARGET(Arithmetic::set_value);
  REGISTER_CATCHER_TARGET(Arithmetic::invert_bool);
}


Arithmetic::~Arithmetic()
{
}


BaseTool::propagation_state_e
Arithmetic::increment(event_handle_t &event)
{
  ASSERT(dynamic_cast<Skinner::Signal *>(event.get_rep()));
  Skinner::Signal *signal = (Skinner::Signal *)(event.get_rep());
  switch (signal->get_vars()->get_type_e("variable")) {
  case Variables::INT_E: {
    Var<int> val(signal->get_vars(), "variable", 0);
    val = val() + 1;
  } break;

  case Variables::DOUBLE_E: {
    Var<double> val(signal->get_vars(), "variable", 0.0);
    val = val() + 1.0;
  } break;

  default:
    throw "Bad type";
    break;
  }

  return BaseTool::CONTINUE_E;
}


BaseTool::propagation_state_e
Arithmetic::decrement(event_handle_t &event)
{
  ASSERT(dynamic_cast<Skinner::Signal *>(event.get_rep()));
  Skinner::Signal *signal = (Skinner::Signal *)(event.get_rep());

  switch (signal->get_vars()->get_type_e("variable")) {
  case Variables::INT_E: {
    Var<int> val(signal->get_vars(), "variable", 0);
    val = val() - 1;
  } break;

  case Variables::DOUBLE_E: {
    Var<double> val(signal->get_vars(), "variable", 0.0);
    val = val() - 1.0;
  } break;

  default:
    throw "Bad type";
    break;
  }

  return BaseTool::CONTINUE_E;
}


BaseTool::propagation_state_e
Arithmetic::add(event_handle_t &event)
{
  ASSERT(dynamic_cast<Skinner::Signal *>(event.get_rep()));
  Skinner::Signal *signal = (Skinner::Signal *)(event.get_rep());
  switch (signal->get_vars()->get_type_e("variable")) {
  case Variables::INT_E: {
    Var<int> val(signal->get_vars(), "variable", 0);
    Var<int> amount(signal->get_vars(), "amount", 0);
    val = val + amount;
  } break;

  case Variables::DOUBLE_E: {
    Var<double> val(signal->get_vars(), "variable", 0.0);
    Var<double> amount(signal->get_vars(), "amount", 0.0);
    val = val + amount;
  } break;

  default:
    throw "Bad type";
    break;
  }

  return BaseTool::CONTINUE_E;
}


BaseTool::propagation_state_e
Arithmetic::subtract(event_handle_t &event)
{
  ASSERT(dynamic_cast<Skinner::Signal *>(event.get_rep()));
  Skinner::Signal *signal = (Skinner::Signal *)(event.get_rep());
  switch (signal->get_vars()->get_type_e("variable")) {
  case Variables::INT_E: {
    Var<int> val(signal->get_vars(), "variable", 0);
    Var<int> amount(signal->get_vars(), "amount", 0);
    val = val - amount;
  } break;

  case Variables::DOUBLE_E: {
    Var<double> val(signal->get_vars(), "variable", 0.0);
    Var<double> amount(signal->get_vars(), "amount", 0.0);
    val = val - amount;
  } break;

  default:
    throw "Bad type";
    break;
  }

  return BaseTool::CONTINUE_E;
}


BaseTool::propagation_state_e
Arithmetic::set_value(event_handle_t &event)
{
  ASSERT(dynamic_cast<Skinner::Signal *>(event.get_rep()));
  Skinner::Signal *signal = (Skinner::Signal *)(event.get_rep());
  signal->get_vars()->copy_var("value", "variable");
  return BaseTool::CONTINUE_E;
}


BaseTool::propagation_state_e
Arithmetic::invert_bool(event_handle_t &event)
{
  ASSERT(dynamic_cast<Skinner::Signal *>(event.get_rep()));
  Skinner::Signal *signal = (Skinner::Signal *)(event.get_rep());
  Var<bool> var(signal->get_vars(), "variable");
  var = !var();
  return BaseTool::CONTINUE_E;
}


BaseTool::propagation_state_e
Arithmetic::continue_if(event_handle_t &event)
{
  ASSERT(dynamic_cast<Skinner::Signal *>(event.get_rep()));
  Skinner::Signal *signal = (Skinner::Signal *)(event.get_rep());
  Var<bool> value(signal->get_vars(), "value",1);
  Var<bool> invert(signal->get_vars(), "invert",0);
  if (invert)
    return value() ? STOP_E : CONTINUE_E;
  else
    return value() ? CONTINUE_E : STOP_E;
}


}
}
