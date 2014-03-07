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
//    File   : SliderButton.cc
//    Author : McKay Davis
//    Date   : Tue Feb  6 17:22:56 2007

#include <Core/Skinner/SliderButton.h>
#include <Core/Skinner/Slider.h>
#include <Core/Skinner/Variables.h>
#include <Core/Events/EventManager.h>

namespace SCIRun {
namespace Skinner {


SliderButton::SliderButton(Variables *vars) :
  Parent(vars),
  slider_(0),
  width_(vars, "SliderButton::width", 10.0),
  x_(0),
  y_(0)
{
  REGISTER_CATCHER_TARGET(SliderButton::do_PointerEvent);
}


SliderButton::~SliderButton()
{
}


BaseTool::propagation_state_e
SliderButton::do_PointerEvent(event_handle_t &event)
{
  PointerSignal *signal = dynamic_cast<PointerSignal *>(event.get_rep());
  PointerEvent *pointer = signal->get_pointer_event();
  x_ = pointer->get_x();
  y_ = pointer->get_y();
  return CONTINUE_E;
}


BaseTool::propagation_state_e
SliderButton::process_event(event_handle_t &event)
{
  set_region(slider_->get_subregion());
  return Parent::process_event(event);
}


}
}
