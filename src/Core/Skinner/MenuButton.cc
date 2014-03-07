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
//    File   : MenuButton.cc
//    Author : McKay Davis
//    Date   : Mon Oct  2 18:22:42 2006

#include <Core/Skinner/MenuButton.h>
#include <Core/Skinner/Variables.h>

namespace SCIRun {
namespace Skinner {


MenuButton::MenuButton(Variables *vars) :
  Parent(vars),
  x1_(vars,"MenuButton::x1",0),
  x2_(vars,"MenuButton::x2",0),
  y1_(vars,"MenuButton::y1",0),
  y2_(vars,"MenuButton::y2",0)
{
  REGISTER_CATCHER_TARGET(MenuButton::redraw);
}


MenuButton::~MenuButton()
{
}


BaseTool::propagation_state_e
MenuButton::redraw(event_handle_t& /*event*/)
{
  x1_ = get_region().x1();
  x2_ = get_region().x2();
  y1_ = get_region().y1();
  y2_ = get_region().y2();
  return CONTINUE_E;
}


}
}
