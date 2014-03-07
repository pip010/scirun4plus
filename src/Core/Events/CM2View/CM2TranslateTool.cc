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
//    File   : CM2TranslateTool.cc
//    Author : Martin Cole
//    Date   : Thu Sep 20 16:22:16 2007

#include <Core/Events/CM2View/CM2TranslateTool.h>

namespace SCIRun {


CM2TranslateTool::CM2TranslateTool(const std::string& name, CM2View* v) :
  BaseTool(name),
  PointerTool(name),
  view_(v),
  start_x_(-1),
  start_y_(-1)
{
}


CM2TranslateTool::~CM2TranslateTool()
{
}


BaseTool::propagation_state_e
CM2TranslateTool::pointer_down(int which, int x, int y,
			       unsigned int mod, int /*time*/)
{
  // shift button 2 is reset.
  if (which == 2 && mod & EventModifiers::SHIFT_E) {
    view_->set_scale(1.0);
    view_->set_pan_x(0.0);
    view_->set_pan_y(0.0);
    return STOP_E;
  }
  if (!(which == 1 && mod & EventModifiers::SHIFT_E)) return CONTINUE_E;
  start_x_ = x;
  start_y_ = y;
  spx_ = view_->pan_x();
  spy_ = view_->pan_y();
  return STOP_E;
}


BaseTool::propagation_state_e
CM2TranslateTool::pointer_motion(int which, int x, int y,
                                  unsigned int mod, int /*time*/)
{
  if (!(which == 1 && mod & EventModifiers::SHIFT_E)) return CONTINUE_E;
  const double dx = double(start_x_ - x) / view_->width();
  const double dy = double(y - start_y_) / view_->height();

  view_->set_pan_x(spx_ + dx / view_->scale());
  view_->set_pan_y(spy_ + dy / view_->scale());
  return STOP_E;
}


BaseTool::propagation_state_e
CM2TranslateTool::pointer_up(int which, int /*x*/, int /*y*/,
                              unsigned int mod, int /*time*/)
{
  if (!(which == 1 && mod & EventModifiers::SHIFT_E)) return CONTINUE_E;
  start_x_ = -1;
  start_y_ = -1;
  return STOP_E;
}


} // namespace SCIRun
