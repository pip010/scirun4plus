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
//    File   : ViewScaleTool.cc
//    Author : Martin Cole
//    Date   : Tue Jun  6 12:31:54 2006

#include <Core/Events/OGLView/ViewScaleTool.h>
#include <Core/Math/MiscMath.h>
#include <sstream>

namespace SCIRun {


ViewScaleTool::ViewScaleTool(const std::string& name, ViewToolInterface* i) :
  BaseTool(name),
  PointerTool(name),
  scene_interface_(i),
  start_x_(-1),
  start_y_(-1),
  start_view_()
{
}


ViewScaleTool::~ViewScaleTool()
{
}


BaseTool::propagation_state_e
ViewScaleTool::pointer_down(int which, int x, int y, unsigned int, int /*time*/)
{
  if (which != 3) return CONTINUE_E;
  scene_interface_->update_mode_string("scale: ");
  start_x_ = x;
  start_y_ = y;
  start_view_ = scene_interface_->view_;
  return STOP_E;
}


BaseTool::propagation_state_e
ViewScaleTool::pointer_motion(int which, int x, int y, unsigned int, int /*time*/)
{
  if (which != 3 || start_x_ < 0 || start_y_ < 0) return CONTINUE_E;
  const int delta = start_y_ + start_x_ - x - y;
  double scale = Max(0.001, Pow(1.01, delta));

  scene_interface_->view_.eyep
    (start_view_.lookat() + scale*(start_view_.eyep()-start_view_.lookat()));
  scene_interface_->need_redraw();
  std::ostringstream str;
  str << "scale: " << 100.0 / scale << "%";
  scene_interface_->update_mode_string(str.str());
  return STOP_E;
}


BaseTool::propagation_state_e
ViewScaleTool::pointer_up(int which, int /*x*/, int /*y*/, unsigned int, int /*time*/)
{
  if (which != 3 || start_x_ < 0 || start_y_ < 0) return CONTINUE_E;
  start_x_ = -1;
  start_y_ = -1;
  scene_interface_->update_mode_string("");
  scene_interface_->need_redraw();
  return STOP_E;
}


} // namespace SCIRun
