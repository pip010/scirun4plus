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
//    File   : CLUTLevelsTool.cc
//    Author : McKay Davis
//    Date   : Sat Oct 14 16:06:12 2006

#include <Applications/Seg3D/CLUTLevelsTool.h>
#include <Applications/Seg3D/Painter.h>

#include <Applications/Seg3D/Seg3DwxGuiUtils.h>

namespace SCIRun {

CLUTLevelsTool::CLUTLevelsTool(Painter *painter) : 
  PointerTool("Color Lookup Table"),
  painter_(painter),
  scale_(1.0), 
  ww_(0), 
  wl_(1.0),
  x_(0),
  y_(0)
{
}


BaseTool::propagation_state_e
CLUTLevelsTool::pointer_down(int b, int x, int y, unsigned int m, int t)
{
  // Find the topmost visible data layer.
  volume_ = 0;
  for (int i = painter_->volumes_.size()-1; i >= 0; i--)
  {
    if (!painter_->volumes_[i]->label_ && painter_->volumes_[i]->visible())
    {
      volume_ = painter_->volumes_[i];
      break;
    }
  }

  if (!volume_.get_rep() || !painter_->cur_window_) {
    return CONTINUE_E;
  }

  ww_ = volume_->clut_max_ - volume_->clut_min_;
  wl_ = volume_->clut_min_ + ww_ / 2.0;
  x_ = x;
  y_ = y;
  
  const double w = painter_->cur_window_->get_region().width();
  const double h = painter_->cur_window_->get_region().height();
  scale_ = (volume_->data_max_ - volume_->data_min_) / sqrt(w*w+h*h);

  return STOP_E;
}


BaseTool::propagation_state_e
CLUTLevelsTool::pointer_motion(int b, int x, int y, unsigned int m, int t)
{
  if (!volume_.get_rep()) { return CONTINUE_E; }

  const float ww = ww_ + scale_ * (y_ - y);  
  const float wl = wl_ + scale_ * (x_ - x);

  volume_->clut_min_ = wl - ww / 2.0;
  volume_->clut_max_ = wl + ww / 2.0;

  std::ostringstream ostrm;

  ostrm << "Window " << ww << "  Level " << wl << "   ";
  ostrm << "Window min " << volume_->clut_min_ << " Window max " << volume_->clut_max_;

  painter_->set_status(ostrm.str());

  volume_->set_slices_dirty();
  painter_->redraw_all();

  return STOP_E;
}


BaseTool::propagation_state_e
CLUTLevelsTool::pointer_up(int b, int x, int y, unsigned int m, int t)
{
  volume_ = 0;
  return QUIT_AND_STOP_E;
}


}  
