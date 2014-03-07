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
//    File   : KeyToolSelectorTool.cc
//    Author : McKay Davis
//    Date   : Sun Oct 15 11:43:45 2006

#include <Applications/Seg3D/KeyToolSelectorTool.h>
#include <Applications/Seg3D/Painter.h>
#include <Applications/Seg3D/AutoviewTool.h>
#include <Applications/Seg3D/BrushTool.h>
#include <Applications/Seg3D/CLUTLevelsTool.h>
#include <Applications/Seg3D/FloodfillTool.h>
#include <Applications/Seg3D/PanTool.h>
#include <Applications/Seg3D/ProbeTool.h>
#include <Applications/Seg3D/StatisticsTool.h>
#include <Applications/Seg3D/ZoomTool.h>

namespace SCIRun {


SliceWindowKeyToolSelectorTool::SliceWindowKeyToolSelectorTool(Painter *p) :
  KeyTool("Painter SliceWindowKeyToolSelectorTool"),
  painter_(p),
  tm_(p->tm_)
{
}
  

SliceWindowKeyToolSelectorTool::~SliceWindowKeyToolSelectorTool()
{
}


BaseTool::propagation_state_e
SliceWindowKeyToolSelectorTool::key_press(string, int keyval,
                                          unsigned int m, unsigned int)
{
  if (!painter_->cur_window_) return STOP_E;

  if (sci_getenv_p("SCI_DEBUG"))
    cerr << "keyval: " << keyval << std::endl;

  SliceWindow &window = *painter_->cur_window_;
  event_handle_t event;
  switch (keyval)
  {
  case SCIRun_equal:    window.zoom_in(event); break;
  case SCIRun_minus:    window.zoom_out(event); break;
  case SCIRun_comma:    window.move_slice(-1); break;
  case SCIRun_period:   window.move_slice(1); break;
  case SCIRun_p:        window.punch_current_slice(event); break;
  case SCIRun_e:        window.erase_current_slice(event); break;

  case SCIRun_c:
    Painter::ThrowSkinnerSignal("Painter::SetLayer", false);
    break;

  case SCIRun_space:
    painter_->toggle_current_volume_visibility();
    break;
    
  case SCIRun_Up:
    if (m & EventModifiers::SHIFT_E)
    {
      window.copy_current_slice_up(event);
    }
    window.move_slice(1);
    break;

  case SCIRun_Down:
    if (m & EventModifiers::SHIFT_E)
    {
      window.copy_current_slice_down(event);
    }
    window.move_slice(-1);
    break;
  }

  return CONTINUE_E;
}  


GlobalKeyToolSelectorTool::GlobalKeyToolSelectorTool(Painter *p) :
  KeyTool("Painter GlobalKeyToolSelectorTool"),
  painter_(p),
  tm_(p->tm_)
{
}
  

GlobalKeyToolSelectorTool::~GlobalKeyToolSelectorTool()
{
}


BaseTool::propagation_state_e
GlobalKeyToolSelectorTool::key_press(string, int keyval,
                                          unsigned int m, unsigned int)
{
  if (sci_getenv_p("SCI_DEBUG"))
    cerr << "keyval: " << keyval << std::endl;

  event_handle_t event;
  switch (keyval) {
  case SCIRun_Left:
    if (m & EventModifiers::SHIFT_E)
    {
      painter_->move_layer_down();
    }
    else
    {
      painter_->current_layer_down();
    }
    return STOP_E;

  case SCIRun_Right:
    if (m & EventModifiers::SHIFT_E)
    {
      painter_->move_layer_up();
    }
    else
    {
      painter_->current_layer_up();
    }
    return STOP_E;
  }

#if 0
  switch (keyval)
  {
  case SCIRun_c:        painter_->CopyLayer(event); break;
  case SCIRun_x:        painter_->DeleteLayer(event); break;
  case SCIRun_v:        painter_->NewLayer(event);break;
  case SCIRun_f:        
  if (painter_->cur_window_) {
    painter_->cur_window_->autoview(painter_->current_volume_);
  } break;

  // Reset CLUT
  case SCIRun_r: { 
    if (painter_->current_volume_.get_rep()) {
      painter_->current_volume_->reset_clut();
    } 
  } break;

  case SCIRun_p:        painter_->opacity_up();break;
  case SCIRun_o:        painter_->opacity_down();break;

  case SCIRun_u:
    if (painter_->current_volume_.get_rep()) {
      painter_->current_volume_->colormap_ = 
        Max(0,painter_->current_volume_->colormap_-1);
      painter_->set_all_slices_tex_dirty();
      painter_->redraw_all();
    } break;
  case SCIRun_i:
    if (painter_->current_volume_.get_rep()) {
      painter_->current_volume_->colormap_ = 
        Min(int(painter_->colormaps_.size()), 
            painter_->current_volume_->colormap_+1);
      painter_->set_all_slices_tex_dirty();
      painter_->redraw_all();
    } break;    
  }
#endif

  return CONTINUE_E;
}  


} // End namespace SCIRun
