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
//    File   : CM2CursorTool.cc
//    Author : Martin Cole
//    Date   : Thu Sep 20 16:22:16 2007

#include <Core/Events/CM2View/CM2CursorTool.h>
#include <Core/Events/CursorChangeEvent.h>

namespace SCIRun {

CM2CursorTool::CM2CursorTool(const std::string& name, CM2View* v) :
  BaseTool(name),
  PointerTool(name),
  view_(v)
{
}


CM2CursorTool::~CM2CursorTool()
{
}


BaseTool::propagation_state_e
CM2CursorTool::pointer_motion(int which, int x, int y,
                              unsigned int /*mod*/, int /*time*/)
{
  if (which != 0) return CONTINUE_E;

  const int old_wid_id = widget_id_;
  const int old_pick_id = pick_id_;
  widget_pick(x,y, false);

  const std::vector<SLIVR::CM2Widget*> &widgets = view_->widgets();

  if (old_wid_id != widget_id_ || old_pick_id != pick_id_)
  {

    std::string cstr("crosshair");
    if (widget_id_ != -1) {
      //!\todo {cursor changes should be abstracted like other gui events}
      cstr = widgets[widget_id_]->tk_cursorname(pick_id_);
    }
    CursorChangeEvent *cce = new CursorChangeEvent(cstr, view_->parent_id());
    event_handle_t event = cce;
    EventManager::add_event(event);
  }
  return STOP_E;
}


void
CM2CursorTool::widget_pick(int x, int y, bool right_button_down)
{
  const std::vector<SLIVR::CM2Widget*> &widgets = view_->widgets();
  const int w = view_->width();
  const int h = view_->height();
  if (!right_button_down)
    for (widget_id_ = widgets.size() - 1; widget_id_>=0; widget_id_--) {
      if (widgets[widget_id_]->get_onState() &&
	  (pick_id_ = widgets[widget_id_]->pick1(x, h - 1 - y, w, h))) {
	break; //found
      }
    }

  if (!pick_id_) {
    for (widget_id_ = widgets.size() - 1; widget_id_>=0; widget_id_--) {
      if (widgets[widget_id_]->get_onState() &&
	  (pick_id_ = widgets[widget_id_]->pick2(x, h - 1 - y, w, h,
						 right_button_down))) {
	break; //found
      }
    }
  }
}


} // namespace SCIRun
