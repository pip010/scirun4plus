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
//    File   : CM2WidgetTransformTool.cc
//    Author : Martin Cole
//    Date   : Thu Sep 20 16:22:16 2007

#include <Core/Events/CM2View/CM2WidgetTransformTool.h>
#include <Core/Events/CM2View/CM2UndoTool.h>

namespace SCIRun {


CM2WidgetTransformTool::CM2WidgetTransformTool(const std::string& name, CM2View* v) :
  CM2CursorTool(name, v),
  first_motion_(false)
{
}


CM2WidgetTransformTool::~CM2WidgetTransformTool()
{
}


BaseTool::propagation_state_e
CM2WidgetTransformTool::pointer_down(int which, int x, int y,
				     unsigned int mod, int /*time*/)
{
  if (which != 1 || mod & EventModifiers::SHIFT_E) return CONTINUE_E;
  cached_interval_ = view_->get_timer_interval();
  view_->set_timer_interval(30);
  widget_pick(x, y, false);

  if (widget_id_ != -1) {
    const std::vector<SLIVR::CM2Widget*> &widgets = view_->widgets();
    for (size_t i = 0; i < widgets.size(); ++i) {
      widgets[i]->select(i == (size_t)widget_id_ ? pick_id_ : 0);
      view_->widget_changed_notify();
    }
  }
  first_motion_ = true;
  return STOP_E;
}


BaseTool::propagation_state_e
CM2WidgetTransformTool::pointer_motion(int which, int x, int y,
                                  unsigned int mod, int /*time*/)
{
  if (which != 1 || mod & EventModifiers::SHIFT_E) return CONTINUE_E;

  const int w = view_->width();
  const int h = view_->height();

  const std::vector<SLIVR::CM2Widget*> &widgets = view_->widgets();

  if (widget_id_ != -1) {
    if (first_motion_) {
      // record the widget state.
      event_handle_t e = new CM2UndoActionEvent(CM2UndoActionEvent::CHANGE_E,
					     widget_id_,
					     widgets[widget_id_]->duplicate(),
					     view_->id());
      EventManager::add_event(e);
      first_motion_ = false;
    }

    widgets[widget_id_]->move(x, h - 1 - y, w, h);
    view_->widget_changed_notify();
  }

  return STOP_E;
}


BaseTool::propagation_state_e
CM2WidgetTransformTool::pointer_up(int which, int /*x*/, int /*y*/,
				   unsigned int mod, int /*time*/)
{
  if (which != 1 || mod & EventModifiers::SHIFT_E) return CONTINUE_E;
  view_->set_timer_interval(cached_interval_);
  view_->widget_changed_notify(false);
  return STOP_E;
}


} // namespace SCIRun
