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
//    File   : CM2UndoTool.cc
//    Author : Martin Cole
//    Date   : Mon Oct  1 13:59:20 2007

#include <Core/Events/CM2View/CM2UndoTool.h>

namespace SCIRun {


CM2UndoTool::CM2UndoTool(const std::string& name, CM2View* v) :
  BaseTool(name),
  view_(v),
  undo_stack_()
{
}


CM2UndoTool::~CM2UndoTool()
{
}


BaseTool::propagation_state_e
CM2UndoTool::process_event(event_handle_t event)
{
  CM2UndoEvent *ue = dynamic_cast<CM2UndoEvent*>(event.get_rep());
  CM2UndoActionEvent *uae = dynamic_cast<CM2UndoActionEvent*>(event.get_rep());
  if (ue) {
    undo();
    return STOP_E;
  } else if (uae) {
    if (uae->action_ == CM2UndoActionEvent::CLEAR_E) {
      while(!undo_stack_.empty()) {
	UndoItem &item = undo_stack_.top();
	delete item.widget_;
	undo_stack_.pop();
      }
    } else {
      undo_stack_.push(UndoItem(uae->action_, uae->widget_idx_, uae->widget_));
    }
    return STOP_E;
  }
  return CONTINUE_E;
}


void
CM2UndoTool::undo()
{
  if (!undo_stack_.empty())
  {
    const UndoItem &item = undo_stack_.top();
    std::vector<SLIVR::CM2Widget*> &widgets = view_->widgets();

    switch (item.action_)
    {
    case CM2UndoActionEvent::CHANGE_E:
      widgets[item.selected_] = item.widget_;
      break;

    case CM2UndoActionEvent::ADD_E:
      widgets.erase(widgets.begin() + item.selected_);
      break;

    case CM2UndoActionEvent::DELETE_E:
      widgets.insert(widgets.begin() + item.selected_, item.widget_);
      break;
    }
    undo_stack_.pop();
    static std::vector<SLIVR::CM2Widget*> add; //no additional widgets...
    SetWidgetsEvent *e = new SetWidgetsEvent(view_->widgets(), add,
					     view_->parent_id());
    e->updating_ = false;
    event_handle_t event = e;
    EventManager::add_event(event);
  }
}


} // namespace SCIRun
