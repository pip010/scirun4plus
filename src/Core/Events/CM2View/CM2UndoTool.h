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
//    File   : CM2UndoTool.h
//    Author : Martin Cole
//    Date   : Mon Oct  1 13:53:40 2007

#if !defined(CM2UndoTool_h)
#define CM2UndoTool_h

#include <Core/Events/CM2View/CM2View.h>
#include <Core/Events/BaseTool.h>
#include <Core/Events/CM2View/share.h>
#include <stack>

namespace SCIRun {

//! undo last widget change.
class SCISHARE CM2UndoEvent : public BaseEvent
{
public:
  CM2UndoEvent(const std::string& target = "") :
    BaseEvent(target)
  {}

  CM2UndoEvent* clone() { return new CM2UndoEvent(*this); }
};

//! record an undo action
class SCISHARE CM2UndoActionEvent : public BaseEvent
{
public:
  enum action_e { CHANGE_E, ADD_E, DELETE_E, CLEAR_E};
  CM2UndoActionEvent(action_e a, int idx, SLIVR::CM2Widget *w,
    const std::string& target = "") :
    BaseEvent(target),
    action_(a),
    widget_idx_(idx),
    widget_(w)
  {}

  CM2UndoActionEvent* clone() { return new CM2UndoActionEvent(*this); }

  action_e          action_;
  int               widget_idx_;
  SLIVR::CM2Widget *widget_;
};

class SCISHARE CM2UndoTool : public BaseTool
{
public:
  struct UndoItem
  {

    int				action_;
    int				selected_;
    SLIVR::CM2Widget*		widget_;
    UndoItem(int a, int s, SLIVR::CM2Widget* w)
      : action_(a), selected_(s), widget_(w) {}
  };
  CM2UndoTool(const std::string& name, CM2View *v);
  virtual ~CM2UndoTool();


  virtual propagation_state_e process_event(event_handle_t event);

protected:
  virtual void              undo();

  CM2View                   *view_;
  std::stack<UndoItem>	     undo_stack_;
};

} // namespace SCIRun

#endif //CM2UndoTool_h
