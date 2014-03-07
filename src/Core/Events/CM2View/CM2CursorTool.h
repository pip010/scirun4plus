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
//    File   : CM2CursorTool.h
//    Author : Martin Cole
//    Date   : Thu Sep 20 16:19:46 2007

#if !defined(CM2CursorTool_h)
#define CM2CursorTool_h

#include <Core/Events/CM2View/CM2View.h>
#include <Core/Events/BaseTool.h>
#include <Core/Events/CM2View/share.h>

namespace SCIRun {

class SCISHARE CM2CursorTool : public PointerTool
{
public:
  CM2CursorTool(const std::string& name, CM2View *v);
  virtual ~CM2CursorTool();


  //! which == button number, x,y in window at event time.
  //! Spawns cursor change events depending on where the cursor is relative
  //! to the widgets.
  virtual propagation_state_e pointer_motion(int which, int x, int y,
                                             unsigned int modifiers,
                                             int time);
protected:
  virtual void widget_pick(int x, int y, bool right_button_down);
  CM2View                 *view_;
  int                      widget_id_; //! id of widget under cursor.
  int                      pick_id_;   //! id of picked object within widget.
};

} // namespace SCIRun

#endif //CM2CursorTool_h
