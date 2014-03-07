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
//    File   : CM2WidgetTransformTool.h
//    Author : Martin Cole
//    Date   : Thu Sep 20 16:19:46 2007

#if !defined(CM2WidgetTransformTool_h)
#define CM2WidgetTransformTool_h

#include <Core/Events/CM2View/CM2View.h>
#include <Core/Events/CM2View/CM2CursorTool.h>
#include <Core/Events/CM2View/share.h>

namespace SCIRun {

class SCISHARE CM2WidgetTransformTool : public CM2CursorTool
{
public:
  CM2WidgetTransformTool(const std::string& name, CM2View *v);
  virtual ~CM2WidgetTransformTool();


  //! which == button number, x,y in window at event time
  virtual propagation_state_e pointer_down(int which, int x, int y,
                                           unsigned int modifiers,
                                           int time);
  //! which == button number, x,y in window at event time
  virtual propagation_state_e pointer_motion(int which, int x, int y,
                                             unsigned int modifiers,
                                             int time);
  //! which == button number, x,y in window at event time
  virtual propagation_state_e pointer_up(int which, int x, int y,
                                         unsigned int modifiers,
                                         int time);
private:
  void widget_changed_notify(bool updating = true);

  bool first_motion_;
  int  cached_interval_;
};

} // namespace SCIRun

#endif //CM2WidgetTransformTool_h
