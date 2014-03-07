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
//    File   : Slider.h
//    Author : McKay Davis
//    Date   : Tue Feb  6 17:26:53 2007

#ifndef SKINNER_SLIDER_H
#define SKINNER_SLIDER_H

#include <Core/Skinner/Parent.h>
#include <map>
#include <set>
#include <vector>

namespace SCIRun {
namespace Skinner {

class SliderButton;


class Slider : public Parent {
public:
  Slider (Variables *variables);
  virtual ~Slider();
  virtual int                       get_signal_id(const string &) const;
  CatcherFunction_t                 process_event;

protected:
  friend class SliderButton;
  RectRegion                        get_subregion();
  CatcherFunction_t                 SliderButton_Maker;
  CatcherFunction_t                 do_PointerEvent;
  bool                              update_value();

  SliderButton *                    button_;
  Var<double>                       min_;
  Var<double>                       max_;
  Var<double>                       resolution_;
  Var<double>                       value_;
  Var<bool>                         vertical_;
  int                               sx_;
  int                               sy_;
  bool                              drag_;
  double                            drag_value_;
};


}
}

#endif
