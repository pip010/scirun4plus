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
//    File   : ColorPicker.h
//    Author : McKay Davis
//    Date   : Thu Feb  8 21:22:15 2007


#ifndef SKINNER_COLORPICKER_H
#define SKINNER_COLORPICKER_H

#include <Core/Skinner/Drawable.h>
#include <Core/Datatypes/Color.h>

namespace SCIRun {
namespace Skinner {

class ColorPicker : public Drawable {
public:
  ColorPicker(Variables *variables);

private:
  CatcherFunction_t            redraw;
  CatcherFunction_t            do_PointerEvent;
  CatcherFunction_t            process_event;
  CatcherFunction_t            hsv_changed;
  CatcherFunction_t            rgba_changed;

  void                         draw_glyph(double, double, double);
  void                         update_from_cache();

  Var<Color>                   color_;
  Var<double>                  h_;
  Var<double>                  s_;
  Var<double>                  v_;
  Var<double>                  r_;
  Var<double>                  g_;
  Var<double>                  b_;
  Var<double>                  a_;
  bool                         modifying_hue_;
  bool                         modifying_satval_;
  SCIRun::Color                cached_rgb_;
  SCIRun::HSVColor             cached_hsv_;
  SCIRun::Skinner::Color       cached_color_;
};


}
} // End namespace SCIRun

#endif
