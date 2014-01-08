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
//    File   : TransferFunction2D.h
//    Author : McKay Davis
//    Date   : Thu Feb  8 15:33:45 2007

#ifndef SKINNER_TRANSFERFUNCTION2D_H
#define SKINNER_TRANSFERFUNCTION2D_H

#include <Core/Skinner/Parent.h>
#include <Core/Volume/ColorMap2.h>
#include <slivr/ShaderProgramARB.h>
#include <slivr/CM2Shader.h>

namespace SCIRun {
class TextureObj;

namespace Skinner {


class TransferFunction2D : public Parent {
public:
  TransferFunction2D (Variables *variables);
  virtual ~TransferFunction2D();
  CatcherFunction_t                 process_event;

private:
  CatcherFunction_t                 do_PointerEvent;
  CatcherFunction_t                 create_rectangle;
  CatcherFunction_t                 create_triangle;
  CatcherFunction_t                 delete_selected_widget;
  CatcherFunction_t                 enable_disable_selected_widget;
  CatcherFunction_t                 redraw;

  TextureObj *                      texture_;
  Var<Color>                        color_;
  Var<bool>                         enabled_;
  SLIVR::CM2ShaderFactory*          shader_factory_;
  ColorMap2Handle                   colormap2_;
  int                               pick1_;
  int                               pick2_;
  int                               widget_;
};


}
}

#endif
