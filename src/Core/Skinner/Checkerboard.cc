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
//    File   : Checkerboard.cc
//    Author : McKay Davis
//    Date   : Fri Feb  9 12:34:48 2007

#include <Core/Skinner/Checkerboard.h>
#include <Core/Math/MiscMath.h>
#include <sci_gl.h>

namespace SCIRun {


Skinner::Checkerboard::Checkerboard(Variables *vars)
  : Skinner::Drawable(vars),
    color1_(vars, "Checkerboard::color1", Color(1.,1.,1.,1.)),
    color2_(vars, "Checkerboard::color2", Color(0.,0.,0.,1.)),
    width_(vars, "Checkerboard::width", 10.0),
    height_(vars, "Checkerboard::height", 10.0)
{
  REGISTER_CATCHER_TARGET(Checkerboard::redraw);
}


BaseTool::propagation_state_e
Skinner::Checkerboard::redraw(event_handle_t &)
{
  const RectRegion &region = get_region();
  const double width = region.width();
  const double height = region.height();
  const double x1 = region.x1();
  const double y1 = region.y1();
  const double x2 = region.x2();
  const double y2 = region.y2();

  glShadeModel(GL_FLAT);
  glBegin(GL_QUADS);

  const double dx = width_;
  const double dy = height_;
  const double nw = width / dx;
  const double nh = height / dy;

  for (int x = 0; x < nw; ++x) {
    for (int y = 0; y < nh; ++y) {
      if ((x + y) % 2) {
        glColor4dv(&(color1_().r));
      } else {
        glColor4dv(&(color2_().r));
      }
      glVertex3d(x1 + dx * x, y2 - dy * y, 0.0);
      glVertex3d(Clamp(x1 + dx * x + dx, x1, x2),
                 y2 - dy * y, 0.0);
      glVertex3d(Clamp(x1 + dx * x + dx, x1, x2),
                 Clamp(y2 - dy * y - dy, y1, y2), 0.0);
      glVertex3d(x1 + dx * x,
                 Clamp(y2 - dy * y - dy, y1, y2), 0.0);
    }
  }

  glEnd();

  CHECK_OPENGL_ERROR();
  return CONTINUE_E;
}


} // end namespace SCIRun
