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
//    File   : ColorPicker.cc
//    Author : McKay Davis
//    Date   : Thu Feb  8 21:22:27 2007

#include <Core/Util/StringUtil.h>
#include <Core/Skinner/Variables.h>
#include <Core/Skinner/ColorPicker.h>
#include <Core/Events/EventManager.h>
#include <Core/Datatypes/Color.h>
#include <Core/Geometry/Point.h>
#include <Core/Geometry/Ray.h>
#include <Core/Geometry/Vector.h>
#include <Core/Geometry/Plane.h>
#include <Core/Math/MiscMath.h>
#include <Core/Math/MinMax.h>
#include <sci_gl.h>

namespace SCIRun {
namespace Skinner{
using namespace SLIVR;


ColorPicker::ColorPicker(Variables *vars)
  : Drawable(vars),
    color_(vars, "ColorPicker::color", Color(1.0, 1.0, 1.0, 1.0)),
    h_(vars, "ColorPicker::hue", 0.0),
    s_(vars, "ColorPicker::saturation", 0.0),
    v_(vars, "ColorPicker::value", 1.0),
    r_(vars, "ColorPicker::red", 1.0),
    g_(vars, "ColorPicker::green", 1.0),
    b_(vars, "ColorPicker::blue", 1.0),
    a_(vars, "ColorPicker::alpha", 1.0),
    modifying_hue_(0),
    modifying_satval_(0)

{
  REGISTER_CATCHER_TARGET(ColorPicker::rgba_changed);
  REGISTER_CATCHER_TARGET(ColorPicker::hsv_changed);

  REGISTER_CATCHER_TARGET(ColorPicker::redraw);
  REGISTER_CATCHER_TARGET(ColorPicker::do_PointerEvent);
}


BaseTool::propagation_state_e
ColorPicker::redraw(event_handle_t &)
{
  const RectRegion &region = get_region();
  const double width = region.width();
  const double height = region.height();
  const double cx = region.x1() + width/2.0;
  const double cy = region.y1() + height/2.0;
  const double outside_wid = Min(width,height)/2.0;
  const double outside_hei = Min(width,height)/2.0;
  const double inside_wid = 0.8 * outside_wid;
  const double inside_hei = 0.8 * outside_hei;

  glShadeModel(GL_SMOOTH);
  glBegin(GL_QUAD_STRIP);
  const double deg2rad = 2.0 * M_PI / 360.0;
  SCIRun::Color color;
  double rad;
  for (int deg = 0; deg <= 360; deg++) {
    rad = deg * deg2rad;
    color = SCIRun::Color(HSVColor(deg, 1.0, 1.0));
    glColor4d(color.r(), color.g(), color.b(), 1.0);
    glVertex3d(cx + sin(rad) * inside_wid, cy + cos(rad) * inside_hei, 0.0);
    glColor4d(color.r(), color.g(), color.b(), 1.0);
    glVertex3d(cx + sin(rad) * outside_wid, cy + cos(rad) * outside_hei, 0.0);
  }
  glEnd();

  glBegin(GL_TRIANGLES);

  double hrad = h_ * deg2rad;
  Point hpoint(cx + sin(hrad) * inside_wid, cy + cos(hrad) * inside_hei, 0.0);
  double zrad = (h_ + 240) * deg2rad;
  Point zpoint(cx + sin(zrad) * inside_wid, cy + cos(zrad) * inside_hei, 0.0);
  double vrad = (h_ + 120) * deg2rad;
  Point vpoint(cx + sin(vrad) * inside_wid, cy + cos(vrad) * inside_hei, 0.0);

  Vector satvec = hpoint - vpoint;
  Vector valvec = vpoint - zpoint;


  color = SCIRun::Color(HSVColor(h_, 1.0, 1.0));

  glColor4d(color.r(), color.g(), color.b(), 1.0);
  glVertex3dv(&hpoint(0));


  color = SCIRun::Color(HSVColor(h_, 0.0, 0.0));

  glColor4d(color.r(), color.g(), color.b(), 1.0);
  glVertex3dv(&zpoint(0));

  color = SCIRun::Color(HSVColor(h_, 0.0, 1.0));
  glColor4d(color.r(), color.g(), color.b(), 1.0);
  glVertex3dv(&vpoint(0));
  glEnd();


  double mcx = cx + sin(hrad) * 0.9 * outside_wid;
  double mcy = cy + cos(hrad) * 0.9 * outside_hei;
  draw_glyph(mcx, mcy, 3.0);

  double s = s_/100.0;
  double v = v_/100.0;
  Point gpoint = zpoint + valvec*v + satvec*v*s;
  draw_glyph(gpoint.x(), gpoint.y(), 6.0);

  CHECK_OPENGL_ERROR();
  return CONTINUE_E;
}


void
ColorPicker::draw_glyph(double mcx, double mcy, double mr)
{
  const double deg2rad = 2.0 * M_PI / 360.0;
  double rad;
  glBegin(GL_TRIANGLE_FAN);
  glColor4d(1.0, 1.0, 1.0, 1.0);
  glColor4d(0.3, 0.3, 0.3, 1.0);
  glVertex3d(mcx, mcy, 0.0);
  for (int deg = 0; deg <= 360; deg+=36) {
    rad = deg * deg2rad;
    glVertex3d(mcx + cos(rad) * mr, mcy + sin(rad) * mr, 0.0);
  }
  glEnd();

  glBegin(GL_TRIANGLE_FAN);
  glColor4d(0.8, 0.8, 0.8, 0.8);
  glColor4d(1.0, 1.0, 1.0, 1.0);
  glVertex3d(mcx-1.0, mcy+1.0, 0.0);
  glColor4d(0.8, 0.8, 0.8, 0.2);
  glColor4d(1.0, 1.0, 1.0, 0.1);
  for (int deg = 0; deg <= 360; deg+=36) {
    rad = deg * deg2rad;
    glVertex3d(mcx + cos(rad) * mr, mcy + sin(rad) * mr, 0.0);
  }
  glEnd();
}


BaseTool::propagation_state_e
ColorPicker::do_PointerEvent(event_handle_t &event)
{
  PointerSignal *signal = dynamic_cast<PointerSignal *>(event.get_rep());
  PointerEvent *pointer = signal->get_pointer_event();

  bool pressing =
    pointer->get_pointer_state() & PointerEvent::BUTTON_1_E;

  if (!pressing) {
    modifying_hue_ = false;
    modifying_satval_ = false;
    return CONTINUE_E;
  }

  bool pressed =
    pointer->get_pointer_state() & PointerEvent::BUTTON_PRESS_E;

  const int x = pointer->get_x();
  const int y = pointer->get_y();
  const RectRegion &region = get_region();
  const double width = region.width();
  const double height = region.height();
  const double ratio = Min(width, height);
  const double cx = region.x1() + width/2.0;
  const double cy = region.y1() + height/2.0;
  const double xx = (x - cx) * 2.0 / ratio;
  const double yy = (y - cy) * 2.0 / ratio;

  if (pressed && pressing && region.inside(x, y)) {
    double dist = sqrt(xx*xx+yy*yy);
    if (dist >= 0.8) {
      modifying_hue_ = true;
    } else {
      modifying_satval_ = true;
    }
  }

  if (!modifying_hue_ && !modifying_satval_) {
    return CONTINUE_E;
  }


  if (modifying_hue_) {
    Vector v1(0.0, 1.0, 0.0);
    Vector v2(xx, yy, 0.0);
    v2.normalize();

    double d = Dot(v2, v1);
    double ad = acos(d) * 180 / M_PI;

    h_ = (xx > 0.0) ? ad : 360-ad;
  }

  if (modifying_satval_) {
    const double inside_wid = 0.8*Min(width,height)/2.0;
    const double inside_hei = 0.8*Min(width,height)/2.0;

    const double deg2rad = 2.0 * M_PI / 360.0;

    double hrad = h_ * deg2rad;
    Point hpoint(cx + sin(hrad) * inside_wid,
                 cy + cos(hrad) * inside_hei, 0.0);
    double zrad = (h_ + 240) * deg2rad;
    Point zpoint(cx + sin(zrad) * inside_wid,
                 cy + cos(zrad) * inside_hei, 0.0);
    double vrad = (h_ + 120) * deg2rad;
    Point vpoint(cx + sin(vrad) * inside_wid,
                 cy + cos(vrad) * inside_hei, 0.0);

    Point cpoint(x,y,0);
    Vector satvec = hpoint - zpoint;
    Vector valvec = vpoint - zpoint;

    Vector normal = (hpoint + vpoint)/2.0 - zpoint;
    normal = -normal;

    Ray vray(zpoint, valvec);
    Ray sray(zpoint, satvec);

    double t1, t2;
    vray.planeIntersectParameter(normal, cpoint, t1);
    sray.planeIntersectParameter(normal, cpoint, t2);

    Point vproj = (zpoint + valvec * t1);

    Point sproj = (zpoint + satvec * t2);

    Vector opp = sproj - vproj;

    v_ = t1 * 100;
    if (Plane(vproj, opp).eval_point(cpoint) > 0.0) {
      s_ = ((cpoint - vproj).length() / opp.length()) * 100.0;
    } else {
      s_ = 0.0;
    }
  }

  hsv_changed(event);

  EventManager::add_event(new WindowEvent(WindowEvent::REDRAW_E));

  return CONTINUE_E;

}


void
ColorPicker::update_from_cache()
{
  h_ = Clamp(cached_hsv_.hue(), 0.0, 360.0);
  s_ = Clamp(cached_hsv_.sat()*100.0, 0.0, 100.0);
  v_ = Clamp(cached_hsv_.val()*100.0, 0.0, 100.0);

  r_ = Clamp(cached_rgb_.r()*255, 0.0, 255.0);
  g_ = Clamp(cached_rgb_.g()*255, 0.0, 255.0);
  b_ = Clamp(cached_rgb_.b()*255, 0.0, 255.0);

  color_ = cached_color_;
}


BaseTool::propagation_state_e
ColorPicker::rgba_changed(event_handle_t& /*event*/)
{
  cached_rgb_ = SCIRun::Color(r_/255.0, g_/255.0, b_/255.0);
  cached_hsv_ = SCIRun::HSVColor(cached_rgb_);
  cached_color_ = Color(cached_rgb_.r(),cached_rgb_.g(),cached_rgb_.b(), a_);

  update_from_cache();

  return CONTINUE_E;
}


BaseTool::propagation_state_e
ColorPicker::hsv_changed(event_handle_t& /*event*/)
{
  cached_hsv_ = SCIRun::HSVColor(h_, s_/100.0, v_/100.0);
  cached_rgb_ = SCIRun::Color(cached_hsv_);
  cached_color_ = Color(cached_rgb_.r(),cached_rgb_.g(),cached_rgb_.b(), a_);

  update_from_cache();

  return CONTINUE_E;
}



BaseTool::propagation_state_e
ColorPicker::process_event(event_handle_t &event)
{
 if (Abs(color_().r - cached_color_.r) > 0.001 ||
     Abs(color_().g - cached_color_.g) > 0.001 ||
     Abs(color_().b - cached_color_.b) > 0.001 ||
     Abs(color_().a - cached_color_.a) > 0.001)
 {
   cached_color_ = color_;
   a_ = color_().a;
   cached_rgb_ =
     SCIRun::Color(cached_color_.r,cached_color_.g,cached_color_.b);
   cached_hsv_ = SCIRun::HSVColor(cached_rgb_);
   update_from_cache();
 }
 else if (Abs(h_       - cached_hsv_.hue()) > 0.01 ||
          Abs(s_/100.0 - cached_hsv_.sat()) > 0.001 ||
          Abs(v_/100.0 - cached_hsv_.val()) > 0.001)
 {
   hsv_changed(event);
 } else if (Abs(r_/255.0 - cached_rgb_.r()) > 0.001 ||
            Abs(g_/255.0 - cached_rgb_.g()) > 0.001 ||
            Abs(b_/255.0 - cached_rgb_.b()) > 0.001 ||
            Abs(a_ - color_().a) > 0.001)
 {
   rgba_changed(event);
 }

  return Drawable::process_event(event);
}


}
} // end namespace SCIRun
