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
//    File   : TransferFunction2D.cc
//    Author : McKay Davis
//    Date   : Thu Feb  8 15:38:33 2007

#include <Core/Skinner/TransferFunction2D.h>
#include <Core/Events/EventManager.h>
//#include <Core/Algorithms/Visualization/NrrdTextureBuilderAlgo.h>
#include <Core/Geom/TextureObj.h>
#include <Core/Events/SceneGraphEvent.h>
#include <Core/Events/ColorMap2Event.h>
#include <Core/Util/Environment.h>
#include <Core/Util/FileUtils.h>
#include <Core/Volume/ColorMap2.h>
#include <Core/Volume/CM2Widget.h>
#include <slivr/CM2Shader.h>
#include <slivr/TextureRenderer.h>
#include <slivr/VolumeRenderer.h>
#include <slivr/ColorMap.h>
#include <Core/Math/MiscMath.h>

namespace SCIRun {
namespace Skinner {
using namespace SLIVR;


TransferFunction2D::TransferFunction2D(Variables *vars) :
  Parent(vars),
  texture_(0),
  color_(vars, "TransferFunction2D::color", Color(1., 1., 1., 1.)),
  enabled_(vars, "TransferFunction2D::enabled", 0),
  shader_factory_(0),
  colormap2_(new ColorMap2()),
  pick1_(0),
  pick2_(0),
  widget_(-1)
{
  REGISTER_CATCHER_TARGET(TransferFunction2D::do_PointerEvent);
  REGISTER_CATCHER_TARGET(TransferFunction2D::redraw);
  REGISTER_CATCHER_TARGET(TransferFunction2D::create_rectangle);
  REGISTER_CATCHER_TARGET(TransferFunction2D::create_triangle);
  REGISTER_CATCHER_TARGET(TransferFunction2D::delete_selected_widget);
  REGISTER_CATCHER_TARGET(TransferFunction2D::enable_disable_selected_widget);

  NrrdDataHandle nrrdh = new NrrdData();

  size_t sizes[NRRD_DIM_MAX];
  sizes[0] = 4;
  sizes[1] = 256;
  sizes[2] = 64;

  nrrdAlloc_nva(nrrdh->nrrd_, nrrdTypeFloat, 3, sizes);
  memset(nrrdh->nrrd_->data, 0, sizes[0]*sizes[1]*sizes[2]);

  texture_ = new TextureObj(nrrdh);
}


TransferFunction2D::~TransferFunction2D()
{
}


void
splat(Nrrd *nrrd, double radius, int x0, int y0,
      double r=1, double g=1, double b=1, double a=1)
{
  //  if (splatcount++ % splatmod) return;
  vector<int> index(3,0);
  const unsigned int wid = Round(radius);
  const unsigned int datawid = (unsigned int)(nrrd->axis[1].size)*4;
  float *data = (float *)nrrd->data;
  for (int y = y0-wid; y <= int(y0+wid); ++y)
    for (int x = x0-wid; x <= int(x0+wid); ++x)
      if (x >= 0 && x < int(nrrd->axis[1].size) &&
          y >= 0 && y < int(nrrd->axis[2].size))
      {
        index[1] = x;
        index[2] = y;

        float dist = sqrt(double(x0-x)*(x0-x)+(y0-y)*(y0-y))/radius;
        if (dist <= 1.0) {
          data[y * datawid + x * 4 + 0] = r;
          data[y * datawid + x * 4 + 1] = g;
          data[y * datawid + x * 4 + 2] = b;
          data[y * datawid + x * 4 + 3] = a;
        }
      }
}


void
line(Nrrd *nrrd, double radius,
     int x0, int y0, int x1, int y1, bool first)
{
  //  splatmod = Ceil(radius*0.3);

  if (x0 < 0 || x0 >= (int) nrrd->axis[1].size ||
      x1 < 0 || x1 >= (int) nrrd->axis[1].size ||
      y0 < 0 || y0 >= (int) nrrd->axis[2].size ||
      y1 < 0 || y1 >= (int) nrrd->axis[2].size) return;
  int dy = y1 - y0;
  int dx = x1 - x0;
  int sx = 1;
  int sy = 1;
  int frac = 0;
  if (dy < 0) {
    dy = -dy;
    sy = -1;
  }
  if (dx < 0) {
    dx = -dx;
    sx = -1;
  }
  dy <<= 1;
  dx <<= 1;
  if (first) splat(nrrd, radius, x0, y0);
  if (dx > dy) {
    frac = dy - (dx >> 1);
    while (x0 != x1) {
      if (frac >= 0) {
        y0 += sy;
        frac -= dx;
      }
      x0 += sx;
      frac += dy;
      splat(nrrd, radius, x0, y0);
    }
  } else {
    frac = dx - (dy >> 1);
    while (y0 != y1) {
      if (frac >= 0) {
        x0 += sx;
        frac -= dy;
      }
      y0 += sy;
      frac += dx;
      splat(nrrd, radius, x0, y0);
    }
  }
}


BaseTool::propagation_state_e
TransferFunction2D::process_event(event_handle_t &event)
{
  ColorMap2Event *cm2e = dynamic_cast<ColorMap2Event *>(event.get_rep());
  if (cm2e) {
    const bool different = (colormap2_ != cm2e->get_data().back());

    enabled_ = false;
    colormap2_ = cm2e->get_data().back();
    
    if (different)
    {
      if (!colormap2_.get_rep())
      {
        widget_ = -1;
        pick1_ = 0;
        pick2_ = 0;
      }
      else
      {
        // Update the widget pointer so that one is always selected.
        widget_ = 0;
        pick1_ = 0;
        pick2_ = 0;

        for (int i = 0; i < (int)colormap2_->widgets().size(); i++)
        {
          //int pick = pick1_ ? pick1_ : pick2_;
          colormap2_->widgets()[i]->select(i == widget_);
          if (i == widget_) {
            Skinner::Color newcolor;
            newcolor.r = colormap2_->widgets()[i]->get_color().r();
            newcolor.g = colormap2_->widgets()[i]->get_color().g();
            newcolor.b = colormap2_->widgets()[i]->get_color().b();
            newcolor.a = colormap2_->widgets()[i]->get_alpha();
            color_ = newcolor;
          }
        }
      }
    }
    EventManager::add_event(new WindowEvent(WindowEvent::REDRAW_E));
  }
  else if (widget_ >= 0)
  {
    SLIVR::CM2Widget *widget = colormap2_->widgets()[widget_];

    enabled_ = widget->get_onState();

    SLIVR::Color wcolor = widget->get_color();
    Color oldcolor(wcolor.r(), wcolor.g(), wcolor.b(),widget->get_alpha());

    if ((Abs(oldcolor.r - color_().r) > 0.001) ||
        (Abs(oldcolor.g - color_().g) > 0.001) ||
        (Abs(oldcolor.b - color_().b) > 0.001) ||
        (Abs(oldcolor.a - color_().a) > 0.001)) {
      widget->set_color(SLIVR::Color(color_().r, color_().g, color_().b));
      widget->set_alpha(color_().a);

      ColorMap2Event *cm2e = new ColorMap2Event(colormap2_);
      EventManager::add_event(cm2e);
    }
  }

  return Parent::process_event(event);
}


BaseTool::propagation_state_e
TransferFunction2D::create_rectangle(event_handle_t& /*event*/)
{
  RectangleCM2Widget *widget = new RectangleCM2Widget();
  widget->set_faux(false);
  colormap2_->widgets().push_back(widget);
  ColorMap2Event *cm2e = new ColorMap2Event(colormap2_);
  EventManager::add_event(cm2e);
  return CONTINUE_E;
}



BaseTool::propagation_state_e
TransferFunction2D::create_triangle(event_handle_t& /*event*/)
{
  TriangleCM2Widget *widget = new TriangleCM2Widget();
  widget->set_faux(false);
  colormap2_->widgets().push_back(widget);
  ColorMap2Event *cm2e = new ColorMap2Event(colormap2_);
  EventManager::add_event(cm2e);
  return CONTINUE_E;
}


BaseTool::propagation_state_e
TransferFunction2D::delete_selected_widget(event_handle_t& /*event*/)
{
  vector<SLIVR::CM2Widget*> &widgets = colormap2_->widgets();

  if (widget_ < 0 || widget_ >= (int)widgets.size()) return STOP_E;

  widgets.erase(widgets.begin() + widget_);

  widget_ = -1;

  ColorMap2Event *cm2e = new ColorMap2Event(colormap2_);
  EventManager::add_event(cm2e);

  return CONTINUE_E;
}


BaseTool::propagation_state_e
TransferFunction2D::enable_disable_selected_widget(event_handle_t& /*event*/)
{
  SLIVR::CM2Widget* widget = colormap2_->widgets()[widget_];

  widget->set_onState(!widget->get_onState());
  enabled_ = !enabled_;

  ColorMap2Event *cm2e = new ColorMap2Event(colormap2_);
  EventManager::add_event(cm2e);

  return CONTINUE_E;
}



BaseTool::propagation_state_e
TransferFunction2D::do_PointerEvent(event_handle_t &event)
{
  PointerSignal *signal = dynamic_cast<PointerSignal *>(event.get_rep());
  PointerEvent *pointer = signal->get_pointer_event();

  bool button1 =
    pointer->get_pointer_state() & PointerEvent::BUTTON_1_E;

  bool pressed =
    pointer->get_pointer_state() & PointerEvent::BUTTON_PRESS_E;

  bool released =
    pointer->get_pointer_state() & PointerEvent::BUTTON_RELEASE_E;

  const RectRegion &region = get_region();
  int x = pointer->get_x();
  int y = pointer->get_y();

  int xx = Round(x - region.x1());
  int yy = Round(y - region.y1());
  int wid = Round(region.width());
  int hei = Round(region.height());

  vector<SLIVR::CM2Widget*> &widgets = colormap2_->widgets();

  if (button1 && pressed && region.inside(x, y))
  {
    for (int i = widgets.size()-1; i >= 0; i--) {
      int pick = widgets[i]->pick1(xx, yy, wid, hei);
      if (pick) {
        widget_ = i;
        pick1_ = pick;
        break;
      }
    }

    if (!pick1_) {
      for (int i = widgets.size()-1; i >= 0; i--) {
        int pick = widgets[i]->pick2(xx, yy, wid, hei, 0);
        if (pick) {
          widget_ = i;
          pick2_ = pick;
          break;
        }
      }
    }

    for (int i = 0; i < (int)widgets.size(); i++) {
      int pick = pick1_ ? pick1_ : pick2_;
      widgets[i]->select(i == widget_ ? pick : 0);
      if (i == widget_) {
        Skinner::Color newcolor;
        newcolor.r = widgets[i]->get_color().r();
        newcolor.g = widgets[i]->get_color().g();
        newcolor.b = widgets[i]->get_color().b();
        newcolor.a = widgets[i]->get_alpha();
        color_ = newcolor;
      }
    }
  }

  if (button1 && (pick1_ || pick2_)) {
    widgets[widget_]->move(xx,yy,wid,hei);
    EventManager::add_event(new WindowEvent(WindowEvent::REDRAW_E));
    ColorMap2Event *cm2e = new ColorMap2Event(colormap2_);
    EventManager::add_event(cm2e);
  }

  if (button1 && released) {
    pick1_ = 0;
    pick2_ = 0;
    EventManager::add_event(new WindowEvent(WindowEvent::REDRAW_E));
  }

  return CONTINUE_E;
}


BaseTool::propagation_state_e
TransferFunction2D::redraw(event_handle_t &)
{
  if (!shader_factory_) {
    SLIVR::ShaderProgramARB::init_shaders_supported();
    if (!SLIVR::ShaderProgramARB::shaders_supported()) {
      cerr << "Cannot init_shader_factory.\n"
           << " 2D widget colors will not be drawn\n";
      shader_factory_ = 0;
    } else {
      shader_factory_ = new SLIVR::CM2ShaderFactory();
    }
  }

  glMatrixMode(GL_TEXTURE);
  glPushMatrix();
  glLoadIdentity();
  //CHECK_OPENGL_ERROR( "TransferFunction2D glLoadIdentity" );

  Point coords[4];
  const RectRegion &region = get_region();


  float tw = 1.0;
  float th = 1.0;

  float tex_coords[8] = {0, Floor(th)-th,
                         tw, Floor(th)-th,
                         tw, Floor(th),
                         0.0, Floor(th) };

  coords[0] = Point(region.x1(), region.y1(), 0.0);
  coords[1] = Point(region.x2(), region.y1(), 0.0);
  coords[2] = Point(region.x2(), region.y2(), 0.0);
  coords[3] = Point(region.x1(), region.y2(), 0.0);

  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  texture_->set_color(1.0, 1.0, 1.0, 1.0);
  texture_->draw(4, coords, tex_coords);
  glMatrixMode(GL_TEXTURE);
  glPopMatrix();

  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glTranslatef(region.x1(), region.y1(), 0);
  glScalef(region.width(), region.height(), 1.0);


  vector<SLIVR::CM2Widget*> &widgets = colormap2_->widgets();
  for(unsigned int i=0; i < widgets.size(); i++) {
    widgets[i]->rasterize(*shader_factory_);
  }

  for(unsigned int i=0; i < widgets.size(); i++) {
    widgets[i]->draw();
  }

  glEnable(GL_BLEND);

  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();

  return CONTINUE_E;
}


} // namespace Skinner
} // namespace SCIRun
