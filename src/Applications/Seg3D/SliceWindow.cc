//  
//  For more information, please see: http://software.sci.utah.edu
//  
//  The MIT License
//  
//  Copyright (c) 2009 Scientific Computing and Imaging Institute,
//  University of Utah.
//  
//  
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
//    File   : SliceWindow.cc
//    Author : McKay Davis
//    Date   : Fri Oct 13 15:08:39 2006


#include <Applications/Seg3D/Painter.h>

#include <sci_comp_warn_fixes.h>

#include <sci_gl.h>

#include <algorithm>

#include <Core/Containers/Array3.h>
#include <Core/Datatypes/Field.h> 
#include <Core/Exceptions/GuiException.h>
#include <Core/Geom/GeomMaterial.h>
#include <Core/Geom/ColorMappedNrrdTextureObj.h>
#include <Core/Geom/GeomSwitch.h>
#include <Core/Skinner/GeomSkinnerVarSwitch.h>
#include <Core/Geom/GeomCull.h>
#include <Core/Geom/GeomGroup.h>
#include <Core/Geom/TexSquare.h>
#include <Core/Geom/FreeType.h>

#include <Core/Math/MiscMath.h>
#include <Core/Math/MinMax.h>
#include <Core/Thread/CleanupManager.h>
#include <Core/Thread/Runnable.h>
#include <Core/Thread/Mutex.h>
#include <Core/Util/Environment.h>
#include <Core/Volume/CM2Widget.h>
#include <Core/Geom/TextRenderer.h>
#include <Core/Geom/FontManager.h>
#include <Core/Geom/GeomColorMappedNrrdTextureObj.h>
#include <Core/Skinner/Variables.h>
#include <Core/Skinner/Signals.h>
#include <Core/Events/EventManager.h>
#include <Core/Events/SceneGraphEvent.h>
#include <Core/Util/FileUtils.h>

#include <stdlib.h>
#include <math.h>

#include <map>
#include <typeinfo>
#include <iostream>

#ifdef _WIN32
#  define snprintf _snprintf
#endif

namespace SCIRun {

SliceWindow::SliceWindow(Skinner::Variables *variables,
                                  Painter *painter) :  
  Skinner::Parent(variables), 
  painter_(painter),
  slices_(),
  recompute_slices_(false),
  center_(0,0,0),
  normal_(0,0,0),
  axis_(2),
  zoom_(variables, "SliceWindow::zoom", 100.0),
  flip_updown_(variables, "SliceWindow::flip_updown",0),
  flip_leftright_(variables, "SliceWindow::flip_leftright",0),
  slab_min_(0),
  slab_max_(0),
  show_guidelines_(1),
  pdown_(0),
  color_(variables, "SliceWindow::color", Skinner::Color(0,0,0,0)),
  show_grid_(variables, "SliceWindow::GridVisible",1),
  grid_font_size_(variables, "SliceWindow::GridFontSize",15),
  show_slices_(variables, "SliceWindow::SlicesVisible",1),
  show_slices_cache_(true),
  groupname_(variables, "SliceWindow::Group","default"),
  geom_switch_(0),
  geom_group_(0),
  selected_slice_(0),
  upper_range_slice_(-1),
  lower_range_slice_(-1)
{
  Skinner::Var<int> axis(variables, "axis", 2);
  set_axis(axis());
  REGISTER_CATCHER_TARGET(SliceWindow::redraw);
  REGISTER_CATCHER_TARGET(SliceWindow::do_PointerEvent);
  REGISTER_CATCHER_TARGET(SliceWindow::Autoview);
  REGISTER_CATCHER_TARGET(SliceWindow::zoom_in);
  REGISTER_CATCHER_TARGET(SliceWindow::zoom_out);
  REGISTER_CATCHER_TARGET(SliceWindow::flip_updown);
  REGISTER_CATCHER_TARGET(SliceWindow::flip_leftright);
}


void
SliceWindow::mark_redraw()
{
  throw_signal("SliceWindow::mark_redraw");
}


void
SliceWindow::render_guide_lines(Point mouse)
{
  if (!show_guidelines_) return;

  //  GLdouble yellow[4] = { 1.0, 0.76, 0.1, 0.8 };
  GLdouble white[4] = { 1.0, 1.0, 1.0, 1.0 };

  double vpw = get_region().width();
  double vph = get_region().height();

  // Push 2D gl view
  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();
  glScaled(2.0, 2.0, 2.0);
  glTranslated(-.5, -.5, -.5);
  glScaled(1.0/vpw, 1.0/vph, 1.0);
  CHECK_OPENGL_ERROR();

  glColor4dv(white);
  glBegin(GL_LINES); 
  glVertex3d(0, mouse.y(), mouse.z());
  glVertex3d(vpw, mouse.y(), mouse.z());
  glVertex3d(mouse.x(), 0, mouse.z());
  glVertex3d(mouse.x(), vph, mouse.z());
  glEnd();
  CHECK_OPENGL_ERROR();

  // Pop 2D gl view
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  CHECK_OPENGL_ERROR();
}



// Renders vertical and horizontal bars that represent
// selected slices in other dimensions.
void
SliceWindow::render_slice_lines(SliceWindows &windows)
{
  NrrdVolumeHandle vol = painter_->current_volume_;
  if (!vol.get_rep()) return;
  double upp = 100.0 / zoom_;    // World space units per one pixel

  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  vector<int> zero_idx(vol->nrrd_handle_->nrrd_->dim, 0);
  Point origin = vol->index_to_world(zero_idx);
  for (unsigned int i = 0; i < windows.size(); ++i) {
    SliceWindow *win = windows[i];
    if (!win->visible_()) continue;
    if (win == this) continue;
    if (win->groupname_() != this->groupname_()) continue;
    if (Dot(normal_, win->normal_) > 0.999) continue;

    // The span vector spans the volume edge to edge
    Vector span = Cross(normal_, win->normal_);
    vector<double> span_index = vol->vector_to_index(span);
    int span_axis = max_vector_magnitude_index(span_index);

    Vector dir1 = span;
    dir1[span_axis-1] = Abs(dir1[span_axis-1]);
    Vector dir2 = win->normal_;
    dir1.normalize();
    dir2.normalize();

    vector<int> pos_idx = vol->world_to_index(win->center_);

    // The lower left corner of the quad, assuming horizontal span
    vector<int> min_idx = pos_idx;
    min_idx[span_axis] = 0;
    const Point min = vol->index_to_world(min_idx);

    // The lower right corner of the quad, assuming horizontal span
    vector<int> max_idx = pos_idx;
    max_idx[span_axis] = vol->nrrd_handle_->nrrd_->axis[span_axis].size;
    const Point max = vol->index_to_world(max_idx);

    vector<int> one_idx = zero_idx;
    one_idx[win->axis_+1] = 1;
    double scale = (vol->index_to_world(one_idx) - origin).length();
    Vector wid = win->normal_;
    wid.normalize();

    wid *= Max(upp, scale);
    const Point min2 = min + wid;
    const Point max2 = max + wid;

    wid.normalize();
    wid *= upp;

    Skinner::Color color = win->color_;
    color.a = 1.0;
    glColor4dv(&(color.r));
    glBegin(GL_QUADS);    

    glVertex3d(min.x(), min.y(), min.z());
    glVertex3d(max.x(), max.y(), max.z());
    glVertex3dv(&(max-wid)(0));
    glVertex3dv(&(min-wid)(0));

    glVertex3d(min2.x(), min2.y(), min2.z());
    glVertex3d(max2.x(), max2.y(), max2.z());
    glVertex3dv(&(max2+wid)(0));
    glVertex3dv(&(min2+wid)(0));

    // Min Arrow Pointer
    const double ascale = upp * 8;
    const Point ula = min2 - ascale * dir1 + ascale * dir2;
    const Point lla = min  - ascale * dir1 - ascale * dir2;

    glVertex3d(min.x(), min.y(), min.z());
    glVertex3d(min2.x(), min2.y(), min2.z());
    glVertex3d(ula.x(), ula.y(), ula.z());
    glVertex3d(lla.x(), lla.y(), lla.z());

    // Max Arrow Pointer 
    const Point lra = max  + ascale * dir1 - ascale * dir2;
    const Point ura = max2 + ascale * dir1 + ascale * dir2;
    glVertex3d(max.x(), max.y(), max.z());
    glVertex3d(lra.x(), lra.y(), lra.z());
    glVertex3d(ura.x(), ura.y(), ura.z());
    glVertex3d(max2.x(), max2.y(), max2.z());

    glEnd();
  }
}



static double
div_d(double dividend, double divisor)
{
  return Floor(dividend/divisor);
}


// Renders vertical and horizontal bars that represent
// selected slices in other dimensions.
void
SliceWindow::render_grid()
{
  double units = zoom_ / 100.0;   // Pixels per world space unit
  const double pixels = 100.0;    // Number of target pixels for grid gap

  std::vector<double> gaps(1, 1.0);
  gaps.push_back(1/2.0);
  gaps.push_back(1/5.0);

  double realdiff = 10000;
  int selected = 0;
  for (unsigned int i = 0; i < gaps.size(); ++i) 
  {
    bool done = false;
    double diff = fabs(gaps[i]*units-pixels);
    while (!done) 
    {
      if (fabs((gaps[i]*10.0)*units - pixels) < diff) 
        gaps[i] *= 10.0;
      else if (fabs((gaps[i]/10.0)*units - pixels) < diff) 
        gaps[i] /= 10.0;
      else
        done = true;
      diff = fabs(gaps[i]*units-pixels);
    }

    if (diff < realdiff) 
    {
      realdiff = diff;
      selected = i;
    }
  }
  double gap = gaps[selected];

  const Skinner::RectRegion &region = get_region();
  const int vw = Ceil(region.width());
  const int vh = Ceil(region.height());

  double grid_color = 0.25;

  glDisable(GL_TEXTURE_2D);
  CHECK_OPENGL_ERROR();

  Point min = screen_to_world(Floor(region.x1()),Floor(region.y1()));
  Point max = screen_to_world(Ceil(region.x2())-1, Ceil(region.y2())-1);

  int xax = x_axis();
  int yax = y_axis();

  if( flip_updown_ ) 
  {
    double tmp = min(yax);
    min(yax) = max(yax);
    max(yax) = tmp;
  }
  
  if( flip_leftright_ ) 
  {
    double tmp = min(xax);
    min(xax) = max(xax);
    max(xax) = tmp;
  }

  min(xax) = div_d(min(xax), gap)*gap;
  min(yax) = div_d(min(yax), gap)*gap;

  std::vector<std::string> lab;
  lab.push_back("X: ");
  lab.push_back("Y: ");
  lab.push_back("Z: ");

  int num = 0;
  Point linemin = min;
  Point linemax = min;
  linemax(yax) = max(yax);
  std::string str;
  TextRenderer *renderer = FontManager::get_renderer(grid_font_size_);
  renderer->set_color(1,1,1,1);
  renderer->set_shadow_color(0,0,0,1);
  renderer->set_shadow_offset(1, -1);
  
  while (linemin(xax) < max(xax)) 
  {
    linemin(xax) = linemax(xax) = min(xax) + gap*num;
    glEnable(GL_BLEND);
    glBlendFunc(GL_ONE_MINUS_DST_COLOR, GL_ZERO);

    glColor4d(grid_color, grid_color, grid_color, 0.25);
    glBegin(GL_LINES);
    glVertex3dv(&linemin(0));
    glVertex3dv(&linemax(0));
    glEnd();

    str = to_string(linemin(xax));

    Point pos = world_to_screen(linemin);

    renderer->render(str, pos.x()+1, 2, 
                     TextRenderer::SHADOW | 
                     TextRenderer:: SW | 
                     TextRenderer::REVERSE);
    
    renderer->render(str, pos.x()+1, vh-2, 
                     TextRenderer::SHADOW | 
                     TextRenderer:: NW | 
                     TextRenderer::REVERSE);
    num++;
  }

  int vert_offset = 10;
/*
  if (lower_range_slice_ >= 0)
  {
    renderer->set_color(1.0,0.5,0.5,1);

    std::string destination_slice;
    if (lower_range_slice_ == upper_range_slice_)
      destination_slice = "Destination slice: " + to_string(lower_range_slice_);
    else
      destination_slice = "Destination range: " + to_string(lower_range_slice_) + ":" + to_string(upper_range_slice_);

    renderer->render(destination_slice, 10, vert_offset, 
                     TextRenderer::SHADOW | 
                     TextRenderer:: SW | 
                     TextRenderer::REVERSE);  
    vert_offset += 25;
  }

  if (selected_slice_ >= 0)
  {
    renderer->set_color(0.5,1,0.5,1);

    std::string source_slice = "Source slice: " + to_string(selected_slice_);
    renderer->render(source_slice, 10, vert_offset, 
                     TextRenderer::SHADOW | 
                     TextRenderer:: SW | 
                     TextRenderer::REVERSE);  
    vert_offset += 25;
  }
*/

  renderer->set_color(1,1,1,1);

  
  num = 0;
  linemin = linemax = min;
  linemax(xax) = max(xax);
  while (linemin(yax) < max(yax)) {
    linemin(yax) = linemax(yax) = min(yax) + gap*num;
    glEnable(GL_BLEND);
    glBlendFunc(GL_ONE_MINUS_DST_COLOR, GL_ZERO);

    glColor4d(grid_color, grid_color, grid_color, 0.25);
    glBegin(GL_LINES);
    glVertex3dv(&linemin(0));
    glVertex3dv(&linemax(0));
    glEnd();

    str = to_string(linemin(yax));
    Point pos = world_to_screen(linemin);
    renderer->render(str, 2, pos.y(), 
                     TextRenderer::SHADOW | 
                     TextRenderer::NW | 
                     TextRenderer::VERTICAL | 
                     TextRenderer::REVERSE);

    renderer->render(str, vw-2, pos.y(), 
                     TextRenderer::SHADOW | 
                     TextRenderer::NE | 
                     TextRenderer::VERTICAL | 
                     TextRenderer::REVERSE);
    num++;
  }

  CHECK_OPENGL_ERROR();
}


// Returns the index to the axis coordinate that is most parallel and 
// in the direction of X in the screen.  
// 0 for x, 1 for y, and 2 for z
int
SliceWindow::x_axis()
{
  Vector adir = Abs(x_dir());
  if ((adir[0] > adir[1]) && (adir[0] > adir[2])) return 0;
  if ((adir[1] > adir[0]) && (adir[1] > adir[2])) return 1;
  return 2;
}


// Returns the index to the axis coordinate that is most parallel and 
// in the direction of Y in the screen.  
// 0 for x, 1 for y, and 2 for z
int
SliceWindow::y_axis()
{
  Vector adir = Abs(y_dir());
  if ((adir[0] > adir[1]) && (adir[0] > adir[2])) return 0;
  if ((adir[1] > adir[0]) && (adir[1] > adir[2])) return 1;
  return 2;
}


void
SliceWindow::setup_gl_view()
{
  NrrdVolumeHandle &vol = painter_->current_volume_;
  if (!vol.get_rep()) return;

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  if (axis_ == 0) {        // Sagittal  screen +X -> +Y, screen +Y -> +Z
    glRotated(-90,0.,1.,0.);
    glRotated(-90,1.,0.,0.);

    if( flip_updown_ )
      glRotated(180,0.,1.,0.);
    
    if( flip_leftright_ )
      glRotated(180,0.,0.,1.);

  } else if (axis_ == 1) { // Coronal   screen +X -> +X, screen +Y -> +Z

    glRotated(-90,1.,0.,0.);

    if( flip_updown_ )
      glRotated(180,1.,0.,0.);

    if( flip_leftright_ )
      glRotated(180,0.,0.,1.);

  } else if (axis_ == 2) { // Axial     screen +X -> +X, screen +Y -> +Y

    if( flip_updown_ ) {
      glRotated(180,1.,0.,0.);
    }

    if( flip_leftright_ ) {
      glRotated(180,0.,1.,0.);
    }
  }

  CHECK_OPENGL_ERROR();
  
  // Do this here because x_axis and y_axis functions use these matrices
  glGetIntegerv(GL_VIEWPORT, gl_viewport_);
  glGetDoublev(GL_MODELVIEW_MATRIX, gl_modelview_matrix_);
  glGetDoublev(GL_PROJECTION_MATRIX, gl_projection_matrix_);
  CHECK_OPENGL_ERROR();

  double hwid = get_region().width()  * 50.0 / zoom_;
  double hhei = get_region().height() * 50.0 / zoom_;

  double cx = center_(x_axis());
  double cy = center_(y_axis());

  double diagonal = hwid*hwid+hhei*hhei;

  double maxz = center_(axis_) + diagonal * zoom_ / 100.0;
  double minz = center_(axis_) - diagonal * zoom_ / 100.0;

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();

  if( flip_updown_ ) {
    cy = -cy;
  }
  
  if( flip_leftright_ ) {
    cx = -cx;
  }

  glOrtho(cx - hwid, cx + hwid, 
	  cy - hhei, cy + hhei, 
	  minz, maxz);

  glGetIntegerv(GL_VIEWPORT, gl_viewport_);
  glGetDoublev(GL_MODELVIEW_MATRIX, gl_modelview_matrix_);
  glGetDoublev(GL_PROJECTION_MATRIX, gl_projection_matrix_);
  CHECK_OPENGL_ERROR();
}


void
SliceWindow::render_slices()
{
  for (unsigned int s = 0; s < slices_.size(); ++s) {
    if (slices_[s]->volume_->dirty_) {
      recompute_slices_ = true;
      break;
    }
  }

  if (show_slices_ != show_slices_cache_) {
    recompute_slices_ = true;
    painter_->redraw_all();
  }
  show_slices_cache_ = show_slices_;

  if (recompute_slices_) {
    recompute_slices_ = false;
    NrrdVolumes volumes;
    painter_->build_volume_list(volumes);
    
    slices_.clear();
    double offset = 0.0;
        
    for (unsigned int i = 0; i < volumes.size(); ++i)
    {
      NrrdVolumeHandle volume = volumes[i];
      volume->purge_unused_slices();
      VolumeSliceHandle slice = 
        volume->get_volume_slice(SLIVR::Plane(SLIVR::Point(center_.x(),
        center_.y(),center_.z()), SLIVR::Vector(normal_.x(),
        normal_.y(),normal_.z() )));
      slices_.push_back(slice);
      int objid = axis_*(volume->label_+1);
      GeomIndexedGroup *ggroup = volume->get_geom_group();
      GeomHandle oldgeom = ggroup->getObj(objid);
      if (oldgeom.get_rep()) {
        ggroup->delObj(objid);
      }
      
      GeomHandle newgeom = slice->get_geom_texture();
      if (newgeom.get_rep()) {
        GeomColorMappedNrrdTextureObj *gcmnto = 
          dynamic_cast<GeomColorMappedNrrdTextureObj *>(newgeom.get_rep());

        gcmnto->set_offset(offset);
        offset += 1.0;
        Skinner::Var<bool> 
          all_slices(painter_->get_vars(),"Painter::slices_visible", true);
        GeomSkinnerVarSwitch *windowswitch = new
          GeomSkinnerVarSwitch(gcmnto, show_slices_);
        GeomSkinnerVarSwitch *allswitch = new
          GeomSkinnerVarSwitch(windowswitch, all_slices);


        ggroup->addObj(allswitch, objid);
      }
    }
  }

  for (unsigned int s = 0; s < slices_.size(); ++s) {
    slices_[s]->draw();
  }
}


void
SliceWindow::extract_slices()
{
  recompute_slices_ = true;
}


void
SliceWindow::set_axis(unsigned int axis)
{
  axis_ = axis % 3;
  normal_ = Vector(axis == 0 ? 1 : 0,
                   axis == 1 ? 1 : 0,
                   axis == 2 ? 1 : 0);
  recompute_slices_ = true;
}


void
SliceWindow::set_probe()
{
  if (painter_->cur_window_ == this) return;
  NrrdVolumeHandle &vol = painter_->current_volume_;
  if (!vol.get_rep()) return;
  Point newcenter = center_;
  newcenter(axis_) = painter_->pointer_pos_(axis_);
  vector<double> nindex = vol->point_to_index(newcenter);
  if (nindex[axis_+1] >= 0 && nindex[axis_+1] < vol->max_index(axis_+1)) {
    center_ = newcenter;
    recompute_slices_ = true;
  }
}


void
SliceWindow::move_slice(int amount)
{
  NrrdVolumeHandle &vol = painter_->current_volume_;
  if (!vol.get_rep()) return;
  
  // transfrom the original center
  std::vector<double> nindex = vol->point_to_index(center_);
  nindex[axis_+1] = Floor(nindex[axis_+1]) + static_cast<double>(amount) + 0.5;

  Point newcenter = vol->index_to_point(nindex);
    
  if (nindex[axis_+1] >= 0 && nindex[axis_+1] < vol->max_index(axis_+1)) 
  {
    Vector delta = newcenter-center_;
    center_ = newcenter;
    recompute_slices_ = true;
    painter_->set_pointer_position(painter_->pointer_pos_ +delta);
    painter_->redraw_all();
  }
}


BaseTool::propagation_state_e
SliceWindow::zoom_in(event_handle_t &)
{
  zoom_ = zoom_*1.1;
  mark_redraw();
  return CONTINUE_E;
}


BaseTool::propagation_state_e
SliceWindow::zoom_out(event_handle_t &)
{
  zoom_ = zoom_/1.1;
  mark_redraw();
  return CONTINUE_E;
}

BaseTool::propagation_state_e
SliceWindow::flip_updown(event_handle_t &)
{
  flip_updown_ = !flip_updown_;
  mark_redraw();
  return CONTINUE_E;
}


BaseTool::propagation_state_e
SliceWindow::flip_leftright(event_handle_t &)
{
  flip_leftright_ = !flip_leftright_;
  mark_redraw();
  return CONTINUE_E;
}

  
Point
SliceWindow::screen_to_world(unsigned int x, unsigned int y)
{
  GLdouble xyz[3];
  gluUnProject(double(x)+0.5, double(y)+0.5, 0,
	       gl_modelview_matrix_, 
	       gl_projection_matrix_,
	       gl_viewport_,
	       xyz+0, xyz+1, xyz+2);
  xyz[axis_] = center_(axis_);
  return Point(xyz[0], xyz[1], xyz[2]);
}


Point
SliceWindow::world_to_screen(const Point &world)
{
  GLdouble xyz[3];
  gluProject(world(0), world(1), world(2),
             gl_modelview_matrix_, 
             gl_projection_matrix_,
	     gl_viewport_,
             xyz+0, xyz+1, xyz+2);
  xyz[0] -= gl_viewport_[0];
  xyz[1] -= gl_viewport_[1];

  return Point(xyz[0], xyz[1], xyz[2]);
}


Vector
SliceWindow::x_dir()
{
  return screen_to_world(1,0) - screen_to_world(0,0);
}


Vector
SliceWindow::y_dir()
{
  return screen_to_world(0,1) - screen_to_world(0,0);
}



BaseTool::propagation_state_e
SliceWindow::Autoview(event_handle_t &)
{
  if (painter_->current_volume_.get_rep())
  {
    autoview(painter_->current_volume_);
  }
  return CONTINUE_E;
}


void
SliceWindow::autoview(NrrdVolumeHandle &volume, double offset)
{
  double wid = get_region().width() -  2*offset;
  double hei = get_region().height() - 2*offset;

  int xax = x_axis();
  int yax = y_axis();

  if (volume.get_rep())
  {
    vector<int> zero(volume->nrrd_handle_->nrrd_->dim, 0);
    vector<int> index = zero;
    index[xax+1] = volume->nrrd_handle_->nrrd_->axis[xax+1].size;
    double w_wid = (volume->index_to_world(index) - 
                    volume->index_to_world(zero)).length();
    double w_ratio = wid/w_wid;
    
    index = zero;
    index[yax+1] = volume->nrrd_handle_->nrrd_->axis[yax+1].size;
    double w_hei = (volume->index_to_world(index) - 
                    volume->index_to_world(zero)).length();
    double h_ratio = hei/w_hei;
    
    zoom_ = Min(w_ratio*100.0, h_ratio*100.0);
    if (zoom_ < 1.0) zoom_ = 100.0; // ??
    center_(xax) = volume->center()(xax);
    center_(yax) = volume->center()(yax);
  } else {
    center_ = Point(0,0,0);
    zoom_ = 100;
  }
  mark_redraw();
}


BaseTool::propagation_state_e
SliceWindow::do_PointerEvent(event_handle_t &event)
{
  ASSERT(dynamic_cast<Skinner::PointerSignal *>(event.get_rep()));
  Skinner::PointerSignal *ps = (Skinner::PointerSignal *)(event.get_rep());
  PointerEvent *pointer = ps->get_pointer_event();
  ASSERT(pointer);

  if (pointer->get_pointer_state() & PointerEvent::MOTION_E) {
    throw_signal("SliceWindow::pointer_motion");
  }

  if (get_region().inside(pointer->get_x(), pointer->get_y()))
  {
    painter_->cur_window_ = this;
    painter_->set_pointer_position(screen_to_world(pointer->get_x(), 
                                                   pointer->get_y()));
  } 
  else if (pointer->get_pointer_state() & PointerEvent::BUTTON_PRESS_E)
  {
    ps->set_signal_result(STOP_E);
    return STOP_E;
  }

  return CONTINUE_E;
}


BaseTool::propagation_state_e
SliceWindow::process_event(event_handle_t &event)
{
  BaseTool::propagation_state_e state = Parent::process_event(event);

  if (state == CONTINUE_E && painter_->cur_window_ == this) {
    painter_->tm_.propagate_event(event);
  }

  return state;
}
   

GeomIndexedGroup *
SliceWindow::get_geom_group()
{
  if (!geom_group_) {
    geom_group_ = new GeomIndexedGroup();
    geom_switch_ = new GeomSkinnerVarSwitch(geom_group_, show_slices_);
    event_handle_t add_geom_switch_event = 
      new SceneGraphEvent(geom_switch_, get_id()+" Transparent");
    EventManager::add_event(add_geom_switch_event);
  }
  return geom_group_;
}


string
double_to_string(double val)
{
  char s[50];
  snprintf(s, 49, "%1.2f", val);
  return string(s);
}


int
SliceWindow::get_signal_id(const string &signalname) const
{
  if (signalname == "SliceWindow::mark_redraw") return 1;
  if (signalname == "SliceWindow::pointer_motion") return 1;
  return 0;
}


BaseTool::propagation_state_e
SliceWindow::redraw(event_handle_t &)
{
  const Skinner::RectRegion &region = get_region();
  if (region.width() <= 0 || region.height() <= 0) return CONTINUE_E;
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  GLint viewport[4];
  glGetIntegerv(GL_VIEWPORT, viewport);    
  glViewport(Floor(region.x1()), Floor(region.y1()), 
             Ceil(region.width()), Ceil(region.height()));

  setup_gl_view();
  CHECK_OPENGL_ERROR();

  // Render the individual slices
  render_slices();

  if (show_grid_())
  {
    render_grid();
  }

  render_slice_lines(painter_->windows_);

  event_handle_t redraw_window = new RedrawSliceWindowEvent(*this);
  painter_->tm_.propagate_event(redraw_window);
  CHECK_OPENGL_ERROR();

  glViewport(viewport[0],viewport[1],viewport[2],viewport[3]);
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  CHECK_OPENGL_ERROR();

  return CONTINUE_E;
}


BaseTool::propagation_state_e
SliceWindow::copy_current_slice(int dir)
{
  if (!painter_->check_for_active_label_volume("Punch current slice hotkey"))
  {
    return STOP_E;
  }

  painter_->volume_lock_.lock();

  // Get the current volume.
  NrrdVolumeHandle &vol = painter_->current_volume_;

  vol->lock.lock();

  if (dir == 0)
  {
    // Undo for punch.
    NrrdDataHandle volnrrd = vol->nrrd_handle_->clone();
    UndoHandle undo =
      new UndoReplaceVolume(painter_, "Undo Punch", vol, volnrrd);
    painter_->push_undo(undo);
  }

  // Get the original slice.
  SLIVR::Plane plane(SLIVR::Point(center_.x(), center_.y(), center_.z()), 
    SLIVR::Vector(normal_.x(),normal_.y(), normal_.z()));
  VolumeSliceHandle slice = vol->get_volume_slice(plane);
  unsigned int clabel = vol->label_;
  const int axis = slice->get_axis();

  size_t start = 0;
  size_t end = vol->nrrd_handle_->nrrd_->axis[axis].size;
  if (dir < 0)
  {
    vector<int> sindex = vol->world_to_index(plane.project(SLIVR::Point(0,0,0)));
    int index = sindex[axis];
    if (index <= 0)
    {
      slice = 0;
      vol->lock.unlock();
      painter_->volume_lock_.unlock();
      painter_->set_status("Already at bottom slice.");
      return QUIT_AND_STOP_E;
    }
    start = index-1;
    end = start+1;
  }
  else if (dir > 0)
  {
    vector<int> sindex = vol->world_to_index(plane.project(SLIVR::Point(0,0,0)));
    int index = sindex[axis];
    if (index >= (int)vol->nrrd_handle_->nrrd_->axis[axis].size-1)
    {
      slice = 0;
      vol->lock.unlock();
      painter_->volume_lock_.unlock();
      painter_->set_status("Already at top slice.");
      return QUIT_AND_STOP_E;
    }
    start = index+1;
    end = start+1;
  }

  NrrdDataHandle curslice = new NrrdData();

  for (size_t i = start; i < end; i++)
  {
    // Pull out the current slice.
    nrrdSlice(curslice->nrrd_, vol->nrrd_handle_->nrrd_, axis, i);

    if (dir != 0)
    {
      NrrdDataHandle undoslice = curslice->clone();
      UndoHandle undo =
        new UndoReplaceSlice(painter_, "Undo Slice Copy",
                             vol, undoslice, axis, i);
      painter_->push_undo(undo);
    }

    // Copy the bitplane from the visible slice.
    VolumeOps::bit_copy(curslice, clabel, slice->nrrd_handle_, clabel);

    // Clear the content for nrrdSplice
    if (vol->nrrd_handle_->nrrd_->content)
    {
      vol->nrrd_handle_->nrrd_->content[0] = 0;
    }

    // Put the results back.
    if (nrrdSplice(vol->nrrd_handle_->nrrd_,
                   vol->nrrd_handle_->nrrd_,
                   curslice->nrrd_,
                   axis, i))
    {
      vol->lock.unlock();
      painter_->volume_lock_.unlock();
      char *err = biffGetDone(NRRD);
      
      cerr << string("Error on line #") 
           << to_string(__LINE__)
           << string(" executing nrrd command: nrrdSplice \n")
           << string("Message: ") 
           << err
           << std::endl;

      free(err);
      return QUIT_AND_STOP_E;
    }
  }

  // Clear the slice pointer.
  slice = 0;

  vol->set_dirty();
  vol->lock.unlock();
  painter_->volume_lock_.unlock();

  painter_->extract_all_window_slices();
  painter_->redraw_all();

  return STOP_E;
}


BaseTool::propagation_state_e
SliceWindow::punch_current_slice(event_handle_t &)
{
  if (!painter_->check_for_active_label_volume("Punch current slice hotkey"))
  {
    return STOP_E;
  }
  
  painter_->set_status("Punching the current slice through the volume.");

  return copy_current_slice(0);
}


BaseTool::propagation_state_e
SliceWindow::copy_current_slice_up(event_handle_t &)
{
  if (!painter_->check_for_active_label_volume("Copy current slice up hotkey"))
  {
    return STOP_E;
  }
  
  painter_->set_status("Copying the current slice upward.");

  return copy_current_slice(1);
}


BaseTool::propagation_state_e
SliceWindow::copy_current_slice_down(event_handle_t &)
{
  if (!painter_->check_for_active_label_volume("Copy current slice down hotkey"))
  {
    return STOP_E;
  }
  
  painter_->set_status("Copying the current slice downward.");

  return copy_current_slice(-1);
}


BaseTool::propagation_state_e
SliceWindow::erase_current_slice(event_handle_t &)
{
  if (!painter_->check_for_active_label_volume("Erase slice hotkey"))
  {
    return STOP_E;
  }
  
  painter_->set_status("Erasing the current slice.");

  painter_->volume_lock_.lock();

  // Get the current volume.
  NrrdVolumeHandle &vol = painter_->current_volume_;

  vol->lock.lock();

  // Get the original slice.
  SLIVR::Plane plane(SLIVR::Point(center_.x(),center_.y(),center_.z()), 
                    SLIVR::Vector(normal_.x(),normal_.y(),normal_.z()));
  VolumeSliceHandle slice = vol->get_volume_slice(plane);
  unsigned int clabel = vol->label_;
  const int axis = slice->get_axis();
  const vector<int> window_center = vol->world_to_index(center_);
  const int coord = window_center[axis];

  NrrdDataHandle curslice = new NrrdData();

  NrrdDataHandle undoslice = slice->nrrd_handle_->clone();
  UndoHandle undo =
    new UndoReplaceSlice(painter_, "Undo Erase Slice",
                         vol, undoslice, axis, coord);
  painter_->push_undo(undo);

  // Copy the bitplane from the visible slice.
  VolumeOps::bit_clear(slice->nrrd_handle_, clabel);

  // Clear the content for nrrdSplice
  if (vol->nrrd_handle_->nrrd_->content)
  {
    vol->nrrd_handle_->nrrd_->content[0] = 0;
  }

  // Put the results back.
  if (nrrdSplice(vol->nrrd_handle_->nrrd_,
                 vol->nrrd_handle_->nrrd_,
                 slice->nrrd_handle_->nrrd_,
                 axis, coord))
  {
    vol->lock.unlock();
    painter_->volume_lock_.unlock();
    char *err = biffGetDone(NRRD);
      
    cerr << string("Error on line #") 
         << to_string(__LINE__)
         << string(" executing nrrd command: nrrdSplice \n")
         << string("Message: ") 
         << err
         << std::endl;

    free(err);
    return QUIT_AND_STOP_E;
  }

  // Clear the slice pointer.
  slice = 0;

  vol->set_dirty();
  vol->lock.unlock();
  painter_->volume_lock_.unlock();

  painter_->extract_all_window_slices();
  painter_->redraw_all();

  return STOP_E;
}


}
