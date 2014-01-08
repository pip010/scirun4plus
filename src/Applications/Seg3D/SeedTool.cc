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
//    File   : SeedTool.cc
//    Author : McKay Davis
//    Date   : Tue Sep 26 18:44:34 2006

#include <Applications/Seg3D/Painter.h>

#include <Applications/Seg3D/SeedTool.h>
#include <Core/Events/EventManager.h>
#include <Core/Events/SceneGraphEvent.h>
#include <Core/Geom/GeomCylinder.h>
#include <Core/Geom/GeomMaterial.h>
#include <sci_gl.h>


namespace SCIRun {

SeedTool::SeedTool(const string &name, Painter *painter) :
  BaseTool(name),
  PointerTool(name),
  painter_(painter),
  seed_lock_(name.c_str()),
  just_one_seed_mode_(false),
  button_(1),
  last_seed_window_(NULL),
  last_seed_index_(-1)
{
}


SeedTool::~SeedTool()
{
  // Clear the 3D seeds.
  event_handle_t scene_event = new SceneGraphEvent(0, "Seeds");
  EventManager::add_event(scene_event);
  painter_->redraw_all();
}


BaseTool::propagation_state_e
SeedTool::pointer_down(int b, int x, int y, unsigned int m, int t)
{
  if (m & EventModifiers::SHIFT_E) { return CONTINUE_E; }

  last_seed_index_ = -1;

  if (b == button_)
  {
    // Add a new seed point.
    seed_lock_.lock();
    last_seed_window_ = painter_->cur_window_;
    if (just_one_seed_mode_)
    {
      seeds_.clear();
    }
    seeds_.push_back(painter_->pointer_pos_);
    last_seed_index_ = seeds_.size() - 1;
    seed_change_callback();

    seed_lock_.unlock();
    painter_->redraw_all();
    return STOP_E;
  }
  else if (!just_one_seed_mode_ && b == 3 && seeds_.size())
  {
    // Change an existing seed point.
    seed_lock_.lock();
    last_seed_window_ = painter_->cur_window_;
    const double max_dist2_to_check = 25.0;
    double mind2 = max_dist2_to_check + 1.0;
    const Point sppos =
      last_seed_window_->world_to_screen(painter_->pointer_pos_);
    for (size_t i = 0; i < seeds_.size(); i++)
    {
      const double d2 =
        (last_seed_window_->world_to_screen(seeds_[i]) - sppos).length2();
      if (d2 < mind2 && d2 < max_dist2_to_check)
      {
        mind2 = d2;
        last_seed_index_ = i;
      }
    }

    if (last_seed_index_ != -1)
    {
      seeds_[last_seed_index_] = painter_->pointer_pos_;
      seed_change_callback();
    }

    seed_lock_.unlock();
    painter_->redraw_all();
    return STOP_E;
  }
  else if (!just_one_seed_mode_ && b == 2 && seeds_.size())
  {
    // Erase an existing seed point.
    seed_lock_.lock();
    last_seed_window_ = painter_->cur_window_;
    const double max_dist2_to_check = 25.0;
    double mind2 = max_dist2_to_check + 1.0;
    const Point sppos =
      last_seed_window_->world_to_screen(painter_->pointer_pos_);
    for (size_t i = 0; i < seeds_.size(); i++)
    {
      const double d2 =
        (last_seed_window_->world_to_screen(seeds_[i]) - sppos).length2();
      if (d2 < mind2 && d2 < max_dist2_to_check)
      {
        mind2 = d2;
        last_seed_index_ = i;
      }
    }

    if (last_seed_index_ != -1)
    {
      seeds_.erase(seeds_.begin() + last_seed_index_);
      last_seed_index_ = -1;
      seed_change_callback();
    }

    seed_lock_.unlock();
    painter_->redraw_all();
    return STOP_E;
  }

  return CONTINUE_E;;
}


BaseTool::propagation_state_e
SeedTool::pointer_motion(int b, int x, int y, unsigned int m, int t)
{
  if (last_seed_index_ != -1)
  {
    seed_lock_.lock();
    seeds_[last_seed_index_] = painter_->pointer_pos_;
    seed_change_callback();
    seed_lock_.unlock();
    painter_->redraw_all();
    return STOP_E;
  }

  return CONTINUE_E;
}


BaseTool::propagation_state_e
SeedTool::pointer_up(int b, int x, int y, unsigned int m, int t)
{
  last_seed_index_ = -1;
  return CONTINUE_E;
}


BaseTool::propagation_state_e 
SeedTool::process_event(event_handle_t event)
{
  RedrawSliceWindowEvent *redraw = 
    dynamic_cast<RedrawSliceWindowEvent *>(event.get_rep());
  if (redraw) {
    draw_gl(redraw->get_window());
  }

  if (dynamic_cast<SetLayerEvent *>(event.get_rep())) {
    seed_lock_.lock();
    seeds_.clear();
    seed_change_callback();
    seed_lock_.unlock();
    painter_->redraw_all();
  }

  if (dynamic_cast<FinishEvent *>(event.get_rep())) {
    run_filter();
  }

  return CONTINUE_E;
}
  

void
SeedTool::draw_gl(SliceWindow &window)
{
  seed_lock_.lock();
  for (size_t i = 0; i < seeds_.size(); ++i) {
    draw_seed(seeds_[i], window);
  }
  seed_lock_.unlock();
}


void
SeedTool::draw_seed(const Point &seed, SliceWindow &window)
{
  Vector left = window.x_dir();
  Vector up = window.y_dir();
  Point center = seed;
  Point p;

  NrrdVolumeHandle curvol = painter_->current_volume_;
  bool in = false;
  if (curvol.get_rep())
  {
    vector<int> pidx = curvol->world_to_index(center);
    vector<int> widx = curvol->world_to_index(window.center_);
    in = (pidx[window.axis_+1] == widx[window.axis_+1]);
  }

  const double scene_scale = painter_->scene_scale();

  const double scale = in ? 7.0 : 4.0;

  // Pixels per world space unit times scene scale
  const double units = scene_scale * window.zoom_ / 100.0;
  const double units_over_2 = units/2.0;
  const double e =
    units_over_2 + Clamp( units_over_2, scale, Max(units, scale));

  for (int pass = 0; pass < 3; ++pass) {
    glLineWidth(Max(scale - pass*2.0, 1.0));
    if (pass == 0)
      glColor4d(0.0, 0.0, 0.0, 1.0);
    else if (pass == 1)
      glColor4d(in ? 0.0 : 0.0, in ? 1.0 : 0.6, 0.0, 1.0);
    else
      glColor4d(in ? 0.5 : 0.5, in ? 1.0 : 0.9, 0.5, 1.0);

    glBegin(GL_LINES);    
    p = center + units_over_2 * up;
    glVertex3dv(&p(0));
    p = center + e * up;
    glVertex3dv(&p(0));
    
    p = center - units_over_2 * up;
    glVertex3dv(&p(0));
    p = center - e * up;
    glVertex3dv(&p(0));
    
    p = center + units_over_2 * left;
    glVertex3dv(&p(0));
    p = center + e * left;
    glVertex3dv(&p(0));
    
    p = center - units_over_2 * left;
    glVertex3dv(&p(0));
    p = center - e * left;
    glVertex3dv(&p(0));
    glEnd();
    CHECK_OPENGL_ERROR();
  }

  glLineWidth(1.0);
}


void
SeedTool::convert_seeds_to_indices(vector<vector<int> > &iseeds,
                                   NrrdVolumeHandle &vol)
{
  iseeds.clear();
  for (size_t i = 0; i < seeds_.size(); i++)
  {
    vector<int> newseed = vol->world_to_index(seeds_[i]);
    if (vol->index_valid(newseed))
    {
      iseeds.push_back(newseed);
    }
  }
}

void
SeedTool::seed_change_callback()
{
  recreate_seed_geom();
}


void
SeedTool::recreate_seed_geom()
{
  const double sscale = painter_->scene_scale();
  const Vector xoff(sscale * 5.0, 0.0, 0.0);
  const Vector yoff(0.0, sscale * 5.0, 0.0);
  const Vector zoff(0.0, 0.0, sscale * 5.0);

  GeomCylinders *lines = new GeomCylinders();
  lines->set_radius(2.0 * sscale);
  for (size_t i = 0; i < seeds_.size(); i++)
  {
    lines->add(seeds_[i]-xoff, seeds_[i]+xoff);
    lines->add(seeds_[i]-yoff, seeds_[i]+yoff);
    lines->add(seeds_[i]-zoff, seeds_[i]+zoff);
  }             

  MaterialHandle green = new Material(Color(0.0, 1.0, 0.0));
  GeomHandle seedgeom = new GeomMaterial(lines, green);

  event_handle_t scene_event =
    new SceneGraphEvent(seedgeom, "Seeds");
  EventManager::add_event(scene_event);
}

}

