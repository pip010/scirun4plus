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
//    File   : PolylineTool.cc
//    Author : McKay Davis
//    Date   : Tue Sep 26 18:44:34 2006

#include <Applications/Seg3D/Painter.h>

#include <Applications/Seg3D/PolylineTool.h>
#include <Core/Geometry/CompGeom.h>
#include <Core/Events/keysyms.h>
#include <sci_gl.h>


namespace SCIRun {

PolylineTool::
PolylineTool(const string &name, Painter *painter) :
  BaseTool(name),
  PointerTool(name),
  painter_(painter),
  seeds_(),
  seed_lock_("PolylineTool"),
  last_seed_window_(NULL),
  last_seed_index_(-1)
{
}


BaseTool::propagation_state_e
PolylineTool::pointer_down(int b, int x, int y, unsigned int m, int t)
{
  if (m & EventModifiers::SHIFT_E) { return CONTINUE_E; }

  last_seed_index_ = -1;

  if (m & EventModifiers::CONTROL_E)
  {
    // Drag all the points.
    seed_lock_.lock();
    drag_seeds_ = seeds_;
    seed_lock_.unlock();
    drag_point_ = painter_->pointer_pos_;
    return STOP_E;
  }
  else if (b == 1)
  {
    // Add a new seed point.
    seed_lock_.lock();
    last_seed_window_ = painter_->cur_window_;
    seeds_.push_back(painter_->pointer_pos_);
    last_seed_index_ = seeds_.size() - 1;
    seed_change_callback(last_seed_index_);

    seed_lock_.unlock();
    painter_->redraw_all();
    return STOP_E;
  }
  else if (b == 3 && seeds_.size())
  {
    // Change an existing seed point.
    seed_lock_.lock();
    last_seed_window_ = painter_->cur_window_;
    const double max_dist2_to_check = 25.0;
    double mind2 = max_dist2_to_check + 1.0;
    const int axis = last_seed_window_->axis_;
    Point curpos = painter_->pointer_pos_;
    curpos(axis) = 0.0;
    const Point sppos = last_seed_window_->world_to_screen(curpos);
    for (size_t i = 0; i < seeds_.size(); i++)
    {
      Point seed = seeds_[i];
      seed(axis) = 0.0;
      const double d2 =
        (last_seed_window_->world_to_screen(seed) - sppos).length2();
      if (d2 < mind2 && d2 < max_dist2_to_check)
      {
        mind2 = d2;
        last_seed_index_ = i;
      }
    }

    if (last_seed_index_ != -1)
    {
      const double axval = seeds_[last_seed_index_](axis);
      seeds_[last_seed_index_] = painter_->pointer_pos_;
      seeds_[last_seed_index_](axis) = axval;
      seed_change_callback(last_seed_index_);
    }
    else if (seeds_.size() > 1)
    {
      last_seed_index_ = compute_nearest_segment_index(curpos, axis);

      if (last_seed_index_ != -1)
      {
        seeds_.insert(seeds_.begin() + last_seed_index_,
                      painter_->pointer_pos_);
        seed_change_callback(-2);
      }
    }

    seed_lock_.unlock();
    painter_->redraw_all();
    return STOP_E;
  }
  else if (b == 2 && seeds_.size())
  {
    // Erase an existing seed point.
    seed_lock_.lock();
    last_seed_window_ = painter_->cur_window_;
    const double max_dist2_to_check = 25.0;
    double mind2 = max_dist2_to_check + 1.0;
    const int axis = last_seed_window_->axis_;
    Point curpos = painter_->pointer_pos_;
    curpos(axis) = 0.0;
    const Point sppos = last_seed_window_->world_to_screen(curpos);
    for (size_t i = 0; i < seeds_.size(); i++)
    {
      Point seed = seeds_[i];
      seed(axis) = 0.0;
      const double d2 =
        (last_seed_window_->world_to_screen(seed) - sppos).length2();
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
      seed_change_callback(-1);
    }

    seed_lock_.unlock();
    painter_->redraw_all();
    return STOP_E;
  }

  return CONTINUE_E;;
}


int
PolylineTool::compute_nearest_segment_index(const Point &curpos, int axis)
{
  // No seed point near.  Check for insertion of new point into an
  // existing line segment.
  int index = -1;
  const double max_dist2_to_check = 25.0;
  double mind2 = max_dist2_to_check + 1.0;
  const Point sppos = last_seed_window_->world_to_screen(curpos);

  for (size_t i = 0; i < seeds_.size(); i++)
  {
    const size_t j = (i + 1) % seeds_.size();
    Point seedi = seeds_[i];
    Point seedj = seeds_[j];
    seedi(axis) = 0.0;
    seedj(axis) = 0.0;
    const Point sseedi = last_seed_window_->world_to_screen(seedi);
    const Point sseedj = last_seed_window_->world_to_screen(seedj);
    const double d2 = distance_to_line2(sppos, sseedi, sseedj);

    if (d2 < mind2 && d2 < max_dist2_to_check)
    {
      mind2 = d2;
      index = i+1;
    }
  }
  return index;
}


BaseTool::propagation_state_e
PolylineTool::pointer_motion(int b, int x, int y, unsigned int m, int t)
{
  if (last_seed_index_ != -1)
  {
    seed_lock_.lock();
    const int axis = last_seed_window_->axis_;
    const double axval = seeds_[last_seed_index_](axis);
    seeds_[last_seed_index_] = painter_->pointer_pos_;
    seeds_[last_seed_index_](axis) = axval;
    seed_change_callback(last_seed_index_);
    seed_lock_.unlock();
    painter_->redraw_all();
    return STOP_E;
  }
  else if (drag_seeds_.size())
  {
    seed_lock_.lock();
    Vector offset = painter_->pointer_pos_ - drag_point_;
    for (size_t i = 0; i < seeds_.size(); i++)
    {
      seeds_[i] = drag_seeds_[i] + offset;
      seed_change_callback(-3);
    }
    seed_lock_.unlock();
    painter_->redraw_all();
    return STOP_E;
  }

  return CONTINUE_E;
}


BaseTool::propagation_state_e
PolylineTool::pointer_up(int b, int x, int y, unsigned int m, int t)
{
  last_seed_index_ = -1;
  drag_seeds_.clear();
  return CONTINUE_E;
}


BaseTool::propagation_state_e 
PolylineTool::process_event(event_handle_t event)
{
  RedrawSliceWindowEvent *redraw = 
    dynamic_cast<RedrawSliceWindowEvent *>(event.get_rep());
  if (redraw) {
    draw_gl(redraw->get_window());
  }

  if (dynamic_cast<SetLayerEvent *>(event.get_rep())) {
    seed_lock_.lock();
    seeds_.clear();
    seed_change_callback(-4);
    seed_lock_.unlock();
    painter_->redraw_all();
  }

  if (dynamic_cast<FinishEvent *>(event.get_rep())) {
    const bool erase =
      painter_->get_vars()->get_bool("Painter::polylinetool_erase");
    return run_filter(erase);
  }

  KeyEvent *keyevent = dynamic_cast<KeyEvent *>(event.get_rep());
  if (keyevent && keyevent->get_keyval() == SCIRun_f)
  {
    return run_filter(false);
  }
  if (keyevent && keyevent->get_keyval() == SCIRun_BackSpace)
  {
    // delete last point placed
    seed_lock_.lock();
    last_seed_window_ = painter_->cur_window_;

    int sz = seeds_.size();
    if (sz > 0)
    {
      last_seed_index_ = sz - 1;
      seeds_.erase(seeds_.begin() + last_seed_index_);
      seed_change_callback(last_seed_index_);
      last_seed_index_ = last_seed_index_ - 1;
    }
    seed_lock_.unlock();
    painter_->redraw_all();
    return STOP_E;
  }
  return CONTINUE_E;
}


void
PolylineTool::draw_gl(SliceWindow &window)
{
  size_t i;

  // Transform the seeds into window space.
  seed_lock_.lock();
  vector<Point> points = seeds_;
  seed_lock_.unlock();

  // Set the color
  glColor4f(0.3, 1.0, 1.0, 1.0);

  // Draw the line loop.
  glLineWidth(2.0);
  glBegin(GL_LINES);
  {
    for (i = 0; i < points.size(); i++)
    {
      glVertex3d(points[i].x(), points[i].y(), points[i].z());
      if (i > 0 && i < points.size()-1)
      {
        glVertex3d(points[i].x(), points[i].y(), points[i].z());
      }
    }
  }
  glEnd();

  glLineWidth(1.0);

  // Draw the implicit close with a different appearance.
  if (points.size() > 2)
  {
    // TODO: Stipple this one.
    glBegin(GL_LINES);
    glVertex3d(points[points.size()-1].x(), points[points.size()-1].y(), points[points.size()-1].z());
    glVertex3d(points[0].x(), points[0].y(), points[0].z());
    glEnd();
  }

  // Draw the control points.
  glPointSize(4.0);

  glBegin(GL_POINTS);
  for (i = 0; i < points.size(); i++)
  {
    glVertex3d(points[i].x(), points[i].y(), points[i].z());
  }
  glEnd();

  glPointSize(1.0);

  CHECK_OPENGL_ERROR();
}


void
PolylineTool::seed_change_callback(int index)
{
}


BaseTool::propagation_state_e
PolylineTool::run_filter(bool erase)
{
  if (!painter_->check_for_active_label_volume("Polyline fill"))
  {
    return STOP_E;
  }

  if (seeds_.size() == 0)
  {
    painter_->set_status("No polyline to rasterize.");
    return STOP_E;
  }

  painter_->volume_lock_.lock();

  // Get the current volume.
  NrrdVolumeHandle &vol = painter_->current_volume_;
  vol->lock.lock();

  if (last_seed_window_ == NULL) return STOP_E;

  // Get the slice.
  SliceWindow *window = last_seed_window_;
  SLIVR::Plane plane(SLIVR::Point(window->center_.x(), window->center_.y(),
    window->center_.z()), SLIVR::Vector(window->normal_.x(),
    window->normal_.y(),window->normal_.z()));
  VolumeSliceHandle slice = vol->get_volume_slice(plane);
  NrrdDataHandle cnrrd = slice->nrrd_handle_;
  unsigned int clabel = vol->label_;

  // Get the mask slice if it's set.
  VolumeSliceHandle mask_slice;
  unsigned int mlabel = 0;
  NrrdDataHandle mnrrd;
  if (painter_->check_for_valid_mask("Polyline fill"))
  {
    mask_slice = painter_->mask_volume_->get_volume_slice(plane);
    mnrrd = mask_slice->nrrd_handle_;
    mlabel = painter_->mask_volume_->label_;
  }

  vol->lock.unlock();
  painter_->volume_lock_.unlock();

  // No slice, nothing to do.
  if (!slice.get_rep()) {
    return CONTINUE_E;
  }

  const int axis = slice->get_axis();

  vector<vector<float> > seeds;
  for (size_t i = 0; i < seeds_.size(); i++)
  {
    vector<double> newseed = vol->point_to_index(seeds_[i]);
    vector<float> newseedf;
    for (size_t j = 0; j < newseed.size(); j++)
    {
      if (j != axis) { newseedf.push_back(newseed[j] - 0.5); }
    }
    seeds.push_back(newseedf);
  }

  // Rasterize the polyline here.
  rasterize(cnrrd, clabel, mnrrd, mlabel, seeds, erase);

  // Put this slice back.
  painter_->volume_lock_.lock();
  vol->lock.lock();

  // Clear the content for nrrdSplice
  if (vol->nrrd_handle_->nrrd_->content) {
    vol->nrrd_handle_->nrrd_->content[0] = 0;
  }

  const vector<int> window_center = vol->world_to_index(window->center_);

  NrrdDataHandle undoslice = new NrrdData();
  nrrdSlice(undoslice->nrrd_, vol->nrrd_handle_->nrrd_,
            axis, window_center[axis]);
  
  if (nrrdSplice(vol->nrrd_handle_->nrrd_,
                 vol->nrrd_handle_->nrrd_,
                 slice->nrrd_handle_->nrrd_,
                 axis, window_center[axis])) {
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
  
  UndoHandle undo =
    new UndoReplaceSlice(painter_,
                         erase?"Undo Polyline Erase":"Undo Polyline Fill",
                         vol, undoslice,
                         (int)axis, window_center[axis]);
  painter_->push_undo(undo);

  // Clear the slice pointers.
  slice = 0;
  mask_slice = 0;

  vol->set_dirty();
  vol->lock.unlock();
  painter_->volume_lock_.unlock();

  painter_->extract_all_window_slices();
  painter_->redraw_all();

  return STOP_E;
}


// Boundary conditions
bool
PolylineTool::check_on_boundary(const vector<vector<float> > &seeds,
                                float x, float y)
{
  const float EPSILON = 1.e-6;

  for (size_t i0 = 0; i0 < seeds.size(); i0++)
  {
    const size_t i1 = (i0+1)%seeds.size();

    const float x0 = seeds[i0][1];
    const float y0 = seeds[i0][2];
    const float x1 = seeds[i1][1];
    const float y1 = seeds[i1][2];

    const float d = y0 - y1;
    if (fabs(d) < EPSILON)
    {
      // parallel
      if (fabs(y0-y) < EPSILON)
      {
        if (x0 < x1)
        {
          if (x > x0-EPSILON && x < x1+EPSILON) return true;
        }
        else
        {
          if (x > x1-EPSILON && x < x0+EPSILON) return true;
        }
      }
      continue;
    }
    else if (y0 < y1)
    {
      if (y < y0-EPSILON || y > y1+EPSILON) continue;
    }
    else
    {
      if (y < y1-EPSILON || y > y0+EPSILON) continue;
    }

    // On wrong side.
    const float u = (x1 - x0) * (y0 - y) - (y1 - y0) * (x0 - x);
    if (fabs(u) < 1.e-6) return true;
  }

  return false;
}


bool
PolylineTool::check_crossings(const vector<vector<float> > &seeds,
                              float x, float y)
{
  const float EPSILON = 1.e-6;

  int inside = 0;
  for (int i = 0; i < 7; i++)
  {
    int crosses = 0;
    const float x3 = x;
    const float y3 = y;
    const float x4 = x + rand()%17 - 8;
    const float y4 = y + rand()%17 - 8;

    for (size_t i0 = 0; i0 < seeds.size(); i0++)
    {
      const size_t i1 = (i0+1)%seeds.size();
      
      const float x1 = seeds[i0][1];
      const float y1 = seeds[i0][2];
      const float x2 = seeds[i1][1];
      const float y2 = seeds[i1][2];

      const float d = (y4 - y3) * (x2 - x1) - (x4 - x3) * (y2 - y1);
      if (fabs(d) < EPSILON) continue; // parallel

      const float u = ((x4 - x3) * (y1 - y3) - (y4 - y3) * (x1 - x3)) / d;
      const float v = ((x2 - x1) * (y1 - y3) - (y2 - y1) * (x1 - x3)) / d;

      if (u < -EPSILON || u > 1+EPSILON) continue; // not on segment
      if (v < -EPSILON) continue; // left
      
      crosses++;
    }
    if (crosses & 1) inside++;
  }

  return (inside > 2);
}


void
PolylineTool::compute_bbox(const vector<vector<float> > &seeds,
                           float &x0, float &y0, float &x1, float &y1)
{
  for (size_t i = 0; i < seeds.size(); i++)
  {
    if (seeds[i][1] < x0) x0 = seeds[i][1];
    if (seeds[i][1] > x1) x1 = seeds[i][1];
    if (seeds[i][2] < y0) y0 = seeds[i][2];
    if (seeds[i][2] > y1) y1 = seeds[i][2];
  }
}


void
PolylineTool::rasterize(NrrdDataHandle dnrrd, unsigned int dlabel,
                        NrrdDataHandle mnrrd, unsigned int mlabel,
                        const vector<vector<float> > &seeds,
                        bool erase)
{
  // Set up the data pointers.
  label_type *dstdata = (label_type *)dnrrd->nrrd_->data;
  const bool masking = mnrrd.get_rep() != NULL;
  label_type *mskdata = 0;
  if (masking) { mskdata = (label_type *)mnrrd->nrrd_->data; }

  const size_t isize = dnrrd->nrrd_->axis[1].size;
  const size_t jsize = dnrrd->nrrd_->axis[2].size;

  float x0 = (float)isize;
  float y0 = (float)jsize;
  float x1 = 0.0;
  float y1 = 0.0;
  compute_bbox(seeds, x0, y0, x1, y1);
  int i0 = (int)floor(x0);
  int j0 = (int)floor(y0);
  int i1 = (int)ceil(x1) + 1;
  int j1 = (int)ceil(y1) + 1;
  if (i0 < 0) i0 = 0;
  if (j0 < 0) j0 = 0;
  if (i1 > (int)isize) i1 = (int)isize;
  if (j1 > (int)jsize) j1 = (int)jsize;
  for (int j = j0; j < j1; j++)
  {
    for (int i = i0; i < i1; i++)
    {
      if (!masking || (mskdata[j * isize + i] & mlabel))
      {
        // We're inside the polygon if we cross an odd number of line
        // segments to get out.
        if (check_crossings(seeds, (float)i, (float)j) ||
            check_on_boundary(seeds, (float)i, (float)j))
        {
          if (erase)
          {
            dstdata[j * isize + i] &= ~dlabel;
          }
          else
          {
            dstdata[j * isize + i] |= dlabel;
          }
        }
      }
    }
  }
}


}

