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
//    File   : BrushTool.cc
//    Author : McKay Davis
//    Date   : Oct 2006

#include <Applications/Seg3D/Painter.h>
#include <Applications/Seg3D/Seg3DFrame.h>
#include <Applications/Seg3D/BrushTool.h>

#include <Applications/Seg3D/GuiCode/seg3devents.h>

#include <sci_comp_warn_fixes.h>

#include <stdlib.h>
#include <math.h>

#include <sci_gl.h>

#include <map>
#include <typeinfo>
#include <iostream>

namespace SCIRun {


// Table driven brush radius from abstract 'size' parameter.  Using
// the size directly as the radius results in a large aliased circle
// (fixed by adding an experimentally determined 0.15 to the larger
// sizes.  It also doesn't cover small touchup brushes well so more
// granularity is added at the smaller sizes.
static double static_radius_table_[41] =
  { 0.0,

    0.5,
    1.0,
    1.5,
    2.0,
    2.6,

    3.05,
    3.6,
    4.05,
    4.4,
    5.15,

    6.15,
    7.15,
    8.15,
    9.15,
    10.15,

    11.15,
    12.15,
    13.15,
    14.15,
    15.15,

    16.15,
    17.15,
    18.15,
    19.15,
    20.15,

    21.15,
    22.15,
    23.15,
    24.15,
    25.15,

    26.15,
    27.15,
    28.15,
    29.15,
    30.15,

    31.15,
    32.15,
    33.15,
    34.15,
    35.15,
  };
    

BrushTool::BrushTool(Painter *painter) :
  BaseTool("Paint Brush"),
  PointerTool("Paint Brush"),
  painter_(painter),
  window_(0),
  slice_(0),
  mask_slice_(0),
  axis_(0),
  value_(0),
  mask_value_(0),
  last_index_(),
  radius_(-1.0),
  radius_table_index_(0),
  draw_cursor_(1),
  erasing_(false),
  force_erasing_(false)
{
}


BrushTool::~BrushTool()
{
}


BaseTool::propagation_state_e
BrushTool::process_event(event_handle_t event)
{
  UpdateBrushRadiusEvent *rchange =
    dynamic_cast<UpdateBrushRadiusEvent *>(event.get_rep());
  if (rchange)
  {
    const int nradius = (int)rchange->get_radius();
    if (nradius != radius_table_index_)
    {
      set_radius(nradius);
    }
    force_erasing_ = rchange->get_force_erasing();
  }
    
  RedrawSliceWindowEvent *redraw = 
    dynamic_cast<RedrawSliceWindowEvent *>(event.get_rep());
  if (redraw) {
    draw_gl(redraw->get_window());
  }

  return CONTINUE_E;
}


BaseTool::propagation_state_e
BrushTool::pointer_down(int b, int x, int y, unsigned int m, int t)
{
  if (m & EventModifiers::SHIFT_E) { return CONTINUE_E; }

  if ((b == 1 || b == 3) &&
      !painter_->check_for_active_label_volume("Brush tool"))
  {
    return STOP_E;
  }

  NrrdVolumeHandle &vol = painter_->current_volume_;
  if (vol->button_ && !vol->button_->layer_visible_())
  {
    painter_->set_status("Painting on non-visible label.  Making visible.");
    vol->button_->layer_visible_ = true;
  }

  if (b == 1 || (b == 3 && vol->label_)) {
    window_ = painter_->cur_window_;
    slice_ = 0;
    painter_->volume_lock_.lock();

    SLIVR::Plane plane(SLIVR::Point(window_->center_.x(),window_->center_.y(),
      window_->center_.z() ), SLIVR::Vector(window_->normal_.x(),
      window_->normal_.y(),window_->normal_.z()));
    slice_ = vol->get_volume_slice(plane);

    mask_slice_ = 0;
    mask_value_ = 0;
    if (painter_->check_for_valid_mask("Brush tool"))
    {
      mask_slice_ = painter_->mask_volume_->get_volume_slice(plane);
      mask_value_ = painter_->mask_volume_->label_;
    }

    painter_->volume_lock_.unlock();

    if (!slice_.get_rep()) {
      return CONTINUE_E;
    }

    axis_ = slice_->get_axis();
    if ((b == 1 && !force_erasing_) ||
        (b == 3 && force_erasing_))
    {
      // Paint
      value_ = vol->label_;
      erasing_ = false;
    }
    else
    {
      // Erase
      value_ = vol->label_;
      erasing_ = true;
    }

    last_index_ = vol->world_to_index(painter_->pointer_pos_);
    last_index_.erase(last_index_.begin()+axis_);
    splat(slice_->nrrd_handle_->nrrd_, radius_, 
          last_index_[1], last_index_[2]);

    slice_->set_region_dirty(last_index_[1], last_index_[2],
                             last_index_[1]+1, last_index_[2]+1,
                             Ceil(radius_));

    painter_->redraw_all();    
    return STOP_E;
  } else if (b == 4) {
    set_radius(Min(radius_table_index_ + 1, 40));
    wxPanel *panel = Painter::global_seg3dframe_pointer_->CurrentToolPanel();
    if (panel)
    {
      wxCommandEvent wxevent(wxEVT_BRUSH_RADIUS_CHANGE, wxID_ANY);
      wxevent.SetInt(radius_table_index_);
      wxPostEvent(panel, wxevent);
    }
    draw_cursor_ = true;
    painter_->redraw_all();    
    return STOP_E;
  } else if (b == 5) {
    set_radius(Max(radius_table_index_ - 1, 1));
    wxPanel *panel = Painter::global_seg3dframe_pointer_->CurrentToolPanel();
    if (panel)
    {
      wxCommandEvent wxevent(wxEVT_BRUSH_RADIUS_CHANGE, wxID_ANY);
      wxevent.SetInt(radius_table_index_);
      wxPostEvent(panel, wxevent);
    }
    draw_cursor_ = true;
    painter_->redraw_all();
    return STOP_E;
  }

  return CONTINUE_E;
}

BaseTool::propagation_state_e
BrushTool::pointer_up(int b, int x, int y, unsigned int m, int t)
{  
  if (!window_ || !slice_.get_rep()) {
    return CONTINUE_E;
  }

  painter_->volume_lock_.lock();
  NrrdVolumeHandle vol = slice_->volume_;  
  vol->lock.lock();

  if (b == 1 || (slice_.get_rep() && b == 3 && vol->label_))
  {
    const vector<int> window_center = vol->world_to_index(window_->center_);

    if (vol->nrrd_handle_->nrrd_->content) {
      vol->nrrd_handle_->nrrd_->content[0] = 0;
    }

    NrrdDataHandle undoslice = new NrrdData();
    if( nrrdSlice(undoslice->nrrd_, vol->nrrd_handle_->nrrd_,
		  axis_, window_center[axis_]) ) {
      vol->lock.unlock();
      painter_->volume_lock_.unlock();
      char *err = biffGetDone(NRRD);

      cerr << to_string(__FILE__)
	   << string("Error on line #") 
           << to_string(__LINE__)
           << string(" executing nrrd command: nrrdSlice \n")
           << string("Message: ") 
           << err
           << std::endl;

      painter_->set_status( string("Arithmetic tool error performing nrrd operation failed: ") + err);
      free(err);
      return STOP_E;
    }

    if (nrrdSplice(vol->nrrd_handle_->nrrd_,
                   vol->nrrd_handle_->nrrd_,
                   slice_->nrrd_handle_->nrrd_,
                   axis_, window_center[axis_])) {
      vol->lock.unlock();
      painter_->volume_lock_.unlock();
      char *err = biffGetDone(NRRD);

      cerr << string("Error on line #") 
           << to_string(__LINE__)
           << string(" executing nrrd command: nrrdSplice \n")
           << string("Message: ") 
           << err
           << std::endl;

      painter_->set_status( string("Arithmetic tool error performing nrrd operation failed: ") + err);
      free(err);
      return STOP_E;
    }

    UndoHandle undo =
      new UndoReplaceSlice(painter_,
                           erasing_?"Undo Erase Stroke":"Undo Brush Stroke",
                           vol, undoslice,
                           (int)axis_, window_center[axis_]);
    painter_->push_undo(undo);

    slice_ = 0;
    mask_slice_ = 0;
    vol->set_dirty();
    vol->lock.unlock();
    painter_->volume_lock_.unlock();
    painter_->extract_all_window_slices();
    // Only need to do this if brush replaces.
    //painter_->set_all_slices_tex_dirty();

    painter_->redraw_all();
    return STOP_E;
  }

  vol->lock.unlock();
  painter_->volume_lock_.unlock();
 
  return CONTINUE_E;
}


BaseTool::propagation_state_e
BrushTool::pointer_motion(int b, int x, int y, unsigned int m, int t)
{
  draw_cursor_ = !(m & EventModifiers::SHIFT_E);
  painter_->redraw_all();
  if (!window_ || !slice_.get_rep()) {
    return CONTINUE_E;
  }
  
  NrrdVolumeHandle vol = slice_->volume_;
  if (!vol.get_rep()) {
    return CONTINUE_E;
  }

  vol->lock.lock();

  if (b == 1 || (slice_.get_rep() && b == 3 && vol->label_))
  {
    vector<int> index = vol->world_to_index(painter_->pointer_pos_);
    index.erase(index.begin()+axis_);
    if (line(slice_->nrrd_handle_->nrrd_, radius_,
             last_index_[1], last_index_[2],
             index[1], index[2], true))
    {
      slice_->set_region_dirty(last_index_[1], last_index_[2], 
                               index[1], index[2], Ceil(radius_));

      painter_->redraw_all();
    }
    last_index_ = index;
    vol->lock.unlock();
    return STOP_E;
  }

  vol->lock.unlock();

  return CONTINUE_E;
}


void
BrushTool::set_radius(int r)
{
  radius_table_index_ = r;
  radius_ = static_radius_table_[radius_table_index_];
  cursor_cache_.clear();

  const float r2 = radius_ * radius_;
  const int wid = Ceil(radius_);
  for (int y = -wid; y <= wid; ++y) {
    for (int x = -wid; x <= wid; ++x) {
      const float dist = x*x+y*y;
      if (dist > r2) continue;
      // Right
      if ((x+1)*(x+1)+(y+0)*(y+0) > r2) {
        cursor_cache_.push_back(x+1);
        cursor_cache_.push_back(y);
        cursor_cache_.push_back(x+1);
        cursor_cache_.push_back(y+1);
      }
      // Top
      if ((x+0)*(x+0)+(y+1)*(y+1) > r2) {
        cursor_cache_.push_back(x);
        cursor_cache_.push_back(y+1);
        cursor_cache_.push_back(x+1);
        cursor_cache_.push_back(y+1);
      }
      // Left 
      if ((x-1)*(x-1)+(y+0)*(y+0) > r2) {
        cursor_cache_.push_back(x);
        cursor_cache_.push_back(y);
        cursor_cache_.push_back(x);
        cursor_cache_.push_back(y+1);
      }
      // Bottom
      if ((x+0)*(x+0)+(y-1)*(y-1) > r2) {
        cursor_cache_.push_back(x);
        cursor_cache_.push_back(y);
        cursor_cache_.push_back(x+1);
        cursor_cache_.push_back(y);
      }
    }
  }
}


void
BrushTool::draw_gl(SliceWindow &window)
{
  if (!draw_cursor_) return;
  if (&window != painter_->cur_window_) return;
  NrrdVolumeHandle &vol = painter_->current_volume_;
  if (!vol.get_rep()) return;

  Color rgb(1.0, 0.0, 0.0);
  if (vol->label_)
  {
    Color c = vol->get_label_color();
    const double m = 0.7;
    rgb = Color(c.r() * m, c.g() * m, c.b() * m);
  }
  glColor4f(rgb.r(), rgb.g(), rgb.b(), 1.0);

  vector<double> upv = vol->vector_to_index(window.y_dir());
  upv[max_vector_magnitude_index(upv)] /=
    Abs(upv[max_vector_magnitude_index(upv)]);
  const Vector up = vol->index_to_vector(upv);

  vector<double> rightv = vol->vector_to_index(window.x_dir());
  rightv[max_vector_magnitude_index(rightv)] /= 
    Abs(rightv[max_vector_magnitude_index(rightv)]);
  const Vector right = vol->index_to_vector(rightv);

  const Point center =
    vol->index_to_world(vol->world_to_index(painter_->pointer_pos_));

  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();

  double mat[16];
  mat[0] = right.x();
  mat[1] = right.y();
  mat[2] = right.z();
  mat[3] = 0.0;
  mat[4] = up.x();
  mat[5] = up.y();
  mat[6] = up.z();
  mat[7] = 0.0;
  mat[8] = 0.0;
  mat[9] = 0.0;
  mat[10] = 1.0;
  mat[11] = 0.0;
  mat[12] = center.x();
  mat[13] = center.y();
  mat[14] = center.z();
  mat[15] = 1.0;

  glMultMatrixd(mat);

  glLineWidth(2.0);

  glVertexPointer(2, GL_INT, 0, &(cursor_cache_.front()));
  glEnableClientState(GL_VERTEX_ARRAY);
  glDrawArrays(GL_LINES, 0, cursor_cache_.size()/2);

  glLineWidth(1.0);

  glPopMatrix();
}


void
BrushTool::splat(Nrrd *nrrd, double radius, int x, int y)
{
  if (mask_slice_.get_rep())
  {
    splat_mask(nrrd, radius, x, y);
  }
  else
  {
    splat_nomask(nrrd, radius, x, y);
  }
}


// Non-masking splat
void
BrushTool::splat_nomask(Nrrd *nrrd, double radius, int x0, int y0)
{ 
  ASSERT(nrrd->type == LabelNrrdType);
  label_type *slicedata = (label_type *)nrrd->data;

  const unsigned int value = erasing_ ? (~value_) : value_;

  const unsigned int wid = Ceil(radius);
  const double r2 = radius * radius;
  for (int y = y0-wid; y <= int(y0+wid); ++y)
  {
    for (int x = x0-wid; x <= int(x0+wid); ++x)
    {
      if (x >= 0 && x < int(nrrd->axis[1].size) &&
          y >= 0 && y < int(nrrd->axis[2].size)) 
      {
        const int dx = x0-x;
        const int dy = y0-y;
        const double dist2 = double(dx*dx + dy*dy);
        if (dist2 <= r2)
        {
          size_t position = (y * nrrd->axis[1].size + x) * nrrd->axis[0].size;
          if (erasing_)
          {
            slicedata[position] &= value;
          }
          else
          {
            slicedata[position] |= value;
          }
        }
      }
    }
  }
}


// Masking splat
void
BrushTool::splat_mask(Nrrd *nrrd, double radius, int x0, int y0)
{ 
  ASSERT(nrrd->type == LabelNrrdType);
  label_type *slicedata = (label_type *)nrrd->data;

  const unsigned int value = erasing_ ? (~value_) : value_;

  Nrrd *mask = mask_slice_->nrrd_handle_->nrrd_;
  label_type *maskdata = (label_type *)mask->data;
  const unsigned int maskvalue = mask_value_;

  const unsigned int wid = Ceil(radius);
  const double r2 = radius * radius;
  for (int y = y0-wid; y <= int(y0+wid); ++y)
  {
    for (int x = x0-wid; x <= int(x0+wid); ++x)
    {
      if (x >= 0 && x < int(nrrd->axis[1].size) &&
          y >= 0 && y < int(nrrd->axis[2].size)) 
      {
        const int dx = x0-x;
        const int dy = y0-y;
        const double dist2 = double(dx*dx + dy*dy);
        if (dist2 <= r2)
        {
          size_t position = (y * nrrd->axis[1].size + x) * nrrd->axis[0].size;
          if (maskdata[position] & maskvalue)
          {
            if (erasing_)
            {
              slicedata[position] &= value;
            }
            else
            {
              slicedata[position] |= value;
            }
          }
        }
      }
    }
  }
}


bool
BrushTool::line(Nrrd *nrrd, double radius,
                int x0, int y0, int x1, int y1, bool first)
{
  int dx = x1 - x0;
  int dy = y1 - y0;
  int sx = 1;
  int sy = 1;
  int frac = 0;
  if (dx < 0) { 
    dx = -dx;  
    sx = -1;
  } 
  if (dy < 0) { 
    dy = -dy;
    sy = -1; 
  } 

  // Our motion got rounded away, nothing to do.
  if (dx == 0 && dy == 0) return false;

  dy <<= 1;
  dx <<= 1;
  if (first) splat(nrrd, radius, x0, y0);
  if (dx > dy)
  {
    frac = dy - (dx >> 1);
    while (x0 != x1)
    {
      if (frac >= 0) {
        y0 += sy;
        frac -= dx;
      }
      x0 += sx;
      frac += dy;
      splat(nrrd, radius, x0, y0);
    }
  }
  else
  {
    frac = dx - (dy >> 1);
    while (y0 != y1)
    {
      if (frac >= 0) {
        x0 += sx;
        frac -= dy;
      }
      y0 += sy;
      frac += dx;
      splat(nrrd, radius, x0, y0);
    }
  }
  return true;
}


}
