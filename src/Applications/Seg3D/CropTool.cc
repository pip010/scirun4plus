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
//    File   : CropTool.cc
//    Author : McKay Davis
//    Date   : Sun Oct  1 23:22:04 2006



#include <Applications/Seg3D/Painter.h>
#include <Applications/Seg3D/CropTool.h>
#include <Core/Events/EventManager.h>
#include <Core/Events/LoadVolumeEvent.h>
#include <sci_gl.h>

#include <Applications/Seg3D/GuiCode/cropvolume.h>
#include <Applications/Seg3D/GuiCode/seg3devents.h>

#ifdef _WIN32
#  undef SCISHARE
#  define SCISHARE __declspec(dllimport)
#else
#  define SCISHARE
#endif

namespace SCIRun {


CropTool::CropTool(Painter *painter) : 
  BaseTool("Crop"),
  PointerTool("Crop"),
  painter_(painter),
  pick_(0)
{
  if (painter_->current_volume_.get_rep()) {
    minmax_[1] = painter_->current_volume_->max_index();
  } else {
    minmax_[1] = vector<int>(4, 0);
  }

  minmax_[0] = vector<int>(minmax_[1].size(), 0);

  pick_minmax_[0] = minmax_[0];
  pick_minmax_[1] = minmax_[1];
  gui_minmax_[0].resize(minmax_[1].size());
  gui_minmax_[1].resize(minmax_[1].size());
  
  Skinner::Variables *vars = painter_->get_vars();
 
  gui_minmax_[0][1] = Skinner::Var<double>(vars, "Painter::crop::min::x", minmax_[0][1]);
  gui_minmax_[0][2] = Skinner::Var<double>(vars, "Painter::crop::min::y", minmax_[0][2]);
  gui_minmax_[0][3] = Skinner::Var<double>(vars, "Painter::crop::min::z", minmax_[0][3]);
  gui_minmax_[1][1] = Skinner::Var<double>(vars, "Painter::crop::max::x", minmax_[1][1]);
  gui_minmax_[1][2] = Skinner::Var<double>(vars, "Painter::crop::max::y", minmax_[1][2]);
  gui_minmax_[1][3] = Skinner::Var<double>(vars, "Painter::crop::max::z", minmax_[1][3]);

  update_to_gui();

  // Set the range to match this layer size.
  //
  // TODO: If the user selects a differently sized data layer after
  // starting a crop the range will be incorrect as it will match the
  // first one.  Need to reupdate the range on layer selection.
  wxPanel *panel = Painter::global_seg3dframe_pointer_->CurrentToolPanel();
  if (panel)
  {
    CropSetRangeStruct *csrs = new CropSetRangeStruct;

    // Set the max range of each spinner crtl.
    if (painter_->current_volume_.get_rep())
    {
      csrs->min_x = Max(minmax_[1][1]-1, 0);
      csrs->min_y = Max(minmax_[1][2]-1, 0);
      csrs->min_z = Max(minmax_[1][3]-1, 0);
      csrs->max_x = minmax_[1][1];
      csrs->max_y = minmax_[1][2];
      csrs->max_z = minmax_[1][3];
    }
    else
    {
      // No layer selected, no range to set.  Use big rather than 0.
      csrs->min_x = 0x7fffffff;
      csrs->min_y = 0x7fffffff;
      csrs->min_z = 0x7fffffff;
      csrs->max_x = 0x7fffffff;
      csrs->max_y = 0x7fffffff;
      csrs->max_z = 0x7fffffff;
    }
    wxCommandEvent event(wxEVT_CROP_SET_RANGE, wxID_ANY);
    event.SetClientData((void *)csrs);
    wxPostEvent(panel, event);
  }

  painter_->redraw_all();
}


CropTool::~CropTool()
{
  painter_->redraw_all();
}


BaseTool::propagation_state_e
CropTool::pointer_motion(int b, int x, int y, unsigned int m, int t) 
{
  if (!painter_->current_volume_.get_rep() || !pick_)
  {
    return CONTINUE_E;
  }

  SliceWindow *window = painter_->cur_window_;
  ASSERT(window);

  vector<int> max_index = painter_->current_volume_->max_index();

  vector<double> idx = 
    painter_->current_volume_->point_to_index(painter_->pointer_pos_);

  if (pick_ == 1)  // Clicked inside crop box
  {
    for (unsigned int i=0; i<2; ++i)
    {
      for (unsigned int a=1; a<4; ++a)
      {
	// Skip the axis that the pointer is on. NOTE: the volume
	// indices have four dimensions.
	if (a-1 == window->axis_) continue;

	// Set delta such the box size does not change when it is
	// moved to the boundary.
	double delta = Clamp(idx[a]-Floor(pick_index_[a]),
			     -double(pick_minmax_[0][a]), 
			     double(max_index[a]-pick_minmax_[1][a]));

        minmax_[i][a] = (int) Floor(pick_minmax_[i][a]+delta);
      }
    }
  }
  else  // Clicked on crop box boundary
  {
    for (unsigned int i=0; i<2; ++i)
    {
      for (unsigned int a=1; a<4; ++a)
      {
	// Skip the axis that the pointer is on. NOTE: the volume
	// indices have four dimensions.
	if (a-1 == window->axis_) continue;

	int newval = Clamp(Round(idx[a]), 0, max_index[a]);

	// Make sure the point was picked and is not the same value as
	// the opposite corner.
	if (Abs(pick_dist_[i][a-1]) < 5.0 && newval != minmax_[(i+1)%2][a])
	{
	  minmax_[i][a] = newval;
	}
      }
    }
  }
  update_to_gui();
  painter_->redraw_all();
  return STOP_E;
}


BaseTool::propagation_state_e
CropTool::pointer_down(int b, int x, int y, unsigned int m, int t) 
{
  if (!painter_->current_volume_.get_rep())
  {
    painter_->set_status("No active volume to crop.");
    return CONTINUE_E;
  }

  if (b == 1 && !m)
  {
    SliceWindow *window = painter_->cur_window_;
    ASSERT(window);
    
    double units = 100.0 / window->zoom_; // world space units per pixel
    pick_ = 1;  // Assume the pointer is inside the crop box.

    for (int i = 0; i < 2; ++i)
    {
      // Save the initial crop boundary.
      pick_minmax_[i] = minmax_[i];

      Point p = painter_->current_volume_->index_to_world(minmax_[i]);

      for (int a = 0; a < 3; ++a)
      {
	if( a == window->axis_)
	{
	  // Skip when the plane and axis align as the distance can be
	  // zero but points with in the plane should be changed.
	  pick_dist_[i][a] = 999;
	}
	else
	  {
	    Vector n(a==0?1:0, a==1?1:0, a==2?1:0);

	    if (i)
	      n *= -1;

	    Plane plane(p, n);

	    pick_dist_[i][a] = plane.eval_point(painter_->pointer_pos_)/units;

	    // Pointer is near a boundary, may be on either side of it.
	    if (Abs(pick_dist_[i][a]) < 5.0)
	      pick_ |= 2;  // Set the second bit.

	    // Pointer is outside the crop boundary, may still be near the
	    // boundary though.
	    if (pick_dist_[i][a] < 0.0)
	      pick_ &= ~1;  // Clear the first bit   ~1 = 011111...
	  }
      }
    }

    pick_index_ = 
      painter_->current_volume_->point_to_index(painter_->pointer_pos_);

    return STOP_E;
  }

  return CONTINUE_E;
}



BaseTool::propagation_state_e
CropTool::pointer_up(int b, int x, int y, unsigned int m, int t) 
{
  if (pick_) {
    for (unsigned int a = 0; a < minmax_[0].size(); ++a)
      if (minmax_[0][a] > minmax_[1][a])
        std::swap(minmax_[0][a],minmax_[1][a]);
    
    pick_ = 0;
    update_to_gui();
    return STOP_E;
  }
  return CONTINUE_E;
}


BaseTool::propagation_state_e 
CropTool::process_event(event_handle_t event)
{
  RedrawSliceWindowEvent *redraw = 
    dynamic_cast<RedrawSliceWindowEvent *>(event.get_rep());
  if (redraw) {
    draw_gl(redraw->get_window());
  }

  if (dynamic_cast<FinishEvent *>(event.get_rep())) {
    finish();
  }

  return CONTINUE_E;
}


void
CropTool::update_to_gui()
{
  for (int i = 0; i < 2; ++i)
    for (int j = 1; j < 4; ++j)
    {
      gui_minmax_[i][j] = minmax_[i][j];
    }

  CropSetRangeStruct *csrs = new CropSetRangeStruct;
  csrs->min_x = minmax_[0][1];
  csrs->min_y = minmax_[0][2];
  csrs->min_z = minmax_[0][3];
  csrs->max_x = minmax_[1][1];
  csrs->max_y = minmax_[1][2];
  csrs->max_z = minmax_[1][3];

  wxCommandEvent event(wxEVT_CROP_SET_MINMAX, wxID_ANY);
  event.SetClientData((void *)csrs);
  wxPanel *panel = Painter::global_seg3dframe_pointer_->CurrentToolPanel();
  if (panel) { wxPostEvent(panel, event); }
}


void
CropTool::update_from_gui()
{
  // Need to delete the for loops for the wx version.
  for (int i = 0; i < 2; ++i)
    for (int j = 1; j < 4; ++j)
    {
      minmax_[i][j] = Round(gui_minmax_[i][j]());
    }
}


int
CropTool::draw_gl(SliceWindow &window)
{
  if (!painter_->current_volume_.get_rep()) return 0;
  update_from_gui();
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  Point ll = painter_->current_volume_->index_to_world(minmax_[0]);
  Point ur = painter_->current_volume_->index_to_world(minmax_[1]);

  Vector dia = ur - ll;

  Vector right = window.x_dir();
  right.normalize();
  right = right*dia;

  // When flipped right is left and left is right.
  if( window.flip_leftright_ ) {
    right *= -1.0;
  }

  Vector up = window.y_dir();
  up.normalize();
  up = up*dia;

  // When flipped up is down and down is up.
  if( window.flip_updown_ ) {
    up *= -1.0;
  }

  Point lr = ll+right;
  Point ul = ll+up;
  ur = ll+right+up;
  
  //  GLdouble blue[4] = { 0.1, 0.4, 1.0, 0.8 };
  //  GLdouble green[4] = { 0.5, 1.0, 0.1, 0.7 };
  //  GLdouble lt_green[4] = { 0.5, 1.0, 0.1, 0.4 };
  //  GLdouble red[4] = { 0.8, 0.2, 0.4, 0.9 };
  GLdouble grey[4] = { 0.6, 0.6, 0.6, 0.6 }; 
  GLdouble white[4] = { 1.0, 1.0, 1.0, 1.0 }; 
  GLdouble black[4] = { 0.0, 0.0, 0.0, 1.0 }; 
  GLdouble yellow[4] = { 1.0, 0.76, 0.1, 1.0 };
  GLdouble lt_yellow[4] = { 0.8, 0.5, 0.1, 1.0 };  

  GLdouble *colors[5] = { lt_yellow, yellow, black, grey, white };
  GLdouble widths[5] = { 11, 9.0, 7.0, 5.0, 1.0 }; 

  glEnable(GL_LINE_SMOOTH);
  for (int pass = 2; pass < 5; ++pass) {
    glColor4dv(colors[pass]);
    glLineWidth(widths[pass]);    

    glBegin(GL_LINE_LOOP);
    {
      glVertex3dv(&ll(0));
      glVertex3dv(&lr(0));
      glVertex3dv(&ur(0));
      glVertex3dv(&ul(0));
    }
    glEnd();
  }
  glLineWidth(1.0);
  glDisable(GL_LINE_SMOOTH);

  widths[0] = 10.0;
  widths[1] = 6.0;
  widths[2] = 2.0;
  
  glEnable(GL_POINT_SMOOTH);
  for (int pass = 0; pass < 5; ++pass) {
    glColor4dv(colors[pass]);
    glPointSize(widths[pass]);
    glBegin(GL_POINTS);
    {
      glVertex3dv(&ll(0));
      glVertex3dv(&lr(0));
      glVertex3dv(&ur(0));
      glVertex3dv(&ul(0));
    }
    glEnd();
  }

  glDisable(GL_POINT_SMOOTH);

  CHECK_OPENGL_ERROR();
  return 0; 
}


BaseTool::propagation_state_e
CropTool::finish()
{
  if (!painter_->current_volume_.get_rep())
  {
    painter_->set_status("Crop tool requires an active volume.");
    return STOP_E;
  }

  // Check for invalid crop volume, creatable with entry boxes.
  // Punt out if volume is not fixable, otherwise just fix it.
  for (unsigned int a = 0; a < minmax_[0].size(); ++a)
  {
    if (minmax_[0][a] > minmax_[1][a])
    {
      std::swap(minmax_[0][a], minmax_[1][a]);
    }
    if (minmax_[0][a] == minmax_[1][a])
    {
      painter_->set_status("Invalid crop, no volume would be left.");
      return STOP_E;
    }
  }      

  size_t *minmax[2] = { new size_t[minmax_[0].size()], 
                        new size_t[minmax_[1].size()] };

  for (int i = 0; i < 2; ++i) {
    for (unsigned int a = 0; a < minmax_[0].size(); ++a) {
      minmax[i][a] = minmax_[i][a]-(i==1?1:0);
    }
  }

  if (minmax[0][2] > minmax[1][2])
  {
    minmax[1][2] = minmax[0][2]+ 10;
  }

  NrrdDataHandle nout_handle = new NrrdData();
  if (nrrdCrop(nout_handle->nrrd_,
               painter_->current_volume_->nrrd_handle_->nrrd_,
               minmax[0], minmax[1])) {
    char *err = biffGetDone(NRRD);
    string str = string("nrrdcrop: ") + err;
    free(err);
    throw str;
  }
  
  delete[] minmax[0];
  delete[] minmax[1];

  // Make a new cropped volume similar to the old volume.
  NrrdVolumeHandle croppedvol =
    new NrrdVolume(painter_, painter_->current_volume_->name_,
                   nout_handle, painter_->current_volume_->label_);
  if (croppedvol->label_)
  {
    croppedvol->set_label_color(painter_->current_volume_->get_label_color());
  }
  
  // Replace the old volume with the new volume.
  NrrdVolumeHandle oldvol = painter_->current_volume_;
  size_t loc = painter_->remove_volume(painter_->current_volume_, false);
  loc = loc + 1; // We're off by one for some reason.
  painter_->insert_volume(croppedvol, loc);
  painter_->current_volume_ = croppedvol;

  // Reset the volume rendering.
  if (painter_->current_volume_ == painter_->current_vrender_target_ &&
      !painter_->current_vrender_target_deferred_)
  {
    LoadVolumeEvent *lve =
      new LoadVolumeEvent(painter_->current_vrender_target_->nrrd_handle_,
                          "", false);
    EventManager::add_event(lve);
  }

  // Push this undo item.
  UndoHandle undo = new UndoReplaceLayer(painter_, "Undo Crop",
                                         oldvol, croppedvol, loc);
  painter_->push_undo(undo);

  // Center the probe on the new cropped volume.
  const Point center = painter_->current_volume_->center();
  for (unsigned int i = 0; i < painter_->windows_.size(); ++ i) {
    painter_->windows_[i]->center_ = center;
  }

  // Reset the slices to match the new volume.
  painter_->rebuild_layer_buttons();
  painter_->extract_all_window_slices();
  painter_->redraw_all();

  // Autoview to the cropped volume.
  event_handle_t empty;
  painter_->Autoview(empty);

  painter_->hide_tool_panel();

  return CONTINUE_E;
}


void
CropTool::set_window_cursor(SliceWindow &window, int cursor) 
{
  return;
  //  if (painter_->event_.window_ != &window || 
  //      window.cursor_pixmap_ == cursor) return;
  //  window.cursor_pixmap_ = cursor;
  string cursor_name;
  switch (cursor) {
  case 1: cursor_name = "bottom_left_corner"; break;
  case 2: cursor_name = "bottom_right_corner"; break;
  case 3: cursor_name = "top_right_corner"; break;
  case 4: cursor_name = "top_left_corner"; break;
  case 5: cursor_name = "sb_h_double_arrow"; break;
  case 6: cursor_name = "sb_v_double_arrow"; break;
  case 7: cursor_name = "sb_h_double_arrow"; break;
  case 8: cursor_name = "sb_v_double_arrow"; break;
  case 9: cursor_name = "fleur"; break;
  case 0: cursor_name = "crosshair"; break;
  }
}


}
