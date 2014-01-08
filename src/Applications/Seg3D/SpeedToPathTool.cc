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
//    File   : SpeedToPathTool.cc
//    Author : David Brayford
//    Date   : May 2008

#include <Applications/Seg3D/Painter.h>

#include <Applications/Seg3D/SpeedToPathTool.h>
#include <Core/Geometry/CompGeom.h>
#include <Core/Events/keysyms.h>

#include <sci_defs/insight_defs.h>
#include <itkDiscreteGaussianImageFilter.h>
#include <itkGradientMagnitudeImageFilter.h>
#include <itkRescaleIntensityImageFilter.h>

namespace SCIRun {


SpeedToPathTool::SpeedToPathTool(const string &name, Painter *painter) :
  PolylineTool(name, painter)
{
}


SpeedToPathTool::~SpeedToPathTool()
{
}


BaseTool::propagation_state_e
SpeedToPathTool::process_event(event_handle_t event)
{
  if (dynamic_cast<SetDataLayerEvent *>(event.get_rep())) 
  {
    generate_speed_function();
    source_volume_ = painter_->current_volume_;

    return STOP_E;
  }

  if (dynamic_cast<FinishEvent *>(event.get_rep())) {
    return run_filter(false);
  }

  return PolylineTool::process_event(event);
}


void
SpeedToPathTool::draw_gl(SliceWindow &window)
{
  seed_lock_.lock();

  if (source_volume_.get_rep())
  {
    const vector<int> window_center =
      source_volume_->world_to_index(window.center_);

    const int sliceloc = window_center[window.axis_+1];
    window_slice_map_t::iterator loc = window_slice_map_.find(&window);
    if (loc == window_slice_map_.end())
    {
      window_slice_map_[&window] = sliceloc;
    }
    else if ((*loc).second != sliceloc)
    {
      window_slice_map_[&window] = sliceloc;
      seed_change_callback(-6);
    }
  }

  // Set the color
  glColor4f(0.21, 0.7, 0.7, 1.0);

  // Draw the closed path.
  glLineWidth(2.0);
  glBegin(GL_LINE_LOOP);
  {
    for (size_t i = 0; i < paths_.size(); i++)
    {
      for (size_t j = 0; j < paths_[i].size(); j++)
      {
        glVertex3d(paths_[i][j].x(), paths_[i][j].y(), paths_[i][j].z());
      }
    } 
  }
  glEnd();
  glLineWidth(1.0);

  size_t i;

  // Transform the seeds into window space.

  vector<Point> points = seeds_;

  seed_lock_.unlock();

  // Draw the control points.
  glColor4f(0.3, 1.0, 1.0, 1.0);
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


int
SpeedToPathTool::compute_nearest_segment_index(const Point &curpos, int axis)
{
  // No seed point near.  Check for insertion of new point into an
  // existing line segment.
  int index = -1;
  const double max_dist2_to_check = 25.0;
  double mind2 = max_dist2_to_check + 1.0;
  const Point sppos = last_seed_window_->world_to_screen(curpos);

  for (size_t i = 0; i < paths_.size(); i++)
  {
    for (size_t j = 0; j < paths_[i].size(); j++)
    {
      Point seedi0 = paths_[i][j];
      seedi0(axis) = 0.0;

      int i1 = i;
      int j1 = (j + 1) % paths_[i].size();
      if (j1 == 0)
      {
        i1 = (i+1) % paths_.size();
      }
      Point seedi1 = paths_[i1][j1];
      seedi1(axis) = 0.0;
      
      const Point sseedi0 = last_seed_window_->world_to_screen(seedi0);
      const Point sseedi1 = last_seed_window_->world_to_screen(seedi1);
      const double d2 = distance_to_line2(sppos, sseedi0, sseedi1);
      if (d2 < mind2 && d2 < max_dist2_to_check)
      {
        mind2 = d2;
        index = i+1;
      }
    }
  }

  return index;
}


void 
SpeedToPathTool::seed_change_callback(int index)
{ 
  if (!source_volume_.get_rep())
  {
    if (painter_->current_volume_.get_rep() &&
        !painter_->current_volume_->label_ &&
        painter_->current_volume_->data_min_ > 0.0 &&
        painter_->current_volume_->data_max_ <= 1.0)
    {
      source_volume_ = painter_->current_volume_;
    }
    else
    {
      painter_->set_status("No speed function.");
      return;
    }
  }

  if (index < 0 || paths_.size() != seeds_.size())
  {
    paths_.clear();
    paths_.resize(seeds_.size());

    // No path.
    if ( seeds_.size() < 2 )
    {
      return;
    }

    for (size_t i = 0; i < seeds_.size(); i++)
    {
      compute_path(i);
    }
  }
  else
  {
    int prev = index-1;
    if (prev < 0) prev = paths_.size()-1;
    compute_path(prev);
    compute_path(index);
  }
}


void
SpeedToPathTool::compute_path(size_t seed_index)
{
  // Polyline.
  paths_[seed_index].push_back(seeds_[seed_index]);
}


BaseTool::propagation_state_e
SpeedToPathTool::run_filter(bool erase)
{
  if (!painter_->check_for_active_label_volume("Speedline fill"))
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
  SLIVR::Plane plane(SLIVR::Point(window->center_.x(),window->center_.y(),
    window->center_.z()), SLIVR::Vector(window->normal_.x(),
    window->normal_.y(),window->normal_.z()));
  VolumeSliceHandle slice = vol->get_volume_slice(plane);
  NrrdDataHandle cnrrd = slice->nrrd_handle_;
  unsigned int clabel = vol->label_;

  // Get the mask slice if it's set.
  VolumeSliceHandle mask_slice;
  unsigned int mlabel = 0;
  NrrdDataHandle mnrrd;
  if (painter_->check_for_valid_mask("Speedline fill"))
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

  vector<vector<float> > path;
  for (size_t i = 0; i < paths_.size(); i++)
  {
    for (size_t j = 0; j < paths_[i].size(); j++)
    {
      {
        vector<double> newseed = vol->point_to_index(paths_[i][j]);
        vector<float> newseedf;
        for (size_t j = 0; j < newseed.size(); j++)
        {
          if (j != axis) { newseedf.push_back(newseed[j] - 0.5); }
        }
        path.push_back(newseedf);
      }
    }
  }

  // Rasterize the polyline here.
  rasterize(cnrrd, clabel, mnrrd, mlabel, path, erase);

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
                         erase?"Undo Speedline Erase":"Undo Speedline Fill",
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


void
SpeedToPathTool::generate_speed_function()
{
  painter_->set_status("Creating Speed Function.");

  // Gaussian blur.
  typedef itk::DiscreteGaussianImageFilter< ITKImageFloat3D, ITKImageFloat3D > FilterTypeGuassian;
  VolumeFilter<FilterTypeGuassian> blur_filter;

  blur_filter.set_volume( painter_->copy_current_layer(" Speed Function") );

  // Set the filter parameters.
  const double variance = 6.0;
  const int max_kernel_width = 40;

  blur_filter->SetVariance(variance);
  blur_filter->SetUseImageSpacingOff();
  if (max_kernel_width > 0)
  {
    blur_filter->SetMaximumKernelWidth(max_kernel_width);
  }

  painter_->start_progress();
  blur_filter();

  // Gradient magnitude
  typedef itk::GradientMagnitudeImageFilter< ITKImageFloat3D, ITKImageFloat3D > FilterType;
  VolumeFilter<FilterType> gm_filter;

  gm_filter.set_volume(painter_->current_volume_);
  gm_filter();
  painter_->current_volume_->reset_data_range();

  ITKDatatypeHandle image_handle =
    nrrd_to_itk_image(painter_->current_volume_->nrrd_handle_);

  // Replace InputImageType with actual image type non typedef
  ITKImageFloat3D::Pointer pImage =
    dynamic_cast<ITKImageFloat3D *>(image_handle->data_.GetPointer());
  
  if (!pImage)
  {
    return;
  } 
  // compute normalization
  // epsilon value used to avoid zero this value might need modifying
  double epsilon = 1.0e-7;
  double max_value = 1.0;

  typedef itk::RescaleIntensityImageFilter < ITKImageFloat3D, ITKImageFloat3D > RescaleFilterType;

  RescaleFilterType::Pointer rescale_filter = RescaleFilterType::New();
  rescale_filter->SetInput( pImage );
  rescale_filter->SetOutputMinimum(epsilon);
  rescale_filter->SetOutputMaximum(max_value);
  rescale_filter->Update();

  SCIRun::ITKDatatypeHandle output_img = new SCIRun::ITKDatatype();
  output_img->data_ = rescale_filter->GetOutput();
  painter_->current_volume_->nrrd_handle_ = 
    itk_image_to_nrrd<ITKImageFloat3D::PixelType, 3>(output_img);

  painter_->current_volume_->reset_data_range();

  painter_->finish_progress();

  UndoHandle undo =
    new UndoReplaceLayer(painter_, "Undo Speed Function", 0,
                         painter_->current_volume_, 0);
  painter_->push_undo(undo);

  painter_->extract_all_window_slices();
  painter_->redraw_all();
}


}
