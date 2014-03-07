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
//    File   : ThresholdTool.cc
//    Author : McKay Davis
//    Date   : Tue Sep 26 18:44:34 2006

#include <Applications/Seg3D/Painter.h>

#include <Applications/Seg3D/ThresholdTool.h>
#include <Applications/Seg3D/Seg3DwxGuiUtils.h>

#include <sci_gl.h>

namespace SCIRun {

ThresholdTool::ThresholdTool(Painter *painter) :
  SeedTool("ThresholdTool", painter),
  lower_threshold_(0.0),
  upper_threshold_(0.0),
  cached_volume_pointer_(NULL)
{
  label_color_ = NrrdVolume::get_next_label_color();

  // TODO:  Pick top visible by default if none are selected.
  painter_->check_for_active_data_volume("Threshold Tool");

  wxPanel *panel = Painter::global_seg3dframe_pointer_->CurrentToolPanel();
  if (panel && painter_->current_volume_.get_rep() &&
      !painter_->current_volume_->label_)
  {
    painter_->current_volume_->reset_data_range();

    lower_threshold_ = painter_->current_volume_->data_min_ +
      (painter_->current_volume_->data_max_ -
       painter_->current_volume_->data_min_) / 4.0;
    upper_threshold_ = painter_->current_volume_->data_max_ -
      (painter_->current_volume_->data_max_ -
       painter_->current_volume_->data_min_) / 4.0;

    ThresholdToolChangeStruct *info = new ThresholdToolChangeStruct;
    info->minval = painter_->current_volume_->data_min_;
    info->maxval = painter_->current_volume_->data_max_;
    info->lower = lower_threshold_;
    info->upper = upper_threshold_;

    // This code appears to cause Windows builds to hang, cause unknown at the moment.
#ifndef _WIN32
    Painter::global_seg3dframe_pointer_->
      SetStatusText(std2wx("Creating a histogram for thresholding"));
    wxBusyCursor(); // Busy cursor until this leaves scope.
    get_histogram(info->histogram);

    Painter::global_seg3dframe_pointer_->SetStatusText(std2wx(""));
#endif

    wxCommandEvent wxevent(wxEVT_THRESHOLDTOOL_CHANGE, wxID_ANY);
    wxevent.SetClientData((void *)info);
    wxPostEvent(panel, wxevent);
  }

  painter_->redraw_all();  
}


ThresholdTool::~ThresholdTool()
{
  texture_cache_map_t::iterator loc = texture_cache_.begin();
  while (loc != texture_cache_.end())
  {
    delete (*loc).second;
    ++loc;
  }
  texture_cache_.clear();
  painter_->redraw_all();  
}


BaseTool::propagation_state_e 
ThresholdTool::process_event(event_handle_t event)
{
  UpdateThresholdToolEvent *update =
    dynamic_cast<UpdateThresholdToolEvent *>(event.get_rep());
  if (update)
  {
    lower_threshold_ = update->get_lower();
    upper_threshold_ = update->get_upper();

    // TODO:  Decide if we want to swap these when lower > upper.
     if (lower_threshold_ > upper_threshold_)
    {
      const double tmp = lower_threshold_;
      lower_threshold_ = upper_threshold_;
      upper_threshold_ = tmp;
    }

    painter_->redraw_all();
    return CONTINUE_E;
  }

  return SeedTool::process_event(event);
}


void
ThresholdTool::seed_change_callback()
{
  SeedTool::seed_change_callback();

  if (!(painter_->current_volume_.get_rep() &&
        !painter_->current_volume_->label_))
  {
    return;
  }

  // Get the data pointers.
  NrrdDataHandle snrrd = painter_->current_volume_->nrrd_handle_;
  float *srcdata = (float *)snrrd->nrrd_->data;

  // Set the filter parameters.
  vector<vector<int> > iseeds;
  convert_seeds_to_indices(iseeds, painter_->current_volume_);

  if (iseeds.empty()) return;

  // Add the seeds.
  float lower = 0, upper = 0;
  for (size_t i = 0; i < iseeds.size(); ++i)
  {
    const float fillval =
      srcdata[VolumeOps::index_to_offset(snrrd->nrrd_, iseeds[i])];

    if (i == 0)
    {
      lower = upper = fillval;
    }
    else
    {
      lower = Min(lower, fillval);
      upper = Max(upper, fillval);
    }
  }

  lower_threshold_ = lower;
  upper_threshold_ = upper;

  wxPanel *panel = Painter::global_seg3dframe_pointer_->CurrentToolPanel();
  if (panel)
  {
    ThresholdToolChangeStruct *info = new ThresholdToolChangeStruct;
    info->minval = painter_->current_volume_->data_min_;
    info->maxval = painter_->current_volume_->data_max_;
    info->lower = lower_threshold_;
    info->upper = upper_threshold_;
    info->histogram.clear();

    wxCommandEvent wxevent(wxEVT_THRESHOLDTOOL_CHANGE, wxID_ANY);
    wxevent.SetClientData((void *)info);
    wxPostEvent(panel, wxevent);
  }
}


void
ThresholdTool::get_histogram(std::vector <int> & histogram )
{
  Nrrd *nrrd_in = painter_->current_volume_->nrrd_handle_->nrrd_;
  Nrrd *weight = 0;
  Nrrd *nrrd_histogram = nrrdNew();

  NrrdRange *range = NULL;

  nrrdRangeNew(AIR_NAN, AIR_NAN);
  nrrdRangeSafeSet(range, nrrd_in, nrrdBlind8BitRangeState);

  int nbins = 200;
  const unsigned int type = string_to_nrrd_type("nrrdTypeInt");

  if (nrrdHisto(nrrd_histogram, nrrd_in, range, weight, nbins, type))
  {
    char *err = biffGetDone(NRRD);
    painter_->set_status( string("Threshold tool error creating nrrd histogram: ") + err);
    free(err);
  }
  else
  {
    int* nrrd_data = (int *) (nrrd_histogram->data);
    
    histogram.clear();
    
    for( unsigned int i=0; i<nbins; ++i)
      histogram.push_back( nrrd_data[i] );

    painter_->set_status( string("Successfully generated a histogram"));

  }

  nrrdNuke(nrrd_histogram);
}


ColorMappedNrrdTextureObj *
ThresholdTool::get_texture(SliceWindow &window)
{
  if (!painter_->current_volume_.get_rep())
  {
    painter_->set_status("No active volume to threshold.");
    return NULL;
  }

  if (painter_->current_volume_->label_)
  {
    painter_->set_status("Active volume is not data.  Cannot threshold labels.");
    return NULL;
  }

  NrrdVolumeHandle vol = painter_->current_volume_;

  if (vol.get_rep() != cached_volume_pointer_)
  {
    texture_cache_map_t::iterator loc = texture_cache_.begin();
    while (loc != texture_cache_.end())
    {
      delete (*loc).second;
      ++loc;
    }
    texture_cache_.clear();
    cached_volume_pointer_ = vol.get_rep();
  }
  
  SLIVR::Plane plane(SLIVR::Point(window.center_.x(),window.center_.y(),
    window.center_.z()), SLIVR::Vector(window.normal_.x(),
    window.normal_.y(),window.normal_.z()));
  VolumeSliceHandle slice = vol->get_volume_slice(plane);

  ColorMappedNrrdTextureObj *texture = NULL;
  texture_cache_map_t::iterator loc = texture_cache_.find(&window);
  if (loc == texture_cache_.end() ||
      (*loc).second->get_min() != slice->pos_)
  {
    ColorMapHandle cmap = vol->get_colormap();
    texture = new ColorMappedNrrdTextureObj(slice->nrrd_handle_, cmap, 1);
    texture->set_label_color(label_color_);
    texture->set_coords(slice->pos_, slice->xdir_, slice->ydir_);
    //texture->set_opacity(0.75);
    if (loc != texture_cache_.end())
    {
      delete (*loc).second;
    }
    texture_cache_[&window] = texture;
  }
  else
  {
    texture = (*loc).second;
  }
  return texture;
}


void
ThresholdTool::draw_gl(SliceWindow &window)
{
  painter_->volume_lock_.lock();

  ColorMappedNrrdTextureObj *texture = get_texture(window);
  if (texture)
  {
    texture->set_threshold(lower_threshold_, upper_threshold_);

    glDisable(GL_CULL_FACE);
    glEnable(GL_BLEND);
    glDisable(GL_LIGHTING);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); 
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glShadeModel(GL_FLAT);

    texture->draw_quad();
  }
  painter_->volume_lock_.unlock();

  SeedTool::draw_gl(window);
}
  

void
ThresholdTool::run_filter()
{
  painter_->volume_lock_.lock();

  // Save off the source.
  NrrdVolumeHandle source_volume = painter_->current_volume_;

  // Make a new label volume
  const string name = painter_->current_volume_->name_ + " Threshold";
  painter_->create_new_label(painter_->current_volume_, name);

  painter_->rebuild_layer_buttons();
  painter_->volume_lock_.unlock();

  VolumeOps::threshold(painter_->current_volume_->nrrd_handle_,
                       painter_->current_volume_->label_,
                       source_volume->nrrd_handle_,
                       lower_threshold_, upper_threshold_);

  UndoHandle undo =
    new UndoReplaceLayer(painter_, "Undo Threshold",
                         0, painter_->current_volume_, 0);
  painter_->push_undo(undo);

  // Redraw everything after completion.
  painter_->extract_all_window_slices();
  painter_->redraw_all();

  painter_->hide_tool_panel();
}


}
