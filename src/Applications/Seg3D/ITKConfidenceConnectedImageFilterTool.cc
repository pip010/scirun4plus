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
//    File   : ITKConfidenceConnectedImageFilterTool.cc
//    Author : McKay Davis
//    Date   : Tue Sep 26 18:44:34 2006

#include <Applications/Seg3D/Painter.h>

#include <Applications/Seg3D/ITKConfidenceConnectedImageFilterTool.h>
#include <sci_gl.h>


namespace SCIRun {

SeedTool::seeds_t ITKConfidenceConnectedImageFilterTool::seed_cache_;


ITKConfidenceConnectedImageFilterTool::
ITKConfidenceConnectedImageFilterTool(Painter *painter) :
  SeedTool("ITKConfidenceConnectedImageFilterTool", painter),
  prefix_("ITKConfidenceConnectedImageFilterTool::"),
  numberOfIterations_(painter_->get_vars(), prefix_+"numberOfIterations"),
  multiplier_(painter_->get_vars(), prefix_+"multiplier"),
  replaceValue_(painter_->get_vars(), prefix_+"replaceValue"),
  initialNeighborhoodRadius_(painter_->get_vars(), 
                             prefix_+"initialNeighborhoodRadius")
{
  seeds_ = seed_cache_;
  painter_->redraw_all();
}


ITKConfidenceConnectedImageFilterTool::~ITKConfidenceConnectedImageFilterTool()
{
  seed_cache_ = seeds_;
}


void
ITKConfidenceConnectedImageFilterTool::run_filter()
{
  if (!painter_->check_for_active_data_volume("Confidence connected filter"))
  {
    return;
  }

  if (seeds_.empty())
  {
    painter_->set_status("Confidence Connected Filter requires one or more seed points.");
  }

  painter_->volume_lock_.lock();

  NrrdVolumeHandle source_volume = painter_->current_volume_;

  vector<vector<int> > iseeds;
  convert_seeds_to_indices(iseeds, source_volume);

  if (iseeds.empty())
  {
    painter_->set_status("All of the seed points were outside the volume.");
    painter_->volume_lock_.unlock();
    return;
  }

  // Make a new label volume.
  NrrdDataHandle nrrdh =
    VolumeOps::create_clear_nrrd(painter_->current_volume_->nrrd_handle_,
                                 LabelNrrdType);
  const string name =
    painter_->unique_layer_name(painter_->current_volume_->name_ + " Confidence Connected");
  painter_->volumes_.push_back(new NrrdVolume(painter_, name, nrrdh, 1));
  painter_->current_volume_ = painter_->volumes_.back();
  painter_->rebuild_layer_buttons();
  painter_->volume_lock_.unlock();

  string status = "Confidence Connected Filter running with " +
    to_string(iseeds.size()) + " seed point";
  if (iseeds.size() > 1) status += "s";
  status += ".";
  painter_->set_status(status);

  typedef itk::ConfidenceConnectedImageFilter
    < ITKImageFloat3D, ITKImageLabel3D > FilterType;

  for (size_t i = 0; i < iseeds.size(); ++i) {
    FilterType::IndexType seed_point;
    for(unsigned int j = 0; j < seed_point.GetIndexDimension(); j++) {
      seed_point[j] = iseeds[i][j+1];
    }
    filter_->AddSeed(seed_point);
  }

  filter_->SetNumberOfIterations(numberOfIterations_);
  filter_->SetMultiplier(multiplier_);
  filter_->SetReplaceValue(replaceValue_);
  filter_->SetInitialNeighborhoodRadius(initialNeighborhoodRadius_);

  filter_.set_volume(painter_->current_volume_);

  painter_->start_progress();
  filter_.start(source_volume->nrrd_handle_);
  painter_->finish_progress();

  UndoHandle undo =
    new UndoReplaceLayer(painter_, "Undo Confidence Connected",
                         0, painter_->current_volume_, 0);
  painter_->push_undo(undo);

  painter_->extract_all_window_slices();
  painter_->redraw_all();

  painter_->hide_tool_panel();
}


}

