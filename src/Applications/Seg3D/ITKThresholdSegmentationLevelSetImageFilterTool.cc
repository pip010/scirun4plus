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
//    File   : ITKThresholdSegmentationLevelSetImageFilterTool.cc
//    Author : McKay Davis
//    Date   : Sat Oct 14 14:52:24 2006

#include <Applications/Seg3D/ITKThresholdSegmentationLevelSetImageFilterTool.h>

#include <Applications/Seg3D/Seg3DFrame.h>
#include <Applications/Seg3D/BrushTool.h>
#include <Applications/Seg3D/Painter.h>
#include <Applications/Seg3D/VolumeOps.h>
#include <Core/Util/Assert.h>


namespace SCIRun {

#if 0
class TSLSFilterRunnable : public Runnable
{
public:
  TSLSFilterRunnable(ITKThresholdSegmentationLevelSetImageFilterTool *tool)
    : tool_(tool) {}
    
  virtual void run() { tool_->run_filter(); }

  ITKThresholdSegmentationLevelSetImageFilterTool *tool_;
};
#endif

ITKThresholdSegmentationLevelSetImageFilterTool::
ITKThresholdSegmentationLevelSetImageFilterTool(Painter *painter) :
  BaseTool("ITKThresholdSegmentationLevelSetImageFilterTool::"),
  Runnable(false),
  painter_(painter),
  seed_volume_(0),
  source_volume_(0),
  threshold_volume_(0),
  LowerThreshold_(painter->get_vars(), get_name()+"LowerThreshold"),  
  UpperThreshold_(painter->get_vars(), get_name()+"UpperThreshold"),
  filter_(),
  started_(false)
{
}


void
ITKThresholdSegmentationLevelSetImageFilterTool::run()
{
  run_filter();
}


BaseTool::propagation_state_e 
ITKThresholdSegmentationLevelSetImageFilterTool::process_event
(event_handle_t event)
{
  if (dynamic_cast<SetLayerEvent *>(event.get_rep()))
  {
    if (!(painter_->current_volume_.get_rep() &&
          painter_->current_volume_->label_))
    {
      painter_->set_status("No active label volume to set as seed volume.");
      return STOP_E;
    }
    seed_volume_ = painter_->current_volume_;
    painter_->redraw_all();
    return STOP_E;
  }

  if (dynamic_cast<SetDataLayerEvent *>(event.get_rep())) 
  {
    if (!(painter_->current_volume_.get_rep() &&
          !painter_->current_volume_->label_))
    {
      painter_->set_status("No active data volume.");
      return STOP_E;
    }

    source_volume_ = painter_->current_volume_;
    painter_->set_status("Using " + source_volume_->name_ +
                         " as data volume.");
  }

  if (dynamic_cast<SetThresholdLayerEvent *>(event.get_rep())) 
  {
    if (!(painter_->current_volume_.get_rep() &&
          painter_->current_volume_->label_))
    {
      painter_->set_status("No active label volume.");
      return STOP_E;
    }

    threshold_volume_ = painter_->current_volume_;
    painter_->set_status("Using " + threshold_volume_->name_ +
                         " as threshold volume.");
  }


  if (dynamic_cast<FinishEvent *>(event.get_rep()))
  {
    set_vars();
    Thread *thread = new Thread(this, "TSLSFilter thread");
    thread->detach();
    return STOP_E;
  }

  return CONTINUE_E;
}


#define SetFilterVarMacro(name, type) \
  filter_->Set##name(painter_->get_vars()->get_##type(toolname+#name));

void
ITKThresholdSegmentationLevelSetImageFilterTool::set_vars()
{
  const string toolname = get_name();
  SetFilterVarMacro(CurvatureScaling, double);
  SetFilterVarMacro(PropagationScaling, double);
  SetFilterVarMacro(AdvectionScaling, double);
  SetFilterVarMacro(EdgeWeight, double);
  SetFilterVarMacro(NumberOfIterations, int);
  SetFilterVarMacro(MaximumRMSError, double);
  SetFilterVarMacro(IsoSurfaceValue, double);
  SetFilterVarMacro(SmoothingIterations,int);
  SetFilterVarMacro(SmoothingTimeStep, double);
  SetFilterVarMacro(SmoothingConductance, double);

  if (painter_->get_vars()->get_bool(toolname+"ReverseExpansionDirection"))
  {
    filter_->SetReverseExpansionDirection(true);
  }
  else
  {
    filter_->SetReverseExpansionDirection(false);
  }    
}
  

void
ITKThresholdSegmentationLevelSetImageFilterTool::run_filter()
{
#if 0
  if (started_) {
    set_vars();
    filter_->ManualReinitializationOn();
    filter_->Modified();
    filter_.update();
    
    return;
  }
  started_ = true;
#endif

  painter_->volume_lock_.lock();

  if (source_volume_.get_rep() == 0)
  {
    if (painter_->current_volume_.get_rep() &&
        !painter_->current_volume_->label_)
    {
      source_volume_ = painter_->current_volume_;
    }
    else if (!(painter_->current_volume_.get_rep() &&
               !painter_->current_volume_->label_))
    {
      int count = 0;
      NrrdVolumeHandle dvol;
      for (size_t i = 0; i < painter_->volumes_.size(); i++)
      {
        if (painter_->volumes_[i]->label_ == 0)
        {
          count++;
          dvol = painter_->volumes_[i];
        }
      }
      if (count == 1)
      {
        source_volume_ = dvol;
        painter_->set_status("Using "
                             + painter_->current_volume_->name_ +
                             " as source data volume.");
      }
    }
  }

  if (seed_volume_.get_rep() == 0)
  {
    if (painter_->current_volume_.get_rep() &&
        painter_->current_volume_->label_)
    {
      seed_volume_ = painter_->current_volume_;
    }
  }

  if (source_volume_.get_rep() == 0)
  {
    painter_->set_status("No source data volume set.");
    painter_->volume_lock_.unlock();
    return;
  }

  if (seed_volume_.get_rep() == 0)
  {
    painter_->set_status("No label seed volume set.");
    painter_->volume_lock_.unlock();
    return;
  }

  if (threshold_volume_.get_rep() == 0)
  {
    painter_->set_status("No threshold label volume set.");
    painter_->volume_lock_.unlock();
    return;
  }

  painter_->current_volume_ = seed_volume_;
  NrrdVolumeHandle new_layer = painter_->copy_current_layer(" Level Set");

  painter_->volume_lock_.unlock();

  NrrdStats stats =
    VolumeOps::nrrd_statistics(source_volume_->nrrd_handle_,
			       threshold_volume_->nrrd_handle_,
			       1 << threshold_volume_->bit());

  LowerThreshold_ = stats.min();
  UpperThreshold_ = stats.max();
  painter_->set_status("Threshold min: " + to_string(LowerThreshold_) + 
                       " Threshold max: " + to_string(UpperThreshold_));
  filter_->SetLowerThreshold(LowerThreshold_);
  filter_->SetUpperThreshold(UpperThreshold_);    

  ITKDatatypeHandle fimageh = source_volume_->get_itk_image();
  FeatureImg *fimage = dynamic_cast<FeatureImg *>(fimageh->data_.GetPointer());
  filter_->SetFeatureImage(fimage);

  new_layer->change_type_from_bit_to_float(1000.0);
  filter_.set_volume(new_layer);
  filter_->Modified();
  painter_->start_progress();
  filter_.start();    
  painter_->finish_progress();

  painter_->volume_lock_.lock();
  new_layer->change_type_from_float_to_bit();
  painter_->volume_lock_.unlock();

  UndoHandle undo =
    new UndoReplaceLayer(painter_, "Undo Threshold Segmentation Level Set",
                         0, new_layer, 0);
  painter_->push_undo(undo);

  painter_->extract_all_window_slices();
  painter_->redraw_all();

  //painter_->hide_tool_panel();
}


}
