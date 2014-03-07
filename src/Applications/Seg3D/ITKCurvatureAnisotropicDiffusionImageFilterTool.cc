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
//    File   : ITKCurvatureAnisotropicDiffusionImageFilterTool.cc
//    Author : McKay Davis
//    Date   : Tue Sep 26 18:44:34 2006

#include <Applications/Seg3D/Painter.h>

#include <Applications/Seg3D/ITKCurvatureAnisotropicDiffusionImageFilterTool.h>
#include <sci_gl.h>


namespace SCIRun {

ITKCurvatureAnisotropicDiffusionImageFilterTool::
ITKCurvatureAnisotropicDiffusionImageFilterTool(Painter *painter) :
  BaseTool("ITK Curvature Anisotropic Diffusion \nImage Filter"),
  PointerTool("ITK Curvature Anisotropic Diffusion \nImage Filter"),
  painter_(painter),
  prefix_("ITKCurvatureAnisotropicDiffusionImageFilterTool::"),
  filter_()
{
}



void
ITKCurvatureAnisotropicDiffusionImageFilterTool::run_filter()
{
  string name = "ITKCurvatureAnisotropic";

  name = painter_->unique_layer_name(name + " Diffused");
  NrrdDataHandle nrrdh = painter_->current_volume_->nrrd_handle_;
  nrrdh.detach();
  NrrdVolume *vol = new NrrdVolume(painter_, name, nrrdh);
  painter_->volume_lock_.lock();
  painter_->volumes_.push_back(vol);
  painter_->current_volume_ = vol;
  painter_->rebuild_layer_buttons();
  painter_->volume_lock_.unlock();

  const int iterations =
    painter_->get_vars()->get_int(prefix_+"numberOfIterations");
  const double timestep = 
    painter_->get_vars()->get_double(prefix_+"timeStep");
  const double conductance =
    (painter_->get_vars()->get_double(prefix_+"conductanceParameter"));

  filter_->SetNumberOfIterations(iterations);
  filter_->SetTimeStep(timestep);
  filter_->SetConductanceParameter(conductance);

  filter_.set_volume(painter_->current_volume_);

  painter_->set_status("Running curvature anisotropic diffusion filter.");
  painter_->start_progress();
  filter_.start(nrrdh);
  painter_->finish_progress();
  
  UndoHandle undo =
    new UndoReplaceLayer(painter_, "Undo Anisotropic Diffusion",
                         0, painter_->current_volume_, 0);
  painter_->push_undo(undo);

  painter_->current_volume_->set_dirty();
  painter_->extract_all_window_slices();
  painter_->redraw_all();

  painter_->hide_tool_panel();
}




BaseTool::propagation_state_e 
ITKCurvatureAnisotropicDiffusionImageFilterTool::process_event(event_handle_t event)
{
  if (dynamic_cast<FinishEvent *>(event.get_rep()))
  {
    if (!painter_->check_for_active_data_volume("Curvature anisotropic diffusion filter"))
    {
      return STOP_E;
    }

    run_filter();

    if (filter_.stopped())
    {
      return CONTINUE_E;
    }
  }

  return CONTINUE_E;
}


}

