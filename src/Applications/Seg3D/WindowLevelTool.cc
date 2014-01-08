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
//    File   : WindowLevelTool.cc
//    Author : McKay Davis
//    Date   : Tue Sep 26 18:44:34 2006

#include <Applications/Seg3D/Painter.h>

#include <Applications/Seg3D/WindowLevelTool.h>
#include <Applications/Seg3D/Seg3DwxGuiUtils.h>

#include <sci_gl.h>

namespace SCIRun {

WindowLevelTool::WindowLevelTool(Painter *painter) :
  BaseTool("WindowLevelTool"),
  window_(0.0),
  level_(0.0),
  painter_(painter)
{
  // TODO:  Pick top visible by default if none are selected.
  painter_->check_for_active_data_volume("WindowLevel Tool");

  wxPanel *panel = Painter::global_seg3dframe_pointer_->CurrentToolPanel();
  if (panel && painter_->current_volume_.get_rep() &&
      !painter_->current_volume_->label_)
  {
    painter_->current_volume_->reset_data_range();

    window_ = (painter_->current_volume_->data_max_ -
       painter_->current_volume_->data_min_);

    level_ = (painter_->current_volume_->data_max_ +
       painter_->current_volume_->data_min_) / 2.0;

    WindowLevelToolChangeStruct *info = new WindowLevelToolChangeStruct;
    info->minval = painter_->current_volume_->data_min_;
    info->maxval = painter_->current_volume_->data_max_;
    info->window = window_;
    info->level  = level_;

    // This code appears to cause Windows builds to hang, cause unknown at the moment.
#ifndef _WIN32
    Painter::global_seg3dframe_pointer_->
      SetStatusText(std2wx("Creating a histogram for windowing and leveling"));
    wxBusyCursor(); // Busy cursor until this leaves scope.
    get_histogram(info->histogram);

    Painter::global_seg3dframe_pointer_->SetStatusText(std2wx(""));
#endif

    wxCommandEvent wxevent(wxEVT_WINDOWLEVELTOOL_CHANGE, wxID_ANY);
    wxevent.SetClientData((void *)info);
    wxPostEvent(panel, wxevent);
  }

  painter_->redraw_all();  
}


WindowLevelTool::~WindowLevelTool()
{
}


BaseTool::propagation_state_e 
WindowLevelTool::process_event(event_handle_t event)
{
  UpdateWindowLevelToolEvent *update =
    dynamic_cast<UpdateWindowLevelToolEvent *>(event.get_rep());
  if (update)
  {
    window_ = update->get_window();
    level_ = update->get_level();

    painter_->current_volume_->clut_min_ = level_ - window_ / 2.0;
    painter_->current_volume_->clut_max_ = level_ + window_ / 2.0;

    painter_->current_volume_->set_slices_dirty();
    painter_->redraw_all();

    painter_->redraw_all();
    return CONTINUE_E;
  }

  return CONTINUE_E;
}


void
WindowLevelTool::get_histogram(std::vector <int> & histogram )
{
  NrrdVolumeHandle source_volume = painter_->current_volume_;

  Nrrd *nrrd_in = source_volume->nrrd_handle_->nrrd_;
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
    painter_->set_status( string("WindowLevel tool error creating nrrd histogram: ") + err);
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

}
