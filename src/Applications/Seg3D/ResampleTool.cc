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
//    File   : ResampleTool.cc
//    Author : McKay Davis
//    Date   : Sun Oct  1 23:22:04 2006



#include <Applications/Seg3D/Painter.h>
#include <Applications/Seg3D/ResampleTool.h>
#include <Core/Events/EventManager.h>
#include <sci_gl.h>

#include <Applications/Seg3D/GuiCode/resampletool.h>
#include <Applications/Seg3D/GuiCode/seg3devents.h>

#ifdef _WIN32
#  undef SCISHARE
#  define SCISHARE __declspec(dllimport)
#else
#  define SCISHARE
#endif

namespace SCIRun {


ResampleTool::ResampleTool(Painter *painter) : 
  BaseTool("Resample"),
  PointerTool("Resample"),
  painter_(painter)
{
}


ResampleTool::~ResampleTool()
{
}


BaseTool::propagation_state_e 
ResampleTool::process_event(event_handle_t event)
{
  if (dynamic_cast<FinishEvent *>(event.get_rep())) {
    return finish();
  }

  return CONTINUE_E;
}


BaseTool::propagation_state_e
ResampleTool::finish()
{
  if (!painter_->current_volume_.get_rep() ||
      painter_->current_volume_->label_)
  {
    painter_->set_status("Resample requires an active data volume.");
    return STOP_E;
  }

  NrrdResampleInfo *info = nrrdResampleInfoNew();

  NrrdKernel *kern = 0;
  double p[NRRD_KERNEL_PARMS_NUM];
  memset(p, 0, NRRD_KERNEL_PARMS_NUM * sizeof(double));
  p[0] = 1.0;

  const string filter = painter_->get_vars()->get_string("Painter::Resample::filter");

  if (filter == "box") {
    kern = nrrdKernelBox;
  } else if (filter == "tent") {
    kern = nrrdKernelTent;
  } else if (filter == "cubicCR") { 
    kern = nrrdKernelBCCubic; 
    p[1] = 0.0;
    p[2] = 0.5;
  } else if (filter == "cubicBS") { 
    kern = nrrdKernelBCCubic; 
    p[1] = 1.0;
    p[2] = 0.0;
  } else if (filter == "gaussian") { 
    kern = nrrdKernelGaussian; 
  } else { // default is quartic
    kern = nrrdKernelAQuartic;
    p[1] = 0.0834; // most accurate as per Teem documenation
  }

  Nrrd *nin = painter_->current_volume_->nrrd_handle_->nrrd_;

  int samples[4];
  samples[0] = 1;
  samples[1] = atoi( painter_->get_vars()->get_string("Painter::Resample::x").c_str() );
  samples[2] = atoi( painter_->get_vars()->get_string("Painter::Resample::y").c_str() );
  samples[3] = atoi( painter_->get_vars()->get_string("Painter::Resample::z").c_str() );

  for (int i = 0; i < 4; i++)
  {
    if (samples[i] < 1 || samples[i] > 4096)
    {
      samples[i] = nin->axis[i].size;
    }
  }

  for (int a = 0; a < 4; a++)
  {
    if (a == 0)
      info->kernel[a] = 0;
    else 
      info->kernel[a] = kern;

    info->samples[a] = (size_t)samples[a];

    memcpy(info->parm[a], p, NRRD_KERNEL_PARMS_NUM * sizeof(double));

    if (info->kernel[a] && 
    	(!(airExists(nin->axis[a].min) && airExists(nin->axis[a].max))))
    {
      nrrdAxisInfoMinMaxSet(nin, a, nin->axis[a].center ? 
                            nin->axis[a].center : nrrdDefaultCenter);
    }

    info->min[a] = nin->axis[a].min;
    info->max[a] = nin->axis[a].max;
  }    
  info->boundary = nrrdBoundaryBleed;
  info->type = nin->type;
  info->renormalize = AIR_TRUE;

  NrrdDataHandle nrrd_handle = new NrrdData();
  Nrrd *nout = nrrd_handle->nrrd_;
  if (nrrdSpatialResample(nout, nin, info)) {
    char *err = biffGetDone(NRRD);
    string errstr(err);
    free(err);
    throw "Trouble resampling: " + errstr;
  }
  nrrdResampleInfoNix(info); 
  
  // TODO: The spacings are wrong after resample.  It appears to
  // assume node centered data?
  //nout->axis[1].spacing *= (nout->axis[1].size - 1.0) / nout->axis[1].size;
  //nout->axis[2].spacing *= (nout->axis[2].size - 1.0) / nout->axis[2].size;
  //nout->axis[3].spacing *= (nout->axis[3].size - 1.0) / nout->axis[3].size;

  // TODO:  Fix the spacial transform for the new nrrd here.
  const string newname =
    painter_->unique_layer_name(painter_->current_volume_->name_ +
                                " Resampled");
  painter_->make_layer(newname, nrrd_handle, 0);

  UndoHandle undo =
    new UndoReplaceLayer(painter_, "Undo Resample",
                         0, painter_->current_volume_, 0);
  painter_->push_undo(undo);

  painter_->extract_all_window_slices();
  painter_->redraw_all();

  // Autoview to the cropped volume.
  event_handle_t empty;
  painter_->Autoview(empty);

  return CONTINUE_E;
}


}
