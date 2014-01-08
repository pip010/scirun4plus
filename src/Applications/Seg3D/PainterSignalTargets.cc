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
//    File   : PainterSignalTargets.cc
//    Author : McKay Davis
//    Date   : Mon Jul  3 21:47:22 2006


#include <Applications/Seg3D/Painter.h>
#include <Applications/Seg3D/CropTool.h>
#include <Applications/Seg3D/ArithmeticTool.h>
#include <Applications/Seg3D/BrushTool.h>
#include <Applications/Seg3D/ResampleTool.h>
#include <Applications/Seg3D/InhomogeneityCorrectionFilter.h>
#include <Applications/Seg3D/ITKConfidenceConnectedImageFilterTool.h>
#include <Applications/Seg3D/ITKCurvatureAnisotropicDiffusionImageFilterTool.h>
#include <Applications/Seg3D/ITKNeighborhoodConnectedImageFilterTool.h>
#include <Applications/Seg3D/ITKThresholdSegmentationLevelSetImageFilterTool.h>
#include <Applications/Seg3D/ITKDiscreteGaussianImageFilterTool.h>
#include <Applications/Seg3D/ITKCannyEdgeDetectionImageFilterTool.h>
#include <Applications/Seg3D/ITKOtsuImageFilterTool.h>
#include <Applications/Seg3D/ITKGradientMagnitudeImageFilterTool.h>
#include <Applications/Seg3D/ITKHoleFillImageFilterTool.h>
#include <Applications/Seg3D/ITKBinaryDilateErodeImageFilterTool.h>
#include <Applications/Seg3D/BrushFloodFill.h>
#include <Applications/Seg3D/FloodfillTool.h>
#include <Applications/Seg3D/FloodFillCopyTool.h>
#include <Applications/Seg3D/LineTool.h>
#include <Applications/Seg3D/PolylineTool.h>
#include <Applications/Seg3D/ThresholdTool.h>
#include <Applications/Seg3D/WindowLevelTool.h>
#include <Applications/Seg3D/LineTool.h>
#include <Applications/Seg3D/SessionReader.h>
#include <Applications/Seg3D/SessionWriter.h>
#include <Applications/Seg3D/VolumeOps.h>
#include <Applications/Seg3D/VolumeFilter.h>

#include <Applications/Seg3D/ITKSpeedToPathGradientDescentImageFilterTool.h>
#include <Applications/Seg3D/ITKSpeedToPathRegularStepGradientDescentImageFilterTool.h>
#include <Applications/Seg3D/ITKSpeedToPathIterateNeighborhoodImageFilterTool.h>

#include <sci_comp_warn_fixes.h>
#include <iostream>
#include <fstream>
#include <sci_gl.h>
#include <Core/Datatypes/Field.h> 
#include <Core/Geom/GeomMaterial.h>
#include <Core/Geom/ColorMappedNrrdTextureObj.h>
#include <Core/Geom/GeomSwitch.h>
#include <Core/Geom/GeomCull.h>
#include <Core/Geom/GeomGroup.h>
#include <Core/Geom/IndexedGroup.h>
#include <Core/Geom/GeomTransform.h>
#include <Core/Geom/TexSquare.h>
#include <Core/Geom/OpenGLViewport.h>
#include <Core/Geom/FreeType.h>
#include <Core/Math/MiscMath.h>
#include <Core/Math/MinMax.h>
#include <Core/Thread/CleanupManager.h>
#include <Core/Thread/Runnable.h>
#include <Core/Thread/Mutex.h>
#include <Core/Util/Environment.h>
#include <Core/Volume/CM2Widget.h>
#include <Core/Geom/TextRenderer.h>
#include <Core/Geom/FontManager.h>
#include <Core/Skinner/Variables.h>
#include <Core/Events/EventManager.h>
#include <Core/Events/LoadVolumeEvent.h>
#include <Core/Events/BaseTool.h>
#include <Core/Util/FileUtils.h>
#include <Core/Events/SceneGraphEvent.h>
#include <slivr/ShaderProgramARB.h>
#include <Core/Volume/ColorMap2.h>
#include <Core/Volume/VolumeRenderer.h>
#include <Core/Skinner/GeomSkinnerVarSwitch.h>

#include <Core/Matlab/matlabfile.h>
#include <Core/Matlab/matlabarray.h>
#include <Core/Matlab/matlabconverter.h>

#include <teem/nrrd.h>

#include <string.h>
#include <stdio.h>

#ifndef _WIN32
#  include <sys/mman.h>
#else
#  include <Core/Util/Rand.h>
#  include <io.h>
#endif

#include <Applications/Seg3D/Seg3DwxGuiUtils.h>

namespace SCIRun {

// The maximum number of segmentations that the importer is allowed to create.
#define MAX_IMPORT_SEGMENTATIONS 64

BaseTool::propagation_state_e 
Painter::InitializeSignalCatcherTargets(event_handle_t &)
{
  REGISTER_CATCHER_TARGET(Painter::SliceWindow_Maker);
  REGISTER_CATCHER_TARGET(Painter::LayerButton_Maker);
    
  REGISTER_CATCHER_TARGET(Painter::StartBrushTool);
  REGISTER_CATCHER_TARGET(Painter::StartCropTool);
  REGISTER_CATCHER_TARGET(Painter::StartArithmeticTool);
  REGISTER_CATCHER_TARGET(Painter::VolumeInformation);
  REGISTER_CATCHER_TARGET(Painter::StartPolylineTool);
  REGISTER_CATCHER_TARGET(Painter::StartResampleTool);
  REGISTER_CATCHER_TARGET(Painter::StartThresholdTool);
  REGISTER_CATCHER_TARGET(Painter::StartWindowLevelTool);
  REGISTER_CATCHER_TARGET(Painter::StartMeasurementTool);

  REGISTER_CATCHER_TARGET(Painter::Autoview);
  REGISTER_CATCHER_TARGET(Painter::CopyLabel);
  REGISTER_CATCHER_TARGET(Painter::DeleteLayer);
  REGISTER_CATCHER_TARGET(Painter::DeleteLayer2);

  REGISTER_CATCHER_TARGET(Painter::MoveLayerUp);
  REGISTER_CATCHER_TARGET(Painter::MoveLayerDown);

  REGISTER_CATCHER_TARGET(Painter::FinishTool);
  REGISTER_CATCHER_TARGET(Painter::ClearTools);  
  REGISTER_CATCHER_TARGET(Painter::SetLayer);
  REGISTER_CATCHER_TARGET(Painter::SetDataLayer);
  REGISTER_CATCHER_TARGET(Painter::SetThresholdLayer);

  REGISTER_CATCHER_TARGET(Painter::ITKBinaryDilate);
  REGISTER_CATCHER_TARGET(Painter::ITKBinaryErode);  

  REGISTER_CATCHER_TARGET(Painter::FlipVolume);
  REGISTER_CATCHER_TARGET(Painter::InvertVolume);
  REGISTER_CATCHER_TARGET(Painter::HistoEqFilter);

  REGISTER_CATCHER_TARGET(Painter::ITKGradientMagnitude);
  REGISTER_CATCHER_TARGET(Painter::InhomogeneityCorrectionFilter);
  REGISTER_CATCHER_TARGET(Painter::ITKBinaryDilateErode);
  REGISTER_CATCHER_TARGET(Painter::ITKCurvatureAnisotropic);
  REGISTER_CATCHER_TARGET(Painter::ITKBinaryDilateErodeImageFilter);
  REGISTER_CATCHER_TARGET(Painter::start_ITKConfidenceConnectedImageFilterTool);
  REGISTER_CATCHER_TARGET(Painter::start_ITKNeighborhoodConnectedImageFilterTool);
  REGISTER_CATCHER_TARGET(Painter::start_ITKThresholdSegmentationLevelSetImageFilterTool);

  REGISTER_CATCHER_TARGET(Painter::start_ITKDiscreteGaussianImageFilterTool);
  REGISTER_CATCHER_TARGET(Painter::start_ITKCannyEdgeDetectionImageFilterTool);

  REGISTER_CATCHER_TARGET(Painter::start_ITKSpeedToPathGradientDescentImageFilterTool);
  REGISTER_CATCHER_TARGET(Painter::start_ITKSpeedToPathRegularStepGradientDescentImageFilterTool);
  REGISTER_CATCHER_TARGET(Painter::start_ITKSpeedToPathIterateNeighborhoodImageFilterTool);
  

  REGISTER_CATCHER_TARGET(Painter::ShowVolumeRendering);
  REGISTER_CATCHER_TARGET(Painter::ShowVolumeRendering2);
  REGISTER_CATCHER_TARGET(Painter::ShowOneIsosurface);
  REGISTER_CATCHER_TARGET(Painter::ShowAllIsosurfaces);

  REGISTER_CATCHER_TARGET(Painter::AbortFilterOn);

  REGISTER_CATCHER_TARGET(Painter::LoadVolume);
  REGISTER_CATCHER_TARGET(Painter::SaveVolume);
  REGISTER_CATCHER_TARGET(Painter::SaveVolumeWithHeader);
  REGISTER_CATCHER_TARGET(Painter::LoadSession);
  REGISTER_CATCHER_TARGET(Painter::SaveSession);
  REGISTER_CATCHER_TARGET(Painter::ImportSegmentation);
  REGISTER_CATCHER_TARGET(Painter::ExportSegmentation1);
  REGISTER_CATCHER_TARGET(Painter::ExportSegmentation2);
  REGISTER_CATCHER_TARGET(Painter::ExportSegmentationWithHeader1);
  REGISTER_CATCHER_TARGET(Painter::ExportSegmentationWithHeader2);
  REGISTER_CATCHER_TARGET(Painter::ExportStatistics);
  REGISTER_CATCHER_TARGET(Painter::ExportStatistics2);

  REGISTER_CATCHER_TARGET(Painter::Undo);

  REGISTER_CATCHER_TARGET(Painter::MedianFilterVolume);

  REGISTER_CATCHER_TARGET(Painter::CreateLabelVolume);  

  REGISTER_CATCHER_TARGET(Painter::RebuildLayers);

  REGISTER_CATCHER_TARGET(Painter::ITKOtsuFilter);
  REGISTER_CATCHER_TARGET(Painter::ITKHoleFillFilter);
  REGISTER_CATCHER_TARGET(Painter::FloodFillCopyFilter);
  REGISTER_CATCHER_TARGET(Painter::LabelInvertFilter);

  REGISTER_CATCHER_TARGET(Painter::FloodFillLayer);
  REGISTER_CATCHER_TARGET(Painter::FloodFillClear);

  REGISTER_CATCHER_TARGET(Painter::RasterPolyline);

  REGISTER_CATCHER_TARGET(Painter::ResetCLUT);
  REGISTER_CATCHER_TARGET(Painter::UpdateBrushRadius);
  REGISTER_CATCHER_TARGET(Painter::UpdateThresholdTool);
  REGISTER_CATCHER_TARGET(Painter::UpdateWindowLevelTool);

  REGISTER_CATCHER_TARGET(Painter::OpenLabelColorDialog);
  REGISTER_CATCHER_TARGET(Painter::SetLabelColor);

  REGISTER_CATCHER_TARGET(Painter::OpenOpacityDialog);
  REGISTER_CATCHER_TARGET(Painter::SetOpacityValue);

  REGISTER_CATCHER_TARGET(Painter::SetMaskLayer);
  REGISTER_CATCHER_TARGET(Painter::ClearMaskLayer);
  REGISTER_CATCHER_TARGET(Painter::MaskDataFilter);
  REGISTER_CATCHER_TARGET(Painter::MaskLabelFilter);
  REGISTER_CATCHER_TARGET(Painter::CombineMaskAnd);
  REGISTER_CATCHER_TARGET(Painter::CombineMaskRemove);
  REGISTER_CATCHER_TARGET(Painter::CombineMaskOr);
  REGISTER_CATCHER_TARGET(Painter::CombineMaskXor);

  REGISTER_CATCHER_TARGET(Painter::MoveScaleLayerUpdateGUI);
  REGISTER_CATCHER_TARGET(Painter::MoveScaleLayer);
  REGISTER_CATCHER_TARGET(Painter::MoveScaleAllLayers);

  REGISTER_CATCHER_TARGET(Painter::ShowVisibleItem);
  REGISTER_CATCHER_TARGET(Painter::RedrawAll);

  return STOP_E;
}

// This builds the slice windows. There are 4 view windows with 3
// orthogonal views in each. This accounts for 12 views - however
// there are 27 slicewidows created.

BaseTool::propagation_state_e 
Painter::SliceWindow_Maker(event_handle_t &event)
{
  Skinner::MakerSignal *maker_signal = 
    dynamic_cast<Skinner::MakerSignal *>(event.get_rep());
  ASSERT(maker_signal);

  SliceWindow *window = new SliceWindow(maker_signal->get_vars(), this);
  windows_.push_back(window);
  maker_signal->set_signal_thrower(window);
  maker_signal->set_signal_name(maker_signal->get_signal_name()+"_Done");

  return MODIFIED_E;
}


BaseTool::propagation_state_e 
Painter::LayerButton_Maker(event_handle_t &event)
{
  Skinner::MakerSignal *maker_signal = 
    dynamic_cast<Skinner::MakerSignal *>(event.get_rep());
  ASSERT(maker_signal);

  LayerButton *lb = new LayerButton(maker_signal->get_vars(), this);
  layer_buttons_.push_back(lb);
  maker_signal->set_signal_thrower(lb);
  maker_signal->set_signal_name(maker_signal->get_signal_name()+"_Done");
  return MODIFIED_E;
}


BaseTool::propagation_state_e 
Painter::StartBrushTool(event_handle_t &event)
{
  tm_.set_tool(new BrushTool(this), 25);
  tm_.set_tool(new BrushFloodFill(this), 24);
 return CONTINUE_E;
}
  

BaseTool::propagation_state_e 
Painter::StartCropTool(event_handle_t &event)
{
  tm_.set_tool(new CropTool(this), 25);
  return CONTINUE_E;
}


BaseTool::propagation_state_e 
Painter::StartArithmeticTool(event_handle_t &event)
{
  tm_.set_tool(new ArithmeticTool(this), 25);
  return CONTINUE_E;
}


BaseTool::propagation_state_e
Painter::VolumeInformation(event_handle_t &event)
{
  if (!current_volume_.get_rep())
  {
    set_status("No volume currently selected.");
    return STOP_E;
  }
  
  NrrdVolumeHandle &vol = current_volume_;
  Nrrd *nrrd = vol->nrrd_handle_->nrrd_;

  VolumeInfoPanelStruct *info = new VolumeInfoPanelStruct;
  info->name = vol->name_;
  info->origin_x = nrrd->axis[1].min;
  info->origin_y = nrrd->axis[2].min;
  info->origin_z = nrrd->axis[3].min;
  info->spacing_x = nrrd->axis[1].spacing;
  info->spacing_y = nrrd->axis[2].spacing;
  info->spacing_z = nrrd->axis[3].spacing;
  info->size_x = nrrd->axis[1].size;
  info->size_y = nrrd->axis[2].size;
  info->size_z = nrrd->axis[3].size;
  info->volume =
    info->size_x * info->spacing_x *
    info->size_y * info->spacing_y *
    info->size_z * info->spacing_z;
  info->is_label = vol->label_ != 0;
  info->is_mask = info->is_label && mask_volume_.get_rep()
    && mask_volume_.get_rep() != current_volume_.get_rep();

  switch( nrrd->space )
  {
  case nrrdSpaceUnknown:
    info->space = string("Unknown");
      break;
  case nrrdSpaceRightAnteriorSuperior:
    info->space = string("RightAnteriorSuperior");
      break;
  case nrrdSpaceLeftAnteriorSuperior:
    info->space = string("LeftAnteriorSuperior");
      break;
  case nrrdSpaceLeftPosteriorSuperior:
    info->space = string("LeftPosteriorSuperior");
      break;
  case nrrdSpaceRightAnteriorSuperiorTime:
    info->space = string("RightAnteriorSuperiorTime");
      break;
  case nrrdSpaceLeftAnteriorSuperiorTime:
    info->space = string("LeftAnteriorSuperiorTime");
      break;
  case nrrdSpaceLeftPosteriorSuperiorTime:
    info->space = string("LeftPosteriorSuperiorTime");
      break;
  case nrrdSpaceScannerXYZ:
    info->space = string("ScannerXYZ");
      break;
  case nrrdSpaceScannerXYZTime:
    info->space = string("ScannerXYZTime");
      break;
  case nrrdSpace3DRightHanded:
    info->space = string("3DRightHanded");
      break;
  case nrrdSpace3DLeftHanded:
    info->space = string("3DLeftHanded");
      break;
  case nrrdSpace3DRightHandedTime:
    info->space = string("3DRightHandedTime");
      break;
  case nrrdSpace3DLeftHandedTime:
    info->space = string("3DLeftHandedTime");
      break;
  }

  if (info->is_label)
  {
    info->labelcount = VolumeOps::bit_count(vol->nrrd_handle_, vol->label_);
    info->labelvolume = info->labelcount *
      info->spacing_x * info->spacing_y * info->spacing_z;
    if (info->is_mask)
    {
      info->maskcount = VolumeOps::bit_count(mask_volume_->nrrd_handle_,
                                             mask_volume_->label_);
      info->maskvolume = info->maskcount *
        info->spacing_x * info->spacing_y * info->spacing_z;
    }
  }
  else
  {
    info->datamin = vol->data_min_;
    info->datamax = vol->data_max_;
  }
  wxCommandEvent wxevent(wxEVT_VOLUME_INFO_PANEL, wxID_ANY);
  wxevent.SetClientData((void *)info);
  wxPostEvent(global_seg3dframe_pointer_, wxevent);

  return CONTINUE_E;
}


BaseTool::propagation_state_e
Painter::StartPolylineTool(event_handle_t &event)
{
  tm_.set_tool(new PolylineTool("PolylineTool", this), 25);
  return CONTINUE_E;
}


BaseTool::propagation_state_e 
Painter::StartResampleTool(event_handle_t &event)
{
  tm_.set_tool(new ResampleTool(this), 25);
  return CONTINUE_E;
}


BaseTool::propagation_state_e
Painter::StartThresholdTool(event_handle_t &event)
{
  tm_.set_tool(new ThresholdTool(this), 25);
  return CONTINUE_E;
}


BaseTool::propagation_state_e
Painter::StartWindowLevelTool(event_handle_t &event)
{
  tm_.set_tool(new WindowLevelTool(this), 25);
  return CONTINUE_E;
}


BaseTool::propagation_state_e
Painter::StartMeasurementTool(event_handle_t &event)
{
  tm_.set_tool(new LineTool(this), 25);
  return CONTINUE_E;
}


BaseTool::propagation_state_e 
Painter::Autoview(event_handle_t &)
{
  if (current_volume_.get_rep())
  {
    for (size_t i = 0; i < windows_.size(); i++)
    {
      windows_[i]->autoview(current_volume_);
    }
  }
  EventManager::add_event(new AutoviewEvent());
  EventManager::add_event(new WindowEvent(WindowEvent::REDRAW_E));
  return CONTINUE_E;
}


BaseTool::propagation_state_e 
Painter::CopyLabel(event_handle_t &)
{
  if (!check_for_active_label_volume("Copy label"))
  {
    return STOP_E;
  }
  
  volume_lock_.lock();

  NrrdDataHandle cnrrd = current_volume_->nrrd_handle_;
  const label_type clabel = current_volume_->label_;

  // Make a new child volume.
  const string name = current_volume_->name_ + " Copy";
  create_new_label(current_volume_, name);
  NrrdDataHandle dnrrd = current_volume_->nrrd_handle_;
  const label_type dlabel = current_volume_->label_;

  volume_lock_.unlock();

  // Copy what was in the current volume into the new volume.
  VolumeOps::bit_copy(dnrrd, dlabel, cnrrd, clabel);

  UndoHandle undo = new UndoReplaceLayer(this, "Undo Copy Layer",
					 0, current_volume_, 0);
  push_undo(undo);

  extract_all_window_slices();
  rebuild_layer_buttons();
  redraw_all();

  hide_tool_panel();

  return CONTINUE_E;
}


BaseTool::propagation_state_e 
Painter::DeleteLayer(event_handle_t &event)
{
  NrrdVolumeHandle layer = current_volume_;
  if (event.get_rep()) {
    Skinner::Signal *signal = dynamic_cast<Skinner::Signal *>(event.get_rep());
    ASSERT(signal);
    Skinner::Variables *vars = signal->get_vars();
    const string name = vars->get_string("LayerButton::name");
    layer = find_volume_by_name(name);
    current_volume_ = layer;
  }

  if (!layer.get_rep()) return STOP_E;

  // Send to wx side for confirmation dialog.
  wxCommandEvent wxevent(wxEVT_COMMAND_LAYER_DELETE_DIALOG, wxID_ANY);
  wxevent.SetString(std2wx("Delete layer " + layer->name_ + "?"));
  wxPostEvent(global_seg3dframe_pointer_, wxevent);

  return CONTINUE_E;
}


BaseTool::propagation_state_e 
Painter::DeleteLayer2(event_handle_t &event)
{
  NrrdVolumeHandle &layer = current_volume_;
  if (!layer.get_rep()) return STOP_E;

  // Clear the volume rendering.
  if (layer == current_vrender_target_)
  {
    // If it's been deferred then no need to update the scene graph.
    if (!current_vrender_target_deferred_)
    {
      NrrdDataHandle empty(0);
      LoadVolumeEvent *lve = new LoadVolumeEvent(empty, "");
      EventManager::add_event(lve);
    }
    current_vrender_target_ = 0;
    current_vrender_target_deferred_ = false;
  }
  
  volume_lock_.lock();
  try {

    // Clear the mask layer if we're deleting it.
    if (current_volume_ == mask_volume_)
    {
      mask_volume_ = 0;
    }

    // Make the new volume list with the deleted volume missing.
    // TODO:  Put the old layer in an undo list instead of deleting.
    int loc = -1;
    NrrdVolumes newvolumes;
    for (size_t i = 0; i < volumes_.size(); i++)
    {
      if (volumes_[i] != layer)
      {
        newvolumes.push_back(volumes_[i]);
      }
      else
      {
        loc = i;
      }
    }
    
    if (loc != -1)
    {
      UndoHandle undo = new UndoReplaceLayer(this, "Undo Layer Delete",
                                             volumes_[loc], 0, loc);
      push_undo(undo);
    }

    volumes_ = newvolumes;

    // Clear the 3D graphics associated with this layer.
    EventManager::add_event(new SceneGraphEvent(0, layer->name_));

    // Remove the label allocation internal bookkeeping for this layer.
    layer->unparent();

    layer = 0;

    // Set the current layer to be near the deleted layer.
    if (loc == -1 || volumes_.size() == 0)
    {
      current_volume_ = 0;
    }
    else
    {
      loc = Min((int)volumes_.size()-1, loc);
      current_volume_ = volumes_[loc];
    }
  }
  catch(...) {
    volume_lock_.unlock();
    throw;
  }
  volume_lock_.unlock();

  rebuild_layer_buttons();
  extract_all_window_slices();
  redraw_all();

  return CONTINUE_E;
}


BaseTool::propagation_state_e
Painter::MoveLayerUp(event_handle_t &event)
{
  if (!current_volume_.get_rep())
  {
    set_status("No active volume to move.");
    return STOP_E;
  }
  move_layer_up();
  return CONTINUE_E;
}


BaseTool::propagation_state_e
Painter::MoveLayerDown(event_handle_t &event)
{
  if (!current_volume_.get_rep())
  {
    set_status("No active volume to move.");
    return STOP_E;
  }
  move_layer_down();
  return CONTINUE_E;
}


BaseTool::propagation_state_e
Painter::ITKBinaryDilateErodeImageFilter(event_handle_t &event)
{
  tm_.set_tool(new ITKBinaryDilateErodeImageFilterTool(this), 25);
  return CONTINUE_E;
}


BaseTool::propagation_state_e 
Painter::ITKBinaryDilate(event_handle_t &event)
{
 tm_.propagate_event(new DilateEvent());
 return CONTINUE_E;
}


BaseTool::propagation_state_e 
Painter::ITKBinaryErode(event_handle_t &event)
{
  tm_.propagate_event(new ErodeEvent());
  return CONTINUE_E;
}


BaseTool::propagation_state_e 
Painter::ITKBinaryDilateErode(event_handle_t &event)
{
  tm_.propagate_event(new DilateErodeEvent());
  return CONTINUE_E;
}


BaseTool::propagation_state_e 
Painter::LoadVolume(const string &filename)
{
  int status;
  NrrdVolumeHandle volume = load_volume<float>(filename, status);
  if (status == 1)
  {
    set_status_dialog("Load of " + filename + " failed!",
                      "Possible reasons include bad file permissions or insufficient memory.");
    return STOP_E;
  }
  else if (status == 2)
  {
    // TODO: This is just a workaround for not just offering to load
    // the session instead anyway.  Session loading needs a confirm
    // dialog because it's destructive and these two things need to be
    // consistent.
    const char *msg = "This file has an unsupported format.";
    if (ends_with(string_tolower(filename), ".ses"))
    {
      msg = "This appears to be a session file and not a volume.";
    }
    set_status_dialog("Load of " + filename + " failed!", msg);
    return STOP_E;
  }
  else if (status == 3)
  {
    set_status_dialog("Load of " + filename + " failed!", "Not a valid file.");
    return STOP_E;
  }

  // If we loaded an image and it appears to be part of a series
  // go ahead and load the whole series in.
  // TODO:  Pop up UI and ask for series parameters
  //        (start, stop, increment, dicom ordering)
  if (!volume.get_rep() || volume->nrrd_handle_->nrrd_->axis[3].size == 1)
  {
    const vector<string> filenames = get_filename_series(filename);
    if (filenames.size() > 1)
    {
      NrrdVolumeHandle volume2 = load_image_series<float>(filenames);
      if (volume2.get_rep()) {
        volume = volume2;
      }
    }
  }
  
  if (!volume.get_rep())
  {
    set_status_dialog("Load of " + filename + " failed!",
                      "No reason given.");
    return STOP_E;
  }
  
  volume_lock_.lock();
  const bool first = volumes_.size() == 0;
  volumes_.push_back(volume);
  current_volume_ = volume;

  for (size_t i = 0; i < windows_.size(); ++ i)
  {
    windows_[i]->center_ = volume->center();
  }

  extract_all_window_slices();
  rebuild_layer_buttons();
  redraw_all();
  set_status("Successfully loaded file: " + filename);
  volume_lock_.unlock();

  // Set the volume rendering target.
  if (!current_vrender_target_.get_rep())
  {
    // If volume rendering is on display this volume now.
    if (get_vars()->get_bool("Painter::volume_visible"))
    {
      event_handle_t empty;
      ShowVolumeRendering(empty);
    }
    else
    {
      // Delay computing the volume rendering until it's turned on.
      current_vrender_target_ = current_volume_;
      current_vrender_target_deferred_ = true;
    }
  }

  if (first)
  {
    // Force the creation of slice SceneGraphEvents before we do the
    // autoview.  Otherwise they aren't created until the first redraw
    // which can happen after the autoview is completed resulting in a
    // failed autoview.
    current_volume_->get_geom_group();

    event_handle_t unused;
    Autoview(unused);
  }

  return CONTINUE_E;
}


BaseTool::propagation_state_e 
Painter::LoadVolume(event_handle_t &event)
{
  const string filename =
    get_vars()->get_string("Painter::LoadVolume::filename");
  return LoadVolume(filename);
}


#define NRRD_EXEC(__nrrd_command__) \
  if (__nrrd_command__) { \
    char *err = biffGetDone(NRRD); \
    cerr << string("Error on line #")+to_string(__LINE__) + \
	    string(" executing nrrd command: \n")+ \
            string(#__nrrd_command__)+string("\n")+ \
            string("Message: ")+string(err); \
    free(err); \
    return STOP_E; \
  }


BaseTool::propagation_state_e 
Painter::FlipVolume(event_handle_t &event)
{
  if (!(current_volume_.get_rep()))
  {
    set_status("No active volume to flip.");
    return STOP_E;
  }

  if (current_volume_->label_)
  {
    set_status("Flipping top level label and all children.");
  }

  Skinner::Var<int> axis(get_vars(), "Painter::flip::axis", 3);

  NrrdDataHandle nout = new NrrdData();
  NRRD_EXEC(nrrdFlip(nout->nrrd_, current_volume_->nrrd_handle_->nrrd_, axis));

  nrrdKeyValueCopy(nout->nrrd_, current_volume_->nrrd_handle_->nrrd_);

  current_volume_->set_nrrd(nout);
  
  string undostr;
  if (axis == 1) undostr = "Undo Sagital Flip";
  else if (axis == 2) undostr = "Undo Coronal Flip";
  else undostr = "Undo Axial Flip";
  UndoHandle undo = new UndoFlipAxis(this, undostr, current_volume_, axis);
  push_undo(undo);

  extract_all_window_slices();
  redraw_all();  
  
  if (current_volume_ == current_vrender_target_ &&
      !current_vrender_target_deferred_)
  {
    LoadVolumeEvent *lve =
      new LoadVolumeEvent(current_vrender_target_->nrrd_handle_, "", false);
    EventManager::add_event(lve);
  }

  return CONTINUE_E;
}


BaseTool::propagation_state_e 
Painter::InvertVolume(event_handle_t &event)
{
  if (!(current_volume_.get_rep()))
  {
    set_status("No active volume to invert.");
    return STOP_E;
  }

  if (current_volume_->label_)
  {
    set_status("Can not invert labeled volumes.");
    return STOP_E;
  }

  Nrrd *nin1 = current_volume_->nrrd_handle_->nrrd_;

  NrrdIter *in1 = nrrdIterNew();
  NrrdIter *in2 = nrrdIterNew();

  nrrdIterSetOwnNrrd(in1, nin1);
  nrrdIterSetValue(in2, current_volume_->data_max_ + current_volume_->data_min_);

  NrrdDataHandle nout_data_handle = new NrrdData();
  if (nrrdArithIterBinaryOp(nout_data_handle->nrrd_, nrrdBinaryOpSubtract, in2, in1)) 
  {
    char *err = biffGetDone(NRRD);
    set_status( string("Arithmetic tool error performing nrrd operation failed: ") + err);
    free(err);
    return STOP_E;
  }

  nrrdKeyValueCopy(nout_data_handle->nrrd_, current_volume_->nrrd_handle_->nrrd_);

  current_volume_->set_nrrd(nout_data_handle);
  
  UndoHandle undo = new UndoInvert(this, "Undo Invert", current_volume_);
  push_undo(undo);

  extract_all_window_slices();
  redraw_all();  
  
  if (current_volume_ == current_vrender_target_ &&
      !current_vrender_target_deferred_)
  {
    LoadVolumeEvent *lve =
      new LoadVolumeEvent(current_vrender_target_->nrrd_handle_, "", false);
    EventManager::add_event(lve);
  }

  return CONTINUE_E;
}


// Unused
BaseTool::propagation_state_e 
Painter::HistoEqFilter(event_handle_t &event)
{
  if (!check_for_active_data_volume("Histogram equalization"))
  {
    return STOP_E;
  }

  int bins = get_vars()->get_int("Painter::HistoEq::bins");
  double alpha = get_vars()->get_double("Painter::HistoEq::alpha");

  if (alpha < 0.0 || alpha > 1.0)
  {
    set_status("Histogram equalization alpha must be between 0.0 and 1.0.  Using 1.0.");
    alpha = 1.0;
  }

  NrrdVolumeHandle src = current_volume_;
  copy_current_layer(" Histogram Equalized");

  NRRD_EXEC(nrrdHistoEq(current_volume_->nrrd_handle_->nrrd_,
                        src->nrrd_handle_->nrrd_,
                        NULL, 3000, bins, alpha));

  current_volume_->reset_data_range();

  UndoHandle undo = new UndoReplaceLayer(this, "Undo Histogram Equalization",
					 0, current_volume_, 0);
  push_undo(undo);

  extract_all_window_slices();
  redraw_all();

  hide_tool_panel();

  return CONTINUE_E;
}


BaseTool::propagation_state_e 
Painter::ITKOtsuFilter(event_handle_t &event)
{
  tm_.set_tool( new ITKOtsuImageFilterTool(this), 25 );
  return CONTINUE_E;
}


BaseTool::propagation_state_e 
Painter::ITKHoleFillFilter(event_handle_t &event)
{
  tm_.set_tool( new ITKHoleFillImageFilterTool(this), 25 );
  return CONTINUE_E;
}


BaseTool::propagation_state_e 
Painter::FloodFillCopyFilter(event_handle_t &event)
{
  tm_.set_tool( new FloodFillCopyTool(this), 25 );
  return CONTINUE_E;
}


BaseTool::propagation_state_e 
Painter::ITKGradientMagnitude(event_handle_t &)
{
  tm_.set_tool( new ITKGradientMagnitudeImageFilterTool(this), 25 );
  return CONTINUE_E;
}


BaseTool::propagation_state_e 
Painter::InhomogeneityCorrectionFilter(event_handle_t &)
{
  tm_.set_tool( new InhomogeneityCorrectionFilterTool(this), 25 );
  return CONTINUE_E;
}


BaseTool::propagation_state_e 
Painter::ITKCurvatureAnisotropic(event_handle_t &event)
{
  tm_.set_tool(new ITKCurvatureAnisotropicDiffusionImageFilterTool(this), 25);
  return CONTINUE_E;
}


BaseTool::propagation_state_e 
Painter::start_ITKNeighborhoodConnectedImageFilterTool(event_handle_t &event)
{
  tm_.set_tool(new ITKNeighborhoodConnectedImageFilterTool(this), 25);
  return CONTINUE_E;
}


BaseTool::propagation_state_e 
Painter::FloodFillClear(event_handle_t &event)
{
  BrushFloodFill *tool =
    dynamic_cast<BrushFloodFill *>(tm_.get_tool(24).get_rep());
  
  if (tool) { tool->flood_fill_slice(true); }
  
  return CONTINUE_E;
}


BaseTool::propagation_state_e 
Painter::FloodFillLayer(event_handle_t &event)
{
  BrushFloodFill *tool =
    dynamic_cast<BrushFloodFill *>(tm_.get_tool(24).get_rep());
  
  if (tool) { tool->flood_fill_slice(false); }

  return CONTINUE_E;
}


BaseTool::propagation_state_e 
Painter::RasterPolyline(event_handle_t &event)
{
  PolylineTool *tool =
    dynamic_cast<PolylineTool *>(tm_.get_tool(25).get_rep());
  
  if (tool)
  {
    const bool erase = get_vars()->get_bool("Painter::polylinetool_erase");
    tool->run_filter(erase);
  }

  return CONTINUE_E;
}



BaseTool::propagation_state_e
Painter::ResetCLUT(event_handle_t &event)
{
  // Find the topmost visible data layer.
  NrrdVolumeHandle volume = 0;
  for (int i = volumes_.size()-1; i >= 0; i--)
  {
    if (!volumes_[i]->label_ && volumes_[i]->visible())
    {
      volume = volumes_[i];
      break;
    }
  }

  if (volume.get_rep() == 0)
  {
    set_status("Reset brightness/contrast requires a visible data volume.");
    return STOP_E;
  }

  volume->reset_clut();
  redraw_all();

  return STOP_E;
}


BaseTool::propagation_state_e
Painter::UpdateBrushRadius(event_handle_t &event)
{
  const double radius = get_vars()->get_double("Painter::brush_radius");
  const bool force_erasing =
    get_vars()->get_bool("Painter::brush_force_erasing");
  event_handle_t sevent =
    new SCIRun::UpdateBrushRadiusEvent(radius, force_erasing);
  tm_.propagate_event(sevent);
  return CONTINUE_E;
}


BaseTool::propagation_state_e
Painter::UpdateThresholdTool(event_handle_t &event)
{
  const double lower = get_vars()->get_double("Painter::thresholdtool_lower");
  const double upper = get_vars()->get_double("Painter::thresholdtool_upper");
  event_handle_t sevent =  new SCIRun::UpdateThresholdToolEvent(lower, upper);
  tm_.propagate_event(sevent);
  return CONTINUE_E;
}


BaseTool::propagation_state_e
Painter::UpdateWindowLevelTool(event_handle_t &event)
{
  const double window = get_vars()->get_double("Painter::windowleveltool_window");
  const double level = get_vars()->get_double("Painter::windowleveltool_level");
  event_handle_t sevent =  new SCIRun::UpdateWindowLevelToolEvent(window, level);
  tm_.propagate_event(sevent);
  return CONTINUE_E;
}


BaseTool::propagation_state_e
Painter::OpenLabelColorDialog(event_handle_t &event)
{
  wxCommandEvent wxevent(wxEVT_COMMAND_COLOUR_PICKER, wxID_ANY);
  colour_picker_data_t *cpd = new colour_picker_data_t;
  Color c = current_volume_->get_label_color();
  cpd->r = c.r();
  cpd->g = c.g();
  cpd->b = c.b();
  wxevent.SetClientData(cpd);
  wxPostEvent(global_seg3dframe_pointer_, wxevent);

  return CONTINUE_E;
}


BaseTool::propagation_state_e
Painter::SetLabelColor(event_handle_t &event)
{
  const double r = get_vars()->get_double("Painter::SetLabelColor::r");
  const double g = get_vars()->get_double("Painter::SetLabelColor::g");
  const double b = get_vars()->get_double("Painter::SetLabelColor::b");

  Color c(r, g, b);
  current_volume_->set_label_color(c);

  rebuild_layer_buttons();
  extract_all_window_slices();
  redraw_all();

  return CONTINUE_E;
}


// Get the layer button event to open the opacity dialog.
BaseTool::propagation_state_e
Painter::OpenOpacityDialog(event_handle_t &event)
{
  wxCommandEvent wxevent(wxEVT_OPACITY_PANEL, wxID_ANY);

  // Pass the opacity of the current volume to the dialog.
  double * opacity = new double;
  *opacity = current_volume_->get_opacity();

  wxevent.SetClientData(opacity);
  wxPostEvent(global_seg3dframe_pointer_, wxevent);

  return CONTINUE_E;
}


// Get the set opacity event from the opacity dialog.
BaseTool::propagation_state_e
Painter::SetOpacityValue(event_handle_t &event)
{
  // Get the new opacity value set by the dialog and set the current
  // volume.
  const double opacity =
    atof( get_vars()->get_string("Painter::SetOpacityValue::value").c_str() );

  current_volume_->set_opacity(opacity);

  // Do this to update the opacity of the volume colormaps.
  if( current_volume_->label_ )
  {
    Color c = current_volume_->get_label_color();
    current_volume_->set_label_color(c);
  }

  rebuild_layer_buttons();
  extract_all_window_slices();
  redraw_all();

  return CONTINUE_E;
}


BaseTool::propagation_state_e
Painter::LabelInvertFilter(event_handle_t &)
{
  if (!check_for_active_label_volume("Label invert"))
  {
    return STOP_E;
  }

  VolumeOps::bit_invert(current_volume_->nrrd_handle_,
                        current_volume_->label_,
                        current_volume_->nrrd_handle_,
                        current_volume_->label_);

  UndoHandle undo = new UndoLabelInvertFilter(this, current_volume_);
  push_undo(undo);

  extract_all_window_slices();
  redraw_all();

  hide_tool_panel();

  return CONTINUE_E;
}


BaseTool::propagation_state_e
Painter::SetMaskLayer(event_handle_t &event)
{
  if (!check_for_active_label_volume("Set mask label"))
  {
    return STOP_E;
  }

  mask_volume_ = current_volume_;
  rebuild_layer_buttons();
  redraw_all();
  return CONTINUE_E;
}


BaseTool::propagation_state_e
Painter::ClearMaskLayer(event_handle_t &event)
{
  mask_volume_ = 0;
  rebuild_layer_buttons();
  redraw_all();
  return CONTINUE_E;
}


BaseTool::propagation_state_e
Painter::MaskDataFilter(event_handle_t &)
{
  if (!check_for_active_data_volume("Mask data"))
  {
    return STOP_E;
  }

  if (!check_for_valid_mask("Mask data", true))
  {
    return STOP_E;
  }

  volume_lock_.lock();

  NrrdDataHandle mnrrd = mask_volume_->nrrd_handle_;
  set_status("Using " + mask_volume_->name_ + " as mask.");

  NrrdDataHandle snrrd = current_volume_->nrrd_handle_;

  current_volume_ = copy_current_layer(" Masked");
  NrrdDataHandle dnrrd = current_volume_->nrrd_handle_;

  volume_lock_.unlock();

  const label_type mlabel = mask_volume_->label_;

  unsigned int maskValue = 
    atoi(get_vars()->get_string("MaskFilter::MaskValue").c_str());

  // Set the new value based on the GUI selection.
  float newval;

  // When using the new masked image minimum or maximum initally set
  // the masked value to the opposite extrema. I.e. if the minimum is
  // needed set the masked value to the maximum so that the new
  // minimum can be found without interfernce from any masked out
  // values.

  switch( maskValue )
  {
  case 0: // current miniumu
  case 3: // new maximum
    newval = current_volume_->data_min_;
    break;

  case 1: // new minimum
  case 4: // current maximum
    newval = current_volume_->data_max_;
    break;

  case 2: // zero
  default:
    newval = 0;
    break;
  }

  VolumeOps::mask_data(dnrrd, snrrd, mnrrd, mlabel, newval);

  current_volume_->reset_data_range();

  // For the new min or max from the mask get the new min max values
  // and remask.
  if( maskValue == 1 || maskValue == 3 )
  {
    NrrdRange range;
    nrrdRangeSet(&range, dnrrd->nrrd_, 0);

    if( maskValue == 1)
	newval = range.min;
    else if( maskValue == 3 )
	newval = range.max;
    
    VolumeOps::mask_data(dnrrd, snrrd, mnrrd, mlabel, newval);
  }

  current_volume_->reset_data_range();

  UndoHandle undo = new UndoReplaceLayer(this, "Undo Mask Data",
					 0, current_volume_, 0);
  push_undo(undo);

  extract_all_window_slices();
  redraw_all();

  hide_tool_panel();

  return CONTINUE_E;
}


BaseTool::propagation_state_e
Painter::MaskLabelFilter(event_handle_t &)
{
  if (!check_for_active_label_volume("Mask label"))
  {
    return STOP_E;
  }

  if (!(mask_volume_.get_rep()))
  {
    set_status("Mask label requires a mask to be selected.");
    return STOP_E;
  }

  volume_lock_.lock();

  NrrdDataHandle mnrrd = mask_volume_->nrrd_handle_;
  const label_type mlabel = mask_volume_->label_;
  set_status("Using " + mask_volume_->name_ + " as mask.");

  NrrdDataHandle cnrrd = current_volume_->nrrd_handle_;
  const label_type clabel = current_volume_->label_;

  // Make a new child volume.
  NrrdDataHandle dnrrd =
    VolumeOps::create_clear_nrrd(current_volume_->nrrd_handle_, LabelNrrdType);
  const label_type dlabel = 1;
  const string name = unique_layer_name(current_volume_->name_ + " Masked");
  volumes_.push_back(new NrrdVolume(this, name, dnrrd, dlabel));
  current_volume_ = volumes_.back();
  rebuild_layer_buttons();
  volume_lock_.unlock();

  // Copy what was in the current volume into the new volume.
  VolumeOps::bit_copy(dnrrd, dlabel, cnrrd, clabel);

  // Mask the new data.
  const size_t dsize = VolumeOps::nrrd_elem_count(dnrrd);
  const size_t msize = VolumeOps::nrrd_elem_count(mnrrd);
  ASSERT(dsize == msize);
  ASSERT(dnrrd->nrrd_->type == LabelNrrdType);
  ASSERT(mnrrd->nrrd_->type == LabelNrrdType);

  label_type *dst = (label_type *)dnrrd->nrrd_->data;
  label_type *mask = (label_type *)mnrrd->nrrd_->data;
  for (size_t i = 0; i < dsize; ++i, ++dst, ++mask) {
    if (!(*mask & mlabel)) { *dst &= ~dlabel; }
  }

  extract_all_window_slices();
  redraw_all();

  hide_tool_panel();

  return CONTINUE_E;
}


BaseTool::propagation_state_e
Painter::CombineMaskAnd(event_handle_t &)
{
  if (!check_for_active_label_volume("Combine mask with logical and"))
  {
    return STOP_E;
  }

  if (!check_for_valid_mask("Combine mask with logical and", true))
  {
    return STOP_E;
  }

  volume_lock_.lock();

  NrrdDataHandle mnrrd = mask_volume_->nrrd_handle_;
  const label_type mlabel = mask_volume_->label_;
  set_status("Using " + mask_volume_->name_ + " as mask.");

  NrrdDataHandle cnrrd = current_volume_->nrrd_handle_;
  const label_type clabel = current_volume_->label_;

  // Make a new child volume.
  const string name = current_volume_->name_ + " and " + mask_volume_->name_;
  current_volume_ = create_new_label(current_volume_, name);
  NrrdDataHandle dnrrd = current_volume_->nrrd_handle_;
  const label_type dlabel = current_volume_->label_;

  volume_lock_.unlock();

  // Copy what was in the current volume into the new volume.
  VolumeOps::bit_and(dnrrd, dlabel, cnrrd, clabel, mnrrd, mlabel);

  UndoHandle undo = new UndoReplaceLayer(this, "Undo Logical And",
					 0, current_volume_, 0);
  push_undo(undo);

  extract_all_window_slices();
  rebuild_layer_buttons();
  redraw_all();

  hide_tool_panel();

  return CONTINUE_E;
}


BaseTool::propagation_state_e
Painter::CombineMaskRemove(event_handle_t &)
{
  if (!check_for_active_label_volume("Combine mask with logical remove"))
  {
    return STOP_E;
  }

  if (!check_for_valid_mask("Combine mask with logical remove", true))
  {
    return STOP_E;
  }

  volume_lock_.lock();

  NrrdDataHandle mnrrd = mask_volume_->nrrd_handle_;
  const label_type mlabel = mask_volume_->label_;
  set_status("Using " + mask_volume_->name_ + " as mask.");

  NrrdDataHandle cnrrd = current_volume_->nrrd_handle_;
  const label_type clabel = current_volume_->label_;

  // Make a new child volume.
  const string name = current_volume_->name_ + " minus " + mask_volume_->name_;
  current_volume_ = create_new_label(current_volume_, name);
  NrrdDataHandle dnrrd = current_volume_->nrrd_handle_;
  const label_type dlabel = current_volume_->label_;

  volume_lock_.unlock();

  // Copy what was in the current volume into the new volume.
  VolumeOps::bit_andnot(dnrrd, dlabel, cnrrd, clabel, mnrrd, mlabel);

  UndoHandle undo = new UndoReplaceLayer(this, "Undo Logical Remove",
					 0, current_volume_, 0);
  push_undo(undo);

  extract_all_window_slices();
  rebuild_layer_buttons();
  redraw_all();

  hide_tool_panel();

  return CONTINUE_E;
}


BaseTool::propagation_state_e
Painter::CombineMaskOr(event_handle_t &)
{
  if (!check_for_active_label_volume("Combine mask with logical or"))
  {
    return STOP_E;
  }

  if (!check_for_valid_mask("Combine mask with logical or", true))
  {
    return STOP_E;
  }

  volume_lock_.lock();

  NrrdDataHandle mnrrd = mask_volume_->nrrd_handle_;
  const label_type mlabel = mask_volume_->label_;
  set_status("Using " + mask_volume_->name_ + " as mask.");

  NrrdDataHandle cnrrd = current_volume_->nrrd_handle_;
  const label_type clabel = current_volume_->label_;

  // Make a new child volume.
  const string name = current_volume_->name_ + " or " + mask_volume_->name_;
  current_volume_ = create_new_label(current_volume_, name);
  NrrdDataHandle dnrrd = current_volume_->nrrd_handle_;
  const label_type dlabel = current_volume_->label_;

  volume_lock_.unlock();

  // Copy what was in the current volume into the new volume.
  VolumeOps::bit_or(dnrrd, dlabel, cnrrd, clabel, mnrrd, mlabel);

  UndoHandle undo = new UndoReplaceLayer(this, "Undo Logical Or",
					 0, current_volume_, 0);
  push_undo(undo);

  extract_all_window_slices();
  rebuild_layer_buttons();
  redraw_all();

  hide_tool_panel();

  return CONTINUE_E;
}


BaseTool::propagation_state_e
Painter::CombineMaskXor(event_handle_t &)
{
  if (!check_for_active_label_volume("Combine mask with logical xor"))
  {
    return STOP_E;
  }

  if (!check_for_valid_mask("Combine mask with logical xor", true))
  {
    return STOP_E;
  }

  volume_lock_.lock();

  NrrdDataHandle mnrrd = mask_volume_->nrrd_handle_;
  const label_type mlabel = mask_volume_->label_;
  set_status("Using " + mask_volume_->name_ + " as mask.");

  NrrdDataHandle cnrrd = current_volume_->nrrd_handle_;
  const label_type clabel = current_volume_->label_;

  // Make a new child volume.
  const string name = current_volume_->name_ + " xor " + mask_volume_->name_;
  current_volume_ = create_new_label(current_volume_, name);
  NrrdDataHandle dnrrd = current_volume_->nrrd_handle_;
  const label_type dlabel = current_volume_->label_;

  volume_lock_.unlock();

  // Copy what was in the current volume into the new volume.
  VolumeOps::bit_xor(dnrrd, dlabel, cnrrd, clabel, mnrrd, mlabel);

  UndoHandle undo = new UndoReplaceLayer(this, "Undo Logical Xor",
					 0, current_volume_, 0);
  push_undo(undo);

  extract_all_window_slices();
  rebuild_layer_buttons();
  redraw_all();

  hide_tool_panel();

  return CONTINUE_E;
}


BaseTool::propagation_state_e
Painter::MoveScaleLayerUpdateGUI(event_handle_t &event)
{
  if (!current_volume_.get_rep())
  {
    set_status("No volume selected.");
    return STOP_E;
  }

  wxPanel *panel = Painter::global_seg3dframe_pointer_->CurrentToolPanel();
  if (panel)
  {
    MoveScaleToolChangeStruct *change = new MoveScaleToolChangeStruct;
    change->spacingx = current_volume_->nrrd_handle_->nrrd_->axis[1].spacing;
    change->spacingy = current_volume_->nrrd_handle_->nrrd_->axis[2].spacing;
    change->spacingz = current_volume_->nrrd_handle_->nrrd_->axis[3].spacing;

    change->originx = current_volume_->nrrd_handle_->nrrd_->axis[1].min;
    change->originy = current_volume_->nrrd_handle_->nrrd_->axis[2].min;
    change->originz = current_volume_->nrrd_handle_->nrrd_->axis[3].min;
    wxCommandEvent wxevent(wxEVT_MOVESCALETOOL_CHANGE, wxID_ANY);
    wxevent.SetClientData((void *)change);
    wxPostEvent(panel, wxevent);
  }

  return STOP_E;
}


BaseTool::propagation_state_e
Painter::MoveScaleLayer(event_handle_t &event)
{
  // TODO:  This won't work on labels, shared nrrd moves more than one.
  if (!current_volume_.get_rep())
  {
    set_status("No volume selected.");
    return STOP_E;
  }

  vector<int> cindex;
  for (size_t i = 0; i < windows_.size(); i++)
  {
    if (windows_[i]->visible())
    {
      const Point &center = windows_[i]->center_;
      cindex = current_volume_->world_to_index(center);
      break;
    }
  }
  
  double xspacing =
    get_vars()->get_double("Painter::movescale::spacing::x");
  double yspacing =
    get_vars()->get_double("Painter::movescale::spacing::y");
  double zspacing =
    get_vars()->get_double("Painter::movescale::spacing::z");

  // TODO: Validate spacings here

  double xorigin =
    get_vars()->get_double("Painter::movescale::origin::x");
  double yorigin =
    get_vars()->get_double("Painter::movescale::origin::y");
  double zorigin =
    get_vars()->get_double("Painter::movescale::origin::z");

  // Validate origin here?

  Nrrd *nrrd = current_volume_->nrrd_handle_->nrrd_;
  
  nrrd->axis[1].spacing = xspacing;
  nrrd->axis[2].spacing = yspacing;
  nrrd->axis[3].spacing = zspacing;

  nrrd->axis[1].min = xorigin;
  nrrd->axis[2].min = yorigin;
  nrrd->axis[3].min = zorigin;

  nrrd->axis[1].max =
    nrrd->axis[1].min + (nrrd->axis[1].size-1) * nrrd->axis[1].spacing;
  nrrd->axis[2].max =
    nrrd->axis[2].min + (nrrd->axis[2].size-1) * nrrd->axis[2].spacing;
  nrrd->axis[3].max =
    nrrd->axis[3].min + (nrrd->axis[3].size-1) * nrrd->axis[3].spacing;

  current_volume_->rebuild_transform();
  current_volume_->set_dirty();
  current_volume_->get_geom_group();

  if (current_volume_ == current_vrender_target_ &&
      !current_vrender_target_deferred_)
  {
    LoadVolumeEvent *lve =
      new LoadVolumeEvent(current_vrender_target_->nrrd_handle_, "", false);
    EventManager::add_event(lve);
  }

  if (cindex.size())
  {
    const Point center = current_volume_->index_to_world(cindex);
    for (size_t i = 0; i < windows_.size(); i++)
    {
      windows_[i]->center_ = center;
    }
  }

  extract_all_window_slices();
  Autoview(event);

  return CONTINUE_E;
}


BaseTool::propagation_state_e
Painter::MoveScaleAllLayers(event_handle_t &event)
{
  if (volumes_.empty()) { return CONTINUE_E; }

  const NrrdVolumeHandle &volume = volumes_[0];

  vector<int> cindex;
  for (size_t i = 0; i < windows_.size(); i++)
  {
    if (windows_[i]->visible())
    {
      const Point &center = windows_[i]->center_;
      cindex = volume->world_to_index(center);
      break;
    }
  }

  double xspacing =
    get_vars()->get_double("Painter::movescale::spacing::x");
  double yspacing =
    get_vars()->get_double("Painter::movescale::spacing::y");
  double zspacing =
    get_vars()->get_double("Painter::movescale::spacing::z");

  // TODO: Validate spacings here

  double xorigin =
    get_vars()->get_double("Painter::movescale::origin::x");
  double yorigin =
    get_vars()->get_double("Painter::movescale::origin::y");
  double zorigin =
    get_vars()->get_double("Painter::movescale::origin::z");

  // Validate origin here?

  for (size_t i = 0; i < volumes_.size(); i++)
  {
    const NrrdVolumeHandle &volume = volumes_[i];
    Nrrd *nrrd = volume->nrrd_handle_->nrrd_;
  
    nrrd->axis[1].spacing = xspacing;
    nrrd->axis[2].spacing = yspacing;
    nrrd->axis[3].spacing = zspacing;

    nrrd->axis[1].min = xorigin;
    nrrd->axis[2].min = yorigin;
    nrrd->axis[3].min = zorigin;

    nrrd->axis[1].max =
      nrrd->axis[1].min + (nrrd->axis[1].size-1) * nrrd->axis[1].spacing;
    nrrd->axis[2].max =
      nrrd->axis[2].min + (nrrd->axis[2].size-1) * nrrd->axis[2].spacing;
    nrrd->axis[3].max =
      nrrd->axis[3].min + (nrrd->axis[3].size-1) * nrrd->axis[3].spacing;

    volume->rebuild_transform();
    volume->set_dirty();
    volume->get_geom_group();
  }

  if (cindex.size())
  {
    const Point center = volume->index_to_world(cindex);
    for (size_t i = 0; i < windows_.size(); i++)
    {
      windows_[i]->center_ = center;
    }
  }

  if (!current_vrender_target_deferred_)
  {
    LoadVolumeEvent *lve =
      new LoadVolumeEvent(current_vrender_target_->nrrd_handle_, "", false);
    EventManager::add_event(lve);
  }

  extract_all_window_slices();
  Autoview(event);

  return CONTINUE_E;
}


BaseTool::propagation_state_e
Painter::ShowVisibleItem(event_handle_t &event)
{
  const string id = get_vars()->get_string("Painter::ShowVisibleItem::id");
  const string group =
    get_vars()->get_string("Painter::ShowVisibleItem::group");
  show_visible_item(id, group);
  return STOP_E;
}


BaseTool::propagation_state_e
Painter::RedrawAll(event_handle_t &event)
{
  // The variables are set by the ThrowSkinnerSignalEvent.
  redraw_all();
  return STOP_E;
}


BaseTool::propagation_state_e 
Painter::start_ITKConfidenceConnectedImageFilterTool(event_handle_t &event)
{
  tm_.set_tool(new ITKConfidenceConnectedImageFilterTool(this), 25);
  return CONTINUE_E;
}


BaseTool::propagation_state_e 
Painter::start_ITKThresholdSegmentationLevelSetImageFilterTool(event_handle_t &event)
{
  tm_.set_tool(new BrushTool(this), 25);
  tm_.set_tool(new ITKThresholdSegmentationLevelSetImageFilterTool(this), 26);
  return CONTINUE_E;
}

BaseTool::propagation_state_e 
Painter::start_ITKSpeedToPathGradientDescentImageFilterTool(event_handle_t &event)
{
  tm_.set_tool(new ITKSpeedToPathGradientDescentImageFilterTool(this), 25);
  return CONTINUE_E;
}

BaseTool::propagation_state_e 
Painter::start_ITKSpeedToPathRegularStepGradientDescentImageFilterTool(event_handle_t &event)
{
  tm_.set_tool(new ITKSpeedToPathRegularStepGradientDescentImageFilterTool(this), 25);
  return CONTINUE_E;
}

BaseTool::propagation_state_e 
Painter::start_ITKSpeedToPathIterateNeighborhoodImageFilterTool(event_handle_t &event)
{
  tm_.set_tool(new ITKSpeedToPathIterateNeighborhoodImageFilterTool(this), 25);
  return CONTINUE_E;
}


BaseTool::propagation_state_e 
Painter::start_ITKDiscreteGaussianImageFilterTool(event_handle_t &event)
{
  tm_.set_tool(new ITKDiscreteGaussianImageFilterTool(this), 25);
  return CONTINUE_E;
}

BaseTool::propagation_state_e 
Painter::start_ITKCannyEdgeDetectionImageFilterTool(event_handle_t &event)
{
  tm_.set_tool(new ITKCannyEdgeDetectionImageFilterTool(this), 25);
  return CONTINUE_E;
}

BaseTool::propagation_state_e 
Painter::FinishTool(event_handle_t &event)
{
  tm_.propagate_event(new FinishEvent());
  return CONTINUE_E;
}

BaseTool::propagation_state_e 
Painter::ClearTools(event_handle_t &event)
{
  // None of the tools currently do anything interesting on QuitEvent
  // except remove themselves.  Just do it directly to exactly which
  // tool slots we care about (rather than do bookkeeping to see how
  // many quit events we need to get rid of the right tools.
  tm_.rm_tool(24);
  tm_.rm_tool(25);
  tm_.rm_tool(26);
  return CONTINUE_E;
}


BaseTool::propagation_state_e 
Painter::SetLayer(event_handle_t &event)
{
  tm_.propagate_event(new SetLayerEvent());
  return CONTINUE_E;
}


BaseTool::propagation_state_e 
Painter::SetDataLayer(event_handle_t &event)
{
  tm_.propagate_event(new SetDataLayerEvent());
  return CONTINUE_E;
}

BaseTool::propagation_state_e 
Painter::SetThresholdLayer(event_handle_t &event)
{
  tm_.propagate_event(new SetThresholdLayerEvent());
  return CONTINUE_E;
}


int
Painter::isosurface_label_volumes(NrrdVolumes &volumes)
{
  Skinner::Var<bool> 
    isosurface_visible(get_vars(), "Painter::isosurface_visible", 1);
  int isocount = 0;
  for (unsigned int i = 0; i < volumes.size(); ++i)
  {
    NrrdVolumeHandle &volume = volumes[i];
    volume->lock.lock();
    if (volume->label_)
    {
      isocount++;

      GeomHandle isosurface;
      NrrdDataHandle nrrdh = volume->nrrd_handle_;
      isosurface = volume->isosurface_label();
      
      Matrix &tmat = volume->transform_;
      Transform transform(volume->min(), 
                          Vector(tmat.get(1,1), 0, 0),
                          Vector(0, tmat.get(2,2), 0),
                          Vector(0, 0, tmat.get(3,3)));

      MaterialHandle mat = volume->get_label_material();
      mat->transparency = 0.5;
      mat->specular = Color(0.1, 0.1, 0.1);
      mat->shininess = 2;
      GeomMaterial *colored_isosurface = new GeomMaterial(isosurface, mat);

      GeomTransform *color_trans_isosurface = \
        new GeomTransform(colored_isosurface, transform);

      GeomSkinnerVarSwitch *isosurface_switch = 
        new GeomSkinnerVarSwitch(color_trans_isosurface, isosurface_visible);
        
      GeomIndexedGroup *volume_group = volume->get_geom_group();
      if (volume_group) {
        volume_group->lock.lock();
        volume_group->addObj(isosurface_switch, -1);
        volume_group->lock.unlock();
      }
    }

    volume->lock.unlock();
  }

  return isocount;
}


BaseTool::propagation_state_e 
Painter::ShowOneIsosurface(event_handle_t &event)
{
  if (!check_for_active_label_volume("Isosurface"))
  {
    return STOP_E;
  }

  volume_lock_.lock();
  NrrdVolumes volumes;
  volumes.push_back(current_volume_);
  isosurface_label_volumes(volumes);
  volume_lock_.unlock();

  return CONTINUE_E;
}


BaseTool::propagation_state_e 
Painter::ShowAllIsosurfaces(event_handle_t &event)
{
  volume_lock_.lock();
  const int count = isosurface_label_volumes(volumes_);
  volume_lock_.unlock();

  if (count == 0)
  {
    set_status("No label volumes to isosurface.");
  }

  return CONTINUE_E;
}


BaseTool::propagation_state_e 
Painter::ShowVolumeRendering(event_handle_t &event)
{
  if (!(current_volume_.get_rep() && !current_volume_->label_))
  {
    set_status("No active data volume to render.");
    return STOP_E;
  }
   
  if (current_volume_ != current_vrender_target_)
  {
    current_vrender_target_ = current_volume_;
    current_vrender_target_deferred_ = true;
  }
  
  get_vars()->set_by_string("Painter::volume_visible", "1");
  ShowVolumeRendering2(event);

  return CONTINUE_E;
}


BaseTool::propagation_state_e 
Painter::ShowVolumeRendering2(event_handle_t &event)
{
  if (!SLIVR::ShaderProgramARB::shaders_supported())
  {
    set_status("Volume rendering is not supported by this hardware.");
    return STOP_E;
  }

  if (current_vrender_target_deferred_ &&
      current_vrender_target_.get_rep())
  {
    current_vrender_target_deferred_ = false;
    set_status("Setting up volume rendering.");
    LoadVolumeEvent *lve =
      new LoadVolumeEvent(current_vrender_target_->nrrd_handle_,
                          current_vrender_target_->full_path_);
    EventManager::add_event(lve);
    rebuild_layer_buttons();
    redraw_all();
  }

  return CONTINUE_E;
}


Nrrd*
Painter::do_CMedian_filter(Nrrd *nin, int radius)
{
  Nrrd *ntmp, *nout;
  nout = nrrdNew();
  ntmp = nin;
	
  if (nrrdCheapMedian(nout, ntmp, 0, 0, radius, 1.0, 2048)) {
    char *err = biffGetDone(NRRD);
    string errstr(err);
    free(err);
    throw "Error doing cheap median: " + errstr;
    return 0;
  }
	
  return nout;
}


BaseTool::propagation_state_e 
Painter::MedianFilterVolume(event_handle_t &event)
{
  if (!check_for_active_data_volume("Median filter"))
  {
    return STOP_E;
  }

  const int radius = get_vars()->get_int("MedianFilterVolume::radius");
  NrrdDataHandle nrrd_handle = current_volume_->get_nrrd_strip_subaxis();
	
  Nrrd *nin = nrrd_handle->nrrd_;
	
  // loop over the tuple axis, and do the median filtering for 
  // each scalar set independently, copy non scalar sets unchanged.
  //vector<string> elems;
  //nrrd_handle->get_tuple_indecies(elems);
	
	
  int min[NRRD_DIM_MAX], max[NRRD_DIM_MAX];
  for (unsigned int i = 0; i < nrrd_handle->nrrd_->dim; i++)
  {
    min[i] = 0;
    max[i] = nrrd_handle->nrrd_->axis[i].size - 1;
  }
	
	
  //! Slice a scalar out of the tuple axis and filter it. So for Vectors
  //! and Tensors, a component wise filtering occurs.
	
  if (nrrdKindSize(nrrd_handle->nrrd_->axis[0].kind) > 1)
  {
    vector<Nrrd*> out;
    for (unsigned int i = 0; i < nrrd_handle->nrrd_->axis[0].size; i++) 
    { 
      Nrrd *sliced = nrrdNew();
      if (nrrdSlice(sliced, nin, 0, i)) {
        char *err = biffGetDone(NRRD);
        cerr<<(string("Trouble with slice: ") + err);
        free(err);
      }
			
      Nrrd *nout_filtered;
      nout_filtered = do_CMedian_filter(sliced, radius);
      if (!nout_filtered) {
        throw "Error filtering, returning";
      }
      out.push_back(nout_filtered);
    }
    // Join the filtered nrrds along the first axis
    NrrdData *nrrd_joined = new NrrdData();
    if (nrrdJoin(nrrd_joined->nrrd_, &out[0], out.size(), 0, 1)) {
      char *err = biffGetDone(NRRD);
      string errstr(err);
      free(err);
      throw "Join Error: " + errstr;
    }
    NrrdDataHandle ntmp(nrrd_joined);
    nrrd_handle = ntmp;
  }
  else
  {
    Nrrd *nout_filtered;
    nout_filtered = do_CMedian_filter(nrrd_handle->nrrd_, radius);
    if (!nout_filtered)
    {
      throw "Error filtering, returning";
    }
    NrrdDataHandle ntmp(new NrrdData(nout_filtered));
    nrrd_handle = ntmp;
  }
	
  //  current_volume_->nrrd_handle_ = nrrd_handle;
    
  //  NrrdDataHandle nrrd_handle = new NrrdData();
  string newname = current_volume_->name_+" Median";
  make_layer(newname, nrrd_handle);

  UndoHandle undo = new UndoReplaceLayer(this, "Undo Median",
					 0, current_volume_, 0);
  push_undo(undo);
  
  extract_all_window_slices();
  rebuild_layer_buttons();
  redraw_all();

  hide_tool_panel();

  return CONTINUE_E;
}


BaseTool::propagation_state_e 
Painter::AbortFilterOn(event_handle_t &event)
{
  for (Filters_t::iterator f = filters_.begin(); 
       f != filters_.end(); ++f) {
    (*f)->stop();
  }
  return CONTINUE_E;
}


BaseTool::propagation_state_e 
Painter::CreateLabelVolume(event_handle_t &event)
{
  if (!current_volume_.get_rep()) return STOP_E;

  // Cache where our current volume is.
  size_t index = volumes_.size()-1;
  size_t i;
  for (i = 0; i < volumes_.size(); i++)
  {
    if (current_volume_ == volumes_[i])
    {
      index = i;
      break;
    }
  }

  create_new_label(current_volume_, current_volume_->name_);

  // Move the volume from the back to right after where the current
  // volume was.
  index++;
  for (i = volumes_.size()-1; i > index; i--)
  {
    volumes_[i] = volumes_[i-1];
  }
  volumes_[index] = current_volume_;

  UndoHandle undo = new UndoReplaceLayer(this, "Undo Create Label",
					 0, current_volume_, 0);
  push_undo(undo);

  extract_all_window_slices();
  rebuild_layer_buttons();
  redraw_all();

  return CONTINUE_E;
}


// David Brayford
BaseTool::propagation_state_e
Painter::LoadSession(const string &filename)
{
  bool success = false;
  string msg = "Possible reasons include bad file permissions or insufficient memory.";

  try {
    SessionReader reader(this);
    success = reader.load_session(filename);
  }
  catch (const char *e)
  {
    msg = string(e);
  }
  catch (string e)
  {
    msg = e;
  }
  catch (...)
  {
  }
  
  if (success)
  {
    set_status("Successfully loaded session: " + filename);
  }
  else
  {
    set_status_dialog("Load of " + filename + " failed!", msg);
    return STOP_E;
  }

  return CONTINUE_E;
}
  

BaseTool::propagation_state_e 
Painter::LoadSession(event_handle_t &event)
{
  const string filename =
    get_vars()->get_string("Painter::LoadSession::filename");
  return LoadSession(filename);
}
  

BaseTool::propagation_state_e
Painter::SaveSession(const string &filename)
{
  volume_lock_.lock();
  const bool success = SessionWriter::write_session(filename, this);
  volume_lock_.unlock();

  if (success) {
    set_status("Successfully saved session: " + filename);
  } else {
    set_status_dialog("Save of " + filename + " failed!",
                      "Possible reasons include bad file permissions or insufficient space.");
    return STOP_E;
  }

  return CONTINUE_E;
}
  

BaseTool::propagation_state_e 
Painter::SaveSession(event_handle_t &event)
{
  const string filename =
    get_vars()->get_string("Painter::SaveSession::filename");
  return SaveSession(filename);
}
  

BaseTool::propagation_state_e
Painter::SaveVolume(const string &filename, const  bool quantize)
{
  if (!current_volume_.get_rep()) {
    set_status("No volume selected.");
    return STOP_E;
  }

  NrrdVolumeHandle volume = current_volume_;

#if 0
  if (volume->label_)
  {
    if (volume->only_one_label())
    {
      if (volume->label_ != 1)
      {
        cout << "MOVING LABEL TO 1\n";
        VolumeOps::bit_copy(volume->nrrd_handle_, 1,
                            volume->nrrd_handle_, volume->label_);
        VolumeOps::bit_clear(volume->nrrd_handle_, volume->label_);
        volume->label_ = 1;
        // TODO:  This kills the outline geometry for some reason.
        volume->set_slices_dirty();
        extract_all_window_slices();
      }
      else
      {
        cout << "DOING NOTHING\n";
      }
    }
    else
    {
      cout << "COPYING BITS TO NEW VOLUME\n";
      NrrdDataHandle justone =
        VolumeOps::create_clear_nrrd(volume->nrrd_handle_, nrrdTypeUChar);
      VolumeOps::bit_copy(justone, 1, volume->nrrd_handle_, volume->label_);
      volume = new NrrdVolume(this, "Temporary Copy", justone);
    }
  }
#else
  if (volume->label_ && !(volume->only_one_label() && volume->label_ == 1))
  {
    NrrdDataHandle justone =
      VolumeOps::create_clear_nrrd(volume->nrrd_handle_, nrrdTypeUChar);
    VolumeOps::bit_copy(justone, 1, volume->nrrd_handle_, volume->label_);
    volume = new NrrdVolume(this, "Temporary Copy", justone);
  }
#endif

  if (volume->write(filename, quantize)) {
    set_status("Successfully saved volume: " + filename);
  }
  else
  {
    set_status_dialog("Save of " + filename + " failed!",
                      "Possible reasons include unknown file extension, bad file permissions or insufficient space.");
    return STOP_E;
  }
  
  return CONTINUE_E;
}

BaseTool::propagation_state_e
Painter::SaveVolumeWithHeader(const string &filename, const string &hdr_filename,
			      const string& series_desc, const  bool quantize)
{
  cout << "In Painter::SaveVolumeWithHeader" << endl;
  if (!current_volume_.get_rep()) {
    set_status("No volume selected.");
    return STOP_E;
  }

  NrrdVolumeHandle volume = current_volume_;

#if 0
  if (volume->label_)
  {
    if (volume->only_one_label())
    {
      if (volume->label_ != 1)
      {
        cout << "MOVING LABEL TO 1\n";
        VolumeOps::bit_copy(volume->nrrd_handle_, 1,
                            volume->nrrd_handle_, volume->label_);
        VolumeOps::bit_clear(volume->nrrd_handle_, volume->label_);
        volume->label_ = 1;
        // TODO:  This kills the outline geometry for some reason.
        volume->set_slices_dirty();
        extract_all_window_slices();
      }
      else
      {
        cout << "DOING NOTHING\n";
      }
    }
    else
    {
      cout << "COPYING BITS TO NEW VOLUME\n";
      NrrdDataHandle justone =
        VolumeOps::create_clear_nrrd(volume->nrrd_handle_, nrrdTypeUChar);
      VolumeOps::bit_copy(justone, 1, volume->nrrd_handle_, volume->label_);
      volume = new NrrdVolume(this, "Temporary Copy", justone);
    }
  }
#else
  if (volume->label_ && !(volume->only_one_label() && volume->label_ == 1))
  {
    NrrdDataHandle justone =
      VolumeOps::create_clear_nrrd(volume->nrrd_handle_, nrrdTypeUChar);
    VolumeOps::bit_copy(justone, 1, volume->nrrd_handle_, volume->label_);
    volume = new NrrdVolume(this, "Temporary Copy", justone);
  }
#endif
  cout << "about to call NrrdVolume::writeWithHeader" << endl;
  if (volume->writeWithHeader(filename, hdr_filename, series_desc, quantize)) {
    set_status("Successfully saved volume: " + filename + " with DICOM header: " + hdr_filename);
  }
  else
  {
    set_status_dialog("Save of " + filename + " with DICOM header failed!",
                      "Possible reasons include unknown file extension, bad file permissions or insufficient space.");
    return STOP_E;
  }
  
  return CONTINUE_E;
}


BaseTool::propagation_state_e 
Painter::SaveVolume(event_handle_t &event)
{
  const string filename =
    get_vars()->get_string("Painter::SaveVolume::filename");
  const bool quantize =
    atoi(get_vars()->get_string("Painter::SaveVolume::quantize").c_str());
  return SaveVolume(filename, quantize);
}

BaseTool::propagation_state_e 
Painter::SaveVolumeWithHeader(event_handle_t &event)
{
  const string filename =
    get_vars()->get_string("Painter::SaveVolumeWithHeader::filename");
  const string hdr_filename =
    get_vars()->get_string("Painter::SaveVolumeWithHeader::hdr_filename");
  const bool quantize =
    atoi(get_vars()->get_string("Painter::SaveVolumeWithHeader::quantize").c_str());
  const string series_desc =
    get_vars()->get_string("Painter::SaveVolumeWithHeader::description");
  return SaveVolumeWithHeader(filename, hdr_filename, series_desc, quantize);
}


struct valset_less_than
{
  bool operator()(const unsigned int s1, const unsigned int s2) const
  {
    return s1 < s2;
  }
};


typedef std::set<unsigned int, valset_less_than > valset_t;

static bool
compute_segmentation_values(valset_t &valset, NrrdVolumeHandle &vol)
{
  const size_t size = VolumeOps::nrrd_elem_count(vol->nrrd_handle_);

  unsigned int *ptr = (unsigned int *) vol->nrrd_handle_->nrrd_->data;
  
  for (size_t i = 0; i < size; i++)
  {
    valset.insert(ptr[i]);

    if ((i & 0xff) == 0 && valset.size() > MAX_IMPORT_SEGMENTATIONS)
    {
      return false;
    }
  }
  
  if (valset.size() > MAX_IMPORT_SEGMENTATIONS)
    return false;
  else
    return true;
}


BaseTool::propagation_state_e
Painter::ImportSegmentation(event_handle_t &event)
{
  const string filename =
    get_vars()->get_string("Painter::ImportSegmentation::filename");
  
  int status;
  NrrdVolumeHandle volume = load_volume<unsigned int>(filename, status);

  if (status == 1)
  {
    set_status_dialog("Import of " + filename + " failed!",
                      "Possible reasons include bad file permissions or insufficient memory.");
    return STOP_E;
  }
  else if (status == 2)
  {
    set_status_dialog("Import of " + filename + " failed!",
                      "This file has an unsupported format.");
    return STOP_E;
  }
  else if (status == 3)
  {
    set_status_dialog("Import of " + filename + " failed!",
                      "This is not a valid file.");
    return STOP_E;
  }

  // If we loaded an image and it appears to be part of a series
  // go ahead and load the whole series in.
  // TODO:  Pop up UI and ask for series parameters
  //        (start, stop, increment, dicom ordering)
  if (!volume.get_rep() || volume->nrrd_handle_->nrrd_->axis[3].size == 1)
  {
    const vector<string> filenames = get_filename_series(filename);
    if (filenames.size() > 1)
    {
      NrrdVolumeHandle volume2 = load_image_series<unsigned int>(filenames);
      if (volume2.get_rep()) {
        volume = volume2;
      }
    }
  }

  // Though a unsigned int was asked for in the load volume that only
  // happens if the ITK reader was used. If reading a mat lab file or
  // vff this may not be the case.
  volume->nrrd_handle_ =
    VolumeOps::convert_to(volume->nrrd_handle_, nrrdTypeUInt );

  valset_t segvals;
  if (!compute_segmentation_values(segvals, volume))
  {
    set_status_dialog("Import of " + filename + " failed!",
                      "This file contains more than " + to_string(MAX_IMPORT_SEGMENTATIONS) + " segmentations.\nIt's probably a data file.");
    return STOP_E;
  }

  for (valset_t::iterator itr = segvals.begin(); itr != segvals.end(); ++itr)
  {
    create_new_label(volume, volume->name_+" Label "+to_string(*itr));
    VolumeOps::extract_label_from_segmentation(current_volume_->nrrd_handle_,
                                               current_volume_->label_,
                                               volume->nrrd_handle_,
                                               *itr);
    extract_all_window_slices();
    rebuild_layer_buttons();
    redraw_all();
  }

  return CONTINUE_E;
}


static export_label_selection_data_t export_segmentation_data;


BaseTool::propagation_state_e
Painter::ExportSegmentation1(event_handle_t &event)
{
  // TODO: Turn this into a merge filter so the user can preview the
  // data and then do a volume->save rather than a blind export.
  // Requires discrete data volume support.
  size_t i;

  // Get a list of all the labels.
  // TODO:  Group the labels by size compatability.
  NrrdVolumes prelabels;
  for (i = 0; i < volumes_.size(); i++)
  {
    if (volumes_[i]->label_)
    {
      prelabels.push_back(volumes_[i]);
    }
  }

  if (prelabels.size() == 0)
  {
    set_status_dialog("Export failed.", "There are no label volumes to export.");
    return STOP_E;
  }

  // Send the list off to the user to determine which labels they want
  // and what numbers they should have.
  export_segmentation_data.nchoices_ = prelabels.size();
  if (export_segmentation_data.choices_)
  {
    delete [] export_segmentation_data.choices_;
  }
  export_segmentation_data.choices_ = new wxString[prelabels.size()];
  export_segmentation_data.selections_.clear();
  wxArrayInt selections;
  for (i = 0; i < prelabels.size(); i++)
  {
    export_segmentation_data.choices_[prelabels.size() - 1 - i] =
      std2wx(prelabels[i]->name_);
    if (prelabels[i]->visible())
    {
      export_segmentation_data.selections_.push_back(prelabels.size() - 1 - i);
    }
  }
  
  // Fire off event
  wxCommandEvent wxevent(wxEVT_COMMAND_EXPORT_LABEL_SELECTION, wxID_ANY);
  wxevent.SetClientData((void *)&export_segmentation_data);
  wxPostEvent(global_seg3dframe_pointer_, wxevent);

  return CONTINUE_E;
}

BaseTool::propagation_state_e
Painter::ExportSegmentationWithHeader1(event_handle_t &event)
{
  // TODO: Turn this into a merge filter so the user can preview the
  // data and then do a volume->save rather than a blind export.
  // Requires discrete data volume support.
  size_t i;

  // Get a list of all the labels.
  // TODO:  Group the labels by size compatability.
  NrrdVolumes prelabels;
  for (i = 0; i < volumes_.size(); i++)
  {
    if (volumes_[i]->label_)
    {
      prelabels.push_back(volumes_[i]);
    }
  }

  if (prelabels.size() == 0)
  {
    set_status_dialog("Export failed.", "There are no label volumes to export.");
    return STOP_E;
  }

  // Send the list off to the user to determine which labels they want
  // and what numbers they should have.
  export_segmentation_data.nchoices_ = prelabels.size();
  if (export_segmentation_data.choices_)
  {
    delete [] export_segmentation_data.choices_;
  }
  export_segmentation_data.choices_ = new wxString[prelabels.size()];
  export_segmentation_data.selections_.clear();
  wxArrayInt selections;
  for (i = 0; i < prelabels.size(); i++)
  {
    export_segmentation_data.choices_[prelabels.size() - 1 - i] =
      std2wx(prelabels[i]->name_);
    if (prelabels[i]->visible())
    {
      export_segmentation_data.selections_.push_back(prelabels.size() - 1 - i);
    }
  }
  
  // Fire off event
  wxCommandEvent wxevent(wxEVT_COMMAND_EXPORT_LABEL_SELECTION_WITH_HEADER, wxID_ANY);
  wxevent.SetClientData((void *)&export_segmentation_data);
  wxPostEvent(global_seg3dframe_pointer_, wxevent);

  return CONTINUE_E;
}


static export_statistics_selection_data_t export_statistics_data;

BaseTool::propagation_state_e
Painter::ExportStatistics(event_handle_t &event)
{
  // TODO: Turn this into a merge filter so the user can preview the
  // data and then do a volume->save rather than a blind export.
  // Requires discrete data volume support.
  size_t i;

  // Get a list of all the labels.
  // TODO:  Group the labels by size compatability.
  NrrdVolumes prelabels;
  for (i = 0; i < volumes_.size(); i++)
  {
    if (volumes_[i]->label_)
    {
      prelabels.push_back(volumes_[i]);
    }
  }

  if (prelabels.size() == 0)
  {
    set_status_dialog("Export failed.", "There are no label volumes to use for the statistices computation.");
    return STOP_E;
  }

  std::string stat_volume;
  for (i = 0; i < volumes_.size(); i++)
  {  
    if (volumes_[i] == current_volume_ && !(volumes_[i]->label_))
    {
      stat_volume = volumes_[i]->name_;
    }
  }
  
  if (stat_volume.empty())
  {
    set_status_dialog("Export failed.", "You need to select a data volume to compute the statistics on.");
    return STOP_E;
  }
  
  export_statistics_data.data_volume_ = stat_volume;
  
  // Send the list off to the user to determine which labels they want
  // and what numbers they should have.
  export_statistics_data.nchoices_ = prelabels.size();
  if (export_statistics_data.choices_)
  {
    delete [] export_statistics_data.choices_;
  }
  export_statistics_data.choices_ = new wxString[prelabels.size()];
  export_statistics_data.selections_.clear();
  wxArrayInt selections;
  for (i = 0; i < prelabels.size(); i++)
  {
    export_statistics_data.choices_[prelabels.size() - 1 - i] =
      std2wx(prelabels[i]->name_);
    if (prelabels[i]->visible())
    {
      export_statistics_data.selections_.push_back(prelabels.size() - 1 - i);
    }
  }
  
  // Fire off event
  wxCommandEvent wxevent(wxEVT_COMMAND_EXPORT_STATISTICS_SELECTION, wxID_ANY);
  wxevent.SetClientData((void *)&export_statistics_data);
  wxPostEvent(global_seg3dframe_pointer_, wxevent);

  return CONTINUE_E;
}


BaseTool::propagation_state_e
Painter::ExportSegmentation2(event_handle_t &event)
{
  size_t i;
  
  // Get a list of all the labels.
  // TODO:  Group the labels by size compatability.
  NrrdVolumes prelabels;
  for (i = 0; i < volumes_.size(); i++)
  {
    if (volumes_[i]->label_)
    {
      prelabels.push_back(volumes_[i]);
    }
  }

  // Something went wrong, labels different now.
  ASSERT(prelabels.size() == export_segmentation_data.nchoices_);

  NrrdVolumes labels;
  for (i = 0; i < export_segmentation_data.selections_.size(); i++)
  {
    labels.push_back(prelabels[prelabels.size() - 1 -
                               export_segmentation_data.selections_[i]]);
  }

  // Create a nrrd.
  NrrdDataHandle segmentation =
    VolumeOps::create_clear_nrrd(labels[0]->nrrd_handle_, nrrdTypeUChar);

  // Put the labels in the nrrd.
  // Merge them in back to front order so that the results match what
  // Seg3D displays.
  unsigned int counter = 1;
  for (int j = ((int)labels.size())-1; j >= 0; j--)
  {
    VolumeOps::merge_label_into_segmentation(segmentation, counter++,
                                             labels[j]->nrrd_handle_,
                                             labels[j]->label_);
  }

  // Write out the segmentation
  const string filename =
    get_vars()->get_string("Painter::ExportSegmentation::filename");

  NrrdVolumeHandle segvol =
    new NrrdVolume(this, "Temporary Segmentation", segmentation);

  if (!segvol->write(filename))
  {
    set_status_dialog("Save of " + filename + " failed!",
                      "Possible reasons include unknown file extension, bad file permissions or insufficient space.");
    return STOP_E;
  } 
  else 
  {
    set_status("Successfully saved " + filename);
  }

  return CONTINUE_E;
}

BaseTool::propagation_state_e
Painter::ExportSegmentationWithHeader2(event_handle_t &event)
{
  size_t i;
  
  // Get a list of all the labels.
  // TODO:  Group the labels by size compatability.
  NrrdVolumes prelabels;
  for (i = 0; i < volumes_.size(); i++)
  {
    if (volumes_[i]->label_)
    {
      prelabels.push_back(volumes_[i]);
    }
  }

  // Something went wrong, labels different now.
  ASSERT(prelabels.size() == export_segmentation_data.nchoices_);

  NrrdVolumes labels;
  // for now use the dicom header info from the first volume selected
  for (i = 0; i < export_segmentation_data.selections_.size(); i++)
  {
    labels.push_back(prelabels[prelabels.size() - 1 -
                               export_segmentation_data.selections_[i]]);
  }

  // Create a nrrd.
  NrrdDataHandle segmentation =
    VolumeOps::create_clear_nrrd(labels[0]->nrrd_handle_, nrrdTypeUChar);

  // Put the labels in the nrrd.
  // Merge them in back to front order so that the results match what
  // Seg3D displays.
  unsigned int counter = 1;
  for (int j = ((int)labels.size())-1; j >= 0; j--)
  {
    VolumeOps::merge_label_into_segmentation(segmentation, counter++,
                                             labels[j]->nrrd_handle_,
                                             labels[j]->label_);
  }

  // Write out the segmentation
  const string filename =
    get_vars()->get_string("Painter::ExportSegmentationWithHeader::filename");
  const string hdr_filename =
    get_vars()->get_string("Painter::ExportSegmentationWithHeader::hdr_filename");
  const string series_desc = 
    get_vars()->get_string("Painter::ExportSegmentationWithHeader::description");

  NrrdVolumeHandle segvol =
    new NrrdVolume(this, "Temporary Segmentation", segmentation);
  // which

  if (!segvol->writeWithHeader(filename, hdr_filename, series_desc))
  {
    set_status_dialog("Save of " + filename + " with DICOM header failed!",
                      "Possible reasons include unknown file extension, bad file permissions or insufficient space.");
    return STOP_E;
  } 
  else 
  {
    set_status("Successfully saved " + filename + " with DICOM header");
  }

  return CONTINUE_E;
}


BaseTool::propagation_state_e
Painter::ExportStatistics2(event_handle_t &event)
{
  size_t i;
  
  // Get a list of all the labels.
  // TODO:  Group the labels by size compatability.
  NrrdVolumes prelabels;
  
  NrrdVolumeHandle DataVolume; 
  for (i = 0; i < volumes_.size(); i++)
  {
    if (volumes_[i]->label_)
    {
      prelabels.push_back(volumes_[i]);
    }
    if (volumes_[i]->name_ == export_statistics_data.data_volume_)
    {
      DataVolume = volumes_[i];
    }
  }

  NrrdVolumes labels;
  for (i = 0; i < export_statistics_data.selections_.size(); i++)
  {
    labels.push_back(prelabels[prelabels.size() - 1 -
                               export_statistics_data.selections_[i]]);
  }


  // Something went wrong, labels different now.
  ASSERT(prelabels.size() == export_segmentation_data.nchoices_);
  ASSERT(DataVolume.get_rep())

  const std::string filename =
    get_vars()->get_string("Painter::ExportStatistics::filename");

  std::ofstream statfile;
  try
  {
    statfile.open(filename.c_str());
  }
  catch(...)
  {
    set_status_dialog("Save of " + filename + " failed!",
                      "Possible reasons include bad file permissions or insufficient space.");
    return STOP_E;
  }

  statfile << "Statistics of "<< DataVolume->name_<< std::endl;
  statfile << "labelname : min max mean median stddev volume " << std::endl;

  for (i = 0; i < labels.size(); i++)
  {
    if (labels[i]->label_)
    {
      NrrdStats stats = VolumeOps::nrrd_statistics(DataVolume->nrrd_handle_,labels[i]->nrrd_handle_,(1 << labels[i]->bit()) );
      statfile << labels[i]->name_ << " : " << stats.min() << " " <<
            stats.max() << " " << stats.mean() << " " << stats.stddev() << " " <<
            stats.volume() << std::endl;
    }
  }

  return CONTINUE_E;
}



BaseTool::propagation_state_e
Painter::Undo(event_handle_t &)
{
  undo_manager_->undo_last();

  return CONTINUE_E;
}


// Unused
BaseTool::propagation_state_e 
Painter::RebuildLayers(event_handle_t &)
{
  rebuild_layer_buttons();
  return CONTINUE_E;
}


NrrdDataHandle
Painter::MatlabNrrd_reader(const string &filename)
{
  NrrdDataHandle mh = 0;
  MatlabIO::matlabfile mf;
  MatlabIO::matlabconverter mc(0);
  MatlabIO::matlabarray ma;
  int numarrays;
  std::string dummytext;

  try
  {
      mf.open(filename, "r");
      numarrays = mf.getnummatlabarrays();
      for (int p=0; p<numarrays; p++)
      {
        ma = mf.getmatlabarrayinfo(p);
        if (mc.sciNrrdDataCompatible(ma, dummytext)) 
        { 
          ma = mf.getmatlabarray(p);
          mc.mlArrayTOsciNrrdData(ma, mh); break; 
        }
      }
      mf.close();
  }
  catch (...)
  {
    mh = 0;
  }

  return mh;
}


NrrdDataHandle
Painter::vff_reader(const string &filename)
{
  // Read the header from the vff file.
  int nrrdSkip = 0;
  int dim = -1;
  size_t size[3];
  int nbits = 16;
  ifstream vffFileStream(filename.c_str(), ios::binary);
	

  if (! vffFileStream)
  {
    throw "VFF Reader:: Could not open file " + filename;
  }
	
  int lineCount = 0;
  int foundFormFeed = 0;	
  char temp[1025];
  while (! vffFileStream.eof() && !foundFormFeed)
  {
    vffFileStream.getline(temp,1024);
    lineCount++;
		
    if (temp[0] == '\f')
    {
      nrrdSkip = lineCount;
      foundFormFeed = 1;
    }
		
    if (strncmp(temp,"rank=",5) == 0)
    {
      istringstream rankLine(&temp[5]);
      rankLine >> dim;
    }
		
    if (strncmp(temp,"bits=",5) == 0)
    {
      istringstream bitsLine(&temp[5]);
      bitsLine >> nbits;
    }
		
    if (strncmp(temp,"size=",5) == 0)
    {
      istringstream sizeLine(&temp[5]);
      sizeLine >> size[0] >> size[1] >> size[2];	    
    }
  }
	
  vffFileStream.close();
	
  // Now create the nrrd from the information contained in the vff header.
  NrrdDataHandle nout = new NrrdData();
  NrrdIoState *nio = nrrdIoStateNew();

  // Set first axes kind to be what user specified.
  nout->nrrd_->axis[0].kind = nrrdKindUnknown;
	
  nout->nrrd_->type = nrrdTypeShort;
  if (nbits == 8){
    nout->nrrd_->type=nrrdTypeUChar;
  }
    
  nout->nrrd_->dim = dim;
	
  // Set the rest of the axes to have a kind of nrrdKindDomain.
  for (int i=1; i<dim; i++) {
    nout->nrrd_->axis[i].kind = nrrdKindDomain;
  }
    
  // Set the nrrd's labels, sizes, and spacings that were
  // parsed in the parse_gui_vars function.
  for(int i=0; i<dim; i++) {
    nout->nrrd_->axis[i].size = size[i];
    nout->nrrd_->axis[i].spacing = 1;
    nrrdAxisInfoMinMaxSet(nout->nrrd_, i, nrrdCenterUnknown);
  }
    
  nrrdIoStateInit(nio);

  nio->lineSkip = nrrdSkip;
  nio->encoding = nrrdEncodingArray[nrrdEncodingTypeRaw];
    
  // Assume only reading in a single file
  FILE *dataFile;
  dataFile = airFopen(filename.c_str(), stdin, "rb");
  nrrdLineSkip(dataFile, nio);
	
  nout->nrrd_->data = calloc(nrrdElementNumber(nout->nrrd_), 
                              nrrdElementSize(nout->nrrd_));
	
  nio->encoding->read(dataFile, nout->nrrd_->data, 
                       nrrdElementNumber(nout->nrrd_), 
                       nout->nrrd_, nio);
    
  airFclose(dataFile);

  // 16 bit VFF is always big endian.
  if (nrrdElementSize(nout->nrrd_) > 1
      && AIR_ENDIAN != airEndianBig)
  {
    nrrdSwapEndian(nout->nrrd_);
  }

  NrrdDataHandle nrrd_float = new NrrdData();
  if (nrrdConvert(nrrd_float->nrrd_, nout->nrrd_, nrrdTypeFloat))
  {
    char *err = biffGetDone(NRRD);
    string errstr(err);
    free(err);
    throw "VFF Reader converting to float: " + errstr;
  }

  return nrrd_float;
}


bool
Painter::MatlabNrrd_writer(NrrdDataHandle mh, const char *filename)
{
  MatlabIO::matlabfile mf;
  MatlabIO::matlabconverter mc(0);
  MatlabIO::matlabarray ma;
  std::string name;
 
  try
  {
      mc.converttostructmatrix();
      mc.sciNrrdDataTOmlArray(mh,ma);
      mh->get_property("name", name);
      if ((name=="")||(!mc.isvalidmatrixname(name))) name = "scirunnrrd";
      mf.open(std::string(filename),"w");
      mf.putmatlabarray(ma,name);
      mf.close();
  }
  catch (...)
  {
    return(false);
  }
  return(true);
}


bool
Painter::check_for_active_data_volume(const char *msg)
{
  if (!(current_volume_.get_rep() && !current_volume_->label_))
  {
    int count = 0;
    NrrdVolumeHandle dvol;
    for (size_t i = 0; i < volumes_.size(); i++)
    {
      if (volumes_[i]->label_ == 0)
      {
        count++;
        dvol = volumes_[i];
      }
    }
    if (count == 1)
    {
      current_volume_ = dvol;
      set_status(string(msg) + " using " + current_volume_->name_ +
                 " as active data volume.");
      return true;
    }
    else
    {
      set_status(string(msg) + " requires an active data volume.");
      return false;
    }
  }
  return true;
}


bool
Painter::check_for_active_label_volume(const char *msg)
{
  if (!(current_volume_.get_rep() && current_volume_->label_))
  {
    set_status(string(msg) + " requires an active label volume.");
    return false;
  }
  return true;
}


bool
Painter::check_for_valid_mask(const char *msg, bool strict)
{
  if (mask_volume_.get_rep() && current_volume_.get_rep())
  {
    if (VolumeOps::compare_nrrd_info(current_volume_->nrrd_handle_,
                                     mask_volume_->nrrd_handle_))
    {
      return true;
    }
    if (strict)
    {
      set_status(string(msg) + " the mask and volume must be the same size.");
    }
    else
    {
      set_status(string(msg) + " mask wrong size.  Not using it.");
    }
    return false;
  }
  if (strict)
  {
    set_status(string(msg) + " requires a mask to be selected.");
  }
  return false;
}


void
Painter::set_status_dialog(const string &s1, const string &s2)
{
  const string s = s1 + "\n" + s2;
  wxString wxstr(s.c_str(), wxConvUTF8);
  wxCommandEvent event(wxEVT_COMMAND_OK_DIALOG, wxID_ANY);
  event.SetString(wxstr);
  wxPostEvent(global_seg3dframe_pointer_, event);
}


void
Painter::set_status(const string &s)
{
  wxString wxstr(s.c_str(), wxConvUTF8);
  wxCommandEvent event(wxEVT_COMMAND_STATUS_BAR_TEXT_CHANGE, wxID_ANY);
  event.SetString(wxstr);
  wxPostEvent(global_seg3dframe_pointer_, event);
}


void
Painter::hide_tool_panel()
{
  wxCommandEvent event(wxEVT_COMMAND_HIDE_TOOL, wxID_ANY);
  wxPostEvent(global_seg3dframe_pointer_, event);
}


void
Painter::clear_all_volumes()
{
  volume_lock_.lock();

  // Clear the volume rendering.
  if (!current_vrender_target_deferred_)
  {
    NrrdDataHandle empty(0);
    LoadVolumeEvent *lve = new LoadVolumeEvent(empty, "");
    EventManager::add_event(lve);
  }
  current_vrender_target_ = 0;
  current_vrender_target_deferred_ = false;

  // Clear the mask layer.
  mask_volume_ = 0;

  // Clear the current layer.
  current_volume_ = 0;

  for (size_t i = 0; i < volumes_.size(); i++)
  {
    // Clear the 3D graphics associated with this layer.
    EventManager::add_event(new SceneGraphEvent(0, volumes_[i]->name_));

    // Remove the label allocation internal bookkeeping for this layer.
    volumes_[i]->unparent();
  }

  volumes_.clear();

  volume_lock_.unlock();

  // Clear the undo stack.
  undo_manager_->clear_undo();

  rebuild_layer_buttons();
  extract_all_window_slices();
  redraw_all();
}


size_t
Painter::remove_volume(NrrdVolumeHandle &layer, bool undo)
{
  size_t result = 0;

  // Clear the volume rendering.
  if (layer == current_vrender_target_)
  {
    // If it's been deferred then no need to update the scene graph.
    if (!current_vrender_target_deferred_)
    {
      NrrdDataHandle empty(0);
      LoadVolumeEvent *lve = new LoadVolumeEvent(empty, "");
      EventManager::add_event(lve);
    }
    current_vrender_target_ = 0;
    current_vrender_target_deferred_ = false;
  }
  
  volume_lock_.lock();
  try {

    // Clear the mask layer if we're deleting it.
    if (layer == mask_volume_)
    {
      mask_volume_ = 0;
    }

    // Make the new volume list with the deleted volume missing.
    // TODO:  Put the old layer in an undo list instead of deleting.
    int loc = -1;
    NrrdVolumes newvolumes;
    for (size_t i = 0; i < volumes_.size(); i++)
    {
      if (volumes_[i] != layer)
      {
        newvolumes.push_back(volumes_[i]);
      }
      else
      {
        loc = i;
      }
    }
    
    if (undo && loc != -1)
    {
      result = loc;
      UndoHandle undo = new UndoReplaceLayer(this, "Undo Layer Delete",
                                             volumes_[loc], 0, loc);
      undo_manager_->push_undo(undo);
    }

    volumes_ = newvolumes;

    // Clear the 3D graphics associated with this layer.
    EventManager::add_event(new SceneGraphEvent(0, layer->name_));

    if (layer == current_volume_)
    {
      // Set the current layer to be near the deleted layer.
      if (loc == -1 || volumes_.size() == 0)
      {
        current_volume_ = 0;
      }
      else
      {
        loc = Min((int)volumes_.size()-1, loc);
        current_volume_ = volumes_[loc];
      }
    }

    // TODO: Unparent?
  }
  catch(...) {
    volume_lock_.unlock();
    throw;
  }
  volume_lock_.unlock();

  return result;
}


void
Painter::insert_volume(NrrdVolumeHandle &volume, size_t index)
{
  if (index < volumes_.size())
  {
    volumes_.insert(volumes_.begin()+index, volume);
  }
  else
  {
    volumes_.push_back(volume);
  }
}


void
Painter::push_undo(UndoHandle &undo)
{
  undo_manager_->push_undo(undo);
}


} // end namespace SCIRun
