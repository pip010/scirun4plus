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
//    File   : ITKNeighborhoodConnectedImageFilterTool.cc
//    Author : McKay Davis
//    Date   : Tue Sep 26 18:44:34 2006

#include <Applications/Seg3D/Painter.h>

#include <Applications/Seg3D/ITKNeighborhoodConnectedImageFilterTool.h>
#include <sci_gl.h>

#include <itkImageToImageFilter.h>
#include <itkCommand.h>
#include <itkNeighborhoodConnectedImageFilter.h>
#include <itkConnectedThresholdImageFilter.h>


namespace SCIRun {

SeedTool::seeds_t ITKNeighborhoodConnectedImageFilterTool::seed_cache_;

ITKNeighborhoodConnectedImageFilterTool::
ITKNeighborhoodConnectedImageFilterTool(Painter *painter) :
  SeedTool("ITKNeighborhoodConnectedImageFilterTool", painter)
{
  seeds_ = seed_cache_;
  painter_->redraw_all();
}


ITKNeighborhoodConnectedImageFilterTool::~ITKNeighborhoodConnectedImageFilterTool()
{
  seed_cache_ = seeds_;
}


void
ITKNeighborhoodConnectedImageFilterTool::run_filter()
{
  if (!painter_->check_for_active_data_volume("Neighborhood connected image filter"))
  {
    return;
  }

  if (seeds_.empty())
  {
    painter_->set_status("No seed points were selected.");
    return;
  }

  // Save the source.
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
  const string name = painter_->current_volume_->name_ + " Neighborhood Connected";
  painter_->create_new_label(painter_->current_volume_, name);

  painter_->rebuild_layer_buttons();
  painter_->volume_lock_.unlock();

  painter_->start_progress();

  label_type value = painter_->current_volume_->label_;

  // Get the data pointers.
  NrrdDataHandle snrrd = source_volume->nrrd_handle_;
  NrrdDataHandle dnrrd = painter_->current_volume_->nrrd_handle_;
  ASSERT(snrrd->nrrd_->type == nrrdTypeFloat);
  ASSERT(dnrrd->nrrd_->type == LabelNrrdType);

  float *srcdata = (float *)snrrd->nrrd_->data;
  label_type *dstdata = (label_type *)dnrrd->nrrd_->data;

  // Array to hold which indices to visit next
  vector<vector<int> > todo[2];
  int current = 0;

  // Allocated a nrrd to mark where the flood fill has visited
  NrrdDataHandle done_handle = new NrrdData();
  size_t size[NRRD_DIM_MAX];
  size[0] = source_volume->nrrd_handle_->nrrd_->axis[0].size;
  size[1] = source_volume->nrrd_handle_->nrrd_->axis[1].size;
  size[2] = source_volume->nrrd_handle_->nrrd_->axis[2].size;
  size[3] = source_volume->nrrd_handle_->nrrd_->axis[3].size;
  nrrdAlloc_nva(done_handle->nrrd_, nrrdTypeUChar, 4, size);

  // Set the visited nrrd to empty.
  memset(done_handle->nrrd_->data, 0, 
         source_volume->nrrd_handle_->nrrd_->axis[0].size *
         source_volume->nrrd_handle_->nrrd_->axis[1].size *
         source_volume->nrrd_handle_->nrrd_->axis[2].size * 
         source_volume->nrrd_handle_->nrrd_->axis[3].size);

  // Copy the mask into the visited array so that it won't be visited.
  if (painter_->check_for_valid_mask("Neighborhood connected image filter"))
  {
    VolumeOps::bit_copy(done_handle, 1,
                        painter_->mask_volume_->nrrd_handle_,
                        painter_->mask_volume_->label_);
  }

  unsigned char *done_data = (unsigned char *)done_handle->nrrd_->data;

  // Add the seeds.
  float min = 0, max = 0;
  unsigned int axes = 0;
  for (size_t i = 0; i < iseeds.size(); ++i)
  {
    if (!axes) { axes = iseeds[i].size(); }

    const size_t offset =
      VolumeOps::index_to_offset(snrrd->nrrd_, iseeds[i]);

    const float fillval = srcdata[offset];
    if (i == 0)
    {
      min = max = fillval;
    }
    else
    {
      min = Min(min, fillval);
      max = Max(max, fillval);
    }

    // Mask this seed point.
    if (done_data[offset]) continue;
    
    // Push back the seed point
    todo[current].push_back(iseeds[i]);
  }

  // Set up the progress counters.
  size_t pcounter = 0;
  size_t pcounter_size = 0;
  for (unsigned int a = 1; a < axes; ++a)
  {
    pcounter_size += dnrrd->nrrd_->axis[a].size;
  }

  // Expand the leading edge until there is none.
  while (!todo[current].empty())
  {
    // Update progress.
    if ((++pcounter & 0xf) == 0)
    {
      painter_->update_progress(pcounter * 100 / pcounter_size);
    }
        
    current = !current;
    todo[current].clear();
    for (unsigned int i = 0; i < todo[!current].size(); ++i)
    {
      dstdata[VolumeOps::index_to_offset(dnrrd->nrrd_,
                                         todo[!current][i])] |= value;
    }
    
    // For each axis
    for (unsigned int i = 0; i < todo[!current].size(); ++i)
    {
      for (unsigned int a = 1; a < axes; ++a)
      {
        // Visit the previous and next neighbor indices along this axis
        for (int dir = -1; dir < 2; dir +=2) {

          // Neighbor index starts as current index
          vector<int> neighbor_index = todo[!current][i];

          // Index adjusted in direction we're looking at along axis
          neighbor_index[a] = neighbor_index[a] + dir;

          // Bail if this index is outside the volume
          if (!source_volume->index_valid(neighbor_index)) continue;
          
          const size_t offset =
            VolumeOps::index_to_offset(snrrd->nrrd_, neighbor_index);
          
          // Check to see if flood fill has already been here
          // Bail if the voxel has been visited
          if (done_data[offset]) continue;
          
          // Now check to see if this pixel is a candidate to be filled
          const float nbrval = srcdata[offset];
          
          // Bail if the voxel is outside the flood fill range
          if (nbrval < min || nbrval > max) continue;

          // Mark this voxel as visited
          done_data[offset] = 1;
          
          todo[current].push_back(neighbor_index);
        }
      }
    }
  }
  painter_->finish_progress();

  UndoHandle undo =
    new UndoReplaceLayer(painter_, "Undo Neighborhood Connected",
                         0, painter_->current_volume_, 0);
  painter_->push_undo(undo);

  painter_->extract_all_window_slices();
  painter_->redraw_all();

  painter_->hide_tool_panel();
}


}
