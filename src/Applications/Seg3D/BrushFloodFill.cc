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
//    File   : BrushFloodFill.cc
//    Author : Michael Callahan
//    Date   : March 2008

#include <Applications/Seg3D/Painter.h>

#include <Applications/Seg3D/BrushFloodFill.h>
#include <Core/Events/keysyms.h>
#include <sci_gl.h>
 

namespace SCIRun {

BrushFloodFill::
BrushFloodFill(Painter *painter) :
  SeedTool("BrushFloodFill", painter),
  value_(0),
  erase_(false)
{
  set_just_one_seed_mode(true);
  set_button(2);
}


void
BrushFloodFill::run_filter()
{
  // This function should not be run directly.  We don't finish it,
  // rather it is called multiple times interactively while painting.
  cout << "BrushFloodFill should not be run directly.\n";
}


BaseTool::propagation_state_e 
BrushFloodFill::process_event(event_handle_t event)
{
  KeyEvent *keyevent = dynamic_cast<KeyEvent *>(event.get_rep());
  if (keyevent && keyevent->get_keyval() == SCIRun_f)
  {
    // TODO: Set the seed point at the current pointer position.
    flood_fill_slice(false);
    return CONTINUE_E;
  }

  return SeedTool::process_event(event);
}


void
BrushFloodFill::flood_fill(Painter *painter,
                           NrrdDataHandle dnrrd, unsigned int dlabel,
                           NrrdDataHandle snrrd, unsigned int slabel,
                           NrrdDataHandle mnrrd, unsigned int mlabel,
                           bool erase, const vector<vector<int> > &seeds)
{
  // Depth first search using a double buffered queue (vector) for the
  // expanding edge.

  // Set up the data pointers.
  label_type *srcdata = (label_type *)snrrd->nrrd_->data;
  label_type *dstdata = (label_type *)dnrrd->nrrd_->data;
  const bool masking = mnrrd.get_rep() != NULL;
  label_type *mskdata = 0;
  if (masking) { mskdata = (label_type *)mnrrd->nrrd_->data; }

  // Two vectors to hold which indices to visit next.
  vector<vector<int> > todo[2];
  int current = 0;

  unsigned int axes = 0;

  // Add all the seed points to the queue.
  for (size_t i = 0; i < seeds.size(); ++i)
  {
    if (!axes) axes = seeds[i].size();

    // Bail if this index is outside the volume
    for (unsigned int a = 1; a < axes; ++a)
    {
      if (seeds[i][a] < 0 ||
          seeds[i][a] >= (int)(dnrrd->nrrd_->axis[a].size)) continue;
    }

    const size_t offset =
      VolumeOps::index_to_offset(snrrd->nrrd_, seeds[i]);

    if (masking && (mskdata[offset] & mlabel)) continue;

    if (erase)
    {
      if (!(srcdata[offset] & slabel)) continue;
      if (!(dstdata[offset] & dlabel)) continue;
      dstdata[offset] &= ~dlabel;
      todo[current].push_back(seeds[i]);
    }
    else
    {
      if (srcdata[offset] & slabel) continue;
      if (dstdata[offset] & dlabel) continue;
      dstdata[offset] |= dlabel;
      todo[current].push_back(seeds[i]);
    }
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
      painter->update_progress(pcounter * 100 / pcounter_size);
    }

    current = !current;
    todo[current].clear();

    // For each axis
    for (unsigned int i = 0; i < todo[!current].size(); ++i)
    {
      for (unsigned int a = 1; a < axes; ++a)
      {
        // Visit the previous and next neighbor indices along this axis
        for (int dir = -1; dir < 2; dir +=2)
        {
          // Neighbor index starts as current index
          vector<int> neighbor_index = todo[!current][i];

          // Index adjusted in direction we're looking at along axis
          neighbor_index[a] = neighbor_index[a] + dir;

          // Bail if this index is outside the volume
          if (neighbor_index[a] < 0 ||
              neighbor_index[a] >= (int)(dnrrd->nrrd_->axis[a].size)) continue;
          
          const size_t offset =
            VolumeOps::index_to_offset(snrrd->nrrd_, neighbor_index);
          
          if (masking && (mskdata[offset] & mlabel)) continue;
          
          if (erase)
          {
            if (!(srcdata[offset] & slabel)) continue;
            if (!(dstdata[offset] & dlabel)) continue;
            dstdata[offset] &= ~dlabel;
          }
          else
          {
            if (srcdata[offset] & slabel) continue;
            if (dstdata[offset] & dlabel) continue;
            dstdata[offset] |= dlabel;
          }

          todo[current].push_back(neighbor_index);
        }
      }
    }
  }
}


BaseTool::propagation_state_e
BrushFloodFill::flood_fill_slice(bool erase)
{
  if (!painter_->check_for_active_label_volume(erase?"Brush flood erase":"Brush flood fill"))
  {
    return STOP_E;
  }

  if (seeds_.empty())
  {
    painter_->set_status("Use the middle mouse button to set the seed point.");
    return STOP_E;
  }

  if (last_seed_window_ == NULL) return STOP_E;

  painter_->volume_lock_.lock();

  // Get the current volume.
  NrrdVolumeHandle &vol = painter_->current_volume_;

  vector<vector<int> > seeds;
  convert_seeds_to_indices(seeds, vol);
  if (seeds.empty())
  {
    painter_->set_status("The seed point is outside of the volume.");
    painter_->volume_lock_.unlock();
    return STOP_E;
  }

  vol->lock.lock();

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
  if (painter_->check_for_valid_mask(erase?"Brush flood erase":"Brush flood fill"))
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

  // Put this slice back
  const int axis = slice->get_axis();
  for (size_t i = 0; i < seeds.size(); i++)
  {
    seeds[i].erase(seeds[i].begin() + axis);
  }

  BrushFloodFill::flood_fill(painter_,
                             cnrrd, clabel, cnrrd, clabel,
                             mnrrd, mlabel, erase, seeds);

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
                         erase?"Undo Erase Fill":"Undo Brush Fill",
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

#if 0
BaseTool::propagation_state_e
BrushFloodFill::flood_fill_volume(bool erase)
{
  NrrdDataHandle mnrrd;
  unsigned int mlabel = 0;
  if (painter_->check_for_valid_mask("Flood fill volume"))
  {
    mnrrd = painter_->mask_volume_->nrrd_handle_;
    mlabel = painter_->mask_volume_->label_;
  }

  flood_fill(painter_,
             painter_->current_volume_->nrrd_handle_,
             painter_->current_volume_->label_,
             painter_->current_volume_->nrrd_handle_,
             painter_->current_volume_->label_,
             mnrrd, mlabel, erase, seeds_);

  // Redraw everything after completion.
  painter_->extract_all_window_slices();
  painter_->redraw_all();

  return CONTINUE_E;
}
#endif

}
