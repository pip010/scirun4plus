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
//    File   : FloodFillCopyTool.cc
//    Author : Michael Callahan
//    Date   : May 2008

#include <Applications/Seg3D/FloodFillCopyTool.h>
#include <Applications/Seg3D/Seg3DFrame.h>
#include <Applications/Seg3D/Painter.h>
#include <Applications/Seg3D/VolumeOps.h>
#include <Core/Util/Assert.h>

namespace SCIRun {

SeedTool::seeds_t FloodFillCopyTool::seed_cache_;


FloodFillCopyTool::FloodFillCopyTool(Painter *painter)
  : SeedTool("FloodFillCopyTool::", painter)
{
  seeds_ = seed_cache_;
  painter_->redraw_all();
}


FloodFillCopyTool::~FloodFillCopyTool()
{
  seed_cache_ = seeds_;
}

void
FloodFillCopyTool::run_filter()
{
  if (!painter_->check_for_active_label_volume("Connected component filter"))
  {
    return;
  }

  if (seeds_.empty())
  {
    painter_->set_status("Connected component filter requires at least one seed point.");
    return;
  }

  painter_->volume_lock_.lock();

  // Save off the source.
  NrrdVolumeHandle source_volume = painter_->current_volume_;

  vector<vector<int> > iseeds;
  convert_seeds_to_indices(iseeds, source_volume);

  if (iseeds.empty())
  {
    painter_->set_status("All of the seed points were outside the volume.");
    painter_->volume_lock_.unlock();
    return;
  }

  // Make a new label volume
  const string name = painter_->current_volume_->name_ + " Connected Component";
  painter_->create_new_label(painter_->current_volume_, name);

  painter_->rebuild_layer_buttons();
  painter_->volume_lock_.unlock();

  NrrdDataHandle mnrrd;
  unsigned int mlabel = 0;
  if (painter_->check_for_valid_mask("Connected component filter"))
  {
    mnrrd = painter_->mask_volume_->nrrd_handle_;
    mlabel = painter_->mask_volume_->label_;
  }

  // TODO:  This is inefficient but easy to implement.  It shouldn't be
  // hard to redo this without the two inverts (just changes the inside
  // test in flood_fill.

  painter_->start_progress();

  // Invert the source so that the selected thing is a hole.
  VolumeOps::bit_invert(source_volume->nrrd_handle_,
                        source_volume->label_,
                        source_volume->nrrd_handle_,
                        source_volume->label_);

  // Flood fill the source 'hole' into the destination.
  BrushFloodFill::flood_fill(painter_,
                             painter_->current_volume_->nrrd_handle_,
                             painter_->current_volume_->label_,
                             source_volume->nrrd_handle_,
                             source_volume->label_,
                             mnrrd, mlabel, false, iseeds);

  painter_->update_progress(90);

  // Reinvert the source back to it's original.
  VolumeOps::bit_invert(source_volume->nrrd_handle_,
                        source_volume->label_,
                        source_volume->nrrd_handle_,
                        source_volume->label_);

  painter_->finish_progress();

  UndoHandle undo =
    new UndoReplaceLayer(painter_, "Undo Connected Component",
                         0, painter_->current_volume_, 0);
  painter_->push_undo(undo);

  painter_->extract_all_window_slices();
  painter_->redraw_all();

  painter_->hide_tool_panel();
}


}
