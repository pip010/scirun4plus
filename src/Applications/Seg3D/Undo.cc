//
//  For more information, please see: http://software.sci.utah.edu
//
//  The MIT License
//
//  Copyright (c) 2009 Scientific Computing and Imaging Institute,
//  University of Utah.
//
//  License for the specific language governing rights and limitations under
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
//    File   : Undo.cc
//    Author : Michael Callahan
//    Date   : June 2008


#include <Applications/Seg3D/Undo.h>
#include <Applications/Seg3D/Painter.h>
#include <Core/Events/EventManager.h>
#include <Core/Events/LoadVolumeEvent.h>

namespace SCIRun {

#define UNDO_MEM_TINY 0
#define UNDO_MEM_SLICE 1
#define UNDO_MEM_VOLUME 2

#define UNDO_SLICE_VOLUME_RATIO 50
#define UNDO_SIZE_IN_VOLUMES 2


Undo::Undo(Painter *painter, const string &label)
  : lock("Undo handle lock"),
    ref_cnt(0),
    painter_(painter),
    label_(label)
{
}


Undo::~Undo()
{
}


UndoManager::UndoManager()
{
}


UndoManager::~UndoManager()
{
}


void
UndoManager::push_undo(const UndoHandle &obj)
{
  // Push this item onto the undo stack.
  undo_stack_.push_front(obj);

  // Keep only the items that we have memory for.
  list<UndoHandle>::iterator itr = undo_stack_.begin();
  int volcount = 0;
  int slcount = 0;
  while (itr != undo_stack_.end())
  {
    const int me = (*itr)->mem_est();
    if (me == UNDO_MEM_VOLUME)
    {
      volcount++;
    }
    else if (me == UNDO_MEM_SLICE)
    {
      slcount++;
      if (slcount > UNDO_SLICE_VOLUME_RATIO)
      {
        volcount++;
        slcount -= UNDO_SLICE_VOLUME_RATIO;
      }
    }
    
    if (volcount > UNDO_SIZE_IN_VOLUMES)
    {
      break;
    }
    ++itr;
  }

  if (itr != undo_stack_.end())
  {
    undo_stack_.erase(itr, undo_stack_.end());
  }
  
  update_gui();
}


void
UndoManager::update_gui()
{
  update_undo_state_t *uus = new update_undo_state_t;
  
  if (undo_stack_.empty())
  {
    uus->enable_ = false;
  }
  else
  {
    uus->enable_ = true;
    uus->label_ = undo_stack_.front()->get_label();
  }

  wxCommandEvent event(wxEVT_COMMAND_UPDATE_UNDO_STATE, wxID_ANY);
  event.SetClientData((void *)uus);
  wxPostEvent(SCIRun::Painter::global_seg3dframe_pointer_, event);
}


void
UndoManager::undo_last()
{
  if (!undo_stack_.empty())
  {
    undo_stack_.front()->undo();
    undo_stack_.pop_front();
  }
  update_gui();
}


void
UndoManager::clear_undo()
{
  undo_stack_.clear();
  update_gui();
}


UndoReplaceLayer::UndoReplaceLayer(Painter *painter,
                                   const string &label,
                                   const NrrdVolumeHandle &oldvol,
                                   const NrrdVolumeHandle &newvol,
                                   size_t loc)
  : Undo(painter, label),
    oldvol_(oldvol),
    newvol_(newvol),
    loc_(loc)
{
}


UndoReplaceLayer::~UndoReplaceLayer()
{
}

void
UndoReplaceLayer::undo()
{
  if (newvol_.get_rep())
  {
    // Delete the new volume
    painter_->remove_volume(newvol_, false);
    newvol_->unparent();
  }
  if (oldvol_.get_rep())
  {
    // Reinsert the old volume into painter_->volumes at loc.
    painter_->insert_volume(oldvol_, loc_);
  }

  painter_->rebuild_layer_buttons();
  painter_->extract_all_window_slices();
  painter_->redraw_all();
}


int
UndoReplaceLayer::mem_est()
{
  // Big if there is a replacement, small if not.
  if (oldvol_.get_rep()) return UNDO_MEM_VOLUME;
  return UNDO_MEM_TINY;
}


UndoLabelInvertFilter::UndoLabelInvertFilter(Painter *painter,
                                             const NrrdVolumeHandle &v)
  : Undo(painter, "Undo Layer Invert"),
    volume_(v)
{
}

UndoLabelInvertFilter::~UndoLabelInvertFilter()
{
}

void
UndoLabelInvertFilter::undo()
{
  VolumeOps::bit_invert(volume_->nrrd_handle_,
                        volume_->label_,
                        volume_->nrrd_handle_,
                        volume_->label_);

  painter_->extract_all_window_slices();
  painter_->redraw_all();
}


int
UndoLabelInvertFilter::mem_est()
{
  return UNDO_MEM_TINY;
}



UndoReplaceSlice::UndoReplaceSlice(Painter *painter,
                                   const string &label,
                                   NrrdVolumeHandle &vol,
                                   NrrdDataHandle &slice,
                                   int axis, int coord)
  : Undo(painter, label),
    volume_(vol),
    slice_(slice),
    axis_(axis),
    coord_(coord)
{
}


UndoReplaceSlice::~UndoReplaceSlice()
{
}
  

void
UndoReplaceSlice::undo()
{
  volume_->lock.lock();
  if (volume_->nrrd_handle_->nrrd_->content)
  {
    volume_->nrrd_handle_->nrrd_->content[0] = 0;
  }

  if (nrrdSplice(volume_->nrrd_handle_->nrrd_,
                 volume_->nrrd_handle_->nrrd_,
                 slice_->nrrd_, axis_, coord_))
  {
    char *err = biffGetDone(NRRD);
    cerr << string("Error on line #") 
         << to_string(__LINE__)
         << string(" executing nrrd command: nrrdSplice \n")
         << string("Message: ") 
         << err
         << std::endl;

    free(err);
  }
  volume_->set_dirty();
  volume_->lock.unlock();

  painter_->extract_all_window_slices();
  painter_->redraw_all();
}


int
UndoReplaceSlice::mem_est()
{
  return UNDO_MEM_SLICE;
}


UndoReplaceVolume::UndoReplaceVolume(Painter *painter,
                                     const string &label,
                                     NrrdVolumeHandle &vol,
                                     NrrdDataHandle &data)
  : Undo(painter, label),
    volume_(vol),
    data_(data)
{
}


UndoReplaceVolume::~UndoReplaceVolume()
{
}
  

void
UndoReplaceVolume::undo()
{
  volume_->lock.lock();
  // TODO: We just replace the layer a the moment.  Should replace it
  // all so as to be more general (non label volumes).
  VolumeOps::bit_copy(volume_->nrrd_handle_, volume_->label_,
                      data_, volume_->label_);
  volume_->lock.unlock();

  painter_->extract_all_window_slices();
  painter_->redraw_all();
}


int
UndoReplaceVolume::mem_est()
{
  return UNDO_MEM_VOLUME;
}


UndoFlipAxis::UndoFlipAxis(Painter *painter,
                           const string &label,
                           NrrdVolumeHandle &vol,
                           int axis)
  : Undo(painter, label),
    volume_(vol),
    axis_(axis)
{
}


UndoFlipAxis::~UndoFlipAxis()
{
}
  

void
UndoFlipAxis::undo()
{
  NrrdDataHandle nout = new NrrdData();
  nrrdFlip(nout->nrrd_, volume_->nrrd_handle_->nrrd_, axis_);

  nrrdKeyValueCopy(nout->nrrd_, volume_->nrrd_handle_->nrrd_);

  volume_->nrrd_handle_ = nout;
  volume_->set_dirty();

  painter_->extract_all_window_slices();
  painter_->redraw_all();  
  
  if (volume_ == painter_->current_vrender_target_ &&
      !painter_->current_vrender_target_deferred_)
  {
    LoadVolumeEvent *lve =
      new LoadVolumeEvent(painter_->current_vrender_target_->nrrd_handle_,
                          "", false);
    EventManager::add_event(lve);
  }

}


int
UndoFlipAxis::mem_est()
{
  return UNDO_MEM_TINY;
}


UndoInvert::UndoInvert(Painter *painter,
		      const string &label,
		      NrrdVolumeHandle &vol)
  : Undo(painter, label),
    volume_(vol)
{
}


UndoInvert::~UndoInvert()
{
}
  

void
UndoInvert::undo()
{
  Nrrd *nin1 = volume_->nrrd_handle_->nrrd_;

  NrrdIter *in1 = nrrdIterNew();
  NrrdIter *in2 = nrrdIterNew();

  nrrdIterSetOwnNrrd(in1, nin1);
  nrrdIterSetValue(in2, volume_->data_max_ + volume_->data_min_);

  NrrdDataHandle nout_data_handle = new NrrdData();
  nrrdArithIterBinaryOp(nout_data_handle->nrrd_,  nrrdBinaryOpSubtract, in2, in1);

  nrrdKeyValueCopy(nout_data_handle->nrrd_, volume_->nrrd_handle_->nrrd_);

  volume_->nrrd_handle_ = nout_data_handle;
  volume_->set_dirty();

  painter_->extract_all_window_slices();
  painter_->redraw_all();  
  
  if (volume_ == painter_->current_vrender_target_ &&
      !painter_->current_vrender_target_deferred_)
  {
    LoadVolumeEvent *lve =
      new LoadVolumeEvent(painter_->current_vrender_target_->nrrd_handle_,
                          "", false);
    EventManager::add_event(lve);
  }

}


int
UndoInvert::mem_est()
{
  return UNDO_MEM_TINY;
}


} // namespace SCIRun

