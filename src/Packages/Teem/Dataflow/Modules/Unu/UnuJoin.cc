/*
   For more information, please see: http://software.sci.utah.edu

   The MIT License

   Copyright (c) 2009 Scientific Computing and Imaging Institute,
   University of Utah.

   
   Permission is hereby granted, free of charge, to any person obtaining a
   copy of this software and associated documentation files (the "Software"),
   to deal in the Software without restriction, including without limitation
   the rights to use, copy, modify, merge, publish, distribute, sublicense,
   and/or sell copies of the Software, and to permit persons to whom the
   Software is furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included
   in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
   OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
   THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
   DEALINGS IN THE SOFTWARE.
*/

//    File   : UnuJoin.cc
//    Author : Martin Cole
//    Date   : Wed Jan 15 10:59:57 2003

#include <Dataflow/Network/Module.h>
#include <Dataflow/GuiInterface/GuiVar.h>
#include <Core/Util/StringUtil.h>
#include <Core/Containers/Array1.h>
#include <Dataflow/Network/Ports/NrrdPort.h>

#include <iostream>
using std::endl;
#include <stdio.h>

namespace SCITeem {
using namespace SCIRun;

class UnuJoin : public Module {
  public:
    UnuJoin(SCIRun::GuiContext *ctx);
    virtual ~UnuJoin() {}
    virtual void execute();

  private:
    NrrdDataHandle    onrrd_handle_;  //! the cached output nrrd handle.
    std::vector<int>       in_generation_; //! all input generation nums.
    int               onrrd_type_;    //! target type for output nrrd.

    GuiInt            join_axis_;
    GuiInt            incr_dim_;
    GuiInt            dim_;
    int               old_axis_;
    int               old_incr_dim_;
};

} // End namespace SCITeem

using namespace SCITeem;

DECLARE_MAKER(UnuJoin)

UnuJoin::UnuJoin(SCIRun::GuiContext *ctx) : 
  Module("UnuJoin", ctx, Filter, "UnuAtoM", "Teem"), 
  onrrd_handle_(0),
  in_generation_(0),
  onrrd_type_(nrrdTypeLast),
  join_axis_(get_ctx()->subVar("join-axis"), 0),
  incr_dim_(get_ctx()->subVar("incr-dim"), 0),
  dim_(get_ctx()->subVar("dim"), 0),
  old_axis_(0),
  old_incr_dim_(0)
{
}

void 
UnuJoin::execute()
{
  std::vector<NrrdDataHandle> nrrd_handles_temp;
  std::vector<NrrdDataHandle> nrrd_handles;
  
  get_dynamic_input_handles("Nrrds",nrrd_handles_temp,true);
  if (nrrd_handles_temp.size() == 0) { return; }

  size_t i = 0;
  bool do_join = false;
  unsigned int max_dim = 0;

  for (i=0; i < nrrd_handles_temp.size(); i++)
  {
    // Although we should not have NrrdData objects with no nrrd inside
    // These appear to occur. Hence ignore those.
    if (nrrd_handles_temp[i]->nrrd_ == 0) continue;
    NrrdDataHandle nrrd_handle = nrrd_handles_temp[i];
    
    // check to see if we need to do the join or can output the cached onrrd.
    if (in_generation_.size() <= i)
    {
      // this is a new input, never been joined.
      do_join = true;
      in_generation_.push_back(nrrd_handle->generation);
      onrrd_type_ = nrrdTypeLast;
    }
    else if (in_generation_[i] != nrrd_handle->generation)
    {
      // different input than last execution
      do_join = true;
      in_generation_[i] = nrrd_handle->generation;
      onrrd_type_ = nrrdTypeLast;
    }

    // the output nrrd must be of one type, so find the type that accomodates
    // all of the nrrds we have as input.
    if (onrrd_type_ == nrrdTypeLast)
    {
      // first time this value is set
      onrrd_type_ = nrrd_handle->nrrd_->type;
    }
    if ((onrrd_type_ != nrrd_handle->nrrd_->type) && 
        (onrrd_type_ != nrrdTypeDouble))
    {
      //! promote to the biggest type
      if (nrrdTypeSize[nrrd_handle->nrrd_->type] > nrrdTypeSize[onrrd_type_])
      {
        onrrd_type_ = nrrd_handle->nrrd_->type;
      }
    }

    if (nrrd_handle->nrrd_->dim > max_dim)
    { 
      max_dim = nrrd_handle->nrrd_->dim;
    }
    nrrd_handles.push_back(nrrd_handle);
  }

  dim_.reset();
  if (max_dim != (unsigned int) dim_.get())
  {
    dim_.set(max_dim);
    dim_.reset();
  }

  // re-join if old axis is different from new
  // axis or incr_dim has changed
  if (old_axis_ != join_axis_.get())
  {
    do_join = true;
    old_axis_ = join_axis_.get();
  }
  if (old_incr_dim_ != incr_dim_.get())
  {
    do_join = true;
    old_incr_dim_ = incr_dim_.get();
  }

  std::vector<Nrrd*> arr(nrrd_handles.size());
  
  if (do_join || !oport_cached("JoinedNrrd"))
  {
    NrrdDataHandle onrrdH;
    
    int i = 0;
    std::string new_label("");
    std::vector<NrrdDataHandle>::iterator iter = nrrd_handles.begin();
    while(iter != nrrd_handles.end())
    {
      NrrdDataHandle nh = *iter;
      ++iter;

      NrrdData* cur_nrrd = nh.get_rep();
      // Does it need conversion to the bigger type?
      if (cur_nrrd->nrrd_->type != onrrd_type_)
      {
        Nrrd* new_nrrd = nrrdNew();
        if (nrrdConvert(new_nrrd, cur_nrrd->nrrd_, onrrd_type_))
        {
          char *err = biffGetDone(NRRD);
          error(std::string("Conversion Error: ") +  err);
          free(err);
          return;
        }
        arr[i] = new_nrrd;
      }
      else
      {
        arr[i] = cur_nrrd->nrrd_;
      }
      ++i;
    }

    join_axis_.reset();
    incr_dim_.reset();

    onrrdH = new NrrdData();
    
    if (nrrdJoin(onrrdH->nrrd_, &arr[0], nrrd_handles.size(),
		 join_axis_.get(), incr_dim_.get()))
    {
      char *err = biffGetDone(NRRD);
      error(std::string("Join Error: ") +  err);
      free(err);
      return;
    }
    send_output_handle("JoinedNrrd", onrrdH, true);

  }

}

