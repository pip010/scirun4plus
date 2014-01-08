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

//    File   : UnuCmedian.cc
//    Author : Martin Cole
//    Date   : Mon Aug 25 10:13:16 2003

#include <Dataflow/Network/Module.h>

#include <Dataflow/GuiInterface/GuiVar.h>
#include <Dataflow/Network/Ports/NrrdPort.h>

#include <sstream>
#include <iostream>
using std::endl;
#include <stdio.h>

namespace SCITeem {

using namespace SCIRun;

class UnuCmedian : public Module {
  public:
    UnuCmedian(SCIRun::GuiContext *ctx);
    virtual ~UnuCmedian() {}
    virtual void execute();

  private:
    Nrrd* do_filter(Nrrd *nin);
    bool is_scalar(const std::string& s) const;

    GuiInt          mode_;
    GuiInt          radius_;
    GuiDouble       weight_;
    GuiInt          bins_;
    GuiInt          pad_;
};

DECLARE_MAKER(UnuCmedian)

UnuCmedian::UnuCmedian(SCIRun::GuiContext *ctx) : 
  Module("UnuCmedian", ctx, Filter, "UnuAtoM", "Teem"), 
  mode_(get_ctx()->subVar("mode")),
  radius_(get_ctx()->subVar("radius")),
  weight_(get_ctx()->subVar("weight")),
  bins_(get_ctx()->subVar("bins")),
  pad_(get_ctx()->subVar("pad"))
{
}


Nrrd*
UnuCmedian::do_filter(Nrrd *nin)
{
  reset_vars();
  Nrrd *ntmp, *nout;
  nout = nrrdNew();
  ntmp = nin;
//   if (pad_.get()) {
//     ntmp = nrrdNew();
//     if (nrrdSimplePad(ntmp, nin, radius_.get(), nrrdBoundaryBleed)) {
//       char *err = biffGetDone(NRRD);
//       error(std::string("Error padding: ") + err);
//       free(err);
//       return 0;
//     }
//   }
//   else {
//     ntmp = nin;
//   }

  if (nrrdCheapMedian(nout, ntmp, pad_.get(), mode_.get(), radius_.get(), 
		      weight_.get(), bins_.get())) 
  {
    char *err = biffGetDone(NRRD);
    error(std::string("Error doing cheap median: ") + err);
    free(err);
    return 0;
  }
  
//   if (pad_.get()) {
//     if (nrrdSimpleCrop(ntmp, nout, radius_.get())) {
//       char *err = biffGetDone(NRRD);
//       error(std::string("Error cropping: ") + err);
//       free(err);
//       return 0;
//     }
//     nrrdNuke(nout);
//     return ntmp;
//   }
  return nout;
}

bool
UnuCmedian::is_scalar(const std::string& s) const
{
  const std::string scalar(":Scalar");
  return (s.find(scalar) < s.size());
}


void
UnuCmedian::execute()
{
  update_state(NeedData);

  NrrdDataHandle nrrd_handle;
  get_input_handle("nin", nrrd_handle,true);

  if (inputs_changed_ || mode_.changed() || radius_.changed() ||
      weight_.changed() || bins_.changed() || pad_.changed() ||
      oport_cached("nout"))
  {
    // Force Teem to be locked befoer calling the Teem library
    NrrdGuard nrrd_guard;

    // Inform module that execution started
    update_state(Executing);  
   
    Nrrd *nin = nrrd_handle->nrrd_;
    NrrdDataHandle nsend(0);

    // loop over the tuple axis, and do the median filtering for 
    // each scalar set independently, copy non scalar sets unchanged.
    //vector<std::string> elems;
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
      std::vector<Nrrd*> out;
      for (unsigned int i = 0; i < nrrd_handle->nrrd_->axis[0].size; i++) 
      { 
        Nrrd *sliced = nrrdNew();
        if (nrrdSlice(sliced, nin, 0, i)) 
        {
          char *err = biffGetDone(NRRD);
          error(std::string("Trouble with slice: ") + err);
          free(err);
        }
        
        Nrrd *nout_filtered;
        nout_filtered = do_filter(sliced);
        if (!nout_filtered) 
        {
          error("Error filtering, returning");
          return;
        }
        out.push_back(nout_filtered);
      }
      // Join the filtered nrrds along the first axis
      NrrdData *nrrd_joined = new NrrdData;
      if (nrrdJoin(nrrd_joined->nrrd_, &out[0], out.size(), 0, 1)) {
        char *err = biffGetDone(NRRD);
        error(std::string("Join Error: ") +  err);
        free(err);
        return;
      }
      NrrdDataHandle ntmp(nrrd_joined);
      send_output_handle("nout", ntmp);
    }
    else
    {
      Nrrd *nout_filtered;
      nout_filtered = do_filter(nrrd_handle->nrrd_);
      if (!nout_filtered)
      {
        error("Error filtering, returning");
        return;
      }
      NrrdDataHandle ntmp(new NrrdData(nout_filtered));
      send_output_handle("nout", ntmp);
    }
  }
}


} // End namespace SCITeem
