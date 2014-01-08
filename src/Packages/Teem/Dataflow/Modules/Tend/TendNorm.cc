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

//    File   : TendNorm.cc
//    Author : Martin Cole
//    Date   : Mon Sep  8 09:46:49 2003

#include <Dataflow/Network/Module.h>

#include <Dataflow/GuiInterface/GuiVar.h>
#include <Dataflow/Network/Ports/NrrdPort.h>
#include <teem/ten.h>

#include <sstream>
#include <iostream>
#include <stdio.h>

namespace SCITeem {

using namespace SCIRun;

class TendNorm : public Module {
  public:
    TendNorm(SCIRun::GuiContext *ctx);
    virtual ~TendNorm() {}
    virtual void execute();

  private:
    GuiDouble       major_weight_;
    GuiDouble       medium_weight_;
    GuiDouble       minor_weight_;
    GuiDouble       amount_;
    GuiDouble       target_;
};

DECLARE_MAKER(TendNorm)

TendNorm::TendNorm(SCIRun::GuiContext *ctx) : 
  Module("TendNorm", ctx, Filter, "Tend", "Teem"), 
  major_weight_(get_ctx()->subVar("major-weight"), 1.0),
  medium_weight_(get_ctx()->subVar("medium-weight"), 1.0),
  minor_weight_(get_ctx()->subVar("minor-weight"), 1.0),
  amount_(get_ctx()->subVar("amount"), 1.0),
  target_(get_ctx()->subVar("target"), 1.0)
{
}


void 
TendNorm::execute()
{
  update_state(NeedData);

  NrrdDataHandle nrrd_handle;
  get_input_handle("nin", nrrd_handle);

  if (inputs_changed_ || major_weight_.changed() || medium_weight_.changed() ||
      minor_weight_.changed() || amount_.changed() ||  target_.changed() ||
      !oport_cached("nout"))
  {
    // Force Teem to be locked befoer calling the Teem library
    NrrdGuard nrrd_guard;

    // Inform module that execution started
    update_state(Executing);
      
    Nrrd *nin = nrrd_handle->nrrd_;
    Nrrd *nout = nrrdNew();

    double weights[3];
    weights[0]=major_weight_.get();
    weights[1]=medium_weight_.get();
    weights[2]=minor_weight_.get();

    if (tenSizeNormalize(nout, nin, weights,
                         (float)amount_.get(), target_.get()))
    {
      char *err = biffGetDone(TEN);
      error(std::string("Error in TendNorm: ") + err);
      free(err);
      return;
    }

    NrrdDataHandle ntmp(new NrrdData(nout));

    send_output_handle("nout", ntmp);
  }
}


} // End namespace SCITeem
